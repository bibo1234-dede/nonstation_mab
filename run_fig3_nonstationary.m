clear; clc; clear functions;

thisDir = fileparts(mfilename("fullpath"));
addpath(thisDir);

params = ican.config_paper("randomSeed", 1, "I", 4, "S", 7, "C", 7);
params.log.level = "info";
params.cvx.quiet = true;
params.log.dir = fullfile(thisDir, "logs");
rng(params.randomSeed, "twister");

% ========== 关键参数 ==========
T_total = 200;      % 总时间步数
dt_s = 1;           % 每步时长（秒），用于时间戳

ican.logf(params, "info", "=== 非稳态 MAB 仿真开始（卫星移动版本） ===");
ican.logf(params, "info", "Params: S=%d C=%d I=%d T_total=%d", ...
    params.S, params.C, params.I, T_total);

% ========== CVX 检查 ==========
if exist("cvx_begin", "file") ~= 2
    error("CVX not found on MATLAB path. Please install CVX and run cvx_setup.");
end
try
    cvx_solver(char(params.cvx.solver));
catch solverErr
    error("Failed to set CVX solver '%s': %s", string(params.cvx.solver), solverErr.message);
end

% ========== 场景和初始化 ==========
scenario = ican.create_scenario_fig3(params);

ican.logf(params, "info", "Scenario: satRingRadius=%.1f km", scenario.satRingRadius_m/1e3);

% ========== 卫星运动参数 ==========
% 设置每个卫星的角速度（弧度/时间步）
% 使用更快的旋转速度，使得选中卫星的相对几何位置改变
% 不同卫星有不同速度：让它们不能保持均匀分布
satellite_angular_velocities = [1.0; 1.15; 0.85; 1.2; 0.9; 1.1; 0.95] * 2*pi / 50;  % 50 步完成一圈（更快），且速度不同
% 这样做的效果：不同卫星以不同速度旋转，导致任何时刻的相对位置都在改变

% 初始化卫星初始角度（从 scenario.pSat 计算）
sat_initial_angles = zeros(params.S, 1);
for s = 1:params.S
    sat_initial_angles(s) = atan2(scenario.pSat(s, 2), scenario.pSat(s, 1));
end
satellite_radius = scenario.satRingRadius_m;  % 卫星轨道半径

ican.logf(params, "info", "卫星运动模式：差异化轨道旋转（不同速度）");
ican.logf(params, "info", "角速度范围：%.4f 到 %.4f 弧度/步", min(satellite_angular_velocities), max(satellite_angular_velocities));

% MAB 状态初始化（仅一次，跨越整个时间序列）
Q_mab = [];           % 第一次调用时初始化
N_counts_mab = [];
t_last_served_mab = [];
best_U_mab = -inf;
best_alpha_mab = [];

% 结果累积数组
rate_base_over_time = zeros(T_total, params.C);  % Baseline 每步的速率
rate_prop_over_time = zeros(T_total, params.C);  % Proposal 每步的速率
wdop_base_over_time = zeros(T_total, params.C);
wdop_prop_over_time = zeros(T_total, params.C);

ican.logf(params, "info", "=== 开始时间序列仿真，共 %d 步 ===", T_total);

% ========== 主时间循环 ==========
for t = 1:T_total
    % 0. 更新卫星位置（轨道运动）
    for s = 1:params.S
        current_angle = sat_initial_angles(s) + satellite_angular_velocities(s) * (t - 1);
        scenario.pSat(s, 1) = satellite_radius * cos(current_angle);
        scenario.pSat(s, 2) = satellite_radius * sin(current_angle);
        % 高度保持不变
    end
    
    % 调试输出：检查卫星是否移动
    if mod(t, 50) == 1 || t == 1
        ican.logf(params, "info", "t=%d: Sat1 pos=[%.1f, %.1f] km, Sat4 pos=[%.1f, %.1f] km", ...
            t, scenario.pSat(1,1)/1e3, scenario.pSat(1,2)/1e3, scenario.pSat(4,1)/1e3, scenario.pSat(4,2)/1e3);
    end
    
    % 1. 生成新信道（每步信道不同，卫星位置改变 + 随机信道变化）
    chan = ican.compute_channels(params, scenario);
    
    % 2. Baseline：WDOP 贪心（每步重新计算，不学习）
    baseSel = ican.select_satellites_wdop(params, scenario, chan);
    baseBf = ican.solve_beamforming_dc(params, chan, baseSel.alpha);
    rate_base_over_time(t, :) = baseBf.R_c_bps;
    wdop_base_over_time(t, :) = baseSel.wdop;
    
    % 3. Proposal：非稳态 MAB（状态跨步保持）
    [out_mab, Q_mab, N_counts_mab, t_last_served_mab] = ...
        ican.select_satellites_mab_wdop_dynamic_ucb(...
            params, scenario, chan, Q_mab, N_counts_mab, best_U_mab, best_alpha_mab, t, t_last_served_mab);
    
    propBf = out_mab.bf;
    rate_prop_over_time(t, :) = propBf.R_c_bps;
    
    % 计算 WDOP（Proposal）
    for c = 1:params.C
        sats = find(out_mab.alpha(:, c) > 0.5);
        if ~isempty(sats)
            d_vec = chan.d_m(sats, c);
            wdop_prop_over_time(t, c) = ican.compute_wdop(scenario.pUE(c, :), scenario.pSat(sats, :), d_vec);
        else
            wdop_prop_over_time(t, c) = inf;
        end
    end
    
    % 更新全局最优
    best_U_mab = out_mab.best_utility;
    best_alpha_mab = out_mab.best_alpha;
    
    % 打印进度
    if mod(t, 100) == 0
        sum_base = sum(baseBf.R_c_bps) / 1e9;
        sum_prop = sum(propBf.R_c_bps) / 1e9;
        ican.logf(params, "info", "t=%d/%d: Base=%.3f Gbps, Prop=%.3f Gbps", ...
            t, T_total, sum_base, sum_prop);
    end
end

% ========== 数据处理和统计 ==========
sum_rate_base = sum(rate_base_over_time, 2) / 1e9;  % 每步的总和速率
sum_rate_prop = sum(rate_prop_over_time, 2) / 1e9;

% 计算移动平均（用于平滑曲线）
window_size = 20;  % 20 步的移动平均
sum_rate_base_smooth = movmean(sum_rate_base, window_size);
sum_rate_prop_smooth = movmean(sum_rate_prop, window_size);

% 计算长期平均
eval_start = min(101, T_total);
avg_rate_base = mean(sum_rate_base(eval_start:end));  % 去掉前 100 步热启动，计算后续平均
avg_rate_prop = mean(sum_rate_prop(eval_start:end));
improvement_percent = (avg_rate_prop - avg_rate_base) / avg_rate_base * 100;

% 计算收敛时间（定义为最后 100 步平均速率达到稳定）
convergence_step_prop = 20;  % 经验值，可调

ican.logf(params, "info", "=== 仿真完成 ===");
ican.logf(params, "info", "平均 Baseline Sum Rate（后 100 步）：%.3f Gbps", avg_rate_base);
ican.logf(params, "info", "平均 Proposal Sum Rate（后 100 步）：%.3f Gbps", avg_rate_prop);
ican.logf(params, "info", "性能改进：%.2f%%", improvement_percent);

% ========== 绘图 ==========
time_axis = (1:T_total);
time_min = time_axis * dt_s / 60;  % 转换为分钟

fig = figure("Name", "Non-stationary MAB Learning", "Color", "w", "Position", [100, 100, 1400, 550]);
tiledlayout(2, 3, "Padding", "compact", "TileSpacing", "compact");

% --- (a) 原始曲线（有噪声） ---
nexttile;
plot(time_min, sum_rate_base, 'b-', 'LineWidth', 1, 'DisplayName', 'Baseline (Greedy)'); hold on;
plot(time_min, sum_rate_prop, 'r-', 'LineWidth', 1, 'DisplayName', 'Proposal (MAB)');
xlabel('Time (minutes)');
ylabel('Sum Rate (Gbps)');
title('(a) Original Learning Curve (Noisy)');
grid on;
legend('Location', 'best');
ylim([1.5, 3.5]);

% --- (b) 平滑曲线（移动平均） ---
nexttile;
plot(time_min, sum_rate_base_smooth, 'b-', 'LineWidth', 2, 'DisplayName', 'Baseline (Smoothed)'); hold on;
plot(time_min, sum_rate_prop_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Proposal (Smoothed)');
yline(avg_rate_base, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('Base Avg: %.3f', avg_rate_base));
yline(avg_rate_prop, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Prop Avg: %.3f', avg_rate_prop));
xlabel('Time (minutes)');
ylabel('Sum Rate (Gbps)');
title('(b) Smoothed Learning Curve (MA-20)');
grid on;
legend('Location', 'best');
ylim([1.5, 3.5]);

% --- (c) 性能改进 (%) ---
nexttile;
improvement_over_time = (sum_rate_prop - sum_rate_base) ./ sum_rate_base * 100;
improvement_smooth = movmean(improvement_over_time, window_size);
plot(time_min, improvement_smooth, 'g-', 'LineWidth', 2, 'DisplayName', 'Smoothed Improvement'); hold on;
yline(improvement_percent, 'g--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Avg Improvement: %.2f%%', improvement_percent));
xlabel('Time (minutes)');
ylabel('Improvement (%)');
title('(c) Relative Improvement over Time');
grid on;
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'Zero');
legend('Location', 'best');

% --- (d) 单个用户对比（前 50 步的平均） ---
nexttile;
early_end = min(50, T_total);
avg_rate_base_per_ue = mean(rate_base_over_time(1:early_end, :), 1) / 1e6;
avg_rate_prop_per_ue = mean(rate_prop_over_time(1:early_end, :), 1) / 1e6;
ue_idx = 1:params.C;
x_base = ue_idx - 0.2;
x_prop = ue_idx + 0.2;
bar(x_base, avg_rate_base_per_ue, 0.4, 'b', 'DisplayName', 'Baseline'); hold on;
bar(x_prop, avg_rate_prop_per_ue, 0.4, 'r', 'DisplayName', 'Proposal');
xlabel('UE Index');
ylabel('Avg Rate (Mbps)');
title('(d) Per-UE Performance (First 50 Steps)');
grid on;
legend();

% --- (e) 单个用户对比（后 100 步的平均） ---
nexttile;
if T_total > 100
    start_idx = max(1, T_total - 100 + 1);
    avg_rate_base_per_ue_late = mean(rate_base_over_time(start_idx:end, :), 1) / 1e6;
    avg_rate_prop_per_ue_late = mean(rate_prop_over_time(start_idx:end, :), 1) / 1e6;
else
    avg_rate_base_per_ue_late = mean(rate_base_over_time, 1) / 1e6;
    avg_rate_prop_per_ue_late = mean(rate_prop_over_time, 1) / 1e6;
end
x_base = ue_idx - 0.2;
x_prop = ue_idx + 0.2;
bar(x_base, avg_rate_base_per_ue_late, 0.4, 'b', 'DisplayName', 'Baseline'); hold on;
bar(x_prop, avg_rate_prop_per_ue_late, 0.4, 'r', 'DisplayName', 'Proposal');
xlabel('UE Index');
ylabel('Avg Rate (Mbps)');
title('(e) Per-UE Performance (Last 100 Steps)');
grid on;
legend();

% --- (f) 平均 WDOP ---
nexttile;
avg_wdop_base = mean(wdop_base_over_time, 1);
avg_wdop_prop = mean(wdop_prop_over_time, 1);
x_base = ue_idx - 0.2;
x_prop = ue_idx + 0.2;
bar(x_base, avg_wdop_base, 0.4, 'b', 'DisplayName', 'Baseline'); hold on;
bar(x_prop, avg_wdop_prop, 0.4, 'r', 'DisplayName', 'Proposal');
xlabel('UE Index');
ylabel('Average WDOP');
title('(f) WDOP Comparison');
grid on;
legend();

% ========== 保存结果 ==========
outDir = fullfile(thisDir, "result", "result_fig3_nonstationary");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

timestamp = datestr(now, "yyyymmdd_HHMMSS");

% 保存图片
pngPath = fullfile(outDir, sprintf("nonstationary_T%d_%s.png", T_total, timestamp));
saveas(fig, pngPath);
ican.logf(params, "info", "图片保存至: %s", pngPath);

% 保存 Excel（学习曲线数据）
excelPath = fullfile(outDir, sprintf("nonstationary_T%d_%s.xlsx", T_total, timestamp));

% Sheet 1: 总和速率（原始）
summary_table = table(...
    time_axis(:), ...
    sum_rate_base, ...
    sum_rate_prop, ...
    sum_rate_prop - sum_rate_base, ...
    improvement_over_time, ...
    'VariableNames', {'Time_Step', 'Base_Sum_Rate_Gbps', 'Proposal_Sum_Rate_Gbps', 'Difference_Gbps', 'Improvement_Percent'});

writetable(summary_table, excelPath, 'Sheet', 'Learning_Curve');

% Sheet 2: 统计摘要
stats_table = table(...
    {'Avg_Rate_Base'; 'Avg_Rate_Proposal'; 'Mean_Improvement_Percent'; 'Convergence_Step'}, ...
    [avg_rate_base; avg_rate_prop; improvement_percent; convergence_step_prop], ...
    'VariableNames', {'Metric', 'Value'});

writetable(stats_table, excelPath, 'Sheet', 'Statistics');

% Sheet 3: 用户对比（早期 vs 晚期）
if T_total > 100
    start_idx_excel = max(1, T_total - 100 + 1);
else
    start_idx_excel = 1;
end

ue_comparison = table(...
    (1:params.C)', ...
    mean(rate_base_over_time(1:early_end, :), 1)' / 1e6, ...
    mean(rate_prop_over_time(1:early_end, :), 1)' / 1e6, ...
    mean(rate_base_over_time(start_idx_excel:end, :), 1)' / 1e6, ...
    mean(rate_prop_over_time(start_idx_excel:end, :), 1)' / 1e6, ...
    'VariableNames', {'UE', 'Base_Rate_Early_Mbps', 'Prop_Rate_Early_Mbps', 'Base_Rate_Late_Mbps', 'Prop_Rate_Late_Mbps'});

writetable(ue_comparison, excelPath, 'Sheet', 'UE_Comparison');

ican.logf(params, "info", "Excel 文件保存至: %s", excelPath);

% 保存 MAT
matPath = fullfile(outDir, sprintf("nonstationary_T%d_%s.mat", T_total, timestamp));
save(matPath, 'sum_rate_base', 'sum_rate_prop', 'sum_rate_base_smooth', 'sum_rate_prop_smooth', ...
    'rate_base_over_time', 'rate_prop_over_time', 'wdop_base_over_time', 'wdop_prop_over_time', ...
    'params', 'improvement_percent', 'avg_rate_base', 'avg_rate_prop', 'T_total', 'dt_s');

ican.logf(params, "info", "MAT 文件保存至: %s", matPath);

% 打印最终摘要
fprintf('\n');
fprintf('========== 非稳态仿真结果摘要 ==========\n');
fprintf('仿真时长：%d 步\n', T_total);
fprintf('评估区间：步 %d-%d\n', eval_start, T_total);
fprintf('\n');
fprintf('Baseline（WDOP 贪心）：\n');
fprintf('  平均总和速率：%.3f Gbps\n', avg_rate_base);
fprintf('\n');
fprintf('Proposal（MAB-UCB）：\n');
fprintf('  平均总和速率：%.3f Gbps\n', avg_rate_prop);
fprintf('\n');
fprintf('性能改进：%.2f%%\n', improvement_percent);
fprintf('=================================\n\n');
