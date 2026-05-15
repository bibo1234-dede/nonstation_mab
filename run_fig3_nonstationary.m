clear; clc; clear functions;
% 非稳态 MAB 主脚本：场景、分组、选星、波束成形与结果导出。

thisDir = fileparts(mfilename("fullpath"));
addpath(thisDir);

% 主流程：参数初始化 -> 场景构建 -> 用户分组 -> 时序仿真 -> 结果导出
params = ican.config_paper("randomSeed", 1, "I", 4, "S", 18, "C", 10);
params.grouping.numGroups = 3;
params.grouping.numLevels = 3;
params.grouping.method = "by_level";      % 选择分组方式："spectral" 或 "by_level"
params.useGrouping = false;               % 启用分组
params.log.level = "info";
params.cvx.quiet = true; % 默认静默，失败时自动启用详细日志重试
params.log.dir = fullfile(thisDir, "logs");
rng(params.randomSeed, "twister");

% 是否运行 Baseline（WDOP 贪心选星 + 波束赋形）
% 设为 false 时，仅运行 Proposal 路径。
runBaseline = true;

T_total = 200;
dt_s = 1;

fprintf("=== 非稳态 MAB 仿真开始 ===\n");
fprintf("参数：S=%d，C=%d，I=%d，T_total=%d\n", params.S, params.C, params.I, T_total);

if exist("cvx_begin", "file") ~= 2
    cvxRoot = "D:\Users\44574\Desktop\瑶的东西\卫星\cvx";
    cvxStartup = fullfile(cvxRoot, "cvx_startup.m");
    if exist(cvxStartup, "file") == 2
        run(cvxStartup);
    end
end
if exist("cvx_begin", "file") ~= 2
    error("未找到 CVX，请先安装 CVX 并执行 cvx_setup。");
end
if exist("cvx_clear", "file") == 2
    cvx_clear;
end

scenario = ican.create_scenario_fig3(params);
fprintf("场景：轨道半径 = %.1f km\n", scenario.R_orbit / 1e3);

% MAB 状态缓存：跨时间步保留，支持折扣更新
Q_mab = [];
N_counts_mab = [];
t_last_served_mab = [];
best_U_mab = -inf;
best_alpha_mab = [];

% 记录每个时刻、每个用户的速率与 WDOP
chan = ican.compute_channels(params, scenario);

% 应用等效信道增益（默认 10 dB，避免数值过大）
gain_amplitude = sqrt(params.effectiveGain_linear);   % 幅度 = sqrt(10^(dB/10))
chan.h = chan.h * gain_amplitude;

if params.log.printArmSpace
    print_full_arm_space(params, scenario, chan);
end

if params.useGrouping
    % 基于优先级的分组
    [user_groups, level_info] = ican.user_grouping_by_level(params, scenario, chan);
    numGroups = numel(user_groups.group_ids);
    fprintf("=== 用户分组完成：%d 组 ===\n", numGroups);
    for k = 1:numGroups
        fprintf("  第 %d 组：%s\n", k, mat2str(user_groups.group_ids{k}(:).'));
    end
else
    user_groups = default_user_groups_all_users(params.C);
    numGroups = 1;
    fprintf("=== 已关闭分组，采用全用户共享的单组占位结构 ===\n");
end

rate_base_over_time = nan(T_total, params.C);
rate_prop_over_time = zeros(T_total, params.C);
wdop_base_over_time = nan(T_total, params.C);
wdop_prop_over_time = zeros(T_total, params.C);
group_rate_base_over_time = nan(T_total, numGroups);
group_rate_prop_over_time = zeros(T_total, numGroups);

fprintf("=== 开始时间序列仿真，共 %d 步 ===\n", T_total);

for t = 1:T_total
    if t == 1 || mod(t, 10) == 0 || t == T_total
        fprintf("[进度] 正在处理 t=%d/%d ...\n", t, T_total);
        drawnow;
    end

    scenario = update_satellite_positions(scenario, t, dt_s);
    chan = ican.compute_channels(params, scenario);

    % 应用等效信道增益（默认 10 dB，避免数值过大）
    gain_amplitude = sqrt(params.effectiveGain_linear);   % 幅度 = sqrt(10^(dB/10))
    chan.h = chan.h * gain_amplitude;

    % 第 1 步：在初始时刻更新用户分组（仅 by_level 模式）
    if params.useGrouping && t == 1
        [user_groups, level_info] = ican.user_grouping_by_level(params, scenario, chan);
        numGroups = numel(user_groups.group_ids);
        group_rate_base_over_time = nan(T_total, numGroups);
        group_rate_prop_over_time = zeros(T_total, numGroups);
        fprintf("=== 用户分组完成：%d 组 ===\n", numGroups);
    end

    % Baseline：可选执行（默认关闭）
    if runBaseline
        baseSel = ican.select_satellites_wdop(params, scenario, chan);
        baseBf = ican.solve_beamforming_dc(params, chan, baseSel.alpha);
        rate_base_over_time(t, :) = baseBf.R_c_bps(:).';
        wdop_base_over_time(t, :) = baseSel.wdop(:).';
        group_rate_base_over_time(t, :) = compute_group_rates(rate_base_over_time(t, :), user_groups);
        % 打印 Baseline 详细信息（每个 UE 的卫星集合与速率、WDOP）
        for cc = 1:params.C
            sats_cc = find(baseSel.alpha(:, cc) > 0.5);
            fprintf('[t=%d][Base] UE%d sats=%s rate=%.3f Mbps WDOP=%.3f\n', ...
                t, cc, mat2str(sats_cc), baseBf.R_c_bps(cc)/1e6, baseSel.wdop(cc));
        end
        % 打印 Baseline 总速率
        total_base_Mbps = sum(baseBf.R_c_bps) / 1e6;
            fprintf('[t=%d][Base] TotalRate=%.3f Mbps\n', t, total_base_Mbps);
    end

    % Proposal：非稳态 MAB + WDOP 软约束 + 可选 Pareto 目标
    [out_mab, Q_mab, N_counts_mab, t_last_served_mab] = ican.select_satellites_mab_wdop_dynamic_ucb( ...
        params, scenario, chan, Q_mab, N_counts_mab, best_U_mab, best_alpha_mab, t, t_last_served_mab, user_groups);

    propBf = out_mab.bf;
    rate_prop_over_time(t, :) = propBf.R_c_bps(:).';
    group_rate_prop_over_time(t, :) = compute_group_rates(rate_prop_over_time(t, :), user_groups);

    % 打印 Proposal 详细信息（每个 UE 的卫星集合与速率、WDOP）
    for cc = 1:params.C
        sats_cc = find(out_mab.alpha(:, cc) > 0.5);
        arm_idx = out_mab.selected_arm_idx(cc);
        fprintf('[t=%d][Prop] UE%d arm#%d sats=%s rate=%.3f Mbps WDOP=%.3f\n', ...
            t, cc, arm_idx, mat2str(sats_cc), propBf.R_c_bps(cc)/1e6, out_mab.wdop_per_user(cc));
    end
    % 打印每颗卫星在 Proposal 下的总速率（Mbps）
    if isfield(propBf, 'satSumRate_bps')
        for s = 1:params.S
            fprintf('[t=%d][Prop] Sat%d totalRate=%.3f Mbps\n', t, s, propBf.satSumRate_bps(s)/1e6);
        end
    end
    % 打印 Proposal 总速率
    total_prop_Mbps = sum(propBf.R_c_bps) / 1e6;
    fprintf('[t=%d][Prop] TotalRate=%.3f Mbps\n', t, total_prop_Mbps);

    for c = 1:params.C
        sats = find(out_mab.alpha(:, c) > 0.5);
        if ~isempty(sats)
            d_vec = chan.d_m(sats, c);
            wdop_prop_over_time(t, c) = ican.compute_wdop(scenario.pUE(c, :), scenario.pSat(sats, :), d_vec);
        else
            wdop_prop_over_time(t, c) = inf;
        end
    end

    best_U_mab = out_mab.best_utility;
    best_alpha_mab = out_mab.best_alpha;

    if mod(t, 20) == 0 || t == 1 || t == T_total
        if runBaseline
            fprintf("[进度] t=%d/%d | 基线=%.3f Mbps | 方案=%.3f Mbps\n", ...
                t, T_total, sum(rate_base_over_time(t, :), "omitnan") / 1e6, sum(rate_prop_over_time(t, :)) / 1e6);
        else
            fprintf("[进度] t=%d/%d | 方案=%.3f Mbps\n", ...
                t, T_total, sum(rate_prop_over_time(t, :)) / 1e6);
        end
        drawnow;
    end
    fprintf("===> 时刻 t=%d 完成 <===\n", t);
end

sum_rate_base = sum(rate_base_over_time, 2, "omitnan") / 1e6;
sum_rate_prop = sum(rate_prop_over_time, 2) / 1e6;
window_size = 20;
sum_rate_base_smooth = movmean(sum_rate_base, window_size);
sum_rate_prop_smooth = movmean(sum_rate_prop, window_size);

avg_wdop_base_time = mean(wdop_base_over_time, 2, "omitnan");
avg_wdop_prop_time = mean(wdop_prop_over_time, 2, "omitnan");

eval_start = max(1, T_total - 99);
avg_rate_base = mean(sum_rate_base(eval_start:end));
avg_rate_prop = mean(sum_rate_prop(eval_start:end));
improvement_percent = (avg_rate_prop - avg_rate_base) / max(avg_rate_base, eps) * 100;

avg_wdop_base = mean(wdop_base_over_time(eval_start:end, :), "all", "omitnan");
avg_wdop_prop = mean(wdop_prop_over_time(eval_start:end, :), "all", "omitnan");

if ~runBaseline
    sum_rate_base(:) = NaN;
    sum_rate_base_smooth(:) = NaN;
    avg_rate_base = NaN;
    avg_wdop_base = NaN;
    improvement_percent = NaN;
end

fprintf("=== 仿真完成 ===\n");
if runBaseline
    fprintf("平均 Baseline Sum Rate（后100步）：%.3f Mbps\n", avg_rate_base);
else
    fprintf("Baseline 已关闭（未执行 WDOP 贪心选星）。\n");
end
fprintf("平均 Proposal Sum Rate（后100步）：%.3f Mbps\n", avg_rate_prop);
if runBaseline
    fprintf("性能改进：%.2f%%\n", improvement_percent);
else
    fprintf("性能改进：N/A（Baseline 关闭）\n");
end

color_urllc = [0.85, 0.20, 0.20];
color_embb = [0.20, 0.55, 0.90];

group_colors = lines(numGroups);
if numGroups >= 1
    group_colors = group_colors(1:numGroups, :);
end

ue_colors = zeros(params.C, 3);
for c = 1:params.C
    if isfield(user_groups, "group_assignment") && numel(user_groups.group_assignment) >= c
        g = user_groups.group_assignment(c);
        if g >= 1 && g <= numGroups
            ue_colors(c, :) = group_colors(g, :);
        end
    end
    if all(ue_colors(c, :) == 0)
        ue_colors(c, :) = [0.5, 0.5, 0.5];
    end
end

fig = figure("Name", "非稳态 MAB 与谱聚类分组", "Color", "w", "Position", [80, 80, 1600, 900]);
tiledlayout(2, 4, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot((1:T_total) * dt_s / 60, sum_rate_base, 'b-', 'LineWidth', 1.2, 'DisplayName', '基线'); hold on;
plot((1:T_total) * dt_s / 60, sum_rate_prop, 'r-', 'LineWidth', 1.2, 'DisplayName', '方案');
xlabel('时间（分钟）');
ylabel('总速率（Mbps）');
title('（a）原始学习曲线');
grid on; legend('Location', 'best');

nexttile;
plot((1:T_total) * dt_s / 60, sum_rate_base_smooth, 'b-', 'LineWidth', 2, 'DisplayName', '基线平滑曲线'); hold on;
plot((1:T_total) * dt_s / 60, sum_rate_prop_smooth, 'r-', 'LineWidth', 2, 'DisplayName', '方案平滑曲线');
yline(avg_rate_base, 'b--', 'LineWidth', 1.2, 'DisplayName', sprintf('基线均值：%.3f', avg_rate_base));
yline(avg_rate_prop, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('方案均值：%.3f', avg_rate_prop));
xlabel('时间（分钟）');
ylabel('总速率（Mbps）');
title('（b）平滑后的学习曲线');
grid on; legend('Location', 'best');

nexttile;
improvement_over_time = (sum_rate_prop - sum_rate_base) ./ max(sum_rate_base, eps) * 100;
plot((1:T_total) * dt_s / 60, movmean(improvement_over_time, window_size), 'g-', 'LineWidth', 2, 'DisplayName', '平滑改进幅度'); hold on;
if runBaseline
    yline(improvement_percent, 'g--', 'LineWidth', 1.2, 'DisplayName', sprintf('平均改进：%.2f%%', improvement_percent));
end
yline(0, 'k--', 'LineWidth', 1.0, 'DisplayName', '零线');
xlabel('时间（分钟）');
ylabel('改进幅度（%）');
title('（c）相对性能改进');
grid on; legend('Location', 'best');

nexttile;
early_end = min(50, T_total);
avg_rate_base_per_ue = mean(rate_base_over_time(1:early_end, :), 1) / 1e6;
avg_rate_prop_per_ue = mean(rate_prop_over_time(1:early_end, :), 1) / 1e6;
for c = 1:params.C
    bar(c - 0.2, avg_rate_base_per_ue(c), 0.35, 'FaceColor', ue_colors(c, :), 'EdgeColor', 'none'); hold on;
    bar(c + 0.2, avg_rate_prop_per_ue(c), 0.35, 'FaceColor', min(1, ue_colors(c, :) + 0.25), 'EdgeColor', 'none');
end
xlabel('用户索引');
ylabel('平均速率（Mbps）');
title('（d）用户级性能（前50步）');
grid on;

nexttile;
if T_total > 100
    late_start = T_total - 99;
else
    late_start = 1;
end
avg_rate_base_per_ue_late = mean(rate_base_over_time(late_start:end, :), 1) / 1e6;
avg_rate_prop_per_ue_late = mean(rate_prop_over_time(late_start:end, :), 1) / 1e6;
for c = 1:params.C
    bar(c - 0.2, avg_rate_base_per_ue_late(c), 0.35, 'FaceColor', ue_colors(c, :), 'EdgeColor', 'none'); hold on;
    bar(c + 0.2, avg_rate_prop_per_ue_late(c), 0.35, 'FaceColor', min(1, ue_colors(c, :) + 0.25), 'EdgeColor', 'none');
end
xlabel('用户索引');
ylabel('平均速率（Mbps）');
title('（e）用户级性能（后100步）');
grid on;

nexttile;
avg_wdop_base_per_ue = mean(wdop_base_over_time(late_start:end, :), 1, "omitnan");
avg_wdop_prop_per_ue = mean(wdop_prop_over_time(late_start:end, :), 1, "omitnan");
for c = 1:params.C
    bar(c - 0.2, avg_wdop_base_per_ue(c), 0.35, 'FaceColor', ue_colors(c, :), 'EdgeColor', 'none'); hold on;
    bar(c + 0.2, avg_wdop_prop_per_ue(c), 0.35, 'FaceColor', min(1, ue_colors(c, :) + 0.25), 'EdgeColor', 'none');
end
yline(params.wdopThreshold, 'k--', 'LineWidth', 1.2, 'DisplayName', 'WDOP 阈值');
xlabel('用户索引');
ylabel('平均 WDOP');
title('（f）WDOP 对比');
grid on;

nexttile;
plot((1:T_total) * dt_s / 60, avg_wdop_base_time, 'b-', 'LineWidth', 1.4, 'DisplayName', '基线平均 WDOP'); hold on;
plot((1:T_total) * dt_s / 60, avg_wdop_prop_time, 'r-', 'LineWidth', 1.4, 'DisplayName', '方案平均 WDOP');
yline(params.wdopThreshold, 'k--', 'LineWidth', 1.2, 'DisplayName', '阈值');
xlabel('时间（分钟）');
ylabel('平均 WDOP');
title('（g）WDOP 随时间变化');
grid on; legend('Location', 'best');

nexttile;
avg_group_base = mean(group_rate_base_over_time(late_start:end, :), 1);
avg_group_prop = mean(group_rate_prop_over_time(late_start:end, :), 1);
bar_data = [avg_group_base(:), avg_group_prop(:)];
b = bar(bar_data, 'grouped');
if numGroups > 0
    b(1).FaceColor = color_urllc;
    b(2).FaceColor = color_embb;
end
xlabel('分组索引');
ylabel('平均速率（Mbps）');
title('（h）分组速率汇总');
grid on;
legend({'Baseline', 'Proposal'}, 'Location', 'best');

if params.useParetoUCB
    fig2 = figure("Name", "帕累托前沿诊断", "Color", "w", "Position", [120, 120, 900, 650]);
    c = 1;
    counts = max(N_counts_mab{c}, eps);
    exploration = params.mabCucb * sqrt(log(max(T_total, 2)) ./ counts);
    ucb1 = Q_mab.q1{c} + exploration;
    ucb2 = Q_mab.q2{c} + exploration;
    front = pareto_front_indices(ucb1, ucb2);
    scatter(ucb1, ucb2, 60, [0.7, 0.7, 0.7], 'filled'); hold on;
    scatter(ucb1(front), ucb2(front), 80, 'r', 'filled');
    xlabel('UCB1（通信）');
    ylabel('UCB2（定位）');
    title('用户 1 的帕累托前沿');
    grid on;
    legend({'All arms', 'Pareto front'}, 'Location', 'best');
end

outDir = fullfile(thisDir, "result", "result_fig3_nonstationary");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

timestamp = datestr(now, "yyyymmdd_HHMMSS");
pngPath = fullfile(outDir, sprintf("nonstationary_T%d_%s.png", T_total, timestamp));
saveas(fig, pngPath);

excelPath = fullfile(outDir, sprintf("nonstationary_T%d_%s.xlsx", T_total, timestamp));
summary_table = table((1:T_total).', sum_rate_base, sum_rate_prop, sum_rate_prop - sum_rate_base, improvement_over_time(:), ...
    'VariableNames', {'Time_Step', 'Base_Sum_Rate_Mbps', 'Proposal_Sum_Rate_Mbps', 'Difference_Mbps', 'Improvement_Percent'});
writetable(summary_table, excelPath, 'Sheet', 'Learning_Curve');

stats_table = table( ...
    {'Avg_Rate_Base'; 'Avg_Rate_Proposal'; 'Mean_Improvement_Percent'; 'Avg_WDOP_Base'; 'Avg_WDOP_Proposal'}, ...
    [avg_rate_base; avg_rate_prop; improvement_percent; avg_wdop_base; avg_wdop_prop], ...
    'VariableNames', {'Metric', 'Value'});
writetable(stats_table, excelPath, 'Sheet', 'Statistics');

ue_comparison = table((1:params.C).', avg_rate_base_per_ue_late(:), avg_rate_prop_per_ue_late(:), avg_wdop_base_per_ue(:), avg_wdop_prop_per_ue(:), ...
    'VariableNames', {'UE', 'Base_Rate_Late_Mbps', 'Prop_Rate_Late_Mbps', 'Base_WDOP_Late', 'Prop_WDOP_Late'});
writetable(ue_comparison, excelPath, 'Sheet', 'UE_Comparison');

group_table = table((1:numGroups).', cellfun(@numel, user_groups.group_ids(:)), avg_group_base(:), avg_group_prop(:), ...
    'VariableNames', {'Group', 'Num_Users', 'Base_Avg_Rate_Mbps', 'Prop_Avg_Rate_Mbps'});
writetable(group_table, excelPath, 'Sheet', 'Group_Comparison');

matPath = fullfile(outDir, sprintf("nonstationary_T%d_%s.mat", T_total, timestamp));
save(matPath, 'sum_rate_base', 'sum_rate_prop', 'sum_rate_base_smooth', 'sum_rate_prop_smooth', ...
    'rate_base_over_time', 'rate_prop_over_time', 'wdop_base_over_time', 'wdop_prop_over_time', ...
    'group_rate_base_over_time', 'group_rate_prop_over_time', 'user_groups', 'similarity_matrix', ...
    'params', 'improvement_percent', 'avg_rate_base', 'avg_rate_prop', 'avg_wdop_base', 'avg_wdop_prop', ...
    'T_total', 'dt_s');

fprintf("图像已保存至：%s\n", pngPath);
fprintf("Excel 已保存至：%s\n", excelPath);
fprintf("MAT 已保存至：%s\n", matPath);

fprintf('\n========== 非稳态仿真结果摘要 ==========');
fprintf('\n仿真时长：%d 步\n', T_total);
fprintf('评估区间：步 %d-%d\n', eval_start, T_total);
if runBaseline
    fprintf('基线平均总和速率：%.3f Mbps\n', avg_rate_base);
else
    fprintf('基线平均总和速率：N/A（Baseline 关闭）\n');
end
fprintf('方案平均总和速率：%.3f Mbps\n', avg_rate_prop);
if runBaseline
    fprintf('性能改进：%.2f%%\n', improvement_percent);
else
    fprintf('性能改进：N/A（Baseline 关闭）\n');
end
fprintf('=================================\n');

function scenario = update_satellite_positions(scenario, t, dt_s)
    if ~isfield(scenario, 'satDynamics')
        return;
    end

    dyn = scenario.satDynamics;
    S = size(scenario.pSat, 1);
    tau = (t - 1) * dt_s;

    for s = 1:S
        r = dyn.radius_m(s);
        Omega = dyn.raan_rad(s);
        inc = dyn.inc_rad(s);
        theta = dyn.phase0_rad(s) + dyn.omega_radps(s) * tau;

        % 圆轨道近似：轨道平面坐标 -> RAAN/倾角旋转到全局坐标
        x_orb = r * cos(theta);
        y_orb = r * sin(theta);

        cO = cos(Omega); sO = sin(Omega);
        ci = cos(inc);   si = sin(inc);

        x = cO * x_orb - sO * ci * y_orb;
        y = sO * x_orb + cO * ci * y_orb;
        z = si * y_orb;

        scenario.pSat(s, :) = [x, y, z];
    end
end

function print_full_arm_space(params, scenario, chan)
    S = params.S;
    C = params.C;
    I = params.I;
    comb = nchoosek(1:S, I);
    fprintf("=== 打印全臂空间：共 %d 个组合/每用户 ===\n", size(comb, 1));
    for c = 1:C
        fprintf("[ArmSpace][UE%d] 开始\n", c);
        for k = 1:size(comb, 1)
            sats = comb(k, :);
            d_vec = chan.d_m(sats, c);
            wdop_val = ican.compute_wdop(scenario.pUE(c, :), scenario.pSat(sats, :), d_vec);
            h_sel = chan.h(:, c, sats);
            rate_proxy = sum(abs(h_sel).^2, "all");
            fprintf("[ArmSpace][UE%d] arm#%d sats=%s WDOP=%.3f rateProxy=%.6g\n", ...
                c, k, mat2str(sats), wdop_val, rate_proxy);
        end
        fprintf("[ArmSpace][UE%d] 结束\n", c);
    end
end

function user_groups = normalize_user_groups(user_groups, C)
    if isfield(user_groups, 'group_ids')
        if ~isfield(user_groups, 'weights') || isempty(user_groups.weights)
            user_groups.weights = ones(numel(user_groups.group_ids), 1);
        end
        if ~isfield(user_groups, 'group_assignment') || isempty(user_groups.group_assignment)
            user_groups.group_assignment = zeros(C, 1);
            for k = 1:numel(user_groups.group_ids)
                ids = user_groups.group_ids{k};
                user_groups.group_assignment(ids(:)) = k;
            end
        end
        user_groups.numGroups = numel(user_groups.group_ids);
        return;
    end
end

function user_groups = default_user_groups_all_users(C)
    user_groups = struct();
    user_groups.group_ids = {(1:C).'};
    user_groups.weights = 1;
    user_groups.group_assignment = ones(C, 1);
    user_groups.numGroups = 1;
end

function rate_per_group = compute_group_rates(rate_row, user_groups)
    numGroups = numel(user_groups.group_ids);
    rate_per_group = zeros(1, numGroups);
    for k = 1:numGroups
        ids = user_groups.group_ids{k};
        ids = ids(ids >= 1 & ids <= numel(rate_row));
        if isempty(ids)
            rate_per_group(k) = 0;
        else
            rate_per_group(k) = mean(rate_row(ids), 'omitnan');
        end
    end
end

function front = pareto_front_indices(x, y)
    n = numel(x);
    isPareto = true(n, 1);
    for i = 1:n
        if ~isPareto(i)
            continue;
        end
        dominated = (x >= x(i) & y >= y(i)) & ((x > x(i)) | (y > y(i)));
        isPareto(dominated) = false;
    end
    front = find(isPareto);
end


