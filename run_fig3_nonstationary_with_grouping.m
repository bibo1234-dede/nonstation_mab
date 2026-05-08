clear; clc; clear functions;

thisDir = fileparts(mfilename("fullpath"));
addpath(thisDir);

params = ican.config_paper("randomSeed", 1, "I", 4, "S", 7, "C", 7);
params.log.level = "info";
params.cvx.quiet = false;  % ← 改成 false 看详细输出，诊断 CVX 问题
params.log.dir = "";  % 输出到屏幕（命令窗口）而不是文件

% ========== Key parameters ==========
T_total = 200;
dt_s = 1;
sigma_disturbance = 0.05;

ican.logf(params, "info", "=== Proposal-only simulation: Non-stationary MAB + User Grouping ===");
ican.logf(params, "info", "Params: S=%d C=%d I=%d T_total=%d WDOP_threshold=%.3f", ...
    params.S, params.C, params.I, T_total, params.wdopThreshold);

% ========== CVX check ==========
if exist("cvx_begin", "file") ~= 2
    error("CVX not found on MATLAB path. Please install CVX and run cvx_setup.");
end
try
    cvx_solver(char(params.cvx.solver));
catch solverErr
    error("Failed to set CVX solver '%s': %s", string(params.cvx.solver), solverErr.message);
end

% ========== Scenario and state initialization ==========
scenario = ican.create_scenario_fig3(params);
ican.logf(params, "info", "Scenario: satRingRadius=%.1f km", scenario.satRingRadius_m/1e3);

Q_mab = [];
N_counts_mab = [];
t_last_served_mab = [];
best_U_mab = -inf;
best_alpha_mab = [];

user_groups = [];
group_stats = [];
user_group_assignment = zeros(params.C, 1);

rate_prop_over_time = zeros(T_total, params.C);
wdop_prop_over_time = zeros(T_total, params.C);

rate_urllc_over_time = zeros(T_total, 1);
rate_embb_over_time = zeros(T_total, 1);
rate_mmtc_over_time = zeros(T_total, 1);

% ========== AR(1) Shadow Fading Initialization ==========
% 一阶高斯-马尔可夫阴影衰落 (时间相关的非稳态模型)
shadow_dB = zeros(params.S, params.C);  % 初值: 0 dB
a_corr = 0.9;  % AR(1) 相关系数 (0.8~0.95, 0.9 表示强相关)
sigma_shadow = 6;  % 对数标准差 (dB), 典型值 4~8 dB
ican.logf(params, "info", "Shadow fading: a_corr=%.2f, sigma=%.1f dB", a_corr, sigma_shadow);

chan_base = ican.compute_channels(params, scenario);

fprintf('\n=== 基础信道计算完成 ===\n'); drawnow;
ican.logf(params, "info", "=== Start time-series simulation, %d steps ===", T_total);

% ========== Main loop ==========
fprintf('=== 进入主循环 ===\n'); drawnow;
for t = 1:T_total
    if mod(t, 10) == 1 || t <= 3 || t == T_total
        fprintf('[Progress] t=%d/%d ...\n', t, T_total);
        drawnow;
    end
    
    if t == 1
        chan = chan_base;
    else
        % AR(1) 阴影衰落更新: shadow_dB(t) = a*shadow_dB(t-1) + sqrt(1-a^2)*sigma*randn()
        shadow_dB = a_corr * shadow_dB + sqrt(1 - a_corr^2) * sigma_shadow * randn(params.S, params.C);
        shadow_linear = 10.^(shadow_dB / 20);  % dB 转换为线性值
        % 阴影衰落作用在信道上 (乘法因子)
        % reshape [S, C] → [1, C, S] 用于隐式广播到 [N, C, S]
        shadow_fading = reshape(shadow_linear, 1, params.C, params.S);  % [1, C, S]
        chan.h = chan.h .* shadow_fading;  % [N, C, S] .* [1, C, S]
        % 距离和方向余弦保持不变 (保留几何结构)
    end

    if t == 1
        [user_groups, group_stats] = ican.user_grouping_strategy(params, chan, scenario);
        fprintf('\n========== User Grouping Completed ==========\n');
        drawnow;

        for c = 1:params.C
            if ismember(c, user_groups.URLLC.user_ids)
                user_group_assignment(c) = 1;
            elseif ismember(c, user_groups.eMBB.user_ids)
                user_group_assignment(c) = 2;
            else
                user_group_assignment(c) = 3;
            end
        end

        fprintf('\n');
        fprintf('========== 用户分组结果 ==========\n');
        drawnow;
        fprintf('URLLC用户 (优先级 %.1f): %s\n', ...
            user_groups.URLLC.priority, mat2str(user_groups.URLLC.user_ids));
        fprintf('  资源分配: %.1f%% | 延迟要求: %d ms | 可靠性: %.3f%%\n', ...
            user_groups.URLLC.resource_ratio*100, user_groups.URLLC.delay_ms, user_groups.URLLC.reliability_percent);
        fprintf('  平均SNR: %.2f dB | 平均信道质量: %.3f\n', ...
            group_stats.avg_SNR_per_group.URLLC, group_stats.avg_channel_quality_per_group.URLLC);

        fprintf('\neMBB用户 (优先级 %.1f): %s\n', ...
            user_groups.eMBB.priority, mat2str(user_groups.eMBB.user_ids));
        fprintf('  资源分配: %.1f%% | 延迟要求: %d ms | 可靠性: %.1f%%\n', ...
            user_groups.eMBB.resource_ratio*100, user_groups.eMBB.delay_ms, user_groups.eMBB.reliability_percent);
        fprintf('  平均SNR: %.2f dB | 平均信道质量: %.3f\n', ...
            group_stats.avg_SNR_per_group.eMBB, group_stats.avg_channel_quality_per_group.eMBB);

        fprintf('\nmMTC用户 (优先级 %.1f): %s\n', ...
            user_groups.mMTC.priority, mat2str(user_groups.mMTC.user_ids));
        fprintf('  资源分配: %.1f%% | 延迟要求: %d ms | 可靠性: %.1f%%\n', ...
            user_groups.mMTC.resource_ratio*100, user_groups.mMTC.delay_ms, user_groups.mMTC.reliability_percent);
        fprintf('  平均SNR: %.2f dB | 平均信道质量: %.3f\n', ...
            group_stats.avg_SNR_per_group.mMTC, group_stats.avg_channel_quality_per_group.mMTC);
        fprintf('===================================\n\n');
        drawnow;
    end

    fprintf('  [t=%d] Calling MAB selection...\n', t);
    drawnow;
    [out_mab, Q_mab, N_counts_mab, t_last_served_mab] = ...
        ican.select_satellites_mab_wdop_dynamic_ucb( ...
            params, scenario, chan, Q_mab, N_counts_mab, ...
            best_U_mab, best_alpha_mab, t, t_last_served_mab, user_groups);
    fprintf('  [t=%d] MAB done, calling beamforming...\n', t);
    drawnow;

    propBf = ican.solve_beamforming_dc(params, chan, out_mab.alpha);
    fprintf('  [t=%d] Beamforming done\n', t);
    drawnow;
    rate_prop_over_time(t, :) = propBf.R_c_bps;

    for c = 1:params.C
        sats = find(out_mab.alpha(:, c) > 0.5);
        if ~isempty(sats)
            d_vec = chan.d_m(sats, c);
            wdop_prop_over_time(t, c) = ican.compute_wdop( ...
                scenario.pUE(c, :), scenario.pSat(sats, :), d_vec);
        else
            wdop_prop_over_time(t, c) = inf;
        end
    end

    rate_urllc_over_time(t) = safe_group_mean(rate_prop_over_time(t, :), user_groups.URLLC.user_ids) / 1e9;
    rate_embb_over_time(t)  = safe_group_mean(rate_prop_over_time(t, :), user_groups.eMBB.user_ids) / 1e9;
    rate_mmtc_over_time(t)  = safe_group_mean(rate_prop_over_time(t, :), user_groups.mMTC.user_ids) / 1e9;

    best_U_mab = out_mab.best_utility;
    best_alpha_mab = out_mab.best_alpha;

    if mod(t, 5) == 0
        sum_prop = sum(propBf.R_c_bps) / 1e9;
        sat_summary = strings(params.C, 1);
        for c = 1:params.C
            sats = find(out_mab.alpha(:, c) > 0.5);
            sat_summary(c) = sprintf("UE%d=%s", c, mat2str(sats(:).'));
        end
        ican.logf(params, "info", "t=%d/%d: Proposal=%.3f Gbps | %s", ...
            t, T_total, sum_prop, strjoin(cellstr(sat_summary), ", "));
    end
end

fprintf('\n========== Main Simulation Loop Completed ==========\n\n');
drawnow;

% ========== Statistics ==========
time_axis = (1:T_total);
time_min = time_axis * dt_s / 60;
window_size = 20;
eval_start = max(1, T_total - 99);

sum_rate_prop = sum(rate_prop_over_time, 2) / 1e9;
sum_rate_prop_smooth = movmean(sum_rate_prop, window_size);

avg_rate_prop = mean(sum_rate_prop(eval_start:end));
avg_rate_urllc = mean(rate_urllc_over_time(eval_start:end));
avg_rate_embb = mean(rate_embb_over_time(eval_start:end));
avg_rate_mmtc = mean(rate_mmtc_over_time(eval_start:end));

avg_wdop_prop = mean(wdop_prop_over_time(eval_start:end, :), "all", "omitnan");
wdop_satisfaction_percent = mean(wdop_prop_over_time(:) <= params.wdopThreshold) * 100;

ican.logf(params, "info", "=== Simulation completed ===");
ican.logf(params, "info", "Average Proposal Sum Rate (last 100 steps): %.3f Gbps", avg_rate_prop);
ican.logf(params, "info", "Average Proposal WDOP (last 100 steps): %.3f", avg_wdop_prop);
ican.logf(params, "info", "Proposal WDOP satisfaction ratio: %.2f%%", wdop_satisfaction_percent);

% ========== Plot ==========
color_urllc = [1, 0, 0];
color_embb = [0.2, 0.6, 1];
color_mmtc = [0, 0.7, 0];

colors_per_ue = zeros(params.C, 3);
for c = 1:params.C
    if user_group_assignment(c) == 1
        colors_per_ue(c, :) = color_urllc;
    elseif user_group_assignment(c) == 2
        colors_per_ue(c, :) = color_embb;
    else
        colors_per_ue(c, :) = color_mmtc;
    end
end

fig = figure("Name", "Proposal Only: Non-stationary MAB + User Grouping", ...
    "Color", "w", "Position", [100, 100, 1400, 850]);
tiledlayout(2, 3, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(time_min, sum_rate_prop, 'r-', 'LineWidth', 1, 'DisplayName', 'Proposal');
xlabel('Time (minutes)');
ylabel('Sum Rate (Gbps)');
title('(a) Proposal Sum Rate');
grid on;
legend('Location', 'best');

nexttile;
plot(time_min, sum_rate_prop_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Proposal MA-20'); hold on;
yline(avg_rate_prop, 'r--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Avg: %.3f Gbps', avg_rate_prop));
xlabel('Time (minutes)');
ylabel('Sum Rate (Gbps)');
title('(b) Smoothed Sum Rate');
grid on;
legend('Location', 'best');

nexttile;
plot(time_min, rate_urllc_over_time, 'Color', color_urllc, 'LineWidth', 2, 'DisplayName', 'URLLC'); hold on;
plot(time_min, rate_embb_over_time, 'Color', color_embb, 'LineWidth', 2, 'DisplayName', 'eMBB');
plot(time_min, rate_mmtc_over_time, 'Color', color_mmtc, 'LineWidth', 2, 'DisplayName', 'mMTC');
xlabel('Time (minutes)');
ylabel('Avg Rate per UE (Gbps)');
title('(c) Group Rate Evolution');
grid on;
legend('Location', 'best');

nexttile;
avg_rate_prop_per_ue_late = mean(rate_prop_over_time(eval_start:end, :), 1) / 1e6;
ue_idx = 1:params.C;
for c = 1:params.C
    bar(ue_idx(c), avg_rate_prop_per_ue_late(c), 0.55, 'FaceColor', colors_per_ue(c, :)); hold on;
end
xlabel('UE Index');
ylabel('Avg Rate (Mbps)');
title('(d) Per-UE Rate (Last 100 Steps)');
grid on;

nexttile;
avg_wdop_prop_per_ue = mean(wdop_prop_over_time(eval_start:end, :), 1, "omitnan");
for c = 1:params.C
    bar(ue_idx(c), avg_wdop_prop_per_ue(c), 0.55, 'FaceColor', colors_per_ue(c, :)); hold on;
end
yline(params.wdopThreshold, 'k--', 'LineWidth', 1.2, 'DisplayName', 'WDOP threshold');
xlabel('UE Index');
ylabel('Average WDOP');
title('(e) Per-UE WDOP (Last 100 Steps)');
grid on;
legend('Location', 'best');

nexttile;
groups = {'URLLC', 'eMBB', 'mMTC'};
bar_data = [avg_rate_urllc, avg_rate_embb, avg_rate_mmtc];
b = bar(bar_data);
b.FaceColor = 'flat';
b.CData = [color_urllc; color_embb; color_mmtc];
set(gca, 'XTickLabel', groups);
ylabel('Average Rate (Gbps)');
title('(f) QoS Group Summary');
grid on;
box off;

% ========== Save results ==========
outDir = fullfile(thisDir, "result", "result_fig3_nonstationary");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

timestamp = datestr(now, "yyyymmdd_HHMMSS");

pngPath = fullfile(outDir, sprintf("proposal_only_nonstationary_T%d_grouping_%s.png", T_total, timestamp));
saveas(fig, pngPath);
ican.logf(params, "info", "Figure saved to: %s", pngPath);

excelPath = fullfile(outDir, sprintf("proposal_only_nonstationary_T%d_grouping_%s.xlsx", T_total, timestamp));

summary_table = table( ...
    time_axis(:), ...
    sum_rate_prop, ...
    rate_urllc_over_time, ...
    rate_embb_over_time, ...
    rate_mmtc_over_time, ...
    mean(wdop_prop_over_time, 2, "omitnan"), ...
    'VariableNames', {'Time_Step', 'Proposal_Sum_Rate_Gbps', ...
                      'URLLC_Avg_Rate_Gbps', 'eMBB_Avg_Rate_Gbps', 'mMTC_Avg_Rate_Gbps', ...
                      'Proposal_Avg_WDOP'});
writetable(summary_table, excelPath, 'Sheet', 'Learning_Curve');

stats_table = table( ...
    {'Avg_Rate_Proposal'; 'Avg_WDOP_Proposal'; 'WDOP_Satisfaction_Percent'; ...
     'URLLC_Avg_Rate'; 'eMBB_Avg_Rate'; 'mMTC_Avg_Rate'}, ...
    [avg_rate_prop; avg_wdop_prop; wdop_satisfaction_percent; ...
     avg_rate_urllc; avg_rate_embb; avg_rate_mmtc], ...
    'VariableNames', {'Metric', 'Value'});
writetable(stats_table, excelPath, 'Sheet', 'Statistics');

qos_types = cell(params.C, 1);
snr_per_user = zeros(params.C, 1);
channel_quality_per_user = zeros(params.C, 1);

for c = 1:params.C
    if user_group_assignment(c) == 1
        qos_types{c} = 'URLLC';
        snr_per_user(c) = group_stats.avg_SNR_per_group.URLLC;
        channel_quality_per_user(c) = group_stats.avg_channel_quality_per_group.URLLC;
    elseif user_group_assignment(c) == 2
        qos_types{c} = 'eMBB';
        snr_per_user(c) = group_stats.avg_SNR_per_group.eMBB;
        channel_quality_per_user(c) = group_stats.avg_channel_quality_per_group.eMBB;
    else
        qos_types{c} = 'mMTC';
        snr_per_user(c) = group_stats.avg_SNR_per_group.mMTC;
        channel_quality_per_user(c) = group_stats.avg_channel_quality_per_group.mMTC;
    end
end

ue_comparison = table( ...
    (1:params.C)', ...
    qos_types, ...
    snr_per_user, ...
    channel_quality_per_user, ...
    mean(rate_prop_over_time(1:min(50, T_total), :), 1)' / 1e6, ...
    mean(rate_prop_over_time(eval_start:end, :), 1)' / 1e6, ...
    mean(wdop_prop_over_time(eval_start:end, :), 1, "omitnan")', ...
    'VariableNames', {'UE', 'QoS_Type', 'SNR_dB', 'Channel_Quality', ...
                      'Proposal_Rate_Early_Mbps', 'Proposal_Rate_Late_Mbps', ...
                      'Proposal_WDOP_Late'});
writetable(ue_comparison, excelPath, 'Sheet', 'UE_Comparison');

group_details = table( ...
    {'URLLC'; 'eMBB'; 'mMTC'}, ...
    [user_groups.URLLC.priority; user_groups.eMBB.priority; user_groups.mMTC.priority], ...
    [user_groups.URLLC.resource_ratio; user_groups.eMBB.resource_ratio; user_groups.mMTC.resource_ratio], ...
    [user_groups.URLLC.delay_ms; user_groups.eMBB.delay_ms; user_groups.mMTC.delay_ms], ...
    [user_groups.URLLC.reliability_percent; user_groups.eMBB.reliability_percent; user_groups.mMTC.reliability_percent], ...
    [group_stats.URLLC_count; group_stats.eMBB_count; group_stats.mMTC_count], ...
    [group_stats.avg_SNR_per_group.URLLC; group_stats.avg_SNR_per_group.eMBB; group_stats.avg_SNR_per_group.mMTC], ...
    'VariableNames', {'QoS_Type', 'Priority', 'Resource_Ratio', 'Delay_ms', ...
                      'Reliability_Percent', 'User_Count', 'Avg_SNR_dB'});
writetable(group_details, excelPath, 'Sheet', 'QoS_Groups');

ican.logf(params, "info", "Excel saved to: %s", excelPath);

matPath = fullfile(outDir, sprintf("proposal_only_nonstationary_T%d_grouping_%s.mat", T_total, timestamp));
save(matPath, 'sum_rate_prop', 'sum_rate_prop_smooth', ...
    'rate_prop_over_time', 'wdop_prop_over_time', ...
    'rate_urllc_over_time', 'rate_embb_over_time', 'rate_mmtc_over_time', ...
    'user_groups', 'group_stats', 'user_group_assignment', ...
    'params', 'avg_rate_prop', 'avg_wdop_prop', 'wdop_satisfaction_percent', ...
    'avg_rate_urllc', 'avg_rate_embb', 'avg_rate_mmtc', 'T_total', 'dt_s', 'sigma_disturbance');

ican.logf(params, "info", "MAT saved to: %s", matPath);

fprintf('\n');
fprintf('========== Proposal-only Simulation Summary ==========\n');
fprintf('Simulation length: %d steps\n', T_total);
fprintf('Evaluation window: steps %d-%d\n', eval_start, T_total);
fprintf('Proposal average sum rate: %.3f Gbps\n', avg_rate_prop);
fprintf('Proposal average WDOP: %.3f\n', avg_wdop_prop);
fprintf('WDOP satisfaction ratio: %.2f%%\n', wdop_satisfaction_percent);
fprintf('\n');
fprintf('Group performance:\n');
fprintf('  URLLC: %.3f Gbps (priority %.1f)\n', avg_rate_urllc, user_groups.URLLC.priority);
fprintf('  eMBB:  %.3f Gbps (priority %.1f)\n', avg_rate_embb, user_groups.eMBB.priority);
fprintf('  mMTC:  %.3f Gbps (priority %.1f)\n', avg_rate_mmtc, user_groups.mMTC.priority);
fprintf('======================================================\n\n');

function val = safe_group_mean(data_row, ids)
    if isempty(ids)
        val = 0;
        return;
    end

    ids = ids(ids <= length(data_row));
    if isempty(ids)
        val = 0;
        return;
    end

    x = data_row(ids);
    if all(isnan(x))
        val = 0;
    else
        val = mean(x(~isnan(x)));
    end
end
