function [user_groups, group_stats] = user_grouping_strategy(params, chan, scenario)
% USER_GROUPING_STRATEGY Group users by QoS class using the actual channel.
%
% Inputs:
%   params   - system parameters
%   chan     - channel struct with fields h and d_m
%   scenario - scenario struct (kept for interface compatibility)
%
% Outputs:
%   user_groups - grouping information
%   group_stats - grouping statistics

if isempty(params) || ~isfield(params, 'C') || ~isfield(params, 'S') || ~isfield(params, 'N')
    error('Params must contain C, S, and N fields.');
end

if isempty(chan) || ~isfield(chan, 'd_m') || ~isfield(chan, 'h')
    error('Chan must contain both d_m and h fields.');
end

C = params.C;
S = params.S;

user_snr_db = zeros(C, 1);
user_channel_quality = zeros(C, 1);
user_path_loss = zeros(C, 1);

for c = 1:C
    distances = chan.d_m(:, c);
    % squeeze 移除单维度 [N, 1, S] → [N, S]，保持内存顺序一致
    channel_matrix = squeeze(chan.h(:, c, :));  % [N, S]
    channel_norms = vecnorm(channel_matrix, 2, 1);
    channel_powers = channel_norms .^ 2;

    % Use the same power/noise quantities as the downstream beamforming code.
    snr_linear = (params.P_W .* channel_powers) ./ max(params.sigma2_W, eps);
    snr_db = 10 * log10(max(snr_linear, realmin));

    [~, best_idx] = sort(snr_db, 'descend');
    num_best = min(3, S);
    strongest_idx = best_idx(1:num_best);

    user_snr_db(c) = mean(snr_db(strongest_idx));
    user_channel_quality(c) = mean(channel_norms(strongest_idx));
    user_path_loss(c) = mean(20 * log10(4 * pi * distances / params.lambda_m));
end

quality_scale = max(user_channel_quality);
if quality_scale > 0
    user_channel_quality = user_channel_quality / quality_scale;
end

SNR_threshold_URLLC = 20;
SNR_threshold_eMBB = 10;
channel_quality_threshold = 0.3;

URLLC_candidates = find(user_snr_db > SNR_threshold_URLLC & ...
                        user_channel_quality > channel_quality_threshold);
eMBB_candidates = find(user_snr_db > SNR_threshold_eMBB & ...
                       user_channel_quality > channel_quality_threshold & ...
                       ~ismember((1:C)', URLLC_candidates));
mMTC_candidates = setdiff(1:C, [URLLC_candidates; eMBB_candidates]);

target_ratio_URLLC = 0.50;
target_ratio_eMBB = 0.30;
target_ratio_mMTC = 0.20;

target_num_URLLC = max(1, round(C * target_ratio_URLLC));
target_num_eMBB = max(1, round(C * target_ratio_eMBB));
target_num_mMTC = max(1, C - target_num_URLLC - target_num_eMBB);

URLLC_users = URLLC_candidates(1:min(numel(URLLC_candidates), target_num_URLLC))';
eMBB_users = eMBB_candidates(1:min(numel(eMBB_candidates), target_num_eMBB))';
mMTC_users = mMTC_candidates(1:min(numel(mMTC_candidates), target_num_mMTC))';

remaining_users = setdiff(1:C, [URLLC_users, eMBB_users, mMTC_users]);
for user_id = remaining_users
    if numel(URLLC_users) < target_num_URLLC
        URLLC_users = [URLLC_users, user_id];
    elseif numel(eMBB_users) < target_num_eMBB
        eMBB_users = [eMBB_users, user_id];
    else
        mMTC_users = [mMTC_users, user_id];
    end
end

URLLC_users = URLLC_users(:);
eMBB_users = eMBB_users(:);
mMTC_users = mMTC_users(:);

user_groups.URLLC.user_ids = URLLC_users;
user_groups.URLLC.priority = 3.5;
user_groups.URLLC.resource_ratio = 0.50;
user_groups.URLLC.delay_ms = 10;
user_groups.URLLC.reliability_percent = 99.999;
user_groups.URLLC.weight = 1.2;
user_groups.URLLC.num_users = numel(URLLC_users);
user_groups.URLLC.avg_snr_db = safe_mean(user_snr_db(URLLC_users));
user_groups.URLLC.avg_channel_quality = safe_mean(user_channel_quality(URLLC_users));
user_groups.URLLC.avg_path_loss_db = safe_mean(user_path_loss(URLLC_users));

user_groups.eMBB.user_ids = eMBB_users;
user_groups.eMBB.priority = 2.5;
user_groups.eMBB.resource_ratio = 0.30;
user_groups.eMBB.delay_ms = 50;
user_groups.eMBB.reliability_percent = 99.0;
user_groups.eMBB.weight = 1.0;
user_groups.eMBB.num_users = numel(eMBB_users);
user_groups.eMBB.avg_snr_db = safe_mean(user_snr_db(eMBB_users));
user_groups.eMBB.avg_channel_quality = safe_mean(user_channel_quality(eMBB_users));
user_groups.eMBB.avg_path_loss_db = safe_mean(user_path_loss(eMBB_users));

user_groups.mMTC.user_ids = mMTC_users;
user_groups.mMTC.priority = 1.5;
user_groups.mMTC.resource_ratio = 0.20;
user_groups.mMTC.delay_ms = 500;
user_groups.mMTC.reliability_percent = 95.0;
user_groups.mMTC.weight = 0.8;
user_groups.mMTC.num_users = numel(mMTC_users);
user_groups.mMTC.avg_snr_db = safe_mean(user_snr_db(mMTC_users));
user_groups.mMTC.avg_channel_quality = safe_mean(user_channel_quality(mMTC_users));
user_groups.mMTC.avg_path_loss_db = safe_mean(user_path_loss(mMTC_users));

group_stats.total_users = C;
group_stats.total_satellites = S;
group_stats.URLLC_count = numel(URLLC_users);
group_stats.eMBB_count = numel(eMBB_users);
group_stats.mMTC_count = numel(mMTC_users);

group_stats.URLLC_percent = (numel(URLLC_users) / C) * 100;
group_stats.eMBB_percent = (numel(eMBB_users) / C) * 100;
group_stats.mMTC_percent = (numel(mMTC_users) / C) * 100;

group_stats.avg_SNR_per_group.URLLC = user_groups.URLLC.avg_snr_db;
group_stats.avg_SNR_per_group.eMBB = user_groups.eMBB.avg_snr_db;
group_stats.avg_SNR_per_group.mMTC = user_groups.mMTC.avg_snr_db;

group_stats.avg_channel_quality_per_group.URLLC = user_groups.URLLC.avg_channel_quality;
group_stats.avg_channel_quality_per_group.eMBB = user_groups.eMBB.avg_channel_quality;
group_stats.avg_channel_quality_per_group.mMTC = user_groups.mMTC.avg_channel_quality;

group_stats.resource_allocation_table = table( ...
    {'URLLC'; 'eMBB'; 'mMTC'}, ...
    [user_groups.URLLC.num_users; user_groups.eMBB.num_users; user_groups.mMTC.num_users], ...
    [group_stats.URLLC_percent; group_stats.eMBB_percent; group_stats.mMTC_percent], ...
    [user_groups.URLLC.resource_ratio*100; user_groups.eMBB.resource_ratio*100; user_groups.mMTC.resource_ratio*100], ...
    [user_groups.URLLC.priority; user_groups.eMBB.priority; user_groups.mMTC.priority], ...
    [user_groups.URLLC.delay_ms; user_groups.eMBB.delay_ms; user_groups.mMTC.delay_ms], ...
    [user_groups.URLLC.reliability_percent; user_groups.eMBB.reliability_percent; user_groups.mMTC.reliability_percent], ...
    'VariableNames', {'QoS_Type', 'Num_Users', 'User_Percent', 'Resource_Percent', 'Priority', 'Delay_ms', 'Reliability_Percent'});

fprintf('\n========== 用户分组详细信息 ==========\n');
fprintf('总用户数: %d | 总卫星数: %d\n\n', C, S);

fprintf('URLLC分组:\n');
fprintf('  用户索引: %s\n', mat2str(URLLC_users'));
fprintf('  用户数: %d (%.1f%%) | 优先级: %.1f\n', ...
    user_groups.URLLC.num_users, group_stats.URLLC_percent, user_groups.URLLC.priority);
fprintf('  资源分配: %.1f%% | 延迟要求: %d ms | 可靠性: %.3f%%\n', ...
    user_groups.URLLC.resource_ratio*100, user_groups.URLLC.delay_ms, user_groups.URLLC.reliability_percent);
fprintf('  平均SNR: %.2f dB | 平均信道质量: %.3f\n\n', ...
    user_groups.URLLC.avg_snr_db, user_groups.URLLC.avg_channel_quality);

fprintf('eMBB分组:\n');
fprintf('  用户索引: %s\n', mat2str(eMBB_users'));
fprintf('  用户数: %d (%.1f%%) | 优先级: %.1f\n', ...
    user_groups.eMBB.num_users, group_stats.eMBB_percent, user_groups.eMBB.priority);
fprintf('  资源分配: %.1f%% | 延迟要求: %d ms | 可靠性: %.1f%%\n', ...
    user_groups.eMBB.resource_ratio*100, user_groups.eMBB.delay_ms, user_groups.eMBB.reliability_percent);
fprintf('  平均SNR: %.2f dB | 平均信道质量: %.3f\n\n', ...
    user_groups.eMBB.avg_snr_db, user_groups.eMBB.avg_channel_quality);

fprintf('mMTC分组:\n');
fprintf('  用户索引: %s\n', mat2str(mMTC_users'));
fprintf('  用户数: %d (%.1f%%) | 优先级: %.1f\n', ...
    user_groups.mMTC.num_users, group_stats.mMTC_percent, user_groups.mMTC.priority);
fprintf('  资源分配: %.1f%% | 延迟要求: %d ms | 可靠性: %.1f%%\n', ...
    user_groups.mMTC.resource_ratio*100, user_groups.mMTC.delay_ms, user_groups.mMTC.reliability_percent);
fprintf('  平均SNR: %.2f dB | 平均信道质量: %.3f\n\n', ...
    user_groups.mMTC.avg_snr_db, user_groups.mMTC.avg_channel_quality);

fprintf('========== 资源分配统计 ==========\n');
disp(group_stats.resource_allocation_table);
fprintf('==================================\n\n');

end

function val = safe_mean(x)
if isempty(x)
    val = 0;
else
    val = mean(x);
end
end
