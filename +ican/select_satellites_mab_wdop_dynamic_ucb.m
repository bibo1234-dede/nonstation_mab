function [out, Q, N_counts, t_last_served] = select_satellites_mab_wdop_dynamic_ucb( ...
    params, scenario, chan, Q, N_counts, best_U, best_alpha, t, t_last_served, user_groups)
% SELECT_SATELLITES_MAB_WDOP_DYNAMIC_UCB  面向 WDOP 的非稳态 MAB 选星器。

S = params.S;
C = params.C;
I = params.I;

if nargin < 10 || isempty(user_groups)
    user_groups = default_user_groups(C);
end

rho = get_param_value(params, "mabRho", 0.98);
cUcb = get_param_value(params, "mabCucb", 1.0);
maxArms = get_param_value(params, "mabMaxArms", 20);
wdopThreshold = get_param_value(params, "wdopThreshold", 6.0);
wdopPenaltyLambda = get_param_value(params, "wdopPenaltyLambda", 0.5);
useParetoUCB = get_param_value(params, "useParetoUCB", false);
paretoAlpha = get_param_value(params, "paretoAlpha", 0.5);

comb = nchoosek(1:S, I);
arm_bank = build_candidate_arms(params, scenario, chan, comb, C, maxArms, wdopThreshold, wdopPenaltyLambda);

if isempty(Q) || ~isstruct(Q) || ~isfield(Q, "q1") || ~isfield(Q, "q2") || ~isfield(Q, "arms")
    Q = struct();
    Q.q1 = cell(C, 1);
    Q.q2 = cell(C, 1);
    Q.arms = cell(C, 1);
    Q.last_t = max(t - 1, 0);
end
if ~isfield(Q, "last_t")
    Q.last_t = max(t - 1, 0);
end

if isempty(N_counts) || ~iscell(N_counts)
    N_counts = cell(C, 1);
end

[Q, N_counts] = align_state_to_arms(Q, N_counts, arm_bank, C);

elapsed = max(0, t - Q.last_t);
if elapsed > 0
    decay = rho^elapsed;
    for c = 1:C
        if ~isempty(N_counts{c})
            N_counts{c} = N_counts{c} * decay;
        end
    end
    Q.last_t = t;
end

user_weights = get_user_weights(user_groups, C);
rate_ref_bps = max(params.bandwidth_Hz * log2(1 + params.P_W / max(params.sigma2_W, eps)), 1);

alpha_t = zeros(S, C);
action_t = zeros(C, 1);
selected_wdop = zeros(C, 1);
user_comm_norm = zeros(C, 1);
user_pos_util = zeros(C, 1);
user_penalty = zeros(C, 1);
user_reward_soft = zeros(C, 1);

for c = 1:C
    Kc = size(arm_bank(c).arms, 1);
    if Kc == 0
        error("select_satellites_mab_wdop_dynamic_ucb:NoArms", "UE%d 没有可选臂。", c);
    end

    counts = max(N_counts{c}, eps);
    q1 = Q.q1{c};
    q2 = Q.q2{c};
    exploration = cUcb * sqrt(log(max(t, 2)) ./ counts);
    ucb1 = q1 + exploration;
    ucb2 = q2 + exploration;

    wdopVals = arm_bank(c).wdop;
    penaltyVals = arm_bank(c).penalty;

    if useParetoUCB
        best_a = select_pareto_arm(ucb1, ucb2, paretoAlpha, penaltyVals, wdopPenaltyLambda);
    else
        selection_score = ucb1 - wdopPenaltyLambda * penaltyVals;
        [~, best_a] = max(selection_score);
    end

    action_t(c) = best_a;
    alpha_t(arm_bank(c).arms(best_a, :), c) = 1;
    selected_wdop(c) = wdopVals(best_a);
end

bf_t = ican.solve_beamforming_dc(params, chan, alpha_t);
R_c_bps = bf_t.R_c_bps(:);

for c = 1:C
    a = action_t(c);
    old_count = N_counts{c}(a);
    new_count = rho * old_count + 1;

    comm_norm = user_weights(c) * max(R_c_bps(c), 0) / rate_ref_bps;
    pos_util = 1 / (1 + selected_wdop(c));
    penalty = max(0, selected_wdop(c) - wdopThreshold)^2;
    soft_reward = user_weights(c) * max(R_c_bps(c), 0) / 1e9 - wdopPenaltyLambda * penalty;

    Q.q1{c}(a) = (rho * old_count * Q.q1{c}(a) + comm_norm) / new_count;
    Q.q2{c}(a) = (rho * old_count * Q.q2{c}(a) + pos_util) / new_count;
    N_counts{c}(a) = new_count;

    user_comm_norm(c) = comm_norm;
    user_pos_util(c) = pos_util;
    user_penalty(c) = penalty;
    user_reward_soft(c) = soft_reward;
end

total_reward = sum(user_reward_soft);
if total_reward > best_U || isempty(best_alpha)
    best_U = total_reward;
    best_alpha = alpha_t;
end

out = struct();
out.alpha = alpha_t;
out.bf = bf_t;
out.best_alpha = best_alpha;
out.best_utility = best_U;
out.utility = total_reward;
out.utility_bps = sum(user_weights .* max(R_c_bps, 0));
out.reward_per_user = user_reward_soft;
out.comm_reward_norm = user_comm_norm;
out.position_reward = user_pos_util;
out.penalty_per_user = user_penalty;
out.wdop_per_user = selected_wdop;
out.selected_arm_idx = action_t;
out.ucb_mode = ternary(useParetoUCB, "pareto", "ucb");
out.arm_bank = arm_bank;

if isempty(t_last_served) || numel(t_last_served) ~= C
    t_last_served = zeros(C, 1);
end
t_last_served(:) = t;

end

function arm_bank = build_candidate_arms(params, scenario, chan, comb, C, maxArms, wdopThreshold, wdopPenaltyLambda)
arm_bank = repmat(struct("arms", [], "wdop", [], "penalty", []), C, 1);

for c = 1:C
    wdopVals = zeros(size(comb, 1), 1);
    for k = 1:size(comb, 1)
        sats = comb(k, :);
        d_vec = chan.d_m(sats, c);
        wdopVals(k) = ican.compute_wdop(scenario.pUE(c, :), scenario.pSat(sats, :), d_vec);
    end

    if wdopPenaltyLambda <= 0
        feasibleMask = wdopVals <= wdopThreshold;
        if any(feasibleMask)
            feasibleComb = comb(feasibleMask, :);
            feasibleWdop = wdopVals(feasibleMask);
        else
            [bestWdop, ~] = min(wdopVals);
            ican.logf(params, "warn", ...
                "UE%d：没有臂满足 WDOP 阈值=%.3f，改用当前最优臂（最佳 WDOP=%.3f）。", ...
                c, wdopThreshold, bestWdop);
            feasibleComb = comb;
            feasibleWdop = wdopVals;
        end
    else
        feasibleComb = comb;
        feasibleWdop = wdopVals;
    end

    [sortedWdop, order] = sort(feasibleWdop, "ascend");
    keepCount = min(maxArms, numel(order));
    keepIdx = order(1:keepCount);
    chosenComb = feasibleComb(keepIdx, :);
    chosenWdop = sortedWdop(1:keepCount);
    chosenPenalty = max(0, chosenWdop - wdopThreshold).^2;

    arm_bank(c).arms = chosenComb;
    arm_bank(c).wdop = chosenWdop(:);
    arm_bank(c).penalty = chosenPenalty(:);
end
end

function [Q, N_counts] = align_state_to_arms(Q, N_counts, arm_bank, C)
for c = 1:C
    newArms = arm_bank(c).arms;
    Kc = size(newArms, 1);
    newQ1 = zeros(Kc, 1);
    newQ2 = zeros(Kc, 1);
    newN = zeros(Kc, 1);

    hasOldState = numel(Q.arms) >= c && ~isempty(Q.arms{c}) && ...
        numel(Q.q1) >= c && ~isempty(Q.q1{c}) && ...
        numel(Q.q2) >= c && ~isempty(Q.q2{c}) && ...
        numel(N_counts) >= c && ~isempty(N_counts{c});

    if hasOldState
        oldArms = Q.arms{c};
        oldQ1 = Q.q1{c};
        oldQ2 = Q.q2{c};
        oldN = N_counts{c};
        [isMatch, oldIdx] = ismember(newArms, oldArms, "rows");
        newQ1(isMatch) = oldQ1(oldIdx(isMatch));
        newQ2(isMatch) = oldQ2(oldIdx(isMatch));
        newN(isMatch) = oldN(oldIdx(isMatch));
    end

    Q.arms{c} = newArms;
    Q.q1{c} = newQ1;
    Q.q2{c} = newQ2;
    N_counts{c} = newN;
end
end

function best_a = select_pareto_arm(ucb1, ucb2, paretoAlpha, penaltyVals, wdopPenaltyLambda)
paretoIdx = pareto_front_indices(ucb1, ucb2);
if numel(paretoIdx) == 1
    best_a = paretoIdx;
    return;
end

score1 = normalize_minmax(ucb1(paretoIdx));
score2 = normalize_minmax(ucb2(paretoIdx));
penalty = penaltyVals(paretoIdx);
score = paretoAlpha * score1 + (1 - paretoAlpha) * score2 - wdopPenaltyLambda * penalty;
[~, loc] = max(score);
best_a = paretoIdx(loc);
end

function idx = pareto_front_indices(x, y)
n = numel(x);
isPareto = true(n, 1);
for i = 1:n
    if ~isPareto(i)
        continue;
    end
    dominated = (x >= x(i) & y >= y(i)) & ((x > x(i)) | (y > y(i)));
    isPareto(dominated) = false;
end
idx = find(isPareto);
end

function y = normalize_minmax(x)
x = x(:);
minX = min(x);
maxX = max(x);
if maxX > minX
    y = (x - minX) / (maxX - minX);
else
    y = zeros(size(x));
end
end

function user_weights = get_user_weights(user_groups, C)
user_weights = ones(C, 1);

if isfield(user_groups, "group_ids") && ~isempty(user_groups.group_ids)
    if isfield(user_groups, "weights") && ~isempty(user_groups.weights)
        weights = user_groups.weights(:);
    else
        weights = ones(numel(user_groups.group_ids), 1);
    end

    for k = 1:numel(user_groups.group_ids)
        ids = user_groups.group_ids{k};
        if isempty(ids)
            continue;
        end
        user_weights(ids(:)) = weights(min(k, numel(weights)));
    end
    return;
end

if isfield(user_groups, "URLLC") && isfield(user_groups, "eMBB") && isfield(user_groups, "mMTC")
    if isfield(user_groups.URLLC, "user_ids")
        user_weights(user_groups.URLLC.user_ids(:)) = get_field_or_default(user_groups.URLLC, "weight", 1.0);
    end
    if isfield(user_groups.eMBB, "user_ids")
        user_weights(user_groups.eMBB.user_ids(:)) = get_field_or_default(user_groups.eMBB, "weight", 1.0);
    end
    if isfield(user_groups.mMTC, "user_ids")
        user_weights(user_groups.mMTC.user_ids(:)) = get_field_or_default(user_groups.mMTC, "weight", 1.0);
    end
end
end

function value = get_param_value(params, fieldName, defaultValue)
value = defaultValue;
if isfield(params, fieldName)
    value = params.(fieldName);
end
end

function value = get_field_or_default(s, fieldName, defaultValue)
value = defaultValue;
if isfield(s, fieldName)
    value = s.(fieldName);
end
end

function user_groups = default_user_groups(C)
user_groups = struct();
user_groups.group_ids = {(1:C).'};
user_groups.weights = 1;
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
