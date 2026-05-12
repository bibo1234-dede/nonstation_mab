function [out, Q, N_counts, t_last_served] = select_satellites_mab_wdop_dynamic_ucb( ...
    params, scenario, chan, Q, N_counts, best_U, best_alpha, t, t_last_served, user_groups)
% SELECT_SATELLITES_MAB_WDOP_DYNAMIC_UCB  面向 WDOP 的非稳态 MAB 选星器。

S = params.S;
C = params.C;
I = params.I;

if nargin < 10 || isempty(user_groups)
    user_groups = default_user_groups(C);
end

useGrouping = get_param_value(params, "useGrouping", false);
rho = get_param_value(params, "mabRho", 0.98);
cUcb = get_param_value(params, "mabCucb", 1.0);
maxArms = get_param_value(params, "mabMaxArms", 20);
wdopThreshold = get_param_value(params, "wdopThreshold", 6.0);
wdopPenaltyLambda = get_param_value(params, "wdopPenaltyLambda", 0.5);
wdopSoftMargin = get_param_value(params, "wdopSoftMargin", 1.5);
satLoadCap = get_param_value(params, "satLoadCap", 4);
satLoadPenaltyLambda = get_param_value(params, "satLoadPenaltyLambda", 0.5);
groupReuseBonusLambda = get_param_value(params, "groupReuseBonusLambda", 0.0);
userReusePenaltyLambda = get_param_value(params, "userReusePenaltyLambda", get_param_value(params, "interGroupReusePenaltyLambda", 0.0));
useParetoUCB = get_param_value(params, "useParetoUCB", false);
paretoAlpha = get_param_value(params, "paretoAlpha", 0.5);
candidatePoolSize = get_param_value(params, "candidatePoolSize", maxArms);
candidateRateWeight = get_param_value(params, "candidateRateWeight", 1.0);
candidateWdopWeight = get_param_value(params, "candidateWdopWeight", 0.25);
debugPrintCandidateArms = get_param_value(params, "debugPrintCandidateArms", false);
debugPrintSelection = get_param_value(params, "debugPrintSelection", false);
debugPrintCandidateLimit = get_param_value(params, "debugPrintCandidateLimit", 12);

if ~useGrouping
    groupReuseBonusLambda = 0;
end

comb = nchoosek(1:S, I);
arm_bank = build_candidate_arms(params, scenario, chan, comb, C, maxArms, candidatePoolSize, candidateRateWeight, candidateWdopWeight, wdopThreshold);

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
if useGrouping
    user_order = get_user_order(user_groups, C);
else
    user_order = mod((0:C-1).' + (t - 1), C) + 1;
end
sat_load = zeros(S, 1);
if useGrouping && isfield(user_groups, "group_ids") && ~isempty(user_groups.group_ids)
    group_sat_load = zeros(max(numel(user_groups.group_ids), 1), S);
else
    group_sat_load = zeros(1, S);
end

alpha_t = zeros(S, C);
action_t = zeros(C, 1);
selected_wdop = zeros(C, 1);
selected_load_penalty = zeros(C, 1);
selected_reuse_penalty = zeros(C, 1);
user_comm_norm = zeros(C, 1);
user_pos_util = zeros(C, 1);
user_penalty = zeros(C, 1);
user_reward_soft = zeros(C, 1);

for orderIdx = 1:numel(user_order)
    c = user_order(orderIdx);
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
    [loadPenaltyVals, feasibleMask] = compute_load_penalty(arm_bank(c).arms, sat_load, satLoadCap);
    if any(~feasibleMask)
        loadPenaltyVals(~feasibleMask) = inf;
    end

    groupBonusVals = zeros(Kc, 1);
    if useGrouping && groupReuseBonusLambda > 0
        groupIdx = get_user_group_index(user_groups, c);
        if groupIdx >= 1 && groupIdx <= size(group_sat_load, 1)
            groupBonusVals = compute_group_bonus(arm_bank(c).arms, group_sat_load(groupIdx, :)', satLoadCap);
        end
    else
        groupIdx = 1;
    end

    reusePenaltyVals = compute_reuse_penalty(arm_bank(c).arms, sat_load, satLoadCap);

    if debugPrintCandidateArms
        print_candidate_arms(params, c, arm_bank(c).arms, arm_bank(c).wdop, arm_bank(c).rateProxy, arm_bank(c).candidateScore, debugPrintCandidateLimit);
    end

    if useParetoUCB
        best_a = select_pareto_arm(ucb1, ucb2, paretoAlpha, penaltyVals, wdopPenaltyLambda, loadPenaltyVals, groupBonusVals, reusePenaltyVals, satLoadPenaltyLambda, groupReuseBonusLambda, userReusePenaltyLambda);
    else
        selection_score = ucb1 - wdopPenaltyLambda * penaltyVals - satLoadPenaltyLambda * loadPenaltyVals - userReusePenaltyLambda * reusePenaltyVals + groupReuseBonusLambda * groupBonusVals;
        if all(~isfinite(selection_score))
            [~, best_a] = min(loadPenaltyVals);
            ican.logf(params, "warn", "UE%d 的所有候选臂都将超过卫星负载上限 %d，改选最小负载惩罚臂。", c, satLoadCap);
        else
            selection_score(~isfinite(selection_score)) = -inf;
            [~, best_a] = max(selection_score);
        end
    end

    action_t(c) = best_a;
    chosen_sats = arm_bank(c).arms(best_a, :);
    alpha_t(chosen_sats, c) = 1;
    selected_wdop(c) = wdopVals(best_a);
    selected_load_penalty(c) = loadPenaltyVals(best_a);
    selected_reuse_penalty(c) = reusePenaltyVals(best_a);
    selected_rate_proxy = arm_bank(c).rateProxy(best_a);
    selected_rate_proxy = selected_rate_proxy(1);

    sat_load(chosen_sats) = sat_load(chosen_sats) + 1;
    if useGrouping && groupIdx >= 1 && groupIdx <= size(group_sat_load, 1)
        group_sat_load(groupIdx, chosen_sats) = group_sat_load(groupIdx, chosen_sats) + 1;
    end

    if debugPrintSelection
        ican.logf(params, "info", "UE%d 选中臂=%s | WDOP=%.3f | rateProxy=%.6g | score=%.6g", ...
            c, mat2str(chosen_sats), selected_wdop(c), selected_rate_proxy, arm_bank(c).candidateScore(best_a));
    end
end

bf_t = ican.solve_beamforming_dc(params, chan, alpha_t);
R_c_bps = bf_t.R_c_bps(:);

for c = 1:C
    a = action_t(c);
    old_count = N_counts{c}(a);
    new_count = rho * old_count + 1;

    comm_norm = user_weights(c) * max(R_c_bps(c), 0) / rate_ref_bps;
    comm_norm = comm_norm(1);
    pos_util = 1 / (1 + selected_wdop(c));
    pos_util = pos_util(1);
    penalty = compute_wdop_soft_penalty(selected_wdop(c), wdopThreshold, wdopSoftMargin);
    penalty = penalty(1);
    soft_reward = user_weights(c) * max(R_c_bps(c), 0) / 1e9 - wdopPenaltyLambda * penalty - satLoadPenaltyLambda * selected_load_penalty(c) - userReusePenaltyLambda * selected_reuse_penalty(c);
    soft_reward = soft_reward(1);

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
out.utility_bps = sum(user_weights .* max(R_c_bps, 0)); %#ok<NASGU>
out.reward_per_user = user_reward_soft;
out.comm_reward_norm = user_comm_norm;
out.position_reward = user_pos_util;
out.penalty_per_user = user_penalty;
out.wdop_per_user = selected_wdop;
out.load_penalty_per_user = selected_load_penalty;
out.sat_load = sat_load;
out.group_sat_load = group_sat_load;
out.selected_arm_idx = action_t;
out.ucb_mode = ternary(useParetoUCB, "pareto", "ucb");
out.arm_bank = arm_bank;

if isempty(t_last_served) || numel(t_last_served) ~= C
    t_last_served = zeros(C, 1);
end
t_last_served(:) = t;

end

function arm_bank = build_candidate_arms(params, scenario, chan, comb, C, maxArms, candidatePoolSize, candidateRateWeight, candidateWdopWeight, wdopThreshold)
arm_bank = repmat(struct("arms", [], "wdop", [], "penalty", [], "rateProxy", [], "candidateScore", []), C, 1);

for c = 1:C
    wdopVals = zeros(size(comb, 1), 1);
    rateProxyVals = zeros(size(comb, 1), 1);
    for k = 1:size(comb, 1)
        sats = comb(k, :);
        d_vec = chan.d_m(sats, c);
        wdopVals(k) = ican.compute_wdop(scenario.pUE(c, :), scenario.pSat(sats, :), d_vec);
        h_sel = chan.h(:, c, sats);
        rateProxyVals(k) = sum(abs(h_sel).^2, "all");
    end

    feasibleMask = wdopVals <= wdopThreshold;
    if ~any(feasibleMask)
        error("select_satellites_mab_wdop_dynamic_ucb:NoWdopFeasibleArm", ...
            "UE%d 没有任何 WDOP <= %.3f 的可行臂。当前最小 WDOP=%.3f。", ...
            c, wdopThreshold, min(wdopVals));
    end

    feasibleComb = comb(feasibleMask, :);
    feasibleWdop = wdopVals(feasibleMask);

    rateScore = normalize_minmax(rateProxyVals(feasibleMask));
    wdopScore = normalize_minmax(feasibleWdop);
    candidateScore = candidateRateWeight * rateScore - candidateWdopWeight * wdopScore;
    [~, order] = sort(candidateScore, "descend");
    keepCount = min(max(candidatePoolSize, maxArms), numel(order));
    keepIdx = order(1:keepCount);
    chosenComb = feasibleComb(keepIdx, :);
    chosenWdop = feasibleWdop(keepIdx);
    chosenRateProxy = rateProxyVals(feasibleMask);
    chosenRateProxy = chosenRateProxy(keepIdx);
    chosenPenalty = compute_wdop_soft_penalty(chosenWdop, wdopThreshold, get_param_value(params, "wdopSoftMargin", 1.5));
    chosenScore = candidateScore(keepIdx);

    arm_bank(c).arms = chosenComb;
    arm_bank(c).wdop = chosenWdop(:);
    arm_bank(c).penalty = chosenPenalty(:);
    arm_bank(c).rateProxy = chosenRateProxy(:);
    arm_bank(c).candidateScore = chosenScore(:);
end
end

function print_candidate_arms(params, c, arms, wdopVals, rateProxyVals, candidateScore, limit)
if isempty(arms)
    ican.logf(params, "info", "UE%d 候选臂为空。", c);
    return;
end

limit = max(1, min(limit, size(arms, 1)));
ican.logf(params, "info", "UE%d 候选臂（前 %d/%d 个）:", c, limit, size(arms, 1));
for k = 1:limit
    ican.logf(params, "info", "  arm#%d sats=%s | WDOP=%.3f | rateProxy=%.6g | score=%.6g", ...
        k, mat2str(arms(k, :)), wdopVals(k), rateProxyVals(k), candidateScore(k));
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

function best_a = select_pareto_arm(ucb1, ucb2, paretoAlpha, penaltyVals, wdopPenaltyLambda, loadPenaltyVals, groupBonusVals, reusePenaltyVals, satLoadPenaltyLambda, groupReuseBonusLambda, userReusePenaltyLambda)
paretoIdx = pareto_front_indices(ucb1, ucb2);
if numel(paretoIdx) == 1
    best_a = paretoIdx;
    return;
end

score1 = normalize_minmax(ucb1(paretoIdx));
score2 = normalize_minmax(ucb2(paretoIdx));
penalty = penaltyVals(paretoIdx);
loadPenalty = loadPenaltyVals(paretoIdx);
groupBonus = groupBonusVals(paretoIdx);
reusePenalty = reusePenaltyVals(paretoIdx);
score = paretoAlpha * score1 + (1 - paretoAlpha) * score2 - wdopPenaltyLambda * penalty - satLoadPenaltyLambda * loadPenalty - userReusePenaltyLambda * reusePenalty + groupReuseBonusLambda * groupBonus;
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

function [loadPenaltyVals, feasibleMask] = compute_load_penalty(arms, loadVec, satLoadCap)
if isempty(arms)
    loadPenaltyVals = zeros(0, 1);
    feasibleMask = true(0, 1);
    return;
end

satLoadCap = max(satLoadCap, 1);
armLoad = loadVec(arms);
projectedLoad = armLoad + 1;
feasibleMask = all(projectedLoad <= satLoadCap, 2);
loadPenaltyVals = mean(projectedLoad / satLoadCap, 2);
end

function groupBonusVals = compute_group_bonus(arms, loadVec, satLoadCap)
if isempty(arms)
    groupBonusVals = zeros(0, 1);
    return;
end

satLoadCap = max(satLoadCap, 1);
armLoad = loadVec(arms);
groupBonusVals = mean(armLoad / satLoadCap, 2);
end

function reusePenaltyVals = compute_reuse_penalty(arms, loadVec, satLoadCap)
if isempty(arms)
    reusePenaltyVals = zeros(0, 1);
    return;
end

satLoadCap = max(satLoadCap, 1);
otherLoad = loadVec(arms);
reusePenaltyVals = mean(max(otherLoad, 0) / satLoadCap, 2);
end

function penalty = compute_wdop_soft_penalty(wdopVals, wdopThreshold, wdopSoftMargin)
wdopSoftMargin = max(wdopSoftMargin, eps);
excess = max(0, wdopVals - wdopThreshold);
penalty = excess / wdopSoftMargin;
end

function groupIdx = get_user_group_index(user_groups, c)
groupIdx = 1;
if isfield(user_groups, "group_assignment") && numel(user_groups.group_assignment) >= c
    groupIdx = user_groups.group_assignment(c);
    if ~isfinite(groupIdx) || groupIdx < 1
        groupIdx = 1;
    end
    return;
end

if isfield(user_groups, "group_ids") && ~isempty(user_groups.group_ids)
    for k = 1:numel(user_groups.group_ids)
        ids = user_groups.group_ids{k};
        if any(ids(:) == c)
            groupIdx = k;
            return;
        end
    end
end
end

function user_order = get_user_order(user_groups, C)
user_order = (1:C).';
if ~isfield(user_groups, "group_ids") || isempty(user_groups.group_ids)
    return;
end

ordered = [];
for k = 1:numel(user_groups.group_ids)
    ids = user_groups.group_ids{k};
    ids = ids(ids >= 1 & ids <= C);
    ordered = [ordered; ids(:)]; %#ok<AGROW>
end

if numel(unique(ordered, "stable")) == C
    user_order = unique(ordered, "stable");
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

end

function value = get_param_value(params, fieldName, defaultValue)
value = defaultValue;
if isfield(params, fieldName)
    value = params.(fieldName);
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
