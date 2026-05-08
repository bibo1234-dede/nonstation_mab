function [out, Q, N_counts, t_last_served] = select_satellites_mab_wdop_dynamic_ucb( ...
    params, scenario, chan, Q, N_counts, best_U, best_alpha, t, t_last_served, user_groups)
% select_satellites_mab_wdop_dynamic_ucb  Non-stationary D-UCB satellite selection.
%
% The feasible arm set is rebuilt at every time step from the current WDOP.
% MAB state is then remapped by the actual satellite-combination rows so Q
% and counts stay attached to the same arm even when the sorted feasible list
% changes. Rewards are kept in log-rate units consistently.

S = params.S;
C = params.C;
I = params.I;

comb = nchoosek(1:S, I);
maxArms = 10;
cUcb = 1.0;
rho = 0.98;

prefComb = build_feasible_arms(params, scenario, chan, comb, C, maxArms);

if isempty(Q) || ~isstruct(Q)
    Q = struct();
    Q.values = cell(C, 1);
    Q.arms = cell(C, 1);
    Q.last_t = max(t - 1, 0);
end
if ~isfield(Q, "values")
    Q.values = cell(C, 1);
end
if ~isfield(Q, "arms")
    Q.arms = cell(C, 1);
end
if ~isfield(Q, "last_t")
    Q.last_t = max(t - 1, 0);
end

if isempty(N_counts) || ~iscell(N_counts)
    N_counts = cell(C, 1);
end

[Q, N_counts] = remap_mab_state(Q, N_counts, prefComb, C);

% D-UCB discounting: decay all historical statistics once per time step.
% This avoids the previous no-op discount caused by refreshing one
% t_last_served(c) value for every user at every step.
elapsed = max(0, t - Q.last_t);
if elapsed > 0
    decay = rho^elapsed;
    for c = 1:C
        Q.values{c} = Q.values{c} * decay;
        N_counts{c} = N_counts{c} * decay;
    end
    Q.last_t = t;
end

alpha_t = zeros(S, C);
action_t = zeros(C, 1);
reward_per_user = zeros(C, 1);

for c = 1:C
    Kc = size(prefComb{c}, 1);
    ucb_values = zeros(Kc, 1);

    for a = 1:Kc
        if N_counts{c}(a) <= eps
            ucb_values(a) = inf;
        else
            exploration = cUcb * sqrt(log(max(t, 2)) / N_counts{c}(a));
            ucb_values(a) = Q.values{c}(a) + exploration;
        end
    end

    [~, best_a] = max(ucb_values);
    action_t(c) = best_a;
    alpha_t(prefComb{c}(best_a, :), c) = 1;
end

bf_t = ican.solve_beamforming_dc(params, chan, alpha_t);
R_c_bps = bf_t.R_c_bps;

for c = 1:C
    w = qos_weight(c, user_groups);
    reward_per_user(c) = w * R_c_bps(c);  % 加权速率，与波束成形目标一致

    a = action_t(c);
    old_count = N_counts{c}(a);
    new_count = old_count + 1;
    Q.values{c}(a) = (Q.values{c}(a) * old_count + reward_per_user(c)) / new_count;
    N_counts{c}(a) = new_count;
end

total_reward = sum(reward_per_user);
if total_reward > best_U || isempty(best_alpha)
    best_U = total_reward;
    best_alpha = alpha_t;
end

out = struct();
out.alpha = alpha_t;
out.utility_bps = total_reward;
out.best_alpha = best_alpha;
out.best_utility = best_U;
out.reward_per_user = reward_per_user;
out.total_reward = total_reward;

% Keep this output meaningful for callers/debug plots without using it for
% discounting. It records when each UE was last updated under the current run.
if isempty(t_last_served) || numel(t_last_served) ~= C
    t_last_served = zeros(C, 1);
end
t_last_served(:) = t;

end

function prefComb = build_feasible_arms(params, scenario, chan, comb, C, maxArms)
prefComb = cell(C, 1);

for c = 1:C
    wdopVals = zeros(size(comb, 1), 1);
    for k = 1:size(comb, 1)
        sats = comb(k, :);
        d_vec = chan.d_m(sats, c);
        wdopVals(k) = ican.compute_wdop( ...
            scenario.pUE(c, :), scenario.pSat(sats, :), d_vec);
    end

    feasibleMask = wdopVals <= params.wdopThreshold;
    if any(feasibleMask)
        feasibleComb = comb(feasibleMask, :);
        feasibleWdop = wdopVals(feasibleMask);
    else
        [bestWdop, ~] = min(wdopVals);
        ican.logf(params, "warn", ...
            "UE%d: no MAB arm satisfies WDOP threshold=%.3f; using best available arms (best WDOP=%.3f).", ...
            c, params.wdopThreshold, bestWdop);
        feasibleComb = comb;
        feasibleWdop = wdopVals;
    end

    [~, order] = sort(feasibleWdop, "ascend");
    keep_idx = order(1:min(maxArms, numel(order)));
    prefComb{c} = feasibleComb(keep_idx, :);
end
end

function [Q, N_counts] = remap_mab_state(Q, N_counts, prefComb, C)
for c = 1:C
    newArms = prefComb{c};
    Kc = size(newArms, 1);
    newQ = zeros(Kc, 1);
    newN = zeros(Kc, 1);

    hasOldState = numel(Q.arms) >= c && ~isempty(Q.arms{c}) && ...
        numel(Q.values) >= c && ~isempty(Q.values{c}) && ...
        numel(N_counts) >= c && ~isempty(N_counts{c});

    if hasOldState
        oldArms = Q.arms{c};
        oldQ = Q.values{c};
        oldN = N_counts{c};
        [isOldArm, oldIdx] = ismember(newArms, oldArms, "rows");
        newQ(isOldArm) = oldQ(oldIdx(isOldArm));
        newN(isOldArm) = oldN(oldIdx(isOldArm));
    end

    Q.arms{c} = newArms;
    Q.values{c} = newQ;
    N_counts{c} = newN;
end
end

function w = qos_weight(c, user_groups)
if ismember(c, user_groups.URLLC.user_ids)
    w = 2.0;
elseif ismember(c, user_groups.eMBB.user_ids)
    w = 1.2;
else
    w = 0.6;
end
end
