function out = select_satellites_mab_wdop(params, scenario, chan)
% select_satellites_mab_wdop  基于WDOP约束的MAB-UCB选星算法
%
% 相比原 select_satellites_mab：
%   将 GDOP 约束替换为 WDOP 约束 (params.wdopThreshold)
%   权重 sigma_i = d_{s,c}，来自 chan.d_m
%   其余 UCB 逻辑、Q表、热启动完全不变

S = params.S; C = params.C; I = params.I;
comb = nchoosek(1:S, I);

ican.logf(params, "info", "MAB-WDOP init: Constructing feasible arms with WDOP constraint...");

% 1. 为每个用户构建候选臂 (WDOP约束替代GDOP约束)
prefComb = cell(C, 1);
prefWdop = cell(C, 1);
Max_Arms = 10;

for c = 1:C
    wdopVals = zeros(size(comb, 1), 1);
    for k = 1:size(comb, 1)
        sats = comb(k, :);
        d_vec = chan.d_m(sats, c);  % 关键：用距离作为权重来源
        wdopVals(k) = ican.compute_wdop( ...
            scenario.pUE(c, :), scenario.pSat(sats, :), d_vec);
    end

    feasibleMask = (wdopVals <= params.wdopThreshold);

    if ~any(feasibleMask)
        error("UE%d has no feasible satellite combination under WDOP threshold=%.3f.", ...
            c, params.wdopThreshold);
    end

    feasibleComb = comb(feasibleMask, :);
    feasibleWdop = wdopVals(feasibleMask);

    % 按WDOP升序排序，截断至前Max_Arms个
    [feasibleWdop, order] = sort(feasibleWdop, "ascend");
    keep_idx = 1:min(Max_Arms, length(feasibleWdop));
    prefComb{c} = feasibleComb(order(keep_idx), :);
    prefWdop{c} = feasibleWdop(keep_idx);
end

% 2. UCB参数初始化 (与原MAB完全相同)
Max_Iters = 25;
c_ucb = 1.0;

Q = cell(C, 1);
N_counts = cell(C, 1);
for c = 1:C
    Kc = size(prefComb{c}, 1);
    Q{c} = zeros(Kc, 1);
    N_counts{c} = zeros(Kc, 1);
end

best_alpha = zeros(S, C);
best_U = -inf;
best_bf = [];

% 初始状态：选WDOP最小的组合（一号臂）
alpha_curr = zeros(S, C);
for c = 1:C
    alpha_curr(prefComb{c}(1, :), c) = 1;
end

bf_init = ican.solve_beamforming_dc(params, chan, alpha_curr);
best_U = bf_init.sumRate_bps;
best_alpha = alpha_curr;
best_bf = bf_init;

for c = 1:C
    Q{c}(1) = best_U;
    N_counts{c}(1) = 1;
end

ican.logf(params, "info", "MAB-WDOP Start: Initial Sum Rate = %.3f Gbps", best_U/1e9);

% 3. UCB迭代探索 (逻辑与原MAB完全相同，仅臂空间基于WDOP筛选)
for t = 1:Max_Iters
    alpha_t = zeros(S, C);
    action_t = zeros(C, 1);

    for c = 1:C
        Kc = size(prefComb{c}, 1);
        ucb_values = zeros(Kc, 1);

        for a = 1:Kc
            if N_counts{c}(a) == 0
                ucb_values(a) = inf;
            else
                q_norm = Q{c}(a) / 1e9;
                explore_term = c_ucb * sqrt(log(t) / N_counts{c}(a));
                ucb_values(a) = q_norm + explore_term;
            end
        end

        [~, best_a] = max(ucb_values);
        action_t(c) = best_a;
        alpha_t(prefComb{c}(best_a, :), c) = 1;
    end

    bf_t = ican.solve_beamforming_dc(params, chan, alpha_t, ...
        "warmStart", best_bf.warmStart);
    reward_t = bf_t.sumRate_bps;

    for c = 1:C
        a = action_t(c);
        N_counts{c}(a) = N_counts{c}(a) + 1;
        Q{c}(a) = Q{c}(a) + (reward_t - Q{c}(a)) / N_counts{c}(a);
    end

    if reward_t > best_U
        best_U = reward_t;
        best_alpha = alpha_t;
        best_bf = bf_t;
        ican.logf(params, "info", ...
            "MAB-WDOP Iter %d: New Best! Sum Rate = %.3f Gbps", t, best_U/1e9);
    else
        ican.logf(params, "info", ...
            "MAB-WDOP Iter %d: Sum Rate = %.3f Gbps (Best = %.3f Gbps)", ...
            t, reward_t/1e9, best_U/1e9);
    end
end

% 4. 计算最终WDOP
wdop = zeros(C, 1);
for c = 1:C
    sats = find(best_alpha(:, c) > 0.5);
    d_vec = chan.d_m(sats, c);
    wdop(c) = ican.compute_wdop(scenario.pUE(c, :), scenario.pSat(sats, :), d_vec);
end

out = struct();
out.alpha = best_alpha;
out.bf = best_bf;
out.gdop = wdop;   % 兼容原字段名
out.wdop = wdop;
out.utility_bps = best_U;
end
