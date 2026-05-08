function sel = select_satellites_wdop(params, scenario, chan)
% select_satellites_wdop  基于WDOP最小化的贪婪选星算法
%
% 相比原 select_satellites_gdop：
%   将 compute_gdop 替换为 compute_wdop
%   权重 sigma_i = d_{s,c}，来自 chan.d_m
%   约束由 params.gdopThreshold 改为 params.wdopThreshold
%
% 输出结构体与原函数完全兼容，字段名保持一致

S = params.S;
C = params.C;
I = params.I;

comb = nchoosek(1:S, I);
nComb = size(comb, 1);

alpha = zeros(S, C);
satIdx = cell(C, 1);
wdop = zeros(C, 1);
wdopAll = zeros(C, nComb);

for c = 1:C
    pc = scenario.pUE(c, :);
    for k = 1:nComb
        sats = comb(k, :);
        d_vec = chan.d_m(sats, c); % 所选卫星到用户c的距离
        wdopAll(c, k) = ican.compute_wdop(pc, scenario.pSat(sats, :), d_vec);
    end

    [wdopMin, idx] = min(wdopAll(c, :));
    bestSats = comb(idx, :);

    if wdopMin > params.wdopThreshold
        ican.logf(params, "warn", ...
            "UE%d: best WDOP %.3f exceeds threshold=%.3f.", ...
            c, wdopMin, params.wdopThreshold);
    else
        ican.logf(params, "info", "UE%d: WDOP-based selection WDOP=%.3f sats=%s", ...
            c, wdopMin, mat2str(bestSats));
    end

    alpha(bestSats, c) = 1;
    satIdx{c} = bestSats;
    wdop(c) = wdopMin;
end

sel = struct();
sel.alpha = alpha;
sel.satIdx = satIdx;
sel.gdop = wdop;    % 保持字段名兼容，存的是WDOP值
sel.wdop = wdop;
sel.wdopAll = wdopAll;
end
