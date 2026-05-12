function sel = select_satellites_wdop(params, scenario, chan)
% SELECT_SATELLITES_WDOP  基于 WDOP 最小化的贪婪选星算法。
%
% 相比旧版 GDOP 选星：
%   1) 将几何指标替换为 WDOP
%   2) 权重来源改为 chan.d_m
%   3) 约束阈值改为 params.wdopThreshold
%
% 输出结构体保持兼容，字段名沿用旧接口。

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
sel.wdop = wdop;
sel.wdopAll = wdopAll;
end
