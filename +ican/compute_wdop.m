function wdop = compute_wdop(pUE_c, pSat_sel, d_vec)
% compute_wdop  基于路径损耗加权的WDOP计算（归一化版本）
%
% 输入:
%   pUE_c    - 用户位置 [1x3]
%   pSat_sel - 所选卫星位置 [Ssel x 3]
%   d_vec    - 所选卫星到该用户的距离 [Ssel x 1]，单位: 米
%
% 输出:
%   wdop     - 加权几何精度衰减因子（无量纲，量级与GDOP相同）
%
% 权重定义:
%   sigma_i = d_{s,c} / d_ref，其中 d_ref = 600km（LEO轨道高度）
%   权重 W_ii = 1 / sigma_i = d_ref / d_{s,c}
%   距离越远，权重越小（信号质量越差）
%
% 对应论文公式:
%   G' = W * G
%   WDOP = sqrt(trace((G'^T G')^{-1}))

pc = pUE_c(:);
Ssel = size(pSat_sel, 1);

if Ssel < 3
    wdop = inf;
    return;
end

if numel(d_vec) ~= Ssel
    error("compute_wdop:SizeMismatch", ...
        "d_vec长度(%d)与卫星数(%d)不匹配.", numel(d_vec), Ssel);
end

% 构建几何矩阵 G (Ssel x 3)
G = zeros(Ssel, 3);
for k = 1:Ssel
    ps = pSat_sel(k, :).';
    d = norm(pc - ps, 2);
    G(k, :) = ((pc - ps) ./ d).';
end

% 构建权重矩阵 W = diag(1/sigma_i)
% sigma_i = d_{s,c} / d_ref，归一化到合理的量级
d_ref = 600e3;  % 参考距离（600 km，LEO轨道高度）
sigma_vec = d_vec(:) / d_ref;  % 归一化权重，量纲为1，范围约 [1.0, 1.5]
sigma_vec = max(sigma_vec, 0.1);  % 防止异常值和除零
W = diag(1 ./ sigma_vec);        % W_ii = d_ref / d_{s,c}

% 加权几何矩阵
G_prime = W * G;

% 计算WDOP
GTG = G_prime.' * G_prime;
if rcond(GTG) < 1e-12
    wdop = inf;
    return;
end

wdop = sqrt(trace(inv(GTG)));
end