function wdop = compute_wdop(pUE_c, pSat_sel, d_vec)
% COMPUTE_WDOP  与 Python 实现一致的 WDOP 计算
% 构造 A 矩阵 (I x 4)：每行为 [-unit_vec, 1]（位置 + 钟差），权重按 1/d^2

pc = pUE_c(:);
Ssel = size(pSat_sel, 1);

if Ssel < 4
    wdop = inf;
    return;
end

if numel(d_vec) ~= Ssel
    error("compute_wdop:SizeMismatch", "d_vec长度(%d)与卫星数(%d)不匹配.", numel(d_vec), Ssel);
end

% 构建 A 矩阵 (I x 4)
A = zeros(Ssel, 4);
dist_list = zeros(Ssel, 1);
for k = 1:Ssel
    ps = pSat_sel(k, :).';
    diff = ps - pc; % from UE to satellite
    d = norm(diff, 2);
    if d <= 1e-12
        wdop = inf;
        return;
    end
    dist_list(k) = d;
    A(k, 1:3) = (-diff ./ d).';
    A(k, 4) = 1.0;
end

% 权重：按距离平方倒数
w = 1.0 ./ max(dist_list(:).^2, 1e-12);

% 归一化权重（按均值），避免数值尺度问题
mean_w = mean(w);
if mean_w <= 0
    wdop = inf;
    return;
end
w = w / mean_w;

% 构造加权信息矩阵 H = A' W A
W = diag(w);
H = A.' * W * A;

% 数值保护
if rcond(H) < 1e-12
    wdop = inf;
    return;
end

wdop = sqrt(trace(inv(H)));
end
