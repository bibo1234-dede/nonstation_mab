function bf = solve_beamforming_zf(params, chan, alpha)
% Baseline 2: Zero-forcing (ZF) beamforming algorithm (Eq. 19 in paper)

    S = params.S; C = params.C; N = params.N;
    w = complex(zeros(N, C, S));
    R_sc_bps = zeros(S, C);
    P = params.P_W;
    sigma2 = params.sigma2_W;
    B = params.bandwidth_Hz;

    for s = 1:S
        served = find(alpha(s, :) > 0.5);
        K = length(served);
        if K == 0, continue; end

        % 1. 构建信道矩阵 H_s (大小为 K x N)
        H_s = zeros(K, N);
        for k = 1:K
            H_s(k, :) = chan.h(:, served(k), s)'; 
        end

        % 2. 计算伪逆 H_s_dagger
        H_s_dagger = pinv(H_s); 

        % 3. 计算功率缩放因子 beta
        beta = sqrt((P * K) / (norm(H_s_dagger, 'fro')^2));

        % 4. 获得波束矩阵 W_s
        W_s = beta * H_s_dagger;

        for k = 1:K
            c = served(k);
            w(:, c, s) = W_s(:, k);
        end

        % 5. 计算 SINR 和速率 (理论上 ZF 的 interf 接近 0)
        for k = 1:K
            c = served(k);
            hi = chan.h(:, c, s);
            desired = abs(hi' * w(:, c, s))^2;
            interf = 0;
            for j = 1:K
                if j == k, continue; end
                interf = interf + abs(hi' * w(:, served(j), s))^2;
            end
            R_sc_bps(s, c) = B * log2(1 + desired / (interf + sigma2));
        end
    end

    bf.w = w;
    bf.R_sc_bps = R_sc_bps;
    bf.R_c_bps = sum(R_sc_bps, 1).';
    bf.sumRate_bps = sum(bf.R_c_bps);
end