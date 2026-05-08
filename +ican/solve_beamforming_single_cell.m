function bf = solve_beamforming_single_cell(params, chan, alpha)
% Baseline 1: Single-cell beamforming algorithm (Eq. 18 in paper)

    S = params.S; C = params.C; N = params.N;
    w = complex(zeros(N, C, S));
    R_sc_bps = zeros(S, C);
    P = params.P_W;
    sigma2 = params.sigma2_W;
    B = params.bandwidth_Hz;

    for s = 1:S
        served = find(alpha(s, :) > 0.5);
        if isempty(served), continue; end
        
        % 1. 计算波束权重 w
        for k = 1:length(served)
            c = served(k);
            h_sc = chan.h(:, c, s);
            % 根据公式 (18) 分配功率
            w(:, c, s) = sqrt(P / (norm(h_sc, 'fro')^2)) * h_sc;
        end
        
        % 2. 计算 SINR 和速率
        for k = 1:length(served)
            c = served(k);
            hi = chan.h(:, c, s);
            desired = abs(hi' * w(:, c, s))^2;
            interf = 0;
            for j = 1:length(served)
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