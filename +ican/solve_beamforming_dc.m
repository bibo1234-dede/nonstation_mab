function bf = solve_beamforming_dc(params, chan, alpha, varargin)
% SOLVE_BEAMFORMING_DC  使用 DC + CVX 求解多卫星波束赋形。

parser = inputParser();
parser.addParameter("warmStart", struct(), @(x) isstruct(x));
parser.addParameter("satIdx", [], @(x) isnumeric(x));
parser.parse(varargin{:});
warm = parser.Results.warmStart;
satIdx = parser.Results.satIdx;

S = params.S;
C = params.C;
N = params.N;

cvx_solver(char(params.cvx.solver));

if isempty(satIdx)
    satIdx = 1:S;
end
satIdx = unique(satIdx(:).');
if any(satIdx < 1) || any(satIdx > S)
    error("solve_beamforming_dc:BadSatIdx", "satIdx must be within [1,%d].", S);
end

w = complex(zeros(N, C, S));
R_sc_bps = zeros(S, C);

sat = repmat(struct( ...
    "servedUeIdx", [], ...
    "dcIters", 0, ...
    "sumRateFromQ_bps", 0, ...
    "cvxStatus", "", ...
    "cvxOptval", NaN), S, 1);

warmOut = struct();
warmOut.Q = cell(S, 1);
warmOut.ueIdx = cell(S, 1);

for s = satIdx
    served = find(alpha(s, :) > 0.5);
    sat(s).servedUeIdx = served;
    if isempty(served)
        continue;
    end

    h_s = chan.h(:, served, s); % N x K

    initQ = [];
    if isfield(warm, "Q") && isfield(warm, "ueIdx") ...
            && numel(warm.Q) >= s && numel(warm.ueIdx) >= s ...
            && ~isempty(warm.Q{s}) && ~isempty(warm.ueIdx{s})
        [initQ, initMapOk] = try_map_warmstart(params, warm.Q{s}, warm.ueIdx{s}, served);
        if ~initMapOk
            initQ = [];
        end
    end

    satRes = solve_dc_one_satellite(params, h_s, initQ);
    sat(s).dcIters = satRes.dcIters;
    sat(s).sumRateFromQ_bps = satRes.sumRateFromQ_bps;
    sat(s).cvxStatus = satRes.cvxStatus;
    sat(s).cvxOptval = satRes.cvxOptval;

    % Map back to global UE indices.
    w(:, served, s) = satRes.w;
    R_sc_bps(s, served) = satRes.R_bps;

    warmOut.Q{s} = satRes.Q;
    warmOut.ueIdx{s} = served;
end

R_c_bps = sum(R_sc_bps, 1).';

bf = struct();
bf.w = w;
bf.R_sc_bps = R_sc_bps;
bf.R_c_bps = R_c_bps;
bf.sumRate_bps = sum(R_c_bps);
bf.satSumRate_bps = sum(R_sc_bps, 2);
bf.warmStart = warmOut;
bf.sat = sat;
end

function satRes = solve_dc_one_satellite(params, h_s, initQ)
% 对单颗卫星执行算法1的求解过程。

B = params.bandwidth_Hz;
sigma2 = params.sigma2_W;
P = params.P_W;

[N, K] = size(h_s);

% 预计算 H_i = h_i h_i^H，便于写成迹形式。
Hi = cell(K, 1);
for i = 1:K
    hi = h_s(:, i);
    Hi{i} = hi * hi';
end

ticSat = tic;

Q_prev = init_random_Q(params, N, K);
if ~isempty(initQ)
    if isequal(size(initQ), [N, N, K])
        Q_prev = initQ;
    end
end

sumRate_prev = sum_rate_from_Q(B, sigma2, Hi, Q_prev);

cvxStatusLast = "";
cvxOptvalLast = NaN;

dcIters = 0;
for iter = 1:params.alg.maxDcIters
    dcIters = iter;

    denomPrev = zeros(K, 1);
    for i = 1:K
        denomPrev(i) = sigma2;
        for k = 1:K
            if k == i
                continue;
            end
            denomPrev(i) = denomPrev(i) + real(trace(Q_prev(:, :, k) * Hi{i}));
        end
        denomPrev(i) = max(denomPrev(i), sigma2); %#ok<AGROW>
    end

    tCvx = tic;
        [Q_new, cvxStatus, cvxOptval] = solve_cvx_dc_step(params, B, sigma2, P, Hi, Q_prev, denomPrev);
        cvxTime_s = toc(tCvx);
        cvxStatusLast = cvxStatus;
        cvxOptvalLast = cvxOptval;

        if ~contains(cvxStatus, "Solved")
            ican.logf(params, "warn", "CVX failed in DC step (status=%s). Falling back to previous Q and continuing. Set params.cvx.quiet=false to view solver output.", cvxStatus);
            % Fall back to previous feasible Q to keep the simulation running.
            Q_new = Q_prev;
            % Record status and continue (exit DC iterations early)
            cvxStatusLast = cvxStatus;
            cvxOptvalLast = cvxOptval;
            break;
        end

    sumRate_new = sum_rate_from_Q(B, sigma2, Hi, Q_new);
    diff_bps = abs(sumRate_new - sumRate_prev);

    ican.logf(params, "debug", ...
        "DC(sat) iter=%d K=%d sumRate=%.3f Mbps diff=%.3f kbps cvx=%s time=%.2fs", ...
        iter, K, sumRate_new/1e6, diff_bps/1e3, cvxStatus, cvxTime_s);

    if diff_bps < params.bfConvThresh_bps
        Q_prev = Q_new;
        sumRate_prev = sumRate_new;
        break;
    end

    Q_prev = Q_new;
    sumRate_prev = sumRate_new;
end

Q_opt = Q_prev;
sumRateFromQ_bps = sumRate_prev;

    % 采用秩1近似恢复波束向量（对应算法1第10步）。
w = complex(zeros(N, K));
for k = 1:K
    Qk = (Q_opt(:, :, k) + Q_opt(:, :, k)')/2;
    [V, D] = eig(Qk);
    [lambdaMax, idx] = max(real(diag(D)));
    lambdaMax = max(lambdaMax, 0);
    w(:, k) = sqrt(lambdaMax) * V(:, idx);
end

R_bps = rates_from_w(B, sigma2, h_s, w);

satRes = struct();
satRes.w = w;
satRes.R_bps = R_bps; % 1xK
satRes.Q = Q_opt;
satRes.sumRateFromQ_bps = sumRateFromQ_bps;
satRes.dcIters = dcIters;
satRes.cvxStatus = cvxStatusLast;
satRes.cvxOptval = cvxOptvalLast;
satRes.time_s = toc(ticSat);
end

function Q = init_random_Q(params, N, K)
%rng(params.randomSeed, "twister");
Q = complex(zeros(N, N, K));
for k = 1:K
    w = randn(N, 1) + 1j*randn(N, 1);
    w = w / norm(w, 2);
    w = sqrt(params.P_W) * w;
    Q(:, :, k) = w * w';
end
end

function [Q_new, cvxStatus, cvxOptval] = solve_cvx_dc_step(params, B, sigma2, P, Hi, Q_prev, denomPrev)
% 进行一次凸化的 DC 步骤（对应问题(15)）。

K = numel(Hi);
N = size(Q_prev, 1);

% 第一次尝试：默认静默求解
if params.cvx.quiet
    cvx_begin sdp quiet
else
    cvx_begin sdp
end
    variable Q(N, N, K) complex

    expression obj
    obj = 0;
    for i = 1:K
        powAll = 0;
        powInterf = 0;
        for k = 1:K
            pk = real(trace(Q(:, :, k) * Hi{i}));
            powAll = powAll + pk;
            if k ~= i
                powInterf = powInterf + pk;
            end
        end
        
        % log(sigma2 + x) = log(sigma2) + log(1 + x/sigma2)。
        powAllN = powAll / sigma2;
        powInterfN = powInterf / sigma2;
        denomPrevN = denomPrev(i) / sigma2;
        obj = obj + (B/log(2))*log(1 + powAllN) - (B/(log(2)*denomPrevN))*powInterfN;
    end
    maximize(obj)
    subject to
        totalPower = 0;
        for k = 1:K
            Q(:, :, k) == hermitian_semidefinite(N);
            totalPower = totalPower + real(trace(Q(:, :, k)));
        end
        totalPower <= P;
cvx_end

Q_new = Q;
cvxStatus = string(cvx_status);
cvxOptval = cvx_optval;
end

function sumRate_bps = sum_rate_from_Q(B, sigma2, Hi, Q)

K = numel(Hi);
rates = zeros(K, 1);
for i = 1:K
    pow = zeros(K, 1);
    for k = 1:K
        pow(k) = real(trace(Q(:, :, k) * Hi{i}));
    end
    num = sigma2 + sum(pow);
    den = sigma2 + (sum(pow) - pow(i));
    rates(i) = B * log2(num/den);
end
sumRate_bps = sum(rates);
end

function R_bps = rates_from_w(B, sigma2, h_s, w)

[~, K] = size(h_s);
R_bps = zeros(1, K);
for i = 1:K
    hi = h_s(:, i);
    desired = abs(hi' * w(:, i))^2;
    interf = 0;
    for k = 1:K
        if k == i
            continue;
        end
        interf = interf + abs(hi' * w(:, k))^2;
    end
    sinr = desired / (interf + sigma2);
    R_bps(i) = B * log2(1 + sinr);
end
end

function [Q_mapped, ok] = try_map_warmstart(params, Q_old, oldUeIdx, newUeIdx)


ok = false;
Q_mapped = [];

oldUeIdx = oldUeIdx(:).';
newUeIdx = newUeIdx(:).';

N = size(Q_old, 1);
K = numel(newUeIdx);
if size(Q_old, 1) ~= size(Q_old, 2)
    return;
end

Q_mapped = complex(zeros(N, N, K));
for k = 1:K
    ue = newUeIdx(k);
    oldPos = find(oldUeIdx == ue, 1);
    if isempty(oldPos)
       
        w = randn(N, 1) + 1j*randn(N, 1);
        w = w / norm(w, 2);
        w = sqrt(params.P_W) * w;
        Q_mapped(:, :, k) = w * w';
    else
        Q_mapped(:, :, k) = Q_old(:, :, oldPos);
    end
end

ok = true;
end
