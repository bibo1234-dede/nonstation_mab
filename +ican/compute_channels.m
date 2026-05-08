function chan = compute_channels(params, scenario)


%rng(params.randomSeed, "twister");

S = params.S;
C = params.C;
N = params.N;
Nx = params.Nx;
Ny = params.Ny;

pUE = scenario.pUE;
pSat = scenario.pSat;

h = complex(zeros(N, C, S));
d_m = zeros(S, C);
dirCos = zeros(S, C, 3);

for s = 1:S
    ps = pSat(s, :).';
    for c = 1:C
        pc = pUE(c, :).';
        diff = pc - ps;
        d = norm(diff, 2);
        d_m(s, c) = d;

        dir = diff ./ d;
        dirCos(s, c, :) = dir;

        
        theta_x = dir(1);
        theta_y = dir(2);

        vx = (1/sqrt(Nx)) * exp(-1j*pi*(0:Nx-1).'*theta_x);
        vy = (1/sqrt(Ny)) * exp(-1j*pi*(0:Ny-1).'*theta_y);
        v = kron(vx, vy); % (Nx*Ny) x 1

        % Pathloss (Eq.(1)).
        g_pl = (params.lambda_m / (4*pi*d))^2;
        g_at = params.atmosAtten;

        % Random phase.
        theta = 2*pi*rand();

        h(:, c, s) = sqrt(g_pl * g_at * N) * exp(-1j*theta) * v;
    end
end

chan = struct();
chan.h = h;
chan.d_m = d_m;
chan.dirCos = dirCos;
end

