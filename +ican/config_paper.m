function params = config_paper(varargin)
 

parser = inputParser();
parser.addParameter("randomSeed", 1, @(x) isnumeric(x) && isscalar(x));
parser.addParameter("I", 4, @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 3);
parser.addParameter("S", 7, @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 1);
parser.addParameter("C", 7, @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 1);
parser.parse(varargin{:});
opts = parser.Results;

params = struct();

% --- Table I ---
params.satHeight_m = 600e3;
params.cellRadius_m = 43.3e3;
params.Nx = 8;
params.Ny = 8;
params.fc_Hz = 4e9;
params.bandwidth_Hz = 50e6; % B in Eq.(5)
params.P_dBw = 26; % per-beam power budget in Eq.(9d)
params.noisePSD_dBmHz = -174;
params.gdopThreshold = 6; % gamma in Eq.(9c)
params.bfConvThresh_bps = 2e6; % delta in Algorithm 1 (0.5 Mbps)

% --- Derived ---
params.c0 = 299792458; % speed of light
params.lambda_m = params.c0 / params.fc_Hz; % wavelength (md->pdf conversion sometimes flips this)
params.N = params.Nx * params.Ny;
params.P_W = 10^(params.P_dBw/10); % dBw -> W
params.noisePSD_WHz = 10^((params.noisePSD_dBmHz - 30)/10);
params.sigma2_W = params.noisePSD_WHz * params.bandwidth_Hz;

% --- Problem size ---
params.S = opts.S;
params.C = opts.C;
params.I = opts.I;

% --- Fig.3 snapshot geometry defaults ---
 
params.scenario = struct();
params.scenario.satRingRadius_m = 130e3;
params.scenario.satLayout = "clustered_ring";
params.scenario.satCoreAngles_deg = [0 90 180 270];
params.scenario.satExtraAngles_deg = [15 30 45];
params.scenario.uePosMode = "cell_center";

% --- Channel model toggles (Eq.(1)) ---
params.atmosAtten = 1.0; % G^{at} in Eq.(1). Paper references 3GPP; not numerically specified.

% --- CVX/MOSEK settings ---
params.cvx = struct();
params.cvx.solver = "mosek";
params.cvx.quiet = true; % set false for full solver logs

% --- Algorithm caps (safety) ---
params.alg = struct();
params.alg.maxDcIters = 5;
params.alg.utilityTol = 1e-9; % accept if U_new >= U_old - tol
params.alg.maxCfgPasses = 10; % allow multiple passes until convergence

% --- Logging ---
params.log = struct();
params.log.level = "info"; % debug|info|warn|error
params.log.toFile = true;
params.log.dir = fullfile("matlab", "logs");

% --- Reproducibility ---
params.randomSeed = opts.randomSeed;

% --- WDOP threshold ---
params.wdopThreshold = 6.0; 
