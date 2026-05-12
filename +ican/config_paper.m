function params = config_paper(varargin)
% CONFIG_PAPER  论文仿真参数统一配置。
 

parser = inputParser();
parser.addParameter("randomSeed", 1, @(x) isnumeric(x) && isscalar(x));
parser.addParameter("I", 6, @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 3);
parser.addParameter("S", 12, @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 1);
parser.addParameter("C", 7, @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 1);
parser.parse(varargin{:});
opts = parser.Results;

params = struct();

% --- 系统参数 ---
params.satHeight_m = 600e3;
% 卫星轨道高度范围（015 仓库使用范围 600km - 1200km）
params.satAltMin_m = 600e3;
params.satAltMax_m = 1200e3;
params.cellRadius_m = 43.3e3;
params.Nx = 4;
params.Ny = 4;
params.fc_Hz = 2e9;
params.bandwidth_Hz = 10e6; % 对应远程复现代码的带宽
% 使用与远程仓库一致的发射功率标度（远程 config: P_max = 30 dBm -> 0 dBW）
params.P_dBw = 0; % dBW
params.noisePSD_dBmHz = -174;
params.gdopThreshold = 6; % 对应公式(9c)中的阈值 gamma
params.bfConvThresh_bps = 2e6; % 对应算法1中的收敛阈值 delta

% --- 派生参数 ---
params.c0 = 299792458; % 光速
params.lambda_m = params.c0 / params.fc_Hz; % 波长
params.N = params.Nx * params.Ny;
params.P_W = 10^(params.P_dBw/10); % dBw -> W
params.noisePSD_WHz = 10^((params.noisePSD_dBmHz - 30)/10);
params.sigma2_W = params.noisePSD_WHz * params.bandwidth_Hz;

% --- 问题规模 ---
params.S = opts.S;
params.C = opts.C;
params.I = opts.I;

% --- Fig.3 几何场景默认值 ---
 
params.scenario = struct();
% 采用与 GitHub 015 仓库一致的默认场景风格：多轨道面卫星 + 局部聚集用户
params.scenario.satLayout = "multi_orbit"; % 可选: multi_orbit / uniform_ring / clustered_ring
params.scenario.satCoreAngles_deg = [0 90 180 270];
params.scenario.satExtraAngles_deg = [15 30 45];
params.scenario.uePosMode = "github_offsets"; % 用户在某一局部服务区附近聚集（如北京附近）
% 可选：使用远端仓库式的经纬度偏移向量来按规则分散用户（单位：度）
params.scenario.geo_lat_deg = 0;   % 参考经度/纬度（deg），远端仓库常用 geo_lat=0, geo_lon=100
params.scenario.geo_lon_deg = 100;
params.scenario.lu_lat_offset = [0, 2, 4, 6];
params.scenario.lu_lon_offset = [0, 0, 0, 0];

% --- 信道模型开关（对应公式(1)） ---
params.atmosAtten = 1.0; % 对应公式(1)中的大气衰减项，论文未给出明确数值

% --- CVX 求解设置 ---
params.cvx = struct();
params.cvx.solver = "sdpt3";
params.cvx.quiet = true; % 设为 false 可查看完整求解日志

% --- 算法上限（安全保护） ---
params.alg = struct();
params.alg.maxDcIters = 5;
params.alg.utilityTol = 1e-9; % 若新效用不低于旧效用减去该容差，则接受
params.alg.maxCfgPasses = 10; % 允许多轮配置迭代，直到收敛

% --- 用户分组配置 ---
params.grouping = struct();
params.grouping.method = "spectral"; % 可选："spectral" 或 "qos"
params.grouping.numGroups = 3;
params.grouping.sigma_d = 50e3;
params.grouping.sigma_h = 0.5;
params.grouping.qosRatios = [0.5, 0.3, 0.2];

% --- MAB 配置 ---
params.mabMaxArms = 20;
params.mabRho = 0.98;
params.mabCucb = 1.0;
params.wdopPenaltyLambda = 0.5;
params.useParetoUCB = false;
params.paretoAlpha = 0.5;

% --- 日志配置 ---
params.log = struct();
params.log.level = "info"; % 可选：debug / info / warn / error
params.log.toFile = true;
params.log.dir = fullfile("matlab", "logs");

% --- 可复现性 ---
params.randomSeed = opts.randomSeed;

% --- WDOP 阈值 ---
params.wdopThreshold = 6.0; 
