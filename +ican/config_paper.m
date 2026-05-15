
function params = config_paper(varargin)
% CONFIG_PAPER  论文仿真参数统一配置。
 

parser = inputParser();
parser.addParameter("randomSeed", 1, @(x) isnumeric(x) && isscalar(x));
parser.addParameter("I", 4, @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 3);
parser.addParameter("S", 18, @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 1);
parser.addParameter("C", 10, @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 1);
parser.parse(varargin{:});
opts = parser.Results;

params = struct();

% --- 系统参数 ---
params.satHeight_m = 600e3;        % 参考卫星轨道高度，单位 m
% 卫星轨道高度范围：用于生成多轨道面卫星时抽样的上下界
params.satAltMin_m = 600e3;        % 最低轨道高度，单位 m
params.satAltMax_m = 1200e3;       % 最高轨道高度，单位 m
params.Nx = 8;                     % 阵列 x 方向天线数
params.Ny = 8;                     % 阵列 y 方向天线数
params.fc_Hz = 2e9;                % 载波频率，单位 Hz
params.bandwidth_Hz = 50e6;        % 系统带宽，单位 Hz
% 发射功率的 dBW 标度；0 dBW = 1 W，20 dBW = 100 W
params.P_dBw = 20;
params.noisePSD_dBmHz = -174;      % 热噪声功率谱密度，单位 dBm/Hz
params.noisePSD_WHz = 10^((params.noisePSD_dBmHz - 30)/10); % dBm/Hz -> W/Hz
params.sigma2_W = params.noisePSD_WHz * params.bandwidth_Hz;  % 接收机噪声功率
params.bfConvThresh_bps = 2e6;     % 波束成形 DC 迭代收敛阈值，单位 bps

% --- 派生参数 ---
params.c0 = 299792458;             % 光速，单位 m/s
params.lambda_m = params.c0 / params.fc_Hz; % 载波波长，单位 m
params.N = params.Nx * params.Ny;  % 单颗卫星阵列总天线数
params.P_W = 10^(params.P_dBw/10); % 发射功率，单位 W

% --- 等效信道增益补偿（保持为适中数值，避免 CVX 数值病态）---
params.effectiveGain_dB = 10;      % 额外信道增益补偿，单位 dB
params.effectiveGain_linear = 10^(params.effectiveGain_dB / 10);  % 对应线性功率倍数

% --- 问题规模 ---
params.S = opts.S;                 % 卫星总数
params.C = opts.C;                 % 用户总数
params.I = opts.I;                 % 每个用户选择的卫星数

% --- Fig.3 几何场景默认值 ---
 
params.scenario = struct();
% 采用与 GitHub 015 仓库一致的默认场景风格：多轨道面卫星 + 局部聚集用户
params.scenario.satLayout = "multi_orbit"; % 卫星布局模式：multi_orbit / uniform_ring / clustered_ring

% --- 信道模型开关（对应公式(1)） ---
params.atmosAtten = 1.0; % 大气衰减系数，1 表示不额外衰减

% --- CVX 求解设置 ---
params.cvx = struct();
params.cvx.solver = "mosek";      % CVX 求解器
params.cvx.quiet = true;          % true 时关闭 CVX 详细输出

% --- 算法上限（安全保护） ---
params.alg = struct();
params.alg.maxDcIters = 5;        % 单次卫星波束成形的 DC 最大迭代次数 (降低以提速 & 减少 MOSEK 数值问题)

% --- 用户分组配置 ---
% 仅支持 by_level（随机优先级分组）模式
params.grouping = struct();
params.useGrouping = false;           % 是否启用用户分组
params.grouping.numLevels = 3;        % 按级别分组时的级别数

% --- MAB 配置 ---
params.mabMaxArms = 50;             % 每个用户保留的最大候选臂数
params.candidatePoolSize = 120;      % 候选臂池大小
params.candidateRateWeight = 1.0;   % 候选臂构建时的速率权重
params.candidateWdopWeight = 0.25;  % 候选臂构建时的 WDOP 权重
params.debugPrintCandidateArms = false; % 是否打印候选臂
params.debugPrintSelection = false;      % 是否打印每轮选择结果
params.debugPrintCandidateLimit = 12;    % 最多打印多少个候选臂
params.mabRho = 0.98;               % 非稳态遗忘因子
params.mabCucb = 1.0;               % UCB 探索系数
params.wdopPenaltyLambda = 0.12;    % WDOP 惩罚权重
params.wdopSoftMargin = 1.5;        % WDOP 软惩罚边际
params.satLoadCap = 6;              % 单颗卫星允许的最大负载（从4降到3，减少竞争）
params.satLoadPenaltyLambda = 2.0;  % 卫星负载惩罚权重（从0.5提到1.0，加强分散激励）
params.groupReuseBonusLambda = 0.10; % 同组卫星重用奖励权重（降低，减弱聚集效应）
params.userReusePenaltyLambda = 0.50; % 跨组卫星重用惩罚权重（提高，鼓励分散）
params.useParetoUCB = false;        % 是否启用 Pareto-UCB
params.paretoAlpha = 0.5;           % Pareto-UCB 中双目标融合系数

% --- 日志配置 ---
params.log = struct();
params.log.level = "info";         % 日志级别：debug / info / warn / error
params.log.printArmSpace = false;   % 是否打印所有臂空间组合

% --- 可复现性 ---
params.randomSeed = opts.randomSeed; % 随机种子

% --- WDOP 阈值 ---
params.wdopThreshold = 6.0;         % WDOP 最大允许阈值
