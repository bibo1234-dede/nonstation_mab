function scenario = create_scenario_fig3(params)
 
if params.S ~= 7 || params.C ~= 7
    error("create_scenario_fig3:SizeMismatch", "This scenario helper assumes S=C=7 for Fig.3.");
end

 

isD_m = sqrt(3) * params.cellRadius_m; % inter-site distance for hex grid

cellCenter = zeros(params.C, 3);
cellCenter(4, :) = [0, 0, 0]; % UE4's cell is the central cell

ringIdx = [1 2 3 5 6 7];
ringAngles_deg = [0 60 120 180 240 300];
for k = 1:numel(ringIdx)
    idx = ringIdx(k);
    ang = deg2rad(ringAngles_deg(k));
    cellCenter(idx, :) = [isD_m*cos(ang), isD_m*sin(ang), 0];
end

pUE = cellCenter;
if isfield(params, "scenario") && isfield(params.scenario, "uePosMode")
    mode = string(params.scenario.uePosMode);
else
    mode = "cell_center";
end

if mode == "random_in_cell"
    if isfield(params, "scenario") && isfield(params.scenario, "ueInCellRadiusFrac")
        radius = params.cellRadius_m * params.scenario.ueInCellRadiusFrac;
    else
        radius = params.cellRadius_m;
    end
    rng(params.randomSeed, "twister");
    for c = 1:params.C
        % Uniform in disk: r = R*sqrt(u), angle = 2pi*v
        u = rand();
        v = rand();
        rr = radius * sqrt(u);
        ang = 2*pi*v;
        offset = [rr*cos(ang), rr*sin(ang), 0];
        pUE(c, :) = cellCenter(c, :) + offset;
    end
end

 

if isfield(params, "scenario") && isfield(params.scenario, "satRingRadius_m")
    satRingRadius_m = params.scenario.satRingRadius_m;
else
    satRingRadius_m = 130e3;
end

if isfield(params, "scenario") && isfield(params.scenario, "satLayout")
    satLayout = string(params.scenario.satLayout);
else
    satLayout = "uniform_ring";
end

pSat = zeros(params.S, 3);

if satLayout == "uniform_ring"
    for s = 1:params.S
        ang = 2*pi*(s-1)/params.S;
        pSat(s, :) = [satRingRadius_m*cos(ang), satRingRadius_m*sin(ang), params.satHeight_m];
    end
elseif satLayout == "clustered_ring"
    if params.S ~= 7
        error("create_scenario_fig3:BadLayout", "satLayout=clustered_ring currently assumes S=7.");
    end
    coreAngles = [0 90 180 270];
    extraAngles = [15 30 45];
    if isfield(params, "scenario") && isfield(params.scenario, "satCoreAngles_deg")
        coreAngles = params.scenario.satCoreAngles_deg;
    end
    if isfield(params, "scenario") && isfield(params.scenario, "satExtraAngles_deg")
        extraAngles = params.scenario.satExtraAngles_deg;
    end
    coreAngles = coreAngles(:).';
    extraAngles = extraAngles(:).';
    if numel(coreAngles) ~= 4 || numel(extraAngles) ~= 3
        error("create_scenario_fig3:BadLayout", "clustered_ring needs 4 core angles and 3 extra angles.");
    end
    angles = [coreAngles, extraAngles];
    for s = 1:params.S
        ang = deg2rad(angles(s));
        pSat(s, :) = [satRingRadius_m*cos(ang), satRingRadius_m*sin(ang), params.satHeight_m];
    end
else
    error("create_scenario_fig3:BadLayout", "Unknown satLayout: %s", satLayout);
end

scenario = struct();
scenario.pUE = pUE;   % Cx3
scenario.pSat = pSat; % Sx3
scenario.satRingRadius_m = satRingRadius_m;
scenario.isd_m = isD_m;
scenario.cellCenter = cellCenter;
scenario.uePosMode = mode;
scenario.satLayout = satLayout;
if satLayout == "clustered_ring"
    scenario.satAngles_deg = angles;
else
    scenario.satAngles_deg = (0:params.S-1) * (360/params.S);
end
end
