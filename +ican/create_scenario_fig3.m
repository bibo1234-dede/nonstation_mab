function scenario = create_scenario_fig3(params)
% CREATE_SCENARIO_FIG3  构建低轨卫星与地面用户的几何场景。

earthRadius_m = 6371e3;
gm_earth = 3.986e14;

if isfield(params, "scenario") && isfield(params.scenario, "satRingRadius_m")
    satRingRadius_m = params.scenario.satRingRadius_m;
else
    satRingRadius_m = 4000e3;
end

if isfield(params, "scenario") && isfield(params.scenario, "satLayout")
    satLayout = string(params.scenario.satLayout);
else
    satLayout = "uniform_ring";
end

isD_m = sqrt(3) * params.cellRadius_m;
layers = build_hex_layers(params.C, params.cellRadius_m);
cellCenter = layers(1:params.C, :);

pUE = cellCenter;
pUE(:, 3) = earthRadius_m;   % 把用户放到地球表面
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
        u = rand();
        v = rand();
        rr = radius * sqrt(u);
        ang = 2*pi*v;
        offset = [rr*cos(ang), rr*sin(ang), 0];
        pUE(c, :) = cellCenter(c, :) + offset;
    end
elseif mode == "clustered_local"
    % 参考远程复现：在北京附近（lat ~40°, lon ~116°）局部随机分布
    lat0 = deg2rad(40);
    lon0 = deg2rad(116);
    rng(params.randomSeed + 1000, "twister");
    for c = 1:params.C
        lat = lat0 + (rand() * 2 - 1) * 0.02; % ±0.02 rad
        lon = lon0 + (rand() * 2 - 1) * 0.02;
        r = earthRadius_m;
        x = r * cos(lat) * cos(lon);
        y = r * cos(lat) * sin(lon);
        z = r * sin(lat);
        pUE(c, :) = [x, y, z];
    end
elseif mode == "github_offsets"
    % 使用远端仓库的经纬度偏移向量按规则分散用户（单位：度）
    if isfield(params.scenario, "lu_lat_offset")
        lat_offsets = params.scenario.lu_lat_offset;
    else
        lat_offsets = zeros(1, params.C);
    end
    if isfield(params.scenario, "lu_lon_offset")
        lon_offsets = params.scenario.lu_lon_offset;
    else
        lon_offsets = zeros(1, params.C);
    end
    if isfield(params.scenario, "geo_lat_deg")
        geo_lat_deg = params.scenario.geo_lat_deg;
    else
        geo_lat_deg = 0;
    end
    if isfield(params.scenario, "geo_lon_deg")
        geo_lon_deg = params.scenario.geo_lon_deg;
    else
        geo_lon_deg = 100;
    end
    for c = 1:params.C
        lat_deg = geo_lat_deg + lat_offsets(mod(c-1, numel(lat_offsets))+1);
        lon_deg = geo_lon_deg + lon_offsets(mod(c-1, numel(lon_offsets))+1);
        lat = deg2rad(lat_deg);
        lon = deg2rad(lon_deg);
        r = earthRadius_m;
        x = r * cos(lat) * cos(lon);
        y = r * cos(lat) * sin(lon);
        z = r * sin(lat);
        pUE(c, :) = [x, y, z];
    end
end

theta0_s = (0:params.S-1).' * 2*pi/params.S;
R_orbit = earthRadius_m + params.satHeight_m;
omega = sqrt(gm_earth / (R_orbit^3));

pSat = zeros(params.S, 3);
if satLayout == "multi_orbit"
    % 参考远程复现的多轨道面生成逻辑：多个轨道面、每面若干颗卫星，并带微小扰动
    num_orbits = max(2, floor(params.S / 4));
    sats_per_orbit = floor(params.S / num_orbits);
    extra = mod(params.S, num_orbits);
    sidx = 1;
    rng(params.randomSeed, "twister");
    for o = 1:num_orbits
        n_sats = sats_per_orbit + (o <= extra);
        raan = 2 * pi * (o - 1) / num_orbits + (rand() - 0.5) * 0.2; % small RAAN jitter
        inc = deg2rad(85 + (rand() - 0.5) * 10); % near-polar ~85° ±5°
        for j = 1:n_sats
            phase = 2 * pi * (j - 1) / n_sats + (rand() - 0.5) * 0.1; % along-track jitter
            if isfield(params, 'satAltMin_m') && isfield(params, 'satAltMax_m')
                alt = params.satAltMin_m + (params.satAltMax_m - params.satAltMin_m) * rand();
            else
                alt = params.satHeight_m;
            end
            r = earthRadius_m + alt;
            x = r * (cos(raan) * cos(phase) - sin(raan) * sin(phase) * cos(inc));
            y = r * (sin(raan) * cos(phase) + cos(raan) * sin(phase) * cos(inc));
            z = r * sin(phase) * sin(inc);
            pSat(sidx, :) = [x, y, z];
            sidx = sidx + 1;
            if sidx > params.S
                break;
            end
        end
        if sidx > params.S
            break;
        end
    end
elseif satLayout == "clustered_ring"
    if params.S ~= 7
        error("create_scenario_fig3:BadLayout", "satLayout=clustered_ring currently assumes S=7.");
    end
    angles_deg = [0 90 180 270 15 30 45];
    for s = 1:params.S
        ang = deg2rad(angles_deg(s));
        pSat(s, :) = [R_orbit*cos(ang), R_orbit*sin(ang), params.satHeight_m];
    end
    theta0_s = deg2rad(angles_deg(:));
else
    error("create_scenario_fig3:BadLayout", "Unknown satLayout: %s", satLayout);
end

scenario = struct();
scenario.pUE = pUE;   % Cx3
scenario.pSat = pSat; % Sx3
scenario.satRingRadius_m = R_orbit;
scenario.isd_m = isD_m;
scenario.cellCenter = cellCenter;
scenario.uePosMode = mode;
scenario.satLayout = satLayout;
scenario.satAngles_deg = rad2deg(theta0_s(:).');
scenario.theta0_s = theta0_s;
scenario.omega = omega;
scenario.R_orbit = R_orbit;
end

function layers = build_hex_layers(C, cellRadius_m)
% 构建以原点为中心的蜂窝式用户布局。

layers = zeros(max(C, 1), 3);
if C <= 0
    return;
end

layers(1, :) = [0, 0, 0];
if C == 1
    return;
end

isD_m = sqrt(3) * cellRadius_m;
ring1_angles = deg2rad(0:60:300);
ring2_angles = deg2rad(0:30:330);
ring3_angles = deg2rad(0:20:340);
ring4_angles = deg2rad(0:15:345);
ring_sets = {
    [ones(6,1), ones(6,1) * isD_m, ring1_angles(:)],
    [2*ones(12,1), 2*ones(12,1) * isD_m, ring2_angles(:)],
    [3*ones(18,1), 3*ones(18,1) * isD_m, ring3_angles(:)],
    [4*ones(24,1), 4*ones(24,1) * isD_m, ring4_angles(:)]
};

idx = 2;
for r = 1:numel(ring_sets)
    ring = ring_sets{r};
    for k = 1:size(ring, 1)
        if idx > C
            return;
        end
        radius = ring(k, 2);
        ang = ring(k, 3);
        layers(idx, :) = [radius*cos(ang), radius*sin(ang), 0];
        idx = idx + 1;
    end
end

while idx <= C
    ring = ceil((idx - 1) / 6);
    ang = 2*pi * mod(idx - 2, 6) / 6;
    radius = ring * sqrt(3) * cellRadius_m;
    layers(idx, :) = [radius*cos(ang), radius*sin(ang), 0];
    idx = idx + 1;
end
end
