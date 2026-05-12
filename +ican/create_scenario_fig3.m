function scenario = create_scenario_fig3(params)
% CREATE_SCENARIO_FIG3  统一为 Python 风格的卫星与 UE 生成
% 仅保留多轨道面卫星生成和基于地理中心的小范围 UE 生成

earthRadius_m = 6371e3;

% 卫星轨道高度范围
if isfield(params, 'satAltMin_m') && isfield(params, 'satAltMax_m')
    alt_min = params.satAltMin_m;
    alt_max = params.satAltMax_m;
else
    alt_min = params.satHeight_m;
    alt_max = params.satHeight_m;
end

S = params.S;
C = params.C;

% 随机性
if isfield(params, 'randomSeed')
    rng(params.randomSeed, 'twister');
else
    rng('shuffle');
end

% 多轨道面卫星分布（参考 Python 实现）
num_orbits = max(2, floor(S / 4));
sats_per_orbit = floor(S / num_orbits);
extra = mod(S, num_orbits);

pSat = zeros(S, 3);
sat_radius = zeros(S, 1);
sat_raan = zeros(S, 1);
sat_inc = zeros(S, 1);
sat_phase0 = zeros(S, 1);
sat_omega = zeros(S, 1);
sats_orbit_id = zeros(S, 1);
sidx = 1;
for o = 1:num_orbits
    n_sats = sats_per_orbit + (o <= extra);
    % RAAN 小扰动
    raan = 2*pi*(o-1)/num_orbits + (rand() - 0.5) * 0.2;
    % 近极地轨道倾角带小扰动
    inc = deg2rad(85 + (rand() - 0.5) * 10);
    for j = 1:n_sats
        phase = 2*pi*(j-1)/n_sats + (rand() - 0.5) * 0.1;
        if alt_max > alt_min
            alt = alt_min + (alt_max - alt_min) * rand();
        else
            alt = alt_min;
        end
        r = earthRadius_m + alt;
        x = r * (cos(raan) * cos(phase) - sin(raan) * sin(phase) * cos(inc));
        y = r * (sin(raan) * cos(phase) + cos(raan) * sin(phase) * cos(inc));
        z = r * sin(phase) * sin(inc);
        pSat(sidx, :) = [x, y, z];

        sat_radius(sidx) = r;
        sat_raan(sidx) = raan;
        sat_inc(sidx) = inc;
        sat_phase0(sidx) = phase;
        sats_orbit_id(sidx) = o;

        % 圆轨道角速度（rad/s），用于时序位置更新
        mu = 3.986004418e14;
        sat_omega(sidx) = sqrt(mu / (r^3));

        sidx = sidx + 1;
        if sidx > S
            break;
        end
    end
    if sidx > S
        break;
    end
end

% UE 位置：以北京附近为中心的小范围随机扰动（与 Python 保持一致）
lat0 = deg2rad(40);
lon0 = deg2rad(116);
rng(params.randomSeed + 1000, 'twister');
pUE = zeros(C, 3);
for c = 1:C
    lat = lat0 + (rand() * 2 - 1) * 0.05; % ±0.05 rad
    lon = lon0 + (rand() * 2 - 1) * 0.05; % ±0.05 rad
    r = earthRadius_m;
    x = r * cos(lat) * cos(lon);
    y = r * cos(lat) * sin(lon);
    z = r * sin(lat);
    pUE(c, :) = [x, y, z];
end

scenario = struct();
scenario.pUE = pUE;   % C x 3
scenario.pSat = pSat; % S x 3
scenario.satLayout = 'multi_orbit';
scenario.satAltMin_m = alt_min;
scenario.satAltMax_m = alt_max;
scenario.R_orbit = mean(sqrt(sum(pSat.^2, 2))); % 近似轨道半径
scenario.geo_center = [rad2deg(lat0), rad2deg(lon0)];

scenario.satDynamics = struct();
scenario.satDynamics.radius_m = sat_radius;
scenario.satDynamics.raan_rad = sat_raan;
scenario.satDynamics.inc_rad = sat_inc;
scenario.satDynamics.phase0_rad = sat_phase0;
scenario.satDynamics.omega_radps = sat_omega;
scenario.satDynamics.orbit_id = sats_orbit_id;
end
