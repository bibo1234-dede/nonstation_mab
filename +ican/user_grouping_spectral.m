function [user_groups, similarity_matrix] = user_grouping_spectral(params, scenario, chan)
% USER_GROUPING_SPECTRAL  基于谱聚类的用户分组。

C = params.C;
K = min(get_grouping_field(params, "numGroups", C), C);
if K < 1
    K = 1;
end

sigma_d = get_grouping_field(params, "sigma_d", 50e3);
sigma_h = get_grouping_field(params, "sigma_h", 0.5);

similarity_matrix = zeros(C, C);
user_gain = zeros(C, 1);
for i = 1:C
    user_gain(i) = mean(abs(chan.h(:, i, :)).^2, "all");
end

for i = 1:C
    for j = 1:C
        d_ij = norm(scenario.pUE(i, :) - scenario.pUE(j, :), 2);
        s_geo = exp(-(d_ij.^2) / (2 * sigma_d^2));

        s_chan = exp(-((user_gain(i) - user_gain(j)).^2) / (2 * sigma_h^2));
        similarity_matrix(i, j) = s_geo * s_chan;
    end
end

similarity_matrix = max(similarity_matrix, similarity_matrix.');
for i = 1:C
    similarity_matrix(i, i) = 1;
end

D = diag(sum(similarity_matrix, 2));
D_inv_sqrt = diag(1 ./ sqrt(max(diag(D), eps)));
L_sym = eye(C) - D_inv_sqrt * similarity_matrix * D_inv_sqrt;
L_sym = (L_sym + L_sym.') / 2;

if K == 1
    group_idx = ones(C, 1);
else
    [V, lambda_vals] = eig(full(L_sym)); %#ok<ASGLU>
    [~, order] = sort(real(diag(lambda_vals)), "ascend");
    U = real(V(:, order(1:K)));
    U = normalize_rows(U);
    group_idx = kmeans(U, K, "Replicates", 5, "MaxIter", 200, "Display", "off");
end

user_groups = struct();
user_groups.method = "spectral";
user_groups.numGroups = K;
user_groups.group_ids = cell(K, 1);
user_groups.weights = ones(K, 1);
user_groups.group_assignment = group_idx(:);
user_groups.group_stats = repmat(struct("num_users", 0), K, 1);
for k = 1:K
    user_groups.group_ids{k} = find(group_idx == k);
    user_groups.group_stats(k).num_users = numel(user_groups.group_ids{k});
end
end

function value = get_grouping_field(params, fieldName, defaultValue)
value = defaultValue;
if isfield(params, "grouping") && isfield(params.grouping, fieldName)
    value = params.grouping.(fieldName);
end
end

function U = normalize_rows(U)
row_norms = sqrt(sum(U.^2, 2));
row_norms(row_norms == 0) = 1;
U = U ./ row_norms;
end

