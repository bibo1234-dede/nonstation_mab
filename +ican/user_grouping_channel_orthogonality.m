function [user_groups, group_info] = user_grouping_channel_orthogonality(params, scenario, chan)
% USER_GROUPING_CHANNEL_ORTHOGONALITY  基于信道正交性的准静态分组。
%
% 目标：让同组用户的信道尽量彼此正交，从而降低同组内复用同一卫星组合时的干扰。

C = params.C;
if isfield(params, "grouping") && isfield(params.grouping, "numGroups")
    num_groups = max(1, min(C, round(params.grouping.numGroups)));
elseif isfield(params, "grouping") && isfield(params.grouping, "numLevels")
    num_groups = max(1, min(C, round(params.grouping.numLevels)));
else
    num_groups = min(3, C);
end

avg_gain = squeeze(mean(mean(abs(chan.h).^2, 1), 3));
avg_gain = avg_gain(:);

similarity_matrix = compute_user_similarity_matrix(chan);
group_ids = cell(num_groups, 1);
group_assignment = zeros(C, 1);

target_sizes = distribute_sizes(C, num_groups);
[~, seed_order] = sort(avg_gain, "descend");

% 先给每个组放一个种子用户，避免空组
for g = 1:num_groups
    seed_user = seed_order(g);
    group_ids{g} = seed_user;
    group_assignment(seed_user) = g;
end

remaining_users = seed_order(num_groups + 1:end);
for idx = 1:numel(remaining_users)
    c = remaining_users(idx);
    best_group = 1;
    best_cost = inf;
    for g = 1:num_groups
        if numel(group_ids{g}) >= target_sizes(g)
            continue;
        end
        members = group_ids{g};
        cost = mean(similarity_matrix(c, members));
        balance_penalty = 0.05 * (numel(group_ids{g}) / max(target_sizes(g), 1));
        total_cost = cost + balance_penalty;
        if total_cost < best_cost
            best_cost = total_cost;
            best_group = g;
        end
    end
    group_ids{best_group} = [group_ids{best_group}, c]; %#ok<AGROW>
    group_assignment(c) = best_group;
end

% 组内按平均增益排序，最高者作为组长
leader_ids = zeros(num_groups, 1);
group_quality = zeros(num_groups, 1);
intra_similarity = zeros(num_groups, 1);
for g = 1:num_groups
    members = group_ids{g};
    if isempty(members)
        continue;
    end
    [~, ord] = sort(avg_gain(members), "descend");
    members = members(ord);
    group_ids{g} = members(:).';
    leader_ids(g) = members(1);
    group_quality(g) = mean(avg_gain(members));
    if numel(members) > 1
        sub_sim = similarity_matrix(members, members);
        intra_similarity(g) = mean(sub_sim(~eye(size(sub_sim, 1))));
    else
        intra_similarity(g) = 0;
    end
end

% 将组按“更正交”的顺序排序，便于 MAB 遍历时优先处理更稳定的组
[~, order] = sort(intra_similarity, "ascend");
group_ids = group_ids(order);
leader_ids = leader_ids(order);
group_quality = group_quality(order);
intra_similarity = intra_similarity(order);

% 重新生成组映射
group_assignment = zeros(C, 1);
for g = 1:num_groups
    group_assignment(group_ids{g}) = g;
end

group_weights = 1 ./ max(intra_similarity(:) + 0.2, eps);
group_weights = group_weights / max(mean(group_weights), eps);

user_groups = struct();
user_groups.method = "orthogonality";
user_groups.numGroups = num_groups;
user_groups.group_ids = group_ids;
user_groups.weights = group_weights;
user_groups.group_assignment = group_assignment;
user_groups.group_stats = repmat(struct("num_users", 0, "leader", 0, "avg_gain", 0, "intra_similarity", 0), num_groups, 1);
user_groups.user_priority_level = ones(C, 1);
user_groups.leader_ids = leader_ids;

for g = 1:num_groups
    user_groups.group_stats(g).num_users = numel(group_ids{g});
    user_groups.group_stats(g).leader = leader_ids(g);
    user_groups.group_stats(g).avg_gain = group_quality(g);
    user_groups.group_stats(g).intra_similarity = intra_similarity(g);
end

group_info = struct();
group_info.similarity_matrix = similarity_matrix;
group_info.avg_gain = avg_gain;
group_info.leader_ids = leader_ids;
group_info.group_quality = group_quality;
group_info.intra_similarity = intra_similarity;

fprintf("=== 基于信道正交性的用户分组 ===\n");
for g = 1:num_groups
    fprintf("组 %d | 组长 UE%d | 用户 %s | 平均增益=%.6f | 组内相似度=%.6f\n", ...
        g, leader_ids(g), mat2str(group_ids{g}), group_quality(g), intra_similarity(g));
end
fprintf("\n");

end

function similarity_matrix = compute_user_similarity_matrix(chan)
C = size(chan.h, 2);
S = size(chan.h, 3);
similarity_matrix = eye(C);
for c1 = 1:C
    h1 = squeeze(chan.h(:, c1, :));
    for c2 = c1 + 1:C
        h2 = squeeze(chan.h(:, c2, :));
        sim_vals = zeros(S, 1);
        valid_count = 0;
        for s = 1:S
            v1 = h1(:, s);
            v2 = h2(:, s);
            denom = norm(v1, 2) * norm(v2, 2);
            if denom > 0
                valid_count = valid_count + 1;
                sim_vals(valid_count) = abs(v1' * v2) / denom;
            end
        end
        if valid_count == 0
            sim = 0;
        else
            sim = mean(sim_vals(1:valid_count));
        end
        similarity_matrix(c1, c2) = sim;
        similarity_matrix(c2, c1) = sim;
    end
end
end

function target_sizes = distribute_sizes(C, G)
base = floor(C / G);
remainder = mod(C, G);
target_sizes = base * ones(G, 1);
target_sizes(1:remainder) = target_sizes(1:remainder) + 1;
end