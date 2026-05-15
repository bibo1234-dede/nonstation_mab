function [user_groups, level_info] = user_grouping_random_priority(params, scenario, chan)
% USER_GROUPING_RANDOM_PRIORITY  基于随机优先级的用户分组。
%
% 将用户随机分为多个优先级：
%   优先级1（Tier 1）：最先处理，选择最优卫星组合
%   优先级2（Tier 2）：次先处理，选择次优卫星组合
%   优先级3（Tier 3）：最后处理，选择剩余卫星
%
% 然后在 MAB 选择时根据优先级来约束候选臂的选择策略。

C = params.C;
num_levels = get_level_field(params, "numLevels", 3);

% 计算每个用户的平均信道增益（用于 MAB 选择时的排序参考）
user_quality = zeros(C, 1);
for c = 1:C
	% 平均信道增益：所有卫星到用户 c 的平均 |h|^2
	user_quality(c) = mean(abs(chan.h(:, c, :)).^2, "all");
end

% 随机分配用户到各优先级
perm = randperm(C);
level_assignment = zeros(C, 1);
users_per_level = ceil(C / num_levels);

for level = 1:num_levels
	start_idx = (level - 1) * users_per_level + 1;
	end_idx = min(level * users_per_level, C);

	user_indices = perm(start_idx:end_idx);
	level_assignment(user_indices) = level;
end

% 构造 user_groups 结构体（兼容现有接口）
user_groups = struct();
user_groups.method = "random_priority";
user_groups.numGroups = num_levels;
user_groups.group_ids = cell(num_levels, 1);
user_groups.weights = ones(num_levels, 1);
user_groups.group_assignment = level_assignment(:);
user_groups.group_stats = repmat(struct("num_users", 0), num_levels, 1);
user_groups.user_priority_level = level_assignment(:);  % 保存优先级

for level = 1:num_levels
	user_groups.group_ids{level} = find(level_assignment == level);
	user_groups.group_stats(level).num_users = numel(user_groups.group_ids{level});
end

% 输出详细信息
level_info = struct();
level_info.user_quality = user_quality;
level_info.level_assignment = level_assignment;
level_info.user_priority_levels = level_assignment;

% 打印分组信息
fprintf("=== 基于随机优先级的用户分组 ===\n");
for level = 1:num_levels
	user_list = user_groups.group_ids{level};
	quality_vals = user_quality(user_list);
	fprintf("优先级 %d（%s）：用户 %s | 平均质量=%.6f\n", ...
		level, level_name(level), mat2str(user_list.'), mean(quality_vals));
end
fprintf("\n");

end

function name = level_name(level)
	names = ["高", "中", "低"];
	if level >= 1 && level <= numel(names)
		name = names(level);
	else
		name = "未知";
	end
end

function value = get_level_field(params, fieldName, defaultValue)
value = defaultValue;
if isfield(params, "grouping") && isfield(params.grouping, fieldName)
	value = params.grouping.(fieldName);
end
end