# 卫星移动的非稳态环境说明

## 什么改变了？

原始版本的"非稳态"是通过**随机信道变化**实现的（每个时间步的随机相位）。

新版本加入了**卫星轨道运动**，这使得非稳态环境更加真实和复杂：
- 卫星位置随时间变化
- 距离 (path loss) 随之改变
- 信道结合了卫星位置变化 + 随机变化

## 核心实现

在 `run_fig3_nonstationary.m` 主循环的**第 0 步**（在计算信道之前）：

```matlab
% 0. 更新卫星位置（轨道运动）
for s = 1:params.S
    current_angle = sat_initial_angles(s) + satellite_angular_velocities(s) * (t - 1);
    scenario.pSat(s, 1) = satellite_radius * cos(current_angle);
    scenario.pSat(s, 2) = satellite_radius * sin(current_angle);
    % 高度保持不变
end
```

**关键参数** (在脚本开头定义)：

| 参数 | 含义 | 默认值 | 建议调整 |
|------|------|--------|---------|
| `satellite_angular_velocities` | 每个卫星的角速度（弧度/时间步） | `2*pi / 100` | 调整数值改变周期 |
| `sat_initial_angles` | 每个卫星的初始角度 | 从 scenario 自动计算 | 一般无需改动 |
| `satellite_radius` | 卫星轨道半径（米） | `scenario.satRingRadius_m` | 一般无需改动 |

## 如何调整运动参数？

### 场景 1：所有卫星以相同角速度旋转（当前默认）

```matlab
satellite_angular_velocities = ones(params.S, 1) * 2*pi / 100;  % 100 步一圈
```

修改参数改变速度/周期：
- `2*pi / 50`  → 50 步完成一圈（更快）
- `2*pi / 200` → 200 步完成一圈（更慢）
- `2*pi / 10`  → 10 步完成一圈（非常快）

### 场景 2：每个卫星有不同的角速度

```matlab
satellite_angular_velocities = [1.0; 1.2; 0.8; 1.0; 1.5; 0.7; 1.1] * 2*pi / 100;
```

这会让卫星以不同速度运动，造成更复杂的相对位置变化。

### 场景 3：分段运动（卫星定期加速/减速）

```matlab
% 前 100 步快速移动，后 100 步慢速移动
if t <= 100
    satellite_angular_velocities = ones(params.S, 1) * 2*pi / 50;
else
    satellite_angular_velocities = ones(params.S, 1) * 2*pi / 200;
end
```

### 场景 4：正逆方向运动

```matlab
% 奇数卫星顺时针，偶数卫星逆时针
satellite_angular_velocities = repmat([1; -1], ceil(params.S/2), 1);
satellite_angular_velocities = satellite_angular_velocities(1:params.S) * 2*pi / 100;
```

## 预期效果

相比原始版本：

| 指标 | 变化 | 原因 |
|------|------|------|
| **卫星可用性** | 大幅波动 | 卫星位置改变导致距离和信噪比变化 |
| **学习难度** | 增加 | MAB 需要适应动态变化的环境 |
| **最优臂变化** | 频繁 | 不同时刻最优的卫星组合不同 |
| **性能差异** | 可能更明显 | 好的 MAB 会学会追踪变化的最优臂 |

## 验证运动

运行脚本后，可以通过以下方式验证卫星确实在移动：

```matlab
% 在主循环中添加调试输出（注释掉或启用）
if mod(t, 50) == 1  % 每 50 步打印一次
    ican.logf(params, "info", "t=%d: Sat 1 位置 = (%.1f, %.1f, %.1f) km", ...
        t, scenario.pSat(1,1)/1e3, scenario.pSat(1,2)/1e3, scenario.pSat(1,3)/1e3);
end
```

## 进一步的改进建议

1. **加入椭圆轨道**：卫星不按圆形路径运动，而是椭圆
   ```matlab
   eccentricity = 0.1;
   current_radius = satellite_radius * (1 - eccentricity*cos(current_angle));
   ```

2. **加入高度变化**：卫星在运动时也改变高度
   ```matlab
   scenario.pSat(s, 3) = params.satHeight_m + 5e3 * sin(current_angle);  % ±5km 高度波动
   ```

3. **不同轨道**：多个卫星星座在不同轨道运动
   ```matlab
   % 分组设置不同半径
   if s <= params.S/2
       current_radius = 120e3;  % 低轨
   else
       current_radius = 140e3;  % 高轨
   end
   ```

4. **阶段性降轨/升轨**：模拟卫星寿命阶段
   ```matlab
   altitude_decay = 1e3 * (t / T_total);  % 逐步降轨
   scenario.pSat(s, 3) = params.satHeight_m - altitude_decay;
   ```

## 性能对比方法

运行两个版本对比：

1. 备份当前文件
2. 注释掉卫星更新代码（或将角速度设为 0）
3. 运行两次获得"静态"和"动态"结果
4. 对比 MAB 学习曲线和最终性能

这会清晰展示卫星运动对 MAB 学习能力的影响。
