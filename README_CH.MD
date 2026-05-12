# Satellite-based-Geolocation

基于 LEO Walker-star 星座，使用四种经典方法估计地面设备位置的端到端脚本：

* **RSS** — 接收信号强度（按历元进行偏移量剖析，无固定参考）
* **TDoA** — 到达时间差（每个历元使用参考卫星）
* **Doppler** — 多普勒频移；单颗卫星随时间观测（全局偏置）
* **AoA** — 到达角，基于卫星处的方位角/俯仰角几何关系

流程：

1. 构建/加载星座与时间；
2. 计算**逐历元可见性**；
3. 一次性为**所有**模态仿真带噪测量值；
4. 运行一个**批量**优化器，可使用**任意组合**的方法，返回单个最优拟合 LLA 点。图中显示 **真实值 vs 估计值 vs 初始值**；
5. 图中显示真实值 vs 估计值 vs 初始值。

---

## 仓库结构（扁平结构）

```text
Main.m                          % 编排脚本：运行完整流程
S101_Parameters.m               % 全局常量、时间窗口、噪声、先验/搜索窗口
S102_CreateSatGeometry.m        % satelliteScenario + Walker-star 星座 + 真实地面站
S103_GeometryAndVisibility.m    % 逐历元可见性（地面与卫星视角）+ 多普勒卫星位置/速度
S201_Measurements.m             % 仿真 RSS、TDoA、多普勒、AoA（真实值 -> 带噪值）
S301_OptimizeEstimate.m         % 联合优化器：通过 Selector 选择方法，可选权重
S401_Plot.m                     % 地理图：真实值 / 估计值 / 初始值 + 摘要文本
README.md
```

## 运行要求

* **MATLAB** R2022b+（推荐 R2025a）
* **工具箱**

  * Satellite Communications Toolbox（`satelliteScenario`、`walkerStar`、`access`、`aer`、`states`）
  * Aerospace Toolbox（`lla2ecef` 等）
  * Mapping Toolbox（用于 `geoplot` 底图）

---

## 快速开始

1. 打开 **`S00_Parameters.m`** 并检查：

   * 仿真窗口：`startTime`、`duration_min`、`sampleTime`
   * 射频与噪声：`Param.f`，IoT/设备噪声标准差（`sigmaRSS`、`sigmaT`、`sigmaf`、`sigmaAz`、`sigmaEl`）
   * 搜索框/先验：`IoTLat`、`IoTLon`、`IoTAlt`、`latMin/Max`、`lonMin/Max`
   * Walker 配置：`walker.NPlanes`、`walker.SatsPerPlane`、`walker.Inc` 等
   * 仰角掩膜：`minElDeg`（默认约 10°）

2. 运行 `Main.m`。默认情况下它会：

   * 构建几何关系；
   * 计算可见性；
   * 仿真测量值；
   * 求解最优拟合位置；
   * 绘制结果。

## 如何选择方法

在 `Main.m` 中设置 `Selector`（每种方法开/关）和可选权重：

```matlab
% [ RSS  TDoA  Doppler  AoA ]
Selector = [ 1    1      0     0 ];                 % 启用 RSS + TDoA
EstimPos = S301_OptimizeEstimate(Selector, sc, LeoSats, InitPos, ...
                                 MeasAll, Param, Const, Device, ...
                                 'Weights', [0.5 1 1 1], ...     % W_RSS, W_TDoA, W_Dop, W_AoA
                                 'FixAltitude', true, ...        % 在 IoTAlt 高度上做 2D 求解
                                 'Display','iter');              % 或 'off'
```

单一方法：将某一个元素设为 1，例如 `[1 0 0 0]` 表示仅使用 RSS。

融合：开启多种方法，并可选择调整 `'Weights'`。

## 每个脚本的作用

#### S102_CreateSatGeometry.m

创建 Walker-star 星座（`satelliteScenario`、`walkerStar`），并在 `IoTPos` 处创建真实地面站。

#### S103_GeometryAndVisibility.m

以两种方式计算方位角/俯仰角/距离历史：

地面视角（`aer(GS, sats)`）：用于 RSS、TDoA、多普勒。可见性掩膜为 `VisMaskG = ELg > minElDeg`。

卫星视角（`aer(sats, GS)`）：用于卫星处的 AoA。可见性掩膜为 `VisMaskS = ELs < -minElDeg`。

同时选择可见历元最多的一颗多普勒卫星，并通过 `states` 提取其 ECEF 位置/速度时间序列。

#### S201_Measurements.m

仅在可见位置生成带噪测量值：

RSS：`RSS = EIRP - FSPL + GRx + noise`。之后会按历元剖析掉常数项。

TDoA：由距离构造 TOA；每个历元的参考卫星为最高仰角卫星；`TDoA = TOA - TOA_ref`。

Doppler（单星）：`fD = -(f/c)·v_r + bias + noise`，其中 `v_r = dot(v_sat, u_LOS)`。

AoA：卫星处的方位角/俯仰角，加入噪声并进行角度环绕处理。

#### S301_OptimizeEstimate.m

联合代价函数 = 各模态特定 RMSE 项的加权和。选项包括：

RSS：每个历元剖析一个偏移量（抵消未知 EIRP 与常数项）。

TDoA：使用同一逐历元可见集合，计算模型值与测量值 TDoA 的 RMSE。

Doppler（单星）：剖析一个全局偏置，例如本振频偏。

AoA：仅基于方位的残差，并将方位角环绕到 `(−180°, 180°]`。

支持 2D 求解（`'FixAltitude', true`）或完整 3D 求解。

#### S401_Plot.m

绘制 `geoplot`：真实值（绿色）、估计值（红色 X）、初始值（蓝色）+ 误差文本。

## 模型片段（供参考）

自由空间路径损耗（dB）：

```text
FSPL = 20·log10(R[m]) + 20·log10(f[Hz]) − 147.55
```

RSS 模型：

```text
RSS ≈ −FSPL + GRx + K_k
```

其中 `K_k` 按历元剖析掉。

多普勒：

```text
f_D ≈ −(f/c)·v_r
v_r = dot(v_sat, u_LOS)
u_LOS = (r_sat − r_dev)/||·||
```

TDoA：

```text
TDOA = TOA − TOA_ref
TOA = range/c
```

## 提示与故障排查

可见性：如果没有卫星可见，可以增大 `duration_min`、移动 `startTime`，或降低 `minElDeg`（仅用于调试）。

初始猜测：`InitPos` 会在 `(IoTLat, IoTLon)` 附近按约 1σ 扰动，并被限制在搜索框内。搜索框越紧，收敛越容易。

高度：固定高度（`'FixAltitude', true`）通常能稳定 RSS/TDoA/Doppler 的求解。

单位：`Param.f` 单位为 Hz；距离为米；速度为 m/s；多普勒为 Hz；角度为度。

噪声真实性：先从较小的标准差开始验证几何关系，然后再扩展到你的实际场景。

权重：融合多种方法时，先从相等权重开始，再进行调节，例如当时钟较精确时提高 TDoA 权重。

## 许可证

MIT

## 致谢

基于 MATLAB Satellite Communications、Aerospace 和 Mapping 工具箱构建。欢迎贡献与提交 pull request！
