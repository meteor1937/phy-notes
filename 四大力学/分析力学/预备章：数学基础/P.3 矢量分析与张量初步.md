---
created: 2026-07-18T03:15:21+08:00
modified: 2026-07-18T13:57:00+08:00
zhihu-title: 矢量分析与张量初步
zhihu-topics:
  - 分析力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061652489871664690
zhihu-toc: true
zhihu-created-at: 2026-07-18 13:57
zhihu-updated-at: 2026-07-21 14:49
---
<!-- zhihu-nav-top:start -->
**上一节：** [P.2 变分法的推广与约束](https://zhuanlan.zhihu.com/p/2061810275808973374)
<!-- zhihu-nav-top:end -->

## P.3 矢量分析与张量初步

分析力学用广义坐标描述系统的位形，用拉格朗日量或哈密顿量概括动力学。要在这套语言中做具体计算，矢量分析和张量代数是日常工具。梯度、散度、旋度在物理空间中描述场的变化；散度定理和斯托克斯定理把体积分与面积分、面积分与线积分联系起来，是变分法中分部积分的高维推广。爱因斯坦求和约定和克罗内克 $\delta$、列维-奇维塔 $\epsilon$ 符号则把多分量的矢量恒等式化为简洁的代数操作。张量的基本概念更进一步：转动惯量、应力、应变、电磁场张量、能量-动量张量，本质上都是张量。这一节把这些数学工具从零梳理一遍，推导完整，建立几何直觉。

---

### 一、标量场与矢量场

**标量场**是空间每一点对应一个实数的函数：

$$
\phi: \mathbb{R}^3 \to \mathbb{R}, \quad \phi(\mathbf{r}) = \phi(x, y, z).
$$

例子：房间里的温度分布、电势、大气的压强。

**矢量场**是空间每一点对应一个矢量的函数：

$$
\mathbf{F}: \mathbb{R}^3 \to \mathbb{R}^3, \quad \mathbf{F}(\mathbf{r}) = (F_x(\mathbf{r}), F_y(\mathbf{r}), F_z(\mathbf{r})).
$$

例子：电场、磁场、流体的速度场、引力场。

这里的“场”字就是“函数”的意思——只是自变量是空间位置。

---

### 二、梯度、散度、旋度的几何直觉

三个一阶微分算子作用在标量场和矢量场上，各有清晰的几何意义。先给直觉，再写分量，最后建立微分恒等式。

#### 2.1 梯度（Gradient）

标量场 $\phi(\mathbf{r})$ 的梯度是矢量场，记作 $\nabla \phi$ 或 $\text{grad}\, \phi$。它的方向是 $\phi$ 增长最快的方向，大小等于该方向上的变化率。

**从方向导数出发推导梯度表达式**。设 $\mathbf{r}(s) = (x(s), y(s), z(s))$ 是一条空间曲线，$\phi$ 沿该曲线的变化率为

$$
\frac{d\phi}{ds} = \frac{\partial \phi}{\partial x}\frac{dx}{ds} + \frac{\partial \phi}{\partial y}\frac{dy}{ds} + \frac{\partial \phi}{\partial z}\frac{dz}{ds}.
$$

令 $\hat{\mathbf{u}} = \frac{d\mathbf{r}}{ds}$ 为曲线的单位切矢量。则

$$
\frac{d\phi}{ds} = \frac{\partial \phi}{\partial x}u_x + \frac{\partial \phi}{\partial y}u_y + \frac{\partial \phi}{\partial z}u_z.
$$

定义矢量

$$
\nabla \phi \equiv \left( \frac{\partial \phi}{\partial x}, \frac{\partial \phi}{\partial y}, \frac{\partial \phi}{\partial z} \right).
$$

则

$$
\frac{d\phi}{ds} = \nabla \phi \cdot \hat{\mathbf{u}}.
$$

当 $\hat{\mathbf{u}}$ 平行于 $\nabla \phi$ 时，方向导数取最大值 $|\nabla \phi|$。所以 $\nabla \phi$ 指向 $\phi$ 增加最快的方向。

在直角坐标系中：

$$
\boxed{\nabla \phi = \frac{\partial \phi}{\partial x} \hat{\mathbf{x}} + \frac{\partial \phi}{\partial y} \hat{\mathbf{y}} + \frac{\partial \phi}{\partial z} \hat{\mathbf{z}}}. \tag{P.3.1}
$$

可以把 $\nabla$ 本身看作一个矢量算子：

$$
\nabla \equiv \hat{\mathbf{x}} \frac{\partial}{\partial x} + \hat{\mathbf{y}} \frac{\partial}{\partial y} + \hat{\mathbf{z}} \frac{\partial}{\partial z}. \tag{P.3.2}
$$

它既有微分性（要求导），又有矢量性（有方向）。这两个性质在后面的运算中会反复出现。

#### 2.2 散度（Divergence）

矢量场 $\mathbf{F}(\mathbf{r})$ 的散度是标量场，记作 $\nabla \cdot \mathbf{F}$ 或 $\text{div}\, \mathbf{F}$。它的几何意义是“从一点流出的净通量”。

考虑一个以 $\mathbf{r}$ 为中心的小立方体，边长分别为 $\Delta x, \Delta y, \Delta z$。计算矢量场 $\mathbf{F}$ 穿出这六个面的净通量。

先算垂直于 $x$ 方向的两个面。在右表面（$x+\Delta x/2$），$F_x$ 向外。在左表面（$x-\Delta x/2$），$F_x$ 向内，通量为负。净贡献为

$$
\left[F_x(x+\Delta x/2, y, z) - F_x(x-\Delta x/2, y, z)\right] \Delta y \Delta z.
$$

泰勒展开 $F_x(x \pm \Delta x/2, y, z) = F_x(x,y,z) \pm \frac{\partial F_x}{\partial x} \frac{\Delta x}{2} + O(\Delta x^2)$，差值给出

$$
\frac{\partial F_x}{\partial x} \Delta x \cdot \Delta y \Delta z = \frac{\partial F_x}{\partial x} \Delta V.
$$

$y, z$ 方向类似。六个面的总净通量为

$$
\left( \frac{\partial F_x}{\partial x} + \frac{\partial F_y}{\partial y} + \frac{\partial F_z}{\partial z} \right) \Delta V.
$$

单位体积的净通量就是散度：

$$
\boxed{\nabla \cdot \mathbf{F} \equiv \frac{\partial F_x}{\partial x} + \frac{\partial F_y}{\partial y} + \frac{\partial F_z}{\partial z}}. \tag{P.3.3}
$$

把 $\nabla$ 视为矢量算子，散度就是 $\nabla$ 与 $\mathbf{F}$ 的“点积”：$\nabla \cdot \mathbf{F}$。

如果 $\nabla \cdot \mathbf{F} > 0$，该点有“源”（流体从这里流出）；$\nabla \cdot \mathbf{F} < 0$，该点有“汇”（流体流入）；$\nabla \cdot \mathbf{F} = 0$，无源无汇，场线只是穿过。静电学中高斯定律 $\nabla \cdot \mathbf{E} = \rho/\epsilon_0$ 就是说电荷是电场的源。

#### 2.3 旋度（Curl）

矢量场 $\mathbf{F}$ 的旋度是矢量场，记作 $\nabla \times \mathbf{F}$ 或 $\text{curl}\, \mathbf{F}$。它的几何意义是“局域旋转”——矢量场绕某点“打转”的趋势。

考虑一个小矩形回路。在 $xy$ 平面中取以 $(x,y)$ 为中心、边长为 $\Delta x, \Delta y$ 的矩形回路，逆时针绕行。计算 $\mathbf{F}$ 沿回路的环量 $\oint \mathbf{F} \cdot d\mathbf{r}$。

沿着四条边的贡献分别为（泰勒展开取一阶）：

- 下边（从 $x-\Delta x/2$ 到 $x+\Delta x/2$，$y-\Delta y/2$）：$+F_x(x, y-\Delta y/2, z) \Delta x$
- 右边（$x+\Delta x/2$，从下到上）：$+F_y(x+\Delta x/2, y, z) \Delta y$
- 上边（反向）：$-F_x(x, y+\Delta y/2, z) \Delta x$
- 左边（反向）：$-F_y(x-\Delta x/2, y, z) \Delta y$

合起来，$x$ 方向两次贡献的差：

$$
-[F_x(x, y+\Delta y/2, z) - F_x(x, y-\Delta y/2, z)] \Delta x \approx -\frac{\partial F_x}{\partial y} \Delta y \Delta x.
$$

$y$ 方向两次贡献的差：

$$
[F_y(x+\Delta x/2, y, z) - F_y(x-\Delta x/2, y, z)] \Delta y \approx \frac{\partial F_y}{\partial x} \Delta x \Delta y.
$$

总环量：

$$
\oint \mathbf{F} \cdot d\mathbf{r} \approx \left( \frac{\partial F_y}{\partial x} - \frac{\partial F_x}{\partial y} \right) \Delta x \Delta y.
$$

$\Delta x \Delta y$ 是回路面积。单位面积的环量定义为旋度的 $z$ 分量：

$$
(\nabla \times \mathbf{F})_z = \frac{\partial F_y}{\partial x} - \frac{\partial F_x}{\partial y}.
$$

对 $yz$ 平面和 $zx$ 平面做同样构造，得到

$$
\boxed{\nabla \times \mathbf{F} \equiv \begin{vmatrix}
\hat{\mathbf{x}} & \hat{\mathbf{y}} & \hat{\mathbf{z}} \\[2pt]
\partial_x & \partial_y & \partial_z \\[2pt]
F_x & F_y & F_z
\end{vmatrix}}. \tag{P.3.4}
$$

展开行列式：

$$
\nabla \times \mathbf{F} = \hat{\mathbf{x}}\left(\frac{\partial F_z}{\partial y} - \frac{\partial F_y}{\partial z}\right) + \hat{\mathbf{y}}\left(\frac{\partial F_x}{\partial z} - \frac{\partial F_z}{\partial x}\right) + \hat{\mathbf{z}}\left(\frac{\partial F_y}{\partial x} - \frac{\partial F_x}{\partial y}\right).
$$

旋度场本身也是矢量场，它在任一点的指向是局部旋转轴的方向（右手定则），大小是单位面积的环量。

#### 2.4 五个重要的微分恒等式

以下是矢量分析中使用频率最高的恒等式。每个都可以从分量定义直接验证。

**(i) 梯度的旋度为零**

$$
\nabla \times (\nabla \phi) = 0. \tag{P.3.5}
$$

**验证**：$(\nabla \times \nabla \phi)_x = \partial_y (\nabla \phi)_z - \partial_z (\nabla \phi)_y = \partial_y \partial_z \phi - \partial_z \partial_y \phi = 0$。混合偏导数可交换次序。

物理意义：保守力（可写为势的负梯度）的旋度为零——无旋场。

**(ii) 旋度的散度为零**

$$
\nabla \cdot (\nabla \times \mathbf{F}) = 0. \tag{P.3.6}
$$

**验证**：
$$
\nabla \cdot (\nabla \times \mathbf{F}) = \partial_x(\partial_y F_z - \partial_z F_y) + \partial_y(\partial_z F_x - \partial_x F_z) + \partial_z(\partial_x F_y - \partial_y F_x).
$$
每一项展开后成对抵消：$\partial_x\partial_y F_z$ 与 $\partial_y\partial_x F_z$ 虽不直接相消，但与 $\partial_y(-\partial_x F_z)$ 合在一起看——逐项配对：$\partial_x\partial_y F_z - \partial_y\partial_x F_z = 0$，等等。所有混合偏导数项均抵消。恒等式成立。

物理意义：磁场的散度为零（$\nabla \cdot \mathbf{B}=0$），正是因为它可以写为矢势的旋度 $\mathbf{B} = \nabla \times \mathbf{A}$。

**(iii) 旋度的旋度——拉普拉斯-梯度-散度关系**

$$
\nabla \times (\nabla \times \mathbf{F}) = \nabla(\nabla \cdot \mathbf{F}) - \nabla^2 \mathbf{F}. \tag{P.3.7}
$$

**验证**：这个推导需要仔细写。取 $x$ 分量：
$$
\begin{aligned}
(\nabla \times (\nabla \times \mathbf{F}))_x &= \partial_y(\nabla \times \mathbf{F})_z - \partial_z(\nabla \times \mathbf{F})_y \\
&= \partial_y(\partial_x F_y - \partial_y F_x) - \partial_z(\partial_z F_x - \partial_x F_z) \\
&= \partial_x\partial_y F_y - \partial_y^2 F_x - \partial_z^2 F_x + \partial_x\partial_z F_z.
\end{aligned}
$$

补上 $\partial_x^2 F_x - \partial_x^2 F_x = 0$，重新组合：

$$
\begin{aligned}
&= (\partial_x^2 F_x + \partial_x\partial_y F_y + \partial_x\partial_z F_z) - (\partial_x^2 + \partial_y^2 + \partial_z^2)F_x \\
&= \partial_x(\partial_x F_x + \partial_y F_y + \partial_z F_z) - \nabla^2 F_x \\
&= \partial_x(\nabla \cdot \mathbf{F}) - \nabla^2 F_x.
\end{aligned}
$$

$y, z$ 分量同理。公式得证。

这是电磁波方程推导中的关键一步。从麦克斯韦方程组出发：
$$
\nabla \times (\nabla \times \mathbf{E}) = \nabla \times \left(-\frac{\partial \mathbf{B}}{\partial t}\right) = -\frac{\partial}{\partial t}(\nabla \times \mathbf{B}) = -\mu_0\epsilon_0 \frac{\partial^2 \mathbf{E}}{\partial t^2}.
$$
左边用恒等式 (iii)：$\nabla(\nabla \cdot \mathbf{E}) - \nabla^2 \mathbf{E}$。真空中 $\nabla \cdot \mathbf{E}=0$，得波动方程 $\nabla^2 \mathbf{E} = \mu_0\epsilon_0 \frac{\partial^2 \mathbf{E}}{\partial t^2}$。

**(iv) 梯度的散度——拉普拉斯算子**

$$
\nabla \cdot (\nabla \phi) = \nabla^2 \phi = \frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} + \frac{\partial^2 \phi}{\partial z^2}. \tag{P.3.8}
$$

这是最基本的二阶标量微分算子。

**(v) 乘积法则**

这三个公式在分部积分中反复出现：
$$
\nabla(\phi\psi) = \phi\nabla\psi + \psi\nabla\phi, \tag{P.3.9}
$$
$$
\nabla \cdot (\phi\mathbf{F}) = \phi\nabla\cdot\mathbf{F} + \mathbf{F}\cdot\nabla\phi, \tag{P.3.10}
$$
$$
\nabla \times (\phi\mathbf{F}) = \phi\nabla\times\mathbf{F} + \nabla\phi \times \mathbf{F}. \tag{P.3.11}
$$

每个都可通过逐分量展开直接验证。

---

### 三、散度定理与斯托克斯定理

这两个定理把积分与微分联系起来，是变分法分部积分的高维形式。

#### 3.1 散度定理（高斯定理）

**定理陈述**：设 $V$ 是三维空间中的有界区域，$\partial V$ 是其闭合边界面（外法向 $\hat{\mathbf{n}}$），$\mathbf{F}$ 是光滑矢量场。则

$$
\boxed{\int_V (\nabla \cdot \mathbf{F}) \, dV = \oint_{\partial V} \mathbf{F} \cdot \hat{\mathbf{n}} \, dS}. \tag{P.3.12}
$$

**推导思路**：把 $V$ 分割成大量小立方体。对每个小立方体 $\Delta V_i$，由散度的定义，$\nabla \cdot \mathbf{F} \, \Delta V_i \approx$ 净通量。对所有小立方体求和，相邻立方体内部的共享面通量互相抵消，只剩下最外边界的面。取极限 $\Delta V_i \to 0$，和变为积分，即得散度定理。

严格证明需要将 $V$ 投影到三个坐标方向分别处理，用微积分基本定理。以 $z$ 方向为例：设 $V$ 在 $xy$ 平面的投影区域为 $D$，$V$ 由下曲面 $z=z_1(x,y)$ 和上曲面 $z=z_2(x,y)$ 界定。则

$$
\int_V \frac{\partial F_z}{\partial z} dV = \iint_D \left[ \int_{z_1(x,y)}^{z_2(x,y)} \frac{\partial F_z}{\partial z} dz \right] dxdy = \iint_D \left[ F_z(x,y,z_2) - F_z(x,y,z_1) \right] dxdy.
$$

上表面的外法向 $z$ 分量为 $n_z \,dS = dxdy$；下表面的外法向向下，$n_z\,dS = -dxdy$。因此上式右侧正是 $\oint_{\partial V} F_z n_z dS$。$x, y$ 方向同理，三者相加即得定理。

**物理含义**：区域内“源”的总和等于穿出边界的总通量。高斯定律 $\oint \mathbf{E}\cdot d\mathbf{S} = Q/\epsilon_0$ 就是这个定理的直接应用。

#### 3.2 斯托克斯定理

**定理陈述**：设 $S$ 是光滑曲面，边界 $\partial S$ 是闭合曲线（方向与曲面法向满足右手定则），$\mathbf{F}$ 是光滑矢量场。则

$$
\boxed{\int_S (\nabla \times \mathbf{F}) \cdot \hat{\mathbf{n}} \, dS = \oint_{\partial S} \mathbf{F} \cdot d\mathbf{r}}. \tag{P.3.13}
$$

**推导思路**：把曲面 $S$ 分成大量小面元 $\Delta S_i$。对每个面元，由旋度的定义，$(\nabla \times \mathbf{F})\cdot \hat{\mathbf{n}} \, \Delta S_i \approx$ 绕该面元边界的环量。对所有面元求和，相邻面元的共享边上线积分方向相反，互相抵消，只剩下最外围边界 $\partial S$ 的线积分。取极限即得定理。

**物理含义**：法拉第定律 $\oint \mathbf{E}\cdot d\mathbf{r} = -\frac{d}{dt}\int \mathbf{B}\cdot d\mathbf{S}$ 是斯托克斯定理的直接体现。保守力满足 $\nabla \times \mathbf{F} = 0$，则由斯托克斯定理，沿任意闭合回路的积分为零——与路径无关。

#### 3.3 分部积分公式的矢量形式

散度定理可以直接给出变分法中所需的矢量分部积分公式。取 $\mathbf{F} = \phi \nabla \psi$，则

$$
\nabla \cdot \mathbf{F} = \nabla \cdot (\phi \nabla \psi) = \phi \nabla^2 \psi + \nabla\phi \cdot \nabla\psi.
$$

散度定理给出

$$
\int_V (\phi \nabla^2 \psi + \nabla\phi \cdot \nabla\psi) dV = \oint_{\partial V} \phi \nabla\psi \cdot \hat{\mathbf{n}} \, dS,
$$

移项得

$$
\int_V \phi \nabla^2 \psi \, dV = \oint_{\partial V} \phi \frac{\partial\psi}{\partial n} dS - \int_V \nabla\phi \cdot \nabla\psi \, dV. \tag{P.3.14}
$$

这是**格林第一恒等式**。交换 $\phi, \psi$ 并相减可得**格林第二恒等式**：

$$
\int_V (\phi \nabla^2 \psi - \psi \nabla^2 \phi) dV = \oint_{\partial V} \left( \phi \frac{\partial\psi}{\partial n} - \psi \frac{\partial\phi}{\partial n} \right) dS. \tag{P.3.15}
$$

这些公式在变分法求解拉普拉斯方程和泊松方程时是基本工具。

---

### 四、曲线坐标中的微分算符

很多物理问题天然适合球坐标或柱坐标。在这些坐标中，梯度、散度、旋度的表达式需要重新推导。

#### 4.1 一般正交曲线坐标

设坐标变换 $(u_1, u_2, u_3)$ 与直角坐标 $(x, y, z)$ 之间一一对应。局部基矢量：

$$
\mathbf{e}_i = \frac{\partial \mathbf{r}}{\partial u_i}, \quad \text{其长度为 } h_i = \left| \frac{\partial \mathbf{r}}{\partial u_i} \right|.
$$

$h_i$ 称为**拉梅系数**（Lamé coefficients）。若三个基矢量相互正交（大多数物理坐标如此），则单位正交基为 $\hat{\mathbf{e}}_i = \frac{1}{h_i} \frac{\partial \mathbf{r}}{\partial u_i}$。线元为

$$
d\mathbf{r} = h_1 du_1 \hat{\mathbf{e}}_1 + h_2 du_2 \hat{\mathbf{e}}_2 + h_3 du_3 \hat{\mathbf{e}}_3.
$$

面积元和体积元分别为

$$
dS_{ij} = h_i h_j \, du_i du_j, \quad dV = h_1 h_2 h_3 \, du_1 du_2 du_3.
$$

在此坐标下的三个基本算子：

**梯度**：由方向导数的定义，$\nabla \phi \cdot d\mathbf{r} = d\phi = \sum \frac{\partial \phi}{\partial u_i} du_i$。又 $\nabla \phi \cdot d\mathbf{r} = \sum (\nabla \phi)_i h_i du_i$。比较得 $(\nabla \phi)_i = \frac{1}{h_i} \frac{\partial \phi}{\partial u_i}$。因此

$$
\boxed{\nabla \phi = \frac{1}{h_1}\frac{\partial \phi}{\partial u_1}\hat{\mathbf{e}}_1 + \frac{1}{h_2}\frac{\partial \phi}{\partial u_2}\hat{\mathbf{e}}_2 + \frac{1}{h_3}\frac{\partial \phi}{\partial u_3}\hat{\mathbf{e}}_3}. \tag{P.3.16}
$$

**散度**：用散度定理的微元分析。考虑以 $(u_1,u_2,u_3)$ 为中心的小体积元。穿过 $u_1$ 为常数的两个面的净通量在 $\Delta u_1 \to 0$ 时为

$$
\frac{\partial}{\partial u_1}(F_1 h_2 h_3) \Delta u_1 \Delta u_2 \Delta u_3,
$$

其中 $F_1 = \mathbf{F} \cdot \hat{\mathbf{e}}_1$。三个方向相加，除以体积元 $h_1h_2h_3\Delta u_1\Delta u_2\Delta u_3$：

$$
\boxed{\nabla \cdot \mathbf{F} = \frac{1}{h_1 h_2 h_3}\left[ \frac{\partial}{\partial u_1}(F_1 h_2 h_3) + \frac{\partial}{\partial u_2}(F_2 h_1 h_3) + \frac{\partial}{\partial u_3}(F_3 h_1 h_2) \right]}. \tag{P.3.17}
$$

**旋度**：用斯托克斯定理的微元分析，或直接类比：

$$
\boxed{\nabla \times \mathbf{F} = \frac{1}{h_1 h_2 h_3}\begin{vmatrix}
h_1\hat{\mathbf{e}}_1 & h_2\hat{\mathbf{e}}_2 & h_3\hat{\mathbf{e}}_3 \\[2pt]
\frac{\partial}{\partial u_1} & \frac{\partial}{\partial u_2} & \frac{\partial}{\partial u_3} \\[2pt]
h_1 F_1 & h_2 F_2 & h_3 F_3
\end{vmatrix}}. \tag{P.3.18}
$$

**拉普拉斯算子**：$\nabla^2 \phi = \nabla \cdot (\nabla \phi)$，将梯度分量代入散度公式：

$$
\boxed{\nabla^2 \phi = \frac{1}{h_1 h_2 h_3}\left[ \frac{\partial}{\partial u_1}\left(\frac{h_2 h_3}{h_1}\frac{\partial \phi}{\partial u_1}\right) + \frac{\partial}{\partial u_2}\left(\frac{h_1 h_3}{h_2}\frac{\partial \phi}{\partial u_2}\right) + \frac{\partial}{\partial u_3}\left(\frac{h_1 h_2}{h_3}\frac{\partial \phi}{\partial u_3}\right) \right]}. \tag{P.3.19}
$$

#### 4.2 柱坐标 $(r, \phi, z)$

坐标变换：$x = r\cos\phi, y = r\sin\phi, z = z$。

$$
h_r = \left|\frac{\partial \mathbf{r}}{\partial r}\right| = |(\cos\phi, \sin\phi, 0)| = 1,
$$
$$
h_\phi = \left|\frac{\partial \mathbf{r}}{\partial \phi}\right| = |(-r\sin\phi, r\cos\phi, 0)| = r,
$$
$$
h_z = 1.
$$

代入一般公式：

**梯度**：
$$
\nabla \psi = \frac{\partial \psi}{\partial r} \hat{\mathbf{r}} + \frac{1}{r}\frac{\partial \psi}{\partial \phi} \hat{\boldsymbol{\phi}} + \frac{\partial \psi}{\partial z} \hat{\mathbf{z}}. \tag{P.3.20}
$$

**散度**（$\mathbf{F} = F_r \hat{\mathbf{r}} + F_\phi \hat{\boldsymbol{\phi}} + F_z \hat{\mathbf{z}}$）：
$$
\nabla \cdot \mathbf{F} = \frac{1}{r}\frac{\partial}{\partial r}(r F_r) + \frac{1}{r}\frac{\partial F_\phi}{\partial \phi} + \frac{\partial F_z}{\partial z}. \tag{P.3.21}
$$

**旋度**：
$$
\nabla \times \mathbf{F} = \left(\frac{1}{r}\frac{\partial F_z}{\partial \phi} - \frac{\partial F_\phi}{\partial z}\right)\hat{\mathbf{r}} + \left(\frac{\partial F_r}{\partial z} - \frac{\partial F_z}{\partial r}\right)\hat{\boldsymbol{\phi}} + \frac{1}{r}\left(\frac{\partial}{\partial r}(r F_\phi) - \frac{\partial F_r}{\partial \phi}\right)\hat{\mathbf{z}}. \tag{P.3.22}
$$

**拉普拉斯算子**：
$$
\nabla^2 \psi = \frac{1}{r}\frac{\partial}{\partial r}\left(r \frac{\partial \psi}{\partial r}\right) + \frac{1}{r^2}\frac{\partial^2 \psi}{\partial \phi^2} + \frac{\partial^2 \psi}{\partial z^2}. \tag{P.3.23}
$$

#### 4.3 球坐标 $(r, \theta, \phi)$

坐标变换：$x = r\sin\theta\cos\phi, y = r\sin\theta\sin\phi, z = r\cos\theta$。

$$
h_r = 1, \quad h_\theta = r, \quad h_\phi = r\sin\theta.
$$

**梯度**：
$$
\nabla \psi = \frac{\partial \psi}{\partial r} \hat{\mathbf{r}} + \frac{1}{r}\frac{\partial \psi}{\partial \theta} \hat{\boldsymbol{\theta}} + \frac{1}{r\sin\theta}\frac{\partial \psi}{\partial \phi} \hat{\boldsymbol{\phi}}. \tag{P.3.24}
$$

**散度**（$\mathbf{F} = F_r \hat{\mathbf{r}} + F_\theta \hat{\boldsymbol{\theta}} + F_\phi \hat{\boldsymbol{\phi}}$）：
$$
\nabla \cdot \mathbf{F} = \frac{1}{r^2}\frac{\partial}{\partial r}(r^2 F_r) + \frac{1}{r\sin\theta}\frac{\partial}{\partial \theta}(\sin\theta \, F_\theta) + \frac{1}{r\sin\theta}\frac{\partial F_\phi}{\partial \phi}. \tag{P.3.25}
$$

**旋度**：
$$
\nabla \times \mathbf{F} = \frac{1}{r\sin\theta}\left[\frac{\partial}{\partial \theta}(\sin\theta F_\phi) - \frac{\partial F_\theta}{\partial \phi}\right]\hat{\mathbf{r}} + \frac{1}{r}\left[\frac{1}{\sin\theta}\frac{\partial F_r}{\partial \phi} - \frac{\partial}{\partial r}(r F_\phi)\right]\hat{\boldsymbol{\theta}} + \frac{1}{r}\left[\frac{\partial}{\partial r}(r F_\theta) - \frac{\partial F_r}{\partial \theta}\right]\hat{\boldsymbol{\phi}}. \tag{P.3.26}
$$

**拉普拉斯算子**：
$$
\nabla^2 \psi = \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2 \frac{\partial \psi}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial \theta}\left(\sin\theta \frac{\partial \psi}{\partial \theta}\right) + \frac{1}{r^2\sin^2\theta}\frac{\partial^2 \psi}{\partial \phi^2}. \tag{P.3.27}
$$

---

### 五、爱因斯坦求和约定

多个指标相乘相加在物理中极其常见。爱因斯坦提出：凡是乘积中重复出现的指标，自动对该指标求和，省略求和号。这就是**爱因斯坦求和约定**。

**基本规则**：
- 重复一次的指标称为“哑指标”（dummy index），自动求和：$a_i b_i \equiv \sum_{i=1}^3 a_i b_i$。
- 不重复的指标称为“自由指标”（free index），等式两边自由指标必须匹配。
- 同一个哑指标不能出现超过两次（若需要更多重复，必须显式写出求和号）。
- 哑指标的符号可以任意替换：$a_i b_i = a_j b_j$。

**示例**：

点积：$\mathbf{a} \cdot \mathbf{b} = a_i b_i$。

矩阵乘法：$(AB)_{ij} = A_{ik} B_{kj}$。（$k$ 求和，$i, j$ 自由）

矩阵的迹：$\text{Tr}(A) = A_{ii}$。

矢量场散度：$\nabla \cdot \mathbf{F} = \partial_i F_i$（其中 $\partial_i \equiv \partial/\partial x_i$）。

矢量的第 $i$ 分量：$(\nabla \phi)_i = \partial_i \phi$。

梯度的散度（拉普拉斯算子）：$\nabla^2 \phi = \partial_i \partial_i \phi$。

旋度的第 $i$ 分量：$(\nabla \times \mathbf{F})_i = \epsilon_{ijk} \partial_j F_k$（$\epsilon$ 符号见下节）。

叉积：$(\mathbf{a} \times \mathbf{b})_i = \epsilon_{ijk} a_j b_k$。

---

### 六、克罗内克 $\delta_{ij}$ 与列维-奇维塔 $\epsilon_{ijk}$ 的恒等式

#### 6.1 克罗内克 $\delta$ 符号

定义：

$$
\delta_{ij} = \begin{cases}
1, & i = j, \\
0, & i \neq j.
\end{cases} \tag{P.3.28}
$$

在代数上，$\delta_{ij}$ 就是单位矩阵的 $(i,j)$ 元素。

**基本性质**：

1. 替换指标：$a_j \delta_{ij} = a_i$。$\delta$ 把指标 $j$ “替换”为 $i$。
2. $\delta_{ii} = 3$（三维空间中迹为 3）。
3. $\delta_{ij}\delta_{jk} = \delta_{ik}$（矩阵乘法的单位元性质）。

#### 6.2 列维-奇维塔 $\epsilon$ 符号

$\epsilon_{ijk}$ 是三维全反对称张量，定义：

$$
\epsilon_{ijk} = \begin{cases}
+1, & (i,j,k) \text{ 是 } (1,2,3) \text{ 的偶排列}, \\
-1, & (i,j,k) \text{ 是 } (1,2,3) \text{ 的奇排列}, \\
0, & \text{任意两个指标相等}.
\end{cases} \tag{P.3.29}
$$

明确写出非零分量：
$$
\epsilon_{123} = \epsilon_{231} = \epsilon_{312} = +1,
$$
$$
\epsilon_{132} = \epsilon_{213} = \epsilon_{321} = -1.
$$

**基本性质**：

1. 全反对称：$\epsilon_{ijk} = -\epsilon_{jik} = -\epsilon_{ikj} = -\epsilon_{kji}$。
2. 交换任意两个指标变号，交换两次恢复原号。
3. $\epsilon_{ijk}$ 给出叉积的分量：$(\mathbf{a} \times \mathbf{b})_i = \epsilon_{ijk} a_j b_k$。
4. 旋度：$(\nabla \times \mathbf{F})_i = \epsilon_{ijk} \partial_j F_k$。
5. 行列式：对于 $3\times 3$ 矩阵 $A$，$\det A = \epsilon_{ijk} A_{1i} A_{2j} A_{3k}$。

#### 6.3 $\epsilon$-$\delta$ 恒等式（最重要）

三维矢量运算中最重要的恒等式是

$$
\boxed{\epsilon_{ijk} \epsilon_{ilm} = \delta_{jl}\delta_{km} - \delta_{jm}\delta_{kl}}. \tag{P.3.30}
$$

**完整推导**：

考虑 $\epsilon_{ijk} \epsilon_{ilm}$。这里 $i$ 是哑指标（求和），$j,k,l,m$ 是自由指标。我们先证明这个恒等式，然后展示它的威力。

**第一步**：构造 $3\times 3$ 矩阵 $M$，使其行列式可以表达这个乘积。注意到 $\epsilon_{ijk}$ 可以理解为：给定 $j,k$ 后唯一剩下一个 $i$ 使 $\epsilon_{ijk} \neq 0$（在三维中）。精确地说，$\sum_i \epsilon_{ijk} \epsilon_{ilm}$ 这个量关于 $(j,k)$ 反对称，关于 $(l,m)$ 反对称，且满足归一化条件。而 $\delta_{jl}\delta_{km} - \delta_{jm}\delta_{kl}$ 也恰好满足完全相同的对称性和归一化。根据唯一性，二者相等。但这样论证需要确认唯一性。我们给出一个更直接的验证。

**第二步**：逐情况验证。由于 $\epsilon_{ijk}$ 有三个指标，$\epsilon_{ijk}\epsilon_{ilm}$ 的自由指标是 $j,k,l,m$，每个可取 $1,2,3$。总共有 $3^4=81$ 种组合。利用对称性简化：

- 如果 $j=k$，左边 $\epsilon_{ijk}=0$（全反对称张量有两个相同指标），所以左边 $=0$。右边：$\delta_{jj}\delta_{km} - \delta_{jm}\delta_{kj} = 3\delta_{km} - \delta_{km} = 2\delta_{km} \neq 0$？等等，这里 $j$ 固定，不是对所有 $j$ 求和！注意 $\epsilon_{ijk}\epsilon_{ilm}$ 中 $i$ 求和但 $j,k,l,m$ 是自由指标，不求和。当 $j=k$ 时，$\epsilon_{ijk} = 0$，所以左边 $=0$。右边 $\delta_{jl}\delta_{km} - \delta_{jm}\delta_{kl}$：若 $j=k$，则第一项 $\delta_{jl}\delta_{jm}$，第二项 $\delta_{jm}\delta_{jl}$，它们相等，相减为零。正确。

- 如果 $j \neq k$，$l=j, m=k$ 则 $\epsilon_{ijk}\epsilon_{ijk} = \sum_i (\epsilon_{ijk})^2$。对于固定的 $j\neq k$，只有一个 $i$ 使 $\epsilon_{ijk} \neq 0$（它是剩下的那个不重复的指标），且 $\epsilon_{ijk}^2=1$。所以和 $=1$。右边：$\delta_{jj}\delta_{kk} - \delta_{jk}\delta_{kj} = 1\cdot 1 - 0 = 1$。正确。

- 如果 $j \neq k$，$l=k, m=j$ 则 $\epsilon_{ijk}\epsilon_{ikj} = -\epsilon_{ijk}\epsilon_{ijk} = -1$（反对称）。右边：$\delta_{jk}\delta_{kj} - \delta_{jj}\delta_{kk} = 0 - 1 = -1$。正确。

- 其余情况通过交换对称性归入上述情形。恒等式成立。

#### 6.4 $\epsilon$-$\delta$ 恒等式的推论

**推论一：两个 $\epsilon$ 的一次求和（缩并一个指标）**

将恒等式中 $j=l$ 并求和（自动进行，因为 $l=j$）：

$$
\epsilon_{ijk} \epsilon_{ijm} = \delta_{jj}\delta_{km} - \delta_{jm}\delta_{kj}.
$$

注意 $\delta_{jj}=3$，$\delta_{jm}\delta_{kj} = \delta_{km}$（当 $j$ 求和时）。所以

$$
\boxed{\epsilon_{ijk} \epsilon_{ijm} = 3\delta_{km} - \delta_{km} = 2\delta_{km}}. \tag{P.3.31}
$$

**推论二：两个 $\epsilon$ 的两次求和（缩并两个指标）**

取 $k=m$ 并求和：

$$
\epsilon_{ijk} \epsilon_{ijk} = 2\delta_{kk} = 2 \times 3 = 6. \tag{P.3.32}
$$

这是 $\epsilon$ 符号所有分量平方和。

#### 6.5 $\epsilon$-$\delta$ 恒等式在矢量运算中的应用

这些恒等式让复杂的矢量恒等式验证变得机械化。

**例一：证明 $\mathbf{a} \times (\mathbf{b} \times \mathbf{c}) = \mathbf{b}(\mathbf{a}\cdot\mathbf{c}) - \mathbf{c}(\mathbf{a}\cdot\mathbf{b})$ **

写第 $i$ 分量：
$$
[\mathbf{a} \times (\mathbf{b} \times \mathbf{c})]_i = \epsilon_{ijk} a_j (\mathbf{b}\times\mathbf{c})_k = \epsilon_{ijk} a_j \epsilon_{klm} b_l c_m.
$$

$\epsilon_{ijk} = \epsilon_{kij}$（偶排列）。把两个 $\epsilon$ 写在一起：
$$
= \epsilon_{kij}\epsilon_{klm} a_j b_l c_m.
$$

用恒等式 (P.3.30)（注意哑指标是 $k$）：
$$
\epsilon_{kij}\epsilon_{klm} = \delta_{il}\delta_{jm} - \delta_{im}\delta_{jl}.
$$

因此
$$
\begin{aligned}
[\mathbf{a} \times (\mathbf{b} \times \mathbf{c})]_i &= (\delta_{il}\delta_{jm} - \delta_{im}\delta_{jl}) a_j b_l c_m \\
&= \delta_{il} b_l \, \delta_{jm} a_j c_m - \delta_{im} c_m \, \delta_{jl} a_j b_l \\
&= b_i (a_j c_j) - c_i (a_j b_j) \\
&= b_i (\mathbf{a}\cdot\mathbf{c}) - c_i (\mathbf{a}\cdot\mathbf{b}).
\end{aligned}
$$

回到矢量形式：$\mathbf{a} \times (\mathbf{b} \times \mathbf{c}) = \mathbf{b}(\mathbf{a}\cdot\mathbf{c}) - \mathbf{c}(\mathbf{a}\cdot\mathbf{b})$。证毕。整个过程没有任何几何直觉，纯代数操作。

**例二：证明旋度的旋度恒等式**

$[\nabla \times (\nabla \times \mathbf{F})]_i = \epsilon_{ijk} \partial_j (\nabla \times \mathbf{F})_k = \epsilon_{ijk} \partial_j (\epsilon_{klm} \partial_l F_m)$。

$\epsilon_{ijk}\epsilon_{klm} = \epsilon_{kij}\epsilon_{klm} = \delta_{il}\delta_{jm} - \delta_{im}\delta_{jl}$。

因此
$$
\begin{aligned}
[\nabla \times (\nabla \times \mathbf{F})]_i &= (\delta_{il}\delta_{jm} - \delta_{im}\delta_{jl}) \partial_j \partial_l F_m \\
&= \partial_i \partial_m F_m - \partial_j \partial_j F_i \\
&= \partial_i(\nabla \cdot \mathbf{F}) - \nabla^2 F_i.
\end{aligned}
$$

矢量形式：$\nabla \times (\nabla \times \mathbf{F}) = \nabla(\nabla \cdot \mathbf{F}) - \nabla^2 \mathbf{F}$。比之前的分量推导简洁得多。

**例三：混合积的恒等式**

$\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c}) = a_i (\mathbf{b}\times\mathbf{c})_i = a_i \epsilon_{ijk} b_j c_k = \epsilon_{ijk} a_i b_j c_k$。

轮换性质来源于 $\epsilon_{ijk}$ 的全反对称：$\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c}) = \mathbf{b} \cdot (\mathbf{c} \times \mathbf{a}) = \mathbf{c} \cdot (\mathbf{a} \times \mathbf{b})$。

---

### 七、张量的基本概念与坐标变换性质

矢量和标量不足以描述所有物理量。转动惯量把一个方向的角速度映射为另一个方向的角动量（$\mathbf{L} = \mathbf{I}\boldsymbol{\omega}$），应力把一个方向的法向映射为力。这些量需要张量。

#### 7.1 从矢量到张量——坐标变换的视角

选择坐标是为了计算方便，物理定律不应依赖于坐标选择。一个好的物理量定义必须包含明确的坐标变换规则。

**标量** $\phi$：在不同坐标系中数值相同。若坐标变换 $\mathbf{r} \to \mathbf{r}'$，则 $\phi'(\mathbf{r}') = \phi(\mathbf{r})$。

**矢量** $\mathbf{v}$（以位置矢量为例）：在直角坐标系中，坐标变换 $x'_i = R_{ij} x_j$（$R$ 是旋转矩阵，正交矩阵 $R^T R = I$），矢量分量按同样规则变换：
$$
v'_i = R_{ij} v_j. \tag{P.3.33}
$$
这就是矢量的定义性特征：它的分量按坐标分量的变换方式变换。

**协变与逆变**：在一般曲线坐标中（不限于直角坐标），需要区分两类矢量。设坐标变换 $x'^i = x'^i(x^1,x^2,x^3)$。

- **逆变矢量**（contravariant）：分量按坐标微分的变换规则变换，$v'^i = \frac{\partial x'^i}{\partial x^j} v^j$。指标在上标位置。位置微元 $dx^i$ 是逆变矢量。
- **协变矢量**（covariant）：分量按梯度的变换规则变换，$w'_i = \frac{\partial x^j}{\partial x'^i} w_j$。指标在下标位置。标量场的梯度 $\partial_i \phi$ 是协变矢量。

在直角坐标中（只涉及旋转），协变和逆变分量一致，但在弯曲空间或曲线坐标中二者不同。分析力学主要使用直角坐标或直接处理广义坐标，所以这里不深入协变与逆变的区别。重要的是理解张量通过坐标变换来定义这一思想。

#### 7.2 张量的定义

**二阶张量**是一个有 $3\times 3=9$ 个分量的量 $T_{ij}$，在坐标变换 $x'_i = R_{ij} x_j$ 下按

$$
T'_{ij} = R_{ik} R_{jl} T_{kl} \tag{P.3.34}
$$

变换。每个指标像矢量分量一样各乘一个旋转矩阵。推广到 $n$ 阶张量：

$$
T'_{i_1 i_2 \ldots i_n} = R_{i_1 j_1} R_{i_2 j_2} \cdots R_{i_n j_n} T_{j_1 j_2 \ldots j_n}. \tag{P.3.35}
$$

**为什么这个定义有用？** 因为它保证了张量方程在坐标变换下形式不变。如果 $T_{ij} = S_{ij}$ 在一个坐标系中成立，乘上 $R$ 矩阵后在新坐标系中也成立。这就是物理定律协变性要求的数学结构。

**重要例子**：

1. **克罗内克 $\delta_{ij}$ 是二阶各向同性张量**：$\delta'_{ij} = R_{ik}R_{jl}\delta_{kl} = R_{ik}R_{jk} = (RR^T)_{ij} = \delta_{ij}$。在任何旋转坐标中分量不变。

2. **列维-奇维塔 $\epsilon_{ijk}$ 是三阶全反对称张量**：
   $$
   \epsilon'_{ijk} = R_{il}R_{jm}R_{kn} \epsilon_{lmn} = (\det R) \epsilon_{ijk}.
   $$
   对于正常旋转 $\det R = +1$，$\epsilon_{ijk}$ 分量不变。它是各向同性张量（在正常旋转下）。

3. **并矢**：由两个矢量构成 $T_{ij} = u_i v_j$。验证变换规则：
   $$
   T'_{ij} = u'_i v'_j = (R_{ik}u_k)(R_{jl}v_l) = R_{ik}R_{jl} u_k v_l = R_{ik}R_{jl} T_{kl}.
   $$
   确实是二阶张量。

任何二阶张量都可以写为并矢的线性组合。

#### 7.3 张量的代数运算

**加法**：$C_{ij} = A_{ij} + B_{ij}$（同阶张量，逐分量相加）。

**数乘**：$B_{ij} = \lambda A_{ij}$。

**张量积**（外积）：$T_{ijkl} = A_{ij} B_{kl}$。阶数相加（$2+2=4$）。

**缩并**（contraction）：对两个指标求和，把张量的阶数降低 2。例如
$$
\text{二阶张量的迹： } T_{ii} \text{（标量，0 阶张量）}.
$$
$$
\text{三阶张量缩并： } S_{ijj} \text{（矢量，1 阶张量）}.
$$

**内积**：先做外积再做缩并。例如矩阵乘法：$(AB)_{ij} = A_{ik}B_{kj}$。这是两个二阶张量的内积，结果还是二阶张量。

**对称化与反对称化**：
$$
T_{(ij)} = \frac{1}{2}(T_{ij} + T_{ji}) \quad \text{（对称部分）},
$$
$$
T_{[ij]} = \frac{1}{2}(T_{ij} - T_{ji}) \quad \text{（反对称部分）}.
$$

任何二阶张量可唯一分解为对称和反对称部分之和：$T_{ij} = T_{(ij)} + T_{[ij]}$。

#### 7.4 二阶张量的物理实例

**转动惯量张量**：角动量与角速度的关系 $\mathbf{L} = \mathbf{I} \boldsymbol{\omega}$，分量形式 $L_i = I_{ij} \omega_j$。$I_{ij}$ 是对称张量（$I_{ij}=I_{ji}$）。在主轴坐标下，$I_{ij}$ 对角化：$I_{ij} = I_i \delta_{ij}$（不对 $i$ 求和）。

**应力张量**：作用在法向为 $\hat{\mathbf{n}}$ 的面元上的力为 $\mathbf{F}_i = \sigma_{ij} n_j$（对 $j$ 求和）。$\sigma_{ij}$ 是对称张量（角动量守恒的要求）。

**电磁场张量**：在相对论电动力学中，电场和磁场统一为一个反对称二阶张量 $F_{\mu\nu}$（$\mu,\nu=0,1,2,3$），$F_{0i} = E_i/c$，$F_{ij} = -\epsilon_{ijk} B_k$。

#### 7.5 张量的商法则

如果一个量 $T_{ijk\ldots}$ 与任意矢量（或张量）内积后得到张量，则 $T$ 本身是张量。例如：如果对任意矢量 $v_j$，$T_{ij} v_j$ 是矢量（即按矢量规则变换），则 $T_{ij}$ 是二阶张量。

**证明思路**：已知 $T'_{ij} v'_j = R_{ik} (T_{kl} v_l)$。又 $v'_j = R_{jl} v_l$，所以 $T'_{ij} R_{jl} v_l = R_{ik} T_{kl} v_l$。两边同乘 $R_{jm}$ 并对 $j$ 求和，用 $R_{jl}R_{jm} = \delta_{lm}$ 得 $T'_{im} = R_{ik} R_{ml} T_{kl}$（此处 $m,l$ 是哑指标的记号调整）。这就是张量变换规则。这一法则在判断一个物理量是否为张量时非常有用。

---

### 八、总结

P.3 节建立了分析力学必需的矢量分析和张量代数基础：

1. **三个微分算子**：梯度 $\nabla\phi$ 指向标量场增长最快的方向；散度 $\nabla\cdot\mathbf{F}$ 测量矢量场的源汇；旋度 $\nabla\times\mathbf{F}$ 测量矢量场的局部旋转。

2. **五个基本恒等式**：$\nabla\times\nabla\phi=0$，$\nabla\cdot(\nabla\times\mathbf{F})=0$，$\nabla\times(\nabla\times\mathbf{F})=\nabla(\nabla\cdot\mathbf{F})-\nabla^2\mathbf{F}$，以及乘积法则。

3. **两个积分定理**：散度定理 $\int_V \nabla\cdot\mathbf{F}\,dV = \oint_{\partial V} \mathbf{F}\cdot\hat{\mathbf{n}}\,dS$，斯托克斯定理 $\int_S \nabla\times\mathbf{F}\cdot\hat{\mathbf{n}}\,dS = \oint_{\partial S}\mathbf{F}\cdot d\mathbf{r}$。它们是变分法中分部积分的高维推广。

4. **曲线坐标**：通过拉梅系数 $h_i$ 给出梯度、散度、旋度、拉普拉斯算子在柱坐标和球坐标中的完整表达式。

5. **爱因斯坦求和约定**：重复指标自动求和，$a_i b_i = \sum_i a_i b_i$。让多分量计算极为简洁。

6. ** $\delta_{ij}$ 和 $\epsilon_{ijk}$ **：$\delta_{ij}$ 是单位张量，$\epsilon_{ijk}$ 是全反对称张量。核心恒等式 $\epsilon_{ijk}\epsilon_{ilm} = \delta_{jl}\delta_{km} - \delta_{jm}\delta_{kl}$ 让所有矢量恒等式可以代数推导。

7. **张量的定义**：按坐标变换规则定义——标量不变，矢量 $v'_i = R_{ij}v_j$，二阶张量 $T'_{ij} = R_{ik}R_{jl}T_{kl}$。这个定义保证了物理方程在不同坐标系中的协变性。

有了这些工具，我们现在可以用最干净的语言处理分析力学中的矢量运算：拉格朗日量中的动能项、刚体的转动惯量张量、电磁场中的洛伦兹力、诺特定理中的守恒流，都可以写得紧凑而精确。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.4 常微分方程与特殊函数](https://zhuanlan.zhihu.com/p/2061652524868821844)

**课程目录：** [分析力学目录](https://zhuanlan.zhihu.com/p/2061652577167717816)
<!-- zhihu-nav-bottom:end -->
