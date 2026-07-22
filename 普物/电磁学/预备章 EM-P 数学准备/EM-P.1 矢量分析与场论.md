---
created: 2026-07-20T19:25:10+08:00
modified: 2026-07-20T20:13:25+08:00
zhihu-title: EM-P.1 矢量分析与场论
zhihu-topics:
  - 电磁学
zhihu-link: https://zhuanlan.zhihu.com/p/2062629242194679720
zhihu-toc: true
zhihu-created-at: 2026-07-20 21:32
zhihu-updated-at: 2026-07-21 14:41
---
## EM-P.1 矢量分析与场论

电磁学中，电场和磁场都是空间位置的矢量函数，描述它们的规律需要一套分析矢量场的数学语言。标量场描述空间中各点的数值（如电势 $V(x,y,z)$），矢量场描述空间中各点的矢量（如电场 $\mathbf{E}(x,y,z)$）。这一节从矢量代数出发，逐步建立标量场的梯度、矢量场的散度和旋度的定义与计算方法，证明高斯散度定理和斯托克斯定理这两个将体积分与面积分、面积分与线积分联系起来的核心公式，引入拉普拉斯算子，阐述亥姆霍兹定理这一矢量分析的基本定理，最后给出柱坐标系和球坐标系中这些微分算子的表达式。

---

### 一、矢量代数基础

#### 1.1 矢量的表示与基本运算

三维空间中的矢量 $\mathbf{A}$ 可以用其沿坐标轴的分量表示。在直角坐标系 $(x,y,z)$ 中，取单位基矢量 $\hat{\mathbf{x}}, \hat{\mathbf{y}}, \hat{\mathbf{z}}$（有时记为 $\mathbf{i}, \mathbf{j}, \mathbf{k}$），则

$$
\mathbf{A} = A_x \hat{\mathbf{x}} + A_y \hat{\mathbf{y}} + A_z \hat{\mathbf{z}}
$$

矢量的模（长度）为

$$
|\mathbf{A}| = A = \sqrt{A_x^2 + A_y^2 + A_z^2}
$$

**标量积（点积）**：两个矢量 $\mathbf{A}$ 和 $\mathbf{B}$ 的点积定义为

$$
\mathbf{A} \cdot \mathbf{B} = A_x B_x + A_y B_y + A_z B_z = |\mathbf{A}| |\mathbf{B}| \cos\theta
$$

其中 $\theta$ 是两矢量间的夹角。点积满足交换律和分配律。它给出一个矢量在另一个矢量方向上的投影乘积，常用于计算功、通量等。

**矢量积（叉积）**：两个矢量 $\mathbf{A}$ 和 $\mathbf{B}$ 的叉积是一个矢量，其方向垂直于 $\mathbf{A}$ 和 $\mathbf{B}$ 所在平面，服从右手定则，大小为 $|\mathbf{A}| |\mathbf{B}| \sin\theta$。在直角坐标系中，

$$
\mathbf{A} \times \mathbf{B} = \begin{vmatrix} \hat{\mathbf{x}} & \hat{\mathbf{y}} & \hat{\mathbf{z}} \\ A_x & A_y & A_z \\ B_x & B_y & B_z \end{vmatrix} = (A_y B_z - A_z B_y)\hat{\mathbf{x}} + (A_z B_x - A_x B_z)\hat{\mathbf{y}} + (A_x B_y - A_y B_x)\hat{\mathbf{z}}
$$

叉积不满足交换律：$\mathbf{A} \times \mathbf{B} = - \mathbf{B} \times \mathbf{A}$。常用于描述力矩、洛伦兹力等。

**混合积**：三个矢量的标量三重积 $\mathbf{A} \cdot (\mathbf{B} \times \mathbf{C})$ 等于它们构成的平行六面体的体积（有正负号），行列式形式为

$$
\mathbf{A} \cdot (\mathbf{B} \times \mathbf{C}) = \begin{vmatrix} A_x & A_y & A_z \\ B_x & B_y & B_z \\ C_x & C_y & C_z \end{vmatrix}
$$

循环轮换：$\mathbf{A} \cdot (\mathbf{B} \times \mathbf{C}) = \mathbf{B} \cdot (\mathbf{C} \times \mathbf{A}) = \mathbf{C} \cdot (\mathbf{A} \times \mathbf{B})$。

**矢量三重积**：$\mathbf{A} \times (\mathbf{B} \times \mathbf{C}) = \mathbf{B}(\mathbf{A} \cdot \mathbf{C}) - \mathbf{C}(\mathbf{A} \cdot \mathbf{B})$，这个恒等式在后续推导旋度的旋度时至关重要。

---

### 二、标量场的梯度

#### 2.1 方向导数与梯度定义

标量场 $f(x,y,z)$ 表示空间中每点的一个数值。我们常需要知道从一点出发沿某个方向 $f$ 的变化快慢。设沿单位矢量 $\hat{\mathbf{u}}$ 的方向移动微小距离 $ds$，函数的变化率为方向导数：

$$
\frac{df}{ds} = \lim_{\Delta s \to 0} \frac{f(\mathbf{r} + \hat{\mathbf{u}} \Delta s) - f(\mathbf{r})}{\Delta s}
$$

由全微分公式，

$$
df = \frac{\partial f}{\partial x} dx + \frac{\partial f}{\partial y} dy + \frac{\partial f}{\partial z} dz = \left( \frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}, \frac{\partial f}{\partial z} \right) \cdot (dx, dy, dz)
$$

若位移沿 $\hat{\mathbf{u}}$ 方向，$(dx,dy,dz) = \hat{\mathbf{u}} ds$，则

$$
\frac{df}{ds} = \left( \frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}, \frac{\partial f}{\partial z} \right) \cdot \hat{\mathbf{u}}
$$

定义矢量

$$
\boxed{\nabla f \equiv \frac{\partial f}{\partial x} \hat{\mathbf{x}} + \frac{\partial f}{\partial y} \hat{\mathbf{y}} + \frac{\partial f}{\partial z} \hat{\mathbf{z}}}
$$

称为标量场 $f$ 的**梯度**。于是方向导数可写为

$$
\frac{df}{ds} = \nabla f \cdot \hat{\mathbf{u}} = |\nabla f| \cos\phi
$$

其中 $\phi$ 是 $\nabla f$ 与 $\hat{\mathbf{u}}$ 的夹角。由此得出梯度的核心性质：
- 梯度指向函数 $f$ **增加最快的方向**，其大小等于该方向上的方向导数最大值。
- 梯度垂直于等值面（$f = \text{常数}$ 的曲面），因为沿等值面方向的方向导数为零。

#### 2.2 梯度的运算法则

梯度算符 $\nabla$（读作“del”或“nabla”）作用于标量场，它本身是一个矢量性微分算符：

$$
\nabla = \hat{\mathbf{x}} \frac{\partial}{\partial x} + \hat{\mathbf{y}} \frac{\partial}{\partial y} + \hat{\mathbf{z}} \frac{\partial}{\partial z}
$$

满足分配律：$\nabla (f + g) = \nabla f + \nabla g$，乘积法则：$\nabla (fg) = f \nabla g + g \nabla f$。

---

### 三、矢量场的散度与高斯散度定理

#### 3.1 矢量场的通量

设矢量场 $\mathbf{F}(x,y,z)$。考虑空间一小块面积元 $d\mathbf{A}$，其方向为法向（通常对于闭合曲面取外法向）。$\mathbf{F}$ 通过该面元的**通量**为 $d\Phi = \mathbf{F} \cdot d\mathbf{A}$。通过整个曲面 $S$ 的通量为

$$
\Phi = \iint_S \mathbf{F} \cdot d\mathbf{A}
$$

#### 3.2 散度的定义

散度描述矢量场在某点“源”或“汇”的强度。定义 $\mathbf{F}$ 在点 $\mathbf{r}$ 的散度为包围该点的小闭合曲面上的通量与所围体积之比，当体积趋于零时的极限：

$$
\boxed{\nabla \cdot \mathbf{F} \equiv \lim_{\Delta V \to 0} \frac{1}{\Delta V} \oiint_{\partial(\Delta V)} \mathbf{F} \cdot d\mathbf{A}}
$$

在直角坐标系中，取以 $(x,y,z)$ 为中心、边长为 $\Delta x, \Delta y, \Delta z$ 的小长方体。分析通过各面的净通量：通过左右两面（与 $x$ 轴垂直）的净通量约为 $\frac{\partial F_x}{\partial x} \Delta x \Delta y \Delta z$，类似可得 $y,z$ 分量。总和除以体积 $\Delta V = \Delta x \Delta y \Delta z$，并取极限得

$$
\boxed{\nabla \cdot \mathbf{F} = \frac{\partial F_x}{\partial x} + \frac{\partial F_y}{\partial y} + \frac{\partial F_z}{\partial z}}
$$

这就是散度在直角坐标下的计算公式。若 $\nabla \cdot \mathbf{F} > 0$，该点有“源”（流体流出）；$\nabla \cdot \mathbf{F} < 0$ 则有“汇”；$\nabla \cdot \mathbf{F} = 0$ 称为无散场。

#### 3.3 高斯散度定理

将任意体积 $V$ 分割成许多小体积元。对每个小体积元，由散度定义近似有 $\mathbf{F}$ 穿过表面的通量 $\approx (\nabla \cdot \mathbf{F}) \Delta V$。将所有体元的通量相加，内部相邻面的通量互相抵消，只剩下外表面 $S = \partial V$ 的通量。取极限即得**高斯散度定理**：

$$
\boxed{\iint_{\partial V} \mathbf{F} \cdot d\mathbf{A} = \iiint_V (\nabla \cdot \mathbf{F}) \, dV}
$$

它将矢量场在闭合曲面上的面积分转化为该面所围体积内的体积分，是电磁学中将积分形式的麦克斯韦方程转化为微分形式的关键工具。

---

### 四、矢量场的旋度与斯托克斯定理

#### 4.1 矢量场的环量

矢量场 $\mathbf{F}$ 沿闭合曲线 $C$ 的线积分 $\oint_C \mathbf{F} \cdot d\mathbf{l}$ 称为**环量**，衡量场的旋转趋势。若环量为零，场为保守场（如静电场）。

#### 4.2 旋度的定义

旋度矢量在某个方向 $\hat{\mathbf{n}}$ 上的投影，定义为围绕该方向的无穷小回路上的环量与回路所围面积之比，当面积趋于零时的极限：

$$
(\nabla \times \mathbf{F}) \cdot \hat{\mathbf{n}} \equiv \lim_{\Delta S \to 0} \frac{1}{\Delta S} \oint_{\partial(\Delta S)} \mathbf{F} \cdot d\mathbf{l}
$$

在直角坐标系中，分别取垂直于 $x,y,z$ 轴的小矩形回路计算环量，可得旋度的各分量。例如，在 $xy$ 平面内取小矩形，环量 $\oint \mathbf{F} \cdot d\mathbf{l} \approx \left( \frac{\partial F_y}{\partial x} - \frac{\partial F_x}{\partial y} \right) \Delta x \Delta y$，该环量对应于旋度的 $z$ 分量。由此得到

$$
\boxed{\nabla \times \mathbf{F} = \begin{vmatrix} \hat{\mathbf{x}} & \hat{\mathbf{y}} & \hat{\mathbf{z}} \\ \frac{\partial}{\partial x} & \frac{\partial}{\partial y} & \frac{\partial}{\partial z} \\ F_x & F_y & F_z \end{vmatrix}}
$$

展开为

$$
\nabla \times \mathbf{F} = \left( \frac{\partial F_z}{\partial y} - \frac{\partial F_y}{\partial z} \right) \hat{\mathbf{x}} + \left( \frac{\partial F_x}{\partial z} - \frac{\partial F_z}{\partial x} \right) \hat{\mathbf{y}} + \left( \frac{\partial F_y}{\partial x} - \frac{\partial F_x}{\partial y} \right) \hat{\mathbf{z}}
$$

旋度不为零的场称为有旋场。

#### 4.3 斯托克斯定理

将任意曲面 $S$ 划分为许多小面元，对每个面元，旋度的法向分量乘以面元面积近似等于绕该面元边界的环量。将所有面元的环量相加，内部相邻边界的线积分相互抵消，只剩下曲面 $S$ 的总边界 $C = \partial S$ 的环量。取极限即得**斯托克斯定理**：

$$
\boxed{\iint_S (\nabla \times \mathbf{F}) \cdot d\mathbf{A} = \oint_{\partial S} \mathbf{F} \cdot d\mathbf{l}}
$$

它将矢量场在闭合曲线上的环量转化为以该曲线为边界的曲面上的旋度面积分，在电磁感应定律和安培定律的分析中不可或缺。

---

### 五、拉普拉斯算子

#### 5.1 标量场的拉普拉斯算子

梯度 $\nabla f$ 是一个矢量场。再取其散度，得到标量场 $f$ 的**拉普拉斯算子**：

$$
\boxed{\nabla^2 f \equiv \nabla \cdot (\nabla f) = \frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2} + \frac{\partial^2 f}{\partial z^2}}
$$

在静电学中，电势满足泊松方程 $\nabla^2 V = -\rho / \varepsilon_0$；无电荷区域满足拉普拉斯方程 $\nabla^2 V = 0$。

#### 5.2 矢量场的拉普拉斯算子

矢量场 $\mathbf{F}$ 的拉普拉斯算子定义为

$$
\nabla^2 \mathbf{F} \equiv \nabla (\nabla \cdot \mathbf{F}) - \nabla \times (\nabla \times \mathbf{F})
$$

在直角坐标系中，这等价于对各分量分别作用标量拉普拉斯算子：

$$
\nabla^2 \mathbf{F} = (\nabla^2 F_x) \hat{\mathbf{x}} + (\nabla^2 F_y) \hat{\mathbf{y}} + (\nabla^2 F_z) \hat{\mathbf{z}}
$$

该恒等式在推导电磁波方程时用到，利用 $\nabla \times (\nabla \times \mathbf{E}) = \nabla (\nabla \cdot \mathbf{E}) - \nabla^2 \mathbf{E}$。

#### 5.3 常用恒等式

以下两个恒等式在电磁学中极为重要：

- 梯度的旋度恒为零：$\nabla \times (\nabla f) = 0$（标量场的梯度场是无旋场）。
- 旋度的散度恒为零：$\nabla \cdot (\nabla \times \mathbf{F}) = 0$（矢量场的旋度场是无散场）。

它们保证了静电势 $\mathbf{E} = -\nabla V$ 自动满足 $\nabla \times \mathbf{E} = 0$，以及磁矢势 $\mathbf{B} = \nabla \times \mathbf{A}$ 自动满足 $\nabla \cdot \mathbf{B} = 0$。

---

### 六、亥姆霍兹定理

**亥姆霍兹定理**是矢量分析的一个基本结论：在无限大空间中（或给定适当边界条件），任何一个足够光滑的矢量场 $\mathbf{F}$ 都可以唯一地分解为一个无旋场 $\mathbf{F}_1$（满足 $\nabla \times \mathbf{F}_1 = 0$）和一个无散场 $\mathbf{F}_2$（满足 $\nabla \cdot \mathbf{F}_2 = 0$）的叠加：

$$
\boxed{\mathbf{F} = -\nabla \phi + \nabla \times \mathbf{A}}
$$

其中标量势 $\phi$ 和矢量势 $\mathbf{A}$ 可以通过 $\mathbf{F}$ 的散度和旋度求得。具体地，令

$$
D = \nabla \cdot \mathbf{F}, \quad \mathbf{C} = \nabla \times \mathbf{F}
$$

则 $\phi$ 和 $\mathbf{A}$ 分别满足泊松方程：

$$
\nabla^2 \phi = -D, \quad \nabla^2 \mathbf{A} = -\mathbf{C} \quad (\text{在库仑规范下})
$$

这就意味着：**一个矢量场由其散度和旋度共同唯一确定**（外加边界条件）。麦克斯韦方程组的核心正是在给定电荷分布（散度源）和电流分布（旋度源）下求解电场和磁场，正是亥姆霍兹定理的直接体现。

---

### 七、曲线坐标系中的梯度、散度、旋度与拉普拉斯算子

电磁学中球对称、柱对称的问题使用相应的曲线坐标更为方便。以下直接给出柱坐标 $(r, \varphi, z)$ 和球坐标 $(r, \theta, \varphi)$ 中的表达式，并简要说明推导思路。

#### 7.1 柱坐标系

坐标变换：$x = r \cos\varphi, \; y = r \sin\varphi, \; z = z$，基矢 $\hat{\mathbf{r}}, \hat{\boldsymbol{\varphi}}, \hat{\mathbf{z}}$。

**梯度**：

$$
\nabla f = \frac{\partial f}{\partial r} \hat{\mathbf{r}} + \frac{1}{r} \frac{\partial f}{\partial \varphi} \hat{\boldsymbol{\varphi}} + \frac{\partial f}{\partial z} \hat{\mathbf{z}}
$$

**散度**：

$$
\nabla \cdot \mathbf{F} = \frac{1}{r} \frac{\partial}{\partial r} (r F_r) + \frac{1}{r} \frac{\partial F_\varphi}{\partial \varphi} + \frac{\partial F_z}{\partial z}
$$

**旋度**：

$$
\nabla \times \mathbf{F} = \left( \frac{1}{r} \frac{\partial F_z}{\partial \varphi} - \frac{\partial F_\varphi}{\partial z} \right) \hat{\mathbf{r}} + \left( \frac{\partial F_r}{\partial z} - \frac{\partial F_z}{\partial r} \right) \hat{\boldsymbol{\varphi}} + \frac{1}{r} \left( \frac{\partial}{\partial r} (r F_\varphi) - \frac{\partial F_r}{\partial \varphi} \right) \hat{\mathbf{z}}
$$

**标量拉普拉斯算子**：

$$
\nabla^2 f = \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial f}{\partial r} \right) + \frac{1}{r^2} \frac{\partial^2 f}{\partial \varphi^2} + \frac{\partial^2 f}{\partial z^2}
$$

#### 7.2 球坐标系

坐标变换：$x = r \sin\theta \cos\varphi, \; y = r \sin\theta \sin\varphi, \; z = r \cos\theta$，基矢 $\hat{\mathbf{r}}, \hat{\boldsymbol{\theta}}, \hat{\boldsymbol{\varphi}}$。

**梯度**：

$$
\nabla f = \frac{\partial f}{\partial r} \hat{\mathbf{r}} + \frac{1}{r} \frac{\partial f}{\partial \theta} \hat{\boldsymbol{\theta}} + \frac{1}{r \sin\theta} \frac{\partial f}{\partial \varphi} \hat{\boldsymbol{\varphi}}
$$

**散度**：

$$
\nabla \cdot \mathbf{F} = \frac{1}{r^2} \frac{\partial}{\partial r} (r^2 F_r) + \frac{1}{r \sin\theta} \frac{\partial}{\partial \theta} (\sin\theta \, F_\theta) + \frac{1}{r \sin\theta} \frac{\partial F_\varphi}{\partial \varphi}
$$

**旋度**：

$$
\nabla \times \mathbf{F} = \frac{1}{r \sin\theta} \left( \frac{\partial}{\partial \theta} (\sin\theta \, F_\varphi) - \frac{\partial F_\theta}{\partial \varphi} \right) \hat{\mathbf{r}} + \frac{1}{r} \left( \frac{1}{\sin\theta} \frac{\partial F_r}{\partial \varphi} - \frac{\partial}{\partial r} (r F_\varphi) \right) \hat{\boldsymbol{\theta}} + \frac{1}{r} \left( \frac{\partial}{\partial r} (r F_\theta) - \frac{\partial F_r}{\partial \theta} \right) \hat{\boldsymbol{\varphi}}
$$

**标量拉普拉斯算子**：

$$
\nabla^2 f = \frac{1}{r^2} \frac{\partial}{\partial r} \left( r^2 \frac{\partial f}{\partial r} \right) + \frac{1}{r^2 \sin\theta} \frac{\partial}{\partial \theta} \left( \sin\theta \, \frac{\partial f}{\partial \theta} \right) + \frac{1}{r^2 \sin^2\theta} \frac{\partial^2 f}{\partial \varphi^2}
$$

这些公式的推导基于坐标变换和度量系数（拉梅系数）。以散度为例，在一般正交曲线坐标系中，$\nabla \cdot \mathbf{F} = \frac{1}{h_1 h_2 h_3} \left[ \frac{\partial}{\partial u_1}(h_2 h_3 F_1) + \frac{\partial}{\partial u_2}(h_3 h_1 F_2) + \frac{\partial}{\partial u_3}(h_1 h_2 F_3) \right]$。柱坐标下 $h_1 = 1, h_2 = r, h_3 = 1$；球坐标下 $h_1 = 1, h_2 = r, h_3 = r\sin\theta$，代入即得上述公式。

---

至此，我们建立了矢量场分析的基本框架：梯度描述标量场的变化率和方向，散度描述矢量场的源汇，旋度描述矢量场的旋转趋势。高斯散度定理和斯托克斯定理将体积分与面积分、面积分与线积分联系起来，是电磁学中由积分方程过渡到微分方程、以及求解具体问题的根本方法。亥姆霍兹定理指明了矢量场由其散度和旋度完全决定，这正是麦克斯韦方程组以 $\nabla \cdot \mathbf{D} = \rho_f$、$\nabla \cdot \mathbf{B} = 0$、$\nabla \times \mathbf{E} = -\partial \mathbf{B}/\partial t$、$\nabla \times \mathbf{H} = \mathbf{J}_f + \partial \mathbf{D}/\partial t$ 的形式完备描述电磁场的数学基础。曲线坐标中的表达式则为求解具有对称性的电磁学问题提供了直接可用的公式。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [EM-P.2 狄拉克 δ 函数与格林函数初步](https://zhuanlan.zhihu.com/p/2062629264005009857)

**课程目录：** [电磁学目录](https://zhuanlan.zhihu.com/p/2062618184079954663)
<!-- zhihu-nav-bottom:end -->
