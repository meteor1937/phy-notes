---
created: 2026-07-17T00:58:38+08:00
modified: 2026-07-17T01:23:26+08:00
zhihu-title: P.3 正交曲线坐标系
zhihu-topics: 电动力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061253464097871213
zhihu-toc: true
zhihu-created-at: 2026-07-17 01:23
review-status: reviewing
physics-checked: true
derivation-checked: true
readability-checked: true
teaching-checked: true
zhihu-layout-checked: false
last-reviewed:
zhihu-updated-at: 2026-07-21 22:07
---
<!-- zhihu-nav-top:start -->
**上一节：** [P.2 亥姆霍兹定理](https://zhuanlan.zhihu.com/p/2061253406543570185)
<!-- zhihu-nav-top:end -->

# P.3 正交曲线坐标系

电磁学边值问题的几何常常是球面、圆柱面或平面。合适的坐标系能让边界成为坐标面，但也会引入随位置变化的基矢和尺度因子。本文从坐标映射的局部几何出发推导正交坐标中的微分算符，再分别检验柱坐标与球坐标公式，并说明坐标奇点和分离变量法的边界。

## 学习目标

完成本篇后，读者应当能够：

1. 从位置映射 $\mathbf r(u_1,u_2,u_3)$ 计算尺度因子、线元、面积元和体积元；
2. 从全微分、微元通量和微元环量推导梯度、散度与旋度；
3. 独立恢复柱坐标和球坐标中的梯度、散度、旋度及标量 Laplacian；
4. 区分普通贝塞尔方程与修正贝塞尔方程对应的分离常数；
5. 判断坐标轴、原点和极轴处的公式奇性是坐标奇性还是物理奇性。

## 先修知识、约定与路线

先修知识包括多元微积分、Jacobian、线积分和面积分、Gauss 定理、Stokes 定理，以及常系数二阶常微分方程。本文采用右手定向的三维欧氏空间；$F_i=\mathbf F\cdot\hat{\mathbf e}_i$ 表示沿单位基矢的物理分量，不是协变或逆变坐标分量。角变量以弧度计，除特别说明外场在所讨论的坐标图内足够光滑。

路线是先建立局部正交坐标的几何量，再由定义推导算符；随后代入柱坐标和球坐标的尺度因子，并用分离变量、量纲、特殊场和积分定理交叉检查。


## 一、为什么需要曲线坐标系

在直角坐标中，三个坐标轴互相垂直，长度单位处处相同——沿着 $x$ 方向走一步，无论你在空间中的哪个位置，这一步的长度就是 $\Delta x$。但在球坐标中，情况不同了：坐标是 $(r, \theta, \phi)$，只有 $r$ 有长度量纲，$\theta$ 和 $\phi$ 是角度。沿 $\theta$ 方向走一步，实际走过的弧长不是 $d\theta$，而是 $r d\theta$ ——取决于你离原点有多远。沿 $\phi$ 方向走一步，弧长是 $r\sin\theta d\phi$ ——既取决于 $r$，又取决于 $\theta$。

更复杂的是，球坐标和柱坐标的基矢量方向随位置变化。在直角坐标中，$\hat{\mathbf{x}}$ 在空间中处处指向同一个方向。但在球坐标中，$\hat{\mathbf{r}}$ 的方向取决于你在哪个点——在北极和赤道，径向单位矢量指向完全不同的空间方向。

我们需要一个系统的方法来处理这些"不走直线、单位长度会变"的坐标系。


## 二、一般正交曲线坐标系

考虑位置映射 $\mathbf r=\mathbf r(u_1,u_2,u_3)$。在 Jacobian 非零的局部区域内，它与直角坐标一一对应；角变量的周期识别和坐标轴、极点等退化位置必须另行处理。若三个坐标切向量 $\partial\mathbf r/\partial u_i$ 两两正交，就称这组局部坐标为**正交曲线坐标系**。

沿坐标 $u_i$ 方向走过微小的坐标变化 $du_i$，空间中的实际位移长度不一定等于 $du_i$，而是乘以一个**尺度因子**（scale factor）$h_i$：

$$
d\mathbf{l} = h_1 du_1 \hat{\mathbf{e}}_1 + h_2 du_2 \hat{\mathbf{e}}_2 + h_3 du_3 \hat{\mathbf{e}}_3
$$

其中 $\hat{\mathbf{e}}_i$ 是第 $i$ 个坐标方向上的单位矢量，三者互相正交：$\hat{\mathbf{e}}_i \cdot \hat{\mathbf{e}}_j = \delta_{ij}$。尺度因子 $h_i$ 是空间位置的函数，由坐标变换的 Jacobian 矩阵决定。具体来说，如果把直角坐标写为 $(u_1, u_2, u_3)$ 的函数 $\mathbf{r}(u_1, u_2, u_3)$，那么：

$$
h_i = \left| \frac{\partial\mathbf{r}}{\partial u_i} \right|
$$

取 $h_i>0$，单位基矢为

$$
\hat{\mathbf e}_i=\frac1{h_i}\frac{\partial\mathbf r}{\partial u_i}.
$$

在右手坐标次序下，Jacobian 的绝对值是 $h_1h_2h_3$。因此线元、三个有向坐标面元的大小和体积元分别为

$$
ds_i=h_i\,du_i,
$$

$$
dS_1=h_2h_3\,du_2du_3,
\quad dS_2=h_3h_1\,du_3du_1,
\quad dS_3=h_1h_2\,du_1du_2,
$$

$$
\boxed{dV=h_1h_2h_3\,du_1du_2du_3}.
$$

这些几何因子是后续微分算符的来源。若某个 $h_i=0$，坐标图在该处退化，含 $1/h_i$ 的分量公式不能逐项直接代入。


## 三、一般正交曲线坐标系中的梯度、散度、旋度

下面始终使用单位基矢分量 $F_i$。基矢随位置变化的效应不需要额外手工求导；它已经包含在沿边长度和坐标面面积的尺度因子中。

### 3.1 梯度

标量场 $\Psi$ 的梯度沿某个方向的分量，等于 $\Psi$ 沿该方向的方向导数。在曲线坐标中，沿 $\hat{\mathbf e}_1$ 方向走过 $h_1du_1$ 的实际距离，标量变化为 $d\Psi=(\partial\Psi/\partial u_1)du_1$。因此该方向的方向导数为 $(1/h_1)\partial\Psi/\partial u_1$。三个方向合起来：

$$
\boxed{\nabla\Psi = \frac{1}{h_1}\frac{\partial\Psi}{\partial u_1}\hat{\mathbf{e}}_1 + \frac{1}{h_2}\frac{\partial\Psi}{\partial u_2}\hat{\mathbf{e}}_2 + \frac{1}{h_3}\frac{\partial\Psi}{\partial u_3}\hat{\mathbf{e}}_3}
$$

### 3.2 散度

散度衡量矢量场从一点向外流出的净通量。取坐标面围成的一个微小曲边长方体——在 $u_1$ 方向的两个面分别位于 $u_1$ 和 $u_1 + du_1$，面积各为 $h_2 h_3 du_2 du_3$。穿过这两个面的通量差除以微元体积 $h_1 h_2 h_3 du_1 du_2 du_3$，取极限，三个方向相加得到：

$$
\boxed{\nabla\cdot\mathbf{F} = \frac{1}{h_1 h_2 h_3}\left[ \frac{\partial}{\partial u_1}(h_2 h_3 F_1) + \frac{\partial}{\partial u_2}(h_3 h_1 F_2) + \frac{\partial}{\partial u_3}(h_1 h_2 F_3) \right]}
$$

注意：散度的表达式里出现了尺度因子乘在矢量分量上再对坐标求导。这和直角坐标中简单的 $\partial F_x/\partial x$ 不同，因为面积元和体积元本身随位置变化。

### 3.3 旋度

旋度的每个分量是单位面积上的有向环量。取法向为 $+\hat{\mathbf e}_1$ 的 $u_2u_3$ 小回路，边界按右手定则取向。沿两条 $u_2$ 边的线积分之差给出 $-\partial_{u_3}(h_2F_2)du_2du_3$，沿两条 $u_3$ 边给出 $+\partial_{u_2}(h_3F_3)du_2du_3$。除以面积 $h_2h_3du_2du_3$，得到

$$
(\nabla\times\mathbf F)_1
=\frac1{h_2h_3}
\left[\frac{\partial(h_3F_3)}{\partial u_2}
-\frac{\partial(h_2F_2)}{\partial u_3}\right].
$$

循环置换指标即可得到其余两项。对右手坐标次序，三个分量可压缩为行列式形式：

$$
\boxed{\nabla\times\mathbf{F} = \frac{1}{h_1 h_2 h_3}\begin{vmatrix} h_1\hat{\mathbf{e}}_1 & h_2\hat{\mathbf{e}}_2 & h_3\hat{\mathbf{e}}_3 \\ \frac{\partial}{\partial u_1} & \frac{\partial}{\partial u_2} & \frac{\partial}{\partial u_3} \\ h_1 F_1 & h_2 F_2 & h_3 F_3 \end{vmatrix}}
$$

这里的行列式只是三个分量公式的记忆方式，微分算符只作用于第三行及其右侧场分量。若坐标次序是左手定向，必须额外改变取向符号，不能直接照搬此式。

### 3.4 拉普拉斯算符

标量 Laplacian 是梯度的散度：$\nabla^2\Psi=\nabla\cdot(\nabla\Psi)$。把梯度分量 $F_i=(1/h_i)\partial\Psi/\partial u_i$ 代入散度公式：

$$
\boxed{\nabla^2\Psi = \frac{1}{h_1 h_2 h_3}\left[ \frac{\partial}{\partial u_1}\!\left( \frac{h_2 h_3}{h_1}\frac{\partial\Psi}{\partial u_1} \right) + \frac{\partial}{\partial u_2}\!\left( \frac{h_3 h_1}{h_2}\frac{\partial\Psi}{\partial u_2} \right) + \frac{\partial}{\partial u_3}\!\left( \frac{h_1 h_2}{h_3}\frac{\partial\Psi}{\partial u_3} \right) \right]}
$$

在具体坐标系中代入尺度因子后，表达式通常会显著简化。


## 四、柱坐标

### 4.1 坐标变换与尺度因子

柱坐标 $(\rho, \phi, z)$ 与直角坐标的关系为：

$$
x = \rho\cos\phi, \quad y = \rho\sin\phi, \quad z = z
$$

位置矢量 $\mathbf{r} = \rho\cos\phi\,\hat{\mathbf{x}} + \rho\sin\phi\,\hat{\mathbf{y}} + z\,\hat{\mathbf{z}}$。计算尺度因子：

- 对 $\rho$ 求偏导：$\frac{\partial\mathbf{r}}{\partial\rho} = \cos\phi\,\hat{\mathbf{x}} + \sin\phi\,\hat{\mathbf{y}}$，长度 $h_\rho = \sqrt{\cos^2\phi + \sin^2\phi} = 1$
- 对 $\phi$ 求偏导：$\frac{\partial\mathbf{r}}{\partial\phi} = -\rho\sin\phi\,\hat{\mathbf{x}} + \rho\cos\phi\,\hat{\mathbf{y}}$，长度 $h_\phi = \rho$
- 对 $z$ 求偏导：$\frac{\partial\mathbf{r}}{\partial z} = \hat{\mathbf{z}}$，长度 $h_z = 1$

$$
h_\rho = 1, \quad h_\phi = \rho, \quad h_z = 1
$$

单位矢量 $\hat{\boldsymbol{\rho}} = \cos\phi\,\hat{\mathbf{x}} + \sin\phi\,\hat{\mathbf{y}}$（径向向外），$\hat{\boldsymbol{\phi}} = -\sin\phi\,\hat{\mathbf{x}} + \cos\phi\,\hat{\mathbf{y}}$（切向，逆时针），$\hat{\mathbf{z}}$ 不变。

### 4.2 柱坐标中的梯度、散度、旋度

代入 $h_\rho=1, h_\phi=\rho, h_z=1$：

**梯度**：

$$
\nabla\Psi = \frac{\partial\Psi}{\partial\rho}\hat{\boldsymbol{\rho}} + \frac{1}{\rho}\frac{\partial\Psi}{\partial\phi}\hat{\boldsymbol{\phi}} + \frac{\partial\Psi}{\partial z}\hat{\mathbf{z}}
$$

**散度**：

$$
\nabla\cdot\mathbf{F} = \frac{1}{\rho}\frac{\partial}{\partial\rho}(\rho F_\rho) + \frac{1}{\rho}\frac{\partial F_\phi}{\partial\phi} + \frac{\partial F_z}{\partial z}
$$

**旋度**：

$$
\nabla\times\mathbf{F} = \frac{1}{\rho}\begin{vmatrix} \hat{\boldsymbol{\rho}} & \rho\hat{\boldsymbol{\phi}} & \hat{\mathbf{z}} \\ \frac{\partial}{\partial\rho} & \frac{\partial}{\partial\phi} & \frac{\partial}{\partial z} \\ F_\rho & \rho F_\phi & F_z \end{vmatrix}
$$

**标量 Laplacian**（把 $\mathbf F=\nabla\Psi$ 代入散度公式）：

$$
\boxed{\nabla^2\Psi = \frac{1}{\rho}\frac{\partial}{\partial\rho}\left( \rho\frac{\partial\Psi}{\partial\rho} \right) + \frac{1}{\rho^2}\frac{\partial^2\Psi}{\partial\phi^2} + \frac{\partial^2\Psi}{\partial z^2}}
$$

柱坐标 Laplacian 的第一项是径向部分。当 $\Psi$ 不依赖方位角 $\phi$ 和 $z$ 时，$\nabla^2\Psi=(1/\rho)d(\rho\,d\Psi/d\rho)/d\rho$。在无源区域令它为零，解为 $\ln\rho$ 和常数；$\ln\rho$ 是无限长线电荷电势的径向形式，但轴线上的源必须另行处理。

### 4.3 柱坐标分离变量与贝塞尔方程

考虑柱坐标中的 Laplace 方程 $\nabla^2\Psi=0$。设 $\Psi(\rho,\phi,z)=R(\rho)\Phi(\phi)Z(z)$。代入后除以 $\Psi$，再乘以 $\rho^2$：

$$
\frac{\rho}{R}\frac{d}{d\rho}\left( \rho\frac{dR}{d\rho} \right) + \frac{1}{\Phi}\frac{d^2\Phi}{d\phi^2} + \rho^2\frac{1}{Z}\frac{d^2Z}{dz^2} = 0
$$

单值性要求 $\Phi(\phi+2\pi)=\Phi(\phi)$，故可取

$$
\Phi''+m^2\Phi=0,
\qquad m\in\mathbb Z.
$$

若轴向方程取 $Z''-k^2Z=0$，则 $Z$ 沿 $z$ 指数变化，径向方程为

$$
\rho\frac{d}{d\rho}\left( \rho\frac{dR}{d\rho} \right) + (k^2\rho^2 - m^2)R = 0
$$

这是 $m$ 阶普通贝塞尔方程，解为 $J_m(k\rho)$ 和 $Y_m(k\rho)$。若轴向方程改取 $Z''+k^2Z=0$，径向方程则为

$$
\rho\frac{d}{d\rho}\left(\rho\frac{dR}{d\rho}\right)
-(k^2\rho^2+m^2)R=0,
$$

其解是修正贝塞尔函数 $I_m(k\rho)$ 和 $K_m(k\rho)$。究竟出现哪一支由原偏微分方程、分离常数符号和边界条件共同决定。若区域包含轴线 $\rho=0$，还要排除在那里发散的 $Y_m$ 或 $K_m$；若区域不含轴线，例如同轴电缆的环形区域，则这些解可以保留。


## 五、球坐标

### 5.1 坐标变换与尺度因子

球坐标 $(r, \theta, \phi)$ 与直角坐标的关系为：

$$
x = r\sin\theta\cos\phi, \quad y = r\sin\theta\sin\phi, \quad z = r\cos\theta
$$

位置矢量 $\mathbf{r} = r\sin\theta\cos\phi\,\hat{\mathbf{x}} + r\sin\theta\sin\phi\,\hat{\mathbf{y}} + r\cos\theta\,\hat{\mathbf{z}}$。计算尺度因子：

- 对 $r$ 求偏导：$\frac{\partial\mathbf{r}}{\partial r} = \sin\theta\cos\phi\,\hat{\mathbf{x}} + \sin\theta\sin\phi\,\hat{\mathbf{y}} + \cos\theta\,\hat{\mathbf{z}}$，长度 $h_r = \sqrt{\sin^2\theta(\cos^2\phi+\sin^2\phi) + \cos^2\theta} = 1$
- 对 $\theta$ 求偏导：$\frac{\partial\mathbf{r}}{\partial\theta} = r\cos\theta\cos\phi\,\hat{\mathbf{x}} + r\cos\theta\sin\phi\,\hat{\mathbf{y}} - r\sin\theta\,\hat{\mathbf{z}}$，长度 $h_\theta = r\sqrt{\cos^2\theta + \sin^2\theta} = r$
- 对 $\phi$ 求偏导：$\frac{\partial\mathbf{r}}{\partial\phi} = -r\sin\theta\sin\phi\,\hat{\mathbf{x}} + r\sin\theta\cos\phi\,\hat{\mathbf{y}}$，长度 $h_\phi = r\sin\theta$

$$
h_r = 1, \quad h_\theta = r, \quad h_\phi = r\sin\theta
$$

### 5.2 球坐标中的梯度、散度、旋度

代入 $h_r=1, h_\theta=r, h_\phi=r\sin\theta$：

**梯度**：

$$
\nabla\Psi = \frac{\partial\Psi}{\partial r}\hat{\mathbf{r}} + \frac{1}{r}\frac{\partial\Psi}{\partial\theta}\hat{\boldsymbol{\theta}} + \frac{1}{r\sin\theta}\frac{\partial\Psi}{\partial\phi}\hat{\boldsymbol{\phi}}
$$

**散度**：

$$
\nabla\cdot\mathbf{F} = \frac{1}{r^2}\frac{\partial}{\partial r}(r^2 F_r) + \frac{1}{r\sin\theta}\frac{\partial}{\partial\theta}(\sin\theta F_\theta) + \frac{1}{r\sin\theta}\frac{\partial F_\phi}{\partial\phi}
$$

**旋度**：

$$
\nabla\times\mathbf{F} = \frac{1}{r^2\sin\theta}\begin{vmatrix} \hat{\mathbf{r}} & r\hat{\boldsymbol{\theta}} & r\sin\theta\hat{\boldsymbol{\phi}} \\ \frac{\partial}{\partial r} & \frac{\partial}{\partial\theta} & \frac{\partial}{\partial\phi} \\ F_r & r F_\theta & r\sin\theta F_\phi \end{vmatrix}
$$

**标量 Laplacian**（把 $\mathbf F=\nabla\Psi$ 代入散度公式）：

$$
\boxed{\nabla^2\Psi = \frac{1}{r^2}\frac{\partial}{\partial r}\left( r^2\frac{\partial\Psi}{\partial r} \right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left( \sin\theta\frac{\partial\Psi}{\partial\theta} \right) + \frac{1}{r^2\sin^2\theta}\frac{\partial^2\Psi}{\partial\phi^2}}
$$

当 $\Psi$ 只依赖 $r$ 时，$\nabla^2\Psi=(1/r^2)d(r^2d\Psi/dr)/dr$。令 $\Psi=u/r$，无源方程化为 $d^2u/dr^2=0$，所以 $\Psi=A+B/r$。$1/r$ 只在 $r>0$ 满足普通 Laplace 方程；原点对应点源奇性。

### 5.3 球坐标分离变量与勒让德方程

考虑球坐标中的 Laplace 方程 $\nabla^2\Psi=0$。设 $\Psi(r,\theta,\phi)=R(r)\Theta(\theta)\Phi(\phi)$。代入后除以 $\Psi$，再乘以 $r^2$：

$$
\frac{1}{R}\frac{d}{dr}\left( r^2\frac{dR}{dr} \right) + \frac{1}{\Theta\sin\theta}\frac{d}{d\theta}\left( \sin\theta\frac{d\Theta}{d\theta} \right) + \frac{1}{\Phi\sin^2\theta}\frac{d^2\Phi}{d\phi^2} = 0
$$

径向部分只依赖于 $r$，角向部分只依赖于 $(\theta,\phi)$。两者之和为零，所以各自必须等于一个常数，且互为相反数。设径向部分等于 $l(l+1)$（这个奇怪的形式是历史选择，但它会给出漂亮的整数解）：

$$
\frac{d}{dr}\left( r^2\frac{dR}{dr} \right) - l(l+1)R = 0
$$

这个方程的解是 $R(r) = Ar^l + Br^{-(l+1)}$。两项分别对应内部解（在原点有限，取 $r^l$）和外部解（在无穷远衰减，取 $r^{-(l+1)}$）。

角向部分为：

$$
\frac{1}{\Theta\sin\theta}\frac{d}{d\theta}\left( \sin\theta\frac{d\Theta}{d\theta} \right) + \frac{1}{\Phi\sin^2\theta}\frac{d^2\Phi}{d\phi^2} = -l(l+1)
$$

再乘以 $\sin^2\theta$，分离 $\Phi(\phi)$ 部分。设 $\Phi(\phi) = e^{im\phi}$（$m$ 为整数以保证 $2\pi$ 周期单值性），代入得到 $\Theta(\theta)$ 满足的方程：

$$
\frac{1}{\sin\theta}\frac{d}{d\theta}\left( \sin\theta\frac{d\Theta}{d\theta} \right) + \left[ l(l+1) - \frac{m^2}{\sin^2\theta} \right]\Theta = 0
$$

令 $x = \cos\theta$，$\Theta(\theta) = P(x)$。$\frac{d}{d\theta} = -\sin\theta\frac{d}{dx}$，$\sin\theta\frac{d}{d\theta} = -(1-x^2)\frac{d}{dx}$。代入整理得到：

$$
\boxed{\frac{d}{dx}\left[ (1-x^2)\frac{dP}{dx} \right] + \left[ l(l+1) - \frac{m^2}{1-x^2} \right]P = 0}
$$

这就是**连带勒让德方程**。当 $m=0$（轴对称情况，$\phi$ 独立），它退化为**勒让德方程**：

$$
\frac{d}{dx}\left[ (1-x^2)\frac{dP}{dx} \right] + l(l+1)P = 0
$$

对轴对称的 $m=0$ 情形，两个线性独立解可取第一类和第二类勒让德函数。要求解在 $x=\pm1$ 都正则，会留下勒让德多项式 $P_l(x)$，并要求 $l=0,1,2,\ldots$；$Q_l(x)$ 在端点发散。一般 $m$ 情形还要求 $|m|\le l$，正则角向解为连带勒让德函数 $P_l^m(x)$。这些整数条件来自球面上的周期性和极点正则性。量子力学角动量使用同一个球面本征值问题，但其物理量子化还依赖量子态与算符的理论结构，不能仅由这里的经典边界条件推出。

勒让德多项式的前几项为：$P_0(x)=1$，$P_1(x)=x$，$P_2(x)=\frac{1}{2}(3x^2-1)$，$P_3(x)=\frac{1}{2}(5x^3-3x)$。它们在 $[-1,1]$ 上正交完备，适合展开只依赖 $\theta$ 的轴对称球面函数；有 $\phi$ 依赖的一般球面函数需要完整的球谐函数基。

连带勒让德函数 $P_l^m(x)$ 与方位因子共同组成球谐函数 $Y_{lm}(\theta,\phi)=N_{lm}P_l^m(\cos\theta)e^{im\phi}$。在给定归一化与相位约定后，$Y_{lm}$ 构成球面平方可积函数的正交完备基，并用于静电多极展开、辐射角分布和量子角动量问题。

### 5.4 坐标退化点不是自动的物理奇点

柱坐标在 $\rho=0$ 处不能定义唯一的 $\phi$ 与 $\hat{\boldsymbol\phi}$；球坐标在 $r=0$ 处全部角坐标退化，在 $\theta=0,\pi$ 处 $\phi$ 和方位基矢退化。因此算符中的 $1/\rho$、$1/r$ 或 $1/\sin\theta$ 不能单独用来判断场发散。正确做法是把场写回直角坐标，或从非退化区域取极限并检查结果是否唯一。

例如常矢量 $\mathbf F=F_0\hat{\mathbf x}$ 在柱坐标中写成

$$
\mathbf F=F_0\cos\phi\,\hat{\boldsymbol\rho}
-F_0\sin\phi\,\hat{\boldsymbol\phi}.
$$

两个分量在轴线上依赖未定义的 $\phi$，但矢量本身处处光滑。这是基矢退化，不是物理奇点。


## 六、三种坐标系的应用场景对比

| 坐标系 | 尺度因子 | 对称性 | 基本方程的解 | 典型应用 |
|--------|---------|--------|------------|---------|
| 直角 | $h_x=h_y=h_z=1$ | 平移 | 三角函数/指数函数 | 矩形波导、无限大平面 |
| 柱坐标 | $h_\rho=1, h_\phi=\rho, h_z=1$ | 柱对称 | 贝塞尔函数（径向）+ 三角函数 | 圆柱波导、同轴电缆、光纤 |
| 球坐标 | $h_r=1, h_\theta=r, h_\phi=r\sin\theta$ | 球面边界或旋转对称 | 球谐函数（角向）+ 幂函数（径向） | 点电荷、球形导体、多极展开 |

优先让问题的边界与坐标面重合。球面边界通常使用球坐标，一般角向解由球谐函数构成，只有轴对称时才退化为勒让德多项式。圆柱面边界通常使用柱坐标，径向解可能是普通或修正贝塞尔函数。平面和矩形边界通常使用直角坐标，分离后的解常为三角函数或指数函数。

坐标系的作用是使边界条件和对称性易于表达，并不保证方程一定可分离。是否能够分离还取决于微分方程、材料参数和边界条件的结构。

## 七、完整例题：径向平方反比场

**已知：** 在 $r>0$ 区域有 $\mathbf F=C\hat{\mathbf r}/r^2$，其中 $C$ 为常数。求普通散度、半径 $R$ 球面的通量，并判断原点是否可以直接代入球坐标散度公式。

对 $r>0$ 使用球坐标散度：

$$
\nabla\cdot\mathbf F
=\frac1{r^2}\frac{d}{dr}\left(r^2\frac C{r^2}\right)=0.
$$

半径 $R$ 的球面上 $\mathbf F\cdot\hat{\mathbf n}=C/R^2$，面积元为 $R^2\sin\theta\,d\theta d\phi$，故

$$
\oint_{S_R}\mathbf F\cdot d\mathbf a
=\frac C{R^2}\,4\pi R^2=4\pi C.
$$

两式不矛盾。第一式只在 $r>0$ 的非退化光滑区域成立；原点既是球坐标退化点，也是该场真正的物理奇点。若在包含原点的全空间用分布表示，则

$$
\nabla\cdot\frac{\hat{\mathbf r}}{r^2}
=4\pi\delta^{(3)}(\mathbf r).
$$

量纲检查也一致：若 $[C]=[F]L^2$，散度单位为 $[F]/L$，而 $C\delta^{(3)}(\mathbf r)$ 同样具有该单位。这个例子说明，分量公式的适用区域必须与积分定理包围的奇点一起检查。

## 八、练习与反馈

### 练习 1：从尺度因子得到测度

写出柱坐标与球坐标的体积元。

**答案：** 由 $dV=h_1h_2h_3du_1du_2du_3$，

$$
dV_{\rm cyl}=\rho\,d\rho\,d\phi\,dz,
\qquad
dV_{\rm sph}=r^2\sin\theta\,dr\,d\theta\,d\phi.
$$

### 练习 2：柱坐标算符

设 $\mathbf F=A\rho\,\hat{\boldsymbol\rho}+Bz\,\hat{\mathbf z}$，其中 $A,B$ 为常数。求散度与旋度。

**答案：**

$$
\nabla\cdot\mathbf F
=\frac1\rho\frac{\partial(A\rho^2)}{\partial\rho}
+\frac{\partial(Bz)}{\partial z}=2A+B,
\qquad
\nabla\times\mathbf F=\mathbf0.
$$

### 练习 3：球坐标径向幂次

设 $\Psi=r^nP_l(\cos\theta)$ 且与方位角无关。在 $r>0$ 代入 Laplace 方程，求允许的 $n$。

**答案：** 利用勒让德方程，角向项给出 $-l(l+1)\Psi/r^2$，径向项给出 $n(n+1)\Psi/r^2$，故

$$
n(n+1)=l(l+1),
\qquad n=l\ \text{或}\ n=-(l+1).
$$

这正对应球内正则支 $r^l$ 与无穷远衰减支 $r^{-(l+1)}$。

### 练习 4：概念诊断

球坐标公式在 $\theta=0$ 含 $1/\sin\theta$，能否据此断言所有场在极轴发散？

**答案：** 不能。极轴处是坐标图退化。必须检查完整矢量或标量的极限；对光滑场，分子中的角向依赖会与尺度因子共同给出有限结果。只有坐标无关的量或积分检验显示发散时，才能判定存在物理奇性。

## 九、主要参考

1. George B. Arfken, Hans J. Weber, Frank E. Harris, *Mathematical Methods for Physicists*, 第 7 版，正交曲线坐标、矢量微分算符、Bessel 函数和球谐函数相关章节。
2. Mary L. Boas, *Mathematical Methods in the Physical Sciences*, 第 3 版，曲线坐标、分离变量和特殊函数相关章节。
3. David J. Griffiths, *Introduction to Electrodynamics*, 第 4 版，第 1 章及附录中的矢量分析公式，以及第 3 章的 Laplace 方程分离变量与多极展开。

这些来源采用的球坐标角符号可能互换。本文固定 $\theta$ 为从 $+z$ 轴量起的极角，$\phi$ 为绕 $z$ 轴的方位角；使用其他教材公式前必须先核对这一约定。

## 十、本篇结论

正交曲线坐标中的几何信息由尺度因子 $h_i$ 编码，Jacobian 为 $h_1h_2h_3$。梯度来自单位弧长上的变化，散度来自单位体积净通量，旋度来自单位面积环量；Laplacian 则由梯度再取散度得到。柱坐标和球坐标公式只在各自非退化坐标图内逐点成立，轴线、原点和极轴必须用极限、直角坐标或分布方法单独检查。

最容易误用的边界有两个：一是忽略分离常数符号而混淆普通与修正贝塞尔函数，二是把坐标分量的发散直接当作物理场发散。后续边值问题将使用这里的算符和正则性条件选择满足具体边界的解支。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.4 狄拉克 δ 函数](https://zhuanlan.zhihu.com/p/2061253558738097791)

**课程目录：** [电动力学目录](https://zhuanlan.zhihu.com/p/2062546254274601559)
<!-- zhihu-nav-bottom:end -->
