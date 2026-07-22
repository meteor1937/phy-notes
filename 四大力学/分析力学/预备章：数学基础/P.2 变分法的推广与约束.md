---
created: 2026-07-18T13:50:58+08:00
modified: 2026-07-18T13:54:53+08:00
zhihu-title: P.2 变分法的推广与约束
zhihu-topics:
  - 分析力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061810275808973374
zhihu-toc: true
zhihu-created-at: 2026-07-18 13:54
zhihu-updated-at: 2026-07-21 14:49
---
<!-- zhihu-nav-top:start -->
**上一节：** [P.1 变分法基础](https://zhuanlan.zhihu.com/p/2061652541293831558)
<!-- zhihu-nav-top:end -->

## P.2 变分法的推广与约束

P.1 节建立了变分法的最基本框架：单个自变量 $x$、单个因变量 $y$、无约束、最低阶导数。这足以处理最简单的变分问题（最短路径、最速降线），但物理世界远比这丰富。分析力学需要处理多个广义坐标（多函数）、多维空间（多自变量）、以及各种约束条件。这一节我们系统地将变分法推广到这些更一般的情形。

整个推广的逻辑链条是：单函数 $\to$ 多函数（方程组）$\to$ 多自变量（场论雏形）$\to$ 带约束（拉格朗日乘子法）$\to$ 等周问题与特征值问题。最后我们讨论变分法与微分方程的等价性——这是分析力学变分原理的数学根基。

---

### 一、多变量泛函（多个自变量）

#### 1.1 问题设定

P.1 中的泛函依赖于一个自变量 $x$ 和一个函数 $y(x)$。现在我们把自变量从一维推广到 $n$ 维。设 $\mathbf{x} = (x_1, x_2, \ldots, x_n) \in \Omega \subset \mathbb{R}^n$，其中 $\Omega$ 是一个有界区域，边界为 $\partial \Omega$。泛函的被积函数依赖于 $n$ 个自变量、函数 $u(\mathbf{x})$ 以及 $u$ 对各个自变量的偏导数 $u_{x_i} \equiv \partial u / \partial x_i$：

$$
J[u] = \int_\Omega F(x_1, \ldots, x_n, u, u_{x_1}, \ldots, u_{x_n}) \, d^n x. \tag{P.2.1}
$$

为书写简洁，记 $\nabla u = (u_{x_1}, \ldots, u_{x_n})$，则

$$
J[u] = \int_\Omega F(\mathbf{x}, u, \nabla u) \, d^n x.
$$

这类泛函出现在静电学（电势 $u(\mathbf{x})$ 使电场能最小）、弹性力学（位移场使弹性能最小）、以及经典场论中。

#### 1.2 变分与分部积分

引入 $u$ 的变分 $\delta u(\mathbf{x})$，它满足边界条件：在 $\partial \Omega$ 上 $\delta u = 0$（固定边界）。泛函的一阶变分为

$$
\delta J = \int_\Omega \left( \frac{\partial F}{\partial u} \delta u + \sum_{i=1}^n \frac{\partial F}{\partial u_{x_i}} \delta u_{x_i} \right) d^n x. \tag{P.2.2}
$$

关键步骤：对含有 $\delta u_{x_i}$ 的每一项做分部积分。注意 $\delta u_{x_i} = \partial (\delta u) / \partial x_i$，因此在第 $i$ 个方向上：

$$
\int_\Omega \frac{\partial F}{\partial u_{x_i}} \frac{\partial (\delta u)}{\partial x_i} \, d^n x.
$$

我们需要高维分部积分公式（即散度定理）。回顾：

$$
\int_\Omega \mathbf{v} \cdot \nabla \phi \, d^n x = \int_{\partial \Omega} \phi \, \mathbf{v} \cdot \mathbf{n} \, dS - \int_\Omega \phi \, \nabla \cdot \mathbf{v} \, d^n x,
$$

其中 $\mathbf{n}$ 是边界外法向单位向量。取 $\mathbf{v} = \frac{\partial F}{\partial u_{x_i}} \mathbf{e}_i$（第 $i$ 个分量为 $F_{u_{x_i}}$，其余为零），$\phi = \delta u$。则 $\mathbf{v} \cdot \nabla (\delta u) = F_{u_{x_i}} \delta u_{x_i}$，$\nabla \cdot \mathbf{v} = \frac{\partial}{\partial x_i} F_{u_{x_i}}$。于是

$$
\int_\Omega F_{u_{x_i}} \delta u_{x_i} \, d^n x
= \int_{\partial \Omega} \delta u \, F_{u_{x_i}} n_i \, dS - \int_\Omega \delta u \, \frac{\partial}{\partial x_i} F_{u_{x_i}} \, d^n x.
$$

对所有 $i=1,\ldots,n$ 求和，得到

$$
\sum_{i=1}^n \int_\Omega F_{u_{x_i}} \delta u_{x_i} \, d^n x
= \int_{\partial \Omega} \delta u \left( \sum_{i=1}^n F_{u_{x_i}} n_i \right) dS - \int_\Omega \delta u \left( \sum_{i=1}^n \frac{\partial}{\partial x_i} F_{u_{x_i}} \right) d^n x.
$$

边界项中 $\sum_i F_{u_{x_i}} n_i = \nabla_{\nabla u} F \cdot \mathbf{n}$，记为法向导数。在固定边界条件下 $\delta u|_{\partial \Omega} = 0$，边界项消失。

将这一结果代入 (P.2.2)：

$$
\delta J = \int_\Omega \left[ \frac{\partial F}{\partial u} - \sum_{i=1}^n \frac{\partial}{\partial x_i} \left( \frac{\partial F}{\partial u_{x_i}} \right) \right] \delta u \, d^n x.
$$

由变分法基本引理（在多维情形下同样成立：若对任意满足边界条件的光滑 $\delta u$ 积分恒为零，则被积函数恒为零），得到**多变量欧拉-拉格朗日方程**：

$$
\boxed{\frac{\partial F}{\partial u} - \sum_{i=1}^n \frac{\partial}{\partial x_i} \left( \frac{\partial F}{\partial u_{x_i}} \right) = 0}. \tag{P.2.3}
$$

#### 1.3 紧凑记号

引入梯度算子 $\nabla$ 和散度算子 $\nabla \cdot$。记 $\nabla u = (u_{x_1}, \ldots, u_{x_n})$，而 $\nabla_{\nabla u} F$ 表示 $F$ 对各个 $u_{x_i}$ 的偏导数构成的向量场。则

$$
\frac{\partial F}{\partial u} - \nabla \cdot \left( \nabla_{\nabla u} F \right) = 0. \tag{P.2.4}
$$

这一形式非常紧凑，且直接推广到场论中的 Euler-Lagrange 方程。

#### 1.4 例：Dirichlet 积分与拉普拉斯方程

取 $F = \frac{1}{2} |\nabla u|^2 = \frac{1}{2} \sum_i u_{x_i}^2$。则

$$
\frac{\partial F}{\partial u} = 0, \qquad
\frac{\partial F}{\partial u_{x_i}} = u_{x_i}.
$$

代入 (P.2.3)：

$$
0 - \sum_{i=1}^n \frac{\partial}{\partial x_i} (u_{x_i}) = -\sum_{i=1}^n u_{x_i x_i} = -\nabla^2 u = 0.
$$

所以使 Dirichlet 积分取极值的函数必须满足**拉普拉斯方程** $\nabla^2 u = 0$。这是静电学中电势满足的基本方程。

#### 1.5 若边界条件为自然边界条件

如果 $\delta u$ 在边界上任意（不固定），则边界项也必须为零：

$$
\int_{\partial \Omega} \delta u \left( \sum_{i=1}^n F_{u_{x_i}} n_i \right) dS = 0, \quad \forall \delta u|_{\partial \Omega}.
$$

由此得到自然边界条件：

$$
\boxed{\left. \sum_{i=1}^n \frac{\partial F}{\partial u_{x_i}} n_i \right|_{\partial \Omega} = 0}, \quad \text{即} \quad \left. \nabla_{\nabla u} F \cdot \mathbf{n} \right|_{\partial \Omega} = 0. \tag{P.2.5}
$$

对于 $F = \frac{1}{2}|\nabla u|^2$，这是 $\nabla u \cdot \mathbf{n} = 0$，即 Neumann 边界条件。

---

### 二、含多个函数的泛函

#### 2.1 问题设定

现在设泛函依赖于 $m$ 个函数 $y_1(x), y_2(x), \ldots, y_m(x)$ 及其一阶导数：

$$
J[y_1, \ldots, y_m] = \int_a^b F(x, y_1, \ldots, y_m, y'_1, \ldots, y'_m) \, dx. \tag{P.2.6}
$$

这是分析力学中最常见的结构：$m$ 个广义坐标 $q_i$，时间 $t$ 作为自变量，$F = L$ 是拉格朗日量。

#### 2.2 变分推导

每个函数独立地变分：$y_i \to y_i + \delta y_i$，其中 $\delta y_i(a) = \delta y_i(b) = 0$（固定端点）。泛函的一阶变分为

$$
\delta J = \int_a^b \sum_{i=1}^m \left( \frac{\partial F}{\partial y_i} \delta y_i + \frac{\partial F}{\partial y'_i} \delta y'_i \right) dx. \tag{P.2.7}
$$

对每一项分别处理。固定 $i$，对 $\delta y'_i$ 项分部积分：

$$
\int_a^b \frac{\partial F}{\partial y'_i} \delta y'_i \, dx
= \left[ \frac{\partial F}{\partial y'_i} \delta y_i \right]_a^b - \int_a^b \frac{d}{dx} \left( \frac{\partial F}{\partial y'_i} \right) \delta y_i \, dx.
$$

端点项因 $\delta y_i(a) = \delta y_i(b) = 0$ 为零。因此

$$
\delta J = \int_a^b \sum_{i=1}^m \left[ \frac{\partial F}{\partial y_i} - \frac{d}{dx} \frac{\partial F}{\partial y'_i} \right] \delta y_i \, dx.
$$

由于各个 $\delta y_i$ 是独立的（可以分别独立地扰动每个函数），要让 $\delta J = 0$ 对任意独立的 $\delta y_i$ 成立，每个括号必须分别恒为零。**具体论证**：取 $\delta y_1$ 任意非零，其余 $\delta y_j = 0$，则第一个括号必须为零。依次可得每个方程。因此得到**欧拉-拉格朗日方程组**：

$$
\boxed{\frac{\partial F}{\partial y_i} - \frac{d}{dx} \frac{\partial F}{\partial y'_i} = 0, \qquad i = 1, 2, \ldots, m.} \tag{P.2.8}
$$

这是分析力学中拉格朗日方程的标准形式——每个广义坐标 $q_i$ 对应一个方程。

#### 2.3 例：平面上的测地线

曲面上两点之间最短路径称为测地线。考虑参数形式 $(x(t), y(t))$，长度泛函为

$$
J[x, y] = \int_{t_0}^{t_1} \sqrt{\dot{x}^2 + \dot{y}^2} \, dt,
$$

其中 $\dot{x} = dx/dt, \dot{y} = dy/dt$。这里 $F = \sqrt{\dot{x}^2 + \dot{y}^2}$，$m=2$，自变量为 $t$。

对 $x$ 的欧拉-拉格朗日方程：

$$
\frac{\partial F}{\partial x} = 0, \qquad
\frac{\partial F}{\partial \dot{x}} = \frac{\dot{x}}{\sqrt{\dot{x}^2 + \dot{y}^2}}.
$$

所以

$$
\frac{d}{dt} \left( \frac{\dot{x}}{\sqrt{\dot{x}^2 + \dot{y}^2}} \right) = 0.
$$

对 $y$ 的方程完全对称。这两个方程给出 $\dot{x} / \sqrt{\dot{x}^2+\dot{y}^2} = C_1$，$\dot{y} / \sqrt{\dot{x}^2+\dot{y}^2} = C_2$，消去 $t$ 得到 $dy/dx = C_2/C_1$ 为常数，即直线。这是平面上的测地线。

#### 2.4 推广：多个自变量 + 多个函数

最一般的情形是 $m$ 个函数 $u^\alpha(\mathbf{x})$（$\alpha=1,\ldots,m$）定义在 $n$ 维区域上。这是经典场论的标准设置（如电磁场的四个势函数 $A^\mu$ 定义在四维时空中）。泛函为

$$
J[u^1,\ldots,u^m] = \int_\Omega F(\mathbf{x}, u^\alpha, \partial_i u^\alpha) \, d^n x, \quad i=1,\ldots,n; \;\alpha=1,\ldots,m. \tag{P.2.9}
$$

结合前两小节的推导，得到**场论中的欧拉-拉格朗日方程组**：

$$
\boxed{\frac{\partial F}{\partial u^\alpha} - \sum_{i=1}^n \frac{\partial}{\partial x_i} \left( \frac{\partial F}{\partial (\partial_i u^\alpha)} \right) = 0, \qquad \alpha = 1, \ldots, m.} \tag{P.2.10}
$$

这是后续章节中连续系统拉格朗日力学的数学基础。

---

### 三、带约束的变分问题——拉格朗日乘子法

物理系统常常受约束。例如悬链线的长度固定，刚体的形状不变，不可压缩流体的密度恒定。在变分法中处理约束的方法，与多元函数极值中的拉格朗日乘子法一脉相承。

#### 3.1 问题的数学形式

考虑泛函

$$
J[y] = \int_a^b F(x, y, y') \, dx,
$$

要求在**约束条件**

$$
\Phi[y] \equiv \int_a^b G(x, y, y') \, dx = C \quad (\text{常数})
$$

下求极值。这类约束称为**等周约束**（isoperimetric constraint），因历史上经典的等周问题（固定周长求最大面积）得名。

更一般地，可以有 $k$ 个约束条件 $\Phi_j[y_1,\ldots,y_m] = C_j$。

#### 3.2 拉格朗日乘子法的动机

回忆多元函数极值中的拉格朗日乘子法：求 $f(x,y)$ 在 $g(x,y)=0$ 下的极值，我们构造 $H(x,y,\lambda) = f(x,y) + \lambda g(x,y)$，然后对 $x,y$ 求无条件极值。为什么这样做？

因为约束 $g=0$ 意味着 $f$ 的梯度与 $g$ 的梯度平行：$\nabla f = -\lambda \nabla g$。将 $\lambda$ 视为独立变量，$\partial H/\partial \lambda = g = 0$ 回到约束。

变分法中的逻辑完全平行：约束 $\Phi[y] = C$ 意味着 $J$ 的“泛函梯度”与 $\Phi$ 的泛函梯度在某一点成比例。因此构造

$$
\mathcal{J}[y] = J[y] - \lambda \Phi[y] = \int_a^b \left[ F(x,y,y') - \lambda G(x,y,y') \right] dx.
$$

或者等价地定义 $F^* = F - \lambda G$。然后对 $F^*$ 写出无约束的欧拉-拉格朗日方程。

#### 3.3 严格推导

我们给出一个严谨的推导。考虑 $J[y]$ 在约束 $\Phi[y] = C$ 下的极值问题。令 $y(x)$ 为满足约束的极值函数。考虑一个依赖于两个小参数 $\epsilon_1, \epsilon_2$ 的扰动族：

$$
y_{\epsilon_1, \epsilon_2}(x) = y(x) + \epsilon_1 \eta_1(x) + \epsilon_2 \eta_2(x),
$$

其中 $\eta_1, \eta_2$ 满足 $\eta_i(a) = \eta_i(b) = 0$。约束条件对任意 $\epsilon_1, \epsilon_2$ 成立：

$$
\Phi[y_{\epsilon_1, \epsilon_2}] = C.
$$

对 $\epsilon_1, \epsilon_2$ 在原点处展开到一阶：

$$
\Phi[y] + \epsilon_1 \delta \Phi[\eta_1] + \epsilon_2 \delta \Phi[\eta_2] + \cdots = C.
$$

但 $\Phi[y] = C$，所以

$$
\epsilon_1 \delta \Phi[\eta_1] + \epsilon_2 \delta \Phi[\eta_2] + O(\epsilon^2) = 0.
$$

这对任意小的 $\epsilon_1, \epsilon_2$ 成立意味着

$$
\delta \Phi[\eta_1] = \delta \Phi[\eta_2] = 0.
$$

即：$\delta \Phi$ 在任意满足端点条件的扰动下为零——这与我们想对 $J$ 建立的极值条件结构相同。但因为约束的存在，$J$ 的极值条件不需要对所有 $\eta$ 成立，只需要对那些**不改变 $\Phi$ 的值**的扰动成立。

因此我们考虑两个线性泛函：

$$
\delta J[\eta] = \int_a^b \left( F_y - \frac{d}{dx}F_{y'} \right) \eta \, dx,
$$
$$
\delta \Phi[\eta] = \int_a^b \left( G_y - \frac{d}{dx}G_{y'} \right) \eta \, dx.
$$

如果存在常数 $\lambda$ 使得对**所有** $\eta$ 有

$$
\delta J[\eta] - \lambda \, \delta \Phi[\eta] = 0,
$$

那么这个条件既包含了约束（对满足 $\delta \Phi[\eta] = 0$ 的 $\eta$ 给出 $\delta J[\eta] = 0$），又自然确定 $\lambda$。这就是拉格朗日乘子。

上述条件等价于对泛函

$$
\mathcal{J}[y] = J[y] - \lambda \Phi[y] = \int_a^b \left( F - \lambda G \right) dx
$$

做无约束变分。因此得到**带等周约束的欧拉-拉格朗日方程**：

$$
\boxed{\frac{\partial (F - \lambda G)}{\partial y} - \frac{d}{dx} \frac{\partial (F - \lambda G)}{\partial y'} = 0}. \tag{P.2.11}
$$

其中 $\lambda$ 是常数乘子，其值由约束条件 $\Phi[y] = C$ 确定。

#### 3.4 推广到多个函数、多个约束

若泛函依赖于 $m$ 个函数 $y_i$，且有 $k$ 个等周约束

$$
\int_a^b G_j(x, y_1, \ldots, y_m, y'_1, \ldots, y'_m) \, dx = C_j, \quad j = 1, \ldots, k,
$$

则构造

$$
F^* = F - \sum_{j=1}^k \lambda_j G_j,
$$

对每个 $y_i$ 写出欧拉-拉格朗日方程，其中 $\lambda_j$ 为常数乘子。

#### 3.5 另一种约束类型：局部约束（逐点约束）

除了积分形式的等周约束，还有**局部约束**（或叫**完整约束**），形式为

$$
g(x, y_1, \ldots, y_m) = 0, \quad \forall x.
$$

这在分析力学中对应完整约束（如刚体约束、质点被限制在曲面上运动）。处理方法是：引入**作为函数的拉格朗日乘子** $\lambda(x)$（注意这里 $\lambda$ 依赖 $x$，不是常数），构造

$$
F^* = F + \lambda(x) g(x, y_1, \ldots, y_m).
$$

然后对 $y_i$ 和 $\lambda$ 写出方程。对 $y_i$ 的欧拉-拉格朗日方程为

$$
\frac{\partial F}{\partial y_i} + \lambda(x) \frac{\partial g}{\partial y_i} - \frac{d}{dx} \frac{\partial F}{\partial y'_i} = 0.
$$

对 $\lambda$ 的方程就是约束本身：$g(x, y_1, \ldots, y_m) = 0$。

这种约束在分析力学正文中会详细展开（达朗贝尔原理、拉格朗日乘子处理非完整约束）。这里我们只需认识到：乘子法的适用范围从积分约束延伸到逐点约束，乘子的性质从常数变为函数。

---

### 四、等周问题与特征值问题

#### 4.1 经典等周问题：固定周长求最大面积

**问题**：在所有长度为 $L$ 的简单闭合曲线中，找到围成面积最大的曲线。

用参数曲线 $(x(t), y(t))$，$t \in [0,1]$，闭合条件 $x(0)=x(1), y(0)=y(1)$。面积由格林公式给出：

$$
A = \frac{1}{2} \int_0^1 (x\dot{y} - y\dot{x}) \, dt.
$$

周长约束：

$$
\int_0^1 \sqrt{\dot{x}^2 + \dot{y}^2} \, dt = L.
$$

构造 $F^* = \frac{1}{2}(x\dot{y} - y\dot{x}) - \lambda \sqrt{\dot{x}^2 + \dot{y}^2}$。写出两个欧拉-拉格朗日方程：

对 $x$：

$$
\frac{\partial F^*}{\partial x} = \frac{1}{2}\dot{y}, \quad
\frac{\partial F^*}{\partial \dot{x}} = -\frac{1}{2}y - \lambda \frac{\dot{x}}{\sqrt{\dot{x}^2+\dot{y}^2}}.
$$

欧拉-拉格朗日方程：

$$
\frac{1}{2}\dot{y} - \frac{d}{dt}\left( -\frac{1}{2}y - \lambda \frac{\dot{x}}{\sqrt{\dot{x}^2+\dot{y}^2}} \right) = 0,
$$

即

$$
\frac{1}{2}\dot{y} + \frac{1}{2}\dot{y} + \frac{d}{dt}\left( \lambda \frac{\dot{x}}{\sqrt{\dot{x}^2+\dot{y}^2}} \right) = 0,
$$

化简为

$$
\dot{y} + \frac{d}{dt}\left( \lambda \frac{\dot{x}}{\sqrt{\dot{x}^2+\dot{y}^2}} \right) = 0.
$$

对 $y$ 的方程类似：

$$
-\dot{x} + \frac{d}{dt}\left( \lambda \frac{\dot{y}}{\sqrt{\dot{x}^2+\dot{y}^2}} \right) = 0.
$$

如果取弧长参数 $s$，使得 $\sqrt{\dot{x}^2+\dot{y}^2} = 1$（即 $ds = dt$），则上述方程大大简化。设 $\dot{x} = dx/ds, \dot{y} = dy/ds$，则 $\dot{x}^2+\dot{y}^2 = 1$。方程变为

$$
\frac{dy}{ds} + \frac{d}{ds}\left( \lambda \frac{dx}{ds} \right) = 0, \quad
-\frac{dx}{ds} + \frac{d}{ds}\left( \lambda \frac{dy}{ds} \right) = 0.
$$

设 $\lambda$ 为常数（通过优化可以证明），积分得

$$
y + \lambda \frac{dx}{ds} = C_1, \quad
-x + \lambda \frac{dy}{ds} = C_2.
$$

消去 $dx/ds, dy/ds$ 并利用 $(dx/ds)^2 + (dy/ds)^2 = 1$，得到

$$
(y - C_1)^2 + (x + C_2)^2 = \lambda^2.
$$

这是一族半径为 $|\lambda|$ 的圆。等周问题的答案是圆——这是变分法历史上最著名的结论之一。

#### 4.2 特征值问题作为带约束的变分问题

特征值问题可以自然地写成带约束变分的形式。考虑泛函

$$
J[y] = \int_a^b (y')^2 \, dx,
$$

在约束条件

$$
\Phi[y] = \int_a^b y^2 \, dx = 1, \qquad y(a) = y(b) = 0
$$

下求极值。构造 $F^* = (y')^2 - \lambda y^2$。欧拉-拉格朗日方程为

$$
\frac{\partial F^*}{\partial y} = -2\lambda y, \quad
\frac{\partial F^*}{\partial y'} = 2y', \quad
\frac{d}{dx}(2y') = 2y''.
$$

所以

$$
-2\lambda y - 2y'' = 0 \quad \Rightarrow \quad y'' + \lambda y = 0.
$$

这是标准的斯特姆-刘维尔特征值问题。边界条件 $y(a)=y(b)=0$ 给出了特征值 $\lambda_n$ 和特征函数 $y_n(x)$。

注意这里 $\lambda$ 的地位：它是乘子，也是特征值。变分原理告诉我们，最小的特征值 $\lambda_1$ 正是 $J[y]$ 在约束 $\int y^2 = 1$ 下的极小值。这被称为**瑞利-里兹原理**：

$$
\lambda_1 = \min_{y: \|y\|=1, \, y(a)=y(b)=0} \int_a^b (y')^2 \, dx.
$$

这一视角在量子力学（基态能量作为能量泛函的极小值）和数值方法（有限元法、里兹法）中极其重要。

---

### 五、变分法与微分方程的等价性

#### 5.1 从微分方程到变分原理

欧拉-拉格朗日方程给了我们一个方向：泛函极值 $\Rightarrow$ 微分方程。反过来，给定一个微分方程，我们能否找到一个泛函，使其欧拉-拉格朗日方程恰好是该微分方程？

答案不总是肯定，但对于许多重要的物理方程，答案是肯定的。这称为微分方程的**变分表述**。

考虑二阶线性微分方程

$$
-\frac{d}{dx}\left( p(x) \frac{dy}{dx} \right) + q(x) y = f(x), \quad x \in [a,b],
$$

带边界条件 $y(a) = y(b) = 0$（或自然边界条件）。我们声称，这个方程的变分泛函是

$$
J[y] = \int_a^b \left[ \frac{1}{2} p(x) (y')^2 + \frac{1}{2} q(x) y^2 - f(x) y \right] dx. \tag{P.2.12}
$$

我们来验证。计算被积函数 $F = \frac{1}{2}p(y')^2 + \frac{1}{2}q y^2 - f y$。

$$
\frac{\partial F}{\partial y} = q y - f, \qquad
\frac{\partial F}{\partial y'} = p y'.
$$

欧拉-拉格朗日方程：

$$
(q y - f) - \frac{d}{dx}(p y') = 0 \quad \Rightarrow \quad -\frac{d}{dx}\left(p \frac{dy}{dx}\right) + q y = f.
$$

确实得到了原微分方程。边界项 $[p y' \delta y]_a^b$：若 $y(a)=y(b)=0$ 则固定边界自动满足；若边界自由，则自然边界条件 $p y'|_{\partial} = 0$（Neumann 条件）。

#### 5.2 自伴随性与变分原理的存在性

为什么有些微分方程可以有变分表述，有些不能？关键条件是**自伴随性**（self-adjointness）。线性微分算子 $\mathcal{L}$ 若满足 $\langle \mathcal{L}u, v \rangle = \langle u, \mathcal{L}v \rangle$（在合适的内积下），则称自伴随。自伴随算子的方程 $\mathcal{L}y = f$ 总可以通过构造对应的“能量”泛函来变分表述。

上述斯特姆-刘维尔算子 $\mathcal{L} = -\frac{d}{dx}p\frac{d}{dx} + q$ 在适当边界条件下是自伴随的，因此有变分表述。物理中几乎所有基本方程——拉普拉斯方程、泊松方程、波动方程、薛定谔方程——都有对应的变分原理。这不是巧合，而是自然界的深刻特征。

#### 5.3 为什么这对分析力学至关重要

分析力学的核心是**哈密顿原理**：力学系统的真实运动，使作用量泛函

$$
S[q] = \int_{t_1}^{t_2} L(q, \dot{q}, t) \, dt
$$

取极值。这个原理的变分直接给出欧拉-拉格朗日方程：

$$
\frac{\partial L}{\partial q_i} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}_i} = 0.
$$

换言之，牛顿第二定律（通常表述为二阶微分方程）等价于一个变分原理。这一等价性赋予了力学全新的视角：我们从“解方程”变成“找极值”，从“局部因果”走向“全局变分”。这不仅是技术上的便利，更是对物理定律结构的深层理解。

---

### 六、总结

P.2 节将 P.1 的基础变分法系统推广到更复杂的结构，这些推广直接服务于分析力学的后续章节：

1. **多变量泛函**：自变量从 $x$ 推广到 $\mathbf{x} \in \mathbb{R}^n$，分部积分利用散度定理，得到 $\frac{\partial F}{\partial u} - \sum_i \partial_i (\partial F/\partial u_{x_i}) = 0$。这是场论欧拉-拉格朗日方程的基础。

2. **多函数泛函**：因变量从单个 $y$ 推广到 $m$ 个函数 $y_i$，得到 $m$ 个独立的欧拉-拉格朗日方程。这是拉格朗日力学的标准框架。

3. **带约束的变分问题**：引入拉格朗日乘子法，处理积分约束（等周约束，乘子为常数）和逐点约束（乘子为函数）。这是处理力学系统中约束问题（达朗贝尔原理、拉格朗日乘子法）的数学基础。

4. **等周问题与特征值问题**：等周问题是带约束变分的历史起源（圆是最优的）；特征值问题可视为带 $L^2$ 约束的变分问题，给出了瑞利-里兹原理。

5. **变分法与微分方程的等价性**：许多物理微分方程（特别是自伴随算子方程）都可以表示为变分原理。哈密顿原理正是这一等价性在力学中的体现——牛顿第二定律等价于作用量泛函的极值条件。

从数学工具的角度，P.1 和 P.2 合在一起，构成了整个分析力学的变分法基础。接下来我们将开始将这些工具用于物理世界，从达朗贝尔原理和哈密顿原理出发，建立拉格朗日力学体系。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.3 矢量分析与张量初步](https://zhuanlan.zhihu.com/p/2061652489871664690)

**课程目录：** [分析力学目录](https://zhuanlan.zhihu.com/p/2061652577167717816)
<!-- zhihu-nav-bottom:end -->
