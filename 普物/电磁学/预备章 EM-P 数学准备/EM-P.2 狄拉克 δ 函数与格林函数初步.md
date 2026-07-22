---
created: 2026-07-20T19:25:52+08:00
modified: 2026-07-20T20:22:15+08:00
zhihu-title: EM-P.2 狄拉克 δ 函数与格林函数初步
zhihu-topics:
  - 电磁学
zhihu-link: https://zhuanlan.zhihu.com/p/2062629264005009857
zhihu-toc: true
zhihu-created-at: 2026-07-20 21:32
zhihu-updated-at: 2026-07-21 14:41
---
<!-- zhihu-nav-top:start -->
**上一节：** [EM-P.1 矢量分析与场论](https://zhuanlan.zhihu.com/p/2062629242194679720)
<!-- zhihu-nav-top:end -->

## EM-P.2 狄拉克 $\delta$ 函数与格林函数初步

在电磁学中，我们经常需要处理“点源”问题——点电荷、点电流元等。这些源的物理尺寸远小于所关心的空间尺度，在数学上可以理想化为一个点。然而，普通的函数语言无法直接描述“在一点处无穷大，在其他处为零，但积分等于有限值”这样的对象。狄拉克 $\delta$ 函数正是为满足这一需求而引入的广义函数。它以极高的密度集中描述点源的物理量，使得我们可以将离散的点源纳入连续的场方程框架中。进而，$\delta$ 函数与格林函数方法结合，为求解泊松方程等非齐次偏微分方程提供了统一而强大的工具——将任意源的场表示为单位点源响应的线性叠加。本节从一维 $\delta$ 函数的定义出发，严格建立其基本性质，推广到三维空间，应用于点电荷密度的描述，证明泊松方程的基本解，并阐述格林函数的基本思想。

---

### 一、一维 $\delta$ 函数的定义与基本性质

#### 1.1 从物理直觉到严格定义

考虑一个沿 $x$ 轴的单位点质量分布：总质量 $m=1$ 全部集中在 $x = x_0$ 这一点上。在经典函数的意义下，质量密度函数 $\rho(x)$ 似乎应该满足：当 $x \neq x_0$ 时 $\rho(x) = 0$，但在 $x = x_0$ 处密度为无穷大，而总质量 $\int_{-\infty}^\infty \rho(x) dx = 1$。这在普通函数的框架下是矛盾的——单个点上的函数值不会影响积分（勒贝格积分意义下零测集不改变积分值）。

狄拉克于 1930 年前后引入了 $\delta$ 函数来处理这类问题。它不是传统意义上的函数，而是一种**广义函数**（或分布），通过与“足够好”的测试函数 $\phi(x)$ 的积分来定义其作用：

$$
\boxed{\int_{-\infty}^\infty \delta(x - x_0) \, \phi(x) \, dx \equiv \phi(x_0)} \tag{P.2.1}
$$

这个定义的含义是：$\delta(x - x_0)$ 是一个“抽取器”——它与任何函数相乘再积分，就精确地抽取出该函数在 $x_0$ 处的值。这一定义完全避开了 $\delta$ 函数在单点上的“值”是什么的问题，而直接规定了它作为积分算子的行为。

取 $\phi(x) = 1$，得到归一化条件

$$
\int_{-\infty}^\infty \delta(x - x_0) \, dx = 1 \tag{P.2.2}
$$

取 $\phi(x) = f(x)$（任意连续函数），则得到**筛选性质**（sifting property）：

$$
\boxed{\int_{-\infty}^\infty f(x) \, \delta(x - x_0) \, dx = f(x_0)} \tag{P.2.3}
$$

#### 1.2 多维空间中的筛选性质

对于任意包含点 $x_0$ 的积分区间 $[a, b]$，有

$$
\int_a^b f(x) \, \delta(x - x_0) \, dx = \begin{cases} f(x_0), & x_0 \in (a, b) \\ 0, & x_0 \notin [a, b] \end{cases} \tag{P.2.4}
$$

如果 $x_0$ 恰好落在区间端点，通常约定 $\delta$ 函数的积分贡献为 $f(x_0)/2$，这一约定来自对 $\delta$ 函数作为某种连续函数极限的对称考虑，但在大多数电磁学应用中积分区间明确包含或不包含源点，端点值并不关键。

#### 1.3 基本代数性质

**偶函数性**：$\delta(-x) = \delta(x)$，即 $\delta(x - x_0) = \delta(x_0 - x)$。这可以从积分变量代换 $x' = -x$ 验证：

$$
\int_{-\infty}^\infty f(x) \, \delta(-x) \, dx = \int_{-\infty}^\infty f(-x') \, \delta(x') \, dx' = f(0) = \int_{-\infty}^\infty f(x) \, \delta(x) \, dx
$$

因此 $\delta(-x)$ 与 $\delta(x)$ 作为广义函数是等价的。

**缩放性质**：对于非零常数 $a$，

$$
\boxed{\delta(a x) = \frac{1}{|a|} \, \delta(x)} \tag{P.2.5}
$$

证明：作变量代换 $y = a x$，$dx = dy/|a|$（绝对值来自积分限的方向），则

$$
\int_{-\infty}^\infty f(x) \, \delta(a x) \, dx = \int_{-\infty}^\infty f\!\left(\frac{y}{a}\right) \, \delta(y) \, \frac{dy}{|a|} = \frac{1}{|a|} f(0)
$$

而 $\int_{-\infty}^\infty f(x) \frac{1}{|a|} \delta(x) dx = \frac{1}{|a|} f(0)$，二者一致。

**与原函数的关系**：$\delta$ 函数可以形式上看作**单位阶跃函数**（Heaviside step function）的导数。阶跃函数定义为

$$
\Theta(x) = \begin{cases} 0, & x < 0 \\ 1, & x > 0 \end{cases} \tag{P.2.6}
$$

（$x=0$ 处的值通常取 $1/2$，但不影响积分。）对于任意测试函数 $\phi(x)$，分部积分给出

$$
\int_{-\infty}^\infty \Theta'(x) \, \phi(x) \, dx = \big[ \Theta(x) \phi(x) \big]_{-\infty}^\infty - \int_{-\infty}^\infty \Theta(x) \, \phi'(x) \, dx = 0 - \int_0^\infty \phi'(x) dx = \phi(0)
$$

与 $\int \delta(x) \phi(x) dx = \phi(0)$ 一致，因此在广义函数意义下

$$
\boxed{\frac{d}{dx} \Theta(x) = \delta(x)} \tag{P.2.7}
$$

这一关系在静电学中处理面电荷分布（极化面电荷、导体表面电荷）时非常有用：体密度与面密度的转换正是通过阶跃函数描述密度的突变。

---

### 二、$\delta$ 函数的常用表示

#### 2.1 矩形函数极限（脉冲函数序列）

考虑一个宽度为 $\varepsilon$、高度为 $1/\varepsilon$ 的矩形脉冲，面积为 $1$：

$$
\delta_\varepsilon(x) = \begin{cases} \dfrac{1}{\varepsilon}, & |x| < \varepsilon/2 \\[8pt] 0, & |x| > \varepsilon/2 \end{cases}
$$

当 $\varepsilon \to 0^+$ 时，该脉冲的宽度趋于零，高度趋于无穷大，但面积始终为 $1$。对任意连续函数 $f(x)$，

$$
\int_{-\infty}^\infty f(x) \, \delta_\varepsilon(x) \, dx = \frac{1}{\varepsilon} \int_{-\varepsilon/2}^{\varepsilon/2} f(x) \, dx \xrightarrow{\varepsilon \to 0} f(0)
$$

因此矩形脉冲序列 $\{\delta_\varepsilon\}$ 弱收敛于 $\delta$ 函数。

#### 2.2 高斯函数极限

高斯函数也具有归一化面积和趋近于零的宽度：

$$
\delta_\sigma(x) = \frac{1}{\sqrt{2\pi}\,\sigma} \, e^{-x^2/(2\sigma^2)}
$$

满足 $\int_{-\infty}^\infty \delta_\sigma(x) dx = 1$。当 $\sigma \to 0^+$ 时，高斯峰无限变窄变高，弱收敛于 $\delta(x)$。高斯函数的平滑性使其在场论和数值计算中非常有用。

#### 2.3 $\operatorname{sinc}$ 函数极限

另一个重要表示来自傅里叶变换：

$$
\delta_\Omega(x) = \frac{1}{\pi} \frac{\sin(\Omega x)}{x} = \frac{\Omega}{\pi} \operatorname{sinc}\!\left(\frac{\Omega x}{\pi}\right)
$$

有 $\int_{-\infty}^\infty \frac{\sin(\Omega x)}{\pi x} dx = 1$（狄利克雷积分）。当 $\Omega \to \infty$ 时，振荡频率趋于无穷，主瓣宽度趋于零，弱收敛于 $\delta(x)$。这一表示直接给出了 $\delta$ 函数的**傅里叶积分表示**：

$$
\boxed{\delta(x) = \frac{1}{2\pi} \int_{-\infty}^\infty e^{i k x} \, dk} \tag{P.2.8}
$$

因为 $\int_{-\Omega}^{\Omega} e^{i k x} dk = \frac{2 \sin(\Omega x)}{x}$，除以 $2\pi$ 即为上述 $\operatorname{sinc}$ 形式。傅里叶表示在电磁波和量子场论中极为重要：动量空间和坐标空间的转换、完备性关系都以 $\delta$ 函数为核心。

---

### 三、三维 $\delta$ 函数与点电荷密度

#### 3.1 三维 $\delta$ 函数的定义

在三维空间中，点源需要三维 $\delta$ 函数来描述。定义

$$
\delta^{(3)}(\mathbf{r} - \mathbf{r}_0) \equiv \delta(x - x_0) \, \delta(y - y_0) \, \delta(z - z_0) \tag{P.2.9}
$$

三维筛选性质为

$$
\boxed{\iiint_{\mathbb{R}^3} f(\mathbf{r}) \, \delta^{(3)}(\mathbf{r} - \mathbf{r}_0) \, dV = f(\mathbf{r}_0)} \tag{P.2.10}
$$

积分区域只要包含 $\mathbf{r}_0$ 即可，否则积分为零。

在正交曲线坐标系中，$\delta$ 函数必须配合体积元的度量系数。例如在球坐标系 $(r,\theta,\varphi)$ 中，

$$
\delta^{(3)}(\mathbf{r} - \mathbf{r}') = \frac{1}{r^2 \sin\theta} \, \delta(r - r') \, \delta(\theta - \theta') \, \delta(\varphi - \varphi') \tag{P.2.11}
$$

这是因为 $dV = r^2 \sin\theta \, dr d\theta d\varphi$，而 $\delta$ 函数本身不携带度量信息，需要因子 $1/(r^2 \sin\theta)$ 来保证在积分 $\int \cdots dV$ 下的筛选性质。

#### 3.2 点电荷的电荷密度

对于位于 $\mathbf{r}_0$ 处的点电荷 $q$，其电荷密度分布为

$$
\boxed{\rho(\mathbf{r}) = q \, \delta^{(3)}(\mathbf{r} - \mathbf{r}_0)} \tag{P.2.12}
$$

验证：在任意包含 $\mathbf{r}_0$ 的体积 $V$ 内积分，

$$
Q_{\text{enc}} = \iiint_V \rho(\mathbf{r}) \, dV = q \iiint_V \delta^{(3)}(\mathbf{r} - \mathbf{r}_0) \, dV = q
$$

如果不包含 $\mathbf{r}_0$，则积分为零。这精确地表达了“所有电荷集中于一点”的物理图像。

对于离散的点电荷系 $\{q_i, \mathbf{r}_i\}$，

$$
\rho(\mathbf{r}) = \sum_i q_i \, \delta^{(3)}(\mathbf{r} - \mathbf{r}_i) \tag{P.2.13}
$$

对于连续分布和点电荷同时存在的情况，只需将连续部分和 $\delta$ 函数部分相加即可，$\delta$ 函数将离散与连续统一于同一个积分框架中。

---

### 四、泊松方程的基本解

#### 4.1 点电荷的电势与基本解

位于原点 $\mathbf{r}=0$ 处的单位正点电荷（$q = \varepsilon_0$，为形式简便，实际点电荷 $q$ 的电势为 $q/(4\pi\varepsilon_0 r)$）在空间中产生的静电势为

$$
\Phi(\mathbf{r}) = \frac{1}{4\pi} \frac{1}{r}, \quad r = |\mathbf{r}| \neq 0 \tag{P.2.14}
$$

（取 $\varepsilon_0 = 1$ 的自然单位制，或理解为乘以适当常数。）该电势满足拉普拉斯方程 $\nabla^2 \Phi = 0$（$r \neq 0$），但在原点 $r=0$ 处有奇异性，那里是源点所在。

我们要证明：$\Phi(\mathbf{r}) = 1/(4\pi r)$ 是泊松方程

$$
\nabla^2 \Phi = -\delta^{(3)}(\mathbf{r}) \tag{P.2.15}
$$

的解。即

$$
\boxed{\nabla^2 \left( \frac{1}{r} \right) = -4\pi \, \delta^{(3)}(\mathbf{r})} \tag{P.2.16}
$$

这是三维空间中标量泊松方程的基本解，是格林函数方法的核心。

#### 4.2 证明 $\nabla^2 (1/r) = -4\pi \delta^{(3)}(\mathbf{r})$

**第一步**：验证 $r \neq 0$ 时 $\nabla^2 (1/r) = 0$。

在球坐标系中，对于只依赖于 $r$ 的函数 $f(r)$，拉普拉斯算子为

$$
\nabla^2 f(r) = \frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{df}{dr} \right)
$$

令 $f(r) = 1/r$，则

$$
\frac{df}{dr} = -\frac{1}{r^2}, \quad r^2 \frac{df}{dr} = -1
$$

$$
\frac{d}{dr} \left( r^2 \frac{df}{dr} \right) = \frac{d}{dr}(-1) = 0
$$

所以

$$
\nabla^2 \left( \frac{1}{r} \right) = \frac{1}{r^2} \cdot 0 = 0, \quad \forall r \neq 0
$$

这验证了除原点外处处满足拉普拉斯方程。

**第二步**：分析原点处的奇异性。

由于 $1/r$ 在 $r=0$ 处不可导（发散），我们需要从积分意义下考察 $\nabla^2 (1/r)$。将 $\nabla^2 (1/r)$ 对任意包含原点的体积 $V$ 积分，利用高斯散度定理：

$$
\iiint_V \nabla^2 \left( \frac{1}{r} \right) dV = \iiint_V \nabla \cdot \left[ \nabla \left( \frac{1}{r} \right) \right] dV = \oiint_{\partial V} \nabla \left( \frac{1}{r} \right) \cdot d\mathbf{A}
$$

其中 $\partial V$ 是 $V$ 的边界闭合曲面。取 $V$ 为以原点为球心、半径 $R$ 的球体 $B_R$，则边界为球面 $S_R$。在球面上，外法向 $\hat{\mathbf{n}} = \hat{\mathbf{r}}$，$d\mathbf{A} = \hat{\mathbf{r}} R^2 \sin\theta d\theta d\varphi = \hat{\mathbf{r}} R^2 d\Omega$（$d\Omega = \sin\theta d\theta d\varphi$ 为立体角元）。而

$$
\nabla \left( \frac{1}{r} \right) = -\frac{1}{r^2} \hat{\mathbf{r}}
$$

在 $r = R$ 处，

$$
\nabla \left( \frac{1}{r} \right) \cdot d\mathbf{A} = \left( -\frac{1}{R^2} \hat{\mathbf{r}} \right) \cdot \left( \hat{\mathbf{r}} R^2 d\Omega \right) = -d\Omega
$$

对全立体角积分：

$$
\oiint_{S_R} \nabla \left( \frac{1}{r} \right) \cdot d\mathbf{A} = -\int_0^{4\pi} d\Omega = -4\pi
$$

该积分值与球半径 $R$ 无关！对于不包含原点的体积，边界曲面上的面积分由于内部无奇点，散度定理给出与体积分一致的值——体积分应为零，面积分也是零（因为 $\nabla^2(1/r) = 0$ 在体积内处处成立）。

因此，函数 $\nabla^2(1/r)$ 的全体积贡献全部来源于原点：在原点以外处处为零，但在包含原点的任意体积上的积分为 $-4\pi$。这正是 $-4\pi \, \delta^{(3)}(\mathbf{r})$ 的定义行为。由此严格地得出

$$
\nabla^2 \left( \frac{1}{r} \right) = -4\pi \, \delta^{(3)}(\mathbf{r}) \tag{P.2.17}
$$

更一般地，若源点不在原点而在 $\mathbf{r}'$，则

$$
\boxed{\nabla^2 \left( \frac{1}{|\mathbf{r} - \mathbf{r}'|} \right) = -4\pi \, \delta^{(3)}(\mathbf{r} - \mathbf{r}')} \tag{P.2.18}
$$

#### 4.3 泊松方程的通解

对于一般的电荷分布 $\rho(\mathbf{r})$，电势满足泊松方程

$$
\nabla^2 V(\mathbf{r}) = -\frac{\rho(\mathbf{r})}{\varepsilon_0}
$$

将 $\rho(\mathbf{r})$ 看作无穷多个点电荷的叠加：$\rho(\mathbf{r}) = \int \rho(\mathbf{r}') \, \delta^{(3)}(\mathbf{r} - \mathbf{r}') \, dV'$。由于泊松方程是线性的，每个点源 $\rho(\mathbf{r}') \delta^{(3)}(\mathbf{r} - \mathbf{r}') dV'$ 产生的电势为

$$
dV(\mathbf{r}) = \frac{\rho(\mathbf{r}') dV'}{4\pi\varepsilon_0} \frac{1}{|\mathbf{r} - \mathbf{r}'|}
$$

叠加所有点源的贡献，得到

$$
\boxed{V(\mathbf{r}) = \frac{1}{4\pi\varepsilon_0} \iiint \frac{\rho(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} \, dV'} \tag{P.2.19}
$$

这正是我们熟知的静电势积分公式。这一推导过程就是格林函数方法的基本逻辑。

---

### 五、格林函数的基本思想

#### 5.1 问题的提法

考虑一个一般的线性非齐次偏微分方程（以泊松方程为例）：

$$
\mathcal{L} u(\mathbf{r}) = f(\mathbf{r}) \tag{P.2.20}
$$

其中 $\mathcal{L}$ 是线性微分算子（静电学中 $\mathcal{L} = \nabla^2$，$f = -\rho/\varepsilon_0$），$f(\mathbf{r})$ 是已知的源分布。我们要求的场 $u(\mathbf{r})$ 满足一定的边界条件。

格林函数 $G(\mathbf{r}, \mathbf{r}')$ 定义为位于 $\mathbf{r}'$ 处的**单位点源**在 $\mathbf{r}$ 处产生的场，即满足

$$
\boxed{\mathcal{L} G(\mathbf{r}, \mathbf{r}') = \delta^{(3)}(\mathbf{r} - \mathbf{r}')} \tag{P.2.21}
$$

（对于非自伴算子，源点所在位置可能是 $\mathbf{r}$ 或 $\mathbf{r}'$，需明确。这里的写法符合静电学惯例：$\mathcal{L}$ 作用于 $\mathbf{r}$。）

#### 5.2 解的叠加表示

一旦求得格林函数 $G(\mathbf{r}, \mathbf{r}')$，由于方程的线性性，任意源分布 $f(\mathbf{r})$ 产生的场是各点单位点源响应的叠加：

$$
\boxed{u(\mathbf{r}) = \iiint G(\mathbf{r}, \mathbf{r}') \, f(\mathbf{r}') \, dV'} \tag{P.2.22}
$$

这可以形式地验证：将 $\mathcal{L}$ 作用于 (P.2.22) 两侧，并利用 (P.2.21)：

$$
\mathcal{L} u(\mathbf{r}) = \iiint \mathcal{L} G(\mathbf{r}, \mathbf{r}') \, f(\mathbf{r}') \, dV' = \iiint \delta^{(3)}(\mathbf{r} - \mathbf{r}') \, f(\mathbf{r}') \, dV' = f(\mathbf{r})
$$

符合原方程。

#### 5.3 自由空间格林函数与边界条件

对于无界自由空间（边界条件为无穷远势场趋于零），泊松方程的格林函数为

$$
G_0(\mathbf{r}, \mathbf{r}') = -\frac{1}{4\pi} \frac{1}{|\mathbf{r} - \mathbf{r}'|} \tag{P.2.23}
$$

因为 $\nabla^2 G_0 = \delta^{(3)}(\mathbf{r} - \mathbf{r}')$（与 (P.2.18) 差一符号，注意 $\mathcal{L} = \nabla^2$ 时 $G_0$ 前系数为 $-1/(4\pi)$）。

对于有界区域或特定边界条件（如导体边界、周期性边界），格林函数需在上述自由空间格林函数的基础上叠加一个齐次方程的解（即 $\mathcal{L} G_h = 0$），使得总 $G$ 满足指定的边界条件。这正是镜像法、本征函数展开法等静边值问题求解方法的数学本质。

格林函数方法不仅限于静电学，在恒定电流、静磁学、电磁波辐射（推迟格林函数）、热传导和量子场论中都是核心求解工具。它以“单位点源响应”为基元，将复杂的场方程求解转化为几何与边界条件的构造问题，体现了线性物理体系最为深刻的叠加原理。

---

狄拉克 $\delta$ 函数为描述点源提供了严密的数学语言，将离散的点电荷自然地纳入连续积分框架。从一维到三维的定义、筛选性质到泊松方程基本解 $\nabla^2 (1/r) = -4\pi \delta^{(3)}(\mathbf{r})$ 的严格证明，$\delta$ 函数将点源的物理直观与偏微分方程的数学结构联系起来。格林函数则将这一思想系统化——把任意源分布产生的场看作点源响应的叠加，为线性场方程的求解提供了统一范式。这些工具在后续的静电势、磁矢势、推迟势和电磁波辐射等章节中将反复出现，是电磁理论大厦的坚实数学地基。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [EM-P.3 复数与复变函数初步](https://zhuanlan.zhihu.com/p/2062629304270431145)

**课程目录：** [电磁学目录](https://zhuanlan.zhihu.com/p/2062618184079954663)
<!-- zhihu-nav-bottom:end -->
