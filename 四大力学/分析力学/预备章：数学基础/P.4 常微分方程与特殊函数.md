---
created: 2026-07-18T03:15:23+08:00
modified: 2026-07-18T13:58:30+08:00
zhihu-title: 常微分方程
zhihu-topics:
  - 分析力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061652524868821844
zhihu-toc: true
zhihu-created-at: 2026-07-18 13:58
zhihu-updated-at: 2026-07-21 14:49
---
<!-- zhihu-nav-top:start -->
**上一节：** [P.3 矢量分析与张量初步](https://zhuanlan.zhihu.com/p/2061652489871664690)
<!-- zhihu-nav-top:end -->

## P.4 常微分方程与特殊函数

分析力学中的运动方程——拉格朗日方程和哈密顿正则方程——几乎全是常微分方程。简谐振子给出常系数线性方程，中心力场中的径向运动引出特殊函数，刚体转动和单摆的非线性振动则需要椭圆函数。这一节从常微分方程的通解结构讲起，建立常系数方程的完整解法，然后进入级数解法这一处理变系数方程的核心工具，接着系统介绍勒让德多项式和贝塞尔函数这两类在物理中最常见的特殊函数，最后简要引入椭圆函数。目标是让每一类方程的来龙去脉、每一步推导都清晰可见。

---

### 一、常微分方程的基本概念

**常微分方程**（ODE）是含有一个自变量 $x$ 和未知函数 $y(x)$ 及其导数的方程。**阶数**是出现的最高阶导数。最一般的 $n$ 阶 ODE 形式为

$$
F(x, y, y', y'', \ldots, y^{(n)}) = 0. \tag{P.4.1}
$$

如果方程对未知函数及其各阶导数是线性的，称为**线性 ODE**：

$$
a_n(x) y^{(n)} + a_{n-1}(x) y^{(n-1)} + \cdots + a_1(x) y' + a_0(x) y = f(x). \tag{P.4.2}
$$

$f(x) \equiv 0$ 时称为**齐次**线性 ODE；$f(x) \not\equiv 0$ 时称为**非齐次**线性 ODE。

在分析力学中，小振动理论给出常系数线性方程组，中心力场问题化为变系数二阶线性 ODE，微扰论经常需要求解非齐次方程。

---

### 二、一阶线性常微分方程

#### 2.1 齐次方程

一阶齐次线性 ODE：

$$
y' + p(x) y = 0. \tag{P.4.3}
$$

这是可分离变量的方程：

$$
\frac{dy}{y} = -p(x) dx.
$$

两边积分：

$$
\ln|y| = -\int p(x) dx + C.
$$

通解为

$$
y(x) = A e^{-\int p(x) dx}, \tag{P.4.4}
$$

其中 $A = \pm e^C$ 是任意常数。给定初始条件 $y(x_0)=y_0$，可唯一确定 $A$。

#### 2.2 非齐次方程——积分因子法

一阶非齐次线性 ODE：

$$
y' + p(x) y = q(x). \tag{P.4.5}
$$

核心方法是**积分因子法**。寻找函数 $\mu(x)$，乘以方程两边后使左边成为全导数。

设积分因子为 $\mu(x)$，则

$$
\mu y' + \mu p y = \mu q.
$$

左边应等于 $(\mu y)' = \mu y' + \mu' y$。比较 $y$ 的系数，要求 $\mu' = \mu p$。这是关于 $\mu$ 的齐次方程，解为

$$
\mu(x) = e^{\int p(x) dx}. \tag{P.4.6}
$$

乘以积分因子后，方程变为

$$
\frac{d}{dx}\left( \mu(x) y \right) = \mu(x) q(x).
$$

积分：

$$
\mu(x) y(x) = \int \mu(x) q(x) dx + C.
$$

通解为

$$
\boxed{y(x) = \frac{1}{\mu(x)} \left[ \int \mu(x) q(x) dx + C \right]}, \quad \mu(x) = e^{\int p(x) dx}. \tag{P.4.7}
$$

**例**：$y' + 2x y = x$。这里 $p(x)=2x$，$\mu(x)=e^{x^2}$。通解：

$$
y(x) = e^{-x^2} \left[ \int x e^{x^2} dx + C \right].
$$

令 $u=x^2$，$du=2x dx$，$\int x e^{x^2} dx = \frac{1}{2}e^{x^2}$。因此

$$
y(x) = e^{-x^2} \left( \frac{1}{2}e^{x^2} + C \right) = \frac{1}{2} + C e^{-x^2}.
$$

解的结构清晰：**通解 = 齐次通解 $C e^{-x^2}$ + 非齐次特解 $1/2$ **。这是一般规律。

---

### 三、二阶线性常微分方程的通解结构

#### 3.1 齐次方程

二阶齐次线性 ODE 的标准形式：

$$
y'' + P(x) y' + Q(x) y = 0. \tag{P.4.8}
$$

**叠加原理**：若 $y_1(x)$ 和 $y_2(x)$ 是解，则 $c_1 y_1(x) + c_2 y_2(x)$ 也是解（直接代入验证，各项线性）。

**线性无关与朗斯基行列式**：两个解线性无关的充要条件是它们的朗斯基（Wronskian）行列式不为零：

$$
W[y_1, y_2](x) = \begin{vmatrix} y_1 & y_2 \\ y_1' & y_2' \end{vmatrix} = y_1 y_2' - y_1' y_2 \neq 0. \tag{P.4.9}
$$

若 $y_1, y_2$ 线性无关，则齐次方程的**通解**为

$$
y_h(x) = c_1 y_1(x) + c_2 y_2(x). \tag{P.4.10}
$$

其中 $c_1, c_2$ 为任意常数。两个线性无关解构成解空间的基。

**朗斯基行列式的 Liouville 公式**：$W(x)$ 满足一阶 ODE $W' + P(x)W = 0$，因此

$$
W(x) = W(x_0) \exp\left(-\int_{x_0}^x P(t) dt\right). \tag{P.4.11}
$$

这意味着 $W(x)$ 要么恒为零（解线性相关），要么永远不为零（解线性无关）。

**降阶法**：若已知一个非零解 $y_1(x)$，可设 $y_2(x) = v(x) y_1(x)$ 代入方程求 $v(x)$。代入后 $v$ 满足一个关于 $v'$ 的一阶方程，可以积分解出 $v$：

$$
v(x) = \int \frac{W(x)}{y_1^2(x)} dx,
$$

其中 $W(x)$ 由 Liouville 公式给出。这是一个非常实用的构造第二个解的方法。

#### 3.2 非齐次方程——常数变易法

二阶非齐次线性 ODE：

$$
y'' + P(x) y' + Q(x) y = f(x). \tag{P.4.12}
$$

已知齐次方程的两个线性无关解 $y_1, y_2$，非齐次方程的**特解**可通过**常数变易法**求得。令

$$
y_p(x) = u_1(x) y_1(x) + u_2(x) y_2(x),
$$

将常数“变易”为函数。代入方程，并施加条件 $u_1' y_1 + u_2' y_2 = 0$（这个条件简化了 $y_p'$ 的表达式），可得方程组：

$$
\begin{cases}
u_1' y_1 + u_2' y_2 = 0, \\
u_1' y_1' + u_2' y_2' = f(x).
\end{cases}
$$

解出

$$
u_1' = -\frac{y_2 f}{W}, \quad u_2' = \frac{y_1 f}{W}. \tag{P.4.13}
$$

积分得 $u_1, u_2$，从而得到特解。非齐次方程的通解为

$$
y(x) = y_h(x) + y_p(x) = c_1 y_1(x) + c_2 y_2(x) + y_p(x). \tag{P.4.14}
$$

---

### 四、常系数线性方程

这是物理中最常见的一类——简谐振子、阻尼振动、受迫振动都属于此。

#### 4.1 齐次方程：特征方程法

$n$ 阶常系数齐次线性 ODE：

$$
a_n y^{(n)} + a_{n-1} y^{(n-1)} + \cdots + a_1 y' + a_0 y = 0, \quad a_n \neq 0. \tag{P.4.15}
$$

设解的形式为 $y = e^{\lambda x}$。代入后导数 $\frac{d^k}{dx^k} e^{\lambda x} = \lambda^k e^{\lambda x}$，提出公因子 $e^{\lambda x} \neq 0$，得到**特征方程**：

$$
a_n \lambda^n + a_{n-1} \lambda^{n-1} + \cdots + a_1 \lambda + a_0 = 0. \tag{P.4.16}
$$

这是 $\lambda$ 的 $n$ 次代数方程，有 $n$ 个根（实或复）。每个根对应一个基解：

- **单实根** $\lambda$：基解 $e^{\lambda x}$。
- ** $k$ 重实根** $\lambda$：基解 $e^{\lambda x}, x e^{\lambda x}, \ldots, x^{k-1} e^{\lambda x}$。
- **单共轭复根** $\alpha \pm i\beta$：基解 $e^{\alpha x}\cos\beta x, e^{\alpha x}\sin\beta x$。
- ** $k$ 重共轭复根** $\alpha \pm i\beta$：基解 $e^{\alpha x}\cos\beta x, e^{\alpha x}\sin\beta x, x e^{\alpha x}\cos\beta x, x e^{\alpha x}\sin\beta x, \ldots, x^{k-1}e^{\alpha x}\cos\beta x, x^{k-1}e^{\alpha x}\sin\beta x$。

通解是这些基解的线性组合。

#### 4.2 二阶齐次方程的详细分类

二阶常系数齐次方程 $a y'' + b y' + c y = 0$（$a \neq 0$），特征方程为 $a\lambda^2 + b\lambda + c = 0$，根为

$$
\lambda = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}.
$$

记判别式 $\Delta = b^2 - 4ac$。

**情况一：$\Delta > 0$ **（两个不等实根 $\lambda_1, \lambda_2$）

$$
y(x) = c_1 e^{\lambda_1 x} + c_2 e^{\lambda_2 x}. \tag{P.4.17}
$$

这是过阻尼振动（阻尼力大于回复力）。

**情况二：$\Delta = 0$ **（重根 $\lambda = -b/(2a)$）

$$
y(x) = (c_1 + c_2 x) e^{\lambda x}. \tag{P.4.18}
$$

这是临界阻尼。

**情况三：$\Delta < 0$ **（共轭复根 $\alpha \pm i\beta$，$\alpha=-b/(2a), \beta=\sqrt{4ac-b^2}/(2a)$）

$$
y(x) = e^{\alpha x}(c_1 \cos\beta x + c_2 \sin\beta x). \tag{P.4.19}
$$

若 $\alpha=0$（无阻尼），就是简谐振动：$y = c_1\cos\beta x + c_2\sin\beta x = A\cos(\beta x - \delta)$。

#### 4.3 非齐次常系数方程——待定系数法

非齐次方程 $a_n y^{(n)} + \cdots + a_0 y = f(x)$。已知齐次通解后，需要求一个特解。当 $f(x)$ 具有特殊形式时（指数函数、多项式、正弦余弦、或它们的乘积），用**待定系数法**：

- 若 $f(x) = P_m(x) e^{\alpha x}$（$P_m$ 为 $m$ 次多项式），设特解为 $x^s Q_m(x) e^{\alpha x}$，其中 $s$ 是 $\alpha$ 作为特征根的重数（$\alpha$ 不是特征根时 $s=0$），$Q_m(x)$ 是待定的 $m$ 次多项式。
- 若 $f(x) = e^{\alpha x}[P_m(x)\cos\beta x + Q_n(x)\sin\beta x]$，设特解为 $x^s e^{\alpha x}[R_k(x)\cos\beta x + S_k(x)\sin\beta x]$，$k=\max(m,n)$，$s$ 是 $\alpha+i\beta$ 作为特征根的重数。

代入原方程，比较同次幂系数，解出待定常数。

**例：受迫谐振子**。$y'' + \omega_0^2 y = F_0 \cos\omega t$。齐次通解 $y_h = c_1\cos\omega_0 t + c_2\sin\omega_0 t$。

若 $\omega \neq \omega_0$（非共振），设 $y_p = A\cos\omega t + B\sin\omega t$。代入：

$$
(-\omega^2 A + \omega_0^2 A)\cos\omega t + (-\omega^2 B + \omega_0^2 B)\sin\omega t = F_0\cos\omega t.
$$

比较系数：$(\omega_0^2 - \omega^2)A = F_0$，$(\omega_0^2 - \omega^2)B = 0$。因此 $A = F_0/(\omega_0^2 - \omega^2), B=0$。特解 $y_p = \frac{F_0}{\omega_0^2 - \omega^2}\cos\omega t$。

若 $\omega = \omega_0$（共振），设 $y_p = t(A\cos\omega_0 t + B\sin\omega_0 t)$。代入求得 $A=0, B=F_0/(2\omega_0)$，特解 $y_p = \frac{F_0}{2\omega_0} t\sin\omega_0 t$。振幅随时间线性增长——这就是共振的数学根源。

---

### 五、级数解法与正则奇点

常系数方程是特例。更一般的变系数方程

$$
y'' + P(x) y' + Q(x) y = 0
$$

往往没有初等函数解。**级数解法**（Frobenius 方法）是处理这类方程的系统工具。

#### 5.1 常点附近的幂级数解法

若 $P(x)$ 和 $Q(x)$ 在 $x_0$ 处解析（即可展开为收敛的泰勒级数），则称 $x_0$ 为方程的**常点**（ordinary point）。在常点附近，方程的解可以展开为幂级数：

$$
y(x) = \sum_{n=0}^{\infty} a_n (x - x_0)^n. \tag{P.4.20}
$$

代入方程，将各项整理为 $(x-x_0)$ 的幂级数，令各次幂系数为零，得到系数 $a_n$ 的递推关系。由此可以逐项求出所有系数，且可证明级数的收敛半径至少为 $x_0$ 到最近奇点的距离。

**例：艾里方程** $y'' - x y = 0$。$P(x)=0, Q(x)=-x$ 处处解析，$x_0=0$ 是常点。设 $y = \sum_{n=0}^{\infty} a_n x^n$。

$$
y'' = \sum_{n=2}^{\infty} n(n-1) a_n x^{n-2} = \sum_{n=0}^{\infty} (n+2)(n+1) a_{n+2} x^n.
$$

代入方程：

$$
\sum_{n=0}^{\infty} (n+2)(n+1) a_{n+2} x^n - \sum_{n=0}^{\infty} a_n x^{n+1} = 0.
$$

第二项中令 $m=n+1$，变为 $\sum_{m=1}^{\infty} a_{m-1} x^m = \sum_{n=1}^{\infty} a_{n-1} x^n$。合并：

$$
2 a_2 + \sum_{n=1}^{\infty} \left[ (n+2)(n+1) a_{n+2} - a_{n-1} \right] x^n = 0.
$$

因此 $a_2 = 0$，且对 $n \geq 1$：

$$
a_{n+2} = \frac{a_{n-1}}{(n+2)(n+1)}.
$$

递推关系跨三步连接。$a_0$ 和 $a_1$ 由初始条件确定（两个自由度）。$a_0 \to a_3 \to a_6 \to \cdots$；$a_1 \to a_4 \to a_7 \to \cdots$；$a_2=0 \to a_5=0 \to a_8=0 \to \cdots$。这样得到两个线性无关的级数解，它们在整复平面上收敛。

#### 5.2 正则奇点与 Frobenius 方法

若 $P(x)$ 或 $Q(x)$ 在 $x_0$ 处有奇点（极点），但 $(x-x_0)P(x)$ 和 $(x-x_0)^2 Q(x)$ 在 $x_0$ 处解析，则称 $x_0$ 为**正则奇点**（regular singular point）。很多物理方程（勒让德方程、贝塞尔方程）都有正则奇点。

设 $x_0=0$ 为正则奇点（不失一般性，可通过平移达到）。Frobenius 方法设解的形式为

$$
y(x) = x^r \sum_{n=0}^{\infty} a_n x^n = \sum_{n=0}^{\infty} a_n x^{n+r}, \quad a_0 \neq 0. \tag{P.4.21}
$$

指数 $r$ 不一定是整数。代入方程，提取最低次幂 $x^{r-2}$ 的系数，得到 $r$ 必须满足的二次方程——**指标方程**（indicial equation）：

$$
r(r-1) + p_0 r + q_0 = 0, \tag{P.4.22}
$$

其中 $p_0 = \lim_{x \to 0} x P(x)$，$q_0 = \lim_{x \to 0} x^2 Q(x)$。

指标方程给出两个根 $r_1, r_2$（设 $\text{Re}(r_1) \geq \text{Re}(r_2)$）。递推关系随后确定所有 $a_n$。根据 $r_1 - r_2$ 是否为整数，解的结构有所不同：

- ** $r_1 - r_2 \notin \mathbb{Z}$ **：两个线性无关解直接对应两个根：
  $$
  y_1(x) = x^{r_1} \sum_{n=0}^{\infty} a_n x^n, \quad y_2(x) = x^{r_2} \sum_{n=0}^{\infty} b_n x^n. \tag{P.4.23}
  $$

- ** $r_1 = r_2$ **：第二个解含有对数项：
  $$
  y_1(x) = x^{r_1} \sum_{n=0}^{\infty} a_n x^n, \quad y_2(x) = y_1(x) \ln x + x^{r_1} \sum_{n=1}^{\infty} b_n x^n. \tag{P.4.24}
  $$

- ** $r_1 - r_2 \in \mathbb{N}^+$ **：第二个解可能含有对数项。

级数解法（Frobenius 方法）是处理分析力学中球坐标和柱坐标下径向方程的核心技术。下面两节展示它在两类最重要的特殊函数上的具体应用。

---

### 六、勒让德多项式

#### 6.1 勒让德方程的来源

在球坐标中求解拉普拉斯方程 $\nabla^2 \Phi = 0$ 时，分离变量 $\Phi(r,\theta,\phi) = R(r)\Theta(\theta)\Phi(\phi)$，极角 $\theta$ 部分满足的方程经过变量代换 $x = \cos\theta$ 后变成

$$
(1 - x^2) \frac{d^2 y}{dx^2} - 2x \frac{dy}{dx} + \nu(\nu+1) y = 0, \quad x \in [-1, 1]. \tag{P.4.25}
$$

这就是**勒让德方程**。改写为标准形式：

$$
y'' - \frac{2x}{1-x^2} y' + \frac{\nu(\nu+1)}{1-x^2} y = 0.
$$

$P(x) = -2x/(1-x^2)$，$Q(x) = \nu(\nu+1)/(1-x^2)$。$x = \pm 1$ 是 $P, Q$ 的一阶极点，是正则奇点。$x=0$ 是常点。

#### 6.2 常点 $x=0$ 附近的级数解

设 $y = \sum_{n=0}^{\infty} a_n x^n$。代入勒让德方程：

$$
(1-x^2)\sum_{n=2}^{\infty} n(n-1) a_n x^{n-2} - 2x \sum_{n=1}^{\infty} n a_n x^{n-1} + \nu(\nu+1) \sum_{n=0}^{\infty} a_n x^n = 0.
$$

把各项调整到同次幂 $x^n$：

第一项：$(1-x^2) \sum n(n-1)a_n x^{n-2} = \sum n(n-1)a_n x^{n-2} - \sum n(n-1)a_n x^n$
$= \sum_{n=0}^{\infty} (n+2)(n+1)a_{n+2} x^n - \sum_{n=0}^{\infty} n(n-1)a_n x^n$。

第二项：$-2x \sum n a_n x^{n-1} = -2 \sum n a_n x^n = -\sum_{n=0}^{\infty} 2n a_n x^n$。

第三项：$\nu(\nu+1) \sum a_n x^n = \sum_{n=0}^{\infty} \nu(\nu+1) a_n x^n$。

合并 $x^n$ 系数（$n \geq 0$）：

$$
(n+2)(n+1) a_{n+2} - n(n-1)a_n - 2n a_n + \nu(\nu+1) a_n = 0.
$$

化简：

$$
(n+2)(n+1) a_{n+2} + [\nu(\nu+1) - n(n+1)] a_n = 0.
$$

注意到 $\nu(\nu+1) - n(n+1) = (\nu - n)(\nu + n + 1)$。因此递推关系为

$$
\boxed{a_{n+2} = -\frac{(\nu - n)(\nu + n + 1)}{(n+2)(n+1)} a_n}. \tag{P.4.26}
$$

$a_0, a_1$ 是两个自由参数。偶次幂由 $a_0$ 生成，奇次幂由 $a_1$ 生成：

$$
y(x) = a_0 \left[ 1 - \frac{\nu(\nu+1)}{2!}x^2 + \frac{\nu(\nu-2)(\nu+1)(\nu+3)}{4!}x^4 - \cdots \right] + a_1 \left[ x - \frac{(\nu-1)(\nu+2)}{3!}x^3 + \cdots \right].
$$

由比值判别法，这两个级数的收敛半径为 $1$（受限于 $x=\pm 1$ 处的正则奇点）。

#### 6.3 勒让德多项式 $P_\ell(x)$

当 $\nu = \ell$ 是非负整数时，递推关系显示 $a_{\ell+2} = 0$（因为 $\nu - \ell = 0$），从而所有更高次幂系数为零。级数在有限项截断，成为一个 $\ell$ 次多项式——这就是**勒让德多项式** $P_\ell(x)$。

通过选择合适的归一化约定，标准定义为 $P_\ell(1) = 1$。由此可以逐个计算出：

$$
\begin{aligned}
P_0(x) &= 1, \\
P_1(x) &= x, \\
P_2(x) &= \frac{1}{2}(3x^2 - 1), \\
P_3(x) &= \frac{1}{2}(5x^3 - 3x), \\
P_4(x) &= \frac{1}{8}(35x^4 - 30x^2 + 3), \\
P_5(x) &= \frac{1}{8}(63x^5 - 70x^3 + 15x).
\end{aligned}
$$

#### 6.4 罗德里格斯公式与生成函数

勒让德多项式有紧凑的解析表达——**罗德里格斯公式**：

$$
\boxed{P_\ell(x) = \frac{1}{2^\ell \ell!} \frac{d^\ell}{dx^\ell} (x^2 - 1)^\ell}. \tag{P.4.27}
$$

**验证**（$\ell=2$）：$\frac{1}{8}\frac{d^2}{dx^2}(x^4 - 2x^2 + 1) = \frac{1}{8}(12x^2 - 4) = \frac{3}{2}x^2 - \frac{1}{2}$，正是 $P_2(x)$。

**生成函数**（母函数）：

$$
\frac{1}{\sqrt{1 - 2x t + t^2}} = \sum_{\ell=0}^{\infty} P_\ell(x) t^\ell, \quad |t| < 1. \tag{P.4.28}
$$

这是静电学中多极展开的基础：将 $1/|\mathbf{r} - \mathbf{r}'|$ 按 $r'/r$ 展开，系数正是勒让德多项式。

#### 6.5 正交归一性

勒让德方程是斯特姆-刘维尔型方程：

$$
\frac{d}{dx}\left[ (1-x^2) \frac{dy}{dx} \right] + \nu(\nu+1) y = 0.
$$

其本征函数（勒让德多项式）在 $[-1, 1]$ 上带权函数 $w(x)=1$ 正交：

$$
\boxed{\int_{-1}^{1} P_\ell(x) P_m(x) \, dx = \frac{2}{2\ell+1} \delta_{\ell m}}. \tag{P.4.29}
$$

**推导**：分别写出 $P_\ell$ 和 $P_m$ 满足的微分方程，前乘后减，积分，利用边界条件即得正交性。归一化常数 $\frac{2}{2\ell+1}$ 可以通过罗德里格斯公式和分部积分求得。

#### 6.6 递推关系

三个重要的递推关系（在数值计算和角动量代数中非常实用）：

$$
(\ell+1) P_{\ell+1}(x) = (2\ell+1) x P_\ell(x) - \ell P_{\ell-1}(x), \tag{P.4.30}
$$
$$
P'_{\ell+1}(x) - P'_{\ell-1}(x) = (2\ell+1) P_\ell(x), \tag{P.4.31}
$$
$$
(1-x^2) P'_\ell(x) = \ell P_{\ell-1}(x) - \ell x P_\ell(x). \tag{P.4.32}
$$

#### 6.7 连带勒让德函数

在球坐标中，若方位角量子数 $m \neq 0$（$\Phi(\phi) = e^{i m\phi}$），极角方程变为连带勒让德方程：

$$
(1-x^2) y'' - 2x y' + \left[ \ell(\ell+1) - \frac{m^2}{1-x^2} \right] y = 0, \quad x=\cos\theta.
$$

解为连带勒让德函数：

$$
P_\ell^m(x) = (-1)^m (1-x^2)^{m/2} \frac{d^m}{dx^m} P_\ell(x), \quad 0 \leq m \leq \ell. \tag{P.4.33}
$$

这是量子力学中角动量本征函数 $\Theta(\theta) \propto P_\ell^m(\cos\theta)$ 的来源。球谐函数就是归一化的连带勒让德函数与 $e^{im\phi}$ 的乘积：

$$
Y_{\ell m}(\theta, \phi) = N_{\ell m} P_\ell^m(\cos\theta) e^{i m\phi}.
$$

---

### 七、贝塞尔函数

#### 7.1 贝塞尔方程的来源

在柱坐标中分离拉普拉斯方程或亥姆霍兹方程 $(\nabla^2 + k^2)\psi = 0$ 时，径向部分满足贝塞尔方程。令 $x = k r$（或直接从圆膜振动问题出发），标准形式为

$$
x^2 y'' + x y' + (x^2 - \nu^2) y = 0, \tag{P.4.34}
$$

或写为

$$
y'' + \frac{1}{x} y' + \left(1 - \frac{\nu^2}{x^2}\right) y = 0.
$$

$\nu \geq 0$ 是常数（在柱对称问题中通常是整数，在球贝塞尔函数中 $\nu = \ell+1/2$ 是半奇数）。$x=0$ 是正则奇点（$(x P(x))|_{x=0} = 1, (x^2 Q(x))|_{x=0} = -\nu^2$ 解析）。

#### 7.2 Frobenius 级数解

设 $y = x^r \sum_{n=0}^{\infty} a_n x^n = \sum_{n=0}^{\infty} a_n x^{n+r}$，$a_0 \neq 0$。代入贝塞尔方程。

先算各项：
$$
y' = \sum_{n=0}^{\infty} (n+r) a_n x^{n+r-1},
$$
$$
y'' = \sum_{n=0}^{\infty} (n+r)(n+r-1) a_n x^{n+r-2}.
$$

代入方程：

$$
\sum_{n=0}^{\infty} (n+r)(n+r-1) a_n x^{n+r} + \sum_{n=0}^{\infty} (n+r) a_n x^{n+r} + \sum_{n=0}^{\infty} a_n x^{n+r+2} - \nu^2 \sum_{n=0}^{\infty} a_n x^{n+r} = 0.
$$

前两项和第四项合并：$x^{n+r}$ 的系数为 $[(n+r)(n+r-1) + (n+r) - \nu^2] a_n = [(n+r)^2 - \nu^2] a_n$。第三项需要移位：$x^{n+r+2}$ 令 $m=n+2$ 变为 $\sum_{m=2}^{\infty} a_{m-2} x^{m+r} = \sum_{n=2}^{\infty} a_{n-2} x^{n+r}$。

总方程：

$$
\sum_{n=0}^{\infty} [(n+r)^2 - \nu^2] a_n x^{n+r} + \sum_{n=2}^{\infty} a_{n-2} x^{n+r} = 0.
$$

**最低次幂 $x^r$（$n=0$）**：系数必须为零：
$$
(r^2 - \nu^2) a_0 = 0.
$$

$a_0 \neq 0$，因此**指标方程**为

$$
r^2 - \nu^2 = 0 \quad \Rightarrow \quad r = \pm \nu. \tag{P.4.35}
$$

** $x^{r+1}$（$n=1$）**：$[(1+r)^2 - \nu^2] a_1 = 0$。对于 $r=\nu \geq 0$，$(1+\nu)^2 - \nu^2 = 1 + 2\nu \neq 0$（除了 $\nu=0$ 的特殊情况，但 $a_1$ 此时也为零），所以 $a_1 = 0$。

** $x^{n+r}$（$n \geq 2$）**：递推关系

$$
[(n+r)^2 - \nu^2] a_n + a_{n-2} = 0.
$$

取 $r = \nu$：

$$
n(n+2\nu) a_n + a_{n-2} = 0 \quad \Rightarrow \quad a_n = -\frac{a_{n-2}}{n(n+2\nu)}. \tag{P.4.36}
$$

由于 $a_1 = 0$，所有奇次系数 $a_3, a_5, \ldots$ 全为零。偶次系数由 $a_0$ 递推。

写 $n=2k$：$a_{2k} = -\frac{a_{2k-2}}{2k(2k+2\nu)} = -\frac{a_{2k-2}}{4k(k+\nu)}$。迭代：

$$
a_{2k} = \frac{(-1)^k a_0}{4^k k! (k+\nu)(k+\nu-1)\cdots(1+\nu)}.
$$

引入伽马函数 $\Gamma(z+1) = z\Gamma(z)$，分母中的乘积可写为 $\Gamma(k+\nu+1)/\Gamma(\nu+1)$。通常取归一化

$$
a_0 = \frac{1}{2^\nu \Gamma(\nu+1)}.
$$

由此定义**第一类贝塞尔函数**（$\nu \geq 0$）：

$$
\boxed{J_\nu(x) = \sum_{k=0}^{\infty} \frac{(-1)^k}{k! \, \Gamma(k+\nu+1)} \left(\frac{x}{2}\right)^{2k+\nu}}. \tag{P.4.37}
$$

当 $\nu = n$ 为非负整数时，$\Gamma(k+n+1) = (k+n)!$：

$$
J_n(x) = \sum_{k=0}^{\infty} \frac{(-1)^k}{k!(k+n)!} \left(\frac{x}{2}\right)^{2k+n}. \tag{P.4.38}
$$

#### 7.3 第二类贝塞尔函数 $Y_\nu(x)$

若 $\nu$ 不是整数，$J_\nu$ 和 $J_{-\nu}$ 线性无关，通解为 $y = c_1 J_\nu(x) + c_2 J_{-\nu}(x)$。

若 $\nu = n$ 为整数，$J_{-n}(x) = (-1)^n J_n(x)$，二者线性相关。此时需要第二类贝塞尔函数（诺伊曼函数）$Y_n(x)$，其定义取极限：

$$
Y_\nu(x) = \frac{J_\nu(x) \cos(\nu\pi) - J_{-\nu}(x)}{\sin(\nu\pi)},
$$
$$
Y_n(x) = \lim_{\nu \to n} Y_\nu(x). \tag{P.4.39}
$$

$Y_n(x)$ 在 $x \to 0$ 处有对数奇异性（来自 $J_{-\nu}$ 中 $\Gamma$ 函数的极点），这在物理问题中对应“内边界条件”要求解在原点有限，从而排除了 $Y_n$。

#### 7.4 贝塞尔函数的性质

**小 $x$ 行为**（$x \to 0$）：
$$
J_\nu(x) \approx \frac{1}{\Gamma(\nu+1)}\left(\frac{x}{2}\right)^\nu, \quad \nu > 0; \qquad J_0(x) \approx 1 - \frac{x^2}{4}. \tag{P.4.40}
$$

**大 $x$ 行为**（$x \to \infty$）：
$$
J_\nu(x) \approx \sqrt{\frac{2}{\pi x}} \cos\left(x - \frac{\nu\pi}{2} - \frac{\pi}{4}\right). \tag{P.4.41}
$$

贝塞尔函数是衰减的正弦波，振幅按 $1/\sqrt{x}$ 减小。这在柱面波传播（如水面波纹、光纤中的光场）中起核心作用。

**递推关系**：
$$
\frac{d}{dx}\left[ x^\nu J_\nu(x) \right] = x^\nu J_{\nu-1}(x), \tag{P.4.42}
$$
$$
\frac{d}{dx}\left[ x^{-\nu} J_\nu(x) \right] = -x^{-\nu} J_{\nu+1}(x), \tag{P.4.43}
$$
$$
J_{\nu-1}(x) + J_{\nu+1}(x) = \frac{2\nu}{x} J_\nu(x), \tag{P.4.44}
$$
$$
J_{\nu-1}(x) - J_{\nu+1}(x) = 2 J'_\nu(x). \tag{P.4.45}
$$

**正交性**（在有限区间 $[0,a]$ 上）：设 $j_{\nu n}$ 是 $J_\nu(x)$ 的第 $n$ 个正零点（$J_\nu(j_{\nu n}) = 0$），则函数系 $\{J_\nu(j_{\nu n} r/a)\}_{n=1}^{\infty}$ 在 $[0,a]$ 上带权 $r$ 正交：

$$
\int_0^a J_\nu\left(\frac{j_{\nu n} r}{a}\right) J_\nu\left(\frac{j_{\nu m} r}{a}\right) r \, dr = \frac{a^2}{2} [J_{\nu+1}(j_{\nu n})]^2 \delta_{nm}. \tag{P.4.46}
$$

这是柱坐标中本征函数展开的基础（如圆膜振动的简正模式、光纤中的传播模式）。

#### 7.5 修正贝塞尔函数与球贝塞尔函数

**修正贝塞尔函数**：若方程为 $x^2 y'' + x y' - (x^2 + \nu^2) y = 0$（虚宗量），解为 $I_\nu(x) = i^{-\nu} J_\nu(ix)$ 和 $K_\nu(x)$。在柱对称的扩散问题和隧道效应中出现。

**球贝塞尔函数**：在球坐标中解亥姆霍兹方程时，径向函数（$r$ 部分）满足球贝塞尔方程（$\nu = \ell+1/2$）。定义

$$
j_\ell(x) = \sqrt{\frac{\pi}{2x}} J_{\ell+1/2}(x), \quad y_\ell(x) = \sqrt{\frac{\pi}{2x}} Y_{\ell+1/2}(x). \tag{P.4.47}
$$

它们有初等表达式。例如：
$$
j_0(x) = \frac{\sin x}{x}, \quad j_1(x) = \frac{\sin x}{x^2} - \frac{\cos x}{x}, \quad j_2(x) = \left(\frac{3}{x^3} - \frac{1}{x}\right)\sin x - \frac{3}{x^2}\cos x.
$$

量子力学中球面波展开（分波法）和散射理论大量使用球贝塞尔函数。

---

### 八、椭圆函数与椭圆积分简介

#### 8.1 从单摆到椭圆积分

单摆的精确运动方程不是简谐振动。质量为 $m$、摆长为 $l$ 的单摆，运动方程为

$$
\ddot{\theta} + \frac{g}{l} \sin\theta = 0.
$$

小角度近似 $\sin\theta \approx \theta$ 给出简谐振动周期 $T_0 = 2\pi\sqrt{l/g}$。但精确的大振幅运动需要用能量积分处理。

能量守恒（设初始角度为 $\theta_0$，从静止释放）：

$$
\frac{1}{2} m l^2 \dot{\theta}^2 - m g l \cos\theta = -m g l \cos\theta_0.
$$

令 $k = \sin(\theta_0/2)$。通过三角恒等式 $\cos\theta = 1 - 2\sin^2(\theta/2)$，可得

$$
\frac{1}{4} \dot{\theta}^2 = \frac{g}{l} \left[ \sin^2\left(\frac{\theta_0}{2}\right) - \sin^2\left(\frac{\theta}{2}\right) \right].
$$

做代换 $\sin(\theta/2) = k \sin\phi$，则 $\dot{\theta} = 2k\cos\phi \, \dot{\phi} / \sqrt{1 - k^2\sin^2\phi}$，代入后积分可得周期：

$$
T = 4 \sqrt{\frac{l}{g}} \int_0^{\pi/2} \frac{d\phi}{\sqrt{1 - k^2 \sin^2\phi}}. \tag{P.4.48}
$$

这个积分不能用初等函数表示。它就是**第一类完全椭圆积分** $K(k)$：

$$
K(k) = \int_0^{\pi/2} \frac{d\phi}{\sqrt{1 - k^2 \sin^2\phi}} = \int_0^1 \frac{dt}{\sqrt{(1-t^2)(1-k^2 t^2)}}. \tag{P.4.49}
$$

单摆的精确周期为

$$
T = 4 \sqrt{\frac{l}{g}} K(k), \quad k = \sin(\theta_0/2).
$$

当 $k \to 0$（小振幅），$K(k) \approx \pi/2 + \pi k^2/8 + \cdots$，$T \approx T_0(1 + k^2/4 + \cdots)$，恢复简谐近似。

#### 8.2 椭圆积分的三种标准类型

**第一类椭圆积分**（Legendre 形式）：

$$
F(\phi, k) = \int_0^\phi \frac{d\theta}{\sqrt{1 - k^2 \sin^2\theta}}. \tag{P.4.50}
$$

完全形式 $K(k) = F(\pi/2, k)$。

**第二类椭圆积分**：

$$
E(\phi, k) = \int_0^\phi \sqrt{1 - k^2 \sin^2\theta} \, d\theta. \tag{P.4.51}
$$

完全形式 $E(k) = E(\pi/2, k)$。椭圆周长公式是 $L = 4a E(e)$（$a$ 为半长轴，$e$ 为离心率）。

**第三类椭圆积分**：

$$
\Pi(\phi, n, k) = \int_0^\phi \frac{d\theta}{(1 - n \sin^2\theta) \sqrt{1 - k^2 \sin^2\theta}}. \tag{P.4.52}
$$

#### 8.3 椭圆函数

椭圆积分的反函数是**椭圆函数**。就像正弦函数是积分 $u = \int_0^y dt/\sqrt{1-t^2}$ 的反函数 $y = \sin u$，椭圆函数来自第一类椭圆积分的反演。

雅可比椭圆函数 $\text{sn}(u, k), \text{cn}(u, k), \text{dn}(u, k)$ 由反演定义：

$$
u = \int_0^{\text{sn}(u,k)} \frac{dt}{\sqrt{(1-t^2)(1-k^2 t^2)}},
$$

且满足恒等式

$$
\text{sn}^2(u,k) + \text{cn}^2(u,k) = 1, \quad k^2 \text{sn}^2(u,k) + \text{dn}^2(u,k) = 1.
$$

当 $k \to 0$，$\text{sn}(u,0) = \sin u, \text{cn}(u,0) = \cos u, \text{dn}(u,0) = 1$。当 $k \to 1$，$\text{sn}(u,1) = \tanh u$。

椭圆函数是**双周期**的（在复平面上有两个独立的周期），这是其区别于三角函数（单周期）的根本特征。

#### 8.4 在分析力学中的出现

椭圆函数在以下经典问题中出现：

1. **单摆的精确运动**：$\sin(\theta/2) = k \, \text{sn}(\sqrt{g/l} \, t, k)$。
2. **不对称陀螺（欧拉陀螺）的转动**：角速度分量由椭圆函数表示。
3. **中心力场中的某些轨道**：如两个引力中心的问题。
4. **可积系统中的有限带势解**（如 KdV 方程的周期解）。

这些问题的共同特征是：运动方程约化为一个包含三次或四次多项式的平方根的积分，这正是椭圆积分的典型形式。

---

### 九、总结

P.4 节从基础到专门，建立了常微分方程和特殊函数的系统知识：

1. **一阶线性 ODE**：通解 $y = \frac{1}{\mu}(\int \mu q + C)$，$\mu = e^{\int p}$。解的结构为齐次通解加非齐次特解。

2. **二阶线性 ODE 的通解结构**：齐次方程的两个线性无关解构成解空间基，朗斯基行列式判断线性无关。非齐次方程用常数变易法求特解。

3. **常系数方程**：特征方程法给出指数形式的基解。重根乘 $x$，复根化为三角函数。待定系数法处理非齐次项。

4. **级数解法和 Frobenius 方法**：常点附近幂级数解，正则奇点附近 Frobenius 级数 $x^r\sum a_n x^n$。指标方程确定指数 $r$，递推关系确定系数。

5. **勒让德多项式**：球坐标拉普拉斯方程的极角部分。$P_\ell(x)$ 满足 $(1-x^2)y'' - 2x y' + \ell(\ell+1)y = 0$，在 $[-1,1]$ 上正交，有罗德里格斯公式和生成函数。连带勒让德函数 $P_\ell^m$ 给出球谐函数的极角部分。

6. **贝塞尔函数**：柱坐标拉普拉斯/亥姆霍兹方程的径向部分。$J_\nu(x)$ 由 Frobenius 级数定义，$Y_\nu(x)$ 在整数阶含对数奇异性。递推关系与三角函数递推类似。球贝塞尔函数 $j_\ell, y_\ell$ 有初等表达式。

7. **椭圆函数与椭圆积分**：来自非线性振动（如单摆精确解）的运动方程。第一类椭圆积分 $K(k)$ 给出单摆的精确周期。雅可比椭圆函数 $\text{sn}, \text{cn}, \text{dn}$ 是椭圆积分的反函数，具有双周期性。

在分析力学的后续章节中，这些工具将反复出现：小振动理论用常系数方程，中心力场用勒让德多项式和球贝塞尔函数，刚体动力学用椭圆函数，连续系统用特殊函数做本征模展开。它们是连接抽象理论框架和具体物理计算的桥梁。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.5 矩阵理论与辛矩阵](https://zhuanlan.zhihu.com/p/2061652505495352454)

**课程目录：** [分析力学目录](https://zhuanlan.zhihu.com/p/2061652577167717816)
<!-- zhihu-nav-bottom:end -->
