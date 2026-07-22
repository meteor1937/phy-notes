---
created: 2026-07-18T03:15:18+08:00
modified: 2026-07-18T13:50:25+08:00
zhihu-title: 变分法
zhihu-topics:
  - 分析力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061652541293831558
zhihu-toc: true
zhihu-created-at: 2026-07-18 13:50
zhihu-updated-at: 2026-07-21 14:49
---
## P.1 变分法基础

分析力学的全部理论大厦，奠基于一个简单而深刻的原理——最小作用量原理：自然界的运动规律，等价于某个被称为“作用量”的量取极值。这句话里，“极值”需要精确的数学定义。微积分教我们求函数的极值，变分法则教我们求**泛函**的极值——泛函是以函数为自变量的映射。

本节从零开始建立变分法的完整框架：先理解泛函是什么，然后定义变分运算 $\delta$，接着一步一步推导出欧拉-拉格朗日方程，最后处理边界条件和高阶导数的推广。

---

### 一、泛函的定义与基本例子

**泛函**是以函数为“自变量”、以实数为输出的映射：

$$
J: \mathcal{F} \to \mathbb{R}, \qquad y(x) \mapsto J[y].
$$

其中 $\mathcal{F}$ 是某个函数空间。微积分关心 $y$ 在不同点 $x$ 处的值；变分法关心整条曲线 $y(x)$ 的形状。泛函 $J[y]$ 的值依赖于 $y$ 在整个区间 $[a,b]$ 上的行为，而非任何单独一点。

下面三个例子展示了变分法处理的问题类型。

**例子一：两点间最短路径**

平面上两点 $(a, y_a)$ 和 $(b, y_b)$，连接它们的曲线长度为

$$
J[y] = \int_a^b \sqrt{1 + (y')^2} \, dx, \qquad y(a) = y_a, \quad y(b) = y_b.
$$

这是一个泛函。直觉告诉我们答案是一条直线，变分法将严格证明这一点。

**例子二：最速降线**

质点从点 $A$ 静止释放，沿光滑曲线滑到点 $B$（$B$ 在 $A$ 斜下方），求下滑时间最短的曲线。能量守恒给出 $v = \sqrt{2gy}$，下滑时间为

$$
J[y] = \int_a^b \frac{\sqrt{1 + (y')^2}}{\sqrt{2gy}} \, dx.
$$

答案是摆线。这是变分法的历史起源问题，由伯努利在 1696 年提出。

**例子三：悬链线**

柔软均匀的链悬挂在两点之间，在重力作用下平衡时，重力势能最小。势能泛函为

$$
J[y] = \int_a^b y \sqrt{1 + (y')^2} \, dx,
$$

约束是链长固定 $\int_a^b \sqrt{1 + (y')^2} \, dx = L$。这是带约束的变分问题。

这三个问题的共同数学结构是

$$
J[y] = \int_a^b F(x, y(x), y'(x)) \, dx. \tag{P.1.1}
$$

被积函数 $F$ 依赖于 $x$、$y(x)$ 和 $y'(x)$。长度泛函对应 $F = \sqrt{1 + (y')^2}$，最速降线对应 $F = \sqrt{1 + (y')^2}/\sqrt{2gy}$，悬链线对应 $F = y\sqrt{1 + (y')^2}$。

问题归结为：给定形如 (P.1.1) 的泛函，如何找到使 $J$ 取极值的函数 $y(x)$？

---

### 二、变分运算 $\delta$ 的定义与规则

微积分中，求函数极值时我们给自变量一个微小增量 $dx$，然后看函数值的增量 $df = f'(x)dx$。变分法中，我们给**函数本身**一个微小扰动 $\delta y(x)$ —— $\delta y$ 是整个函数的改变，在每一 $x$ 处定义。

**形式化定义**：考虑一族函数

$$
y_\epsilon(x) = y(x) + \epsilon \eta(x),
$$

其中 $\eta(x)$ 是一个固定的函数（代表扰动的“形状”），$\epsilon$ 是一个小参数。$\delta y$ 就是这族函数对 $\epsilon$ 的一阶变化：

$$
\delta y = \left. \frac{d}{d\epsilon} \right|_{\epsilon=0} y_\epsilon = \eta(x).
$$

类似地，$y'$ 的变分为

$$
\delta y' = \left. \frac{d}{d\epsilon} \right|_{\epsilon=0} y'_\epsilon = \eta'(x).
$$

**变分运算的基本规则**：

1. **线性性**：
   $$
   \delta(c_1 y_1 + c_2 y_2) = c_1 \delta y_1 + c_2 \delta y_2.
   $$

2. **乘积法则**：
   $$
   \delta(y_1 y_2) = (\delta y_1) y_2 + y_1 (\delta y_2).
   $$

3. **链式法则**：对任意可微函数 $g$，
   $$
   \delta g(y) = g'(y) \delta y.
   $$

4. **变分与求导可交换**——最关键的一条规则：
   $$
   \delta\left( \frac{dy}{dx} \right) = \frac{d}{dx}(\delta y). \tag{P.1.2}
   $$

   **证明**：
   $$
   \delta(y') = \left. \frac{d}{d\epsilon} \right|_{\epsilon=0} \frac{d}{dx} y_\epsilon
   = \frac{d}{dx} \left. \frac{d}{d\epsilon} \right|_{\epsilon=0} y_\epsilon
   = \frac{d}{dx}(\delta y).
   $$
   交换 $\epsilon$ 导数和 $x$ 导数的次序依赖于 $y_\epsilon$ 足够光滑，我们默认满足。

**泛函的变分**定义为 $J$ 对函数扰动的线性响应：

$$
\delta J = \left. \frac{d}{d\epsilon} \right|_{\epsilon=0} J[y + \epsilon \eta].
$$

将 $F(x, y+\epsilon\eta, y'+\epsilon\eta')$ 在 $(y, y')$ 处泰勒展开到一阶：

$$
F(x, y+\epsilon\eta, y'+\epsilon\eta') = F(x, y, y') + \epsilon\left( \frac{\partial F}{\partial y} \eta + \frac{\partial F}{\partial y'} \eta' \right) + O(\epsilon^2).
$$

对 $\epsilon$ 求导并取 $\epsilon=0$，得到

$$
\delta J = \int_a^b \left( \frac{\partial F}{\partial y} \delta y + \frac{\partial F}{\partial y'} \delta y' \right) dx. \tag{P.1.3}
$$

这是变分法的**基本出发点**。它把泛函的变分表示为 $\delta y$ 和 $\delta y'$ 的线性组合。但 $\delta y$ 和 $\delta y'$ 不是独立的—— $\delta y' = (\delta y)'$。下一步我们需要通过分部积分，把 $\delta y'$ 项转化为 $\delta y$，从而提取公因子。

---

### 三、欧拉-拉格朗日方程的完整推导

要求 $J[y]$ 取极值，一阶必要条件为

$$
\delta J = 0.
$$

我们要从这一条件推导出 $y(x)$ 必须满足的微分方程。

**第一步：分部积分**

对 (P.1.3) 中的第二项做分部积分。记 $F_y = \partial F/\partial y$，$F_{y'} = \partial F/\partial y'$。则

$$
\int_a^b F_{y'} \delta y' \, dx.
$$

令 $u = F_{y'}$，$dv = \delta y' \, dx = d(\delta y)$，于是 $du = \frac{d}{dx}F_{y'} \, dx$，$v = \delta y$。分部积分公式 $\int u \, dv = uv - \int v \, du$ 给出

$$
\int_a^b F_{y'} \delta y' \, dx = \Big[ F_{y'} \delta y \Big]_a^b - \int_a^b \frac{d}{dx}\left( F_{y'} \right) \delta y \, dx. \tag{P.1.4}
$$

**第二步：处理边界项**

我们首先讨论**固定端点**的情形：$y(a)$ 和 $y(b)$ 是给定的，因此 $\delta y(a) = \delta y(b) = 0$。边界项 $[F_{y'} \delta y]_a^b$ 消失为零。

将 (P.1.4) 代回 (P.1.3)：

$$
\delta J = \int_a^b \left[ F_y - \frac{d}{dx} F_{y'} \right] \delta y \, dx. \tag{P.1.5}
$$

**第三步：变分法基本引理**

我们现在有

$$
\int_a^b \left[ F_y - \frac{d}{dx} F_{y'} \right] \delta y \, dx = 0,
$$

且这一等式必须对**任意**满足 $\delta y(a) = \delta y(b) = 0$ 的光滑函数 $\delta y(x)$ 成立。

**变分法基本引理**：若 $\phi(x)$ 在 $[a,b]$ 上连续，且对任意满足 $\eta(a)=\eta(b)=0$ 的光滑函数 $\eta(x)$ 均有

$$
\int_a^b \phi(x) \eta(x) \, dx = 0,
$$

则 $\phi(x) \equiv 0$ 在 $(a,b)$ 上恒成立。

**证明**（反证法）：假设存在 $x_0 \in (a,b)$ 使得 $\phi(x_0) \neq 0$。不失一般性，设 $\phi(x_0) > 0$。由 $\phi$ 的连续性，存在 $\epsilon > 0$ 使得在 $(x_0-\epsilon, x_0+\epsilon) \subset (a,b)$ 上 $\phi(x) > 0$。构造一个“隆起函数”：

$$
\eta(x) = \begin{cases}
[(x - x_0 + \epsilon)(x - x_0 - \epsilon)]^2, & x \in (x_0-\epsilon, x_0+\epsilon), \\
0, & \text{否则}.
\end{cases}
$$

这个函数在 $(x_0-\epsilon, x_0+\epsilon)$ 内严格为正，在区间外及端点处为零，且光滑。则

$$
\int_a^b \phi(x)\eta(x) \, dx = \int_{x_0-\epsilon}^{x_0+\epsilon} \phi(x)\eta(x) \, dx > 0,
$$

与假设矛盾。$\square$

应用这个引理，取 $\phi(x) = F_y - \frac{d}{dx}F_{y'}$，$\eta = \delta y$，立即得到

$$
\boxed{\frac{\partial F}{\partial y} - \frac{d}{dx} \frac{\partial F}{\partial y'} = 0}. \tag{P.1.6}
$$

这就是**欧拉-拉格朗日方程**。它把“找使泛函取极值的函数”这个变分问题，转化为了一个微分方程问题。

**检验：最短路径问题**

取 $F = \sqrt{1 + (y')^2}$。则

$$
\frac{\partial F}{\partial y} = 0, \qquad
\frac{\partial F}{\partial y'} = \frac{y'}{\sqrt{1 + (y')^2}}.
$$

欧拉-拉格朗日方程给出

$$
\frac{d}{dx}\left( \frac{y'}{\sqrt{1 + (y')^2}} \right) = 0.
$$

所以 $\frac{y'}{\sqrt{1 + (y')^2}} = C$（常数）。两边平方并整理：$(y')^2 = C^2[1 + (y')^2]$，即 $(y')^2(1-C^2) = C^2$，所以 $y' = C/\sqrt{1-C^2} = \text{常数}$。积分得 $y = kx + b$，是一条直线。变分法严格证明了我们朴素的几何直觉。

**补充：$\frac{d}{dx}$ 是全导数**

注意 (P.1.6) 中 $\frac{d}{dx}$ 是**全导数**，而非偏导数。因为 $F_{y'} = F_{y'}(x, y(x), y'(x))$，其对 $x$ 的全导数为

$$
\frac{d}{dx} F_{y'} = \frac{\partial F_{y'}}{\partial x} + \frac{\partial F_{y'}}{\partial y} y' + \frac{\partial F_{y'}}{\partial y'} y''.
$$

展开后的欧拉-拉格朗日方程是一个二阶常微分方程：

$$
F_y - F_{y'x} - F_{y'y} y' - F_{y'y'} y'' = 0.
$$

---

### 四、边界条件的完整讨论

上面的推导假设了端点固定。现在我们考虑更一般的情况。分部积分后，完整形式为

$$
\delta J = \Big[ F_{y'} \delta y \Big]_a^b + \int_a^b \left( F_y - \frac{d}{dx} F_{y'} \right) \delta y \, dx. \tag{P.1.7}
$$

$\delta J = 0$ 要求内部积分和边界项分别处理。$\delta y$ 在内部任意，欧拉-拉格朗日方程仍然成立。边界项则给出额外条件。

**情况一：两端固定（Dirichlet 边界条件）**

$y(a), y(b)$ 给定 $\Rightarrow$ $\delta y(a) = \delta y(b) = 0$。边界项自动为零。只需满足欧拉-拉格朗日方程。

**情况二：一端固定、一端自由**

设 $y(a)$ 固定，$y(b)$ 自由。则 $\delta y(a) = 0$，但 $\delta y(b)$ 任意。边界项变为

$$
\Big[ F_{y'} \delta y \Big]_a^b = F_{y'}\big|_{x=b} \, \delta y(b) - F_{y'}\big|_{x=a} \, \delta y(a) = F_{y'}\big|_{x=b} \, \delta y(b).
$$

$\delta J = 0$ 要求这一项为零。由于 $\delta y(b)$ 任意，必须

$$
\boxed{\left. \frac{\partial F}{\partial y'} \right|_{x=b} = 0}. \tag{P.1.8}
$$

这称为**自然边界条件**（natural boundary condition）或横截性条件（transversality condition）。

**情况三：两端都自由**

$$
\left. \frac{\partial F}{\partial y'} \right|_{x=a} = 0, \qquad
\left. \frac{\partial F}{\partial y'} \right|_{x=b} = 0.
$$

**实例**：取 $F = \frac{1}{2}(y')^2$，欧拉-拉格朗日方程为 $y'' = 0$，通解 $y = Cx + D$。

- 两端固定 $y(a)=y_a, y(b)=y_b$：$C,D$ 唯一确定。
- $y(a)$ 固定、$y(b)$ 自由：自然边界条件 $F_{y'}|_{x=b} = y'(b) = 0$，得 $C=0$。$y(a)=y_a$ 给出 $D=y_a$。解为 $y(x) = y_a$，一条水平直线。物理直觉：没有端点约束时，使曲线“斜率平方积分”最小的就是不要有斜率。

---

### 五、含高阶导数的泛函

有些物理问题（如弹性梁的弯曲能）的被积函数包含二阶或更高阶导数：

$$
J[y] = \int_a^b F(x, y, y', y'', \ldots, y^{(n)}) \, dx.
$$

我们先推导 $n=2$ 的情况，再推广到一般。

设 $J[y] = \int_a^b F(x, y, y', y'') \, dx$，端点条件为 $y(a), y(b), y'(a), y'(b)$ 均固定，因此 $\delta y$ 和 $\delta y'$ 在端点处均为零。

写出变分：

$$
\delta J = \int_a^b \left( F_y \delta y + F_{y'} \delta y' + F_{y''} \delta y'' \right) dx. \tag{P.1.9}
$$

**第一项** $F_y \delta y$：不需处理。

**第二项** $F_{y'} \delta y'$：分部积分一次。

$$
\int_a^b F_{y'} \delta y' \, dx = \Big[ F_{y'} \delta y \Big]_a^b - \int_a^b \frac{d}{dx}(F_{y'}) \delta y \, dx.
$$

边界项因 $\delta y(a) = \delta y(b) = 0$ 为零。

**第三项** $F_{y''} \delta y''$：需要分部积分**两次**。

第一次：

$$
\int_a^b F_{y''} \delta y'' \, dx = \Big[ F_{y''} \delta y' \Big]_a^b - \int_a^b \frac{d}{dx}(F_{y''}) \delta y' \, dx.
$$

边界项因 $\delta y'(a) = \delta y'(b) = 0$ 为零。对剩余积分做第二次分部积分：

$$
-\int_a^b \frac{d}{dx}(F_{y''}) \delta y' \, dx
= -\Big[ \frac{d}{dx}(F_{y''}) \delta y \Big]_a^b + \int_a^b \frac{d^2}{dx^2}(F_{y''}) \delta y \, dx.
$$

边界项因 $\delta y(a) = \delta y(b) = 0$ 再次为零。

综上，

$$
\int_a^b F_{y''} \delta y'' \, dx = \int_a^b \frac{d^2}{dx^2}(F_{y''}) \delta y \, dx.
$$

**合并三项**：

$$
\delta J = \int_a^b \left( F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} \right) \delta y \, dx.
$$

由变分法基本引理，得到含二阶导数的欧拉-拉格朗日方程：

$$
\boxed{\frac{\partial F}{\partial y} - \frac{d}{dx}\frac{\partial F}{\partial y'} + \frac{d^2}{dx^2}\frac{\partial F}{\partial y''} = 0}. \tag{P.1.10}
$$

**检验**：取 $F = \frac{1}{2}(y'')^2$。则 $F_y = 0$，$F_{y'} = 0$，$F_{y''} = y''$。方程给出 $\frac{d^2}{dx^2}(y'') = 0$，即 $y^{(4)} = 0$。通解为三次多项式 $y = ax^3 + bx^2 + cx + d$，四个积分常数由四个边界条件 $y(a), y(b), y'(a), y'(b)$ 确定。这是合理的：使梁的弯曲能 $\int (y'')^2 dx$ 最小的形状，在没有额外约束时就是三次样条。

**推广到 $n$ 阶**：

观察 $k=0,1,2$ 项的符号规律：$F_y$ 前系数 $+1 = (-1)^0$，$-\frac{d}{dx}F_{y'}$ 前系数 $-1 = (-1)^1$，$+\frac{d^2}{dx^2}F_{y''}$ 前系数 $+1 = (-1)^2$。每次分部积分多出一个负号。对一般 $n$ 阶，结果为

$$
\boxed{\sum_{k=0}^n (-1)^k \frac{d^k}{dx^k} \frac{\partial F}{\partial y^{(k)}} = 0}. \tag{P.1.11}
$$

其中 $y^{(0)} = y$。记忆方法：$y^{(k)}$ 对应项为 $(-1)^k$ 乘以对 $x$ 的 $k$ 阶全导数。

边界条件也随之增加：$n$ 阶导数需要固定 $y, y', \ldots, y^{(n-1)}$ 在两端点，否则会导出相应的自然边界条件。

---

### 六、首次积分与特殊情况

** $F$ 不显含 $y$ **：若 $\partial F/\partial y = 0$，欧拉-拉格朗日方程简化为

$$
\frac{d}{dx} F_{y'} = 0 \quad \Rightarrow \quad F_{y'} = \text{常数}.
$$

这是最短路径例子的情况。常被称为“循环坐标”对应的守恒量——这是分析力学中诺特定理的前身。

** $F$ 不显含 $x$ **：若 $\partial F/\partial x = 0$，存在一个“首次积分”（Beltrami 恒等式）：

$$
F - y' \frac{\partial F}{\partial y'} = \text{常数}. \tag{P.1.12}
$$

**推导**：计算 $\frac{d}{dx}(F - y' F_{y'})$：

$$
\begin{aligned}
\frac{d}{dx}(F - y' F_{y'})
&= \frac{\partial F}{\partial x} + F_y y' + F_{y'} y'' - y'' F_{y'} - y' \frac{d}{dx}F_{y'} \\
&= \frac{\partial F}{\partial x} + y'\left( F_y - \frac{d}{dx}F_{y'} \right).
\end{aligned}
$$

如果 $F$ 不显含 $x$（$\partial F/\partial x = 0$），且 $y$ 满足欧拉-拉格朗日方程（括号内为零），则整个导数为零。因此 $F - y' F_{y'}$ 是常数。这个结果在分析力学中对应于能量守恒。

---

### 七、总结

P.1 节从零建立了变分法的核心框架：

1. **泛函**是以函数为自变量的映射，形如 $J[y] = \int_a^b F(x, y, y') dx$。

2. **变分** $\delta y$ 是函数本身的微小扰动。核心规则：$\delta(y') = (\delta y)'$。

3. **泛函变分**的表达式：
   $$
   \delta J = \int_a^b \left( F_y \delta y + F_{y'} \delta y' \right) dx.
   $$

4. **欧拉-拉格朗日方程**通过分部积分 + 变分法基本引理导出：
   $$
   \frac{\partial F}{\partial y} - \frac{d}{dx}\frac{\partial F}{\partial y'} = 0.
   $$

5. **边界条件**：固定端点给出标准方程；自由端点给出自然边界条件 $F_{y'} = 0$。

6. **高阶导数**推广：
   $$
   \sum_{k=0}^n (-1)^k \frac{d^k}{dx^k} \frac{\partial F}{\partial y^{(k)}} = 0.
   $$

变分法的本质，是把“从函数空间中选一条最优曲线”这个无穷维搜索问题，转化为一个可解的微分方程。这个转化，正是分析力学从牛顿范式走向拉格朗日范式的数学基础。从下一节开始，我们将把这一工具用于物理世界，建立最小作用量原理和拉格朗日力学的完整体系。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.2 变分法的推广与约束](https://zhuanlan.zhihu.com/p/2061810275808973374)

**课程目录：** [分析力学目录](https://zhuanlan.zhihu.com/p/2061652577167717816)
<!-- zhihu-nav-bottom:end -->
