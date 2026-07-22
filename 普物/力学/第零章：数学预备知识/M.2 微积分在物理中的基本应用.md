---
created: 2026-07-20T06:03:34+08:00
modified: 2026-07-20T06:31:58+08:00
zhihu-title: M.2 微积分在物理中的基本应用
zhihu-topics:
  - 力学
zhihu-link: https://zhuanlan.zhihu.com/p/2062423053758363601
zhihu-toc: true
zhihu-created-at: 2026-07-20 10:40
---
<!-- zhihu-nav-top:start -->
**上一节：** [M.1 向量代数](https://zhuanlan.zhihu.com/p/2062422859478196964)
<!-- zhihu-nav-top:end -->

## M.2 微积分在物理中的基本应用

微积分是物理学的语言。没有微积分，物理学只能描述匀速直线运动——世界是静止的、均匀的、平庸的。有了微积分，我们可以描述变化：速度如何随时间改变，力如何在空间中累积，物理量如何在无限小的尺度上定义。这一节我们从物理学的角度重新审视微积分：导数=变化率，积分=累积量，微分方程=物理定律的数学表达。

我们的目标不是重复数学课上的推导，而是建立“看到物理量就想到微积分操作”的本能反应。

---

### 一、导数——变化率

在物理学中，导数无处不在。它的物理意义永远是：**一个物理量随另一个物理量变化的快慢**。

**速度**：位置随时间的变化率。

$$
\mathbf{v}(t) = \frac{d\mathbf{r}}{dt} = \lim_{\Delta t \to 0} \frac{\mathbf{r}(t+\Delta t) - \mathbf{r}(t)}{\Delta t}. \tag{1}
$$

速度是矢量，其方向沿轨迹切线。在直角坐标系中，速度的分量为：

$$
v_x = \frac{dx}{dt}, \quad v_y = \frac{dy}{dt}, \quad v_z = \frac{dz}{dt}. \tag{2}
$$

**加速度**：速度随时间的变化率。

$$
\mathbf{a}(t) = \frac{d\mathbf{v}}{dt} = \frac{d^2\mathbf{r}}{dt^2}. \tag{3}
$$

加速度是矢量。在直角坐标系中：

$$
a_x = \frac{dv_x}{dt} = \frac{d^2x}{dt^2}, \quad a_y = \frac{dv_y}{dt} = \frac{d^2y}{dt^2}, \quad a_z = \frac{dv_z}{dt} = \frac{d^2z}{dt^2}. \tag{4}
$$

**一个关键的变量代换技巧**：

在力学中，我们经常遇到加速度是位置的函数，而不是时间的函数（例如弹簧力 $F = -kx$，加速度 $a = -(k/m)x$，依赖于位置 $x$）。这时我们想把运动方程中的变量从时间 $t$ 变换到位置 $x$。

关键恒等式是：

$$
a = \frac{dv}{dt} = \frac{dv}{dx}\frac{dx}{dt} = v\frac{dv}{dx}. \tag{5}
$$

推导：链式法则的直接应用。第一层，加速度是速度对时间的导数；第二层，速度 $v$ 可以看作位置 $x$ 的函数（通过运动轨迹 $v(x)$），然后用链式法则展开 $\frac{dv}{dt} = \frac{dv}{dx}\frac{dx}{dt}$；第三层，$\frac{dx}{dt} = v$。于是得到 $a = v\frac{dv}{dx}$。

这个技巧非常有用。举例：已知加速度 $a(x) = -\omega^2 x$（简谐运动），我们可以直接积分得到速度与位置的关系，而不需要先求 $x(t)$：

$$
v\frac{dv}{dx} = -\omega^2 x \quad \Rightarrow \quad v\,dv = -\omega^2 x\,dx.
$$

两边积分：

$$
\int v\,dv = \int -\omega^2 x\,dx \quad \Rightarrow \quad \frac{1}{2}v^2 = -\frac{1}{2}\omega^2 x^2 + C.
$$

这就是能量守恒的雏形。

---

### 二、微分——微元的思想

微分在物理中的核心作用是：**把一个整体切割成无限小的部分，在每个小部分上物理量近似均匀，然后求和（积分）得到整体**。

**一维细杆的质量**：

设细杆的线密度为 $\lambda(x)$（单位长度的质量），在位置 $x$ 处取一小段 $dx$，这一小段的质量为 $dm = \lambda(x)dx$。整个细杆的质量：

$$
M = \int dm = \int_{x_1}^{x_2} \lambda(x)\,dx. \tag{6}
$$

如果密度均匀 $\lambda = \text{常数}$，则 $M = \lambda L$，其中 $L = x_2 - x_1$。

**薄圆盘的质量**：

设面密度为 $\sigma(r)$（单位面积的质量），在半径 $r$ 处取一个宽度为 $dr$ 的圆环，圆环面积为 $dA = 2\pi r\,dr$，质量为 $dm = \sigma(r)\,dA = \sigma(r) 2\pi r\,dr$。整个圆盘的质量：

$$
M = \int dm = \int_0^R \sigma(r) 2\pi r\,dr. \tag{7}
$$

如果面密度均匀 $\sigma = \text{常数}$，则 $M = \sigma \pi R^2$。

**三维物体的质量**：

设体密度为 $\rho(\mathbf{r})$（单位体积的质量），取体积微元 $dV$，质量 $dm = \rho(\mathbf{r})\,dV$。总质量：

$$
M = \int_V dm = \int_V \rho(\mathbf{r})\,dV. \tag{8}
$$

在直角坐标系中，$dV = dx\,dy\,dz$。在球坐标系中，$dV = r^2\sin\theta\,dr\,d\theta\,d\phi$。

**电量的微元取法**：

完全类似。若电荷线密度为 $\lambda_e(x)$，则 $dq = \lambda_e(x)dx$；若面密度为 $\sigma_e(x,y)$，则 $dq = \sigma_e dA$；若体密度为 $\rho_e(\mathbf{r})$，则 $dq = \rho_e dV$。

**微元法的核心原则**：在微元尺度上，物理量可以视为常数。这就是微积分处理物理问题的本质——将非均匀问题转化为均匀问题的积分。

---

### 三、积分——累积量

物理学中的许多重要物理量都是某种量的积分（累积）：
- 位移 = 速度对时间的积分
- 功 = 力对位移的积分
- 冲量 = 力对时间的积分
- 质量 = 密度对体积的积分
- 电荷 = 电荷密度对体积的积分

**由加速度求速度、由速度求位置**：

已知加速度 $a(t)$ 和初始条件 $v(0) = v_0$，速度：

$$
v(t) = v_0 + \int_0^t a(t')\,dt'. \tag{9}
$$

已知速度 $v(t)$ 和初始位置 $x(0) = x_0$，位置：

$$
x(t) = x_0 + \int_0^t v(t')\,dt'. \tag{10}
$$

积分是导数的逆运算。给定加速度 $a(t)$，通过两次积分得到位置 $x(t)$，这完整确定了运动。

**变力做功**：

力 $\mathbf{F}(\mathbf{r})$ 沿路径 $C$ 从 $\mathbf{r}_1$ 到 $\mathbf{r}_2$ 做的功为：

$$
W = \int_C \mathbf{F}\cdot d\mathbf{r}. \tag{11}
$$

这个积分是**线积分**（路径积分），它不只是一维积分的简单推广。沿路径 $C$ 的参数化：$\mathbf{r}(s)$，其中 $s$ 是路径参数（通常是弧长），$d\mathbf{r} = \frac{d\mathbf{r}}{ds}ds$，则

$$
W = \int_{s_1}^{s_2} \mathbf{F}(\mathbf{r}(s))\cdot \frac{d\mathbf{r}}{ds}\,ds. \tag{12}
$$

**功的几何意义**：在一维情况下，$W = \int_{x_1}^{x_2} F(x)\,dx$，这就是 $F$-$x$ 曲线下方的面积。在三维中，功是力场沿路径的“累积”——力在路径每个点上的切向分量乘以路径长度，再全部加起来。

**冲量**：

力 $\mathbf{F}(t)$ 在时间间隔 $[t_1, t_2]$ 内的冲量为：

$$
\mathbf{I} = \int_{t_1}^{t_2} \mathbf{F}(t)\,dt. \tag{13}
$$

冲量是力在时间上的累积。冲量定理 $\mathbf{I} = \Delta\mathbf{p}$ 是动量定理的积分形式。

---

### 四、微分方程初步

物理学中，绝大多数定律都是微分方程——它们描述物理量之间的关系，而这个关系涉及导数。

**一阶线性微分方程——积分因子法**：

考虑方程：

$$
\frac{dx}{dt} + a(t)x = b(t). \tag{14}
$$

这是线性方程，因为 $x$ 和 $\frac{dx}{dt}$ 都是一次幂。

求解方法：找积分因子 $\mu(t)$，使左边的式子变为全导数。

推导：我们希望

$$
\frac{d}{dt}[\mu(t)x] = \mu(t)\frac{dx}{dt} + \frac{d\mu}{dt}x.
$$

由原方程 $\frac{dx}{dt} = b - a x$，代入：

$$
\frac{d}{dt}[\mu x] = \mu(b - ax) + \frac{d\mu}{dt}x = \mu b + \left(\frac{d\mu}{dt} - \mu a\right)x.
$$

为了使右边只含 $x$ 的项消失，我们令：

$$
\frac{d\mu}{dt} = \mu a(t). \tag{15}
$$

解这个方程（分离变量）：

$$
\frac{d\mu}{\mu} = a(t)dt \quad \Rightarrow \quad \ln\mu = \int a(t)dt \quad \Rightarrow \quad \mu(t) = e^{\int a(t)dt}. \tag{16}
$$

于是：

$$
\frac{d}{dt}[\mu x] = \mu b(t).
$$

两边积分：

$$
\mu(t)x(t) = \int \mu(t)b(t)dt + C.
$$

所以：

$$
x(t) = \frac{1}{\mu(t)}\left[\int \mu(t)b(t)dt + C\right]. \tag{17}
$$

这就是一阶线性微分方程的通解。

**举例**：阻尼运动 $m\frac{dv}{dt} = -kv$，即 $\frac{dv}{dt} + \frac{k}{m}v = 0$。这里 $a(t) = k/m$（常数），$b(t) = 0$。积分因子 $\mu(t) = e^{(k/m)t}$。通解为：

$$
v(t) = v_0 e^{-(k/m)t}. \tag{18}
$$

**分离变量法**：

当方程可以写成 $g(x)\frac{dx}{dt} = h(t)$ 的形式时，可以直接分离变量：

$$
g(x)\,dx = h(t)\,dt. \tag{19}
$$

两边积分：

$$
\int g(x)\,dx = \int h(t)\,dt + C. \tag{20}
$$

这种方法适用于大量一阶微分方程。例如自由落体 $m\frac{dv}{dt} = mg$，即 $\frac{dv}{dt} = g$，分离变量 $dv = g\,dt$，积分得 $v = gt + C$。

**微分方程的物理意义**：微分方程不是静止的公式，它描述的是“变化如何依赖于当前状态”。给定初始状态，解微分方程就是预测未来。

---

### 五、量纲与单位制

**量纲**是物理量的基本属性，不依赖单位的选择。力学中三个基本量纲是：

- 长度：$L$
- 质量：$M$
- 时间：$T$

其他力学量的量纲可由这三个导出：

| 物理量 | 量纲 |
|--------|------|
| 速度 | $L T^{-1}$ |
| 加速度 | $L T^{-2}$ |
| 力 | $M L T^{-2}$ |
| 功/能量 | $M L^2 T^{-2}$ |
| 功率 | $M L^2 T^{-3}$ |
| 动量 | $M L T^{-1}$ |
| 角动量 | $M L^2 T^{-1}$ |
| 压强 | $M L^{-1} T^{-2}$ |

**量纲齐次性原理**：任何物理方程的两端必须具有相同的量纲。这个简单原则是检验公式正确性的第一道防线。

**举例**：单摆的周期 $T$ 可能依赖于摆长 $l$、重力加速度 $g$ 和质量 $m$。设 $T = k l^\alpha g^\beta m^\gamma$，其中 $k$ 是无量纲常数。量纲方程：

$$
T = L^\alpha (L T^{-2})^\beta M^\gamma = L^{\alpha+\beta} T^{-2\beta} M^\gamma.
$$

比较两边量纲：$T^1 = T^{-2\beta} \Rightarrow \beta = -1/2$；$L^0 = L^{\alpha+\beta} \Rightarrow \alpha = 1/2$；$M^0 = M^\gamma \Rightarrow \gamma = 0$。因此：

$$
T = k \sqrt{\frac{l}{g}}. \tag{21}
$$

这就是单摆周期公式的定性形式。精确计算给出 $k = 2\pi$。量纲分析告诉我们：质量 $m$ 不出现——所有单摆在小角度下周期与质量无关。

**单位制**：国际单位制（SI）规定长度、质量、时间的基本单位分别是米（m）、千克（kg）、秒（s）。由这些导出力的单位——牛顿：$1\,\text{N} = 1\,\text{kg}\cdot\text{m}\cdot\text{s}^{-2}$。

---

### 六、总结

把 M.2 节的全部内容压缩成几条线：

1. **导数 = 变化率**：速度 $v = dx/dt$，加速度 $a = dv/dt$。关键变量代换技巧 $a = v\,dv/dx$，用于将 $a(x)$ 转化为 $v(x)$ 的积分。

2. **微分 = 微元**：把整体切成无限小的小块，每个小块视为常数，然后积分。$dm = \lambda dx$，$dm = \sigma dA$，$dm = \rho dV$。

3. **积分 = 累积量**：速度积分得位移，加速度积分得速度，力对位移积分得功，力对时间积分得冲量，密度对体积积分得质量。

4. **变力做功**：$W = \int \mathbf{F}\cdot d\mathbf{r}$ 是一维积分 $W = \int F(x)dx$ 的推广，是力沿路径的累积。

5. **一阶线性微分方程**：$\frac{dx}{dt} + a(t)x = b(t)$，积分因子法求解，$\mu(t) = e^{\int a(t)dt}$。

6. **分离变量法**：$g(x)dx = h(t)dt$，两边积分得隐式解。

7. **量纲分析**：物理方程两端量纲必须相同，可用于猜测公式形式和检验计算结果。

微积分是物理学的语言。导数告诉我们“变化有多快”，积分告诉我们“总共累积了多少”，微分方程告诉我们“变化如何依赖于状态本身”。这三者共同构成了从牛顿力学到量子场论的一切定量物理的基础。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [M.3 坐标系](https://zhuanlan.zhihu.com/p/2062423889158861169)

**课程目录：** [力学目录](https://zhuanlan.zhihu.com/p/2062487053170906945)
<!-- zhihu-nav-bottom:end -->
