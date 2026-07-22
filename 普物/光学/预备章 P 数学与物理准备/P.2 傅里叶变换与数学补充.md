---
created: 2026-07-20T17:22:30+08:00
modified: 2026-07-20T17:29:27+08:00
zhihu-title: P.2 傅里叶变换与数学补充
zhihu-topics:
  - 光学
zhihu-link: https://zhuanlan.zhihu.com/p/2062589994406089592
zhihu-toc: true
zhihu-created-at: 2026-07-20 17:29
zhihu-updated-at: 2026-07-21 14:37
---
<!-- zhihu-nav-top:start -->
**上一节：** [P.1 振动与波动回顾](https://zhuanlan.zhihu.com/p/2062589956233737211)
<!-- zhihu-nav-top:end -->

## P.2 傅里叶变换与数学补充

这一节的任务是建立处理光波场在空间上变化的核心数学工具——傅里叶分析。我们从周期函数的傅里叶级数出发，严格过渡到非周期函数的傅里叶积分，引入空间频率与角谱的概念。接着推导几个常用的傅里叶变换对：矩形函数、高斯函数和梳状函数。然后，将这些工具与标量衍射理论联系起来，阐明任意光场可以展开为不同方向传播的平面波的叠加（角谱展开）。最后，给出傍轴近似下菲涅耳近似和夫琅禾费近似的数学形式，为后续衍射章节奠定严格基础。

---

### 一、从傅里叶级数到傅里叶积分

#### 1.1 周期函数的傅里叶级数

设函数 $f(x)$ 在区间 $[-L/2, L/2]$ 上有定义，且满足狄利克雷条件，周期延拓为周期 $L$ 的函数。其傅里叶级数展开为

$$
f(x) = \sum_{n=-\infty}^{\infty} c_n \, e^{i 2\pi n x / L} \tag{1}
$$

其中傅里叶系数

$$
c_n = \frac{1}{L} \int_{-L/2}^{L/2} f(x) \, e^{-i 2\pi n x / L} \, dx \tag{2}
$$

这组公式表明：任何一个周期函数都可以看作一系列离散频率 $k_n = 2\pi n/L$ 的复指数函数的线性叠加。

如果采用角波数 $k_n = 2\pi n/L$，相邻频率的间隔为

$$
\Delta k = k_{n+1} - k_n = \frac{2\pi}{L}
$$

则级数可写为

$$
f(x) = \sum_{n=-\infty}^{\infty} c_n \, e^{i k_n x} , \quad
c_n = \frac{1}{L} \int_{-L/2}^{L/2} f(x) \, e^{-i k_n x} \, dx
$$

#### 1.2 向非周期函数的过渡：傅里叶积分

当周期 $L \to \infty$ 时，函数不再周期重复，变为非周期函数。此时，频率间隔 $\Delta k = 2\pi/L \to 0$，离散的 $k_n$ 趋近于连续的 $k$ 变量，求和过渡为积分。

重写傅里叶级数展开式：

$$
f(x) = \sum_{n=-\infty}^{\infty} \left[ \frac{1}{L} \int_{-L/2}^{L/2} f(\xi) \, e^{-i k_n \xi} \, d\xi \right] e^{i k_n x}
$$

注意到 $\frac{1}{L} = \frac{\Delta k}{2\pi}$，所以

$$
f(x) = \sum_{n=-\infty}^{\infty} \frac{\Delta k}{2\pi} \left[ \int_{-L/2}^{L/2} f(\xi) \, e^{-i k_n \xi} \, d\xi \right] e^{i k_n x}
$$

当 $L \to \infty$ 时，$\Delta k \to 0$，求和变成对 $k$ 的积分，积分区间变为 $(-\infty, \infty)$，方括号内的积分变为在整个实轴上的积分：

$$
\boxed{f(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \left[ \int_{-\infty}^{\infty} f(\xi) \, e^{-i k \xi} \, d\xi \right] e^{i k x} \, dk} \tag{3}
$$

定义函数 $f(x)$ 的**傅里叶变换**为

$$
\boxed{F(k) = \int_{-\infty}^{\infty} f(x) \, e^{-i k x} \, dx} \tag{4}
$$

则**傅里叶逆变换**为

$$
\boxed{f(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} F(k) \, e^{i k x} \, dk} \tag{5}
$$

这两个公式构成了傅里叶积分变换对。$F(k)$ 称为 $f(x)$ 的**频谱**或**角谱**（angular spectrum）。物理上，$F(k)$ 表示组成 $f(x)$ 的各种空间频率 $k$ 分量的复振幅。

#### 1.3 对称形式的傅里叶变换

在光学中，常采用空间频率 $\nu_x = k/(2\pi)$ 作为变量，这样正逆变换更具有对称性。定义

$$
\tilde{f}(\nu_x) = \int_{-\infty}^{\infty} f(x) \, e^{-i 2\pi \nu_x x} \, dx \tag{6}
$$

$$
f(x) = \int_{-\infty}^{\infty} \tilde{f}(\nu_x) \, e^{i 2\pi \nu_x x} \, d\nu_x \tag{7}
$$

这里 $\nu_x$ 的量纲是长度的倒数，表示单位长度内的周期数。若 $k$ 是角波数，则 $k = 2\pi \nu_x$。本文主要采用 (4)-(5) 的形式，有时为了方便也会使用 (6)-(7) 的形式，两者本质上等价，只需注意系数 $1/(2\pi)$ 的位置。

---

### 二、空间频率与角谱概念

#### 2.1 从时间信号到空间信号

在振动与波动回顾（P.1）中，我们用复数表示时间简谐信号 $e^{-i\omega t}$。对空间分布的光场，可以类似地引入**空间频率**的概念。

设 $x$ 轴上的光场复振幅分布为 $u(x)$（暂时忽略时间因子 $e^{-i\omega t}$，在单色光近似下可将其分离）。将 $u(x)$ 做傅里叶变换：

$$
U(k_x) = \int_{-\infty}^{\infty} u(x) \, e^{-i k_x x} \, dx \tag{8}
$$

逆变换为

$$
u(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} U(k_x) \, e^{i k_x x} \, dk_x \tag{9}
$$

这里 $k_x$ 是 $x$ 方向的空间角频率（角波数），物理上代表一个沿 $x$ 方向传播的平面波波矢的 $x$ 分量。一个单色平面波具有形式

$$
e^{i (k_x x + k_y y + k_z z - \omega t)}
$$

在 $z = 0$ 平面上，其横向部分为 $e^{i (k_x x + k_y y)}$。因此，(9) 式的物理含义是：任意一个横向光场分布 $u(x,y)$ 都可以分解为无穷多个不同空间频率 $(k_x, k_y)$ 的平面波的叠加，每个平面波的复振幅由 $U(k_x, k_y)$ 给出。这个 $U(k_x, k_y)$ 就称为该光场的**角谱**。

#### 2.2 二维傅里叶变换

推广到二维，设光场在 $z=0$ 平面的复振幅为 $u(x,y)$，其角谱定义为

$$
U(k_x, k_y) = \iint_{-\infty}^{\infty} u(x,y) \, e^{-i (k_x x + k_y y)} \, dx\,dy \tag{10}
$$

逆变换为

$$
u(x,y) = \frac{1}{(2\pi)^2} \iint_{-\infty}^{\infty} U(k_x, k_y) \, e^{i (k_x x + k_y y)} \, dk_x\,dk_y \tag{11}
$$

每个分量 $e^{i(k_x x + k_y y)}$ 代表一个在 $(k_x, k_y)$ 方向传播的平面波在 $z=0$ 面上的相位分布。需注意，并非任意 $(k_x, k_y)$ 都对应一个可传播的平面波，还必须满足波矢大小等于介质中的波数 $k = n \omega/c$ 这一约束，即

$$
k_x^2 + k_y^2 + k_z^2 = k^2 \tag{12}
$$

若 $k_x^2 + k_y^2 \le k^2$，则 $k_z = \sqrt{k^2 - k_x^2 - k_y^2}$ 为实数，对应**传播波**；若 $k_x^2 + k_y^2 > k^2$，则 $k_z = i \sqrt{k_x^2 + k_y^2 - k^2}$ 为纯虚数，对应沿 $z$ 指数衰减的**隐失波**（倏逝波）。

---

### 三、几个基本函数的傅里叶变换对

为了后续处理衍射问题，我们需要熟练掌握几个常用函数的傅里叶变换。这部分推导将非常细致，每个积分都完全展开。

#### 3.1 矩形函数

矩形函数定义为

$$
\text{rect}(x) = \begin{cases}
1, & |x| < 1/2 \\
1/2, & |x| = 1/2 \\
0, & |x| > 1/2
\end{cases}
$$

（积分时边界值不影响结果，可忽略 $x=\pm 1/2$ 的点。）宽度为 $a$ 的矩形函数为 $\text{rect}(x/a)$。

求其傅里叶变换：

$$
F(k) = \int_{-\infty}^{\infty} \text{rect}\!\left(\frac{x}{a}\right) e^{-i k x} \, dx
= \int_{-a/2}^{a/2} e^{-i k x} \, dx
$$

计算积分：

$$
\begin{aligned}
F(k) &= \left. \frac{e^{-i k x}}{-i k} \right|_{-a/2}^{a/2} \\
&= \frac{1}{-i k} \left( e^{-i k a/2} - e^{i k a/2} \right) \\
&= \frac{2}{k} \cdot \frac{e^{i k a/2} - e^{-i k a/2}}{2i} \\
&= \frac{2}{k} \sin\!\left( \frac{k a}{2} \right) \\
&= a \frac{\sin(k a/2)}{k a/2}
\end{aligned}
$$

定义 $\text{sinc}$ 函数

$$
\text{sinc}(\xi) \equiv \frac{\sin(\pi \xi)}{\pi \xi}
$$

则 $k a/2 = \pi (k a / (2\pi)) = \pi \nu_x a$，所以

$$
F(k) = a \, \text{sinc}\!\left( \frac{k a}{2\pi} \right)
$$

或者用角频率表示：

$$
\boxed{\int_{-\infty}^{\infty} \text{rect}\!\left(\frac{x}{a}\right) e^{-i k x} dx = a \, \text{sinc}\!\left( \frac{k a}{2\pi} \right)} \tag{13}
$$

有时也定义 $\text{sinc}(u) = \sin u / u$，则结果为 $a \, \text{sinc}(k a/2)$。在光学中，缝宽为 $a$ 的单缝夫琅禾费衍射的振幅分布正比于 $\text{sinc}$ 函数，其根源就在这里。

#### 3.2 高斯函数

高斯函数定义为

$$
f(x) = e^{-a x^2}, \quad a > 0
$$

求傅里叶变换：

$$
F(k) = \int_{-\infty}^{\infty} e^{-a x^2} e^{-i k x} dx = \int_{-\infty}^{\infty} e^{-a x^2 - i k x} dx
$$

采用配方法：

$$
a x^2 + i k x = a \left( x^2 + \frac{i k}{a} x \right) = a \left[ \left( x + \frac{i k}{2a} \right)^2 + \frac{k^2}{4a^2} \right]
$$

因此

$$
\begin{aligned}
F(k) &= \int_{-\infty}^{\infty} e^{-a \left( x + \frac{i k}{2a} \right)^2} e^{-\frac{k^2}{4a}} dx \\
&= e^{-\frac{k^2}{4a}} \int_{-\infty}^{\infty} e^{-a u^2} du \quad \left( u = x + \frac{i k}{2a} \right)
\end{aligned}
$$

这个积分是经典的高斯积分，通过将积分路径从实轴平移到复平面上的直线 $\text{Im}(u) = k/(2a)$，利用柯西积分定理可知其值等于实轴上的高斯积分：

$$
\int_{-\infty}^{\infty} e^{-a u^2} du = \sqrt{\frac{\pi}{a}}
$$

于是

$$
F(k) = \sqrt{\frac{\pi}{a}} \, e^{-\frac{k^2}{4a}} \tag{14}
$$

这是一个极为优美的结果：高斯函数的傅里叶变换仍然是高斯函数。若空间分布是宽度为 $\sigma$ 的高斯型，即 $f(x) = e^{-x^2/(2\sigma^2)}$，则对应 $a = 1/(2\sigma^2)$，代入 (14) 得

$$
F(k) = \sqrt{2\pi} \, \sigma \, e^{-\sigma^2 k^2 / 2} \tag{15}
$$

或写成

$$
\boxed{\mathcal{F}\!\left[ e^{-x^2/(2\sigma^2)} \right] = \sqrt{2\pi}\,\sigma\, e^{-\sigma^2 k^2/2}}
$$

在光学中，激光束的基模（TEM₀₀）横截面光场就是高斯分布，其角谱也是高斯型，这一性质在激光传输和光学谐振腔理论中至关重要。

#### 3.3 梳状函数（狄拉克梳）

梳状函数 $\text{comb}(x)$ 定义为等间距分布的 δ 函数序列：

$$
\text{comb}(x) = \sum_{n=-\infty}^{\infty} \delta(x - n) \tag{16}
$$

间距为 $a$ 的梳状函数可写为 $\frac{1}{a} \text{comb}(x/a) = \sum_{n} \delta(x - n a)$。

求梳状函数的傅里叶变换：

$$
F(k) = \int_{-\infty}^{\infty} \sum_{n=-\infty}^{\infty} \delta(x - n) \, e^{-i k x} dx = \sum_{n=-\infty}^{\infty} e^{-i k n}
$$

这个无穷级数并不在普通函数意义下收敛。我们需要用分布理论（广义函数）处理。利用泊松求和公式：

$$
\sum_{n=-\infty}^{\infty} \delta(x - n) = \sum_{m=-\infty}^{\infty} e^{i 2\pi m x} \quad \text{（在分布意义下）} \tag{17}
$$

那么其傅里叶变换为

$$
\begin{aligned}
F(k) &= \int_{-\infty}^{\infty} \left[ \sum_{m=-\infty}^{\infty} e^{i 2\pi m x} \right] e^{-i k x} dx \\
&= \sum_{m=-\infty}^{\infty} \int_{-\infty}^{\infty} e^{-i (k - 2\pi m) x} dx \\
&= \sum_{m=-\infty}^{\infty} 2\pi \, \delta(k - 2\pi m)
\end{aligned}
$$

因为 $\int_{-\infty}^{\infty} e^{-i (k - k_0) x} dx = 2\pi \delta(k - k_0)$。所以

$$
\boxed{\mathcal{F}\!\left[ \sum_{n=-\infty}^{\infty} \delta(x - n) \right] = 2\pi \sum_{m=-\infty}^{\infty} \delta(k - 2\pi m)} \tag{18}
$$

对于间距为 $a$ 的梳状函数 $\text{III}_a(x) = \sum_n \delta(x - n a)$，通过变量代换 $x' = x/a$ 可得

$$
\mathcal{F}\!\left[ \sum_{n=-\infty}^{\infty} \delta(x - n a) \right] = \frac{2\pi}{a} \sum_{m=-\infty}^{\infty} \delta\!\left(k - \frac{2\pi m}{a}\right) \tag{19}
$$

这个结果表明：空间域的梳状函数变换到频率域仍然是梳状函数，只是间距变为 $2\pi/a$。这一性质是衍射光栅产生多级离散光谱的数学基础。

---

### 四、标量衍射理论中的平面波展开思想

#### 4.1 亥姆霍兹方程的平面波解

在单色光场中，电场或磁场的一个标量分量 $u(\mathbf{r}, t) = U(\mathbf{r}) e^{-i\omega t}$ 满足亥姆霍兹方程

$$
\nabla^2 U + k^2 U = 0, \quad k = \frac{\omega}{c} = \frac{2\pi}{\lambda} \tag{20}
$$

最简解是平面波：

$$
U(\mathbf{r}) = A e^{i \mathbf{k} \cdot \mathbf{r}} = A e^{i(k_x x + k_y y + k_z z)}
$$

其中波矢 $\mathbf{k} = (k_x, k_y, k_z)$ 满足色散关系

$$
k_x^2 + k_y^2 + k_z^2 = k^2 \tag{21}
$$

由于亥姆霍兹方程是线性的，任意多个平面波的叠加也是解：

$$
U(\mathbf{r}) = \iiint A(\mathbf{k}) \, e^{i \mathbf{k} \cdot \mathbf{r}} \, \delta(k_x^2 + k_y^2 + k_z^2 - k^2) \, dk_x dk_y dk_z
$$

但更常见的做法是，利用 (21) 消去 $k_z$，将解限制在沿 $z$ 正向传播的分量上：

$$
k_z = \sqrt{k^2 - k_x^2 - k_y^2}, \quad (\text{取 } k_z > 0 \text{ 或 } \text{Im}(k_z) > 0)
$$

于是，任意一个在 $z$ 半空间传播的标量波场可以写为

$$
U(x,y,z) = \iint_{-\infty}^{\infty} A(k_x, k_y) \, e^{i(k_x x + k_y y + k_z z)} \, dk_x dk_y \tag{22}
$$

这就是**角谱展开**：平面波的复振幅 $A(k_x, k_y)$ 是 $z=0$ 平面上的场 $U_0(x,y) = U(x,y,0)$ 的二维傅里叶变换。

#### 4.2 传播作为角谱的线性滤波器

在 $z=0$ 平面，由逆傅里叶变换：

$$
U_0(x,y) = \iint_{-\infty}^{\infty} A(k_x, k_y) \, e^{i(k_x x + k_y y)} \, dk_x dk_y
$$

因此

$$
A(k_x, k_y) = \frac{1}{(2\pi)^2} \iint_{-\infty}^{\infty} U_0(x,y) \, e^{-i(k_x x + k_y y)} \, dx dy \tag{23}
$$

将 $A(k_x, k_y)$ 代入 (22) 式，得到传播后的场：

$$
U(x,y,z) = \iint_{-\infty}^{\infty} A(k_x, k_y) \, e^{i k_z z} \, e^{i(k_x x + k_y y)} \, dk_x dk_y \tag{24}
$$

这表明：在角谱域，光场的传播等价于每个平面波分量乘上一个相位因子 $e^{i k_z z}$。这个因子是一个传递函数：

$$
H(k_x, k_y; z) = e^{i k_z z} = e^{i z \sqrt{k^2 - k_x^2 - k_y^2}} \tag{25}
$$

对于 $k_x^2 + k_y^2 \le k^2$ 的分量，$k_z$ 为实数，相位因子为纯相位移动；对于 $k_x^2 + k_y^2 > k^2$ 的分量，$k_z$ 为正虚数，$e^{i k_z z} = e^{- z \sqrt{k_x^2 + k_y^2 - k^2}}$ 随 $z$ 指数衰减，即为隐失波，在远离衍射屏几个波长后就消失殆尽，不参与远场衍射。

---

### 五、傍轴近似与菲涅耳近似、夫琅禾费近似

实际光学系统大多满足**傍轴条件**，即光波的传播方向与光轴（$z$ 轴）夹角很小，这意味着角谱 $A(k_x, k_y)$ 只在 $k_x, k_y \ll k$ 的区域有显著非零值。

#### 5.1 傍轴近似下 $k_z$ 的泰勒展开

傍轴条件：

$$
k_x^2 + k_y^2 \ll k^2
$$

对平方根做泰勒展开：

$$
k_z = k \sqrt{1 - \frac{k_x^2 + k_y^2}{k^2}} = k \left[ 1 - \frac{1}{2} \frac{k_x^2 + k_y^2}{k^2} - \frac{1}{8} \left(\frac{k_x^2 + k_y^2}{k^2}\right)^2 - \cdots \right]
$$

保留到一次项（即 $(k_x^2 + k_y^2)/k^2$ 的一阶项），称为**菲涅耳近似**：

$$
k_z \approx k - \frac{k_x^2 + k_y^2}{2k} \tag{26}
$$

此时传递函数 (25) 变为

$$
H(k_x, k_y; z) \approx e^{i k z} \, \exp\!\left[ -i \frac{z}{2k} (k_x^2 + k_y^2) \right] \tag{27}
$$

若只保留到零阶项，即 $k_z \approx k$，则称为**夫琅禾费近似**，但要注意这个近似需要在什么条件下才能替代菲涅耳近似，实际上夫琅禾费近似对场本身的表达式有更强的要求，我们将马上说明。

#### 5.2 菲涅耳衍射积分

将菲涅耳近似 (26) 代入角谱传播公式 (24)：

$$
\begin{aligned}
U(x,y,z) &\approx e^{i k z} \iint_{-\infty}^{\infty} A(k_x, k_y) \, e^{i(k_x x + k_y y)} \, \exp\!\left[ -i \frac{z}{2k} (k_x^2 + k_y^2) \right] \, dk_x dk_y
\end{aligned}
$$

根据卷积定理，角谱域乘以高斯二次相位因子相当于空间域与一个二次相位函数的卷积。具体地，代入 $A(k_x, k_y)$ 的表达式 (23)：

$$
\begin{aligned}
U(x,y,z) &= e^{i k z} \iint_{-\infty}^{\infty} \left[ \frac{1}{(2\pi)^2} \iint_{-\infty}^{\infty} U_0(\xi,\eta) e^{-i(k_x \xi + k_y \eta)} d\xi d\eta \right] \\
&\qquad \times e^{i(k_x x + k_y y)} e^{-i \frac{z}{2k}(k_x^2 + k_y^2)} dk_x dk_y
\end{aligned}
$$

交换积分次序，整理得

$$
U(x,y,z) = \frac{e^{i k z}}{(2\pi)^2} \iint_{-\infty}^{\infty} U_0(\xi,\eta) \left[ \iint_{-\infty}^{\infty} e^{i[k_x (x-\xi) + k_y (y-\eta)]} e^{-i \frac{z}{2k}(k_x^2 + k_y^2)} dk_x dk_y \right] d\xi d\eta
$$

方括号内的二维积分可分离为两个一维积分的乘积，每个都是高斯型傅里叶逆变换核。计算标准积分：

$$
\int_{-\infty}^{\infty} e^{i p u} e^{-i \alpha u^2} du = \sqrt{\frac{\pi}{i\alpha}} \, e^{-i \frac{p^2}{4\alpha}} \quad (\alpha > 0)
$$

这里 $\alpha = z/(2k)$，对 $k_x$ 积分：

$$
\int_{-\infty}^{\infty} e^{i k_x (x-\xi)} e^{-i \frac{z}{2k} k_x^2} dk_x = \sqrt{\frac{2\pi k}{i z}} \, e^{i \frac{k}{2z} (x-\xi)^2}
$$

同理对 $k_y$ 积分。相乘后，得到

$$
U(x,y,z) = \frac{e^{i k z}}{i \lambda z} \iint_{-\infty}^{\infty} U_0(\xi,\eta) \, \exp\!\left\{ i \frac{k}{2z} \left[ (x-\xi)^2 + (y-\eta)^2 \right] \right\} d\xi d\eta \tag{28}
$$

这个公式就是**菲涅耳衍射积分**。其中 $\lambda = 2\pi/k$ 是波长。因子 $1/(i\lambda z)$ 保证了能量守恒。

菲涅耳衍射积分成立的条件（菲涅耳近似条件）是：在展开 $k_z$ 时忽略的二阶及以上项引起的相位误差远小于 $\pi$。可以证明，它要求

$$
z^3 \gg \frac{\pi}{4\lambda} \left[ (x-\xi)^2 + (y-\eta)^2 \right]_{\text{max}}^2 \tag{29}
$$

这在旁轴且观察距离不太远时容易满足。

#### 5.3 夫琅禾费衍射积分

如果观测距离 $z$ 进一步增大，满足更强的条件

$$
z \gg \frac{k (\xi^2 + \eta^2)_{\text{max}}}{2} \tag{30}
$$

那么菲涅耳积分指数中的二次相位因子

$$
\frac{k}{2z} (\xi^2 + \eta^2)
$$

可以忽略。将菲涅耳积分 (28) 中的平方项展开：

$$
(x-\xi)^2 + (y-\eta)^2 = x^2 + y^2 - 2(x\xi + y\eta) + \xi^2 + \eta^2
$$

于是

$$
U(x,y,z) = \frac{e^{i k z}}{i \lambda z} e^{i \frac{k}{2z}(x^2 + y^2)} \iint_{-\infty}^{\infty} \left[ U_0(\xi,\eta) \, e^{i \frac{k}{2z}(\xi^2 + \eta^2)} \right] e^{-i \frac{k}{z}(x\xi + y\eta)} d\xi d\eta
$$

当 $\frac{k}{2z}(\xi^2 + \eta^2) \ll \pi$，即 $z \gg \frac{\pi}{\lambda} (\xi^2 + \eta^2)_{\text{max}}$ 时，括号中的二次相位因子近似为 $1$。忽略此项，得到

$$
\boxed{U(x,y,z) = \frac{e^{i k z}}{i \lambda z} e^{i \frac{k}{2z}(x^2 + y^2)} \iint_{-\infty}^{\infty} U_0(\xi,\eta) \, e^{-i \frac{2\pi}{\lambda z}(x\xi + y\eta)} d\xi d\eta} \tag{31}
$$

这就是**夫琅禾费衍射公式**。关键在于，除去积分前的相位因子，观察面上的场正比于孔径平面场 $U_0(\xi,\eta)$ 的二维傅里叶变换，频率变量为

$$
k_x' = \frac{k x}{z}, \quad k_y' = \frac{k y}{z}
\quad \text{或空间频率} \quad \nu_x = \frac{x}{\lambda z},\; \nu_y = \frac{y}{\lambda z}
$$

夫琅禾费衍射因此又称为**远场衍射**，它直接将衍射屏的复振幅分布与远场角谱对应起来。透镜可以在焦平面上实现夫琅禾费衍射，而无需极远的距离，这是傅里叶光学的基础。

---

通过以上推导，我们从傅里叶级数出发，严格建立了一维和二维傅里叶变换的数学框架，引入了空间频率和角谱，推导了矩形、高斯、梳状函数的变换对，并将这些工具无缝衔接到标量衍射理论中，详细给出了菲涅耳近似和夫琅禾费近似的数学来源及相应的衍射积分公式。这套数学语言是后续分析光的干涉、衍射、成像和信息处理的基石。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [1.1 费马原理](https://zhuanlan.zhihu.com/p/2062590210819602399)

**课程目录：** [光学目录](https://zhuanlan.zhihu.com/p/2062614039084192499)
<!-- zhihu-nav-bottom:end -->
