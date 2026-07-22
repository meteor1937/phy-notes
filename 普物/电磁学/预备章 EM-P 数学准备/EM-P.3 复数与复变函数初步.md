---
created: 2026-07-20T19:25:52+08:00
modified: 2026-07-20T20:23:01+08:00
zhihu-title: EM-P.3 复数与复变函数初步
zhihu-topics:
  - 电磁学
zhihu-link: https://zhuanlan.zhihu.com/p/2062629304270431145
zhihu-toc: true
zhihu-created-at: 2026-07-20 21:32
zhihu-updated-at: 2026-07-21 14:41
---
<!-- zhihu-nav-top:start -->
**上一节：** [EM-P.2 狄拉克 δ 函数与格林函数初步](https://zhuanlan.zhihu.com/p/2062629264005009857)
<!-- zhihu-nav-top:end -->

## EM-P.3 复数与复变函数初步

电磁学从静电场、恒定电流、静磁场进入时变场后，所面对的方程在时间导数下由实函数变为复指数形式，求解大大简化。复数及其运算成为处理交流电路和时谐电磁波不可或缺的数学语言。正弦量的加减、微积分运算通过复数表示转化为代数运算，阻抗和导纳的概念也因此自然产生。这一节系统回顾复数的三种表示形式及其相互转换，推导并深刻理解欧拉公式，建立时谐量与其相量表示之间的一一对应关系，阐明相量法如何将微积分方程化为代数方程，最后引入阻抗与导纳的初步概念，为后续交流电路分析和电磁波章节奠定计算基础。

---

### 一、复数的三种表示形式及其转换

#### 1.1 代数形式与复平面

复数 $z$ 的代数形式为

$$
z = x + i y
$$

其中 $x, y$ 为实数，$i$ 为虚数单位，满足

$$
i^2 = -1
$$

$x = \operatorname{Re}(z)$ 称为复数的**实部**，$y = \operatorname{Im}(z)$ 称为复数的**虚部**（注意虚部是实数 $y$，而非 $i y$）。

在直角坐标平面上，以横轴为实轴、纵轴为虚轴建立**复平面**（Argand 图）。复数 $z = x + i y$ 对应于点 $(x, y)$，从原点指向该点的矢量长度为复数的**模**：

$$
|z| = \sqrt{x^2 + y^2}
$$

该矢量与正实轴的夹角（以弧度计，逆时针为正）称为复数的**辐角**，记为 $\arg z$，满足

$$
\tan(\arg z) = \frac{y}{x}, \quad x \neq 0
$$

辐角不是唯一的：$\arg z + 2\pi n$（$n$ 为任意整数）对应同一个复数。通常取 $(-\pi, \pi]$ 内的值作为**辐角主值**，记为 $\operatorname{Arg} z$。

#### 1.2 三角形式

由极坐标关系 $x = r \cos\theta$, $y = r \sin\theta$（其中 $r = |z|$, $\theta = \arg z$），代数形式转化为三角形式：

$$
z = r (\cos\theta + i \sin\theta) \tag{P.3.1}
$$

三角形式清晰地分离了模与辐角，便于乘除运算。

#### 1.3 指数形式与欧拉公式

欧拉公式建立了复指数与三角函数之间的深刻联系：

$$
\boxed{e^{i\theta} = \cos\theta + i \sin\theta} \tag{P.3.2}
$$

证明：将指数函数 $e^{i\theta}$ 在 $\theta = 0$ 附近作泰勒展开：

$$
e^{i\theta} = \sum_{n=0}^\infty \frac{(i\theta)^n}{n!} = \sum_{n=0}^\infty \frac{i^n \theta^n}{n!}
$$

将实部和虚部分开。利用 $i^0 = 1$, $i^1 = i$, $i^2 = -1$, $i^3 = -i$, $i^4 = 1$，周期为 4：

$$
\begin{aligned}
e^{i\theta} &= \left( 1 - \frac{\theta^2}{2!} + \frac{\theta^4}{4!} - \frac{\theta^6}{6!} + \cdots \right) + i \left( \theta - \frac{\theta^3}{3!} + \frac{\theta^5}{5!} - \frac{\theta^7}{7!} + \cdots \right) \\
&= \cos\theta + i \sin\theta
\end{aligned}
$$

实部正是 $\cos\theta$ 的泰勒级数，虚部正是 $\sin\theta$ 的泰勒级数。欧拉公式得证。

利用欧拉公式，三角形式立即化为**指数形式**：

$$
\boxed{z = r e^{i\theta}} \tag{P.3.3}
$$

指数形式使乘除、幂运算极为简便。例如 $z_1 z_2 = r_1 r_2 e^{i(\theta_1 + \theta_2)}$：模相乘，辐角相加。

三种形式总结：
- 代数形式 $z = x + i y$ ——便于加减
- 三角形式 $z = r(\cos\theta + i\sin\theta)$ ——体现模与辐角的分离
- 指数形式 $z = r e^{i\theta}$ ——便于乘除与微积分

#### 1.4 共轭复数

$z = x + i y$ 的共轭复数定义为

$$
\bar{z} = z^* = x - i y = r e^{-i\theta}
$$

模不变，辐角反号。基本性质：

$$
|z|^2 = z z^*, \quad \operatorname{Re}(z) = \frac{z + z^*}{2}, \quad \operatorname{Im}(z) = \frac{z - z^*}{2i}
$$

---

### 二、复数的运算

#### 2.1 加减法

$$
z_1 \pm z_2 = (x_1 \pm x_2) + i (y_1 \pm y_2)
$$

用代数形式最为方便。

#### 2.2 乘除法

设 $z_1 = r_1 e^{i\theta_1}$, $z_2 = r_2 e^{i\theta_2}$，

- 乘法：$z_1 z_2 = r_1 r_2 e^{i(\theta_1 + \theta_2)}$
- 除法：$\dfrac{z_1}{z_2} = \dfrac{r_1}{r_2} e^{i(\theta_1 - \theta_2)} \quad (z_2 \neq 0)$

在代数形式下，除法通过乘以分母的共轭实现：

$$
\frac{x_1 + i y_1}{x_2 + i y_2} = \frac{(x_1 + i y_1)(x_2 - i y_2)}{x_2^2 + y_2^2}
$$

#### 2.3 幂与方根

由指数形式，$n$ 次幂为

$$
z^n = r^n e^{i n \theta}
$$

$n$ 次方根有 $n$ 个值（复数域中代数基本定理的体现）：

$$
z^{1/n} = r^{1/n} \, \exp\!\left( i \, \frac{\theta + 2\pi k}{n} \right), \quad k = 0, 1, 2, \dots, n-1
$$

#### 2.4 几个常用恒等式

由 $e^{\pm i\theta} = \cos\theta \pm i \sin\theta$，可反解出

$$
\cos\theta = \frac{e^{i\theta} + e^{-i\theta}}{2}, \quad \sin\theta = \frac{e^{i\theta} - e^{-i\theta}}{2i} \tag{P.3.4}
$$

这就是用复指数表示三角函数的恒等式，在将 $\cos^2\theta$、$\sin\theta\cos\theta$ 等化为倍角时特别方便，也为相量分析中将实数正弦量“提升”为复数指数形式提供了数学基础。

复数恒等式 $\operatorname{Re}(z_1 + z_2) = \operatorname{Re}(z_1) + \operatorname{Re}(z_2)$ 是线性的，但乘法不满足实部可分离：$\operatorname{Re}(z_1 z_2) \neq \operatorname{Re}(z_1) \operatorname{Re}(z_2)$。这在使用相量法计算功率时需要特别注意，将在后续功率计算中详细讨论。

---

### 三、时谐量的相量表示

#### 3.1 时谐实函数与复表示

交流电路和时谐电磁场中最常遇到的激励和响应形式为**正弦（余弦）时谐量**：

$$
a(t) = A_m \cos(\omega t + \phi) \tag{P.3.5}
$$

其中 $A_m$ 为振幅，$\omega$ 为角频率（单位 $\text{rad/s}$），$\phi$ 为初相位。

任何一个余弦函数都可以看作某个复指数函数的实部：

$$
A_m \cos(\omega t + \phi) = \operatorname{Re}\!\big[ A_m e^{i(\omega t + \phi)} \big] = \operatorname{Re}\!\big[ \tilde{A} e^{i\omega t} \big]
$$

其中

$$
\boxed{\tilde{A} \equiv A_m e^{i\phi}} \tag{P.3.6}
$$

称为该时谐量的**相量**（phasor）。相量是一个复数，其模等于时谐量的振幅 $A_m$，其辐角等于时谐量的初相位 $\phi$。

#### 3.2 相量运算的基本规则

在线性系统中，如果所有激励都是同频率的时谐量，且系统已达到稳态，那么所有响应也都是同频率的时谐量。此时，我们可以将全部时谐量用其相量表示，在复数域中进行运算，最后取实部恢复时间函数。

**规则 1（线性叠加）**：若 $a(t) \leftrightarrow \tilde{A}$, $b(t) \leftrightarrow \tilde{B}$，则

$$
\alpha a(t) + \beta b(t) \quad \leftrightarrow \quad \alpha \tilde{A} + \beta \tilde{B}
$$

其中 $\alpha, \beta$ 为实常数。

**规则 2（微分）**：若 $a(t) \leftrightarrow \tilde{A}$，则

$$
\frac{da}{dt} \quad \leftrightarrow \quad i\omega \tilde{A}
$$

因为 $\frac{d}{dt} [ A_m \cos(\omega t + \phi) ] = -\omega A_m \sin(\omega t + \phi) = \omega A_m \cos(\omega t + \phi + \pi/2)$，这恰好对应于相量乘以 $i\omega$（因为 $i = e^{i\pi/2}$，乘 $i$ 使相位增加 $\pi/2$）。严格地：

$$
\frac{d}{dt} \operatorname{Re}(\tilde{A} e^{i\omega t}) = \operatorname{Re}(i\omega \tilde{A} e^{i\omega t})
$$

**规则 3（积分）**：若 $a(t) \leftrightarrow \tilde{A}$，则

$$
\int a(t) dt \quad \leftrightarrow \quad \frac{\tilde{A}}{i\omega} \quad (\text{忽略积分常数，对应稳态解})
$$

这些规则意味着：在相量域中，对时间的微分化为乘以 $i\omega$，对时间的积分化为除以 $i\omega$。微分方程化为代数方程。

#### 3.3 从相量恢复时域函数

计算得到相量 $\tilde{A}$ 后，恢复时域函数的步骤为：

$$
a(t) = \operatorname{Re}\!\big[ \tilde{A} e^{i\omega t} \big]
$$

如果相量是以极坐标形式给出的 $\tilde{A} = A_m e^{i\phi}$，那么

$$
a(t) = A_m \cos(\omega t + \phi)
$$

关于振幅的约定：相量的模通常取为时谐函数的**振幅** $A_m$。但在电力工程中，相量的模常取为**有效值**（RMS 值）$A_{\text{rms}} = A_m / \sqrt{2}$。本课程采用振幅约定，即 $\tilde{A} = A_m e^{i\phi}$，与电磁波方程中的复振幅保持一致。

#### 3.4 相量图

复数可以在复平面上用矢量表示。同频率的各相量在复平面上构成一个“定格”的图像——它们之间的夹角对应相位差，长度对应振幅比。由于所有相量都以同一角速度 $\omega$ 绕原点旋转（乘以 $e^{i\omega t}$），取实部即为物理量在实轴上的投影，相量之间的相对关系在匀速旋转中保持不变。因此，用静止的相量图分析各量间的相位关系和幅值关系十分直观。在分析 RLC 电路中电压与电流的超前与滞后关系时，相量图是极为有效的辅助工具。

---

### 四、相量法在时谐场方程中的应用

#### 4.1 麦克斯韦方程组的相量形式

电磁场中的时变麦克斯韦方程组包含时间导数。在时谐场假设下，所有场量均以 $e^{i\omega t}$ 规律变化，$\partial/\partial t \to i\omega$。麦克斯韦方程组在相量形式下化为不含时间的复代数-微分方程（此处以真空中为例）：

$$
\begin{aligned}
\nabla \times \tilde{\mathbf{E}} &= -i\omega \tilde{\mathbf{B}} \\
\nabla \times \tilde{\mathbf{H}} &= \tilde{\mathbf{J}}_f + i\omega \tilde{\mathbf{D}} \\
\nabla \cdot \tilde{\mathbf{D}} &= \tilde{\rho}_f \\
\nabla \cdot \tilde{\mathbf{B}} &= 0
\end{aligned}
$$

这样，时谐电磁场的求解不再需要处理时间变量 $t$。空间偏微分方程得到后，配以适当的边界条件即可求解复振幅 $\tilde{\mathbf{E}}, \tilde{\mathbf{B}}$ 等，物理场只需最后取实部乘 $e^{i\omega t}$。

#### 4.2 电路方程的相量化

在 RLC 电路中，电感元件的电压-电流关系 $v_L = L \, di/dt$，在相量域中化为 $\tilde{V}_L = i\omega L \tilde{I}$；电容元件 $i_C = C \, dv/dt$ 化为 $\tilde{I}_C = i\omega C \tilde{V}_C$，即 $\tilde{V}_C = \frac{1}{i\omega C} \tilde{I}$。基尔霍夫定律在相量形式下为各相量的代数和为零。整个交流电路的分析因此简化为复数代数方程组的求解，与直流电阻电路在形式上完全类似，只是将电阻推广为复数阻抗。

---

### 五、阻抗与导纳的初步概念

#### 5.1 阻抗的定义

对于交流电路中的一个二端元件（或二端网络），定义其**阻抗**为电压相量与电流相量之比：

$$
\boxed{Z \equiv \frac{\tilde{V}}{\tilde{I}}} \tag{P.3.7}
$$

阻抗是复数，单位为欧姆（$\Omega$）。其模 $|Z| = V_m / I_m$ 表示电压振幅与电流振幅之比，辐角 $\arg Z = \phi_v - \phi_i$ 表示电压超前于电流的相位角。

电阻 $R$、电感 $L$、电容 $C$ 三种基本元件的阻抗分别为

$$
Z_R = R, \quad Z_L = i\omega L, \quad Z_C = \frac{1}{i\omega C} = -i \frac{1}{\omega C} \tag{P.3.8}
$$

#### 5.2 导纳

阻抗的倒数称为**导纳**：

$$
Y \equiv \frac{1}{Z} = G + i B \quad (\text{西门子}, \text{S}) \tag{P.3.9}
$$

实部 $G$ 为电导，虚部 $B$ 为电纳。导纳在并联电路分析中尤为方便：并联元件的导纳直接相加。

#### 5.3 阻抗的串联与并联

相量域中，基尔霍夫定律成立，阻抗的串并联规则与电阻完全相同：

- 串联：$Z_{\text{eq}} = Z_1 + Z_2 + \cdots$
- 并联：$\dfrac{1}{Z_{\text{eq}}} = \dfrac{1}{Z_1} + \dfrac{1}{Z_2} + \cdots$ 或 $Y_{\text{eq}} = Y_1 + Y_2 + \cdots$

利用这些规则，可以将复杂的交流电路化简为等效阻抗，从而求得总电流、各支路电流和各元件电压。

---

复数从代数、三角到指数形式的统一描述，由欧拉公式 $e^{i\theta} = \cos\theta + i\sin\theta$ 架起了连接指数函数与三角函数的桥梁。时谐量的相量表示将时间微积分运算转化为代数运算（$\partial/\partial t \to i\omega$），使得交流电路的稳态分析和时谐电磁场的求解大大简化。阻抗和导纳作为相量电压与相量电流之比，在交流电路中扮演了与直流电路中电阻完全相同的角色，串并联规则完全一致，使得成熟的直流电路分析技术可全套移植到相量域中。这套复数与相量方法是整个电磁学后半部分——从 RLC 电路到麦克斯韦方程组时谐解——的统一计算语言。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [1.1 电荷与库仑定律](https://zhuanlan.zhihu.com/p/2062629338571453016)

**课程目录：** [电磁学目录](https://zhuanlan.zhihu.com/p/2062618184079954663)
<!-- zhihu-nav-bottom:end -->
