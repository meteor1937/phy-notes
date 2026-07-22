---
created: 2026-07-20T17:22:30+08:00
modified: 2026-07-20T17:29:20+08:00
zhihu-title: P.1 振动与波动回顾
zhihu-topics:
  - 光学
zhihu-link: https://zhuanlan.zhihu.com/p/2062589956233737211
zhihu-toc: true
zhihu-created-at: 2026-07-20 17:29
zhihu-updated-at: 2026-07-21 14:37
---
## P.1 振动与波动回顾

这一节的任务是把振动和波动中最基本、最核心的语言建立起来。我们从质点动力学出发，完整推导简谐振动的运动方程和解，引入相位与复数表示；然后将振动传播出去，得到一维简谐波，再从连续介质的微元分析导出波动方程，讨论行波解、叠加原理、干涉、驻波，最后用两列频率相近的波叠加引出相速度和群速度的区别。

---

### 一、简谐振动

#### 1.1 运动方程的建立

最简单的振动系统是一端固定的轻弹簧连接一个质量为 $m$ 的质点，在光滑水平面上运动。取平衡位置为坐标原点 $x=0$，质点偏离平衡位置的位移为 $x(t)$。弹簧的恢复力满足胡克定律 $F = -k x$，其中 $k$ 是劲度系数。根据牛顿第二定律：

$$
m \frac{d^2 x}{dt^2} = -k x
$$

引入常数

$$
\omega^2 \equiv \frac{k}{m} \quad (\omega > 0)
$$

方程化为

$$
\frac{d^2 x}{dt^2} + \omega^2 x = 0 \tag{1}
$$

这是一个二阶常系数线性齐次微分方程。

#### 1.2 方程的求解

设试探解 $x(t) = e^{\lambda t}$，代入 (1) 式：

$$
\frac{d^2}{dt^2} e^{\lambda t} + \omega^2 e^{\lambda t} = \lambda^2 e^{\lambda t} + \omega^2 e^{\lambda t} = (\lambda^2 + \omega^2) e^{\lambda t} = 0
$$

因为 $e^{\lambda t} \neq 0$，得到特征方程

$$
\lambda^2 + \omega^2 = 0 \quad \Rightarrow \quad \lambda = \pm i\omega
$$

通解为两个线性无关解的线性组合：

$$
x(t) = C_1 e^{i\omega t} + C_2 e^{-i\omega t} \tag{2}
$$

其中 $C_1, C_2$ 为任意复常数。因为位移 $x(t)$ 必须是实数，可以令 $C_2 = C_1^*$，或直接化为实数形式。利用欧拉公式

$$
e^{i\theta} = \cos\theta + i\sin\theta, \quad e^{-i\theta} = \cos\theta - i\sin\theta
$$

将 (2) 式展开：

$$
\begin{aligned}
x(t) &= C_1 (\cos\omega t + i\sin\omega t) + C_2 (\cos\omega t - i\sin\omega t) \\
&= (C_1 + C_2) \cos\omega t + i(C_1 - C_2) \sin\omega t
\end{aligned}
$$

令 $A_1 = C_1 + C_2$，$A_2 = i(C_1 - C_2)$，若初始条件为实数，则 $A_1, A_2$ 均为实数，于是通解可以写为

$$
x(t) = A_1 \cos\omega t + A_2 \sin\omega t \tag{3}
$$

更常见的写法是利用三角恒等式合成为一个余弦函数。设

$$
A_1 = A\cos\phi, \quad A_2 = -A\sin\phi
$$

则

$$
\begin{aligned}
x(t) &= A\cos\phi \cos\omega t - A\sin\phi \sin\omega t \\
&= A\cos(\omega t + \phi)
\end{aligned}
$$

其中振幅 $A = \sqrt{A_1^2 + A_2^2}$，初相位 $\phi$ 满足 $\tan\phi = -A_2/A_1$。这样我们得到了简谐振动的标准形式：

$$
x(t) = A \cos(\omega t + \phi) \tag{4}
$$

式中：
- $\omega = \sqrt{k/m}$ 称为角频率（单位 rad/s）；
- $A$ 是振幅，由初始能量决定；
- $\phi$ 是初相位，由初始位移和速度共同决定。

速度与加速度：

$$
v(t) = \frac{dx}{dt} = -\omega A \sin(\omega t + \phi) = \omega A \cos\!\left(\omega t + \phi + \frac{\pi}{2}\right)
$$

$$
a(t) = \frac{d^2 x}{dt^2} = -\omega^2 A \cos(\omega t + \phi) = -\omega^2 x(t)
$$

#### 1.3 复数表示法

在涉及叠加和相位运算时，直接用三角函数比较繁琐。引入复数表示可以极大简化计算。将实数解 $x(t) = A\cos(\omega t + \phi)$ 视为一个复函数的实部：

$$
x(t) = \operatorname{Re}\!\big[ \tilde{x}(t) \big], \quad 
\tilde{x}(t) = A e^{i(\omega t + \phi)} = \tilde{A} e^{i\omega t}
$$

其中复振幅 $\tilde{A} \equiv A e^{i\phi}$ 同时携带了振幅和初相位的信息。

使用复数表示后，微分运算变成简单的乘法：

$$
\frac{d}{dt} \tilde{x}(t) = i\omega \tilde{x}(t), \quad 
\frac{d^2}{dt^2} \tilde{x}(t) = (i\omega)^2 \tilde{x}(t) = -\omega^2 \tilde{x}(t)
$$

这直接恢复了运动方程 (1)。在涉及多个同频率振动叠加时，只需对复振幅进行代数相加，最后取实部即可得到物理结果。

---

### 二、从振动到波动：一维简谐波

#### 2.1 行波的波函数

若振动的质点不是孤立的，而是连续介质的一部分，振动就会以有限速度向周围传播，形成波。我们以沿 $x$ 轴拉紧的弦上的横波为例。设 $t$ 时刻位于 $x$ 处的质点偏离平衡位置的横向位移为 $y(x,t)$。

假定 $x=0$ 处的质点做简谐振动：

$$
y(0,t) = A \cos(\omega t + \phi)
$$

如果这个振动状态以恒定速率 $v$ 沿着 $+x$ 方向无畸变地传播，那么 $x>0$ 处质点的振动形式与原点相同，只是在时间上滞后了 $x/v$。也就是说，$x$ 处 $t$ 时刻的位移等于原点在 $t - x/v$ 时刻的位移：

$$
y(x,t) = y\!\left(0, t - \frac{x}{v}\right) = A \cos\!\left[\omega\left(t - \frac{x}{v}\right) + \phi\right] \tag{5}
$$

这就得到了沿 $+x$ 方向传播的简谐行波的波函数。为书写方便，定义角波数

$$
k \equiv \frac{\omega}{v} \quad (\text{单位 rad/m})
$$

则波函数写成更对称的形式：

$$
y(x,t) = A \cos(kx - \omega t + \phi') \tag{6}
$$

其中 $\phi' = -\phi$（相位可重新定义）。通常将简谐波函数写为

$$
y(x,t) = A \cos(kx - \omega t + \phi_0) \tag{7}
$$

或使用复数表示

$$
\tilde{y}(x,t) = A e^{i(kx - \omega t + \phi_0)} \tag{8}
$$

物理位移为 $y = \operatorname{Re}[\tilde{y}]$。

#### 2.2 波长、周期、频率和波速

波函数中的相位为 $\Phi(x,t) = kx - \omega t + \phi_0$。在某固定时刻 $t$，空间上相邻两同相位点之间的距离称为波长 $\lambda$。即

$$
k(x+\lambda) - \omega t = kx - \omega t + 2\pi \quad \Rightarrow \quad k\lambda = 2\pi
$$

$$
\lambda = \frac{2\pi}{k} \tag{9}
$$

类似地，在某固定位置 $x$，相位变化 $2\pi$ 所需的时间称为周期 $T$：

$$
kx - \omega(t+T) = kx - \omega t - 2\pi \quad \Rightarrow \quad \omega T = 2\pi
$$

$$
T = \frac{2\pi}{\omega}, \quad f = \frac{1}{T} = \frac{\omega}{2\pi} \tag{10}
$$

波速 $v$ 由等相位面的移动速度给出。令 $\Phi(x,t) = \text{常数}$，对 $t$ 求导：

$$
\frac{d\Phi}{dt} = k\frac{dx}{dt} - \omega = 0 \quad \Rightarrow \quad \frac{dx}{dt} = \frac{\omega}{k}
$$

这正是波的传播速度：

$$
v = \frac{\omega}{k} = \lambda f \tag{11}
$$

---

### 三、波动方程的导出

简谐波 (7) 是波动方程的解。波动方程本身可以从连续介质的动力学推导出来。我们仍然用张紧的弦作为模型，导出支配 $y(x,t)$ 的偏微分方程。

#### 3.1 微元受力分析

设弦的线密度为 $\mu$（单位 kg/m），弦上张力大小为 $T$，小振动假设下 $T$ 近似不变，且振动方向沿 $y$ 轴，位移很小，故弦上各点的切线斜率 $\partial y/\partial x \ll 1$。

取弦上一小段，从 $x$ 到 $x+dx$，其质量为 $\mu\,dx$。两端受到大小近似为 $T$ 的张力，方向沿弦的切线。左端张力的 $y$ 分量为

$$
F_y(x) = -T \sin\theta_1 \approx -T \left.\frac{\partial y}{\partial x}\right|_{x}
$$

右端张力的 $y$ 分量为

$$
F_y(x+dx) = T \sin\theta_2 \approx T \left.\frac{\partial y}{\partial x}\right|_{x+dx}
$$

合力为

$$
F_{\text{net}} = T \left( \left.\frac{\partial y}{\partial x}\right|_{x+dx} - \left.\frac{\partial y}{\partial x}\right|_{x} \right)
$$

利用偏导数的差分近似：

$$
\left.\frac{\partial y}{\partial x}\right|_{x+dx} - \left.\frac{\partial y}{\partial x}\right|_{x} \approx \frac{\partial^2 y}{\partial x^2} dx
$$

因此

$$
F_{\text{net}} = T \frac{\partial^2 y}{\partial x^2} dx
$$

#### 3.2 牛顿第二定律

对该微元应用牛顿第二定律（横向）：

$$
\mu\,dx \frac{\partial^2 y}{\partial t^2} = T \frac{\partial^2 y}{\partial x^2} dx
$$

化简得到

$$
\frac{\partial^2 y}{\partial t^2} = \frac{T}{\mu} \frac{\partial^2 y}{\partial x^2} \tag{12}
$$

令

$$
v = \sqrt{\frac{T}{\mu}} \tag{13}
$$

即得经典的一维波动方程：

$$
\boxed{\frac{\partial^2 y}{\partial t^2} = v^2 \frac{\partial^2 y}{\partial x^2}} \tag{14}
$$

#### 3.3 行波解：达朗贝尔解

为了求解 (14)，引入新变量

$$
\xi = x - vt, \quad \eta = x + vt
$$

由链式法则：

$$
\frac{\partial}{\partial x} = \frac{\partial}{\partial \xi} + \frac{\partial}{\partial \eta}, \quad
\frac{\partial}{\partial t} = -v\frac{\partial}{\partial \xi} + v\frac{\partial}{\partial \eta}
$$

二阶导数：

$$
\frac{\partial^2}{\partial x^2} = \frac{\partial^2}{\partial \xi^2} + 2\frac{\partial^2}{\partial \xi \partial \eta} + \frac{\partial^2}{\partial \eta^2}
$$

$$
\frac{\partial^2}{\partial t^2} = v^2 \frac{\partial^2}{\partial \xi^2} - 2v^2 \frac{\partial^2}{\partial \xi \partial \eta} + v^2 \frac{\partial^2}{\partial \eta^2}
$$

代入 (14)：

$$
v^2 \left( \frac{\partial^2 y}{\partial \xi^2} - 2\frac{\partial^2 y}{\partial \xi \partial \eta} + \frac{\partial^2 y}{\partial \eta^2} \right)
= v^2 \left( \frac{\partial^2 y}{\partial \xi^2} + 2\frac{\partial^2 y}{\partial \xi \partial \eta} + \frac{\partial^2 y}{\partial \eta^2} \right)
$$

消去共同因子 $v^2$ 并整理，交叉项抵消，得到

$$
\frac{\partial^2 y}{\partial \xi \partial \eta} = 0 \tag{15}
$$

这表明 $\partial y/\partial \xi$ 与 $\eta$ 无关，只依赖于 $\xi$；再对 $\eta$ 积分一次，$y$ 本身是一个只依赖于 $\xi$ 的函数和一个只依赖于 $\eta$ 的函数之和：

$$
y(x,t) = f(\xi) + g(\eta) = f(x - vt) + g(x + vt) \tag{16}
$$

这就是一维波动方程的达朗贝尔通解。$f(x-vt)$ 表示以速率 $v$ 沿 $+x$ 方向传播的行波，$g(x+vt)$ 表示以速率 $v$ 沿 $-x$ 方向传播的行波。简谐波 (7) 对应取 $f(x-vt) = A\cos[k(x-vt)+\phi_0]$ 的特例。

---

### 四、叠加原理与干涉

波动方程 (14) 是线性的：若 $y_1$ 和 $y_2$ 都是解，则任意线性组合 $c_1 y_1 + c_2 y_2$ 也是解。这导致了波的叠加原理。

考虑两列同频率、同振动方向、在空间中叠加的简谐波。设它们在空间中某点的振动分别为

$$
y_1 = A_1 \cos(\omega t + \phi_1), \quad y_2 = A_2 \cos(\omega t + \phi_2)
$$

用复数表示：

$$
\tilde{y}_1 = A_1 e^{i(\omega t + \phi_1)}, \quad \tilde{y}_2 = A_2 e^{i(\omega t + \phi_2)}
$$

合振动 $\tilde{y} = \tilde{y}_1 + \tilde{y}_2 = (A_1 e^{i\phi_1} + A_2 e^{i\phi_2}) e^{i\omega t}$。令合复振幅

$$
\tilde{A} = A_1 e^{i\phi_1} + A_2 e^{i\phi_2}
$$

其模平方为

$$
\begin{aligned}
|\tilde{A}|^2 &= (A_1 e^{i\phi_1} + A_2 e^{i\phi_2})(A_1 e^{-i\phi_1} + A_2 e^{-i\phi_2}) \\
&= A_1^2 + A_2^2 + A_1 A_2 \big(e^{i(\phi_1-\phi_2)} + e^{-i(\phi_1-\phi_2)}\big) \\
&= A_1^2 + A_2^2 + 2 A_1 A_2 \cos(\phi_1 - \phi_2)
\end{aligned}
$$

因此合振动的振幅 $A$ 满足

$$
A^2 = A_1^2 + A_2^2 + 2 A_1 A_2 \cos\Delta\phi \tag{17}
$$

其中 $\Delta\phi = \phi_1 - \phi_2$ 是两列波在该点的相位差。光强（或平均能流）正比于振幅平方，所以

$$
I = I_1 + I_2 + 2\sqrt{I_1 I_2} \cos\Delta\phi \tag{18}
$$

这就是**干涉**的基本公式。当 $\cos\Delta\phi = 1$，即

$$
\Delta\phi = 2m\pi \quad (m = 0, \pm1, \pm2, \dots)
$$

时，$A = A_1 + A_2$，$I$ 最大，称为**相长干涉**；当 $\cos\Delta\phi = -1$，即

$$
\Delta\phi = (2m+1)\pi
$$

时，$A = |A_1 - A_2|$，$I$ 最小，称为**相消干涉**。

相位差通常来源于两列波从波源到观察点的光程不同。若两列波在介质中的波长均为 $\lambda$，光程差 $\delta$ 与相位差的关系为

$$
\Delta\phi = \frac{2\pi}{\lambda} \delta \tag{19}
$$

---

### 五、驻波

驻波是两列振幅相同、频率相同、传播方向相反的行波叠加形成的特殊干涉图样。取两列波为

$$
y_1 = A \cos(kx - \omega t), \quad y_2 = A \cos(kx + \omega t)
$$

叠加：

$$
y = y_1 + y_2 = A \big[ \cos(kx - \omega t) + \cos(kx + \omega t) \big]
$$

利用三角恒等式 $\cos\alpha + \cos\beta = 2\cos\frac{\alpha+\beta}{2} \cos\frac{\alpha-\beta}{2}$：

$$
\begin{align}
y &= A \cdot 2 \cos\!\left( \frac{(kx - \omega t) + (kx + \omega t)}{2} \right)
\cos\!\left( \frac{(kx - \omega t) - (kx + \omega t)}{2} \right) \\
&= 2A \cos(kx) \cos(\omega t) \tag{20}
\end{align}
$$

这个表达式不再具有行波中“相位沿 $x$ 传播”的形式，而是空间部分和时间部分分离。因子 $\cos(\omega t)$ 表明所有点均以相同相位（同相或反相）振动；因子 $2A \cos(kx)$ 给出各点的振幅。

- 当 $|\cos(kx)| = 1$，即 $kx = m\pi$ 或 $x = m\frac{\lambda}{2}$ 时，振幅最大，为 $2A$，这些点称为**波腹**。
- 当 $\cos(kx) = 0$，即 $kx = (2m+1)\frac{\pi}{2}$ 或 $x = (2m+1)\frac{\lambda}{4}$ 时，振幅恒为零，这些点称为**波节**。

节点始终保持静止，相邻节点之间各点做振幅随位置变化的简谐振动。驻波中没有净能量沿 $x$ 方向传播，能量被局限在各段之间来回转换。

---

### 六、波包、相速度与群速度

前面的讨论都针对单一频率的简谐波。然而实际的信号总是由一群频率不同的傅里叶分量叠加而成，形成一个**波包**。波包的运动速度——群速度——通常与单个单色波的相速度不同。

#### 6.1 两波叠加形成波包

为简单起见，考虑两列振幅相同、频率和波数有微小差异的简谐波，都沿 $+x$ 方向传播：

$$
y_1 = A \cos(k_1 x - \omega_1 t), \quad y_2 = A \cos(k_2 x - \omega_2 t)
$$

其中

$$
k_1 = k_0 - \Delta k,\; k_2 = k_0 + \Delta k, \quad
\omega_1 = \omega_0 - \Delta\omega,\; \omega_2 = \omega_0 + \Delta\omega
$$

并且 $\Delta k \ll k_0$，$\Delta\omega \ll \omega_0$。利用和差化积：

$$
\begin{align}
y &= y_1 + y_2 \\
&= A \big[ \cos(k_1 x - \omega_1 t) + \cos(k_2 x - \omega_2 t) \big] \\
&= 2A \cos\!\left( \frac{k_1+k_2}{2}x - \frac{\omega_1+\omega_2}{2}t \right)
\cos\!\left( \frac{k_1-k_2}{2}x - \frac{\omega_1-\omega_2}{2}t \right) \\
&= 2A \cos\!\big( k_0 x - \omega_0 t \big) \,
\cos\!\big( \Delta k \, x - \Delta\omega \, t \big) \tag{21}
\end{align}
$$

这个结果表明，合振动可以看作一个“载波”

$$
\cos(k_0 x - \omega_0 t)
$$

其振幅受到一个低频包络

$$
A_{\text{env}}(x,t) = 2A \cos(\Delta k \, x - \Delta\omega \, t)
$$

的调制。

#### 6.2 相速度

载波的等相位面移动速度就是**相速度**：

$$
v_p = \frac{\omega_0}{k_0} \tag{22}
$$

这就是前面定义的单色波速度。

#### 6.3 群速度

包络的“波峰”移动速度由 $\Delta k \, x - \Delta\omega \, t = \text{常数}$ 决定。对 $t$ 求导：

$$
\Delta k \frac{dx}{dt} - \Delta\omega = 0 \quad \Rightarrow \quad
v_g \equiv \frac{dx}{dt} = \frac{\Delta\omega}{\Delta k}
$$

当 $\Delta k \to 0$ 时，取极限得到**群速度**：

$$
\boxed{v_g = \frac{d\omega}{dk}} \tag{23}
$$

波包整体的能量以群速度传播。将 $\omega = v_p k$ 代入，可得相速度与群速度的关系：

$$
v_g = \frac{d}{dk}(v_p k) = v_p + k \frac{d v_p}{dk} \tag{24}
$$

若介质无色散，即 $v_p$ 不随 $k$ 改变（$d v_p/dk = 0$），则 $v_g = v_p$。若存在色散，$v_g$ 与 $v_p$ 不同：正常色散区 $d v_p/dk < 0$，$v_g < v_p$；反常色散区则可能 $v_g > v_p$。群速度才是在实验中观察到的信号和能量的传播速度。

---

以上从简谐振子的运动方程出发，完整地建立了描述波动所需的全部基本概念和数学工具，包括复振幅、波动方程、行波解、干涉条件、驻波以及相速度与群速度的严格区分。这些内容是光学中干涉、衍射、偏振等一切现象的共同基础。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.2 傅里叶变换与数学补充](https://zhuanlan.zhihu.com/p/2062589994406089592)

**课程目录：** [光学目录](https://zhuanlan.zhihu.com/p/2062614039084192499)
<!-- zhihu-nav-bottom:end -->
