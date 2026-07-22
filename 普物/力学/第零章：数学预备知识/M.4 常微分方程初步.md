---
created: 2026-07-20T06:03:40+08:00
modified: 2026-07-20T06:35:17+08:00
zhihu-title: M.4 常微分方程初步
zhihu-topics:
  - 力学
zhihu-link: https://zhuanlan.zhihu.com/p/2062423944645415902
zhihu-toc: true
zhihu-created-at: 2026-07-20 10:40
zhihu-updated-at: 2026-07-21 14:24
---
<!-- zhihu-nav-top:start -->
**上一节：** [M.3 坐标系](https://zhuanlan.zhihu.com/p/2062423889158861169)
<!-- zhihu-nav-top:end -->

## M.4 常微分方程初步

力学中，牛顿第二定律 $\mathbf{F} = m\mathbf{a}$ 通常表现为一个微分方程——因为它把力（位置的函数）和加速度（位置的二阶导数）联系了起来。弹簧振子、单摆、行星轨道、阻尼运动、受迫振动……所有这些问题的核心都是“解一个微分方程”。

这一节我们系统梳理力学中最常遇到的几类微分方程：一阶线性方程、二阶常系数齐次方程、二阶常系数非齐次方程。我们的目标不是抽象地讨论微分方程理论，而是为后续各章的运动方程求解提供可以直接使用的工具。

---

### 一、一阶线性微分方程

**标准形式**：

$$
\frac{dx}{dt} + a(t)x = b(t). \tag{1}
$$

其中 $a(t)$ 和 $b(t)$ 是已知函数。这个方程称为一阶线性微分方程，因为它包含未知函数 $x$ 和它的一阶导数，且 $x$ 和 $\frac{dx}{dt}$ 都是一次幂。如果 $b(t) = 0$，称为齐次方程；如果 $b(t) \neq 0$，称为非齐次方程。

#### 1.1 齐次方程

先看最简单的情况 $b(t) = 0$：

$$
\frac{dx}{dt} + a(t)x = 0. \tag{2}
$$

分离变量：

$$
\frac{dx}{x} = -a(t)dt. \tag{3}
$$

两边积分：

$$
\ln|x| = -\int a(t)dt + C. \tag{4}
$$

取指数：

$$
x(t) = C e^{-\int a(t)dt}. \tag{5}
$$

其中 $C$ 是任意常数，由初始条件 $x(0) = x_0$ 确定。

**例子**：阻尼运动 $\frac{dv}{dt} + \frac{k}{m}v = 0$。这里 $a(t) = k/m$ 是常数。通解：

$$
v(t) = C e^{-(k/m)t}. \tag{6}
$$

若初始速度为 $v(0) = v_0$，则 $C = v_0$，得到：

$$
v(t) = v_0 e^{-(k/m)t}. \tag{7}
$$

#### 1.2 非齐次方程——积分因子法

对于非齐次方程(1)，分离变量法不再适用。标准解法是**积分因子法**。

寻找一个函数 $\mu(t)$（称为积分因子），使得方程两边乘以 $\mu(t)$ 后，左边恰好是一个“乘积的导数”：

$$
\mu(t)\frac{dx}{dt} + \mu(t)a(t)x = \frac{d}{dt}[\mu(t)x]. \tag{8}
$$

展开右边的导数：

$$
\frac{d}{dt}[\mu(t)x] = \mu(t)\frac{dx}{dt} + \frac{d\mu}{dt}x. \tag{9}
$$

比较(8)和(9)，要求：

$$
\mu(t)a(t) = \frac{d\mu}{dt}. \tag{10}
$$

整理：

$$
\frac{d\mu}{dt} = \mu(t)a(t). \tag{11}
$$

分离变量：

$$
\frac{d\mu}{\mu} = a(t)dt. \tag{12}
$$

积分：

$$
\ln|\mu| = \int a(t)dt \quad \Rightarrow \quad \mu(t) = e^{\int a(t)dt}. \tag{13}
$$

（积分常数可以取 $0$，因为我们只需要一个特定的 $\mu$ 即可。）

现在将(1)两边乘以 $\mu(t)$：

$$
\mu\frac{dx}{dt} + \mu a x = \mu b. \tag{14}
$$

左边正好是 $\frac{d}{dt}(\mu x)$，因此：

$$
\frac{d}{dt}(\mu x) = \mu b. \tag{15}
$$

两边积分：

$$
\mu(t)x(t) = \int \mu(t)b(t)dt + C. \tag{16}
$$

于是：

$$
x(t) = \frac{1}{\mu(t)}\left[\int \mu(t)b(t)dt + C\right]. \tag{17}
$$

代入 $\mu(t) = e^{\int a(t)dt}$：

$$
\boxed{x(t) = e^{-\int a(t)dt}\left[\int e^{\int a(t)dt} b(t)dt + C\right]}. \tag{18}
$$

这就是一阶线性微分方程的通解公式。

**例子**：电路中的 $RL$ 回路（仅作类比，不用深究物理）$\frac{dI}{dt} + \frac{R}{L}I = \frac{V}{L}$。这里 $a = R/L$ 为常数，$b = V/L$ 为常数。$\int a\,dt = (R/L)t$，$\mu(t) = e^{(R/L)t}$。代入(18)：

$$
I(t) = e^{-(R/L)t}\left[\int e^{(R/L)t}\frac{V}{L}dt + C\right] = e^{-(R/L)t}\left[\frac{V}{R}e^{(R/L)t} + C\right] = \frac{V}{R} + C e^{-(R/L)t}.
$$

若初始电流 $I(0) = 0$，则 $C = -V/R$，得到：

$$
I(t) = \frac{V}{R}\left(1 - e^{-(R/L)t}\right). \tag{19}
$$

在力学中，一阶线性方程最常见的场景是速度依赖的阻尼力 $F = -kv$。我们已经在 M.2 中用分离变量法处理过了。积分因子法是它的系统化推广。

---

### 二、二阶线性常系数齐次方程

**标准形式**：

$$
\ddot{x} + p\dot{x} + qx = 0. \tag{20}
$$

其中 $p, q$ 是常数。这是二阶线性常系数齐次方程。“二阶”指最高阶导数为二阶；“常系数”指 $p, q$ 不随时间变化；“齐次”指右边为 $0$。

这类方程在力学中出现频率极高——简谐振动、阻尼振动都归于此。

**求解思路**：指数函数有一个重要性质——求导后仍是指数函数。设试探解为 $x(t) = e^{\lambda t}$，其中 $\lambda$ 是待定常数。代入(20)：

$$
\ddot{x} = \lambda^2 e^{\lambda t}, \quad \dot{x} = \lambda e^{\lambda t}.
$$

代入：

$$
(\lambda^2 + p\lambda + q)e^{\lambda t} = 0.
$$

由于 $e^{\lambda t} \neq 0$，必须有：

$$
\lambda^2 + p\lambda + q = 0. \tag{21}
$$

这称为**特征方程**（characteristic equation）。它是一个关于 $\lambda$ 的一元二次方程，解为：

$$
\lambda_{1,2} = \frac{-p \pm \sqrt{p^2 - 4q}}{2}. \tag{22}
$$

根据判别式 $\Delta = p^2 - 4q$ 的符号，有三种情况。

#### 2.1 情况一：$\Delta > 0$（过阻尼）

两个不同的实根 $\lambda_1, \lambda_2$。通解为：

$$
x(t) = C_1 e^{\lambda_1 t} + C_2 e^{\lambda_2 t}. \tag{23}
$$

由于 $\Delta > 0$，两个根都是负数（因为 $p>0$ 时阻尼项为正，系统耗散能量）。$x(t)$ 随时间指数衰减到零，没有振荡。这对应过阻尼状态。

**物理例子**：粘性很大的流体中的阻尼运动，系统缓慢回到平衡位置。

#### 2.2 情况二：$\Delta = 0$（临界阻尼）

重根 $\lambda_1 = \lambda_2 = -p/2$。如果仍用(23)的形式，$C_1 e^{\lambda t} + C_2 e^{\lambda t} = (C_1+C_2)e^{\lambda t}$，无法覆盖所有可能的初始条件。需要寻找第二个线性无关的解。

当特征方程有重根 $\lambda$ 时，方程的通解为：

$$
x(t) = (C_1 + C_2 t)e^{\lambda t}. \tag{24}
$$

验证：设 $x = (C_1 + C_2 t)e^{\lambda t}$，代入方程可验证满足。这是二阶常系数齐次方程重根情况下的标准解法。

**物理意义**：临界阻尼系统恢复到平衡位置的速度最快，且没有过冲。这是许多工程系统追求的状态——例如汽车减震器的理想阻尼。

#### 2.3 情况三：$\Delta < 0$（欠阻尼——简谐振动）

两个共轭复根。设 $p = 2\gamma$，$q = \omega_0^2$（这对应弹簧-阻尼系统的标准记法，$\gamma$ 是阻尼系数，$\omega_0$ 是无阻尼固有频率）：

特征方程为：

$$
\lambda^2 + 2\gamma\lambda + \omega_0^2 = 0. \tag{25}
$$

解为：

$$
\lambda_{1,2} = -\gamma \pm i\sqrt{\omega_0^2 - \gamma^2}. \tag{26}
$$

令 $\omega = \sqrt{\omega_0^2 - \gamma^2}$（阻尼振荡频率），则 $\lambda_{1,2} = -\gamma \pm i\omega$。

通解为：

$$
x(t) = C_1 e^{(-\gamma + i\omega)t} + C_2 e^{(-\gamma - i\omega)t} = e^{-\gamma t}\left(C_1 e^{i\omega t} + C_2 e^{-i\omega t}\right). \tag{27}
$$

利用欧拉公式 $e^{\pm i\omega t} = \cos\omega t \pm i\sin\omega t$，可以将通解改写为实函数形式。设 $C_1 = \frac{A}{2}e^{i\phi}$，$C_2 = \frac{A}{2}e^{-i\phi}$（这是为了把两个复数常数合并为两个实常数 $A$ 和 $\phi$），则：

$$
x(t) = A e^{-\gamma t}\cos(\omega t + \phi). \tag{28}
$$

其中 $A$ 和 $\phi$ 由初始条件确定。

**特殊情况 $\gamma = 0$ **（无阻尼）：

这就是**简谐振动**的方程：

$$
\ddot{x} + \omega_0^2 x = 0. \tag{29}
$$

通解为：

$$
x(t) = A\cos(\omega_0 t + \phi). \tag{30}
$$

或者等价地写为：

$$
x(t) = A_1\cos\omega_0 t + A_2\sin\omega_0 t. \tag{31}
$$

两种形式等价，只是常数 $A, \phi$ 和 $A_1, A_2$ 的对应关系为：

$$
A = \sqrt{A_1^2 + A_2^2}, \quad \tan\phi = -\frac{A_2}{A_1}. \tag{32}
$$

**初始条件确定常数**：设 $t=0$ 时 $x(0) = x_0$，$\dot{x}(0) = v_0$。利用形式(31)：

$$
x(0) = A_1 = x_0, \quad \dot{x}(0) = A_2\omega_0 = v_0 \Rightarrow A_2 = \frac{v_0}{\omega_0}.
$$

于是：

$$
x(t) = x_0\cos\omega_0 t + \frac{v_0}{\omega_0}\sin\omega_0 t. \tag{33}
$$

这就是由初始位置和初始速度唯一确定的运动。

---

### 三、二阶线性常系数非齐次方程

**标准形式**：

$$
\ddot{x} + p\dot{x} + qx = f(t). \tag{34}
$$

其中 $f(t)$ 是已知函数（称为驱动项或外源项）。方程的通解由两部分组成：

$$
x(t) = x_h(t) + x_p(t). \tag{35}
$$

其中 $x_h(t)$ 是齐次方程(20)的通解，$x_p(t)$ 是任意一个特解。

这一节我们只讨论一种情况——简谐驱动力：

$$
\ddot{x} + 2\gamma\dot{x} + \omega_0^2 x = \frac{F_0}{m}\cos\omega t. \tag{36}
$$

这是受迫振动的标准方程。$\omega$ 是驱动力的角频率，$\omega_0$ 是系统的固有频率，$\gamma$ 是阻尼系数。

**求特解**：由于驱动力是 $\cos\omega t$，特解应取相同频率的振荡。设：

$$
x_p(t) = A\cos\omega t + B\sin\omega t. \tag{37}
$$

代入(36)（这里不展开全部代数，只给结果）。更方便的处理是用复数表示法。令驱动项为 $\frac{F_0}{m}e^{i\omega t}$，设特解为 $x_p(t) = \text{Re}[\tilde{A} e^{i\omega t}]$。代入后得到复振幅：

$$
\tilde{A} = \frac{F_0/m}{\omega_0^2 - \omega^2 + 2i\gamma\omega}. \tag{38}
$$

振幅（实振幅）为：

$$
A_0 = |\tilde{A}| = \frac{F_0/m}{\sqrt{(\omega_0^2 - \omega^2)^2 + (2\gamma\omega)^2}}. \tag{39}
$$

相位为：

$$
\delta = \arctan\left(\frac{2\gamma\omega}{\omega_0^2 - \omega^2}\right). \tag{40}
$$

因此特解为：

$$
x_p(t) = A_0\cos(\omega t - \delta). \tag{41}
$$

完整的解为齐次解（衰减项）加特解（稳态项）：

$$
x(t) = A_h e^{-\gamma t}\cos(\omega_d t + \phi) + \frac{F_0/m}{\sqrt{(\omega_0^2 - \omega^2)^2 + (2\gamma\omega)^2}}\cos(\omega t - \delta). \tag{42}
$$

其中 $\omega_d = \sqrt{\omega_0^2 - \gamma^2}$。

**共振**：当驱动频率 $\omega$ 接近固有频率 $\omega_0$ 时，振幅(39)的分母取得最小值，振幅急剧增大。在 $\omega = \sqrt{\omega_0^2 - 2\gamma^2}$ 处振幅取最大值。当阻尼 $\gamma$ 很小时，共振频率几乎等于 $\omega_0$。

---

### 四、初始条件确定积分常数——一个完整的例子

用二阶非齐次方程的完整求解来展示“初始条件如何确定常数”。

设无阻尼受迫振动：

$$
\ddot{x} + \omega_0^2 x = \frac{F_0}{m}\cos\omega t. \tag{43}
$$

初始条件为 $x(0) = 0$，$\dot{x}(0) = 0$（从静止开始）。

齐次解：$x_h(t) = C_1\cos\omega_0 t + C_2\sin\omega_0 t$。

特解（由(39)和(41)，这里 $\gamma=0$）：$x_p(t) = \frac{F_0/m}{\omega_0^2 - \omega^2}\cos\omega t$。（如果 $\omega = \omega_0$，则需另作处理，这是共振情况。）

因此通解：

$$
x(t) = C_1\cos\omega_0 t + C_2\sin\omega_0 t + \frac{F_0/m}{\omega_0^2 - \omega^2}\cos\omega t.
$$

代入 $x(0) = 0$：

$$
C_1 + \frac{F_0/m}{\omega_0^2 - \omega^2} = 0 \quad \Rightarrow \quad C_1 = -\frac{F_0/m}{\omega_0^2 - \omega^2}.
$$

$\dot{x}(t) = -C_1\omega_0\sin\omega_0 t + C_2\omega_0\cos\omega_0 t - \frac{F_0/m}{\omega_0^2 - \omega^2}\omega\sin\omega t$。

代入 $\dot{x}(0) = 0$：

$$
C_2\omega_0 = 0 \quad \Rightarrow \quad C_2 = 0.
$$

因此：

$$
x(t) = \frac{F_0/m}{\omega_0^2 - \omega^2}(\cos\omega t - \cos\omega_0 t). \tag{44}
$$

这就是从静止出发的受迫振动响应。物理图像：拍频现象，当 $\omega$ 接近 $\omega_0$ 时，振幅随时间缓慢变化，出现“拍”。

---

### 五、总结

把 M.4 节的全部内容压缩成几条线：

1. **一阶线性方程** $\dot{x} + a(t)x = b(t)$：积分因子法，通解为 $x(t) = e^{-\int a dt}[\int e^{\int a dt} b\,dt + C]$。

2. **二阶常系数齐次方程** $\ddot{x} + p\dot{x} + qx = 0$：特征方程 $\lambda^2 + p\lambda + q = 0$。根据判别式 $\Delta = p^2 - 4q$ 分三种情况：
   - $\Delta > 0$：过阻尼，$x = C_1 e^{\lambda_1 t} + C_2 e^{\lambda_2 t}$
   - $\Delta = 0$：临界阻尼，$x = (C_1 + C_2 t)e^{\lambda t}$
   - $\Delta < 0$：欠阻尼（简谐振动），$x = A e^{-\gamma t}\cos(\omega t + \phi)$

3. **二阶常系数非齐次方程** $\ddot{x} + p\dot{x} + qx = f(t)$：通解 = 齐次解 + 特解。对简谐驱动力 $f(t) = F_0\cos\omega t$，特解振幅为 $\frac{F_0/m}{\sqrt{(\omega_0^2-\omega^2)^2 + (2\gamma\omega)^2}}$，在 $\omega \approx \omega_0$ 时发生共振。

4. **初始条件确定常数**：常数由 $x(0)$ 和 $\dot{x}(0)$ 唯一确定。$n$ 阶微分方程需要 $n$ 个初始条件。

微分方程是物理定律的数学形式。掌握了这些基本方程的解法，就掌握了从牛顿定律出发预测运动轨迹的核心工具。后续各章中，弹簧振子、单摆、阻尼运动、受迫振动、电路振荡——所有这些现象的解都来自这一节的方法。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [M.5 泰勒展开与近似方法](https://zhuanlan.zhihu.com/p/2062423971132355087)

**课程目录：** [力学目录](https://zhuanlan.zhihu.com/p/2062487053170906945)
<!-- zhihu-nav-bottom:end -->
