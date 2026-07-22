---
created: 2026-07-20T06:03:37+08:00
modified: 2026-07-20T06:34:33+08:00
zhihu-title: M.3 坐标系
zhihu-topics:
  - 力学
zhihu-link: https://zhuanlan.zhihu.com/p/2062423889158861169
zhihu-toc: true
zhihu-created-at: 2026-07-20 10:40
zhihu-updated-at: 2026-07-21 14:24
---
<!-- zhihu-nav-top:start -->
**上一节：** [M.2 微积分在物理中的基本应用](https://zhuanlan.zhihu.com/p/2062423053758363601)
<!-- zhihu-nav-top:end -->

## M.3 坐标系

坐标系的选取是描述物理运动的第一步。同一个物理现象——比如行星绕太阳运动、小球沿斜面下滑——可以在不同的坐标系中描述。但不同坐标系带来的数学复杂度天差地别。选对了坐标系，运动方程几乎可以“读出来”；选错了，你会在代数迷宫中浪费大量时间。

这一节我们系统梳理力学中最常用的三种坐标系：直角坐标系、自然坐标系、极坐标系。重点放在极坐标系——因为它是描述圆周运动和中心力场的天然语言，也是后续天体力学和刚体力学的基础。我们不仅要给出公式，还要把每个公式的来龙去脉推导清楚，特别是“基向量也随时间变化”这个初学者最容易忽略的关键点。

---

### 一、直角坐标系

直角坐标系 $(x, y, z)$ 是我们在初等数学中最熟悉的坐标系。它的基向量 $\hat{i}, \hat{j}, \hat{k}$ 有三个核心特征：两两垂直、长度均为 $1$、方向固定不随时间变化。

位置矢量（也叫位矢）为：

$$
\mathbf{r}(t) = x(t)\hat{i} + y(t)\hat{j} + z(t)\hat{k}. \tag{1}
$$

由于基向量不随时间变化，速度只需对位置坐标求导：

$$
\mathbf{v}(t) = \frac{d\mathbf{r}}{dt} = \dot{x}\hat{i} + \dot{y}\hat{j} + \dot{z}\hat{k}. \tag{2}
$$

加速度再对速度求导：

$$
\mathbf{a}(t) = \frac{d\mathbf{v}}{dt} = \ddot{x}\hat{i} + \ddot{y}\hat{j} + \ddot{z}\hat{k}. \tag{3}
$$

点号表示对时间求导，这是力学中的标准简写：$\dot{x} = \frac{dx}{dt}$，$\ddot{x} = \frac{d^2x}{dt^2}$。

**优点**：基向量固定，求导简单直接，不会出现任何额外项。**缺点**：对于曲线运动，直角坐标往往需要耦合三个方向的运动方程，物理图像不够直观。例如圆周运动用 $x(t) = R\cos\omega t$，$y(t) = R\sin\omega t$ 描述，虽然正确，但没有直接体现“半径不变、角度匀速增加”的几何本质。

---

### 二、自然坐标系

自然坐标系（又称 Frenet-Serret 坐标系）是沿着轨迹本身建立的坐标架。它不依赖于固定的坐标轴，而是紧贴着质点的运动轨迹。这个坐标系在处理“已知轨迹形状、求运动参数”的问题时非常自然，比如过山车在曲线轨道上的受力分析，或者赛车在弯道中的加速度分析。

**建立方式**：

设质点沿一条已知曲线运动。我们用弧长 $s(t)$ 标记质点的位置——即从轨迹上的某个参考点开始，沿轨迹量取的距离。在轨迹上每一点，定义两个相互垂直的单位向量：

- 切向单位向量 $\hat{\tau}$：沿轨迹切线方向，指向 $s$ 增大的方向（即运动方向）。
- 法向单位向量 $\hat{n}$：沿轨迹法线方向，指向轨迹的凹侧（即曲率中心方向）。

这两个向量都随位置变化——质点运动到轨迹上的不同点时，$\hat{\tau}$ 和 $\hat{n}$ 的方向都会改变。

位置矢量虽然没有简单的分量形式，但弧长 $s$ 是直接可测的。速度总是沿切线方向：

$$
\mathbf{v} = \frac{d\mathbf{r}}{dt} = \frac{d\mathbf{r}}{ds}\frac{ds}{dt} = \frac{ds}{dt}\hat{\tau} = \dot{s}\hat{\tau} = v\hat{\tau}. \tag{4}
$$

其中 $v = \dot{s}$ 是瞬时速率。

加速度需要同时对 $\hat{\tau}$ 求导：

$$
\mathbf{a} = \frac{d\mathbf{v}}{dt} = \frac{d}{dt}(v\hat{\tau}) = \dot{v}\hat{\tau} + v\frac{d\hat{\tau}}{dt}. \tag{5}
$$

现在计算 $\frac{d\hat{\tau}}{dt}$。$\hat{\tau}$ 沿轨迹变化，其变化率可以通过弧长来连接：

$$
\frac{d\hat{\tau}}{dt} = \frac{d\hat{\tau}}{ds}\frac{ds}{dt} = v\frac{d\hat{\tau}}{ds}. \tag{6}
$$

在微分几何中，$\frac{d\hat{\tau}}{ds} = \frac{1}{\rho}\hat{n}$，其中 $\rho$ 是该点的曲率半径。这个关系可以通过考察单位切向量在无穷小弧长上的转动角度来证明。于是：

$$
\frac{d\hat{\tau}}{dt} = \frac{v}{\rho}\hat{n}. \tag{7}
$$

代入(5)：

$$
\mathbf{a} = \dot{v}\hat{\tau} + \frac{v^2}{\rho}\hat{n}. \tag{8}
$$

这就是自然坐标系中的加速度表达式。两项都有明确的物理意义：
- 切向加速度 $\dot{v}$：改变速率的大小。
- 法向加速度 $v^2/\rho$：改变速度的方向，指向曲率中心。

对于匀速圆周运动，$v$ 为常数，$\dot{v}=0$，只有法向加速度 $v^2/R$，指向圆心。这是我们都熟悉的结果。

---

### 三、极坐标系

极坐标系是处理绕固定点运动（如天体运动、圆周运动）最自然的坐标系。

**建立方式**：在平面内取一个固定点 $O$（极点）和一条固定射线（极轴）。质点的位置由两个量决定：到极点的距离 $r$（径向坐标）和与极轴的夹角 $\theta$（角坐标）。在每一点定义两个单位向量：

- 径向单位向量 $\hat{r}$：沿径向向外，即从极点 $O$ 指向质点 $P$ 的方向。
- 横向单位向量 $\hat{\theta}$：沿角度增大的方向，即逆时针方向，垂直于 $\hat{r}$。

位置矢量用极坐标表达非常简单：

$$
\mathbf{r} = r\hat{r}. \tag{9}
$$

但关键点是：** $\hat{r}$ 和 $\hat{\theta}$ 的方向依赖于 $\theta$，而 $\theta$ 随时间变化。因此这两个基向量也随时间变化！** 这是极坐标与直角坐标最本质的区别，也是初学者最容易犯错的地方。

**基向量随时间的变化**：

把 $\hat{r}$ 和 $\hat{\theta}$ 用直角坐标的基向量 $\hat{i}, \hat{j}$ 表示：

$$
\hat{r} = \cos\theta\,\hat{i} + \sin\theta\,\hat{j}, \tag{10}
$$

$$
\hat{\theta} = -\sin\theta\,\hat{i} + \cos\theta\,\hat{j}. \tag{11}
$$

验证：当 $\theta = 0$ 时，$\hat{r} = \hat{i}$，$\hat{\theta} = \hat{j}$，符合直觉——此时径向沿 $x$ 轴，横向沿 $y$ 轴。当 $\theta = \pi/2$ 时，$\hat{r} = \hat{j}$，$\hat{\theta} = -\hat{i}$，也正确。

现在对时间求导。注意 $\hat{i}, \hat{j}$ 是固定不动的，唯一变化的是 $\cos\theta$ 和 $\sin\theta$ 中的 $\theta(t)$。

对 $\hat{r}$ 求导：

$$
\begin{aligned}
\frac{d\hat{r}}{dt} &= \frac{d}{dt}(\cos\theta\,\hat{i} + \sin\theta\,\hat{j}) \\
&= (-\sin\theta\,\dot{\theta})\hat{i} + (\cos\theta\,\dot{\theta})\hat{j} \\
&= \dot{\theta}(-\sin\theta\,\hat{i} + \cos\theta\,\hat{j}).
\end{aligned}
$$

括号内正好是 $\hat{\theta}$ 的表达式(11)。因此：

$$
\dot{\hat{r}} = \dot{\theta}\,\hat{\theta}. \tag{12}
$$

对 $\hat{\theta}$ 求导：

$$
\begin{aligned}
\frac{d\hat{\theta}}{dt} &= \frac{d}{dt}(-\sin\theta\,\hat{i} + \cos\theta\,\hat{j}) \\
&= (-\cos\theta\,\dot{\theta})\hat{i} + (-\sin\theta\,\dot{\theta})\hat{j} \\
&= -\dot{\theta}(\cos\theta\,\hat{i} + \sin\theta\,\hat{j}).
\end{aligned}
$$

括号内正好是 $\hat{r}$ 的表达式(10)。因此：

$$
\dot{\hat{\theta}} = -\dot{\theta}\,\hat{r}. \tag{13}
$$

这两个公式(12)和(13)是极坐标中一切运动学推导的基石。

**速度的推导**：

从位置矢量(9)出发，用乘法法则求导：

$$
\mathbf{v} = \frac{d\mathbf{r}}{dt} = \frac{d}{dt}(r\hat{r}) = \dot{r}\hat{r} + r\dot{\hat{r}}.
$$

代入(12)：

$$
\mathbf{v} = \dot{r}\hat{r} + r\dot{\theta}\hat{\theta}. \tag{14}
$$

两项的物理意义极其清晰：
- 径向速度：$\dot{r}\hat{r}$，沿径向向外或向内的运动。
- 横向速度：$r\dot{\theta}\hat{\theta}$，绕极点转动的速度，其大小等于角速度乘以半径。

如果 $r$ 不变（圆周运动），径向速度为零，速度只有横向分量 $v = r\dot{\theta}$，这正是线速度与角速度的关系。

**加速度的推导**：

对速度(14)继续求导：

$$
\mathbf{a} = \frac{d\mathbf{v}}{dt} = \frac{d}{dt}(\dot{r}\hat{r} + r\dot{\theta}\hat{\theta}).
$$

逐项展开。第一项：

$$
\frac{d}{dt}(\dot{r}\hat{r}) = \ddot{r}\hat{r} + \dot{r}\dot{\hat{r}} = \ddot{r}\hat{r} + \dot{r}\dot{\theta}\hat{\theta}.
$$

第二项用乘法法则展开，注意 $r$、$\dot{\theta}$、$\hat{\theta}$ 三个量都随时间变化：

$$
\frac{d}{dt}(r\dot{\theta}\hat{\theta}) = \dot{r}\dot{\theta}\hat{\theta} + r\ddot{\theta}\hat{\theta} + r\dot{\theta}\dot{\hat{\theta}}.
$$

代入(13)中的 $\dot{\hat{\theta}} = -\dot{\theta}\hat{r}$：

$$
= \dot{r}\dot{\theta}\hat{\theta} + r\ddot{\theta}\hat{\theta} + r\dot{\theta}(-\dot{\theta}\hat{r}) = \dot{r}\dot{\theta}\hat{\theta} + r\ddot{\theta}\hat{\theta} - r\dot{\theta}^2\hat{r}.
$$

把两项加起来，$\hat{r}$ 方向的系数：$\ddot{r} - r\dot{\theta}^2$。$\hat{\theta}$ 方向的系数：$\dot{r}\dot{\theta} + \dot{r}\dot{\theta} + r\ddot{\theta} = 2\dot{r}\dot{\theta} + r\ddot{\theta}$。

因此：

$$
\mathbf{a} = (\ddot{r} - r\dot{\theta}^2)\hat{r} + (r\ddot{\theta} + 2\dot{r}\dot{\theta})\hat{\theta}. \tag{15}
$$

这是一般平面曲线运动在极坐标中的完整加速度表达式。我们将在后续天体力学中反复使用它。

**特殊情况**：圆周运动 $r = R$ 为常数。则 $\dot{r} = 0$，$\ddot{r} = 0$。代入(15)：

$$
\mathbf{a} = -R\dot{\theta}^2\hat{r} + R\ddot{\theta}\hat{\theta}. \tag{16}
$$

如果角速度恒定 $\dot{\theta} = \omega$，$\ddot{\theta} = 0$，则：

$$
\mathbf{a} = -R\omega^2\hat{r}. \tag{17}
$$

负号表示加速度指向 $- \hat{r}$ 方向，即指向圆心。这正是向心加速度 $a = R\omega^2 = v^2/R$。

---

### 四、柱坐标与球坐标简介

**柱坐标**：在极坐标基础上增加一个竖直方向 $z$。位置矢量为：

$$
\mathbf{r} = \rho\hat{\rho} + z\hat{k}. \tag{18}
$$

其中 $\rho$ 是到 $z$ 轴的距离，$\hat{\rho}$ 是径向基向量（与极坐标的 $\hat{r}$ 相同），$z$ 是高度，$\hat{k}$ 是竖直方向的固定单位向量。

速度为：

$$
\mathbf{v} = \dot{\rho}\hat{\rho} + \rho\dot{\phi}\hat{\phi} + \dot{z}\hat{k}. \tag{19}
$$

加速度为：

$$
\mathbf{a} = (\ddot{\rho} - \rho\dot{\phi}^2)\hat{\rho} + (\rho\ddot{\phi} + 2\dot{\rho}\dot{\phi})\hat{\phi} + \ddot{z}\hat{k}. \tag{20}
$$

柱坐标适用于圆柱对称问题，如带电圆柱周围的电场、圆管中的流动。

**球坐标**：用 $(r, \theta, \phi)$ 三个量描述空间点，其中 $r$ 是到原点的距离，$\theta$ 是极角（与 $z$ 轴的夹角），$\phi$ 是方位角（在 $xy$ 平面上的角度）。球坐标的基向量为 $\hat{r}, \hat{\theta}, \hat{\phi}$，位置矢量：

$$
\mathbf{r} = r\hat{r}. \tag{21}
$$

速度在球坐标中的表达式较复杂：

$$
\mathbf{v} = \dot{r}\hat{r} + r\dot{\theta}\hat{\theta} + r\sin\theta\,\dot{\phi}\hat{\phi}. \tag{22}
$$

加速度更为复杂，在后续课程（如电动力学中的偶极辐射、量子力学中的氢原子）中会用到，这里仅作了解。

---

### 五、总结

把 M.3 节的全部推导压缩成几条线：

1. **直角坐标系**：基向量固定，求导最简单，但曲线运动需要耦合三个坐标，物理图像不够直观。

2. **自然坐标系**：沿轨迹的切向和法向建立坐标。速度 $v = \dot{s}\hat{\tau}$，加速度 $\mathbf{a} = \dot{v}\hat{\tau} + \frac{v^2}{\rho}\hat{n}$。特别适合已知轨迹形状的运动分析。

3. **极坐标的核心推导**：
   - 基向量随时间变化：$\dot{\hat{r}} = \dot{\theta}\hat{\theta}$，$\dot{\hat{\theta}} = -\dot{\theta}\hat{r}$
   - 速度：$\mathbf{v} = \dot{r}\hat{r} + r\dot{\theta}\hat{\theta}$
   - 加速度：$\mathbf{a} = (\ddot{r} - r\dot{\theta}^2)\hat{r} + (r\ddot{\theta} + 2\dot{r}\dot{\theta})\hat{\theta}$

4. **关键区别**：直角坐标系中基向量固定，极坐标系中基向量随 $\theta$ 旋转。这个区别导致了极坐标加速度中出现了 $-r\dot{\theta}^2$（向心项）和 $2\dot{r}\dot{\theta}$（科里奥利项）等额外项。这些“额外项”不是数学的繁冗，而是物理的真实——它们对应着曲线运动中速度方向改变的几何效应。

5. **选择坐标系的原则**：尽量让约束条件在坐标中表现为常数的坐标系。圆周运动选极坐标（$r$ 固定），抛体运动选直角坐标（$a_x=0, a_y=-g$ 最简单），已知轨道形状选自然坐标系。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [M.4 常微分方程初步](https://zhuanlan.zhihu.com/p/2062423944645415902)

**课程目录：** [力学目录](https://zhuanlan.zhihu.com/p/2062487053170906945)
<!-- zhihu-nav-bottom:end -->
