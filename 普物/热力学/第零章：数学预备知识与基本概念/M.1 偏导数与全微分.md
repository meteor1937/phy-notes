---
created: 2026-07-20T15:24:01+08:00
modified: 2026-07-20T15:32:19+08:00
zhihu-title: M.1 偏导数与全微分
zhihu-topics:
  - 热力学
zhihu-link: https://zhuanlan.zhihu.com/p/2062559913633698674
zhihu-toc: true
zhihu-created-at: 2026-07-20 16:36
zhihu-updated-at: 2026-07-21 14:28
---
## M.1 偏导数与全微分

在进入热力学之前，我们必须先掌握一门语言——偏导数和全微分。热力学处理的物理量通常是多个变量的函数：压强 $p$ 依赖于温度 $T$ 和体积 $V$，内能 $U$ 依赖于温度 $T$ 和体积 $V$，熵 $S$ 依赖于温度 $T$ 和压强 $p$。当我们改变一个变量而保持另一个不变时，物理量如何变化？这正是偏导数回答的问题。

更重要的是，热力学中有一类量（如内能、熵）只依赖于系统的状态，与如何到达该状态无关——它们叫做**状态函数**。另一类量（如功、热量）依赖于具体的过程路径——它们叫做**过程量**。偏导数和全微分提供了精确区分这两类量的数学工具：状态函数的微分是**恰当微分**，过程量的微分是**非恰当微分**。

这一节我们从偏导数的定义出发，建立全微分、链式法则、循环关系等一系列工具。这些工具在热力学的后续推导中会反复出现——从麦克斯韦关系到克拉佩龙方程，无一例外。

---

### 一、偏导数的定义与物理含义

**偏导数**的定义与普通导数完全相同，只是函数依赖于多个变量，我们在求导时只让一个变量变化，其余变量保持固定。

设 $f$ 是 $x$ 和 $y$ 的二元函数：$f = f(x, y)$。$f$ 对 $x$ 的偏导数定义为：

$$
\left( \frac{\partial f}{\partial x} \right)_y = \lim_{\Delta x \to 0} \frac{f(x+\Delta x, y) - f(x, y)}{\Delta x}. \tag{1}
$$

下标 $y$ 的含义是：在取极限的过程中，$y$ 保持固定不变。同理，$f$ 对 $y$ 的偏导数为：

$$
\left( \frac{\partial f}{\partial y} \right)_x = \lim_{\Delta y \to 0} \frac{f(x, y+\Delta y) - f(x, y)}{\Delta y}. \tag{2}
$$

**下标的重要性**：在热力学中，变量的选取不是唯一的。同一个物理量可以用不同的自变量组合来表示。例如，内能 $U$ 可以看成 $U(T,V)$，也可以看成 $U(S,V)$，还可以看成 $U(T,p)$。在不同的变量组合下，同一个物理量对同一个自变量的偏导数可能不同。因此，下标必须明确写出，指明哪些变量在求导过程中保持不变——否则表达式就没有确定的意义。

**物理含义**：偏导数 $\left( \frac{\partial f}{\partial x} \right)_y$ 是在 $y$ 不变的条件下，$f$ 随 $x$ 的变化率。在热力学中，这对应着特定的物理过程：
- $\left( \frac{\partial U}{\partial T} \right)_V$：体积不变时内能随温度的变化率（定容热容）
- $\left( \frac{\partial U}{\partial T} \right)_p$：压强不变时内能随温度的变化率（不同于定容热容）
- $\left( \frac{\partial p}{\partial T} \right)_V$：体积不变时压强随温度的变化率

这两个偏导数不同，因为一个是 $V$ 固定，一个是 $p$ 固定。下标正是为了区分这种情况。

---

### 二、全微分

如果函数 $f = f(x, y)$ 可微，那么当 $x$ 和 $y$ 都发生微小变化时，$f$ 的总变化量可以写成：

$$
df = \left( \frac{\partial f}{\partial x} \right)_y dx + \left( \frac{\partial f}{\partial y} \right)_x dy. \tag{3}
$$

这就是**全微分**。它告诉我们：$f$ 的微小变化等于每个自变量单独变化引起的贡献之和。

**几何意义**：在三维空间中，$f = f(x, y)$ 是一个曲面。全微分 $df$ 是曲面上某一点的切平面的高度变化。

**推广到更多变量**：如果 $f = f(x_1, x_2, \ldots, x_n)$，则：

$$
df = \sum_{i=1}^n \left( \frac{\partial f}{\partial x_i} \right)_{x_j} dx_i. \tag{4}
$$

其中下标 $x_j$ 表示除了 $x_i$ 之外的所有变量保持固定。

**热力学中的一个关键区分**：全微分 $df$ 是一个**无穷小量**，它描述的是 $f$ 在状态空间中从 $(x,y)$ 到 $(x+dx, y+dy)$ 的变化。如果 $f$ 是状态函数（如内能、熵），那么 $df$ 是一个恰当微分——它沿着任何路径的积分只取决于起点和终点。如果 $f$ 不是状态函数（如功、热量），那么 $\delta f$ 不是恰当微分——它的积分依赖于路径。我们将在下一节详细区分。

**计算示例**：设 $f(x,y) = x^2 y + y^3$。

$$
\left( \frac{\partial f}{\partial x} \right)_y = 2xy, \quad \left( \frac{\partial f}{\partial y} \right)_x = x^2 + 3y^2.
$$

全微分为：

$$
df = 2xy \, dx + (x^2 + 3y^2) \, dy.
$$

---

### 三、恰当微分与非恰当微分——状态函数与过程量

这是偏微分在热力学中最重要的应用。我们需要严格区分两种微分。

**恰当微分**：如果存在一个函数 $f(x, y)$，使得微分 $df$ 可以写成(3)的形式，则称 $df$ 是恰当微分。物理上，$f$ 对应一个**状态函数**——它的值只依赖于系统的状态，与如何到达该状态无关。

**非恰当微分**：如果不存在这样的函数 $f$，则称该微分是非恰当微分。物理上，它对应一个**过程量**——它的值依赖于具体的路径，不是状态的函数。

**判据**：对于二元函数 $f(x, y)$，微分形式

$$
df = P(x, y) \, dx + Q(x, y) \, dy \tag{5}
$$

是恰当微分的**充要条件**是：

$$
\left( \frac{\partial P}{\partial y} \right)_x = \left( \frac{\partial Q}{\partial x} \right)_y. \tag{6}
$$

这个条件来自混合偏导数的可交换性，我们将在下一节证明。

**物理例子**：
- 内能 $U$：$dU$ 是恰当微分——内能是状态函数。
- 熵 $S$：$dS = \delta Q_{\text{rev}}/T$ 是恰当微分——熵是状态函数。
- 功 $\delta W = p\,dV$：非恰当微分——功是过程量，依赖路径。
- 热量 $\delta Q$：非恰当微分——热量是过程量，依赖路径。

**记号约定**：在热力学中，我们通常用 $d$ 表示恰当微分（状态函数），用 $\delta$ 表示非恰当微分（过程量）。这个记号区分不是数学上的必然，而是物理上的惯例，用来提醒我们这两类量本质不同。

---

### 四、混合偏导数的可交换性

如果一个函数的二阶偏导数连续，则混合偏导数的求导次序可以交换：

$$
\frac{\partial^2 f}{\partial x \partial y} = \frac{\partial^2 f}{\partial y \partial x}. \tag{7}
$$

这是微积分中的一个基本定理。它在热力学中的核心应用是：从全微分 $df = P\,dx + Q\,dy$ 出发，如果 $df$ 是恰当微分（即存在 $f$ 使得 $P = (\partial f/\partial x)_y$，$Q = (\partial f/\partial y)_x$），那么：

$$
\left( \frac{\partial P}{\partial y} \right)_x = \frac{\partial^2 f}{\partial y \partial x}, \quad \left( \frac{\partial Q}{\partial x} \right)_y = \frac{\partial^2 f}{\partial x \partial y}.
$$

由混合偏导数的可交换性，两者相等，从而得到恰当微分的判据(6)。

**热力学中的典型应用**：麦克斯韦关系就是混合偏导数可交换性在热力学势上的直接应用。例如，从 $dF = -S\,dT - p\,dV$ 出发，$P = -S$，$Q = -p$，恰当微分条件给出：

$$
\left( \frac{\partial S}{\partial V} \right)_T = \left( \frac{\partial p}{\partial T} \right)_V.
$$

这是四个麦克斯韦关系之一。没有混合偏导数的可交换性，就没有麦克斯韦关系。

---

### 五、链式法则——热力学中变量之间的隐式关系

在热力学中，物理量之间的关系通常是隐式的。例如，状态方程 $f(p, V, T) = 0$ 给出了 $p, V, T$ 之间的约束关系——它们不是独立的，只有两个独立变量。我们需要工具来处理这种隐式关系下的导数关系。

**链式法则的标准形式**：如果 $f = f(x, y)$，而 $x = x(t)$，$y = y(t)$，则：

$$
\frac{df}{dt} = \left( \frac{\partial f}{\partial x} \right)_y \frac{dx}{dt} + \left( \frac{\partial f}{\partial y} \right)_x \frac{dy}{dt}. \tag{8}
$$

这是基本的链式法则。热力学中的更复杂情况是：三个变量 $x, y, z$ 之间存在一个约束关系，只有两个是独立的。这时我们需要推导变量之间的偏导数关系。

**三个变量之间的偏导数关系**：

设 $x, y, z$ 满足一个隐式关系 $F(x, y, z) = 0$，因此只有两个独立变量。固定 $z$ 时，$x$ 和 $y$ 之间的变化率由 $\left( \frac{\partial x}{\partial y} \right)_z$ 描述。我们可以通过全微分来推导这些量之间的关系。

对 $F(x, y, z) = 0$ 取全微分：

$$
\left( \frac{\partial F}{\partial x} \right)_{y,z} dx + \left( \frac{\partial F}{\partial y} \right)_{x,z} dy + \left( \frac{\partial F}{\partial z} \right)_{x,y} dz = 0. \tag{9}
$$

固定 $z$（$dz = 0$），可得：

$$
\left( \frac{\partial F}{\partial x} \right)_{y,z} dx + \left( \frac{\partial F}{\partial y} \right)_{x,z} dy = 0.
$$

所以：

$$
\left( \frac{\partial x}{\partial y} \right)_z = - \frac{\left( \frac{\partial F}{\partial y} \right)_{x,z}}{\left( \frac{\partial F}{\partial x} \right)_{y,z}}. \tag{10}
$$

同理，固定 $x$ 或 $y$ 可以得到其他偏导数关系。这是处理隐式约束时求偏导的标准方法。

---

### 六、循环关系

这是热力学变量变换中最常用、也最容易被忽略的关系。它告诉我们：如果三个变量 $x, y, z$ 之间存在一个函数关系，使得其中任意两个可以确定第三个，那么：

$$
\boxed{ \left( \frac{\partial x}{\partial y} \right)_z \left( \frac{\partial y}{\partial z} \right)_x \left( \frac{\partial z}{\partial x} \right)_y = -1 }. \tag{11}
$$

这个关系看起来反直觉——三个偏导数相乘为什么是 $-1$，而不是 $+1$？我们现在完整推导它。

**推导**：设 $x = x(y, z)$，$y = y(x, z)$，$z = z(x, y)$。这个关系的前提是三个变量中只有两个独立——它们之间存在一个约束。

首先，写出 $dx$ 的全微分，把 $x$ 看成 $y$ 和 $z$ 的函数：

$$
dx = \left( \frac{\partial x}{\partial y} \right)_z dy + \left( \frac{\partial x}{\partial z} \right)_y dz. \tag{12}
$$

现在考虑固定 $x$（$dx = 0$）的情况。在固定 $x$ 的条件下，$dy$ 和 $dz$ 不再是独立的，它们之间有关系。由(12)：

$$
0 = \left( \frac{\partial x}{\partial y} \right)_z dy + \left( \frac{\partial x}{\partial z} \right)_y dz.
$$

解得：

$$
\frac{dy}{dz} \Big|_{x} = - \frac{\left( \frac{\partial x}{\partial z} \right)_y}{\left( \frac{\partial x}{\partial y} \right)_z}. \tag{13}
$$

但左边正是 $\left( \frac{\partial y}{\partial z} \right)_x$。因此：

$$
\left( \frac{\partial y}{\partial z} \right)_x = - \frac{\left( \frac{\partial x}{\partial z} \right)_y}{\left( \frac{\partial x}{\partial y} \right)_z}. \tag{14}
$$

或者等价地：

$$
\left( \frac{\partial y}{\partial z} \right)_x \left( \frac{\partial x}{\partial y} \right)_z = - \left( \frac{\partial x}{\partial z} \right)_y. \tag{15}
$$

现在，$\left( \frac{\partial x}{\partial z} \right)_y$ 的倒数就是 $\left( \frac{\partial z}{\partial x} \right)_y$（因为固定 $y$ 时，$x$ 和 $z$ 之间的函数关系是反函数关系）。因此：

$$
\left( \frac{\partial y}{\partial z} \right)_x \left( \frac{\partial x}{\partial y} \right)_z \left( \frac{\partial z}{\partial x} \right)_y = -1. \tag{16}
$$

这就证明了循环关系。

**一个具体的数值验证**：

取 $pV = nRT$（理想气体状态方程），三个变量 $p, V, T$ 之间存在约束关系。固定 $T$ 时 $p \propto 1/V$，所以 $\left( \frac{\partial p}{\partial V} \right)_T = -\frac{nRT}{V^2}$。

固定 $p$ 时 $V \propto T$，所以 $\left( \frac{\partial V}{\partial T} \right)_p = \frac{nR}{p}$。

固定 $V$ 时 $p \propto T$，所以 $\left( \frac{\partial T}{\partial p} \right)_V = \frac{V}{nR}$。

三者相乘：

$$
\left( \frac{\partial p}{\partial V} \right)_T \left( \frac{\partial V}{\partial T} \right)_p \left( \frac{\partial T}{\partial p} \right)_V = \left( -\frac{nRT}{V^2} \right) \left( \frac{nR}{p} \right) \left( \frac{V}{nR} \right) = -\frac{nRT}{pV} = -1.
$$

循环关系成立。

**物理含义**：循环关系告诉我们，三个热力学变量之间的偏导数不能随意选择——它们必须满足这个乘积为 $-1$ 的约束。这意味着三个变量中只有两个是真正独立的，第三个由状态方程确定。

---

### 七、雅可比行列式——热力学变量变换的系统工具

在热力学中，我们经常需要从一个变量组合 $(x, y)$ 变换到另一个变量组合 $(u, v)$。雅可比行列式提供了系统地进行这种变换的工具。

**定义**：雅可比行列式（Jacobian）定义为：

$$
\frac{\partial (x, y)}{\partial (u, v)} \equiv
\begin{vmatrix}
\left( \frac{\partial x}{\partial u} \right)_v & \left( \frac{\partial x}{\partial v} \right)_u \\[6pt]
\left( \frac{\partial y}{\partial u} \right)_v & \left( \frac{\partial y}{\partial v} \right)_u
\end{vmatrix}
= \left( \frac{\partial x}{\partial u} \right)_v \left( \frac{\partial y}{\partial v} \right)_u - \left( \frac{\partial x}{\partial v} \right)_u \left( \frac{\partial y}{\partial u} \right)_v. \tag{17}
$$

**性质**：
1. 反函数雅可比：$\frac{\partial (x, y)}{\partial (u, v)} \cdot \frac{\partial (u, v)}{\partial (x, y)} = 1$
2. 链式法则：$\frac{\partial (x, y)}{\partial (u, v)} = \frac{\partial (x, y)}{\partial (r, s)} \cdot \frac{\partial (r, s)}{\partial (u, v)}$
3. 交换性：$\frac{\partial (x, y)}{\partial (u, v)} = -\frac{\partial (y, x)}{\partial (u, v)}$

**热力学中的应用**：

当我们想把一个偏导数从一种变量组合转换到另一种时，雅可比行列式提供了标准方法。例如：

$$
\left( \frac{\partial U}{\partial V} \right)_T = \frac{\partial (U, T)}{\partial (V, T)} = \frac{\partial (U, T)}{\partial (p, T)} \cdot \frac{\partial (p, T)}{\partial (V, T)}.
$$

使用雅可比行列式，我们可以系统地处理这类变换，而不需要每次都重新推导链式法则。

**与循环关系的联系**：循环关系实际上是雅可比行列式的一个特例：

$$
\left( \frac{\partial x}{\partial y} \right)_z \left( \frac{\partial y}{\partial z} \right)_x \left( \frac{\partial z}{\partial x} \right)_y = -1.
$$

这个关系可以证明 $\frac{\partial (x, y)}{\partial (y, z)} \cdot \frac{\partial (y, z)}{\partial (z, x)} \cdot \frac{\partial (z, x)}{\partial (x, y)} = -1$，直接来自雅可比行列式的链式法则。

---

### 八、总结

把 M.1 节的全部内容压缩成几条线：

1. **偏导数**：$\left( \frac{\partial f}{\partial x} \right)_y$ 表示在 $y$ 保持不变的条件下 $f$ 对 $x$ 的变化率。下标必须明确写出——热力学中同一物理量在不同变量组合下有不同的偏导数。

2. **全微分**：$df = \left( \frac{\partial f}{\partial x} \right)_y dx + \left( \frac{\partial f}{\partial y} \right)_x dy$。它把函数的总变化分解为各独立变量变化的贡献之和。

3. **恰当微分与非恰当微分**：
   - 恰当微分：存在状态函数 $f$ 使得 $df = Pdx + Qdy$，满足 $\left( \frac{\partial P}{\partial y} \right)_x = \left( \frac{\partial Q}{\partial x} \right)_y$。对应状态函数。
   - 非恰当微分：不存在这样的 $f$，积分依赖路径。对应过程量（功、热量）。

4. **混合偏导数可交换性**：$\frac{\partial^2 f}{\partial x \partial y} = \frac{\partial^2 f}{\partial y \partial x}$。这是麦克斯韦关系的数学基础。

5. **链式法则**：$\frac{df}{dt} = \left( \frac{\partial f}{\partial x} \right)_y \frac{dx}{dt} + \left( \frac{\partial f}{\partial y} \right)_x \frac{dy}{dt}$。用于处理变量随参数变化的情况。

6. **循环关系**：$\left( \frac{\partial x}{\partial y} \right)_z \left( \frac{\partial y}{\partial z} \right)_x \left( \frac{\partial z}{\partial x} \right)_y = -1$。三变量约束下的偏导数乘积关系，推导中使用了固定一个变量时的全微分。

7. **雅可比行列式**：$\frac{\partial (x, y)}{\partial (u, v)}$ 是热力学变量变换的系统工具，用于在变量组合之间转换偏导数。

这些偏微分工具是热力学的数学语言。在接下来的章节中，它们将被反复使用——从热力学第一定律的微分形式 $dU = T\,dS - p\,dV$，到麦克斯韦关系的导出，再到克拉佩龙方程的推导，无一例外。掌握这些工具，就是掌握了阅读和书写热力学的语法。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [M.2 理想气体状态方程](https://zhuanlan.zhihu.com/p/2062559974912528631)

**课程目录：** [热力学目录](https://zhuanlan.zhihu.com/p/2062556088797500605)
<!-- zhihu-nav-bottom:end -->
