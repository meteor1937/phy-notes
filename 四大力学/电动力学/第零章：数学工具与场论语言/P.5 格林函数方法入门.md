---
created: 2026-07-17T00:59:45+08:00
modified: 2026-07-17T01:26:47+08:00
zhihu-title: P.5 格林函数方法入门
zhihu-topics: 电动力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061253748400337330
zhihu-toc: true
zhihu-created-at: 2026-07-17 01:27
review-status: reviewing
physics-checked: true
derivation-checked: true
readability-checked: true
teaching-checked: true
zhihu-layout-checked: false
last-reviewed:
zhihu-updated-at: 2026-07-21 14:52
---
<!-- zhihu-nav-top:start -->
**上一节：** [P.4 狄拉克 δ 函数](https://zhuanlan.zhihu.com/p/2061253558738097791)
<!-- zhihu-nav-top:end -->

# P.5 格林函数方法入门

前面几节我们建立了矢量分析和 $\delta$ 函数的工具。现在把这些工具组合起来，解决电磁学中最核心的数学问题：给定电荷分布，求电势。

这个问题在数学上就是解泊松方程 $\nabla^2\phi = -\rho/\epsilon_0$。如果电荷分布是简单的——比如点电荷、均匀带电球——我们可以用高斯定律或者直接积分来求电势。但如果电荷分布是任意的、没有对称性的，该怎么系统地求解？

格林函数方法的核心思路是：先求线性算符对单位点源的响应，再把一般源写成点源的连续叠加。方法本身不要求源或区域具有对称性，但 Green 函数会依赖区域形状和边界条件；复杂几何中的困难往往正是构造这个核函数。

这一节以 Poisson 方程为例，建立从算符约定、点源归一化、Green 第二恒等式到边界积分表示的完整逻辑链。

## 学习目标

完成本篇后，读者应当能够：

1. 在固定符号约定下定义 Poisson 算符及其 Green 函数；
2. 验证 $1/(4\pi R)$ 的点源归一化，并从卷积复算 Poisson 方程；
3. 用 Green 第二恒等式推导含边界项的积分表示；
4. 区分 Dirichlet 与 Neumann Green 函数，并说明 Neumann 问题的兼容条件与常数零模；
5. 用点电荷极限、量纲和远场行为检查所得电势。

## 先修知识、约定与路线

需要掌握 [[P.1 矢量分析]]、[[P.4 狄拉克 δ 函数]] 中的散度定理和分布恒等式。本文采用 SI 制。场点记为 $\mathbf r$，源点和积分变量记为 $\mathbf r'$，并定义

$$
R=|\mathbf r-\mathbf r'|,
\qquad
d^3r'=d\tau'.
$$

全文统一使用正算符

$$
L=-\nabla^2,
\qquad
L\phi=\frac{\rho}{\epsilon_0},
\qquad
L_{\mathbf r}G(\mathbf r,\mathbf r')=\delta^{(3)}(\mathbf r-\mathbf r').
$$

等价地，$\nabla^2G=-\delta^{(3)}$。推导先在全空间得到自由空间核，再用 Green 第二恒等式恢复边界项，最后讨论有限区域的 Dirichlet 与 Neumann 条件。假定源足够局域、积分收敛；涉及点源时按分布理解等式。

---

## 一、线性系统的叠加原理与点源响应

格林函数方法的物理基础是叠加原理。静电学是线性的：如果电荷分布 $\rho_1$ 产生电势 $\phi_1$，电荷分布 $\rho_2$ 产生电势 $\phi_2$，那么 $\rho_1 + \rho_2$ 产生的电势就是 $\phi_1 + \phi_2$。数学上，这是因为泊松方程 $\nabla^2\phi = -\rho/\epsilon_0$ 是线性的——拉普拉斯算符是线性算符，不包含 $\phi$ 的非线性项。

叠加原理意味着，如果我们能把任意电荷分布分解为一些“基本单元”的叠加，那么总电势就是这些基本单元各自产生的电势的叠加。那么，最自然的基本单元是什么？是点电荷。

一个位于 $\mathbf{r}'$、电荷量为 $q$ 的点电荷，其电荷密度为 $\rho(\mathbf{r}) = q\,\delta^{(3)}(\mathbf{r} - \mathbf{r}')$。任意电荷分布 $\rho(\mathbf{r})$ 可以看作是连续分布的点电荷的叠加——每个空间点 $\mathbf{r}'$ 上的电荷量 $\rho(\mathbf{r}') d\tau'$ 就是一个“微小的点电荷”，位于 $\mathbf{r}'$。那么，总电势就应该是所有 $\mathbf{r}'$ 上这些微小点电荷产生的电势的叠加。

一个位于 $\mathbf{r}'$ 的单位点电荷（$q=1$）产生的电势，我们在物理上已经知道是 $1/(4\pi\epsilon_0|\mathbf{r} - \mathbf{r}'|)$。因此，位于 $\mathbf{r}'$、电荷量为 $\rho(\mathbf{r}') d\tau'$ 的微小点电荷产生的电势就是 $\frac{\rho(\mathbf{r}') d\tau'}{4\pi\epsilon_0|\mathbf{r} - \mathbf{r}'|}$。把所有 $\mathbf{r}'$ 上的贡献叠加起来（积分），就得到了总电势：

$$
\phi(\mathbf{r}) = \frac{1}{4\pi\epsilon_0} \int \frac{\rho(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} d\tau'
$$

这就是库仑积分的电势形式。这一步只是库仑定律与叠加原理给出的物理预期，还没有检查积分是否满足 Poisson 方程。严格地说，$1/(4\pi R)$ 是算符 $L=-\nabla^2$ 的 Green 函数；乘上 $1/\epsilon_0$ 后才是单位点电荷的电势核。下面逐项验证符号与归一化。


## 二、格林函数的严格定义

给定线性微分算符 $\mathcal L$、定义域和边界条件，Green 函数定义为

$$
\mathcal{L}_{\mathbf{r}} G(\mathbf{r}, \mathbf{r}') = \delta^{(3)}(\mathbf{r} - \mathbf{r}')
$$

这里 $\mathcal{L}_{\mathbf{r}}$ 表示算符作用于 $G$ 的第一个变量 $\mathbf{r}$，$\mathbf{r}'$ 是点源的位置参数。格林函数 $G(\mathbf{r}, \mathbf{r}')$ 的物理意义是：位于 $\mathbf{r}'$ 处的单位强度的点源，在观测点 $\mathbf{r}$ 处产生的场。

对于泊松方程，算符是 $\mathcal{L} = -\nabla^2$（注意这里的符号约定——泊松方程为 $-\nabla^2\phi = \rho/\epsilon_0$，所以格林函数满足 $-\nabla^2 G = \delta$ 或者等价地 $\nabla^2 G = -\delta$）。我们采用约定：

$$
\boxed{\nabla^2 G(\mathbf{r}, \mathbf{r}') = -\delta^{(3)}(\mathbf{r} - \mathbf{r}')}
$$

由 [[P.4 狄拉克 δ 函数]] 中的分布恒等式 $\nabla^2(1/R)=-4\pi\delta^{(3)}(\mathbf r-\mathbf r')$，可知自由空间中满足无穷远趋零条件的解为

$$
\boxed{G(\mathbf{r}, \mathbf{r}') = \frac{1}{4\pi|\mathbf{r} - \mathbf{r}'|}}
$$

归一化可由小球通量独立检查。在以 $\mathbf r'$ 为中心、半径为 $a$ 的球上，

$$
-\oint_{R=a}\nabla G\cdot d\mathbf a
=-\left(-\frac{1}{4\pi a^2}\right)4\pi a^2=1.
$$

因此 $-\nabla^2G$ 在原点的总源强确为 1。平移不变性使自由空间核只依赖 $R$；在有边界区域中，这一性质通常不再成立。


## 三、从格林函数到一般解：卷积

有了格林函数，任意源分布 $f(\mathbf{r})$ 下方程 $\mathcal{L}\phi(\mathbf{r}) = f(\mathbf{r})$ 的解可以通过叠加原理得到。将源 $f(\mathbf{r})$ 写为点源的连续叠加：

$$
f(\mathbf{r}) = \int \delta^{(3)}(\mathbf{r} - \mathbf{r}') f(\mathbf{r}') d\tau'
$$

对于每个点源 $\delta^{(3)}(\mathbf{r} - \mathbf{r}') f(\mathbf{r}') d\tau'$，其产生的响应是 $G(\mathbf{r}, \mathbf{r}') f(\mathbf{r}') d\tau'$。将所有点源的响应线性叠加，得到通解：

$$
\boxed{\phi(\mathbf{r}) = \int G(\mathbf{r}, \mathbf{r}') f(\mathbf{r}') d\tau'}
$$

只要 $G$ 与 $\phi$ 采用相同边界条件，并且没有未固定的齐次零模，这个积分就给出相应解。边界条件并未消失，而是被编码进 $G$；构造满足复杂边界的 $G$ 往往是主要难点。

对于本文约定的方程 $L\phi=\rho/\epsilon_0$，应取 $f=\rho/\epsilon_0$。代入 $G=1/(4\pi R)$ 得

$$
\begin{aligned}
\phi(\mathbf r)
&=\int G(\mathbf r,\mathbf r')\frac{\rho(\mathbf r')}{\epsilon_0}\,d^3r'\\
&=\frac{1}{4\pi\epsilon_0}
\int\frac{\rho(\mathbf r')}{|\mathbf r-\mathbf r'|}\,d^3r'.
\end{aligned}
$$

这里没有可任意“吸收”的负号：电势的参考自由度只允许加常数，不能把 $\phi$ 整体乘以 $-1$。直接让 $L=-\nabla^2$ 作用在积分上，得到

$$
L\phi
=\frac1{\epsilon_0}\int
\delta^{(3)}(\mathbf r-\mathbf r')\rho(\mathbf r')\,d^3r'
=\frac{\rho(\mathbf r)}{\epsilon_0},
$$

所以卷积的符号和归一化均正确。该全空间公式还假定电势在无穷远选为零；无限长、总电荷不收敛的源不满足此前提，必须改用与几何相配的边界条件和 Green 函数。


## 四、用 Green 第二恒等式恢复边界项

上面的推导利用了“点源叠加”的物理直觉。为了让这个推导在数学上更严密，我们可以借助格林恒等式（Green's identities）来展示格林函数和偏微分方程解的精确关系。

Green 第二恒等式是散度定理的直接推论。对于两个在体积 $V$ 内二次可微的标量函数 $\phi$ 和 $\psi$，有

$$
\int_V \left( \phi \nabla^2\psi - \psi \nabla^2\phi \right) d\tau = \oint_S \left( \phi \nabla\psi - \psi \nabla\phi \right) \cdot d\mathbf{a}
$$

对 Laplace 算符配合同一组齐次 Dirichlet 或适当归一化的 Neumann 条件，Green 第二恒等式还给出互易性 $G(\mathbf r,\mathbf r')=G(\mathbf r',\mathbf r)$。因此，定义中对场变量成立的 $-\nabla^2G=\delta^{(3)}$ 也可对源变量写成 $-\nabla'^2G=\delta^{(3)}$。

令 $\psi(\mathbf r')=G(\mathbf r,\mathbf r')$，其中 $\mathbf r$ 固定、$\mathbf r'$ 是积分变量。由于 $G$ 在 $\mathbf r'=\mathbf r$ 奇异，严格做法是先从 $V$ 中挖去一个以 $\mathbf r$ 为中心的小球，再令球半径趋零；这个小球面通量正好产生下文的 $-\phi(\mathbf r)$ 项。

将这两个函数代入格林恒等式（积分变量为 $\mathbf{r}'$，作用在 $\mathbf{r}'$ 上的拉普拉斯记作 $\nabla'^2$）：

$$
\begin{aligned}
\int_V\left(\phi\nabla'^2G-G\nabla'^2\phi\right)\,d^3r'
=\oint_S\left(\phi\nabla'G-G\nabla'\phi\right)\cdot d\mathbf a'
\end{aligned}
$$

代入 $\nabla'^2G=-\delta^{(3)}$ 与 $\nabla'^2\phi=-\rho/\epsilon_0$，左边成为

$$
-\phi(\mathbf r)
+\frac1{\epsilon_0}\int_VG(\mathbf r,\mathbf r')\rho(\mathbf r')\,d^3r'.
$$

移项后得到完整表示

$$
\boxed{
\begin{aligned}
\phi(\mathbf r)={}&\frac1{\epsilon_0}\int_VG\rho\,d^3r'\\
&-\oint_S\phi\frac{\partial G}{\partial n'}\,da'
+\oint_SG\frac{\partial\phi}{\partial n'}\,da'.
\end{aligned}}
$$

这里 $\partial/\partial n'$ 沿区域 $V$ 的外法向。全空间中若 $\phi$、$G$ 及其梯度使无穷远面积分趋零，两个边界项消失，公式退化为库仑积分。有限区域中不能丢掉它们。


## 五、Green 函数是线性系统的脉冲响应

从信号与系统的角度来看，格林函数方法就是线性系统理论在偏微分方程中的推广。一个线性系统由输入 $f$ 和输出 $\phi$ 来描述，$\phi = \mathcal{L}^{-1} f$。格林函数 $G$ 就是这个系统的“脉冲响应”——当输入是一个 $\delta$ 脉冲（点源）时，系统的输出就是 $G$。对于任意的输入 $f$，输出就是脉冲响应与输入的卷积：$\phi = G * f$。

在静电学中，输入是 $f=\rho/\epsilon_0$，输出是电势，算符的脉冲响应是 $G=1/(4\pi R)$。若直接把输入称作 $\rho$，则相应的电势核应写成 $G/\epsilon_0$，两种说法不能混用。

同一方法也适用于 Helmholtz、波动和扩散算符，但必须同时指定定义域、边界或初始条件。例如 $(\nabla^2+k^2)G=-\delta^{(3)}$ 的自由空间外行 Green 函数是 $e^{ikR}/(4\pi R)$；若不加 Sommerfeld 辐射条件，入行解或齐次解同样允许。只有在线性算符和边界条件都具有平移不变性时，$G$ 才只依赖 $\mathbf r-\mathbf r'$，积分才是通常意义的卷积。


## 六、例题：均匀带电实心球

**已知条件：** 半径为 $a$ 的球内电荷密度为常数 $\rho_0$，球外为零，并取 $\phi(\infty)=0$。

**待求量：** 用自由空间 Green 函数求球内外电势，并检查连续性与远场极限。

球对称积分所需的角向恒等式是

$$
\int\frac{d\Omega'}{|\mathbf r-\mathbf r'|}
=\frac{4\pi}{\max(r,r')}.
$$

它可由 $1/|\mathbf r-\mathbf r'|$ 的 Legendre 展开得到；角积分只保留 $\ell=0$ 项。对于 $r\ge a$，始终有 $r\ge r'$，所以

$$
\phi_{\mathrm{out}}(r)
=\frac{\rho_0}{\epsilon_0r}\int_0^a r'^2\,dr'
=\frac{\rho_0a^3}{3\epsilon_0r}
=\frac{Q}{4\pi\epsilon_0r},
$$

其中 $Q=4\pi a^3\rho_0/3$。对于 $r<a$，在 $r'=r$ 处分开径向积分：

$$
\begin{aligned}
\phi_{\mathrm{in}}(r)
&=\frac{\rho_0}{\epsilon_0}
\left[
\frac1r\int_0^r r'^2\,dr'
+\int_r^a r'\,dr'
\right]\\
&=\frac{\rho_0}{6\epsilon_0}(3a^2-r^2).
\end{aligned}
$$

在 $r=a$ 处两式都给出 $\rho_0a^2/(3\epsilon_0)$，且法向导数连续，因为表面没有额外面电荷。远处结果只依赖总电荷 $Q$，回到点电荷势；量纲 $\rho_0a^2/\epsilon_0$ 为伏特。这三项分别检查了边界匹配、远场极限和单位。


## 七、Dirichlet 与 Neumann 边界条件

全空间核 $G_0=1/(4\pi R)$ 对应无穷远趋零条件。有限区域的 Dirichlet 问题给定 $\phi|_S$，可选 $G_D|_S=0$，使完整表示中的最后一个边界项消失：

$$
\phi(\mathbf r)
=\frac1{\epsilon_0}\int_VG_D\rho\,d^3r'
-\oint_S\phi\frac{\partial G_D}{\partial n'}\,da'.
$$

Neumann 问题给定 $\partial\phi/\partial n$。常数满足 $L(1)=0$，所以解只能确定到加法常数；相应 Green 函数不能在有限区域同时满足 $-\nabla^2G_N=\delta^{(3)}$ 与齐次法向导数。对该方程在 $V$ 积分会得到左边通量为零、右边为 1 的矛盾。常用修正是

$$
-\nabla'^2G_N(\mathbf r,\mathbf r')
=\delta^{(3)}(\mathbf r-\mathbf r')-\frac1{|V|},
\qquad
\frac{\partial G_N}{\partial n'}\bigg|_S=0.
$$

原 Poisson 问题还必须满足由散度定理得到的兼容条件

$$
\oint_S\frac{\partial\phi}{\partial n}\,da
=-\frac1{\epsilon_0}\int_V\rho\,d^3r.
$$

指定一个参考值或令体积平均电势为零，才可消除常数零模。

镜像法可理解为构造 $G_D$ 的一种技巧：在自由空间核上加入关于源变量调和的修正，使合成核满足指定边界条件。镜像电荷位于物理区域之外，不是额外真实电荷。

---

## 八、常见误区

1. **把 $G$ 直接称为单位点电荷的电势。** 本文的 $G$ 是 $L=-\nabla^2$ 的点源响应；电势核是 $G/\epsilon_0$。
2. **在 $L=-\nabla^2$ 与 $\nabla^2$ 两种约定间漏掉负号。** 应从 $LG=\delta$ 与 $L\phi=f$ 同时出发，不能靠改变电势零点修正符号。
3. **认为 Green 函数让边界条件消失。** 边界条件决定 $G$；若使用自由空间核解决有限区域问题，就会漏掉边界积分。
4. **把所有线性问题都写成平移不变卷积。** 有边界或空间变系数时，$G(\mathbf r,\mathbf r')$ 一般不能简化成 $G(\mathbf r-\mathbf r')$。
5. **给 Neumann Green 函数同时施加齐次法向导数和单个净单位源。** 两者违反积分兼容条件，必须减去零模并另外固定电势常数。
6. **把 Helmholtz 核 $e^{ikR}/(4\pi R)$ 当成无条件唯一。** 外行、入行或其他齐次部分由辐射及边界条件选择。

## 九、练习与反馈

### 练习 1：符号检查

若改用 $\mathcal L=\nabla^2$，使 $\mathcal L\phi=-\rho/\epsilon_0$，相应 Green 函数和卷积中的源应怎样写？

**答案：** 可取 $\mathcal G=-1/(4\pi R)$，满足 $\nabla^2\mathcal G=\delta^{(3)}$；再写

$$
\phi=\int\mathcal G\left(-\frac{\rho}{\epsilon_0}\right)d^3r'
=\frac1{4\pi\epsilon_0}\int\frac{\rho}{R}\,d^3r'.
$$

两个负号相消，与本文约定得到同一电势。

### 练习 2：Dirichlet 边界项

接地边界满足 $\phi|_S=0$。使用 $G_D|_S=0$ 时，完整表示中的两个边界积分是否都为零？

**答案：** 是。含 $\phi\,\partial_{n'}G_D$ 的项因 $\phi=0$ 消失，含 $G_D\,\partial_{n'}\phi$ 的项因 $G_D=0$ 消失。因此接地区域内的势只由 $G_D$ 与体电荷积分给出。

### 练习 3：Neumann 兼容性

区域内没有体电荷，却在整个边界指定恒定且非零的外法向导数 $\partial_n\phi=b$。在有限体积区域中是否总有解？

**答案：** 不一定。兼容条件要求 $b|S|=0$，所以当 $b\ne0$ 时无解。物理上，边界总电通量不能在内部净电荷为零时非零。

---

## 十、本篇结论

在算符约定

$$
L=-\nabla^2,
\qquad
LG=\delta^{(3)}
$$

下，自由空间 Green 函数为 $G=1/(4\pi R)$，局域电荷分布产生的电势为

$$
\boxed{
\phi(\mathbf r)
=\frac1{4\pi\epsilon_0}
\int_{\mathbb R^3}\frac{\rho(\mathbf r')}{|\mathbf r-\mathbf r'|}\,d^3r'
}.
$$

这个公式成立的关键不是“做一次卷积”本身，而是点源归一化、算符符号和无穷远边界条件三者一致。有限区域必须保留 Green 第二恒等式给出的边界项；Neumann 问题还要处理积分兼容条件与常数零模。下一阶段可在具体导体几何中通过镜像法或本征函数展开构造相应 Green 函数。

## 参考资料

1. David J. Griffiths, *Introduction to Electrodynamics*, 4th ed., Cambridge University Press, 2017, Chapter 3。用于静电唯一性、边值问题与镜像法的教材背景。
2. John David Jackson, *Classical Electrodynamics*, 3rd ed., Wiley, 1998, §§1.7–1.10。用于 Dirichlet/Neumann 边值问题、Green 函数和镜像法。
3. George B. Arfken, Hans J. Weber, and Frank E. Harris, *Mathematical Methods for Physicists*, 7th ed., Academic Press, 2013, Chapter 9。用于 Green 函数作为逆算符核及不同微分算符的推广。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.6 张量初步](https://zhuanlan.zhihu.com/p/2061253800107701070)

**课程目录：** [电动力学目录](https://zhuanlan.zhihu.com/p/2062546254274601559)
<!-- zhihu-nav-bottom:end -->
