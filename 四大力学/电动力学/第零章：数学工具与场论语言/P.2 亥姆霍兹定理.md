---
created: 2026-07-17T00:58:24+08:00
modified: 2026-07-17T01:21:51+08:00
zhihu-title: P.2 亥姆霍兹定理
zhihu-topics: 电动力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061253406543570185
zhihu-toc: true
zhihu-created-at: 2026-07-17 01:21
review-status: reviewing
physics-checked: true
derivation-checked: true
readability-checked: true
teaching-checked: true
zhihu-layout-checked: false
last-reviewed:
zhihu-updated-at: 2026-07-21 22:07
---
<!-- zhihu-nav-top:start -->
**上一节：** [P.1 矢量分析](https://zhuanlan.zhihu.com/p/2061253023981212089)
<!-- zhihu-nav-top:end -->

# P.2 Helmholtz 分解与重建定理

矢量场的散度给出局部净通量密度，旋度给出局部环量密度。Helmholtz 定理回答逆问题：已知这两类局部数据和适当边界信息，能否重建原矢量场？答案是可以，但“适当边界信息”和区域拓扑不能省略。

## 学习目标

完成本篇后，读者应当能够：

1. 陈述全空间 Helmholtz 定理的光滑性、衰减和兼容条件；
2. 从散度和旋度构造标量势、矢量势，并逐步验证重建结果；
3. 证明给定边界行为后重建的唯一性；
4. 在 Fourier 空间写出纵向、横向投影算符；
5. 解释该定理与 Maxwell 方程的关系，同时避免把纵横分解等同于近场远场分解。

## 先修知识与记号

需要掌握 [[P.1 矢量分析]]、[[P.4 狄拉克 δ 函数]] 和自由空间 Poisson Green 函数

$$
\nabla^2\frac1{|\mathbf r-\mathbf r'|}
=-4\pi\delta^{(3)}(\mathbf r-\mathbf r').
$$

令

$$
D(\mathbf r)=\nabla\cdot\mathbf F,
\qquad
\mathbf C(\mathbf r)=\nabla\times\mathbf F.
$$

本文在三维欧氏空间、笛卡尔坐标和 SI 制下讨论实或复矢量场。为避免把技术条件藏在“足够光滑”中，以下取一组便于检验的充分条件：$\mathbf F\in C^2(\mathbb R^3)$，$D$、$\mathbf C$ 连续且按绝对值可积，定义势的两个积分及其所需导数收敛；此外 $\mathbf F(\mathbf r)\to0$，并且下文分部积分出现的无穷远面积分趋于零。紧支撑的 $D$、$\mathbf C$ 满足积分方面的要求；更弱的衰减条件也可以，但必须逐项保证相同的收敛与边界项结论。

推导路线是：先由 $D$、$\mathbf C$ 构造两个 Poisson 势，再分别检查候选场的散度和旋度，最后用无穷远条件排除差场中的调和自由度。Fourier 表述和电磁学例子只解释结果，不替代这三步证明。

---

## 一、为什么还需要边界条件

若 $\mathbf F_1$ 和 $\mathbf F_2$ 有相同散度和旋度，则差

$$
\mathbf W=\mathbf F_1-\mathbf F_2
$$

满足

$$
\nabla\cdot\mathbf W=0,
\qquad
\nabla\times\mathbf W=0.
$$

常矢量场就是非零例子，所以只给散度和旋度并不能在没有边界条件时保证唯一。有限区域中还可能存在由区域拓扑支持的调和场。全空间衰减条件或有限区域边界数据的作用，就是排除这些额外自由度。这里的“调和场”指同时满足 $\nabla\cdot\mathbf H=0$ 与 $\nabla\times\mathbf H=0$ 的场，不只是每个分量满足 Laplace 方程这一局部性质。

旋度数据也不是任意矢量场。恒等式

$$
\nabla\cdot(\nabla\times\mathbf F)=0
$$

要求

$$
\boxed{\nabla\cdot\mathbf C=0}.
$$

这是重建问题的必要兼容条件。

---

## 二、构造纵场和横场

定义标量势

$$
\boxed{
\Phi(\mathbf r)
=\frac1{4\pi}
\int_{\mathbb R^3}
\frac{D(\mathbf r')}{|\mathbf r-\mathbf r'|}
\,d^3r'
}
$$

和矢量势

$$
\boxed{
\mathbf A(\mathbf r)
=\frac1{4\pi}
\int_{\mathbb R^3}
\frac{\mathbf C(\mathbf r')}{|\mathbf r-\mathbf r'|}
\,d^3r'
}.
$$

候选重建场为

$$
\boxed{
\mathbf F_{\mathrm{rec}}
=-\nabla\Phi+\nabla\times\mathbf A
}.
$$

第一项记为纵场

$$
\mathbf F_\parallel=-\nabla\Phi,
$$

第二项记为横场

$$
\mathbf F_\perp=\nabla\times\mathbf A.
$$

名称来自 Fourier 空间中它们分别平行、垂直于波矢，而不是来自实空间场线外观。

---

## 三、验证散度

让 Laplacian 对场点 $\mathbf r$ 作用：

$$
\begin{aligned}
\nabla^2\Phi
&=\frac1{4\pi}
\int D(\mathbf r')
\nabla^2\frac1{|\mathbf r-\mathbf r'|}
\,d^3r'\\
&=-\int D(\mathbf r')
\delta^{(3)}(\mathbf r-\mathbf r')\,d^3r'\\
&=-D(\mathbf r).
\end{aligned}
$$

因此

$$
\begin{aligned}
\nabla\cdot\mathbf F_{\mathrm{rec}}
&=-\nabla^2\Phi
+\nabla\cdot(\nabla\times\mathbf A)\\
&=D.
\end{aligned}
$$

重建场具有要求的散度。

---

## 四、验证旋度

首先证明所构造的 $\mathbf A$ 满足 Coulomb 条件。记 $R=|\mathbf r-\mathbf r'|$，利用

$$
\nabla\frac1R=-\nabla'\frac1R,
$$

得到

$$
\begin{aligned}
\nabla\cdot\mathbf A
&=\frac1{4\pi}
\int\mathbf C(\mathbf r')\cdot
\nabla\frac1R\,d^3r'\\
&=-\frac1{4\pi}
\int\mathbf C\cdot\nabla'\frac1R\,d^3r'.
\end{aligned}
$$

对源坐标分部积分：

$$
\int\mathbf C\cdot\nabla'\frac1R\,d^3r'
=\oint_{S_\infty}
\frac{\mathbf C\cdot d\mathbf a'}R
-\int\frac{\nabla'\cdot\mathbf C}{R}\,d^3r'.
$$

第一项由衰减条件消失，第二项由 $\nabla'\cdot\mathbf C=0$ 消失，所以

$$
\boxed{\nabla\cdot\mathbf A=0}.
$$

与标量势同理，

$$
\nabla^2\mathbf A=-\mathbf C.
$$

于是

$$
\begin{aligned}
\nabla\times\mathbf F_{\mathrm{rec}}
&=\nabla\times(\nabla\times\mathbf A)\\
&=\nabla(\nabla\cdot\mathbf A)-\nabla^2\mathbf A\\
&=\mathbf C.
\end{aligned}
$$

重建场也具有要求的旋度。至此存在性得到显式证明。

---

## 五、唯一性证明

设另一个满足相同 $D$、$\mathbf C$ 和无穷远边界条件的场为 $\mathbf F_2$。差场 $\mathbf W$ 满足

$$
\nabla\cdot\mathbf W=0,
\qquad
\nabla\times\mathbf W=0.
$$

使用恒等式

$$
\nabla^2\mathbf W
=\nabla(\nabla\cdot\mathbf W)
-\nabla\times(\nabla\times\mathbf W),
$$

可得

$$
\nabla^2\mathbf W=0.
$$

所以每个直角分量都是全空间调和函数。取任意分量 $w$ 和包含观察点的球 $B_R$。最大值原理分别用于 $w$ 与 $-w$，给出

$$
|w(\mathbf r)|\leq \max_{|\mathbf x|=R}|w(\mathbf x)|,
\qquad \mathbf r\in B_R.
$$

令 $R\to\infty$；由于 $\mathbf W$ 在无穷远一致趋零，右端趋零，故每个实分量都为零。对于复值场，分别对实部和虚部使用同一论证，也得到每个分量都为零：

$$
\mathbf W=0.
$$

因此

$$
\boxed{\mathbf F_1=\mathbf F_2}.
$$

有限区域中不能直接使用这个结论；还要指定足以排除调和场的边界数据。多连通区域也可能存在既无散又无旋但环量非零的场。例如去掉 $z$ 轴的空间中，$\mathbf H=\hat{\boldsymbol\phi}/s$ 在区域内无散、无旋，但沿环绕 $z$ 轴一周的闭合曲线有 $\oint\mathbf H\cdot d\boldsymbol\ell=2\pi$。这个例子说明局部微分条件不能替代拓扑与边界信息。

有限区域的重建还会保留边界项。由 Green 恒等式或对 $\nabla'\cdot(\mathbf F/R)$、$\nabla'\times(\mathbf F/R)$ 分部积分，可得到体积分之外依赖边界法向量和场迹的贡献。因此，不能把本篇全空间公式原样截断到有限区域；必须使用与给定法向或切向边界数据相配套的 Green 函数。

---

## 六、全空间 Helmholtz 定理

在上述光滑、收敛、兼容和边界条件下，

$$
\boxed{
\begin{aligned}
\mathbf F(\mathbf r)
={}&-\nabla
\left[
\frac1{4\pi}
\int\frac{\nabla'\cdot\mathbf F(\mathbf r')}
{|\mathbf r-\mathbf r'|}\,d^3r'
\right]\\
&+\nabla\times
\left[
\frac1{4\pi}
\int\frac{\nabla'\times\mathbf F(\mathbf r')}
{|\mathbf r-\mathbf r'|}\,d^3r'
\right].
\end{aligned}
}
$$

第一项无旋并携带全部散度；第二项无散并携带全部旋度。分解是非局域的：某一点的纵场和横场一般都需要全空间源数据才能确定。

---

## 七、Fourier 空间中的投影

采用

$$
\mathbf F(\mathbf r)
=\int\frac{d^3k}{(2\pi)^3}
\widetilde{\mathbf F}(\mathbf k)
e^{i\mathbf k\cdot\mathbf r}.
$$

对 $\mathbf k\ne0$，定义单位矢量

$$
\hat{\mathbf k}=\frac{\mathbf k}{k}.
$$

纵向分量是

$$
\boxed{
\widetilde{\mathbf F}_\parallel
=\hat{\mathbf k}
(\hat{\mathbf k}\cdot\widetilde{\mathbf F})
},
$$

横向分量是

$$
\boxed{
\widetilde{\mathbf F}_\perp
=\widetilde{\mathbf F}
-\hat{\mathbf k}
(\hat{\mathbf k}\cdot\widetilde{\mathbf F})
}.
$$

用投影算符写成

$$
P^\parallel_{ij}=\frac{k_ik_j}{k^2},
\qquad
P^\perp_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}.
$$

它们满足

$$
(P^\parallel)^2=P^\parallel,
\quad
(P^\perp)^2=P^\perp,
\quad
P^\parallel P^\perp=0,
\quad
P^\parallel+P^\perp=I.
$$

由于

$$
i\mathbf k\cdot\widetilde{\mathbf F}_\perp=0,
\qquad
i\mathbf k\times\widetilde{\mathbf F}_\parallel=0,
$$

“纵”和“横”的名称在这里完全明确。$\mathbf k=0$ 模由平均值或边界条件单独决定，不能除以 $k^2$。

### 平面波例题

**已知条件：** 设单一非零 Fourier 模为

$$
\mathbf F=\mathbf F_0e^{i\mathbf k\cdot\mathbf r},
\qquad
\mathbf F_0=(F_x,0,F_z),
\qquad
\mathbf k=k\hat{\mathbf z}.
$$

**待求量：** 求它相对于 $\mathbf k$ 的纵向与横向分量，并验证两部分分别无旋、无散。

**解题路线：** 因为 $\hat{\mathbf k}=\hat{\mathbf z}$，先把振幅 $\mathbf F_0$ 投影到 $z$ 方向，再从总振幅中减去纵向投影。于是

$$
\mathbf F_\parallel
=F_z\hat{\mathbf z}e^{ikz},
\qquad
\mathbf F_\perp
=F_x\hat{\mathbf x}e^{ikz}.
$$

对平面波，$\nabla$ 作用等价于乘 $i\mathbf k$。因此可直接检查：

$$
\nabla\times\mathbf F_\parallel=0,
\qquad
\nabla\cdot\mathbf F_\perp=0.
$$

该结果没有单位问题：两个投影算符均无量纲，分解前后的场具有相同单位。易错点是把实空间坐标轴误当成固定的“纵向”；纵向始终由当前非零波矢 $\mathbf k$ 决定。

---

## 八、与 Maxwell 方程的关系

Maxwell 方程在每个时刻给出

$$
\nabla\cdot\mathbf E=\frac\rho{\epsilon_0},
\qquad
\nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t},
$$

$$
\nabla\cdot\mathbf B=0,
\qquad
\nabla\times\mathbf B
=\mu_0\mathbf J
+\mu_0\epsilon_0\frac{\partial\mathbf E}{\partial t}.
$$

给定某一时刻的散度、旋度和空间边界行为，Helmholtz 定理决定该时刻的场。但时变问题中的旋度数据本身含有未知场的时间导数，所以该定理不替代动力学初值。完整演化还需要初始 $\mathbf E,\mathbf B$、源的时间行为和边界条件。

静电情形 $\nabla\times\mathbf E=0$，电场是纯纵向：

$$
\mathbf E=-\nabla
\left[
\frac1{4\pi\epsilon_0}
\int\frac{\rho(\mathbf r')}{R}\,d^3r'
\right].
$$

静磁情形 $\nabla\cdot\mathbf B=0$，磁场是纯横向：

$$
\mathbf B=\nabla\times
\left[
\frac{\mu_0}{4\pi}
\int\frac{\mathbf J(\mathbf r')}{R}\,d^3r'
\right],
$$

其中假定电荷和稳恒电流局域分布、相应积分收敛，并有 $\nabla\cdot\mathbf J=0$。势的量纲也提供交叉检查：$\Phi$ 的量纲为 $[\mathbf F]L$，$\mathbf A$ 同样为 $[\mathbf F]L$，再作用一次空间导数才回到 $[\mathbf F]$。

---

## 九、不能混同的三个概念

1. **纵向/横向是空间投影。** 它由 $\mathbf k$ 方向或散度、旋度定义。
2. **近场/远场是距离渐近分类。** 它按 $1/R^n$ 的主导阶次划分。
3. **束缚场/辐射场是能量传播性质。** 是否把净能量带到无穷远还要检查 Poynting 通量。

横场可以包含近场项，不能把“横向”无条件等同于“辐射”；纵向分量在不同规范和表述中的因果结构也需要结合完整电磁场判断。

---

## 十、常见误区

1. **只知道散度和旋度就无条件唯一。** 还需要边界行为和拓扑条件。
2. **任意 $\mathbf C$ 都能作为旋度。** 必须满足 $\nabla\cdot\mathbf C=0$。
3. **无散场总能在任意区域写成全局旋度。** 多连通区域可能需要额外调和部分。
4. **纵场就是近场，横场就是远场。** 这是两种不同分类。
5. **Helmholtz 定理替代 Maxwell 初值问题。** 它是空间重建定理，不给出时间演化。
6. **矢量势唯一。** $\mathbf A\to\mathbf A+\nabla\chi$ 不改变旋度。
7. **$k=0$ 模可使用 $k_ik_j/k^2$。** 零模必须由平均值或边界条件处理。

---

## 十一、练习与反馈

### 练习 1

把常矢量 $\mathbf a$ 相对于非零波矢 $\mathbf k$ 分解为纵向和横向分量。

**答案：**

$$
\mathbf a_\parallel
=\frac{\mathbf k(\mathbf k\cdot\mathbf a)}{k^2},
\qquad
\mathbf a_\perp
=\mathbf a-\mathbf a_\parallel.
$$

### 练习 2

某全空间光滑场散度和旋度均为零，但不在无穷远衰减。能否由 Helmholtz 唯一性断言它为零？

**答案：** 不能。任意非零常矢量场就是反例。

### 练习 3

为什么稳恒电流必须满足 $\nabla\cdot\mathbf J=0$，才能直接作为磁场重建中的旋度源？

**答案：** 因为 $\nabla\times\mathbf B=\mu_0\mathbf J$ 两边取散度后左边恒为零，所以右边也必须无散；这与稳恒连续性方程一致。

### 练习 4

已知 $\widetilde{\mathbf F}=(1,2,3)$，$\mathbf k=(0,0,k)$。求纵横分量。

**答案：**

$$
\widetilde{\mathbf F}_\parallel=(0,0,3),
\qquad
\widetilde{\mathbf F}_\perp=(1,2,0).
$$

---

## 十二、本篇结论

在适当的光滑、收敛、兼容和边界条件下，矢量场由散度、旋度和边界行为唯一确定，并分解为

$$
\mathbf F=\mathbf F_\parallel+\mathbf F_\perp,
\qquad
\nabla\times\mathbf F_\parallel=0,
\qquad
\nabla\cdot\mathbf F_\perp=0.
$$

该分解在实空间表现为 Poisson Green 函数卷积，在 Fourier 空间表现为相对于 $\mathbf k$ 的正交投影。最容易误用的边界是忽略无穷远条件或有限区域的调和部分：相同散度和旋度并不自动意味着相同矢量场。该定理解释了 Maxwell 方程为什么同时规定散度和旋度，但不能替代时变电磁场的初始边值问题。后续 [[P.5 格林函数方法入门]] 将把这里的自由空间构造推广到具有具体边界条件的区域。

## 参考资料

1. David J. Griffiths, *Introduction to Electrodynamics*, 4th ed., Cambridge University Press, 2017, §1.6 “The Theory of Vector Fields”。本文的全空间重建公式、$1/(4\pi)$ 归一化和唯一性讨论以此为主要教材依据。
2. Jason Cantarella, Dennis DeTurck, and Herman Gluck, “Vector Calculus and the Topology of Domains in 3-Space,” *The American Mathematical Monthly* **109** (2002), 409–442, DOI: 10.1080/00029890.2002.11919871。用于有限区域、调和部分和拓扑限制的讨论。
3. John David Jackson, *Classical Electrodynamics*, 3rd ed., Wiley, 1998, Chapter 1。用于 Green 函数、边值问题与静电唯一性的背景；本篇不把有限区域的结论简化为无边界项的全空间公式。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.3 正交曲线坐标系](https://zhuanlan.zhihu.com/p/2061253464097871213)

**课程目录：** [电动力学目录](https://zhuanlan.zhihu.com/p/2062546254274601559)
<!-- zhihu-nav-bottom:end -->
