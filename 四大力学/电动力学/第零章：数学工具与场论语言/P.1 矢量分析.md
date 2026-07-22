---
created: 2026-07-17T00:55:57+08:00
modified: 2026-07-17T01:19:49+08:00
zhihu-title: P.1 矢量分析
zhihu-topics: 电动力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061253023981212089
zhihu-toc: true
zhihu-created-at: 2026-07-17 01:19
review-status: reviewing
physics-checked: true
derivation-checked: true
readability-checked: true
teaching-checked: true
zhihu-layout-checked: false
last-reviewed:
zhihu-updated-at: 2026-07-21 22:07
---
# P.1 矢量分析

电动力学同时使用场的局部变化和区域上的总效应。梯度、散度、旋度描述局部变化；线积分、面积分、体积分描述有限区域；Gauss 定理和 Stokes 定理把两种描述联系起来。本篇只在三维欧氏空间的直角坐标系中推导公式。柱坐标和球坐标中的表达式将在 [[P.3 正交曲线坐标系]] 中从尺度因子重新推导，不能把直角坐标公式机械替换变量后使用。

## 学习目标

完成本篇后，读者应当能够：

1. 从全微分推导梯度和方向导数，并解释梯度与等值面的关系；
2. 从小长方体净通量推导散度，从小矩形环量推导旋度；
3. 正确使用 Gauss 定理、Stokes 定理和梯度线积分定理，并检查取向；
4. 逐分量证明电动力学常用的矢量恒等式，说明它们成立所需的光滑性；
5. 用量纲、积分定理和简单场例交叉检查散度、旋度及其物理解释。

## 先修知识、约定与路线

先修知识包括一元和多元微积分、偏导数、链式法则、二重与三重积分，以及矢量的点积和叉积。

本文采用：

- 右手直角坐标系 $(x,y,z)$；
- 正交单位基矢 $(\hat{\mathbf x},\hat{\mathbf y},\hat{\mathbf z})$；
- 位置矢量 $\mathbf r=x\hat{\mathbf x}+y\hat{\mathbf y}+z\hat{\mathbf z}$；
- 标量场记为 $\phi(\mathbf r)$，矢量场记为 $\mathbf F(\mathbf r)$；
- 重复指标求和只在明确引入指标记号后使用；
- 除非特别说明，场至少具有公式所需阶数的连续偏导。

推导路线是：先从场的微小变化定义梯度，再从通量和环量定义散度与旋度；随后证明三个积分定理和常用恒等式，最后通过完整例题与误区检查建立可操作的判断方法。

---

## 一、场、微分与算符

### 1. 标量场和矢量场

若区域 $D\subset\mathbb R^3$ 中每一点 $\mathbf r$ 对应一个实数，得到标量场

$$
\phi:D\to\mathbb R.
$$

温度 $T(\mathbf r)$、电势 $\phi(\mathbf r)$ 和质量密度 $\rho(\mathbf r)$ 都是标量场。若每一点对应一个三维矢量，得到矢量场

$$
\mathbf F:D\to\mathbb R^3,
\qquad
\mathbf F=F_x\hat{\mathbf x}+F_y\hat{\mathbf y}+F_z\hat{\mathbf z}.
$$

电场、磁场、流体速度和电流密度都是矢量场。场的单位由具体物理量决定；微分算符会额外带来长度的倒数。例如，若 $\phi$ 的单位是伏特，则 $\nabla\phi$ 的单位是伏特每米。

### 2. nabla 算符

在直角坐标系中定义

$$
\nabla
=\hat{\mathbf x}\frac{\partial}{\partial x}
+\hat{\mathbf y}\frac{\partial}{\partial y}
+\hat{\mathbf z}\frac{\partial}{\partial z}.
$$

$\nabla$ 不是普通矢量，而是矢量形式的微分算符。它与标量场、矢量场做形式上的乘法时分别产生梯度、散度和旋度：

$$
\nabla\phi,\qquad
\nabla\cdot\mathbf F,\qquad
\nabla\times\mathbf F.
$$

这些公式在直角坐标中简洁，是因为基矢不随位置改变。曲线坐标的基矢和尺度因子依赖位置，必须重新推导。

---

## 二、梯度：标量场的一阶变化

### 1. 从全微分出发

设 $\phi(x,y,z)$ 在点 $\mathbf r$ 可微。位移

$$
d\mathbf r=dx\,\hat{\mathbf x}+dy\,\hat{\mathbf y}+dz\,\hat{\mathbf z}
$$

引起的标量变化为全微分

$$
d\phi
=\frac{\partial\phi}{\partial x}dx
+\frac{\partial\phi}{\partial y}dy
+\frac{\partial\phi}{\partial z}dz.
$$

定义梯度

$$
\boxed{
\nabla\phi
=\frac{\partial\phi}{\partial x}\hat{\mathbf x}
+\frac{\partial\phi}{\partial y}\hat{\mathbf y}
+\frac{\partial\phi}{\partial z}\hat{\mathbf z}
},
$$

则全微分可以写成坐标无关的点积形式

$$
\boxed{d\phi=\nabla\phi\cdot d\mathbf r}.
$$

这不是记忆技巧，而是梯度的基本特征：它是唯一能够把任意微小位移线性映射为一阶标量变化的矢量。

### 2. 方向导数

沿单位矢量 $\hat{\mathbf u}$ 取曲线

$$
\mathbf r(s)=\mathbf r_0+s\hat{\mathbf u}.
$$

链式法则给出

$$
\frac{d}{ds}\phi(\mathbf r_0+s\hat{\mathbf u})
=\nabla\phi\cdot\hat{\mathbf u}.
$$

因此点 $\mathbf r_0$ 处沿 $\hat{\mathbf u}$ 的方向导数是

$$
\boxed{D_{\hat{\mathbf u}}\phi=\nabla\phi\cdot\hat{\mathbf u}}.
$$

由 Cauchy--Schwarz 不等式，

$$
|D_{\hat{\mathbf u}}\phi|
\le |\nabla\phi|\,|\hat{\mathbf u}|
=|\nabla\phi|.
$$

当 $\hat{\mathbf u}$ 与 $\nabla\phi$ 同向时取最大值 $|\nabla\phi|$；反向时取最小值 $-|\nabla\phi|$。

### 3. 梯度与等值面

设等值面由 $\phi(\mathbf r)=C$ 给出，$d\mathbf r_t$ 是面内切向位移。沿等值面移动时 $d\phi=0$，所以

$$
\nabla\phi\cdot d\mathbf r_t=0.
$$

只要 $\nabla\phi\ne0$，梯度就垂直于等值面。若 $\nabla\phi=0$，该点是驻点，不能仅凭梯度确定等值面的法向。

### 4. 完整例题：方向导数

取

$$
\phi(x,y,z)=x^2y+3z.
$$

梯度为

$$
\nabla\phi
=2xy\,\hat{\mathbf x}
+x^2\,\hat{\mathbf y}
+3\,\hat{\mathbf z}.
$$

在点 $(1,2,0)$，

$$
\nabla\phi=(4,1,3).
$$

沿方向 $(1,2,2)$ 的单位矢量是 $\hat{\mathbf u}=(1,2,2)/3$，所以

$$
D_{\hat{\mathbf u}}\phi
=(4,1,3)\cdot\frac{(1,2,2)}{3}
=4.
$$

最大方向导数为 $|\nabla\phi|=\sqrt{26}$。由于 $4<\sqrt{26}$，结果符合方向导数上界。

**学习检查点。** 若把方向矢量 $(1,2,2)$ 直接代入点积而不先单位化，算出的量为什么不是方向导数？因为曲线参数每增加一个单位时，实际位移长度不再是一个单位；结果同时包含了方向和参数尺度。方向导数按单位长度定义，因此必须使用单位矢量。

---

## 三、散度：单位体积的净流出

### 1. 通量

有向面元写成

$$
d\mathbf a=\hat{\mathbf n}\,dS,
$$

其中 $\hat{\mathbf n}$ 是所选法向。矢量场穿过曲面 $S$ 的通量定义为

$$
\Phi_S=\int_S\mathbf F\cdot d\mathbf a.
$$

闭合曲面的法向约定为向外，记为

$$
\oint_{\partial V}\mathbf F\cdot d\mathbf a.
$$

点积只保留法向分量；切向分量不穿过曲面。

### 2. 坐标无关定义

以点 $\mathbf r_0$ 为中心取逐渐缩小的规则体积 $\Delta V$，缩小时保持各方向尺度同阶，排除越来越细长或扁平的退化体积。若极限存在，定义

$$
\boxed{
\nabla\cdot\mathbf F(\mathbf r_0)
=\lim_{\Delta V\to0}
\frac{1}{\Delta V}
\oint_{\partial(\Delta V)}
\mathbf F\cdot d\mathbf a
}.
$$

散度是标量。它描述净流出通量的一阶体密度，而不是场线图上“箭头是否看起来张开”。

### 3. 直角坐标表达式的推导

取边长为 $\Delta x,\Delta y,\Delta z$ 的小长方体，并让三条边以有界的长宽比同时趋于零。两张垂直于 $x$ 轴的面给出的通量必须先在面上积分：

$$
\begin{aligned}
\Delta\Phi_x
={}&\int_{y-\Delta y/2}^{y+\Delta y/2}
\int_{z-\Delta z/2}^{z+\Delta z/2}\\
&\left[
F_x\left(x+\frac{\Delta x}{2},y',z'\right)
-F_x\left(x-\frac{\Delta x}{2},y',z'\right)
\right]dz'\,dy'.
\end{aligned}
$$

在长方体中心对被积函数做一阶 Taylor 展开。切向的一阶奇函数项在对称积分区间内抵消，保留法向差的主导项：

$$
\Delta\Phi_x
=\frac{\partial F_x}{\partial x}
\Delta x\Delta y\Delta z
+o(\Delta V).
$$

另外两对面同理。因此

$$
\Delta\Phi
=\left(
\frac{\partial F_x}{\partial x}
+\frac{\partial F_y}{\partial y}
+\frac{\partial F_z}{\partial z}
\right)\Delta V+o(\Delta V).
$$

除以 $\Delta V$ 并取极限：

$$
\boxed{
\nabla\cdot\mathbf F
=\frac{\partial F_x}{\partial x}
+\frac{\partial F_y}{\partial y}
+\frac{\partial F_z}{\partial z}
}.
$$

若 $\mathbf F$ 的单位是 $[F]$，散度单位为 $[F]/\mathrm m$。

### 4. 正负号的含义

- $\nabla\cdot\mathbf F>0$：小体积的净流出为正；
- $\nabla\cdot\mathbf F<0$：净流入为正；
- $\nabla\cdot\mathbf F=0$：局部净通量的一阶项为零。

零散度不表示 $\mathbf F=0$，也不表示场没有空间变化。例如均匀场和刚体旋转速度场都可以散度为零。

---

## 四、旋度：单位面积的局部环量

### 1. 环量与取向

沿有向闭合曲线 $C$ 的环量定义为

$$
\Gamma_C=\oint_C\mathbf F\cdot d\mathbf l.
$$

曲线方向与张成曲面的法向由右手定则关联：右手四指沿曲线正向弯曲，拇指指向正法向。

### 2. 坐标无关定义

取法向为 $\hat{\mathbf n}$ 的规则小面元 $\Delta S$，边界为 $\partial(\Delta S)$；缩小时保持面元各方向尺度同阶。旋度在该法向的分量定义为

$$
\boxed{
(\nabla\times\mathbf F)\cdot\hat{\mathbf n}
=\lim_{\Delta S\to0}
\frac{1}{\Delta S}
\oint_{\partial(\Delta S)}
\mathbf F\cdot d\mathbf l
}.
$$

对三个独立法向求得分量后，就确定矢量 $\nabla\times\mathbf F$。

### 3. $z$ 分量的推导

在 $xy$ 平面取边长 $\Delta x,\Delta y$ 的小矩形，沿从 $+z$ 方向看为逆时针的方向积分。上下两边的总贡献为

$$
\begin{aligned}
\Delta\Gamma_x
={}&\int_{x-\Delta x/2}^{x+\Delta x/2}
\left[
F_x\left(x',y-\frac{\Delta y}{2},z\right)
-F_x\left(x',y+\frac{\Delta y}{2},z\right)
\right]dx'\\
={}&-\frac{\partial F_x}{\partial y}
\Delta x\Delta y+o(\Delta S).
\end{aligned}
$$

第一项来自下边沿 $+x$ 方向的积分，第二项来自上边沿 $-x$ 方向的积分。左右两边分别沿 $-y$ 和 $+y$ 方向，它们的总贡献为

$$
\begin{aligned}
\Delta\Gamma_y
={}&\int_{y-\Delta y/2}^{y+\Delta y/2}
\left[
F_y\left(x+\frac{\Delta x}{2},y',z\right)
-F_y\left(x-\frac{\Delta x}{2},y',z\right)
\right]dy'\\
={}&\frac{\partial F_y}{\partial x}
\Delta x\Delta y+o(\Delta S).
\end{aligned}
$$

将两组贡献相加，除以 $\Delta S=\Delta x\Delta y$ 并取极限：

$$
(\nabla\times\mathbf F)_z
=\frac{\partial F_y}{\partial x}
-\frac{\partial F_x}{\partial y}.
$$

循环置换坐标后得到

$$
\boxed{
\begin{aligned}
\nabla\times\mathbf F
={}&
\left(\frac{\partial F_z}{\partial y}
-\frac{\partial F_y}{\partial z}\right)\hat{\mathbf x}\\
&+\left(\frac{\partial F_x}{\partial z}
-\frac{\partial F_z}{\partial x}\right)\hat{\mathbf y}\\
&+\left(\frac{\partial F_y}{\partial x}
-\frac{\partial F_x}{\partial y}\right)\hat{\mathbf z}.
\end{aligned}
}
$$

旋度单位同样是 $[F]/\mathrm m$。行列式写法可用于记忆，但真正定义来自单位面积环量。

### 4. 旋度与角速度

刚体绕 $z$ 轴以角速度 $\boldsymbol\Omega=\Omega\hat{\mathbf z}$ 转动时，

$$
\mathbf v=\boldsymbol\Omega\times\mathbf r
=(-\Omega y,\Omega x,0).
$$

直接求导：

$$
\nabla\cdot\mathbf v=0,
\qquad
\nabla\times\mathbf v=2\boldsymbol\Omega.
$$

这个例子也体现了可微速度场的一般局部结论：速度梯度的反对称部分描述局部刚体转动，其角速度是 $\nabla\times\mathbf v/2$。剪切和形变还由速度梯度的对称部分描述，不能由旋度单独确定。小桨轮图像对速度场有用，但对一般抽象矢量场，旋度的可靠定义仍是环量密度。

**学习检查点。** 散度和旋度都含一阶空间导数，但输出类型不同：散度把矢量场映为标量，检验小体积的净流出；旋度把矢量场映为矢量，其法向分量检验小回路的环量面密度。判断一个具体问题该算哪一个量时，应先看问题问的是“穿出边界多少”还是“沿边界绕行多少”。

---

## 五、Laplacian：二阶空间变化

标量场的 Laplacian 定义为梯度的散度：

$$
\boxed{
\nabla^2\phi
=\nabla\cdot(\nabla\phi)
=\frac{\partial^2\phi}{\partial x^2}
+\frac{\partial^2\phi}{\partial y^2}
+\frac{\partial^2\phi}{\partial z^2}
}.
$$

对直角坐标中的矢量场，矢量 Laplacian 按分量定义：

$$
\nabla^2\mathbf F
=(\nabla^2F_x)\hat{\mathbf x}
+(\nabla^2F_y)\hat{\mathbf y}
+(\nabla^2F_z)\hat{\mathbf z}.
$$

静电势在无电荷区域满足 Laplace 方程 $\nabla^2\phi=0$；有电荷时满足 Poisson 方程 $\nabla^2\phi=-\rho/\epsilon_0$。Laplacian 含两个空间导数，单位为原场单位除以长度平方。

---

## 六、三个积分定理

### 1. 梯度线积分定理

设光滑曲线 $C$ 由 $\mathbf r(t)$ 参数化，起点 $A$、终点 $B$。链式法则给出

$$
\frac{d}{dt}\phi(\mathbf r(t))
=\nabla\phi\cdot\frac{d\mathbf r}{dt}.
$$

对 $t$ 积分：

$$
\boxed{
\int_C\nabla\phi\cdot d\mathbf l
=\phi(B)-\phi(A)
}.
$$

因此梯度场的线积分只依赖端点。对闭合曲线，

$$
\oint_C\nabla\phi\cdot d\mathbf l=0.
$$

### 2. Gauss 散度定理

设 $V$ 是分片光滑有界体积，边界 $\partial V$ 取外法向，$\mathbf F$ 在包含 $V$ 的区域上具有连续一阶偏导，则

$$
\boxed{
\oint_{\partial V}\mathbf F\cdot d\mathbf a
=\int_V\nabla\cdot\mathbf F\,dV
}.
$$

证明思路是把 $V$ 分割成小体积。每个内部公共面在相邻两块中法向相反，通量严格抵消，只留下外边界；小体积通量除以体积在极限中趋于散度。

**交叉检查。** 取 $\mathbf F=a\mathbf r$，有 $\nabla\cdot\mathbf F=3a$。半径 $R$ 的球面上 $\mathbf F\cdot\hat{\mathbf n}=aR$，所以

$$
\oint_S\mathbf F\cdot d\mathbf a
=aR\cdot4\pi R^2=4\pi aR^3.
$$

体积分为

$$
\int_V3a\,dV
=3a\cdot\frac{4\pi R^3}{3}
=4\pi aR^3.
$$

两边相等。

### 3. Stokes 定理

设 $S$ 是分片光滑有向曲面，边界 $\partial S=C$ 按右手定则取向，$\mathbf F$ 在邻域内具有连续一阶偏导，则

$$
\boxed{
\oint_C\mathbf F\cdot d\mathbf l
=\int_S(\nabla\times\mathbf F)\cdot d\mathbf a
}.
$$

把曲面划分为小面元时，内部公共边被相反方向走过而抵消，只留下外边界。

对刚体旋转场，在 $xy$ 平面半径 $R$ 的圆周上，

$$
\oint_C\mathbf v\cdot d\mathbf l
=\Omega R\cdot2\pi R
=2\pi\Omega R^2.
$$

而

$$
\int_S(\nabla\times\mathbf v)\cdot d\mathbf a
=2\Omega\cdot\pi R^2.
$$

两种计算一致，也检查了旋度中的因子 2。

**学习检查点。** 使用 Gauss 定理时必须同时给出闭合曲面和外法向；使用 Stokes 定理时必须同时给出曲面法向及与之匹配的边界方向。若把 Stokes 定理中的边界方向反转而保持法向不变，线积分变号，等式将不再成立。

---

## 七、常用恒等式及成立条件

以下恒等式假设相关场至少具有连续二阶偏导。若场在点电荷、面电荷或材料界面处不光滑，需要分区使用或改用分布形式。

### 1. 梯度的旋度为零

以 $z$ 分量为例：

$$
[\nabla\times(\nabla\phi)]_z
=\frac{\partial^2\phi}{\partial x\partial y}
-\frac{\partial^2\phi}{\partial y\partial x}=0.
$$

其余分量同理，因此

$$
\boxed{\nabla\times(\nabla\phi)=0}.
$$

### 2. 旋度的散度为零

展开后六个混合偏导两两抵消：

$$
\boxed{\nabla\cdot(\nabla\times\mathbf F)=0}.
$$

所以用 $\mathbf B=\nabla\times\mathbf A$ 定义磁场时，$\nabla\cdot\mathbf B=0$ 自动满足。

### 3. 旋度的旋度

逐分量展开，或使用 Levi--Civita 符号恒等式

$$
\epsilon_{ijk}\epsilon_{klm}
=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl},
$$

可得

$$
\boxed{
\nabla\times(\nabla\times\mathbf F)
=\nabla(\nabla\cdot\mathbf F)-\nabla^2\mathbf F
}.
$$

以指标记号验证：

$$
\begin{aligned}
[\nabla\times(\nabla\times\mathbf F)]_i
&=\epsilon_{ijk}\partial_j
(\epsilon_{klm}\partial_lF_m)\\
&=(\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl})
\partial_j\partial_lF_m\\
&=\partial_i(\partial_mF_m)-\partial_j\partial_jF_i.
\end{aligned}
$$

右边分别是 $\nabla(\nabla\cdot\mathbf F)$ 和 $\nabla^2\mathbf F$ 的第 $i$ 个分量。

### 4. 乘积法则

由普通微积分的 Leibniz 法则逐分量可得：

$$
\boxed{
\nabla\cdot(\phi\mathbf F)
=\nabla\phi\cdot\mathbf F
+\phi\nabla\cdot\mathbf F
},
$$

$$
\boxed{
\nabla\times(\phi\mathbf F)
=\nabla\phi\times\mathbf F
+\phi\nabla\times\mathbf F
},
$$

$$
\boxed{
\nabla\cdot(\mathbf F\times\mathbf G)
=\mathbf G\cdot(\nabla\times\mathbf F)
-\mathbf F\cdot(\nabla\times\mathbf G)
},
$$

$$
\boxed{
\begin{aligned}
\nabla(\mathbf F\cdot\mathbf G)
={}&(\mathbf F\cdot\nabla)\mathbf G
+(\mathbf G\cdot\nabla)\mathbf F\\
&+\mathbf F\times(\nabla\times\mathbf G)
+\mathbf G\times(\nabla\times\mathbf F).
\end{aligned}
}
$$

这里

$$
(\mathbf F\cdot\nabla)
=F_x\frac{\partial}{\partial x}
+F_y\frac{\partial}{\partial y}
+F_z\frac{\partial}{\partial z}
$$

是方向微分算符，不是点积后的普通数。

**学习检查点。** 在旋度的旋度恒等式中，若已知 $\nabla\cdot\mathbf F=0$，可以删去的是 $\nabla(\nabla\cdot\mathbf F)$，而不是 $\nabla^2\mathbf F$。这一点正是从无源 Maxwell 方程推导波动方程时的关键化简。

---

## 八、局部结论与整体结论

恒等式 $\nabla\times\nabla\phi=0$ 总是给出“梯度场无旋”。反方向则需要区域条件：

> 若 $\nabla\times\mathbf F=0$ 且定义域单连通，则存在标量势 $\phi$，使 $\mathbf F=\nabla\phi$。

“单连通”排除了绕洞而不能连续收缩为一点的闭合路径。若定义域有洞，即使处处旋度为零，闭合环量也可能不为零。

同样，$\nabla\cdot(\nabla\times\mathbf A)=0$ 保证旋度场无散；从无散场构造全局矢势还需要适当的区域拓扑和边界条件。[[P.2 亥姆霍兹定理]] 将系统说明这些存在性与唯一性条件。

---

## 九、电动力学中的两个完整应用

### 1. 从点电荷电势恢复电场

对 $r\ne0$，

$$
\phi(\mathbf r)
=\frac{q}{4\pi\epsilon_0r},
\qquad
r=\sqrt{x^2+y^2+z^2}.
$$

先计算

$$
\frac{\partial r}{\partial x}=\frac{x}{r},
\qquad
\frac{\partial}{\partial x}\frac1r
=-\frac{x}{r^3}.
$$

循环置换后，

$$
\nabla\frac1r
=-\frac{x\hat{\mathbf x}+y\hat{\mathbf y}+z\hat{\mathbf z}}{r^3}
=-\frac{\hat{\mathbf r}}{r^2}.
$$

所以

$$
\boxed{
\mathbf E=-\nabla\phi
=\frac{q}{4\pi\epsilon_0}
\frac{\hat{\mathbf r}}{r^2}
}.
$$

量纲为伏特每米，也等于牛顿每库仑。对 $r\ne0$ 直接计算 $\nabla\cdot\mathbf E=0$，但包围原点的球面通量是 $q/\epsilon_0$。两者并不矛盾：原点是奇点，普通导数推导不适用。[[P.4 狄拉克 δ 函数]] 将证明

$$
\nabla\cdot\frac{\hat{\mathbf r}}{r^2}
=4\pi\delta^{(3)}(\mathbf r).
$$

### 2. 从 Maxwell 旋度方程得到波动方程

无源真空中，

$$
\nabla\times\mathbf E
=-\frac{\partial\mathbf B}{\partial t},
\qquad
\nabla\times\mathbf B
=\mu_0\epsilon_0
\frac{\partial\mathbf E}{\partial t},
\qquad
\nabla\cdot\mathbf E=0.
$$

对 Faraday 定律取旋度：

$$
\nabla\times(\nabla\times\mathbf E)
=-\frac{\partial}{\partial t}
(\nabla\times\mathbf B).
$$

这里交换了空间旋度与时间偏导，要求 $\mathbf B$ 的相关混合偏导连续。

左边用旋度的旋度恒等式，右边代入 Ampere--Maxwell 定律：

$$
\nabla(\nabla\cdot\mathbf E)-\nabla^2\mathbf E
=-\mu_0\epsilon_0
\frac{\partial^2\mathbf E}{\partial t^2}.
$$

无源条件使第一项为零，得到

$$
\boxed{
\nabla^2\mathbf E
-\mu_0\epsilon_0
\frac{\partial^2\mathbf E}{\partial t^2}=0
}.
$$

与标准波动方程比较，波速为

$$
c=\frac1{\sqrt{\mu_0\epsilon_0}}.
$$

量纲也与波速一致：

$$
[\mu_0\epsilon_0]=\mathrm{s^2\,m^{-2}},
\qquad
\left[\frac1{\sqrt{\mu_0\epsilon_0}}\right]
=\mathrm{m\,s^{-1}}.
$$

每一步都能定位到一个 Maxwell 方程、无源条件或矢量恒等式；没有在 Maxwell 方程之外增加介质本构假设。

---

## 十、常见误区

1. **散度为零等于场为零。** 错。均匀场可以非零但散度为零。
2. **旋度为零必然存在全局单值势。** 只在适当的单连通区域成立。
3. **场线弯曲就表示旋度非零。** 场线外观不能替代局部环量计算。
4. **旋度就是角速度。** 对刚体速度场，局部角速度是旋度的一半；对一般场不能直接这样解释。
5. **直角坐标分量公式可直接用于球坐标。** 错。曲线坐标的基矢和尺度因子随位置改变。
6. **Gauss 定理和 Stokes 定理不需要光滑条件。** 奇点、界面和分布源需要分区或广义函数处理。
7. **改变曲面不影响 Stokes 积分。** 只有共享同一边界、场在两曲面之间光滑且不存在奇点时才能这样比较。
8. **$\nabla$ 可以完全当作普通矢量移动。** 它是微分算符，乘积法则决定它作用于右侧哪些因子。

---

## 十一、练习与反馈

### 练习 1：梯度与等值面

对

$$
\phi=x^2+y^2-2z
$$

求点 $(1,-1,0)$ 处的梯度、最大方向导数，并写出该点等值面的法向。

**答案：**

$$
\nabla\phi=(2x,2y,-2),
$$

所以该点为 $(2,-2,-2)$，最大方向导数为 $2\sqrt3$，等值面法向可取 $(1,-1,-1)$。

### 练习 2：散度与旋度

给定

$$
\mathbf F=(xy,yz,zx),
$$

计算散度和旋度。

**答案：**

$$
\nabla\cdot\mathbf F=y+z+x,
$$

$$
\nabla\times\mathbf F
=(-y,-z,-x).
$$

例如 $x$ 分量是 $\partial_y(zx)-\partial_z(yz)=0-y=-y$。

### 练习 3：判断推理

某区域内 $\nabla\times\mathbf F=0$。能否立即断言任意闭合曲线的环量为零？

**答案：** 不能。还需确认区域的拓扑和场的光滑性。若区域单连通且场足够光滑，可以由 Stokes 定理和势函数存在性得到闭合环量为零；若区域含洞或排除了奇点，结论可能失败。

### 练习 4：Gauss 定理

用 Gauss 定理计算 $\mathbf F=(x,y,z)$ 穿过边长 $L$、中心在原点的立方体的外向通量。

**答案：**

$$
\nabla\cdot\mathbf F=3,
$$

立方体体积为 $L^3$，所以通量为 $3L^3$。直接逐面计算应得到相同结果。

---

## 十二、参考来源

本文采用的定义、积分定理、矢量恒等式与电磁学符号可在下列标准教材中交叉核对：

1. David J. Griffiths, *Introduction to Electrodynamics*, 第 4 版，第 1 章 “Vector Analysis”。本文的梯度、散度、旋度、积分定理以及点电荷场例与该章约定一致。
2. George B. Arfken, Hans J. Weber, Frank E. Harris, *Mathematical Methods for Physicists*, 第 7 版，矢量分析相关章节。本文用指标记号证明矢量恒等式所需的 Levi--Civita 符号和 Kronecker delta 恒等式可由此核对。
3. Jerrold E. Marsden, Anthony J. Tromba, *Vector Calculus*, 第 6 版，关于方向导数、曲线与曲面积分、Gauss 定理和 Stokes 定理的章节。本文对取向和光滑条件的表述以这些定理的标准假设为依据。

这些来源使用的势函数正负号可能随应用而变：一般数学陈述可写 $\mathbf F=\nabla\phi$；静电学按 $\mathbf E=-\nabla\phi$ 定义电势。两种写法服务于不同约定，并不矛盾。

---

## 十三、本篇结论

直角坐标中的核心定义是

$$
d\phi=\nabla\phi\cdot d\mathbf r,
$$

$$
\nabla\cdot\mathbf F
=\lim_{\Delta V\to0}
\frac1{\Delta V}
\oint\mathbf F\cdot d\mathbf a,
$$

$$
(\nabla\times\mathbf F)\cdot\hat{\mathbf n}
=\lim_{\Delta S\to0}
\frac1{\Delta S}
\oint\mathbf F\cdot d\mathbf l.
$$

Gauss 定理把散度的体积分变成闭合面通量，Stokes 定理把旋度的面积分变成边界环量。梯度无旋、旋度无散和旋度的旋度恒等式建立了势函数、矢势和电磁波方程的数学基础。

最容易误用的边界是：局部微分条件不自动等价于全局势表示；奇点和非单连通区域必须单独处理。下一篇 [[P.2 亥姆霍兹定理]] 将回答一个更完整的问题：给定矢量场的散度、旋度和边界行为，能否重建该矢量场，以及重建在什么条件下唯一。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.2 亥姆霍兹定理](https://zhuanlan.zhihu.com/p/2061253406543570185)

**课程目录：** [电动力学目录](https://zhuanlan.zhihu.com/p/2062546254274601559)
<!-- zhihu-nav-bottom:end -->
