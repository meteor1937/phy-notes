---
created: 2026-07-17T00:59:00+08:00
modified: 2026-07-17T15:09:36+08:00
zhihu-title: P.4 狄拉克 δ 函数
zhihu-topics: 电动力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061253558738097791
zhihu-toc: true
zhihu-created-at: 2026-07-17 01:25
zhihu-updated-at: 2026-07-21 14:52
review-status: reviewing
physics-checked: true
derivation-checked: true
readability-checked: true
teaching-checked: true
zhihu-layout-checked: false
last-reviewed:
---
<!-- zhihu-nav-top:start -->
**上一节：** [P.3 正交曲线坐标系](https://zhuanlan.zhihu.com/p/2061253464097871213)
<!-- zhihu-nav-top:end -->

# P.4 Dirac $\delta$ 分布

点电荷把有限电荷集中在零体积中，普通函数无法同时表达“除一点外处处为零”和“全空间积分有限”。Dirac $\delta$ 的正确对象不是逐点函数值，而是它对测试函数的作用。本篇从分布定义出发，推导平移、缩放、复合函数和导数规则，再严格证明

$$
\nabla^2\frac1r=-4\pi\delta^{(3)}(\mathbf r).
$$

## 学习目标

读者完成本篇后应能够：

1. 用测试函数定义 $\delta$ 分布，解释为什么“$\delta(0)=\infty$”不是严格定义；
2. 推导 $\delta(ax)$、$\delta[g(x)]$ 和 $\delta'$ 的作用规则；
3. 在直角、柱和球坐标中正确写出点、线、面源；
4. 用挖去小球和测试函数的方法证明 $\nabla^2(1/r)$ 的分布恒等式；
5. 处理移动点源、点偶极子和 Green 函数中的 Jacobian。

## 先修知识与约定

需要掌握换元积分、分部积分、多元微积分、Gauss 定理和 [[P.1 矢量分析]]。本文在 SI 制下讨论电荷密度；位置矢量为 $\mathbf r$，$r=|\mathbf r|$。测试函数 $f$ 属于 $C_c^\infty$，即无限次可微且具有紧支集；紧支集保证无穷远边界项消失。等号若涉及 $\delta$，默认表示分布意义下相等，即两边作用于任意测试函数都给出相同结果。

---

## 一、从普通函数到分布

### 1. 点源的要求

位于 $x=a$ 的单位点源应满足

$$
\int_{-\infty}^{\infty}f(x)\delta(x-a)\,dx=f(a).
$$

这条筛选公式就是定义。更抽象地写，

$$
\boxed{\langle\delta_a,f\rangle=f(a)}.
$$

$\delta_a$ 是把测试函数映射到数 $f(a)$ 的线性泛函。它没有必须指定的普通函数值，因此“原点无穷大、其他地方为零”只能作为直观图像，不能用于代数运算。

取 $f=1$（或在源附近等于 1 的紧支测试函数）得到归一化

$$
\int_{-\infty}^{\infty}\delta(x-a)\,dx=1.
$$

### 2. 逼近恒等族

高斯函数族

$$
\delta_\sigma(x)
=\frac1{\sqrt{2\pi}\sigma}
e^{-x^2/(2\sigma^2)}
$$

满足积分为 1。令 $x=\sigma u$：

$$
\int f(x)\delta_\sigma(x)\,dx
=\int\frac{du}{\sqrt{2\pi}}
e^{-u^2/2}f(\sigma u).
$$

当 $\sigma\to0^+$ 时，$f(\sigma u)\to f(0)$，在适当支配条件下可把极限移入积分，得到 $f(0)$。所以

$$
\delta_\sigma\longrightarrow\delta
$$

是分布收敛，不是逐点收敛。不同的归一化窄峰可以趋于同一个 $\delta$ 分布。

---

## 二、一维 $\delta$ 的运算规则

### 1. 平移与偶性

由定义，

$$
\int f(x)\delta(x-a)\,dx=f(a).
$$

令 $a=0$ 并作 $x\mapsto-x$ 换元，可得

$$
\boxed{\delta(-x)=\delta(x)}.
$$

### 2. 缩放

对非零常数 $a$，

$$
\begin{aligned}
\int f(x)\delta(ax)\,dx
&=\frac1{|a|}
\int f\left(\frac u a\right)\delta(u)\,du\\
&=\frac{f(0)}{|a|}.
\end{aligned}
$$

因此

$$
\boxed{\delta(ax)=\frac1{|a|}\delta(x)}.
$$

绝对值来自 $a<0$ 时积分方向反转，不能写成 $1/a$。

### 3. 复合函数

设 $g(x)$ 的零点 $x_i$ 都是单根，即

$$
g(x_i)=0,qquad g'(x_i)\ne0.
$$

在每个零点附近，

$$
g(x)\simeq g'(x_i)(x-x_i).
$$

分区换元并使用缩放规则：

$$
\boxed{
\delta[g(x)]
=\sum_i\frac{\delta(x-x_i)}{|g'(x_i)|}
}.
$$

例如

$$
\delta(x^2-a^2)
=\frac{\delta(x-a)+\delta(x+a)}{2|a|}
\quad(a\ne0).
$$

若零点不是单根，公式失效，必须回到调节或积分定义。

### 4. 分布导数

普通分部积分启发定义

$$
\boxed{
\langle\delta'_a,f\rangle
=-\langle\delta_a,f'\rangle=-f'(a)
}.
$$

一般地，

$$
\langle\delta_a^{(n)},f\rangle
=(-1)^nf^{(n)}(a).
$$

负号不是约定装饰，而是把微分从分布转移到测试函数时分部积分产生的。

### 5. 常用恒等式

对任意光滑 $f$，

$$
f(x)\delta(x-a)=f(a)\delta(x-a).
$$

取 $f(x)=x$、$a=0$：

$$
x\delta(x)=0.
$$

对它求分布导数：

$$
\delta(x)+x\delta'(x)=0,
$$

所以

$$
\boxed{x\delta'(x)=-\delta(x)}.
$$

---

## 三、多维 $\delta$ 与坐标变换

### 1. 直角坐标

三维点源定义为

$$
\boxed{
\delta^{(3)}(\mathbf r-\mathbf r_0)
=\delta(x-x_0)\delta(y-y_0)\delta(z-z_0)
},
$$

并满足

$$
\int d^3r\,
f(\mathbf r)
\delta^{(3)}(\mathbf r-\mathbf r_0)
=f(\mathbf r_0).
$$

点电荷和移动点电荷分别写成

$$
\rho(\mathbf r)=q\delta^{(3)}(\mathbf r-\mathbf r_0),
$$

$$
\rho(\mathbf r,t)
=q\delta^{(3)}[\mathbf r-\mathbf r_q(t)],
\qquad
\mathbf J=\rho\,\mathbf v_q(t).
$$

### 2. 一般坐标的 Jacobian

若 $\mathbf r=\mathbf r(u^1,u^2,u^3)$，并且坐标变换在源点 $u_0$ 附近一一对应、Jacobian 不为零，体积元为

$$
d^3r=|J(u)|\,du^1du^2du^3,
$$

则点源必须带逆 Jacobian：

$$
\delta^{(3)}(\mathbf r-\mathbf r_0)
=\frac{
\delta(u^1-u_0^1)
\delta(u^2-u_0^2)
\delta(u^3-u_0^3)}
{|J(u_0)|}.
$$

这里 $J=\det[\partial(x,y,z)/\partial(u^1,u^2,u^3)]$。分母在 $u_0$ 取值，是因为三个一维 $\delta$ 已经把积分变量筛选到源点。该公式不能直接用于坐标奇点。

例如在柱坐标 $(s,\varphi,z)$ 的正则点 $s_0>0$，

$$
\delta^{(3)}(\mathbf r-\mathbf r_0)
=\frac{\delta(s-s_0)\delta(\varphi-\varphi_0)\delta(z-z_0)}{s_0}.
$$

对位于 $z$ 轴上的点，角坐标不唯一，更适合写成 $\delta(x)\delta(y)\delta(z-z_0)$，或在轴对称积分中使用 $\delta_+(s)\delta(z-z_0)/(2\pi s)$，其中 $s\in[0,\infty)$。

球坐标的 Jacobian 为 $r^2\sin\theta$。原点的两个角坐标都不唯一，因此不能把上一条正则坐标公式机械地代入 $r_0=0$。对原点处的各向同性点源，常用径向分布表示

$$
\boxed{
\delta^{(3)}(\mathbf r)
=\frac{\delta_+(r)}{4\pi r^2}
}
$$

其中 $r$ 的积分范围是 $[0,\infty)$，半轴上的 $\delta_+$ 定义为

$$
\int_0^\infty\delta_+(r)f(r)\,dr=f(0).
$$

这个记号避免了把整条实轴上的偶分布 $\delta(r)$ 限制到半轴时常见的 $1/2$ 约定歧义。验证时，体积元 $r^2\,dr\,d\Omega$ 消去分母，角积分再给出 $4\pi$，最终得到 $f(0)$。

### 3. 线源和面源

$z=0$ 平面上的面电荷密度 $\sigma(x,y)$ 对应

$$
\rho(x,y,z)=\sigma(x,y)\delta(z).
$$

$z$ 轴上的线电荷密度 $\lambda(z)$ 对应

$$
\rho(x,y,z)=\lambda(z)\delta(x)\delta(y).
$$

量纲检查：三维、二维和一维 $\delta$ 的单位分别是 $\mathrm m^{-3}$、$\mathrm m^{-2}$、$\mathrm m^{-1}$，从而 $\rho$ 始终是 $\mathrm C\,\mathrm m^{-3}$。

---

## 四、严格证明 $\nabla^2(1/r)=-4\pi\delta^{(3)}(\mathbf r)$

### 1. 原点之外

对 $r>0$，球对称 Laplacian 给出

$$
\nabla^2\frac1r
=\frac1{r^2}\frac{d}{dr}
\left(r^2\frac{d}{dr}\frac1r\right)
=0.
$$

但这不能决定原点的分布贡献，因为 $1/r$ 在原点不光滑。

### 2. 测试函数与挖球区域

要证明分布恒等式，必须证明对任意测试函数 $f$，

$$
\left\langle\nabla^2\frac1r,f\right\rangle
=-4\pi f(0).
$$

按分布导数定义，

$$
\left\langle\nabla^2\frac1r,f\right\rangle
=\int_{\mathbb R^3}\frac1r\nabla^2f\,d^3r.
$$

挖去半径 $\epsilon$ 的小球，记

$$
\Omega_\epsilon
=\mathbb R^3\setminus B_\epsilon.
$$

在 $\Omega_\epsilon$ 中 $\nabla^2(1/r)=0$。Green 第二恒等式给出

$$
\int_{\Omega_\epsilon}
\left[
\frac1r\nabla^2f
-f\nabla^2\frac1r
\right]d^3r
=\oint_{\partial\Omega_\epsilon}
\left[
\frac1r\partial_nf
-f\partial_n\frac1r
\right]dS.
$$

测试函数有紧支集，外边界贡献为零。小球内边界对区域 $\Omega_\epsilon$ 的外法向是 $\mathbf n=-\hat{\mathbf r}$，所以

$$
\partial_nf=-\partial_rf,
\qquad
\partial_n\frac1r=\frac1{\epsilon^2}.
$$

在 $r=\epsilon$ 上 $dS=\epsilon^2d\Omega$，故

$$
\int_{\Omega_\epsilon}
\frac1r\nabla^2f\,d^3r
=-\epsilon\int d\Omega\,\partial_rf
-\int d\Omega\,f(\epsilon,\Omega).
$$

当 $\epsilon\to0$，第一项趋零，第二项趋于

$$
-\int d\Omega\,f(0)=-4\pi f(0).
$$

因此

$$
\boxed{
\nabla^2\frac1r
=-4\pi\delta^{(3)}(\mathbf r)
}.
$$

这个证明没有跨越奇点直接使用普通 Gauss 定理，而是在光滑区域计算后取极限。

### 3. 平移后的形式

平移不改变 Laplacian：

$$
\boxed{
\nabla^2\frac1{|\mathbf r-\mathbf r'|}
=-4\pi\delta^{(3)}(\mathbf r-\mathbf r')
}.
$$

所以自由空间 Poisson Green 函数

$$
G(\mathbf r,\mathbf r')
=\frac1{4\pi|\mathbf r-\mathbf r'|}
$$

满足

$$
\nabla^2G=-\delta^{(3)}(\mathbf r-\mathbf r').
$$

---

## 五、从点源到电磁模型

### 1. 点电荷的 Poisson 方程

取

$$
\phi(\mathbf r)
=\frac{q}{4\pi\varepsilon_0r}.
$$

使用上面的恒等式：

$$
\nabla^2\phi
=-\frac q{\varepsilon_0}\delta^{(3)}(\mathbf r)
=-\frac\rho{\varepsilon_0}.
$$

这同时解释了为什么 $r>0$ 处 Laplacian 为零，而包围原点的电通量不为零。

### 2. 理想点偶极子

电偶极矩为 $\mathbf p$、位于原点的理想偶极子可写成

$$
\boxed{
\rho(\mathbf r)
=-\mathbf p\cdot\nabla
\delta^{(3)}(\mathbf r)
}.
$$

检查总电荷：

$$
\int\rho\,d^3r=0.
$$

检查偶极矩第 $i$ 分量：

$$
\begin{aligned}
\int r_i\rho(\mathbf r)\,d^3r
&=-p_j\int r_i\partial_j\delta^{(3)}(\mathbf r)\,d^3r\\
&=p_j\partial_jr_i|_{\mathbf r=0}
=p_j\delta_{ij}=p_i.
\end{aligned}
$$

因此符号和归一化都正确。

### 3. 推迟时间中的 Jacobian

若 $g(t')=0$ 定义推迟时刻且根为单根，

$$
\delta[g(t')]
=\frac{\delta(t'-t_R)}{|g'(t_R)|}.
$$

对移动点电荷，

$$
g(t')=t-t'
-\frac{|\mathbf r-\mathbf r_q(t')|}{c}.
$$

记 $\mathbf R(t')=\mathbf r-\mathbf r_q(t')$、$R=|\mathbf R|$、$\hat{\mathbf R}=\mathbf R/R$ 和 $\boldsymbol\beta=\mathbf v_q/c$。在保持观测事件 $(\mathbf r,t)$ 不变时，

$$
\frac{dR}{dt'}
=\hat{\mathbf R}\cdot\frac{d\mathbf R}{dt'}
=-\hat{\mathbf R}\cdot\mathbf v_q(t'),
$$

因此

$$
g'(t')
=-1+\hat{\mathbf R}\cdot\boldsymbol\beta,
$$

对亚光速点电荷，$|\hat{\mathbf R}\cdot\boldsymbol\beta|<1$，故 $g'(t_R)<0$ 且

$$
|g'(t_R)|=1-\hat{\mathbf R}\cdot\boldsymbol\beta.
$$

所以筛选积分产生

$$
\frac1{1-\hat{\mathbf R}\cdot\boldsymbol\beta}.
$$

这个因子是 Lienard--Wiechert 势的一部分，不能遗漏。

---

## 六、Fourier 表示

采用正、反变换约定

$$
\widetilde f(k)=\int_{-\infty}^{\infty}f(x)e^{-ikx}\,dx,
$$

$$
f(x)=\int\frac{dk}{2\pi}\widetilde f(k)e^{ikx},
$$

$k$ 空间中的常数 $1$ 经过反变换给出

$$
\boxed{
\delta(x)=\int_{-\infty}^{\infty}
\frac{dk}{2\pi}e^{ikx}
}.
$$

验证时把它放入测试函数积分：

$$
\int dx\,f(x)
\int\frac{dk}{2\pi}e^{ik(x-a)}
=f(a),
$$

其中最后一步使用 Fourier 反演公式。三维形式为

$$
\delta^{(3)}(\mathbf r)
=\int\frac{d^3k}{(2\pi)^3}
e^{i\mathbf k\cdot\mathbf r}.
$$

不同 Fourier 约定会重新分配 $2\pi$ 因子，必须保持正反变换一致。

---

## 七、常见误区

1. **把 $\delta(0)=\infty$ 当作可代数运算的数。** $\delta$ 由对测试函数的作用定义。
2. **把窄峰极限理解为逐点极限。** 正确概念是分布收敛。
3. **缩放公式遗漏绝对值。** $\delta(ax)=\delta(x)/|a|$。
4. **复合函数公式用于重根。** $g'(x_i)=0$ 时不能使用单根公式。
5. **在球坐标只写 $\delta(r)$。** 必须包含体积元 Jacobian、角向归一化，并说明半轴 $\delta_+$ 的约定。
6. **跨越 $r=0$ 直接使用普通微分。** 奇点贡献必须用分布或挖球极限处理。
7. **移动点源只做位置替换。** 隐式根会产生 $1/|g'|$。
8. **忽略 $\delta'$ 分部积分的负号。** 负号决定偶极矩方向。

---

## 八、练习与反馈

### 练习 1

化简 $\delta(3x-6)$。

**答案：**

$$
\delta(3x-6)=\frac13\delta(x-2).
$$

### 练习 2

计算

$$
\int_{-\infty}^{\infty}
(x^2+1)\delta'(x-2)\,dx.
$$

**答案：**

$$
-\frac d{dx}(x^2+1)\bigg|_{x=2}=-4.
$$

### 练习 3

验证 $z=0$ 面源 $\rho=\sigma(x,y)\delta(z)$ 对任意柱体积分得到面上的总电荷。

**答案：** 先对 $z$ 积分得到 1，剩余积分为

$$
Q=\int_A\sigma(x,y)\,dx\,dy.
$$

### 练习 4

为什么只证明 $\nabla^2(1/r)=0$ 对 $r>0$ 成立，不能推出它在全空间为零？

**答案：** 原点是奇点，普通导数结论没有覆盖原点；分布导数还包含由任意小球面通量固定的 $-4\pi\delta^{(3)}(\mathbf r)$。

### 学习检查点

不看正文，尝试回答以下问题；答案都能在相邻章节中定位：

1. 为什么 $\delta[g(x)]$ 的单根公式必须含 $|g'(x_i)|$，而不是 $g'(x_i)$？
2. 在球坐标原点，为什么不能直接把 $J=r^2\sin\theta$ 代入一般坐标公式？
3. 挖球证明中，小球面法向为什么是 $-\hat{\mathbf r}$？这个符号如何决定最终的 $-4\pi$？
4. 推迟时间积分中的 Jacobian 在什么条件下等于 $1-\hat{\mathbf R}\cdot\boldsymbol\beta$？

---

## 九、本篇结论

$\delta$ 是分布，不是具有无穷大点值的普通函数。它的核心规则是筛选性质、换元 Jacobian 和分布分部积分。三维恒等式

$$
\nabla^2\frac1r=-4\pi\delta^{(3)}(\mathbf r)
$$

把点源、Poisson 方程和自由空间 Green 函数统一起来。最容易出错的地方是坐标 Jacobian、复合函数的 $1/|g'|$ 和奇点处普通微积分的失效。下一篇 [[P.5 格林函数方法入门]] 将把点源响应叠加为一般边值问题的解。

## 参考资料

1. L. Schwartz, *Théorie des distributions*, Hermann, 1950--1951：分布与分布导数的奠基性定义。
2. G. B. Arfken, H. J. Weber, F. E. Harris, *Mathematical Methods for Physicists*, 7th ed., Academic Press, 2013，第 1 章：Dirac $\delta$、坐标变换与 Fourier 表示。
3. J. D. Jackson, *Classical Electrodynamics*, 3rd ed., Wiley, 1998，第 1、6、14 章：点源、Green 函数与推迟势中的 $\delta$ 函数。
4. D. J. Griffiths, *Introduction to Electrodynamics*, 4th ed., Pearson, 2013，第 1、3、10 章：三维 $\delta$、Poisson 方程和 Liénard--Wiechert 势。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [P.5 格林函数方法入门](https://zhuanlan.zhihu.com/p/2061253748400337330)

**课程目录：** [电动力学目录](https://zhuanlan.zhihu.com/p/2062546254274601559)
<!-- zhihu-nav-bottom:end -->
