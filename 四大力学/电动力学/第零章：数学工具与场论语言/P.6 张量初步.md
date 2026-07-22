---
created: 2026-07-17T00:59:59+08:00
modified: 2026-07-17T15:26:57+08:00
zhihu-title: P.6 张量初步
zhihu-topics: 电动力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061253800107701070
zhihu-toc: true
zhihu-created-at: 2026-07-17 01:28
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
**上一节：** [P.5 格林函数方法入门](https://zhuanlan.zhihu.com/p/2061253748400337330)
<!-- zhihu-nav-top:end -->

## P.6 张量初步

梯度、散度、旋度、点乘和叉乘足以紧凑地写出四条麦克斯韦方程，但场的动量守恒还需要描述“某个方向的动量穿过某个法向的面”。这同时涉及动量分量和面法向两个方向，因而需要二阶张量。麦克斯韦应力张量 $T_{ij}$ 的第一个指标 $i$ 标记动量分量，第二个指标 $j$ 标记所穿过面的法向；$T_{ij}n_j$ 是法向为 $\mathbf n$ 的面上的牵引力第 $i$ 分量。

本篇只讨论三维欧几里得空间中的笛卡尔正交坐标。这个限制使指标可以用 $\delta_{ij}$ 升降而不改变数值，也使普通偏导数足够使用；一般曲线坐标还需要度规、协变导数和体积元，不能直接照搬所有公式。

### 学习目标

完成本篇后，读者应能：

1. 区分自由指标和哑指标，并检查一条指标等式是否合法；
2. 用正交变换律判断标量、矢量和二阶张量；
3. 推导 $\epsilon_{ijk}\epsilon_{imn}=\delta_{jm}\delta_{kn}-\delta_{jn}\delta_{km}$ 并用它证明叉乘恒等式；
4. 区分并积 $A_iB_j$、转置并积 $B_iA_j$ 和缩并 $A_iB_i$；
5. 从麦克斯韦方程逐项推出场的局域动量守恒式，并解释其中的符号约定。

### 先修知识与约定

需要掌握 [[P.1 矢量分析]] 中的梯度、散度、旋度和 Gauss 定理。空间指标 $i,j,k,\ldots$ 取 $1,2,3$；重复指标默认求和。本文使用 SI 制，$\varepsilon_0\mu_0=c^{-2}$，电荷和电流密度分别为 $\rho$、$\mathbf J$。为简化记号，二阶张量既写作 $\mathbf T$，也写其分量 $T_{ij}$。


### 一、为什么要引入张量指标记号

用矢量记号写出来的公式很紧凑，但在需要分量操作时——比如证明一个含叉乘的恒等式，或者计算两个矢量的并积——矢量记号就显得不够用了。我们想要一种既能显式写出分量、又能紧凑地求和的记号。

**爱因斯坦求和约定**解决了这个矛盾。在本文的笛卡尔欧几里得约定中，同一项里重复两次的指标默认从 $1$ 到 $3$ 求和；这种被求和的指标称为哑指标，可以一致地改名。只出现一次的指标称为自由指标，它必须在等号两边逐项出现且位置一致。一个指标在同一项里出现三次通常是非法写法。

比如，矢量 $\mathbf{A}$ 和 $\mathbf{B}$ 的点乘可以写为：

$$
\mathbf{A}\cdot\mathbf{B} = \sum_{i=1}^{3} A_i B_i \quad \rightarrow \quad A_i B_i
$$

第二个写法里，指标 $i$ 出现了两次，按照求和约定，自动对所有 $i$ 求和。$\sum$ 符号被省略了，但意思没变。这看起来只是省了点笔墨，但到后面处理多个矢量和复杂恒等式时，这种紧凑性会带来实质性的简化。

矢量 $\mathbf{A}$ 的散度在指标记号下写为 $\partial_i A_i$（$\partial_i \equiv \partial/\partial x_i$）。旋度的第 $i$ 个分量呢？它需要一种能表达“交错”结构的符号。

指标记号的关键不是省略 $\sum$，而是显式编码变换律。坐标轴作正交变换 $x'_i=R_{ij}x_j$，其中 $R_{ik}R_{jk}=\delta_{ij}$，则矢量和二阶张量分别满足

$$
A'_i=R_{ij}A_j,
\qquad
T'_{ij}=R_{ip}R_{jq}T_{pq}.
$$

缩并后的 $A_iB_i$ 在正交变换下不变：

$$
A'_iB'_i
=R_{ij}R_{ik}A_jB_k
=\delta_{jk}A_jB_k
=A_jB_j.
$$

这说明“张量”不是任意带很多下标的数组，而是按规定方式随坐标基变换的几何对象。


### 二、克罗内克 $\delta$ 符号

**克罗内克 $\delta$ 符号** $\delta_{ij}$ 定义为：

$$
\boxed{\delta_{ij} = \begin{cases} 1, & i = j \\ 0, & i \neq j \end{cases}}
$$

它的两个指标都写成下标（在本文的正交笛卡尔基中，上下标数值没有区别）。$\delta_{ij}$ 就是 $3\times3$ 单位矩阵的 $(i,j)$ 分量。

$\delta_{ij}$ 在求和中的核心作用是“换指标”。比如：

$$
\delta_{ij} A_j = A_i
$$

因为当 $j$ 被求和时，只有 $j=i$ 那一项活下来。这看起来是平凡的操作，但在复杂的恒等式化简中，$\delta$ 符号能帮我们消去多余的求和，把表达式大幅度简化。

两个 $\delta$ 符号的乘积也有直观的几何意义。$\delta_{ij}\delta_{ji}$ 中 $i,j$ 都是哑指标。把求和暂时显式写出，先固定 $i$ 对 $j$ 求和，再对 $i$ 求和：

$$
\sum_i\sum_j\delta_{ij}\delta_{ji}
=\sum_i\delta_{ii}
=3.
$$

在爱因斯坦记号中直接写作 $\delta_{ij}\delta_{ji}=\delta_{ii}=3$；最后一个 $\delta_{ii}$ 已经包含对 $i$ 的求和。它就是三维单位矩阵的迹。


### 三、列维-奇维塔 $\epsilon$ 符号

**列维-奇维塔符号** $\epsilon_{ijk}$ 在选定右手正交基后由以下规则定义：

- $\epsilon_{123} = 1$
- 交换任意两个指标，符号变号：$\epsilon_{ijk} = -\epsilon_{jik} = -\epsilon_{ikj}$
- 如果有两个指标相同，则 $\epsilon_{ijk} = 0$（因为交换相同指标时，$\epsilon$ 必须同时等于它自己和它的相反数）

因此，三个指标的所有 $3^3 = 27$ 种组合中，只有 $6$ 种非零：$\epsilon_{123} = \epsilon_{231} = \epsilon_{312} = 1$，$\epsilon_{132} = \epsilon_{213} = \epsilon_{321} = -1$，其余 $21$ 种全为零。

$\epsilon_{ijk}$ 是表达“反对称结构”的完美工具。矢量叉乘可以用它写出来：$\mathbf{C} = \mathbf{A}\times\mathbf{B}$ 的第 $i$ 个分量为：

$$
C_i = \epsilon_{ijk} A_j B_k
$$

验证一下：$C_1 = \epsilon_{1jk} A_j B_k = \epsilon_{123}A_2 B_3 + \epsilon_{132}A_3 B_2 = A_2 B_3 - A_3 B_2$，确实是叉乘的 $x$ 分量。

旋度也同样可以用 $\epsilon$ 符号写出：$\nabla\times\mathbf{A}$ 的第 $i$ 个分量为 $\epsilon_{ijk} \partial_j A_k$（$\partial_j \equiv \partial/\partial x_j$）。

必须注意反射变换。对正交矩阵 $R$，

$$
R_{ip}R_{jq}R_{kr}\epsilon_{pqr}
=(\det R)\epsilon_{ijk}.
$$

因此 $\epsilon_{ijk}$ 在保持取向的旋转（$\det R=1$）下像三阶张量，在反射（$\det R=-1$）下多出一个负号；它更准确地说是赝张量分量。叉乘也因此产生轴矢量，而不是普通极矢量。


### 四、核心恒等式：$\epsilon_{ijk}\epsilon_{imn}$ 的缩并

这是整个张量代数中最重要的一个恒等式。它把两个 $\epsilon$ 符号的乘积（对公共指标求和）转化为 $\delta$ 符号的组合：

$$
\boxed{\epsilon_{ijk} \epsilon_{imn} = \delta_{jm}\delta_{kn} - \delta_{jn}\delta_{km}}
$$

这里指标 $i$ 是公共的求和指标，$j,k$ 和 $m,n$ 分别来自两个 $\epsilon$ 符号。

**逐情形证明。** 若 $j=k$ 或 $m=n$，左边至少一个 $\epsilon$ 为零；右边两个 $\delta$ 乘积也相等，故两边都为零。以下只需考虑 $j\ne k$ 且 $m\ne n$。

- 若有序对 $(m,n)=(j,k)$，只有使 $(i,j,k)$ 成为 $1,2,3$ 排列的那个 $i$ 有贡献，左边为 $\epsilon_{ijk}^2=1$；右边为 $1-0=1$。
- 若 $(m,n)=(k,j)$，同一个非零项变号，左边为 $-1$；右边为 $0-1=-1$。
- 其余情况下，集合 $\{m,n\}$ 与 $\{j,k\}$ 不同。两个 $\epsilon$ 若要同时非零，各自的三个指标都必须包含 $1,2,3$，这要求两组有序对包含同样两个元素，产生矛盾；故左边为零。右边的两个乘积也都为零。

所有指标组合已经穷尽，因此恒等式成立。它也可以记成二维子行列式

$$
\epsilon_{ijk}\epsilon_{imn}
=\begin{vmatrix}
\delta_{jm}&\delta_{jn}\\
\delta_{km}&\delta_{kn}
\end{vmatrix}.
$$

进一步缩并会得到更简单的形式。令 $k=n$ 并对 $k$ 求和：

$$
\epsilon_{ijk} \epsilon_{imk} = \delta_{jm}\delta_{kk} - \delta_{jk}\delta_{km}
$$

$\delta_{kk}=3$，$\delta_{jk}\delta_{km} = \delta_{jm}$（对 $k$ 求和时，$\delta_{jk}$ 把 $k$ 换成 $j$）。所以：

$$
\boxed{\epsilon_{ijk} \epsilon_{imk} = 3\delta_{jm} - \delta_{jm} = 2\delta_{jm}}
$$

再缩并一次：令 $j=m$ 并对 $j$ 求和，$\epsilon_{ijk} \epsilon_{ijk} = 2\delta_{jj} = 6$。这个结果是预期之中的：三个非零指标的全排列有 $3! = 6$ 种，每种贡献 $1$，共 $6$。


### 五、用张量语言重写矢量恒等式

有了 $\epsilon_{ijk}\epsilon_{imn}$ 恒等式，以前需要逐分量验证的矢量恒等式现在可以用代数简化一步到位。

#### 5.1 矢量三重积：$\mathbf{A}\times(\mathbf{B}\times\mathbf{C})$

用指标记号写出来：$[\mathbf{A}\times(\mathbf{B}\times\mathbf{C})]_i = \epsilon_{ijk} A_j (\mathbf{B}\times\mathbf{C})_k = \epsilon_{ijk} A_j (\epsilon_{kmn} B_m C_n)$。这里 $\epsilon_{kmn}$ 中的 $k$ 是公共指标。把 $\epsilon_{ijk}$ 和 $\epsilon_{kmn}$ 的指标重新排列，利用 $\epsilon$ 的循环对称性 $\epsilon_{ijk} = \epsilon_{kij}$，我们得到 $\epsilon_{kij}\epsilon_{kmn}$。应用恒等式（把 $k$ 作为公共求和指标，相当于原公式中的 $i$）：

$$
\epsilon_{kij} \epsilon_{kmn} = \delta_{im}\delta_{jn} - \delta_{in}\delta_{jm}
$$

代回去：

$$
\begin{aligned}
[\mathbf{A}\times(\mathbf{B}\times\mathbf{C})]_i &= (\delta_{im}\delta_{jn} - \delta_{in}\delta_{jm}) A_j B_m C_n \\
&= B_i (A_j C_j) - C_i (A_j B_j) \\
&= B_i (\mathbf{A}\cdot\mathbf{C}) - C_i (\mathbf{A}\cdot\mathbf{B})
\end{aligned}
$$

这就完成了证明：$\mathbf{A}\times(\mathbf{B}\times\mathbf{C}) = \mathbf{B}(\mathbf{A}\cdot\mathbf{C}) - \mathbf{C}(\mathbf{A}\cdot\mathbf{B})$。逐分量方法需要分别展开三个分量再合并，指标法则把共同结构集中在一次 $\epsilon$ 缩并中。

#### 5.2 旋度之旋度：$\nabla\times(\nabla\times\mathbf{F})$

$[\nabla\times(\nabla\times\mathbf{F})]_i = \epsilon_{ijk} \partial_j (\epsilon_{kmn} \partial_m F_n) = \epsilon_{kij}\epsilon_{kmn} \partial_j \partial_m F_n$。同样用恒等式，$\epsilon_{kij}\epsilon_{kmn} = \delta_{im}\delta_{jn} - \delta_{in}\delta_{jm}$。代回：

$$
\begin{aligned}
[\nabla\times(\nabla\times\mathbf{F})]_i &= (\delta_{im}\delta_{jn} - \delta_{in}\delta_{jm}) \partial_j \partial_m F_n \\
&= \partial_i (\partial_j F_j) - \partial_j \partial_j F_i \\
&= [\nabla(\nabla\cdot\mathbf{F}) - \nabla^2\mathbf{F}]_i
\end{aligned}
$$

所以 $\nabla\times(\nabla\times\mathbf{F}) = \nabla(\nabla\cdot\mathbf{F}) - \nabla^2\mathbf{F}$。这是电磁波推导中最关键的矢量恒等式，用 $\epsilon$ 符号证明它比逐分量展开简洁得多。

#### 5.3 两个叉乘的点乘

坡印廷定理使用的恒等式是

$$
\nabla\cdot(\mathbf E\times\mathbf B)
=\mathbf B\cdot(\nabla\times\mathbf E)
-\mathbf E\cdot(\nabla\times\mathbf B).
$$

从左边出发并使用乘积法则：

$$
\begin{aligned}
\partial_i(\epsilon_{ijk}E_jB_k)
&=\epsilon_{ijk}(\partial_iE_j)B_k
+\epsilon_{ijk}E_j\partial_iB_k\\
&=B_k\epsilon_{kij}\partial_iE_j
-E_j\epsilon_{jik}\partial_iB_k\\
&=\mathbf B\cdot(\nabla\times\mathbf E)
-\mathbf E\cdot(\nabla\times\mathbf B).
\end{aligned}
$$

第二行第一项用了循环置换 $\epsilon_{ijk}=\epsilon_{kij}$，第二项用了交换前两指标后变号。原式本身等于正的 $\nabla\cdot(\mathbf E\times\mathbf B)$；在坡印廷能量守恒方程中出现负散度，是因为移项，而不是这个恒等式自带负号。


### 六、并矢与二阶张量

两个矢量 $\mathbf{A}$ 和 $\mathbf{B}$ 的**并积**（或张量积）记作 $\mathbf A\mathbf B$，它是一个二阶张量，分量为 $(\mathbf A\mathbf B)_{ij}=A_iB_j$。转置并积 $\mathbf B\mathbf A$ 的同一 $(i,j)$ 分量是 $B_iA_j$，一般不等于 $A_iB_j$。注意 $B_jA_i=A_iB_j$ 只是同两个标量因子的交换，不能把它误认成转置分量。

单位张量（二阶）就是 $\delta_{ij}$。任意二阶张量 $\overset{\leftrightarrow}{\mathbf{T}}$ 可以分解为对称部分和反对称部分：$T_{ij} = \frac{1}{2}(T_{ij} + T_{ji}) + \frac{1}{2}(T_{ij} - T_{ji})$。

两个矢量的点乘有两种推广方式。矢量 $\mathbf{A}$ 点乘张量 $\overset{\leftrightarrow}{\mathbf{T}}$，得到一个矢量：$(\mathbf{A}\cdot\overset{\leftrightarrow}{\mathbf{T}})_j = A_i T_{ij}$（对第一个指标求和）。反过来，$(\overset{\leftrightarrow}{\mathbf{T}}\cdot\mathbf{A})_i = T_{ij} A_j$（对第二个指标求和）。这两个结果一般不同，除非 $\overset{\leftrightarrow}{\mathbf{T}}$ 是对称的。

张量散度的指标位置存在两种教材约定。本文为了与牵引力 $(\mathbf T\cdot\mathbf n)_i=T_{ij}n_j$ 一致，定义

$$
(\nabla\cdot\mathbf T)_i\equiv\partial_jT_{ij}.
$$

有些教材把 $\partial_iT_{ij}$ 记作 $\nabla\cdot\mathbf T$；它等于本文的 $\nabla\cdot\mathbf T^{\mathsf T}$。对称张量两种写法数值相同，但对一般张量必须先确认约定。


### 七、麦克斯韦应力张量的指标形式

真空中含电荷与电流时，麦克斯韦应力张量定义为

$$
T_{ij}
=\varepsilon_0\left(E_iE_j-\frac12\delta_{ij}E^2\right)
+\frac1{\mu_0}\left(B_iB_j-\frac12\delta_{ij}B^2\right).
$$

它是对称张量。$T_{ij}$ 的单位是 $\mathrm{N\,m^{-2}}$，既可解释为应力，也可解释为带符号的动量通量。场动量密度与洛伦兹力密度分别为

$$
\mathbf g=\varepsilon_0\mathbf E\times\mathbf B
=\frac{\mathbf S}{c^2},
\qquad
\mathbf f=\rho\mathbf E+\mathbf J\times\mathbf B,
$$

其中 $\mathbf S=\mu_0^{-1}\mathbf E\times\mathbf B$。目标是证明

$$
\boxed{
\partial_tg_i-\partial_jT_{ij}=-f_i
}
$$

或等价地 $\partial_jT_{ij}=f_i+\partial_tg_i$。注意这里是减去 $\partial_jT_{ij}$；若把真正“向外的场动量通量”定义为 $\Pi_{ij}=-T_{ij}$，守恒式才写成 $\partial_tg_i+\partial_j\Pi_{ij}=-f_i$。

#### 7.1 两个代数准备

对任意矢量场 $\mathbf A$，由两个 $\epsilon$ 的缩并可得

$$
[\mathbf A\times(\nabla\times\mathbf A)]_i
=A_k\partial_iA_k-A_j\partial_jA_i.
$$

因此

$$
(\partial_jA_i)A_j-A_k\partial_iA_k
=-[\mathbf A\times(\nabla\times\mathbf A)]_i.
$$

另一方面，$\partial_i(A_kA_k/2)=A_k\partial_iA_k$，这里使用了乘积法则并把重复指标 $k$ 求和。

#### 7.2 电场部分

对电场部分求散度：

$$
\begin{aligned}
\partial_jT^{(E)}_{ij}
&=\varepsilon_0\partial_j
\left(E_iE_j-\frac12\delta_{ij}E_kE_k\right)\\
&=\varepsilon_0\left[
(\partial_jE_i)E_j+E_i\partial_jE_j-E_k\partial_iE_k
\right]\\
&=\rho E_i
-\varepsilon_0[\mathbf E\times(\nabla\times\mathbf E)]_i\\
&=\rho E_i
+\varepsilon_0[\mathbf E\times\partial_t\mathbf B]_i.
\end{aligned}
$$

第三行用了 Gauss 定律 $\partial_jE_j=\rho/\varepsilon_0$，第四行用了 Faraday 定律 $\nabla\times\mathbf E=-\partial_t\mathbf B$。

#### 7.3 磁场部分

同样地，利用 $\partial_jB_j=0$，

$$
\begin{aligned}
\partial_jT^{(B)}_{ij}
&=\frac1{\mu_0}
\left[(\partial_jB_i)B_j-B_k\partial_iB_k\right]\\
&=-\frac1{\mu_0}[\mathbf B\times(\nabla\times\mathbf B)]_i\\
&=(\mathbf J\times\mathbf B)_i
+\varepsilon_0[(\partial_t\mathbf E)\times\mathbf B]_i.
\end{aligned}
$$

最后一行代入 Ampère--Maxwell 定律

$$
\nabla\times\mathbf B
=\mu_0\mathbf J+\mu_0\varepsilon_0\partial_t\mathbf E
$$

并使用 $-\mathbf B\times\mathbf J=\mathbf J\times\mathbf B$。

#### 7.4 合并与积分形式

相加后，两项时间导数恰好组成叉乘的乘积法则：

$$
\begin{aligned}
\partial_jT_{ij}
&=\rho E_i+(\mathbf J\times\mathbf B)_i\\
&\quad+\varepsilon_0
[\mathbf E\times\partial_t\mathbf B
+(\partial_t\mathbf E)\times\mathbf B]_i\\
&=f_i+\partial_tg_i.
\end{aligned}
$$

这就得到局域守恒式。对固定体积 $V$ 积分并使用 Gauss 定理，

$$
\frac{d}{dt}\int_Vg_i\,dV
=\oint_{\partial V}T_{ij}n_j\,dA
-\int_Vf_i\,dV.
$$

右边第一项是边界应力给体积内系统带来的净动量，第二项是场交给物质的动量。量纲上，$\partial_tg_i$、$\partial_jT_{ij}$ 与 $f_i$ 都是 $\mathrm{N\,m^{-3}}$，为推导提供独立检查。

---

### 八、常见误区

1. **自由指标在等号两边不一致。** 例如 $A_i=B_j$ 没有说明 $i,j$ 的关系；合法张量等式必须逐个自由分量成立。
2. **同一项中一个指标出现三次。** $A_iB_iC_i$ 不符合标准求和约定，应显式写求和或改写所需对象。
3. **把 $\epsilon_{ijk}$ 当作反射下不变的普通张量。** 反射会带来 $\det R=-1$，叉乘结果是轴矢量。
4. **混淆 $A_iB_j$ 与 $B_iA_j$。** 后者是转置并积的分量；$B_jA_i$ 则与 $A_iB_j$ 完全相同。
5. **忘记应力张量与动量通量的负号。** 按本文定义，$T_{ij}$ 满足 $\partial_tg_i-\partial_jT_{ij}=-f_i$；若采用通量 $\Pi_{ij}=-T_{ij}$，散度项才带正号。
6. **把笛卡尔偏导公式直接搬到曲线坐标。** 球坐标和柱坐标需要协变导数或包含尺度因子的专用公式。

### 九、练习与反馈

#### 练习 1：检查指标合法性

判断下列式子是否符合求和约定：$A_i=B_{ij}C_j$、$A_iB_i=C_jD_j$、$A_iB_iC_i$。

**答案。** 第一式合法，自由指标为 $i$、哑指标为 $j$；第二式合法，两边都是标量，左右哑指标可以取不同字母；第三式中 $i$ 出现三次，不符合标准约定。

#### 练习 2：叉乘模长

用 $\epsilon$ 缩并证明

$$
|\mathbf A\times\mathbf B|^2
=A^2B^2-(\mathbf A\cdot\mathbf B)^2.
$$

**答案。**

$$
\begin{aligned}
(\mathbf A\times\mathbf B)_i(\mathbf A\times\mathbf B)_i
&=\epsilon_{ijk}\epsilon_{imn}A_jB_kA_mB_n\\
&=(\delta_{jm}\delta_{kn}-\delta_{jn}\delta_{km})
A_jB_kA_mB_n\\
&=A_jA_jB_kB_k-A_jB_jA_kB_k.
\end{aligned}
$$

#### 练习 3：静电场中的表面力

真空导体表面外侧电场为 $\mathbf E=E_n\mathbf n$，磁场为零。求 $\mathbf T\cdot\mathbf n$。

**答案。** 代入 $E_i=E_nn_i$ 和 $n_jn_j=1$：

$$
T_{ij}n_j
=\varepsilon_0\left(E_n^2n_i-\frac12E_n^2n_i\right)
=\frac12\varepsilon_0E_n^2n_i.
$$

它是沿外法向的电磁压力，单位为 $\mathrm{Pa}$。

---

### 十、本篇结论

指标记号的约束来自变换律，而不只是省略求和号。$\delta_{ij}$ 表示欧几里得度规与单位映射，$\epsilon_{ijk}$ 编码取向和反对称结构；核心缩并式把叉乘计算化为 $\delta$ 的代数。二阶张量 $T_{ij}$ 的两个自由指标分别记录两个方向信息，缩并一个指标后才得到矢量。

在本文的应力约定下，电磁场动量守恒为

$$
\boxed{
\partial_tg_i-\partial_jT_{ij}
=-(\rho\mathbf E+\mathbf J\times\mathbf B)_i
}.
$$

这个式子的适用范围是 SI 制、真空宏观麦克斯韦方程和固定笛卡尔坐标。介质中的应力与动量分配还涉及 Abraham--Minkowski 等约定，不能未经说明直接沿用真空结论。

### 参考资料

1. D. J. Griffiths, *Introduction to Electrodynamics*, 4th ed., Pearson, 2013，第 1、8 章：指标记号、Levi-Civita 符号和电磁动量守恒。
2. J. D. Jackson, *Classical Electrodynamics*, 3rd ed., Wiley, 1998，第 1、6 章：张量代数、麦克斯韦应力张量与守恒律。
3. G. B. Arfken, H. J. Weber, F. E. Harris, *Mathematical Methods for Physicists*, 7th ed., Academic Press, 2013，第 2 章：Cartesian tensors、$\delta$ 与 $\epsilon$ 恒等式。
4. L. D. Landau and E. M. Lifshitz, *The Classical Theory of Fields*, 4th ed., Butterworth-Heinemann, 1975，第 4 章：电磁场能量-动量张量。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [1.1 库仑定律与场概念的建立](https://zhuanlan.zhihu.com/p/2061254017586553481)

**课程目录：** [电动力学目录](https://zhuanlan.zhihu.com/p/2062546254274601559)
<!-- zhihu-nav-bottom:end -->
