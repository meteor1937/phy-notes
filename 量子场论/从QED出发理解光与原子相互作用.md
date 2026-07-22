---
created: 2026-07-21T01:03:43+08:00
modified: 2026-07-21T01:14:50+08:00
zhihu-title: 从QED出发理解光与原子相互作用
zhihu-topics:
  - 量子场论
zhihu-link: https://zhuanlan.zhihu.com/p/2062704363005539458
zhihu-toc: true
zhihu-created-at: 2026-07-21 01:14
---
# 光与原子相互作用的量子场论图景
——从 QED 第一原理推导爱因斯坦系数及其本质启示

## 引言

1917 年，爱因斯坦以天才的热力学洞察力，将黑体辐射的普朗克公式与原子的玻尔兹曼分布相结合，唯象地引入了自发辐射、受激辐射与受激吸收三个系数 $A_{21}, B_{21}, B_{12}$，并证明它们满足两个确定的关系。然而，爱因斯坦的推导存在根本的认知局限：他无法说明为什么三个系数不是相互独立的；他必须借助热平衡的细致平衡条件，而无法触及这些系数本身的微观起源；他也无法解释自发辐射的驱动力究竟来自何处；更关键的是，在他的框架里，不同原子对光子的吸收与辐射被假定为相互独立，这只是一种为凑出平衡分布而引入的假设，无法从第一性原理加以证明。

这些问题，只有在量子电动力学（QED）建立之后才得到彻底解决。本文将从头开始，严格地量子化自由电磁场，建立原子与辐射的偶极相互作用哈密顿量，然后运用含时微扰论精确计算三种跃迁过程的速率。我们将亲眼看到，**受激吸收与受激辐射不过是同一个电子‑光子相互作用顶角在光子数表象下的两种表现形式，它们的系数相等源自相互作用哈密顿量的厄米性所蕴含的时间反演对称性，是一种与热平衡条件无关的微观可逆性**；而**自发辐射则毫无神秘之处，它就是真空电磁场的零点涨落所诱导的受激辐射，其速率严格等于受激辐射系数乘以真空模式密度——那多出来的“相空间积分”，正是对真空里所有可接纳光子模式的求和**。三种现象共享同一项相互作用，爱因斯坦系数关系是量子场论结构内蕴的必然，而非热力学凑巧的拼凑。

## 1. 自由电磁场的量子化

一切从没有电荷、没有原子的纯辐射场开始。在库仑规范（$\nabla\cdot\mathbf{A}=0$，标量势 $\phi=0$）下，电磁场由矢势 $\mathbf{A}(\mathbf{r},t)$ 完全描述，满足波动方程

$$
\nabla^2 \mathbf{A} - \frac{1}{c^2}\frac{\partial^2 \mathbf{A}}{\partial t^2} = 0.
$$

为赋予场粒子性，我们将其限制在一个边长为 $L$ 的立方体腔内，施加周期性边界条件。此时平面波构成完备集：

$$
\mathbf{A}(\mathbf{r},t) = \sum_{\mathbf{k},\lambda} \sqrt{\frac{\hbar}{2\omega_k \varepsilon_0 V}} \left[ a_{\mathbf{k},\lambda} \boldsymbol{\varepsilon}_{\mathbf{k},\lambda} e^{i(\mathbf{k}\cdot\mathbf{r} - \omega_k t)} + a_{\mathbf{k},\lambda}^\dagger \boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^* e^{-i(\mathbf{k}\cdot\mathbf{r} - \omega_k t)} \right],
$$

其中 $V = L^3$ 为量子化体积，$\omega_k = c|\mathbf{k}|$，$\lambda=1,2$ 标记两个独立的横偏振方向。单位偏振矢量满足横波条件 $\mathbf{k}\cdot\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}=0$，以及正交归一条件 $\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^* \cdot \boldsymbol{\varepsilon}_{\mathbf{k},\lambda'} = \delta_{\lambda\lambda'}$。系数 $a_{\mathbf{k},\lambda}$ 与 $a_{\mathbf{k},\lambda}^\dagger$ 被提升为光子湮灭和产生算符，满足玻色对易关系

$$
[a_{\mathbf{k},\lambda}, a_{\mathbf{k}',\lambda'}^\dagger] = \delta_{\mathbf{k}\mathbf{k}'}\delta_{\lambda\lambda'}, \qquad [a_{\mathbf{k},\lambda}, a_{\mathbf{k}',\lambda'}] = [a_{\mathbf{k},\lambda}^\dagger, a_{\mathbf{k}',\lambda'}^\dagger] = 0.
$$

光子数态 $|n_{\mathbf{k},\lambda}\rangle$ 是粒子数算符的本征态，满足

$$
a^\dagger a |n\rangle = n |n\rangle, \qquad a|n\rangle = \sqrt{n}\,|n-1\rangle, \qquad a^\dagger|n\rangle = \sqrt{n+1}\,|n+1\rangle.
$$

注意到一个对后文至关重要的因子：**$\sqrt{n+1}$**。它包含两项：$n$ 部分与初态光子数成正比，将给出受激行为；常数 $1$ 部分即使光子数为零也依然存在，这正是真空涨落效应的代数根源。

电场算符由 $\mathbf{E} = -\partial\mathbf{A}/\partial t$ 得到：

$$
\mathbf{E}(\mathbf{r},t) = i\sum_{\mathbf{k},\lambda} \sqrt{\frac{\hbar\omega_k}{2\varepsilon_0 V}} \left[ a_{\mathbf{k},\lambda} \boldsymbol{\varepsilon}_{\mathbf{k},\lambda} e^{i(\mathbf{k}\cdot\mathbf{r} - \omega_k t)} - a_{\mathbf{k},\lambda}^\dagger \boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^* e^{-i(\mathbf{k}\cdot\mathbf{r} - \omega_k t)} \right]. \tag{1}
$$

这就是我们即将用来与原子耦合的量子化电场。它同时包含了光子的湮灭项（$a$）和产生项（$a^\dagger$），为吸收和辐射过程提供了统一的数学结构。

## 2. 原子模型与偶极相互作用

真实原子是复杂的多电子系统，但只需考虑与近共振辐射耦合的两个能级，便足以揭示三种跃迁过程的本质。设基态 $|g\rangle$ 能量为 $E_g$，激发态 $|e\rangle$ 能量为 $E_e$，跃迁频率 $\omega_0 = (E_e-E_g)/\hbar$。原子自由哈密顿量写为

$$
H_A = E_g |g\rangle\langle g| + E_e |e\rangle\langle e|.
$$

在可见光至紫外波段，辐射波长（$\sim 10^3\text{ Å}$）远大于原子尺度（$\sim 1\text{ Å}$），因此电子感受到的电场在整个原子体积内高度均匀。我们可以略去场随空间变化的相位，直接取原子质心（设为原点）处的电场，这就是电偶极近似。在该近似下，相互作用哈密顿量为

$$
H_I = -\mathbf{d} \cdot \mathbf{E}(0,t),
$$

其中 $\mathbf{d} = -e\,\mathbf{r}$ 是原子电偶极算符。对于无永久偶极矩的原子，$\mathbf{d}$ 仅有非对角元，我们记为

$$
\mathbf{d} = \mathbf{d}_{eg} |e\rangle\langle g| + \mathbf{d}_{ge} |g\rangle\langle e|,
$$

且 $\mathbf{d}_{ge} = \mathbf{d}_{eg}^*$。

于是总哈密顿量为

$$
H = H_0 + H_I, \qquad H_0 = H_A + H_{\text{rad}},
$$

其中 $H_{\text{rad}} = \sum_{\mathbf{k},\lambda} \hbar\omega_k \, a_{\mathbf{k},\lambda}^\dagger a_{\mathbf{k},\lambda}$。

**关键点一：唯一顶角的显现。**  
整个理论只有一个相互作用项 $-\mathbf{d}\cdot\mathbf{E}(0,t)$。它既包含 $a$（湮灭），又包含 $a^\dagger$（产生）。当我们将它作用于不同光子数背景的初态时，将分别给出吸收、受激辐射和自发辐射。三种现象的本质差异，仅在于初态光子数的不同，以及由此决定的 $a$ 与 $a^\dagger$ 谁起作用、以及 Bose 增强因子是 $n$ 还是 $n+1$。这一顶角，在量子场论的费曼图中，就是唯一的电子‑光子顶点。所谓“顶角的交叉对称性和时间反演对称性”，正是指这个顶点在互换初末态（并将吸收光子变为发射光子）时，振幅的模方不变。

## 3. 含时微扰论与费米黄金规则

为计算跃迁速率，我们转入相互作用绘景。令

$$
|\psi_I(t)\rangle = e^{iH_0 t/\hbar} |\psi_S(t)\rangle, \qquad H_I(t) = e^{iH_0 t/\hbar} H_I e^{-iH_0 t/\hbar}.
$$

利用玻色算符和原子投影算符在 $H_0$ 下的简单演化，

$$
a_{\mathbf{k},\lambda}(t) = a_{\mathbf{k},\lambda} e^{-i\omega_k t}, \quad a^\dagger_{\mathbf{k},\lambda}(t) = a^\dagger_{\mathbf{k},\lambda} e^{i\omega_k t},
$$
$$
|e\rangle\langle g|(t) = |e\rangle\langle g| e^{i\omega_0 t}, \quad |g\rangle\langle e|(t) = |g\rangle\langle e| e^{-i\omega_0 t},
$$

我们得到相互作用绘景中的哈密顿量：

$$
H_I(t) = -\mathbf{d}(t)\cdot\mathbf{E}(0,t),
$$

其中

$$
\mathbf{d}(t) = \mathbf{d}_{eg} e^{i\omega_0 t} |e\rangle\langle g| + \mathbf{d}_{ge} e^{-i\omega_0 t} |g\rangle\langle e|,
$$
$$
\mathbf{E}(0,t) = i\sum_{\mathbf{k},\lambda} C_k \left( a_{\mathbf{k},\lambda} \boldsymbol{\varepsilon}_{\mathbf{k},\lambda} e^{-i\omega_k t} - a_{\mathbf{k},\lambda}^\dagger \boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^* e^{i\omega_k t} \right),
$$
且 $C_k \equiv \sqrt{\hbar\omega_k/(2\varepsilon_0 V)}$。

设 $t=0$ 时系统处于 $H_0$ 的本征态 $|i\rangle$，则由一阶含时微扰论，跃迁到另一本征态 $|f\rangle$ 的振幅为

$$
c_f(t) = -\frac{i}{\hbar} \int_0^t dt' \, \langle f|H_I(t')|i\rangle.
$$

在长时极限下，跃迁概率呈时间线性增长，单位时间跃迁速率由费米黄金规则给出：

$$
W_{i\to f} = \frac{2\pi}{\hbar} \bigl| \langle f| H_I(0) | i \rangle \bigr|^2 \, \delta(E_f - E_i). \tag{2}
$$

其中 $H_I(0)$ 是剥离所有时间振荡因子后的“骨架”算符：

$$
H_I(0) = -\mathbf{d}(0)\cdot\mathbf{E}(0,0),
$$
$$
\mathbf{d}(0) = \mathbf{d}_{eg}|e\rangle\langle g| + \mathbf{d}_{ge}|g\rangle\langle e|,
$$
$$
\mathbf{E}(0,0) = i\sum_{\mathbf{k},\lambda} C_k \left( a_{\mathbf{k},\lambda} \boldsymbol{\varepsilon}_{\mathbf{k},\lambda} - a_{\mathbf{k},\lambda}^\dagger \boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^* \right).
$$

下面，我们将同一个 $H_I(0)$ 作用于三种不同的初态，得出三种过程的速率。它们将共享相同的耦合常数 $|\mathbf{d}_{eg}|$ 和相同的相空间结构，系数之间的关联由此诞生。

## 4. 受激吸收

**初态：** 原子处于基态 $|g\rangle$，某一特定模式 $(\mathbf{k},\lambda)$ 中有 $n_{\mathbf{k},\lambda}$ 个光子。  
$|i\rangle = |g, n_{\mathbf{k},\lambda}\rangle$，能量 $E_i = E_g + n\hbar\omega_k$。

**末态：** 原子跃迁到激发态 $|e\rangle$，该模式光子数减 1。  
$|f\rangle = |e, n_{\mathbf{k},\lambda}-1\rangle$，能量 $E_f = E_e + (n-1)\hbar\omega_k$。

能量守恒的 $\delta$ 函数要求 $\hbar\omega_k = E_e - E_g = \hbar\omega_0$，即 $\omega_k = \omega_0$，共振吸收。

要在矩阵元中得到非零贡献，必须从 $H_I(0)$ 中取出能够 **湮灭一个光子并将原子从 $|g\rangle$ 提升到 $|e\rangle$** 的项。审视 $\mathbf{E}(0,0)$：湮灭算符 $a_{\mathbf{k},\lambda}$ 的项带有系数 $i C_k \boldsymbol{\varepsilon}_{\mathbf{k},\lambda}$。同时，$\mathbf{d}(0)$ 中能使原子上升的是 $\mathbf{d}_{eg}|e\rangle\langle g|$。二者结合，给出该过程的等效算符：

$$
H_I^{\text{(abs)}}(0) = -\mathbf{d}_{eg}|e\rangle\langle g| \cdot \bigl(i C_k \boldsymbol{\varepsilon}_{\mathbf{k},\lambda} a_{\mathbf{k},\lambda}\bigr) = -i C_k (\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda})\, |e\rangle\langle g|\, a_{\mathbf{k},\lambda}.
$$

计算矩阵元，利用 $a|n\rangle = \sqrt{n}\,|n-1\rangle$，得到

$$
\langle f| H_I^{\text{(abs)}}(0) |i\rangle = -i C_k (\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}) \sqrt{n}.
$$

代入黄金规则 (2)，并注意 $C_k^2 = \hbar\omega_k/(2\varepsilon_0 V)$，以及 $\delta(E_f-E_i) = \delta(\hbar\omega_k - \hbar\omega_0) = \hbar^{-1}\delta(\omega_k - \omega_0)$，得

$$
W_{g\to e}^{\text{abs},(\mathbf{k},\lambda)} = \frac{2\pi}{\hbar} \cdot \frac{\hbar\omega_k}{2\varepsilon_0 V} |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}|^2 \, n_{\mathbf{k},\lambda} \cdot \frac{1}{\hbar}\delta(\omega_k - \omega_0)
= \frac{\pi\omega_k}{\varepsilon_0 V \hbar} |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}|^2 \, n_{\mathbf{k},\lambda} \, \delta(\omega_k - \omega_0). \tag{3}
$$

这正是单个模式贡献的吸收速率。我们暂时停在这里，先不去做模式求和，而是转向辐射过程，以便进行直接对比。

## 5. 受激辐射与自发辐射——同一顶角的另一面

现在考虑原子最初处于激发态 $|e\rangle$，同一模式中有 $n_{\mathbf{k},\lambda}$ 个光子。

**初态：** $|i\rangle = |e, n_{\mathbf{k},\lambda}\rangle$，能量 $E_i = E_e + n\hbar\omega_k$。  
**末态：** $|f\rangle = |g, n_{\mathbf{k},\lambda}+1\rangle$，能量 $E_f = E_g + (n+1)\hbar\omega_k$。

能量守恒再次要求 $\omega_k = \omega_0$。

这次，我们需要**产生一个光子并让原子退激发**。从 $\mathbf{E}(0,0)$ 中挑选含 $a^\dagger$ 的项： $-i C_k \boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^* a_{\mathbf{k},\lambda}^\dagger$。而 $\mathbf{d}(0)$ 中使原子下降的是 $\mathbf{d}_{ge}|g\rangle\langle e|$。因此相关的部分为

$$
H_I^{\text{(em)}}(0) = -\mathbf{d}_{ge}|g\rangle\langle e| \cdot \bigl(-i C_k \boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^* a_{\mathbf{k},\lambda}^\dagger\bigr) = i C_k (\mathbf{d}_{ge}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^*)\, |g\rangle\langle e|\, a_{\mathbf{k},\lambda}^\dagger.
$$

矩阵元为

$$
\langle f| H_I^{\text{(em)}}(0) |i\rangle = i C_k (\mathbf{d}_{ge}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^*) \langle n+1|a^\dagger|n\rangle = i C_k (\mathbf{d}_{ge}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}^*) \sqrt{n+1}.
$$

取模方，利用 $|\mathbf{d}_{ge}\!\cdot\!\boldsymbol{\varepsilon}^*| = |\mathbf{d}_{eg}^*\!\cdot\!\boldsymbol{\varepsilon}^*| = |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}|$，得到

$$
|\langle f| H_I(0) |i\rangle|^2 = C_k^2 \, |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}|^2 \, (n_{\mathbf{k},\lambda} + 1).
$$

**这便是整个故事中最具启示性的表达式。** 因子 $n+1$ 将辐射速率自然切割为两部分：与 $n$ 成正比的部分，是受激辐射；而与 $n$ 无关的常数 $1$，是自发辐射。两者来自同一矩阵元、同一顶角，区别仅在于初态光子数是否为零。

代入黄金规则，单模辐射速率为

$$
W_{e\to g}^{\text{em},(\mathbf{k},\lambda)} = \frac{2\pi}{\hbar} C_k^2 |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}|^2 (n_{\mathbf{k},\lambda} + 1) \, \delta(\hbar\omega_k - \hbar\omega_0)
= \frac{\pi\omega_k}{\varepsilon_0 V \hbar} |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}|^2 (n_{\mathbf{k},\lambda}+1) \, \delta(\omega_k - \omega_0). \tag{4}
$$

对比 (3) 式与 (4) 式，我们可以看到：除去因子 $n$ 与 $n+1$ 的差异，两个速率的系数**完全一致**。若我们只关心受激部分（即正比于 $n$ 的项），则吸收与受激辐射的单模速率在模密度和偏振权重上完全对称。这就是时间反演对称性在跃迁振幅层面的直接表达：一个光子参与的过程与其逆过程，有相同的概率幅模方。这种对称性，与系统是否处于热平衡毫无关系，是哈密顿量厄米性内禀的结果。

## 6. 模式求和与爱因斯坦系数的诞生

实际原子面对的是连续分布的大量电磁模式。总跃迁速率需对所有模式求和。对于热平衡或宽带辐射场，不同模式间不存在相干，我们只需将各模式速率相加：

$$
W_{g\to e}^{\text{total}} = \sum_{\mathbf{k},\lambda} W_{g\to e}^{\text{abs},(\mathbf{k},\lambda)}, \qquad 
W_{e\to g}^{\text{total}} = \sum_{\mathbf{k},\lambda} W_{e\to g}^{\text{em},(\mathbf{k},\lambda)}.
$$

取热力学极限，将求和转化为积分：

$$
\sum_{\mathbf{k}} \longrightarrow \frac{V}{(2\pi)^3} \int d^3k.
$$

利用球坐标 $d^3k = k^2 dk d\Omega$，且 $k = \omega/c$，$dk = d\omega/c$，故 $d^3k = \frac{\omega^2}{c^3} d\omega d\Omega$。于是

$$
\sum_{\mathbf{k},\lambda} \longrightarrow \frac{V}{(2\pi)^3} \sum_{\lambda=1}^2 \int d\Omega \int_0^\infty d\omega \frac{\omega^2}{c^3}.
$$

先处理吸收。将 (3) 式代入并积分，

$$
W_{g\to e}^{\text{total}} = \frac{V}{(2\pi)^3} \sum_\lambda \int d\Omega \int_0^\infty d\omega \frac{\omega^2}{c^3} \cdot \frac{\pi\omega}{\varepsilon_0 V \hbar} |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}|^2 \, n(\omega) \, \delta(\omega - \omega_0).
$$

体积 $V$ 精确抵消——这是量子化体积作为中间工具不留下痕迹的体现。利用 $\delta$ 函数对频率积分，得

$$
W_{g\to e}^{\text{total}} = \frac{1}{8\pi^2 \varepsilon_0 \hbar c^3} \int d\Omega \sum_\lambda |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}|^2 \;\cdot\; \omega_0^3 \, n(\omega_0).
$$

其中偏振与角度求和的结果是（详见附录）

$$
\int d\Omega \sum_{\lambda=1}^2 |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}|^2 = \frac{8\pi}{3} |\mathbf{d}_{eg}|^2.
$$

因此

$$
W_{g\to e}^{\text{total}} = \frac{1}{8\pi^2 \varepsilon_0 \hbar c^3} \cdot \frac{8\pi}{3} |\mathbf{d}_{eg}|^2 \, \omega_0^3 \, n(\omega_0) = \frac{\omega_0^3 |\mathbf{d}_{eg}|^2}{3\pi\varepsilon_0 \hbar c^3} \, n(\omega_0). \tag{5}
$$

该式表明，吸收速率正比于光子数 $n(\omega_0)$。回想光场能量密度的定义：在单位体积、单位角频率间隔内，电磁模式数（含两个偏振）为

$$
g(\omega) = \frac{\omega^2}{\pi^2 c^3},
$$

每个模式平均能量为 $\hbar\omega\, n(\omega)$，故能量密度

$$
\rho(\omega) = \hbar\omega \, g(\omega) \, n(\omega) = \frac{\hbar\omega^3}{\pi^2 c^3} n(\omega). \tag{6}
$$

由此解出 $n(\omega_0) = \frac{\pi^2 c^3}{\hbar\omega_0^3} \rho(\omega_0)$，代入 (5) 式：

$$
W_{g\to e}^{\text{total}} = \frac{\omega_0^3 |\mathbf{d}_{eg}|^2}{3\pi\varepsilon_0 \hbar c^3} \cdot \frac{\pi^2 c^3}{\hbar\omega_0^3} \rho(\omega_0) = \frac{\pi |\mathbf{d}_{eg}|^2}{3\varepsilon_0 \hbar^2} \, \rho(\omega_0).
$$

依照爱因斯坦的记号，吸收速率写为 $W_{g\to e} = B_{12} \, \rho(\omega_0)$，于是我们读出了受激吸收系数：

$$
B_{12} = \frac{\pi |\mathbf{d}_{eg}|^2}{3\varepsilon_0 \hbar^2}. \tag{7}
$$

现在处理辐射过程。将 (4) 式作同样的模式求和：

$$
W_{e\to g}^{\text{total}} = \frac{V}{(2\pi)^3} \sum_\lambda \int d\Omega \int d\omega \frac{\omega^2}{c^3} \cdot \frac{\pi\omega}{\varepsilon_0 V \hbar} |\mathbf{d}_{eg}\!\cdot\!\boldsymbol{\varepsilon}|^2 (n(\omega)+1) \delta(\omega-\omega_0)
= \frac{\omega_0^3 |\mathbf{d}_{eg}|^2}{3\pi\varepsilon_0 \hbar c^3} \bigl( n(\omega_0) + 1 \bigr).
$$

结果清晰地分为两项：

$$
W_{e\to g}^{\text{total}} = \underbrace{\frac{\omega_0^3 |\mathbf{d}_{eg}|^2}{3\pi\varepsilon_0 \hbar c^3} n(\omega_0)}_{\text{受激辐射}} \;+\; \underbrace{\frac{\omega_0^3 |\mathbf{d}_{eg}|^2}{3\pi\varepsilon_0 \hbar c^3}}_{\text{自发辐射}}.
$$

第一项用 (6) 式转换为能量密度，得到受激辐射速率 $B_{21} \rho(\omega_0)$，其中

$$
B_{21} = \frac{\pi |\mathbf{d}_{eg}|^2}{3\varepsilon_0 \hbar^2} = B_{12}. \tag{8}
$$

第二项即为自发辐射速率，记作

$$
A_{21} = \frac{\omega_0^3 |\mathbf{d}_{eg}|^2}{3\pi\varepsilon_0 \hbar c^3}. \tag{9}
$$

对比 $A_{21}$ 与 $B_{21}$，立即得到爱因斯坦第一关系：

$$
\frac{A_{21}}{B_{21}} = \frac{\hbar\omega_0^3}{\pi^2 c^3}. \tag{10}
$$

若改用频率 $\nu = \omega/2\pi$，并利用 $h = 2\pi\hbar$，则上式化为教科书中的经典形式：

$$
\frac{A_{21}}{B_{21}} = \frac{8\pi h \nu^3}{c^3}.
$$

至此，三个爱因斯坦系数全部由单一的原子参数——跃迁偶极矩阵元 $|\mathbf{d}_{eg}|$ 表达出来，且之间的关系成为必然，而非假定。

## 7. 为什么系数不独立：三项本质相同

上述数学推导中蕴含着深刻的物理，我们将其提炼为三个核心观点，逐一批注。

### 7.1 同一顶角与微观可逆性

整个推导的起点与终点都只有 **一个** 相互作用哈密顿量：

$$
H_I = -\mathbf{d}\cdot\mathbf{E}.
$$

我们从未针对吸收、受激辐射或自发辐射引入不同的耦合项。吸收过程抽取的是 $\mathbf{E}$ 中的湮灭算符 $a$，受激辐射和自发辐射抽取的是产生算符 $a^\dagger$。这两个算符恰为厄米共轭，而原子部分的矩阵元在两种过程（上升与下降）中互为复共轭。因此，

$$
|\langle f|H_I|i\rangle|_{\text{吸收}} = |\langle i|H_I|f\rangle|_{\text{辐射}}.
$$

跃迁概率的模方严格相等，导致 $B_{12} = B_{21}$（非简并情形）。**这是一种内蕴于量子力学的微观可逆性，它不依赖于系统是否达到热平衡，也不要求细致平衡条件。** 在相对论性 QED 中，这进一步表述为交叉对称性：将吸收过程的入射光子“交叉”到末态成为出射光子，振幅解析延拓后模方不变。爱因斯坦当初是通过热平衡下与普朗克公式对比才“推导”出这个关系，而在量子场论中，它是理论结构的内在对称性，热力学不过是其宏观展现。

### 7.2 自发辐射：真空涨落诱导的受激辐射

(4) 式中的因子 $n+1$ 是量子电动力学最漂亮的馈赠之一。它告诉我们，即使在 $n=0$ 的完全黑暗环境中，激发态原子仍然会衰变——速率正比于那个孤零零的“1”。这一项的起源，是产生算符对真空态的作用：$a^\dagger|0\rangle = \sqrt{1}|1\rangle$。

我们究竟该如何理解这个“1”？在经典图像里，真空就是什么都没有，何来驱动？量子场论的回答是：**真空不空**。每一个电磁模式都存在不可消除的零点振荡，其能量为 $\frac{1}{2}\hbar\omega$。这种零点涨落构成一个遍在的、各向同性的“噪声场”。当原子处于激发态时，这个真空噪声便扮演了与真实光子完全相同的角色，诱导原子向下跃迁并释放出一个真实光子。用公式表达，自发辐射速率正是受激辐射系数乘以真空的有效模式密度：

$$
A_{21} = B_{21} \times \frac{\hbar\omega_0^3}{\pi^2 c^3}.
$$

而因子 $\frac{\hbar\omega_0^3}{\pi^2 c^3}$ 恰好是 $\hbar\omega_0 \, g(\omega_0)$，其中 $g(\omega_0) = \frac{\omega_0^2}{\pi^2 c^3}$ 是单位角频率间隔内的模式数（每单位体积）。从计算过程看，这个因子正是来自**对真空里所有可以接纳光子的模式进行相空间积分**——即“末态求和”。自发辐射多出的那个相空间积分，就是将光子发射到任意方向、任意偏振的所有可能末态进行遍历求和的结果。换句话说，自发辐射不是某种神秘的独立机制，它就是**在所有真空模式中同时发生的、由真空涨落诱导的受激辐射**。

### 7.3 爱因斯坦唯象理论中的独立性问题

在爱因斯坦 1917 年的推导中，他做出了一个关键假设：处于不同能级的原子，它们对辐射的吸收和发射过程是彼此独立的，因此总反应速率可以简单地正比于各能级的原子数。这相当于将原子与辐射场之间的相互作用视为泊松过程，不计原子间的关联，也不考虑原子对辐射场的反作用导致的相干效应。

在量子场论视角下，这一假设对应于什么？它对应于我们处理辐射场时，将不同电磁模式视为自由的、无相互作用的玻色子气体，并且忽略了原子通过辐射场媒介而产生的有效相互作用（例如超辐射、偶极-偶极耦合）。在我们的推导中，我们正是这样做的：我们把初态和末态均取为光子数态（Fock 态），假设各模式独立演化，原子跃迁只改变单个模式的光子数，然后将所有模式的贡献线性叠加。这就**严格证明了热平衡（或宽带非相干场）条件下，原子吸收和发射的独立性确实成立**，但它并非原理，而是一阶微扰和自由场近似的结果。一旦辐射场具有高度相干性，或原子密度足够高，这种独立性就会瓦解，我们就会看到拉比振荡、集体自发辐射等非独立现象。爱因斯坦的唯象框架之所以能成功，恰恰是因为普通热光源满足这种独立性条件。

因此，三个爱因斯坦系数看似是经验输入，实则是量子电动力学中同一个顶角在不同光子统计背景下的三种表现。它们的非独立性，不是巧合，而是“物质-光只有一种基本相互作用”这一事实的同义反复。

## 8. 推广至简并情形

上述推导假设了上下能级非简并。若低能级有 $g_1$ 重简并，高能级有 $g_2$ 重简并，我们需要对各子能级的跃迁求和。吸收系数 $B_{12}$ 通常定义为对所有低能级子能级平均、高能级求和；受激辐射系数 $B_{21}$ 则相反。对于每一对具体子能级，其跃迁矩阵元仍然满足 $|\mathbf{d}_{\alpha\beta}| = |\mathbf{d}_{\beta\alpha}|$。将所有通道的贡献相加，不难得出

$$
g_1 B_{12} = g_2 B_{21}.
$$

此关系与热平衡的玻尔兹曼因子完全兼容，但其根源依然是量子力学基本顶角的对称性，与温度无关。

## 9. 结语

我们从真空中电磁场的量子化出发，沿着最朴素的量子力学含时微扰论路径，亲手算出了原子与光相互作用的三条基本通道的速率。三条通道最后汇入完全相同的偶极矩阵元 $|\mathbf{d}_{eg}|$ 和相同的结构因子。

---

## 附录：偏振求和与角向积分

我们需要计算

$$
I = \int d\Omega \sum_{\lambda=1}^2 |\mathbf{d}\cdot\boldsymbol{\varepsilon}_{\mathbf{k},\lambda}|^2.
$$

对固定波矢 $\mathbf{k}$，两个偏振矢量 $\boldsymbol{\varepsilon}_{\mathbf{k},1}$、$\boldsymbol{\varepsilon}_{\mathbf{k},2}$ 和单位矢量 $\hat{\mathbf{k}} = \mathbf{k}/|\mathbf{k}|$ 构成一组正交归一基底。投影算符 $\sum_\lambda \varepsilon_\lambda^i \varepsilon_\lambda^j = \delta^{ij} - \hat{k}^i\hat{k}^j$。因此

$$
\sum_\lambda |\mathbf{d}\cdot\boldsymbol{\varepsilon}|^2 = d_i^* d_j (\delta^{ij} - \hat{k}^i\hat{k}^j) = |\mathbf{d}|^2 - |\mathbf{d}\cdot\hat{\mathbf{k}}|^2.
$$

取 $\mathbf{d}$ 沿 $z$ 轴方向，大小为 $d$，则 $|\mathbf{d}\cdot\hat{\mathbf{k}}|^2 = d^2 \cos^2\theta$，其中 $\theta$ 为极角。于是

$$
I = \int_0^{2\pi} d\phi \int_0^\pi \sin\theta d\theta \; d^2 (1 - \cos^2\theta)
= 2\pi d^2 \int_{-1}^1 (1 - u^2) du = 2\pi d^2 \left[ u - \frac{u^3}{3} \right]_{-1}^1
= 2\pi d^2 \left(2 - \frac{2}{3}\right) = \frac{8\pi}{3} |\mathbf{d}|^2.
$$

这一干净的结果直接用于文中总速率的化简。