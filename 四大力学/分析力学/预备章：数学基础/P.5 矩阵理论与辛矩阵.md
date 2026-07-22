---
created: 2026-07-18T03:15:26+08:00
modified: 2026-07-18T13:59:53+08:00
zhihu-title: 矩阵理论
zhihu-topics:
  - 分析力学
zhihu-link: https://zhuanlan.zhihu.com/p/2061652505495352454
zhihu-toc: true
zhihu-created-at: 2026-07-18 13:59
zhihu-updated-at: 2026-07-21 14:49
---
<!-- zhihu-nav-top:start -->
**上一节：** [P.4 常微分方程与特殊函数](https://zhuanlan.zhihu.com/p/2061652524868821844)
<!-- zhihu-nav-top:end -->

## P.5 矩阵理论与辛矩阵

分析力学从拉格朗日力学过渡到哈密顿力学后，相空间取代了位形空间，正则方程取代了二阶拉格朗日方程。在这个新框架中，矩阵理论——特别是辛矩阵——扮演着核心角色。小振动理论需要本征值问题和对角化来求简正频率和简正模式；二次型描述动能和势能；而辛矩阵则刻画了哈密顿相流和正则变换的本质结构。这一节从本征值问题的几何意义讲起，建立矩阵对角化与二次型的关系，然后引入辛矩阵的定义、性质及其与哈密顿力学的深刻联系。所有推导从零开始，完整写出。

---

### 一、线性变换的几何视角

矩阵不是一堆数字的随意排列。一个 $n \times n$ 实矩阵 $A$ 代表 $\mathbb{R}^n$ 上的一个**线性变换**：

$$
\mathbf{x} \mapsto \mathbf{y} = A\mathbf{x}, \quad y_i = \sum_{j=1}^n A_{ij} x_j.
$$

线性意味着两条：$A(\mathbf{x}_1 + \mathbf{x}_2) = A\mathbf{x}_1 + A\mathbf{x}_2$，$A(c\mathbf{x}) = c A\mathbf{x}$。

在几何上，这个变换可以拉伸、旋转、反射、剪切向量。理解一个线性变换的核心问题是：是否存在某些特殊方向，变换后方向不变，只被拉伸或压缩？这就是本征值问题。

---

### 二、本征值问题与本征向量

#### 2.1 定义与特征方程

设 $A$ 是 $n \times n$ 矩阵（实或复）。若存在非零向量 $\mathbf{v} \neq \mathbf{0}$ 和标量 $\lambda$，使得

$$
A \mathbf{v} = \lambda \mathbf{v}, \tag{P.5.1}
$$

则称 $\lambda$ 为 $A$ 的**本征值**（eigenvalue），$\mathbf{v}$ 为对应的**本征向量**（eigenvector）。

几何意义：在 $\mathbf{v}$ 方向上，线性变换 $A$ 等同于数乘 $\lambda$。这个方向是 $A$ 的“不变方向”。

改写为

$$
(A - \lambda I) \mathbf{v} = \mathbf{0}.
$$

这是一个齐次线性方程组。它有非零解 $\mathbf{v} \neq \mathbf{0}$ 的充要条件是系数矩阵不可逆，即行列式为零：

$$
\boxed{\det(A - \lambda I) = 0}. \tag{P.5.2}
$$

这就是**特征方程**。左边是一个关于 $\lambda$ 的 $n$ 次多项式——**特征多项式**：

$$
p(\lambda) = \det(A - \lambda I) = (-1)^n \lambda^n + c_{n-1} \lambda^{n-1} + \cdots + c_1 \lambda + c_0.
$$

代数基本定理保证在复数域中有 $n$ 个根（计重数），这就是 $A$ 的本征值。记本征值集合为 $\{\lambda_1, \lambda_2, \ldots, \lambda_n\}$（称为谱）。

#### 2.2 一个完整计算的例子

取 $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。

特征多项式：
$$
\det\begin{pmatrix} 2-\lambda & 1 \\ 1 & 2-\lambda \end{pmatrix} = (2-\lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda - 1)(\lambda - 3).
$$

本征值 $\lambda_1 = 1, \lambda_2 = 3$。

对于 $\lambda_1 = 1$，解 $(A - I)\mathbf{v} = \mathbf{0}$：
$$
\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \Rightarrow v_1 + v_2 = 0.
$$
取 $\mathbf{v}_1 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$（或任意非零倍数）。

对于 $\lambda_2 = 3$，解 $(A - 3I)\mathbf{v} = \mathbf{0}$：
$$
\begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix} \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \Rightarrow -v_1 + v_2 = 0.
$$
取 $\mathbf{v}_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。

几何上，$A$ 在 $(1, -1)$ 方向上保持分量不变（本征值 $1$），在 $(1, 1)$ 方向上拉伸 $3$ 倍（本征值 $3$）。

#### 2.3 本征值的性质

**迹与本征值之和**：
$$
\text{Tr}(A) = \sum_{i=1}^n A_{ii} = \sum_{i=1}^n \lambda_i. \tag{P.5.3}
$$

证明思路：特征多项式展开中 $\lambda^{n-1}$ 的系数等于 $-\text{Tr}(A)$，同时也等于 $-\sum \lambda_i$。

**行列式与本征值之积**：
$$
\det(A) = \prod_{i=1}^n \lambda_i. \tag{P.5.4}
$$

证明：在 $\det(A - \lambda I)$ 中令 $\lambda = 0$，得 $\det(A) = p(0) = \prod (-\lambda_i) \times (-1)^n = \prod \lambda_i$。

**实矩阵的复本征值成对共轭**：若 $A$ 是实矩阵，$\lambda$ 是本征值且 $\mathbf{v}$ 是本征向量，则取复共轭得 $A\mathbf{v}^* = \lambda^* \mathbf{v}^*$，所以 $\lambda^*$ 也是本征值，$\mathbf{v}^*$ 是相应的本征向量。

#### 2.4 相似变换与本征值不变性

若 $P$ 是可逆矩阵，矩阵 $B = P^{-1} A P$ 称为 $A$ 的**相似变换**。相似矩阵有相同的本征值：

$$
\begin{aligned}
\det(B - \lambda I) &= \det(P^{-1} A P - \lambda P^{-1} P) = \det(P^{-1}(A - \lambda I) P) \\
&= \det(P^{-1}) \det(A - \lambda I) \det(P) = \det(A - \lambda I).
\end{aligned}
$$

本征值是相似不变量。对角化就是找一个相似变换，把 $A$ 变成最简单的形式——对角矩阵。

---

### 三、矩阵对角化

#### 3.1 对角化的条件

如果 $A$ 有 $n$ 个线性无关的本征向量 $\mathbf{v}_1, \ldots, \mathbf{v}_n$（对应的本征值为 $\lambda_1, \ldots, \lambda_n$），将这 $n$ 个本征向量排成列构成矩阵 $P$：

$$
P = \begin{pmatrix} | & | & & | \\ \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\ | & | & & | \end{pmatrix}.
$$

则

$$
A P = A \begin{pmatrix} \mathbf{v}_1 & \cdots & \mathbf{v}_n \end{pmatrix} = \begin{pmatrix} A\mathbf{v}_1 & \cdots & A\mathbf{v}_n \end{pmatrix} = \begin{pmatrix} \lambda_1 \mathbf{v}_1 & \cdots & \lambda_n \mathbf{v}_n \end{pmatrix}.
$$

另一方面，

$$
P \begin{pmatrix} \lambda_1 & & \\ & \ddots & \\ & & \lambda_n \end{pmatrix} = \begin{pmatrix} \lambda_1 \mathbf{v}_1 & \cdots & \lambda_n \mathbf{v}_n \end{pmatrix}.
$$

因此

$$
A P = P \Lambda, \quad \Lambda = \text{diag}(\lambda_1, \ldots, \lambda_n).
$$

由于 $\mathbf{v}_i$ 线性无关，$P$ 可逆。两边同时右乘 $P^{-1}$：

$$
\boxed{A = P \Lambda P^{-1}}, \quad \text{或} \quad \boxed{P^{-1} A P = \Lambda}. \tag{P.5.5}
$$

对角化的充要条件是：$A$ 有 $n$ 个线性无关的本征向量。

**充分条件**（更易验证）：如果 $A$ 的所有 $n$ 个本征值互不相同，则对应的 $n$ 个本征向量自动线性无关—— $A$ 可对角化。

**注意**：即使有重本征值，仍可能对角化（如单位矩阵 $I$，本征值全为 $1$，但任何非零向量都是本征向量，存在完备的本征向量系）。也可能不可对角化——此时只能化为若尔当标准型（见后）。

#### 3.2 对称矩阵的对角化——正交对角化

物理中最重要的矩阵类之一是对称矩阵 $A = A^T$。实对称矩阵有两个核心定理：

**定理一**：实对称矩阵的所有本征值都是实数。

**证明**：设 $A\mathbf{v} = \lambda \mathbf{v}$，$\mathbf{v} \neq \mathbf{0}$。取内积：
$$
\lambda \|\mathbf{v}\|^2 = \langle \mathbf{v}, A\mathbf{v} \rangle = \langle A\mathbf{v}, \mathbf{v} \rangle = \lambda^* \|\mathbf{v}\|^2.
$$
因为 $A$ 对称，$\langle \mathbf{v}, A\mathbf{v} \rangle = \mathbf{v}^\dagger A \mathbf{v} = (A\mathbf{v})^\dagger \mathbf{v} = \langle A\mathbf{v}, \mathbf{v} \rangle$。所以 $\lambda = \lambda^*$，即 $\lambda \in \mathbb{R}$。

**定理二**：实对称矩阵不同本征值对应的本征向量正交。

**证明**：设 $A\mathbf{v}_1 = \lambda_1 \mathbf{v}_1$，$A\mathbf{v}_2 = \lambda_2 \mathbf{v}_2$，$\lambda_1 \neq \lambda_2$。则
$$
\lambda_1 \langle \mathbf{v}_2, \mathbf{v}_1 \rangle = \langle \mathbf{v}_2, A\mathbf{v}_1 \rangle = \langle A\mathbf{v}_2, \mathbf{v}_1 \rangle = \lambda_2 \langle \mathbf{v}_2, \mathbf{v}_1 \rangle.
$$
因此 $(\lambda_1 - \lambda_2) \langle \mathbf{v}_2, \mathbf{v}_1 \rangle = 0$，而 $\lambda_1 \neq \lambda_2$，故 $\langle \mathbf{v}_2, \mathbf{v}_1 \rangle = 0$。

结合这两个定理：实对称矩阵可以找到一组**正交归一**的本征向量基 $\{\mathbf{u}_1, \ldots, \mathbf{u}_n\}$（对于重本征值，在特征子空间中做 Gram-Schmidt 正交化）。将这些正交归一基排成正交矩阵 $Q = (\mathbf{u}_1, \ldots, \mathbf{u}_n)$（$Q^T Q = I$，即 $Q^{-1} = Q^T$），则

$$
\boxed{A = Q \Lambda Q^T}, \quad Q^T Q = I. \tag{P.5.6}
$$

这就是**正交对角化**——实对称矩阵可以通过正交变换对角化。

#### 3.3 若尔当标准型简介（不可对角化的情况）

若 $A$ 没有足够的线性无关本征向量，则无法对角化。此时最简形式是**若尔当标准型**（Jordan normal form）。存在可逆矩阵 $P$ 使得

$$
P^{-1} A P = J = \begin{pmatrix}
J_{n_1}(\lambda_1) & & \\
& \ddots & \\
& & J_{n_k}(\lambda_k)
\end{pmatrix},
$$

其中每个若尔当块为

$$
J_m(\lambda) = \begin{pmatrix}
\lambda & 1 & & \\
& \lambda & \ddots & \\
& & \ddots & 1 \\
& & & \lambda
\end{pmatrix}_{m \times m}.
$$

$1$ 在上次对角线的位置表示本征向量的缺失——对应的是广义本征向量。在分析力学的小振动理论中，我们通常处理的是对称矩阵（动能和势能矩阵），因此总是可以对角化。若尔当型在耗散系统和某些微扰问题中可能出现。

---

### 四、二次型

#### 4.1 定义与矩阵表示

**二次型**是 $n$ 个变量的二次齐次多项式：

$$
Q(\mathbf{x}) = \sum_{i,j=1}^n A_{ij} x_i x_j,
$$

其中 $A_{ij}$ 可以取为对称的（因为 $x_i x_j = x_j x_i$，总可令 $A_{ij} = A_{ji}$）。用矩阵表示：

$$
\boxed{Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}}, \quad A = A^T. \tag{P.5.7}
$$

在分析力学中，动能就是广义速度的二次型：$T = \frac{1}{2} \dot{\mathbf{q}}^T M \dot{\mathbf{q}}$（$M$ 是质量矩阵，对称正定）。势能在平衡点附近展开也是二次型：$V \approx \frac{1}{2} \mathbf{q}^T K \mathbf{q}$（$K$ 是刚度矩阵，对称）。

#### 4.2 通过坐标变换化简二次型

对坐标做线性变换 $\mathbf{x} = P \mathbf{y}$（$P$ 可逆），则

$$
Q = \mathbf{x}^T A \mathbf{x} = (P\mathbf{y})^T A (P\mathbf{y}) = \mathbf{y}^T (P^T A P) \mathbf{y}.
$$

因此二次型的矩阵按**合同变换** $A \to P^T A P$ 变化。合同变换保持对称性。

**目标**：找一个 $P$，使 $P^T A P$ 尽可能简单——对角形。由于 $A$ 是实对称矩阵，它可以正交对角化：$Q^T A Q = \Lambda = \text{diag}(\lambda_1, \ldots, \lambda_n)$。令 $\mathbf{x} = Q \mathbf{y}$，则

$$
Q = \mathbf{x}^T A \mathbf{x} = \mathbf{y}^T \Lambda \mathbf{y} = \sum_{i=1}^n \lambda_i y_i^2. \tag{P.5.8}
$$

这就是**主轴定理**（principal axis theorem）：任何实二次型都可以通过正交变换化为平方和（无交叉项），系数正是 $A$ 的本征值。

#### 4.3 惯性定律与正定性

**西尔维斯特惯性定律**：虽然合同变换 $P^T A P$ 的具体形式依赖于 $P$ 的选择，但正本征值的个数 $n_+$、负本征值的个数 $n_-$ 和零本征值的个数 $n_0$ 是合同的**不变量**。三元组 $(n_+, n_-, n_0)$ 称为二次型的惯性指数。

**正定性**：若对所有 $\mathbf{x} \neq \mathbf{0}$，$\mathbf{x}^T A \mathbf{x} > 0$，则 $A$ 是**正定矩阵**（positive definite）。等价条件：
- 所有本征值 $\lambda_i > 0$。
- 所有顺序主子式 $> 0$（西尔维斯特判据）。
- 可做 Cholesky 分解 $A = L L^T$（$L$ 下三角可逆）。

类似定义半正定（$\geq 0$）、负定（$<0$）、不定（有正有负）。动能矩阵 $M$ 是正定的，因为它正比于 $\dot{\mathbf{q}}^2$ 的期望值。

#### 4.4 同时对角化两个二次型

小振动理论的核心问题是同时对角化 $T$ 和 $V$：

$$
T = \frac{1}{2} \dot{\mathbf{q}}^T M \dot{\mathbf{q}}, \quad V = \frac{1}{2} \mathbf{q}^T K \mathbf{q}.
$$

我们需要找一个坐标变换，使 $T$ 变为单位矩阵（各模动能解耦）而 $V$ 变为对角矩阵（各模势能解耦）。这等价于求解广义本征值问题：

$$
\boxed{K \mathbf{v} = \lambda M \mathbf{v}}. \tag{P.5.9}
$$

**推导**：由于 $M$ 正定，可做 Cholesky 分解 $M = L L^T$。令 $\mathbf{q} = (L^T)^{-1} \mathbf{u}$，则

$$
T = \frac{1}{2} \dot{\mathbf{u}}^T \dot{\mathbf{u}}, \quad V = \frac{1}{2} \mathbf{u}^T (L^{-1} K (L^T)^{-1}) \mathbf{u} \equiv \frac{1}{2} \mathbf{u}^T \tilde{K} \mathbf{u}.
$$

$\tilde{K}$ 是实对称矩阵，可正交对角化：$Q^T \tilde{K} Q = \Omega^2 = \text{diag}(\omega_1^2, \ldots, \omega_n^2)$。令 $\mathbf{u} = Q \boldsymbol{\eta}$，则

$$
T = \frac{1}{2} \dot{\boldsymbol{\eta}}^T \dot{\boldsymbol{\eta}}, \quad V = \frac{1}{2} \boldsymbol{\eta}^T \Omega^2 \boldsymbol{\eta}.
$$

在简正坐标 $\boldsymbol{\eta}$ 下，$n$ 个自由度完全解耦，每个满足 $\ddot{\eta}_i + \omega_i^2 \eta_i = 0$。这一整套步骤在分析力学的小振动章节中会完整展开，这里先建立数学框架。

---

### 五、辛矩阵的定义与性质

哈密顿力学的相空间是偶数维的：$2n$ 维，由 $n$ 个广义坐标 $q_i$ 和 $n$ 个广义动量 $p_i$ 组成。这个空间的几何不是欧几里得的，而是**辛的**（symplectic）。辛矩阵是保持这个辛几何结构的线性变换。

#### 5.1 辛内积与辛矩阵的定义

将 $2n$ 维相空间向量记作 $\boldsymbol{\xi} = (q_1, \ldots, q_n, p_1, \ldots, p_n)^T$。定义**辛矩阵** $J$：

$$
\boxed{J = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}}, \tag{P.5.10}
$$

其中 $I_n$ 是 $n \times n$ 单位矩阵，$0$ 是 $n \times n$ 零矩阵。$J$ 是 $2n \times 2n$ 矩阵。

$J$ 的基本性质：
$$
J^T = -J, \quad J^2 = -I_{2n}, \quad J^{-1} = -J = J^T. \tag{P.5.11}
$$

$J$ 定义了一个反对称双线性型——**辛内积**：

$$
\omega(\boldsymbol{\xi}, \boldsymbol{\eta}) = \boldsymbol{\xi}^T J \boldsymbol{\eta} = \sum_{i=1}^n (q_i^\xi p_i^\eta - p_i^\xi q_i^\eta). \tag{P.5.12}
$$

对于单个向量，它可定义一个“辛面积”（二维相面积之和）。

**辛矩阵的定义**：一个 $2n \times 2n$ 实矩阵 $M$ 称为**辛矩阵**（symplectic matrix），当且仅当

$$
\boxed{M^T J M = J}. \tag{P.5.13}
$$

这个条件等价于：$M$ 保持辛内积不变，即对任意 $\boldsymbol{\xi}, \boldsymbol{\eta}$，
$$
\omega(M\boldsymbol{\xi}, M\boldsymbol{\eta}) = \omega(\boldsymbol{\xi}, \boldsymbol{\eta}).
$$

因为 $(M\boldsymbol{\xi})^T J (M\boldsymbol{\eta}) = \boldsymbol{\xi}^T (M^T J M) \boldsymbol{\eta} = \boldsymbol{\xi}^T J \boldsymbol{\eta}$。

#### 5.2 辛矩阵的基本性质

**性质一：辛矩阵构成群**

所有 $2n \times 2n$ 辛矩阵的集合记为 $\text{Sp}(2n, \mathbb{R})$（实辛群）。

- 单位元 $I_{2n}$ 是辛矩阵：$I^T J I = J$。
- 若 $M_1, M_2$ 是辛矩阵，则 $M_1 M_2$ 也是：$(M_1 M_2)^T J (M_1 M_2) = M_2^T (M_1^T J M_1) M_2 = M_2^T J M_2 = J$。
- 若 $M$ 是辛矩阵，则 $M^{-1}$ 也是：由 $M^T J M = J$，两边左乘 $(M^T)^{-1}$、右乘 $M^{-1}$ 得 $J = (M^T)^{-1} J M^{-1}$。注意 $(M^{-1})^T = (M^T)^{-1}$，所以 $(M^{-1})^T J M^{-1} = J$。

因此 $\text{Sp}(2n, \mathbb{R})$ 是 $\text{GL}(2n, \mathbb{R})$（一般线性群）的子群。

**性质二：行列式为 $+1$ **

$$
\boxed{\det M = 1}. \tag{P.5.14}
$$

**证明**：由 $M^T J M = J$，取行列式：
$$
\det(M^T) \det(J) \det(M) = \det(J).
$$
$\det(M^T) = \det(M)$，所以 $(\det M)^2 \det J = \det J$。只要证明 $\det J \neq 0$（确实 $\det J = (\det I_n)^2 = 1$，因为 $J$ 通过置换行列可化为 $\begin{pmatrix} I_n & 0 \\ 0 & I_n \end{pmatrix}$ 加一个符号——准确说，$\det J = 1$，因为 $J = \begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}$ 经过 $n$ 次交换行可将 $-I$ 移到左上，$I$ 移到右下，符号为 $(-1)^n$，同时 $-I$ 贡献 $(-1)^n$，总乘积为 $1$）。因此 $(\det M)^2 = 1$，$\det M = \pm 1$。进一步可证 $\det M = +1$（辛矩阵恒为正行列式，这是辛群连通性的结果—— $\text{Sp}(2n,\mathbb{R})$ 是连通的，而 $\det$ 连续，在 $I$ 处为 $1$，故整体为 $1$）。

**性质三：本征值的倒数对称性**

若 $\lambda$ 是辛矩阵 $M$ 的本征值，则 $1/\lambda$ 也是本征值。

**证明**：由 $M^T J M = J$ 及 $J^{-1} = -J$，有 $M^{-1} = -J M^T J$。特征多项式：
$$
\det(M - \lambda I) = \det( -J M^T J - \lambda I) \quad (\text{注意这需要谨慎使用})。
$$

更直接的方法：从 $M^T J M = J$ 可得 $M$ 与 $(M^T)^{-1}$ 相似（用 $J$ 做相似变换）。实际上：
$$
J M J^{-1} = -J M J = -J M (-J^{-1})? 
$$
让我们仔细推导：$M^T J M = J \Rightarrow M^T = J M^{-1} J^{-1}$。取逆：$(M^T)^{-1} = (J M^{-1} J^{-1})^{-1} = J M J^{-1}$。但 $(M^T)^{-1} = (M^{-1})^T$。所以 $M^{-1}$ 的转置相似于 $M$。而 $M^{-1}$ 与本征值 $\lambda^{-1}$ 对应，且 $M^T$ 与 $M$ 有相同的本征值（特征多项式相同）。因此若 $\lambda$ 是本征值，$\lambda^{-1}$ 也是本征值。

更精确的结论：特征多项式是**回文**的（palindromic）：$p(\lambda) = \lambda^{2n} p(1/\lambda)$（对于适当的归一化）。

由此推出：本征值以四元组 $(\lambda, \lambda^{-1}, \bar{\lambda}, \bar{\lambda}^{-1})$ 出现（考虑实矩阵的复共轭对称性）。

**性质四：辛矩阵的分块表示**

将 $M$ 写为 $n \times n$ 分块：

$$
M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}.
$$

辛条件 $M^T J M = J$ 展开：

$$
\begin{pmatrix} A^T & C^T \\ B^T & D^T \end{pmatrix}
\begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}
\begin{pmatrix} A & B \\ C & D \end{pmatrix} =
\begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}.
$$

先算后两个矩阵的乘积：
$$
\begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}
\begin{pmatrix} A & B \\ C & D \end{pmatrix} =
\begin{pmatrix} C & D \\ -A & -B \end{pmatrix}.
$$

再乘以 $(M^T)$：
$$
\begin{pmatrix} A^T & C^T \\ B^T & D^T \end{pmatrix}
\begin{pmatrix} C & D \\ -A & -B \end{pmatrix} =
\begin{pmatrix} A^T C - C^T A & A^T D - C^T B \\ B^T C - D^T A & B^T D - D^T B \end{pmatrix}.
$$

要求等于 $\begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}$。因此：

$$
\boxed{
\begin{aligned}
A^T C - C^T A &= 0 \quad (\text{即 } A^T C \text{ 对称}), \\
B^T D - D^T B &= 0 \quad (\text{即 } B^T D \text{ 对称}), \\
A^T D - C^T B &= I_n, \\
B^T C - D^T A &= -I_n \quad (\text{即 } D^T A - B^T C = I_n).
\end{aligned}} \tag{P.5.15}
$$

这四条是辛条件的分块形式。当 $M$ 是正则变换的雅可比矩阵时，这些条件直接对应于泊松括号的保持。

---

### 六、辛矩阵与哈密顿力学的联系

辛矩阵不是抽象的数学概念——它是哈密顿力学几何结构的代数化身。

#### 6.1 哈密顿正则方程的辛形式

将相空间坐标记作 $\boldsymbol{\xi} = (q_1, \ldots, q_n, p_1, \ldots, p_n)^T$。哈密顿函数 $H(\boldsymbol{\xi}, t)$。哈密顿正则方程为

$$
\dot{q}_i = \frac{\partial H}{\partial p_i}, \quad \dot{p}_i = -\frac{\partial H}{\partial q_i}.
$$

用 $J$ 和梯度向量 $\nabla_{\boldsymbol{\xi}} H = (\partial H/\partial q_1, \ldots, \partial H/\partial q_n, \partial H/\partial p_1, \ldots, \partial H/\partial p_n)^T$，可以紧凑地写成

$$
\boxed{\dot{\boldsymbol{\xi}} = J \nabla_{\boldsymbol{\xi}} H}. \tag{P.5.16}
$$

验证：$(J \nabla H)_i = \sum_{j=1}^{2 n} J_{ij} (\nabla H)_j$。对 $i \leq n$（$q$ 部分），$J_{i, n+i}=1$，其余为零，所以 $(J \nabla H)_i = \partial H/\partial p_i = \dot{q}_i$。对 $i > n$（$p$ 部分），$J_{i, i-n} = -1$，所以 $(J \nabla H)_i = -\partial H/\partial q_{i-n} = \dot{p}_{i-n}$。这正是正则方程。

#### 6.2 线性哈密顿系统的相流是辛矩阵

考虑**线性**哈密顿系统：$H = \frac{1}{2} \boldsymbol{\xi}^T S \boldsymbol{\xi}$，其中 $S$ 是 $2 n \times 2 n$ 实对称矩阵。则 $\nabla H = S \boldsymbol{\xi}$，正则方程为

$$
\dot{\boldsymbol{\xi}} = J S \boldsymbol{\xi} \equiv K \boldsymbol{\xi}, \quad K = J S. \tag{P.5.17}
$$

$K$ 满足 $J K = -K^T J$（因为 $(J S)^T = S^T J^T = -S J$，而 $J K = J (J S) = -S$，$K^T J = (J S)^T J = -S J J = S$，等等——关键性质是 $K$ 是**无穷小辛矩阵**：$J K + K^T J = 0$）。

解为 $\boldsymbol{\xi}(t) = e^{K t} \boldsymbol{\xi}(0)$。**演化矩阵 $M(t) = e^{K t}$ 是辛矩阵**：

**证明**：需验证 $M(t)^T J M(t) = J$。对 $t$ 求导：
$$
\frac{d}{dt}(M^T J M) = \dot{M}^T J M + M^T J \dot{M}.
$$
$\dot{M} = K M$，$\dot{M}^T = M^T K^T$。所以
$$
\frac{d}{dt}(M^T J M) = M^T K^T J M + M^T J K M = M^T (K^T J + J K) M.
$$
由 $K = J S$，$S^T = S$：
$$
K^T J + J K = (J S)^T J + J (J S) = S^T J^T J + J^2 S = S (-J^2) - S? 
$$
仔细算：$J^T = -J$，所以 $(J S)^T = S^T J^T = -S J$。$K^T J = -S J J = S$（因为 $J^2 = -I$）。$J K = J (J S) = J^2 S = -S$。所以 $K^T J + J K = S - S = 0$。因此 $M^T J M$ 是常数。$t=0$ 时 $M(0)=I$，$I^T J I = J$。所以对所有 $t$，$M(t)^T J M(t) = J$。

**结论**：线性哈密顿系统的相流是单参数辛矩阵群。这直接对应了刘维尔定理——相空间体积（更准确地说，辛体积形式）在演化下不变。

#### 6.3 正则变换的雅可比矩阵是辛矩阵

更一般地，坐标变换 $\boldsymbol{\xi} \to \boldsymbol{\eta} = \boldsymbol{\eta}(\boldsymbol{\xi})$ 如果是**正则变换**（保持哈密顿方程形式），则其雅可比矩阵

$$
M_{ij} = \frac{\partial \eta_i}{\partial \xi_j}
$$

在每一点都是辛矩阵：

$$
M^T J M = J.
$$

**推导概要**：正则变换保持泊松括号：$\{\eta_i, \eta_j\}_{\boldsymbol{\xi}} = J_{ij}$。而

$$
\{\eta_i, \eta_j\} = \sum_{k,l} \frac{\partial \eta_i}{\partial \xi_k} J_{kl} \frac{\partial \eta_j}{\partial \xi_l} = (M J M^T)_{ij}.
$$

要求等于 $J_{ij}$，即 $M J M^T = J$，或等价地 $M^T J M = J$（利用 $J^{-1} = -J$ 和 $M^T J M = J \iff M J M^T = J$，因为 $M^T J M = J \Rightarrow M^{-T} J M^{-1} = J \Rightarrow$ 等等——实际上两者等价，因为辛群在逆和转置下封闭）。

这是哈密顿力学的核心几何事实：**正则变换就是相空间上保持辛形式的微分同胚**。辛条件 $M^T J M = J$ 是正则变换在局部（线性近似或线性变换）上的体现。

#### 6.4 辛本征值与稳定性

对于 $2 n \times 2 n$ 维哈密顿矩阵 $K = J S$（无穷小辛矩阵），若 $\lambda$ 是本征值，则 $-\lambda$ 也是（因为 $K$ 无迹且特征多项式只含偶次项或某种对称性——实际上，从 $J K + K^T J = 0$ 可推出 $K$ 与 $-K^T$ 相似，从而特征值成对 $\pm \lambda$）。

结合起来：线性哈密顿系统的本征值以四元组 $(\lambda, -\lambda, \bar{\lambda}, -\bar{\lambda})$ 出现。这决定了平衡点的稳定性分类：
- 若存在 $\text{Re}(\lambda) > 0$，平衡点不稳定（双曲型）。
- 若所有本征值纯虚数且对角化，平衡点椭圆型（稳定，小振动）。
- 若出现重本征值而无法对角化（若尔当块），可能出现多项式增长的共振不稳定性（Krein 碰撞）。

这套分类在分析力学的小振动和稳定性分析（如 KAM 定理的前提条件）中至关重要。

---

### 七、总结

P.5 节建立了分析力学必需的矩阵理论基础：

1. **本征值与本征向量**：$A\mathbf{v} = \lambda \mathbf{v}$。本征值由特征方程 $\det(A - \lambda I) = 0$ 确定。实对称矩阵的本征值全为实数，本征向量可正交归一化。

2. **矩阵对角化**：若存在 $n$ 个线性无关本征向量，$A = P \Lambda P^{-1}$。实对称矩阵可正交对角化：$A = Q \Lambda Q^T$。不可对角化时化为若尔当标准型。

3. **二次型**：$Q = \mathbf{x}^T A \mathbf{x}$（$A$ 对称）。通过正交变换可化为平方和 $\sum \lambda_i y_i^2$（主轴定理）。正定性由本征值的正负号决定。小振动理论中同时对角化 $T$ 和 $V$ 等价于求解广义本征值问题 $K\mathbf{v} = \lambda M \mathbf{v}$。

4. **辛矩阵**：$M^T J M = J$，其中 $J = \begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}$。辛矩阵行列式为 $1$，本征值成对 $(\lambda, 1/\lambda)$，构成辛群 $\text{Sp}(2 n,\mathbb{R})$。

5. **与哈密顿力学的联系**：正则方程为 $\dot{\boldsymbol{\xi}} = J \nabla H$。线性哈密顿系统的演化矩阵是辛矩阵。正则变换的雅可比矩阵处处是辛矩阵。辛本征值结构决定了平衡点的稳定性。

矩阵理论不是孤立的数学技巧，而是相空间几何的代数语言。本征值问题对角化动能和势能，辛矩阵保证正则变换保持相空间的辛结构，二者共同构成了分析力学从计算到几何的完整数学基础。至此，预备章 P.1–P.5 的数学工具全部建立，下一步可以正式进入分析力学的基本原理。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [1.1 牛顿力学的回顾与局限](https://zhuanlan.zhihu.com/p/2061652458804343618)

**课程目录：** [分析力学目录](https://zhuanlan.zhihu.com/p/2061652577167717816)
<!-- zhihu-nav-bottom:end -->
