---
created: 2026-07-20T06:03:31+08:00
modified: 2026-07-20T06:25:24+08:00
zhihu-title: M.1 向量代数
zhihu-topics:
  - 力学
zhihu-link: https://zhuanlan.zhihu.com/p/2062422859478196964
zhihu-toc: true
zhihu-created-at: 2026-07-20 10:40
zhihu-updated-at: 2026-07-21 14:24
---
## M.1 向量代数

在进入力学之前，我们必须先掌握描述空间关系的语言。物理世界是三维的，力、速度、加速度这些物理量不仅有大小，还有方向——它们都是向量。标量（如质量、温度、能量）用一个数就可以完整描述；向量（如位移、速度、力）需要同时给出大小和方向。

这一节我们从最基础的几何直觉出发，建立向量的完整代数体系。每一步都给出详细的推导过程，不跳过任何代数步骤。

---

### 一、向量的几何表示与分量表示

**向量的几何定义**：向量是既有大小又有方向的量。在几何上，我们用一条有向线段来表示向量——线段的长度代表大小，箭头的指向代表方向。

向量通常用粗体字母或带箭头的字母表示：$\mathbf{A}$、$\vec{A}$ 或 $\mathbf{v}$、$\vec{v}$。手写时写作 $\vec{A}$。

**向量的相等**：两个向量相等，当且仅当它们的大小相等且方向相同。向量与它在空间中的位置无关——平移一个向量不改变它。这是向量的“自由”性质，在物理中极其重要：力的作用点虽然有意义，但力本身的向量表示只关心大小和方向。

**向量的分量表示**：在三维空间中选定一个直角坐标系 $(x,y,z)$，基向量为 $\hat{i}, \hat{j}, \hat{k}$（或 $\hat{e}_x, \hat{e}_y, \hat{e}_z$）。任意向量 $\mathbf{A}$ 可以唯一地分解为三个基向量的线性组合：

$$
\mathbf{A} = A_x \hat{i} + A_y \hat{j} + A_z \hat{k}. \tag{1}
$$

其中 $A_x, A_y, A_z$ 是 $\mathbf{A}$ 在三个坐标轴上的投影——它们是实数，称为向量的分量。

向量的大小（模长）由勾股定理给出：

$$
|\mathbf{A}| = \sqrt{A_x^2 + A_y^2 + A_z^2}. \tag{2}
$$

分量表示与几何表示完全等价——给定分量可以画出向量，给定向量可以读出分量。分量表示的优势在于，它将几何运算转化为数的运算。

**方向余弦**：向量 $\mathbf{A}$ 与三个坐标轴的夹角 $\alpha, \beta, \gamma$ 满足：

$$
\cos\alpha = \frac{A_x}{|\mathbf{A}|}, \quad \cos\beta = \frac{A_y}{|\mathbf{A}|}, \quad \cos\gamma = \frac{A_z}{|\mathbf{A}|}. \tag{3}
$$

且 $\cos^2\alpha + \cos^2\beta + \cos^2\gamma = 1$。这三个余弦完全确定了向量的方向。

**单位向量**：一个向量除以它自己的模长，得到该方向的单位向量：

$$
\hat{\mathbf{A}} = \frac{\mathbf{A}}{|\mathbf{A}|}, \quad |\hat{\mathbf{A}}| = 1. \tag{4}
$$

单位向量只表示方向，不包含大小信息。

---

### 二、向量加减法

**向量加法**：两个向量 $\mathbf{A}$ 和 $\mathbf{B}$ 相加，得到一个新的向量 $\mathbf{C} = \mathbf{A} + \mathbf{B}$。

几何上有两种等价方法：

**平行四边形法则**：将 $\mathbf{A}$ 和 $\mathbf{B}$ 的起点放在同一点，以它们为邻边作平行四边形，从公共起点出发的对角线就是 $\mathbf{A} + \mathbf{B}$。

**三角形法则**：将 $\mathbf{B}$ 的起点放在 $\mathbf{A}$ 的终点，从 $\mathbf{A}$ 的起点到 $\mathbf{B}$ 的终点的连线就是 $\mathbf{A} + \mathbf{B}$。

这两种方法等价——平行四边形的一条对角线正是三角形法则中连接起点和终点的线段。

**分量形式**：设 $\mathbf{A} = A_x\hat{i} + A_y\hat{j} + A_z\hat{k}$，$\mathbf{B} = B_x\hat{i} + B_y\hat{j} + B_z\hat{k}$。则

$$
\mathbf{A} + \mathbf{B} = (A_x + B_x)\hat{i} + (A_y + B_y)\hat{j} + (A_z + B_z)\hat{k}. \tag{5}
$$

推导：直接合并同类项即可，因为基向量相同。

向量加法满足：
- 交换律：$\mathbf{A} + \mathbf{B} = \mathbf{B} + \mathbf{A}$
- 结合律：$(\mathbf{A} + \mathbf{B}) + \mathbf{C} = \mathbf{A} + (\mathbf{B} + \mathbf{C})$
- 零向量 $\mathbf{0}$ 是加法单位元：$\mathbf{A} + \mathbf{0} = \mathbf{A}$
- 每个向量有逆向量：$\mathbf{A} + (-\mathbf{A}) = \mathbf{0}$

**向量减法**：$\mathbf{A} - \mathbf{B} = \mathbf{A} + (-\mathbf{B})$，其中 $-\mathbf{B}$ 是与 $\mathbf{B}$ 大小相同、方向相反的向量。在几何上，$\mathbf{A} - \mathbf{B}$ 是从 $\mathbf{B}$ 的终点指向 $\mathbf{A}$ 的终点的向量（当两者起点重合时）。

分量形式：

$$
\mathbf{A} - \mathbf{B} = (A_x - B_x)\hat{i} + (A_y - B_y)\hat{j} + (A_z - B_z)\hat{k}. \tag{6}
$$

---

### 三、数乘向量

一个实数 $\lambda$（标量）与向量 $\mathbf{A}$ 相乘，得到新向量 $\lambda\mathbf{A}$：

- 大小：$|\lambda\mathbf{A}| = |\lambda| \cdot |\mathbf{A}|$
- 方向：若 $\lambda > 0$，方向与 $\mathbf{A}$ 相同；若 $\lambda < 0$，方向与 $\mathbf{A}$ 相反；若 $\lambda = 0$，得到零向量。

数乘满足：
- 分配律（对标量加法）：$(\lambda + \mu)\mathbf{A} = \lambda\mathbf{A} + \mu\mathbf{A}$
- 分配律（对向量加法）：$\lambda(\mathbf{A} + \mathbf{B}) = \lambda\mathbf{A} + \lambda\mathbf{B}$
- 结合律：$\lambda(\mu\mathbf{A}) = (\lambda\mu)\mathbf{A}$

分量形式：

$$
\lambda\mathbf{A} = (\lambda A_x)\hat{i} + (\lambda A_y)\hat{j} + (\lambda A_z)\hat{k}. \tag{7}
$$

---

### 四、点乘（标量积）

两个向量的点乘（也称内积、标量积）定义为：

$$
\mathbf{A}\cdot\mathbf{B} = |\mathbf{A}|\,|\mathbf{B}|\cos\theta, \tag{8}
$$

其中 $\theta$ 是 $\mathbf{A}$ 与 $\mathbf{B}$ 之间的夹角（$0 \le \theta \le \pi$）。

点乘的结果是一个标量（实数），不是向量。

**点乘的几何意义**：$\mathbf{A}\cdot\mathbf{B}$ 等于 $\mathbf{A}$ 的模长乘以 $\mathbf{B}$ 在 $\mathbf{A}$ 方向上的投影（即 $B\cos\theta$），也等于 $\mathbf{B}$ 的模长乘以 $\mathbf{A}$ 在 $\mathbf{B}$ 方向上的投影（即 $A\cos\theta$）。这就是点乘在物理中频繁出现的原因——它天然地表达了“一个量在另一个量方向上的分量”。

点乘满足：
- 交换律：$\mathbf{A}\cdot\mathbf{B} = \mathbf{B}\cdot\mathbf{A}$（因为 $\cos\theta$ 对称）
- 分配律：$\mathbf{A}\cdot(\mathbf{B} + \mathbf{C}) = \mathbf{A}\cdot\mathbf{B} + \mathbf{A}\cdot\mathbf{C}$
- 与数乘结合：$(\lambda\mathbf{A})\cdot\mathbf{B} = \lambda(\mathbf{A}\cdot\mathbf{B})$

**基向量的点乘**：$\hat{i}, \hat{j}, \hat{k}$ 两两垂直，且各自长度为 $1$：

$$
\hat{i}\cdot\hat{i} = 1, \quad \hat{j}\cdot\hat{j} = 1, \quad \hat{k}\cdot\hat{k} = 1,
$$
$$
\hat{i}\cdot\hat{j} = 0, \quad \hat{j}\cdot\hat{k} = 0, \quad \hat{k}\cdot\hat{i} = 0. \tag{9}
$$

用克罗内克符号 $\delta_{ij}$ 可统一写为 $\hat{e}_i \cdot \hat{e}_j = \delta_{ij}$，其中 $\delta_{ij} = 1$ 当 $i=j$，否则为 $0$。

**分量形式的推导**：设 $\mathbf{A} = A_x\hat{i} + A_y\hat{j} + A_z\hat{k}$，$\mathbf{B} = B_x\hat{i} + B_y\hat{j} + B_z\hat{k}$。利用分配律展开：

$$
\begin{aligned}
\mathbf{A}\cdot\mathbf{B} =& (A_x\hat{i} + A_y\hat{j} + A_z\hat{k}) \cdot (B_x\hat{i} + B_y\hat{j} + B_z\hat{k}) \\
=& A_xB_x(\hat{i}\cdot\hat{i}) + A_xB_y(\hat{i}\cdot\hat{j}) + A_xB_z(\hat{i}\cdot\hat{k}) \\
&+ A_yB_x(\hat{j}\cdot\hat{i}) + A_yB_y(\hat{j}\cdot\hat{j}) + A_yB_z(\hat{j}\cdot\hat{k}) \\
&+ A_zB_x(\hat{k}\cdot\hat{i}) + A_zB_y(\hat{k}\cdot\hat{j}) + A_zB_z(\hat{k}\cdot\hat{k}).
\end{aligned}
$$

代入(9)中的基向量点乘结果，非对角项均为 $0$，对角项均为 $1$，得到：

$$
\mathbf{A}\cdot\mathbf{B} = A_xB_x + A_yB_y + A_zB_z. \tag{10}
$$

这是点乘最常用的计算形式。

**点乘的物理应用**：
- 功：$W = \mathbf{F}\cdot\Delta\mathbf{r}$（力在位移方向上的分量乘以位移）
- 功率：$P = \mathbf{F}\cdot\mathbf{v}$
- 磁通量：$\Phi = \mathbf{B}\cdot\mathbf{S}$
- 判断垂直：若 $\mathbf{A}\cdot\mathbf{B} = 0$，则 $\mathbf{A} \perp \mathbf{B}$（或其中一个是零向量）

---

### 五、叉乘（矢量积）

两个向量的叉乘（也称外积、矢量积）定义为：

$$
\mathbf{A}\times\mathbf{B} = (|\mathbf{A}|\,|\mathbf{B}|\sin\theta)\,\hat{\mathbf{n}}, \tag{11}
$$

其中 $\theta$ 是 $\mathbf{A}$ 与 $\mathbf{B}$ 之间的夹角（$0 \le \theta \le \pi$），$\hat{\mathbf{n}}$ 是垂直于 $\mathbf{A}$ 和 $\mathbf{B}$ 所在平面的单位向量，方向由**右手定则**确定：右手四指从 $\mathbf{A}$ 弯向 $\mathbf{B}$（经过较小的角度），大拇指指向 $\hat{\mathbf{n}}$ 的方向。

叉乘的结果是一个向量，其大小等于 $\mathbf{A}$ 和 $\mathbf{B}$ 张成的平行四边形的面积（$|\mathbf{A}||\mathbf{B}|\sin\theta$），方向垂直于该平行四边形所在的平面。

**叉乘的核心性质**：

- 反交换律：$\mathbf{A}\times\mathbf{B} = -\mathbf{B}\times\mathbf{A}$
  证明：由定义，$\mathbf{B}\times\mathbf{A}$ 的大小相同，但方向由从 $\mathbf{B}$ 到 $\mathbf{A}$ 的右手定则决定，这与 $\mathbf{A}\times\mathbf{B}$ 的方向相反，因此差一个负号。
  
- 分配律：$\mathbf{A}\times(\mathbf{B} + \mathbf{C}) = \mathbf{A}\times\mathbf{B} + \mathbf{A}\times\mathbf{C}$
- 与数乘结合：$(\lambda\mathbf{A})\times\mathbf{B} = \lambda(\mathbf{A}\times\mathbf{B})$

- **不满足结合律**：$(\mathbf{A}\times\mathbf{B})\times\mathbf{C} \neq \mathbf{A}\times(\mathbf{B}\times\mathbf{C})$。具体关系将在后面推导。

**基向量的叉乘**：

根据右手定则和正交性：

$$
\hat{i}\times\hat{i} = 0, \quad \hat{j}\times\hat{j} = 0, \quad \hat{k}\times\hat{k} = 0,
$$
$$
\hat{i}\times\hat{j} = \hat{k}, \quad \hat{j}\times\hat{k} = \hat{i}, \quad \hat{k}\times\hat{i} = \hat{j},
$$
$$
\hat{j}\times\hat{i} = -\hat{k}, \quad \hat{k}\times\hat{j} = -\hat{i}, \quad \hat{i}\times\hat{k} = -\hat{j}. \tag{12}
$$

这些关系可以统一用列维-奇维塔符号 $\epsilon_{ijk}$ 写为 $\hat{e}_i \times \hat{e}_j = \sum_{k} \epsilon_{ijk} \hat{e}_k$，其中 $\epsilon_{ijk}$ 是全反对称张量，$\epsilon_{123}=1$，交换任意两个指标改变符号，重复指标则为 $0$。

**分量形式的推导**：设 $\mathbf{A} = A_x\hat{i} + A_y\hat{j} + A_z\hat{k}$，$\mathbf{B} = B_x\hat{i} + B_y\hat{j} + B_z\hat{k}$。利用分配律展开：

$$
\begin{aligned}
\mathbf{A}\times\mathbf{B} =& (A_x\hat{i} + A_y\hat{j} + A_z\hat{k}) \times (B_x\hat{i} + B_y\hat{j} + B_z\hat{k}) \\
=& A_xB_x(\hat{i}\times\hat{i}) + A_xB_y(\hat{i}\times\hat{j}) + A_xB_z(\hat{i}\times\hat{k}) \\
&+ A_yB_x(\hat{j}\times\hat{i}) + A_yB_y(\hat{j}\times\hat{j}) + A_yB_z(\hat{j}\times\hat{k}) \\
&+ A_zB_x(\hat{k}\times\hat{i}) + A_zB_y(\hat{k}\times\hat{j}) + A_zB_z(\hat{k}\times\hat{k}).
\end{aligned}
$$

代入(12)中的基向量叉乘结果：

$$
\begin{aligned}
\mathbf{A}\times\mathbf{B} =& A_xB_y\hat{k} + A_xB_z(-\hat{j}) + A_yB_x(-\hat{k}) + A_yB_z\hat{i} + A_zB_x\hat{j} + A_zB_y(-\hat{i}) \\
=& (A_yB_z - A_zB_y)\hat{i} + (A_zB_x - A_xB_z)\hat{j} + (A_xB_y - A_yB_x)\hat{k}.
\end{aligned}
$$

因此：

$$
\mathbf{A}\times\mathbf{B} =
\begin{vmatrix}
\hat{i} & \hat{j} & \hat{k} \\
A_x & A_y & A_z \\
B_x & B_y & B_z
\end{vmatrix}. \tag{13}
$$

展开这个行列式（按第一行展开）即可得到上式：

$$
\mathbf{A}\times\mathbf{B} = (A_yB_z - A_zB_y)\hat{i} - (A_xB_z - A_zB_x)\hat{j} + (A_xB_y - A_yB_x)\hat{k}.
$$

注意第二项的负号来自行列式展开的符号规则（$-$ 对应 $\hat{j}$ 项）。这等价于上面写的形式。

**行列式记法是最容易记忆的**——第一行是基向量，第二行是 $\mathbf{A}$ 的分量，第三行是 $\mathbf{B}$ 的分量。

**叉乘的物理应用**：
- 力矩：$\boldsymbol{\tau} = \mathbf{r}\times\mathbf{F}$（力臂乘以力，方向沿转轴）
- 角动量：$\mathbf{L} = \mathbf{r}\times\mathbf{p}$
- 洛伦兹力：$\mathbf{F} = q\mathbf{v}\times\mathbf{B}$
- 安培力：$d\mathbf{F} = I\,d\mathbf{l}\times\mathbf{B}$

---

### 六、混合积（标量三重积）

三个向量 $\mathbf{A}$、$\mathbf{B}$、$\mathbf{C}$ 的混合积定义为：

$$
\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C}). \tag{14}
$$

先做叉乘得到一个向量，再与第三个向量做点乘，结果是标量。

**几何意义**：混合积的绝对值等于以 $\mathbf{A}$、$\mathbf{B}$、$\mathbf{C}$ 为三条邻边的平行六面体的体积。

为什么？$\mathbf{B}\times\mathbf{C}$ 的大小等于 $\mathbf{B}$ 和 $\mathbf{C}$ 张成的平行四边形的面积，方向垂直于该平面。$\mathbf{A}$ 与这个向量的点乘等于 $\mathbf{A}$ 在该垂直方向上的投影（即高）乘以底面积——这正是平行六面体的体积。

**分量形式**：

由(10)和(13)：

$$
\begin{aligned}
\mathbf{B}\times\mathbf{C} &= (B_yC_z - B_zC_y)\hat{i} + (B_zC_x - B_xC_z)\hat{j} + (B_xC_y - B_yC_x)\hat{k}, \\
\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C}) &= A_x(B_yC_z - B_zC_y) + A_y(B_zC_x - B_xC_z) + A_z(B_xC_y - B_yC_x).
\end{aligned}
$$

这正是行列式：

$$
\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C}) =
\begin{vmatrix}
A_x & A_y & A_z \\
B_x & B_y & B_z \\
C_x & C_y & C_z
\end{vmatrix}. \tag{15}
$$

**混合积的循环性质**：

$$
\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C}) = \mathbf{B}\cdot(\mathbf{C}\times\mathbf{A}) = \mathbf{C}\cdot(\mathbf{A}\times\mathbf{B}). \tag{16}
$$

证明：在行列式(15)中，将行做循环置换 $(A\to B\to C\to A)$，行列式不变（因为三次置换相当于偶置换）。交换任意两个向量（即交换两行）则行列式变号。

**几何应用**：若 $\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C}) = 0$，则三向量共面（平行六面体体积为零），即它们线性相关。

---

### 七、矢量恒等式：三重叉积（BAC-CAB）

三个向量的三重叉积 $\mathbf{A}\times(\mathbf{B}\times\mathbf{C})$ 是一个向量。它满足以下恒等式（称为“BAC-CAB 恒等式”）：

$$
\mathbf{A}\times(\mathbf{B}\times\mathbf{C}) = \mathbf{B}(\mathbf{A}\cdot\mathbf{C}) - \mathbf{C}(\mathbf{A}\cdot\mathbf{B}). \tag{17}
$$

记忆口诀：**BAC 减 CAB**——括号外的向量 $\mathbf{A}$ 依次与括号内的两个向量点乘，分别乘以括号内另一个向量。

下面给出完整推导，不省略任何代数步骤。

设 $\mathbf{A} = (A_x, A_y, A_z)$，$\mathbf{B} = (B_x, B_y, B_z)$，$\mathbf{C} = (C_x, C_y, C_z)$。

先计算 $\mathbf{D} = \mathbf{B}\times\mathbf{C}$，其分量为：

$$
D_x = B_yC_z - B_zC_y, \quad D_y = B_zC_x - B_xC_z, \quad D_z = B_xC_y - B_yC_x. \tag{18}
$$

然后计算 $\mathbf{E} = \mathbf{A}\times\mathbf{D}$。利用叉乘分量公式(13)：

$$
E_x = A_yD_z - A_zD_y, \quad E_y = A_zD_x - A_xD_z, \quad E_z = A_xD_y - A_yD_x. \tag{19}
$$

我们只详细推导 $x$ 分量，$y$ 和 $z$ 分量完全对称。

将(18)代入 $E_x$：

$$
\begin{aligned}
E_x &= A_y(B_xC_y - B_yC_x) - A_z(B_zC_x - B_xC_z) \\
&= A_yB_xC_y - A_yB_yC_x - A_zB_zC_x + A_zB_xC_z.
\end{aligned}
$$

重新组合，将含有 $B_x$ 的项放在一起，含有 $C_x$ 的项放在一起：

$$
E_x = B_x(A_yC_y + A_zC_z) - C_x(A_yB_y + A_zB_z). \tag{20}
$$

现在，我们要验证它等于 $B_x(\mathbf{A}\cdot\mathbf{C}) - C_x(\mathbf{A}\cdot\mathbf{B})$。

计算右边的 $x$ 分量：

$$
B_x(\mathbf{A}\cdot\mathbf{C}) - C_x(\mathbf{A}\cdot\mathbf{B}) = B_x(A_xC_x + A_yC_y + A_zC_z) - C_x(A_xB_x + A_yB_y + A_zB_z).
$$

展开：

$$
= B_xA_xC_x + B_x(A_yC_y + A_zC_z) - C_xA_xB_x - C_x(A_yB_y + A_zB_z).
$$

第一项 $B_xA_xC_x$ 和第三项 $C_xA_xB_x$ 相等（因为乘法交换），相互抵消。于是得到：

$$
= B_x(A_yC_y + A_zC_z) - C_x(A_yB_y + A_zB_z).
$$

这与(20)完全一致。因此 $E_x = B_x(\mathbf{A}\cdot\mathbf{C}) - C_x(\mathbf{A}\cdot\mathbf{B})$。

同理可证 $y$ 和 $z$ 分量也成立。因此(17)得证。

**注意**：括号的位置不能挪动。$(\mathbf{A}\times\mathbf{B})\times\mathbf{C}$ 给出不同的结果：

$$
(\mathbf{A}\times\mathbf{B})\times\mathbf{C} = \mathbf{B}(\mathbf{A}\cdot\mathbf{C}) - \mathbf{A}(\mathbf{B}\cdot\mathbf{C}). \tag{21}
$$

这个公式可以这样记忆：最外侧的向量 $\mathbf{C}$ 与最近的向量 $\mathbf{B}$ 点乘，乘以内侧的另一个向量 $\mathbf{A}$，再减去另一项。

(17)和(21)的区别在于哪个向量在点乘中与括号外的向量配对——注意括号位置决定配对规则。

---

### 八、总结

把 M.1 节的全部内容压缩成几条线：

1. **向量**是既有大小又有方向的量。分量表示 $\mathbf{A} = A_x\hat{i} + A_y\hat{j} + A_z\hat{k}$ 将几何运算转化为数的运算，模长 $|\mathbf{A}| = \sqrt{A_x^2 + A_y^2 + A_z^2}$。

2. **向量加法**：平行四边形法则/三角形法则，分量对应相加：$\mathbf{A} + \mathbf{B} = (A_x+B_x)\hat{i} + (A_y+B_y)\hat{j} + (A_z+B_z)\hat{k}$。加法满足交换律和结合律。

3. **数乘**：改变向量大小（$\lambda<0$ 时反向），分量对应相乘：$\lambda\mathbf{A} = (\lambda A_x)\hat{i} + (\lambda A_y)\hat{j} + (\lambda A_z)\hat{k}$。

4. **点乘** $\mathbf{A}\cdot\mathbf{B} = |\mathbf{A}||\mathbf{B}|\cos\theta$，结果是标量。分量形式 $A_xB_x + A_yB_y + A_zB_z$。交换律、分配律成立。物理中表达“投影”和“功”。

5. **叉乘** $\mathbf{A}\times\mathbf{B} = |\mathbf{A}||\mathbf{B}|\sin\theta\,\hat{\mathbf{n}}$，结果是向量，方向由右手定则确定。分量形式为行列式。反交换律成立，不满足结合律。物理中表达“力矩”和“角动量”。

6. **混合积** $\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C})$ 是标量，绝对值等于平行六面体体积。分量形式为 $3\times3$ 行列式。循环置换不变，交换两向量变号。判断共面的判据。

7. **三重叉积**：$\mathbf{A}\times(\mathbf{B}\times\mathbf{C}) = \mathbf{B}(\mathbf{A}\cdot\mathbf{C}) - \mathbf{C}(\mathbf{A}\cdot\mathbf{B})$（BAC-CAB 恒等式），证明已完整展示。括号位置不同导致公式不同。

这些向量代数的规则是力学乃至整个物理学的语言基础。点乘、叉乘、混合积、三重叉积在后续章节中会反复出现——从功、力矩、角动量到刚体运动、电磁学，无一例外。掌握这些运算，就是掌握了物理学的语法。

---

<!-- zhihu-nav-bottom:start -->
**下一节：** [M.2 微积分在物理中的基本应用](https://zhuanlan.zhihu.com/p/2062423053758363601)

**课程目录：** [力学目录](https://zhuanlan.zhihu.com/p/2062487053170906945)
<!-- zhihu-nav-bottom:end -->
