# Singular Value Decomposition (SVD)




奇异值分解(SVD)是一种矩阵因子分解方法，在线性代数中，被广泛应用。  
奇异值分解也是一种矩阵近似的方法，这个近似是在弗罗贝尼乌斯范数（Frobenius norm) 意义下的近似。  
奇异值分解是在平方损失(弗罗贝尼乌斯范数)意义下对矩阵的最优近似，即数据压缩。

## 定义

### 1. 奇异值分解

将一个 $m \times n$ 的实矩阵$A$，$A \in \mathbf{R}^{m \times n}$ 表示为以下三个实矩阵乘积形式的运算，即矩阵因子分解：

$$
A=U \Sigma V^{\mathrm{T}}
$$

其中， 
- $U$是$m$阶正交矩阵: $UU^T = I$
- $V$为$n$阶正交矩阵: $VV^T = I$
- $\Sigma$是由降序排列的非负的对角线元素组成的$m \times n$的对角矩阵: 
  - $\Sigma = diag(\sigma_1, \sigma_2, \cdots, \sigma_p)$
  - $\sigma_{1} \geqslant \sigma_{2} \geqslant \cdots \geqslant \sigma_{p} \geqslant 0$
  - $p = \min(m, n)$

那么，称

- $U \Sigma V^{\mathrm{T}}$ : 矩阵A的奇异值分解
- $\sigma_i$ 为A的奇异值(singluar value)
- $U$ 的列向量为左奇异向量(left singular vector)
- $V$ 的列向量为右奇异向量(right singular vector)

⚠️注意，矩阵的奇异值分解不唯一。

![SVD](/images/stats/SVD.png)


**实际常用的是：**

### 2. 紧奇异值分解(compact singular value decomposition)： 无损压缩  

设$m \times n$的实矩阵$A$，其秩为$\operatorname{rank}(A)=r$, $r \leqslant \min (m, n)$, 则紧奇异值分解为：

$$
A=U_r \Sigma_r V^{\mathrm{T}}_r
$$

其中，

- $U_r$ 是 $m \times r$ 矩阵
- $V_r$ 是 $n \times r$ 矩阵
- $\Sigma_r$ 是 $r$ 阶对角阵
- $r = \operatorname{rank}(A)$


### 3. 截断奇异值分解(truncated singular value decomposition)： 有损压缩

一般讲奇异值分解，实际上多指截断奇异值分解

在奇异值分解中， 只取最大的$k$个奇异值($k < r, r= \operatorname{rank}(A)$ )对应的部分，就得到截断奇异值分解

设 $m \times n$ 的实矩阵 $A$，其秩为 $\operatorname{rank}(A)=r$, 且 $0 < k < r$, 则截断奇异值分解为：

$$
A \approx U_{k} \Sigma_{k} V_{k}^{\mathrm{T}}
$$

其中，

- $U_k$ 是$m \times k$ 矩阵(前 $k$列)
- $V_k$ 是$n \times k$ 矩阵(前 $k$列)
- $\Sigma_k$ 是 $k$ 阶对角阵(前 $k$个)
- $r = \operatorname{rank}(A)$


## 几何解释

从线性变换的角度理解奇异值分解：
将 $m \times n$ 的实矩阵 $A$表示为从 $n$ 维空间 $R_n$ 到 $m$ 维空间 $R_m$的一个线性变换：

$$
T: x \rightarrow A x
$$

$x$, $Ax$为各自空间的向量。

那么**线性变换**可以理解为：

1. 一个坐标系的旋转或反射变换
2. 一个坐标轴的缩放变换
3. 另一个坐标系的旋转或反射变换

对A进行奇异值分解，U和V都是正交矩阵

- V的列向量构成Rn空间的一组标准正交基，表示Rn空间中正交坐标系的旋转或反射变换
- U的列向量都成Rm空间的一组标准正交基，表示Rm空间中正交坐标系的旋转或反射变换
- $\Sigma$的对角元素是一组非负实数，表示Rn中的原始正交坐标系坐标轴的缩放变换

![SVD](/images/stats/1200px-Singular-Value-Decomposition.svg.png)

## 奇异值计算

矩阵$A$的奇异值分解可以通过求对阵矩阵$A^TA$的特征值和特征向量得到。

- $A^TA$的特征向量构成正交矩阵$V$的列
- $A^TA$的特征值$\lambda_j$的平方根为奇异值$\sigma_i$，即 

$$
\sigma_{j}=\sqrt{\lambda_{j}}, \quad j=1,2, \cdots, n
$$

- 对$\sigma_i$从大到小排列，得到对角矩阵 $\Sigma$
- 求正奇异值对应的左奇异向量，再扩充的 $A^T$ 的标准正交基，构成正交矩阵 $U$ 的列


### 求值过程

1. 求对阵矩阵$A^TA$的特征值和特征向量

计算对称矩阵 $W= A^TA$  
求解特征方程 $(W - \lambda I)x = 0$
得到特征值$\lambda_i$，并将之降序排列  

$$
\lambda_{1} \geqslant \lambda_{2} \geqslant \cdots \geqslant \lambda_{n} \geqslant 0
$$

将特征值 $\lambda_i$ 代入特征方程求的对应特征向量

2. 求 $n$ 阶正交矩阵 $V$

将特征向量单位化， 得到单位特征向量构成 $n$ 阶正交矩阵V

$$
V=\left[\begin{array}{llll}
v_{1} & v_{2} & \cdots & v_{n}
\end{array}\right]
$$

3. 求 $m \times n$ 对角矩阵 $\Sigma$

计算A的奇异值

$$
\sigma_{j}=\sqrt{\lambda_{j}}, \quad j=1,2, \cdots, n
$$

构造 $m \times n$ 矩阵对角矩阵 $\Sigma$， 主对角线元素是奇异值， 其余元素为 0

$$
\Sigma=\operatorname{diag}\left(\sigma_{1}, \sigma_{2}, \cdots, \sigma_{n}\right)
$$

4. 求 $m$ 阶正交矩阵 $U$

对 $A$ 的前 $r$个正奇异值， 令

$$
u_{j}=\frac{1}{\sigma_{j}} A v_{j}, \quad j=1,2, \cdots, r
$$

得到

$$
U_{1}=\left[\begin{array}{llll}
u_{1} & u_{2} & \cdots & u_{r}
\end{array}\right]
$$

求 $A^T$ 的零空间的一组标准正交基

$$
\lbrace u_{r+1}, u_{r+2}, \cdots, u_{m} \rbrace
$$

令

$$
U_{2}=\left[\begin{array}{llll}
u_{r+1} & u_{r+2} & \cdots & u_{m}
\end{array}\right]
$$

且令

$$
U=\left[\begin{array}{ll}
U_{1} & U_{2}
\end{array}\right]
$$

5. 得到奇异值分解

$$
A=U \Sigma V^{\mathrm{T}}
$$


### 举例

求: 矩阵 $A$ 的奇异值分解

$$
A=\left[\begin{array}{ll}
1 & 1 \cr
2 & 2 \cr
0 & 0
\end{array}\right]
$$

解：  
1. 求 $A^TA$ 的特征值和特征向量  

$$
A^{\mathrm{T}} A=\left[\begin{array}{lll}
1 & 2 & 0 \cr
1 & 2 & 0
\end{array}\right]\left[\begin{array}{ll}
1 & 1 \cr
2 & 2 \cr
0 & 0
\end{array}\right]=\left[\begin{array}{ll}
5 & 5 \cr
5 & 5
\end{array}\right]
$$  

满足特征方程

$$
\left(A^{\mathrm{T}} A-\lambda I\right) x=0
$$

得到齐次线性方程组

$$
\begin{cases}
(5-\lambda) x_{1} + & 5 x_{2}=0 \cr
5 x_{1} + & (5-\lambda) x_{2}=0
\end{cases}
$$

该方程组有非零解的充要条件是

$$
\left|\begin{array}{cc}
5-\lambda & 5 \cr
5 & 5-\lambda
\end{array}\right|=0
$$

即

$$
\lambda^{2}-10 \lambda=0
$$

得到 $\lambda_1 = 10 $, $\lambda_2 = 0$. 

特征值代入线性方程组，分别得到对应单位特征向量

$$
v_{1}=\left[\begin{array}{c}
\frac{1}{\sqrt{2}} \cr
\frac{1}{\sqrt{2}}
\end{array}\right], v_{2}=\left[\begin{array}{c}
\frac{1}{\sqrt{2}} \cr
-\frac{1}{\sqrt{2}}
\end{array}\right] 
$$


2. 求正交矩阵 $V$

构造正交矩阵

$$
V=\left[\begin{array}{cc}
\frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \cr
\frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}}
\end{array}\right]
$$

3. 求对角矩阵 $\Sigma$

奇异值为 $\sigma_1 = \sqrt{\lambda_1} = \sqrt{10}$, $\sigma_2 = 0$, 那么

$$
\Sigma=\left[\begin{array}{cc}
\sqrt{10} & 0 \cr
0 & 0 \cr
0 & 0
\end{array}\right]
$$

⚠️注意： 为了 $\Sigma$ 能与$U$和$V$进行矩阵乘法， $\Sigma$要加上零行向量

4. 求正交矩阵 $U$

基于 $A$ 的奇异值计算得到列向量 $u_1$

$$
u_{1}=\frac{1}{\sigma_{1}} A v_{1}=\frac{1}{\sqrt{10}}\left[\begin{array}{cc}
1 & 1 \cr
2 & 2 \cr
0 & 0
\end{array}\right]\left[\begin{array}{c}
\frac{1}{\sqrt{2}} \cr
\frac{1}{\sqrt{2}}
\end{array}\right]=\left[\begin{array}{c}
\frac{1}{\sqrt{5}} \cr
\frac{2}{\sqrt{5}} \cr
0
\end{array}\right]
$$

列向量$u_2$， $u_3$是 $A^T$的零空间 $N(A^T)$ 的一组标准正交基，故而求解以下线性方程组

$$
A^{\mathrm{T}} x=\left[\begin{array}{lll}
1 & 2 & 0 \cr
1 & 2 & 0
\end{array}\right]\left[\begin{array}{l}
x_{1} \cr
x_{2} \cr
x_{3}
\end{array}\right]=\left[\begin{array}{l}
0 \cr
0
\end{array}\right]
$$

即

$$
\begin{array}{c}
x_{1}+2 x_{2}+0 x_{3}=0 \cr
x_{1}=-2 x_{2}+0 x_{3}
\end{array}
$$

分别取 $x_2$，$x_3$ 为 $(1,0)$ 和 $(0,1)$ 得到 $N(A^T)$的基

$$
(-2,1,0)^{\mathrm{T}}, \quad(0,0,1)^{\mathrm{T}}
$$

得 $N(A^T)$ 的一组标准正交基是

$$
u_{2}=\left(-\frac{2}{\sqrt{5}}, \frac{1}{\sqrt{5}}, 0\right)^{\mathrm{T}}, \quad u_{3}=(0,0,1)^{\mathrm{T}}
$$


最后，构造正交矩阵 $U$

$$
U=\left[\begin{array}{ccc}
\frac{1}{\sqrt{5}} & -\frac{2}{\sqrt{5}} & 0 \cr
\frac{2}{\sqrt{5}} & \frac{1}{\sqrt{5}} & 0 \cr
0 & 0 & 1
\end{array}\right]
$$


5. 矩阵 $A$ 的奇异值分解


$$
A=U \Sigma V^{\mathrm{T}}=\left[\begin{array}{ccc}
\frac{1}{\sqrt{5}} & -\frac{2}{\sqrt{5}} & 0 \cr
\frac{2}{\sqrt{5}} & \frac{1}{\sqrt{5}} & 0 \cr
0 & 0 & 1
\end{array}\right]\left[\begin{array}{cc}
\sqrt{10} & 0 \cr
0 & 0 \cr
0 & 0
\end{array}\right]\left[\begin{array}{cc}
\frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \cr
\frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}}
\end{array}\right]
$$





## 矩阵的外积展开式表示

将A的奇异值分解看成矩阵 $U\Sigma$ 和 $V^T$ 的乘积，  
将 $U\Sigma$ 按列向量分块， 

$$
U \Sigma=\left[\begin{array}{llll}
\sigma_{1} u_{1} & \sigma_{2} u_{2} & \cdots & \sigma_{n} u_{n}
\end{array}\right]
$$


$V^T$ 按行向量分块

$$
V^{\mathrm{T}}=\left[\begin{array}{c}
v_{1}^{\mathrm{T}} \cr
v_{2}^{\mathrm{T}} \cr
\vdots \cr
v_{n}^{\mathrm{T}}
\end{array}\right]
$$

那么外积展开式为

$$
A=\sigma_{1} u_{1} v_{1}^{\mathrm{T}}+\sigma_{2} u_{2} v_{2}^{\mathrm{T}}+\cdots+\sigma_{n} u_{n} v_{n}^{\mathrm{T}}
$$

或者

$$
A=\sum_{k=1}^{n} A_{k}=\sum_{k=1}^{n} \sigma_{k} u_{k} v_{k}^{\mathrm{T}}
$$

其中 $A_{k}=\sigma_{k} u_{k} v_{k}^{\mathrm{T}}$ 是 $m \times n$ 矩阵

而


$$
u_{i} v_{j}^{\mathrm{T}}=\left[\begin{array}{c}
u_{1 i} \cr
u_{2 i} \cr
\vdots \cr
u_{m i}
\end{array}\right]\left[\begin{array}{cccc}
v_{1 j} & v_{2 j} & \cdots & v_{n j}
\end{array}\right] 
= \left[\begin{array}{cccc}
u_{1 i} v_{1 j} & u_{1 i} v_{2 j} & \cdots & u_{1 i} v_{n j} \cr
u_{2 i} v_{1 j} & u_{2 i} v_{2 j} & \cdots & u_{2 i} v_{n j} \cr
\vdots & \vdots & & \vdots \cr
u_{m i} v_{1 j} & u_{m i} v_{2 j} & \cdots & u_{m i} v_{n j}
\end{array}\right]
$$

## 补充：

### Determinant 行列式

Determinant:  
The determinant of a `square matrix` is a `scalar` that provides information about the matrix. e.g. `Invertibility`  

Geometrically, it can be viewed as the `volume scaling factor` of the `linear transformation` described by the matrix. 

- 一个矩阵的行列式就是一个超平行多面体的（有向的）面积/体积，这个多面体的每条边对应着对应矩阵的列；
- 矩阵 $A$ 的行列式 $det(A)$ 就是线性变换 $A$ 下的图形面积或体积的伸缩因子。
- 矩阵的行列式的几何意义是矩阵对应的线性变换前后的面积比

Property:

1. $Det(I) = 1$
2. Exchanging rows only reverse the sign of det
3. Determinant is "linear" for each row
4. $det(A) \neq 0$, $A$ is invertible
5. Cramer's rule: $A^{-1} = \frac{1}{det(A)}C^T$
  - $det(A)$: scalar
  - $C$: cofactors of $A$ (C has the same size as $A$)
  - $C^T$ is adjugate of $A$ (伴随矩阵)

### Eigenvalue and Eigenvector

Eigen (German word): "unique to" or "belonging to"

if $Av = \lambda v$ ($v$ is a vector , $\lambda$ is a scalar)

- $A$ must be square
- $v$ is an eigenvector of $A$, exluding zero vector
- $\lambda$ is an eigenvalue of $A$ that correponds to $v$

$T$ is a linear operator if $T(v) = \lambda v$ ( $v$ is a vector, $\lambda$ is a scalar)

- $v$ is an eigenvector of $T$, exluding zero vector
- $\lambda$ is an eigenvalue of $T$ that correponds to $v$

An eigenvector of A corresponds to a unique eigenvalue.  
An eigenvalue of A has infinitely many eigenvectors.  

=> how to find eigenvalues t :

$$
det(A - tI_n) = 0
$$

参考： 李航《统计学习方法》
