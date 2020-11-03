# Canonical Correlation Analysis (CCA)


典型相关分析(CCA) ，一种常用降维算法，也可以用于多个线性空间相关性计算。比如同一对象的多模态数据，单细胞多组学数据整合等。

## 原理

### 定义
假设有两个特征空间， $S1 = x_1 \in R^{d1}$, $S2 = x_2 \in R^{d2}$, 将两个特征向量合并，有

$$
\mathbf{x} = [\begin{array}{c}
\mathbf{x}_1 \cr
\mathbf{x}_2
\end{array} ] 
$$

$$
E(\mathbf{x})=[ \begin{array}{c}
\mu_{1} \cr
\mu_{2}
\end{array} ]
$$

$$
\Sigma=\bigg[ \begin{array}{cc}
\Sigma_{11} & \Sigma_{12} \cr
\Sigma_{21} & \Sigma_{22}
\end{array} \bigg]
$$

其中， $\Sigma$ 为协方差矩阵， 且$\Sigma_{12} = \Sigma_{21}^T$。

### 计算相关系数

令 $\mathbf{a}$, $\mathbf{b}$, 且满足

$$
u = \mathbf{a}^T \mathbf{x}_1 \\\ 
v = \mathbf{b}^T \mathbf{x}_2
$$

计算 u, v 的方差和协方差

$$
\operatorname{var}(u) = a^T \Sigma_{11} a \\\ 
\operatorname{var}(v) = b^T \Sigma_{22}b  \\\ 
\operatorname{cov}(u, v) = a^T \Sigma_{12} b
$$

因此，u 和 v的相关系数为：

$$
\operatorname{Corr}(u, v)=\frac{\operatorname{cov}(u, v)}{\sqrt{\operatorname{var}(u)} \sqrt{\operatorname{var}(v)}}
$$

代入 a， b

$$
\operatorname{Corr}(u, v)=\frac{\mathbf{a}^{T} \Sigma_{12} \mathbf{b}}{\sqrt{\mathbf{a}^{T} \Sigma_{11} \mathbf{a}} \sqrt{\mathbf{b}^{T} \Sigma_{22} \mathbf{b}}}
$$


## 求解

目标： 最大化 $\operatorname{Corr}(u,v)$

方法：

### 1. 固定分母，最大化分子

$$
\max_{\mathbf{a}, \mathbf{b}} \mathbf{a}^{T} \Sigma_{12} \mathbf{b}
$$

$$
\text { s.t. } \quad \mathbf{a}^{T} \Sigma_{11} \mathbf{a}=1, \quad \mathbf{b}^{T} \Sigma_{22} \mathbf{b}=1
$$

### 2. 构造拉格朗日等式

$$
L=\mathbf{a}^{T} \Sigma_{12} \mathbf{b}-\frac{\lambda_{1}}{2}\left(\mathbf{a}^{T} \Sigma_{11} \mathbf{a}-1\right)-\frac{\lambda_{2}}{2}\left(\mathbf{b}^{T} \Sigma_{22} \mathbf{b}-1\right)
$$

### 3. 求导

$$
\begin{array}{l}
\frac{\partial L}{\partial \mathbf{a}}=\Sigma_{12} \mathbf{b}-\lambda_{1} \Sigma_{11} \mathbf{a}=0 \\\ 
\frac{\partial L}{\partial \mathbf{b}}=\Sigma_{21} \mathbf{a}-\lambda_{2} \Sigma_{22} \mathbf{b}=0
\end{array}
$$

得： $\lambda_{1}=\lambda_{2}=\mathbf{a}^{T} \Sigma_{12} \mathbf{b}$

### 4. 求得最大相关系数

令 $\lambda = \lambda_1 = \lambda_2$， 有：

$$
\begin{array}{l}
\Sigma_{11}^{-1} \Sigma_{12} \mathbf{b}=\lambda \mathbf{a} \\\ 
\Sigma_{22}^{-1} \Sigma_{21} \mathbf{a}=\lambda \mathbf{b}
\end{array}
$$

即

$$
\left(\begin{array}{cc}
\Sigma_{11}^{-1} & 0 \\\ 
0 & \Sigma_{22}^{-1}
\end{array}\right)\left(\begin{array}{cc}
0 & \Sigma_{12} \\\ 
\Sigma_{21} & 0
\end{array}\right)\left(\begin{array}{l}
\mathbf{a} \\\ 
\mathbf{b}
\end{array}\right)=\lambda\left(\begin{array}{l}
\mathbf{a} \\\ 
\mathbf{b}
\end{array}\right)
$$

令

$$
B=\left(\begin{array}{cc}
\Sigma_{11} & 0 \\\ 
0 & \Sigma_{22}
\end{array}\right), \quad A=\left(\begin{array}{cc}
0 & \Sigma_{12} \\\ 
\Sigma_{21} & 0
\end{array}\right) \quad \mathbf{w}=\left(\begin{array}{l}
\mathbf{a} \\\ 
\mathbf{b}
\end{array}\right)
$$

上式表示为：

$$B^{-1}A\mathbf{w} = \lambda \mathbf{w}$$

因此，求解 $B^{-1}A$ 的特征值和特征向量，将原来的$x_1$, $x_2$做映射。

$\lambda$ 就是相关系数u和v的相关系数， u和v 就是一对典型变量(canonical variables).

### 求解方法2
从偏导数等式得到

$$
\Sigma_{11}^{-1} \Sigma_{12} \Sigma_{22}^{-1} \Sigma_{21} \mathbf{a}=\lambda^{2} \mathbf{a}
$$

因此，求出$\lambda$ 和 $\mathbf{a}$, 再代入，求$\mathbf{b}$.


## 参考

[典型关联分析](https://www.cnblogs.com/jerrylead/archive/2011/06/20/2085491.html)
