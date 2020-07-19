# Boosting





提升（Boosting）方法： 通过改变训练样本的权重（概率分布），学习n个分类器，并将这些分类器线性组合，提高分类性能。

## 1. AdaBoost
AdaBoost通过提高被前一轮弱分类器错误分类样本的权值，从而降低被正确分类样本的权值，并采取加权多数表决的方法达到分类目的。

输入：训练数据集$T=\{(x_1, y_1), (x_2, y_2), \cdots, (x_N, y_N)\}$, $\mathcal{Y} = \{-1,+1\}$;  
输出：分类器$G(x)$

1). 初始化训练数据权值分布

$$D_1 = (w_{11}, \cdots, w_{1i}, \cdots, w_{1N}), w_{1i} = \frac{1}{N}, i = 1,2,\cdots,N$$  

2). 对 $m = 1，2，\cdots, M$  
a.对权值分布$D_m$的训练数据集学习，得到基本分类器 

$$
G_{m}(x): \mathcal{X} \rightarrow\{-1,+1\}
$$

b.计算$G(x)$在训练数据集上的分类误差率

$$
e_{m}=\sum_{i=1}^{N} P\left(G_{m}\left(x_{i}\right) \neq y_{i}\right)=\sum_{i=1}^{N} w_{m i} I\left(G_{m}\left(x_{i}\right) \neq y_{i}\right)
$$

c. 计算$G(x)$的系数

$$
\alpha_{m}=\frac{1}{2} \log \frac{1-e_{m}}{e_{m}}
$$

d. 更新训练数据的权值分布

$$
D_{m+1}=\left(w_{m+1,1}, \cdots, w_{m+1, i}, \cdots, w_{m+1, N}\right)
$$


$$
w_{m+1, i} = \frac{w_{m i}}{Z_{m}} \exp \left(-\alpha_{m} y_{i} G_{m}\left(x_{i}\right)\right), \quad i=1,2, \cdots, N
$$


其中，

$$
Z_{m}=\sum_{i=1}^{N} w_{m i} \exp \left(-\alpha_{m} y_{i} G_{m}\left(x_{i}\right)\right)
$$

3）构建基本线性分类器组合

$$
f(x)=\sum_{m=1}^{M} \alpha_{m} G_{m}(x)
$$

得到最终分类器

$$
\begin{aligned}
G(x) &=\operatorname{sign}(f(x)) \cr
&=\operatorname{sign}\left(\sum_{m=1}^{M} \alpha_{m} G_{m}(x)\right)
\end{aligned}
$$


### 1.1 AdaBoost算法误差
AdaBoost算法最终分类器训练误差界为

$$
\frac{1}{N} \sum_{i=1}^{N} I\left(G\left(x_{i}\right) \neq y_{i}\right) \leqslant \frac{1}{N} \sum_{i} \exp \left(-y_{i} f\left(x_{i}\right)\right)=\prod_{m} Z_{m}
$$

这一定理说明，每一轮选取适当的$G_m$使$Z_m$最小，从而使训练误差下降最快。

对于二分类问题：

$$
\begin{aligned}
\prod_{m=1}^{M} Z_{m} &=\prod_{m=1}^{M}[2 \sqrt{e_{m}\left(1-e_{m}\right)}] \cr
&=\prod_{m=1}^{M} \sqrt{\left(1-4 \gamma_{m}^{2}\right)} \cr
& \leqslant \exp \left(-2 \sum_{m=1}^{M} \gamma_{m}^{2}\right)
\end{aligned}
$$

其中， $\gamma_{m}=\frac{1}{2}-e_{m}$

### 1.2 AdaBoost算法解释
AdaBooost可以认为：模型为加法模型，损失函数为指数函数，学习算法为前向分布算法的二分类学习方法

#### 1.2.1 前向分步算法

考虑加法模型（additive model）

$$
f(x)=\sum_{m=1}^{M} \beta_{m} b\left(x ; \gamma_{m}\right)
$$

其中，$b(x; \gamma_m)$为基函数，$gamma_m$为参数， $\beta_m$为系数。

在给定训练集和损失函数$L(y,f(x))$的条件下，学习加法模型$f(x)$成为经验风险极小化（损失函数极小化）问题：

$$
\min_{\beta_{m}, \gamma_{m}} \sum_{i=1}^{N} L\left(y_{i}, \sum_{m=1}^{M} \beta_{m} b\left(x_{i} ; \gamma_{m}\right)\right)
$$

**前向分布算法思想**是： 从前向后，每一步只学一个基函数及其系数，逐步逼近优化目标函数，达到优化步骤简化的目的。

因此，每一步只需优化如下损失函数：

$$
\min_{\beta, \gamma} \sum_{i=1}^{N} L\left(y_{i}, \beta b\left(x_{i} ; \gamma\right)\right)
$$

**算法步骤**

输入：训练数据集$T=\lbrace (x_1, y_1), (x_2, y_2), \cdots, (x_N, y_N)\rbrace$, 损失函数$L(y,f(x))$;基函数集$\lbrace b(x;\gamma) \rbrace$;  
输出：加法模型$f(x)$

1）初始化$f_0(x) = 0$  
2) 对$m = 1,2,\cdots, M$  
a.极小化损失函数

$$
\left(\beta_{m}, \gamma_{m}\right)=\arg \min _{\beta, \gamma} \sum_{i=1}^{N} L\left(y_{i}, f_{m-1}\left(x_{i}\right)+\beta b\left(x_{i} ; \gamma\right)\right)
$$

得到参数$\beta_m$, $\gamma_m$。

b.更新

$$
f_{m}(x)=f_{m-1}(x)+\beta_{m} b\left(x ; \gamma_{m}\right)
$$

3）得到加法模型

$$
f(x)=f_{M}(x)=\sum_{m=1}^{M} \beta_{m} b\left(x ; \gamma_{m}\right)
$$

## 2. Boosting Tree

提升树🌲是以决策树为基本分类器的提升方法

### 2.1 提升树模型
采用加法模型（基函数的线性组合）与前向分布算法：

$$
f_{M}(x)=\sum_{m=1}^{M} T\left(x ; \Theta_{m}\right)
$$

其中 $T\left(x ; \Theta_{m}\right)$表示决策树，$\Theta_{m}$决策树参数， $M$为树的个数

### 2.2 提升树算法
采用**加法模型**和**前向分布算法**实现学习优化的过程。

首先确定提升树$f_{0}(x)=0$， 第$m$步的模型是

$$
f_{m}(x)=f_{m-1}(x)+T\left(x ; \Theta_{m}\right)
$$

其中， $f_{m-1}(x)$为当前模型，通过经验风险极小化确定下一刻决策树的参数$\Theta_{m}$：

$$
\hat \Theta_m = \arg \min_{\Theta_{m}} \sum_{i=1}^{N} L(y_{i}, f_{m-1} (x_{i})+T (x_{i} ; \Theta_{m} ))
$$

#### 2.2.1 回归问题提升树

训练数据集:
$T=\lbrace (x_1, y_1), (x_2, y_2), \cdots, (x_N, y_N)\rbrace$, $x_{i} \in \mathcal{X} \subseteq \mathbf{R}^{n}$, $\mathcal{X}$为输入空间， $\mathcal{Y} \subseteq \mathbf{R}$;  

将输入空间划分为$J$个互不相交的区域$R1，R2, \cdots, R_J$， 并且每个区域上确定输出的常量$c_j$，那么树可以表示为：

$$
T(x ; \Theta)=\sum_{j=1}^{J} c_{j} I\left(x \in R_{j}\right)
$$

其中， 

$$
\Theta=\lbrace \left(R_{1}, c_{1}\right),\left(R_{2}, c_{2}\right), \cdots,\left(R_{J}, c_{J}\right)\rbrace
$$
表示树的却与划分和各个取悦是那个的常数。


采用一下前向分布算法

$$
\begin{aligned}
&f_{0}(x)=0\cr
&\begin{array}{l}
f_{m}(x)=f_{m-1}(x)+T\left(x ; \Theta_{m}\right), \quad m=1,2, \cdots, M \cr
f_{M}(x)=\sum_{m=1}^{M} T\left(x ; \Theta_{m}\right)
\end{array}
\end{aligned}
$$

求解$\hat \Theta_{m}$，
若用平方误差损失函数：

$$
L(y, f(x))=(y-f(x))^{2}
$$

则损失函数为：

$$
\begin{aligned}
L\left(y, f_{m-1}(x)+T\left(x ; \Theta_{m}\right)\right) &=\left[y-f_{m-1}(x)-T\left(x ; \Theta_{m}\right)\right]^{2} \cr
&=\left[r-T\left(x ; \Theta_{m}\right)\right]^{2}
\end{aligned}
$$

这里， 

$$
r=y-f_{m-1}(x)
$$

是当前模型拟合数据的残差（residual）。因此对于回归问题提升树，只需拟合当前模型残差。得到$T\left(x ; \Theta_{m}\right)$，更新模型，得到$f_m(x)$。

## 3. 梯度提升
当损失函数不是简单的平方损失、指数损失时，提升树的优化就很难。梯度提升算法利用最速下降法的近似方法，计算损失函数的负梯度在当前模型的值

$$
-\left[\frac{\partial L\left(y, f\left(x_{i}\right)\right)}{\partial f\left(x_{i}\right)}\right]_{f(x)=f_{m-1}(x)}
$$

并将其作为回归问题提升树算法中的残差近似值，拟合一个回归树。


输入： 训练数据集$T=\lbrace (x_1, y_1), (x_2, y_2), \cdots, (x_N, y_N)\rbrace$, $x_{i} \in \mathcal{X} \subseteq \mathbf{R}^{n}$,$\mathcal{X}$为输入空间， $\mathcal{Y} \subseteq \mathbf{R}$; 损失函数$L(y,f(x))$  
 输出： 回归树$\hat f(x)$

1) 初始化

$$
f_{0}(x)=\arg \min _{c} \sum_{i=1}^{N} L\left(y_{i}, c\right)
$$

2) 对 $m=1，2，\cdots, M$

  (1) 对 $i=1，2，\cdots, N$计算

$$
r_{m i}=-\left[\frac{\partial L\left(y_{i}, f\left(x_{i}\right)\right)}{\partial f\left(x_{i}\right)}\right]_{f(x)=f_{m-1}(x)}
$$

  (2) 对$r_{mi}$拟合一个回归树，得到第$m$颗树的节点区域$R_{mj}$

  (3) 对$j=1,2,\cdots, J$, 计算

$$
c_{m j}=\arg \min _{c} \sum_{x_{i} \in R_{m j}} L\left(y_{i}, f_{m-1}\left(x_{i}\right)+c\right)
$$

  (4)更新

$$
f_{m}(x)=f_{m-1}(x)+\sum_{j=1}^{J} c_{m j} I\left(x \in R_{m j}\right)
$$

3) 得到回归树

$$
\hat{f}(x)=f_{M}(x)=\sum_{m=1}^{M} \sum_{j=1}^{J} c_{m j} I\left(x \in R_{m j}\right)
$$




参考： 李航《统计学习方法》
