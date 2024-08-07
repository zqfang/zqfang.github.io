---
title: 'Markov Chain Mento Carlo (MCMC)'
date: 2020-06-06
categories: ["Machine Learning"]
comments: true
tags: ['Monte Carlo', "Statistical Learning"]
#markup: "mmark"
math: true
---


蒙特卡罗方法，又称统计模拟方法(statistical simulation method), 通过概率模型的随机抽样进行进行近似数值计算的方法。  
马可夫蒙特卡罗法（Markov Chain Monte Carlo, MCMC）则是以马可夫链为概率模型的蒙特卡罗方法。

Metropolis-Hastings算法是最基本的MCMC。
Gibbs sampling是更简单、使用更广泛的MCMC。


强烈建议先阅读这里（解释最清楚）： https://zhuanlan.zhihu.com/p/37121528


# Markov Chain Monte Carlo (MCMC)

## 为什么需要MCMC ？
1. 概率密度函数未知，累计分布函数没有反函数
2. 高维随机变量

## 蒙特卡罗法（Monte Carlo）
蒙特卡罗法要解决的问题是，假设概率分布的定义己知，**通过抽样获得概率分布的随机样本**，并通过得到的随机样本对概率分布的特征进行分析

蒙特卡罗方法于20世纪40年代美国在第二次世界大战中研制原子弹的“曼哈顿计划”计划时首先提出，为保密选择用赌城摩纳哥的Monte Carlo作为代号。

### 1. 蒙特卡罗方法的核心  

蒙特卡罗方法的核心是随机抽样(random sampling)
- 直接抽样
- 接受-拒绝抽样： 适用于概率密度函数复杂，不能直接抽样的情况
- 重要性抽样： 适用于概率密度函数复杂，不能直接抽样的情况


接受-拒绝抽样思想：找一个可以直接抽样的建议分布（proposal distribution），其概率密度函数为$q(x)$, 并且$q(x)$的$c$倍一定大于$p(x)$， 其中$c > 0$,按照$q(x)$进行抽样，假设得到结果$x^\ast$， 再按照$\frac{p(x^\ast)}{cq(x^\ast )}$的比例随机决定是否接受$x^\ast$。落到$p(x)$范围内的就接受，落到$p(x)$范围外的就拒绝❌。

![Sampling](/images/stats/rejsmp.png)

这些抽样方法的**缺点**： 
- 抽样效率低， 比如 $p(x)$ 占 $cq(x)$ 涵盖体积比例很低
- 当x为高维数据时，很难寻找合适的建议分布

一个解决办法就是MCMC.

### MCMC

MCMC是一种对高维随机向量抽样的方法， 此方法模拟一个马氏链， 使马氏链的平稳分布（(Stationary Distribution)）为目标分布， 由此产生大量的近似服从目标分布的样本， 但样本不是相互独立的。 MCMC的目标分布密度函数或概率函数可以只计算到差一个常数倍的值。

### 2. 数学期望估计(Estimation of mathematical expectation)

按照概率分布 $p(x)$ 独立抽取n个样本后计算函数的样本均值

$$
\hat f_{n}=\frac{1}{n} \sum_{i=1}^{n} f\left(x_{i}\right)
$$

作为数学期望的近似值。

根据大数定律可知，当样本容量增大是，样本均值以概率1收敛性于数学期望

$$
\hat f_{n} \rightarrow E_{p(x)}[f(x)], \quad n \rightarrow \infty
$$

于是，得到数学期望的近似计算方法

$$
E_{p(x)}[f(x)] \approx \frac{1}{n} \sum_{i=1}^{n} f\left(x_{i}\right)
$$


### 3. 蒙特卡罗积分（Monte carlo intergration）

计算函数 $h(x)$ 积分

$$
\int_{\mathcal{X}} h(x) \mathrm{d} x
$$

将 $h(x)$ 分解成 $f(x)$ 和概率密度函数 $p(x)$ 的乘积，即函数 $h(x)$ 的积分可以表示为函数 $f(x)$ 关于概率密度函数 $p(x)$ 的数学期望：

$$
\int_{\mathcal{X}} h(x) \mathrm{d} x=\int_{\mathcal{X}} f(x) p(x) \mathrm{d} x=E_{p(x)}[f(x)]
$$

因此，可利用样本均值计算近似积分：

$$
\int_{\mathcal{X}} h(x) \mathrm{d} x=E_{p(x)}[f(x)] \approx \frac{1}{n} \sum_{i=1}^{n} f\left(x_{i}\right)
$$

更进一步

$$
\begin{aligned}
E_{p(z)}[f(z)] &= \int f(z) p(z) dz \cr
&= \int \underbrace{f(z) \frac{p(z)}{q(z)}}_{new  \tilde{f} (z)} q(z) dz \cr 
& \approx \frac{1}{N} \sum_{n=1}^{N} f(z^{i}) \frac{p(z^{i})}{q(z^{i})}
\end{aligned}
$$


## Markov Chain
### 定义
#### 离散状态马可夫链
马可夫性：
随机变量$X_t$只依赖$X_{t-1}$，而不依赖过去的随机变量 $\lbrace X_{0}, X_{1}, \cdots, X_{t-2} \rbrace$。即

$$
P\left(X_{t} | X_{0}, X_{1}, \cdots, X_{t-1}\right)=P\left(X_{t} | X_{t-1}\right), \quad t=1,2, \cdots
$$

马可夫链或马可夫过程（markov process）指：
具有马可夫性的随机序列 $X=\lbrace X_{0}, X_{1}, \cdots, X_{t}, \cdots \rbrace$。

马可夫链的转移条件概率分布为 $P(X_t | X_{t-1})$ 。转移概率分布决定马可夫链的特性。


时间齐次马可夫链（time homogenous Markov Chain）是指转移状态分布于t无关的马可夫链

状态转移矩阵: $$P_{m \times m}$$

平稳分布：


马可夫链 $X$， 其状态空间为$\mathcal{S}$， 转移矩阵为 $P = (p_{ij})$， 如果存在状态空间 $\mathcal{S}$ 上的一个分布 

$$
\pi = \left[\begin{array}{c}
\pi_1 \cr
\pi_2 \cr
\vdots
\end{array}\right]
$$

使得 $\pi=P\pi$, 则称$\pi$为马可夫链$X = \{X_0, X_1, \cdots, X_t, \cdots \}$ 的平稳分布

#### 连续状态马可夫链
定义在连续状态空间，转移概率分布有概率转移核（trainsition kernel）表示

$$
P(x, A) = \int_{A} p(x, y) dy
$$

转移核$P(x, A)$表示转移概率

$$
P (X_t = A | X_{t-1} = x) = P (x, A)
$$

#### 马可夫链的性质

1. 不可约 (irreducible): 时刻 0 从状态 $j$ 出发，时刻 $t$ 到达状态 $i$ 的概率大于 0 ，则称此马尔可夫链 $X$ 是不可约的

$$
P(X_t = i | X_0 = j) > 0
$$

2. 非周期：不纯在一个状态，使得再返回到这个状态所经历的时间长呈周期性

3. 正常返(positive recurrent): 任意一个状态$i$，从其他任意状态 $j$ 出发，当时间趋近无穷时，首次转移到这个状态$i$的概率 $p^t_{ij}$ 不为0

$$
\lim_{t \rightarrow \infty} p^t_{ij} > 0
$$

4. 遍历定理：满足相应条件的马尔可夫链，当时间趋于无穷时，马尔可 夫链的状态分布趋近于平稳分布，随机变量的函数的样本均值以概率 1 收敛于该函数 的数学期望

马可夫链 $X$， 其状态空间为$\mathcal{S}$， 若马可夫链 $X$ 不可约、非周期且正常返， 则马可夫链有唯一的平稳分布 $\pi = (\pi_1, \pi_2, \cdots)^T$， 并且转移概率的极限分布是马可夫链的平稳分布

$$
\lim_{t \rightarrow \infty} P(X_t = i | X_0 = j) = \pi_i, i = 1,2, \cdots ; j = 1,2,\cdots
$$

若 $f(X)$是定义在状态空间上的函数 $E_{pi}[ | f(X) | ] <  \infty$, 则

$$
P \{ \hat{f_t} \rightarrow E_{pi}[  f(X) ] \} = 1
$$

这里，

$$
\hat{f_t} = \frac{1}{t} \sum^t_{s=1} f(x_s)
$$


关于平稳分布 $\pi = (\pi_1, \pi_2, \cdots)^T$ 的数学期望 $E_{pi}[f(X)] = \sum f(i)\pi_i$, 有

$$
 \hat{f_t} \rightarrow E_{pi}[  f(X) ], t \rightarrow \infty 
$$

处处成立或以概率1成立。

## Markov Chain Monte Carlo

马可夫蒙特卡罗法更适合随机变量是多元的、密度函数是非标准形式的、随机变量各分量不独立等情况

基本思想：

在随机变量$x$的状态空间$\mathcal{S}$上定一个满足遍历定理的马可夫链，使其平稳分布就是抽样的目标分布 $p(x)$， 然后在这个马可夫链上进行随机游走，每个时刻得到一个样本。根据遍历定理，当时间趋于无穷是，样本的分布趋近平稳分布，样本函数均值趋近函数的数学期望

$$
\hat{E}f = \frac{1}{n-m} \sum^{n}_{i=m+1} f(x_i)
$$


马尔可夫链蒙特卡罗法的关键是如何构建转移核函数或转移矩阵， 包括： Metropolis-Hastings 和 Gibbs sampling。

马尔可夫链蒙特卡罗法中得到的样本序列，**相邻的样本点是相关的**，而不是独立的  
马尔可夫链蒙特卡罗法的**收敛性的判断**通常是**经验性**的
 
基本步骤：

1) 在随机变量$x$的状态空间$\mathcal{S}$上构建一个满足遍历定理的马可夫链，使其平稳分布为目标分布 $p(x)$
2) 从状态空间的某一点 $X_0$ 出发，用构造的马可夫链进行随机游走，产生样本序列 $\{x_0, x_1, \cdots, x_t, \cdots \}$。

3) 应用马可夫链的遍历定理， 确定正整数 $m$ 和 $n$， $m < n$， 得到样本集合 $\{x_{m+1}, x_{m+2}， \cdots, x_n \}$ 求得函数f的（遍历）均值

$$
\hat{E}f = \frac{1}{n-m} \sum^{n}_{i=m+1} f(x_i)
$$

几个重要问题需要注意：
- 如何定义马可夫链，保证MCMC成立
- 如何确定收敛步数m，保证抽样无偏性
- 如何确定迭代步数n， 保证遍历均值的计算精度

### Metropolis-Hastings
Metropolis-Hasting是马尔可夫链蒙特卡罗法的代表算法

可以对多元变量的每一变量的条件分布依次分别进行抽样， 从而实现对整个多元变量的一次抽样，这就是单分量 Metropolis- Hastings (singlecomponent Metropolis- Hastings) 算法。


#### M-H采样python实现  
https://zhuanlan.zhihu.com/p/37121528

假设目标平稳分布是一个均值3，标准差2的正态分布，而选择的马尔可夫链状态转移矩阵 $Q(i,j)$ 的条件转移概率是以 $i$ 为均值,方差1的正态分布在位置 $j$ 的值。

```python

from scipy.stats import norm


def norm_dist_prob(theta):
    y = norm.pdf(theta, loc=3, scale=2)
    return y


T = 5000
pi = [0 for i in range(T)]
sigma = 1
t = 0
while t < T - 1:
    t = t + 1
    pi_star = norm.rvs(loc=pi[t - 1], scale=sigma, size=1,
                       random_state=None)  #状态转移进行随机抽样
    alpha = min(
        1, (norm_dist_prob(pi_star[0]) / norm_dist_prob(pi[t - 1])))  #alpha值

    u = random.uniform(0, 1)
    if u < alpha:
        pi[t] = pi_star[0]
    else:
        pi[t] = pi[t - 1]

plt.scatter(pi, norm.pdf(pi, loc=3, scale=2), label='Target Distribution')
num_bins = 50
plt.hist(pi,
         num_bins,
         density=1,
         facecolor='red',
         alpha=0.7,
         label='Samples Distribution')
plt.legend()
plt.show()

```

### Gibbs Sampling
吉布斯抽样，可以认为是 Metropolis-Hastings 算法的特殊情况，但是更容易实现，因而被广泛使用。

吉布斯抽样用于多元变量**联合分布**的抽样和估计。 其基本做法是，从联合概率分布定义满条件概率分布，依次对满条件概率分布进行抽样，得到样本的序列。


`吉布斯抽样`适合于**满条件概率分布** `容易抽样` 的情况，而`单分量MetropolisHastings` 算法适合于满条件概率分布`不容易抽样`的情况，这时使用容易抽样的条件分布作建议分布。



#### 二维Gibbs采样实例python实现  
(阅读：https://zhuanlan.zhihu.com/p/37121528)
假设我们要采样的是一个二维正态分布 $N(\mu, \Sigma)$ ，其中： $\mu=(\mu_{1}, \mu_{2})= (5, -1)$ , $\Sigma = \begin{pmatrix}
\sigma^{2}_{1} &   \rho \sigma_{1}\sigma_{2}b\rho \sigma_{2}& 
\sigma^{2}_{2}\end{pmatrix} = \begin{pmatrix}
 1& 1b1 & 
4\end{pmatrix}$;

而采样过程中的需要的状态转移条件分布为：

$P(x_{1}|x_{2}) = N(\mu_{1}+ \rho \sigma_{1}/\sigma_{2}(x_{2} - \mu_{2}), (1 - \rho^{2})\sigma^{2}_{1})$

$P(x_{2}|x_{1}) = N(\mu_{2}+ \rho \sigma_{2}/\sigma_{1}(x_{1} - \mu_{1}), (1 - \rho^{2})\sigma^{2}_{2})$



```python
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import multivariate_normal

samplesource = multivariate_normal(mean=[5,-1], cov=[[1,0.5],[0.5,2]])

def p_ygivenx(x, m1, m2, s1, s2):
    return (random.normalvariate(m2 + rho * s2 / s1 * (x - m1), math.sqrt(1 - rho ** 2) * s2))

def p_xgiveny(y, m1, m2, s1, s2):
    return (random.normalvariate(m1 + rho * s1 / s2 * (y - m2), math.sqrt(1 - rho ** 2) * s1))

N = 5000
K = 20
x_res = []
y_res = []
z_res = []
m1 = 5
m2 = -1
s1 = 1
s2 = 2

rho = 0.5
y = m2

for i in range(N):
    for j in range(K):
        x = p_xgiveny(y, m1, m2, s1, s2)   #y给定得到x的采样
        y = p_ygivenx(x, m1, m2, s1, s2)   #x给定得到y的采样
        z = samplesource.pdf([x,y])
        x_res.append(x)
        y_res.append(y)
        z_res.append(z)

num_bins = 50
plt.hist(x_res, num_bins,density=1, facecolor='green', alpha=0.5,label='x')
plt.hist(y_res, num_bins, density=1, facecolor='red', alpha=0.5,label='y')
plt.title('Histogram')
plt.legend()
plt.show()
```



参考：
1. 李航《统计学习方法》
2. https://zhuanlan.zhihu.com/p/37121528