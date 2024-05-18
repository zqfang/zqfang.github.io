---
title: 'Probability, P-value, Likelihood'
date: 2020-04-29
categories: ["Statistic"]
tags: ["P-Value"]
comments: true
#markup: "mmark"
math: true
---


# Probability, P-value, Likelihood

1. **Probability**: the level of possibility of something happening or being true.
   - We determine the possibility of an event.
   - We know the parameters associated with the event and assume them to be trustworthy.

2. **Likelihood**: the chance that something will happen.
   - We have some observations.
   - We have an explanation (or parameters).
   - Likelihood helps us quantify whether the explanation is trustworthy.

If you notice closely, “likelihood” is the only synonym of “probability”.

## Probability and likelihood 

[likehood & maximum likehood](http://fangs.in/post/thinkstats/likelihood/)

在非正式场合似然（likelihood）和概率（Probability）几乎是一对同义词，但是在统计学中似然和概率却是两个不同的概念。

**概率**: 在特定环境下某件事情发生的可能性，也就是结果没有产生之前依据环境所对应的参数来预测某件事情发生的可能性。  
>比如抛硬币，抛之前我们不知道最后是哪一面朝上，但是根据硬币的性质我们可以推测任何一面朝上的可能性均为50%，这个概率只有在抛硬币之前才是有意义的，抛完硬币后的结果便是确定的；

**似然**: 刚好相反，是在确定的结果下去推测产生这个结果的可能环境（参数）。  
>假设随机抛掷一枚硬币1,000次，结果500次人头朝上，500次数字朝上，那么两面朝上的概率均为50%。<span style="color: red">运用出现的结果来判断这个事情本身的性质（参数）</span>，也就是似然。

**当结果和参数相互对应，似然和概率在数值上相等**。 用 θ 表示环境对应的参数，x 表示结果，那么概率可以表示为：

$$P(x | \theta )$$  

$p(x \vert θ)$ 是条件概率的表示方法。θ 是前置条件，理解为在 θ 的前提下，事件 x 发生的概率，相对应的似然可以表示为:  

$$\mathcal{L}(\theta | x)$$  

可以理解为已知结果为 x ，参数为 θ (似然函数里 θ 是变量，这里说的参数和变量是相对与概率而言的)对应的概率，即：  

$$\mathcal{L}(\theta | x)=P(x | \theta)$$

两者在数值上相等，但是意义并不相同, $\mathcal{L}$ 是关于 θ 的函数，而 P 则是关于 x 的函数。

## Probability and P-value

A `p-value` is the probability that random chance generated the data, or something else that is equal or rarer.

A p-value is composed of three parts:

1. The probability random chance would result in the observation.  
2. The probability of observing something else that is equally rare.
3. The probability of observing something rarer or more extreme. 

But `probability`

$$ \text{Probalibility} = \frac{ \text{Number of outcomes of interest}} { \text{The total number of outcomes}}$$


In `hypothesis testing`, p-values are numbers, between 0 and 1, that, how small does a p-value have to be before we are confident that interested A is different from B.

Statquest: `P Values, clearly explained`  

{{< youtube 5Z9OIYA8He8 >}}


### PDF (probability density function)

PDF：概率密度函数（probability density function）, 连续型随机变量的概率密度函数是一个描述某个确定的取值点附近的可能性的函数。  
数学表示：用PDF在某一区间上的积分来刻画随机变量落在这个区间中的概率

$$
\operatorname{Pr}(a \leq X \leq b)=\int_{a}^{b} f_{X}(x) d x
$$

### PMF (probability mass function)

PMF : 概率质量函数（probability mass function), 在概率论中，概率质量函数是离散随机变量在各特定取值上的概率。  
数学表示： PMF其实就是高中所学的离散型随机变量的分布律。  

$$
f_{X}(x)=\operatorname{Pr}(X=x)
$$

### CDF (cumulative distribution function)

CDF : 累积分布函数 (cumulative distribution function)，是概率密度函数的积分，能完整描述一个实随机变量X的概率分布。

CDF是PDF的（从负无穷$-\infty$到当前值的）`积分`，PDF是CDF的`导数`．（为了便于概率的计算，引入CDF的概念）  
CDF相当于其左侧的`面积`，也相当于小于该值的概率，负无穷的CDF值为０，正无穷的CDF值总为１．
对于连续变量，有

$$
F_{X}(x)=\operatorname{Pr}(X \leq x)=\int_{-\infty}^{x} f_{X}(t) dt
$$

对于离散型变量，有如

$$
F_{X}(x)=\operatorname{Pr}(X \leq x)= 
\begin{cases}
0 \text { if } x<0 \cr
\frac{1}{2} \text { if } 0 \leq x<1 \cr
1 \text { if } x \geq 1
\end{cases}
$$


### Central Limit Theorem

中心极限定理（Central Limit Theorem）

给定一个任意分布的总体，
每次从这些总体中随机抽取 n 个抽样，一共抽 m 次，
然后把这 m 组抽样分别求出平均值，
当m足够大时，这m次的平均值的分布（称为抽样分布）接近正态分布。

**独立同分布**的中心极限定理

$$
\lim_{n \rightarrow \infty} F_{\mathcal{X}}(x) = \lim_{n \rightarrow \infty} P \Bigg\lbrace \frac{\sum_{k=1}^{n}X_k - n\mu}{\sqrt{n}\sigma} \leq x \Bigg\rbrace = \int_{-\infty}^{x} \frac{1}{\sqrt{2\pi}}e^{-\frac{t^2}{2}}dt
$$

德莫佛－拉普拉斯定理: 设随机变量序列$\lbrace \eta_1, \eta_2,\cdots, \eta_n \rbrace$ 服从参数为$n, p (0 < p < 1)$ 的**二项分布**


$$
\lim_{n \rightarrow +\infty} P \Bigg\lbrace \frac{\eta_n - np }{\sqrt{np(1-p)}} \leq x \Bigg\rbrace = \int_{-\infty}^x \frac{1}{\sqrt{2\pi}}e^{-\frac{t^2}{2}}dt
$$




### Law of Large Numbers 
当样本数据无限大时，样本均值趋于总体均值

$$
\bar{X} = \frac{1}{n} \sum_{k=1}^n X_k
\xrightarrow{p} \mu
$$

大数定律告诉我们能用频率近似代替概率；能用样本均值近似代替总体均值。


**辛钦大数定律:**

设$X_1, X_2, \cdots, X_n$是相互独立且服从同分布的随机变量序列， 具有数学期望$E(X_k) = \mu$, $k=1,2,3,\cdots$。对于任意 $\epsilon > 0$, 有

$$
\lim_{n \rightarrow +\infty} P \bigg\lbrace | \frac{1}{n} \sum_{k=1}^n X_k - \mu | < \epsilon \bigg\rbrace = 1
$$

**切比雪夫大数定律:** 随机变量序列$X$具有相同期望和方差, 样本均值依概率$p$收敛于 $\mu$

$$
\frac{1}{n} \sum_{k=1}^n X_k
\xrightarrow{p} \mu
$$

**伯努利大数定律:** $n_A$是n次独立重复试验中事件A发生的次数， p是事件A在每次试验中发生的概率，任意$\epsilon > 0$

$$
\lim_{n \rightarrow +\infty} P \bigg\lbrace | \frac{n_A}{n} - p | < \epsilon  \bigg\rbrace = 1
$$


比较

| 定律 | 分布 | 期望 | 方差	| 结论 |  
| --- | --- | --- | --- | --- |  
| 辛钦大数定律 | 相互独立且同分布 | 存在 | | 估算期望 |  
| 切比雪夫大数定律 | 相互独立 | 相同 | 相同 | 估算期望 |  
| 伯努利大数定律 | 二项分布 | 相同 | 相同 | 频率=概率 |  
| 相同点：$n \rightarrow +\infty$, 依概率趋近 | 条件组件变得严格 |

### Confidence interval

置信区间（confidence interval）  
置信区间是指由样本统计量所构造的总体参数的估计区间。  
置信区间展现的是这个参数的真实值落在测量值（推测值）的周围的可信程度。

## Additional: StatQuest 

How to calculate `P-value`

{{< youtube JQc3yx0-Q9E >}}
