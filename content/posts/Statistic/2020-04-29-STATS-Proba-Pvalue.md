---
title: 'Probability and P-value'
date: 2020-04-29
categories: ["Statistic"]
tags: ["P-Value"]
comments: true
#markup: "mmark"
math: true
---

# Probability and P-value

## 1. PDF (probability density function)
1. PDF：概率密度函数（probability density function）, 连续型随机变量的概率密度函数是一个描述某个确定的取值点附近的可能性的函数。  
数学表示：用PDF在某一区间上的积分来刻画随机变量落在这个区间中的概率

$$
\operatorname{Pr}(a \leq X \leq b)=\int_{a}^{b} f_{X}(x) d x
$$

## 2.PMF (probability mass function)
PMF : 概率质量函数（probability mass function), 在概率论中，概率质量函数是离散随机变量在各特定取值上的概率。  
数学表示： PMF其实就是高中所学的离散型随机变量的分布律。  

$$
f_{X}(x)=\operatorname{Pr}(X=x)
$$

## 3. CDF (cumulative distribution function)
CDF : 累积分布函数 (cumulative distribution function)，是概率密度函数的积分，能完整描述一个实随机变量X的概率分布。

CDF是PDF的（从负无穷-oo到当前值的）积分，PDF是CDF的导数．（为了便于概率的计算，引入CDF的概念）  
CDF相当于其左侧的面积，也相当于小于该值的概率，负无穷的CDF值为０，正无穷的CDF值总为１．
对于连续变量，有

$$
F_{X}(x)=\operatorname{Pr}(X \leq x)=\int_{-\infty}^{x} f_{X}(t) d t
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


## 4. Central Limit Theorem
中心极限定理（Central Limit Theorem）

给定一个任意分布的总体，
每次从这些总体中随机抽取 n 个抽样，一共抽 m 次，
然后把这 m 组抽样分别求出平均值，
当m足够大时，这m次的平均值的分布（称为抽样分布）接近正态分布。

## 5. Confidence interval
置信区间（confidence interval）
置信区间是指由样本统计量所构造的总体参数的估计区间。
置信区间展现的是这个参数的真实值落在测量值（推测值）的周围的可信程度。

## 6. Probability and P-value

A `p-value` is the probability that random chance generated the data, or something else that is equal or rarer.

A p-value is composed of three parts:
1. The probability random chance would result in the observation.  
2. The probability of observing something else that is equally rare.
3. The probability of observing something rarer or more extreme. 


But `probability`

$$ Probalibility = \frac{Number of outcomes of interest} {The total number of outcomes}$$


In `hypothesis testing`, p-values are numbers, between 0 and 1, that, how small does a p-value have to be before we are confident that interested A is different from B.

## 7. Additional resources 

Statquest: `P Values, clearly explained`
{{< youtube 5Z9OIYA8He8 >}}


and `How to calculate P-value`

{{< youtube JQc3yx0-Q9E >}}
