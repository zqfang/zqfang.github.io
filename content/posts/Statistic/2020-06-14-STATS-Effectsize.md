---
title: 'Power, Effect size, Sample size'
date: 2020-06-14
categories: ["Statistic"]
comments: true
math: true
draft: false
---

样本量、效应量、显著水平和统计功效的统计原理和计算  
效应量通常用三种方式来衡量：标准均差（standardized mean difference），几率（odd ratio），(3) 相关系数（correlation coefficient）。

## 统计原理
可视化样本量、效应量、α和统计功效的关系  
[Significance](https://rpsychologist.com/d3/nhst/)

The `probability` of a `type I error` is usually denoted by $\alpha$ and is commonly referred to as the signiﬁcance level of a test.  
\The `probability` of a `type II error` is usually denoted by $\beta$.


### 显著水平 α

| Type Error | 定义 |  表示 | 举例|
| ---  | --- | --- | --- |
| Ⅰ 型错误 | 拒绝实际上成立的$H_0$ | Ⅰ 型错误的概率用显著水平 α 表示, | 假阳性、误诊 |
| II 型错误 | 拒绝（“接受”）实际上不成立的$H_0$ | II 型错误概率用 β 表示 |  假阴性、漏诊 | 

 - α 常取值0.05、0.01，α可以取单尾、双尾。需假设检验前预先设定。
 - β 只取单尾

形象化理解  

![type-I-II-error](/images/stats/type-erros.jpeg)

### 功效（power）
功效（power）：正确拒绝原假设的概率，记作1-β, 即

$$
1 - \beta = \operatorname{Pr} ( \text{rejecting } H_0 | H_1 true)
$$

假设检验的功效受以下三个因素影响：

样本量 (n)：其他条件保持不变，样本量越大，功效就越大。  
显著性水平 (α)： 其他条件保持不变，显著性水平越低，功效就越小。  
两总体之间的差异：其他条件保持不变，总体参数的真实值和估计值之间的差异越大，功效就越大。也可以说，效应量（effect size）越大，功效就越大。

### 效应量（effect size）

效应量： 样本间差异或相关程度的量化指标。效应量越大，表示两个总体重叠的程度越小，效应越明显。


效应量通常用三种方式来衡量：(1) 标准均差（standardized mean difference），(2) 几率（odd ratio），(3) 相关系数（correlation coefficient）。

#### Difference family: Effect sizes based on differences between means

- 标准均差（standardized mean difference）
Standardized mean difference: 基于总体均值和方差， 效应量为

$$
\theta = \frac{\mu_1 - \mu_2}{\sigma}
$$

- Cohen’s d : 两总体均值之间的标准差异。适用于两组样本的样本量和方差相似的情况。

$$
d = \frac{ \bar{x}_1 - \bar{x}_2}{s} = \frac{\mu_1 - \mu_2}{\sigma}
$$

s是样本方差

$$
s = \sqrt{\frac{(n_1 - 1)s^2_1 + (n_2 - 1)s^2_2}{n_1 + n_2}}
$$

$$
s_1 = \frac{1}{n_1 -1} \sum_{i=1}^{n_1} (x_{1,i} - \bar{x}_1)^2
$$ 

d = 0.01 to 2.0


| Effect size | d | Reference|
| --- | --- | --- |
| Very small | 0.01 | Sawilowsky, 2009
| Small |0.20| Cohen, 1988
| Medium |0.50 | Cohen, 1988
| Large |0.80 | Cohen, 1988
| Very large | 1.20 | Sawilowsky, 2009
| Huge |2.0| Sawilowsky, 2009


- Hedges’ g: 是cohen的方法的改进，适用于两组样本的样本量不同的情况。

$$
g = \frac{ \bar{x}_1 - \bar{x}_2}{s^*}
$$

而$s^{*}$是

$$
s^* = \sqrt{\frac{(n_1 - 1)s^2_1 + (n_2 - 1)s^2_2}{n_1 + n_2 -2 }}
$$

- Glass’s Δ （delta）: 和cohen的方法类似，但是只除以控制组(control)的标准差。适用于两组样本的方差不同的情况。

$$
\Delta = \frac{ \bar{x}_1 - \bar{x}_2}{s_2}
$$


#### Categorical family: Effect sizes for associations among categorical variables

-  Odd ratio (OR)  
The odds ratio is the odds of success in the treatment group relative to the odds of success in the control group.  适用于binary数据。

- Relative risk (RR) or risk ratio
the risk (probability) of an event relative to some independent variable.

- Risk difference or absolute risk reduction  
the difference in risk (probability) of an event between two groups

- Cramer’s φ (Phi) or Cramer’s V: 用于测算类别型数据 (nominal data) 的效应量。当类别型变量包含2个类别时，使用Cramer’s phi，如果超过2个类别，那么使用Cramer’s V。

- Cohen's w 

...

#### Correlation family: Effect sizes based on "variance explained"
- Pearson r correlation  

|Effect size| r |
| --- | --- |
| small | ~ 0.1 |  
| medium |  ~ 0.3 |
| large | r > 0.5 |  


- Cohen’s $f^2$: 用于测算方差分析ANOVA，多元回归之类的效应量。

多元回归的效应量

$$
f^2 = \frac{R^2}{1-R^2}
$$

where $R^2$ is the squared multiple correlation ([Coefficient of determination](https://en.wikipedia.org/wiki/Coefficient_of_determination))


## 功效、效应量和样本量计算

### 计算样本量

[determining-sample-size](https://www.datasciencecentral.com/profiles/blogs/determining-sample-size-in-one-picture)
![determining-sample-size](/images/stats/sample-size-determination.png)


### 计算效应量

用`statsmodels`库计算功效，效应量和样本量的函数都是同一个，只要把需要计算的那个值仍然设为None，把其他想要达到的数值填上即可


单样本t检验：
```python
statsmodels.stats.power.tt_solve_power(effect_size=None, 
                                        nobs=None, 
                                        alpha=None, 
                                        power=None, 
                                        alternative='two-sided')
```

独立样本t检验：
```python
statsmodels.stats.power.tt_ind_solve_power(effect_size=None, 
                                            nobs1=None, 
                                            alpha=None, 
                                            power=None, 
                                            ratio=1.0, 
                                            alternative='two-sided')
```
卡方拟合优度检验：
```python
statsmodels.stats.power.GofChisquarePower.solve_power(effect_size=None,
                                                    nobs=None, 
                                                    alpha=None, 
                                                    power=None, 
                                                    n_bins=2)
```
F方差齐性检验：
```python
statsmodels.stats.power.FTestPower.solve_power(effect_size=None, 
                                                df_num=None, 
                                                df_denom=None, 
                                                nobs=None, 
                                                alpha=None, 
                                                power=None, 
                                                ncc=1)
```
方差分析：
```python
statsmodels.stats.power.FTestAnovaPower.solve_power(effect_size=None, 
                                                    nobs=None, 
                                                    alpha=None, 
                                                    power=None, 
                                                    k_groups=2)
```

## StatQuest

Sample size

{{< youtube 67zCIqdeXpo >}}


Power Analysis

{{< youtube VX_M3tIyiYk >}}


参考： 
1. https://en.wikipedia.org/wiki/Effect_size
2. https://www.cnblogs.com/HuZihu/p/12009535.html

