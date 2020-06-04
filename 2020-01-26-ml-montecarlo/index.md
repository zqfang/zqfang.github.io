# Mento Carlos



蒙特卡罗方法，又称统计模拟方法(statistical simulation method), 通过概率模型的随机抽样进行进行近似数值计算的方法。

### 1. 蒙特卡罗方法的核心  

蒙特卡罗方法的核心是随机抽样(random sampling)

### 2. Monte Carlo intergration:

计算h(x)积分

$$
\int_{\mathcal{X}} h(x) \mathrm{d} x
$$

将h(x)分解成f(x)和概率密度函数p(x)的乘积，即函数h(x)的积分可以表示为函数f(x)关于概率密度函数p(x)的数学期望：

$$
\int_{\mathcal{X}} h(x) \mathrm{d} x=\int_{\mathcal{X}} f(x) p(x) \mathrm{d} x=E_{p(x)}[f(x)]
$$

因此，可利用样本均值计算近似积分：

$$
\int_{\mathcal{X}} h(x) \mathrm{d} x=E_{p(x)}[f(x)] \approx \frac{1}{n} \sum_{i=1}^{n} f\left(x_{i}\right)
$$

变形

$$
\begin{aligned}
\mathrm{E}_{p(z)}[f(z)] &=\int f(z) p(z) dz \\
&=\int \underbrace{f(z) \frac{p(z)}{q(z)}}_{new \tilde{f}(z)} q(z) dz \\ 
& \approx \frac{1}{N} \sum_{n=1}^{N} f(z^{i}) \frac{p(z^{i})}{q(z^{i})}
\end{aligned}
$$



### 3. Markov Chain Monte Carlo

Markov chain or markov process:

$$
P\left(X_{t} | X_{0}, X_{1}, \cdots, X_{t-1}\right)=P\left(X_{t} | X_{t-1}\right), \quad t=1,2, \cdots
$$

 


参考： 李航《统计学习方法》
