# Fourier transform for biologist


A biologist's way to learn Fourier transform

## Visual intuition in 3D

This is an awesome introduction

{{< youtube r18Gi8lSkfM >}}


## Fourier Series

Discrete Fourier transform (DFT)


A Fourier series is a periodic function composed of harmonically related sinusoids, combined by a weighted summation.

周期性函数可以变换为正余弦函数的和

$$
\begin{aligned}
f(t) &= \frac {a_0}{2} + \sum_n a_n \sin(n\omega t + \varphi_n) \cr
&=  \frac {a_0}{2} +  \sum_n a_n \sin(n\omega t) + \sum_n b_n \cos(n\omega t)
\end{aligned}
$$


Note: 3 orthogonal bases ( 1, sin, cos )

## Fourier Transform

### Euler's formula

For any real number $x$,

$$
e^{i\varphi} = \cos \varphi + i \sin \varphi
$$

$i$ is the imagenary unit. let $\varphi = \omega t$, get complex exponentials

$$
e^{i\omega t} = \cos (\omega t) + i \sin (\omega t)
$$

clockwise roation: $e^{i\omega t}$  
anti-clockwise roation: $e^{ - i\omega t}$

when $\varphi = \pi$, `Eluer's identity`

$$
e^{i\pi} + 1 = 0
$$



### Fourier transform for non-periodic function

For non-periodic function $f(t)$, we multiply $e^{ - i\omega t}$ to get signals only present in $e^{ - i\omega t}$, then

$$
\int_{-\infty}^{+\infty} f(t) e^{ - j\omega t} dt \rightarrow \begin{cases} = 0, \text{without } \omega \cr 
\neq 0, \text{with } \omega 
\end{cases}
$$

Note: the bases of $e^{ - j\omega t}$ are **orthogonal**, if a signal is orthogonal to these bases, their dot product equal to 0. 

- if a signal in $f(t)$ not present in $e^{ - j\omega t}$: $$\int_{-\infty}^{+\infty} f(t) e^{ - j\omega t} dt = 0$$
- if a signal in $f(t)$ present in $e^{ - j\omega t}$: $$\int_{-\infty}^{+\infty} f(t) e^{ - j\omega t} dt \neq 0$$

Now, the **Fourier transform definition**:

$$
F(\omega) = \int_{-\infty}^{+\infty} f(t) e^{ - j\omega t} dt
$$

and **inverse Fourier transform**:

$$
f(t) = \int_{-\infty}^{+\infty} F(\omega) e^{ - j\omega t}d \omega
$$

### 2D Fourier transform

$$
F(u, v) = \int_{-\infty}^{+\infty} \int_{-\infty}^{+\infty} f(x,y) e^{ - j(ux + vy)} dxdy
$$

Now, the **orthogonal bases** become (1, $\sin(ux + vy)$, $\cos(ux + vy)$), which could be further decomposited into: 
- $\sin(ux)\sin(vy)$
- $\sin(ux)\cos(vy)$
- $\cos(ux)\sin(vy)$
- $\cos(ux)\cos(vy)$


### Common application

1. audio: $t$ is the time domain
2. image: now $t$ become the location in the image
   - low freq: contour
   - high freq: detail

图像频率特性分析:

频谱图上的每一个像素点都代表一个频率值，幅值由像素点亮度变码而得。对于一幅图像，图像信号的频率特性如下：

* 直流分量: 表示预想的平均灰度
* 低频分量: 代表了大面积背景区域和缓慢变化部分
* 高频分量: 代表了它的边缘、细节、跳跃部分以及颗粒噪声
* 振幅: 描述了图像灰度的亮度
* 相位: 决定了图像是什么样子


## Laplace Transform 

Euler's Formula:

$$
e^{i\omega t} = \cos (\omega t) + i \sin (\omega t)
$$


we get 

$$
\cos(\omega t) = \frac{1}{2} (e^{j\omega t} + e^{-j\omega t})
$$


Fourier transform:

$$
F(\omega) = \int_{-\infty}^{+\infty} f(t) e^{ - j\omega t} dt
$$

The **problem** of Fourier transformation is that each component of sinusoids **keep constant magnitude** while oscillating.  

For $f(t)$ like $y = x^2$, when $x \rightarrow \infty$, FT do not perform well

To solve this problem, we could multiply a decay fuction $e^{-\sigma t}$, $\sigma > 0$. 


$$
\begin{aligned}
F(\omega) &= \int_{-\infty}^{+\infty} f(t) e^{-\sigma t} e^{ - j\omega t} dt \cr
&= \int_{-\infty}^{+\infty} f(t) e^{-t(\sigma + j\omega)} dt
\end{aligned}
$$

With a decay fuction, FT perform well when $x \rightarrow \infty$.

Let complex number $s = \sigma + j\omega$, then `Laplace Transform` is

$$
F(s) = \int_{-\infty}^{+\infty} f(t) e^{-st} dt
$$

**Laplace Transform**: 

- A generalized Fourier transform
- The **magnitude keep increasing/decreasing** as oscillation continue.
- if the real componet of $s$ is 0 ( $\sigma = 0$ ), then the **magnitude will stay constant**


**Inverse Laplace transform**:

$$
f(t) = \frac {1}{2 \pi i} \int_{c - i\infty}^{c + i\infty} F(s)e^{st} ds
$$


A visual intution of Laplace Transform:

{{<youtube 6MXMDrs6ZmA>}}

## Fourier Transform and Convolution

Convolution:
$$
h(x) = (f \star g)(x) = \int_{-\infty}^{+\infty} f(u)(g(x-u))du
$$

discrete form:

$$
h[n] = (f \star g)[n] = \sum_{u = -\infty}^{+\infty} f[u]g[n-u]
$$

- 时域卷积定理：时域上的卷积对应频域上的乘积 
  
  $$
  F[f(t) \star g(t)] = F_f(\omega) \cdot F_g (\omega)
  $$ 

- 频域卷积定理：频域内的卷积对应时域内的乘积

  $$
  F[f(t) \cdot g(t)] = \frac {1} {2\pi} F_f(\omega) \star F_g (\omega)
  $$ 
  

Examlple：

```python
import numpy as np
f=np.array([1,1,1])
g=np.array([2,3,2,6])
fg=np.convolve(f,g) # convolute
n = len(f)+len(g)-1
N = 2**(int(np.log2(n))+1)
a=np.fft.rfft(f,N) # make f and g have equal length N, and tranform
b=np.fft.rfft(g,N) 
c=a*b   # note here: just element wise multiplication after FT
fft_fg=np.fft.irfft(c)[:n] # inverse FT, only top n are valided.
print(fg)
print(fft_fg)
```
the output is

```
[2 5 7 11 8 6]
[2. 5. 7. 11. 8. 6.]
```

**Summary**:

Weighted sum => convolution => multiplication after fourier tranform.


A short animation explained what convolution is 

{{<youtube f0t-OCG79-U>}}


## Reference

1. [傅立叶变换](https://www.youtube.com/watch?v=0LuyxzqI3Hk)（李永乐老师）  
2. [Fourier Series](https://en.wikipedia.org/wiki/Fourier_series)
3. [Euler's formula](https://en.wikipedia.org/wiki/Euler%27s_formula)

![Euler](/images/ml/Euler_formula.svg.png)
