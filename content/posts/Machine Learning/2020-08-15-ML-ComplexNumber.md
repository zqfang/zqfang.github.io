---
title: "Complex number for biologist"
date: 2020-10-01
categories: ["Machine Learning"]
tags: ["Deep Learning", "Complex"]
comments: true
math: true
---

A biologist like me might have never had a numerical computing training. I don't even known what a complex number really means. Here are some useful basics to keep in mind.

## Complex number

complex number $a+bi$ lives in a 2d `complex plane`, including  

- real axis: $a$
- imagnary axis: $i$
  - `orthognal` to real axis 
  - $i \rightarrow 90 \degree \text{rotation}$


## 2 ways of representation 

- $z = a + bi$
- $z = r \cos(\phi) + r \sin(\phi) i = r e^{i \phi}$

![representation](/images/ml/complexnumber.png)

## 3 Facts about Multiplication

- $z \cdot 1 = z$
- $z \cdot i = \operatorname{Rot90}(z)$
  -  e.g. $i \cdot i = -1$
- $z \cdot ( c + di) = c \cdot z + d \cdot (zi)$
  - e.g. $(2+i)(2-i) = 2 \cdot 2 + 2i -2i - i^2 = 5 + 0i$


## expotential

form:
$$
\exp (i \theta)=1+i \theta+\frac{(i \theta)^{2}}{2}+\frac{(i \theta)^{3}}{6}+\frac{(i \theta)^{4}}{24}+\cdots
$$

derivative:

$$
\frac{d}{d t} e^{i t}=i \cdot e^{ i \cdot t}
$$


$$
i^n \cdot i^k = i^{n+k}
$$


## expotential form to find complex roots

Example:  
![representation2](/images/ml/complexnumber2.png)


## Reference

[e](https://en.wikipedia.org/wiki/E_(mathematical_constant))