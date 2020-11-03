---
title: "Complex number fundamentals for biologist"
date: 2020-08-01
categories: ["Machine Learning"]
tags: ["Deep Learning", "Complex"]
comments: true
math: true
draft: true
---


## Complex number

complex number lives in a 2d `complex plane`, where $a+bi$

- real axis: $a$
- imagnary axis: $i$
  - `orthognal` to real axis 
  - $i \rightarrow 90 \degree \text{rotation}$

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



