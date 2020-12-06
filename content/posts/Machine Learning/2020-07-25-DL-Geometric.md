---
title: "Geometric Deep Learning"
date: 2020-07-25
categories: ["Machine Learning"]
tags: ["Deep Learning","Graph"]
comments: true
math: true
---

Introduction of Graph Neural Networks

## Data

- Eculidean Structure Data: image, video, voice ...
  * easy to find adjacent neighbors
  * easy to define `distance`
- Non-Eculidean data: Graph, Manifold
  * hard to define adjacent neighbors
  * or the numbers of adjacent nodes varies.
  * means hard to define `distance`, `convolution` ...  

## Embed (project) Non-Eculidean Data into Eculidean Space

using geometric deep learning

## Graph Neural Network

### Common tasks

1. graph classification: classcify graphs according to its topology. each graph has a label.
   - most common
   - definition:  graph $G = (A, F)$ 
     - Adjacency matrix (of G): $A \in \lbrace 0,1 \rbrace ^{n \times n}$
     - Feature matrix (of nodes): $F \in R^{n \times d}$, with n nodes, d features
     - Given $\mathcal{D} = \lbrace (G_1, y_1), \cdots, (G_n, y_n) \rbrace$, learn
  
    $$
    f: \mathcal{G} \rightarrow \mathcal{Y}
    $$

2. node classification: each node has a label
3. generative graph (models): e.g. virtual drug screen

### Algebra Presentation of Graphs

1. Adjacency matrix

$$
A_{i j}= \begin{cases}
1 & \text { if }\lbrace v_{i}, v_{j}\rbrace \in E \text { and } i \neq j \cr
0 & \text { otherwise }
\end{cases}
$$

2. Degree matrix: D is a diagonal matrix, where

$$
D_{ii} = d(v_i)
$$

3. Laplacian matrix: if we consider all edges in graph $G$ to be undirected, then Laplacian matrix $L$ could be defined as

$$
L = D-A
$$

Thus, we have the elements:

$$
L_{i j}=\begin{cases}
d\left(v_{i}\right) & \text { if } i=j \cr
-1 & \text { if }\lbrace v_{i}, v_{j}\rbrace \in E \text { and } i \neq j \cr
0 & \text { otherwise. }
\end{cases}.
$$

What and why [Laplacian matrix](https://www.zhihu.com/question/54504471/answer/630639025) 

The key to understand Laplacian is PDE or Heat equation

$$
\frac{\partial T}{\partial t}(x, t)=\alpha \cdot \frac{\partial^{2} T}{\partial x^{2}}(x, t)
$$


4. Symmetric normalized Laplacian

the symmetric normalized Laplacian is define as:

$$
\begin{aligned}
L^{sym} &=D^{-\frac{1}{2}} L D^{-\frac{1}{2}} \cr
&=I-D^{-\frac{1}{2}} A D^{-\frac{1}{2}}
\end{aligned}
$$

The elements are given by:

$$
L_{i j}^{s y m}=\begin{cases}
1 & \text { if } i=j \text { and } d\left(v_{i}\right) \neq 0 \cr
-\frac{1}{\sqrt{d\left(v_{i}\right) d\left(v_{j}\right)}} & \text { if } \lbrace v_{i}, v_{j} \rbrace \in E \text { and } i \neq j \cr
0 & \text { otherwise. }
\end{cases}
$$

5. Random walk nmormalized Laplacian
  
$$
L^{rw} = D^{-1}L = I - D^{-1}A
$$

The elements can be computed by:

$$
L_{i j}^{r w}= \begin{cases}
1 & \text { if } i=j \text { and } d\left(v_{i}\right) \neq 0 \cr
-\frac{1}{d\left(v_{i}\right)} & \text { if }\lbrace v_{i}, v_{j}\rbrace \in E \text { and } i \neq j \cr
0 & \text { otherwise }
\end{cases}
$$

6. Incidence Matrix

$$
M_{i j}= \begin{cases}
1 & \text { if } \exists k \text { s.t } e_{j}=\lbrace v_{i}, v_{k} \rbrace \cr
-1 & \text { if } \exists k \text { s.t } e_{j}=\lbrace v_{k}, v_{i} \rbrace \cr
0 & \text { otherwise. }
\end{cases}
$$

for undrected graph, the corresponding incidence matrix statisfies that

$$
M_{i j}=\begin{cases}
1 & \text { if } \exists k \text { s.t } e_{j}=\lbrace v_{i}, v_{k}\rbrace \cr
0 & \text { otherwise. }
\end{cases}.
$$


## Convolution on Spectral

The key to understand graph convolution: 
- [Laplacian matrix](//www.zhihu.com/question/54504471/answer/630639025)   
- [Laplacian matrix wikepedia](https://en.wikipedia.org/wiki/Laplacian_matrix)  
- [Newton's law of cooling](https://en.wikipedia.org/wiki/Newton%27s_law_of_cooling)
- [Heat equation](https://en.wikipedia.org/wiki/Heat_equation) 

Now, the convolution operation is defined in the Fourier domain by computing the eigendecomposition of the  Laplacian Matrix

$$
L = U \Lambda U^{-1} = U \Lambda U^T
$$

Note: $U$ is an orthognal matrix, $U^{-1} = U^T$

Then, given $x \in R^n$, 

- the fourier transform: $\hat{x} = U^{T} x$ 
- reverse fourier transform: $x = U \hat{x}$ 

Finally, given signal $x$ and kernel $y$, the graph fourier transform ($*_{\mathcal{g}}$) is

$$
x *_{\mathcal{G}} y=U\left(\left(U^{T} x\right) \odot\left(U^{T} y\right)\right)
$$

$\odot$: element-wise multiplication

As we have a kernel $g_{\theta}(\sdot)$,

$$
y=g_{\theta}(L)(x)=g_{\theta}\left(U \Lambda U^{T}\right) x=U g_{\theta}(\Lambda) U^{T} x 
$$

where 

$$
g_{\theta}(\Lambda)=\operatorname{diag}(\theta)=\left[\begin{array}{ccc}
\theta_{1} & \cdots & 0 \cr
\vdots & \ddots & \vdots \cr
0 & \cdots & \theta_{n-1}
\end{array}\right]
$$

The learned parameters are in $\operatorname{diag}(\theta)$

The **problems**:

1. lost local connectivity on space (e.g. CNN on images preserve locality)
2. computational complexity $O(n)$, not well generalized on large scale Graphs 

Need more knowledge of the Chebyshev ploynomials to get deeper. see my next post about GNN.


## Convolution on Spatial

Another way to understand graph convolution: `Message Passing`.

### Message passing

Message passing: node $\mathcal{S}_1$ and its neigbor $\mathcal{N}$ (B1, B2, B3), aggregate $\mathcal{N}$'s message to $\mathcal{S}_1$.

for example, aggreate (sum) each nodes's features $H^{(l)} \in R^d$,

$$
\sum_{u \in \mathcal{N}(v)} H^{(l)}(u) \in \mathbb{R}^{d_{i}}
$$

node $v$'s neigbors: $\mathcal{N(v)}$, layer: $l$

![message passing](/images/ml/20190519220830445.png)

Generally, we add a linear transform matrix $W^{(l)} \in R^{d_i \times d_o}$ to change the feature dimension.

$$
\left(\sum_{u \in \mathcal{N}(v)} H^{(l)}(u)\right) W^{(l)} \in \mathbb{R}^{d_{o}}
$$

After add activate function, get a more compact equation

$$
f(H^{(l)}, A) = \sigma ( A H^{(l)}W^{(l)})
$$

$A$: Adjacency Matrix

### Example
Given graph
![graph](/images/ml/20190519223501518.png)

- input feature (10 dim): $f_{in} \in R^{10}$
- output feature (20 dim): $f_{out} \in R^{20}$
- each node's feature: $H^{(l)} \in R^{6 \times 10}$ 
- weight: $W^{(l)} \in R^{10 \times 20}$
- adjcency matrix: $A \in R^{6 \times 6}$
  
Message passing step:

1. feature dimension change: $HW \in R^{6 \times 20}$
2. select the neigborhood nodes: $AHW$

The **problems**:

- Each node have different degree, make the scale of output feature map will completely change (see each row of adjcency matrix). So, we have to `normalize laplacian matrix`.
- Each node did not include information from itself. So need to make a `self connection`.
  
### Definition

- **Adjcency Matrix**:
   
   $$
   \tilde{A} = A + I_n
   $$

- **Degree Matrix**:
   
  $$
  \tilde D_{ii} = \sum_{j} \tilde{A}_{ij}
  $$

- **Random Walk Normalization** of A: make row sum equal to 1
  
  $$
  \tilde{A} = D^{-1}A
  $$

  ![rwnorm](/images/ml/20190519225607797.png)

- **Symmetric Normalization**: used more in practice, more dynamic.

  $$
  A = D^{-\frac{1}{2}}AD^{-\frac{1}{2}}
  $$

  ![symmnorm](/images/ml/20190519225933458.png)

- **Laplacian matrix normliaztion**:
  
  $$
  \begin{aligned}
  L^{sym} &= D^{-\frac{1}{2}}LD^{-\frac{1}{2}} \cr 
  &= D^{-\frac{1}{2}}(D-A)D^{-\frac{1}{2}} \cr 
  &= I_n - D^{-\frac{1}{2}}AD^{-\frac{1}{2}} 
  \end{aligned}
  $$

Finally, we have

$$
H^{(l+1)} = \sigma ( \tilde{D}^{-\frac{1}{2}}\tilde{A}\tilde{D}^{-\frac{1}{2}}H^{(l)}W^{(l)})
$$

- $\tilde{A} = A+ I_n$
- $\tilde D_ii = \sum_j \tilde A_{ij}$