# Graph: GNN review


More about Graph Neural Network

## Algebra presentation of Graphs

### 1. Adjacency matrix

$$
A_{i j}= \begin{cases}
1 & \text { if }\lbrace v_{i}, v_{j}\rbrace \in E \text { and } i \neq j \cr
0 & \text { otherwise }
\end{cases}
$$

### 2. Degree matrix: D is a diagonal matrix, where

$$
D_{ii} = d(v_i)
$$

### 3. Laplacian matrix

What and why [Laplacian matrix](https://www.zhihu.com/question/54504471/answer/630639025)

if we consider all edges in graph $G$ to be undirected, then Laplacian matrix $L$ could be defined as

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

### 4. Symmetric normalized Laplacian

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

### 5. Random walk nmormalized Laplacian
  
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

### 6. Incidence Matrix

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


## Vanilla Graph Neural Networks

A node is defined by its features and related nodes in the graph. The aim of GNN is to lean a state embedding $h_v \in R^s$, which encodes the information of the neighborhood, for each node. The state embedding $h_v$ is used to produce an output $O_v$, such as the distribution of the predicted node lable.

### Model


$$
h_v = f(X_v, X_{co[v]}, h_{ne[v]}, X_{ne[v]})
$$

$$
o_v = g(h_v, X_v)
$$


- $f$: local transition function, shared amoing all nodes
- $g$: local output function
- $X$: the input feature
- $h$: hidden state
- $co[v]$: the set of edges connected to node v
- $ne[v]$: the set of neighbors of node $v$

Now, we have a compact form as 

$$
H = F(H,X) \\\ O = G(H, X_N)
$$

- $F$: the global transition function
- $G$: the global output funciton

GNN use the classic iterative scheme to compute the state

$$
H^{t+1} = F(H^t, X)
$$

Next question is how to learn parameters of $f$ and $g$. The loss can be written as 

$$
loss = \sum_{i = 1}^p(t_i - o_i)
$$

where 

- $p$: the number of supervised nodes
- the state $h_{v}^{t}$ are iteratively updated until a time step $T$.

## Graph Convolutional Networks

Spectral approaches and spatial approaches. Four calssic models (Spectral Network, ChebNet, GCN, and AGCN)

### Spectral approaches 
#### Spectral Network

the convolution operation is defined in the Fourier domain by computing the eigendecomposition of the graph Laplacian. The operation can be defined as multiplication of a signal x (a scalar for each node) with a filter $g_{\theta} = \text{diag}(\theta)$:



$$
\mathbf{g_{\theta}}  \star \mathbf{x}=\mathbf{U g}_{\theta}(\Lambda) \mathbf{U}^{T} \mathbf{x}
$$

- $\mathbf{U}$: the matrix of eigenvectors of the normalized graph Laplacian $L = I_N - D^{-\frac{1}{2}}AD^{-\frac{1}{2}} = \mathbf{U}\Lambda\mathbf{U}^T$
- This operation results in potentially intense computations and non-spatially localized filters

#### CHEBNET

the operation is 

$$
\mathbf{g_\theta} \star \mathbf{x} \approx \sum_{k=0}^{K} \boldsymbol{\theta}_{k} \mathbf{T}_{k}(\tilde{\mathbf{L}}) \mathbf{x}
$$

- $\tilde{\mathbf{L}} = \frac{2}{\lambda_{max}}L - I_N$
- $\lambda_{max}$: the largest eigenvalue of $L$.
- $\theta$: a vector of Chebyshev coefficients
- $T_{k}(x)$: the Chebyshev ploynomials $T_k(x) = 2xT_{k-1}(x) -T_{k-2}(x)$. $T_0(x) = 1$, $T_1(x) = x$.

#### GCN

$$
\mathbf{g_\theta^{\prime}} \star \mathbf{x} \approx \theta_{0}^{\prime} \mathbf{x}+\theta_{1}^{\prime}\left(\mathbf{L}-\mathbf{I}_{N}\right) \mathbf{x}=\theta_{0}^{\prime} \mathbf{x}-\theta_{1}^{\prime} \mathbf{D}^{-\frac{1}{2}} \mathbf{A} \mathbf{D}^{-\frac{1}{2}} \mathbf{x}
$$

then, constraining the number of parameters with $\theta = \theta_{0}^{\prime} = - \theta_{1}^{\prime}$, get

$$
\mathbf{g_\theta} \star \mathbf{x} \approx \theta\left(\mathbf{I}_{N}+\mathbf{D}^{-\frac{1}{2}} \mathbf{A} \mathbf{D}^{-\frac{1}{2}}\right) \mathbf{x}
$$

this operator could lead to numberical instabilities and exploding/vanishing gradients. Finally, introudce the *renormalization trick*: $I_N + \mathbf{D}^{-\frac{1}{2}} \mathbf{A} \mathbf{D}^{-\frac{1}{2}} \rightarrow \tilde{\mathbf{D}}^{-\frac{1}{2}} \tilde{\mathbf{A}} \tilde{\mathbf{D}}^{-\frac{1}{2}}$, with $\tilde{A} = A + I_N$, and $\tilde{D} = \sum_j\tilde{A}_{ij}$. Finally

$$
\mathbf{Z}=\tilde{\mathbf{D}}^{-\frac{1}{2}} \tilde{\mathbf{A}} \tilde{\mathbf{D}}^{-\frac{1}{2}} \mathbf{X} \Theta
$$

- $\Theta \in \mathbb{R}^{C \times F}$: a matrix of filter paramters
-  $\mathbf{Z} \in \mathbb{R}^{N \times F}$: the convolved signal matrix

#### AGCN
Adaptive Graph Convolution Network (AGCN) is proposed to learn the underlying relations. AGCN learns a "residual" graph Laplacian $\mathbf{L}_{res}$ and add it to the original Lapalcian matrix

$$
\widehat{\mathbf{L}}=\mathbf{L}+\alpha \mathbf{L}_{r e s}
$$

- $\mathbf{L}_{res}$ is computed by learned graph adjacency matrix $\widehat{\mathbf{A}}$

$$
\begin{aligned}
\mathbf{L_{res}} &=\mathbf{I}-\widehat{\mathbf{D}}^{-\frac{1}{2}} \widehat{\mathbf{A}} \widehat{\mathbf{D}}^{-\frac{1}{2}} \cr 
\widehat{\mathbf{D}} &=\operatorname{degree}(\widehat{\mathbf{A}})
\end{aligned}
$$

- $\widehat{\mathbf{A}}$ is computed via a learned metric


The idea behind the adaptive metric is that Euclidean distance is not suitatble for graph structured data and the metric should be adaptive the the task and input features. ACGN use the generalized Mahalanobis distance

$$
D(\mathbf x_{i}, \mathbf x_{j})= \sqrt{ (\mathbf x_{i}-\mathbf x_{j})^{T} \mathbf{M}(\mathbf x_{i}-\mathbf x_{j})}
$$

where M is a learned prameter taht statisfies $M = W_d W_{d}^T$.


### Spatial Methods
spatial approches defined convolutions directly on the graph, operating on spatially close neighbors. The major challenge of spatial approches is defining the convolution operation with differently sized neighborhoods and maintaining the local invariance of CNNs.

#### Nerual FPS

Use different weight matris for nodes with different degrees

$$
\begin{aligned}
\mathbf{x} &=\mathbf{h}_{v}^{t-1}+\sum_{i=1}^{\left|N_{v}\right|} \mathbf{h}_{i}^{t-1} \\
\mathbf{h}_{v}^{t} &=\sigma\left(\mathbf{x} \mathbf{W}_{t}^{\left|N_{v}\right|}\right)
\end{aligned}
$$

#### Patchy-SAN
first selects and normalized exactly k neigbors for each node, then normalized neighhorhos servers as the receptive filed and the convolutional operation is applied. 

1. Node sequence selection
2. Neighorhood assembly
3. Graph normalization
4. Convolution architecture.

![patchy-san](/images/ml/patchy-san.png)
   
#### DCNN

the diffusion-convolution neural network (DCNN)

for node classification

$$
H = \sigma(W^c \odot P^*X)
$$

- P*: N x K x N tensor, contains the power series $\lbrace P, P^2, \cdots, P^K \rbrace$
- X: N x F tensor of input features

for graph classification

$$
H = \sigma(W^c \odot 1^T_N P^* X /N)
$$

$1_N$: N x 1 vector of ones.

#### DGCN
dual grpah convolutional network, which jointly consider the local consistency and global consitency on graphs. It use 2 convolutional networks to capture the local/global consistency and adopts an unsupervised loss to ensemble them.

first network:

$$
\mathbf{Z}=\tilde{\mathbf{D}}^{-\frac{1}{2}} \tilde{\mathbf{A}} \tilde{\mathbf{D}}^{-\frac{1}{2}} \mathbf{X} \Theta
$$

second network replace the adjacency matrix with positive pointwise mutual information (PPMI) matrix

$$
\mathbf H^{\prime}=\sigma(\mathbf D_{P}^{-\frac{1}{2}} \mathbf{X}_{P} \mathbf D_{P}^{-\frac{1}{2}} \mathbf H \Theta)
$$

Assemble these two convolutions via the final loss function

$$
L=L_{0}\left(\operatorname{Conv}_{A}\right)+\lambda(t) L_{r e g}\left(\operatorname{Conv}_{A}, \operatorname{Conv}_{P}\right)
$$

$\lambda(t)$ is the dynamic weight to blance the importance of these two loss functions.


$$
L_{0}\left(\operatorname{Conv}_{A}\right)=-\frac{1}{\left|y_{L}\right|} \sum_{l \in y_{L}} \sum_{i=1}^{c} Y_{l, i} \ln \left(\widehat{Z}_{l, i}^{A}\right)
$$

$$
L_{r e g}\left(\operatorname{Conv}_{A}, \operatorname{Conv}_{P}\right)=\frac{1}{n} \sum_{i=1}^{n}\left\|\widehat{Z}_{i,:}^{P}-\widehat{Z}_{i,:}^{A}\right\|^{2}
$$

#### LGCN
learnable graph conovlutional networks



#### MONET
A propose spatial-domain model (MoNet) on non-Euclidean domains which could generalized serveal previous network. GCNN, ACNN, GCN, DCNN
#### GRAPHSAGE
a general inductive framework.

$$
\begin{aligned}
\mathbf h_{N_{v}}^{t} &=\text { AGGREGATE }_{t}\left(\lbrace \mathbf{h}_{u}^{t-1}, \forall u \in N_{v}\rbrace \right) \cr
\mathbf{h}_{v}^{t} &=\sigma\left(\mathbf{W}^{t} \cdot\left[\mathbf{h}_{v}^{t-1} \| \mathbf{h}_{N_{v}}^{t}\right]\right)
\end{aligned}
$$


1. Mean aggregator
   
$$
\mathbf h_{v}^{t}=\sigma\left(\mathbf{W} \cdot \operatorname{MEAN}\left(\lbrace \mathbf{h}_{v}^{t-1}\rbrace \cup\lbrace \mathbf{h}_{u}^{t-1}, \forall u \in N_{v}\rbrace \right)\right.
$$

2. LSTM aggregator
   
3. Pooling aggregator
   
$$
\mathbf h_{N_{v}}^{t} =\max \left(\lbrace \sigma\left(\mathbf{W}_{\mathrm{pool}} \mathbf{h}_{u}^{t-1}+\mathbf{b}\right), \forall u \in N_{v}\rbrace \right)
$$

an unsupervised loss function which encouage nearby nodes to have similar representations while distant nodes have different representations:

$$
J_{G}\left(\mathbf{z}_{u}\right)=-\log \left(\sigma\left(\mathbf{z}_{u}^{T} \mathbf{z}_{v}\right)\right)-Q \cdot E_{v_{n} \sim P_{n}(v)} \log \left(\sigma\left(-\mathbf{z}_{u}^{T} \mathbf{z}_{v_{n}}\right)\right)
$$

