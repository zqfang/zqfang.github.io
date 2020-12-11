# Graph: Node2Vec


Node Embedings are learnt in the same way as `word2vec` (skip-gram model)

However, graphs could be (un)directed, (un)weighted, (a)cyclic and are basically much more complex than the strucure of a sequence...

So how do we generate "corpus" from a graph ?

## Random walk on the graph

Given a graph and a starting point, we **select a neighbor of it at random**; then we select a neigbor of this point at random, and move to it, etc.

**The (random) sequence of points selected** this way is a random walk on the graph

### Why Random walks

1. Expressivity: Flexible stochastic definition of node similarity that incorporates both **local** and **higher-order** neighborhood information
2. Efficiency: Do not need to consider all node pairs when training; only need to consider pairs that co-occur on random walks

## Sampling strategy

![node2vec2](/images/ml/node2vec2.png)

Node2vec's sampling strategy accepts 4 argument:
- **Number of walks**: number of random walks to be generated from each node in the graph
- **Walk length**: how many nodes are in each random walk
- **P**: return hyperparameter
  - return back to the previous node
- **Q**: Inout hyperparameter 
  - Moving outwards : (DFS biased or BFS baised control)
  - intuitively, **q** is the "ratio" of BFS vs. DFS

Also, the standard skip-gram parameters
- context window size
- number of iterations
- etc.

### Node2Vec: Biased Walks

![biasedwalk](/images/ml/randwalk.png)


## Principal

Idea: use flexible, biased random walks that can trade off between local and global views of the network.

![node2vec1](/images/ml/node2vec1.png)

Consider you are on the random walk, and have just transitioned from node $t$ to node $v$ in the above diagram

the probability to transition from $v$ to any neighbors is edge $\alpha$, where $\alpha$ is depened on the hyperparameters. 

- $\bold{P}$: controls the probability to go back to $t$ after visiting $v$
- $\bold{Q}$: controls the probability to go explore undiscovered parts of the graphs.

So, the final travel probability is a function of 

$$
\alpha_{p q}(t, x)=\begin{cases} 
\frac{1}{p} & \text { if } d_{t x}=0 \cr 
1 & \text { if } d_{t x}=1 \cr 
\frac{1}{q} & \text { if } d_{t x}=2
\end{cases}
$$

where,
- $d_{tx}$: the length of shortest path of $t$ and $v$.
  - 0: $x$ is $t$
  - 1: $x$ is neighbor to $t$
  - 2: $x$ and $t$ not connected  

Using the sampling strategy, node2vec will generate "senences" (the directed subgraphs) which are will be used for embedding just like text sentences are used in word2vec.

### Alias Sampling
Preproccsing of transition probabilities for guiding random walk

 Alias Sampling: sampling of nodes while simulating the random walk can be done efficiently in $O(1)$.

[Understand Alias Method](http://shomy.top/2017/05/09/alias-method-sampling/)

[What's really going on ?](https://zhuanlan.zhihu.com/p/56136631)


## Applications

- Clustering/community detection
- Node classification
- Link prediction: predict edge $(i,j)$ based on $f(z_i,z_j)$
  - where we can:
    - concatenate: $f(z_i,z_j) = g([z_i, z_j])$ 
    - Hadamard: $f(z_i,z_j) = g( z_i \star z_j)$ (per coordinate product)
    - Sum/Avg: $f(z_i,z_j) = g(z_i + z_j)$
    - Distance: $f(z_i,z_j) = g( || z_i - z_j || _2$

## Reference

[Jure Leskovec, Stanford CS224W: Machine Learning with Graphs](http://cs224w.stanford.edu)  
[node2vec: Scalable Feature Learning for Networks](https://cs.stanford.edu/~jure/pubs/node2vec-kdd16.pdf)  
[Embeddings for Graph Data](https://towardsdatascience.com/node2vec-embeddings-for-graph-data-32a866340fef)
