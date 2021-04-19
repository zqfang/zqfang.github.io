---
title: "Graph: GraphRNN"
date: 2020-12-12
categories: ["Machine Learning with Graphs"]
comments: true
tags: ["Deep Learning", "Graph"]
math: true
---

## Why is it interesting

1. Drug discovery
  - discovery highly drug-like molecules 
  - complete an existing molecule to optimize a desired property

2. Discovering novel structures 
3. Network science

## Why is it hard
1. Large and variable output
2. Non-unique representations
   - $n$-node graph can be represented in $n!$ ways
   - Hard to compute/optimize objective functions
3. Complex dependencies 
   - edge fprmation has long-range dependencies 

## Graph Generative Model

Given: Graphs sampled from $p_{data}(G)$  
Goal:   
  - learn the distribution $p_{model}(G)$
  - sample from $p_{model}(G)$

Setup:
  - Assume we try to learn a generative model from a set of points (i.e., graphs) $\lbrace x_i \rbrace$
    - $p_{data}(x)$ is the data distribution, which is never known to us, but we have sampled $
\boldsymbol{x}_{i} \sim p_{data}(\boldsymbol{x})$.
    - $p_{model}(\boldsymbol{x}; \theta)$ it the model, parametrized by $\theta$, that we use to approximate $p_{data}(x)$.

Auto-regressive models

$p_{model}(\boldsymbol{x}; \theta)$ is used for both density estimation and sampling (from the probability density)
  - Apply chain rule: Joint distritbution is a product of conditional distribution

$$
p_{\text {model}}(\boldsymbol{x} ; \theta)=\prod_{t=1}^{n} p_{\text {model}}\left(x_{t} \mid x_{1}, \ldots, x_{t-1} ; \theta\right)
$$

- $\boldsymbol{x}$ is a vector, $x_t$ it the $t$-th dimension. E.g. $\boldsymbol{x}$ is a sentence, $x_t$ is the $t$-th word.
- For graph generation,$x_t$ will be the $t$-th action (add node, add edge)

## GraphRNN: Generating Graph
Idea: Generating graphs via sequentially adding nodes and edges. 

### Model Graphs as Sequences 
Graphs $G$ with node ordering $\pi$ can be uniquely mapped into a sequence of node and edge additions $S^{\pi}$.

The sequence $S^{\pi}$ has two levels:
- Node-level: add nodes, one at a time
- Edge-level: add edges between existing nodes

![grnn1](/images/ml/grnn1.png)


We transformed graph generation problem into a sequence generation problem.  
Need to model 2 processes
- generate a state for a new node (node-level sequence)
- Generate edges for the new node based on its state (Edge-level sequence)

![grnn2](/images/ml/grnn2.png)

![grnn3](/images/ml/grnn3.png)
### GraphRNN
Relationship between **node-level RNN** and **edge-level RNN**
- Node-level RNN generate the initial state for edge-level RNN
- Edge-level RNN generates edges for the new node, then updates node-level RNN state using generated results

![grnn4](/images/ml/grnn4.png)

### Issue: Tractability

- Any node can connect to any prior node
- Too many step for edge generation
  - need to generate full adjacency matrix
  - complex too-long edge dependencies

Solution: **Tractablity via BFS**

- Breadth-First Search node ordering
![grnn5](/images/ml/grnn5.png)

- Benefits:
  - Reduce possible node orderings
    - From $O(n!)$ to number of distinct BFS orderings
  - Reduce steps for edge generation
    - reducing nuber of previous nodes to look at

![grnn6](/images/ml/grnn6.png)

## Evaluating generated graphs

**Challege**: There is no efficient Graph isomorphism test that can be applied to **any class** of graphs  

**Solution**:
  - Visual similarity
  - Graph statistics similarity 


## Reference

[Jure Leskovec, Stanford CS224W: Machine Learning with Graphs](http://cs224w.stanford.edu)