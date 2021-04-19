---
title: "Graph: Semi-supervised Node Classification"
date: 2020-12-10
categories: ["Machine Learning with Graphs"]
comments: true
tags: ["Deep Learning", "Graph"]
math: true
---


Problems: Given a network with labels on some nodes, how do we assign labels to all other nodes in the network?

classification label of an object $O$ in network may depend on:
- Features of $O$
- Labels of the objects in $O$'s neighborhood
- Features of objects in $O$'s neigborhood

## Collective classification models

- Reational clasifiers
- Iterative classifications
- Loopy belief propagation
   
### Intuition 

Simultaneous classification of interlinked nodes using correlations

### Applications

- Document classification
- Part of speech tagging
- Link prediction
- Optical character recognition
- Entity resolution in sensor networks
- Image/3D data segmentation
- Spam and fraud detection

### Markov Assumption
The label ${Y_i}$ of one onde $i$ depends on the labels of its neighbors ${N_i}$.

$$P(Y_i \vert i ) = P (Y_i \vert N_i)$$

### Steps

1. Local Classifier: Assign initial labels
   - predicts label based on node attributes/  features
   - standard classification task
   - Does not use network information 
2. Retional Classifier: Capture correlations based on the network 
   - learns a classifer to label one node based on the labels and /or attributes of tis neighbors
   - This is where network information is used
3. Collective Inference: Propagate correlations through networks 
   - Apply relational classifier to each node iteratively
   - Iterate until the inconsistency between neighboring labels is minimized 
   - Network structure substantially affects the final prediction 


## Relational classifers
1. Basic idea:  
Class probability of ${Y_i}$ is a weighted average of class probabilities of its neigbors.

2. Iteration
   - labeled nodes: initialize with ground-truth Y labels
   - unlabeled nodes: initialize Y uniformly 
   - Update all nodes in a random order until convergence or until maximum number of iteration is reached 

3. Repeat for eahc node $i$ and label $c$

$$
P( Y_i = c) = \frac{1}{\sum_{(i,j) \in E}W(i,j)}\sum_{(i,j) \in E}W(i,j)P(Y_j=c)
$$

   - $W(i,j)$ is the edge strength from $i$ to $j$

However,
   - Model cannot use node feature information
   - Convergence is not guaranteed 

## Iterative classification

### Main idea
Classify node $i$ based on its attributes as well as labels of neigbor set $N_i$.

  - Create a flat vector $a_i$ for each node $i$
  - Train a classifier to classify using $a_i$
  - **Aggregate** various numbers of neighbors using: count, mode, proportion, mean, exists, etc.

### Basic architecture of iterative classifers

1. Bootstrap phase
    - Convert each node $i$ to a flat vector $a_i$
    - Use local classifier $f(a_i)$ to compute best value for $Y_i$
2. Iteration phase: iterate till convergence
    - Repeat for each node $i$
      - Update node vector $a_i$
      - Update label $Y_i$ to $f(a_i)$. This is a hard assignment
    - Iterate until class labels stabilized or max number of iterations is reached

However,  
  Convergence is not guaranteed. Run for max number of iterations

### Applications

fake reviewer/review detction

## Belief propagation

**Belief propagation** is a dynamic programming approach to answering conditional probability queries in a graphical model

Notation
- Label-label potential matrix $\psi$: Dependency between a node and its neigbor
- $\phi(Y_i, Y_j)$ equals the probability of a node $j$ being in state $Y_j$ given that it has a $i$ neigbbor in state $Y_i$
- Prior belief $\phi$: Probability $\phi _i (Y_i)$ of node $i$ being in state ${Y_i}$
- $m_{i \rightarrow j}(Y_j)$ is $i$'s estimate of $j$ being in state $Y_j$ 
- $\mathcal{L}$ is the set of all states

![Loopy-BP1](/images/ml/loopy-BP1.png)

![Loopy-BP1](/images/ml/loopy-BP1.png)

## Reference

[Jure Leskovec, Stanford CS224W: Machine Learning with Graphs](http://cs224w.stanford.edu)