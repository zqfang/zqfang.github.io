# Graph: GCN and GAT


Graph Convolutional Network and Graph Attention

## Why deep graph encoder ?

Limitations of Shallow Encoders (e.g. node2vec)

- $O( | V | )$ parameters are needed:
  - No sharing of parameters between nodes
  - Every node has its own unique embedding
- Inherently "transductive":
  - Can not generate embeddings for nodes that are not seen during training
- Do not incorporate node features
  - Many graphs have features that we can and should leverage

## Graph Convolutional Network

Could get embedding for unseen nodes!!!

**Aggreate Neighbors**: Generate node embeddings based on local network neighborhoods.  
1. Intuition: Nodes aggregate information from their neighors using neural networks.
![GCN1](/images/ml/gcn1.png)

1. Computation graph: defined by networkneigborhood
![GCN2](/images/ml/gcn2.png)

1. Layers: Model can be of arbitary depth
   - Nodes have embeddings at each layer
   - Layer-0 embedding:  Node $u$'s input feature $x_u$
   - Layer-K embedding gets information from nodes that are K hops away
  
![GCN3](/images/ml/gcn3.png)

### Training
1. Unsuperised training
   - use only the graph structure
   - "similar: nodes have similar embeddings
   - unspuervise loss function could based on: 
     - Random walks (node2vec, DeepWalk, struc2vec)
     - Graph factorization
     - Node proximity in the graph
2. Supervised training
   - directly train the mode for a supervised task (e.g. node classification) 

### GraphSAGE: Generalized neigborhood aggregation

![GCN4](/images/ml/gcn4.png)

where,
- $W_k$, ${B_k}$ is learnable weighted matrices.
- $h_v^0 = x_v$: initial 0-th layer embeddings are equal to node features.
- $\mathbf{h}_{v}^{k-1}$: previous layer embedding of $v$.
- $\mathbf{z} _v = \mathbf{h} _{v}^{k}$: Embedding after $k$ layers of neigborhood aggregation
- $\sigma$: Non-linearity, e.g., ReLU 

AGG:
- Mean: take a weighted average of neighbors 
  $$\text{AGG} = \sum _{u \in N(v)} \frac{ \mathbf{h} _{u}^{k-1} } { | N(v) | }$$

- Pool: Transform neighbor vectors and apply symmetric vector function
  $$\text{AGG} = \gamma ([\mathbf{Q}\mathbf{h}_{u}^{k-1}, \forall u \in (N(v))])$$  
    - $\gamma$: element-wise mean/max

- LSTM: Apply LSTM to reshuffled of neighbors
  $$\text{AGG} = \text{LSTM} ([\mathbf{h}_{u}^{k-1}, \forall u \in \pi(N(v))])$$

## Graph Attention Network

### Simple Neighorhood Aggregation

The formula

$$
\mathbf{h} _{v}^{k} =  \sigma (\mathbf{W} _{k}  \sum _{u \in N(v)} \frac{ \mathbf{h} _{u}^{k-1}}{ | N(v) | } + \mathbf{B} _{k} \mathbf{h} _{v}^{k-1} )
$$

Equivalently rewritten in vector form:

$$
\mathbf{H} ^{(l+1)} = \sigma ( \mathbf{H} ^{(l)} \mathbf{W} _{0}^{l} + \tilde{\mathbf{A}} \mathbf{H} ^{(l)} \mathbf{W} _{0}^{l})
$$

with $\tilde{A} = D^{- \frac{1}{2}} A D^{- \frac{1}{2}}$.

### Graph convolutional operator

- Aggregates messages across neighborhoods. $N(v)$
- $\alpha _{vu} = \frac{1}{ | N(v) |}$ is the **weighting factor** of node $u$'s message to node $v$
- $\alpha _{vu}$ id defined explicitly based on the structural properties of the graph
- All neighbors $u \in N(v$ are equally important to node $v$)

### Attention strategy

Allows for (implicitly) specifying different importance values ($\alpha_{vu}$) to different neighbors

Compute embedding $\mathbf{h}_{v}^{k}$ of each node in the graph following:
- Nodes attend over their neighorhoods' message
- Implicitly specifying different weights to different nodes in a neighborhood


1. Attention Mechanism
   - Compute attention coefficients $e_{vu}$ across pairs of nodes $u$, $v$ based on their messages:
   $$e_{vu} = a (\mathbf{W}_{k} \mathbf{h}_{u}^{k-1}, \mathbf{W}_{k} \mathbf{h}_{v}^{k-1})$$
   - Normalize coefficients using the softmax function in order to comparable across different neighorhoods:
   
    $$
    \begin{aligned}
    \alpha_{vu} &=\frac{ \exp (e_{v u} )} {\sum_{k \in N(v)} \exp (e_{v k} )} \cr 
    \boldsymbol{h}_{v}^{k} &=\sigma (\sum_{u \in N(v)} \alpha_{v u} \boldsymbol{W}_{k} \boldsymbol{h}_{u}^{k-1} )
    \end{aligned}
    $$

2. Multi-head attention
   - Attention operations in a given layer are independently replicated R times (each replica with different parameters)
   - Outputs are aggregated (by concatenating or adding)



## Reference

[Jure Leskovec, Stanford CS224W: Machine Learning with Graphs](http://cs224w.stanford.edu) 
