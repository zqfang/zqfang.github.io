# Graph: Foundation


Basics.

## Definition

1. Graph: $G(V, E)$
2. Adjacency Matrix: $A$
3. Degree: $D$, the number of nodes that are adjacent to $v$.
4. Neighbors: $N$, the number of $N_{v(i)}$ is equal to $D_{v(i)}$.

### Connectivity
5. Walk
   1. A walk on a graph is an alternating sequence of nodes and edges, starting with a node and ending with a node where each edge is incident with the nodes immediately preceding and following it.
   2. A walk starting at node $u$ and ending at node $v$ is called a $u$-$v$ walk.
   3. The length of a **walk**: the number of **edges** in this walk.
   4. **Trail**: A trail is a walk whose **edges** are distinct
   5. **Path**: A path is a walk whose **nodes** are distinct
6. Subgraph
7. Connected Component: Given a graph $G(V,E)$, a subgraph $G \prime (V \prime, E \prime)$ is said to be a connected component if there is at least one path between any pair of nodes in the graph and the nodes in $V \prime$ are not adjacent to any vertices in $V/V\prime$.
8. Connected Graph
9. Shortest Path: The shortest path between node $v_s$ and node $v_t$ is defined as:

$$
p^{sp}_{st} = \arg \min _{ p \in paths} \vert p \vert
$$


10.  Diameter: the diameter of a graph is defined as the length of the longest shortest path in the graph.

### Centrality

In a graph, the centrality of a node measures the importance of the node in the
graph.

11. Degree Centrality: measure the centrality of a given node based on its degree
12. The eigenvector centrality:  (Bonacich, 1972, 2007) defines the centrality score of a given node $v_i$ by considering the centrality scores of its neighboring nodes as:

$$
\boldsymbol c_{e} (v_i) = \frac{1} {\lambda} \sum^{N}_{j=1} A_{i,j} \cdot \boldsymbol{c}_{e} {v_j}
$$

or reform:


$$
\lambda \boldsymbol{c}_e = \boldsymbol{A} \cdot \boldsymbol{c}_e
$$


Clearly, $\boldsymbol{c}_e$ is an eigenvector of the matrix $\boldsymbol{A}$ with its corresponding eigenvalue $\lambda$.

13. Katz centrality:  It's a variant of the eigenvector centrality, which not only considers the centrality scores of the neighbors but also includes a small constant for the central node itself.

$$
c_{k} (v_i) = \alpha \sum^{N}_{j=1} A_{i,j} c_k (v_j) + \beta
$$

where $\beta$ is a constant. the Katz centrality is reformed as:

$$
\bold{c}_k = \alpha \bold{A} \bold{c}_k + \boldsymbol{\beta}
$$

$$
(\bold I - \alpha \cdot A )\bold{c}_k = \boldsymbol{\beta}
$$

14. Betweenness Centrality: Another way to measure the importance of a node is to check whether it is at an important position in the graph. Specifically, if there are many paths passing through a node, it is at an important position in the graph.

### Spectral Graph Theory

15. Laplacian Matrix: $L = D -A$
16. Normalized Laplacian Matrix: 

$$
\begin{aligned}
L &=D^{-\frac{1}{2}} (D - A) D^{-\frac{1}{2}} \cr
&=I-D^{-\frac{1}{2}} A D^{-\frac{1}{2}}
\end{aligned}
$$



