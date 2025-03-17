---
title: "Graph: Augmentation"
date: 2021-08-12
categories: ["Machine Learning with Graphs"]
comments: true
tags: ["Deep Learning", "Graph"]
math: false
---


Graph Augmentations using PyG


## Graph Structure Agumentation


### Half-Hop

HalfHop adds a “slow node” to all edges with some probability p . Note that these slow nodes have averaged features from the parent nodes, and additionally are undirected.

### Virtual Node

VirtualNode (Gilmer 2017) appends a virtual node to the given homogeneous graph that is connected to all other nodes. This allows information to travel long distances during the propagation phase, especially if such paths are normally sparse.

## Node Feature Transformation

### Random Walk Positional Encoding (RWPE)

this method is based on the random walk diffusion process. 
This method concatenates the extra random-walk features to each node’s features.

### Laplacian Eigenvector Positional Encoding (LapPE)

This takes the first lappe_k eigenvectors of the graph’s laplacian matrix and adds it to the node feature matrix. This is not only unique, but is distance-sensitive with respect to the Euclidean norm. One thing to be careful of is that since lappe_k is the number of eigenvectors to look at, this shouldn’t exceed the number of nodes.

## Training Agumentation

### Dropout Edge

PyG’s dropout_edge function randomly removes edges from the graph

### Mask Features

PyG’s mask_feature function randomly masks parts of node features. 

## References

[CS224W](https://medium.com/stanford-cs224w/graph-augmentations-using-pyg-6e8d1e093450)