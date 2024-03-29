---
title: "Graph: Train, valid, and test dataset split for link prediction"
date: 2021-08-12
categories: ["Machine Learning with Graphs"]
comments: true
tags: ["Deep Learning", "Graph"]
math: true
---

## Link Prediction

- Link prediction is a common task in knowledgegraph's link completeion. 
- Link prediction is usually an **unsupervised or self-supervised task**, which means that sometimes we need to split the dataset and create corresponding labels on our own.

## How to prepare train, valid, test datasets ?  

For link prediction, we will split edges twice 

- Step 1: Assign 2 types of edges in the original graph
    - Message edges: Used for GNN message passing 
    - Supervision edges: Use for computing objectives 
-  After step 1:
    - Only message edges will remain in the graph 
    - Supervision edges are used as supervision for edge predictions made by the model, will not be fed into GNN!
- Step 2: Split edges into train / validation / test

### Option 1: Inductive setting

- training / validation / test sets are on **different graphs**
- The dataset consists of multiple graphs
- Each split can only observe the graph(s) within the split. A successful model should generalize to unseen graphs
- Applicable to node / edge / graph tasks

![inductive](/images/ml/data.split.edge.inductive.png)

### Option 2: Transductive 

- training / validation / test sets are on the **same graph**
- The dataset consists of one graph
- The entire graph can be observed in all dataset splits, we only split the labels
- Only applicable to node / edge prediction tasks

![inductive](/images/ml/data.split.edge.transductive.png)



### Code
Option 1: PyG's `RandomLinkSplit`

```python
from torch_geometric.transforms import RandomLinkSplit, RandomNodeSplit

## designed for transductive learning
tfs = RandomLinkSplit(is_undirected=True, 
                      add_negative_train_samples=True,
                      neg_sampling_ratio=1.0,
                      key = "edge_label", # supervision label
                      disjoint_train_ratio=0,# disjoint mode if > 0
                      # edge_types=None, # for heteroData
                      # rev_edge_types=None, # for heteroData
                      )
train_data, val_data, test_data = tfs(data)
# Here, *_data.edge_index denotes the graph structure used for message passing,
# *_data.edge_label_index and *_data.edge_label denote the training/evaluation edges and their corresponding labels. 

## if inductive learning, need subgraph. e.g 
from torch_geometric.utils import subgraph
train_mask = torch.rand(data.num_nodes) < 0.5
test_mask = ~train_mask

train_data = copy.copy(data)
train_data.edge_index, _ = subgraph(train_mask, data.edge_index, relabel_nodes=True)
train_data.x = data.x[train_mask]

test_data = copy.copy(data)
test_data.edge_index, _ = subgraph(test_mask, data.edge_index, relabel_nodes=True)
test_data.x = data.x[test_mask]
```

Option 2: `deepsnap`'s `GraphDataset`

- The `GraphDataset` is compatible with `Pytorch Geometric` ! 
```python
from deepsnap.dataset import GraphDataset
from deepsnap.hetero_graph import HeteroGraph

hetero = HeteroGraph(H)
dataset = GraphDataset([hetero], 
                        task='link_pred', 
                        edge_train_mode="disjoint", 
                        edge_message_ratio=0.8, 
                        edge_negative_sampling_ratio=2)
dataset_train, dataset_val, dataset_test = dataset.split(transductive=True, 
                                                         split_ratio=[0.8, 0.1, 0.1])
# dataset could be use for PyG or deepsnap's high-level API
```

Check the all docs [here](https://snap.stanford.edu/deepsnap/notes/colab.html)

The content blew is almost the same as in `colab notebooks`.  It's just for easy and quick viewing in any devices.


### General rules
In general, edges in the graph will be splitted to two types: 
1. `message passing` edges: used for GNN message passing 
2. `supervision` edges: used in loss function for backpropagation.
    - Need to include `negative sampling edges`, the edges **not existed** in the original graph.


DeepSNAP's `GraphDataset` will automatically generate labels for all edges. 
- Negative edges: label 0. 
- Positive supervision edges: usually label 1. 
    - If the original edges already have labels (started from 0), all the labels will be added by 1. 

In addition to edges split and negative edge sampling, edges in each of the train, validation and test sets usually need to be **disjoint**.


## Transductive Link Prediction Split

`DeepSNAP` link prediction contains two main split modes (edge_train_mode: all, disjoin)

### Split Mode: All

The figure blew shows the supervision edges in train (blue), validation (red) and test (green) sets. Notice that all original edges in `all` mode will be included in the supervision edges.



**To be more specific**:

* At `training` time: the training supervision edges are **same** with the training message passing edges.
    - The $\text{training supervision edges} == \text{training message passing edges}$
* At `validation` time: the message passing edges are the training message passing edges and training supervision edges (still the training message passing edges in this case).
    -  The $\text{validation supervision edges} \notin \text{training supervision edges}$:  `disjoint`  with `training supervision edges`: 
    -  The $\text{validation message passing edges} = \text{training message passing edges} + \text{training supervision edges}$
* At `test` time: the message passing edges are the union of training message passing edges, training supervision edges, and validation supervision edges. 
    - The $\text{test supervision edges} \notin \lbrace \text{training supervision edges},  \text{valid supervision edges} \rbrace$:  `disjoint` with `training supervision` edges and `validation supervision` edges.
    - The $\text{test message passing edges} = \text{validation supervison edges} + \text{training message passing edges} + \text{training supervision edges}$


### Split Mode: Disjoint

The figure blow shows the supervision edges in train (blue), validation (red), test (green) sets and the training message passing edges (grey). Notice that not all original edges in `disjoint` mode will be included in the supervision edges.


**To be more specific**:

* At `training` time: the training supervision edges are disjoint with the training message passing edges.
    - The $\text{training supervision edges} \notin \text{training message passing edges}$
* At `validation` time: the message passing edges are the union of training message passing edges and training supervision edges. Notice that the validation supervision edges are disjoint with training supervision edges.
    -  The $\text{validation message passing edges} = \text{training message passing edges} + \text{training supervision edges}$
    -  The $\text{validation supervision edges} \notin \text{training supervision edges}$  
* At `test` time: the message passing edges are the training message passing edges, training supervision edges, and validation supervision edges. The test supervision edges are disjoint with training supervision edges and validation supervision edges.
    -  The $\text{test message passing edges} = \text{validation supervison edges} + \text{training message passing edges} + \text{training supervision edges}$
    -  The $\text{validation supervison edges} \notin \lbrace \text{training supervision edges},  \text{valid supervision edges} \rbrace$


## Inductive Link Prediction Split

For inductive link prediction in DeepSNAP, graphs will be splitted to different (train, validation and test) sets. Each graph in the same set will have message passing edges and supervision edges (which are same in this case). But supervision and message passing edges in each graph in different sets are disjoint.


## Negative Sampling Ratio and Resampling

For `link_pred` task, DeepSNAP will automatically and randomly sample negative edges when:
* The dataset is splitted to several datasets, such that one dataset is splitted to train, validation and test.
* The `Batch` of the graph is called or used (this will resample all negative edges).

The number or ratio of negative edges can be controlled by specifying the `edge_negative_sampling_ratio`, which has the default value 1. The resampling can be disabled by setting `resample_negatives` to be False. The example below shows how to set different number or ratio of negative edges.

Training set negative edges will be resampled when the Batch object is called
However, to reduce the computation cost, the validation and test sets negative edges will not be resampled.
```python
dataset = GraphDataset([dg], task=task)
dataset_train, dataset_val, dataset_test = dataset.split(
            transductive=True, split_ratio=[0.8, 0.1, 0.1])
dataloaders = {
    "train": DataLoader(
        dataset_train, collate_fn=Batch.collate([]), shuffle=True),
    "val": DataLoader(
        dataset_val, collate_fn=Batch.collate([]), shuffle=True),
    "test": DataLoader(
        dataset_test, collate_fn=Batch.collate([]), shuffle=True),
}
```

## Message Passing Ratio

Here is an example of adjusting the number of message passing edges and supervision edges in `disjoint` mode. We can control the number of edges by adjusting the `edge_message_ratio`, which defines the ratio between message-passing edges and supervision edges in the training set.


## Node Split
See also dataset split for node classification


![NodeSplit](/images/ml/data.split.node.png)


## Reference

1. [DeepSnap](https://snap.stanford.edu/deepsnap/modules/dataset.html#deepsnap-graphdataset)
2. [Jure Leskovec, Stanford CS224W: Machine Learning with Graphs](http://cs224w.stanford.edu) 