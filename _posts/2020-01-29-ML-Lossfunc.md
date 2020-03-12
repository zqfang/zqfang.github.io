---
title: 'Loss function for multi-label classification'
date: 2020-01-29
categories: ["Machine Learning"]
comments: true
---

Multi-label classification, tasks commonly be seen on health record data (multi symptoms).

Loss function design:  

1. Multi binary cross-entropy  
each class has a binary output

2. Label smoothing, another regularization technique 

It’s designed to make the model a little bit less certain of it’s decision by changing a little bit its target: instead of wanting to predict 1 for the correct class and 0 for all the others, we ask it to predict 1-ε for the correct class and ε for all the others, with ε a (small) positive number and N the number of classes. This can be written as:

$$
\text {loss}=(1-\varepsilon) c e(i)+\varepsilon \sum c e(j) / N
$$

where ce(x) is cross-entropy of x (i.e. −log(px)), and i is the correct class.



finally, for multi-label loss function:

$$
(1-\epsilon) \sum_{i}\left(-\frac{\log p_{i}}{n}\right)+\frac{\epsilon}{N} \sum\left(-\log p_{i}\right)
$$

See the fastai implementation here:
[LabelSmoothingCrossEntropy](https://github.com/fastai/fastai2/blob/master/fastai2/layers.py#L285)

about line 285:
```python
class LabelSmoothingCrossEntropy(Module):
    y_int = True
    def __init__(self, eps:float=0.1, reduction='mean'): self.eps,self.reduction = eps,reduction

    def forward(self, output, target):
        c = output.size()[-1]
        log_preds = F.log_softmax(output, dim=-1)
        if self.reduction=='sum': loss = -log_preds.sum()
        else:
            loss = -log_preds.sum(dim=-1)
            if self.reduction=='mean':  loss = loss.mean()
        return loss*self.eps/c + (1-self.eps) * F.nll_loss(log_preds, target.long(), reduction=self.reduction)

    def activation(self, out): return F.softmax(out, dim=-1)
    def decodes(self, out):    return out.argmax(dim=-1)
```