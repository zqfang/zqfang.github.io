
---
title: "SwiGLU"
date: 2024-07-29
categories: ["Nature Language Processing"]
comments: true
tags: ["Deep Learning", "LLM", "NLP"]
math: true
---


The most commonly used activation function in LLM.

![SwiGLU](/images/nlp/swish.jpg)



```python
class SwiGLU(nn.Module):
    def __init__(self, w1, w2, w3):
        super.__init__()
        self.w1 = w1
        self.w2 = w2
        slef.w3 = w3
    def forward(self, x):
        x1 = F.linear(x, self.w1.weight)
        x2 = F.linear(x, self.w2.weight)
        hidden = F.silu(x1)* x2
        return F.linear(hidden, self.w3.weight)

```