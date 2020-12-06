---
title: "NLP: Word2Vec"
date: 2020-12-05
categories: ["Nature Language Processing"]
comments: true
tags: ["Deep Learning", "NLP"]
math: true
---

Word2Vec

## CBOW
Continuous Bag of Words Model (CBOW)

When trainning, use `N-gram` language model. That's for a target word, select $m$ (window) words before and after.

Model
![parsing](/images/nlp/cbow_network_arch.png)

1. one-hot encoding get $2m$ vectors:
$$X = (x^{c-m}, \cdots, x^{c-1}, x^{c+1}, \cdots, x^{c+m})$$

2. Embeding Vector $\mathcal{V} \in R^{n \times \mathcal{V}}$, 

$$
\left(v_{(c-m)}=\mathcal{V} x^{(c-m)}, v_{(c-m+1)}=\mathcal{V} x^{(c-m+1)}, \ldots, v_{(c+m)}=\mathcal{V} x^{(c+m)}\right)
$$

3. average

$$
\hat{v}=\frac{v_{c-m}+v_{c-m+1}+\ldots+v_{c-m}}{2 m}
$$

4. multiplut output layer matrix $\mathcal{U} \in R^{n \times \mathcal{V}}$,

$$
z = \mathcal{U} \hat{v}
$$

5. then $\hat{y}$,

$$ \hat{y} = \operatorname{softmax}(z)$$

6. optimization: cross-entropy

$$
\begin{aligned}
\operatorname{minimize} \mathcal{J} &=-\log P (w_{c} \mid w_{c-m}, \cdots, w_{c-1}, w_{c+1}, \cdots, w_{c+m}) \cr 
&=-\log P\left(u_{c} \mid \hat{v}\right) \cr 
&=-\log \frac{\exp \left(u_{c}^{T} \hat{v}\right)}{\sum_{j=1}^{|V|} \exp \left(u_{j}^{T} \hat{v}\right)} \cr 
&=-u_{c}^{T} \hat{v}+\log \sum_{j=1}^{|V|} \exp \left(u_{j}^{T} \hat{v}\right)
\end{aligned}
$$

## Skip-gram

it's on the opposite of CBOW

![parsing](/images/nlp/skip_gram_nn.jpg)

1. generate one-hot encoding for $x$
2. multiply embeding

$$v_c = \mathcal{V}x$$

3. multiply output matrix $\mathcal{U}$, get $2m$ vectors.

$$
u = \mathcal{U}v_c = u_{c-m}, \cdots, u_{c-1}, u_{c+1, \cdots, u_{c+m}}
$$ 

4. for each vector, apply `softmax`, get 

$$
y^{(c-m)}, \cdots, y^{(c-1)}, y^{(c+1)}, \cdots, y^{(c+m)}
$$

5. loss 

$$
\begin{aligned}
\text { minimize } J &=-\log P (w_{c-m}, \ldots, w_{c-1}, w_{c+1}, \ldots, w_{c+m} ) \cr 
&=-\log \prod_{j=0, j \neq m} P\left(w_{c-m+j} \mid w_{c}\right) \cr 
&=-\log \prod_{j=0, j \neq m}^{2 m} P\left(u_{c-m+j} \mid v_{c}\right) \cr 
&=-\log \prod_{j=0, j \neq m} \frac{2 m}{\sum_{k=1}^{|V|} \exp \left(u_{k}^{T} v_{c}\right)} \cr 
&=-\sum_{j=0, j \neq m}^{2 m} u_{c-m+j}^{T} v_{c}+2 m \log \sum_{k=1}^{|V|} \exp \left(u_{k}^{T} v_{c}\right)
\end{aligned}
$$


## Others

### Sub-sampling
use probability $P$ to random delete words, e.g. "the", "a".

$$
P = 1 - \sqrt{\frac{\text{sample}}{\text{freq}(w)}}
$$

### Neagtive sampling

When training, update postive word and partial of negative words.

### Hiearchical Softmax
build `Huffman Tree` acroding to the word freqencies.
the higer freq of words, the higher word levels, then the learning become easier and faster. 


## Reference 

[Word2Vec- The Skip-Gram Model](http://mccormickml.com/2016/04/19/word2vec-tutorial-the-skip-gram-model/)
[word2vec](https://blog.razrlele.com/p/2455)