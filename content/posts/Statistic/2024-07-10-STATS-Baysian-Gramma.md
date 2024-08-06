---
title: 'Bayesian grammar'
date: 2024-07-16
categories: ["Machine Learning"]
tags: ["Bayes"]
comments: true
#markup: "mmark"
math: true
draft: true
---



## Basics



## Bayesian Inference

## Building a Bayesian model for events

Bayes' Rule is given by:
$$
\begin{equation}
P(H | E) = \frac{P(E | H) \cdot P(H)}{P(E)}
\end{equation}
$$

Where:

$P(H | E)$ is the `posterior probability`: the probability of the hypothesis $H$ given the evidence $E$.

$P(E | H)$ is the `likelihood`: the probability of the evidence $E$ given that the hypothesis $H$ is true.

$P(H)$ is the `prior probability`: the initial probability of the hypothesis $H$ before observing the evidence $E$.

$P(E)$ is the `marginal likelihood` or `evidence`: the total probability of the evidence $E$ under all possible hypotheses.




### Prior probability model

As a first step in our Bayesian analysis, we’ll formalize our prior understanding (probability) of a event happens

As a valid probability model must: 
1. it accounts for all possible events (all articles must be fake or real); 
2. it assigns prior probabilities to each event; 
3. these probabilities sum to one.

### Conditional probability & likelihood

In the second step of our Bayesian analysis, we’ll summarize the insights from the new data we collected. 

the `conditional probability `of  A given  B,  $$P(A | B)$$, measures the probability of observing  Ain light of the information that  B occurred. For example, the certainty of an event  
A might increase in light of new data  B.

Conditional probabilities are fundamental to Bayesian analyses.

###  Normalizing constants

To be continue


### Posterior probability model via Bayes’ Rule!

Bayes's Rule

`posterior  = ( prior * likelihood) / normalization_constant`









## Reference

1. [Bayes Rules!](https://www.bayesrulesbook.com/chapter-2#building-a-bayesian-model-for-events)
2. [Zhang Zhenhu](https://www.zhangzhenhu.com/glm/source/%E8%B4%9D%E5%8F%B6%E6%96%AF%E4%BC%B0%E8%AE%A1/content.html)