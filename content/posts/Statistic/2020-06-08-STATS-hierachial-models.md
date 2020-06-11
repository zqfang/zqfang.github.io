---
title: 'Multilevel (Hierachical) Models'
date: 2020-06-13
categories: ["Statistic"]
comments: true
math: true
---
What's Multilevel models, and how to deal with it

## What is multilevel model
Multilevel model AKA: 
- multilevel Models
- random-effects models 
- hierarchical models 
- variance-components models
- random-coefficient models
- mixed models


Many kinds of data, including observational data collected in the human and biological sciences, have **a hierarchical or clustered structure, or non-hierarchical structures.**.


## A Simple Example

Given a set of repeated measures data giving growth patterns for a sample of 26 boys in Oxford, England. The height of each boy is measured on nine different occasions.

![height](/images/stats/hm.1.png)

We could try modelling the growth pattern with a simple linear regression

$$
H = \beta_0 + \beta_1 A + \epsilon
$$ (1)

where H and A represent height and age and $\epsilon$ represents the variation in height that
cannot be explained by the linear relationship with age.

However,  If we try to use Model 1 for the complete set of data, the fit will be very poor


To make model more realistic, we allow the intercept in Model 1 to vary form subject to subject.

$$
H_{ij} = \beta_{0j} + \beta_1 A_{ji} + \epsilon_{ij}
$$ (2)

Now, assume that the individual intercepts follow a normal distribution with variance $\tau_{0}$,

$$
\beta_{0j} = \beta_{0} + \mu_{0j}
$$ (3)

where, $\mu_{0j} \text{\textasciitilde} N (0, \tau_0)$

Model 3 accounts for the variation from one subject to another

Fitting the multilevel model to the data, and obtain much better predictions
![height](/images/stats/hm.2.png)



## How to deal with hierachical structures

In simple linear regression, where we try to fit our data to a straight line, the effects of a cluster can be addressed by allowing multiple levels of random effects, or residuals specific to the clusters.

We can also find the correlation coefficient between members of clusters, and that can be used to set up variables that show the ‘fixed effects’ that are taking place between related data points. These correlations will be represented by coefficients in a modeling equation.


1: Hierarchical structures : model all levels simultaneously
2: Non- Hierarchical structures
a) cross-classified structure
b) multiple membership with weights



## How do multilevel models differ from regression models?

$$
H_{ij} = (\beta_0 + \mu_{0j}) + \beta_1 A_{ij} + \epsilon_{ij} = \beta_0 + \beta_1 A_{ij} + \mu_{0j} + \epsilon_{ij}
$$


The feature that distinguishes this model from an ordinary regression model is the presence
of two random variables 
- the measurement level random variable $\epsilon_{ij}$
- the subject level 
random variable $\mu_{0j}$


Because multilevel models contain a mix of fixed effects and random effects, they are sometimes known as mixed-effects models.


## Benefits of multilevel modelling

- Generalize to a wider population
- Fewer parameters are needed
- Information can be shared between groups


参考：

1. https://en.wikipedia.org/wiki/Multilevel_model
2. http://www.statstutor.ac.uk/resources/uploaded/multilevelmodelling.pdf