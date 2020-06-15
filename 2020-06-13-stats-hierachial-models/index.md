# Multilevel (Hierachical) Models

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
$$ 

where 
- H: height
- A: age 
- $\epsilon$: the variation in height that
cannot be explained by the linear relationship with age.

However,  If we try to use the model above for the complete set of data, the fit will be very poor (see figure above)


To make model more realistic, we allow the intercept in the model above to vary from subject to subject. 

New multilevel model:

$$
H_{ij} = \beta_{0j} + \beta_1 A_{ji} + \epsilon_{ij}
$$

Now, assume that the individual intercepts follow a normal distribution with variance $\tau_{0}$,

$$
\beta_{0j} = \beta_{0} + \mu_{0j}
$$ 

where, $\mu_{0j} \sim \mathcal{N} (0, \tau_0)$, $\mu_{0j}$ accounts for the variation from one subject to another

Fitting the multilevel model to the data, and obtain much better predictions
![height](/images/stats/hm.2.png)



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


## How to deal with hierachical structures

1. Hierarchical structures : model all levels simultaneously  
2. Non- Hierarchical structures
  -  cross-classified structure
  -  multiple membership with weights

...

参考：

1. https://en.wikipedia.org/wiki/Multilevel_model
2. http://www.statstutor.ac.uk/resources/uploaded/multilevelmodelling.pdf
3. https://www.cs.princeton.edu/courses/archive/fall11/cos597C/lectures/hierarchical-models.pdf
