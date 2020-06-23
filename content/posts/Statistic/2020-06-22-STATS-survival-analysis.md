---
title: "Survival Analysis"
date: 2020-06-22
categories: ["Statistic"]
comments: true
math: true
---

## Censoring

Censoring

- Surivial without Censoring
- Surivial with Censoring

## Kaplan Meier Curve

More individual in each group, better sepration of the group, better p-value

* Takes censoring into account
* Estimates probabilitu of "survival" on a given day
* Conditional probability of surviving on a given day:

$$
\frac {N_{ \text{"alive" day before}} - N_{ \text{"dying" nextday}}} { \text{"alive" day before}}
$$


Kaplan-Meier survival curve

* Survival times $t_1 \leq t_2 \leq \cdots \leq t_n$
* The proportion of subjects, $S(t)$, surviving beyoind any follow up time $t$ is estimated by (conditional probability):

$$
S(t) = \frac {r_1 - d_1}{r_1} \times \frac {r_2 - d_2}{r_2} \times \cdots \times \frac{r_p - d_p}{r_p}
$$

where

- $t_p$ is the largest survival time less han or equal to $t$
- $r_i$ is the number of subjects alive just before time $t_i$
- $d_i$ = numebr who died at time $t_i$
- for censored obeservations $d_i = 0$

## Statistic

### Log Rank Test

- Compares survival times of two independent groups.
- Assumes that the relative risk of event (e.g. death) between the two groups is constant (proportional hazards)
- Ranks the survial times combined and compared observed and expected rates


**Null hypothesis**: the rates of events (death) in the two groups are equal

under $H_0$, 

$$
X^2 = \frac { (O_A - E_A)^2}{E_A} + \frac { (O_B - E_B)^2}{E_B} \sim \chi^2
$$


- $O_A$: observed events in group A
- $E_A$: expected events in gorup A under null hypohesis
   
expect = (proportion in risk set) * (# of failures over both groups)


$$
e_{1j} = ( \frac{ n_{1j}}{ n_{1j} + n_{2j}}) \times ( m_{1j} + m_{2j})
$$

$$
e_{2j} = ( \frac{ n_{2j}}{ n_{1j} + n_{2j}}) \times ( m_{1j} + m_{2j})
$$




### Cox Regression

- Extends comparison of survial times to allow different predictors (estimate k variable together)
- Models the `hazard`: probability of dying at a point in time, given survival to that point in time

$$
H(t) = H_0(t) \times \exp(b_1X_1 + b_2X_2 + \cdots + b_kX_k )
$$

- Model links to a baseline hazard, $H_0(t)$



* Can accomodate many variables, both discrete and continuous measures of event times
* Proportional hazards assumption: the hazard for any individual is a fixed proportion of the hazard for any other individual

Hazard ratio

- Exp(B) give the hazard ratio (or relative hazard/risk)
