# Annova


## one-way ANOVA from scratch

1. Calculate the Sum of Squares Total (SST):

$$
SS_{total} = \sum_{j=1}^k \sum_{i=1}^l (X_{ij} - \bar{X})^2
$$

2. Calculate the Sum of Squares Within Groups (SSW): 

$$
SS_{within} = \sum_{j=1}^k \sum_{i=1}^l (X_{ij} - \bar{X_j})^2
$$


3. Calculate the Sum of Squares Between Groups (SSB):

$$
SS_{between} = \sum_{j=1}^k n_j ( \bar X_{j} - \bar{X}) ^2
$$

$n_j$: numbers of individual point in group j.

Verify that  

$$
SS_{total} = SS_{between} + SS_{within}
$$

4. Calculate the Degrees of Freedom (df)  

Calculate the Degrees of Freedom Total (DFT)

$$
df_{total} = n -1
$$


Calculate the Degrees Between k Groups (DFB)

$$
df_{bewteen} = k -1
$$

Calculate the Degrees of Freedom Within Groups (DFW)

$$
df_{within} = n - k
$$



Verify that 
$$ 
df_t = df_w + df_b
$$

1. Calculate the Mean Squares


Calculate the Mean Squares Between (MSB)

$$
MS_{between} = \frac{SS_{between}}{df_{between}}
$$


Calculate the Mean Squares Within (MSW)

$$
MS_{within} = \frac{SS_{within}}{df_{within}}
$$

6. Calculate the F Statistic

$$
F = \frac{ MS_{between}}{MS_{within}}
$$

7. get pvalue

```python
import scipy.stats as stat
pvalue = stat.f.sf(F, dfb, dfw) # sf: pvalue = 1 - stat.f.cdf() 
```

## ANOVA Effect size
Omega squared (ω2) is a measure of effect size, or the degree of association for a population. It is an estimate of how much variance in the response variables are accounted for by the explanatory variables. Omega squared is widely viewed as a lesser biased alternative to eta-squared, especially when sample sizes are small.


MSerror: mean square error SSE/df(error)


Formula
$$
\omega^2  = \frac {SS_{Effect} - df_{Effect} MS_{error}}{SS_{total} + MS_{error}}
$$

for multi-factor, completely randomized design,

Formula
$$
\omega^2  = \frac {SS_{Effect} - df_{Effect} MS_{errors}} {SS_{Effect} + (N-df_{Effect}) MS_{error}}
$$

Interpreting Results  

* ω2 can have values between ± 1.
* Zero indicates no effect.
* If the observed F is less than one, ω2 will be negative.

## ANOVA Post-hoc comparison

ANOVA does not tell which group are significantly different from each other. To know the pairs of significant different groups, we could perform multiple pairwise comparison (Post-hoc comparison) analysis using Tukey HSD test.

```python
from statsmodels.stats.multicomp import pairwise_tukeyhsd
# perform multiple pairwise comparison (Tukey HSD)
m_comp = pairwise_tukeyhsd(endog=d_melt['value'], groups=d_melt['treatments'], alpha=0.05)
```

## Two-way (two factor) ANOVA

example

```python
# load packages
import statsmodels.api as sm
from statsmodels.formula.api import ols
# Ordinary Least Squares (OLS) model
# C(Genotype):C(years) represent interaction term
model = ols('value ~ C(Genotype) + C(years) + C(Genotype):C(years)', data=d_melt).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
anova_table
```

see more about ANOVA in python [here](https://reneshbedre.github.io/blog/anova.html)



$$
H_0: \beta_1 = 0 \\
H_A: \beta_1 \neq 0
$$

$$
H_0: \frac {p_{new}} {p_{old}} =  \frac {p_{new}} {1 - p_{new}} = 1 \\
H_A:  \frac {p_{new}} {1 - p_{new}} \neq 1
$$


