# Quantile normalization


Quantile normalization is frequently used in microarray data analysis. It was introduced as `quantile standardization` and then renamed as `quantile normalization`.


## Assumptions of quantile normalization
- The roughly same distribution of values across samples
- Most genes are not differentially expressed


## How Q-normalization work

row: genes  
column: samples/Arrays  

Procedure:
1. order values within each sample
2. determine a rank from lowest to highest and record the order within each sample
3. Average across rows and substitute value with average
4. re-order averaged values in the original order recorded in 2.

Tied rank entries ?  
Average the tied rank entries' mean values and substitute.

![Q-normalization](/images/stats/quantile-normalization.png)




## When NOT to normalize
Consider a dilution experiment. In which distributions are supposed to decrease (left plot), Q-normalization does the totally wrong thing (right plot). When you expect a real difference in distributions, Q-normalization will create weird artifacts.

![Q-normalization](/images/stats/quantile-normalization-2.png)




## StatQuest: Quantile Normalization
{{< youtube ecjN6Xpv6SE >}}


## reference

1. https://en.wikipedia.org/wiki/Quantile_normalization
2. BIOMEDIN 245: Statistical and Machine Learning Methods for Genomics, Stanford 
