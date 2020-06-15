# Quantile normalization


Quantile normalization is frequently used in microarray data analysis. It was introduced as `quantile standardization` and then renamed as `quantile normalization`.


## Quantile, quartile, percentile ???

`Quantiles` are just the lines that divide data into equally sized groups.  
`percentiles` are just quantiles that divide the data into 100 equally sized groups

Example:

0 quartile = 0 quantile = 0th percentile

1 quartile = 0.25 quantile = 25th percentile

2 quartile = .5 quantile = 50th percentile (median)

3 quartile = .75 quantile = 75th percentile

4 quartile = 1 quantile = 100th percentile


## Quantile normalization
Quantile normalization transform the statistical distributions across samples to be the same.
### Assumptions
- The roughly same distribution of values across samples
- Most genes are not differentially expressed

Assume global differences in the distribution are induced by only technical variation!

### How Q-normalization work

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




### When NOT to normalize
Consider a dilution experiment. In which distributions are supposed to decrease (left plot), Q-normalization does the totally wrong thing (right plot). When you expect a real difference in distributions, Q-normalization will create weird artifacts.

![Q-normalization](/images/stats/quantile-normalization-2.png)




## Smooth quantile ormalizaiton 

### Assumptions
the statistical distribution of each sample should be the same ( or have the same distributional shape) `within biological groups or conditions`, but `allowing` that they may `differ between groups`

### How to 

At each quantile, a weight is computed comparing the variability between groups relative to the total variability between and within groups



Let gene(g) denote the $g^{th}$ row after sorting each column in the data. For each row, gene(g), we compute the weight $w(g)$ ∈ [0,1], where a weight of 0 implies quantile normalization within groups is applied and a weight of 1 indicates quantile normalization is applied. The weight at each row depends on the between group sum of squares SSB(g) and total sum of squares SST(g), as follows:

$$
w_{(g)} = \operatorname{median} \bigg\lbrace 1- \frac{SSB_{(i)}}{SST_{(i)}} \bigg\rbrace \text{for } j = g -k, \cdots, g, \cdots, g+k
$$

where $k$ = floor(Total number of genes * 0.05). The number 0.05 is a flexible parameter that can be altered to change the window of the number of genes considered. 

![Q-normalization](/images/stats/smoothqn.jpeg)


## StatQuest: Quantile Normalization
{{< youtube ecjN6Xpv6SE >}}


## reference

1. https://en.wikipedia.org/wiki/Quantile_normalization
2. BIOMEDIN 245: Statistical and Machine Learning Methods for Genomics, Stanford 
3. Hicks SC, Okrah K, Paulson JN, Quackenbush J, Irizarry RA, Bravo HC. Smooth quantile normalization. Biostatistics. 2018;19(2):185‐198. doi:10.1093/biostatistics/kxx028
