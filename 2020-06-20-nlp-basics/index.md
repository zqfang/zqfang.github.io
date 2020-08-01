# NLP: A short guide for biologist


NLP Basics for the newbies like me

## Languwage model

Models that **assigns probabilities to sequences of words** are called languwage models.

## Count-based Representation

### 1. one-hot representation
### 2. BoW: Bag of words
Blow describes the `occurrence of words` within a document. including

- A Vocabulary of known words
- A measure of the presence of known words, e.g. count

### 3. TF or TF-IDF representation: Term Frequency Inverse Document Frequency

 - TF: the sum of the one-hot representation of a phrase, sentence or document's constituent words

 $$
 TF (w) = \frac { \text{ Number of the term w appears in the document }} { \text{Number of terms in the document}}
 $$

 - IDF: penalizes common tokens and rewards rare tokens

 $$
 IDF(w) = \log \frac{\text{Number of documents}}{\text{Number of documents with term w}}
 $$

 - TF-IDF: $TF(w) \times IDF(w)$

### 4. Positive pointwise Mutual Information (PPMI)

An alternative weighting function to tf-idf

PPMI is a measure of how often two events x and y occur, compared with what we would expect if they were independent:

$$
I(x,y) = \log_2 \frac{P(x, y)}{{P(x)}{P(y)}}
$$

The pointwise mutual information between a target word w and a context word c is then defined as:

$$
PMI(w,c) = \log_2 \frac{P(w, c)}{{P(w)}{P(c)}}
$$

PMI is a useful tool whenever we need to find words that are strongly associated.

`PPMI` replaces all negative PMI values with zeros:

$$
\operatorname{PPMI}(w,c) = \max (\log_2 \frac{P(w, c)}{{P(w)}{P(c)}}, 0 )
$$

However, **PMI has the problem**: very rare words tend to have very high PMI values. So, use a different function $P_{\alpha}(c)$ that raise the probability of the context word to the power of $\alpha$:

$$
\operatorname{PPMI_{\alpha}}(w,c) = \max (\log_2 \frac{P(w, c)}{{P(w)}{P_{\alpha}(c)}}, 0 )
$$

## Word Embedding

### word2vec

The intuition here is taht we could just use running text as implicitly supervised training data for such a classifer: a word $s$ that occurs near the target word $apricot$ acts as gold 'correct answer' to the question "Is word $w$ likely to show up near $apricot$?"

This advoids the need for any sort of hand-labeled superivsion signal.

Alogrithm:

- `Skip-gram with negative sampling`, aslo called SGNS
- Skip-gram trains a probabilistic classifier that given a test targe word $t$ and its context window of $k$ workds $c_{1:k}$, assigns a probability based on how similar this contxt window is to the target word.

$$
\begin{aligned}
P\left(+\mid t, c_{1: k}\right) &=\prod_{i=1}^{k} \frac{1}{1+e^{-t \cdot c_{i}}} \cr
\log P\left(+\mid t, c_{1: k}\right) &=\sum_{i=1}^{k} \log \frac{1}{1+e^{-t \cdot c_{i}}}
\end{aligned}
$$

- Note: similarity between embeddings (Cosine)

$$
Similarity(t,c) \approx t \cdot c
$$
- Skip-gram makes the strong assumption that all contxt words are independent


The intuition of skip-gram is:

1. Treat the target word and a neighboring context word as positive examples;
2. Randomly sample other words in the lexicon to get negative samles;
3. Use logistic regression to train a classifer to distinguid those two case;
4. use the regression weights as the embeddings.


### GloVe

## Concepts

### Corpora, Tokens, and Types

1. corpus (plural: corpora): a text dataset
2. tokens (English): words and numeric sequences separated by white-spaces characters or punctuation
3. instance or data point: the text along with its metatdata
4. dataset: a collection of instances
5. types: unique tokens present in a corpus.
6. vocabulary or lexicon: the set of all types in a corpus

![corpus](/images/nlp/corpus.png)

### Unigram, Bigrams, Trigrams, ... , N-grams

N-grams are fixed-length consecutive token sequence occurring in the text

### Lemmas and Stems

`Lemmas` are the root forms of words.  e.g. the root form of the word fly, can be inflected into other words -- flow, flew, flies, flown, flowing ...

Lemmas, also called `citation form`.  

Stemming: use handcrafted rules to strip endings of words to reduce them to a common form called `stems`

### Word Senses and Semantics

Senses: the different meanings of a word

### Categorizing Words: POS Tagging

Part-of-speech (POS): also known as word classes, or syntactic categories

POS divided into two broad supercateogries: 

- `closed class` types: `function words` *like of, it, and*, or *you*
- `open class` types: nouns, verbs, adjectives, adverbs

POS tagging: labeling individual words or tokens

Common algorithms to do tagging

- HMM: Hidden Markov Models
- MEMM: maximum Entropy Markov Models

### Categorizing Spans: Chunking and Named Entity Recognition

a span of text: a contiguous multi-token boundary.

`chunking` or `shallow parsing`: identify the noun phrases (NP) and verb pharses (VP) in a span of text.

A `named entity` is a string mention of a real-world concept like a person, location, organization, drug name, et. al.

### Coreference

The task of deciding whether two strings refer to same entity

### Minimumn Edit distance

A way to quantify string similarity.

### Perplexity

The **perplexity** (sometimes called PP) of a language model on a test set is the inverse probability of the test set, normalised by the number of words.

$$\begin{aligned} 
\mathrm{PP}(W) &=P(w_{1} w_{2} \ldots w_{N})^{-\frac{1}{N}} \cr &=\sqrt[N]{\frac{1}{P(w_{1} w_{2} \ldots w_{N})}} \cr
&= \sqrt[N]{\prod_{i=1}^{N} \frac{1}{P(w_{i} \mid w_{1} \ldots w_{i-1})}} 
\end{aligned}
$$

Another way to hink about perplexity: as **the weighted average branching factor** of a language.


### Entropy

Entropy is a measure of information. The entropy of the random variable X is:

$$
H(X)=-\sum_{x \in \chi} p(x) \log _{2} p(x)
$$

the log can be computed in any base. If we use log base 2, the resulting value of entropy will mesured in **bits**.

One intuitive way to think about entorpy is as a **lower bound** on the number of bits it would take to encode a certain desision or piece of information in the optimal coding scheme.

### Cross-entropy
The cross-entropy is useful when we don't know the actual probability distribution p that generated some data.

The cross-entropy of m (a model of p) on p is defined by

$$
H(p,m) = \lim_{n \rightarrow \infty} - \frac{1}{n}\sum_{W \in L} p (w_1,\cdots,w_n) \log m (w_1, cdots, w_n)
$$


the cross-entropy $H(p,m)$ is an upper bound on the entropy $H(p)$. For any model m:

$$
H(p) \leq H(p,m)
$$

The more accurate m is, the closer the cross-entropy $H(p,m)$ will be to the true entropy $H(p)$.


### The relation of perplexity and cross-entropy

The approximation to the cross-entropy of a model $M = P(w_i | W_{i-N+1} \cdots W_{i-1})$ on a sequence of words W is 

$$
H(W) = - \frac{1}{N} \log P(w_1 w_2 \cdots w_N)
$$

The perplexity of a model P on a seqence of words W is defined as exp of this cross-entropy

$$
\begin{aligned}
\operatorname{Perplexity}(W) &=2^{H(W)} \cr
&=P\left(w_{1} w_{2} \ldots w_{N}\right)^{-\frac{1}{N}} \cr
&=\sqrt[N]{\frac{1}{P\left(w_{1} w_{2} \ldots w_{N}\right)}} \cr
&=\sqrt[N]{\prod_{i=1}^{N} \frac{1}{P\left(w_{i} \mid w_{1} \ldots w_{i-1}\right)}}
\end{aligned}
$$

## statistical testing

The approach to computing p-values(x) in NLP is to use non-parametric tests. e.g.

- bootstrap test
- approximate randomization

**bootstrapping** refers to repeated drawing large numbers of smaller samples with replacement from an orignial larger sample.

- the intuition of the bootstrap test is that we can create many virtual test sets from an observed test set by repeated sampling from it.
- the method only maks the assumption that sample is representative of the population

