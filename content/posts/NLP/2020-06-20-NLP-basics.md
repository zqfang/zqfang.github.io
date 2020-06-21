---
title: 'NLP basic concepts'
date: 2020-6-20
categories: ["Nature Language Processing"]
comments: true
tags: ["Statistical Learning", "NLP"]
math: true
draft: true
---

Just some basic concepts
## Traditional NLP
### Count-based Representation
1. one-hot representation
2. BoW: Bag of words
Bow describes the `occurrence of words` within a document. including 
  - A Vocabulary of known words
  - A measure of the presence of known words, e.g. count


3. TF or TF-IDF representation: Term Frequency Inverse Document Frequency 
  - TF: the sum of the one-hot representation of a phrase, sentence or document's constituent words
    
    $$
    TF (w) = \frac { \text{ Number of the term w appears in the document }} { \text{Number of terms in the document}} 
    $$

  - IDF: penalizes common tokens and rewards rare tokens
  
    $$
    IDF(w) = \log \frac{\text{Number of documents}}{\text{Number of documents with term w}}
    $$

  - TF-IDF: $TF(w) * IDF(w)$




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

Stemming: use handcrafted rules to strip endings of words to reduce them to a common form called `stems`

### Categorizing Words: POS Tagging
part-of-speech (POS) tagging: labeling individual words or tokens

### Categorizing Spans: Chunking and Named Entity Recognition

a span of text: a contiguous multitoken boundary.

`chunking` or `shallow parsing`: identify the noun phrases (NP) and verb pharses (VP) in a span of text.

A `named entity` is a string mention of a real world concept like a person, location, organization, drug name ... 


### Structure of Sentences: Parse trees

Shallow parsing identifies `phrasal units`, the task of identifying the relationship between them is called `parsing`.

`Parse trees` indicate how different `grammatical units` in a sentence are related hiearachically. (constituent parse)

`dependency parsing`

### Word Senses and Semantics
Senses: the different meanings of a word



