---
title: 'NLP concepts'
date: 2020-6-20
categories: ["Nature Language Processing"]
comments: true
tags: ["NLP"]
math: true
draft: true
---

Basic concepts in NLP

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

   - TF-IDF: $TF(w) \times IDF(w)$

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

### Word Senses and Semantics

Senses: the different meanings of a word

### Categorizing Words: POS Tagging

part-of-speech (POS) tagging: labeling individual words or tokens

### Categorizing Spans: Chunking and Named Entity Recognition

a span of text: a contiguous multi-token boundary.

`chunking` or `shallow parsing`: identify the noun phrases (NP) and verb pharses (VP) in a span of text.

A `named entity` is a string mention of a real-world concept like a person, location, organization, drug name, et. al.

### Structure of Sentences: Parse trees

Shallow parsing identifies `phrasal units`, the task of identifying the relationship between them is called `parsing`.

![parsing](/images/nlp/parsing.png)

1. `Parse trees` indicate how different `grammatical units` in a sentence are related hierachically. (aslo refer to `constituent parse`, chart-based )

2. `dependency parsing`: directed graph (graph-based)
   - node -> word
   - edge -> relation
   - all the words have one incoming edge, except ROOT
   - there is a unique path from each word to ROOT

![parsing](/images/nlp/parsing.classifier.png)

maximum spanning tree

SyntaxNet:  
![syntaxnet](/images/nlp/parsingtree.gif)

## Reference

Hung-yi Lee: [Deep Learning for human language processing](https://www.youtube.com/watch?v=9erBrs-VIqc)
