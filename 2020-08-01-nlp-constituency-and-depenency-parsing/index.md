# NLP: Parse trees



## Structure of Sentences: Parse trees

Shallow parsing identifies `phrasal units`, the task of identifying the relationship between them is called `parsing`.

![parsing](/images/nlp/parsing.png)

1. `Parse trees` indicate how different `grammatical units` in a sentence are related hierachically. (aslo refer to `constituent parse`, chart-based )

2. `dependency parsing`: directed graph (graph-based)
   - node -> word
   - edge -> relation
   - all the words have one incoming edge, except ROOT
   - there is a unique path from each word to ROOT

![parsing](/images/nlp/parsing.classifier.png)


Contex-free grammars

common application:

- grammer checking
- semantic analysis
- question anwsering
- information extraction

## Constituency parsing (Syntatic parsing)

The task of recognizing a sentence and assigning a syntactic structure to it. 

### CKY Parsing

The most widely used dynamic-porgramming based approach to parsing. See also Earley algorithm and chart parsing

## Dependency parsing

Dependency grammars `grammatical relation` provides the basis for the binary relations that comprise dependency structures.

Dependency grammars allow to classify the kinds of grammatical relations, or mammatical function

## Dependency Formalisms

dependency structures are simply directed graphs: $G = (V,A)$, which refer to as arcs

**Dependency tree** is a directed graph that statisfies the following constrains:

1. there is a single designated root node that has no incoming arcs
2. with teh exception of the root node, each vertex has exactly one incoming arc.
3. there is a unique path from the root node to each vertex in $V$

## Transition-Based Dependency parsing

Shift-reduce parsing

SyntaxNet:  
![syntaxnet](/images/nlp/parsingtree.gif)

## Graph-Based Dependency parsing

motivations:

- capable of producing non-projective trees
- parsing accuarcy, particularly with respect to longer dependencies.

### Maximun spanning tree (MST)

1. construct a fully-connected, weighted, directed, rooted graph where the vertices are input words and the directed edges represent *all possible* head-dependent assignments. The weights reflect the score for each possible head-dependent relation.
   - every vertex in a spanning tree has exactly one incoming edge
   - absolute values of the edge scores are not critical to determining its maximum spanning tree. But relative weights of the edges entering each vertex that matters


## Reference

Hung-yi Lee: [Deep Learning for human language processing](https://www.youtube.com/watch?v=9erBrs-VIqc)
