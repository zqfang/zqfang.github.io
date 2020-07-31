# NLP: RNN and Self-attension


## Backpropagation Through Time

## Long Short-Term Memory

1. Delete information from the context that is no longer needed: `Forget Gate` f

$$
f_t = \sigma (U_f h_{t-1} + W_f X_t)
$$

$$
k_t = c_{t-1} \odot f_t
$$

2. Compute the actual information we need to extract from the previous hidden stat and current inputs

$$
g_t = \tanh (U_g h_{t-1} + W_g x_t)
$$

3. Select the information to add to the current context: `Add Gate` i 

$$
i_t = \sigma (U_i h_{t-1} + W_i X_t)
$$

$$
j_t = g_{t} \odot i_t
$$

4. Get new context vector

$$
c_t = j_t + k_t
$$

5. `Output Gate` o: decide what information is required for the current hiddent state (as opposed to what information need to be preseved for future decicions) 

$$
o_t = \sigma (U_o h_{t-1} + W_o x_t)
$$

$$
h_t = o_t \odot \tanh (c_t)
$$


## Gated Recurrent Units

GRU ease the tranning burden by dispensing with the use of a separate context vector, and by reducing the number of gates to 2: 

1. a reset gate, $r$: decide which aspects of the previous hidden state are relevant to the current context and what can be ignored. 

$$
r_t = \sigma (U_r h_{t-1} + W_r x_t)
$$

Then computing an intermediate representation for the new hidden stat at time $t$

$$
\tilde h_t = \tanh (U(r_t \odot h_{t-1}) + Wx_t)
$$


2. an update gate, $z$: detemine which aspects of the new intermedicate representation will be used directly and which aspects of the previous stat need to be preseverd for future use
   
$$
z_t = \sigma (U_z h_{t-1} + W_z x_t)
$$

$$
h_t = (1- z_t)h_{t-1} + z_t \tilde h_t
$$


## Attention mechanism

Consider Encoder to Deconder Network, a decoder

$$
h_i^d = g(\hat y_{i-1}, h_{i-1}^d, c_i)
$$

1. computing context vector $c_i$:  a vector of scores that capture the relevance of each encoder hidden state to the decoder state captured in $h_{i-1}^d$. That's, at each state $i$ during decoding, we'll compute $score(h_{i-1}^d, h_j^e)$ for each encoder state $j$. Recall similarity

$$
score(h_{i-1}^d, h_j^e) = h_{i-1}^d \cdot h_j^e
$$

2. make a more robust similarity score by adding a learnable weights, $W_s$:

$$
score(h_{i-1}^d, h_j^e) = h_{i-1}^d W_s h_j^e
$$

3. normalize the scores

$$
\begin{aligned}
\alpha_{ij} &= \operatorname{softmax} (score(h_{i-1}^d, h_j^e)) \cr 
&= \frac {\exp (score(h_{i-1}^d, h_j^e))} {\sum_k score(h_{i-1}^d, h_j^e)} 
\end{aligned}
$$

4. finally, give $\alpha$,

$$
c_i = \sum_j \alpha_{ij}h_j^e
$$

![attention](/images/nlp/encoder-decoder.png)


## Self-Attention

1. create 3 vectors from each of the encoders' input vector, a `Query` vector, a `Key` vector and a `Value` vector, then multiplying the embedding (of word) `X`

    ![attention1](/images/nlp/self-attention-1.png)

2. calculate a score by taking the dot product of the `Query` with the `Key`
3. divide the scores by the square root of the dimension of the key vector (a more stable gradients), then pass the result grought a softmax.

4. multiply each `Value` vector by the softmax score 
5. sum up the weighted value vectors

![attention2](/images/nlp/self-attention-2.png)

6. **multi-head attention**: to focus on different region and give "representation subspace"

![multiheads](/images/nlp/self-attention-3.png)



## Transformer

A transformer of two stacked encoder and decoder looks like this

![transformer](/images/nlp/transformers.png)

1. transformer use `positoinal encoding` vector to representing the order of the sequence. It follows a specific pattern that the model learns, which helps it determine the position of each word, or the distance between different words in the sequence.

2. transformer use `LayerNorm`


## Reference

[transformer](https://jalammar.github.io/illustrated-transformer/)
