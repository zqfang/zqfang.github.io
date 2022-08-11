---
title: "Rust Notes: NDarray"
date: 2022-02-10
categories: ["Coding"]
tags: ["Rust"]
comments: true
description: ""
hiddenFromHomePage: true
draft: true
---



## Vec to Array

1. simple case
```rust
use ndarray::{concatenate, s, Axis, Array, Array2, Array3, arr2, arr3};

let vec1 = vec![[1,2,3,4], [1,2,3,4]];
let array = Array::from_vec(vec1); // ndim = 1,  Nd will flatten into 1d by using from_vec.
let array2 = Array2<i64> = arr2(&vec1); // ndim = 2
```

2. nested Vecs/Arrays

see docs [here](https://docs.rs/ndarray/latest/ndarray/struct.ArrayBase.html)
```rust
// you know ahead-of-time the shape of the Array
let mut arr = Array2::zeros((2, 3));
for (i, mut row) in arr.axis_iter_mut(Axis(0)).enumerate() {
    // Perform calculations and assign to `row`; this is a trivial example:
    row.fill(i);
}
assert_eq!(arr, array![[0, 0, 0], [1, 1, 1]]);

// you don't know ahead-of-time the shape of the Array
// append data to a flat Vec, then conert it using ::from_shape_vec()

let ncols = 3;
let mut data = Vec::new();
let mut nrows = 0;
let arr = Array2::from_shape_vec((nrows, ncols), data)?;

// If neither of these options works for you
// using Iterator::flatten() then ::from_shape_vec()
let nested: Vec<Array2<i32>> = vec![
    array![[1, 2, 3], [4, 5, 6]],
    array![[7, 8, 9], [10, 11, 12]],
];
let inner_shape = nested[0].dim();
let shape = (nested.len(), inner_shape.0, inner_shape.1);
let flat: Vec<i32> = nested.iter().flatten().cloned().collect();
let arr = Array3::from_shape_vec(shape, flat)?;
assert_eq!(arr, array![
    [[1, 2, 3], [4, 5, 6]],
    [[7, 8, 9], [10, 11, 12]],
]);
```


## indexing and slicing


```rust
println!("slice)

```
