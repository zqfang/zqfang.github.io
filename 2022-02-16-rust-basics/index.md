# Rus: Advanced Trait


Some codesnape for the usage of Trait
## Trait
A trait defines functionality a particular type has and can share with other types.

### Defined and Implement a Trait

```rust
pub trait Summary {
    fn summarize(&self) -> String;
}


pub struct Tweet {
    pub username: String,
    pub content: String,
    pub reply: bool,
    pub retweet: bool,
}

impl Summary for Tweet {
    fn summarize(&self) -> String {
        format!("{}: {}", self.username, self.content)
    }
}

```


### Default Implementataions
```rust
pub trait Summary {
    fn summarize(&self) -> String {
        String::from("(Read more...)")
    }
}

```
You could overwirte trait methods when implement structs!

### Trait as Parameters


```rust
pub fn notify(item: &impl Summary) {
    println!("Breaking news! {}", item.summarize());
}
```

### Trait Bound Synatx

- `<T: Trait1 + Trait2>`
- `function<T, U>() -> T where T: Trait1, U: Trait2`

```rust
pub fn notify<T: Summary>(item: &T) {
    println!("Breaking news! {}", item.summarize());
}

// multiple trait bound
pub fn notify<T: Summary + Display>(item: &T) {}

// where
fn some_function<T, U>(t: &T, u: &U) -> i32
where
    T: Display + Clone,
    U: Clone + Debug,
{
    unimplemented!()
}

```

### Returning Instances that implement Traits

- `-> impl SomeTrait `

```rust
pub trait Summary {
    fn summarize(&self) -> String;
}

pub struct NewsArticle {
    pub headline: String,
    pub location: String,
    pub author: String,
    pub content: String,
}

impl Summary for NewsArticle {
    fn summarize(&self) -> String {
        format!("{}, by {} ({})", self.headline, self.author, self.location)
    }
}

pub struct Tweet {
    pub username: String,
    pub content: String,
    pub reply: bool,
    pub retweet: bool,
}

impl Summary for Tweet {
    fn summarize(&self) -> String {
        format!("{}: {}", self.username, self.content)
    }
}

fn returns_summarizable() -> impl Summary {
    Tweet {
        username: String::from("horse_ebooks"),
        content: String::from(
            "of course, as you probably already know, people",
        ),
        reply: false,
        retweet: false,
    }
}

```


### Trait Object
Define
- `Box<dyn Draw>`
- `&dyn Draw`
- `impl<T> SomeStruct<T> where T: Draw `

```rust
trait Draw {
    fn draw(&self) -> String;
}

// NOTE: Box<dyn Draw> and &dyn Draw are both worked 
fn draw1(x: Box<dyn Draw>) {
    x.draw(); // Deref to .
}

fn draw2(x: &dyn Draw) {
    x.draw();
}

// NOTE: used in a vec
pub struct Screen {
    pub components: Vec<Box<dyn Draw>>,
}
impl Screen {
    pub fn run(&self) {
        for component in self.components.iter() {
            component.draw();
        }
    }
}
```

use generic 
```rust
pub struct Screen<T: Draw> {
    pub components: Vec<T>,
}
impl<T> Screen<T>
    where T: Draw {
    pub fn run(&self) {
        for component in self.components.iter() {
            component.draw();
        }
    }
}
```

### Associated Type
Note: Associated Type has nothing to do with associated function
- define inside trait
- has to be assign a type in `impl`
- `type Ithem`

```rust
// Item has to be defined in the impl
pub trait Iterator {
    type Item;
    fn next(&mut self) -> Option<Self::Item>;
}

```

### Default generic type's parameter


```rust
// RHS=Self -> default to its own type 
trait Add<RHS=Self> {
    type Output;
    fn add(self, rhs: RHS) -> Self::Output;
}

struct Point {
    x: i32,
    y: i32,
}

impl Add for Point {
    type Output = Point; // defined here

    fn add(self, other: Point) -> Point {
        Point {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

fn main() {
    assert_eq!(Point { x: 1, y: 0 } + Point { x: 2, y: 3 },
               Point { x: 3, y: 3 });
}
```

### Call method with same name

Struct's member function has priority!

```rust
trait Pilot {
    fn fly(&self);
}

trait Wizard {
    fn fly(&self);
}

struct Human;

impl Pilot for Human {
    fn fly(&self) {
        println!("This is your captain speaking.");
    }
}

impl Wizard for Human {
    fn fly(&self) {
        println!("Up!");
    }
}

impl Human {
    fn fly(&self) {
        println!("*waving arms furiously*");
    }
}
```

example
```rust
fn main() {
    let person = Human;
    Pilot::fly(&person); // call Pilot's fly 
    Wizard::fly(&person); // call Wizard's fly 
    person.fly(); // call self.fly.  Human's fly 
}
```

### Call method (without &self) with same name 

use `as` to limit!!!
- `<Dog as Animal>::baby_name()`

```rust
trait Animal {
    fn baby_name() -> String;
}

struct Dog;

impl Dog {
    fn baby_name() -> String {
        String::from("Spot")
    }
}

impl Animal for Dog {
    fn baby_name() -> String {
        String::from("puppy")
    }
}
```

```rust
fn main() {
    println!("A baby dog is called a {}", <Dog as Animal>::baby_name());
}
```








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

