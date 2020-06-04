# What is Big O




### What on earth is **Big O**?  
Time complexity and space complexity


### Time complexity
O(f(n)): number of commands need to execute. proportional to f(n).   
表示运行算法所需要执行的指令数，和f(n)成正。 
  
严格来讲，O(f(n))表示算法执行的上界。业界默认为算法执行的最低上界(最坏情况）。   
n represents the data scale 数据规模

when n is a large number, the constant is usually ignored.

| algorithm  | n of cmd   |
|:---------- |:---------- |
| Binary reserach $O(logn)$  | $a*logn$ |
| Max/min in an array $O(n)$ | b*n|
| merge sort $O(nlogn)$      | $c*nlogn$|
| select sort $O(n^2)$     | $d*n^2$ |
| quick sort $O(nlogn)$      | e*nlogn|
| adjacent graph           | $O(V+E)$ |
| Lazy Prim                | $O(ElogE)$ |
| Prim                     | $O(ElogV)$ |
| Kruskal                  | $O(ElogE)$ | 
| Dijkstra                 | $O(ElogV)$ |
| Bellman-Ford             | $O(EV)$ | 

minimum span tree  
Shortest path tree (Single source shortest path)

### Space complexity

| cmd | complexity|
|------------|--------|
| new an array | $O(n)$  
| new 2d array | $O(n^2)$  
| new an constant space | $O(1)$

recursive function: the depth (n) of a recursive function, the extra space need $O(n)$.

### Make sense of n  
If you want to solve the problem in 1 second,  
then an algorithm of  

| complexity | cmds n|
|------------|--------|
| $O(n^2)$   | could exec cmd n = $10^4$ |  
| $O(n)$     | could exec cmd n = $10^8$ | 
| $O(nlogn)$ | could exec cmd n = $10^7$ |


### Example

1. binarySearch

```
    find from n element  
    find from $n/2$ element  
    find from $n/4$  element
    ...
    find from 1 
```

That's, need how many steps of search when n = 1?  $log_{2}n = O(logn)$.

2. int2string. Set num > 0

```cpp
string int2string(int num) {
  string s="";
  while(num) {
    s += '0' + num%10;
    num /= 10;
  }
  reverse(s); // O(n) 
  return s;
}
```
That is, how many "/10" steps when num = 0? $log_{10}n = O(logn)$.

3. Case: two nested for loop, not always $O(n^2)$

```cpp
void hello(int n){
  for (int sz =1; sz < n; sz ++ sz) // logn here
      for( int i=1; i < n;; i++)  //n
          cout<<"hello, complex" <<endl;
}
```
So, should be $O(nlogn)$

4. isPrime: $O(\sqrt{n})$  
```cpp
// set n > 1
bool isPrime(int n){
  for( int x =2; x*x <= n; x++){
    if( n%x == 0) return false;
    return true;
  }
}
```
5. recursive function

* single recursive function call
```cpp
int binarySearch(int arr[], int l, int r, int target) {
  if (l>r) return -1;
  int mid = l +(r-l)/2;
  if (arr[mid] == target) return mid;
  else if (arr[mid] > target)
      return binarySearch(arr, ;, mid-1, target);
  else
      return binarySearch(arr, mid+1, r, target);
}
```
each step need O(1), so overall complexity depend on recursive exec depth.  
That is, if each function call needs time T, then time complexity: O(T*depth) -> O(n).

Another example: recursion depth logn, them time complexity O(logn).
```cpp
double pow( double x, int n){
  assert(n >=0);
  if (n==0) return 1.0;
  double t = pow(x, n/2);
  if( n%2) return x*t*t;
  return t*t;
}
```


* multi recursive exec
how many exec step? 
```cpp
int f(int n) {
    assert(n >=0);
    if(n == 0) return 1;
    return f(n-1) + f(n);
}
```
that's, count how many nodes on a full binary tree. $2^{n+1} -1 = O(2^n)$

how to think about this?
```cpp
void mergeSort(int arr[]. int l. int r){
  if (l >=r) return;
  int mid = (l+r) /2;
  mergeSort(arr, l, mid);
  mergeSort(arr, mid+1, r);
  merge(arr, l, mid, r);
}
```
For binary tree, complexity for each level O(n), while tree depth O(logn). Overall, O(nlogn)


6. Amortized time  
i.e. dynamic vector/stack/deque

```cpp
template<T>
class MyVector{
private:
    T* data;
    int capacity;
    int size;
    //O(n)
    void resize(int newCapacity){
      assert(newCapacity >= size);
      T* newData = new T[newCapacity];
      for(int i = 0; i < size; i++ ){
        newData[i] = data[i];
      }
      delete[] data;
      data = newData;
      capacity = newCapacity; 
    }
public: 
    MyVector() {
      data = new T[10];
      capacity = 10;
      size = 0;
    }
    ~MyVector() {
      delete[] data;
    }
    // Average: O(1)
    void push_back(T e){
      //assert(size < capacity)
      if (size == capacity)
          resize (2 *capacity);
      data[size++] = e;    
    }
    // Average O(1)
    T pop_back(){
      assert (size >0);
      T ret = data[size-1];
      size --;
      // note the denominator here. To advoid ossilation of space complexity 
      if(size == capacity / 4)
          resize(capacity /2); 
      return ret;
    }
};
```



