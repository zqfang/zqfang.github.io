# C++ Notes: Return Pointer or Reference



## Explaination

A pointer or reference could not be return if they point/refer to a local variable stored in stack inside a function (local variable stored in stack will be destoried automatically when return, and the pointer become `wild`)

Situations when a function **could return pointer or reference**

* variable defined outside a function scope
* global variable
* local static variable
* local variable stored in heap ( `new` opterator, `malloc()`)

Other process could not access the memory of variable stored in heap until it is released. That's why.


## Example

```cpp
#include <iostream>
#include <string.h>
#include <stdlib.h>
 
using namespace std;
 
string& f1(const string &s)
{
    static string result = s;
    return result;
}
 
string &f2(const string &s)
{
    string *p = new string;
    *p = s;
    return *p;
}
   
int *f3()
{
   int *a = (int *)malloc(sizeof(int) * 10);
   *a = 10;
   *(a + 1) = 11;
   return a;
}
 
int &f4()
{
   int *a = (int *)malloc(sizeof(int) * 10);
   *a = 10;
   *(a + 1) = 11;
   return *a;
}
 
int main()
{
    int *a = &f7();
    cout<<(*(a + 1))<<endl;
    free(a); // free the memory when done.
    return 0;
}
```

Return `*this`, or alrealy exist objets
```cpp
// ref
A& A::operator++()
{
    count++;
    return *this; // already existed object, created outside
}
```

## What's the defference between `new` and `malloc()`

`malloc()` is a function that **takes a number (of bytes)** as its argument; it __returns a void*__ pointing to unitialized storage.  
`new` is an operator that **takes a type** and (optionally) a set of initializers for that type as its arguments; it **returns a pointer to an** (optionally) initialized **object of its type**.

The difference is most obvious when you want to allocate an object of a user-defined type with non-trivial initialization semantics
