---
title: "C++ Notes"
date: 2020-02-10
categories: ["Coding"]
tags: ["C++"]
comments: true
description: "Get answers for C/C++ within ? s"
---

# C++ Notes
Just some C/C++ code snippets to keep in mind. C/C++ is tremendous complicated, but it's still the most powerful programming language. 

### Table of Contents
* [char to Int](#char-to-int)  
* [Pointer](#pointer)  
* [Return Reference/Pointer](#return-reference)  
* [Array as Argument](#array) 
* [Operator that can't be overloaded](#operator-could-not-be-overloaded) 
* [Object Instantization](#object-instantization)  
* [Object Relationship](#object-relationships)  
* [Virtual Function and Ploymorphism](#virtual-functions-and-runtime-polymorphism)
* [Const](#const)  
* [Friend](#friend)  
* [Design Pattern](#design-pattern)  
  - [Singleton](#singleton)  
  - [Delegate](#delegate)



## Char to Int

#### 1. `Char` is `ASCII`  
C store `Char` as `ASCII` (Int) by default. So, `Char` is equal to ASCII code.

```cpp
char a = 'A'; // 65
int c = a; //c = 65
```

#### 2. `Char` to `Int`, to `String`, and vice versa

- `char` and `int`
    ```cpp
    char c='5'
    int res = c -'0' ; // 5

    int i=5;
    char res = I + '0'; // '5'
    ```

- `char*`,`int` and `string`
    ```cpp
    // char * to string
    const char * str_c =  "hello";
    std::string str = str_c;

    // string to char*
    str.c_str(); // return const char*

    // int to string
    std::to_string()

    // string to int
    std::stoi()
    ```

#### 3. `Char` as subscript of array: legal!  
This is useful when create `hashmap`. e.g. counting chars
```cpp
int test[200] = {0};
test['A'] = 1; // legal
test['b'] = 2;  // legal
```

[Back to top](#table-of-contents)

## Operator Could not be Overloaded

| operator | memo|
|---|---|
| .  | Member access or dot operator |
| .* | Pointer-to-member Operator |
| :: | Scope Resolution operator |
| ? : | conditional operator, ternary operator |
| sizeof | object size operator, built-in operations |
| typeid | object type operator, built-in operations |

OK, What's `.*` ?

```cpp
//we have a class
struct X
{
   void f() {}
   void g() {}
};

typedef void (X::*pointer)();
//ok, let's take a pointer and assign f to it.
pointer somePointer = &X::f;
//now I want to call somePointer. But for that, I need an object
X x;
//now I call the member function on x like this
(x.*somePointer)(); //will call x.f()
//now, suppose x is not an object but a pointer to object
X* px = new X;
//I want to call the memfun pointer on px. I use ->*
(px ->* somePointer)(); //will call px->f();
```
Now, you can't use x.somePointer(), or px->somePointer() because there is no such member in class X. 

see [here](https://stackoverflow.com/questions/6586205/what-are-the-pointer-to-member-and-operators-in-c)


[Back to top](#table-of-contents)

## Pointer
Pointer syntax  
**Rule:** read from right to left
```cpp
int a; // an int 
int *a; // a pointer point to int
int **a; // secondary int pointer, point to another int pointer 
int a[10]; // int array 
int *a[10]; // a poiter array, point to int 
int (*a)[10]; // a int pointer point to an int array 
int (*a)(int); // a pointer point to a function, will return an int
int (*a[10])(int); //  a poiter array, point to a function，will return an int
```

Declare two pointers
```cpp
int* a, b; // equal to int* a; int b;
int *a, *b; // correct way
```

[Back to top](#table-of-contents)

## Return Reference
Return Reference when define function
- return ref/pointer only if array or given object's ref/pointer. 
- return value only if temp object created 

```cpp
// ref
RMB& RMB::operator++()
{
    yuan++;
    return *this; // already existed object, created outside
}
// value
RMB RMB::operator++(int)
{
    RMB temp(yuan);//create object
    yuan++;
    return temp;//return a copy of temp
}
```
[Back to top](#table-of-contents)

## Array
### Array as `formal arguments`  
An Array could not copy to anther Array (**Copy pointer is not allowed!**), so `call-by-value` is not allowed.  
So, use array pointer:
```cpp
//these are same 
void print(const int*);
void print(const int[]);
void print(const int[5]);
```
multi-dimension array
```cpp
void print(const int(*p)[3], int rowsize);
void print(const int p[][5], int rowsize);
```
When use pointer to an Array, the dimension is unknown. So, need an extra argument to specify it explicitly.

Example:
```cpp
void print1(int (*p)[3]) {
    cout<<p[1][1]<<endl;
}
void print2(int p[][3]) {
    cout<<p[0][0]<<endl;
}

int a[2][3]={ {1,2},{3,4} };
print1(a); // 4
print2(a); // 1

int b[2][4]={ {1,2,5,6},{3,4,7,8} };
print1(b); // error
```
[Back to top](#table-of-contents)

## Object Instantization

### 1. without new

stack
```cpp
ClassName object(param); //  A a(1);  
ClassName object2 =  ClassName(param); //  A b = A(1);
```
### 2. with new

heap
```cpp
ClassName *object = new ClassName(param);//A *a = new A();
delete object;
```
### 3. copy constructor

```cpp
// 
```

### 4. Smart Pointer

```cpp 
std::unique_ptr<ClassName> object (new ClassName(param));
// recommend this way of instantization
std::unique_ptr<ClassName> object = std::make_unique<ClassName>(param);
```
[Back to top](#table-of-contents)

## Friend
The `friend` declaration appears in a `class body` and grants a `function` or `another class` **access** to `private` and `protected` **members** of the class where the friend declaration appears.

### 1. `friend function`

Declare anywhere inside a class, but define outside
```cpp
// function
friend <type> <Name>(<arguments>);
```

Example:
```cpp
class A
{
public:
	A(int _a):a(_a){};
    // non-member function
	friend int getA_a(A &_classA);
private:
	int a;
};

// without the friend keyword 
int getA_a(A &_classA)
{   //access member by formal arguments
	return _classA.a; 
}

A _classA(3);
std::cout<<getA_a(_classA); // 3
```

[Back to top](#table-of-contents)
### 2. `friend class` 

Delare inside class, define outside
```cpp
// class
friend class <Name>;
```
**Note:** `friend class X {};` is an error

Example:
```cpp
class B
{
public:
	B(int _b):b(_b){};
	friend class C; // friend class
private:
	int b;
};
 
class C
{
public:
	int getB_b(B _classB){
        //access member by formal arguments
		return _classB.b;
	};
};

B _classB(3);
C _classC; // an instance of a friend class
_classC.getB_b(_classB);
```

[Back to top](#table-of-contents)
### 3. Others: friend ostream, friend template ...

```cpp
class Y {
    int data; 
    // the non-member function operator<< will have access to Y's private members
    friend std::ostream& operator<<(std::ostream& out, const Y& o);
    friend char* X::foo(int); // members of other classes can be friends too
    friend X::X(char), X::~X(); // constructors and destructors can be friends
};

// this operator<< still needs to be defined, as a non-member
std::ostream& operator<<(std::ostream& out, const Y& y)
{   // can access private member Y::data
    return out << y.data; 
}
```

[Back to top](#table-of-contents)


## Object Relationships

relation types
  - "is-a"
  - "has-a"
  - "uses-a"
  - "depends-on"


| Property | Composition |	Aggregation	| Association |
|---|---|---|---|---|
| Relationship type | Whole/part | Whole/part | Otherwise unrelated |
| Members can belong to multiple classes | No | Yes | Yes |
| Members existence managed by class | Yes | No | No |
| Directionality |Unidirectional |Unidirectional | Unidirectional or bidirectional |
|Relationship verb | Part-of | Has-a | Uses-a |



### Composition: `has a data member`
Building complex objects from simpler ones is called **object composition** .

object composition models a __“has-a”__ relationship between two objects.
In C++, It means structs and classes can have data members of various types.

```cpp
class A；
class B 
{
public:
    B(){}
    ~B(){}
private:
    A a;
    int b;
}；
```

Summary:
1. Typically use normal member variables
2. Can use pointer members if the class handles object allocation/deallocation itself
3. Responsible for creation/destruction of parts

### Aggregation: "has a"
Unlike a composition, parts can belong to more than one object at a time, and the whole object is not responsible for the existence and lifespan of the parts.

Summary:
1. Typically use pointer or reference members that point to or reference objects that live outside the scope of the aggregate class
2. Not responsible for creating/destroying parts


### Association: "uses a"
Association models as “uses-a” relationship. The doctor “uses” the patient (to earn income). The patient uses the doctor (for whatever health purposes they need).


### Delegate: "has a"
or called `pImpl(Pointer to IMPLementation)`

Delegate: Composition by reference

has a pointer of another object

```cpp
class A；
class B 
{
public:
    B(){}
    ~B(){}
private:
    A *a;
    int b;
}；
```


### Inheritance: "is a"
public, protected, private

```cpp
class A
{
public:
    A(){}
    virtual ~A(){}
}
class B : public A
{
};
```


[Back to top](#table-of-contents)

## Virtual Functions and Runtime Polymorphism

1. Declare: vitrual keyword
    ```cpp
    class TestA {
    public:
        virtual void func() {  cout << "virtual function" << endl; }
    };

    class Test : public TestA {
    public:
        virtual void func() {  //  virtual could be omited 
            cout << "Test virtual function" << endl;
        }
        ~Test() { }
    };

    TestA* t = new Test; // parent pointer point to child (Ploymorphism)
    t->func();          // Results：Test virtual function
    delete t;

    ```

2. Member Could not be virtual
    * inline function
    * constructor
    * non-member function 
    * static function: only one copy of all objects.
    * friend function: it's non-member function
    * member function template !

3. Pure virtual function

  - declare
    ```cpp
    virtual void fun() = 0; 
    ```
  - class with pure virtual functoin could not be instantized!
  - a derived class of virtual class has to define pure virtual function. then the derived class could be instantized.
  - `abstract class`: class with pure virtual function
  


4. virtual deconstrutor
  - A parent pointer point to it's child. When delete the parent pointer, only parent constuctor is called. if declared a virtual deconstuctor, child's deconstuctor is called first, then the parent deconsturctor.
  - virtural keyword could be omited if a parent deconstructor is declared.
  - `delete` a pointer will only called object's deconstructor where the pointer point to.     


[Back to top](#table-of-contents)

## Const

### 1. `const` before or behind type/class, the syntax semantic are same
```cpp
// they are same 
const int x;  // (int x) is const/inmutable 
int const x;  // (const x) has type int 
```

### 2. Pointer with const: `const int* p`, `int const* p` and `int *const p`

__Dirty trick:__ use `*` as a separator, `const` restrict the type according to the side where it belong to    

**point to const:** These two expression are same 
```cpp
// -> (const int) | p;  p : a mutable pointer points to a const/immutable int
const int * p;

// -> (int const) | p; p2: a mutable pointer points a const which has type int
int const * p2;
```

**const pointer:** But these two not the same
```cpp
// -> int | (const p);  p3: a const pointer, point to an mutable int
int * const p3;

// -> (const int) | const p; p4: a const pointer, pointing to an immutable/const int
const int * const p4; 
```
### 3. class member `func` with const: `() const`  
a. `const object` 
- could not change class variable
- could not call `non-const` function

```cpp
class Number
{
public:
    void set(int num) { number = num; }
    int get() { return number; }
    int get2() const {return number;}
    int number = 0;
};

// Example
const Number n;
n.number = 1; // Error, n is const
n.set(1); // Error, n is const, non-const `set()`
n.get(); // Error, non-const `get()`
n.get2(); // OK
```

b. `() const`  
- could not change class variable, except static 
- could get variable

```cpp
class Number {
private:
  int a;
  static int b;
  const int c = 20;
public:
	void set() {  
         a = 10; // error when `this` argument has type 'const'
    void set2() const {
        b = 20; // OK
    }         	
	int get() const {  // OK
		return a; // did not change a
	}
};

const Number n;
n.set(); // Error
n.set2(); // OK
n.get(); // OK
```

Easy to understand, when pointer `this` is const
```cpp
void Number::set(const Number *const this, int num) { number = num; } // illegal -> const this
```

c. `() const` overloading

[Back to top](#table-of-contents)








## Design Pattern
### Singleton
Define
```cpp
class Singleton
{
    private:
        /* Here will be the instance stored. */
        static Singleton* instance;
        /* Private constructor to prevent instancing. */
        Singleton() {};
    public:
        /* Static access method. */
        static Singleton* getInstance()
        {
            if (instance == 0)
                instance = new Singleton();
            return instance;
        }
};

/* NULL, because instance will be initialized on demand. */
Singleton* Singleton::instance = 0;
```

Usage
```cpp
#include <iostream>
int main()
{
    //new Singleton(); // Won't work
    Singleton* s = Singleton::getInstance(); // Ok
    Singleton* r = Singleton::getInstance();

    /* The addresses will be the same. */
    std::cout << s << std::endl;
    std::cout << r << std::endl;
}
```

[Back to top](#table-of-contents)
### Delegate
Don't confuse with delegate constructor!!!


A delegate is a class that wraps a `pointer` or `reference` to an `object instance`, a member method of that object's class to be called on that object instance, and `provides a method to trigger` that call.

Example 1
```cpp
#include <iostream>
using namespace std;

class RealPrinter {
public:
    void print() { std::cout << "real-printer" << std::endl; }
};

class Printer {
public:
    Printer() : p(RealPrinter()) {}
    void print() { p.print(); }
private:
    RealPrinter p;
};

int main()
{
    Printer* printer = new Printer();
    printer->print();
}
```

Example 2:
```cpp
#include <iostream>
class I //interface {
public:
    virtual void f() = 0;
    virtual void g() = 0;
};

class A : public I {
public:
    void f(){std::cout << "A::f()" << std::endl;}
    void g(){std::cout << "A::g()" << std::endl;}
};

class B : public I {
public:
    void f(){std::cout << "B::f()" << std::endl;}
    void g(){std::cout << "B::g()" << std::endl;}
};


class C : public I {
public:
    C() { m_i = new A();/*delegation*/ }

    void f(){ m_i->f(); }
    void g(){ m_i->g(); }

    // normal attributes
    void toA(){ m_i = new A(); }
    void toB(){ m_i = new B(); }

private:
    I* m_i;
}

int main()
{
    C cc = C();
    cc.f();     // output: A::f()
    cc.g();     // output: A::g()

    cc.toB();
    cc.f();     // output: B::f()
    cc.g();     // output: B::g()
}
```
[Back to top](#table-of-contents)

### Composite
Composite is a structural design pattern that allows composing objects into a tree-like structure and work with the it as if it was a singular object.



[Back to top](#table-of-contents)