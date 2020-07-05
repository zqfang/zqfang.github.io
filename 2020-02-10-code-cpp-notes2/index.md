# C++ Notes: Advanced 1


Just some advanced C/C++ code snippets to keep in mind.

## Header file naming

Never used some header file name with `std`. Sometimes, compiler could not find the std headers.!!!

## `default` and `delete` in class

special class member:

- default constructor 
- deconstructor
- copy constructor
- operater `=`

when use default and delete

1. default

```cpp
class X {
public:
   X()=default; // with this, you could declare like this: X x;
   X(int){};
};
X x; // works
```

2. delete: prohibit func call marked by `delete`

```cpp
class X {
public: 
    X(); 
    X(const X&) = delete;
    X& operator = (const X &) = delete;
}; 
// example 
X x1; 
X x2=x1;   // Error, copy constructor is prohibited
```


## Pointer and Smart Pointer cast

```cpp
#include <memory> // smart pointer header

// [pomter to smart pointer
struct Base {};
struct Derived: Base {};

// pointer to smart
Base *p1 = new Derived(); // upcast
std::shared_ptr<Base> sp(p1);

//  a polymorphic type
Base *p = new Derived(); // upcast, dynamic_cast is unnecessary
Derived* dp = dynamic_cast<Derived*> (p); // downcast

// smart pointer convert to pointer
std::shared_ptr<Base> smart = std::make_shared<Derived>();
Base* p2 = smart.get(); // .get()

// smart pointer cast
// downcast
std::shared_ptr<Derived> dsmart = std::dynamic_pointer_cast<Derived>(smart);

// upcast
// case 1
std::shared_ptr<Base> foo(new Derived());
// case 2
std::shared_ptr<Derived> bar = std::make_shared<Base>();
std::shared_ptr<Base> foo = std::dynamic_pointer_cast<A>(bar);
```

## FileIO

The simplest example

```cpp
#include <iostream>
#include <fstream>

// output file
std::ofstream output;
output.open("test.compact.txt");
// read input file
string line;
std::ifstream input("test.chrX.vcf");
if (input.is_open()) {
    while (getline(input, line))
        output << line <<'\n';
}
input.close();
output.close();

```

## Strings: split and strip

split string by delimiter

```cpp
std::vector<std::string> split(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}
```

strip strings

```cpp
std::string trim(const std::string& str, const std::string delimiter = " \n\r\t")
{
//    std::string s;
//    s.erase(s.find_last_not_of(" \n\r\t")+1);
    size_t first = str.find_first_not_of(delimiter);
    if (std::string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(delimiter);
    return str.substr(first, (last - first + 1));
}
```

## Containor

### remove duplicated elements: O(nlogn)

```cpp
sort(nums.begin(),nums.end()); // inplace
// unique do not change vector size, only put dup elements to end of containor
// and return a iter which points to the first non-uniqdup element 
vector<int>::iterator iter = unique(nums.begin(), nums.end());
nums.erase(iter, nums.end()); //remove duplciates inplace
// nums.resize( std::distance(nums.begin(),iter) );
```

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

### Composite

Composite is a structural design pattern that allows composing objects into a tree-like structure and work with the it as if it was a singular object.

