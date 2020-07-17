# C++ Notes:  STL


Just some advanced C/C++ code snippets to keep in mind.

## Header naming

Never used some header file name with `std`. Sometimes, compiler could not find the std headers.!!!

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

### Remove elements

1. Vector or Deque: algorithm remove() followed by erase()
   
```cpp
vector<int> vec = {1,1,2,3,4,4,6};
auto itr = remove(vec.begin(), vec.end(), 4);
vec.erase(iter, vec.end());
vec.shrink_to_fit(); // reduce capacity 
```

2. List: member function .remove() 
3. Associative Container or Unordered Container: .erase()


## Lambda

1. Syntax
```cpp
auto basicLambda = [] { cout << "Hello, world!" << endl; };
basicLambda();   // Hello, world!


// return type
auto add = [](int a, int b) -> int { return a + b; };
// inference return type
auto multiply = [](int a, int b) { return a * b; };
```

2. Capture parameters
   
```cpp
int x = 10;   
auto add_x = [x](int a) { return a + x; };  // copy capture x
auto multiply_x = [&x](int a) { return a * x; };  // ref capture x
```

- []：默认不捕获任何变量；
- [=]：默认以值捕获所有变量；
- [&]：默认以引用捕获所有变量；
- [ x ]：仅以值捕获x，其它变量不捕获；
- [&x]：仅以引用捕获x，其它变量不捕获；
- [=, &x]：默认以值捕获所有变量，但是x是例外，通过引用捕获；
- [&, x]：默认以引用捕获所有变量，但是x是例外，通过值捕获；
- [this]：通过引用捕获当前对象（其实是复制指针）；
- [*this]：通过传值方式捕获当前对象；

1. capture expression

```cpp
// capture by expression
int x = 4;
auto y = [&r = x, x = x + 1] { r += 2; return x * x; }();
// x = 6，y = 25

// initialize directly
auto z = [str = "string"]{ return str; }();
// z: const char *
```

2. generic: `auto`

```cpp
auto add = [](auto x, auto y) { return x + y; };

int x = add(2, 3);   // 5
double y = add(2.5, 3.5);  // 6.0
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

