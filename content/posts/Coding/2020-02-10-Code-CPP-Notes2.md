---
title: "C++ Notes 2"
date: 2020-02-10
categories: ["Coding"]
tags: ["C++"]
comments: true
draft: true
description: "Get answers for C/C++ within ? s"
hiddenFromHomePage: true
---

Just some useful code snippets.

## Pointer and Smart Pointer cast

```cpp
#include <memory>

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

## extern, static

global variable

- defined outside all functions and available to all functions.
- unaffected by scopes and are always available (exists until the program ends)

extern

- declare a `global variable` (exists on the whole project): variable could be used in `multi- .cpp` files
- `extern "C" {/* c code */}`: compile c code.

static

- declare a `local global variable` (file scope): only be accessed in its translation unit or `.o file`, that's, in the file where it is created.
- declare a static class member (class scope): initialization should be `outside class body`
  - static data member
  - static function:
    - no `this` pointer: only access to other static member/function
    - could declare as private

## constexpr

The constexpr specifier `declares` that it is possible to `evaluate the value` of the function or variable at `compile time`.

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

## strings: split and strip

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
