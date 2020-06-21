---
title: "C++ Notes 2"
date: 2020-02-10
categories: ["Coding"]
tags: ["C++"]
comments: true
draft: true
description: "Get answers for C/C++ within ? s"
---



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



## extern, static
global variable
  - defined outside all functions and available to all functions.
  - unaffected by scopes and are always available (exists until the program ends)

extern
  - declare a `global variable` (exists on the whole project): variable could be used in `multi- .cpp` files
  - extern "C" {/* c code*/}: compile c code.

static
  - declare a `local global variable` (file scope): only be accessed in its translation unit or `.o file`, that's, in the file where it is created.
  - declare a static class member (class scope): initialization should be `outside class body`
    - static data member
    - static function:
      - no `this` pointer: only access to other static member/function
      - could declare as private

## constexpr

The constexpr specifier `declares` that it is possible to `evaluate the value` of the function or variable at `compile time`.