---
title: 'Pointer in UnSafe Swift'
date: 2020-05-30
categories: ["Coding"]
draft: false
tags: ["Swift"]
comments: true
description: "Get answers for Swift within 10s"
---


## Pointers
Unsafe Swift pointers use a predictable naming scheme:

<strong>Unsafe</strong> <span style="color: red">[Mutable]</span><span style="color: green">[Raw]</span><span style="color: turquoise">[Buffer]</span><strong>Pointer</strong><span style="color: magenta">[\<T>]</span> 



Explain:  

<b>Pointers</b> are just memory addresses.   
Direct memory access is <b>Unsafe</b>.  
<span style="color: red">Mutable</span> means you can write to it.  
<span style="color: green">Raw</span> means it points to a blob of bytes.  
<span style="color: turquoise">Buffer</span> means that is works like a collection.  
Generic <span style="color: magenta">[\<T>]</span> pointers are typed.



## Working with Pointers

| C Pointer | Swift Type |
| --- | --- |
| int * | UnsaftMutablePointer<Int32> |
| const int * | UnsafePointer<Int32> |
| NSDate ** | AutoreleasingUnsafeMutablePointer<NSDate> |
| struct UnknowType * | OpaquePointer |
| void * | UnsafeMutableRawPointer |
| const void * | UnsafeRawPointer



Explain: see [here](https://swift.gg/2016/12/13/swift-and-c-everything-you-need-to-know/)


## Usage

see [here](https://www.raywenderlich.com/7181017-unsafe-swift-using-pointers-and-interacting-with-c) and [here](https://swift.gg/2016/12/13/swift-and-c-everything-you-need-to-know/)
