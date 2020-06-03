# Swift Cheat Sheet


# Swift Cheat Sheet
Stolen from [iwasrobbed](https://github.com/iwasrobbed/Swift-CheatSheet).   
I simplify it and add some more. 

It's a high level and a quick reference to Swift.
The purpose of this cheat sheet is to teach myself and get answers within 10s. 


### Table of Contents

* [Code Document](#code-documentation)
* [Data Types](#data-types)
* [Operators](#operators)
* [Operator Overloading](#operator-overloading)
* [Declaring Classes](#declaring-classes)
* [Declarations](#declarations)
* [Lazy Property](#lazy-property)
* [Property Observer](#property-observer)
* [Literals](#literals)
* [Functions](#functions)
* [Constants and Variables](#constants-and-variables)
* [Naming Conventions](#naming-conventions)
* [Closures](#closures)
* [Generics](#generics)
* [Control Statements](#control-statements)
* [Extension](#extension)
* [Protocol](#protocol)
* [Protocol Extension](#protocol-extension)
* [Error Handling](#error-handling)
* [Passing Information](#passing-information)
* [User Defaults](#user-defaults)
* [Common Patterns](#common-patterns)
* [Unicode Support](#unicode-support)
* [File IO](#fileio)

## Code Documentation
Two ways of commenting: 
  - //
  - /* ... */

Two ways of documenting with markdown (Reconigzed by xcode):
  - ///
  - /**  ... */

### Markdown
a few keywords that xcode can recognized automatically, with the format like **- <Keyword>**.  

The most common:  Prameters, Throws, Returns
```swift
/**
- Prameters: 
  - argument1: This is arg1
  - argument2: This is arg2 
- Returns: The results string.
- Throws: `Error` if nil 
*/
```

Other keywords
```swift
/**
- Precondition:
- Postcondition: 
- Requires: All the information in the object should be sorted
- Invariant: The object will maintained sorted
- Complexity: O(n^2)
- Important:
- Warning: Very computation consuming
- Attention: Same as Warning
- Note: something to keep in mind
- Remark: Same as note
*/
```

Metadata
```swift
/**
- Author: 
- Authors:
- Copyright:
- Date:
- Since:
- Version:
*/
```

[Back to top](#table-of-contents)

### MARK

Using `MARK` to organize your code:

```swift
// MARK: - Use mark to logically organize your code
// Declare some functions or variables here
// MARK: - They also show up nicely in the properties/functions list in Xcode
// Declare some more functions or variables here
```

### FIXME

Using `FIXME` to remember to fix your code:

```swift
// Some broken code might be here
// FIXME: Use fixme to create a reminder to fix broken code later
```

`FIXME` works a lot like `MARK` because it makes organizing code easier, but it's used exclusively when you need to remember to fix something.

### TODO

Using `TODO` to remember to add, delete, or generally refactor your code:

````swift
// Some incomplete code might be here
// TODO: Use todo to create a reminder to finish things up later
````

`TODO` is very similar to `FIXME` and `MARK`, but it's used exclusively when you need to remember to add, delete, or change your code later.

**Auto-generating method documentation:**
In a method's preceding line, press `⌥ Option + ⌘ Command + /` to automatically generate a documentation stub for your method.

[Back to top](#table-of-contents)

## Data Types

### Size

Permissible sizes of data types are determined by how many bytes of memory are allocated for that specific type and whether it's a 32-bit or 64-bit environment.  In a 32-bit environment, `long` is given 4 bytes, which equates to a total range of `2^(4*8)` (with 8 bits in a byte) or `4294967295`.  In a 64-bit environment, `long` is given 8 bytes, which equates to `2^(8*8)` or `1.84467440737096e19`.

For a complete guide to 64-bit changes, please [see the transition document](https://developer.apple.com/library/mac/documentation/Darwin/Conceptual/64bitPorting/transition/transition.html#//apple_ref/doc/uid/TP40001064-CH207-TPXREF101).

### C Primitives

Unless you have a good reason to use C primitives, you should just use the Swift types to ensure compability going foward.

In fact, Swift just aliases C types to a Swift equivalent:

```swift
// C char is aliased as an Int8 and unsigned as UInt8
let aChar = CChar()
let anUnsignedChar = CUnsignedChar()
print("C char size: \(MemoryLayout.size(ofValue: aChar)) with min: \(Int8.min) and max: \(Int8.max)")
// C char size: 1 with min: -128 and max: 127
print("C unsigned char size: \(MemoryLayout.size(ofValue: anUnsignedChar)) with min: \(UInt8.min) and max: \(UInt8.max)")
// C unsigned char size: 1 with min: 0 and max: 255

// C short is aliased as an Int16 and unsigned as UInt16
let aShort = CShort()
let unsignedShort = CUnsignedShort()
print("C short size: \(MemoryLayout.size(ofValue: aShort)) with min: \(Int16.min) and max: \(Int16.max)")
// C short size: 2 with min: -32768 and max: 32767
print("C unsigned short size: \(MemoryLayout.size(ofValue: unsignedShort)) with min: \(UInt16.min) and max: \(UInt16.max)")
// C unsigned short size: 2 with min: 0 and max: 65535

// C int is aliased as an Int32 and unsigned as UInt32
let anInt = CInt()
let unsignedInt = CUnsignedInt()
print("C int size: \(MemoryLayout.size(ofValue: anInt)) with min: \(Int32.min) and max: \(Int32.max)")
// C int size: 4 with min: -2147483648 and max: 2147483647
print("C unsigned int size: \(MemoryLayout.size(ofValue: unsignedInt)) with min: \(UInt32.min) and max: \(UInt32.max)")
// C unsigned int size: 4 with min: 0 and max: 4294967295

// C long is aliased as an Int and unsigned as UInt
let aLong = CLong()
let unsignedLong = CUnsignedLong()
print("C long size: \(MemoryLayout.size(ofValue: aLong)) with min: \(Int.min) and max: \(Int.max)")
// C long size: 8 with min: -9223372036854775808 and max: 9223372036854775807
print("C unsigned long size: \(MemoryLayout.size(ofValue: unsignedLong)) with min: \(UInt.min) and max: \(UInt.max)")
// C unsigned long size: 8 with min: 0 and max: 18446744073709551615

// C long long is aliased as an Int64 and unsigned as UInt64
let aLongLong = CLongLong()
let unsignedLongLong = CUnsignedLongLong()
print("C long long size: \(MemoryLayout.size(ofValue: aLongLong)) with min: \(Int64.min) and max: \(Int64.max)")
// C long long size: 8 with min: -9223372036854775808 and max: 9223372036854775807
print("C unsigned long long size: \(MemoryLayout.size(ofValue: unsignedLongLong)) with min: \(UInt64.min) and max: \(UInt64.max)")
// C unsigned long long size: 8 with min: 0 and max: 18446744073709551615
```

From the [docs](https://developer.apple.com/library/ios/documentation/swift/conceptual/buildingcocoaapps/InteractingWithCAPIs.html):

| C Type | Swift Type
| :---: | :---:
| bool | CBool
| char, signed char | CChar
| unsigned char | CUnsignedChar
| short | CShort
| unsigned short | CUnsignedShort
| int | CInt
| unsigned int | CUnsignedInt
| long | CLong
| unsigned long | CUnsignedLong
| long long | CLongLong
| unsigned long long | CUnsignedLongLong
| wchar_t | CWideChar
| char16_t | CChar16
| char32_t | CChar32
| float | CFloat
| double | CDouble

#### Integers

Integers can be signed or unsigned.  When signed, they can be either positive or negative and when unsigned, they can only be positive.

Apple states: _Unless you need to work with a specific size of integer, always use `Int` for integer values in your code. This aids code consistency and interoperability. Even on 32-bit platforms, `Int` [...] is large enough for many integer ranges._

**Fixed width integer types with their accompanying byte sizes as the variable names:**

```swift
// Exact integer types
let aOneByteInt: Int8 = 127
let aOneByteUnsignedInt: UInt8 = 255
let aTwoByteInt: Int16 = 32767
let aTwoByteUnsignedInt: UInt16 = 65535
let aFourByteInt: Int32 = 2147483647
let aFourByteUnsignedInt: UInt32 = 4294967295
let anEightByteInt: Int64 = 9223372036854775807
let anEightByteUnsignedInt: UInt64 = 18446744073709551615

// Minimum integer types
let aTinyInt: Int8 = 127
let aTinyUnsignedInt: UInt8 = 255
let aMediumInt: Int16 = 32767
let aMediumUnsignedInt: UInt16  = 65535
let aNormalInt: Int32  = 2147483647
let aNormalUnsignedInt: UInt32 = 4294967295
let aBigInt: Int64 = 9223372036854775807
let aBigUnsignedInt: UInt64 = 18446744073709551615

// The largest supported integer type
let theBiggestInt: IntMax = 9223372036854775807
let theBiggestUnsignedInt: UIntMax = 18446744073709551615
```

#### Floating Point

Floats cannot be signed or unsigned.

```swift
// Single precision (32-bit) floating-point. Use it when floating-point values do not require 64-bit precision.
let aFloat = Float()
print("Float size: \(MemoryLayout.size(ofValue: aFloat))")
// Float size: 4

// Double precision (64-bit) floating-point. Use it when floating-point values must be very large or particularly precise.
let aDouble = Double()
print("Double size: \(MemoryLayout.size(ofValue: aDouble))")
// Double size: 8
```

#### Boolean

```swift
// Boolean
let isBool: Bool = true // Or false
```

In Objective-C comparative statements, `0` and `nil` were considered `false` and any non-zero/non-nil values were considered `true`. However, this is not the case in Swift. Instead, you'll need to directly check their value such as `if x == 0` or `if object != nil`

#### Primitives

**nil** : Used to specify a null object pointer.  When classes are first initialized, all properties of the class point to `nil`.

### Enum & Bitmask Types

Enumeration types can be defined as follows:

```swift
// Specifying a typed enum with a name (recommended way)
enum UITableViewCellStyle: Int {
    case default, valueOne, valueTwo, subtitle
}

// Accessing it:
let cellStyle: UITableViewCellStyle = .default
```

As of Swift 3, all enum options should be named in lowerCamelCased.

#### Working with Bitmasks

Newer Swift versions have a nice substitute for the old `NS_OPTIONS` macro for creating bitmasks to compare to.

An example for posterity:

```swift
struct Options: OptionSet {
    let rawValue: Int

    init(rawValue: Int) {
        self.rawValue = rawValue
    }

    init(number: Int) {
        self.init(rawValue: 1 << number)
    }

    static let OptionOne = Options(number: 0)
    static let OptionTwo = Options(number: 1)
    static let OptionThree = Options(number: 2)
}

let options: Options = [.OptionOne, .OptionTwo]

options.contains(.OptionOne) // true
options.contains(.OptionThree) // false
```

### Type Casting

Sometimes it is necessary to cast an object into a specific class or data type.  Examples of this would be casting from a `Float` to an `Int` or from a `UITableViewCell` to a subclass such as `RPTableViewCell`.

#### Checking Types

Swift uses `is` and `as` both for checking object types as well as conformance to a given protocol.

#### Operator: is

Checking object type using `is`:

```swift
if item is Movie {
    movieCount += 1
    print("It is a movie.")
} else if item is Song {
    songCount += 1
    print("It is a song.")
}

```

The `is` operator returns `true `if an instance is of that object type, or conforms to the specified protocol, and returns `false` if it does not.

#### Operators: as? and as!

If you want to be able to easily access the data during one of these checks, you can use `as?` to optionally (or `as!` to force) unwrap the object when necessary:

```swift
for item in library {
    if let movie = item as? Movie {
        print("Director: \(movie.director)")
    } else if let song = item as? Song {
        print("Artist: \(song.artist)")
    }
}
```

The `as?` version of the downcast operator returns an optional value of the object or protocol's type, and this value is `nil` if the downcast fails or this instance does not conform to the specified protocol.

The `as!` version of the downcast operator forces the downcast to the specified object or protocol type and triggers a runtime error if the downcast does not succeed.

#### Casting from Generic Types

If you're working with `AnyObject` objects given from the Cocoa API, you can use:

```swift
for movie in someObjects as! [Movie] {
    // do stuff
}
```

If given an array with `Any` objects, you can use a `switch` statement with the type defined for each `case`:

```swift
var things = [Any]()

for thing in things {
    switch thing {
    case 0 as Int:
        print("Zero as an Int")
    case let someString as! String:
        print("S string value of \"\(someString)\"")
    case let (x, y) as! (Double, Double):
        print("An (x, y) point at \(x), \(y)")
    case let movie as! Movie:
        print("A movie called '\(movie.name)' by director \(movie.director)")
    default:
        print("Didn't match any of the cases specified")
    }
}
```

#### Basic Casting

Swift also offers some simple methods of casting between it's given data types.

```swift
// Example 1:
let aDifferentDataType: Float = 3.14
let anInt: Int = Int(aDifferentDataType)

// Example 2:
let aString: String = String(anInt)
```

[Back to top](#table-of-contents)

## Operators

Swift supports most standard C operators and improves several capabilities to eliminate common coding errors. The assignment operator `=` does not return a value, to prevent it from being mistakenly used when the equal to operator `==` is intended.

Arithmetic operators (`+`, `-`, `*`, `/`, `%`) detect and disallow value overflow, to avoid unexpected results when working with numbers that become larger or smaller than the allowed value range of the type that stores them.

#### Arithmetic Operators

Operator | Purpose
| :---: | ---
| + | Addition
| - | Subtraction
| * | Multiplication
| / | Division
| % | Remainder

#### Comparative Operators

| Operator | Purpose
| :---: | ---
| == | Equal to
| === | Identical to
| != | Not equal to
| !== | Not identical to
| ~= | Pattern match
| > | Greater than
| < | Less than
| >= | Greater than or equal to
| <= | Less than or equal to

#### Assignment Operators

| Operator | Purpose
| :---: | ---
| = | Assign
| += | Addition
| -= | Subtraction
| *= | Multiplication
| /= | Division
| %= | Remainder
| &= | Bitwise AND
| &#124;= | Bitwise Inclusive OR
| ^= | Exclusive OR
| <<= | Shift Left
| >>= | Shift Right

#### Logical Operators

| Operator | Purpose
| :---: | ---
| ! | NOT
| && | Logical AND
| &#124;&#124; | Logical OR

#### Range Operators

| Operator | Purpose
| :---: | ---
| ..< | Half-open range
| ... | Closed range

#### Bitwise Operators

| Operator | Purpose
| :---: | ---
| & | Bitwise AND
| &#124; | Bitwise Inclusive OR
| ^ | Exclusive OR
| ~ | Unary complement (bit inversion)
| << | Shift Left
| >> | Shift Right

#### Overflow and Underflow Operators

Typically, assigning or incrementing an integer, float, or double past it's range would result in a runtime error. However, if you'd instead prefer to safely truncate the number of available bits, you can opt-in to have the variable overflow or underflow using the following operators:

| Operator | Purpose
| :---: | ---
| &+ | Addition
| &- | Subtraction
| &* | Multiplication

Example for unsigned integers (works similarly for signed):

```swift
var willOverflow = UInt8.max
// willOverflow equals 255, which is the largest value a UInt8 can hold
willOverflow = willOverflow &+ 1
// willOverflow is now equal to 0

var willUnderflow = UInt8.min
// willUnderflow equals 0, which is the smallest value a UInt8 can hold
willUnderflow = willUnderflow &- 1
// willUnderflow is now equal to 255
```

#### Other Operators

| Operator | Purpose
| :---: | ---
| ?? | Nil coalescing
| ?: | Ternary conditional
| ! | Force unwrap object value
| ? | Safely unwrap object value

[Back to top](#table-of-contents)

## Operator Overloading

Swift allows you to overwrite existing operators or define new operators for existing or custom types. For example, this is why in Swift you can join strings using the `+` operator, even though it is typically used for math.

Operator overloading is limited to the following symbols, `/ = - + * % < > ! & | ^ . ~`, however you cannot overload the `=` operator by itself (it must be combined with another symbol).

Operators can be specified as:
* `prefix`: goes before an object such as `-negativeNumber`
* `infix`: goes between two objects, such as `a + b`
* `postfix`: goes after an object, such as `unwrapMe!`

### Custom operators

* `associativity`: defines how operators of the same precedence are grouped together (left, right)
* `precedence`: gives some operators higher priority than others; these operators are applied first.

Refer [Operator Declarations](https://developer.apple.com/documentation/swift/swift_standard_library/operator_declarations) to see full details about operator `associativity` and `precedence`.

**Example:**  

DefaultPrecedence group
```swift
// declare first and set rules with a precedence group
infix operator ** // use DefaultPrecedence group
```

Custom Precedence group
```swift
// define a custom precedence group
precedencegroup ExponentiationPrecedence {
    higherThan: MultiplicationPrecedence
    associativity: right // none, left, right
    //assignment: false
}

// now, replace original declaration of ** with 
infix operator **: ExponentiationPrecedence
```
That's it.

```swift
// impelment 
infix func ** (x: Double, p: Double) -> Double {
    return pow(x, p)
} 
2**3  // 8
2**3**2 // 512
1+2**3**2 // 513
5*2**3**2  // 2560
```

see also [docs](https://docs.swift.org/swift-book/LanguageGuide/AdvancedOperators.html#//apple_ref/doc/uid/TP40014097-CH27-ID41)

[Back to top](#table-of-contents)

## Declaring Classes

Classes are typically declared using separate `.swift` files, but multiple classes can also be created within the same file if you'd like to organize it that way.

Unlike Objective-C, there's no need for an interface file (`.h`) in Swift.

The implementation file should contain (in this order):
* Any needed `import` statements
* A `class` declaration which contains any constants or variables necessary for the class
* All public and private functions

Example:

MyClass.swift

```swift
import UIKit

class MyClass {
    // Declare any constants or variables at the top
    let kRPErrorDomain = "com.myIncredibleApp.errors"
    var x: Int, y: Int
    // MARK: - Class Methods, e.g. MyClass.functionName()
    class func alert() {
        print("This is a class function.")
    }
    // MARK: - Instance Methods, e.g. myClass.functionName()
    init(x: Int, y: Int) {
        self.x = x
        self.y = y
    }
    // MARK: - Private Methods
    private func pointLocation() -> String {
        return "x: \(x), y: \(y)"
    }
}
```

#### Instantiation

When you want to create a new instance of a class, you use the syntax:

```swift
let myClass = MyClass(x: 1, y: 2)
```

where `x` and `y` are variables that are passed in at the time of instantiation.

[Back to top](#table-of-contents)

## Declarations

More info [here in the docs](https://developer.apple.com/library/ios/documentation/swift/conceptual/Swift_Programming_Language/Declarations.html#//apple_ref/doc/uid/TP40014097-CH34-XID_704).

#### Preprocessor

Swift doesn't come with a preprocessor so it only supports a limited number of statements for build time. Things like `#define` have been replaced with global constants defined outside of a class.

| Directive | Purpose
| :---: | ---
| #if | An `if` conditional statement
| #elif | An `else if` conditional statement
| #else | An `else` conditional statement
| #endif | An `end if` conditional statement

#### Imports

| Directive | Purpose
| :---: | ---
| import | Imports a framework

#### Constants & Variables

| Directive | Purpose
| :---: | ---
| let | Declares local or global constant
| var | Declares a local or global variable
| class | Declares a class-level constant or variable
| static | Declares a static type

#### Classes, Structure, Functions and Protocols

| Directive | Purpose
| :---: | ---
| typealias | Introduces a named alias of an existing type
| enum | Introduces a named enumeration
| struct | Introduces a named structure
| class | Begins the declaration of a class
| init | Introduces an initializer for a class, struct or enum
| init? | Produces an optional instance or an implicitly unwrapped optional instance; can return `nil`
| deinit | Declares a function called automatically when there are no longer any references to a class object, just before the class object is deallocated
| func | Begins the declaration of a function
| protocol | Begins the declaration of a formal protocol
| static | Defines as type-level within struct or enum
| convenience | Delegate the init process to another initializer or to one of the class’s designated initializers
| extension | Extend the behavior of class, struct, or enum
| subscript | Adds subscripting support for objects of a particular type, normally for providing a convenient syntax for accessing elements in a collective, list or sequence
| override | Marks overriden initializers

#### Operators

| Directive | Purpose
| :---: | ---
| operator | Introduces a new infix, prefix, or postfix operator

#### Declaration Modifiers

| Directive | Purpose
| :---: | ---
| dynamic | Marks a member declaration so that access is always dynamically dispatched using the Objective-C runtime and never inlined or devirtualized by the compiler
| final | Specifies that a class can’t be subclassed, or that a property, function, or subscript of a class can’t be overridden in any subclass
| lazy | Indicates that the property’s initial value is calculated and stored at most once, when the property is first accessed
| optional | Specifies that a protocol’s property, function, or subscript isn’t required to be implemented by conforming members
| required | Marks the initializer so that every subclass must implement it
| weak | Indicates that the variable or property has a weak reference to the object stored as its value

#### Access Control

| Directive | Purpose
| :---: | ---
| open | Can be subclassed outside of its own module and its methods overridden as well; truly open to modification by others and useful for framework builders
| public | Can only be subclassed by its own module or have its methods overridden by others within the same module
| internal | (Default) Indicates the entities are only available to the entire module that includes the definition, e.g. an app or framework target
| fileprivate | Indicates the entities are available only from within the source file where they are defined
| private | Indicates the entities are available only from within the declaring scope within the file where they are defined (e.g. within the `{ }` brackets only)


```swift
public class AccessLevelsShowCase {
    
    // Property accessible for other modules
    public var somePublicProperty = 0

    // Property accessible from the module 
    var someInternelProperty = 1

    // Property accessible from its own defining source file
    fileprivate func someFilePrivateMethod() {}

    // Property accessible fro its enclosing declaration
    private func somePrivateMethod() {}
}
```
[Back to top](#table-of-contents)

## Literals

Literals are compiler directives which provide a shorthand notation for creating common objects.

| Syntax | What it does
| :---: | ---
| `"string"` | Returns a `String` object
| `28` | Returns an `Int`
| `3.14`, `0xFp2`, `1.25e2` | Returns a `Double` object
| `true`, `false` | Returns a `Bool` object
| `[]` | Returns an `Array` object
| `[keyName:value]` | Returns a `Dictionary` object
| `0b` | Returns a binary digit
| `0o` | Returns an octal digit
| `0x` | Returns a hexadecimal digit

#### Strings

Special characters can be included:

* Null Character: `\0`
* Backslash: `\\` (can be used to escape a double quote)
* Horizontal Tab: `\t`
* Line Feed: `\n`
* Carriage Return: `\r`
* Double Quote: `\"`
* Single Quote: `\'`
* Unicode scalar: `\u{n}` where n is between one and eight hexadecimal digits

Multiline string literal
```swift
let json = """
{ "username": "David", "loginCount": 2}
"""
```
#### Array Access Syntax

```swift
let example = [ "hi", "there", 23, true ]
print("item at index 0: \(example[0])")
```

#### Dictionary Access Syntax

```swift
let example = [ "hi" : "there", "iOS" : "people" ]
if let value = example["hi"] {
    print("hi \(value)")
}
```

#### Mutability

For mutable literals, declare it with `var`; immutable with `let`.

[Back to top](#table-of-contents)

## Functions

#### Declaration Syntax

Functions without a return type use this format:

```swift
// Does not return anything or take any arguments
func doWork() {
    // Code
}
```

`class` precedes declarations of class functions:

```swift
// Call on a class, e.g. MyClass.someClassFunction()
class func someClassFunction() {
    // Code
}
```

`static` is similar to class functions where you don't need an instance of the class or struct in order to call a method on it:

```swift
// Call on a class/struct, e.g. MyStruct.someStaticFunction()
static func someStaticFunction() {
    // Code
}
```

Declare instance functions:

```swift
// Called on an instance of a class, e.g. myClass.someInstanceFunction()
func doMoreWork() {
    // Code
}
```

Function arguments are declared within the parentheses:

```swift
// Draws a point
func draw(point: CGPoint)
```

Return types are declared as follows:

```swift
// Returns a String object for the given String argument
func sayHelloToMyLilFriend(lilFriendsName: String) -> String {
    return "Oh hello, \(lilFriendsName). Cup of tea?"
}
```

You can have multiple return values, referred to as a tuple:

```swift
// Returns multiple objects
func sayHelloToMyLilFriend(lilFriendsName: String) -> (msg: String, nameLength: Int) {
    return ("Oh hello, \(lilFriendsName). Cup of tea?", countElements(lilFriendsName))
}

var hello = sayHelloToMyLilFriend("Rob")
print(hello.msg) // "Oh hello, Rob. Cup of tea?"
print(hello.nameLength) // 3
```

And those multiple return values can be optional:

```swift
func sayHelloToMyLilFriend(lilFriendsName: String) -> (msg: String, nameLength: Int)?
```

By default, external parameter names are given when you call the function, but you can specify that one or more are not shown in the method signature by putting a `_` symbol in front of the parameter name:

```swift
func sayHelloToMyLilFriend(_ lilFriendsName: String) {
    // Code
}

sayHelloToMyLilFriend("Rob")
```

or you can rename the variable once within the method scope:

```swift
func sayHelloToMyLilFriend(friendsName lilFriendsName: String) {
    // Code
}

sayHelloToMyLilFriend(friendsName: "Rob") // and local variable is `lilFriendsName`
```


You can also specify default values for the parameters:

```swift
func sayHelloToMyLilFriend(_ lilFriendsName: String = "Rob") {
    // Code
}

sayHelloToMyLilFriend() // "Oh hello, Rob. Cup of tea?"
sayHelloToMyLilFriend("Jimbob") // "Oh hello, Jimbob. Cup of tea?"
```

Swift also supports variadic parameters so you can have an open-ended number of parameters passed in:

```swift
func sayHelloToMyLilFriends(_ lilFriendsName: String...) {
    // Code
}

sayHelloToMyLilFriends("Rob", "Jimbob", "Cletus")
// "Oh hello, Rob, Jimbob and Cletus. Cup of tea?"
```

And lastly, you can also use a prefix to declare input parameters as `inout`.

An in-out parameter has a value that is passed in to the function, is modified by the function, and is passed back out of the function to replace the original value.

You may remember `inout` parameters from Objective-C where you had to sometimes pass in an `&error` parameter to certain methods, where the `&` symbol specifies that you're actually passing in a pointer to the object instead of the object itself. The same applies to Swift's `inout` parameters now as well.

#### Calling Functions

Functions are called using dot syntax: `myClass.doWork()` or `self.sayHelloToMyLilFriend("Rob Phillips")`

`self` is a reference to the function's containing class.

At times, it is necessary to call a function in the superclass using `super.someMethod()`.

[Back to top](#table-of-contents)

## Constants and Variables

Declaring a constant or variable allows you to maintain a reference to an object within a class or to pass objects between classes.

Constants are defined with `let` and variables with `var`. By nature, constants are obviously immutable (i.e. cannot be changed once they are instantiated) and variables are mutable.

```swift
class MyClass {
	let text = "Hello" // Constant
	var isComplete: Bool // Variable
}
```

There are many ways to declare properties in Swift, so here are a few examples:

```swift
var myInt = 1 // inferred type
var myExplicitInt: Int = 1 // explicit type
var x = 1, y = 2, z = 3 // declare multiple variables

let (a,b) = (1,2) // declare multiple constants
```

#### Getters and Setters

In Objective-C, variables were backed by getters, setters, and private instance variables created at build time. However, in Swift getters and setters are only used for computed properties and constants actually don't have a getter or setter at all.

The getter is used to read the value, and the setter is used to write the value. The setter clause is optional, and when only a getter is needed, you can omit both clauses and simply return the requested value directly. However, if you provide a setter clause, you must also provide a getter clause.

You can overrride the getter and setter of a property to create the illusion of the Objective-C property behavior, but you'd need to store them as a private property with a different name (not recommended for most scenarios):

```swift
private var _x: Int = 0

var x: Int {
    get {
        print("Accessing x...")
        return _x
    }
    set {
        print("Setting x...")
        _x = newValue
    }
}
```

#### Property Observer
Swift also has callbacks for when a property will be or was set using `willSet` and `didSet` shown below:

* `willset`: before assignment
* `didSet`: after assignment

```swift
class LightBulb {
    static let maxCurrent = 30
    var current = 0 {
        willSet(newCurrent) { // do something before value assignment
          // newValue -> newCurrent
          print("Current value changed, the change is \(abs(current- newCurrent))")
        }
        didSet { // do somthing afther value assignment
          if current == LightBulb.maxCurrent { 
              print("current get to maximum point")
              }
            // oldValue  
        }
    }
}
let bulb = LightBulb()
bulb.current = 20
bulb.current = 30
bulb.current = 40
```
[Back to top](#table-of-contents)


#### Lazy Property
* `lazy`: only compute once and remember the value, won't re-compute if called again.

```swift
class ClosedRange {
    let start: Int
    let end: Int
    var width: Int {
        return end - start +1
    }
    // note the =
    lazy var sum: Int = {
        var res = 0
        print("run")
        for i in self.start...self.end{
            res += 1
        }
        return
    }() // don't forget ()

    init?(start: Int, end: Int){
        if start > end {
            return nil
        }
        self.start = start
        self.end = end 
    }
}

// example
if let range = ClosedRange(start: 0, end: 10_000) {
    range.width //1001
    range.sum // will print out "run"
    range.sum // 
    range.sum // 
}
```
[Back to top](#table-of-contents)

#### Accessing

#### Local Variables

Local variables and constants only exist within the scope of a function.

```swift
func doWork() {
    let localStringVariable = "Some local string variable."
    self.doSomething(string: localStringVariable)
}
```

[Back to top](#table-of-contents)

## Naming Conventions

The general rule of thumb: Clarity and brevity are both important, but clarity should never be sacrificed for brevity.

#### Functions and Properties

These both use `camelCase` where the first letter of the first word is lowercase and the first letter of each additional word is capitalized.

#### Class names and Protocols

These both use `CapitalCase` where the first letter of every word is capitalized.

### Enums

The options in an enum should be `lowerCamelCased`

#### Functions

These should use verbs if they perform some action (e.g. `performInBackground`).  You should be able to infer what is happening, what arguments a function takes, or what is being returned just by reading a function signature.

Example:

```swift
// Correct
func move(from start: Point, to end: Point) {}

// Incorrect (likely too expressive, but arguable)
func moveBetweenPoints(from start: Point, to end: Point) {}

// Incorrect (not expressive enough and lacking argument clarity)
func move(x: Point, y: Point) {}
```

[Back to top](#table-of-contents)

## Closures

Closures in Swift are similar to blocks in Objective-C and are essentially chunks of code, typically organized within a `{}` clause, that are passed between functions or to execute code as a callback within a function. Swift's `func` functions are actually just a special case of a closure in use.

#### Syntax

```swift
{ (params) -> returnType in
    statements
}
```

#### Examples

```swift
// Map just iterates over the array and performs whatever is in the closure on each item
let people = ["Rob", "Jimbob", "Cletus"]
people.map({
    (person: String) -> String in
    "Oh hai, \(person)..."
})
// Oh hai, Rob
// Oh hai, Jimbob
// Oh hai, Cletus

// Closure for alphabetically reversing an array of names, where sorted is a Swift library function
let names = ["Francesca", "Joe", "Bill", "Sally", ]
var reversed = names.sorted { (s1: String, s2: String) -> Bool in
    return s1 > s2
}
// Or on a single line:
reversed = names.sorted{ (s1: String, s2: String) -> Bool in return s1 > s2 }
// Or because Swift can infer the Bool type:
reversed = names.sorted { s1, s2 in return s1 > s2 }
// Or because the return statement is implied:
reversed = names.sorted { s1, s2 in s1 > s2 }
// Or even shorter using shorthand argument names, such as $0, $1, $2, etc.:
reversed = names.sorted { $0 > $1 }
// Or just ridiculously short because Swift's String greater-than operator implementation exactly matches this function definition:
reversed = names.sorted(by: >)
```

If the closure is the last parameter to the function, you can also use the trailing closure pattern. This is especially useful when the closure code is especially long and you'd like some extra space to organize it:

```swift
func someFunctionThatTakesAClosure(closure: () -> ()) {
    // function body goes here
}

// Instead of calling like this:
someFunctionThatTakesAClosure({
    // closure's body goes here
})

// You can use trailing closure like this:
someFunctionThatTakesAClosure() {
    // trailing closure's body goes here
}
```
#### Capturing Values

A closure can capture constants and variables from the surrounding context in which it is defined. The closure can then refer to and modify the values of those constants and variables from within its body, even if the original scope that defined the constants and variables no longer exists.

In Swift, the simplest form of a closure that can capture values is a nested function, written within the body of another function. A nested function can capture any of its outer function’s arguments and can also capture any constants and variables defined within the outer function.

```swift
func makeIncrementor(forIncrement amount: Int) -> () -> Int {
    var runningTotal = 0
    func incrementor() -> Int {
        runningTotal += amount
        return runningTotal
    }
    return incrementor
}
```

Swift determines what should be captured by reference and what should be copied by value. You don’t need to annotate a variable to say that they can be used within the nested function. Swift also handles all memory management involved in disposing of variables when they are no longer needed by the function.

#### Capturing Self

If you create a closure that references `self.*` it will capture `self` and retain a strong reference to it. This is sometimes the intended behavior, but often could lead to retain cycles where both objects won't get deallocated at the end of their lifecycles.

The two best options are to use `unowned` or `weak`. This might look a bit messy, but saves a lot of headache.

Use `unowned` when you know the closure will only be called if `self` still exists, but you don't want to create a strong (retain) reference.

Use `weak` if there is a chance that `self` will not exist, or if the closure is not dependent upon `self` and will run without it. If you do use `weak` also remember that `self` will be an optional variable and should be checked for existence.

```swift
typealias SomeClosureType = (_ value: String) -> ()

class SomeClass {
    fileprivate var currentValue = ""

    init() {
        someMethod { (value) in // Retained self
            self.currentValue = value
        }

        someMethod { [unowned self] (value) in // Not retained, but expected to exist
            self.currentValue = value

        }

        someMethod { [weak self] value in // Not retained, not expected to exist
            // Or, alternatively you could do
            guard let sSelf = self else { return }

            // Or, alternatively use `self?` without the guard
            sSelf.currentValue = value
        }
    }

    func someMethod(closure: SomeClosureType) {
        closure("Hai")
    }
}
```

Reference: [Apple: Automatic Reference Counting](https://developer.apple.com/library/ios/documentation/Swift/Conceptual/Swift_Programming_Language/AutomaticReferenceCounting.html)

[Back to top](#table-of-contents)

## Generics

Coming soon...

[Back to top](#table-of-contents)

## Control Statements

Swift uses all of the same control statements that other languages have:

#### If-Else If-Else

```swift
if someTestCondition {
    // Code to execute if the condition is true
} else if someOtherTestCondition {
    // Code to execute if the other test condition is true
} else {
    // Code to execute if the prior conditions are false
}
```

As you can see, parentheses are optional.

#### Ternary Operators

The shorthand notation for an `if-else` statement is a ternary operator of the form: `someTestCondition ? doIfTrue : doIfFalse`

Example:

```swift
func stringForTrueOrFalse(trueOrFalse: Bool) -> String {
    return trueOrFalse ? "True" : "False"
}
```

#### Nil Coalescing Operators

In Swift, we need to consider the use of `optional` values. One very basic way to handle `nil` cases is with an `if-else` statement:

```swift
func stringForOptionalExistence(optionalValue: String?) -> String {
  if optionalValue != nil {
    return optionalValue
  } else {
    return "Empty"
  }
}
```

In this particular case, we are returning `optionalValue` if it is not `nil`, and `"Empty"` if `optionalValue` is `nil`. The shorthand notation for this type of `if(!=nil)-else` statement is a nil coalescing operator of the form: `optionalValue ?? nonOptionalValue`

Example:

```swift
func stringForOptionalExistence(optionalValue: String?) -> String {
  return optionalValue ?? "Empty"
}
```

#### For Loops

Swift enables you to use ranges inside of `for` loops now:

```swift
for index in 1...5 {
    print("\(index) times 5 is \(index * 5)")
}

// Or if you don't need the value of the index
let base = 3, power = 10
var answer = 1
for _ in 1...power {
    answer *= base
}
print("\(base) to the power of \(power) is \(answer)")
// prints "3 to the power of 10 is 59049"
```


#### Enumerating arrays & dictionaries

```swift
// We explicitly cast to the Movie class from AnyObject class
for movie in someObjects as [Movie] {
    // Code to execute each time
}

// Enumerating simple array
let names = ["Anna", "Alex", "Brian", "Jack"]
for name in names {
    print("Hello, \(name)!")
}

// Enumerating simple dictionary
let numberOfLegs = ["spider": 8, "ant": 6, "cat": 4]
for (animalName, legCount) in numberOfLegs {
    print("\(animalName)s have \(legCount) legs")
}
```

If you need to cast to a certain object type, see the earlier discussion about the `as!` and `as?` keywords.

#### While Loop

```swift
while someTestCondition {
   // Code to execute while the condition is true
}
```

#### Repeat While Loop

```swift
repeat {
    // Code to execute while the condition is true
} while someTestCondition
```

#### Switch

Switch statements are often used in place of `if` statements if there is a need to test if a certain variable matches the value of another constant or variable.  For example, you may want to test if an error code integer you received matches an existing constant value or if it's a new error code.

```swift
switch errorStatusCode {
    case .network:
        // Code to execute if it matches

     case .wifi:
        // Code to execute if it matches

     default:
        // Code to execute if nothing else matched
}
```

Switch statements in Swift do not fall through the bottom of each case and into the next one by default. Instead, the entire switch statement finishes its execution as soon as the first matching switch case is completed, without requiring an explicit `break` statement. This makes the switch statement safer and easier to use than in C, and avoids executing more than one switch case by mistake.

#### Exiting Loops

Although `break` is not required in Swift, you can still use a `break` statement to match and ignore a particular case, or to break out of a matched case before that case has completed its execution.

* `return` : Stops execution and returns to the calling function.  It can also be used to return a value from a function.
* `break` : Used to stop the execution of a loop.

[Back to top](#table-of-contents)

## Extension
Extensions add new functionality to an existing class, structure, enumeration or protocol type
```swift
extension String {
    // Extending String type to calculate if a String instance is truthy of falsy
    var boolValue:Bool {
        if self == "1"
        return true
    }
    return false
}
let isTrue = "0".boolValue
```
[Back to top](#table-of-contents)
## Protocol
Define
```swift
protocol Codable {
    // definitions
    var description: String 
    var mustBeSettable: Int { get set }
    var doesNotNeedToBeSettable: Int { get }
    func dance () -> Double
    static func someTypeMethod()
    mutating func toggle() // modify (or mutate) the instance it belongs to
    init(someParameter: Int) // require specific initializers 
}
```

Usage
```swift
import Foundation
struct UserInfo: Codable {
    let username: String
    let loginCount: Int
}

extension UserInfo: CustomStringConvertible {
    var description: String {
        return "\(username) has tried to login \(loginCount) time(s)"
    }
}
```
### Protocol Extension
Protocols can be extended to provide method, initializer, subscript, and computed property implementations to conforming types. 

**Very Important and Useful:**  `Implementation` to any `method` or computed `property` requirement of that `protocol` can only be in `extension`

```swift
extension RandomNumberGenerator {
    func randomBool() -> Bool {
        return random() > 0.5
    }
}
```
By creating an extension on the protocol, all conforming types automatically gain this method implementation without any additional modification.

You can use protocol extensions to provide a default implementation to any method or computed property requirement of that protocol. 
```swift
extension PrettyTextRepresentable  {
    var prettyTextualDescription: String {
        return textualDescription
    }
}
```

[Back to top](#table-of-contents)
## Error Handling

Representing an Error
```swift
enum BeverageMachineError: Error {
    case invalidSelection
    case insufficientFunds
    case outOfStock
}
func selectBeverage (_ selection: Int) throws -> String {
    // do something
    return "Waiting for beverage..."
}

// us do...catch to handle error throwed by func
let message:String
do {
    message = try selectBeverage(20)
} catch BeverageMachineError.invalidSelection {
    print("Invalid selection")
} catch BeverageMachineError.insufficientFunds {
    print("Insufficient Funds")
} catch BeverageMachineError.outOfStock {
    print("Out of Stock")
} catch {
    print ("Generic error")
}

// if throw error, return nil
let message = try? selectBeverage(10)
// if throw error, get a runtime error
let message = try! selectBeverage(10)
```

[Back to top](#table-of-contents)

## Passing Information

Coming soon...

[Back to top](#table-of-contents)

## User Defaults

User defaults are basically a way of storing simple preference values which can be saved and restored across app launches.  It is not meant to be used as a data storage layer, like Core Data or sqlite.

### Storing Values

```swift
let userDefaults = UserDefaults.standard
userDefaults.setValue("Some Value", forKey: "RPSomeUserPreference")
```

### Retrieving Values

```swift
let userDefaults = UserDefaults.standard
let someValue = userDefaults.value(forKey: "RPSomeUserPreference") as AnyObject?
```

There are also other convenience functions on `UserDefaults` instances such as `bool(forKey:...)`, `string(forKey:...)`, etc.

[Back to top](#table-of-contents)

## Common Patterns

For a comprehensive list of design patterns, as established by the Gang of Four, look here: [Design Patterns in Swift](https://github.com/ochococo/Design-Patterns-In-Swift)

### Singletons

Singleton's are a special kind of class where only one instance of the class exists for the current process. They are a convenient way to share data between different parts of an app without creating global variables or having to pass the data around manually, but they should be used sparingly since they often create tighter coupling between classes.

To turn a class into a singleton, you use the following implementation where the function name is prefixed with `shared` plus another word which best describes your class.  For example, if the class is a network or location manager, you would name the function `sharedManager` instead of `sharedInstance`.

```swift
class MyClass {

    // MARK: - Instantiation

    // Naming convention:
    // sharedInstance, sharedManager, sharedController, etc.
    // depending on the class type
    static let sharedInstance = MyClass()

    // This prevents others from using the default '()' initializer for this class.
    private init() {}

    var isReady = true

    // More class code here
}
```

**Explanation**: The static constant `sharedInstance` is run as `dispatch_once` the first time that variable is accessed to make sure the initialization is atomic. This ensures it is thread safe, fast, lazy, and also bridged to ObjC for free. More from [here](http://krakendev.io/blog/the-right-way-to-write-a-singleton).

**Usage**: You would get a reference to that singleton class in another class with the following code:

```swift
// Now you could do
let myClass = MyClass.sharedInstance
let answer = myClass.isReady ? "Yep!" : "Nope!"
print("Are you ready to rock and roll? \(answer)")
```

[Back to top](#table-of-contents)

## Unicode Support

Although I don't recommend this, Swift will compile even if you use emoji's in your code since it offers Unicode support.

More info from Apple [here](https://developer.apple.com/library/ios/documentation/swift/conceptual/Swift_Programming_Language/StringsAndCharacters.html)

[Back to top](#table-of-contents)


## FileIO
### C style FileIO

```swift
let fd = fopen("aFile.txt", "w")
fwrite("Hello Swift!", 12, 1, fd)

let res = fclose(file)
if res != 0 {
    print(strerror(errno))
}

let fd = fopen("aFile.txt", "r")
var array = [Int8](count: 13, repeatedValue: 0)
fread(&array, 12, 1, fd)
fclose(fd)

let str = String.fromCString(array)
print(str) // Hello Swift!
```

### Swift Style FileIO
```swift
let path = Bundle.main.path(forResource:"test", ofType: "txt")
// read
let lines = try? String(contentsOfFile: path!)
                     .split{$0 == "\n"}
                     .map(String.init)

// write
do {
    let lines = self._outlines.joined(separator: "\n")
    try lines.write(to: url, atomically: false, encoding: .utf8)
} catch{}
```
[Back to top](#table-of-contents)
