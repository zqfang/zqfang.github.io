# C++ Notes: Concurrency


Learn C++11 thread library.  Code snippets from [Concurrent Programming with C++11](https://www.youtube.com/playlist?list=PL5jc9xFGsL8E12so1wlMS0r0hTQoJL74M)

## Process vs. Threads
 
![concurrency](/images/cpp-concurrency.png)

## Usage Summary

A short summary of thread library in STL

1. thread and async
```cpp
/* thread */
std::thread t1(factorial, 6); // create a new thread
std::this_thread::sleep_for(chrono::milliseconds(3));
chrono::steady_clock::time_point tp = chrono::steady_clock::now() + chrono::microseconds(4);
std::this_thread::sleep_until(tp);

/* async() */
std::future<int> fu = async(factorial, 6); // create a new thread
```

2. mutex
```cpp
/* Mutex */
std::mutex mu;
std::lock_guard<mutex> locker(mu);
std::unique_lock<mutex> ulocker(mu);
ulocker.try_lock();
ulocker.try_lock_for(chrono::nanoseconds(500));
ulocker.try_lock_until(tp);
```

3. condition variable
```cpp
/* Condition Variable */
std:condition_variable cond;
cond.wait_for(ulocker, chrono::microseconds(2));
cond.wait_until(ulocker, tp);
```

4. future and promise
```cpp
/* Future and Promise */
std::promise<int> p; 
std::future<int> f = p.get_future();
f.get();
f.wait();
f.wait_for(chrono::milliseconds(2));
f.wait_until(tp);
```

1. Packaged task
```cpp
 /* Packaged Task */
 std::packaged_task<int(int)> t(factorial);
 std::future<int> fu2 = t.get_future();
 t(6);
```


## Cases

### thread

two way to create new thread
```cpp
std::thread t1(func);
std::async(std::launch::async, func);
```

exmample of thread

```cpp
#incldue <thread>
void function1() {
    std::cout <<"hello"<<std::endl;
}

std::tread t1(function1); // t1 start running
// t1.join();  // main thread wait for t1 to finish
t1.detach(); // t1 will freely on its own -- deamon process``

/// once detach, forever detach. 

if (t1.joinable())
  t1.join();  // if detached, this line crashed.
```

### thread managment

```cpp
class Fctor {
public:
    void operator()(std::string & s) {
        std::cout<<"this is thread"<<std::endl;
    }
};
std::string = "string int"
std::thread t1((Fctor()), s); // alway pass by value
std::thread t2((Fctor()), std::ref(s)); // pass by ref

std::thread t3((Fctor()), std::move(s)); // move s from main to thread
std::thread t4 = std::move(t3); // thread could not be copy, only move

try {
    std::cout<<"this is main"<<std::endl;
} catch (...) {
    t1.join();
    t2.jion();
    throw;
}

t1.join();
t2.join();
t4.join();
```

if oversubscription, limit threads with maximum cpu cores
```cpp
std::thread::hardware_concurrency(); // indication, number of cpu cores
```

### race condition and mutex 

```cpp
#include <thread>
#include <mutex>


std::mutex mu;
void shared_print(std::string msg, int id) {
    std::lock_guard<std::mutex> guard(mu); // RAII
    //mu.lock(); // +> safeguard no two threads using cout at the same time, 
    std::cout<<msg << id << std::endl; 
    //mu.unlock(); // if error thrown between .lock() and .unlock(), generate zombie process
}

void function_1() {
    for(int i=0; i < 100; i++)
        shared_print("from functino t1", i);
}

int main() {
    std::thread t1(function_1);
    for (int i=0; i < 100; i++)
        shared_print("from main", i);
    t1.join()
}

```

more practical example

```cpp
#include <thread>
#include <mutex>

class LogFile {
    std::mutex mu;
    std::ofstream f;
public:
    LogFile() {
        f.open("log.txt");
    }
    void shared_print(std::string msg, int id) {
        std::lock_guard<std::mutex> guard(mu); // RAII
        f << msg << id << std::endl; 
    }
    // never return f to outside wworld
    // never pass f as an argument for user
};
void function_1(LogFile& log) {
    for(int i=0; i < 100; i++)
        log.shared_print("from functino t1", i);
}

int main() {
    LogFile log;
    std::thread t1(function_1, std::ref(log));
    for (int i=0; i < 100; i++)
        log.shared_print("from main", i);
    t1.join()
}

```

### Advoid Deadlock 

1. prefer locking single mutex
2. Advoid locking a mutex and then calling a user provded function
3. use std::lock() to lock more than one mutex
4. lock the mutexs in same order.


### unique_lock and lazy initialization

```cpp

class LogFile {
    std::mutex mu;
    std::ofstream f;
    std::once_flag flag;
public:
    LogFile() {
        f.open("log.txt");
    }
    void shared_print(std::string msg, int id) {
        // if you need to check whether a file is open in each call, use once_flag to rescue
        // std::call_once(flag, [&](){ f.open("log.txt");}) // file only open once

        std::unique_lock<std::mutex> locker(mu, std::defer_lock); // note here
        // do something else
        locker.lock();
        f << msg << id << std::endl; 
        locker.unlock();

        // call again
        locker.lock();
        // do something ...

        locker.unlock();

        std::unique_lock<std::mutex> locker2 = std::move(locker); // could change ownership
    }  

};
```

### condition variable

Condition variable is to synchronize the execution order of threads

```cpp
std::condition_variable cond;

// usage 1
std::unique_lock<std::mutex> locker(mu);
cond.wait(locker); // spurious wake

cond.wait(locker, [](){return !q.empty();}) //

cond.notify_one(); // notify one waiting thread
cond.notify_all(); // 

```

### future and promise

Future and promise provide a convenience way to communicate between threads.

e.g. return value to main thread

```cpp
int factorial(int n) {
    int res = 1;
    return res+ n;
}
 
int x;

std::future<int>  fu = std::async(function, 4); // future, get something in future
x = fu.get();
// fu.get(); //crash

std::future<int> fu2 = std::async(std::launch::deferred, factorial, 4); // means not excuate unitl call .get()
x = fu2.get(); // only excuate fu2 when called get

std::future<int> fu3 = std::async(std::launch::async | std::launch::deferred , factorial, 4); // create new thread by calling async or not 
x = fu3.get(); // only excuate fu2 when called get
```

usage of promise
```cpp
int factorial(std::future<int> &f) {
    int res = 1;
    int N = f.get(); // note here
    return res + N;
}
int x;
std::promise<int> p;
std::future<int> f = p.get_future();
std::future<int> fu4 = std::async(std::launch::async, factorial, std::ref(f));
// do something else  ...
//// if p not set, throw error
// p.set_exception(std::make_exception_ptr)(std::runtime_error("To err is human"));
// set p
p.set_value(4);

// get from child 
x = fu4.get();
```

`shared_future` for multi-threads
```cpp
int factorial(std::shared_future<int> &f) {
    int res = 1;
    int N = f.get(); // note here
    return res + N;
}
int x;
std::promise<int> p;
std::future<int> f = p.get_future();
std::shared_future<int> sf = f.shared();

std::future<int> fu5 = std::async(std::launch::async, factorial, sf);
std::future<int> fu6 = std::async(std::launch::async, factorial, sf);
std::future<int> fu7 = std::async(std::launch::async, factorial, sf);
p.set_value(4);
// get from child 
x = fu4.get();
```



### using callable object

```cpp
class A {
public:
    void f(int, char) {}
    long g(double x) {return 0;}
    int operator()(int n) {return 0;}
};

A a;
std::thread t1(a, 6); // copy of a() in a different thread
std::thread t2(std::ref(a), 6) // a() in a different thread
std::thread t3(A(), 6); // temp A
std::thread t4([](int x){return x*x;}, 6);

std::thread t5(&A::f, a, 6, 'w'); // copy of a.f(6,'w') in a different thread
// these feature could be used in
// std::bind, std::async
```

### packagee tasks

packaged_task provides a way to implement a task pool.  It can conveniently convey the returned value from a task to a different thread

```cpp
std::thread t1(factorial, 6); // could pass args
std::packaged_task<int(int)> t(factorial); // could not pass additional args
std::packaged_task<int()> t2(std::bind(factorial, 6)); // now we could pass args using std::bind
// ...
t(6); // in a different context, t alwaly return void, so
int x = t.get_future().get(); // get value

// call t2 by
t2();
```

```cpp
int factorial(int N) {
    int res =1;
    for (int i=N; i > 1; i --)
        res *= i;
    return res;
}

std::deque<std::packaged_task<int()>> task_q;
std::mutex mu; 
std::condition_variable cond;

void thread_1() {
    std::packaged_task<int()> t;
    {  
        //std::lock_guard<std::mutex> locker(mu); // advoid data race
        std::unique_lock<std::mutex> locker(mu);
        cond.wait(locker, [](){return !taks_q.empty();})
        t = std::move(task_q.front());
        task_q.pop_front();
    }
    t();
}

int main() {
    std::thread t1(thread_1);  // so, task_q run in t1;
    std::packaged_task<int()> t(std::bind(factorical, 6));
    std::future<int> fu = t.get_future(); // get returned value to main thread

    {
        std::lock_guard<std::mutex> locker(mu);
        task_q.push_bask(std::move(t));
    }
    std::cout<<fu.get();
    t1.join();
    return 0;
}

```


Summary: 3 method to get a future

- promise::get_future()
- packaged_task::get_future()
- async() returns a future



