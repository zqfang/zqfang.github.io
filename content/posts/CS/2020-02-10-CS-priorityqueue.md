---
title: 'Priority Queque'
date: 2020-02-09
categories: ["Algorithm and data structure"]
tags: ["Python", "C++"]
comments: true
---


## min-heap
```python

import heapq
# create a priorty queque in python, and heapq will ensure the list maintains the heap property.
priority_queue = [] 
# add element. The smallest element will always be at the root, i.e., priority_queue[0].
heapq.heappush(priority_queue, (priority, item))
## Use heapq.heappop to remove and return the smallest element from the priority queue.
item = heapq.heappop(priority_queue)
```
### min-heap example
```python
import heapq

class Task:
    def __init__(self, priority, description):
        self.priority = priority
        self.description = description

    def __lt__(self, other):
        return self.priority < other.priority

priority_queue = []

# Add tasks to the priority queue
heapq.heappush(priority_queue, Task(1, "task1"))
heapq.heappush(priority_queue, Task(3, "task3"))
heapq.heappush(priority_queue, Task(2, "task2"))
heapq.heappush(priority_queue, Task(0, "task0"))

# Pop tasks from the priority queue
while priority_queue:
    task = heapq.heappop(priority_queue)
    print(f"Processing {task.description} with priority {task.priority}")
```


## Max-heap
use the min-heap implementation by storing negative values of the priorities. This way, the largest value (in terms of original priority) is always the smallest (in terms of negated priority).

```python
import heapq

# Initialize an empty heap
max_heap = []

# Function to push an item onto the max-heap
def heappush_max(heap, item):
    heapq.heappush(heap, -item)

# Function to pop the largest item from the max-heap
def heappop_max(heap):
    return -heapq.heappop(heap)

# Function to peek at the largest item in the max-heap
def heap_max(heap):
    return -heap[0]

# Push items onto the max-heap
heappush_max(max_heap, 10)
heappush_max(max_heap, 30)
heappush_max(max_heap, 20)
heappush_max(max_heap, 40)

# Peek at the largest item
print("Max item:", heap_max(max_heap))  # Output: Max item: 40

# Pop items from the max-heap
print("Popped item:", heappop_max(max_heap))  # Output: Popped item: 40
print("Popped item:", heappop_max(max_heap))  # Output: Popped item: 30
print("Popped item:", heappop_max(max_heap))  # Output: Popped item: 20
print("Popped item:", heappop_max(max_heap))  # Output: Popped item: 10

```

## Priority Queue in C++

### Max-Heap
```cpp
#include <iostream>
#include <queue>

int main() {
    // Create a max-heap priority queue
    std::priority_queue<int> pq;

    // Push elements into the priority queue
    pq.push(10);
    pq.push(20);
    pq.push(15);

    // Display the top element
    std::cout << "Top element: " << pq.top() << std::endl;  // Output: 20

    // Pop the top element
    pq.pop();

    // Display the top element again
    std::cout << "Top element after pop: " << pq.top() << std::endl;  // Output: 15

    return 0;
}
```

### Min-Heap
```cpp
#include <iostream>
#include <queue>
#include <vector>

struct Task {
    int priority;
    std::string description;

    Task(int p, std::string d) : priority(p), description(d) {}
};

struct CompareTask {
    bool operator()(Task const& t1, Task const& t2) {
        // Higher priority comes first
        return t1.priority < t2.priority; // max-heap
        // return t1.priority > t2.priority;; // Invert the comparison to create a min-heap
    }
};

int main() {
    // Create a priority queue for Tasks
    std::priority_queue<Task, std::vector<Task>, CompareTask> taskQueue;

    // Push tasks into the priority queue
    taskQueue.push(Task(1, "Low priority task"));
    taskQueue.push(Task(3, "High priority task"));
    taskQueue.push(Task(2, "Medium priority task"));

    // Process tasks by priority
    while (!taskQueue.empty()) {
        Task t = taskQueue.top();
        std::cout << "Processing task: " << t.description << " with priority " << t.priority << std::endl;
        taskQueue.pop();
    }

    return 0;
}
```