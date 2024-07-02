---
title: 'Dijkstra'
date: 2020-02-12
categories: ["Algorithm and data structure"]
tags: ["Python", "C++"]
comments: true
---

给出一个有向图，一个起点，一个终点，问起点到终点的最短路径。

## Dijkstra

* dijkstra 算法可以同时求 起点到所有节点的最短路径
* 权值不能为负数
* minDist数组 用来记录 每一个节点距离源点的最小距离。


```python
import heapq

def dijkstra(graph, start):
    # Create a priority queue
    priority_queue = []
    # Initialize the distances with infinity
    distances = {vertex: float('infinity') for vertex in graph}
    # Distance to the start node is zero
    distances[start] = 0
    # Push the start node to the priority queue
    heapq.heappush(priority_queue, (0, start))

    while priority_queue:
        current_distance, current_vertex = heapq.heappop(priority_queue)

        # Nodes can be added multiple times to the priority queue, we only process nodes with the smallest distance
        if current_distance > distances[current_vertex]:
            continue

        for neighbor, weight in graph[current_vertex].items():
            distance = current_distance + weight

            # Only consider this new path if it's better
            if distance < distances[neighbor]:
                distances[neighbor] = distance
                heapq.heappush(priority_queue, (distance, neighbor))

    return distances

# Example graph as an adjacency list
graph = {
    'A': {'B': 1, 'C': 4},
    'B': {'A': 1, 'C': 2, 'D': 5},
    'C': {'A': 4, 'B': 2, 'D': 1},
    'D': {'B': 5, 'C': 1}
}

# Running the algorithm
start_node = 'A'
distances = dijkstra(graph, start_node)
print("Shortest distances from node A:", distances)
```