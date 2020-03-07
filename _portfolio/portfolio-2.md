---
layout: post
title: "YOLOv3 Inference Framework C++"
excerpt: "Torch (C++) implementation of YOLOv3, which works on Windows, Mac, Linux.<br/><img src='/images/yolo_arch1.png'>"
collection: portfolio
date: 2019-08-12
permalink: /portfolio/portfolio-2
---

A good example to learn deep learning model which is coded in C++.

## Performance
test code:
```
yolov3 models/yolov3.cfg models/yolov3.weights images
```

Results:

1. tested with CPU: `Core i9`
    - Windows
    - MAC: average time (682 ms/image).
    
2. tested with GPU: `Tesla V100`
    - Linux: average time (22 ms/image).


## Features
Supports  
- NMS
- Soft NMS
- Weighted NMS

## TODO
- Support training
- ...



## Credits

This repo is created based on the implementations below:  
[weixu000](https://github.com/weixu000/libtorch-yolov3-deepsort)  
[PyTorch-YOLOv3](https://github.com/eriklindernoren/PyTorch-YOLOv3)  
[YOLO_v3_tutorial_from_scratch](https://github.com/ayooshkathuria/YOLO_v3_tutorial_from_scratch)
