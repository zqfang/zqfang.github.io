# YOLOv3 Inference Framework in C++

<br/><img src='/images/yolo_arch1.png'>
My journey to object detection began with YOLOv3. I think it's really a good starting point for someone like me without computer vison background to understand what's going on behind the scence. While learning object detection, I made a simple modified C++ version (with LibTorch) based on others' work [^1] [^2] [^3].  

See code [here](https://github.com/zqfang/YOLOv3_CPP)

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

`YOLOv3_CPP` repo is created based on the implementations below:  
 
[^1]: [weixu000](https://github.com/weixu000/libtorch-yolov3-deepsort).
[^2]: [PyTorch-YOLOv3](https://github.com/eriklindernoren/PyTorch-YOLOv3).
[^3]: [YOLO_v3_tutorial_from_scratch](https://github.com/ayooshkathuria/YOLO_v3_tutorial_from_scratch).

