---
layout: page
title: About Me
permalink: /about/
search_exclude: true
---

My name is Zhuoqing Fang  
I'm now a postdoc working on computational genomics at Stanford.  
I like traveling, hiking, sleeping.  


### Education

* Ph.D., Developmental Biology, Shanghai Institutes for Biological Sciences, University of Chinese Academy of Science
  * Date: June 2019
* Visiting Graduate, Computational Biology, University of California Los Angeles
  * Date: April to July 2015
* B.S., Life Science, China Agricultural University
  * Date: June 2012

### Work experience

* Nov 2019 - : Postdoc
  * Stanford University
  * Machine learning in genomics
  * Advisor: Prof. Gary Peltz

* May 2019 - Oct 2019: Deep Learning Researcher
  * Fosun Aitrox
  * Duties included: 
    - Build models for digital pathology images.
      - Pytorch, and Tensorflow
      - Multi-Instance Learning for digital pathology images
      - Supervised learning: Faster R-CNN, YOLOv3, SSD
    - Deploy models onto microscopy
      - Build inference framework (C++) for cloud and local computers

  
### Skills

* Programming 
  * Python, R, C++


* Deep learning
  * Computer vision
    - Image registration
    - Object detection, segmentation, classification
  * Framework
    - Pytorch
    - Tensorflow 2.0
    - S4TF: Swift for Tensorflow. The next generation of DL framework I'd like to try.
 
* Computational genomics
  * Single Cell RNA-seq
  * ChIP-seq
  * Other bioinfo stuffs

* Stem cell biology
  * CRISPR/cas9 gene editing
  * Molecular bio skills


### Publications

<ul>{% for post in site.publications reversed %}
  {% include archive-single-cv.html %}
{% endfor %}</ul>