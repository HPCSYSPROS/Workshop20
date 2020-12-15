caffe
=====

This caffe benchmark is training a image classifier using the ResNet-50
architecture on a 100-category ImageNet dataset.

Input
=====

Input is too large to include here.

The online dataset has 1K categories. For benchmarking we used the first 100
categories.  The data can be obtained from:

http://www.image-net.org/challenges/LSVRC/2012/nnoupb/ILSVRC2012_img_train.tar

and

http://www.image-net.org/challenges/LSVRC/2012/nnoupb/ILSVRC2012_img_val.tar

Code modifications
==================

There are no code modifications.

Optimisations
=============

We use the -xCOMMON-AVX512 flag to build it.

Third-party dependencies
========================

Caffe depends on the following third-party dependencies:

1. Boost            1.65.1
2. MKL              17.0.4
3. gflags           2.2.0
4. glog             0.3.4
5. HDF5             1.8.16
6. leveldb          1.19
7. lmdb             2.8
8. mlsl             l_mlsl_2017.1.016
9. opencv           2.4.13
10. protocol-buffer 3.0.0
11. snappy          1.1.3

Timing
======

Timing output is in output/log-resnet50-1024skx-mkl2017-1.log
