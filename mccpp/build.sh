#!/bin/bash

g++ -fPIC --shared -std=c++11 -Wall pairwise.cpp  export.cpp \
-I${CONDA_PREFIX}/include \
-I${CONDA_PREFIX}/include/python3.6m \
-I${CONDA_PREFIX}/include/eigen3 \
-L${CONDA_PREFIX}/lib -lpython3.6m \
-o pairwise.so
