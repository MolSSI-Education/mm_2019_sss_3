#!/bin/bash

g++ -fPIC --shared -Wall box.cpp export.cpp \
    -I ${CONDA_PREFIX}/include \
    -I ${CONDA_PREFIX}/include/python3.7m \
    -I ${CONDA_PREFIX}/include/eigen3 \
    -I ${CONDA_PREFIX}/lib/python3.7m -lpython3.7m \
    -o box_cpp.so
