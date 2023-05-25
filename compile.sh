#!/bin/bash

cmake -S . -B build
cmake --build build --target supra
#cmake --build build --target sym
