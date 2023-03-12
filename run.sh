#!/bin/bash

exe=./exe

if [ -f "$exe" ]; then
    rm $exe
fi

make

if [ -f "$exe" ]; then
    time $exe
fi