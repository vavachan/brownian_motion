#!/bin/bash
for m in  4.0 5.0 1.0 2.0 3.0 
do
echo $m
for n in {1..5} 
do
echo $n
python manleap.py $n $m
done
done
