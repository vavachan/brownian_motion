#!/bin/sh
for i in $(ls -d *Mean_square*)
do
gnuplot << EOF
set print "slopes.txt" append
f(x)=m*x+c
fit [6:10][100:300] f(x) '$i' u 1:2 via m,c
print m
EOF
done
