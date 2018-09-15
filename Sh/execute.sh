#!/bin/sh

gfortran ../Module/Src/lorenz96.f90 -o lorenz96.out

./lorenz96.out

rm lorenz96.out

mv output.grd /mnt/d/Data/Assimilation/Result/
mv out.ctl /mnt/d/Data/Assimilation/Result/
