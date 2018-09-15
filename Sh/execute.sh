#!/bin/sh

pdir=/home/jtr/Assimilation
ddir=/mnt/d/Data/Assimilation

gfortran ${pdir}/Module/Src/lorenz96.f90 -o ${pdir}/Module/Exe/lorenz96.out

cd ${pdir}/Work

${pdir}/Module/Exe/lorenz96.out
#rm ${pdir}/Module/Exe/lorenz96.out

mv output.grd ${ddir}/Result/
mv out.ctl ${ddir}/Result/

cd ${pdir}
