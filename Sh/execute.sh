#!/bin/sh

# setting for calculation
# model
nx=40
nt=15000
dt=0.01d0
forcing=16.0d0
# option
imem=1
pmem=1
ptb=0.01d0

# path to directory
pdir=/home/jtr/Assimilation
wdir=${pdir}/Work
sdir=${pdir}/Module/Src
edir=${pdir}/Module/Exe
ddir=/mnt/d/Data/Assimilation

# compile
gfortran ${sdir}/lorenz96.f90 -o ${edir}/lorenz96.out

# execution start
cd ${wdir}

# make namelist
cat << EOF > namelist
&model
nx=${nx}
nt=${nt}
dt=${dt}
forcing=${forcing}
/
&option
imem=${imem}
pmem=${pmem}
ptb=${ptb}
/
EOF

${edir}/lorenz96.out

mv output.grd ${ddir}/Result/
mv out.ctl ${ddir}/Result/

# execution end
cd ${pdir}
