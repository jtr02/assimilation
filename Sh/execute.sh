#!/bin/sh

# setting for calculation
# model
nx=40
nt=15000
dt=0.01d0
forcing=8.0d0
# option
imem=10
pmem=100
# unit
ufc=51
ueg=52

# path to directory
pdir=/home/jtr/Assimilation
wdir=${pdir}/Work
sdir=${pdir}/Module/Src
edir=${pdir}/Module/Exe
ddir=/mnt/d/Data/Assimilation

# compile
gfortran ${sdir}/lorenz96.f90 -o ${edir}/lorenz96.out

# execution start
rm -r ${wdir}
mkdir ${wdir}
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
/
&unit
ufc=${ufc}
ueg=${ueg}
/
EOF

${edir}/lorenz96.out

mv out${ufc}.grd ${ddir}/Result/
mv out${ufc}.ctl ${ddir}/Result/fc.ctl
mv out${ueg}.grd ${ddir}/Result/
mv out${ueg}.ctl ${ddir}/Result/eg.ctl

# execution end
cd ${pdir}
