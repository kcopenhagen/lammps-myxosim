#!/bin/bash

eps=$1
epspath=$(printf "eps%03d/" $eps)
d=$2
dpath=$(printf "d%04d/" `echo "$d*100/1" | bc`)
l=$3
lpath=$(printf "l%04d/" `echo "$l*100/1" | bc`)
rho=$4
rhopath=$(printf "rho%04d/" `echo "$rho*100/1" | bc`)
v=$5
vpath=$(printf "v%04d/" `echo "$v*100/1" | bc`)
Kagar=$6
Kagarpath=$(printf "Kagar%03d/" $Kagar)
Kstiff=$7
Kstiffpath=$(printf "Kstiff%03d/" $Kstiff)
rev=$8
revpath=$(printf "rev%03d/" $rev)
max_n=$9

sed s/@e/$eps/ in.spfr3Dall > in.spfr3Dt
sed s/@d/$d/ in.spfr3Dt > in.spfr3Du
sed s/@d/$d/ Initial_all.m > Initial_allt.m 
req=`echo "$d*0.6" | bc`
sed s/@r/$req/ in.spfr3Du > in.spfr3Dt
s=`echo "$d/(e(l(2)*(1/6.0)))" | bc -l`
sed s/@s/$s/ in.spfr3Dt > in.spfr3Du
sed s/@L/$l/ Initial_allt.m > Initial_allu.m
sed s/@n/$rho/ Initial_allu.m > Initial_allt.m
sed s/@v/$v/ in.spfr3Du > in.spfr3Dt
sed s/@A/$Kagar/ in.spfr3Dt > in.spfr3Du
w=`echo "$Kagar*0.3" | bc`
sed s/@w/$w/ in.spfr3Du > in.spfr3Dt
sed s/@K/$Kstiff/ in.spfr3Dt > in.spfr3Du
sed s/@R/$rev/ in.spfr3Du > in.spfr3Dt

fpath1="/scratch/gpfs/kc32/lammps-myxo-data/PoissonRevs/"

f=$fpath1$epspath$dpath$lpath$rhopath$vpath$Kagarpath$Kstiffpath$revpath
i=0
j=$(printf "run_%03d" $i)

fname=$f$j
if [ -d "$fname" ]
then
while [ -d "$fname" ]
do
	let i=$i+1
	if [ $i -gt $max_n ]
	then
		echo "Maximum runs reached"
		exit 1
	fi
	j=$(printf "run_%03d" $i)
	fname=$f$j
done
fi

echo $fname

mkdir -p $fname/cells
cp run.slurm3Dall_exp $fname
cp in.spfr3Dt $fname/in.spfr3Dall_exp
cp Initial_allt.m $fname/Initial_all.m
cd $fname

sbatch run.slurm3Dall_exp
