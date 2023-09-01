#!/bin/bash

IFS='	'
re='^[0-9]+([.][0-9]+)?$'

while read -ra line
do
	temp=$(echo ${line[@]:1})
	nn=${line[0]}_s
	declare $nn="$temp"
done < allparams.txt

IFS=' '
read -ra Epss <<< "$eps_s"
read -ra Ds <<< "$d_s"
read -ra Ls <<< "$l_s"
read -ra Rhos <<< "$rho_s"
read -ra Vs <<< "$v_s"
read -ra Kagars <<< "$Kagar_s"
read -ra Kstiffs <<< "$Kstiff_s"
read -ra Revs <<< "$rev_s"
read -ra max_n <<< "$max_n_s"


for eps in "${Epss[@]}"
do
	for d in "${Ds[@]}"
	do
		for l in "${Ls[@]}"
		do
			for rho in "${Rhos[@]}"
			do
				for v in "${Vs[@]}"
				do
					for Kagar in "${Kagars[@]}"
					do
						for Kstiff in "${Kstiffs[@]}"
						do
							for rev in "${Revs[@]}"
							do
								./run_slurm3Dallexp.sh $eps $d $l $rho $v $Kagar $Kstiff $rev $max_n
							done
						done
					done
				done
			done
		done
	done
done
