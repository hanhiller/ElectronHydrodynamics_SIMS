#!/bin/bash

echo "Lowest Temperature, Highest Tempature, Spacing:"
read lowT highT dT

echo "Particle Density:"
read partDensity

echo "Probe Tip Location (x y):"
read x y

echo "Number of CPUs:"
read nCPU

echo "Time Steps:"
read timeSteps

for T in $(seq ${lowT} ${dT} ${highT})
do
	file="./fullCircle_T${T}_Probe${x},${y}.txt"

	echo "[DEFAULT]" > $file
	echo "" >> $file
	echo "Temperature = $T" >> $file
	echo "Particle Density = $partDensity" >> $file
	echo "probeCenterX = $x" >> $file
	echo "probeCenterY = $y" >> $file
	echo "Number of CPUs = $nCPU" >> $file
	echo "time steps = $timeSteps" >> $file
	echo "" >> $file
	echo "Constriction Width = 0.3" >> $file
	echo "Scattering Probability = 0" >> $file
	echo "Source Drain Ratio = 1.2" >> $file
	echo "save interval = 10000" >> $file
	echo "diameter = 10" >> $file
	echo "injector height and width = 1.5,0.6" >> $file
	echo "" >> $file
	echo "base output path = ./" >> $file
    	echo "initCondFile = fullCircle_T${T}_Probe${x},${y}" >> $file
done
