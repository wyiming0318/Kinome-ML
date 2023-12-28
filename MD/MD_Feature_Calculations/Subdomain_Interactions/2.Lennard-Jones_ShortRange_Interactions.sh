
# 8. Calculate desired interations and output to .xvg
echo "###########################################" >> ../calc.log
for d in */
do {
cd ./$d
name="$(basename "$PWD")"
if [ -f "LJ.xvg" ]
then
	echo "Calculating interaction data $name ... Already exist!" >> ../../calc.log
fi
if [ ! -f "LJ.xvg" ]
then
gmx energy -f ener.edr -s energy_calc.tpr -o LJ.xvg -b 50000 -e 100000 <<EOD
19
23
27
31
39
43
47
55
59
67
EOD
	 if [ -f "LJ.xvg" ]
	 then
	 	echo "Calculating interaction data $name ... Done!" >> ../../calc.log
	 fi
	 if [ ! -f "LJ.xvg" ]
	 then
	 	echo "Calculating interaction data $name ... Failed" >> ../../calc.log
	 fi
cp ./LJ.xvg ./LJ.csv
fi
cd ../
} done
echo "###########################################" >> ../calc.log
echo >> ../calc.log
echo >> ../calc.log
