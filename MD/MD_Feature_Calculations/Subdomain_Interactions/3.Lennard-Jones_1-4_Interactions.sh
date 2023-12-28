
# 8. Calculate desired interations and output to .xvg
echo "###########################################" >> ../calc.log
for d in */
do {
cd ./$d
name="$(basename "$PWD")"
if [ -f "LJ14.xvg" ]
then
	echo "Calculating interaction data $name ... Already exist!" >> ../../calc.log
fi
if [ ! -f "LJ14.xvg" ]
then
gmx energy -f ener.edr -s energy_calc.tpr -o LJ14.xvg -b 50000 -e 100000 <<EOD
21
25
29
33
41
45
49
57
61
69
EOD
	 if [ -f "LJ14.xvg" ]
	 then
	 	echo "Calculating interaction data $name ... Done!" >> ../../calc.log
	 fi
	 if [ ! -f "LJ14.xvg" ]
	 then
	 	echo "Calculating interaction data $name ... Failed" >> ../../calc.log
	 fi
cp ./LJ14.xvg ./LJ14.csv
fi
cd ../
} done
echo "###########################################" >> ../calc.log
echo >> ../calc.log
echo >> ../calc.log
