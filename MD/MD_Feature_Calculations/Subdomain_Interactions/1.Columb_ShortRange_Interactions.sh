# define variables
dir_data="Coul"

# 1. Create a folder to hold gmx calculation data for simulations
echo "###########################################" >> calc.log
echo "Creating data folder '$dir_data'..." >> calc.log

if [ -d "$dir_data" ]; then
	echo "Folder already exist." >> calc.log
else
	mkdir $dir_data
	echo "Data folder created!" >> calc.log
fi
echo "###########################################" >> calc.log
echo >> calc.log
echo >> calc.log
# 2. Create folders for all mutants and WT
echo "###########################################" >> calc.log
echo "Creating subfolders under '$dir_data'..." >> calc.log
for d in */
do {

cd ./$d
name="$(basename "$PWD")"
DIR="../$dir_data/$name"
if [ $name != $dir_data ] && [ -d "$DIR" ] 
then
	echo "$name ... Already exist!" >> ../calc.log
fi
if [ $name != $dir_data ] && [ ! -d "$DIR" ] 
then
	mkdir $DIR
	echo "$name ... Done!" >> ../calc.log
fi

cd ../
} done
echo "###########################################" >> calc.log
echo >> calc.log
echo >> calc.log
# 3. Copy md2.mdp from simulation folder to calculation folder
echo "###########################################" >> calc.log
cd ./$dir_data
for d in */
do {
cd ./$d
name="$(basename "$PWD")"
if [ -f "md2.mdp" ] 
then
	echo "Copying md2.mdp from simulation folder to $name ... Already exist!" >> ../../calc.log
fi
if [ ! -f "md2.mdp" ] 
then
	cp ../../$name/md2.mdp ./
	echo "Copying md2.mdp from simulation folder to $name ... Done!" >> ../../calc.log
fi
cd ../
} done
echo "###########################################" >> ../calc.log
echo >> ../calc.log
echo >> ../calc.log
# 4. Generate index file for energy groups
echo "###########################################"  >> ../calc.log
for d in */
do {
cd ./$d
name="$(basename "$PWD")"
if [ -f "energy_groups.ndx" ] 
then
	echo "Writing index file for energy groups of $name ... Already exist!"  >> ../../calc.log
fi
if [ ! -f "energy_groups.ndx" ]
then

gmx make_ndx -f ../../$name/proteinOnly.gro -o energy_groups.ndx <<EOD
1 & r 727-732
1 & r 760-778
1 & r 839-850
1 & r 861-884
name 10 PL
name 11 ACH
name 12 CL
name 13 AL
q
EOD
	if [ -f "energy_groups.ndx" ]
	then
		echo "Writing index file for energy groups of $name ... Done!" >> ../../calc.log
	fi
	if [ ! -f "energy_groups.ndx" ]
	then
		echo "Writing index file for energy groups of $name ... Failed!" >> ../../calc.log
	fi
fi
cd ../
} done
echo "For failed ones, please ensure the name of the .gro file in simulation folder is correct" >> ../calc.log
echo "###########################################" >> ../calc.log
echo >> ../calc.log
echo >> ../calc.log
# 5. Add energy groups to the new .mdp file
echo "###########################################"  >> ../calc.log
for d in */
do {
cd ./$d
name="$(basename "$PWD")"
echo $name
if [ -f "energy_calc.mdp" ]
then
	echo "Writing energy groups to 'energy_calc.mdp' $name ... Already exist!" >> ../../calc.log
fi
if [ ! -f "energy_calc.mdp" ]
then
	cp ./md2.mdp ./energy_calc.mdp
	echo "; ENERGY GROUPS " >> energy_calc.mdp
	echo "energygrps = PL ACH CL AL" >> energy_calc.mdp
	echo "verlet-buffer-tolerance = -1" >> energy_calc.mdp
	echo "Writing energy groups to 'energy_calc.mdp' $name ... Done!" >> ../../calc.log
	sed -i 's/rlist                    = 1/rlist                    = 4.4/' energy_calc.mdp
	sed -i 's/rcoulomb                 = 1/rcoulomb                 = 4.4/' energy_calc.mdp
	sed -i 's/rvdw                     = 1/rvdw                     = 4.4/' energy_calc.mdp
	sed -i 's/coulombtype              = PME/coulombtype              = Cut-off/' energy_calc.mdp
	sed -i 's/tc-grps                  = Protein Non-Protein/tc-grps                  = Protein/' energy_calc.mdp
	sed -i 's/tau-t                    = 0.2 0.2/tau-t                    = 0.2/' energy_calc.mdp
	sed -i 's/ref-t                    = 300 300/ref-t                    = 300/' energy_calc.mdp
	sed -i 's/nstlist                  = 10/nstlist                  = 1/' energy_calc.mdp
	
fi
cd ../
} done
echo "###########################################" >> ../calc.log
echo >> ../calc.log
echo >> ../calc.log
# 6. Use GROMACS 'gmx grompp' to generate new .tpr file
echo "###########################################" >> ../calc.log
for d in */
do {
cd ./$d
name="$(basename "$PWD")"
if [ -f "energy_calc.tpr" ]
then
	echo "Generating new .tpr file for energy calculation $name ... Already exist!" >> ../../calc.log
fi
if [ ! -f "energy_calc.tpr" ]
then
gmx grompp -maxwarn 1 -f energy_calc.mdp -c ../../$name/proteinOnly.gro -p ../../$name/proteinOnly.top -n 		energy_groups.ndx -o energy_calc.tpr
	if [ -f "energy_calc.tpr" ]
	then
		echo "Generating new .tpr file for energy calculation $name ... Done!" >> ../../calc.log
	fi
	if [ ! -f "energy_calc.tpr" ]
	then
		echo "Generating new .tpr file for energy calculation $name ... Failed!" >> ../../calc.log
	fi 
fi

cd ../
} done
echo "###########################################" >> ../calc.log
echo >> ../calc.log
echo >> ../calc.log
# 7. Apply GROMACS "gmx mdrun -rerun" to get new .edr file
echo "###########################################" >> ../calc.log
for d in */
do {
cd ./$d
name="$(basename "$PWD")"
if [ -f "ener.edr" ]
then
	echo "Generating new .edr file $name ... Already exist!" >> ../../calc.log
fi
if [ ! -f "ener.edr" ]
then
gmx mdrun -rerun ../../$name/proteinOnly.xtc -s energy_calc.tpr
	if [ -f "ener.edr" ]
	then
		echo "Generating new .edr file $name ... Done!" >> ../../calc.log
	fi
	if [ ! -f "ener.edr" ]
	then
		echo "Generating new .edr file $name ... Failed!" >> ../../calc.log
	fi
fi
rm ./traj.trr
cd ../
} done
echo "###########################################" >> ../calc.log
echo >> ../calc.log
echo >> ../calc.log
# 8. Calculate desired interations and output to .xvg
echo "###########################################" >> ../calc.log
for d in */
do {
cd ./$d
name="$(basename "$PWD")"
if [ -f "energy.xvg" ]
then
	echo "Calculating interaction data $name ... Already exist!" >> ../../calc.log
fi
if [ ! -f "energy.xvg" ]
then
gmx energy -f ener.edr -s energy_calc.tpr -b 50000 -e 100000 <<EOD
17
19
21
23
25
27
29
31
37
39
41
43
45
47
53
55
57
59
65
67
EOD
	 if [ -f "energy.xvg" ]
	 then
	 	echo "Calculating interaction data $name ... Done!" >> ../../calc.log
	 fi
	 if [ ! -f "energy.xvg" ]
	 then
	 	echo "Calculating interaction data $name ... Failed" >> ../../calc.log
	 fi
cp ./energy.xvg ./energy.csv
fi
cd ../
} done
echo "###########################################" >> ../calc.log
echo >> ../calc.log
echo >> ../calc.log
