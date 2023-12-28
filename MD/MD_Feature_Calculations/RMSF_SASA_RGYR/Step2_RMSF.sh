for d in */

do {
	
cd ./$d
	
name="$(basename "$PWD")" # get mutation ID
name=$(basename "$PWD")
name_RMSF1="$(basename "$PWD")_RMSF_P-loop"
name_RMSF2="$(basename "$PWD")_RMSF_aC-Helix"
name_RMSF3="$(basename "$PWD")_RMSF_C-loop"
name_RMSF4="$(basename "$PWD")_RMSF_A-loop"
name_RMSF5="$(basename "$PWD")_RMSF_rest"

gmx rmsf -f *.trr -s *.tpr -n $name.ndx -o $name_RMSF1.xvg -b 50000 <<EOD
6
EOD
gmx rmsf -f *.trr -s *.tpr -n $name.ndx -o $name_RMSF2.xvg -b 50000 <<EOD
7
EOD
gmx rmsf -f *.trr -s *.tpr -n $name.ndx -o $name_RMSF3.xvg -b 50000 <<EOD
8
EOD
gmx rmsf -f *.trr -s *.tpr -n $name.ndx -o $name_RMSF4.xvg -b 50000 <<EOD
9
EOD
gmx rmsf -f *.trr -s *.tpr -n $name.ndx -o $name_RMSF5.xvg -b 50000 <<EOD
10
EOD
cd ../
}
done
