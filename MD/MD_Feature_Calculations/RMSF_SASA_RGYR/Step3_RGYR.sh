for d in */

do {
	
cd ./$d
	
name="$(basename "$PWD")" # get mutation ID
name=$(basename "$PWD")
name_RGYR1="$(basename "$PWD")_RGYR_P-loop"
name_RGYR2="$(basename "$PWD")_RGYR_aC-Helix"
name_RGYR3="$(basename "$PWD")_RGYR_C-loop"
name_RGYR4="$(basename "$PWD")_RGYR_A-loop"
name_RGYR5="$(basename "$PWD")_RGYR_rest"

gmx gyrate -f *.trr -s *.tpr -n $name.ndx -o $name_RGYR1.xvg -b 50000 <<EOD
1
EOD
gmx gyrate -f *.trr -s *.tpr -n $name.ndx -o $name_RGYR2.xvg -b 50000 <<EOD
2
EOD
gmx gyrate -f *.trr -s *.tpr -n $name.ndx -o $name_RGYR3.xvg -b 50000 <<EOD
3
EOD
gmx gyrate -f *.trr -s *.tpr -n $name.ndx -o $name_RGYR4.xvg -b 50000 <<EOD
4
EOD
gmx gyrate -f *.trr -s *.tpr -n $name.ndx -o $name_RGYR5.xvg -b 50000 <<EOD
5
EOD

cd ../
}
done
