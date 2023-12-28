for d in */

do {
	
cd ./$d
	
name="$(basename "$PWD")" # get mutation ID
name=$(basename "$PWD")
name_SASA="$(basename "$PWD")_SASA"

gmx sasa -f *.trr -s *.tpr -n $name.ndx -o $name_SASA.xvg -b 50 -tu ns -output 1 2 3 4 5 -surface 0

cd ../
}
done
