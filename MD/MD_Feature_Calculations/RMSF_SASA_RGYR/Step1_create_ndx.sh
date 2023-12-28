for d in */

do {
	
cd ./$d
	
name="$(basename "$PWD")" # get mutation ID
gmx make_ndx -f *.tpr -o $name.ndx <<EOD
3
1
del 0-18
r 727-732
r 760-778
r 839-850
r 861-884
1 & !2 & !3 & !4 & !5
0&2
0&3
0&4
0&5
0&6
del 0
name 1 P-loop
name 2 aC-Helix
name 3 C-loop
name 4 A-loop
name 5 rest
name 6 backbone_P-loop
name 7 backbone_aC-Helix
name 8 backbone_C-loop
name 9 backbone_A-loop
name 10 backbone_rest
q
EOD
cd ../
}
done
