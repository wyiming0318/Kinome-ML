for d in */

do {
	
cd ./$d
	
ID="$(basename "$PWD")" # get mutation ID
name=$(basename "$PWD")
name="$(basename "$PWD")_SASA"
start=19

############################### P-loop ################################
end=$(wc -l < $name.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $3; totalsq+=($3*$3) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $3;} END { print (total/NR)}')
# write to file
echo $ID $avg >> ../data_SASA_P-loop.txt

############################### aC-Helix ################################
end=$(wc -l < $name.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $4; totalsq+=($4*$4) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $4;} END { print (total/NR)}')
# write to file
echo $ID $avg >> ../data_SASA_aC-Helix.txt

############################### C-loop ################################
end=$(wc -l < $name.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $5; totalsq+=($5*$5) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $5;} END { print (total/NR)}')
# write to file
echo $ID $avg >> ../data_SASA_C-loop.txt

############################### A-loop ################################
end=$(wc -l < $name.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $6; totalsq+=($6*$6) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $6;} END { print (total/NR)}')
# write to file
echo $ID $avg >> ../data_SASA_A-loop.txt

############################### rest ################################
end=$(wc -l < $name.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $7; totalsq+=($7*$7) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name.xvg | awk '{ total += $7;} END { print (total/NR)}')
# write to file
echo $ID $avg >> ../data_SASA_rest.txt
cd ../
}
done
