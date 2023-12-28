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
start=18

############################### P-loop ################################
end=$(wc -l < $name_RMSF1.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF1.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF1.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RMSF_P-loop.txt

############################### aC-Helix ################################
end=$(wc -l < $name_RMSF2.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF2.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF2.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RMSF_aC-Helix.txt

############################### C-loop ################################
end=$(wc -l < $name_RMSF3.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF3.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF3.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RMSF_C-loop.txt

############################### A-loop ################################
end=$(wc -l < $name_RMSF4.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF4.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF4.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RMSF_A-loop.txt

############################### rest ################################
end=$(wc -l < $name_RMSF5.xvg)
# calculate std
#std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF5.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RMSF5.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RMSF_rest.txt
cd ../
}
done
