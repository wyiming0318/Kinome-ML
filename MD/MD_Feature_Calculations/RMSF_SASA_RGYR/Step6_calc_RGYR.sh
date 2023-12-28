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
start=28

############################### P-loop ################################
end=$(wc -l < $name_RGYR1.xvg)
# calculate std
#####std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR1.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR1.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RGYR_P-loop.txt

############################### aC-Helix ################################
end=$(wc -l < $name_RGYR2.xvg)
# calculate std
#####std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR2.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR2.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RGYR_aC-Helix.txt

############################### C-loop ################################
end=$(wc -l < $name_RGYR3.xvg)
# calculate std
#####std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR3.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR3.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RGYR_C-loop.txt

############################### A-loop ################################
end=$(wc -l < $name_RGYR4.xvg)
# calculate std
#####std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR4.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR4.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RGYR_A-loop.txt

############################### rest ################################
end=$(wc -l < $name_RGYR5.xvg)
# calculate std
#####std=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR5.xvg | awk '{ total += $2; totalsq+=($2*$2) } END { print sqrt(((totalsq/NR)-(total/NR)*(total/NR))) }')
# calculate avg
avg=$(nawk -vs="$start" -ve="$end" 'NR>e{exit}NR>=s' $name_RGYR5.xvg | awk '{ total += $2;} END { print (total/NR)}')
# write to file
echo $name $avg >> ../data_RGYR_rest.txt

cd ../
}
done
