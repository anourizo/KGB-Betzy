#!/bin/bash

i=0

while read c1; do w_kess_arr[i++]=$c1; done < w0.txt


for ((i = 0 ; i < 12 ; i++)); do

w_kessence=${w_kess_arr[$i]}
echo $w_kessence

# Calculate using Python
par_smg=$(python -c "print(3 * ((1 + ($w_kessence)) / 1e-07))")

output_dir='output_'$i''
#echo $output_dir;

#########
[ ! -d "./outputs/$output_dir" ] && mkdir -p "./outputs/$output_dir"


sed -e 's/expansion_smg = 6.878623e-01, -9.000000e-01, 0.000000e+00/expansion_smg = 6.878623e-01, '$w_kessence', 0.000000e+00/g' -e 's/parameters_smg = 3.000000e+06, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00/parameters_smg = '$par_smg', 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00/g' -e 's#root = output/file_class#root = 'outputs/$output_dir/file_class'#g' file_classparameters.ini> file_classparameters_$output_dir.ini

sed 's/file_classparameters.ini/file_classparameters_'$output_dir'.ini/g' Run_hiclass.sh> Run_number$i.sh
#
bash Run_number$i.sh

mv file_classparameters_$output_dir.ini ./outputs/$output_dir
mv Run_number$i.sh ./outputs/$output_dir

done
