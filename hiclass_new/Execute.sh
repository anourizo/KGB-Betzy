#!/bin/bash

# Loop through folders
for ((i=0; i<200; i++))
do
    folder_name="output_$i"
    file_path="../Desktop/SimulationBasedEmulator/kev-hicalss/$folder_name/file_classparameters.ini"

    if [ -f "$file_path" ]; then
        ./class "$file_path"
        echo "Executed command in $folder_name"
    else
        echo "File $file_path not found"
    fi
done

