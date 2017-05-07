#!/bin/sh
#./main.py

#genereate
#pack_alg="elrand"
#file_name="foam_el1_3"
#file_path="results/$file_name"
#nvol=9
#./foam_generate.py -c input.json -o "$file_path.geo" -p $pack_alg -g
#./foam_convert.py -i "$file_path.geo" -o "$file_path" -n $nvol -c geo2fe
#./foam_relax.py -i "$file_path.fe" -o $file_path"_wet" -s 0.6 -p 0.92 -m wet

##analyze foam
#json_files=$file_path*.json
#for json in $json_files
#do
#echo $json
#./foam_analyze.py -i $json --all -f auto
#done


#foam_testing_dry_unrelaxed.json
# foam_testing_wet_relaxed.geo.json
#manual relax
#output foam_testing_relaxed.fe
#./foam_analyzator.py -i results/foam_testing_wet_relaxed.geo.json --all -f wetjson
#./main.py -c batch/exp_sphere.json -gta
#./main.py -c batch/exp_dispsphere.json -gta
#./main.py -c batch/exp_regell.json -gta
#./main.py -c batch/exp_dispell.json -gta
#Generate
#./par_study_annealing.py -c json_config/exp_densesphere.json -f se_cmd/relax.cmd -g -n 10
#./par_study_annealing.py -c json_config/exp_densesphere.json -f se_cmd/relax.cmd -g -n 10
#./par_study_annealing.py -c json_config/exp_densesphere.json -f se_cmd/relax.cmd -g -n 10

cmd_dir="se_cmd"
cfg_dir="json_config"
relax_cmds=(relaxNE1.cmd) # relaxNE0.cmd relaxNE4.cmd   relaxNE0.cmd relaxNE1.cmd
#(relax0.cmd relax1.cmd relax2.cmd relax3.cmd)
config_files=(exp_dispell5.json exp_dispell6.json exp_dispell7.json) # exp_dispell2.json exp_dispell3.json exp_dispell4.json exp_dispell5.json exp_dispell6.json exp_dispell7.json

if [ "$1" = "-g" -o "$2" = "-g" ]
then
    echo "Generating structures..."
    for cfg in "${config_files[@]}"
	do
		./par_study_annealing.py -c $cfg_dir/$cfg -f None -g -n 10
	done
fi

if [ "$1" = "-ra" -o "$2" = "-ra" ]
then
    echo "Relaxing structures and analyzing..."
    for cfg in "${config_files[@]}"
	do
		for cmd in "${relax_cmds[@]}"
		do
			./par_study_annealing.py -c $cfg_dir/$cfg -f $cmd_dir/$cmd -ra -n 10
		done
	done
elif [ "$1" = "-rsa" -o "$2" = "-rsa" ]
then
	echo "Relaxing structures, optimizing strut content and analyzing..."
	for cfg in "${config_files[@]}"
	do
		./par_study_annealing.py -c $cfg_dir/$cfg -f $cmd_dir/$cmd -rsa -n 10
	done
elif [ "$1" = "-a" -o "$2" = "-a" ]
then
	echo "Analyzing..."
	for cfg in "${config_files[@]}"
	do
		./par_study_annealing.py -c $cfg_dir/$cfg -f "nonerelaxed" -a
	done
fi









