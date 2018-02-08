#!/bin/sh
#generate all diff exp

cmd_dir="se_cmd"
cfg_dir="json_config"
relax_cmd=relaxNE5.cmd

#config_files=(diff_sp02.json diff_sp01.json diff_sp03.json diff_sp04.json diff_sp05.json)
config_files=(diff_resf01.json diff_resf02.json diff_resf03.json diff_resf04.json diff_resf05.json diff_resf06.json diff_resf07.json)
#config_files=(diff_sp05.json)

if [ "$1" = "-ra" -o "$2" = "-ra" ]
then
    echo "Relaxing structures and analyzing..."
    for cfg in "${config_files[@]}"
	do
	./par_study_annealing.py -c $cfg_dir/$cfg -f $cmd_dir/$relax_cmd -sa -n 1
	done
fi

if [ "$1" = "-cd" -o "$2" = "-cd" ]
then
	echo "Optimizing strut content and porosity, analyzing and generate output for diffusion model..."
	for cfg in "${config_files[@]}"
	do
		./par_study_annealing.py -c $cfg_dir/$cfg -f $cmd_dir/$relax_cmd -spda --const-res -n 1
	done
fi

if [ "$1" = "-sd" -o "$2" = "-sd" ]
then
	echo "Optimizing strut content and porosity, analyzing and generate output for diffusion model..."
	for cfg in "${config_files[@]}"
	do
		./par_study_annealing.py -c $cfg_dir/$cfg -f $cmd_dir/$relax_cmd -spda -n 1
	done
fi

