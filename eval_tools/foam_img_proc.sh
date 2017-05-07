#!/bin/sh
#ristretto relax0before*.png&
#ristretto relax1*.png&
#ristretto relax2*.png&
#ristretto relax3*.png&

relax_cmds=(relaxNE0.cmd relaxNE1.cmd) #relaxNE1 relaxNE2 relaxNE3
#(relax0 relax1 relax2 relax3)
for cmd in "${relax_cmds[@]}"
do
	convert $cmd"before"*.png -append "result_"$cmd"before_dry.png"
	convert $cmd"after-dry"*.png -append "result_"$cmd"after_dry.png"
	#convert $cmd"after-wet"*.png -append "result_"$cmd"after_wet.png"
	ristretto "result_"$cmd*.png&
done



#convert relax0before*.png -append relax0before_dry.png
#convert relax0(after-dry*.png -append relax0(after_dry.png
#convert relax0(after-wet*.png -append relax0(after_wet.png
