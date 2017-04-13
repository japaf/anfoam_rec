/*
Print volume and surface of border body created by wetfoam2.cmd
print_stat
*/
function real comp_area(integer maxid)
{
	tot_area:=0.0;
	foreach body bb do {
		foreach bb.facet ff where ff.frontbody<maxid and ff.backbody<maxid do {
			tot_area:=tot_area+ff.area;
		};
	};
	return tot_area/2;
}

print_stat:={
	printf "{\"area\": %f,", comp_area(max(body,id));
	printf "\"volume\": %f}\n", bodies[max(body,id)].volume;
}

