/*
Print volume and surface of border body created by wetfoam2.cmd
print_stat
*/
function real comp_area(integer maxid)
{
	foamface_mark; //functional, excludes struts
	facet_area:=0;
	maxfmark:=foamface_count-max(bodies,id)+1;
	foreach facet ff where fmark > 0 and fmark <= maxfmark do
		  { 
			facet_area += ff.area;
		  };
/*	
tot_area:=0.0;
	foreach body bb do {
		foreach bb.facet ff where ff.frontbody<maxid and ff.backbody<maxid do {
			tot_area:=tot_area+ff.area;
		};
	};
*/
	return facet_area;
}

print_stat:={
	printf "{\"area\": %f,", comp_area(max(body,id));
	printf "\"volume\": %f}\n", bodies[max(body,id)].volume;
}

