//Jiri Kolar
//This file depends on foamface.cmd
//Load this foamface.cmd before
//for example:
read "../se_cmd/foamface.cmd"

//print analysis of wet foam in json format

an_wet_json:= {

	printf "{\"border-body\": %d,\n",border_body;
	//foam_signature; not functional for wet foam
	foamface_mark; //functional, excludes struts
	foamedge_mark; // functional but mark whole face as one edge = perimeter
	
	vertices_count:=max(vertices,id);
	//print foamface_count;
	//foamface_count, foamedge_count defined in foamsignature
	define faces_area real[foamface_count];
	define used_faces integer[foamface_count];
	define used_vertices integer[vertices_count]; 
	face_edge_counts := 0;
	define face_perimeter real[foamedge_count];
	define used_edges integer[foamedge_count];
	edge_counted := 0;
	foreach edge ee where emark > 0 do
	  { 
		face_perimeter[ee.emark] += ee.length;
	  };
	 foreach facet ff where fmark > 0 do
	  { 
		faces_area[ff.fmark] += ff.area;
	  };
	  printf "\"cells\": [";
	  first_cell_output:=1;
	  foreach body bb do{
		if first_cell_output=1 then {
			printf "{";
			first_cell_output:=0;
		}else{
			printf ",{";
		};
		  //print id
		  printf "\"id\": %d,",bb.id;
		  //print volume
		  printf "\"volume\": %f,",bb.volume;
		  //mark which edges,faces,vertices print
		  used_edges:=0;
		  used_faces:=0;
		  used_vertices:=0;
		  foreach bb.facets ff do{
		  foreach ff.edges ee do {
			if ee.emark>0 then
				used_edges[ee.emark]:=1;
			used_vertices[ee.vertex[1].id]:=1;
			used_vertices[ee.vertex[2].id]:=1;
			};
			used_faces[ff.fmark] := 1;
		  };
		  //print edge length
		  printf "\"face-perimeter\": ";
		  printf "[";
		  first_output:=1;
		  nedges:=0;
		  for ( inx := 1 ; inx <= foamedge_count; inx += 1) { 
			if used_edges[inx]=1 then {
				if first_output=1 then {
					printf "%f",face_perimeter[inx];
					first_output:=0;
				} else {
					printf ",%f",face_perimeter[inx];
				};
			 };			
		   };
			printf "],";
		  //print facet area
		  printf "\"facet-area\": ";
		  printf "[";
		  first_output:=1;
		  nedges:=0;
		  for ( inx := 1 ; inx <= foamface_count; inx += 1) { 
			if used_faces[inx]=1 then {
				if first_output=1 then {
					printf "%f",faces_area[inx];
					first_output:=0;
				} else {
					printf ",%f",faces_area[inx];
				};
			 };			
		   };
			printf "],";
		
		  //vertices
		  printf "\"vertices\": ";
		  printf "[";
		  first_output:=1;
		  nedges:=0;
		  for ( inx := 1 ; inx <= vertices_count; inx += 1) { 
			if used_vertices[inx]=1 then {
				if first_output=1 then {
					printf "[%f,%f,%f]",vertices[inx].X,vertices[inx].Y,vertices[inx].Z;
					first_output:=0;
				} else {
					printf ",[%f,%f,%f]",vertices[inx].X,vertices[inx].Y,vertices[inx].Z;
				};
			 };			
		   };
			printf "]}";
	  };
	  printf "]}\n";
}

/*
usage: an_wet_json >>> "an.json"
*/
