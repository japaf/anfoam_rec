/*
Jiri Kolar
This file depends on foamface.cmd (Load together with this script)
For example: 
	read "../se_cmd/foamface.cmd"

To analyze wet foam
usage: an_wet_json >>> "anwet.json"
To analyze dry foam
usage: an_dry_json >>> "andry.json"
*/

border_body:=0;

//print analysis of dry foam in json format
an_dry_facet_area:= {
	foamface_mark; //functional, excludes struts
	facet_area:=0;
	foreach facet ff where fmark > 0 do
		  { 
			facet_area += ff.area;
		  };
	printf "%f, %f \n", facet_area, total_area;
}

an_wet_facet_area:= {
	foamface_mark; //functional, excludes struts
	facet_area:=0;
	maxfmark:=foamface_count-max(bodies,id)+1;
	foreach facet ff where fmark > 0 and fmark <= maxfmark do
		  { 
			facet_area += ff.area;
		  };
	printf "%f, %f \n", facet_area, total_area;
}


an_dry_json:= {
	printf "{\"cell-types\": ";
	foam_signature;
	vertices_count:=max(vertices,id);
	//foamface_count, foamedge_count defined in foamsignature
	define faces_area real[foamface_count];
	define used_faces integer[foamface_count];
	define used_vertices integer[vertices_count];	
	define edges_length real[foamedge_count];
	define used_edges integer[foamedge_count];
	face_edge_counts := 0;
	edge_counted := 0;
	edges_length:=0;
	faces_area:=0;
	foreach edge ee where emark > 0 do
	  { 
		edges_length[ee.emark] += ee.length;
	  };
	 foreach facet ff where fmark > 0 do
	  { 
		faces_area[ff.fmark] += ff.area;
	  };
	  printf ",\n";
	  
	  printf "\"facet-area\": [";
	  first_output:=1;
	  for ( inx := 1 ; inx <= foamface_count; inx += 1) { 
			if first_output=1 then {
				printf "%f",faces_area[inx];
				first_output:=0;
			} else {
				printf ",%f",faces_area[inx];
			};		
	   };
	   printf "],";  
	   
	  printf "\"angles\": ";
	  printf "[";
	  vvfirst:=1;
	  anglefirst:=1;
	  foreach vertex vv do {
		if sum(vv.edge eee,eee.emark>0)==4 then {
			if vvfirst==1 then {
				printf "[";
			} else {
				printf ",[";
			};
			vvfirst:=0;
			anglefirst:=1;
			foreach vv.edge ee0 where emark>0 do {
			foreach vv.edge ee1 where emark>0 and id>ee0.id do {
			if ee0.id!=ee1.id then {
			dir0:=1;
			dir1:=1;
			if ee0.vertices[1].id!=vv.id then
				dir0:=dir0*-1;
			if ee1.vertices[1].id!=vv.id then
				dir1:=dir1*-1;
			norm0:=sqrt(ee0.x**2+ee0.y**2+ee0.z**2);
			norm1:=sqrt(ee1.x**2+ee1.y**2+ee1.z**2);
			cosa:=dir0*ee0.edge_vector*dir1*ee1.edge_vector/norm0/norm1;
			tetrahedral_angle:=acos(cosa)*180/pi;
			if anglefirst==1 then {
				printf "%f",tetrahedral_angle;
			}else{
				printf ",%f",tetrahedral_angle;
			};
			anglefirst:=0;
			};
			};
			};
			printf "]";
		}
	  };
	  printf "]";
	  printf ",\n";
	  
	  //cell descriptors
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
		  printf "\"edge-length\": ";
		  printf "[";
		  first_output:=1;
		  nedges:=0;
		  for ( inx := 1 ; inx <= foamedge_count; inx += 1) { 
			if used_edges[inx]=1 then {
				if first_output=1 then {
					printf "%f",edges_length[inx];
					first_output:=0;
				} else {
					printf ",%f",edges_length[inx];
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
	faces_area:=0;
	face_perimeter:=0;
	edge_counted := 0;
		
	foreach edge ee where emark > 0 do
	  { 
		face_perimeter[ee.emark] += ee.length;
	  };
	 foreach facet ff where fmark > 0 do
	  { 
		faces_area[ff.fmark] += ff.area;
	  };
	  
	  printf "\"facet-area\": [";
		  first_output:=1;
		  for ( inx := 1 ; inx <= foamface_count; inx += 1) { 
				if first_output=1 then {
					printf "%f",faces_area[inx];
					first_output:=0;
				} else {
					printf ",%f",faces_area[inx];
				};		
		   };
			printf "],";
			
	 //cell info
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
To analyze wet foam
usage: an_wet_json >> "anwet.json"
To analyze dry foam
usage: an_dry_json >>> "andry.json"
*/
