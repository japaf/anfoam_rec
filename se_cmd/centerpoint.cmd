define vertex attribute vmark integer

procedure compute_centerpoint(integer vol_id)
{
 set vertex vmark 0;
 foreach facet ff where ff.frontbody=vol_id or ff.backbody=vol_id do {
	 for ( inx := 1 ; inx <= ff.valence; inx += 1) {
	  ff.vertices[inx].vmark:=1;
	 };
 };
 cx:=0;
 cy:=0;
 cz:=0;
 nn:=0;
 foreach vertex vv where vmark=1 do {
 nn:=nn+1;
 cx:=cx+vv.x;
 cy:=cy+vv.y;
 cz:=cz+vv.z;
 };
 printf "%f %f %f\n", cx/nn,cy/nn,cz/nn;
}

compute_all_centerpoints:= {
 detorus;
 foreach body bb do {
	compute_centerpoint(bb.id); 
 };
}
