
            // FOR_ALL(elem_gid, 0, mesh.num_elems, {

            //     for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                        
            //         int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);

                    
            //         double zone_coords[3]; 
            //         zone_coords[0] = 0.0;
            //         zone_coords[1] = 0.0;
            //         zone_coords[2] = 0.0;

            //         // get the coordinates of the zone center
            //         for (int node_lid = 0; node_lid < mesh.num_nodes_in_zone; node_lid++){
                        
            //             zone_coords[0] += node_coords(stage, mesh.nodes_in_zone(zone_gid, node_lid), 0);
            //             zone_coords[1] += node_coords(stage, mesh.nodes_in_zone(zone_gid, node_lid), 1);
            //             zone_coords[2] += node_coords(stage, mesh.nodes_in_zone(zone_gid, node_lid), 2);
                        
            //         } // end loop over nodes in element
            //         zone_coords[0] = zone_coords[0]/mesh.num_nodes_in_zone;
            //         zone_coords[1] = zone_coords[1]/mesh.num_nodes_in_zone;
            //         zone_coords[2] = zone_coords[2]/mesh.num_nodes_in_zone;

            //         // p = rho*ie*(gamma - 1)
            //         double gamma = elem_state_vars(elem_gid,0); // gamma value

            //         double temp_pres = 0.0;
            //         temp_pres = 0.25*( cos(2.0*PI*zone_coords[0]) + cos(2.0*PI*zone_coords[1]) ) + 1.0;
                    
            //         zone_sie(stage, zone_gid) = temp_pres/((gamma - 1.0));
                    
            //     }// end loop over zones
            // });// end fill sie with TG exact
            
			// Kokkos::fence();
			//
			//
			//
			//

            // FOR_ALL(node_gid, 0, mesh.num_nodes, {
            //     node_vel(stage, node_gid, 0) = sin( PI * node_coords(1, node_gid, 0) ) * cos(PI * node_coords(1, node_gid, 1)); 
            //     node_vel(stage, node_gid, 1) =  -1.0*cos(PI * node_coords(1, node_gid, 0)) * sin(PI * node_coords(1, node_gid, 1)); 
            //     node_vel(stage, node_gid, 2) = 0.0;
            // });// end fill vel with TG exact
            // Kokkos::fence();
        
