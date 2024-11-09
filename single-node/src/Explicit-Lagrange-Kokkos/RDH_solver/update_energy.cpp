#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include "ref_surf_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;

void update_energy(DViewCArrayKokkos <double> &zone_sie,
                   const DViewCArrayKokkos <double> &PSI,
                   const DViewCArrayKokkos <double> &m,
                   const mesh_t &mesh,
                   const int stage){

        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);

                zone_sie(1, zone_gid) = zone_sie(stage, zone_gid) - PSI(stage, elem_gid, zone_lid)/m(zone_gid);
                
            }
        });

}// end update energy