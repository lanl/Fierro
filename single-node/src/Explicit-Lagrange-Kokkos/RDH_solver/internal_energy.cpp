#include <Teuchos_RCP.hpp>
#include <MueLu.hpp>

#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void update_internal_energy(DViewCArrayKokkos <double> &zone_sie,
                     const size_t stage,
                     const mesh_t &mesh,
                     const CArrayKokkos <double> &A1,
                     const CArrayKokkos <double> &lumped_mass){
    
    
    FOR_ALL(zone_gid, 0, mesh.num_zones,{

        zone_sie( 1, zone_gid ) = zone_sie(stage, zone_gid) - A1(zone_gid)/lumped_mass(zone_gid);

    });//end for all
    Kokkos::fence();
    

}// end update momentum