#ifndef GEOMETRY_H
#define GEOMETRY_H 

#include "elements.h"
#include "swage.h"

void get_gauss_pt_jacobian(swage::mesh_t& mesh, elements::ref_element& ref_elem);

void get_vol_jacobi(swage::mesh_t& mesh, elements::ref_element& ref_elem);

void get_vol_hex(swage::mesh_t& mesh, elements::ref_element& ref_elem);





#endif // GEOMETRY_H

