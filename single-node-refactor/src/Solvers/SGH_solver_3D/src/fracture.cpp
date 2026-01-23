#include <stdio.h>
#include </home/alexholmes814/MATAR/src/include/matar.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "mesh.h"
#include "state.h"
#include "fracture.h"
#include "fracture_stress_bc.h"

using namespace mtr; // matar namespace

cohesive_zones_t::cohesive_zones_t() {
// constructor for cohesive zones
}

// initialize the identification of cohesive zones
// this is an algorithim for identifying cohesive zones in a mesh
// it loops over all nodal coordinates and identifies overlapping nodal coordinate pairs
void cohesive_zones_t::initialize(Mesh_t& mesh, State_t& State){
    // the following code counts the number of boundary nodes and checks for node overlaps (2 nodes with the same coordinates)
    // this is the beginning step to setting up cohesive zones for fracture
                   
    // counting the number of boundary nodes
    size_t num_bdy_nodes = mesh.num_bdy_nodes;
    //std::cout << "Number of boundary nodes: " << num_bdy_nodes << std::endl;
    printf("Total boundary nodes: %zu\n", mesh.num_bdy_nodes);

    // total number of boundary nodes across all sets
    //size_t total_bdy_nodes = 0;
    //for (size_t i = 0; i < mesh.num_bdy_sets; ++i) {
    //    std::cout << "Boundary nodes in set " << i << ": " << mesh.num_bdy_nodes_in_set(i) << std::endl;
    
    const double tol = 1e-8; //0.000001; //e-3; // adjust as needed; added just in case coordinate pairs are close but not exactly equal
    size_t overlap_index = 0; // counts unique overlapping nodes (2 unique overlapping nodes = 1 overlapping node pair)
    size_t pair_count = 0; // counts how many overlapping node pairs exist
    
    
    // count unique overlapping nodes
    for (size_t i = 0; i < num_bdy_nodes; ++i) {
        size_t node_i = mesh.bdy_nodes(i);
        for (size_t j = i + 1; j < num_bdy_nodes; ++j) {
            size_t node_j = mesh.bdy_nodes(j);

            bool overlap = true;
            for (size_t k = 0; k < 3; ++k) {
                if (std::abs(State.node.coords(node_i, k) - State.node.coords(node_j, k)) > tol) {
                    overlap = false;
                    break;
                }
            }

            if (overlap) {
                ++pair_count;
            }
        }
    }

    
    // allocate only the size of overlapping nodes 
    //CArrayKokkos<size_t> overlapping_node_gids(pair_count, 2, "overlapping_node_gids");
    overlapping_node_gids = CArrayKokkos<size_t> (pair_count, 2, "overlapping_node_gids");

    // second pass: store actual overlapping node pairs
    size_t pair_index = 0; // fills the rows (pairs) that are added to 2D overlapping_node_gids array 
    
    // store node IDs in the array
    for (size_t i = 0; i < num_bdy_nodes; ++i) {
        size_t node_i = mesh.bdy_nodes(i);
        for (size_t j = i + 1; j < num_bdy_nodes; ++j) {
            size_t node_j = mesh.bdy_nodes(j);

            bool overlap = true;
            for (size_t k = 0; k < 3; ++k) {
                if (std::abs(State.node.coords(node_i, k) - State.node.coords(node_j, k)) > tol) {
                    overlap = false;
                    break;
                }
            }

            if (overlap) {
                //++pair_count;
                printf("Overlap (cohesive zone) found between node %zu and node %zu\n", node_i, node_j);
                overlapping_node_gids(pair_index, 0) = node_i;
                overlapping_node_gids(pair_index, 1) = node_j;
                ++pair_index;
               
            }
        }
    }

    printf("Total overlapping node pairs: %zu\n", pair_count);

    // print overlapping node coordinates
    for (size_t i = 0; i < pair_index; ++i) {
        size_t node_i = overlapping_node_gids(i, 0);
        size_t node_j = overlapping_node_gids(i, 1);

        printf("Overlapping Pair: %zu <-> %zu\n", node_i, node_j);

        printf("    Node %zu coords: ", node_i);
        for (size_t k = 0; k < 3; ++k) {
            printf("%g ", State.node.coords(node_i, k));
        }
        printf("\n");

        printf("    Node %zu coords: ", node_j);
            for (size_t k = 0; k < 3; ++k) {
                printf("%g ", State.node.coords(node_j, k));
            }
            printf("\n");
    }

    // ======================== test for function for cohesive_zone_elem_count in fracture.cpp: which finds the max number of elements that any cohesive zone node is part of ========================
    //size_t max_elem_in_cohesive_zone = cohesive_zone_elem_count(overlapping_node_gids, mesh.elems_in_node, mesh);
    max_elem_in_cohesive_zone = cohesive_zone_elem_count(overlapping_node_gids, mesh.elems_in_node, mesh);
    printf("Max elements connected to any cohesive zone node: %zu\n", max_elem_in_cohesive_zone);
    // ======================== END test for function in cohesive_zone_elem_count fracture.cpp: which finds the max number of elements that any cohesive zone node is part of ========================
 
    
    // ======================== face-by-face cross-check debug (compute_face_geometry) ========================
{
    // quick sanity check
    if (mesh.num_nodes_in_elem != 8 || mesh.num_dims != 3) {
        printf("[debug] face-geometry check only implemented for HEX8/3D\n");
    } else {
        printf("======================== face-by-face cross-check ========================\n");

        // print the nodes in each element
        for (size_t elem = 0; elem < mesh.num_elems; ++elem) {
            printf("Element %zu nodes: ", elem);
            for (size_t ln = 0; ln < mesh.num_nodes_in_elem; ++ln) {
                printf("%zu ", mesh.nodes_in_elem(elem, ln));
            }
            printf("\n");
        
        // // loop faces
        // // for each surf, look up the patch_id
        // // patch_id = mesh.patches_in_elem(elem, surf)
        // // then fetch four node IDs: g[a] = mesh.nodes_in_patch(patch_id, a)
        // // so node order is whatever nodes_in_patch says for that patch_id
         for (size_t surf = 0; surf < 6; ++surf) {


        // mesh patch mapping (Fierro nodal indexing convention)
        const size_t patch_id = mesh.patches_in_elem(elem, surf);
        size_t face_gid[4];
        for (size_t a = 0; a < 4; a++){
            face_gid[a] = mesh.nodes_in_patch(patch_id, a);
        } 
                // print face node IDs
                printf("  Face %zu node IDs: %zu %zu %zu %zu\n",
                        surf, face_gid[0], face_gid[1], face_gid[2], face_gid[3]);

                DCArrayKokkos<double> &X = State.node.coords; // (num_nodes x 3)

                printf("    coords[gid0]: %g %g %g\n", X(face_gid[0],0), X(face_gid[0],1), X(face_gid[0],2));
                printf("    coords[gid1]: %g %g %g\n", X(face_gid[1],0), X(face_gid[1],1), X(face_gid[1],2));
                printf("    coords[gid2]: %g %g %g\n", X(face_gid[2],0), X(face_gid[2],1), X(face_gid[2],2));
                printf("    coords[gid3]: %g %g %g\n", X(face_gid[3],0), X(face_gid[3],1), X(face_gid[3],2));

                // simple centroid check: average of 4 face nodes
                double cx_simple = 0.25 * (X(face_gid[0],0) + X(face_gid[1],0) + X(face_gid[2],0) + X(face_gid[3],0));
                double cy_simple = 0.25 * (X(face_gid[0],1) + X(face_gid[1],1) + X(face_gid[2],1) + X(face_gid[3],1));
                double cz_simple = 0.25 * (X(face_gid[0],2) + X(face_gid[1],2) + X(face_gid[2],2) + X(face_gid[3],2));

                
                // stack buffers + ViewCArrayKokkos wrappers 
                // stack buffer = temporary memory that only exists in this function scope (temp raw storage)
                double n_buf[3], r_buf[3], s_buf[3], cen_buf[3];
                ViewCArrayKokkos<double> n(&n_buf[0], 3);
                ViewCArrayKokkos<double> r(&r_buf[0], 3);
                ViewCArrayKokkos<double> s(&s_buf[0], 3);
                ViewCArrayKokkos<double> cenface(&cen_buf[0], 3);

                // call compute_face_geometry and compare
                compute_face_geometry(
                    State.node.coords,   
                    mesh,
                    State.node.coords,   
                    mesh.nodes_in_elem,  
                    surf,
                    elem,
                    n, r, s, cenface
                );

                // clean math before printing (prevent -0.0s in output)
                auto pz = [](double v) { return (std::fabs(v) < 1e-13) ? 0.0 : v; };

                // print calculated results
                printf("    centroid(simple avg): %g %g %g\n", cx_simple, cy_simple, cz_simple);
                printf("    centroid(computed):   %g %g %g\n", cenface(0), cenface(1), cenface(2));
                printf("    r: %g %g %g\n", r(0), r(1), r(2));
                printf("    s: %g %g %g\n", s(0), s(1), s(2));
                printf("    n: %g %g %g\n", n(0), n(1), n(2));
            }
        }

        printf("==========================================================================\n");
    }
}
//    ======================== END face-by-face cross-check debug (compute_face_geometry) ========================

//    ======================== cohesive_zone_info debug ========================

//    CArrayKokkos<int> cz_info = build_cohesive_zone_info(
//        mesh,
//        State,
//       overlapping_node_gids,
//        max_elem_in_cohesive_zone,
//        tol
//    );
    cz_info = build_cohesive_zone_info(
        mesh,
        State,
        overlapping_node_gids,
        max_elem_in_cohesive_zone,
        tol
    );

{

    printf("\n================== cohesive_zone_info debug ==================\n");
    printf("num_overlapping_node_pairs=%zu  max_elem_in_cohesive_zone=%zu\n",
           overlapping_node_gids.dims(0), max_elem_in_cohesive_zone);

    

    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
        if (i != 0) continue; // ONLY DO FIRST PAIR FOR TEST AGAINST TIME
        const size_t nodeA = overlapping_node_gids(i,0);
        const size_t nodeB = overlapping_node_gids(i,1);

        printf("\n-- Pair %zu  (A (node gid) = %zu, B (node gid) = %zu) --\n", i, nodeA, nodeB);

        // print A-side elems
        //printf("  A elems (IDs): ", mesh.elems_in_node.stride(nodeA));
        // unused argument: mesh.elems_in_node.stride(nodeA)
        printf("  A elems (IDs): ");
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {
            const int eA = cz_info(i, 0*max_elem_in_cohesive_zone + j);
            if (eA >= 0) printf("%d ", eA);
        }
        printf("\n");

        // print B-side elems
        //printf("  A elems (IDs): ", mesh.elems_in_node.stride(nodeB));
        // unused argument: mesh.elems_in_node.stride(nodeB)
        printf("  B elems (IDs): ");
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {
            const int eB = cz_info(i, 1*max_elem_in_cohesive_zone + j);
            if (eB >= 0) printf("%d ", eB);
        }
        printf("\n");

        // print stored local corners
        printf("  kA (local corner index): ");
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {
            const int kA = cz_info(i, 4*max_elem_in_cohesive_zone + j);
            if (kA >= 0) printf("%d ", kA);
        }
        printf("\n");

        printf("  kB (local corner index): ");
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {
            const int kB = cz_info(i, 5*max_elem_in_cohesive_zone + j);
            if (kB >= 0) printf("%d ", kB);
        }
        printf("\n");

        // print everything thats in cohesive_zone_info: check what was stored (-1 means empty slot)
        printf(" checking everything thats in cohesive_zone_info: ");
        for (size_t j = 0; j < 6*max_elem_in_cohesive_zone; ++j){
            printf("%d ", cz_info(i, j));
        }
        printf("\n");

        
        // ---------- re-run the *same* candidate-face search (ABS distance) and print the true first match ----------
        // small stack buffers + views
        double nA_buf[3], rA_buf[3], sA_buf[3], cA_buf[3];
        double nB_buf[3], rB_buf[3], sB_buf[3], cB_buf[3];
        ViewCArrayKokkos<double> nA(&nA_buf[0],3), rA(&rA_buf[0],3), sA(&sA_buf[0],3), cA(&cA_buf[0],3);
        ViewCArrayKokkos<double> nB(&nB_buf[0],3), rB(&rB_buf[0],3), sB(&sB_buf[0],3), cB(&cB_buf[0],3);

        printf("  Stored matches (all):\n");

        int printed = 0;

        // reuse buffers 
        double nA_buf2[3], rA_buf2[3], sA_buf2[3], cA_buf2[3];
        double nB_buf2[3], rB_buf2[3], sB_buf2[3], cB_buf2[3];
        ViewCArrayKokkos<double> nA2(&nA_buf2[0],3), rA2(&rA_buf2[0],3), sA2(&sA_buf2[0],3), cA2(&cA_buf2[0],3);
        ViewCArrayKokkos<double> nB2(&nB_buf2[0],3), rB2(&rB_buf2[0],3), sB2(&sB_buf2[0],3), cB2(&cB_buf2[0],3);

        for (size_t slotA = 0; slotA < max_elem_in_cohesive_zone; ++slotA) {
            const int eA = cz_info(i, 0*max_elem_in_cohesive_zone + slotA);
            const int fA = cz_info(i, 2*max_elem_in_cohesive_zone + slotA); // face chosen for this A-slot
            if (eA < 0 || fA < 0) continue; // no stored match in this A-slot

            // geometry for stored A face
            compute_face_geometry(State.node.coords, mesh,
                                 State.node.coords, mesh.nodes_in_elem,
                                 (size_t)fA, (size_t)eA,
                                 nA2, rA2, sA2, cA2);

            // find which B-slot was paired (the builder filled exactly one B-slot face)
            int partnerB = -1, eB_hit = -1, fB_hit = -1;
            double dist_hit = 0.0, dot_hit = 0.0;

            for (size_t slotB = 0; slotB < max_elem_in_cohesive_zone; ++slotB) {
                const int eB = cz_info(i, 1*max_elem_in_cohesive_zone + slotB);
                const int fB = cz_info(i, 3*max_elem_in_cohesive_zone + slotB); // face chosen for this B-slot
                if (eB < 0 || fB < 0) continue;

                // geometry for stored B face
                compute_face_geometry(State.node.coords, mesh,
                                      State.node.coords, mesh.nodes_in_elem,
                                      (size_t)fB, (size_t)eB,
                                       nB2, rB2, sB2, cB2);

                const double dx = cA2(0) - cB2(0);
                const double dy = cA2(1) - cB2(1);
                const double dz = cA2(2) - cB2(2);
                const double dist = sqrt(dx*dx + dy*dy + dz*dz);
                const double dot  = nA2(0)*nB2(0) + nA2(1)*nB2(1) + nA2(2)*nB2(2);

                if (dist <= tol && dot <= -1.0 + tol) {
                    partnerB = (int)slotB;
                    eB_hit = eB; fB_hit = fB;
                    dist_hit = dist; dot_hit = dot;
                    break; // exactly one B-slot should partner this A-slot
                }
            }

            ++printed;
            if (partnerB >= 0) {
                printf("    [%d] A(slot=%zu): elem=%d face=%d  <->  B(slot=%d): elem=%d face=%d\n",
                        printed, slotA, eA, fA, partnerB, eB_hit, fB_hit);
                printf("         cen(A)=(%.6g, %.6g, %.6g) nA=(%.6g, %.6g, %.6g)\n",
                        cA2(0), cA2(1), cA2(2), nA2(0), nA2(1), nA2(2));
                printf("         cen(B)=(%.6g, %.6g, %.6g) nB=(%.6g, %.6g, %.6g)\n",
                        cB2(0), cB2(1), cB2(2), nB2(0), nB2(1), nB2(2));
                printf("    checks: |dcentroid|=%.6g  (tol=%.6g)   dot(nA,nB)=%.6g\n",
                        dist_hit, tol, dot_hit);
            } else {
                     // shouldn't happen if builder filled both sides; helpful sanity check
                     printf("    [%d] A(slot=%zu): elem=%d face=%d  <->  B: NOT FOUND (unexpected)\n",
                     printed, slotA, eA, fA);
                    }
    }

    if (printed == 0) printf("   No matched faces found.\n");
}
    // oriented debug
// // THROTTLE OUTPUT COMMENT OUT //
//     printf("[CZ::init] this=%p  pairs=%zu  maxcz=%zu  info_rows=%zu\n",
//            (void*)this,
//             overlapping_node_gids.dims(0),
//             max_elem_in_cohesive_zone,
//             cz_info.dims(0));
// // THROTTLE OUTPUT COMMENT OUT //
    // end oriented debug

    printf("\n==============================================================\n");
} // end cohesive_zone_info debug
    // ======================== END cohesive_zone_info debug ========================
} // end cohesive_zones_t::initialize()

    // ======================== cohesive_zone_orientation debug ========================
// COMMENT OUT HERE TO STOP FUNCTION DEBUGGING PRINTS
void cohesive_zones_t::debug_oriented(Mesh_t& mesh,
                                      State_t& State,
                                      CArrayKokkos<size_t>& overlap,
                                      CArrayKokkos<int>& info,
                                      size_t maxcz,
                                      double tol) {
    //    ======================== cohesive_zone_orientation debug ========================

    // debug oriented
// // THROTTLE OUTPUT COMMENT OUT //
//     printf("[CZ::debug] this=%p  pairs(in)=%zu  maxcz(in)=%zu  info_rows(in)=%zu\n",
//           (void*)this, overlap.dims(0), maxcz, info.dims(0));
// // THROTTLE OUTPUT COMMENT OUT //
    // end debug oriented


    // call oriented() to compute cohesive zone normals at t and t+dt 
    CArrayKokkos<double> cohesive_zone_orientation(overlap.dims(0), 6, "cohesive_zone_orientation");

    // initialize to zero
    //cohesive_zone_orientation.set_values(0.0);
    //overlapping_node_gids = CArrayKokkos<size_t> (overlapping_node_gids.dims(0), 2, "overlapping_node_gids");
    //cz_info = CArrayKokkos<int> cz_info;
    //tol = double tol;

{
// // THROTTLE OUTPUT COMMENT OUT //
//     printf("\n================== cohesive_zone_orientation debug ==================\n");
// // THROTTLE OUTPUT COMMENT OUT //
    // checking execute() to make sure inputs are correct
    //printf("execute() check: num_overlapping_node_pairs=%zu overlapping_node_gids=(%zu x %zu) max_elem_in_cohesive_zone=%zu\n", 
    //       overlapping_node_gids.dims(0), overlapping_node_gids.dims(1),
    //       max_elem_in_cohesive_zone, cz_info.dims(0));

    // exiting loop if overlapping node pairs = 0
    //if (overlapping_node_gids.dims(0) == 0){
    //    printf("execute() check: no overlapping node pairs found, skipping cohesive_zone_orientation debug\n");
    //    printf("\n======================================================\n");
    //    return;
    //}

    // assigning pos to State.node.coords
    DCArrayKokkos<double> &pos   = State.node.coords;//State.node.coords_n0; // nodes + ut
    //DCArrayKokkos<double> &pos = State.node.coords; // nodes + ut + us

    // debug: checking State.node.coords_n0 for each cycle
    //for( int i = 0; i < mesh.num_nodes; i++) {
    //    for (int j = 0; j < 3; j++) {
    //        printf("%f  State.node.coords_n0 for each cycle\n", State.node.coords_n0(i,j));
    //    }
    //    printf("\n");
    //}   

    // call oriented()
    //oriented(mesh, pos, pos, overlap, info, maxcz, tol, cohesive_zone_orientation);
    oriented(mesh, pos, overlap, info, maxcz, tol, cohesive_zone_orientation);


    // loop over overlapping node pairs
    for (size_t i = 0; i < overlap.dims(0); ++i) {
        if (i != 0) continue; // ONLY DO FIRST PAIR FOR TEST AGAINST TIME
        const size_t nodeA = overlap(i,0);
        const size_t nodeB = overlap(i,1);
// // THROTTLE OUTPUT COMMENT OUT //
//         printf("\n-- Pair %zu  (A gid=%zu, B gid=%zu) --\n", i, nodeA, nodeB);
// // THROTTLE OUTPUT COMMENT OUT //
        // accumulators exactly like oriented()
        double sum_t [3] = {0.0, 0.0, 0.0};
        double sum_dt[3] = {0.0, 0.0, 0.0};
        int cnt = 0;

        // temp views for compute_face_geometry
        double nA_t_buf[3], rA_t_buf[3], sA_t_buf[3], cA_t_buf[3];
        double nB_t_buf[3], rB_t_buf[3], sB_t_buf[3], cB_t_buf[3];
        double nA_dt_buf[3], rA_dt_buf[3], sA_dt_buf[3], cA_dt_buf[3];
        ViewCArrayKokkos<double> nA_t (&nA_t_buf[0], 3),  rA_t (&rA_t_buf[0], 3),  sA_t (&sA_t_buf[0], 3),  cA_t (&cA_t_buf[0], 3);
        ViewCArrayKokkos<double> nB_t (&nB_t_buf[0], 3),  rB_t (&rB_t_buf[0], 3),  sB_t (&sB_t_buf[0], 3),  cB_t (&cB_t_buf[0], 3);
        ViewCArrayKokkos<double> nA_dt(&nA_dt_buf[0], 3), rA_dt(&rA_dt_buf[0], 3), sA_dt(&sA_dt_buf[0], 3), cA_dt(&cA_dt_buf[0], 3);

        // contributors header
// // THROTTLE OUTPUT COMMENT OUT //
//         printf("  contributors (A-side matched faces over all slots):\n");  
// // THROTTLE OUTPUT COMMENT OUT //        
        // walk over all slot-keyed A-side matches and accumulate (blocks [0] elems, [2] faces)
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {
            const int eA = info(i, 0*max_elem_in_cohesive_zone + j); // A elem at slot j
            const int fA = info(i, 2*max_elem_in_cohesive_zone + j); // A face at slot j
            const int kA = info(i, 4*max_elem_in_cohesive_zone + j); // A local corner slot j
            if (eA < 0 || fA < 0) continue; // skip if -1        

        // geometry: A at t, A at t+dt, B at t (oriented uses A for orientation; B is sanity check)
        compute_face_geometry(pos,   mesh, pos,   mesh.nodes_in_elem,
                              static_cast<size_t>(fA), static_cast<size_t>(eA),
                              nA_t, rA_t, sA_t, cA_t);
        compute_face_geometry(pos, mesh, pos, mesh.nodes_in_elem,
                              static_cast<size_t>(fA), static_cast<size_t>(eA),
                              nA_dt, rA_dt, sA_dt, cA_dt);
        //compute_face_geometry(pos,   mesh, pos,   mesh.nodes_in_elem,
        //                      static_cast<size_t>(fB), static_cast<size_t>(eB),
        //                      nB_t, rB_t, sB_t, cB_t);

            // accumulate like oriented()
            sum_t [0] += nA_t (0); sum_t [1] += nA_t (1); sum_t [2] += nA_t (2);
            sum_dt[0] += nA_dt(0); sum_dt[1] += nA_dt(1); sum_dt[2] += nA_dt(2);
            ++cnt;
// // THROTTLE OUTPUT COMMENT OUT //
//             // per-face print
//             printf("  A[j=%zu]: eA elem=%d fA face=%d  kA local corner=%d "
//                    "cen_t=(%.6g, %.6g, %.6g) n_t=(%.6g, %.6g, %.6g)  |  "
//                    "cen_tdt=(%.6g, %.6g, %.6g) n_tdt=(%.6g, %.6g, %.6g)\n",
//                    j, eA, fA, kA,
//                    cA_t(0),  cA_t(1),  cA_t(2),  nA_t(0),  nA_t(1),  nA_t(2),
//                    cA_dt(0), cA_dt(1), cA_dt(2), nA_dt(0), nA_dt(1), nA_dt(2));
// // THROTTLE OUTPUT COMMENT OUT //
            // accumulate normals (exactly like oriented())
            //sum_t [0] += nA_t (0);  sum_t [1] += nA_t (1);  sum_t [2] += nA_t (2);
            //sum_dt[0] += nA_dt(0);  sum_dt[1] += nA_dt(1);  sum_dt[2] += nA_dt(2);
            //cnt += 1;
        }

        if (cnt == 0) {
            printf("  (no contributing A-side faces found in block [2])\n");
            printf("  stored cohesive_zone_orientation: t=(%.6g, %.6g, %.6g)  tdt=(%.6g, %.6g, %.6g)\n",
                   cohesive_zone_orientation(i,0), cohesive_zone_orientation(i,1), cohesive_zone_orientation(i,2),
                   cohesive_zone_orientation(i,3), cohesive_zone_orientation(i,4), cohesive_zone_orientation(i,5));
            continue;
        }

        // alignment test (exactly like oriented()): flip dt sum if needed
        const double dot_align = sum_t[0]*sum_dt[0] + sum_t[1]*sum_dt[1] + sum_t[2]*sum_dt[2];
        if (dot_align < 0.0) {
            sum_dt[0] = -sum_dt[0];
            sum_dt[1] = -sum_dt[1];
            sum_dt[2] = -sum_dt[2];
        }

        // normalize both sums
        double mag_t  = sqrt(sum_t[0]*sum_t[0] + sum_t[1]*sum_t[1] + sum_t[2]*sum_t[2]);
        double mag_dt = sqrt(sum_dt[0]*sum_dt[0] + sum_dt[1]*sum_dt[1] + sum_dt[2]*sum_dt[2]);

        double current_norm[3] = {0.0,0.0,0.0};
        double next_norm[3] = {0.0,0.0,0.0};
        if (mag_t  > 0.0) { current_norm[0] = sum_t [0]/mag_t;  current_norm[1] = sum_t [1]/mag_t;  current_norm[2] = sum_t [2]/mag_t; }
        if (mag_dt > 0.0) { next_norm[0] = sum_dt[0]/mag_dt; next_norm[1] = sum_dt[1]/mag_dt; next_norm[2] = sum_dt[2]/mag_dt; }
// // THROTTLE OUTPUT COMMENT OUT //
//         printf("  averaged (pre-norm)  t=(%.6g, %.6g, %.6g)  tdt=(%.6g, %.6g, %.6g)  cnt=%d\n",
//                sum_t[0], sum_t[1], sum_t[2], sum_dt[0], sum_dt[1], sum_dt[2], cnt);
//         printf("  align: dot(sum_t, sum_tdt)=%.6g\n", dot_align);
//         printf("  averaged (unit)      t=(%.6g, %.6g, %.6g)  tdt=(%.6g, %.6g, %.6g)\n",
//                current_norm[0], current_norm[1], current_norm[2], next_norm[0], next_norm[1], next_norm[2]);

//         // compare to oriented() output
//         printf("  stored cohesive_zone_orientation:   t=(%.6g, %.6g, %.6g)  tdt=(%.6g, %.6g, %.6g)\n",
//                cohesive_zone_orientation(i,0), cohesive_zone_orientation(i,1), cohesive_zone_orientation(i,2),
//                cohesive_zone_orientation(i,3), cohesive_zone_orientation(i,4), cohesive_zone_orientation(i,5));

//         printf("  diff vs stored:      t=(%.6g, %.6g, %.6g)  tdt=(%.6g, %.6g, %.6g)\n",
//                current_norm[0]-cohesive_zone_orientation(i,0), current_norm[1]-cohesive_zone_orientation(i,1), current_norm[2]-cohesive_zone_orientation(i,2),
//                next_norm[0]-cohesive_zone_orientation(i,3), next_norm[1]-cohesive_zone_orientation(i,4), next_norm[2]-cohesive_zone_orientation(i,5));        
// // THROTTLE OUTPUT COMMENT OUT //
    }
// // THROTTLE OUTPUT COMMENT OUT //
//     printf("\n======================================================\n");
// // THROTTLE OUTPUT COMMENT OUT //
    }
    
} // end cohesive_zones_t::debug_orient()
//COMMENT OUT HERE TO STOP FUNCTION DEBUGGING PRINTS
//    ======================== END cohesive_zone_orientation debug ========================

//    ======================== ucmap debug ========================
//COMMENT OUT HERE TO STOP FUNCTION DEBUGGING PRINTS
void cohesive_zones_t::debug_ucmap(
    const DCArrayKokkos<double>& pos,                 // coords (positions at t)
    const DCArrayKokkos<double>& vel,                 // vel   (velocities at t)
    double dt,                                        // same dt passed to ucmap()
    const CArrayKokkos<double>&  cohesive_zone_orientation,          //
    const CArrayKokkos<size_t>&  overlapping_node_gids,
    const CArrayKokkos<double>&  local_opening              // what ucmap() already wrote
){
// // THROTTLE OUTPUT COMMENT OUT //
//     printf("\n==================== ucmap debug ====================\n");
// // THROTTLE OUTPUT COMMENT OUT //

    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
        if (i != 0) continue; // ONLY DO FIRST PAIR FOR TEST AGAINST TIME
        const size_t A = overlapping_node_gids(i,0);
        const size_t B = overlapping_node_gids(i,1);

        // calculate displacement between overlapping node pairs in global frame at time t
        const double u_rel_x_t = pos(B,0) - pos(A,0);
        const double u_rel_y_t = pos(B,1) - pos(A,1);
        const double u_rel_z_t = pos(B,2) - pos(A,2);

        // calculate velocity difference between overlapping node pairs in global frame at time t
        const double v_rel_x_t = vel(B,0) - vel(A,0);
        const double v_rel_y_t = vel(B,1) - vel(A,1);
        const double v_rel_z_t = vel(B,2) - vel(A,2);

        // calculate displacement at t+dt  (forward euler predicted)
        const double u_x_tdt = u_rel_x_t + dt*v_rel_x_t;
        const double u_y_tdt = u_rel_y_t + dt*v_rel_y_t;
        const double u_z_tdt = u_rel_z_t + dt*v_rel_z_t;

        // normal at time t from oriented() (already unitized)
        const double current_norm_x = cohesive_zone_orientation(i,0);
        const double current_norm_y = cohesive_zone_orientation(i,1);
        const double current_norm_z = cohesive_zone_orientation(i,2);

        // replicating ucmap() exactly

        // dotting with the normal vector to get the normal component of displacement magnitude at time t
        const double u_norm_mag_t = (u_rel_x_t*current_norm_x + u_rel_y_t*current_norm_y + u_rel_z_t*current_norm_z);     

        // tangential components at time t
        const double u_tan_x_t = u_rel_x_t - u_norm_mag_t*current_norm_x;
        const double u_tan_y_t = u_rel_y_t - u_norm_mag_t*current_norm_y;
        const double u_tan_z_t = u_rel_z_t - u_norm_mag_t*current_norm_z;

        // tangential magnitude at time t assuiming that us* == ur*
        const double u_tan_mag_t = sqrt(u_tan_x_t*u_tan_x_t + u_tan_y_t*u_tan_y_t + u_tan_z_t*u_tan_z_t); 

        // velocity components dotted with normal vector to get normal rate of velocity
        const double v_norm_t = v_rel_x_t*current_norm_x + v_rel_y_t*current_norm_y + v_rel_z_t*current_norm_z;    

        // tangential rates of velocity
        const double v_tan_x_t  = v_rel_x_t - v_norm_t*current_norm_x;           
        const double v_tan_y_t  = v_rel_y_t - v_norm_t*current_norm_y;
        const double v_tan_z_t  = v_rel_z_t - v_norm_t*current_norm_z;

        // Forward Euler update
        // normal mangitude at time t+dt
        const double u_norm_mag_tdt = u_norm_mag_t + dt*v_norm_t;

        // tangential components of displacement at time t+dt
        const double u_tan_x_tdt = u_tan_x_t + dt*v_tan_x_t;
        const double u_tan_y_tdt = u_tan_y_t + dt*v_tan_y_t;
        const double u_tan_z_tdt = u_tan_z_t + dt*v_tan_z_t;

        // tangential magnitude of displacement at time t+dt
        const double u_tan_mag_tdt = sqrt(u_tan_x_tdt*u_tan_x_tdt + u_tan_y_tdt*u_tan_y_tdt + u_tan_z_tdt*u_tan_z_tdt);
// // THROTTLE OUTPUT COMMENT OUT //
//         // prints
//         printf("\n-- Pair %zu (A gid=%zu, B gid=%zu) --\n", i, A, B);
//         printf("  current_norm (at t)=(%.6g, %.6g, %.6g), ||current_norm||=%.6g\n",
//                current_norm_x, current_norm_y, current_norm_z, sqrt(current_norm_x*current_norm_x + current_norm_y*current_norm_y + current_norm_z*current_norm_z));
//         printf("  Nodal Displacement at t: u_t=(%.6g, %.6g, %.6g)\n", u_rel_x_t, u_rel_y_t, u_rel_z_t);
//         printf("  Nodal Velocities at t: v_t=(%.6g, %.6g, %.6g)\n", v_rel_x_t, v_rel_y_t, v_rel_z_t);
//         printf("  dt=%.6g\n", dt);
//         //printf("  Forward Euler Method:\n");
//         printf("  Predicted Nodal Displacement at t+dt: u_tdt=(%.6g, %.6g, %.6g)\n",
//                u_x_tdt, u_y_tdt, u_z_tdt);
//         printf("    Normal Crack Opening Magnitude at t: u_norm_mag_t=%.9g\n", u_norm_mag_t);
//         printf("    Tangential Crack Opening Magnitude at t: u_tan_mag_t=%.9g\n", u_tan_mag_t);
//         printf("    -> Forward Euler Predicted Normal Crack Opening at t+dt: u_norm_mag_tdt=%.9g\n", u_norm_mag_tdt);
//         printf("    -> Forward Euler Predicted Tangential Crack Opening Magnitude at t+dt: u_tan_mag_tdt=%.9g\n", u_tan_mag_tdt);
//         printf("  stored local_opening: [u_norm_mag_t=%.9g, u_tan_mag_t=%.9g, u_norm_mag_tdt=%.9g, u_tan_mag_tdt=%.9g]\n",
//                local_opening(i,0), local_opening(i,1), local_opening(i,2), local_opening(i,3));

//         printf("  diff vs stored: "
//                "d_un_t=%.3e d_utan_t=%.3e d_un_tdt=%.3e d_utan_tdt=%.3e\n",
//                u_norm_mag_t - local_opening(i,0),
//                u_tan_mag_t - local_opening(i,1),
//                u_norm_mag_tdt - local_opening(i,2),
//                u_tan_mag_tdt - local_opening(i,3));
//     }

//     printf("\n==============================================================\n");
// // THROTTLE OUTPUT COMMENT OUT //
    
}
}
// ======================== END ucmap debug ========================    

// ======================== cohesive_zone_var_update debug ========================
// ======================== cohesive_zone_var_update DEBUG ========================

KOKKOS_FUNCTION
void cohesive_zones_t::debug_cohesive_zone_var_update(
    const CArrayKokkos<double>& local_opening,
    const double dt,
    const double time_value, // ADDED IN FOR DEBUGGING
    const CArrayKokkos<size_t>& overlapping_node_gids,
    const RaggedRightArrayKokkos<double>& stress_bc_global_vars,
    const int bdy_set,
    const ViewCArrayKokkos<double>& internal_vars,       
    const ViewCArrayKokkos<double>& delta_internal_vars  
) {

    // read cohesive zone input parameters
    const double E_inf   = stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::E_inf);
    const double a1      = stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::a1);
    const double n_exp   = stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::n_exp);
    const double u_n_star= stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::u_n_star);
    const double u_t_star= stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::u_t_star);
    const int    num_prony_terms     = (int)(stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::num_prony_terms) + 0.5);

    // E_dt
    double E_dt = E_inf;
    for (int j = 0; j < num_prony_terms; ++j) {
        const int    base = fractureStressBC::BCVars::prony_base + 2*j;
        const double Ej   = stress_bc_global_vars(bdy_set, base);
        const double tauj = stress_bc_global_vars(bdy_set, base + 1);
        const double tau_eff = (tauj > 0.0) ? tauj : std::numeric_limits<double>::min();
        const double one_minus_exp = 1.0 - exp(-dt / tau_eff);
        E_dt += Ej * tau_eff * (one_minus_exp / dt);
    }
// // THROTTLE OUTPUT COMMENT OUT //
//     printf("================== cohesive_zone_var_update debug ==================\n");
// // THROTTLE OUTPUT COMMENT OUT //

// // THROTTLE OUTPUT COMMENT OUT //
//     printf("cohesive zone input parameters (bdy_set=%d) --\n", bdy_set);
//     printf("  E_inf=%.9g  a1=%.9g  n_exp=%.9g  u_n*=%.9g  u_t*=%.9g  num_prony_terms=%d  dt=%.9g\n",
//            E_inf, a1, n_exp, u_n_star, u_t_star, num_prony_terms, dt);
//     printf("  E_dt=%.9g\n", E_dt);
// // THROTTLE OUTPUT COMMENT OUT //
    // loop over each cohesive zone node pair
    for (size_t i = 0; i < overlapping_node_gids.dims(0); i++){
        // uncomment below for when using excel for alpha, lambda, traction comparison 
        if (i != 0) continue; // ONLY DO FIRST PAIR FOR TEST AGAINST TIME 

        const size_t gidA = overlapping_node_gids(i,0);
        const size_t gidB = overlapping_node_gids(i,1);
// // THROTTLE OUTPUT COMMENT OUT //
//         printf("\n-- Pair %zu (A gid=%zu, B gid=%zu) --\n", i, gidA, gidB);
// // THROTTLE OUTPUT COMMENT OUT //
        // inputs
        const double u_norm_mag_t   = local_opening(i,0);
        const double u_tan_mag_t   = local_opening(i,1);
        const double u_norm_mag_tdt = local_opening(i,2);
        const double u_tan_mag_tdt = local_opening(i,3);

        // recompute lambdas for reporting (mirrors update)
        double lambda_t   = sqrt((u_norm_mag_t   / u_n_star)*(u_norm_mag_t   / u_n_star) + (u_tan_mag_t   / u_t_star)*(u_tan_mag_t   / u_t_star));
        double lambda_tdt = sqrt((u_norm_mag_tdt / u_n_star)*(u_norm_mag_tdt / u_n_star) + (u_tan_mag_tdt / u_t_star)*(u_tan_mag_tdt / u_t_star));
        const double lam_dot = (lambda_tdt - lambda_t) / dt;

        // damage values (before/after clamp) as stored
        const double d_alpha = delta_internal_vars(i,1);
        const double alpha_t = internal_vars(i,1);
        const double alpha_tdt = alpha_t + d_alpha;

        // prony new stresses are stored directly at [4 + j]
        double sigma_sum = 0.0, sigma_sum_exp = 0.0;

// // THROTTLE OUTPUT COMMENT OUT //        printf("  local_opening(t):     u_norm_mag_t=%.9g  u_tan_mag_t=%.9g\n", u_norm_mag_t,   u_tan_mag_t);
//         printf("  local_opening(tdt):  u_norm_mag_tdt=%.9g  u_tan_mag_tdt=%.9g\n", u_norm_mag_tdt, u_tan_mag_tdt);
//         printf("  lambda_t=%.9g  lambda_tdt=%.9g  lambda_dot_t=%.9g  (stored lambda_dot_t=%.9g)\n",
//                 lambda_t, lambda_tdt, lam_dot, delta_internal_vars(i,0));
// // THROTTLE OUTPUT COMMENT OUT //

        // show damage
        double dadt_report;
        if (lam_dot > 0.0) {
            const double lam_mid = 0.5*(lambda_t + lambda_tdt);
            dadt_report = a1 * pow(lam_mid, n_exp);
        } else {
            dadt_report = 0.0;
          }
        //printf("  alpha(t)=%.9g  d_alpha(dt)=%.9g  (expected dadt*dt=%.9g)  alpha(t+dt)=%.9g\n",
                //alpha_t, d_alpha, dadt_report*dt, alpha_tdt);
// // THROTTLE OUTPUT COMMENT OUT //
//         printf("  alpha(t)=%.9g  d_alpha(dt)=%.9g  alpha(t+dt)=%.9g\n",
//                 alpha_t, d_alpha, alpha_tdt);                
// // THROTTLE OUTPUT COMMENT OUT //

// // THROTTLE OUTPUT COMMENT OUT //
//           printf("  alpha(t)=%.9g\n", alpha_t);
// // THROTTLE OUTPUT COMMENT OUT //
        // per Prony term details and sums
        for (int j = 0; j < num_prony_terms; ++j) {
            const int    base   = fractureStressBC::BCVars::prony_base + 2*j;
            const double Ej     = stress_bc_global_vars(bdy_set, base);
            const double tauj   = stress_bc_global_vars(bdy_set, base + 1);
            const double tau_eff= (tauj > 0.0) ? tauj : std::numeric_limits<double>::min();
            const double aexp   = exp(-dt / tau_eff);
            const double sigma_current = internal_vars(i, 4 + j);          
            const double sigma_next = delta_internal_vars(i, 4 + j);
// // THROTTLE OUTPUT COMMENT OUT //    
//             printf("  Prony[%d]: E=%.9g  tau=%.9g  a=exp(-dt/tau)=%.9g  sigma_current=%.9g  sigma_next=%.9g\n",
//                     j, Ej, tauj, aexp, sigma_current, sigma_next);
// // THROTTLE OUTPUT COMMENT OUT //
            sigma_sum   += sigma_next;
            sigma_sum_exp += (1.0 - exp(-dt / tau_eff)) * sigma_next;
        }
// // THROTTLE OUTPUT COMMENT OUT //
//         printf("  sigma_sum=%.9g  sigma_sum_exp=%.9g\n", sigma_sum, sigma_sum_exp);
// // THROTTLE OUTPUT COMMENT OUT //
        // scalar pieces in traction increment formula
        const double deltaE_term  = E_dt * lam_dot * dt;
        const double elastic_term = E_inf * lambda_t + sigma_sum;   
        const double damp_term    = -sigma_sum_exp;

        // inverses from production function
        //const double inv_uns_lambda_tdt = (u_n_star > 0.0 && lambda_tdt > 0.0) ? 1.0 / (u_n_star * lambda_tdt) : 0.0;
        //const double inv_uns_lambda_t   = (u_n_star > 0.0 && lambda_t   > 0.0) ? 1.0 / (u_n_star * lambda_t  ) : 0.0;
        //const double inv_uts_lambda_tdt = (u_t_star > 0.0 && lambda_tdt > 0.0) ? 1.0 / (u_t_star * lambda_tdt) : 0.0;
        //const double inv_uts_lambda_t   = (u_t_star > 0.0 && lambda_t   > 0.0) ? 1.0 / (u_t_star * lambda_t  ) : 0.0;

        // inverses (guarded)
        const double inv_uns_tdt = (u_n_star > 0.0 && lambda_tdt > 0.0) ? 1.0 / (u_n_star * lambda_tdt) : 0.0;
        const double inv_uns_t   = (u_n_star > 0.0 && lambda_t   > 0.0) ? 1.0 / (u_n_star * lambda_t  ) : 0.0;
        const double inv_uts_tdt = (u_t_star > 0.0 && lambda_tdt > 0.0) ? 1.0 / (u_t_star * lambda_tdt) : 0.0;
        const double inv_uts_t   = (u_t_star > 0.0 && lambda_t   > 0.0) ? 1.0 / (u_t_star * lambda_t  ) : 0.0;

// // THROTTLE OUTPUT COMMENT OUT // 
//         printf("  terms: deltaE_term=%.9g  elastic_term=%.9g  damp_term=%.9g\n", deltaE_term, elastic_term, damp_term);
// // THROTTLE OUTPUT COMMENT OUT // 
        // traction increment terms breakdown
        const double one_m_a_tdt = (1.0 - alpha_tdt);
        const double one_m_a_t   = (1.0 - alpha_t);
        const double termN1 = u_norm_mag_tdt * inv_uns_tdt * one_m_a_tdt * deltaE_term; // normal traction viscous term t+dt
        const double termN2 = u_norm_mag_tdt * inv_uns_tdt * one_m_a_tdt * elastic_term; // normal traction elastic term t+dt
        const double termN3 = -u_norm_mag_t   * inv_uns_t   * one_m_a_t   * elastic_term; // normal traction elastic term t
        const double termN4 = u_norm_mag_tdt * inv_uns_tdt * one_m_a_tdt * damp_term; // normal traction damping term t+dt
        const double termT1 = u_tan_mag_tdt * inv_uts_tdt * one_m_a_tdt * deltaE_term; // tangential traction viscous term t+dt
        const double termT2 = u_tan_mag_tdt * inv_uts_tdt * one_m_a_tdt * elastic_term; // tangential traction elastic term t+dt
        const double termT3 = -u_tan_mag_t   * inv_uts_t   * one_m_a_t   * elastic_term; // tangential traction elastic term t
        const double termT4 = u_tan_mag_tdt * inv_uts_tdt * one_m_a_tdt * damp_term; // tangential traction damping term t+dt
// // THROTTLE OUTPUT COMMENT OUT // 
//         printf("  traction_norm terms: [%.9g, %.9g, %.9g, %.9g]\n", termN1, termN2, termN3, termN4);
//         printf("  traction_tan terms: [%.9g, %.9g, %.9g, %.9g]\n", termT1, termT2, termT3, termT4);
// // THROTTLE OUTPUT COMMENT OUT // 

// // THROTTLE OUTPUT COMMENT OUT //
//         printf("  Tn(t) (normal traction at time t) = %.9g\n", fabs(termN3));
// // THROTTLE OUTPUT COMMENT OUT //

// // THROTTLE OUTPUT COMMENT OUT //
//         //printf(" normal traction at time t: Tn_t=%.9g\n", termN3);
//         printf(" normal traction at time t+dt: Tn_tdt=%.9g\n", termN2);
        

    
//         printf("  STORED: lambda_dot_t=%.9g, delta_a=%.9g, traction_norm(normal traction increment)=%.9g, traction_tan(tangential traction increment)=%.9g\n",
//                 delta_internal_vars(i,0), delta_internal_vars(i,1),
//                 delta_internal_vars(i,2), delta_internal_vars(i,3));
// // THROTTLE OUTPUT COMMENT OUT // 

// // THROTTLE OUTPUT COMMENT OUT // 
//         printf("  STORED: lambda_dot_t=%.9g\n", delta_internal_vars(i,0));
//         printf("  STORED: delta_a=%.9g\n", delta_internal_vars(i,1));
//         printf("  STORED: traction_norm(normal increment)=%.9g\n", delta_internal_vars(i,2));
//         printf("  STORED: traction_tan(tangential increment)=%.9g\n", delta_internal_vars(i,3));
// // THROTTLE OUTPUT COMMENT OUT // 
        //if (num_prony_terms > 0) {
//
         //   printf("  STORED : ");
         //   for (int j = 0; j < num_prony_terms; ++j) {
        //        printf("%s%.9g", (j==0?"[":", "), delta_internal_vars(i, 4 + j));
        //    }
        //    printf("]\n");

        // uncomment below for when using excel for alpha, lambda, traction comparison 

// THROTTLE OUTPUT COMMENT OUT //
       // FOR CONSTANT RATE OF OPENING BEFORE COHESIVE ZONE LOADS TEST:
       //CSV STYLE: time, alpha, lambda_t, Tn_t
       if (i == 0){
       printf("%.9e, %.9e, %.9e, %.9e\n",
              time_value,
              alpha_t,
              lambda_t,
              fabs(termN3));
       }
// THROTTLE OUTPUT COMMENT OUT //

// // THROTTLE OUTPUT COMMENT OUT //
//        // FOR CONSTANT RATE OF OPENING AFTER COHESIVE ZONE LOADS CONTRIBUTION TEST:
//        //CSV STYLE: time, alpha, lambda_t, Tn_t
//        if (i == 0){
//            const double Tn_t = internal_vars(i,2); // traction normal at time t after cohesive_zone_loads contribution
//            const double dTn   = delta_internal_vars(i,2); // traction normal increment after cohesive_zone_loads contribution
//            const double Tn_tdt = Tn_t + dTn; // traction normal at time t+dt after cohesive_zone_loads contribution
//         printf("%.9e, %.9e, %.9e, %.9e\n",
//                  time_value,
//                  alpha_t,
//                  lambda_t,
//                  fabs(termN3));
//        }
// // THROTTLE OUTPUT COMMENT OUT //
    }
// // THROTTLE OUTPUT COMMENT OUT //
//     printf("\n================ end cohesive_zone_var_update debug ================\n");
// // THROTTLE OUTPUT COMMENT OUT //
}    
// ======================== END cohesive_zone_var_update debug ======================== 

// ======================== cohesive_zone_loads debug ========================
KOKKOS_FUNCTION
void cohesive_zones_t::debug_cohesive_zone_loads(
    Mesh_t &mesh,
    const DCArrayKokkos<double> &pos,
    const CArrayKokkos<size_t> &overlapping_node_gids,
    const CArrayKokkos<double> &cohesive_zone_orientation,
    //const CArrayKokkos<int> &cz_info,
    //const size_t max_elem_in_cohesive_zone,
    const ViewCArrayKokkos<double> &internal_vars,
    const ViewCArrayKokkos<double> &delta_internal_vars,
    const CArrayKokkos<double> &pair_area,
    const ViewCArrayKokkos<double> &F_cz
)
{

// // THROTTLE OUTPUT COMMENT OUT //
//     printf("\n================= debug_cohesive_zone_loads =================\n");
//     printf(" number of cohesive zone node pairs = %zu\n", overlapping_node_gids.dims(0));
//     //printf(" num_nodes = %zu\n", mesh.num_nodes);
// // THROTTLE OUTPUT COMMENT OUT //  
    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
        if (i != 0) continue; // ONLY DO FIRST PAIR FOR TEST AGAINST TIME

        const size_t gidA = overlapping_node_gids(i, 0);
        const size_t gidB = overlapping_node_gids(i, 1);

        // traction normal at time t
        const double Tn_t = internal_vars(i, 2);

        // traction normal increment (dt)
        const double dTn = delta_internal_vars(i, 2);

        // normal traction at t+dt (traction normal at time t plus increment)
        const double Tn_tdt = Tn_t + dTn;

        // traction tangential at time t
        const double Tt_t = internal_vars(i, 3);

        // traction tangential increment (dt)
        const double dTt = delta_internal_vars(i, 3);

        // traction tangential at t+dt (traction tangential at time t plus increment)
        const double Tt_tdt = Tt_t + dTt;

        // normal traction increment at t+dt
        //const double Tn_tdt = internal_vars(i,2) + delta_internal_vars(i,2);

        // tangential traction increment at t+dt
        //const double Tt_tdt = internal_vars(i,3) + delta_internal_vars(i,3);

        // orientation at t+dt
        const double nx = cohesive_zone_orientation(i, 3);
        const double ny = cohesive_zone_orientation(i, 4);
        const double nz = cohesive_zone_orientation(i, 5);

        // effective cohesive zone node pair area
        const double Aeff = pair_area(i);

        // forces on the cohesive zone node pair
        const double FxA = F_cz(3*gidA    );
        const double FyA = F_cz(3*gidA + 1);
        const double FzA = F_cz(3*gidA + 2);

        const double FxB = F_cz(3*gidB    );
        const double FyB = F_cz(3*gidB + 1);
        const double FzB = F_cz(3*gidB + 2);

        // checking sum of the forces (should be equal and opposite == 0)
        const double Fsum_x = FxA + FxB;
        const double Fsum_y = FyA + FyB;
        const double Fsum_z = FzA + FzB;
// // THROTTLE OUTPUT COMMENT OUT //
//         printf("pair %zu: (gidA=%zu, gidB=%zu)\n", i, gidA, gidB);
//         printf("  area_total = %.15e\n", Aeff);
//         printf("  Tn_t   = %+ .6e,  dTn   = %+ .6e,  Tn_tdt = %+ .6e\n",
//                Tn_t, dTn, Tn_tdt);
//         printf("  Tt_t   = %+ .6e,  dTt   = %+ .6e,  Tt_tdt = %+ .6e\n",
//                Tt_t, dTt, Tt_tdt);
//        //printf("  Tn_tdt = %+ .6e,  Tt_tdt = %+ .6e\n", Tn_tdt, Tt_tdt);        
//         printf("  n_tdt  = ( %+ .6e, %+ .6e, %+ .6e )\n",
//                nx, ny, nz);
//         printf("  F_A    = ( %+ .6e, %+ .6e, %+ .6e )\n",
//                FxA, FyA, FzA);
//         printf("  F_B    = ( %+ .6e, %+ .6e, %+ .6e )\n",
//                FxB, FyB, FzB);
//         printf("  F_A + F_B = ( %+ .6e, %+ .6e, %+ .6e )  [should == 0]\n\n",
//                Fsum_x, Fsum_y, Fsum_z);
// // THROTTLE OUTPUT COMMENT OUT //               
    }
// // THROTTLE OUTPUT COMMENT OUT //
//     printf("==============================================================\n\n");
// // THROTTLE OUTPUT COMMENT OUT //
}
// ======================== END cohesive_zone_loads debug ========================

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn cohesive_zone_elem_count
/// \brief Returns the maximum number of elements connected to any node in the cohesive zone overlapping node pairs
/// This value is used to size data structures that depend on the maximum connectivity per node
/// \param overlapping_node_gids 2D array (num_pairs x 2) containing node pairs involved in cohesive zones
/// \param elems_in_node RaggedRightArray mapping each node to the elements it belongs to
/// \param mesh Reference to the mesh containing connectivity information
/// \return Maximum number of elements connected to any node in any cohesive pair
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t cohesive_zones_t::cohesive_zone_elem_count(const CArrayKokkos<size_t>& overlapping_node_gids,
               const RaggedRightArrayKokkos<size_t>& elems_in_node, const Mesh_t& mesh) {

    size_t max_elem_in_cohesive_zone = 0;
    FOR_REDUCE_MAX(i, 0, overlapping_node_gids.dims(0),
                   j, 0, overlapping_node_gids.dims(1), max_elem_in_cohesive_zone, {
        if (max_elem_in_cohesive_zone < mesh.elems_in_node.stride(overlapping_node_gids(i,j))) {
            max_elem_in_cohesive_zone = mesh.elems_in_node.stride(overlapping_node_gids(i,j));
        }
    }, max_elem_in_cohesive_zone);

    return max_elem_in_cohesive_zone;
}

/// \brief Computes face geometry vectors and centroid for a given element surface
///
/// This function computes the geometric properties of a specified surface (face) — 
/// also referred to as a "patch" per the nodal indexing convention in mesh.h — 
/// for a first-order hexahedral element.
/// Specifically, it calculates the orthonormal in-plane basis vectors r and s, 
/// the outward unit normal vector n, and the centroid cenface of the face in physical space
///
/// \param nodes Global nodal coordinates array (num_nodes x 3) from the mesh
/// \param conn Element-to-node connectivity array (num_elems x nodes_in_elem) from the mesh
/// \param surf Local surface (patch) ID [0–5] corresponding to a face of a hex element 
///             (per the face-node ordering in mesh.h) (which face)
/// \param elem Index of the element from which the surface is extracted (whcihc element)
/// \param n Output normal vector to the face (length 3, unit magnitude)
/// \param r Output in-plane direction vector r (length 3, unit magnitude)
/// \param s Output in-plane direction vector s (length 3, unit magnitude)
/// \param cenface Output centroid of the face in physical space (length 3)
///
/// \note This function assumes first-order hexahedral elements (nodes_in_elem = 8)
///
KOKKOS_FUNCTION
void cohesive_zones_t::compute_face_geometry(const DCArrayKokkos<double> &nodes, // unused
                            const Mesh_t &mesh,
                            const DCArrayKokkos<double> &node_coords,
                            const DCArrayKokkos<size_t> &conn, // unused
                            const size_t surf,
                            const size_t elem,
                            ViewCArrayKokkos<double> &n,
                            ViewCArrayKokkos<double> &r,
                            ViewCArrayKokkos<double> &s,
                            ViewCArrayKokkos<double> &cenface)
                            const {
 
    // building face-to-global node id mapping for HEX8 from mesh.h
    size_t face_gid[4];
    switch (surf) {
        case 0: // [0,4,6,2]
            face_gid[0] = mesh.nodes_in_elem(elem, 0);
            face_gid[1] = mesh.nodes_in_elem(elem, 4);
            face_gid[2] = mesh.nodes_in_elem(elem, 6);
            face_gid[3] = mesh.nodes_in_elem(elem, 2);
            break;
        case 1: // [1,3,7,5]
            face_gid[0] = mesh.nodes_in_elem(elem, 1);
            face_gid[1] = mesh.nodes_in_elem(elem, 3);
            face_gid[2] = mesh.nodes_in_elem(elem, 7);
            face_gid[3] = mesh.nodes_in_elem(elem, 5);
            break;
        case 2: // [0,1,5,4]
            face_gid[0] = mesh.nodes_in_elem(elem, 0);
            face_gid[1] = mesh.nodes_in_elem(elem, 1);
            face_gid[2] = mesh.nodes_in_elem(elem, 5);
            face_gid[3] = mesh.nodes_in_elem(elem, 4);
            break;
        case 3: // [3,2,6,7]
            face_gid[0] = mesh.nodes_in_elem(elem, 3);
            face_gid[1] = mesh.nodes_in_elem(elem, 2);
            face_gid[2] = mesh.nodes_in_elem(elem, 6);
            face_gid[3] = mesh.nodes_in_elem(elem, 7);
            break;
        case 4: // [0,2,3,1]
            face_gid[0] = mesh.nodes_in_elem(elem, 0);
            face_gid[1] = mesh.nodes_in_elem(elem, 2);
            face_gid[2] = mesh.nodes_in_elem(elem, 3);
            face_gid[3] = mesh.nodes_in_elem(elem, 1);
            break;
        case 5: // [4,5,7,6]
            face_gid[0] = mesh.nodes_in_elem(elem, 4);
            face_gid[1] = mesh.nodes_in_elem(elem, 5);
            face_gid[2] = mesh.nodes_in_elem(elem, 7);
            face_gid[3] = mesh.nodes_in_elem(elem, 6);
            break;
        default:
            // shouldn’t happen for HEX8; zero out and return
            for (int j = 0; j < 3; ++j) { n(j)=0; r(j)=0; s(j)=0; cenface(j)=0; }
            return;
    }

    // shape function derivatives at face center
    double xi = 0.0, eta = 0.0;
    double dN_dxi[4], dN_deta[4];

    dN_dxi[0]  = -0.25 * (1.0 - eta);
    dN_dxi[1]  =  0.25 * (1.0 - eta);
    dN_dxi[2]  =  0.25 * (1.0 + eta);
    dN_dxi[3]  = -0.25 * (1.0 + eta);

    dN_deta[0] = -0.25 * (1.0 - xi);
    dN_deta[1] = -0.25 * (1.0 + xi);
    dN_deta[2] =  0.25 * (1.0 + xi);
    dN_deta[3] =  0.25 * (1.0 - xi);

    // zero out accumulators
    // r = (rx, ry, rz)
    // s = (sx, sy, sz)
    // orthogonal in-plane vectors
    double rx = 0.0, ry = 0.0, rz = 0.0;
    double sx = 0.0, sy = 0.0, sz = 0.0;

    for (int j = 0; j < 3; ++j) cenface(j) = 0.0;

    for (int a = 0; a < 4; ++a) {
        size_t node_id = face_gid[a];

        double x = node_coords(node_id, 0);
        double y = node_coords(node_id, 1);
        double z = node_coords(node_id, 2);

        // centroid
        cenface(0) += 0.25 * x;
        cenface(1) += 0.25 * y;
        cenface(2) += 0.25 * z;

        // derivatives
        rx += dN_dxi[a]  * x;
        ry += dN_dxi[a]  * y;
        rz += dN_dxi[a]  * z;

        sx += dN_deta[a] * x;
        sy += dN_deta[a] * y;
        sz += dN_deta[a] * z;
    }

    // normalize r
    double mag_r = sqrt(rx*rx + ry*ry + rz*rz);
    r(0) = rx / mag_r;
    r(1) = ry / mag_r;
    r(2) = rz / mag_r;

    // normalize s
    double mag_s = sqrt(sx*sx + sy*sy + sz*sz);
    s(0) = sx / mag_s;
    s(1) = sy / mag_s;
    s(2) = sz / mag_s;

    // cross product n = r x s
    double nx = r(1)*s(2) - r(2)*s(1);
    double ny = r(2)*s(0) - r(0)*s(2);
    double nz = r(0)*s(1) - r(1)*s(0);

    // normalize n
    double mag_n = sqrt(nx*nx + ny*ny + nz*nz);
    n(0) = nx / mag_n;
    n(1) = ny / mag_n;
    n(2) = nz / mag_n;
                            
    // final cleanup of the -0.0s in the output vectors
    auto zap0 = [](double &v){ if (fabs(v) < 1e-13) v = 0.0; };
    zap0(r(0)); zap0(r(1)); zap0(r(2));
    zap0(s(0)); zap0(s(1)); zap0(s(2));
    zap0(n(0)); zap0(n(1)); zap0(n(2));
}

// this array stores the releveant elements and surfaces for each cohesive zone
// essentially, it makes a map to grab mesh info
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn cohesive_zones_t::build_cohesive_zone_info
/// \brief Build per–cohesive zone node-pair lookup tables of incident elements, corner indices, and the matched (opposing) faces.
///
/// This routine assembles, for every overlapping_node_gids.dims(0) overlapping node pair (A,B), the mesh connectivity that a cohesive zone needs:
/// 1) which elements touch node A and node B, 
/// 2) the local-corner index (k) of A/B inside each such element, and
/// 3) for each element “slot”, which face on the A-side opposes which face on the B-side (centroid coincidence within
///    a tolerance and nearly opposite unit normals). The result is a compact integer table used later to orient and
///    apply cohesive-zone physics.
///
/// \param mesh  Reference to the mesh (element–node connectivity and elems_in_node ragged map are consumed).
/// \param state Reference to the state (node coordinates are used to compute face geometry, normals, and centroids).
/// \param overlapping_node_gids 2D array (num_pairs x 2) of global node IDs, one row per cohesive pair: [A_gid, B_gid].
/// \param max_elem_in_cohesive_zone Upper bound on the number of elements incident to any node in any pair
///                                  (typically from cohesive_zone_elem_count); sizes all per-pair “slot” columns.
/// \param tol Centroid-coincidence tolerance used when declaring two faces to be a match; normals must also be opposite
///            within a dot-product check (dot <= -1 + tol).
/// \return CArrayKokkos<int> cohesive_zone_info table with shape (num_pairs, 6 * max_elem_in_cohesive_zone).
///
/// \details
/// The returned table is organized in 6 contiguous column blocks, each of length max_elem_in_cohesive_zone:
///   [0*max .. 1*max-1] : A-side element IDs          (elements incident to node A), -1 if slot empty
///   [1*max .. 2*max-1] : B-side element IDs          (elements incident to node B), -1 if slot empty
///   [2*max .. 3*max-1] : A-side matched face IDs     (per A element-slot; face index in that A element), -1 if none
///   [3*max .. 4*max-1] : B-side matched face IDs     (per B element-slot; face index in that B element), -1 if none
///   [4*max .. 5*max-1] : kA local-corner indices     (node A’s local corner 0..7 inside each A element), -1 if none
///   [5*max .. 6*max-1] : kB local-corner indices     (node B’s local corner 0..7 inside each B element), -1 if none
///
/// Internally, for each element-slot we derive up to three candidate faces incident to its local corner (k). The face
/// IDs use the code’s hexahedral face numbering convention {0..5}. For each A-slot we search B-slots to find one
/// opposing/coincident face pair (centroids within tol; dot(nA,nB) <= -1+tol). When a match is found, we record:
///   A-face -> block [2] at the same A-slot, and B-face -> block [3] at the same B-slot.
/// This allows multiple distinct A-slots (and B-slots) in a pair to record separate matches (e.g., two CZ faces),
/// while still enforcing “at most one match per slot”.
///
/// Notes:
///  - All outputs are initialized to -1.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CArrayKokkos<int> cohesive_zones_t::build_cohesive_zone_info(
    const Mesh_t& mesh,
    const State_t& state,
    const CArrayKokkos<size_t>& overlapping_node_gids,   // (overlapping_node_gids.dims(0) [number of overlapping node pairs] x 2)
    const size_t max_elem_in_cohesive_zone,              // from cohesive_zone_elem_count()
    const double tol                                      // centroid coincidence tolerance
) {
    // output: (rows = #pairs, cols = 6 * max_elem_in_cohesive_zone)
    // column blocks (each of length max_elem_in_cohesive_zone):
    // [0]   elems A-side: stores elements incident to nodeA (incident meaning all elements that have nodeA in their connectivity)
    // [1]   elems B-side: stores elements incident to nodeB (incident meaning all elements that have nodeB in their connectivity)
    // [2]   matched face ids A-side (filled later)
    // [3]   matched face ids B-side (filled later)
    // [4]   local-corner index in element for nodeA (filled when discover k)
    // [5]   local-corner index in element for nodeB (filled when discover k)
    CArrayKokkos<int> cohesive_zone_info(
        overlapping_node_gids.dims(0),
        6 * max_elem_in_cohesive_zone,
        "cohesive_zone_info"
    );
    cohesive_zone_info.set_values(-1);

    // intermediate faces table (same shape as original vczfaces)
    // for each row i:
    //   slots [0 .. 3*max-1]     : up to 3 faces for each A-side element
    //   slots [3*max .. 6*max-1] : up to 3 faces for each B-side element
    // max 3 faces per element corner
    CArrayKokkos<int> cohesive_zone_faces(
        overlapping_node_gids.dims(0),
        6 * max_elem_in_cohesive_zone,
        "cohesive_zone_faces"
    );
    cohesive_zone_faces.set_values(-1);


    // fill the first two blocks of cohesive_zone_info with incident element lists taken from mesh.elems_in_node (A- and B-side)

    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
        const size_t nodeA = overlapping_node_gids(i, 0);
        const size_t nodeB = overlapping_node_gids(i, 1);

        const size_t degA = mesh.elems_in_node.stride(nodeA);
        for (size_t j = 0; j < max_elem_in_cohesive_zone && j < degA; ++j) {
            cohesive_zone_info(i, 0 + j) = static_cast<int>( mesh.elems_in_node(nodeA, j) );
        }

        const size_t degB = mesh.elems_in_node.stride(nodeB);
        for (size_t j = 0; j < max_elem_in_cohesive_zone && j < degB; ++j) {
            cohesive_zone_info(i, max_elem_in_cohesive_zone + j) = static_cast<int>( mesh.elems_in_node(nodeB, j) );
        }
    }

    // build candidate faces + store local corner k slot-keyed
    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
        // Walk element slots for this pair
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {

            // A-side
            {
                const int elemA = cohesive_zone_info(i, 0 + j);
                if (elemA != -1) {
                    // find local corner kA of nodeA in elemA
                    int kA = -1;
                    for (int k = 0; k < 8; ++k) {
                        if (mesh.nodes_in_elem(static_cast<size_t>(elemA), static_cast<size_t>(k))
                            == overlapping_node_gids(i, 0)) { kA = k; break; }
                    }
                    // store kA slot-keyed in block [4]
                    cohesive_zone_info(i, 4*max_elem_in_cohesive_zone + j) = kA;

                    // store 3 face candidates for A-slot j
                    if (kA >= 0) {
                        switch (kA) { // three faces incident to each local corner k
                            case 0:
                                cohesive_zone_faces(i, 3*j + 0) = 0;
                                cohesive_zone_faces(i, 3*j + 1) = 2;
                                cohesive_zone_faces(i, 3*j + 2) = 4; break;
                            case 1:
                                cohesive_zone_faces(i, 3*j + 0) = 1;
                                cohesive_zone_faces(i, 3*j + 1) = 2;
                                cohesive_zone_faces(i, 3*j + 2) = 4; break;
                            case 2:
                                cohesive_zone_faces(i, 3*j + 0) = 0;
                                cohesive_zone_faces(i, 3*j + 1) = 3;
                                cohesive_zone_faces(i, 3*j + 2) = 4; break;
                            case 3:
                                cohesive_zone_faces(i, 3*j + 0) = 1;
                                cohesive_zone_faces(i, 3*j + 1) = 3;
                                cohesive_zone_faces(i, 3*j + 2) = 4; break;
                            case 4:
                                cohesive_zone_faces(i, 3*j + 0) = 0;
                                cohesive_zone_faces(i, 3*j + 1) = 2;
                                cohesive_zone_faces(i, 3*j + 2) = 5; break;
                            case 5:
                                cohesive_zone_faces(i, 3*j + 0) = 1;
                                cohesive_zone_faces(i, 3*j + 1) = 2;
                                cohesive_zone_faces(i, 3*j + 2) = 5; break;
                            case 6:
                                cohesive_zone_faces(i, 3*j + 0) = 0;
                                cohesive_zone_faces(i, 3*j + 1) = 3;
                                cohesive_zone_faces(i, 3*j + 2) = 5; break;
                            case 7:
                                cohesive_zone_faces(i, 3*j + 0) = 1;
                                cohesive_zone_faces(i, 3*j + 1) = 3;
                                cohesive_zone_faces(i, 3*j + 2) = 5; break;
                            default: break;
                        }
                    }
                }
            }

            // B-side
            {
                const int elemB = cohesive_zone_info(i, max_elem_in_cohesive_zone + j);
                if (elemB != -1) {
                    // find local corner kB of nodeB in elemB
                    int kB = -1;
                    for (int k = 0; k < 8; ++k) {
                        if (mesh.nodes_in_elem(static_cast<size_t>(elemB), static_cast<size_t>(k))
                            == overlapping_node_gids(i, 1)) { kB = k; break; }
                    }
                    // store kB slot-keyed in block [5]
                    cohesive_zone_info(i, 5*max_elem_in_cohesive_zone + j) = kB;

                    // store 3 face candidates for B-slot j (upper half offset = 3*max)
                    const size_t base = 3*max_elem_in_cohesive_zone + 3*j;
                    if (kB >= 0) {
                        switch (kB) { // three faces incident to each local corner k
                            case 0:
                                cohesive_zone_faces(i, base + 0) = 0;
                                cohesive_zone_faces(i, base + 1) = 2;
                                cohesive_zone_faces(i, base + 2) = 4; break;
                            case 1:
                                cohesive_zone_faces(i, base + 0) = 1;
                                cohesive_zone_faces(i, base + 1) = 2;
                                cohesive_zone_faces(i, base + 2) = 4; break;
                            case 2:
                                cohesive_zone_faces(i, base + 0) = 0;
                                cohesive_zone_faces(i, base + 1) = 3;
                                cohesive_zone_faces(i, base + 2) = 4; break;
                            case 3:
                                cohesive_zone_faces(i, base + 0) = 1;
                                cohesive_zone_faces(i, base + 1) = 3;
                                cohesive_zone_faces(i, base + 2) = 4; break;
                            case 4:
                                cohesive_zone_faces(i, base + 0) = 0;
                                cohesive_zone_faces(i, base + 1) = 2;
                                cohesive_zone_faces(i, base + 2) = 5; break;
                            case 5:
                                cohesive_zone_faces(i, base + 0) = 1;
                                cohesive_zone_faces(i, base + 1) = 2;
                                cohesive_zone_faces(i, base + 2) = 5; break;
                            case 6:
                                cohesive_zone_faces(i, base + 0) = 0;
                                cohesive_zone_faces(i, base + 1) = 3;
                                cohesive_zone_faces(i, base + 2) = 5; break;
                            case 7:
                                cohesive_zone_faces(i, base + 0) = 1;
                                cohesive_zone_faces(i, base + 1) = 3;
                                cohesive_zone_faces(i, base + 2) = 5; break;
                            default: break;
                        }
                    }
                }
            }
        }
    }

    // find ALL opposing/coincident face matches (one per element slot)
    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {

    // small stack buffers + Views
    double nA[3], rA[3], sA[3], cA[3];
    double nB[3], rB[3], sB[3], cB[3];
    ViewCArrayKokkos<double> nAj(&nA[0],3), rAj(&rA[0],3), sAj(&sA[0],3), cAj(&cA[0],3);
    ViewCArrayKokkos<double> nBk(&nB[0],3), rBk(&rB[0],3), sBk(&sB[0],3), cBk(&cB[0],3);

    // IDs of three faces per corner k 
    auto push_three_faces = [](int k, int (&out)[3]) {
        switch (k){
        case 0: out[0]=0; out[1]=2; out[2]=4; break;
        case 1: out[0]=1; out[1]=2; out[2]=4; break;
        case 2: out[0]=0; out[1]=3; out[2]=4; break;
        case 3: out[0]=1; out[1]=3; out[2]=4; break;
        case 4: out[0]=0; out[1]=2; out[2]=5; break;
        case 5: out[0]=1; out[1]=2; out[2]=5; break;
        case 6: out[0]=0; out[1]=3; out[2]=5; break;
        case 7: out[0]=1; out[1]=3; out[2]=5; break;
        default: out[0]=out[1]=out[2]=-1; break;
        }
    };

    // loop over A element slots; fill at most one match per A slot
    for (size_t slotA = 0; slotA < max_elem_in_cohesive_zone; ++slotA) {

        // element A
        const int eA = cohesive_zone_info(i, 0*max_elem_in_cohesive_zone + slotA);
        if (eA < 0) continue;

        // skip if this A slot already has a matched face (filled earlier)
        if (cohesive_zone_info(i, 2*max_elem_in_cohesive_zone + slotA) >= 0) continue;

        // corner kA
        const int kA = cohesive_zone_info(i, 4*max_elem_in_cohesive_zone + slotA);
        if (kA < 0) continue;

        // candidate faces for A slot
        int fA_cand[3]; push_three_faces(kA, fA_cand);
        bool matched_this_A_slot = false;

        // search A side
        for (int tA = 0; tA < 3 && !matched_this_A_slot; ++tA) {
            const int fA = fA_cand[tA];
            if (fA < 0) continue;

            // compute A face geometry
            compute_face_geometry(state.node.coords, mesh,
                                  state.node.coords, mesh.nodes_in_elem,
                                  static_cast<size_t>(fA), static_cast<size_t>(eA),
                                  nAj, rAj, sAj, cAj);

            // search B side
            for (size_t slotB = 0; slotB < max_elem_in_cohesive_zone && !matched_this_A_slot; ++slotB) {
                const int eB = cohesive_zone_info(i, 1*max_elem_in_cohesive_zone + slotB);
                if (eB < 0) continue;

                // skip B slot if already filled
                if (cohesive_zone_info(i, 3*max_elem_in_cohesive_zone + slotB) >= 0) continue;

                // corner kB
                const int kB = cohesive_zone_info(i, 5*max_elem_in_cohesive_zone + slotB);
                if (kB < 0) continue;

                // candidate faces for B slot
                int fB_cand[3]; push_three_faces(kB, fB_cand);

                // search B candidate faces
                for (int tB = 0; tB < 3 && !matched_this_A_slot; ++tB) {
                    const int fB = fB_cand[tB];
                    if (fB < 0) continue;

                    // compute B face geometry
                    compute_face_geometry(state.node.coords, mesh,
                                          state.node.coords, mesh.nodes_in_elem,
                                          static_cast<size_t>(fB), static_cast<size_t>(eB),
                                          nBk, rBk, sBk, cBk);

                    // ABS centroid distance + opposite normals
                    const double dx = cAj(0) - cBk(0);
                    const double dy = cAj(1) - cBk(1);
                    const double dz = cAj(2) - cBk(2);
                    const double dist = sqrt(dx*dx + dy*dy + dz*dz);
                    const double dot  = nAj(0)*nBk(0) + nAj(1)*nBk(1) + nAj(2)*nBk(2);

                    // check that match is within tolerance

                    if (dist <= tol && dot <= -1.0 + tol) {

                        // record the match in the following slots
                        cohesive_zone_info(i, 2*max_elem_in_cohesive_zone + slotA) = fA; // A-face for A slot
                        cohesive_zone_info(i, 3*max_elem_in_cohesive_zone + slotB) = fB; // B-face for B slot
                        matched_this_A_slot = true; // done with this A slot; move to next A slot
                    }
                }
            }
        }
    }
}
    return cohesive_zone_info;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn cohesive_zones_t::oriented
/// \brief Average and orient cohesive-zone face normals for each overlapping_node_gids.dims(0) overlapping node pair (reference and current configs).
///
/// For every cohesive zone node pair (A,B), this routine scans the per-pair “slots” produced by build_cohesive_zone_info(),
/// gathers all *matched* A-side faces (block [2]) together with their parent A-elements (block [0]), computes the
/// unit face normals at time t (pos) and at time t+dt (pos), sums them, enforces a consistent sign over time by
/// aligning the t+dt sum to the t sum, normalizes both sums, and writes the result to cohesive_zone_orientation:
///     cohesive_zone_orientation(i,:) = [ n_ref_x, n_ref_y, n_ref_z,  n_cur_x, n_cur_y, n_cur_z ].
/// If a pair has no contributing matched A-faces, the stored orientation remains zero.
///
/// Inputs are assumed to be in the *current mesh connectivity* (hexes with face IDs {0..5}) and with cz_info already
/// populated by build_cohesive_zone_info(). Face geometry (centroid, in-plane directions, and outward unit normal)
/// is computed on the fly via compute_face_geometry() using pos and pos.
///
/// \param mesh  Mesh object (provides nodes_in_elem and element/face topology used by compute_face_geometry).
/// \param pos   Node coordinates at reference time t            (num_nodes x 3).
/// \param pos Node coordinates at current time t + dt         (num_nodes x 3).
/// \param overlapping_node_gids 2D array of cohesive pairs      (num_pairs x 2) with global node IDs [A_gid, B_gid].
/// \param cz_info Integer table from build_cohesive_zone_info()  (num_pairs x 6*max); blocks are:
///                 [0] A-elements per slot, [1] B-elements per slot,
///                 [2] matched A-faces per A-slot, [3] matched B-faces per B-slot,
///                 [4] kA local corner per A-slot, [5] kB local corner per B-slot.
/// \param max_elem_in_cohesive_zone Slot count per pair (same value used to size the cz_info blocks).
/// \param tol Centroid-coincidence tolerance used during matching.
/// \return cohesive_zone_orientation Output (num_pairs x 6): per-pair unit normals at t and t+dt:
///                   columns 0..2 → current_norm (from pos), columns 3..5 → next_norm (from pos).
///
/// \details
/// Algorithm per pair i:
///   1) Initialize sums sum_t = 0, sum_dt = 0, cnt = 0.
///   2) For each A-slot j = 0..max-1:
///        -Read eA = cz_info(i, [0] + j) and fA = cz_info(i, [2] + j).
///        -If both are valid (>= 0), call compute_face_geometry(pos,   eA, fA) → nA_t, and
///                                     compute_face_geometry(pos, eA, fA) → nA_dt.
///        -Accumulate: sum_t  += nA_t;  sum_dt += nA_dt;  ++cnt.
///   3) If cnt == 0 -> leave zeros for this pair and continue.
///   4) Temporal sign consistency: if dot(sum_t, sum_dt) < 0, flip sum_dt = -sum_dt.
///   5) Normalize: current_norm = sum_t / ||sum_t||,  next_norm = sum_dt / ||sum_dt|| (guarded against zero magnitude).
///   6) Store into cohesive_zone_orientation(i,0..5).
///
/// Notes:
///  -This averages *all* matched A-side Cohesive Zone faces recorded for the pair (not just one), so if there are multiple CZ faces
///    between the same A/B regions, both contribute to the average.
///  -The B-side matches are not needed here once the A-side matches have been established; orientation uses A-faces.
///  -compute_face_geometry() is assumed to return outward unit normals consistent with the element’s local face ordering.
///  -If a face degenerates (nearly zero area), its normal magnitude can be ill-defined; after summation, zero-magnitude
///    checks ensure we do not divide by zero.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void cohesive_zones_t::oriented(
    Mesh_t& mesh,
    DCArrayKokkos<double>& pos,      // current  coords (num_nodes x 3)
    //DCArrayKokkos<double>& pos,    // current ("t+dt") coords (num_nodes x 3) 
    CArrayKokkos<size_t>& overlapping_node_gids, // (nvcz x 2): A and B node ids per cohesive pair
    CArrayKokkos<int>& cz_info,      // from build_cohesive_zone_info()
    size_t max_elem_in_cohesive_zone,
    double tol,                 // centroid coincidence tolerance (ABS distance)
    CArrayKokkos<double>& cohesive_zone_orientation       // (overlapping_node_gids.dims(0) x 6): [nx_t,ny_t,nz_t, nx_tdt,ny_tdt,nz_tdt]
) 
{
    // zero out output array
    cohesive_zone_orientation.set_values(0.0);

    // temp views for compute_face_geometry
    double nA_t_buf[3], rA_t_buf[3], sA_t_buf[3], cA_t_buf[3];
    double nB_t_buf[3], rB_t_buf[3], sB_t_buf[3], cB_t_buf[3];
    double nA_dt_buf[3], rA_dt_buf[3], sA_dt_buf[3], cA_dt_buf[3];
    double nB_dt_buf[3], rB_dt_buf[3], sB_dt_buf[3], cB_dt_buf[3];
    ViewCArrayKokkos<double> nA_t (&nA_t_buf[0], 3),  rA_t (&rA_t_buf[0], 3),  sA_t (&sA_t_buf[0], 3),  cA_t (&cA_t_buf[0], 3);
    ViewCArrayKokkos<double> nB_t (&nB_t_buf[0], 3),  rB_t (&rB_t_buf[0], 3),  sB_t (&sB_t_buf[0], 3),  cB_t (&cB_t_buf[0], 3);
    ViewCArrayKokkos<double> nA_dt(&nA_dt_buf[0], 3), rA_dt(&rA_dt_buf[0], 3), sA_dt(&sA_dt_buf[0], 3), cA_dt(&cA_dt_buf[0], 3);
    ViewCArrayKokkos<double> nB_dt(&nB_dt_buf[0], 3), rB_dt(&rB_dt_buf[0], 3), sB_dt(&sB_dt_buf[0], 3), cB_dt(&cB_dt_buf[0], 3);
    
    // pull the single matched faces that build_cohesive_zone_info() wrote:
    // A-faces are in block [2], B-faces are in block [3]
    // A-elems in block [0], B-elems in block [1]
    // basically, find the first non -1 on each side
    // find first true A/B face match (abs centroid distance <= tol and opposite normals)
    // A-side element slots are in block #0; their local corners are in block #4

    // looping through each cohesive zone pair
    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {

        // accumulators for averaving normals of cohesive zone faces
        double sum_t [3] = {0.0, 0.0, 0.0};
        double sum_dt [3] = {0.0, 0.0, 0.0};
        int cnt = 0;

        // find first matched face on A side (block[2]) and B side (block[3]) that is greater than or equal to zero
        // A side
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {
            //const int f_try = cz_info(i, 2*max_elem_in_cohesive_zone + j);
            //if (f_try >= 0) { jA = static_cast<int>(j); fA = f_try; break; }
            const int eA = cz_info(i, 0*max_elem_in_cohesive_zone + j); // A elem slot j
            const int fA = cz_info(i, 2*max_elem_in_cohesive_zone + j); // A matched face slot j
            if (eA < 0 || fA < 0) {
                continue; // no contributing face in this slot
            }            
        
        
        // compute A-side normal at t (reference) and t+dt (current)
        compute_face_geometry(pos,   mesh, pos,   mesh.nodes_in_elem,
                              static_cast<size_t>(fA), static_cast<size_t>(eA),
                              nA_t, rA_t, sA_t, cA_t);
        compute_face_geometry(pos, mesh, pos, mesh.nodes_in_elem,
                              static_cast<size_t>(fA), static_cast<size_t>(eA),
                              nA_dt, rA_dt, sA_dt, cA_dt);

        // optional: compute B at t for a quick sanity check (opposing normals)            
        //compute_face_geometry(pos,   mesh, pos,   mesh.nodes_in_elem,
        //                      static_cast<size_t>(fB), static_cast<size_t>(eB),
        //                      nB_t, rB_t, sB_t, cB_t);



        // accumulate normals
        // reference normal = nA_t (reference)
        // current normal = nA_dt (current)        
        sum_t [0] += nA_t (0);  sum_t [1] += nA_t (1);  sum_t [2] += nA_t (2);
        sum_dt[0] += nA_dt(0);  sum_dt[1] += nA_dt(1);  sum_dt[2] += nA_dt(2);
        cnt += 1;        
        } // end for j

        if (cnt == 0) {
            // no matched A-faces found for this VCZ row; leave zeros
            continue;
        }        

        // align t+dt sum to t sum (keep a consistent sign across time)
        double dot_align = sum_t[0]*sum_dt[0] + sum_t[1]*sum_dt[1] + sum_t[2]*sum_dt[2];
        if (dot_align < 0.0) {
            sum_dt[0] = -sum_dt[0];
            sum_dt[1] = -sum_dt[1];
            sum_dt[2] = -sum_dt[2];
        }

        // normalize
        const double mag_t  = sqrt(sum_t [0]*sum_t [0] + sum_t [1]*sum_t [1] + sum_t [2]*sum_t [2]);
        const double mag_dt = sqrt(sum_dt[0]*sum_dt[0] + sum_dt[1]*sum_dt[1] + sum_dt[2]*sum_dt[2]);

        double current_norm[3] = {0.0,0.0,0.0};
        double next_norm[3] = {0.0,0.0,0.0};
        if (mag_t  > 0.0) { current_norm[0] = sum_t [0]/mag_t;  current_norm[1] = sum_t [1]/mag_t;  current_norm[2] = sum_t [2]/mag_t; }
        if (mag_dt > 0.0) { next_norm[0] = sum_dt[0]/mag_dt; next_norm[1] = sum_dt[1]/mag_dt; next_norm[2] = sum_dt[2]/mag_dt; }

        // store
        cohesive_zone_orientation(i,0) = current_norm[0]; // nx_t (current)
        cohesive_zone_orientation(i,1) = current_norm[1]; // ny_t (current)
        cohesive_zone_orientation(i,2) = current_norm[2]; // nz_t (current)
        cohesive_zone_orientation(i,3) = next_norm[0]; // nx_tdt (next)
        cohesive_zone_orientation(i,4) = next_norm[1]; // ny_tdt (next)
        cohesive_zone_orientation(i,5) = next_norm[2]; // nz_tdt (next)
    }
} // end oriented()


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// NEED TO REWRITE THIS FUNCTION HEADER
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void cohesive_zones_t::ucmap(
    const DCArrayKokkos<double>& pos, // State.node.coords // rename as pos
    const DCArrayKokkos<double>& vel, // State.node.vel //rename as vel
    const CArrayKokkos<double>& cohesive_zone_orientation,
    const CArrayKokkos<size_t>& overlapping_node_gids,
    const double dt, // from sgh_execute.cpp timestep driver
    CArrayKokkos<double>& local_opening    // (overlapping_node_gids.dims(0) x 4): [un_t, utan_t, un_tdt, utan_tdt]
)
{
    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
        const size_t NodeA = overlapping_node_gids(i,0);
        const size_t NodeB = overlapping_node_gids(i,1);

        // calulate displacement between overlapping node pairs in global frame at time t
        const double u_rel_x_t = pos(NodeB,0) - pos(NodeA,0);
        const double u_rel_y_t = pos(NodeB,1) - pos(NodeA,1);
        const double u_rel_z_t = pos(NodeB,2) - pos(NodeA,2);

        // calculate velocity difference between overlapping node pairs in global frame at time t
        const double v_rel_x_t = vel(NodeB,0) - vel(NodeA,0);
        const double v_rel_y_t = vel(NodeB,1) - vel(NodeA,1);
        const double v_rel_z_t = vel(NodeB,2) - vel(NodeA,2);

        // normal at time t from oriented() (already unitized)
        const double current_norm_x = cohesive_zone_orientation(i,0); // current normal at time t x direx
        const double current_norm_y = cohesive_zone_orientation(i,1); // current normal at time t y direc
        const double current_norm_z = cohesive_zone_orientation(i,2); // current normal at time t z direc

        // dotting with the normal vector to get the normal component of displacement at time t
        const double u_norm_mag_t = (u_rel_x_t*current_norm_x + u_rel_y_t*current_norm_y + u_rel_z_t*current_norm_z);

        // tangential components of dispacement at time t
        const double u_tan_x_t = u_rel_x_t - u_norm_mag_t*current_norm_x;
        const double u_tan_y_t = u_rel_y_t - u_norm_mag_t*current_norm_y;
        const double u_tan_z_t = u_rel_z_t - u_norm_mag_t*current_norm_z;

        // tangential magnitude of displacement at time t assuming that us* == ur*
        const double u_tan_mag_t = sqrt(u_tan_x_t*u_tan_x_t + u_tan_y_t*u_tan_y_t + u_tan_z_t*u_tan_z_t);

        // relative velocity (velocity difference) vel components projected onto normal/tangent directions

        // velocity components dotted with normal vector to get normal rate of velocity
        const double v_norm_t = v_rel_x_t*current_norm_x + v_rel_y_t*current_norm_y + v_rel_z_t*current_norm_z;

        // tangential rates of velocity
        const double v_tan_x_t  = v_rel_x_t - v_norm_t*current_norm_x;
        const double v_tan_y_t  = v_rel_y_t - v_norm_t*current_norm_y;
        const double v_tan_z_t  = v_rel_z_t - v_norm_t*current_norm_z;

        // Forward-Euler update 

        // normal: u_norm_mag_tdt = u_norm_mag_t + dt*v_norm_t
        const double u_norm_mag_tdt   = u_norm_mag_t + dt*v_norm_t;

        // tangential: 
        const double u_tan_x_tdt   = u_tan_x_t + dt*v_tan_x_t;
        const double u_tan_y_tdt   = u_tan_y_t + dt*v_tan_y_t;
        const double u_tan_z_tdt   = u_tan_z_t + dt*v_tan_z_t;

        // tangential magnitude at time t+dt
        const double u_tan_mag_tdt = sqrt(u_tan_x_tdt*u_tan_x_tdt + u_tan_y_tdt*u_tan_y_tdt + u_tan_z_tdt*u_tan_z_tdt);
        
        // store
        local_opening(i,0) = u_norm_mag_t; // normal crack opening magnitude at time t
        local_opening(i,1) = u_tan_mag_t; // tangential crack opening magnitude at time t
        local_opening(i,2) = u_norm_mag_tdt; // forward eueler predicted normal crack opening magnitude at time t+dt
        local_opening(i,3) = u_tan_mag_tdt; // forward euler predicted tangential crack opening magnitude at time t+dt
    }
}

KOKKOS_FUNCTION
void cohesive_zones_t::cohesive_zone_var_update(
    const CArrayKokkos<double>& local_opening,
    const double dt,
    const double time_value, // ADDED IN FOR DEBUGGING
    const CArrayKokkos<size_t>& overlapping_node_gids,
    const RaggedRightArrayKokkos<double>& stress_bc_global_vars, // BC parameters per boundary set for fractureStressBC
    const int bdy_set,
    const ViewCArrayKokkos<double>& internal_vars,      // (overlapping_node_gids.dims(0), 4 + num_prony_terms)
    const ViewCArrayKokkos<double>& delta_internal_vars // (overlapping_node_gids.dims(0), 4 + num_prony_terms) 
                                                        // lambda_dot_t, d_alpha
)
{
    // read cohesive zone parameters from stress_bc_global_vars for this boundary set
    const double E_inf = stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::E_inf);
    const double a1   = stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::a1);
    const double n_exp = stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::n_exp);
    const double u_n_star  = stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::u_n_star);
    const double u_t_star  = stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::u_t_star);
    const int    num_prony_terms  = (int)(stress_bc_global_vars(bdy_set, fractureStressBC::BCVars::num_prony_terms) + 0.5); // 0.5 for rounding to the nearest int

    // calculating E_dt
    // for j: E_dt += E_j * tau_j * (1 - exp(-dt/tau_j)) / dt
    double E_dt = E_inf;
    for (int j = 0; j < num_prony_terms; ++j) {
        const int prony_base = fractureStressBC::BCVars::prony_base + 2*j;
        const double E_j  = stress_bc_global_vars(bdy_set, prony_base);
        const double tau_j = stress_bc_global_vars(bdy_set, prony_base + 1);
        const double tau_eff       = (tau_j > 0.0) ? tau_j : std::numeric_limits<double>::min(); // same logic as Gavin's code to avoid div by zero
        const double one_minus_exp = 1.0 - exp(-dt / tau_eff);
        E_dt += E_j * tau_eff * (one_minus_exp / dt);
    }

    // loop over each cohesive zone node pair
    for (size_t i = 0; i < overlapping_node_gids.dims(0); i++){

        // reading in local openings (normal and tangential displacements) at t and t+dt
        const double u_norm_mag_t = local_opening(i,0);
        const double u_tan_mag_t = local_opening(i,1);
        const double u_norm_mag_tdt = local_opening(i,2);
        const double u_tan_mag_tdt = local_opening(i,3);

        // calculating lambda_t and lambda _tdt values
        double lambda_t = sqrt((u_norm_mag_t / u_n_star) * (u_norm_mag_t / u_n_star) + (u_tan_mag_t / u_t_star) * (u_tan_mag_t / u_t_star));
        double lambda_tdt = sqrt((u_norm_mag_tdt / u_n_star) * (u_norm_mag_tdt / u_n_star) + (u_tan_mag_tdt / u_t_star) * (u_tan_mag_tdt / u_t_star));

        // calculating lambda_dot_t
        const double lambda_dot_t = (lambda_tdt - lambda_t) / dt;
        delta_internal_vars(i,0) = lambda_dot_t; // lambda rate at t

        // d_alpha_dt (damage growth/increment) over this step
        double d_alpha_dt;
        if (lambda_dot_t > 0.0){
            const double lambda_mid = 0.5*(lambda_tdt + lambda_t); // forward euler mid-point; growth when loading > 0
            d_alpha_dt = a1 * pow(lambda_mid, n_exp);
        }else {
            d_alpha_dt = 0.0;
        } 
        delta_internal_vars(i,1) = d_alpha_dt * dt; // damage variable, increment over the step

        // updating delta prony stresses for prony terms
       for (int j = 0; j < num_prony_terms; ++j) {
            const int    prony_base    = fractureStressBC::BCVars::prony_base + 2*j;
            const double E_j     = stress_bc_global_vars(bdy_set, prony_base); // in Gavin's code, this is Eandrhom(j,0)
            const double tau_j   = stress_bc_global_vars(bdy_set, prony_base + 1); // in Gavin's code, this is Eandrhom(j,1)
            const double tau_eff = (tau_j > 0.0) ? tau_j : std::numeric_limits<double>::min(); // same logic as Gavin's code to avoid div by zero
            const double a       = exp(-dt / tau_eff);
            delta_internal_vars(i, 4 + j) = a * internal_vars(i, 4 + j) + E_j * tau_eff * lambda_dot_t * (1.0 - a); // prony branch stresses 4 columns 
        }

        // calculating sigma sums and sigma product sums (deltaE_term in the residual traction)
        double sigma_sum = 0.0;
        double sigma_sum_exp = 0.0;
        for (int j = 0; j < num_prony_terms; ++j) {
            const int    prony_base    = fractureStressBC::BCVars::prony_base + 2*j;
            const double tau_j   = stress_bc_global_vars(bdy_set, prony_base + 1); // in Gavin's code, this is Eandrhom(j,1)
            const double tau_eff = (tau_j > 0.0) ? tau_j : std::numeric_limits<double>::min(); // same logic as Gavin's code to avoid div by zero
            const double sigma_j = delta_internal_vars(i, 4 + j); // used to update prony stresses
            sigma_sum     += sigma_j;
            sigma_sum_exp += (1.0 - exp(-dt / tau_eff)) * sigma_j;
        }

        // enforcing alpha domain limitations (clamp to 0 or 1)
        const double alpha_t  = internal_vars(i,1); // damage at time t beginning of step
        double       delta_a  = delta_internal_vars(i,1);
        if (alpha_t + delta_a > 1.0) {
            delta_a = 1.0 - alpha_t;
            delta_internal_vars(i,1) = delta_a;
        }
        // damage at the end of the step
        const double alpha_tdt = alpha_t + delta_a;

        // tractions at t+dt
        // scalar terms
        const double deltaE_term  = E_dt * lambda_dot_t * dt;     
        const double elastic_term = E_inf * lambda_t + sigma_sum;  
        const double damp_term    = -sigma_sum_exp;

        // inverses gaurded
        const double inv_uns_lambda_tdt = (u_n_star > 0.0 && lambda_tdt > 0.0) ? 1.0 / (u_n_star * lambda_tdt) : 0.0; // if lambda = 0, set inv to 0 to avoid div by zero (lambda == 0, traction == 0)
        const double inv_uns_lambda_t   = (u_n_star > 0.0 && lambda_t   > 0.0) ? 1.0 / (u_n_star * lambda_t  ) : 0.0;
        const double inv_uts_lambda_tdt = (u_t_star > 0.0 && lambda_tdt > 0.0) ? 1.0 / (u_t_star * lambda_tdt) : 0.0;
        const double inv_uts_lambda_t   = (u_t_star > 0.0 && lambda_t   > 0.0) ? 1.0 / (u_t_star * lambda_t  ) : 0.0;



        delta_internal_vars(i,2) = // normal traction increment
              u_norm_mag_tdt * inv_uns_lambda_tdt * (1.0 - alpha_tdt) * deltaE_term
            + u_norm_mag_tdt * inv_uns_lambda_tdt * (1.0 - alpha_tdt) * elastic_term
            - u_norm_mag_t   * inv_uns_lambda_t   * (1.0 - alpha_t  ) * elastic_term
            + u_norm_mag_tdt * inv_uns_lambda_tdt * (1.0 - alpha_tdt) * damp_term;

        delta_internal_vars(i,3) = // tangential traction increment
              u_tan_mag_tdt * inv_uts_lambda_tdt * (1.0 - alpha_tdt) * deltaE_term
            + u_tan_mag_tdt * inv_uts_lambda_tdt * (1.0 - alpha_tdt) * elastic_term
            - u_tan_mag_t   * inv_uts_lambda_t   * (1.0 - alpha_t  ) * elastic_term
            + u_tan_mag_tdt * inv_uts_lambda_tdt * (1.0 - alpha_tdt) * damp_term;

        // everything stored:
        // delta_internal_vars(i,0) : lambda_dot_t
        // delta_internal_vars(i,1) : delta_a
        // delta_internal_vars(i,2) : normal traction increment
        // delta_internal_vars(i,3) : tangential traction increment
        // delta_internal_vars(i, 4 + j) : prony internal variables 
    }
}        

KOKKOS_FUNCTION
void cohesive_zones_t::cohesive_zone_loads(
    Mesh_t &mesh,
    const DCArrayKokkos<double> &pos,
    const CArrayKokkos<size_t> &overlapping_node_gids,
    const CArrayKokkos<double> &cohesive_zone_orientation,
    const CArrayKokkos<int> &cz_info,
    const size_t max_elem_in_cohesive_zone,
    const ViewCArrayKokkos<double> &internal_vars,
    const ViewCArrayKokkos<double> &delta_internal_vars,
    CArrayKokkos<double> &pair_area,
    const ViewCArrayKokkos<double> &F_cz  
)
{
    // zero out the cohesive zone force vector
    for (size_t a = 0; a < F_cz.dims(0); a++){
        F_cz(a) = 0.0;
    }

    // temp views for geometry (temp storage to return geometry from compute_face_geometry)
    // plain C array wrapped in ViewCArrayKokkos instead of allocating a new CArrayKokkos each time
    double nA_buf[3], rA_buf[3], sA_buf[3], cA_buf[3];
    ViewCArrayKokkos<double> nA(&nA_buf[0],3); // face normal vector
    ViewCArrayKokkos<double> rA(&rA_buf[0],3); // orthogonal vector in face plane
    ViewCArrayKokkos<double> sA(&sA_buf[0],3); // orthogonal vector in face plane
    ViewCArrayKokkos<double> cA(&cA_buf[0],3); // face centroid

    // looping over cohesive zone node pairs 
    for (size_t i = 0; i < overlapping_node_gids.dims(0); i++){

        // global node IDs for the cohesive zone node pairs
        const size_t gidA = overlapping_node_gids(i,0);
        const size_t gidB = overlapping_node_gids(i,1);

        // tractions at t+dt (normal and tangential)
        const double Tn_tdt = internal_vars(i,2) + delta_internal_vars(i,2);  // normal traction
        const double Tt_tdt = internal_vars(i,3) + delta_internal_vars(i,3);  // tangential traction

        // normal direction from cohesive zone orientation
        const double nx = cohesive_zone_orientation(i,3);
        const double ny = cohesive_zone_orientation(i,4);
        const double nz = cohesive_zone_orientation(i,5);

        // effective area for this cohesive zone node pair
        double area_total = 0.0;

        // from cohesive zone info, loop over elements connected to this cohesive zone node pair
        for (size_t j = 0; j < max_elem_in_cohesive_zone; ++j) {
            const int eA = cz_info(i, 0*max_elem_in_cohesive_zone + j); // A elem id
            const int fA = cz_info(i, 2*max_elem_in_cohesive_zone + j); // A face id

            // guard: if no element/face, skip
            if (eA < 0 || fA < 0) continue;

            // gp = 1/sqrt(3) for 2-point gauss quadrature for integrating over a quad HEX8 face
            const double gp = 0.5773502691896257;

            // building 4 GPs on the patch fA (Fierro patch numbering)
            double gps[4][3]; // gps[k][0] = xi
                              // gps[k][1] = eta
                              // gps[k][2] = zeta

            switch (fA) {
            case 0: // xi-minus face
                gps[0][0] = -1.0; gps[0][1] = -gp; gps[0][2] = -gp;
                gps[1][0] = -1.0; gps[1][1] =  gp; gps[1][2] = -gp;
                gps[2][0] = -1.0; gps[2][1] = -gp; gps[2][2] =  gp;
                gps[3][0] = -1.0; gps[3][1] =  gp; gps[3][2] =  gp;
                break;
            case 1: // xi-plus face
                gps[0][0] =  1.0; gps[0][1] = -gp; gps[0][2] = -gp;
                gps[1][0] =  1.0; gps[1][1] =  gp; gps[1][2] = -gp;
                gps[2][0] =  1.0; gps[2][1] = -gp; gps[2][2] =  gp;
                gps[3][0] =  1.0; gps[3][1] =  gp; gps[3][2] =  gp;
                break;
            case 2: // eta-minus face
                gps[0][1] = -1.0; gps[0][0] = -gp; gps[0][2] = -gp;
                gps[1][1] = -1.0; gps[1][0] =  gp; gps[1][2] = -gp;
                gps[2][1] = -1.0; gps[2][0] = -gp; gps[2][2] =  gp;
                gps[3][1] = -1.0; gps[3][0] =  gp; gps[3][2] =  gp;
                break;
            case 3: // eta-plus face
                gps[0][1] =  1.0; gps[0][0] = -gp; gps[0][2] = -gp;
                gps[1][1] =  1.0; gps[1][0] =  gp; gps[1][2] = -gp;
                gps[2][1] =  1.0; gps[2][0] = -gp; gps[2][2] =  gp;
                gps[3][1] =  1.0; gps[3][0] =  gp; gps[3][2] =  gp;
                break;
            case 4: // zeta-minus face
                gps[0][2] = -1.0; gps[0][0] = -gp; gps[0][1] = -gp;
                gps[1][2] = -1.0; gps[1][0] =  gp; gps[1][1] = -gp;
                gps[2][2] = -1.0; gps[2][0] = -gp; gps[2][1] =  gp;
                gps[3][2] = -1.0; gps[3][0] =  gp; gps[3][1] =  gp;
                break;
            case 5: // zeta-plus face
                gps[0][2] =  1.0; gps[0][0] = -gp; gps[0][1] = -gp;
                gps[1][2] =  1.0; gps[1][0] =  gp; gps[1][1] = -gp;
                gps[2][2] =  1.0; gps[2][0] = -gp; gps[2][1] =  gp;
                gps[3][2] =  1.0; gps[3][0] =  gp; gps[3][1] =  gp;
                break;
            default:
                // bad face id
                break;
            }

            // local array for element nodal coords (pos) 8 nodes x 3 coords
            double x[8][3];
            for (int a = 0; a < 8; ++a) {
                // element connectivity
                x[a][0] = pos(mesh.nodes_in_elem((size_t)eA, (size_t)a),0);
                x[a][1] = pos(mesh.nodes_in_elem((size_t)eA, (size_t)a),1);
                x[a][2] = pos(mesh.nodes_in_elem((size_t)eA, (size_t)a),2);
            }

            // signs at each node for Fierro HEX8 ordering:
            // node 0(-,-,-) node 1(+,-,-) node 2(-,+,-) node 3(+,+,-) node 4(-,-,+) node 5(+,-,+) node 6(-,+,+) node 7(+,+,+)
            const double sign_xi[8] = {-1, +1, -1, +1, -1, +1, -1, +1};
            const double sign_eta[8] = {-1, -1, +1, +1, -1, -1, +1, +1};
            const double sign_zeta[8] = {-1, -1, -1, -1, +1, +1, +1, +1};

            double area_face = 0.0;

            // loop over 4 surface gauss points to compute area
            for (int k = 0; k < 4; ++k) {
                const double xi   = gps[k][0];
                const double eta  = gps[k][1];
                const double zeta = gps[k][2];

                // initialize jacobian to zero
                // jacobian J(m,o): m=0(xi),1(eta),2(zeta); o=0(x),1(y),2(z)
                double J[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

                // compute jacobian by summing over element nodes
                for (int a = 0; a < 8; ++a) {
                    // shape function derivatives
                    const double dN_dxi   = 0.125 * sign_xi[a] * (1.0 + sign_eta[a]*eta)  * (1.0 + sign_zeta[a]*zeta);
                    const double dN_deta  = 0.125 * sign_eta[a] * (1.0 + sign_xi[a]*xi)   * (1.0 + sign_zeta[a]*zeta);
                    const double dN_dzeta = 0.125 * sign_zeta[a] * (1.0 + sign_xi[a]*xi)   * (1.0 + sign_eta[a]*eta);

                    for (int o = 0; o < 3; ++o) {
                    J[0][o] += x[a][o] * dN_dxi;
                    J[1][o] += x[a][o] * dN_deta;
                    J[2][o] += x[a][o] * dN_dzeta;
                    }
                }

                // choosing the two surface tangents (rows of J) based on which param is fixed
                int a_row = 0, b_row = 1;
                if (fA == 0 || fA == 1) { a_row = 1; b_row = 2; } // xi fixed = eta,zeta
                if (fA == 2 || fA == 3) { a_row = 0; b_row = 2; } // eta fixed = xi,zeta
                if (fA == 4 || fA == 5) { a_row = 0; b_row = 1; } // zeta fixed = xi,eta

                // extract tangent vectors 1 and 2 (tangent vectors on the face)
                const double tan_1_x = J[a_row][0], tan_1_y = J[a_row][1], tan_1_z = J[a_row][2];
                const double tan_2_x = J[b_row][0], tan_2_y = J[b_row][1], tan_2_z = J[b_row][2];

                // crossing tangent vectors to get normal vector (normal to the surface at that gauss point)
                const double cross_x = tan_1_y*tan_2_z - tan_1_z*tan_2_y;
                const double cross_y = tan_1_z*tan_2_x - tan_1_x*tan_2_z;
                const double cross_z = tan_1_x*tan_2_y - tan_1_y*tan_2_x;

                // are contribution from this gauss point
                area_face += sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z); 
            }

            // area lumping because of 4 corners of face (.25 each corner)
            area_total += 0.25 * area_face;

        }

        // store effective area for this cohesive zone node pair (debugging)
        pair_area(i) = area_total;

        // guard: if there are no faces contributing to area, skip
        if (area_total <= 0.0) continue;

        // normal force vectors
        double Fn_x = Tn_tdt * area_total * nx;
        double Fn_y = Tn_tdt * area_total * ny;
        double Fn_z = Tn_tdt * area_total * nz;

        // tangential unit vector for B-A seperation
        const double dx = pos(gidB,0) - pos(gidA,0);
        const double dy = pos(gidB,1) - pos(gidA,1);
        const double dz = pos(gidB,2) - pos(gidA,2);

        // component of sepration in normal direction
        const double udotn = dx*nx + dy*ny + dz*nz;

        // project out normal component to get tangential component
        double tx = dx - udotn * nx;
        double ty = dy - udotn * ny;
        double tz = dz - udotn * nz;

        // normalize tangential vector
        double tmag = sqrt(tx*tx + ty*ty + tz*tz);

        // if tmag > 0 divide to get unit vector
        if (tmag > 0.0) {
            tx /= tmag;
            ty /= tmag;
            tz /= tmag;

        } else {
            // if no tangential opening, just skip tangential traction
            tx = ty = tz = 0.0;
        }

        // tangential force vectors
        double Ft_x = Tt_tdt * area_total * tx;
        double Ft_y = Tt_tdt * area_total * ty;
        double Ft_z = Tt_tdt * area_total * tz;

        // total cohesive zone force on the pair (equal and opposite)
        const double Fx = Fn_x + Ft_x;
        const double Fy = Fn_y + Ft_y;
        const double Fz = Fn_z + Ft_z;

        // global scale; apply as equal and opposite nodal forces
        F_cz(3*gidA    ) += Fx;
        F_cz(3*gidA + 1) += Fy;
        F_cz(3*gidA + 2) += Fz;
        F_cz(3*gidB    ) -= Fx;
        F_cz(3*gidB + 1) -= Fy;
        F_cz(3*gidB + 2) -= Fz;
    }
}
