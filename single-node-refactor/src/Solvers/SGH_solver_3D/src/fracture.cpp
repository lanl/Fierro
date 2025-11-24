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

    //     auto push_three_faces = [](int k, int (&out)[3]) {
    //         // three faces incident to each local corner k  (matching the mapping used in build_cohesive_zone_info)
    //         switch (k) {
    //             case 0: out[0]=0; out[1]=2; out[2]=4; break;
    //             case 1: out[0]=1; out[1]=2; out[2]=4; break;
    //             case 2: out[0]=0; out[1]=3; out[2]=4; break;
    //             case 3: out[0]=1; out[1]=3; out[2]=4; break;
    //             case 4: out[0]=0; out[1]=2; out[2]=5; break;
    //             case 5: out[0]=1; out[1]=2; out[2]=5; break;
    //             case 6: out[0]=0; out[1]=3; out[2]=5; break;
    //             case 7: out[0]=1; out[1]=3; out[2]=5; break;
    //             default: out[0]=out[1]=out[2]=-1; break;
    //         }
    //     };

    //     // find the local corner index k of a given global node in element e (or -1)
    //     auto find_k = [&](int e, size_t gid)->int {
    //         if (e < 0) return -1;
    //         for (int k = 0; k < 8; ++k) {
    //             if (mesh.nodes_in_elem(static_cast<size_t>(e), static_cast<size_t>(k)) == gid) return k;
    //         }
    //         return -1;
    //     };

    //     bool found = false;
    //     int eA_hit=-1, fA_hit=-1, eB_hit=-1, fB_hit=-1;
    //     double dist_hit = 0.0, dot_hit = 0.0;

    //     // build and test candidates exactly as in the face matcher
    //     for (size_t slotA = 0; slotA < max_elem_in_cohesive_zone && !found; ++slotA) {
    //         const int eA = cz_info(i, 0*max_elem_in_cohesive_zone + slotA);
    //         if (eA < 0) continue;

    //         const int kA = find_k(eA, nodeA);
    //         if (kA < 0) continue;

    //         int fA_cand[3]; push_three_faces(kA, fA_cand);

    //         for (int tA = 0; tA < 3 && !found; ++tA) {
    //             const int fA = fA_cand[tA];
    //             if (fA < 0) continue;

    //             // geometry for A face
    //             compute_face_geometry(State.node.coords, mesh,
    //                                   State.node.coords, mesh.nodes_in_elem,
    //                                   static_cast<size_t>(fA), static_cast<size_t>(eA),
    //                                   nA, rA, sA, cA);

    //             for (size_t slotB = 0; slotB < max_elem_in_cohesive_zone && !found; ++slotB) {
    //                 const int eB = cz_info(i, 1*max_elem_in_cohesive_zone + slotB);
    //                 if (eB < 0) continue;

    //                 const int kB = find_k(eB, nodeB);
    //                 if (kB < 0) continue;

    //                 int fB_cand[3]; push_three_faces(kB, fB_cand);

    //                 for (int tB = 0; tB < 3 && !found; ++tB) {
    //                     const int fB = fB_cand[tB];
    //                     if (fB < 0) continue;

    //                     // geometry for B face
    //                     compute_face_geometry(State.node.coords, mesh,
    //                                           State.node.coords, mesh.nodes_in_elem,
    //                                           static_cast<size_t>(fB), static_cast<size_t>(eB),
    //                                           nB, rB, sB, cB);

    //                     // ABS centroid distance + opposite normals
    //                     const double dx = cA(0) - cB(0);
    //                     const double dy = cA(1) - cB(1);
    //                     const double dz = cA(2) - cB(2);
    //                     const double dist = sqrt(dx*dx + dy*dy + dz*dz);
    //                     const double dot  = nA(0)*nB(0) + nA(1)*nB(1) + nA(2)*nB(2);

    //                     if (dist <= tol && dot <= -1.0 + tol) {
    //                         found = true;
    //                         eA_hit = eA; fA_hit = fA;
    //                         eB_hit = eB; fB_hit = fB;
    //                         dist_hit = dist; dot_hit = dot;
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     if (!found) {
    //         printf("  No matched faces found.\n");
    //     } else {
    //         printf("  Matched faces (first/true): A(elem=%d, face=%d)  B(elem=%d, face=%d)\n",
    //                eA_hit, fA_hit, eB_hit, fB_hit);
    //         printf("    centroid(A)=(%.6g, %.6g, %.6g)  nA=(%.6g, %.6g, %.6g)\n",
    //                cA(0), cA(1), cA(2), nA(0), nA(1), nA(2));
    //         printf("    centroid(B)=(%.6g, %.6g, %.6g)  nB=(%.6g, %.6g, %.6g)\n",
    //                cB(0), cB(1), cB(2), nB(0), nB(1), nB(2));
    //         printf("    checks: |dcentroid|=%.6g  (tol=%.6g)   dot(nA,nB)=%.6g\n",
    //                dist_hit, tol, dot_hit);
    //     }
    // }
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
    printf("[CZ::init] this=%p  pairs=%zu  maxcz=%zu  info_rows=%zu\n",
           (void*)this,
            overlapping_node_gids.dims(0),
            max_elem_in_cohesive_zone,
            cz_info.dims(0));
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
    printf("[CZ::debug] this=%p  pairs(in)=%zu  maxcz(in)=%zu  info_rows(in)=%zu\n",
          (void*)this, overlap.dims(0), maxcz, info.dims(0));
    // end debug oriented


    // call oriented() to compute cohesive zone normals at t and t+dt 
    CArrayKokkos<double> cohesive_zone_orientation(overlap.dims(0), 6, "cohesive_zone_orientation");

    // initialize to zero
    //cohesive_zone_orientation.set_values(0.0);
    //overlapping_node_gids = CArrayKokkos<size_t> (overlapping_node_gids.dims(0), 2, "overlapping_node_gids");
    //cz_info = CArrayKokkos<int> cz_info;
    //tol = double tol;

{
    printf("\n================== cohesive_zone_orientation debug ==================\n");

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
        const size_t nodeA = overlap(i,0);
        const size_t nodeB = overlap(i,1);

        printf("\n-- Pair %zu  (A gid=%zu, B gid=%zu) --\n", i, nodeA, nodeB);

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
        printf("  contributors (A-side matched faces over all slots):\n");  
        
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

            // per-face print
            printf("  A[j=%zu]: eA elem=%d fA face=%d  kA local corner=%d "
                   "cen_t=(%.6g, %.6g, %.6g) n_t=(%.6g, %.6g, %.6g)  |  "
                   "cen_tdt=(%.6g, %.6g, %.6g) n_tdt=(%.6g, %.6g, %.6g)\n",
                   j, eA, fA, kA,
                   cA_t(0),  cA_t(1),  cA_t(2),  nA_t(0),  nA_t(1),  nA_t(2),
                   cA_dt(0), cA_dt(1), cA_dt(2), nA_dt(0), nA_dt(1), nA_dt(2));

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

        printf("  averaged (pre-norm)  t=(%.6g, %.6g, %.6g)  tdt=(%.6g, %.6g, %.6g)  cnt=%d\n",
               sum_t[0], sum_t[1], sum_t[2], sum_dt[0], sum_dt[1], sum_dt[2], cnt);
        printf("  align: dot(sum_t, sum_tdt)=%.6g\n", dot_align);
        printf("  averaged (unit)      t=(%.6g, %.6g, %.6g)  tdt=(%.6g, %.6g, %.6g)\n",
               current_norm[0], current_norm[1], current_norm[2], next_norm[0], next_norm[1], next_norm[2]);

        // compare to oriented() output
        printf("  stored cohesive_zone_orientation:   t=(%.6g, %.6g, %.6g)  tdt=(%.6g, %.6g, %.6g)\n",
               cohesive_zone_orientation(i,0), cohesive_zone_orientation(i,1), cohesive_zone_orientation(i,2),
               cohesive_zone_orientation(i,3), cohesive_zone_orientation(i,4), cohesive_zone_orientation(i,5));

        printf("  diff vs stored:      t=(%.6g, %.6g, %.6g)  tdt=(%.6g, %.6g, %.6g)\n",
               current_norm[0]-cohesive_zone_orientation(i,0), current_norm[1]-cohesive_zone_orientation(i,1), current_norm[2]-cohesive_zone_orientation(i,2),
               next_norm[0]-cohesive_zone_orientation(i,3), next_norm[1]-cohesive_zone_orientation(i,4), next_norm[2]-cohesive_zone_orientation(i,5));        

    }
    
    printf("\n======================================================\n");
    }
    
} // end cohesive_zones_t::debug_orient()
// COMMENT OUT HERE TO STOP FUNCTION DEBUGGING PRINTS
    // ======================== END cohesive_zone_orientation debug ========================

    // ======================== ucmap debug ========================
// COMMENT OUT HERE TO STOP FUNCTION DEBUGGING PRINTS
void cohesive_zones_t::debug_ucmap(
    const DCArrayKokkos<double>& pos,                 // coords (positions at t)
    const DCArrayKokkos<double>& vel,                 // vel   (velocities at t)
    double dt,                                        // same dt passed to ucmap()
    const CArrayKokkos<double>&  cohesive_zone_orientation,          //
    const CArrayKokkos<size_t>&  overlapping_node_gids,
    const CArrayKokkos<double>&  local_opening              // what ucmap() already wrote
){
    printf("\n==================== ucmap debug ====================\n");

    for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
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

        // prints
        printf("\n-- Pair %zu (A gid=%zu, B gid=%zu) --\n", i, A, B);
        printf("  current_norm (at t)=(%.6g, %.6g, %.6g), ||current_norm||=%.6g\n",
               current_norm_x, current_norm_y, current_norm_z, sqrt(current_norm_x*current_norm_x + current_norm_y*current_norm_y + current_norm_z*current_norm_z));
        printf("  Nodal Displacement at t: u_t=(%.6g, %.6g, %.6g)\n", u_rel_x_t, u_rel_y_t, u_rel_z_t);
        printf("  Nodal Velocities at t: v_t=(%.6g, %.6g, %.6g)\n", v_rel_x_t, v_rel_y_t, v_rel_z_t);
        printf("  dt=%.6g\n", dt);
        //printf("  Forward Euler Method:\n");
        printf("  Predicted Nodal Displacement at t+dt: u_tdt=(%.6g, %.6g, %.6g)\n",
               u_x_tdt, u_y_tdt, u_z_tdt);
        printf("    Normal Crack Opening Magnitude at t: u_norm_mag_t=%.9g\n", u_norm_mag_t);
        printf("    Tangential Crack Opening Magnitude at t: u_tan_mag_t=%.9g\n", u_tan_mag_t);
        printf("    -> Forward Euler Predicted Normal Crack Opening at t+dt: u_norm_mag_tdt=%.9g\n", u_norm_mag_tdt);
        printf("    -> Forward Euler Predicted Tangential Crack Opening Magnitude at t+dt: u_tan_mag_tdt=%.9g\n", u_tan_mag_tdt);
        printf("  stored local_opening: [u_norm_mag_t=%.9g, u_tan_mag_t=%.9g, u_norm_mag_tdt=%.9g, u_tan_mag_tdt=%.9g]\n",
               local_opening(i,0), local_opening(i,1), local_opening(i,2), local_opening(i,3));

        printf("  diff vs stored: "
               "d_un_t=%.3e d_utan_t=%.3e d_un_tdt=%.3e d_utan_tdt=%.3e\n",
               u_norm_mag_t - local_opening(i,0),
               u_tan_mag_t - local_opening(i,1),
               u_norm_mag_tdt - local_opening(i,2),
               u_tan_mag_tdt - local_opening(i,3));
    }

    printf("\n==============================================================\n");
}

// ======================== END ucmap debug ========================    

// ======================== cohesive_zone_var_update debug ========================
// ======================== cohesive_zone_var_update DEBUG ========================

KOKKOS_FUNCTION
void cohesive_zones_t::debug_cohesive_zone_var_update(
    const CArrayKokkos<double>& local_opening,
    const double dt,
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

    printf("================== cohesive_zone_var_update debug ==================\n");
    printf("cohesive zone input parameters (bdy_set=%d) --\n", bdy_set);
    printf("  E_inf=%.9g  a1=%.9g  n_exp=%.9g  u_n*=%.9g  u_t*=%.9g  num_prony_terms=%d  dt=%.9g\n",
           E_inf, a1, n_exp, u_n_star, u_t_star, num_prony_terms, dt);
    printf("  E_dt=%.9g\n", E_dt);

    // loop over each cohesive zone node pair
    for (size_t i = 0; i < overlapping_node_gids.dims(0); i++){

        const size_t gidA = overlapping_node_gids(i,0);
        const size_t gidB = overlapping_node_gids(i,1);
        printf("\n-- Pair %zu (A gid=%zu, B gid=%zu) --\n", i, gidA, gidB);

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

        printf("  local_opening(t):     u_norm_mag_t=%.9g  u_tan_mag_t=%.9g\n", u_norm_mag_t,   u_tan_mag_t);
        printf("  local_opening(tdt):  u_norm_mag_tdt=%.9g  u_tan_mag_tdt=%.9g\n", u_norm_mag_tdt, u_tan_mag_tdt);
        printf("  lambda_t=%.9g  lambda_tdt=%.9g  lambda_dot_t=%.9g  (stored lambda_dot_t=%.9g)\n",
                lambda_t, lambda_tdt, lam_dot, delta_internal_vars(i,0));

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
        printf("  alpha(t)=%.9g  d_alpha(dt)=%.9g  alpha(t+dt)=%.9g\n",
                alpha_t, d_alpha, alpha_tdt);                

        // per Prony term details and sums
        for (int j = 0; j < num_prony_terms; ++j) {
            const int    base   = fractureStressBC::BCVars::prony_base + 2*j;
            const double Ej     = stress_bc_global_vars(bdy_set, base);
            const double tauj   = stress_bc_global_vars(bdy_set, base + 1);
            const double tau_eff= (tauj > 0.0) ? tauj : std::numeric_limits<double>::min();
            const double aexp   = exp(-dt / tau_eff);
            const double sigma_current = internal_vars(i, 4 + j);          
            const double sigma_next = delta_internal_vars(i, 4 + j);    
            printf("  Prony[%d]: E=%.9g  tau=%.9g  a=exp(-dt/tau)=%.9g  sigma_current=%.9g  sigma_next=%.9g\n",
                    j, Ej, tauj, aexp, sigma_current, sigma_next);
            sigma_sum   += sigma_next;
            sigma_sum_exp += (1.0 - exp(-dt / tau_eff)) * sigma_next;
        }
        printf("  sigma_sum=%.9g  sigma_sum_exp=%.9g\n", sigma_sum, sigma_sum_exp);

        // scalar pieces in traction increment formula
        const double deltaE_term  = E_dt * lam_dot * dt;
        const double elastic_term = E_inf * lambda_t + sigma_sum;   
        const double damp_term    = -sigma_sum_exp;

        // inverses (guarded)
        const double inv_uns_tdt = 1.0 / (u_n_star * ((lambda_tdt > 0.0) ? lambda_tdt : std::numeric_limits<double>::min()));
        const double inv_uns_t   = 1.0 / (u_n_star * ((lambda_t   > 0.0) ? lambda_t   : std::numeric_limits<double>::min()));
        const double inv_uts_tdt = 1.0 / (u_t_star * ((lambda_tdt > 0.0) ? lambda_tdt : std::numeric_limits<double>::min()));
        const double inv_uts_t   = 1.0 / (u_t_star * ((lambda_t   > 0.0) ? lambda_t   : std::numeric_limits<double>::min()));
        printf("  terms: deltaE_term=%.9g  elastic_term=%.9g  damp_term=%.9g\n", deltaE_term, elastic_term, damp_term);

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
        printf("  traction_norm terms: [%.9g, %.9g, %.9g, %.9g]\n", termN1, termN2, termN3, termN4);
        printf("  traction_tan terms: [%.9g, %.9g, %.9g, %.9g]\n", termT1, termT2, termT3, termT4);
        printf("  Tn(t) (normal traction at time t) = %.9g\n", fabs(termN3));
        //printf(" normal traction at time t: Tn_t=%.9g\n", termN3);
        //printf(" normal traction at time t+dt: Tn_tdt=%.9g\n", termN2);
        

    
        //printf("  STORED: lambda_dot_t=%.9g, delta_a=%.9g, traction_norm(normal traction increment)=%.9g, traction_tan(tangential traction increment)=%.9g\n",
                //delta_internal_vars(i,0), delta_internal_vars(i,1),
                //delta_internal_vars(i,2), delta_internal_vars(i,3));

        printf("  STORED: lambda_dot_t=%.9g\n", delta_internal_vars(i,0));
        printf("  STORED: delta_a=%.9g\n", delta_internal_vars(i,1));
        printf("  STORED: traction_norm(normal traction increment)=%.9g\n", delta_internal_vars(i,2));
        printf("  STORED: traction_tan(tangential traction increment)=%.9g\n", delta_internal_vars(i,3));

        //if (num_prony_terms > 0) {
//
         //   printf("  STORED : ");
         //   for (int j = 0; j < num_prony_terms; ++j) {
        //        printf("%s%.9g", (j==0?"[":", "), delta_internal_vars(i, 4 + j));
        //    }
        //    printf("]\n");

       // }

    }
    printf("\n================ end cohesive_zone_var_update debug ================\n");
}    
// ======================== END cohesive_zone_var_update debug ======================== 

// **************************************************************** Fierro Conversion **************************************************************** 
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
// **************************************************************** Fierro Conversion **************************************************************** 

// **************************************************************** Fierro Conversion **************************************************************** 
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

// **************************************************************** Fierro Conversion **************************************************************** 

// **************************************************************** Fierro Conversion **************************************************************** 
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
// // ORIGINAL MATCHED FACE ALGORITHM (COMMENT TO TRY NEW ONE)
//     // find FIRST opposing/coincident face match, write at slots
//     for (size_t i = 0; i < overlapping_node_gids.dims(0); ++i) {
//         bool found = false;

//         // small stack buffers + Views
//         double nA[3], rA[3], sA[3], cA[3];
//         double nB[3], rB[3], sB[3], cB[3];
//         ViewCArrayKokkos<double> nAj(&nA[0],3), rAj(&rA[0],3), sAj(&sA[0],3), cAj(&cA[0],3);
//         ViewCArrayKokkos<double> nBk(&nB[0],3), rBk(&rB[0],3), sBk(&sB[0],3), cBk(&cB[0],3);

//         // A candidates span indices [0 .. 3*max-1]
//         for (int j = 0; j < static_cast<int>(3 * max_elem_in_cohesive_zone) && !found; ++j) {
//             const int fA = cohesive_zone_faces(i, j);
//             const size_t slotA = static_cast<size_t>(j / 3);
//             const int eA = cohesive_zone_info(i, slotA);
//             if (fA < 0 || eA < 0) continue;

//             compute_face_geometry(
//                 state.node.coords, mesh,
//                 state.node.coords, mesh.nodes_in_elem,
//                 static_cast<size_t>(fA), static_cast<size_t>(eA),
//                 nAj, rAj, sAj, cAj
//             );

//             // B candidates span indices [0 .. 3*max-1] but stored with offset +3*max
//             for (int k = 0; k < static_cast<int>(3 * max_elem_in_cohesive_zone) && !found; ++k) {
//                 const int fB = cohesive_zone_faces(i, k + 3*max_elem_in_cohesive_zone);
//                 const size_t slotB = static_cast<size_t>(k / 3);
//                 const int eB = cohesive_zone_info(i, slotB + max_elem_in_cohesive_zone);
//                 if (fB < 0 || eB < 0) continue;

//                 compute_face_geometry(
//                     state.node.coords, mesh,
//                     state.node.coords, mesh.nodes_in_elem,
//                     static_cast<size_t>(fB), static_cast<size_t>(eB),
//                     nBk, rBk, sBk, cBk
//                 );

//                 // ABS centroid distance + opposite normals
//                 const double dx = cAj(0) - cBk(0);
//                 const double dy = cAj(1) - cBk(1);
//                 const double dz = cAj(2) - cBk(2);
//                 const double dist = sqrt(dx*dx + dy*dy + dz*dz);
//                 const double dot  = nAj(0)*nBk(0) + nAj(1)*nBk(1) + nAj(2)*nBk(2);

//                 if (dist <= tol && dot <= -1.0 + 1.0e-8) {
//                     // write at the slots that produced the match
//                     cohesive_zone_info(i, 2*max_elem_in_cohesive_zone + slotA) = fA; // A-face at A-slot
//                     cohesive_zone_info(i, 3*max_elem_in_cohesive_zone + slotB) = fB; // B-face at B-slot
//                     found = true;
//                 }
//             }
//         }
//     }
// // ORIGINAL MATCHED FACE ALGORITHM (COMMENT TO TRY NEW ONE)

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

// **************************************************************** Fierro Conversion **************************************************************** 

// **************************************************************** Fierro Conversion **************************************************************** 
// averaging all the normals of each individual cohesive zone
// averages based on the faces that are contributing to the cohesive zone

// variable map:
// Gavin --> Fierro:
// nodes --> State.node.coords
// nodes + ut --> pos (reference + total displacemetn up to time t)
// nodes + ut + us --> pos (reference + total displacement up to time t + this step)
// conn --> mesh.nodes_in_elem
// nvcz --> overlapping_node_gids.dims(0)
// vczconn --> overlapping_node_gids
// vczinfo --> cohesive_zone_info
// interpvals / MasterShapes() --> compute_face_geometry()
// el0, el1 --> eA, eB (incident elements from cohesive_zone_info block 0 and block 1)
// lp0, lp1 --> fA, fB (matched faces from cohesive_zone_info block 2 and block 3)
// ln0, ln1 --> kA, kB (local corner indices from cohesive_zone_info block 4 and block 5)
// currel0t, currel0tdt --> compute_face_geometry() 

// compute normals on the fly with compute_face_geometry()
// orient needs cohesive_zone_info A side elems, B side elems, matched A faces, matched B faces, local kA, local kB (corner indices)
// for each row i (cohesive zone pair) = (overlapping_node_gids.dims(0)), loop j and when both matched faces are >= 0, call compute_face_geometry()...........
//... for the faces on the current config, sum normals, normalize, store. (same as Gavin's average and normalize but using compute_face_geometry()

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

// **************************************************************** Fierro Conversion **************************************************************** 

// **************************************************************** FROM GAVIN'S CODE **************************************************************** 
// // calculates the opening of the crack in global coordinates and then maps it to local coordinates using normal vectors in
// // vczorient array
// // inputs: nodes, ut, us, vczorient, vczconn, nvcz
// // output: ulocvcz (4 columns: unt untdt utt utdt)
// void ucmap(CArray<double>& nodes, CArray<double>& ut, CArray<double>& us, CArray<double>& vczorient, CArray<int>& vczconn, int nvcz, CArray<double>& ulocvcz) {
//     // intermediate variables
//     auto uglobt = CArray <double> (3);
//     auto uglobtdt = CArray <double> (3);
//     double umagt = 0;
//     double umagtdt = 0;
    
//     // looping through each vcz
//     for (int i = 0; i < nvcz; i++) {
//         // calculating displacement between vcz nodes in global frame
//         uglobt(0) = (nodes(vczconn(i,1),0) + ut(3*vczconn(i,1))) - (nodes(vczconn(i,0),0) + ut(3*vczconn(i,0)));
//         uglobt(1) = (nodes(vczconn(i,1),1) + ut(3*vczconn(i,1)+1)) - (nodes(vczconn(i,0),1) + ut(3*vczconn(i,0)+1));
//         uglobt(2) = (nodes(vczconn(i,1),2) + ut(3*vczconn(i,1)+2)) - (nodes(vczconn(i,0),2) + ut(3*vczconn(i,0)+2));
//         uglobtdt(0) = (nodes(vczconn(i,1),0) + ut(3*vczconn(i,1)) + us(3*vczconn(i,1))) - (nodes(vczconn(i,0),0) + ut(3*vczconn(i,0)) + us(3*vczconn(i,0)));
//         uglobtdt(1) = (nodes(vczconn(i,1),1) + ut(3*vczconn(i,1)+1) + us(3*vczconn(i,1)+1)) - (nodes(vczconn(i,0),1) + ut(3*vczconn(i,0)+1) + us(3*vczconn(i,0)+1));
//         uglobtdt(2) = (nodes(vczconn(i,1),2) + ut(3*vczconn(i,1)+2) + us(3*vczconn(i,1)+2)) - (nodes(vczconn(i,0),2) + ut(3*vczconn(i,0)+2) + us(3*vczconn(i,0)+2));
        
//         // dotting with normal vector to get normal component of displacement
//         ulocvcz(i,0) = uglobt(0)*vczorient(i,0) + uglobt(1)*vczorient(i,1) + uglobt(2)*vczorient(i,2);
//         ulocvcz(i,2) = uglobtdt(0)*vczorient(i,3) + uglobtdt(1)*vczorient(i,4) + uglobtdt(2)*vczorient(i,5);

//         // calcuting magnitude of global u vectors
//         umagt = sqrt(uglobt(0)*uglobt(0) + uglobt(1)*uglobt(1) + uglobt(2)*uglobt(2));
//         umagtdt = sqrt(uglobtdt(0)*uglobtdt(0) + uglobtdt(1)*uglobtdt(1) + uglobtdt(2)*uglobtdt(2));
        
//         // calculating tangent component of displacement assuming that us* == ur*
//         ulocvcz(i,1) = sqrt(abs(umagt*umagt - ulocvcz(i,0)*ulocvcz(i,0)));
//         ulocvcz(i,3) = sqrt(abs(umagtdt*umagtdt - ulocvcz(i,2)*ulocvcz(i,2)));
//     }
    
// }
// **************************************************************** FROM GAVIN'S CODE **************************************************************** 

// **************************************************************** Fierro Conversion **************************************************************** 
// variable map:
// Gavin --> Fierro:
// nodes --> State.node.coords
// nodes + ut --> pos (reference + total displacemetn up to time t)
// nodes + ut + us --> pos (reference + total displacement up to time t + this step)
// vczorient --> cohesive_zone_orientation
// vczconn --> overlapping_node_gids
// nvcz --> overlapping_node_gids.dims(0)

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

// **************************************************************** Fierro Conversion **************************************************************** 

// **************************************************************** FROM GAVIN'S CODE **************************************************************** 
// // This function calculates delta_delta_internal variables (ddinvars)
// // Inputs: cohezive zone local displacments (uvcz [un(t), ut(t), un(t)+dun, ut(t)+dut]x3), time step size (delt), cohesive zone properties (Einf, a1, n, Eandrhom, uns, uts),
// // NPT (number of prony terms), nvcz (number of cohesive zones), internal variables [lambdadot(t-delt) , alpha , <tn,tt> , sig_1 through sig_m] (invars),
// // Output: delta_internal variables [lambdadot(t), d_alpha, d_<tn,tt>, d_sig_1 through d_sig_m]xNCZ (dinvars)
// void internalvars(CArray<double> uvcz, double delt, double Einf, double a1, double n, CArray<double> Eandrhom, double uns, double uts, int npt, int nvcz, CArray<double> invars, CArray<double> dinvars) {
    
//     // Initializing variables for calculations
//     double lambdadt = 0;
//     double lambdat = 0;
//     double dadt = 0;
//     double Edelt = 0;
//     double sumsig = 0;
//     double sumsigexp = 0;

//     // Calculating E(delt) NOTE: This is a redundant calculation so long as delt is uniform across time steps that ideally should be moved elsewhere to avoid repeating it unless delt is made to vary across steps
//     Edelt += Einf;
//     for (int i = 0; i < npt; i++) {
//         Edelt += Eandrhom(i,0) * Eandrhom(i,1) * (1 - exp(-delt / Eandrhom(i,1))) / delt;
//     }

//     // Looping over each unique pair
//     for (int i = 0; i < nvcz; i++) {

//         // Calculating the lambda(t)+dlambda and lambda(t) values
//         lambdadt = sqrt((uvcz(i,2) / uns) * (uvcz(i,2) / uns) + (uvcz(i,3) / uts) * (uvcz(i,3) / uts));
//         lambdat = sqrt((uvcz(i,0) / uns) * (uvcz(i,0) / uns) + (uvcz(i,1) / uts) * (uvcz(i,1) / uts));

//         // Initializing and calculating lambda rate for current load step (lambdadot(t))
//         dinvars(i,0) = (lambdadt - lambdat) / delt;

//         // Calculating dalpha/dt for current load step
//         if (dinvars(i,0) > 0) {
//             dadt = a1 * pow(((lambdadt + lambdat) / 2), n);
//         }
//         else {
//             dadt = 0;
//         }

//         // Preventing divide by zero errors in final calculation: the value is set equal to the minimum possible positive value that can be represented by a double in C++
//         if (lambdat == 0) {
//             lambdat = std::numeric_limits<double>::min();
//         }
//         if (lambdadt == 0) {
//             lambdadt = std::numeric_limits<double>::min();
//         }
        
//         // Updating dalpha
//         dinvars(i,1) = dadt * delt;

//         // Updating dsigma values for prony terms
//         for (int j = 0; j < npt; j++) {
//             dinvars(i,4 + j) = exp(-delt / Eandrhom(j,1)) * invars(i,4 + j) + Eandrhom(j,0) * Eandrhom(j,1) * invars(i,0) * (1 - exp(-delt / Eandrhom(j,1)));
//         }

//         // Resetting sumsig and sumsigexp on each iteration of the loop prior to calculating them
//         sumsig = 0;
//         sumsigexp = 0;

//         // Calculating sigma sums and sigma product sums (the deltaE term in the residual traction) in the residual traction terms
//         for (int j = 0; j < npt; j++) {
//             sumsig += dinvars(i,4 + j);
//             sumsigexp += (1 - exp(-delt / Eandrhom(j,1))) * dinvars(i,4 + j);
//         }

//         // Enforcing alpha domain limitations
//         if ((invars(i,1) + dinvars(i,1)) > 1) {
//             dinvars(i,1) = 1 - invars(i,1);
//         }

//         // Updating dtraction values for each cohesive zone along the edge
//         dinvars(i,2) = uvcz(i,2) / (uns * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * Edelt * dinvars(i,0) * delt
//             + uvcz(i,2) / (uns * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * (Einf * lambdat + sumsig)
//             - uvcz(i,0) / (uns * lambdat) * (1 - (invars(i,1))) * (Einf * lambdat + sumsig)
//             + uvcz(i,2) / (uns * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * -sumsigexp;
//         dinvars(i,3) = uvcz(i,3) / (uts * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * Edelt * dinvars(i,0) * delt
//             + uvcz(i,3) / (uts * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * (Einf * lambdat + sumsig)
//             - uvcz(i,1) / (uts * lambdat) * (1 - (invars(i,1))) * (Einf * lambdat + sumsig)
//             + uvcz(i,3) / (uts * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * -sumsigexp;
//     }
// }
// **************************************************************** FROM GAVIN'S CODE **************************************************************** 

// **************************************************************** Fierro Conversion **************************************************************** 
KOKKOS_FUNCTION
void cohesive_zones_t::cohesive_zone_var_update(
    const CArrayKokkos<double>& local_opening,
    const double dt,
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


        const double inv_uns_lambda_tdt = 1.0 / (u_n_star * ((lambda_tdt>0.0)?lambda_tdt:std::numeric_limits<double>::min())); // avoid div by zero
        const double inv_uns_lambda_t   = 1.0 / (u_n_star * ((lambda_t >0.0)?lambda_t :std::numeric_limits<double>::min()));
        const double inv_uts_lambda_tdt = 1.0 / (u_t_star * ((lambda_tdt>0.0)?lambda_tdt:std::numeric_limits<double>::min()));
        const double inv_uts_lambda_t   = 1.0 / (u_t_star * ((lambda_t >0.0)?lambda_t :std::numeric_limits<double>::min()));



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
// **************************************************************** Fierro Conversion **************************************************************** 

// **************************************************************** FROM GAVIN'S CODE **************************************************************** 
// // This function calculates the loads due to the cohesive zone tractions in the global frame of reference
// // inputs: vczorient, vczconn, nvcz, invars, dinvars, vczinfo, nodes, conn, ut, us
// // output: Fvcz (global load vector from cohesive zones), KVCZ (global stiffness contributions from cohesive zones)
// void vczloads(CArray<double> vczorient, CArray<int> vczconn, double nvcz, CArray<double> invars, CArray<double> dinvars, CArray<int> vczinfo, CArray<double> nodes, CArray<int> conn, CArray<double> ut, CArray<double> us, CArray<double> Eandrhom, int npt, double delt, double Einf, double uns, double uts, CArray<double> Fvcz,CArray<double> Kvcz) {
//     // reset Fvcz and Kvcz
//     Fvcz.set_values(0);
//     Kvcz.set_values(0);
    
//     // intermediate variables
//     auto curr = CArray <double> (8,3); // current configuration of an element
//     double gp = 0.5773502691896257;
//     auto gps = CArray <double> (4,3);
//     auto J = CArray <double> (3,3);
//     double area;
//     auto interpvals = FArray <double> (8,4);
//     auto Fn = CArray <double> (3);
//     auto Ft = CArray <double> (3);
//     double tanmag;
//     double udotn;
//     auto uglobtdt = CArray <double> (3);
//     auto tanvec = CArray <double> (3);
//     double kn;
//     double kt;
//     double alpha;
//     double beta;
//     auto kvcz = CArray <double> (6,6);

//     // Edelt calculator for stiffness of truss element values
//     double Edelt = 0;
//     Edelt += Einf;
//     for (int i = 0; i < npt; i++) {
//         Edelt += Eandrhom(i,0) * Eandrhom(i,1) * (1 - exp(-delt / Eandrhom(i,1))) / delt;
//     }

//     // looping over each cohesive zone
//     for (int i = 0; i < nvcz; i++) {
//         // looping over each element for the cohesive zones
//         for (int j = 0; j < vczinfo.dims(1)/6; j++) {
//             // checking that the element index in vczinfo is populated with an element number
//             if (vczinfo(i,j) != -1) {

//                 // calculate current configuration of the element in the index
//                 for (int k = 0; k < 8; k++) {
//                     for (int m = 0; m < 3; m++) {
//                         curr(k,m) = nodes(conn(vczinfo(i,j),k),m) + ut(3 * conn(vczinfo(i,j),k) + m) + us(3 * conn(vczinfo(i,j),k) + m);
//                     }
//                 }
                
//                 // finding gauss points to loop through based on patch number
//                 switch(vczinfo(i,2*vczinfo.dims(1)/6+j)){
//                 case 0:
//                     gps(0,0) = -1;
//                     gps(1,0) = -1;
//                     gps(2,0) = -1;
//                     gps(3,0) = -1;
//                     gps(0,1) = -gp;
//                     gps(1,1) = gp;
//                     gps(2,1) = -gp;
//                     gps(3,1) = gp;
//                     gps(0,2) = -gp;
//                     gps(1,2) = -gp;
//                     gps(2,2) = gp;
//                     gps(3,2) = gp;
//                     break;
//                 case 1:
//                     gps(0,0) = 1;
//                     gps(1,0) = 1;
//                     gps(2,0) = 1;
//                     gps(3,0) = 1;
//                     gps(0,1) = -gp;
//                     gps(1,1) = gp;
//                     gps(2,1) = -gp;
//                     gps(3,1) = gp;
//                     gps(0,2) = -gp;
//                     gps(1,2) = -gp;
//                     gps(2,2) = gp;
//                     gps(3,2) = gp;
//                     break;
//                 case 2:
//                     gps(0,1) = -1;
//                     gps(1,1) = -1;
//                     gps(2,1) = -1;
//                     gps(3,1) = -1;
//                     gps(0,0) = -gp;
//                     gps(1,0) = gp;
//                     gps(2,0) = -gp;
//                     gps(3,0) = gp;
//                     gps(0,2) = -gp;
//                     gps(1,2) = -gp;
//                     gps(2,2) = gp;
//                     gps(3,2) = gp;
//                     break;
//                 case 3:
//                     gps(0,1) = 1;
//                     gps(1,1) = 1;
//                     gps(2,1) = 1;
//                     gps(3,1) = 1;
//                     gps(0,0) = -gp;
//                     gps(1,0) = gp;
//                     gps(2,0) = -gp;
//                     gps(3,0) = gp;
//                     gps(0,2) = -gp;
//                     gps(1,2) = -gp;
//                     gps(2,2) = gp;
//                     gps(3,2) = gp;
//                     break;
//                 case 4:
//                     gps(0,2) = -1;
//                     gps(1,2) = -1;
//                     gps(2,2) = -1;
//                     gps(3,2) = -1;
//                     gps(0,0) = -gp;
//                     gps(1,0) = gp;
//                     gps(2,0) = -gp;
//                     gps(3,0) = gp;
//                     gps(0,1) = -gp;
//                     gps(1,1) = -gp;
//                     gps(2,1) = gp;
//                     gps(3,1) = gp;
//                     break;
//                 case 5:
//                     gps(0,2) = 1;
//                     gps(1,2) = 1;
//                     gps(2,2) = 1;
//                     gps(3,2) = 1;
//                     gps(0,0) = -gp;
//                     gps(1,0) = gp;
//                     gps(2,0) = -gp;
//                     gps(3,0) = gp;
//                     gps(0,1) = -gp;
//                     gps(1,1) = -gp;
//                     gps(2,1) = gp;
//                     gps(3,1) = gp;
//                     break;
//                 }
                
//                 // jacobian calculations for the patch area
//                 area = 0;
//                 for (int k = 0; k < 4; k++) {

//                     // calculating jacobian
//                     MasterShapes(interpvals, gps(k,0),gps(k,1),gps(k,2));
//                     J.set_values(0);
//                     for (int m = 0; m < 3; m++) {
//                         for (int o = 0; o < 3; o++) {
//                             for (int p = 0; p < 8; p++) {
//                                 J(m,o) +=  curr(p,o)*interpvals(p,m+1);
//                             }
//                         }
//                     }
                    
//                     // cross product of jacobian on surface to find area
//                     switch(2*vczinfo.dims(1)/6+j) {
//                     case 0:
//                         area += sqrt(pow((J(1,1)*J(2,2) - J(2,1)*J(1,2)),2) + pow((J(1,0)*J(2,2) - J(2,0)*J(1,2)),2) + pow((J(1,0)*J(2,1) - J(2,0)*J(1,1)),2));
//                         break;
//                     case 1:
//                         area += sqrt(pow((J(1,1)*J(2,2) - J(2,1)*J(1,2)),2) + pow((J(1,0)*J(2,2) - J(2,0)*J(1,2)),2) + pow((J(1,0)*J(2,1) - J(2,0)*J(1,1)),2));
//                         break;
//                     case 2:
//                         area += sqrt(pow((J(0,1)*J(2,2) - J(2,1)*J(0,2)),2) + pow((J(0,0)*J(2,2) - J(2,0)*J(0,2)),2) + pow((J(0,0)*J(2,1) - J(2,0)*J(0,1)),2));
//                         break;
//                     case 3:
//                         area += sqrt(pow((J(0,1)*J(2,2) - J(2,1)*J(0,2)),2) + pow((J(0,0)*J(2,2) - J(2,0)*J(0,2)),2) + pow((J(0,0)*J(2,1) - J(2,0)*J(0,1)),2));
//                         break;
//                     case 4:
//                         area += sqrt(pow((J(0,1)*J(1,2) - J(1,1)*J(0,2)),2) + pow((J(0,0)*J(1,2) - J(1,0)*J(0,2)),2) + pow((J(0,0)*J(1,1) - J(1,0)*J(0,1)),2));
//                         break;
//                     case 5:
//                         area += sqrt(pow((J(0,1)*J(1,2) - J(1,1)*J(0,2)),2) + pow((J(0,0)*J(1,2) - J(1,0)*J(0,2)),2) + pow((J(0,0)*J(1,1) - J(1,0)*J(0,1)),2));
//                         break;
//                     }
//                 }
                
//                 // calculating normal force vector
//                 Fn(0) = (invars(i,2)+dinvars(i,2))*(area/4)*vczorient(i,3);
//                 Fn(1) = (invars(i,2)+dinvars(i,2))*(area/4)*vczorient(i,4);
//                 Fn(2) = (invars(i,2)+dinvars(i,2))*(area/4)*vczorient(i,5);
                
//                 // calculating tangent vector     u_t = u - (u dot n)n
//                 uglobtdt(0) = (nodes(vczconn(i,1),0) + ut(3*vczconn(i,1)) + us(3*vczconn(i,1))) - (nodes(vczconn(i,0),0) + ut(3*vczconn(i,0)) + us(3*vczconn(i,0)));
//                 uglobtdt(1) = (nodes(vczconn(i,1),1) + ut(3*vczconn(i,1)+1) + us(3*vczconn(i,1)+1)) - (nodes(vczconn(i,0),1) + ut(3*vczconn(i,0)+1) + us(3*vczconn(i,0)+1));
//                 uglobtdt(2) = (nodes(vczconn(i,1),2) + ut(3*vczconn(i,1)+2) + us(3*vczconn(i,1)+2)) - (nodes(vczconn(i,0),2) + ut(3*vczconn(i,0)+2) + us(3*vczconn(i,0)+2));
//                 udotn = uglobtdt(0)*vczorient(i,3) + uglobtdt(1)*vczorient(i,4) + uglobtdt(2)*vczorient(i,5);
//                 tanvec(0) = uglobtdt(0) - udotn*vczorient(i,3);
//                 tanvec(1) = uglobtdt(1) - udotn*vczorient(i,4);
//                 tanvec(2) = uglobtdt(2) - udotn*vczorient(i,5);
//                 tanmag = sqrt(tanvec(0)*tanvec(0) + tanvec(1)*tanvec(1) + tanvec(2)*tanvec(2));
//                 if (tanmag != 0) {
//                     tanvec(0) /= tanmag;
//                     tanvec(1) /= tanmag;
//                     tanvec(2) /= tanmag;
//                 }

//                 // calculating tangent force vector
//                 Ft(0) = (invars(i,3)+dinvars(i,3))*(area/4)*tanvec(0);
//                 Ft(1) = (invars(i,3)+dinvars(i,3))*(area/4)*tanvec(1);
//                 Ft(2) = (invars(i,3)+dinvars(i,3))*(area/4)*tanvec(2);
//                 //prarr(Ft);
//                 // adding to global force vector
//                 Fvcz(3*vczconn(i,0)) += Fn(0) + Ft(0);
//                 Fvcz(3*vczconn(i,0)+1) += Fn(1) + Ft(1);
//                 Fvcz(3*vczconn(i,0)+2) += Fn(2) + Ft(2);
//                 Fvcz(3*vczconn(i,1)) -= Fn(0) + Ft(0);
//                 Fvcz(3*vczconn(i,1)+1) -= Fn(1) + Ft(1);
//                 Fvcz(3*vczconn(i,1)+2) -= Fn(2) + Ft(2);


// DO NOT WORRY ABOUT THIS THIS IS FOR QUASI STIFFNESS CONTRIBUTION
//                 // calculating rotation angles of truss element
//                 // if u_x is zero it causes divide by zero errors
//                 if (uglobtdt(0) == 0) {
                    
//                     if (uglobtdt(1) == 0) {
//                         // u_x and u_y both zero then no rotation wrt z axis
//                         alpha = 0;
//                     } else if (uglobtdt(1) > 0) {
//                         // u_x is zero and u_y is positive, truss orientation wrt z axis is 90 degrees (straight up)
//                         alpha = M_PI/2;
//                     } else {
//                         // u_x is zero and u_y is negative, truss orientation wrt z axis is 270 degrees (straight down)
//                         alpha = 3*M_PI/2;
//                     }

//                     if (uglobtdt(2) == 0) {
//                         // u_x and u_z both zero then no rotation wrt y axis
//                         beta = 0;
//                     } else if (uglobtdt(2) > 0) {
//                         // u_x is zero and u_y is positive, truss orientation wrt y axis is 90 degrees (straight up)
//                         beta = 3*M_PI/2;
//                     } else {
//                         // u_x is zero and u_y is negative, truss orientation wrt y axis is 270 degrees (straight down)
//                         beta = M_PI/2;
//                     }

//                 } else if (uglobtdt(0) > 0) {
//                     // if u_x is positive then the orientation wrt z and y axes is in either quadrant 1 or 4 of unit circle
//                     alpha = atan(uglobtdt(1)/uglobtdt(0));
//                     beta = atan(-uglobtdt(2)/uglobtdt(0));
//                 } else {
//                     // if u_x is negative then the orientation wrt z and y axes is in either quadrant 2 or 3 of unit circle
//                     // this case requires a correction to account for full 360 rotation due to range limits of arctan
//                     alpha = atan(uglobtdt(1)/uglobtdt(0))+M_PI;
//                     beta = atan(-uglobtdt(2)/uglobtdt(0))+M_PI;
//                 }
//                 /* if (uglobtdt(0) != 0) {
//                     alpha = atan(uglobtdt(1)/uglobtdt(0));
//                     beta = atan(-uglobtdt(2)/uglobtdt(0));
//                 } else {
//                     alpha = 0;
//                     beta = 0;
//                 } */
//                 //printf("a = %f    b = %f\n\n",alpha*180/3.14159,beta*180/3.14159);

//                 // calculating directional stiffness in local frame
//                 // note: only one tangent value needed bc of assumption urs == uss
//                 kn = (1 - (invars(i,1) + dinvars(i,1))) * Edelt * (area/4) / uns;
//                 kt = (1 - (invars(i,1) + dinvars(i,1))) * Edelt * (area/4) / uts;

//                 // populating local stiffness matrix based on +ccw turns
//                 // first rotation about x2 axis (beta), second rotation about x3 axis (alpha)
//                 kvcz(0,0) = cos(alpha)*cos(alpha)*cos(beta)*cos(beta)*kn + (cos(alpha)*cos(alpha)*sin(beta)*sin(beta) + sin(alpha)*sin(alpha))*kt;
//                 kvcz(0,1) = cos(alpha)*cos(beta)*cos(beta)*sin(alpha)*kn + (cos(alpha)*sin(alpha)*sin(beta)*sin(beta) - cos(alpha)*sin(alpha))*kt;
//                 kvcz(0,2) = (-cos(alpha)*cos(beta)*sin(beta))*kn + cos(alpha)*cos(beta)*sin(beta)*kt;
//                 kvcz(0,3) = (-cos(alpha)*cos(alpha)*cos(beta)*cos(beta))*kn + (- cos(alpha)*cos(alpha)*sin(beta)*sin(beta) - sin(alpha)*sin(alpha))*kt;
//                 kvcz(0,4) = (-cos(alpha)*cos(beta)*cos(beta)*sin(alpha))*kn + (cos(alpha)*sin(alpha) - cos(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(0,5) = cos(alpha)*cos(beta)*sin(beta)*kn + (-cos(alpha)*cos(beta)*sin(beta))*kt;
//                 kvcz(1,1) = cos(beta)*cos(beta)*sin(alpha)*sin(alpha)*kn + (cos(alpha)*cos(alpha) + sin(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(1,2) = (-cos(beta)*sin(alpha)*sin(beta))*kn + cos(beta)*sin(alpha)*sin(beta)*kt;
//                 kvcz(1,3) = (-cos(alpha)*cos(beta)*cos(beta)*sin(alpha))*kn + (cos(alpha)*sin(alpha) - cos(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(1,4) = (-cos(beta)*cos(beta)*sin(alpha)*sin(alpha))*kn + (- cos(alpha)*cos(alpha) - sin(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(1,5) = cos(beta)*sin(alpha)*sin(beta)*kn + (-cos(beta)*sin(alpha)*sin(beta))*kt;
//                 kvcz(2,2) = sin(beta)*sin(beta)*kn + cos(beta)*cos(beta)*kt;
//                 kvcz(2,3) = cos(alpha)*cos(beta)*sin(beta)*kn + (-cos(alpha)*cos(beta)*sin(beta))*kt;
//                 kvcz(2,4) = cos(beta)*sin(alpha)*sin(beta)*kn + (-cos(beta)*sin(alpha)*sin(beta))*kt;
//                 kvcz(2,5) = (-sin(beta)*sin(beta))*kn + (-cos(beta)*cos(beta))*kt;
//                 kvcz(3,3) = cos(alpha)*cos(alpha)*cos(beta)*cos(beta)*kn + (cos(alpha)*cos(alpha)*sin(beta)*sin(beta) + sin(alpha)*sin(alpha))*kt;
//                 kvcz(3,4) = cos(alpha)*cos(beta)*cos(beta)*sin(alpha)*kn + (cos(alpha)*sin(alpha)*sin(beta)*sin(beta) - cos(alpha)*sin(alpha))*kt;
//                 kvcz(3,5) = (-cos(alpha)*cos(beta)*sin(beta))*kn + cos(alpha)*cos(beta)*sin(beta)*kt;
//                 kvcz(4,4) = cos(beta)*cos(beta)*sin(alpha)*sin(alpha)*kn + (cos(alpha)*cos(alpha) + sin(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(4,5) = (-cos(beta)*sin(alpha)*sin(beta))*kn + cos(beta)*sin(alpha)*sin(beta)*kt;
//                 kvcz(5,5) = sin(beta)*sin(beta)*kn + cos(beta)*cos(beta)*kt;

//                 // symmetric fill
//                 kvcz(1,0) = kvcz(0,1);
//                 kvcz(2,0) = kvcz(0,2);
//                 kvcz(2,1) = kvcz(1,2);
//                 kvcz(3,0) = kvcz(0,3);
//                 kvcz(3,1) = kvcz(1,3);
//                 kvcz(3,2) = kvcz(2,3);
//                 kvcz(4,0) = kvcz(0,4);
//                 kvcz(4,1) = kvcz(1,4);
//                 kvcz(4,2) = kvcz(2,4);
//                 kvcz(4,3) = kvcz(3,4);
//                 kvcz(5,0) = kvcz(0,5);
//                 kvcz(5,1) = kvcz(1,5);
//                 kvcz(5,2) = kvcz(2,5);
//                 kvcz(5,3) = kvcz(3,5);
//                 kvcz(5,4) = kvcz(4,5);

//                 // populating Kvcz based on vcz connectivity
//                 for (int m = 0; m < 2; m++) {
//                     for (int o = 0; o < 2; o++) {
//                         for (int p = 0; p < 3; p++) {
//                             for (int q = 0; q < 3; q++) {
//                                 Kvcz(3*vczconn(i,m)+p,3*vczconn(i,o)+q) += kvcz(3*m + p, 3*o + q); 
//                             }
//                         }
//                     }
//                 }

//             }
//         }
//     }

// }
// **************************************************************** FROM GAVIN'S CODE **************************************************************** 

// **************************************************************** FROM GAVIN'S CODE **************************************************************** 
// // This function calls all necessary functions to update the vcz behavior based on a displacement field
// // inputs: nodes, conn, ne, ut, us, vczconn, vczinfo, nvcz, interpvals, vczorient, ulocvcz, npt, vcz material parameters (Einf, a1, n, Eandrhom, uns, uts), delt, invars, dinvars
// // outputs: Fvcz (global load vector due to vczs) and Kvcz (global stiffness matrix to characterize vczs)
// void vczupdate(CArray<double> nodes, CArray<int> conn, int ne, CArray<double> ut, CArray<double> us, CArray<int> vczconn, CArray<int> vczinfo, int nvcz, FArray<double> interpvals, CArray<double> vczorient, CArray<double> ulocvcz, int npt, double Einf, double a1, double n, CArray<double> Eandrhom, double uns, double uts, double delt, CArray<double> invars, CArray<double> dinvars, CArray<double> Fvcz, CArray<double> Kvcz) {
//     // find global orientation normal vectors
//     oriented(nodes,conn,ne,nvcz,ut,us,vczinfo,interpvals,vczorient);
//     // map displacement into local normal and tangent frame of reference;
//     ucmap(nodes,ut,us,vczorient,vczconn,nvcz,ulocvcz);
//     printf(" ");
//     // updating interal variables of cohesive zones
//     internalvars(ulocvcz,delt,Einf,a1,n,Eandrhom,uns,uts,npt,nvcz,invars,dinvars);
//     // populating global force vector and stiffness matrix
//     vczloads(vczorient,vczconn,nvcz,invars,dinvars,vczinfo,nodes,conn,ut,us,Eandrhom,npt,delt,Einf,uns,uts,Fvcz,Kvcz);
// }

// int main()
// {
//     // Initializing input variables from first 2 lines of input file then populating them with readFirstLines()
//     int NNODE;
//     int NEL;
//     int NDBC;
//     int NPL;
//     int NTL;
//     int BCFLAG;
//     int NLS;
//     double E;
//     double nu;
//     double t;
//     double tol;
//     int NUP;
//     double a1;
//     double n;
//     double Einf;
//     double delt;
//     int NPT;
//     double uns;
//     double uts;
//     std::string filename = "Input.txt";
//     readFirstLines(filename, NNODE, NEL, NDBC, NPL, NTL, BCFLAG, NLS, E, nu, t, tol, NUP, a1, n, Einf, delt, NPT, uns, uts);

//     // Initializing input variables based upon sizes from readFirstLines output then populating them with readTheRest()
//     auto NODES = CArray <double> (NNODE,3);
//     auto CONN = CArray <int> (NEL,8);
//     auto Eandrhom = CArray <double> (NPT,2);
//     auto UPs = CArray <int> (NUP,2);
//     int DBCcols;
//     int PLScols;
//     int TLScols;
//     if (BCFLAG == 0) {
//         DBCcols = 3;
//         PLScols = 3;
//         TLScols = 5;
//     }
//     else {
//         DBCcols = 2+NLS;
//         PLScols = 2+NLS;
//         TLScols = 2+(3*NLS);
//     }
//     auto DBCS = CArray <double> (NDBC, DBCcols);
//     auto PLS = CArray <double> (NPL, PLScols);
//     auto TLS = CArray <double> (NTL, TLScols);
//     readTheRest(filename, NUP, NPT, NNODE, NEL, NDBC, NPL, NTL, BCFLAG, NLS, NODES, CONN, DBCS, PLS, TLS, Eandrhom, UPs);
    
//     // Defining gauss point interpolation function values
//     // (gp__ where V means volume, p# means patch number, and the final number is that subset gauss point)
//     // 2 pt gauss quadrature (linear element) so V spans 0-7, p# spans 0-5, each patch spans 0-3, and weights are always 1
//     double gp = 0.5773502691896257;
//     auto gpV0 = FArray <double> (8,4);
//     MasterShapes(gpV0,-gp,-gp,-gp);
//     auto gpV1 = FArray <double> (8,4);
//     MasterShapes(gpV1,-gp,-gp,gp);
//     auto gpV2 = FArray <double> (8,4);
//     MasterShapes(gpV2,gp,-gp,gp);
//     auto gpV3 = FArray <double> (8,4);
//     MasterShapes(gpV3,gp,-gp,-gp);
//     auto gpV4 = FArray <double> (8,4);
//     MasterShapes(gpV4,-gp,gp,-gp);
//     auto gpV5 = FArray <double> (8,4);
//     MasterShapes(gpV5,-gp,gp,gp);
//     auto gpV6 = FArray <double> (8,4);
//     MasterShapes(gpV6,gp,gp,gp);
//     auto gpV7 = FArray <double> (8,4);
//     MasterShapes(gpV7,gp,gp,-gp);
    
//     // patch definition: [0, 1, 2, 3, 4, 5] = [xi1=-1, xi1=+1, xi2=-1, xi2=+1, xi3=-1, xi3=+1]
//     auto gpp00 = FArray <double> (8,4);
//     MasterShapes(gpp00,-1,-gp,-gp);
//     auto gpp01 = FArray <double> (8,4);
//     MasterShapes(gpp01,-1,-gp,gp);
//     auto gpp02 = FArray <double> (8,4);
//     MasterShapes(gpp02,-1,gp,-gp);
//     auto gpp03 = FArray <double> (8,4);
//     MasterShapes(gpp03,-1,gp,gp);

//     auto gpp10 = FArray <double> (8,4);
//     MasterShapes(gpp10,1,-gp,-gp);
//     auto gpp11 = FArray <double> (8,4);
//     MasterShapes(gpp11,1,-gp,gp);
//     auto gpp12 = FArray <double> (8,4);
//     MasterShapes(gpp12,1,gp,-gp);
//     auto gpp13 = FArray <double> (8,4);
//     MasterShapes(gpp13,1,gp,gp);

//     auto gpp20 = FArray <double> (8,4);
//     MasterShapes(gpp20,-gp,-1,-gp);
//     auto gpp21 = FArray <double> (8,4);
//     MasterShapes(gpp21,-gp,-1,gp);
//     auto gpp22 = FArray <double> (8,4);
//     MasterShapes(gpp22,gp,-1,-gp);
//     auto gpp23 = FArray <double> (8,4);
//     MasterShapes(gpp23,gp,-1,gp);

//     auto gpp30 = FArray <double> (8,4);
//     MasterShapes(gpp30,-gp,1,-gp);
//     auto gpp31 = FArray <double> (8,4);
//     MasterShapes(gpp31,-gp,1,gp);
//     auto gpp32 = FArray <double> (8,4);
//     MasterShapes(gpp32,gp,1,-gp);
//     auto gpp33 = FArray <double> (8,4);
//     MasterShapes(gpp33,gp,1,gp);

//     auto gpp40 = FArray <double> (8,4);
//     MasterShapes(gpp40,-gp,-gp,-1);
//     auto gpp41 = FArray <double> (8,4);
//     MasterShapes(gpp41,gp,-gp,-1);
//     auto gpp42 = FArray <double> (8,4);
//     MasterShapes(gpp42,-gp,gp,-1);
//     auto gpp43 = FArray <double> (8,4);
//     MasterShapes(gpp43,gp,gp,-1);

//     auto gpp50 = FArray <double> (8,4);
//     MasterShapes(gpp50,-gp,-gp,1);
//     auto gpp51 = FArray <double> (8,4);
//     MasterShapes(gpp51,gp,-gp,1);
//     auto gpp52 = FArray <double> (8,4);
//     MasterShapes(gpp52,-gp,gp,1);
//     auto gpp53 = FArray <double> (8,4);
//     MasterShapes(gpp53,gp,gp,1);

//     // calculating C matrix for linear elastic bulk elements
//     CArray <double> CMat = CMaterial(E,nu);
    
//     // Initializing necessary arrays and variables to be used during load step -> elements nested loops:
//     // elcoords 8x3, Kg 3NNx3NN, Fg 3NNx1, F02g 3NNx1, F01g 3NNx1, ut 3NNx1, us 3NNx1, dus 3NNx1, Kel 24x24, F01el 24x1,
//     // uel 8x3, S01 6x1, E01 3x3, gradu 3x3, detJ, dpsig 8x3
//     auto elcoords = CArray <double> (8,3);
//     auto Kg = CArray <double> (3*NNODE,3*NNODE);
//     auto Fg = CArray <double> (3*NNODE);
//     auto F02g = CArray <double> (3*NNODE);
//     auto F01g = CArray <double> (3*NNODE);
//     auto ut = CArray <double> (3*NNODE);
//     auto us = CArray <double> (3*NNODE);
//     auto dus = CArray <double> (3*NNODE);
//     auto Kel = CArray <double> (24,24);
//     auto F01el = CArray <double> (24);
//     auto uel = CArray <double> (8,3);
//     auto S01 = CArray <double> (6);
//     auto E01 = CArray <double> (3,3);
//     auto gradu = CArray <double> (3,3);
//     double detJ;
//     auto dpsig = CArray <double> (8,3);

//     // Initialize matrices that see +=
//     Kg.set_values(0);
//     Fg.set_values(0);
//     F02g.set_values(0);
//     F01g.set_values(0);
//     ut.set_values(0);
//     us.set_values(0);
//     dus.set_values(0);
//     Kel.set_values(0);
//     F01el.set_values(0);

//     // calling for mesh preprocessing to handle viscoelastic cohesive zones
//     auto vczconn = CArray <int> (NUP,2);
//     auto interpvals = FArray <double> (8,4);
//     CArray <int> vczinfo = VCZmeshpreprocess(NODES,CONN,NEL,NUP,UPs,interpvals);

//     // allocations for arrays needed to handle cohesive zone behavior
//     auto invars = CArray<double> (NUP,4+NPT);
//     auto dinvars = CArray<double> (NUP,4+NPT);
//     invars.set_values(0);
//     dinvars.set_values(0);
//     auto vczorient = CArray <double> (NUP,6);
//     auto ulocvcz = CArray <double> (NUP,4);
//     auto Fvcz = CArray<double> (3*NNODE);
//     auto Kvcz = CArray<double> (3*NNODE,3*NNODE);

//     // Defining boundary condition stepping arrays
//     int DBsize;
//     int PLsize;
//     int TLsize;
//     // BCFLAG: 0 for uniform stepping, 1 for custom stepping
//     // defining number of columns in step arrays
//     if (BCFLAG == 0) {
//         DBsize = 1;
//         PLsize = 1;
//         TLsize = 3;
//     }
//     else {
//         DBsize = NLS;
//         PLsize = NLS;
//         TLsize = 3*NLS;
//     }
//     // initializing step arrays
//     auto DBstep = CArray <double> (NDBC,DBsize);
//     auto PLstep = CArray <double> (NPL,PLsize);
//     auto TLstep = CArray <double> (NTL,TLsize);
//     // populating step arrays
//     if (BCFLAG == 0) {
//         for (int i = 0; i < NDBC; i++) {
//             DBstep(i,0) = DBCS(i,2) / NLS;
//         }

//         for (int i = 0; i < NPL; i++) {
//             PLstep(i,0) = PLS(i,2) / NLS;
//         }

//         for (int i = 0; i < NTL; i++) {
//             TLstep(i,0) = TLS(i,2) / NLS;
//             TLstep(i,1) = TLS(i,3) / NLS;
//             TLstep(i,2) = TLS(i,4) / NLS;
//         }
//     }
//     else {
//         for (int i = 0; i < NDBC; i++) {
//             for (int j = 0; j < NLS; j++) {
//                 if (j == 0) {
//                     // accounting for first load step assuming that the body is initially undeformed
//                     DBstep(i,j) = DBCS(i,2);
//                 }
//                 else {
//                     DBstep(i,j) = DBCS(i,j + 2) - DBCS(i,j + 1);
//                 }
//             }
//         }

//         for (int i = 0; i < NPL; i++) {
//             for (int j = 0; j < NLS; j++) {
//                 if (j == 0) {
//                     // accounting for first load step assuming that the body is initially undeformed
//                     PLstep(i,j) = PLS(i,2);
//                 }
//                 else {
//                     PLstep(i,j) = PLS(i,j + 2) - PLS(i,j + 1);
//                 }
//             }
//         }

//         for (int i = 0; i < NTL; i++) {
//             for (int j = 0; j < NLS; j++) {
//                 for (int k = 0; k < 3; k++) {
//                     if (j == 0) {
//                         // accounting for first load step assuming that the body is initially undeformed
//                         TLstep(i,k) = TLS(i,2 + k);
//                     }
//                     else {
//                         TLstep(i,3 * j + k) = TLS(i,3 * j + k + 2) - TLS(i,3 * (j - 1) + k + 2);
//                     }
//                 }
//             }
//         }
//     }

//     // Defining arrays for applying boundary conditions on the current load step
//     // PL02 carries current point loads, TL02 carries current traction loads, DB02 carries current prescribed displacement
//     // DBiter carries zeros for maintaining that the displacement step vector satifies DB02 during iteration
//     auto PL02 = CArray <double> (NPL,3);
//     auto TL02 = CArray <double> (NTL,5);
//     auto DB12 = CArray <double> (NDBC,3);
//     auto DBiter = CArray <double> (NDBC,3);
//     PL02.set_values(0);
//     TL02.set_values(0);
//     DB12.set_values(0);
//     DBiter.set_values(0);
//     // pulling node numbers and dofs for point loads and dirichlet, pulling element number and patch number for traction loads
//     for (int i = 0; i < NPL; i++) {
//         for (int j = 0; j < 2; j++) {
//             PL02(i,j) = PLS(i,j);
//         }
//     }
//     for (int i = 0; i < NTL; i++) {
//         for (int j = 0; j < 2; j++) {
//             TL02(i,j) = TLS(i,j);
//         }
//     }
//     for (int i = 0; i < NDBC; i++) {
//         for (int j = 0; j < 2; j++) {
//             DB12(i,j) = DBCS(i,j);
//             DBiter(i,j) = DBCS(i,j);
//         }
//     }

//     // finding global location of gauss points
//     CArray <double> gpcoords = gpglob(NEL,gp,elcoords,uel,NODES,CONN,ut,us,gpV0,gpV1,gpV2,gpV3,gpV4,gpV5,gpV6,gpV7);

//     // Opening Output file
//     std::ofstream outputFile("output.txt");

//     // Writing initial configuration to output file
//     outputFile << "Initial Configuration:" << "\n" << "Node     x     y     z" << "\n";
//     for (int i = 0; i < NNODE; i++) {
//         outputFile << i << "     " << NODES(i,0) << "     " << NODES(i,1) << "     " << NODES(i,2) << "\n";
//     }
//     outputFile << "\n";

//     // Initializing numerator and denominator of enorm equation for simpler calculation
//     double nnorm = 0;
//     double dnorm = 0;

//     // Defining max allowed iterations for a single load step and euclidean norm
//     // Euclidean norm of displacement is initialized at greater than tolerance
//     // so that the first check does not fail and the convergence loop will run
//     int maxiter = 100;
//     double enorm = 2 * tol;

//     // Beginning of load step loop
//     for (int i = 0; i < NLS; i++) {
//         // Updating boundary conditions for the current step from 01 to 02
//         if (BCFLAG == 0) {
//             for (int j = 0; j < NDBC; j++) {
//                 DB12(j,2) = DBstep(j,0);
//             }
//             for (int j = 0; j < NPL; j++) {
//                 PL02(j,2) += PLstep(j,0);
//             }
//             for (int j = 0; j < NTL; j++) {
//                 TL02(j,2) += TLstep(j,0);
//                 TL02(j,3) += TLstep(j,1);
//                 TL02(j,4) += TLstep(j,2);
//             }
//         }
//         else {
//             for (int j = 0; j < NDBC; j++) {
//                 DB12(j,2) = DBstep(j,i);
//             }
//             for (int j = 0; j < NPL; j++) {
//                 PL02(j,2) += PLstep(j,i);
//             }
//             for (int j = 0; j < NTL; j++) {
//                 for (int k = 0; k < 3; k++) {
//                     TL02(j,2 + k) += TLstep(j,3 * i + k);
//                 }
//             }
//         }

//         // Applying the loads for this load step to define F02g for the load step
//         PointLoad(F02g, PL02, NPL);
//         for (int j = 0; j < NTL; j++) {
//             ElemCoords(elcoords, uel, static_cast<int>(TL02(j,0)), NODES, CONN, ut, us);
//             switch(static_cast<int>(TL02(j,1))) {
//             case 0:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp00, gpp01, gpp02, gpp03);
//                 break;
//             case 1:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp10, gpp11, gpp12, gpp13);
//                 break;
//             case 2:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp20, gpp21, gpp22, gpp23);
//                 break;
//             case 3:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp30, gpp31, gpp32, gpp33);
//                 break;
//             case 4:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp40, gpp41, gpp42, gpp43);
//                 break;
//             case 5:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp50, gpp51, gpp52, gpp53);
//                 break;
//             }
//         }

//         // resetting dinvars for the load step
//         dinvars.set_values(0);
        
//         // initializing cohesive zone vector and array for first iteration
//         vczupdate(NODES,CONN,NEL,ut,us,UPs,vczinfo,NUP,interpvals,vczorient,ulocvcz,NPT,Einf,a1,n,Eandrhom,uns,uts,delt,invars,dinvars,Fvcz,Kvcz);
        
//         // enforcing dirichlet conditions on us
//         for (int k = 0; k < NDBC; k++) {
//             us(3 * static_cast<int>(DB12(k,0)) + static_cast<int>(DB12(k,1))) = DB12(k,2);
//         }

//         // Resetting Euclidean norm before initiating convergence check loop so that
//         // the subsequent load steps after the first do not fail to run
//         enorm = 2 * tol;

//         // Beginning of convergence check loop
//         for (int j = 0; j < maxiter; j++) {
//             // Checking for convergence
//             if (enorm < tol) {
//                 // uncomment to print iteration number to console
//                 //printf("j is %i\n", j);
//                 break;
//             }
            
//             //printf("LS is %i\nIT is %i\n\n",i,j);

//             // Resetting global arrays
//             Kg.set_values(0);
//             F01g.set_values(0);
//             Fg.set_values(0);

//             // Resetting enorm variables
//             nnorm = 0;
//             dnorm = 0;

//             // Beginning of element loop
//             for (int k = 0; k < NEL; k++) {
//                 // Pulling element geometry
//                 ElemCoords(elcoords, uel, k, NODES, CONN, ut, us);

//                 // Resetting element arrays
//                 Kel.set_values(0);
//                 F01el.set_values(0);

//                 // Calculating element matrices
//                 // Each gauss point is added to overall element matrices and so doesn't need to be reset between function calls
//                 Gradients(dpsig, gradu, detJ, gpV0, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV1, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV2, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV3, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV4, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV5, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV6, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV7, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
                
//                 // Assembly to global matrices
//                 for (int m = 0; m < 8; m++) {
//                     for (int o = 0; o < 8; o++) {
//                         Kg(3*CONN(k,m),3*CONN(k,o)) += Kel(3*m,3*o);
//                         Kg(3*CONN(k,m),3*CONN(k,o)+1) += Kel(3*m,3*o + 1);
//                         Kg(3*CONN(k,m),3*CONN(k,o)+2) += Kel(3*m,3*o + 2);
//                         Kg(3*CONN(k,m)+1,3*CONN(k,o)) += Kel(3*m+1,3*o);
//                         Kg(3*CONN(k,m)+1,3*CONN(k,o)+1) += Kel(3*m+1,3*o + 1);
//                         Kg(3*CONN(k,m)+1,3*CONN(k,o)+2) += Kel(3*m+1,3*o + 2);
//                         Kg(3*CONN(k,m)+2,3*CONN(k,o)) += Kel(3*m+2,3*o);
//                         Kg(3*CONN(k,m)+2,3*CONN(k,o)+1) += Kel(3*m+2,3*o + 1);
//                         Kg(3*CONN(k,m)+2,3*CONN(k,o)+2) += Kel(3*m+2,3*o + 2);
//                     }
//                 }
//                 for (int m = 0; m < 8; m++) {
//                     F01g(3*CONN(k,m)) += F01el(3*m);
//                     F01g(3*CONN(k,m) + 1) += F01el(3*m + 1);
//                     F01g(3*CONN(k,m) + 2) += F01el(3*m + 2);
//                 }
//             }

//             // adding vcz stiffness contributions
//             for (int k = 0; k < 3*NNODE; k++) {
//                 for (int m = 0; m < 3*NNODE; m++) {
//                     Kg(k,m) += Kvcz(k,m);
//                 }
//             }
            
//             // Forming full Fg
//             for (int k = 0; k < 3 * NNODE; k++) {
//                 Fg(k) = F02g(k) - F01g(k) + Fvcz(k);
//             }
            
//             // Applying Dirichlet boundary conditions
//             Dirichlet(Kg,Fg,DBiter,NDBC,NNODE);

//             // Solving system of equations
//             GaussElim(Kg,Fg,dus);
            
//             // updating us
//             for (int k = 0; k < 3 * NNODE; k++) {
//                 us(k) += dus(k);
//             }

//             // updating vcz for next iteration
//             vczupdate(NODES,CONN,NEL,ut,us,UPs,vczinfo,NUP,interpvals,vczorient,ulocvcz,NPT,Einf,a1,n,Eandrhom,uns,uts,delt,invars,dinvars,Fvcz,Kvcz);

//             // solving for error norm
//             for (int k = 0; k < 3 * NNODE; k++) {
//                 nnorm += dus(k)*dus(k);
//                 dnorm += us(k)*us(k);
//             }
//             enorm = sqrt(nnorm/dnorm);

//         }
//         // updating total displacement
//         for (int k = 0; k < 3 * NNODE; k++) {
//             ut(k) += us(k);
//         }

//         // reset displacement step and force vector that holds tractions and point loads
//         us.set_values(0);
//         F02g.set_values(0);

//         // Updating internal variables post load step convergence
//         for (int j = 0; j < NUP; j++) {
//             invars(j,0) = dinvars(j,0);
//             invars(j,1) += dinvars(j,1);
//             invars(j,2) += dinvars(j,2);
//             invars(j,3) += dinvars(j,3);
//             for (int k = 0; k < NPT; k++) {
//                 invars(j,4 + k) = dinvars(j,4 + k);
//             }
//         }

//         // printing cohesive zone values for output checking (probably should add to output file in so way at some point?)
//         /* for (int j = 0; j < 1; j++) {
//             printf("%f     %f\n",invars(j,1),invars(j,2));
//         } */
//         printf("%f     %f\n",invars(0,1),invars(0,2));
//         //std::cout << invars(0,1) << "     " << invars(0,2) << std::endl;
//         //std::cout << invars.dims(0) << "     " << invars.dims(1) << std::endl;
//         //printf("\n");
//         //prarr(invars);


//         // Post processing and writing to output
//         // Outputting displacement field
//         outputFile << "Load Step " << i << ":" << "\n";
//         outputFile << "Total Displacement:" << "\n" << "Node     x           y              z" << "\n";
//         for (int j = 0; j < NNODE; j++) {
//             outputFile << j << "     " << std::setprecision(2) << std::scientific << ut(3 * j) << "     " << ut(3 * j + 1) << "     " << ut(3 * j + 2) << "\n";
//         }
//         // Outputting strain field
//         outputFile << "\n" << "x          y          z          Exx           Eyy            Ezz            Exy          Exz         Eyz" << "\n";
//         for (int j = 0; j < NEL; j++) {
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV0,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,0) << "     " << gpcoords(j,1) << "     " << gpcoords(j,2) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV1,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,3) << "     " << gpcoords(j,4) << "     " << gpcoords(j,5) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV2,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,6) << "     " << gpcoords(j,7) << "     " << gpcoords(j,8) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV3,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,9) << "     " << gpcoords(j,10) << "     " << gpcoords(j,11) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV4,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,12) << "     " << gpcoords(j,13) << "     " << gpcoords(j,14) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV5,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,15) << "     " << gpcoords(j,16) << "     " << gpcoords(j,17) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV6,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,18) << "     " << gpcoords(j,19) << "     " << gpcoords(j,20) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV7,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,21) << "     " << gpcoords(j,22) << "     " << gpcoords(j,23) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//         }
//         // Outputting stress field
//         outputFile << "\n" << "x          y          z          Sxx           Syy            Szz             Sxy           Sxz          Syz" << "\n";
//         for (int j = 0; j < NEL; j++) {
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV0,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,0) << "     " << gpcoords(j,1) << "     " << gpcoords(j,2) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV1,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,3) << "     " << gpcoords(j,4) << "     " << gpcoords(j,5) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV2,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,6) << "     " << gpcoords(j,7) << "     " << gpcoords(j,8) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV3,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,9) << "     " << gpcoords(j,10) << "     " << gpcoords(j,11) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV4,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,12) << "     " << gpcoords(j,13) << "     " << gpcoords(j,14) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV5,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,15) << "     " << gpcoords(j,16) << "     " << gpcoords(j,17) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV6,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,18) << "     " << gpcoords(j,19) << "     " << gpcoords(j,20) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV7,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,21) << "     " << gpcoords(j,22) << "     " << gpcoords(j,23) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//         }
//     }
// }
// **************************************************************** FROM GAVIN'S CODE **************************************************************** 