// ==============================================================
// Voxelization code for binary stl files
// ==============================================================
// ABOUT: This code is a Matarized version of the original code written by:
// Adam A (2023). Mesh voxelisation (https://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation), MATLAB Central File Exchange. Retrieved June 28, 2023.

#include <stdio.h>
#include <array>
#include <variant>
#include <chrono>
#include <fstream>

#include <math.h>

#include <matar.h>
#include "stl-to-voxelvtk.h"
#include <chrono>

using namespace mtr; // matar namespace

// Functions used within MAIN
std::tuple<CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, unsigned int> binary_stl_reader(std::string stl_file_path); // BINARY STL READER FUNCTION
void main_function(CArray<bool> &gridOUTPUT, int &gridX, int &gridY, int &gridZ, CArray<float> &normal, CArray<float> &v1X, CArray<float> &v1Y, CArray<float> &v1Z, CArray<float> &v2X, CArray<float> &v2Y, CArray<float> &v2Z, CArray<float> &v3X, CArray<float> &v3Y, CArray<float> &v3Z, unsigned int &n_facets, double &voxel_dx, double &voxel_dy, double &voxel_dz); // VOXELIZATION FUNCTION

// ==============================================================
// ---------------------------- MAIN ----------------------------
// ==============================================================
std::tuple<CArray<bool>, double, double, double> Voxelizer::create_voxel_vtk(std::string stl_file_path, std::string vtk_file_path, int gridX, int gridY, int gridZ, double length_x, double length_y, double length_z) {
    
//    Kokkos::initialize(argc, argv);
//    {
        // Start Clock
        auto start = std::chrono::high_resolution_clock::now();
            
        // ***** USER-DEFINED INPUTS *****
        //std::string stl_file_path = argv[1]; // location of BINARY "stl" file
        //std::string vtk_file_path = argv[2];
        //int gridX = atoi(argv[3]); // voxel resolution x-direction
        //int gridY = atoi(argv[4]); // voxel resolution y-direction
        //int gridZ = atoi(argv[5]); // voxel resolution z-direction
        
        // Read the stl file
        auto [normal, v1X, v1Y, v1Z, v2X, v2Y, v2Z, v3X, v3Y, v3Z, n_facets] = binary_stl_reader(stl_file_path);
        
//        v1X.update_device();
//        v1Y.update_device();
//        v1Z.update_device();
//        v2X.update_device();
//        v2Y.update_device();
//        v2Z.update_device();
//        v3X.update_device();
//        v3Y.update_device();
//        v3Z.update_device();
//        Kokkos::fence();
        
        // Voxelize the part in the x,y,z directions
        double voxel_dx;
        double voxel_dy;
        double voxel_dz;
    
        // X-direction voxelization
        CArray<bool> gridOUTPUTX(gridY+2,gridZ+2,gridX+2);
        FOR_LOOP (i, 0, gridY+2,
                 j, 0, gridZ+2,
                 k, 0, gridX+2,{
            gridOUTPUTX(i,j,k) = 0;
        });
//        Kokkos::fence();
        main_function(gridOUTPUTX, gridY, gridZ, gridX, normal, v1Y, v1Z, v1X, v2Y, v2Z, v2X, v3Y, v3Z, v3X, n_facets, voxel_dx, voxel_dy, voxel_dz);
        main_function(gridOUTPUTX, gridY, gridZ, gridX, normal, v1Y, v1Z, v1X, v2Y, v2Z, v2X, v3Y, v3Z, v3X, n_facets, voxel_dx, voxel_dy, voxel_dz);
        // Y-direction voxelization
        CArray<bool> gridOUTPUTY(gridZ+2,gridX+2,gridY+2);
        FOR_LOOP (i, 0, gridZ+2,
                 j, 0, gridX+2,
                 k, 0, gridY+2,{
            gridOUTPUTY(i,j,k) = 0;
        });
//        Kokkos::fence();
        main_function(gridOUTPUTY, gridZ, gridX, gridY, normal, v1Z, v1X, v1Y, v2Z, v2X, v2Y, v3Z, v3X, v3Y, n_facets, voxel_dx, voxel_dy, voxel_dz);
        
        // Z-direction voxelization
        CArray<bool> gridOUTPUTZ(gridX+2,gridY+2,gridZ+2);
        FOR_LOOP (i, 0, gridX+2,
                 j, 0, gridY+2,
                 k, 0, gridZ+2,{
            gridOUTPUTZ(i,j,k) = 0;
        });
//        Kokkos::fence();
        main_function(gridOUTPUTZ, gridX, gridY, gridZ, normal, v1X, v1Y, v1Z, v2X, v2Y, v2Z, v3X, v3Y, v3Z, n_facets, voxel_dx, voxel_dy, voxel_dz);

        // Sum the voxelization in the XYZ-directions
        CArray<int> TOTgrid(gridX+2,gridY+2,gridZ+2);
        FOR_LOOP(i,0,gridX+2,
                j,0,gridY+2,
                k,0,gridZ+2,{
            TOTgrid(i,j,k) = gridOUTPUTX(j,k,i) + gridOUTPUTY(k,i,j) + gridOUTPUTZ(i,j,k);
        });
//        Kokkos::fence();
        CArray<bool> OUTPUTgrid(gridZ,gridY,gridX);
        FOR_LOOP(i,1,gridX+1,
                j,1,gridY+1,
                k,1,gridZ+1,{
            OUTPUTgrid(k-1,j-1,i-1) = TOTgrid(i,j,k)>=2;
        });
//        Kokkos::fence();

        // VTK WRITER
        if (length_x*length_y*length_z > 0) {
            voxel_dx = length_x/gridX;
            voxel_dy = length_y/gridY;
            voxel_dz = length_z/gridZ;
        }
    
        int i,j,k;
        const char* cvtk_file_path = vtk_file_path.c_str(); // convert std::string to C-style string
        auto out=std::fopen(cvtk_file_path,"w"); // open the file

        fprintf(out,"# vtk DataFile Version 3.0\n"); // write the header
        fprintf(out,"Header\n");
        fprintf(out,"ASCII\n");
        fprintf(out,"DATASET RECTILINEAR_GRID\n");
        fprintf(out,"DIMENSIONS %d %d %d\n",gridX+1,gridY+1,gridZ+1);
        
        fprintf(out,"X_COORDINATES %d float\n", gridX+1); // nodal x coordinates
        for (i=0; i<(gridX+1); i++) {
            fprintf(out,"%lf ",double(i)*double(voxel_dx));
        }
        fprintf(out,"\nY_COORDINATES %d float\n", gridY+1); // nodal y coordinates
        for (i=0; i<(gridY+1); i++) {
            fprintf(out,"%lf ",double(i)*double(voxel_dy));
        }
        fprintf(out,"\nZ_COORDINATES %d float\n", gridZ+1); // nodal z coordinates
        for (i=0; i<(gridZ+1); i++) {
            fprintf(out,"%lf ",double(i)*double(voxel_dz));
        }
        
        fprintf(out,"\n\nCELL_DATA %d\n",gridX*gridY*gridZ); // material (1) or void (0) region definition
        fprintf(out,"SCALARS density float 1\n");
        fprintf(out,"LOOKUP_TABLE default\n");
        for (i=0;i<gridZ;i++){
            for (j=0;j<gridY;j++){
                for (k=0;k<gridX;k++){
                    fprintf(out,"%d ",OUTPUTgrid(i,j,k));
                }
            }
        }

        fclose(out);
        
        // END CLOCK
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        std::cout << "Total Time: " << duration.count() << " milliseconds" << std::endl;

//    } // end of kokkos scope

//    Kokkos::finalize();

    printf("\nfinished\n\n");
    
    
    return {OUTPUTgrid, voxel_dx, voxel_dy, voxel_dz};
}

std::tuple<double, double, double> Voxelizer::create_voxel_vtk_GUI(std::string stl_file_path, std::string vtk_file_path, int gridX, int gridY, int gridZ, double origin_x, double origin_y, double origin_z, double length_x, double length_y, double length_z) {
    
//    Kokkos::initialize(argc, argv);
//    {
        // Start Clock
        auto start = std::chrono::high_resolution_clock::now();
            
        // ***** USER-DEFINED INPUTS *****
        //std::string stl_file_path = argv[1]; // location of BINARY "stl" file
        //std::string vtk_file_path = argv[2];
        //int gridX = atoi(argv[3]); // voxel resolution x-direction
        //int gridY = atoi(argv[4]); // voxel resolution y-direction
        //int gridZ = atoi(argv[5]); // voxel resolution z-direction
        
        // Read the stl file
        auto [normal, v1X, v1Y, v1Z, v2X, v2Y, v2Z, v3X, v3Y, v3Z, n_facets] = binary_stl_reader(stl_file_path);
        
//        v1X.update_device();
//        v1Y.update_device();
//        v1Z.update_device();
//        v2X.update_device();
//        v2Y.update_device();
//        v2Z.update_device();
//        v3X.update_device();
//        v3Y.update_device();
//        v3Z.update_device();
//        Kokkos::fence();
        
        // Voxelize the part in the x,y,z directions
        double voxel_dx;
        double voxel_dy;
        double voxel_dz;
    
        // X-direction voxelization
        CArray<bool> gridOUTPUTX(gridY+2,gridZ+2,gridX+2);
        FOR_LOOP (i, 0, gridY+2,
                 j, 0, gridZ+2,
                 k, 0, gridX+2,{
            gridOUTPUTX(i,j,k) = 0;
        });
//        Kokkos::fence();
        main_function(gridOUTPUTX, gridY, gridZ, gridX, normal, v1Y, v1Z, v1X, v2Y, v2Z, v2X, v3Y, v3Z, v3X, n_facets, voxel_dx, voxel_dy, voxel_dz);
        // Y-direction voxelization
        CArray<bool> gridOUTPUTY(gridZ+2,gridX+2,gridY+2);
        FOR_LOOP (i, 0, gridZ+2,
                 j, 0, gridX+2,
                 k, 0, gridY+2,{
            gridOUTPUTY(i,j,k) = 0;
        });
//        Kokkos::fence();
        main_function(gridOUTPUTY, gridZ, gridX, gridY, normal, v1Z, v1X, v1Y, v2Z, v2X, v2Y, v3Z, v3X, v3Y, n_facets, voxel_dx, voxel_dy, voxel_dz);
        
        // Z-direction voxelization
        CArray<bool> gridOUTPUTZ(gridX+2,gridY+2,gridZ+2);
        FOR_LOOP (i, 0, gridX+2,
                 j, 0, gridY+2,
                 k, 0, gridZ+2,{
            gridOUTPUTZ(i,j,k) = 0;
        });
//        Kokkos::fence();
        main_function(gridOUTPUTZ, gridX, gridY, gridZ, normal, v1X, v1Y, v1Z, v2X, v2Y, v2Z, v3X, v3Y, v3Z, n_facets, voxel_dx, voxel_dy, voxel_dz);

        // Sum the voxelization in the XYZ-directions
        CArray<int> TOTgrid(gridX+2,gridY+2,gridZ+2);
        FOR_LOOP(i,0,gridX+2,
                j,0,gridY+2,
                k,0,gridZ+2,{
            TOTgrid(i,j,k) = gridOUTPUTX(j,k,i) + gridOUTPUTY(k,i,j) + gridOUTPUTZ(i,j,k);
        });
//        Kokkos::fence();
        CArray<bool> OUTPUTgrid(gridZ,gridY,gridX);
        FOR_LOOP(i,1,gridX+1,
                j,1,gridY+1,
                k,1,gridZ+1,{
            OUTPUTgrid(k-1,j-1,i-1) = TOTgrid(i,j,k)>=2;
        });
//        Kokkos::fence();

        // VTK WRITER
        if (length_x*length_y*length_z > 0) {
            voxel_dx = length_x/gridX;
            voxel_dy = length_y/gridY;
            voxel_dz = length_z/gridZ;
        }
    
        int i,j,k;
        const char* cvtk_file_path = vtk_file_path.c_str(); // convert std::string to C-style string
        auto out=std::fopen(cvtk_file_path,"w"); // open the file

        fprintf(out,"# vtk DataFile Version 3.0\n"); // write the header
        fprintf(out,"Header\n");
        fprintf(out,"ASCII\n");
        fprintf(out,"DATASET RECTILINEAR_GRID\n");
        fprintf(out,"DIMENSIONS %d %d %d\n",gridX+1,gridY+1,gridZ+1);
        
        fprintf(out,"X_COORDINATES %d float\n", gridX+1); // nodal x coordinates
        for (i=0; i<(gridX+1); i++) {
            fprintf(out,"%lf ",double(i)*double(voxel_dx)+origin_x);
        }
        fprintf(out,"\nY_COORDINATES %d float\n", gridY+1); // nodal y coordinates
        for (i=0; i<(gridY+1); i++) {
            fprintf(out,"%lf ",double(i)*double(voxel_dy)+origin_y);
        }
        fprintf(out,"\nZ_COORDINATES %d float\n", gridZ+1); // nodal z coordinates
        for (i=0; i<(gridZ+1); i++) {
            fprintf(out,"%lf ",double(i)*double(voxel_dz)+origin_z);
        }
        
        fprintf(out,"\n\nCELL_DATA %d\n",gridX*gridY*gridZ); // material (1) or void (0) region definition
        fprintf(out,"SCALARS density float 1\n");
        fprintf(out,"LOOKUP_TABLE default\n");
        for (i=0;i<gridZ;i++){
            for (j=0;j<gridY;j++){
                for (k=0;k<gridX;k++){
                    fprintf(out,"%d ",OUTPUTgrid(i,j,k));
                }
            }
        }

        fclose(out);
        
        // END CLOCK
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        std::cout << "Total Time: " << duration.count() << " milliseconds" << std::endl;

//    } // end of kokkos scope

//    Kokkos::finalize();

    printf("\nfinished\n\n");
    
    
    return {voxel_dx, voxel_dy, voxel_dz};
}





// ==============================================================
// ------------------- VOXELIZATION FUNCTIONS -------------------
// ==============================================================

// BINARY STL READER - (Note: it can ONLY read binary stl files)
std::tuple<CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, unsigned int> binary_stl_reader(std::string stl_file_path){
    // Open the binary file
    std::string filename = stl_file_path;
    auto input = std::ifstream{filename, std::ifstream::in | std::ifstream::binary};
    
    // Check that the file is actually open
    if (input.is_open()) {
        std::cout << "Opening .stl file... \n";
    }
    else {
        std::cout << "WARNING: .stl file is not open \n";
    }
    
    // Initialize variables
    char heading[81];
    unsigned int n_facets;
    
    // Read the file's heading
    char* ptr1 = heading;
    input.read(ptr1, 80);
    heading[80] = '\0'; // to ensure a proper C string
    std::cout << "File heading: " << heading << "\n";
    
    // Read the number of facets in the file
    unsigned int* ptr2 = &n_facets;
    input.read(reinterpret_cast<char*>(ptr2), sizeof(unsigned int));
    std::cout << "Voxelizing " << n_facets << " facets..." << "\n";
    
    // Initialize storage variables
    CArray<float> normal(n_facets,3);
    CArray<float> v1X(n_facets);
    CArray<float> v1Y(n_facets);
    CArray<float> v1Z(n_facets);
    CArray<float> v2X(n_facets);
    CArray<float> v2Y(n_facets);
    CArray<float> v2Z(n_facets);
    CArray<float> v3X(n_facets);
    CArray<float> v3Y(n_facets);
    CArray<float> v3Z(n_facets);
    float normalp[3];
    float v1p[3];
    float v2p[3];
    float v3p[3];
    
    // Pull data from file
    float* ptr3 = normalp;
    float* ptr4 = v1p;
    float* ptr5 = v2p;
    float* ptr6 = v3p;
    for (size_t i = 0; i < n_facets; ++i) {
        size_t n_bytes = 3 * sizeof(float);
        input.read(reinterpret_cast<char*>(ptr3), n_bytes);
        input.read(reinterpret_cast<char*>(ptr4), n_bytes);
        input.read(reinterpret_cast<char*>(ptr5), n_bytes);
        input.read(reinterpret_cast<char*>(ptr6), n_bytes);
        input.seekg(2,input.cur);
        for (int j=0; j<3; j++) {
            normal(i,j) = normalp[j];
        }
        v1X(i) = v1p[0];
        v1Y(i) = v1p[1];
        v1Z(i) = v1p[2];
        v2X(i) = v2p[0];
        v2Y(i) = v2p[1];
        v2Z(i) = v2p[2];
        v3X(i) = v3p[0];
        v3Y(i) = v3p[1];
        v3Z(i) = v3p[2];
//        v1X.host(i) = v1p[0];
//        v1Y.host(i) = v1p[1];
//        v1Z.host(i) = v1p[2];
//        v2X.host(i) = v2p[0];
//        v2Y.host(i) = v2p[1];
//        v2Z.host(i) = v2p[2];
//        v3X.host(i) = v3p[0];
//        v3Y.host(i) = v3p[1];
//        v3Z.host(i) = v3p[2];
    }
    input.close();
    return {normal, v1X, v1Y, v1Z, v2X, v2Y, v2Z, v3X, v3Y, v3Z, n_facets};
}

// VOXELIZATION FUNCTION
void main_function(CArray<bool> &gridOUTPUT, int &gridX, int &gridY, int &gridZ, CArray<float> &normal, CArray<float> &v1X, CArray<float> &v1Y, CArray<float> &v1Z, CArray<float> &v2X, CArray<float> &v2Y, CArray<float> &v2Z, CArray<float> &v3X, CArray<float> &v3Y, CArray<float> &v3Z, unsigned int &n_facets, double &voxel_dx, double &voxel_dy, double &voxel_dz){
    // Find the global maximum and minimum values of the mesh
    float meshXmax;
    float meshXmin;
    float meshYmax;
    float meshYmin;
    float meshZmax;
    float meshZmin;
    
    // Global maximum x-direction
    FOR_REDUCE_MAX(i,0,n_facets,meshXmax, {
        if (v1X(i) > meshXmax | v2X(i) > meshXmax | v3X(i) > meshXmax) {
            if (v1X(i) > v2X(i) && v1X(i) > v3X(i)) {
                meshXmax = v1X(i);
            } else if (v2X(i) > v1X(i) && v2X(i) > v3X(i)) {
                meshXmax = v2X(i);
            } else if (v3X(i) > v1X(i) && v3X(i) > v2X(i)) {
                meshXmax = v3X(i);
            }
        }
    }, meshXmax);
    
    // Global minimum x-direction
    FOR_REDUCE_MIN(i,0,n_facets,meshXmin, {
        if (v1X(i) < meshXmin | v2X(i) < meshXmin | v3X(i) < meshXmin) {
            if (v1X(i) < v2X(i) && v1X(i) < v3X(i)) {
                meshXmin = v1X(i);
            } else if (v2X(i) < v1X(i) && v2X(i) < v3X(i)) {
                meshXmin = v2X(i);
            } else if (v3X(i) < v1X(i) && v3X(i) < v2X(i)) {
                meshXmin = v3X(i);
            }
        }
    }, meshXmin);

    // Global maximum y-direction
    FOR_REDUCE_MAX(i,0,n_facets,meshYmax, {
        if (v1Y(i) > meshYmax | v2Y(i) > meshYmax | v3Y(i) > meshYmax) {
            if (v1Y(i) > v2Y(i) && v1Y(i) > v3Y(i)) {
                meshYmax = v1Y(i);
            } else if (v2Y(i) > v1Y(i) && v2Y(i) > v3Y(i)) {
                meshYmax = v2Y(i);
            } else if (v3Y(i) > v1Y(i) && v3Y(i) > v2Y(i)) {
                meshYmax = v3Y(i);
            }
        }
    }, meshYmax);

    // Global minimum y-direction
    FOR_REDUCE_MIN(i,0,n_facets,meshYmin, {
        if (v1Y(i) < meshYmin | v2Y(i) < meshYmin | v3Y(i) < meshYmin) {
            if (v1Y(i) < v2Y(i) && v1Y(i) < v3Y(i)) {
                meshYmin = v1Y(i);
            } else if (v2Y(i) < v1Y(i) && v2Y(i) < v3Y(i)) {
                meshYmin = v2Y(i);
            } else if (v3Y(i) < v1Y(i) && v3Y(i) < v2Y(i)) {
                meshYmin = v3Y(i);
            }
        }
    }, meshYmin);

    // Global maximum z-direction
    FOR_REDUCE_MAX(i,0,n_facets,meshZmax, {
        if (v1Z(i) > meshZmax | v2Z(i) > meshZmax | v3Z(i) > meshZmax) {
            if (v1Z(i) > v2Z(i) && v1Z(i) > v3Z(i)) {
                meshZmax = v1Z(i);
            } else if (v2Z(i) > v1Z(i) && v2Z(i) > v3Z(i)) {
                meshZmax = v2Z(i);
            } else if (v3Z(i) > v1Z(i) && v3Z(i) > v2Z(i)) {
                meshZmax = v3Z(i);
            }
        }
    }, meshZmax);

    // Global minimum z-direction
    FOR_REDUCE_MIN(i,0,n_facets,meshZmin, {
        if (v1Z(i) < meshZmin | v2Z(i) < meshZmin | v3Z(i) < meshZmin) {
            if (v1Z(i) < v2Z(i) && v1Z(i) < v3Z(i)) {
                meshZmin = v1Z(i);
            } else if (v2Z(i) < v1Z(i) && v2Z(i) < v3Z(i)) {
                meshZmin = v2Z(i);
            } else if (v3Z(i) < v1Z(i) && v3Z(i) < v2Z(i)) {
                meshZmin = v3Z(i);
            }
        }
    }, meshZmin);
//    Kokkos::fence();

    // Define the grid that the voxel mesh will be generated on
//    float voxwidth;
    
    // Voxel grid x-direction
    voxel_dx = (meshXmax-meshXmin)/(gridX);
    int elx = gridX+2;
    CArray<float> gridCOx(elx+1);
    gridCOx(0) = meshXmin;
    gridCOx(elx) = meshXmax;
    FOR_LOOP (i, 0, elx-1, {
        gridCOx(i+1) = (meshXmin+(voxel_dx/2)) + voxel_dx*i;
    });
//    Kokkos::fence();

    // Voxel grid y-direction
    voxel_dy = (meshYmax-meshYmin)/(gridY);
    int ely = gridY+2;
    CArray<float> gridCOy(ely+1);
    gridCOy(0) = meshYmin;
    gridCOy(ely) = meshYmax;
    FOR_LOOP(i, 0, ely-1, {
        gridCOy(i+1) = (meshYmin+(voxel_dy/2)) + voxel_dy*i;
    });
//    Kokkos::fence();

    // Voxel grid z-direction
    voxel_dz = (meshZmax-meshZmin)/(gridZ);
    int elz = gridZ+2;
    CArray<float> gridCOz(elz+1);
    gridCOz(0) = meshZmin;
    gridCOz(elz) = meshZmax;
    FOR_LOOP(i, 0, elz-1, {
        gridCOz(i+1) = (meshZmin+(voxel_dz/2)) + voxel_dz*i;
    });
//    Kokkos::fence();

    // Z-DIRECTIONAL RAY EXECUTION
    // Get the minimum and maximum x,y,z coordinates for each facet
    CArray<float> facetXmin(n_facets);
    CArray<float> facetXmax(n_facets);
    CArray<float> facetYmin(n_facets);
    CArray<float> facetYmax(n_facets);
    CArray<float> facetZmin(n_facets);
    CArray<float> facetZmax(n_facets);
    
    FOR_LOOP (i,0,n_facets, {
        // Facet minimum x-direction
        if (v1X(i) <= v2X(i) && v1X(i) <= v3X(i)) {
            facetXmin(i) = v1X(i);
        } else if (v2X(i) <= v1X(i) && v2X(i) <= v3X(i)) {
            facetXmin(i) = v2X(i);
        } else if (v3X(i) <= v1X(i) && v3X(i) <= v2X(i)) {
            facetXmin(i) = v3X(i);
        }

        // Facet maximum x-direction
        if (v1X(i) >= v2X(i) && v1X(i) >= v3X(i)) {
            facetXmax(i) = v1X(i);
        } else if (v2X(i) >= v1X(i) && v2X(i) >= v3X(i)) {
            facetXmax(i) = v2X(i);
        } else if (v3X(i) >= v1X(i) && v3X(i) >= v2X(i)) {
            facetXmax(i) = v3X(i);
        }

        // Facet minimum y-direction
        if (v1Y(i) <= v2Y(i) && v1Y(i) <= v3Y(i)) {
            facetYmin(i) = v1Y(i);
        } else if (v2Y(i) <= v1Y(i) && v2Y(i) <= v3Y(i)) {
            facetYmin(i) = v2Y(i);
        } else if (v3Y(i) <= v1Y(i) && v3Y(i) <= v2Y(i)) {
            facetYmin(i) = v3Y(i);
        }

        // Facet maximum y-direction
        if (v1Y(i) >= v2Y(i) && v1Y(i) >= v3Y(i)) {
            facetYmax(i) = v1Y(i);
        } else if (v2Y(i) >= v1Y(i) && v2Y(i) >= v3Y(i)) {
            facetYmax(i) = v2Y(i);
        } else if (v3Y(i) >= v1Y(i) && v3Y(i) >= v2Y(i)) {
            facetYmax(i) = v3Y(i);
        }

        // Facet minimum z-direction
        if (v1Z(i) <= v2Z(i) && v1Z(i) <= v3Z(i)) {
            facetZmin(i) = v1Z(i);
        } else if (v2Z(i) <= v1Z(i) && v2Z(i) <= v3Z(i)) {
            facetZmin(i) = v2Z(i);
        } else if (v3Z(i) <= v1Z(i) && v3Z(i) <= v2Z(i)) {
            facetZmin(i) = v3Z(i);
        }

        // Facet maximum z-direction
        if (v1Z(i) >= v2Z(i) && v1Z(i) >= v3Z(i)) {
            facetZmax(i) = v1Z(i);
        } else if (v2Z(i) >= v1Z(i) && v2Z(i) >= v3Z(i)) {
            facetZmax(i) = v2Z(i);
        } else if (v3Z(i) >= v1Z(i) && v3Z(i) >= v2Z(i)) {
            facetZmax(i) = v3Z(i);
        }
    });
//    Kokkos::fence();
    
    // Create a list to record all rays that fail to voxelize
    CArray<int> correctionLIST(ely*elx,2);
    CArray<int> facetCROSSLISTcounter(1);
    CArray<int> correctionLISTcounter(1);
    
    // Find which facets could be crossed by the ray in the y-direction
    CArray<int> possibleCROSSLISTy(ely*n_facets,ely);
    CArray<int> counter(1);
    CArray<int> count1(ely);
    counter(0) = 0;
    FOR_LOOP (loopY,0,ely,
             i,0,n_facets,{
        if (facetYmin(i) <= gridCOy(loopY) && facetYmax(i) >= gridCOy(loopY)) {
            possibleCROSSLISTy(counter(0),loopY) = i;
            counter(0)++;
//            Kokkos::atomic_add(&counter(0),1);
        }
        if (i == n_facets-1) {
            count1(loopY) = counter(0);
            counter(0) = 0;
        }
    });
//    Kokkos::fence();
    
    // Initialize variables
    CArray<int> possibleCROSSLIST(n_facets);
    CArray<int> facetCROSSLIST(n_facets);
    CArray<int> vertexCROSSLIST(n_facets);
    CArray<float> coN(n_facets);
    CArray<float> gridCOzCROSS(n_facets);
    CArray<float> gridCOzCROSSunique(n_facets);
    CArray<bool> voxelsINSIDE(elz);
    
    // Loop through each pixel in the x-y plane by passing rays in the z-direction and finding where they cross the voxelized mesh
    correctionLISTcounter(0) = 0;
    FOR_LOOP (loopY,0,ely,
            loopX,0,elx,{
        facetCROSSLISTcounter(0) = 0;

        // Find which facets could be crossed by the ray in the x-direction
        int count2 = 0;
        for (int i = 0; i < count1(loopY); i++) {
            if (facetXmin(possibleCROSSLISTy(i,loopY)) <= gridCOx(loopX) && facetXmax(possibleCROSSLISTy(i,loopY)) >= gridCOx(loopX)) {
                possibleCROSSLIST(count2) = possibleCROSSLISTy(i,loopY);
                count2++;
            }
        }
        
        // Check if the ray actually passed through the facet
        if (count2 > 0) {
            
            // For each facet, check if the ray actually crossed it and didn't just go close by
            // Check if the ray crossed vertices
            int count3 = 0;
            for (int i=0; i<count2; i++) {
                if (v1X(possibleCROSSLIST(i))==gridCOx(loopX) && v1Y(possibleCROSSLIST(i))==gridCOy(loopY) || v2X(possibleCROSSLIST(i))==gridCOx(loopX) && v2Y(possibleCROSSLIST(i))==gridCOy(loopY) || v3X(possibleCROSSLIST(i))==gridCOx(loopX) && v3Y(possibleCROSSLIST(i))==gridCOy(loopY)) {
                    vertexCROSSLIST(count3) = possibleCROSSLIST(i);
                    count3++;
                }
            }
            
            // If the ray crossed a vertice determine if it needs correction
            if (count3 > 0) {
                float coN_max;
                float coN_min;
                for (int i=0; i<count3; i++) {
                    coN(i) = normal(vertexCROSSLIST(i),2);
                    if (coN(i)>coN_max) {
                        coN_max = coN(i);
                    }
                    if (coN(i)<coN_min) {
                        coN_min = coN(i);
                    }
                }
                if (coN_max<0 || coN_min>0) {
                    for (int i=0; i<count3; i++) {
                        facetCROSSLIST(facetCROSSLISTcounter(0)) = vertexCROSSLIST(i);
                        facetCROSSLISTcounter(0)++;
//                        Kokkos::atomic_add(&facetCROSSLISTcounter(0),1);
                    }
                } else {
                    correctionLIST(correctionLISTcounter(0),0) = loopX;
                    correctionLIST(correctionLISTcounter(0),1) = loopY;
                    correctionLISTcounter(0)++;
//                    Kokkos::atomic_add(&correctionLISTcounter(0),1);;
                }
            }
            
            // Check if the ray crossed facets
            int loopCHECKFACET;
            if (count2 > 0) {
                for (int i=0; i<count2; i++) {
                    loopCHECKFACET = possibleCROSSLIST(i);
                    
                    // Check if the ray crossed a facet
                    float Y1predicted = v2Y(loopCHECKFACET) - ((v2Y(loopCHECKFACET)-v3Y(loopCHECKFACET))*(v2X(loopCHECKFACET)-v1X(loopCHECKFACET))/(v2X(loopCHECKFACET)-v3X(loopCHECKFACET)));
                    float YRpredicted = v2Y(loopCHECKFACET) - ((v2Y(loopCHECKFACET)-v3Y(loopCHECKFACET))*(v2X(loopCHECKFACET)-gridCOx(loopX))/(v2X(loopCHECKFACET)-v3X(loopCHECKFACET)));
                    // Check if the ray is on the same side of the 2-3 edge as the 1st vertex
                    if (Y1predicted > v1Y(loopCHECKFACET) && YRpredicted > gridCOy(loopY) || Y1predicted < v1Y(loopCHECKFACET) && YRpredicted < gridCOy(loopY)) {
                        float Y2predicted = v3Y(loopCHECKFACET) - ((v3Y(loopCHECKFACET)-v1Y(loopCHECKFACET))*(v3X(loopCHECKFACET)-v2X(loopCHECKFACET))/(v3X(loopCHECKFACET)-v1X(loopCHECKFACET)));
                        YRpredicted = v3Y(loopCHECKFACET) - ((v3Y(loopCHECKFACET)-v1Y(loopCHECKFACET))*(v3X(loopCHECKFACET)-gridCOx(loopX))/(v3X(loopCHECKFACET)-v1X(loopCHECKFACET)));
                        // Check if the ray is on the same side of the 3-1 edge as the 2nd vertex
                        if (Y2predicted > v2Y(loopCHECKFACET) && YRpredicted > gridCOy(loopY) || Y2predicted < v2Y(loopCHECKFACET) && YRpredicted < gridCOy(loopY)) {
                            float Y3predicted = v1Y(loopCHECKFACET) - ((v1Y(loopCHECKFACET)-v2Y(loopCHECKFACET))*(v1X(loopCHECKFACET)-v3X(loopCHECKFACET))/(v1X(loopCHECKFACET)-v2X(loopCHECKFACET)));
                            YRpredicted = v1Y(loopCHECKFACET) - ((v1Y(loopCHECKFACET)-v2Y(loopCHECKFACET))*(v1X(loopCHECKFACET)-gridCOx(loopX))/(v1X(loopCHECKFACET)-v2X(loopCHECKFACET)));
                            // Check if the ray is on the same side of the 1-2 edge as the 3rd vertex
                            if (Y3predicted > v3Y(loopCHECKFACET) && YRpredicted > gridCOy(loopY) || Y3predicted < v3Y(loopCHECKFACET) && YRpredicted < gridCOy(loopY)) {
                                facetCROSSLIST(facetCROSSLISTcounter(0)) = loopCHECKFACET;
                                facetCROSSLISTcounter(0)++;
//                                Kokkos::atomic_add(&facetCROSSLISTcounter(0),1);
                            }
                        }
                    }
                }
                
                // Find the z-coordinate of the locations where the ray crosses each facet or vertex
                if (facetCROSSLISTcounter(0) > 0) {
                    int i = 0;
                    int loopFINDZ;
                    for (int j=0; j<facetCROSSLISTcounter(0) ; j++) {
                        loopFINDZ = facetCROSSLIST(j);
                        float planecoA = v1Y(loopFINDZ)*(v2Z(loopFINDZ)-v3Z(loopFINDZ))+v2Y(loopFINDZ)*(v3Z(loopFINDZ)-v1Z(loopFINDZ))+v3Y(loopFINDZ)*(v1Z(loopFINDZ)-v2Z(loopFINDZ));
                        float planecoB = v1Z(loopFINDZ)*(v2X(loopFINDZ)-v3X(loopFINDZ))+v2Z(loopFINDZ)*(v3X(loopFINDZ)-v1X(loopFINDZ))+v3Z(loopFINDZ)*(v1X(loopFINDZ)-v2X(loopFINDZ));
                        float planecoC = v1X(loopFINDZ)*(v2Y(loopFINDZ)-v3Y(loopFINDZ))+v2X(loopFINDZ)*(v3Y(loopFINDZ)-v1Y(loopFINDZ))+v3X(loopFINDZ)*(v1Y(loopFINDZ)-v2Y(loopFINDZ));
                        float planecoD = -v1X(loopFINDZ)*(v2Y(loopFINDZ)*v3Z(loopFINDZ)-v3Y(loopFINDZ)*v2Z(loopFINDZ))-v2X(loopFINDZ)*(v3Y(loopFINDZ)*v1Z(loopFINDZ)-v1Y(loopFINDZ)*v3Z(loopFINDZ))-v3X(loopFINDZ)*(v1Y(loopFINDZ)*v2Z(loopFINDZ)-v2Y(loopFINDZ)*v1Z(loopFINDZ));
                        if (sqrt(pow(planecoC,2)) < 1E-14) {
                            planecoC = 0;
                        }
                        gridCOzCROSS(i) = (-planecoD-planecoA*gridCOx(loopX)-planecoB*gridCOy(loopY))/planecoC;
                        i++;
                    }

                    // Keep only the unique values
                    int flagup;
                    int ucnt = 0;
                    for (int i=0; i<facetCROSSLISTcounter(0); i++) {
                        for (int j=0; j<facetCROSSLISTcounter(0); j++) {
                            if (gridCOzCROSS(i) == gridCOzCROSS(j) && i != j) {
                                flagup = 1;
                            }
                        }
                        if (flagup == 1) {
                            flagup = 0;
                        } else {
                            gridCOzCROSSunique(ucnt) = gridCOzCROSS(i);
                            ucnt++;
                        }
                    }
                    
                    // Label all the voxels that the ray passes through
                    if (ucnt%2 == 0) {
                        for (int loopASSIGN=0; loopASSIGN<ucnt/2; loopASSIGN++) {
                            float temp;
                            for (int j=0; j<ucnt; j++) {
                                for (int i=0; i<ucnt-1; i++) {
                                    if (gridCOzCROSSunique(i)>gridCOzCROSSunique(i+1)) {
                                        temp = gridCOzCROSSunique(i);
                                        gridCOzCROSSunique(i) = gridCOzCROSSunique(i+1);
                                        gridCOzCROSSunique(i+1) = temp;
                                    }
                                }
                            }
                            for (int i=0; i<elz; i++) {
                                voxelsINSIDE(i) = (gridCOz(i)>gridCOzCROSSunique(2*loopASSIGN) && gridCOz(i)<gridCOzCROSSunique(2*loopASSIGN+1));

                                // NOTE: IF THERE ARE EVER ANY ERRORS WITH THE CODE THIS IS LIKELY WHERE THEY ARE OCCURING (ROUNDING ERROR)
                                if (sqrt(pow(gridCOz(i)-gridCOzCROSSunique(2*loopASSIGN),2))/sqrt(pow(gridCOz(i),2)) < 1E-4) {
                                    voxelsINSIDE(i) = 1;
                                }
                                if (sqrt(pow(gridCOz(i)-gridCOzCROSSunique(2*loopASSIGN+1),2))/sqrt(pow(gridCOz(i),2)) < 1E-4) {
                                    voxelsINSIDE(i) = 1;
                                }
                                gridOUTPUT(loopX,loopY,i) += voxelsINSIDE(i);
                            }
                         }
                    } else if (ucnt != 0) {
                        if (count3 == 0) {
                            correctionLIST(correctionLISTcounter(0),0) = loopX;
                            correctionLIST(correctionLISTcounter(0),1) = loopY;
                            correctionLISTcounter(0)++;
//                            Kokkos::atomic_add(&correctionLISTcounter(0),1);
                        }
                    }
                }
            }
        }
    });
//    Kokkos::fence();

    // Use interpolation to fill in the rays that could not be voxelised
    if (correctionLISTcounter(0) > 0) {
        CArray<int> voxelsforcorrection(1);
        CArray<int> c1(1);
        CArray<int> c2(1);
        CArray<int> c3(1);
        CArray<int> c4(1);
        CArray<int> c5(1);
        CArray<int> c6(1);
        CArray<int> c7(1);
        CArray<int> c8(1);
        FOR_LOOP(loopC,0,correctionLISTcounter(0),
                k,0,elz,{
            c1(0) = -1;
            c2(0) = -1;
            c3(0) = -1;
            c4(0) = -1;
            c5(0) = -1;
            c6(0) = -1;
            c7(0) = -1;
            c8(0) = -1;
            if (correctionLIST(loopC,0)-1 == -1) {
                c1(0) = 0;
                c2(0) = 0;
                c3(0) = 0;
            }
            if (correctionLIST(loopC,0)+1 == elx) {
                c6(0) = 0;
                c7(0) = 0;
                c8(0) = 0;
            }
            if (correctionLIST(loopC,1)-1 == -1) {
                c1(0) = 0;
                c4(0) = 0;
                c6(0) = 0;
            }
            if (correctionLIST(loopC,1)+1 == ely) {
                c3(0) = 0;
                c5(0) = 0;
                c8(0) = 0;
            }
            if (c1(0) == -1) {
                c1(0) = gridOUTPUT(correctionLIST(loopC,0)-1,correctionLIST(loopC,1)-1,k);
            }
            if (c2(0) == -1) {
                c2(0) = gridOUTPUT(correctionLIST(loopC,0)-1,correctionLIST(loopC,1),k);
            }
            if (c3(0) == -1) {
                c3(0) = gridOUTPUT(correctionLIST(loopC,0)-1,correctionLIST(loopC,1)+1,k);
            }
            if (c4(0) == -1) {
                c4(0) = gridOUTPUT(correctionLIST(loopC,0),correctionLIST(loopC,1)-1,k);
            }
            if (c5(0) == -1) {
                c5(0) = gridOUTPUT(correctionLIST(loopC,0),correctionLIST(loopC,1)+1,k);
            }
            if (c6(0) == -1) {
                c6(0) = gridOUTPUT(correctionLIST(loopC,0)+1,correctionLIST(loopC,1)-1,k);
            }
            if (c7(0) == -1) {
                c7(0) = gridOUTPUT(correctionLIST(loopC,0)+1,correctionLIST(loopC,1),k);
            }
            if (c8(0) == -1) {
                c8(0) = gridOUTPUT(correctionLIST(loopC,0)+1,correctionLIST(loopC,1)+1,k);
            }
            voxelsforcorrection(0) = c1(0)+c2(0)+c3(0)+c4(0)+c5(0)+c6(0)+c7(0)+c8(0);
            if (voxelsforcorrection(0) >= 4) {
                gridOUTPUT(correctionLIST(loopC,0),correctionLIST(loopC,1),k) = 1;
            }
        });
    }
}





