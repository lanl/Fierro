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

#include <matar.h>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <stl-to-voxelvtk.h>

using namespace mtr; // matar namespace

// Functions used within MAIN
std::tuple<CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, unsigned int> binary_stl_reader(const char* stl_file_path); // BINARY STL READER FUNCTION
void main_function(CArray<bool> &gridOUTPUT, int &gridX, int &gridY, int &gridZ, CArray<float> &normalX, CArray<float> &normalY, CArray<float> &normalZ, CArray<float> &v1X, CArray<float> &v1Y, CArray<float> &v1Z, CArray<float> &v2X, CArray<float> &v2Y, CArray<float> &v2Z, CArray<float> &v3X, CArray<float> &v3Y, CArray<float> &v3Z, unsigned int &n_facets, CArray<float> &nodex, CArray<float> &nodey, CArray<float> &nodez); // VOXELIZATION FUNCTION

double sign (double x1, double y1, double v1x, double v1y, double v2x, double v2y) {
    return (x1 - v2x) * (v1y - v2y) - (v1x - v2x) * (y1 - v2y);
}

void compute_normal(float* normal, float* a, float* b, float* c) {
    float ba[3];
    float ca[3];
    for (size_t i = 0; i < 3; i++) {
        ba[i] = b[i] - a[i];
        ca[i] = c[i] - a[i];
    }
    normal[0] =  (ba[1] * ca[2] - ba[2] * ca[1]);
    normal[1] = -(ba[0] * ca[2] - ba[2] * ca[0]);
    normal[2] =  (ba[0] * ca[1] - ba[1] * ca[0]);
}

void Voxelizer::create_voxel_vtk(
        const char* stl_file_path, const char* vtk_file_path, 
        int gridX, int gridY, int gridZ,
        bool use_index_space) {
    
    // Start Clock
    auto start = std::chrono::high_resolution_clock::now();

    // Read the stl file
    auto [normalX, normalY, normalZ, v1X, v1Y, v1Z, v2X, v2Y, v2Z, v3X, v3Y, v3Z, n_facets] = binary_stl_reader(stl_file_path);
    
    // Voxelize the part in the x,y,z directions
    CArray<float> nodex(gridX+2);
    CArray<float> nodey(gridY+2);
    CArray<float> nodez(gridZ+2);
    
    // X-direction voxelization
    CArray<bool> gridOUTPUTX(gridY+2,gridZ+2,gridX+2);
    FOR_ALL (i, 0, gridY+2,
                j, 0, gridZ+2,
                k, 0, gridX+2,{
        gridOUTPUTX(i,j,k) = 0;
    });
    main_function(gridOUTPUTX, gridY, gridZ, gridX, normalY, normalZ, normalX, v1Y, v1Z, v1X, v2Y, v2Z, v2X, v3Y, v3Z, v3X, n_facets, nodey, nodez, nodex);

    // Y-direction voxelization
    CArray<bool> gridOUTPUTY(gridZ+2,gridX+2,gridY+2);
    FOR_ALL (i, 0, gridZ+2,
                j, 0, gridX+2,
                k, 0, gridY+2,{
        gridOUTPUTY(i,j,k) = 0;
    });
    main_function(gridOUTPUTY, gridZ, gridX, gridY, normalZ, normalX, normalY, v1Z, v1X, v1Y, v2Z, v2X, v2Y, v3Z, v3X, v3Y, n_facets, nodez, nodex, nodey);
    
    // Z-direction voxelization
    CArray<bool> gridOUTPUTZ(gridX+2,gridY+2,gridZ+2);
    FOR_ALL (i, 0, gridX+2,
                j, 0, gridY+2,
                k, 0, gridZ+2,{
        gridOUTPUTZ(i,j,k) = 0;
    });
    main_function(gridOUTPUTZ, gridX, gridY, gridZ, normalX, normalY, normalZ, v1X, v1Y, v1Z, v2X, v2Y, v2Z, v3X, v3Y, v3Z, n_facets, nodex, nodey, nodez);

    // Sum the voxelization in the XYZ-directions
    CArray<int> TOTgrid(gridX+2,gridY+2,gridZ+2);
    FOR_ALL(i,0,gridX+2,
            j,0,gridY+2,
            k,0,gridZ+2,{
        TOTgrid(i,j,k) = gridOUTPUTX(j,k,i) + gridOUTPUTY(k,i,j) + gridOUTPUTZ(i,j,k);
    });
    CArray<bool> OUTPUTgrid(gridZ,gridY,gridX);
    FOR_ALL(i,1,gridX+1,
            j,1,gridY+1,
            k,1,gridZ+1,{
        OUTPUTgrid(k-1,j-1,i-1) = TOTgrid(i,j,k) > 0;
    });

    // VTK WRITER
    int i,j,k;
    auto out=fopen(vtk_file_path,"w"); // open the file

    fprintf(out,"# vtk DataFile Version 3.0\n"); // write the header
    fprintf(out,"Header\n");
    fprintf(out,"ASCII\n");
    fprintf(out,"DATASET RECTILINEAR_GRID\n");
    fprintf(out,"DIMENSIONS %d %d %d\n",gridX+1,gridY+1,gridZ+1);
    
    fprintf(out,"X_COORDINATES %d float\n", gridX+1); // nodal x coordinates
    for (i=0; i<(gridX+1); i++) {
        fprintf(out,"%f ", use_index_space ? i : nodex(i));
    }
    fprintf(out,"\nY_COORDINATES %d float\n", gridY+1); // nodal y coordinates
    for (i=0; i<(gridY+1); i++) {
        fprintf(out,"%f ", use_index_space ? i : nodey(i));
    }
    fprintf(out,"\nZ_COORDINATES %d float\n", gridZ+1); // nodal z coordinates
    for (i=0; i<(gridZ+1); i++) {
        fprintf(out,"%f ", use_index_space ? i : nodez(i));
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

    printf("\nfinished\n\n");
}

// ==============================================================
// ------------------- VOXELIZATION FUNCTIONS -------------------
// ==============================================================

// BINARY STL READER - (Note: it can ONLY read binary stl files)
std::tuple<CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, unsigned int> binary_stl_reader(const char* stl_file_path){
    // Open the binary file
    const char* filename = stl_file_path;
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
    CArray<float> normalX(n_facets);
    CArray<float> normalY(n_facets);
    CArray<float> normalZ(n_facets);
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
        if (std::sqrt(normalp[0] * normalp[0] + normalp[0] * normalp[0] + normalp[0] * normalp[0]) <= (1-1e-6))
            compute_normal(normalp, v1p, v2p, v3p);
        normalX(i) = normalp[0];
        normalY(i) = normalp[1];
        normalZ(i) = normalp[2];
        v1X(i) = v1p[0];
        v1Y(i) = v1p[1];
        v1Z(i) = v1p[2];
        v2X(i) = v2p[0];
        v2Y(i) = v2p[1];
        v2Z(i) = v2p[2];
        v3X(i) = v3p[0];
        v3Y(i) = v3p[1];
        v3Z(i) = v3p[2];
    }
    input.close();
    return {normalX, normalY, normalZ, v1X, v1Y, v1Z, v2X, v2Y, v2Z, v3X, v3Y, v3Z, n_facets};
}

// VOXELIZATION FUNCTION
void main_function(CArray<bool> &gridOUTPUT, int &gridX, int &gridY, int &gridZ, CArray<float> &normalX, CArray<float> &normalY, CArray<float> &normalZ, CArray<float> &v1X, CArray<float> &v1Y, CArray<float> &v1Z, CArray<float> &v2X, CArray<float> &v2Y, CArray<float> &v2Z, CArray<float> &v3X, CArray<float> &v3Y, CArray<float> &v3Z, unsigned int &n_facets, CArray<float> &nodex, CArray<float> &nodey, CArray<float> &nodez){
    // Find the global maximum and minimum values of the mesh
    float meshXmax;
    float meshXmin;
    float meshYmax;
    float meshYmin;
    float meshZmax;
    float meshZmin;
    
    // Global maximum x-direction
    REDUCE_MAX(i,0,n_facets,meshXmax, {
        if (v1X(i) > meshXmax || v2X(i) > meshXmax || v3X(i) > meshXmax) {
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
    REDUCE_MIN(i,0,n_facets,meshXmin, {
        if (v1X(i) < meshXmin || v2X(i) < meshXmin || v3X(i) < meshXmin) {
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
    REDUCE_MAX(i,0,n_facets,meshYmax, {
        if (v1Y(i) > meshYmax || v2Y(i) > meshYmax || v3Y(i) > meshYmax) {
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
    REDUCE_MIN(i,0,n_facets,meshYmin, {
        if (v1Y(i) < meshYmin || v2Y(i) < meshYmin || v3Y(i) < meshYmin) {
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
    REDUCE_MAX(i,0,n_facets,meshZmax, {
        if (v1Z(i) > meshZmax || v2Z(i) > meshZmax || v3Z(i) > meshZmax) {
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
    REDUCE_MIN(i,0,n_facets,meshZmin, {
        if (v1Z(i) < meshZmin || v2Z(i) < meshZmin || v3Z(i) < meshZmin) {
            if (v1Z(i) < v2Z(i) && v1Z(i) < v3Z(i)) {
                meshZmin = v1Z(i);
            } else if (v2Z(i) < v1Z(i) && v2Z(i) < v3Z(i)) {
                meshZmin = v2Z(i);
            } else if (v3Z(i) < v1Z(i) && v3Z(i) < v2Z(i)) {
                meshZmin = v3Z(i);
            }
        }
    }, meshZmin);

    // Define the grid that the voxel mesh will be generated on
    float voxwidth;
    
    // Voxel grid x-direction
    voxwidth = (meshXmax-meshXmin)/(gridX);
    int elx = gridX+2;
    CArray<float> gridCOx(elx);
    FOR_ALL (i, 0, elx, {
        gridCOx(i) = (meshXmin+(voxwidth/2)) + voxwidth*(i-1);
        nodex(i) = gridCOx(i) + voxwidth / 2;
    });

    // Voxel grid y-direction
    voxwidth = (meshYmax-meshYmin)/(gridY);
    int ely = gridY+2;
    CArray<float> gridCOy(ely);
    FOR_ALL(i, 0, ely, {
        gridCOy(i) = (meshYmin+(voxwidth/2)) + voxwidth*(i-1);
        nodey(i) = gridCOy(i) + voxwidth / 2;
    });

    // Voxel grid z-direction
    voxwidth = (meshZmax-meshZmin)/(gridZ);
    int elz = gridZ+2;
    CArray<float> gridCOz(elz);
    FOR_ALL(i, 0, elz, {
        gridCOz(i) = (meshZmin+(voxwidth/2)) + voxwidth*(i-1);
        nodez(i) = gridCOz(i) + voxwidth / 2;
    });
    double voxwidth_z = voxwidth;

    // Z-DIRECTIONAL RAY EXECUTION
    // Get the minimum and maximum x,y,z coordinates for each facet
    CArray<float> facetXmin(n_facets);
    CArray<float> facetXmax(n_facets);
    CArray<float> facetYmin(n_facets);
    CArray<float> facetYmax(n_facets);
    CArray<float> facetZmin(n_facets);
    CArray<float> facetZmax(n_facets);
    
    FOR_ALL (i,0,n_facets, {
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
    
    // Create a list to record all rays that fail to voxelize
    CArray<int> correctionLIST(n_facets,2);
    size_t facetCROSSLISTcounter = 0;
    size_t correctionLISTcounter = 0;
    
    // Find which facets could be crossed by the ray in the y-direction
    CArray<int> possibleCROSSLISTy(n_facets,ely);
    CArray<int> counter(1);
    CArray<int> count1(ely);
    counter(0) = 0;
    FOR_ALL (loopY,0,ely,
             i,0,n_facets,{
        if (facetYmin(i) <= gridCOy(loopY) && facetYmax(i) >= gridCOy(loopY)) {
            possibleCROSSLISTy(counter(0),loopY) = i;
            counter(0)++;
        }
        if (i == n_facets-1) {
            count1(loopY) = counter(0);
            counter(0) = 0;
        }
    });
    

    std::vector<size_t> cross_list;
    std::vector<double> zs;

    for(size_t loopY = 0; loopY < ely; loopY++) {
        double y = gridCOy(loopY);
        
        for(size_t loopX = 0; loopX < elx; loopX++) {
            double x = gridCOx(loopX);

            cross_list.clear();
            zs.clear();
            
            for (size_t i = 0; i < n_facets; i++) {
                if ((facetYmin(i) <= y && y <= facetYmax(i))
                    && (facetXmin(i) <= x && x <= facetXmax(i))) {
                    
                    float d1, d2, d3;
                    bool has_neg, has_pos;
                    float fuzz = 1e-6;

                    d1 = sign(x, y, v1X(i), v1Y(i), v2X(i), v2Y(i));
                    d2 = sign(x, y, v2X(i), v2Y(i), v3X(i), v3Y(i));
                    d3 = sign(x, y, v3X(i), v3Y(i), v1X(i), v1Y(i));

                    has_neg = (d1 < -fuzz) || (d2 < -fuzz) || (d3 < -fuzz);
                    has_pos = (d1 > fuzz) || (d2 > fuzz) || (d3 > fuzz);
                    
                    if (!(has_neg && has_pos)) {
                        cross_list.push_back(i);
                    }
                }
            }

            for (size_t i : cross_list) {
                zs.push_back(-(
                      (x - v1X(i)) * normalX(i)
                    + (y - v1Y(i)) * normalY(i)
                    + (0 - v1Z(i)) * normalZ(i)
                ) / normalZ(i));
            }

            std::sort(zs.begin(), zs.end());

            if (cross_list.size() == 0)
                continue;
            
            std::vector<double> unique_zs;
            for (size_t i = 0; i < zs.size(); i++) 
                if ((i == zs.size() - 1) || (zs[i+1] - zs[i] >= 1e-5))
                    unique_zs.push_back(zs[i]);
            
            for (size_t i = 0; i < unique_zs.size() - 1; i++) {
                if (i % 2 == 1)
                    continue;
                size_t z1 = (unique_zs[i] - meshZmin) / voxwidth_z + 1;
                size_t z2 = (unique_zs[i + 1] - meshZmin) / voxwidth_z + 1;

                for (size_t loopZ = z1; loopZ <= z2; loopZ++) {
                    gridOUTPUT(loopX, loopY, loopZ) = 1;
                }
            }
        }
    }
}




