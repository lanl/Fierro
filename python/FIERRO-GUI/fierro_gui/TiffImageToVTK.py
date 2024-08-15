    
def TiffImageToVTK(folder_dir, out, out_dir, file_format):
    # import the modules

    # opperating system
    import os
    from os import listdir
    from glob import glob

    #import open3d as o3d

    # image
    from PIL import Image
    #from libtiff import TIFF
    import cv2 # type: ignore

    # math
    import numpy as np
    from numpy import asarray

    # from multiprocessing import Pool
    from multiprocessing.pool import ThreadPool as Pool

    #c writer
    import c_writer # type: ignore

    # voxel2stl
    import voxel2stl # type: ignore



    # -------------------------------------------------------
    # This gives the index value of the point or the elem
    # the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
    # the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
    #--------------------------------------------------------
    #
    # Returns a global id for a given i,j,k
    def get_id(i, j, k, num_i, num_j):
        return (i + j*num_i + k*num_i*num_j)


    #  coding to align image
    #def get_alignment_matrix(box: o3d.geometry.OrientedBoundingBox) -> np.ndarray:
    #    box_points = np.asarray(box.get_box_points())
    #    axes = [
    #        box_points[1] - box_points[0],
    #        box_points[2] - box_points[0],
    #        box_points[3] - box_points[0]
    #    ]
    #    axes = [a / np.linalg.norm(a) for a in axes]
    #    return np.row_stack(axes)
    #
    #def get_alignment_matrix_from_points(points: np.ndarray) -> np.ndarray:
    #    cloud = o3d.geometry.PointCloud()
    #    cloud.points = o3d.utility.Vector3dVector(points)
    #    cloud.paint_uniform_color((1, 0, 0))
    #    return get_alignment_matrix(cloud.get_minimal_oriented_bounding_box())
    #
    # Code for aligning a point cloud ^








    # get the path/directory
    # folder_dir = "./data-alot-tiff"

    # count the number of files and save file name and index
    # args[0][0] = id of first image
    # args[0][1] = file name including path
    #args = list(enumerate(glob(f'{folder_dir}/*.jpg')))
    args = list(enumerate(glob(f'{folder_dir}/*.tif')))
    num_files = len(args)



    # determine the dimensions of the image 
    img= Image.open(args[0][1])
    NxImage = img.size[0]
    NyImage = img.size[1]
    NzImage = num_files



    # create an array to store voxel values
    voxelImageValues = np.zeros([NxImage, NyImage, NzImage])
    fileread = np.zeros(NzImage)


    # define worker function before a Pool is instantiated
    def read_data(file_id, file_name):
        
        print("reading file " + str(file_id+1) + " of " + str(num_files))

        # PIL images into NumPy arrays
        img= Image.open(file_name)
            
            
        img_t = img.transpose(Image.ROTATE_270) # numpy reads the transpose of the image
        img_t = img_t.rotate(12, Image.NEAREST)
            
        
            
        # Normalise to range 0..255
        #norm = (img_t.astype(np.float)-img_t.min())*255.0 / (img_t.max()-img_t.min())
        arr = asarray(img_t)
        voxelImageValues[0:NxImage, 0:NyImage, file_id] = ((arr - arr.min()) * (1/(arr.max() - arr.min()) * 255)).astype('uint8')
        #voxelImageValues /= 2712
            
        # PIL package randomonly doesn't open files, so recording that it did read the file
        fileread[file_id] = 1


    # end function




    pool_size = 6  # your "parallelness"
    pool = Pool(pool_size)


    # save the data from the files
    for file_id, file_name in args:

        pool.apply_async(read_data, (file_id, file_name))

    # end for

    pool.close()
    pool.join()

    #checking to see if PIL failed to read image file
    print(" ")
    for file_id, file_name in args:
        if(fileread[file_id]==0):
        
            print("PIL python packaged failure, re-reading image file " + str(file_id+1))
            
            # will try 10 times to read the image file
            for attempt in range(10):
                print("attempt " + str(attempt+1))
                
                # try reading the data again
                read_data(file_id, file_name)
                
                # exit when I read it
                if (fileread[file_id]==1):
                    print("sucessfully read image file")
                    print(" ")
                    break
                # end if
                
                if (attempt == 49):
                    print("ERROR: image file could not be read by the PIL package")
                    exit(1)
                # end if
                
            #end for attempt
            
        # end if on proper reading of data
    # end for file_id



    # coarsen resolution
    Nx = 256
    Ny = 256
    Nz = 256

    # contour value
    contour_level = 127  # half way from 0:255


    Nx = min(Nx, NxImage)
    Ny = min(Ny, NyImage)
    Nz = min(Nz, NzImage)

    reduce_resX = round(NxImage/Nx)
    reduce_resY = round(NyImage/Ny)
    reduce_resZ = round(NzImage/Nz)

    print("reducing res by = " + str(reduce_resX) + " x " + str(reduce_resY) + " x " + str(reduce_resZ))


    voxelValues = np.zeros([Nx,Ny,Nz])


    # a function to reduce resolution
    def reduce_voxel_res(k):
        print("making new layer = " + str(k) + " of " + str(Nz))
        for j in range(Ny):
            for i in range(Nx):
        
                kstart = k*reduce_resZ
                kstop = min(NzImage,kstart+reduce_resZ)
                
                jstart = j*reduce_resY
                jstop = min(NyImage,jstart+reduce_resY)
                
                istart = i*reduce_resX
                istop = min(NxImage,istart+reduce_resX)
                
                #for kcount in range(kstart,kstop):
                #    for jcount in range(jstart, jstop):
                #        for icount in range(istart, istop):
                #            voxelValues[i][j][k] += voxelImageValues[icount][jcount][kcount]
                #end for
                
                # average
                #voxelValues[i][j][k] = np.sum(voxelImageValues[istart:istop, jstart:jstop, kstart:kstop])
                #voxelValues[i][j][k] /= reduce_resX*reduce_resY*reduce_resZ
                
                # min, max, or mean
                voxelValues[i][j][k] = np.mean(voxelImageValues[istart:istop, jstart:jstop, kstart:kstop])
                
            # end for
        # end for

    # end function



    poolRes = Pool(pool_size)
    for k in range(Nz):
        poolRes.apply_async(reduce_voxel_res, (k,))
    # end for

    poolRes.close()
    poolRes.join()

    #voxelValues = voxelImageValues

    print("\n")
    print("writing vtk file")


    num_elems  = Nx*Ny*Nz
    num_points = (Nx+1)*(Ny+1)*(Nz+1)
    Lx = 1
    Ly = 1
    Lz = 1
    Dx = Lx/Nx
    Dy = Ly/Ny
    Dz = Lz/Nz


    if (file_format == ".vtk"):
        # Rectilinear VTK format
        # fname = "VoxelValues.vtk"
        file1 = open(out_dir+'/'+out+'.vtk','w')
        file1.write("# vtk DataFile Version 3.0\n")
        file1.write("Mesh for Fierro\n")
        file1.write("ASCII \n")
        file1.write("DATASET RECTILINEAR_GRID\n\n")
        file1.write("DIMENSIONS " + str(Nx+1) + " " + str(Ny+1) + " " + str(Nz+1) + "\n\n")

        file1.write("X_COORDINATES " + str(Nx+1) + " float\n")
        for i in range(0, Nx+1):
            x_pt = i*Dx
            file1.write(str(x_pt) + " ")
        # end for
        file1.write("\n")
        file1.write("\n")


        file1.write("Y_COORDINATES " + str(Ny+1) + " float\n")
        for j in range(0, Ny+1):
            y_pt = j*Dy
            file1.write(str(y_pt) + " ")
        # end for
        file1.write("\n")
        file1.write("\n")

        file1.write("Z_COORDINATES " + str(Nz+1) + " float\n")
        for k in range(0, Nz+1):
            z_pt = k*Dz
            file1.write(str(z_pt) + " ")
        # end for
        file1.write("\n")
        file1.write("\n")

        # Write the scalar elem variable to file
        file1.write("CELL_DATA " + str(num_elems) + "\n")
        file1.write("SCALARS density float 1 \n")  # the 1 is number of scalar components [1:4]
        file1.write("LOOKUP_TABLE default \n")
        file1.flush()


        # flatten data to 1D following Fortran convention
        print("3D voxel resolution = " + str(Nx) + " x " + str(Ny) + " x " + str(Nz))
        print("total number of voxels = ", voxelValues.size)


        # Closing file
        file1.close()

        # write out values  (95 on grey scale)
        c_writer.write((voxelValues >= contour_level).flatten('F'), out_dir+'/'+out+'.vtk')
        # contour level was 127, half of 255 is 127.5

    # creating and writing STL file
    elif (file_format == ".stl"):
        print(" ")

        xvals = np.zeros(Nx)
        for i in range(0, Nx):
            xvals[i] = i*Dx + Dx*0.5

        yvals = np.zeros(Ny)
        for j in range(0, Ny):
            yvals[j] = j*Dy + Dy*0.5

        # adding padding
        zvals = np.zeros(Nz+2)
        for k in range(0, Nz+2):
            zvals[k] = k*Dz - Dz*0.5  # includes padding


        PattedVoxelValues = np.zeros([Nx, Ny, Nz+2])
        PattedVoxelValues[0:Nx,0:Ny,1:Nz+1] = voxelValues

        # write a STL file
        voxel2stl.write((PattedVoxelValues >= contour_level).flatten('F'), xvals, yvals, zvals, out, out_dir)
        # contour level was 127, half of 255 is 127.5


    # Unstructured VTK format
    '''
    # create a file
    file1 = open("VoxelValues.vtk", "w")
    file1.write("# vtk DataFile Version 2.0\n")
    file1.write("Mesh for Fierro\n")
    file1.write("ASCII \n")
    file1.write("DATASET UNSTRUCTURED_GRID\n\n")

    num_elems  = Nx*Ny*Nz
    num_points = (Nx+1)*(Ny+1)*(Nz+1)
    Lx = 1
    Ly = 1
    Lz = 1
    Dx = Lx/Nx
    Dy = Ly/Ny
    Dz = Lz/Nz

    file1.write("POINTS " + str(num_points) + " float\n")

    # write points
    for k in range(0, Nz+1):
        for j in range(0, Ny+1):
            for i in range(0, Nx+1):
                x_pt = i*Dx
                y_pt = j*Dy
                z_pt = k*Dz
                file1.write(str(x_pt) + " " + str(y_pt) + " " + str(z_pt) + "\n")
    # end for


    # write elem values
    file1.write("\n")
    file1.write("CELLS " + str(num_elems) + " " + str(num_elems+num_elems*8) + "\n")  # size=all printed values


    # write all global point numbers for this elem
    for k in range(0, Nz):
        for j in range(0, Ny):
            for i in range(0, Nx):
            
                file1.write("8 ")  # // num points in this elem
                for kcount in range(k, k+2):
                    for jcount in range(j, j+2):
                        for icount in range(i, i+2):
                            
                            # get the gobal id
                            gid = icount + jcount*(Nx+1) + kcount*(Nx+1)*(Ny+1)#get_id(icount,jcount,kcount,Nx,Ny)
                            
                            print(str(icount) + " " + str(jcount) + " " + str(kcount) + " ")
                            print(gid)
                            
                            file1.write(str(gid) + " ")
                
                # done exporting nodes in the this elem
                print(" \n")
                file1.write("\n")
    # end for



    file1.write("\n")
    file1.write("CELL_TYPES " + str(num_elems) + "\n")
    # elem types:
    # FE convention linear hex = 12, linear quad = 9
    # element types: https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    # element types: https://kitware.github.io/vtk-js/api/Common_DataModel_CellTypes.html
    # vtk format: https://www.kitware.com//modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
    for k in range(0, Nz):
        for j in range(0, Ny):
            for i in range(0, Nx):
            
                file1.write("11 \n") # linear i,j,k hex is type 11
        
    # end for


    # Write the scalar elem variable to file
    file1.write("\n")
    file1.write("CELL_DATA " + str(num_elems) + "\n")
    file1.write("SCALARS elem_var float 1\n")  # the 1 is number of scalar components [1:4]
    file1.write("LOOKUP_TABLE default\n")
    for k in range(0, Nz):
        for j in range(0, Ny):
            for i in range(0, Nx):
            
                file1.write( str(testValues[i][j][k]) + "\n" );
            
    # end for



    # Closing file
    file1.close()
    '''
