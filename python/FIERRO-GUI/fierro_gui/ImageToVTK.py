# ************************ Read in Image Files ************************

#imageDir = 'Path/to/your/images'
#outDir = 'Path/to/your/vtk'
#out = 'Name of your output file'

def ImageToVTK(imageDir, out, outDir, file_type):
    import PIL.Image
    import numpy as np
    import glob, os, csv
    os.chdir(imageDir)

#    # Convert any jpg or jpeg files to png
#    for file in glob.glob("*.jpg") + glob.glob("*.jpeg"):
##        print ("File: " + file)
#        name, file_extension = os.path.splitext(file)
#        if file_extension != ".png":
##            print ("Name: " + name + " file_extension: " + file_extension)
#            im = PIL.Image.open(file)
#            im.save(name + ".png")

    # Get the list of files
    all_files = os.listdir(imageDir)
    fileNames = sorted([os.path.join(imageDir, f) for f in all_files if f.endswith(file_type)])
        
    NumFiles = len(fileNames)
    img = PIL.Image.open(fileNames[0])
    rows, columns = img.size
    imageStack = np.empty((NumFiles, columns, rows), dtype=int)
    for i in range(NumFiles):
        img = PIL.Image.open(fileNames[i])
        if img.mode == '1':
            img_array = np.array(img)
            binary_array = (img_array > 0).astype(np.uint8)
            binary_array = np.rot90(binary_array, k=-1)
            binary_array = binary_array.transpose()
        elif img.mode == 'L':
            img_array = np.array(img)
            threshold = 128
            binary_array = (img_array >= threshold).astype(np.uint8)
        elif img.mode == 'RGB':
            img_array = np.array(img)
            grayscale = np.dot(img_array[...,:3],[0.2989, 0.5870, 0.1140])
            threshold = 128
            binary_array = (grayscale >= threshold).astype(np.uint8)
            binary_array = np.rot90(binary_array, k=-1)
            binary_array = binary_array.transpose()
        imageStack[i, :, :] = binary_array

    # *************** Convert Image Stack to VTK **********************

    img_arr=imageStack.flatten()

    # VTK Header Information
    VTK_FileVersion = '# vtk DataFile Version 3.0'
    SecondLine = 'Cycle in 3D Grid'
    ThirdLine = 'ASCII'
    FourthLine = 'DATASET STRUCTURED_POINTS'
    FifthLine = 'DIMENSIONS'+' '+str(rows)+' '+str(columns)+' '+str(NumFiles)
    SixthLine = 'ORIGIN 0 0 0'
    SeventhLine = 'SPACING 1 1 1'
    EigthLine = 'POINT_DATA'+' '+str(rows*columns*NumFiles)
    NinthLine = 'SCALARS density float'
    TenthLine = 'LOOKUP_TABLE default'

    # Write Header and Data
    print ("Outputting VTK file in " + outDir + '/' + out + '.vtk')
    with open(outDir+'/'+out+'.vtk','w') as VTK_File:
        VTK_File.write(VTK_FileVersion+'\n')
        VTK_File.write(SecondLine+'\n')
        VTK_File.write(ThirdLine+'\n')
        VTK_File.write(FourthLine+'\n')
        VTK_File.write(FifthLine+'\n')
        VTK_File.write(SixthLine+'\n')
        VTK_File.write(SeventhLine+'\n')
        VTK_File.write(EigthLine+'\n')
        VTK_File.write(NinthLine+'\n')
        VTK_File.write(TenthLine+'\n')

        for i in range(rows*columns*NumFiles):
            csv.writer(VTK_File).writerow([img_arr[i]])

    VTK_File.close()
