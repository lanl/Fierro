# ************************ Read in Image Files ************************

#imageDir = 'Path/to/your/images'
#outDir = 'Path/to/your/vtk'
#out = 'Name of your output file'

def ImageToVTK(imageDir, out, outDir):
    import PIL.Image
    import numpy as np
    import glob, os, csv
    os.chdir(imageDir)

    # Convert any jpg or jpeg files to png
    for file in glob.glob("*.jpg") + glob.glob("*.jpeg"):
        print ("File: " + file)
        name, file_extension = os.path.splitext(file)
        if file_extension != ".png":
            print ("Name: " + name + " file_extension: " + file_extension)
            im = PIL.Image.open(file)
            im.save(name + ".png")

    fileNames = []
    for file in glob.glob("*.png"):
        fileNames.append(file)

    NumFiles = len(fileNames)

    img = PIL.Image.open(fileNames[0])
    rows,columns = img.size
    imageStack = np.empty((NumFiles,columns,rows),dtype=int)

    for i in range(NumFiles):
        img = PIL.Image.open(fileNames[i])
        # pattern matching to add channel variable manually to any file that was 
        # converted from jpg (which are read in as 3 channel images)
        match img:
            case [a, b, c]:
                if c is None:
                    imageStack[i,:,3] = img
                else:
                    imageStack[i,:,:] = img
        
    imageStack = imageStack[1:(np.shape(imageStack)[0]-1),:,:]
    NumFiles = np.shape(imageStack)[0]

    # *************** Convert Image Stack to VTK **********************

    img_arr=imageStack.reshape(rows*columns*NumFiles,1)

    Length = len(img_arr)
    n = int(Length/1) # Number of values in chunked vtk file. Change denominator
    # if you want to write chunks of the data out at a time. Helpful for big files.

    if (Length/n)-round(Length/n) > 0:
        ChunkNumber = round(Length/n)+1
        
    elif (Length/n)-round(Length/n) <= 0:
        ChunkNumber = round(Length/n)    

    def splitData(l,n):
        
        for i in range(0,Length,n):
            yield img_arr[i:i+n]
        
    chunkedData = list(splitData(img_arr,n))

    # VTK Header Information
    VTK_FileVersion = '# vtk DataFile Version 3.0'
    SecondLine = 'Cycle in 3D Grid'
    ThirdLine = 'ASCII'
    FourthLine = 'DATASET STRUCTURED_POINTS'
    FifthLine = 'DIMENSIONS'+' '+str(rows)+' '+str(columns)+' '+str(NumFiles)
    SixthLine = 'ORIGIN 0 0 0'
    SeventhLine = 'SPACING 1 1 1'
    EigthLine = 'POINT_DATA'+' '+str(rows*columns*NumFiles)
    NinthLine = 'SCALARS Cycle float'
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

        for i in range(ChunkNumber):
            
            csv.writer(VTK_File,delimiter='\n').writerows(np.reshape(chunkedData[i],(-1,1)))

    VTK_File.close()
