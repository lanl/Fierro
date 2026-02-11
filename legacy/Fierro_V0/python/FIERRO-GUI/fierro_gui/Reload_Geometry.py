from paraview.simple import *
from LAFFT2VTK import *
from PySide6.QtWidgets import QTableWidgetItem

def Upload_Batch_Geometry(self):
    # Import Geometry [.stl]
    if ".stl" in self.file_type:
        # Show the stl part
        self.stl = paraview.simple.STLReader(FileNames = self.in_file_path)
        self.display = paraview.simple.Show(self.stl, self.render_view)
        paraview.simple.ResetCamera()
        
        # Define number of voxels based on previous inputs
        self.INNumberOfVoxelsX.setText(self.TParts.item(0,7).text())
        self.INNumberOfVoxelsY.setText(self.TParts.item(0,8).text())
        self.INNumberOfVoxelsZ.setText(self.TParts.item(0,9).text())
        
        # Define geometry size based on previous inputs
        if self.BStlDimensions.isChecked():
            # Get the bounds of the stl
            self.stl.UpdatePipeline()
            self.bounds = self.stl.GetDataInformation().GetBounds()
            self.stlLx = round(self.bounds[1] - self.bounds[0],4)
            self.stlLy = round(self.bounds[3] - self.bounds[2],4)
            self.stlLz = round(self.bounds[5] - self.bounds[4],4)
            self.stlOx = self.bounds[0]
            self.stlOy = self.bounds[2]
            self.stlOz = self.bounds[4]
            
            # Open voxelization tab and fill out data
            self.INOriginX.setText(str(self.stlOx))
            self.INOriginY.setText(str(self.stlOy))
            self.INOriginZ.setText(str(self.stlOz))
            self.INLengthX.setText(str(self.stlLx))
            self.INLengthY.setText(str(self.stlLy))
            self.INLengthZ.setText(str(self.stlLz))
        elif self.BCustomDimensions.isChecked():
            # Open voxelization tab and fill out data
            self.INOriginX.setText(self.TParts.item(0,1).text())
            self.INOriginY.setText(self.TParts.item(0,2).text())
            self.INOriginZ.setText(self.TParts.item(0,3).text())
            self.INLengthX.setText(self.TParts.item(0,4).text())
            self.INLengthY.setText(self.TParts.item(0,5).text())
            self.INLengthZ.setText(self.TParts.item(0,6).text())
            
        # Voxelize geometry
        self.BVoxelizeGeometry.click()
                
    # Import Image Stack [.png, jpg, .tiff]
    elif ".png" in self.file_type or ".jpg" in self.file_type or ".jpeg" in self.file_type or ".tif" in self.file_type or ".tiff" in self.file_type:
        # Display image stack in Paraview window
        if self.file_type == ".png":
            self.imagestack_reader = paraview.simple.PNGSeriesReader(FileNames=self.is_file_names)
            paraview.simple.SetDisplayProperties(Representation = "Surface")
            text = self.INPartName.text()
            self.variable_name = f"part_{text}"
            setattr(self, self.variable_name, paraview.simple.Threshold(Input = self.imagestack_reader, ThresholdMethod = "Above Upper Threshold", UpperThreshold = 255))
            self.display = paraview.simple.Show(getattr(self, self.variable_name), self.render_view)
            paraview.simple.Hide(self.imagestack_reader)
            self.render_view.ResetCamera()
            self.render_view.StillRender()
        elif self.file_type == ".jpg" or self.file_type == ".jpeg":
            self.imagestack_reader = paraview.simple.JPEGSeriesReader(FileNames=self.is_file_names)
            paraview.simple.UpdatePipeline()
            self.imagestack_reader = paraview.simple.ExtractComponent(Input = self.imagestack_reader)
            paraview.simple.SetDisplayProperties(self.imagestack_reader, Representation = "Surface")
            text = self.INPartName.text()
            self.variable_name = f"part_{text}"
            setattr(self, self.variable_name, paraview.simple.Threshold(Input = self.imagestack_reader, ThresholdMethod = "Above Upper Threshold", UpperThreshold = 128))
            self.display = paraview.simple.Show(getattr(self, self.variable_name), self.render_view)
            paraview.simple.Hide(self.imagestack_reader)
            self.render_view.ResetCamera()
            self.render_view.StillRender()
        elif self.file_type == ".tif" or self.file_type == ".tiff":
            self.imagestack_reader = paraview.simple.TIFFSeriesReader(FileNames=self.is_file_names)
            paraview.simple.SetDisplayProperties(Representation = "Surface")
            text = self.INPartName.text()
            self.variable_name = f"part_{text}"
            setattr(self, self.variable_name, paraview.simple.Threshold(Input = self.imagestack_reader, ThresholdMethod = "Above Upper Threshold", UpperThreshold = 128))
            self.display = paraview.simple.Show(getattr(self, self.variable_name), self.render_view)
            paraview.simple.Hide(self.imagestack_reader)
            self.render_view.ResetCamera()
            self.render_view.StillRender()
            
        # create initial vtk file
        ImageToVTK(self.in_file_path, 'VTK_Geometry_' + self.INPartName.text(), self.voxelizer_dir, self.file_type)
        
        if self.BVTKFileProperties.isChecked():
            # Get the bounds of the image stack
            self.imagestack_reader.UpdatePipeline()
            self.bounds = self.imagestack_reader.GetDataInformation().GetBounds()
            self.vtkLx = round(self.bounds[1] - self.bounds[0],4)
            self.vtkLy = round(self.bounds[3] - self.bounds[2],4)
            self.vtkLz = round(self.bounds[5] - self.bounds[4],4)
            self.vtkOx = self.bounds[0]
            self.vtkOy = self.bounds[2]
            self.vtkOz = self.bounds[4]
            self.extents = self.imagestack_reader.GetDataInformation().GetExtent()
            self.vtkNx = int(self.extents[1] - self.extents[0])+1
            self.vtkNy = int(self.extents[3] - self.extents[2])+1
            self.vtkNz = int(self.extents[5] - self.extents[4])+1
        
            # Set file properties
            self.INvtkvx.setText(str(self.vtkNx))
            self.INvtkvy.setText(str(self.vtkNy))
            self.INvtkvz.setText(str(self.vtkNz))
            self.INvox.setText(str(self.vtkOx))
            self.INvoy.setText(str(self.vtkOy))
            self.INvoz.setText(str(self.vtkOz))
            self.INvlx.setText(str(self.vtkLx))
            self.INvly.setText(str(self.vtkLy))
            self.INvlz.setText(str(self.vtkLz))
            
        elif self.BVTKCustomProperties.isChecked():
            # Open voxelization tab and fill out data
            self.INvox.setText(self.TParts.item(0,1).text())
            self.INvoy.setText(self.TParts.item(0,2).text())
            self.INvoz.setText(self.TParts.item(0,3).text())
            self.INvlx.setText(self.TParts.item(0,4).text())
            self.INvly.setText(self.TParts.item(0,5).text())
            self.INvlz.setText(self.TParts.item(0,6).text())
            self.INvtkvx.setText(str(self.vtkNx))
            self.INvtkvy.setText(str(self.vtkNy))
            self.INvtkvz.setText(str(self.vtkNz))
            
        # Alter vtk file to new specifications
        
        self.BAddVTKGeometry.click()
        
    # Import Data Set [.vtk, .xdmf]
    elif ".vtk" in self.file_type or ".txt" in self.file_type:
        # VTK IMPORT
        if self.file_type == ".vtk":
            # Paraview window
            self.vtk_reader = paraview.simple.LegacyVTKReader(FileNames = self.in_file_path)
            paraview.simple.SetDisplayProperties(Representation = "Surface")
            text = self.INPartName.text()
            self.variable_name = f"part_{text}"
            setattr(self, self.variable_name, paraview.simple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
            self.display = paraview.simple.Show(getattr(self, self.variable_name), self.render_view)
            paraview.simple.Hide(self.vtk_reader)
            self.render_view.ResetCamera()
            self.render_view.StillRender()
            
            # Read what type of vtk file this is
            with open(self.in_file_path, 'r') as file:
                for index, line in enumerate(file):
                    if index == 3:
                        vtk_type = line.strip()
                        break
                        
            # Get the voxels in the vtk
            self.extents = self.vtk_reader.GetDataInformation().GetExtent()
            self.vtkNx = int(self.extents[1] - self.extents[0])
            self.vtkNy = int(self.extents[3] - self.extents[2])
            self.vtkNz = int(self.extents[5] - self.extents[4])
                        
            # STRUCTURED_POINTS vtk file
            if "STRUCTURED_POINTS" in vtk_type:
                # adjust number of voxels
                self.vtkNx += 1
                self.vtkNy += 1
                self.vtkNz += 1
            
            # Define geometry size based on previous inputs
            if self.BVTKFileProperties.isChecked():
                # Get the length of the vtk
                self.vtk_reader.UpdatePipeline()
                self.bounds = self.vtk_reader.GetDataInformation().GetBounds()
                self.vtkLx = round(self.bounds[1] - self.bounds[0],4)
                self.vtkLy = round(self.bounds[3] - self.bounds[2],4)
                self.vtkLz = round(self.bounds[5] - self.bounds[4],4)
                
                # Get the origin of the vtk
                self.vtkOx = self.bounds[0]
                self.vtkOy = self.bounds[2]
                self.vtkOz = self.bounds[4]
                
                # Open vtk tab and fill out data
                self.INvox.setText(str(self.vtkOx))
                self.INvoy.setText(str(self.vtkOy))
                self.INvoz.setText(str(self.vtkOz))
                self.INvlx.setText(str(self.vtkLx))
                self.INvly.setText(str(self.vtkLy))
                self.INvlz.setText(str(self.vtkLz))
                self.INvtkvx.setText(str(self.vtkNx))
                self.INvtkvy.setText(str(self.vtkNy))
                self.INvtkvz.setText(str(self.vtkNz))
            elif self.BVTKCustomProperties.isChecked():
                # Open voxelization tab and fill out data
                self.INvox.setText(self.TParts.item(0,1).text())
                self.INvoy.setText(self.TParts.item(0,2).text())
                self.INvoz.setText(self.TParts.item(0,3).text())
                self.INvlx.setText(self.TParts.item(0,4).text())
                self.INvly.setText(self.TParts.item(0,5).text())
                self.INvlz.setText(self.TParts.item(0,6).text())
                self.INvtkvx.setText(str(self.vtkNx))
                self.INvtkvy.setText(str(self.vtkNy))
                self.INvtkvz.setText(str(self.vtkNz))
            # Alter vtk file to new specifications
            self.BAddVTKGeometry.click()
                
        # TXT IMPORT
        elif self.file_type == ".txt":
            # output file location
            vtk_location = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
            
            # convert text file to vtk file for visualization
            los_alamos_to_vtk(self.in_file_path, vtk_location)
            
            # Paraview window
            self.vtk_reader = paraview.simple.LegacyVTKReader(FileNames = vtk_location)
            paraview.simple.SetDisplayProperties(Representation = "Surface")
            text = self.INPartName.text()
            self.variable_name = f"part_{text}"
            setattr(self, self.variable_name, self.vtk_reader)
            self.display = paraview.simple.Show(getattr(self, self.variable_name), self.render_view)
            paraview.simple.Show(self.txt_reader)
            self.render_view.ResetCamera()
            self.render_view.StillRender()
                        
            # Get the avaliable arrays for coloring
            self.vtk_reader.UpdatePipeline()
            data_maps = get_paraview_variables(self.vtk_reader)
            self.INSelectColorBy.clear()
            self.INSelectColorBy.addItems(data_maps) # add options to combo box
            
            # Get the length of the xdmf
            self.bounds = self.vtk_reader.GetDataInformation().GetBounds()
            self.vtkLx = round(self.bounds[1] - self.bounds[0],4)
            self.vtkLy = round(self.bounds[3] - self.bounds[2],4)
            self.vtkLz = round(self.bounds[5] - self.bounds[4],4)
            
            # Get the origin of the xdmf
            self.vtkOx = self.bounds[0]
            self.vtkOy = self.bounds[2]
            self.vtkOz = self.bounds[4]
            
            # Get the voxels in the xdmf
            self.extents = self.vtk_reader.GetDataInformation().GetExtent()
            self.vtkNx = int(self.extents[1] - self.extents[0])
            self.vtkNy = int(self.extents[3] - self.extents[2])
            self.vtkNz = int(self.extents[5] - self.extents[4])
            
            # Add part to the table
            row = self.TParts.rowCount()
            self.TParts.insertRow(row)
            self.TParts.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
            self.TParts.setItem(row, 1, QTableWidgetItem(str(self.vtkOx)))
            self.TParts.setItem(row, 2, QTableWidgetItem(str(self.vtkOy)))
            self.TParts.setItem(row, 3, QTableWidgetItem(str(self.vtkOz)))
            self.TParts.setItem(row, 4, QTableWidgetItem(str(self.vtkLx)))
            self.TParts.setItem(row, 5, QTableWidgetItem(str(self.vtkLy)))
            self.TParts.setItem(row, 6, QTableWidgetItem(str(self.vtkLz)))
            self.TParts.setItem(row, 7, QTableWidgetItem(str(self.vtkNx)))
            self.TParts.setItem(row, 8, QTableWidgetItem(str(self.vtkNy)))
            self.TParts.setItem(row, 9, QTableWidgetItem(str(self.vtkNz)))
            self.TParts.setItem(row, 10, QTableWidgetItem(self.in_file_path))
            
                
    # Delete previously existing part geometries
    self.TParts.removeRow(0)
    
def Reload_Geometry(self):
    if self.file_type == ".vtk":
        # Paraview window
        self.vtk_reader = paraview.simple.LegacyVTKReader(FileNames = self.in_file_path)
        paraview.simple.SetDisplayProperties(Representation = "Surface")
        text = self.INPartName.text()
        self.variable_name = f"part_{text}"
        setattr(self, self.variable_name, paraview.simple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
        self.display = paraview.simple.Show(getattr(self, self.variable_name), self.render_view)
        paraview.simple.Hide(self.vtk_reader)
        self.render_view.ResetCamera()
        self.render_view.StillRender()
    elif self.file_type == ".txt":
        # input file location
        Data_filename = self.TParts.item(0,10).text()
    
        # output file location
        vtk_location = self.voxelizer_dir + '/VTK_Geometry_' + self.TParts.item(0,0).text() + '.vtk'
        
        # convert text file to vtk file for visualization
        los_alamos_to_vtk(Data_filename, vtk_location)
        
        # Paraview window
        self.txt_reader = paraview.simple.LegacyVTKReader(FileNames = vtk_location)
        paraview.simple.SetDisplayProperties(Representation = "Surface")
        text = self.INPartName.text()
        self.variable_name = f"part_{text}"
        setattr(self, self.variable_name, self.txt_reader)
        self.display = paraview.simple.Show(getattr(self, self.variable_name), self.render_view)
        paraview.simple.Show(self.txt_reader)
        self.render_view.ResetCamera()
        self.render_view.StillRender()
        
        # Open up window to change color map
        self.SAGeometryScrollArea.verticalScrollBar().setValue(0)
        self.GeometryOptions.setCurrentIndex(3)

        # Get the avaliable arrays for coloring
        self.txt_reader.UpdatePipeline()
        data_maps = get_paraview_variables(self.txt_reader)
        self.INSelectColorBy.clear()
        self.INSelectColorBy.addItems(data_maps) # add options to combo box
    else:
        print("ERROR: can only reload visualizations for vtk files at the moment")
        
# Function to get avaliable variables from Paraview
def get_paraview_variables(data_reader):
    output = data_reader.GetClientSideObject().GetOutputDataObject(0)
    point_data = output.GetPointData()
    cell_data = output.GetCellData()
    point_arrays = [point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays())]
    cell_arrays = [cell_data.GetArrayName(i) for i in range(cell_data.GetNumberOfArrays())]
    all_arrays = point_arrays + cell_arrays
    return all_arrays
