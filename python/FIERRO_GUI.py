# To import images using resource file you must conver the .rc file to .py using (in command line):
# Compile Icons: pyside6-rcc ./fierro_gui/Icons/IconResourceFile.qrc -o IconResourceFile_rc.py

# Run Command: python -m fierro_gui.gui

from fierro_gui.ui_FIERRO_GUI import Ui_MainWindow


from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)


from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox, QFormLayout,
    QFrame, QGridLayout, QHBoxLayout, QHeaderView,
    QLabel, QLineEdit, QMainWindow, QMenu,
    QMenuBar, QPlainTextEdit, QProgressBar, QPushButton,
    QSizePolicy, QSpacerItem, QSplitter, QStackedWidget,
    QStatusBar, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget)


from PySide6.QtWidgets import ( QTableWidgetItem, QSplashScreen)

import os, glob#, sys
import time
from pathlib import Path
#import numpy as np
from PySide6.QtWidgets import (QFileDialog, QMessageBox)
#from PySide6.QtCore import (QTimer, QProcess)
import matplotlib
matplotlib.use('qt5agg')
#from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
#from matplotlib.figure import Figure
#from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
#import csv
#import re
import shutil
import paraview.simple as pvsimple
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PySide6.QtCore import *
import PySide6.QtConcurrent
from PySide6.QtGui import QDesktopServices, QMovie, QPixmap
from PySide6.QtCore import QRunnable, Slot, QThreadPool
#import fierro_voxelizer
#import fierro_mesh_builder
#import tempfile
#import time
import subprocess
from importlib import reload

from Geometry import *
from EVPFFT_Lattice import *
from Mesh_Builder_WInput import *
from FIERRO_Setup import *
from ImageToVTK import *
from TiffImageToVTK import *
from Dream3DReader import *
import DeveloperInputs

# Info on multithreading: https://www.pythonguis.com/tutorials/multithreading-pyside6-applications-qthreadpool/
class Worker(QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs

    @Slot()  # QtCore.Slot
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        self.fn(*self.args, **self.kwargs)

class FIERRO_GUI(Ui_MainWindow):

    def __init__(self):
        self.threadpool = QThreadPool()
        print("Multithreading with maximum %d threads" % self.threadpool.maxThreadCount())

    def setupUi(self, MainWindow):
        super().setupUi(MainWindow)

        app = QApplication.instance()

        # Set up multithreading
        self.threadpool = QThreadPool()
        print("Multithreading with maximum %d threads" % self.threadpool.maxThreadCount())

        # Set up tab and tool navigation
        self.NavigationMenu.currentChanged.connect(lambda: self.ToolWindow.setCurrentIndex(self.NavigationMenu.currentIndex()))
        self.INSelectGeometryImport.currentIndexChanged.connect(lambda: self.GeometryOptions.setCurrentIndex(self.INSelectGeometryImport.currentIndex()))
        self.INSelectSolverSettings.currentIndexChanged.connect(lambda: self.SolverSettingsOptions.setCurrentIndex(self.INSelectSolverSettings.currentIndex()))
        self.INSelectDefineMaterials.currentIndexChanged.connect(lambda: self.DefineMaterialsOptions.setCurrentIndex(self.INSelectDefineMaterials.currentIndex()))
        self.INSelectBoundaryConditions.currentIndexChanged.connect(lambda: self.BoundaryConditionsOptions.setCurrentIndex(self.INSelectBoundaryConditions.currentIndex()))
        self.INRunSelection.currentIndexChanged.connect(lambda: self.RunOptions.setCurrentIndex(self.INRunSelection.currentIndex()))
        self.INSelectPostprocessing.currentIndexChanged.connect(lambda: self.PostprocessingOptions.setCurrentIndex(self.INSelectPostprocessing.currentIndex()))
        self.INSelectAssignMaterials.currentIndexChanged.connect(lambda: self.AssignMaterialsOptions.setCurrentIndex(self.INSelectAssignMaterials.currentIndex()))

        # Set up pipeline selection
        selectionComboBoxes = [self.INSelectSolverSettings, self.INSelectBoundaryConditions, self.INSelectDefineMaterials, self.INRunSelection,self.INSelectPostprocessing, self.INSelectAssignMaterials]

        # When starting the gui automatically disable functions until a pipeline is chosen
        for i in range(3, self.NavigationMenu.count()):
            self.NavigationMenu.setTabEnabled(i, False)
        for i in selectionComboBoxes:
            i.setEnabled(False)

        # Function to enable certain tabs/tools based on chosen pipeline
        def SetupPipeline(selection):
            if (selection == 0): # Input/Output only
                for i in range(3, self.NavigationMenu.count()):
                    self.NavigationMenu.setTabEnabled(i, False)
                for i in selectionComboBoxes:
                    i.setEnabled(False)
            else: # Solvers
                for i in range(3, self.NavigationMenu.count()):
                    self.NavigationMenu.setTabEnabled(i, True)
                selection -= 1
                for i in selectionComboBoxes:
                    allItems = [i.model().item(j) for j in range(i.count())]
                    for k in allItems:
                        k.setEnabled(False)
                    i.setEnabled(True)
                    i.model().item(selection).setEnabled(True)
                    i.setCurrentIndex(selection)
        self.INPipelineSelection.currentIndexChanged.connect(lambda: SetupPipeline(self.INPipelineSelection.currentIndex()))

        # Paraview imports
        self.render_view = pvsimple.CreateRenderView()
        self.Paraview = QVTKRenderWindowInteractor(rw=self.render_view.GetRenderWindow(),iren=self.render_view.GetInteractor())
        self.paraviewLayout.addWidget(self.Paraview)

        # Set up loading animation widget
        self.LLoading = QLabel()
        self.LLoading.setScaledContents(1)
        self.LLoading.setFixedSize(330,100)
        self.LLoading.setParent(self.Paraview)
        self.LLoading.lower()
        self.LLoading.hide()
        #self.LLoading.setAttribute(Qt.WA_StyledBackground, True)
        #self.LLoading.setStyleSheet('background-color: #52576E;')
        #self.LLoading.setAutoFillBackground(True)
        movie = QMovie(":/Logos/Logos/FierroLoading.gif")
        movie.setSpeed(150)
        #print("Valid? ", movie.isValid())
        #print("Error: ", movie.lastErrorString())
        #print("Supported: ", movie.supportedFormats())
        self.LLoading.setMovie(movie)
        
        # Open Url in Help Menu
        def openUrl():
            url = QUrl('https://lanl.github.io/Fierro/')
            if not QDesktopServices.openUrl(url):
                QMessageBox.warning(self, 'Open Url', 'Could not open url')
        
        # Help menu
        self.actionManual.triggered.connect(openUrl)
        
        # Fierro Setup Menu
        self.actionChange_Working_Directory.triggered.connect(self.open_fierro_setup_dialog)

        def loadingAnimation():
            print ("Start load animation")
            self.LLoading.raise_()
            self.LLoading.show()
            movie.start()
            for i in range(100000):
                app.processEvents()

        def stopLoadingAnimation():
            print ("End load animation")
            movie.stop()
            self.LLoading.lower()
            self.LLoading.hide()
            
        # Upload Geometry 
        def geometry_upload_click():
            if not self.INPartName.text():
                self.warning_message("Please name the part")
            # Upload .stl or .vtk
            elif self.INSelectGeometryImport.currentIndex() == 0:
                global b3_filename
                b3_filename = QFileDialog.getOpenFileName(
                    filter="Geometry File (*.stl *.vtk)",
                )
                # Paraview window
                self.file_type = b3_filename[0][-4:-1]
                if self.file_type == '.st':
                    # Show the stl part
                    self.stl = pvsimple.STLReader(FileNames = b3_filename)
                    pvsimple.Show(self.stl, self.render_view)
                    pvsimple.ResetCamera(view=None)
                    
                    # Turn on settings
                    self.LSTLVoxelization.setEnabled(True)
                    self.LNumberOfVoxelsX.setEnabled(True)
                    self.INNumberOfVoxelsX.setEnabled(True)
                    self.LNumberOfVoxelsY.setEnabled(True)
                    self.INNumberOfVoxelsY.setEnabled(True)
                    self.LNumberOfVoxelsZ.setEnabled(True)
                    self.INNumberOfVoxelsZ.setEnabled(True)
                    self.BVoxelizeGeometry.setEnabled(True)
                    self.BStlDimensions.setEnabled(True)
                    self.BCustomDimensions.setEnabled(True)
                    self.LOriginX.setEnabled(True)
                    self.LOriginY.setEnabled(True)
                    self.LOriginZ.setEnabled(True)
                    self.INOriginX.setEnabled(True)
                    self.INOriginY.setEnabled(True)
                    self.INOriginZ.setEnabled(True)
                    
                    # Get the bounds of the stl
                    self.stl.UpdatePipeline()
                    self.bounds = self.stl.GetDataInformation().GetBounds()
                    self.stlLx = round(self.bounds[1] - self.bounds[0],4)
                    self.stlLy = round(self.bounds[3] - self.bounds[2],4)
                    self.stlLz = round(self.bounds[5] - self.bounds[4],4)
                    self.INLengthX.setText(str(self.stlLx))
                    self.INLengthY.setText(str(self.stlLy))
                    self.INLengthZ.setText(str(self.stlLz))
                    
                elif self.file_type == '.vt':
                    # Paraview window
                    print (b3_filename)
                    self.vtk_reader = pvsimple.LegacyVTKReader(FileNames = b3_filename)
                    pvsimple.SetDisplayProperties(Representation = "Surface")
                    text = self.INPartName.text()
                    self.variable_name = f"part_{text}"
                    setattr(self, self.variable_name, pvsimple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
                    pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                    pvsimple.Hide(self.vtk_reader)
                    self.render_view.ResetCamera()
                    self.render_view.StillRender()
                    
                    # Turn off settings
                    self.INNumberOfVoxelsX.setText(QCoreApplication.translate("MainWindow", u"32", None))
                    self.INNumberOfVoxelsY.setText(QCoreApplication.translate("MainWindow", u"32", None))
                    self.INNumberOfVoxelsZ.setText(QCoreApplication.translate("MainWindow", u"32", None))
                    self.LSTLVoxelization.setEnabled(False)
                    self.LNumberOfVoxelsX.setEnabled(False)
                    self.INNumberOfVoxelsX.setEnabled(False)
                    self.LNumberOfVoxelsY.setEnabled(False)
                    self.INNumberOfVoxelsY.setEnabled(False)
                    self.LNumberOfVoxelsZ.setEnabled(False)
                    self.INNumberOfVoxelsZ.setEnabled(False)
                    self.BVoxelizeGeometry.setEnabled(False)
                    
                    # Get the length of the vtk
                    self.vtk_reader.UpdatePipeline()
                    self.bounds = self.vtk_reader.GetDataInformation().GetBounds()
                    self.vtkLx = round(self.bounds[1] - self.bounds[0],4)
                    self.vtkLy = round(self.bounds[3] - self.bounds[2],4)
                    self.vtkLz = round(self.bounds[5] - self.bounds[4],4)
                    self.INLengthX.setText(str(self.vtkLx))
                    self.INLengthY.setText(str(self.vtkLy))
                    self.INLengthZ.setText(str(self.vtkLz))
                    
                    # Get the voxels in the vtk
                    self.extents = self.vtk_reader.GetDataInformation().GetExtent()
                    self.vtkNx = int(self.extents[1] - self.extents[0])
                    self.vtkNy = int(self.extents[3] - self.extents[2])
                    self.vtkNz = int(self.extents[5] - self.extents[4])
                    self.INNumberOfVoxelsX.setText(str(self.vtkNx))
                    self.INNumberOfVoxelsY.setText(str(self.vtkNy))
                    self.INNumberOfVoxelsZ.setText(str(self.vtkNz))
                    
                    # Get the origin of the vtk
                    self.INOriginX.setText(str(self.bounds[0]))
                    self.INOriginY.setText(str(self.bounds[2]))
                    self.INOriginZ.setText(str(self.bounds[4]))
                    
                    # Rename the file and save it to directory location
                    new_file_path = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
                    if b3_filename[0] != new_file_path:
                        shutil.copy(b3_filename[0], new_file_path)
                    
                    # Add part to the table
                    row = self.TParts.rowCount()
                    self.TParts.insertRow(row)
                    self.TParts.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
                    self.TParts.setItem(row, 1, QTableWidgetItem(self.INOriginX.text()))
                    self.TParts.setItem(row, 2, QTableWidgetItem(self.INOriginY.text()))
                    self.TParts.setItem(row, 3, QTableWidgetItem(self.INOriginZ.text()))
                    self.TParts.setItem(row, 4, QTableWidgetItem(self.INLengthX.text()))
                    self.TParts.setItem(row, 5, QTableWidgetItem(self.INLengthY.text()))
                    self.TParts.setItem(row, 6, QTableWidgetItem(self.INLengthZ.text()))
                    self.TParts.setItem(row, 7, QTableWidgetItem(self.INNumberOfVoxelsX.text()))
                    self.TParts.setItem(row, 8, QTableWidgetItem(self.INNumberOfVoxelsY.text()))
                    self.TParts.setItem(row, 9, QTableWidgetItem(self.INNumberOfVoxelsZ.text()))
                    self.INPartName.clear()
                    self.INOriginX.clear()
                    self.INOriginY.clear()
                    self.INOriginZ.clear()
                    self.INLengthX.clear()
                    self.INLengthY.clear()
                    self.INLengthZ.clear()
                    self.INNumberOfVoxelsX.clear()
                    self.INNumberOfVoxelsY.clear()
                    self.INNumberOfVoxelsZ.clear()
                    
                    # Add part as an option for material assignment
                    self.INPartMaterial.clear()
                    self.INPartMaterial.addItem("global")
                    for i in range(self.TParts.rowCount()):
                        self.INPartMaterial.addItem(self.TParts.item(i,0).text())
                    for i in range(self.TBasicGeometries.rowCount()):
                        self.INPartMaterial.addItem(self.TBasicGeometries.item(i,0).text())
            # Upload image stacks
            elif self.INSelectGeometryImport.currentIndex() == 1:
                global file_format
                global ct_folder_path
                ct_dialog = QFileDialog
                ct_folder_path = ct_dialog.getExistingDirectory(None, "Select Folder")
                
                # If the image stack is composed of pngs/jpgs
                print (os.listdir(ct_folder_path)[0])
                if os.listdir(ct_folder_path)[0].endswith(".png") or os.listdir(ct_folder_path)[0].endswith(".jpg") or os.listdir(ct_folder_path)[0].endswith(".jpeg"):
                    self.BTiffToStl.setEnabled(False)
                    self.BImageToVTK.setEnabled(True)
                    self.LImageFileFormat.setText("Image File Format: <i>.png/.jpg</i>")
                    file_format = ".png"
                # If the image stack is composed of tifs
                elif os.listdir(ct_folder_path)[0].endswith(".tif") or os.listdir(ct_folder_path)[0].endswith(".tiff"):
                    self.BTiffToStl.setEnabled(True)
                    self.BImageToVTK.setEnabled(True)
                    self.LImageFileFormat.setText("Image File Format: <i>.tif</i>")
                    file_format = ".tif"
                else:
                    self.warning_message("Invalid directory") 
                self.LUploadedDirectory.setText("Directory: <i>" +  os.path.basename(ct_folder_path) + "</i>")
            # Upload polycrystals
            elif self.INSelectGeometryImport.currentIndex() == 2:
                global dream_filename
                dream_filename = QFileDialog.getOpenFileName(
                    filter="Dream File (*.xdmf* *.dream3d)",
                )
                loadingAnimation()
                
                # Get xdmf file from dream3d file
                ext = os.path.splitext(dream_filename[0])
                if (ext[1] == ".dream3d"):
                    print ("File to write: ", dream_filename[0])
                    dream_filename = dream_filename[0]
                    #worker = Worker(self.RunDream3DReader) # Any other args, kwargs are passed to the run function

                    # Execute
                    #self.threadpool.start(worker)

                    Dream3DReader(str(dream_filename), 'DREAM_Geometry_' + self.INPartName.text(), self.voxelizer_dir)
                    dream_filename = self.voxelizer_dir + '/DREAM_Geometry_' + str(self.INPartName.text()) + '.xdmf'

                # Paraview window
                #stl_location = self.voxelizer_dir + '/DREAM_Geometry_' + str(self.INPartName.text()) + '.stl'
                print (dream_filename)
                # Show the stl part
                self.dream = pvsimple.XDMFReader(FileNames = dream_filename)
                if (self.dream):
                    print ("Success")
                else:
                    print ("Failure")
                pvsimple.SetDisplayProperties(Representation = "Surface")
                text = self.INPartName.text()
                self.variable_name = f"part_{text}"
                print ("var name: " + self.variable_name)
                setattr(self, self.variable_name, "")
                #pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                display = pvsimple.Show(self.dream, self.render_view)
                pvsimple.ColorBy(display, 'EulerAngles')
                self.INSelectColorBy.setEnabled(True)
                pvsimple.ResetCamera(view=None)
                row = self.TDream.rowCount()
                self.TDream.insertRow(row)
                self.TDream.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
                self.INPartName.clear()
                stopLoadingAnimation()
                def change_color_by():
                    # Dream Color By Selection
                    #print (self.INSelectColorBy.itemText(self.INSelectColorBy.currentIndex()))
                    color = None
                    match self.INSelectColorBy.itemText(self.INSelectColorBy.currentIndex()):
                        case "Euler Angles":
                            color = 'EulerAngles'
                        case "Feature Ids":
                            color = 'FeatureIds'
                        case "IPF Color":
                            color = 'IPFColor'
                        case "Phases":
                            color = 'Phases'

                    pvsimple.ColorBy(display, color)
                    pvsimple.ResetCamera(view=None)
                self.INSelectColorBy.currentIndexChanged.connect(change_color_by)
            else:
                self.warning_message('ERROR: Incorrect file type')
        self.BUploadGeometryFile.clicked.connect(geometry_upload_click)

        
            

        # Convert image stack to .vtk
        def image_to_VTK():
            print(ct_folder_path)
            if file_format == ".png":
                print ("PNG/JPG")
                # Call ImageToVTK function in python file of same name
                
                worker = Worker(self.RunImageToVTK) # Any other args, kwargs are passed to the run function

                # Execute
                self.threadpool.start(worker)
                loadingAnimation()
                
            elif file_format == ".tif":
                print ("TIFF")
                # Call TiffImageToVTK function in python file of same name, specifying to save only .vtk output
                
                worker = Worker(self.RunTiffImageToVTK) # Any other args, kwargs are passed to the run function

                # Execute
                self.threadpool.start(worker)
                loadingAnimation()
                
                #TiffImageToVTK(ct_folder_path, 'VTK_Geometry_' + self.INPartName.text(), self.voxelizer_dir, '.vtk')

            # Paraview window
            vtk_location = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
            print ("file location: " + vtk_location)
            self.vtk_reader = pvsimple.LegacyVTKReader(FileNames = vtk_location)
            pvsimple.SetDisplayProperties(Representation = "Surface")
            text = self.INPartName.text()
            self.variable_name = f"part_{text}"
            print ("var name: " + self.variable_name)
            setattr(self, self.variable_name, pvsimple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
            pvsimple.Show(getattr(self, self.variable_name), self.render_view)
            pvsimple.Hide(self.vtk_reader)
            self.render_view.ResetCamera()
            self.render_view.StillRender()
            
            self.INNumberOfVoxelsX.setText(QCoreApplication.translate("MainWindow", u"32", None))
            self.INNumberOfVoxelsY.setText(QCoreApplication.translate("MainWindow", u"32", None))
            self.INNumberOfVoxelsZ.setText(QCoreApplication.translate("MainWindow", u"32", None))
            self.LSTLVoxelization.setEnabled(False)
            self.LNumberOfVoxelsX.setEnabled(False)
            self.INNumberOfVoxelsX.setEnabled(False)
            self.LNumberOfVoxelsY.setEnabled(False)
            self.INNumberOfVoxelsY.setEnabled(False)
            self.LNumberOfVoxelsZ.setEnabled(False)
            self.INNumberOfVoxelsZ.setEnabled(False)
            self.BVoxelizeGeometry.setEnabled(False)
            
            # Get the length of the vtk
            self.vtk_reader.UpdatePipeline()
            self.bounds = self.vtk_reader.GetDataInformation().GetBounds()
            self.vtkLx = round(self.bounds[1] - self.bounds[0],4)
            self.vtkLy = round(self.bounds[3] - self.bounds[2],4)
            self.vtkLz = round(self.bounds[5] - self.bounds[4],4)
            self.INLengthX.setText(str(self.vtkLx))
            self.INLengthY.setText(str(self.vtkLy))
            self.INLengthZ.setText(str(self.vtkLz))
            
            # Get the voxels in the vtk
            self.extents = self.vtk_reader.GetDataInformation().GetExtent()
            self.vtkNx = int(self.extents[1] - self.extents[0])
            self.vtkNy = int(self.extents[3] - self.extents[2])
            self.vtkNz = int(self.extents[5] - self.extents[4])
            self.INNumberOfVoxelsX.setText(str(self.vtkNx))
            self.INNumberOfVoxelsY.setText(str(self.vtkNy))
            self.INNumberOfVoxelsZ.setText(str(self.vtkNz))
            
            # Get the origin of the vtk
            self.INOriginX.setText(str(self.bounds[0]))
            self.INOriginY.setText(str(self.bounds[2]))
            self.INOriginZ.setText(str(self.bounds[4]))
            
            # Add part to the table
            row = self.TParts.rowCount()
            self.TParts.insertRow(row)
            self.TParts.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
            self.TParts.setItem(row, 1, QTableWidgetItem(self.INOriginX.text()))
            self.TParts.setItem(row, 2, QTableWidgetItem(self.INOriginY.text()))
            self.TParts.setItem(row, 3, QTableWidgetItem(self.INOriginZ.text()))
            self.TParts.setItem(row, 4, QTableWidgetItem(self.INLengthX.text()))
            self.TParts.setItem(row, 5, QTableWidgetItem(self.INLengthY.text()))
            self.TParts.setItem(row, 6, QTableWidgetItem(self.INLengthZ.text()))
            self.TParts.setItem(row, 7, QTableWidgetItem(self.INNumberOfVoxelsX.text()))
            self.TParts.setItem(row, 8, QTableWidgetItem(self.INNumberOfVoxelsY.text()))
            self.TParts.setItem(row, 9, QTableWidgetItem(self.INNumberOfVoxelsZ.text()))
            self.INPartName.clear()
            self.INOriginX.clear()
            self.INOriginY.clear()
            self.INOriginZ.clear()
            self.INLengthX.clear()
            self.INLengthY.clear()
            self.INLengthZ.clear()
            self.INNumberOfVoxelsX.clear()
            self.INNumberOfVoxelsY.clear()
            self.INNumberOfVoxelsZ.clear()
            
            # Add part as an option for material assignment
            self.INPartMaterial.clear()
            self.INPartMaterial.addItem("global")
            for i in range(self.TParts.rowCount()):
                self.INPartMaterial.addItem(self.TParts.item(i,0).text())
            for i in range(self.TBasicGeometries.rowCount()):
                self.INPartMaterial.addItem(self.TBasicGeometries.item(i,0).text())
            stopLoadingAnimation()
        self.BImageToVTK.clicked.connect(image_to_VTK)

        # Convert tiff image stack to .stl
        def image_to_STL():
            loadingAnimation()
            print(ct_folder_path)
            print ("TIFF")
            # Call TiffImageToVTK function in python file of same name, specifying to save only .stl output
            TiffImageToVTK(ct_folder_path, 'STL_Geometry_' + self.INPartName.text(), self.voxelizer_dir, '.stl')
            
            # Paraview window
            stl_location = self.voxelizer_dir + '/STL_Geometry_' + str(self.INPartName.text()) + '.stl'
            # Show the stl part
            self.stl = pvsimple.STLReader(FileNames = stl_location)
            text = self.INPartName.text()
            self.variable_name = f"part_{text}"
            print ("var name: " + self.variable_name)
            setattr(self, self.variable_name, "")
            pvsimple.Show(self.stl, self.render_view)
            pvsimple.ResetCamera(view=None)
            
            
            # Turn on settings
            self.LSTLVoxelization.setEnabled(True)
            self.LNumberOfVoxelsX.setEnabled(True)
            self.INNumberOfVoxelsX.setEnabled(True)
            self.LNumberOfVoxelsY.setEnabled(True)
            self.INNumberOfVoxelsY.setEnabled(True)
            self.LNumberOfVoxelsZ.setEnabled(True)
            self.INNumberOfVoxelsZ.setEnabled(True)
            self.BVoxelizeGeometry.setEnabled(True)
            self.BStlDimensions.setEnabled(True)
            self.BCustomDimensions.setEnabled(True)
            self.LOriginX.setEnabled(True)
            self.LOriginY.setEnabled(True)
            self.LOriginZ.setEnabled(True)
            self.INOriginX.setEnabled(True)
            self.INOriginY.setEnabled(True)
            self.INOriginZ.setEnabled(True)
            
            # Get the bounds of the stl
            self.stl.UpdatePipeline()
            self.bounds = self.stl.GetDataInformation().GetBounds()
            self.stlLx = round(self.bounds[1] - self.bounds[0],4)
            self.stlLy = round(self.bounds[3] - self.bounds[2],4)
            self.stlLz = round(self.bounds[5] - self.bounds[4],4)
            self.INLengthX.setText(str(self.stlLx))
            self.INLengthY.setText(str(self.stlLy))
            self.INLengthZ.setText(str(self.stlLz))
            
            # Add part to the table
            row = self.TParts.rowCount()
            self.TParts.insertRow(row)
            self.TParts.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
            self.TParts.setItem(row, 1, QTableWidgetItem(self.INOriginX.text()))
            self.TParts.setItem(row, 2, QTableWidgetItem(self.INOriginY.text()))
            self.TParts.setItem(row, 3, QTableWidgetItem(self.INOriginZ.text()))
            self.TParts.setItem(row, 4, QTableWidgetItem(self.INLengthX.text()))
            self.TParts.setItem(row, 5, QTableWidgetItem(self.INLengthY.text()))
            self.TParts.setItem(row, 6, QTableWidgetItem(self.INLengthZ.text()))
            self.TParts.setItem(row, 7, QTableWidgetItem(self.INNumberOfVoxelsX.text()))
            self.TParts.setItem(row, 8, QTableWidgetItem(self.INNumberOfVoxelsY.text()))
            self.TParts.setItem(row, 9, QTableWidgetItem(self.INNumberOfVoxelsZ.text()))
            self.INPartName.clear()
            self.INOriginX.clear()
            self.INOriginY.clear()
            self.INOriginZ.clear()
            self.INLengthX.clear()
            self.INLengthY.clear()
            self.INLengthZ.clear()
            self.INNumberOfVoxelsX.clear()
            self.INNumberOfVoxelsY.clear()
            self.INNumberOfVoxelsZ.clear()
            
            # Add part as an option for material assignment
            self.INPartMaterial.clear()
            self.INPartMaterial.addItem("global")
            for i in range(self.TParts.rowCount()):
                self.INPartMaterial.addItem(self.TParts.item(i,0).text())
            for i in range(self.TBasicGeometries.rowCount()):
                self.INPartMaterial.addItem(self.TBasicGeometries.item(i,0).text())
            stopLoadingAnimation()
        self.BTiffToStl.clicked.connect(image_to_STL)
        
        # Allow for custom dimensions of stl files
        def custom_dimensions_checked():
            if self.BCustomDimensions.isChecked():
                self.LLengthX.setEnabled(True)
                self.LLengthY.setEnabled(True)
                self.LLengthZ.setEnabled(True)
                self.INLengthX.setEnabled(True)
                self.INLengthY.setEnabled(True)
                self.INLengthZ.setEnabled(True)
            else:
                self.LLengthX.setEnabled(False)
                self.LLengthY.setEnabled(False)
                self.LLengthZ.setEnabled(False)
                self.INLengthX.setEnabled(False)
                self.INLengthY.setEnabled(False)
                self.INLengthZ.setEnabled(False)
                self.INLengthX.setText(str(self.stlLx))
                self.INLengthY.setText(str(self.stlLy))
                self.INLengthZ.setText(str(self.stlLz))
        self.BCustomDimensions.clicked.connect(custom_dimensions_checked)
        self.BStlDimensions.clicked.connect(custom_dimensions_checked)
        
        # Voxelize Geometry
        def voxelize_geometry_click():
            if not self.INNumberOfVoxelsX.text() or not self.INNumberOfVoxelsY.text() or not self.INNumberOfVoxelsZ.text() or not self.INLengthX.text() or not self.INLengthY.text() or not self.INLengthZ.text() or not self.INOriginX.text() or not self.INOriginY.text() or not self.INOriginZ.text():
                self.warning_message('ERROR: Number of voxels NOT defined')
            else:
                loadingAnimation()
                # Run voxelization executable
                reload(DeveloperInputs)
                if self.UserConfig == "Developer":
                    executable_path = DeveloperInputs.fierro_voxelizer_exe
                elif self.UserConfig == "User":
                    executable_path = "fierro-voxelizer"
                vtk_location = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
                arguments = [b3_filename[0], vtk_location, self.INNumberOfVoxelsX.text(), self.INNumberOfVoxelsY.text(), self.INNumberOfVoxelsZ.text(), self.INOriginX.text(), self.INOriginY.text(), self.INOriginZ.text(), self.INLengthX.text(), self.INLengthY.text(), self.INLengthZ.text()]
                command = [executable_path] + arguments
                try:
                    process = subprocess.Popen(command)
                    process.wait()
                except Exception as e:
                    self.warning_message("ERROR: fierro-voxelizer executable")
                    return
                
                # Paraview window
                pvsimple.Delete(self.stl)
                self.vtk_reader = pvsimple.LegacyVTKReader(FileNames = vtk_location)
                pvsimple.SetDisplayProperties(Representation = "Surface")
                text = self.INPartName.text()
                self.variable_name = f"part_{text}"
                setattr(self, self.variable_name, pvsimple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
                pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                pvsimple.Hide(self.vtk_reader)
                self.render_view.ResetCamera()
                self.render_view.StillRender()
                
                
                # Add part to the table
                row = self.TParts.rowCount()
                self.TParts.insertRow(row)
                self.TParts.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
                self.TParts.setItem(row, 1, QTableWidgetItem(self.INOriginX.text()))
                self.TParts.setItem(row, 2, QTableWidgetItem(self.INOriginY.text()))
                self.TParts.setItem(row, 3, QTableWidgetItem(self.INOriginZ.text()))
                self.TParts.setItem(row, 4, QTableWidgetItem(self.INLengthX.text()))
                self.TParts.setItem(row, 5, QTableWidgetItem(self.INLengthY.text()))
                self.TParts.setItem(row, 6, QTableWidgetItem(self.INLengthZ.text()))
                self.TParts.setItem(row, 7, QTableWidgetItem(self.INNumberOfVoxelsX.text()))
                self.TParts.setItem(row, 8, QTableWidgetItem(self.INNumberOfVoxelsY.text()))
                self.TParts.setItem(row, 9, QTableWidgetItem(self.INNumberOfVoxelsZ.text()))
                self.INPartName.clear()
                self.INOriginX.clear()
                self.INOriginY.clear()
                self.INOriginZ.clear()
                self.INLengthX.clear()
                self.INLengthY.clear()
                self.INLengthZ.clear()
                self.INNumberOfVoxelsX.clear()
                self.INNumberOfVoxelsY.clear()
                self.INNumberOfVoxelsZ.clear()
                
                # Add part as an option for material assignment
                self.INPartMaterial.clear()
                self.INPartMaterial.addItem("global")
                for i in range(self.TParts.rowCount()):
                    self.INPartMaterial.addItem(self.TParts.item(i,0).text())
                for i in range(self.TBasicGeometries.rowCount()):
                    self.INPartMaterial.addItem(self.TBasicGeometries.item(i,0).text())
            stopLoadingAnimation()
        self.BVoxelizeGeometry.clicked.connect(voxelize_geometry_click)
        
        # Delete any previously loaded geometries from table
        global table
        table = self.TParts
        def delete_part():
            print ("Delete Table: ", table)
            current_row = table.currentRow()
            if current_row < 0:
                return QMessageBox.warning(QMessageBox(),"Warning","Please select a part to delete")

            button = QMessageBox.question(
                QMessageBox(),
                'Confirmation',
                'Are you sure that you want to delete the selected row?',
                QMessageBox.Yes |
                QMessageBox.No
            )
            if button == QMessageBox.StandardButton.Yes:
                self.newvar = "part_" + table.item(current_row,0).text()
                pvsimple.Delete(getattr(self, self.newvar))
                self.render_view.ResetCamera()
                self.render_view.StillRender()
                table.removeRow(current_row)
                
                # delete from material assignment options
                self.INPartMaterial.clear()
                self.INPartMaterial.addItem("global")
                for i in range(table.rowCount()):
                    self.INPartMaterial.addItem(table.item(i,0).text())
                for i in range(self.TBasicGeometries.rowCount()):
                    self.INPartMaterial.addItem(self.TBasicGeometries.item(i,0).text())

        # Function to set what table to delete part from
        def set_table(table_selection):
            global table
            table = table_selection
            print ("Table: ", table)
            delete_part()
        self.BDeleteGeometry.clicked.connect(lambda: set_table(self.TParts))
        self.BDeleteDream.clicked.connect(lambda: set_table(self.TDream))
            
        # Global Mesh Generation
        # Deactivate 2D option for now
        self.INDimension.model().item(1).setEnabled(False)
        def mesh_class():
            if str(self.INCoordinateSystem.currentText()) == 'Rectangular':
                if str(self.INDimension.currentText()) == '2D':
                    self.MeshInputs2.setCurrentIndex(1)
                else:
                    self.MeshInputs2.setCurrentIndex(0)
            else:
                if str(self.INDimension.currentText()) == '2D':
                    self.MeshInputs2.setCurrentIndex(2)
                else:
                    self.MeshInputs2.setCurrentIndex(3)
        self.INCoordinateSystem.currentIndexChanged.connect(mesh_class)
        self.INDimension.currentIndexChanged.connect(mesh_class)
            
        # Write input file for mesh builder
        def global_mesh_click():
            Mesh_Builder_WInput(self, self.GLOBAL_MESH, self.global_mesh_dir)
        self.BGenerateGlobalMesh.clicked.connect(global_mesh_click)
            
        # ======= EVPFFT SOLVER BULK RESPONSE PIPELINE =======
        self.EVPFFTSteps = 3 # Define the number of steps for evpfft
        EVPFFT_Lattice(self)
            
        # ======= EXPLICIT SOLVER SGH PIPELINE =======
        Geometry(self)

    def RunImageToVTK(self):
        print ("RunImageToVTK")
        ImageToVTK(ct_folder_path, 'VTK_Geometry_' + self.INPartName.text(), self.voxelizer_dir) # Any other args, kwargs are passed to the run function

    def RunTiffImageToVTK(self):
        print ("RunTiffImageToVTK")
        TiffImageToVTK(ct_folder_path, 'VTK_Geometry_' + self.INPartName.text(), self.voxelizer_dir, '.vtk')

    def RunDream3DReader(self):
        print ("RunDream3DReader")
        Dream3DReader(str(dream_filename), 'DREAM_Geometry_' + self.INPartName.text(), self.voxelizer_dir)
        
    # ========== WARNING MESSAGE ============
    def warning_message(self, msg):
        message = QMessageBox()
        message.setText(msg)
        message.exec()
        
    # ========== FIERRO SETUP ============
    def open_fierro_setup_dialog(self, gself):
        dialog = FierroSetup(gself)
        dialog.settingChanged.connect(self.handleSettingChanged)
        dialog.setWindowModality(Qt.WindowModal)
        if dialog.exec() == QDialog.Accepted:
            self.directory = dialog.get_directory()
            # Create a temp file for the voxelization
            self.voxelizer_dir = os.path.join(self.directory, 'voxelizer')
            os.makedirs(self.voxelizer_dir, exist_ok=True)
            
            # Create temp files for evpfft
            self.evpfft_dir = os.path.join(self.directory, 'evpfft')
            os.makedirs(self.evpfft_dir, exist_ok=True)
            self.ELASTIC_PARAMETERS_0 = os.path.join(self.evpfft_dir, 'elastic_parameters_0.txt')
            self.ELASTIC_PARAMETERS_1 = os.path.join(self.evpfft_dir, 'elastic_parameters_1.txt')
            self.PLASTIC_PARAMETERS = os.path.join(self.evpfft_dir, 'plastic_parameters.txt')
            self.EVPFFT_INPUT = os.path.join(self.evpfft_dir, 'evpfft_lattice_input.txt')
            
            # Create a temp file for the global mesh
            self.global_mesh_dir = os.path.join(self.directory, 'global_mesh')
            os.makedirs(self.global_mesh_dir, exist_ok=True)
            self.GLOBAL_MESH = os.path.join(self.global_mesh_dir, 'global_mesh.yaml')
            self.GLOBAL_MESH_OUTPUT = os.path.join(self.global_mesh_dir, 'mesh.vtk')
            
            # Create temp files for SGH
            sgh_dir = os.path.join(self.directory, 'sgh')
            os.makedirs(sgh_dir, exist_ok=True)
            self.EXPLICIT_SGH_INPUT = os.path.join(sgh_dir, 'explicit_sgh_input.yaml')
        
        else:
            self.warning_message("ERROR: Working directory was not defined")
    
    # Output which user profile the user selected
    def handleSettingChanged(self, setting):
        if setting == "User":
            self.UserConfig = "User"
        elif setting == "Developer":
            self.UserConfig = "Developer"

