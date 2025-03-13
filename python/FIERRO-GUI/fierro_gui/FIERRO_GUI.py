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
    QVBoxLayout, QWidget, QToolTip)


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
from PySide6.QtGui import QDesktopServices, QMovie, QPixmap, QColor, QBrush, QFont
from PySide6.QtCore import QRunnable, Slot, QThreadPool
#import fierro_voxelizer
#import fierro_mesh_builder
#import tempfile
#import time
import subprocess
from importlib import reload

from SGH import *
from Homogenization import *
from Bulk_Forming import *
from Mesh_Builder_WInput import *
from FIERRO_Setup import *
from ImageToVTK import *
from Save_Load import *
from LAFFT2VTK import *
#from TiffImageToVTK import *
#from Dream3DReader import *
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
        self.INSelectSolverSettings.currentIndexChanged.connect(lambda: self.SolverSettingsOptions.setCurrentIndex(self.INSelectSolverSettings.currentIndex()))
        self.INSelectDefineMaterials.currentIndexChanged.connect(lambda: self.DefineMaterialsOptions.setCurrentIndex(self.INSelectDefineMaterials.currentIndex()))
        self.INSelectBoundaryConditions.currentIndexChanged.connect(lambda: self.BoundaryConditionsOptions.setCurrentIndex(self.INSelectBoundaryConditions.currentIndex()))
        self.INRunSelection.currentIndexChanged.connect(lambda: self.RunOptions.setCurrentIndex(self.INRunSelection.currentIndex()))
        self.INSelectPostprocessing.currentIndexChanged.connect(lambda: self.PostprocessingOptions.setCurrentIndex(self.INSelectPostprocessing.currentIndex()))
        self.INSelectAssignMaterials.currentIndexChanged.connect(lambda: self.AssignMaterialsOptions.setCurrentIndex(self.INSelectAssignMaterials.currentIndex()))
        self.INPipelineSelection.currentIndexChanged.connect(lambda: self.PipelineExtras.setCurrentIndex(self.INPipelineSelection.currentIndex()))

        # Set up pipeline selection
        selectionComboBoxes = [self.INSelectGeometry, self.INSelectSolverSettings, self.INSelectBoundaryConditions, self.INSelectDefineMaterials, self.INRunSelection,self.INSelectPostprocessing, self.INSelectAssignMaterials]

        # When starting the gui automatically disable functions until a pipeline is chosen
        for i in range(2, self.NavigationMenu.count()):
            self.NavigationMenu.setTabEnabled(i, False)
        for i in selectionComboBoxes:
            i.setEnabled(False)
            
        # Set up geometry import flag
        self.geometry_tab = 0;

        # Function to enable certain tabs/tools based on chosen pipeline
        def SetupPipeline(selection):
            # Enable tabs
            for i in range(2, self.NavigationMenu.count()):
                self.NavigationMenu.setTabEnabled(i, True)
            selection -= 1
            for i in selectionComboBoxes:
                allItems = [i.model().item(j) for j in range(i.count())]
                for k in allItems:
                    k.setEnabled(False)
                i.setEnabled(True)
                i.model().item(selection).setEnabled(True)
                i.setCurrentIndex(selection)
                
            # Enable solver specific items
            if (selection == 0): # Explicit SGH Solver
                # Define geometry imports
                self.INSelectGeometryImport.clear()
                self.INSelectGeometryImport.addItem("Import Geometry (.stl, .vtk)")
                self.INSelectGeometryImport.addItem("Import Image Stack (.png)")
#                self.INSelectGeometryImport.addItem("Import Polycrstalline Data Set (.txt)")
                self.INSelectGeometryImport.addItem("Create Basic Part")
            elif (selection == 1): # Homogenization Solver
                # Define geometry imports
                self.INSelectGeometryImport.clear()
                self.INSelectGeometryImport.addItem("Import Geometry (.stl, .vtk)")
                self.INSelectGeometryImport.addItem("Import Image Stack (.png, .jpg, .tif)")
                self.INSelectGeometryImport.addItem("Import Polycrystalline Data Set (.txt)")
                
                # Turn off tabs
                self.NavigationMenu.setTabEnabled(3, False)
                self.NavigationMenu.setTabEnabled(5, False)
                self.MaterialMenu.setTabEnabled(1, False)
            elif (selection == 2): # Bulk Forming
                # Define geometry imports
                self.INSelectGeometryImport.clear()
                self.INSelectGeometryImport.addItem("Import Polycrystalline Data Set (.txt)")
                self.INSelectGeometryImport.addItem("Import Geometry (.stl, .vtk)")
                # Turn off tabs
                self.NavigationMenu.setTabEnabled(3, False)
                self.MaterialMenu.setTabEnabled(1, True)
        self.INPipelineSelection.currentIndexChanged.connect(lambda: SetupPipeline(self.INPipelineSelection.currentIndex()))
        
        # Disable SGH Solver Pipeline **just for now**
        self.INPipelineSelection.model().item(0).setEnabled(False)
        self.INPipelineSelection.model().item(1).setEnabled(False)
        
        # Make information text bigger
        font = QToolTip.font()
        font.setPointSize(16)
        QToolTip.setFont(font)

        # Paraview imports
        self.render_view = pvsimple.CreateRenderView()
        self.Paraview = QVTKRenderWindowInteractor(rw=self.render_view.GetRenderWindow(),iren=self.render_view.GetInteractor())
        self.paraviewLayout.addWidget(self.Paraview)

        # Set up loading animation widget
        self.LLoading = QLabel()
        self.LLoading.setScaledContents(1)
        self.LLoading.setFixedSize(792,240)
        self.LLoading.lower()
        self.LLoading.hide()
        self.paraviewLayout.addWidget(self.LLoading, alignment=Qt.AlignBottom)
        #self.LLoading.setParent(self.Paraview)
        #self.LLoading.setAttribute(Qt.WA_StyledBackground, True)
        #self.LLoading.setStyleSheet('background-color: #52576E;')
        #self.LLoading.setAutoFillBackground(True)
        movie = QMovie(u":/Logos/Logos/FierroLoading.gif")
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
        
        # Save menu
        self.actionSaveAs.triggered.connect(self.open_save_dialog)
        
        # Load menu
        self.actionOpen.triggered.connect(self.open_load_dialog)

        
        def loadingAnimation():
            hi = 1
#            print ("Start load animation")
#            self.LLoading.raise_()
#            self.LLoading.show()
#            movie.start()
#            for i in range(100000):
#                app.processEvents()

        def stopLoadingAnimation():
            hi = 1
#            print ("End load animation")
#            movie.stop()
#            self.LLoading.lower()
#            self.LLoading.hide()

        # Define Units
        self.old_units = self.INUnits.currentText()
        self.new_units = 'new'
        self.length_variables = [self.INOriginX, self.INOriginY, self.INOriginZ, self.INLengthX, self.INLengthY, self.INLengthZ, self.INvox, self.INvoy, self.INvoz, self.INvlx, self.INvly, self.INvlz]
        self.pressure_variables = [self.INYoungsModulus, self.INEip, self.INEop, self.INGop, self.INEx, self.INEy, self.INEz, self.INGxy, self.INGxz, self.INGyz]
        def units():
            self.new_units = self.INUnits.currentText()
            if self.old_units != self.new_units:
                # Convert from m->mm
                if 'm' in self.old_units and 'mm' in self.new_units:
                    j = 1
                    for i in self.length_variables:
                        # Change labels
                        length = getattr(self, f'Ulength{j}', None)
                        if length:
                            length.setText("[mm]")
                        j+=1
                        # Change variable vlaues
                        if i.text():
                            m_mm = float(i.text())*1000
                            i.setText(f'{m_mm}')
                    # Change part geometry definition table
                    for j in range(self.TParts.rowCount()):
                        for i in range(1,7):
                            val = float(self.TParts.item(j,i).text())*1000
                            self.TParts.setItem(j,i,QTableWidgetItem(f'{val}'))
                            
                # Convert from mm->m
                if 'mm' in self.old_units and 'm' in self.new_units:
                    j = 1
                    for i in self.length_variables:
                        # Change labels
                        length = getattr(self, f'Ulength{j}', None)
                        if length:
                            length.setText("[m]")
                        j+=1
                        # Change variable values
                        if i.text():
                            mm_m = float(i.text())/1000
                            i.setText(f'{mm_m}')
                    # Change part geometry definition table
                    for j in range(self.TParts.rowCount()):
                        for i in range(1,7):
                            val = float(self.TParts.item(j,i).text())/1000
                            self.TParts.setItem(j,i,QTableWidgetItem(f'{val}'))
                            
                # Convert from Pa->MPa
                if 'Pa' in self.old_units and 'MPa' in self.new_units:
                    j = 1
                    for i in self.pressure_variables:
                        # Change labels
                        pressure = getattr(self, f'UPressure{j}', None)
                        if pressure:
                            pressure.setText("[MPa]")
                        j+=1
                        # Change variable values
                        if i.text():
                            Pa_MPa = float(i.text())/1000000
                            i.setText(f'{Pa_MPa}')
                    # Change table - material definition values
                    for j in range(self.TMaterials.rowCount()):
                        for i in range(2,23):
                            if self.TMaterials.item(j,i) is not None:
                                val = float(self.TMaterials.item(j,i).text())/1000000
                                self.TMaterials.setItem(j,i,QTableWidgetItem(f'{val}'))
                            if self.TMaterials_2.item(j,i) is not None:
                                val = float(self.TMaterials_2.item(j,i).text())/1000000
                                self.TMaterials_2.setItem(j,i,QTableWidgetItem(f'{val}'))
                    # Change table - homogenized constants values
                    if self.THomogenization.item(3,0) is not None:
                        for i in [0,1,2,6,7,8]:
                            val = float(self.THomogenization.item(i,0).text())/1000000
                            self.THomogenization.setItem(i,0,QTableWidgetItem(f'{val}'))
                    # Homogenization table labels
                    self.THomogenization.setItem(-1,3,QTableWidgetItem("[MPa]"))
                    self.THomogenization.setItem(0,3,QTableWidgetItem("[MPa]"))
                    self.THomogenization.setItem(1,3,QTableWidgetItem("[MPa]"))
                    self.THomogenization.setItem(5,3,QTableWidgetItem("[MPa]"))
                    self.THomogenization.setItem(6,3,QTableWidgetItem("[MPa]"))
                    self.THomogenization.setItem(7,3,QTableWidgetItem("[MPa]"))
                    # Adjust homogenization solver settings
                    self.INErrorTolerance.setText('0.00001')
                    
                # Convert from MPa->Pa
                if 'MPa' in self.old_units and 'Pa' in self.new_units:
                    j = 1
                    for i in self.pressure_variables:
                        # Change labels
                        pressure = getattr(self, f'UPressure{j}', None)
                        if pressure:
                            pressure.setText("[Pa]")
                        j+=1
                        # Change variable values
                        if i.text():
                            MPa_Pa = float(i.text())*1000000
                            i.setText(f'{MPa_Pa}')
                    # Change table - material definition values
                    for j in range(self.TMaterials.rowCount()):
                        for i in range(2,23):
                            val = float(self.TMaterials.item(j,i).text())*1000000
                            self.TMaterials.setItem(j,i,QTableWidgetItem(f'{val}'))
                            val = float(self.TMaterials_2.item(j,i).text())*1000000
                            self.TMaterials_2.setItem(j,i,QTableWidgetItem(f'{val}'))
                    # Change table - homogenized constants values
                    if self.THomogenization.item(3,0) is not None:
                        for i in [0,1,2,6,7,8]:
                            val = float(self.THomogenization.item(i,0).text())*1000000
                            self.THomogenization.setItem(i,0,QTableWidgetItem(f'{val}'))
                    # Homogenization table labels
                    self.THomogenization.setItem(-1,3,QTableWidgetItem("[Pa]"))
                    self.THomogenization.setItem(0,3,QTableWidgetItem("[Pa]"))
                    self.THomogenization.setItem(1,3,QTableWidgetItem("[Pa]"))
                    self.THomogenization.setItem(5,3,QTableWidgetItem("[Pa]"))
                    self.THomogenization.setItem(6,3,QTableWidgetItem("[Pa]"))
                    self.THomogenization.setItem(7,3,QTableWidgetItem("[Pa]"))
                    # Adjust homogenization solver settings
                    self.INErrorTolerance.setText('10')
            self.old_units = self.INUnits.currentText()
        self.INUnits.currentIndexChanged.connect(units)
        self.units = units
            
        # Upload Geometry
        def geometry_upload_click():
            global file_type
            self.geometry_tab = 1
            # Check for part name
            if not self.INPartName.text():
                self.warning_message("Please name the part")
            
            # Import Geometry [.stl, .vtk]
            elif ".stl" in self.INSelectGeometryImport.currentText() or ".vtk" in self.INSelectGeometryImport.currentText():
                # Get filename
                global Geometry_filename
                Geometry_filename, _ = QFileDialog.getOpenFileName(
                    filter="Geometry File (*.stl *.vtk)",
                )
                file_type = os.path.splitext(Geometry_filename)[1].lower()
                
                # STL IMPORT
                if ".stl" in file_type:
                    # Show the stl part
                    self.stl = pvsimple.STLReader(FileNames = Geometry_filename)
                    self.display = pvsimple.Show(self.stl, self.render_view)
                    pvsimple.ResetCamera()
                    
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
                    self.SAGeometryScrollArea.verticalScrollBar().setValue(0)
                    self.GeometryOptions.setCurrentIndex(1)
                    self.INOriginX.setText(str(self.stlOx))
                    self.INOriginY.setText(str(self.stlOy))
                    self.INOriginZ.setText(str(self.stlOz))
                    self.INLengthX.setText(str(self.stlLx))
                    self.INLengthY.setText(str(self.stlLy))
                    self.INLengthZ.setText(str(self.stlLz))
                
                # VTK IMPORT
                elif file_type == ".vtk":
                    # Set initial structured_points flag to 0
                    self.structured_flag = 0
                    
                    # Paraview window
                    self.vtk_reader = pvsimple.LegacyVTKReader(FileNames = Geometry_filename)
                    pvsimple.SetDisplayProperties(Representation = "Surface")
                    text = self.INPartName.text()
                    self.variable_name = f"part_{text}"
                    setattr(self, self.variable_name, pvsimple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
                    self.display = pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                    pvsimple.Hide(self.vtk_reader)
                    self.render_view.ResetCamera()
                    self.render_view.StillRender()
                    
                    # Rename the file and save it to directory location
                    self.new_file_path = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
                    if Geometry_filename != self.new_file_path:
                        shutil.copy(Geometry_filename, self.new_file_path)
                    
                    # Read what type of vtk file this is
                    with open(self.new_file_path, 'r') as file:
                        for index, line in enumerate(file):
                            if index == 3:
                                vtk_type = line.strip()
                                break
                                
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
                        
                        # structured_points flag
                        self.structured_flag = 1
                    
                    # Open up window to alter dimensions
                    self.GeometryOptions.setCurrentIndex(5)
                    
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
                        
            # Import Image Stack [.png, jpg, .tiff]
            elif ".png" in self.INSelectGeometryImport.currentText() or ".jpg" in self.INSelectGeometryImport.currentText() or ".tif" in self.INSelectGeometryImport.currentText():
                # get global variables
                global ct_folder_path
                ct_dialog = QFileDialog
                ct_folder_path = ct_dialog.getExistingDirectory(None, "Select Folder")

                # get file type
                file_type = os.path.splitext(os.listdir(ct_folder_path)[1])[1]
                if file_type == "":
                    warning_message("WARNING: might be hidden files in foler, couldn't determine file type")
                
                # Display image stack in Paraview window
                all_files = os.listdir(ct_folder_path)
                file_names = sorted([os.path.join(ct_folder_path, f) for f in all_files if f.endswith(file_type)])
                if not file_names:
                    self.warning_message("No files found with the image extension in the folder.")
                if file_type == ".png":
                    self.imagestack_reader = pvsimple.PNGSeriesReader(FileNames=file_names)
                    pvsimple.SetDisplayProperties(Representation = "Surface")
                    text = self.INPartName.text()
                    self.variable_name = f"part_{text}"
                    setattr(self, self.variable_name, pvsimple.Threshold(Input = self.imagestack_reader, ThresholdMethod = "Above Upper Threshold", UpperThreshold = 255))
                    self.display = pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                    pvsimple.Hide(self.imagestack_reader)
                    self.render_view.ResetCamera()
                    self.render_view.StillRender()
                elif file_type == ".jpg" or file_type == ".jpeg":
                    self.imagestack_reader = pvsimple.JPEGSeriesReader(FileNames=file_names)
                    pvsimple.UpdatePipeline()
                    self.imagestack_reader = pvsimple.ExtractComponent(Input = self.imagestack_reader)
                    pvsimple.SetDisplayProperties(self.imagestack_reader, Representation = "Surface")
                    text = self.INPartName.text()
                    self.variable_name = f"part_{text}"
                    setattr(self, self.variable_name, pvsimple.Threshold(Input = self.imagestack_reader, ThresholdMethod = "Above Upper Threshold", UpperThreshold = 128))
                    self.display = pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                    pvsimple.Hide(self.imagestack_reader)
                    self.render_view.ResetCamera()
                    self.render_view.StillRender()
                elif file_type == ".tif" or file_type == ".tiff":
                    self.imagestack_reader = pvsimple.TIFFSeriesReader(FileNames=file_names)
                    pvsimple.SetDisplayProperties(Representation = "Surface")
                    text = self.INPartName.text()
                    self.variable_name = f"part_{text}"
                    setattr(self, self.variable_name, pvsimple.Threshold(Input = self.imagestack_reader, ThresholdMethod = "Above Upper Threshold", UpperThreshold = 128))
                    self.display = pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                    pvsimple.Hide(self.imagestack_reader)
                    self.render_view.ResetCamera()
                    self.render_view.StillRender()
                
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
                
                # Voxelize PNG files
                self.new_file_path = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
                image_to_VTK()
                self.structured_flag = 1
                
                # Open up window to alter dimensions
                self.GeometryOptions.setCurrentIndex(5)
                
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
                
            # Import Data Set [.vtk, .xdmf]
            elif ".txt" in self.INSelectGeometryImport.currentText():
                # Get filename
                global Data_filename
                Data_filename, _ = QFileDialog.getOpenFileName(
                    filter="Data Set File (*.txt)",
                )
                file_type = os.path.splitext(Data_filename)[1].lower()
                        
                # XDMF IMPORT
                if file_type == ".xdmf":
                    # Paraview window
                    self.xdmf_reader = pvsimple.XDMFReader(FileNames = Data_filename)
                    pvsimple.SetDisplayProperties(Representation = "Surface")
                    text = self.INPartName.text()
                    self.variable_name = f"part_{text}"
                    setattr(self, self.variable_name, self.xdmf_reader)
                    self.display = pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                    self.render_view.ResetCamera()
                    self.render_view.StillRender()
                    
                    # Open up window to change color map
                    self.SAGeometryScrollArea.verticalScrollBar().setValue(0)
                    self.GeometryOptions.setCurrentIndex(3)
                    
                    # Get the avaliable arrays for coloring
                    self.xdmf_reader.UpdatePipeline()
                    data_maps = get_paraview_variables(self.xdmf_reader)
                    self.INSelectColorBy.clear()
                    self.INSelectColorBy.addItems(data_maps) # add options to combo box
                    
                    # Get the length of the xdmf
                    self.xdmf_reader.UpdatePipeline()
                    self.bounds = self.xdmf_reader.GetDataInformation().GetBounds()
                    self.xdmfLx = round(self.bounds[1] - self.bounds[0],4)
                    self.xdmfLy = round(self.bounds[3] - self.bounds[2],4)
                    self.xdmfLz = round(self.bounds[5] - self.bounds[4],4)
                    
                    # Get the origin of the xdmf
                    self.xdmfOx = self.bounds[0]
                    self.xdmfOy = self.bounds[2]
                    self.xdmfOz = self.bounds[4]
                    
                    # Get the voxels in the xdmf
                    self.extents = self.xdmf_reader.GetDataInformation().GetExtent()
                    self.xdmfNx = int(self.extents[1] - self.extents[0])
                    self.xdmfNy = int(self.extents[3] - self.extents[2])
                    self.xdmfNz = int(self.extents[5] - self.extents[4])
                    
                    # Add part to the table
                    row = self.TParts.rowCount()
                    self.TParts.insertRow(row)
                    self.TParts.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
                    self.TParts.setItem(row, 1, QTableWidgetItem(str(self.xdmfOx)))
                    self.TParts.setItem(row, 2, QTableWidgetItem(str(self.xdmfOy)))
                    self.TParts.setItem(row, 3, QTableWidgetItem(str(self.xdmfOz)))
                    self.TParts.setItem(row, 4, QTableWidgetItem(str(self.xdmfLx)))
                    self.TParts.setItem(row, 5, QTableWidgetItem(str(self.xdmfLy)))
                    self.TParts.setItem(row, 6, QTableWidgetItem(str(self.xdmfLz)))
                    self.TParts.setItem(row, 7, QTableWidgetItem(str(self.xdmfNx)))
                    self.TParts.setItem(row, 8, QTableWidgetItem(str(self.xdmfNy)))
                    self.TParts.setItem(row, 9, QTableWidgetItem(str(self.xdmfNz)))
                    self.INPartName.clear()
                    
                    # Add part as an option for material assignment
                    self.INRegion.clear()
                    self.INRegion_2.clear()
                    self.INRegion_3.clear()
                    for i in range(self.TParts.rowCount()):
                        self.INRegion.addItem(self.TParts.item(i,0).text())
                        self.INRegion_2.addItem(self.TParts.item(i,0).text())
                        self.INRegion_3.addItem(self.TParts.item(i,0).text())
                    for i in range(self.TBasicGeometries.rowCount()):
                        self.INRegion.addItem(self.TBasicGeometries.item(i,0).text())
                    self.INRegion.addItem("global")
                    self.INRegion_2.addItem("global")
                    self.INRegion_3.addItem("global")
                    
                # Text file
                elif file_type == ".txt":
                    # output file location
                    vtk_location = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
                    
                    # convert text file to vtk file for visualization
                    los_alamos_to_vtk(Data_filename, vtk_location)
                    
                    # Paraview window
                    self.txt_reader = pvsimple.LegacyVTKReader(FileNames = vtk_location)
                    pvsimple.SetDisplayProperties(Representation = "Surface")
                    text = self.INPartName.text()
                    self.variable_name = f"part_{text}"
                    setattr(self, self.variable_name, self.txt_reader)
                    self.display = pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                    pvsimple.Show(self.txt_reader)
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
                    
                    # Get the length of the xdmf
                    self.bounds = self.txt_reader.GetDataInformation().GetBounds()
                    self.vtkLx = round(self.bounds[1] - self.bounds[0],4)
                    self.vtkLy = round(self.bounds[3] - self.bounds[2],4)
                    self.vtkLz = round(self.bounds[5] - self.bounds[4],4)
                    
                    # Get the origin of the xdmf
                    self.vtkOx = self.bounds[0]
                    self.vtkOy = self.bounds[2]
                    self.vtkOz = self.bounds[4]
                    
                    # Get the voxels in the xdmf
                    self.extents = self.txt_reader.GetDataInformation().GetExtent()
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
                    self.TParts.setItem(row, 10, QTableWidgetItem(Data_filename))
                    self.INPartName.clear()
                    
                    # Add part as an option for material assignment
                    self.INRegion.clear()
                    self.INRegion_2.clear()
                    self.INRegion_3.clear()
                    for i in range(self.TParts.rowCount()):
                        self.INRegion.addItem(self.TParts.item(i,0).text())
                        self.INRegion_2.addItem(self.TParts.item(i,0).text())
                        self.INRegion_3.addItem(self.TParts.item(i,0).text())
                    for i in range(self.TBasicGeometries.rowCount()):
                        self.INRegion.addItem(self.TBasicGeometries.item(i,0).text())
                    
        self.BUploadGeometryFile.clicked.connect(geometry_upload_click)
        
        # Warn User if no geometry was uploaded
        def no_geometry():
            if self.TParts.rowCount() == 0 and self.geometry_tab == 1:
                self.warning_message('WARNING: Model contains no geometry (the part was not fully imported into the model)')
            self.geometry_tab = 0;
        self.NavigationMenu.currentChanged.connect(no_geometry)
        
        # Function to get avaliable variables from Paraview
        def get_paraview_variables(data_reader):
            output = data_reader.GetClientSideObject().GetOutputDataObject(0)
            point_data = output.GetPointData()
            cell_data = output.GetCellData()
            point_arrays = [point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays())]
            cell_arrays = [cell_data.GetArrayName(i) for i in range(cell_data.GetNumberOfArrays())]
            all_arrays = point_arrays + cell_arrays
            return all_arrays
            
        # Change color map of Paraview variables
        def color_map():
            pvsimple.ColorBy(self.display, self.INSelectColorBy.currentText())
            pvsimple.ResetCamera(view=None)
        self.INSelectColorBy.currentIndexChanged.connect(color_map)

        # Convert image stack to .vtk
        def image_to_VTK():
            # Call ImageToVTK function in python file of same name
            worker = Worker(self.RunImageToVTK) # Any other args, kwargs are passed to the run function

            # Execute
            self.threadpool.start(worker)
            loadingAnimation()
            stopLoadingAnimation()
        
        # Voxelize Geometry
        def voxelize_geometry_click():
            if not self.INNumberOfVoxelsX.text() or not self.INNumberOfVoxelsY.text() or not self.INNumberOfVoxelsZ.text() or not self.INLengthX.text() or not self.INLengthY.text() or not self.INLengthZ.text() or not self.INOriginX.text() or not self.INOriginY.text() or not self.INOriginZ.text():
                self.warning_message('ERROR: Geometry Input Settins Are Missing')
            else:
                # Get the input file location
                if self.INBatchJob.isChecked() and self.INHomogenizationBatchFile.text() != "":
                    global Geometry_filename
                    Geometry_filename = self.INHomogenizationBatchFile.text()
            
                # Get the output vtk location
                vtk_location = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
                
                # Voxelize STL files
                loadingAnimation()
                # Run voxelization executable
                reload(DeveloperInputs)
                if self.UserConfig == "Developer":
                    executable_path = DeveloperInputs.fierro_voxelizer_exe
                elif self.UserConfig == "User":
                    executable_path = "fierro-voxelizer"
                arguments = [Geometry_filename, vtk_location, self.INNumberOfVoxelsX.text(), self.INNumberOfVoxelsY.text(), self.INNumberOfVoxelsZ.text(), self.INOriginX.text(), self.INOriginY.text(), self.INOriginZ.text(), self.INLengthX.text(), self.INLengthY.text(), self.INLengthZ.text()]
                command = [executable_path] + arguments
                try:
                    process = subprocess.Popen(command)
                    process.wait()
                except Exception as e:
                    self.warning_message("ERROR: fierro-voxelizer executable")
                    return
                pvsimple.Delete(self.stl)
                
                # Paraview window
                self.vtk_reader = pvsimple.LegacyVTKReader(FileNames = vtk_location)
                pvsimple.SetDisplayProperties(Representation = "Surface")
                text = self.INPartName.text()
                self.variable_name = f"part_{text}"
                setattr(self, self.variable_name, pvsimple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
                self.display = pvsimple.Show(getattr(self, self.variable_name), self.render_view)
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
                self.TParts.setItem(row, 10, QTableWidgetItem(vtk_location))
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
                self.INRegion.clear()
                self.INRegion_2.clear()
                self.INRegion_3.clear()
                for i in range(self.TParts.rowCount()):
                    self.INRegion.addItem(self.TParts.item(i,0).text())
                    self.INRegion_2.addItem(self.TParts.item(i,0).text())
                    self.INRegion_3.addItem(self.TParts.item(i,0).text())
                for i in range(self.TBasicGeometries.rowCount()):
                    self.INRegion.addItem(self.TBasicGeometries.item(i,0).text())
                self.INRegion.addItem("global")
                self.INRegion_2.addItem("global")
                self.INRegion_3.addItem("global")
                    
                # Leave input page
                self.GeometryOptions.setCurrentIndex(0)
            stopLoadingAnimation()
        self.BVoxelizeGeometry.clicked.connect(voxelize_geometry_click)
        
        # Allow for custom dimensions of stl files
        def custom_dimensions_checked():
            if self.BCustomDimensions.isChecked():
                self.INOriginX.setEnabled(True)
                self.INOriginY.setEnabled(True)
                self.INOriginZ.setEnabled(True)
                self.INLengthX.setEnabled(True)
                self.INLengthY.setEnabled(True)
                self.INLengthZ.setEnabled(True)
                self.LOriginX.setEnabled(True)
                self.LOriginY.setEnabled(True)
                self.LOriginZ.setEnabled(True)
                self.LLengthX.setEnabled(True)
                self.LLengthY.setEnabled(True)
                self.LLengthZ.setEnabled(True)
            else:
                self.INOriginX.setEnabled(False)
                self.INOriginY.setEnabled(False)
                self.INOriginZ.setEnabled(False)
                self.INLengthX.setEnabled(False)
                self.INLengthY.setEnabled(False)
                self.INLengthZ.setEnabled(False)
                self.LOriginX.setEnabled(False)
                self.LOriginY.setEnabled(False)
                self.LOriginZ.setEnabled(False)
                self.LLengthX.setEnabled(False)
                self.LLengthY.setEnabled(False)
                self.LLengthZ.setEnabled(False)
        self.BCustomDimensions.clicked.connect(custom_dimensions_checked)
        self.BStlDimensions.clicked.connect(custom_dimensions_checked)
        
        # Modiffy vtk file
        def modify_vtk():
            # Batch geometry run information
            if self.INBatchJob.isChecked() and self.INHomogenizationBatchFile.text() != "":
                # Rename the file and save it to directory location
                self.new_file_path = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
                Data_filename = self.INHomogenizationBatchFile.text()
                _, file_type = os.path.splitext(Data_filename)
                file_type = file_type[1:]
                print("FIERRO_GUI FILE TYPE IS: ", file_type)
                if Data_filename != self.new_file_path and os.path.isfile(Data_filename) and 'vtk' in file_type:
                    shutil.copy(Data_filename, self.new_file_path)
                # Read what type of vtk file this is
                with open(self.new_file_path, 'r') as file:
                    for index, line in enumerate(file):
                        if index == 3:
                            vtk_type = line.strip()
                            break
                if "STRUCTURED_POINTS" in vtk_type:
                    self.structured_flag = 1
                else:
                    self.structured_flag = 0
                    
            # Fully modify vtk file if it is a structured_points format or the user made changes
            if self.BVTKCustomProperties.isChecked() or self.structured_flag == 1:
                # read lines of original file
                with open(self.new_file_path, 'r') as file:
                    lines = file.readlines()
                    
                # write rectilinear_grid vtk file
                with open(self.new_file_path, "w") as file:
                    # write the header
                    file.write("# vtk DataFile Version 3.0\n")
                    file.write("Header\n")
                    file.write("ASCII\n")
                    file.write("DATASET RECTILINEAR_GRID\n")
                    
                    # number of nodes
                    nodesx = int(self.INvtkvx.text())+1
                    nodesy = int(self.INvtkvy.text())+1
                    nodesz = int(self.INvtkvz.text())+1
#                    nodesx = self.vtkNx+1
#                    nodesy = self.vtkNy+1
#                    nodesz = self.vtkNz+1
                    
                    # spacing of nodes
                    spacex = float(self.INvlx.text())/float(self.INvtkvx.text())
                    spacey = float(self.INvly.text())/float(self.INvtkvy.text())
                    spacez = float(self.INvlz.text())/float(self.INvtkvz.text())
                    
                    # write the node data
                    file.write(f'DIMENSIONS {nodesx} {nodesy} {nodesz}\n')
                    file.write(f'X_COORDINATES {nodesx} float\n')
                    for i in range(0, nodesx):
                        x_pt = i*spacex + float(self.INvox.text())
                        file.write(f'{x_pt} ')
                    file.write("\n")
                    file.write(f'Y_COORDINATES {nodesy} float\n')
                    for i in range(0, nodesy):
                        y_pt = i*spacey + float(self.INvoy.text())
                        file.write(f'{y_pt} ')
                    file.write("\n")
                    file.write(f'Z_COORDINATES {nodesz} float\n')
                    for i in range(0, nodesz):
                        z_pt = i*spacez + float(self.INvoz.text())
                        file.write(f'{z_pt} ')
                    file.write("\n\n")
                    
                    # write the cell data
                    num_elems = int(self.INvtkvx.text())*int(self.INvtkvy.text())*int(self.INvtkvz.text())
                    file.write(f'CELL_DATA {num_elems}\n')
                    file.write("SCALARS density float 1\n")
                    file.write("LOOKUP_TABLE default\n")
                    if self.structured_flag == 1:
                        file.write(' '.join(line.strip() for line in lines[10:]))
                    else:
                        file.writelines(lines[15:])
                    
                # Show new vtk file in paraview
                self.newvar = "part_" + self.INPartName.text()
                pvsimple.Delete(getattr(self, self.newvar))
                self.vtk_reader = pvsimple.LegacyVTKReader(FileNames = self.new_file_path)
                pvsimple.SetDisplayProperties(Representation = "Surface")
                text = self.INPartName.text()
                self.variable_name = f"part_{text}"
                setattr(self, self.variable_name, pvsimple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
                self.display = pvsimple.Show(getattr(self, self.variable_name), self.render_view)
                pvsimple.Hide(self.vtk_reader)
                self.render_view.ResetCamera()
                self.render_view.StillRender()
            
            # Add part to the table
            row = self.TParts.rowCount()
            self.TParts.insertRow(row)
            self.TParts.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
            self.TParts.setItem(row, 1, QTableWidgetItem(self.INvox.text()))
            self.TParts.setItem(row, 2, QTableWidgetItem(self.INvoy.text()))
            self.TParts.setItem(row, 3, QTableWidgetItem(self.INvoz.text()))
            self.TParts.setItem(row, 4, QTableWidgetItem(self.INvlx.text()))
            self.TParts.setItem(row, 5, QTableWidgetItem(self.INvly.text()))
            self.TParts.setItem(row, 6, QTableWidgetItem(self.INvlz.text()))
            self.TParts.setItem(row, 7, QTableWidgetItem(self.INvtkvx.text()))
            self.TParts.setItem(row, 8, QTableWidgetItem(self.INvtkvy.text()))
            self.TParts.setItem(row, 9, QTableWidgetItem(self.INvtkvz.text()))
            self.TParts.setItem(row, 10, QTableWidgetItem(self.new_file_path))
            self.INPartName.clear()
            self.INvox.clear()
            self.INvoy.clear()
            self.INvoz.clear()
            self.INvlx.clear()
            self.INvly.clear()
            self.INvlz.clear()
            self.INPartName.clear()
            
            # Add part as an option for material assignment
            self.INRegion.clear()
            self.INRegion_2.clear()
            self.INRegion_3.clear()
            for i in range(self.TParts.rowCount()):
                self.INRegion.addItem(self.TParts.item(i,0).text())
                self.INRegion_2.addItem(self.TParts.item(i,0).text())
                self.INRegion_3.addItem(self.TParts.item(i,0).text())
            self.INRegion.addItem("global")
            self.INRegion_2.addItem("global")
            self.INRegion_3.addItem("global")

            # Open up window to change color map
            self.SAGeometryScrollArea.verticalScrollBar().setValue(0)
            self.GeometryOptions.setCurrentIndex(3)

            # Get the avaliable arrays for coloring
            self.vtk_reader.UpdatePipeline()
            data_maps = get_paraview_variables(self.vtk_reader)
            self.INSelectColorBy.clear()
            self.INSelectColorBy.addItems(data_maps) # add options to combo box
        self.BAddVTKGeometry.clicked.connect(modify_vtk)
        
        # Allow for custom dimensions of vtk or image stack files
        def custom_dimensions_checked_2():
            if self.BVTKCustomProperties.isChecked():
                self.INvox.setEnabled(True)
                self.INvoy.setEnabled(True)
                self.INvoz.setEnabled(True)
                self.INvlx.setEnabled(True)
                self.INvly.setEnabled(True)
                self.INvlz.setEnabled(True)
                self.Lvox.setEnabled(True)
                self.Lvoy.setEnabled(True)
                self.Lvoz.setEnabled(True)
                self.Lvlx.setEnabled(True)
                self.Lvly.setEnabled(True)
                self.Lvlz.setEnabled(True)
            else:
                self.INvox.setEnabled(False)
                self.INvoy.setEnabled(False)
                self.INvoz.setEnabled(False)
                self.INvlx.setEnabled(False)
                self.INvly.setEnabled(False)
                self.INvlz.setEnabled(False)
                self.Lvox.setEnabled(False)
                self.Lvoy.setEnabled(False)
                self.Lvoz.setEnabled(False)
                self.Lvlx.setEnabled(False)
                self.Lvly.setEnabled(False)
                self.Lvlz.setEnabled(False)
        self.BVTKCustomProperties.clicked.connect(custom_dimensions_checked_2)
        self.BVTKFileProperties.clicked.connect(custom_dimensions_checked_2)
        
        # Delete any previously loaded geometries from table
        global table
        table = self.TParts
        def delete_part():
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
                self.INMaterial.clear()
                for i in range(table.rowCount()):
                    self.INMaterial.addItem(table.item(i,0).text())
                for i in range(self.TBasicGeometries.rowCount()):
                    self.INMaterial.addItem(self.TBasicGeometries.item(i,0).text())
                self.INMaterial.addItem("global")

        # Function to set what table to delete part from
        def set_table(table_selection):
            global table
            table = table_selection
            delete_part()
        self.BDeleteGeometry.clicked.connect(lambda: set_table(self.TParts))
            
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
        
        # Paraview camera settings
        def reset_camera():
            pvsimple.ResetCamera()
        self.BResetCamera.clicked.connect(reset_camera)
            
        # ======= HOMOGENIZATION PIPELINE =======
        Homogenization(self)
        
        # ======== BULK FORMING PIPELINE ========
        Bulk_Forming(self)
            
        # ======= EXPLICIT SOLVER SGH PIPELINE =======
        SGH(self)

    def RunImageToVTK(self):
        print ("RunImageToVTK")
        ImageToVTK(ct_folder_path, 'VTK_Geometry_' + self.INPartName.text(), self.voxelizer_dir, file_type) # Any other args, kwargs are passed to the run function

#    def RunTiffImageToVTK(self):
#        print ("RunTiffImageToVTK")
#        TiffImageToVTK(ct_folder_path, 'VTK_Geometry_' + self.INPartName.text(), self.voxelizer_dir, '.vtk')

#    def RunDream3DReader(self):
#        print ("RunDream3DReader")
#        Dream3DReader(str(dream_filename), 'DREAM_Geometry_' + self.INPartName.text(), self.voxelizer_dir)
        
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
            self.evpfft_dir = os.path.join(self.directory, 'homogenization')
            os.makedirs(self.evpfft_dir, exist_ok=True)
            
            # Create a temp file for the global mesh
            self.global_mesh_dir = os.path.join(self.directory, 'global_mesh')
            os.makedirs(self.global_mesh_dir, exist_ok=True)
            self.GLOBAL_MESH = os.path.join(self.global_mesh_dir, 'global_mesh.yaml')
            self.GLOBAL_MESH_OUTPUT = os.path.join(self.global_mesh_dir, 'mesh.vtk')
            
            # Create temp files for SGH
            sgh_dir = os.path.join(self.directory, 'sgh')
            os.makedirs(sgh_dir, exist_ok=True)
            self.EXPLICIT_SGH_INPUT = os.path.join(sgh_dir, 'explicit_sgh_input.yaml')
            
            # Create temp files for bulk forming
            self.bulk_forming_dir = os.path.join(self.directory, 'bulk_forming')
            os.makedirs(self.bulk_forming_dir, exist_ok=True)
        
        else:
            self.warning_message("ERROR: Working directory was not defined")
            
    # Save file
    def open_save_dialog(self):
        Save(self)
        
    # Load file
    def open_load_dialog(self):
        # Temporairly disconnect functions
        self.INUnits.currentIndexChanged.disconnect()
        # Load GUI
        Load(self)
        # Load visualizations
        if self.TParts.item(0,10):
            self.in_file_path = self.TParts.item(0,10).text()
            self.file_type = os.path.splitext(self.in_file_path)[1].lower()
            Reload_Geometry(self)
        # Load variables [Homogenization Solver]
        if self.INPipelineSelection.currentText() == "Homogenization":
            # Load working directory where job files are saved
            if self.INHomogenizationJobDir.text() != "":
                self.working_directory = self.INHomogenizationJobDir.text()
                self.THomogenization.setEnabled(True)
            # Update material selections
            rows = self.TParts.rowCount()
            if rows > 0:
                self.INRegion_2.clear()
                for i in range(rows):
                    self.INRegion_2.addItem(self.TParts.item(i,0).text())
                self.INRegion_2.addItem("global")
            rows = self.TMaterials.rowCount()
            if rows > 0:
                self.INMaterial.clear()
                for i in range(rows):
                    self.INMaterial.addItem(self.TMaterials.item(i,0).text())
        elif self.INPipelineSelection.currentText() == "Bulk Forming":
            # Load working directory where job files are saved
            if self.INBFJobDir.text() != "":
                self.working_directory = self.INBFJobDir.text()
            # Update material selections
            rows = self.TParts.rowCount()
            if rows > 0:
                self.INRegion_3.clear()
                for i in range(rows):
                    self.INRegion_3.addItem(self.TParts.item(i,0).text())
            rows = self.TMaterials_2.rowCount()
            if rows > 0:
                self.INMaterial_2.clear()
                for i in range(rows):
                    self.INMaterial_2.addItem(self.TMaterials_2.item(i,0).text())
            
        # Reconnect functions
        self.INUnits.currentIndexChanged.connect(self.units)
        self.old_units = self.INUnits.currentText()
        self.new_units = 'new'
    
    # Output which user profile the user selected
    def handleSettingChanged(self, setting):
        if setting == "User":
            self.UserConfig = "User"
        elif setting == "Developer":
            self.UserConfig = "Developer"





#### ZACK'S ADDITIONAL FUNCTIONS
## Convert tiff image stack to .stl
#def image_to_STL():
#    loadingAnimation()
#    print(ct_folder_path)
#    print ("TIFF")
#    # Call TiffImageToVTK function in python file of same name, specifying to save only .stl output
#    TiffImageToVTK(ct_folder_path, 'STL_Geometry_' + self.INPartName.text(), self.voxelizer_dir, '.stl')
#
#    # Paraview window
#    stl_location = self.voxelizer_dir + '/STL_Geometry_' + str(self.INPartName.text()) + '.stl'
#    # Show the stl part
#    self.stl = pvsimple.STLReader(FileNames = stl_location)
#    text = self.INPartName.text()
#    self.variable_name = f"part_{text}"
#    print ("var name: " + self.variable_name)
#    setattr(self, self.variable_name, "")
#    pvsimple.Show(self.stl, self.render_view)
#    pvsimple.ResetCamera(view=None)
#
#    # Get the bounds of the stl
#    self.stl.UpdatePipeline()
#    self.bounds = self.stl.GetDataInformation().GetBounds()
#    self.stlLx = round(self.bounds[1] - self.bounds[0],4)
#    self.stlLy = round(self.bounds[3] - self.bounds[2],4)
#    self.stlLz = round(self.bounds[5] - self.bounds[4],4)
#    self.INLengthX.setText(str(self.stlLx))
#    self.INLengthY.setText(str(self.stlLy))
#    self.INLengthZ.setText(str(self.stlLz))
#
#    # Add part to the table
#    row = self.TParts.rowCount()
#    self.TParts.insertRow(row)
#    self.TParts.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
#    self.TParts.setItem(row, 1, QTableWidgetItem(self.INOriginX.text()))
#    self.TParts.setItem(row, 2, QTableWidgetItem(self.INOriginY.text()))
#    self.TParts.setItem(row, 3, QTableWidgetItem(self.INOriginZ.text()))
#    self.TParts.setItem(row, 4, QTableWidgetItem(self.INLengthX.text()))
#    self.TParts.setItem(row, 5, QTableWidgetItem(self.INLengthY.text()))
#    self.TParts.setItem(row, 6, QTableWidgetItem(self.INLengthZ.text()))
#    self.TParts.setItem(row, 7, QTableWidgetItem(self.INNumberOfVoxelsX.text()))
#    self.TParts.setItem(row, 8, QTableWidgetItem(self.INNumberOfVoxelsY.text()))
#    self.TParts.setItem(row, 9, QTableWidgetItem(self.INNumberOfVoxelsZ.text()))
#    self.INPartName.clear()
#    self.INOriginX.clear()
#    self.INOriginY.clear()
#    self.INOriginZ.clear()
#    self.INLengthX.clear()
#    self.INLengthY.clear()
#    self.INLengthZ.clear()
#    self.INNumberOfVoxelsX.clear()
#    self.INNumberOfVoxelsY.clear()
#    self.INNumberOfVoxelsZ.clear()
#
#    # Add part as an option for material assignment
#    self.INPartMaterial.clear()
#    self.INPartMaterial.addItem("global")
#    for i in range(self.TParts.rowCount()):
#        self.INPartMaterial.addItem(self.TParts.item(i,0).text())
#    for i in range(self.TBasicGeometries.rowCount()):
#        self.INPartMaterial.addItem(self.TBasicGeometries.item(i,0).text())
#    stopLoadingAnimation()
#self.BTiffToStl.clicked.connect(image_to_STL)
#

# Upload polycrystals
#            elif self.INSelectGeometryImport.currentIndex() == 3:
#                global dream_filename
#                dream_filename = QFileDialog.getOpenFileName(
#                    filter="Dream File (*.xdmf* *.dream3d)",
#                )
#                loadingAnimation()
#
#                # Get xdmf file from dream3d file
#                ext = os.path.splitext(dream_filename[0])
#                if (ext[1] == ".dream3d"):
#                    print ("File to write: ", dream_filename[0])
#                    dream_filename = dream_filename[0]
#                    #worker = Worker(self.RunDream3DReader) # Any other args, kwargs are passed to the run function
#
#                    # Execute
#                    #self.threadpool.start(worker)
#
#                    Dream3DReader(str(dream_filename), 'DREAM_Geometry_' + self.INPartName.text(), self.voxelizer_dir)
#                    dream_filename = self.voxelizer_dir + '/DREAM_Geometry_' + str(self.INPartName.text()) + '.xdmf'
#
#                # Paraview window
#                #stl_location = self.voxelizer_dir + '/DREAM_Geometry_' + str(self.INPartName.text()) + '.stl'
#                print (dream_filename)
#                # Show the stl part
#                self.dream = pvsimple.XDMFReader(FileNames = dream_filename)
#                if (self.dream):
#                    print ("Success")
#                else:
#                    print ("Failure")
#                pvsimple.SetDisplayProperties(Representation = "Surface")
#                text = self.INPartName.text()
#                self.variable_name = f"part_{text}"
#                print ("var name: " + self.variable_name)
#                setattr(self, self.variable_name, "")
#                #pvsimple.Show(getattr(self, self.variable_name), self.render_view)
#                display = pvsimple.Show(self.dream, self.render_view)
#                pvsimple.ColorBy(display, 'EulerAngles')
#                self.INSelectColorBy.setEnabled(True)
#                pvsimple.ResetCamera(view=None)
#                row = self.TDream.rowCount()
#                self.TDream.insertRow(row)
#                self.TDream.setItem(row, 0, QTableWidgetItem(self.INPartName.text()))
#                self.INPartName.clear()
#                stopLoadingAnimation()
#                def change_color_by():
#                    # Dream Color By Selection
#                    #print (self.INSelectColorBy.itemText(self.INSelectColorBy.currentIndex()))
#                    color = None
#                    match self.INSelectColorBy.itemText(self.INSelectColorBy.currentIndex()):
#                        case "Euler Angles":
#                            color = 'EulerAngles'
#                        case "Feature Ids":
#                            color = 'FeatureIds'
#                        case "IPF Color":
#                            color = 'IPFColor'
#                        case "Phases":
#                            color = 'Phases'
#
#                    pvsimple.ColorBy(display, color)
#                    pvsimple.ResetCamera(view=None)
#                self.INSelectColorBy.currentIndexChanged.connect(change_color_by)
#            else:
#                self.warning_message('ERROR: Incorrect file type')
