# To import images using resource file you must conver the .rc file to .py using (in command line):
# pyside6-rcc IconResourceFile.qrc -o IconResourceFile.py

from fierro_gui.ui_FIERRO_GUI import Ui_MainWindow

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QAction, QBrush, QColor, QConicalGradient,
    QCursor, QFont, QFontDatabase, QGradient,
    QIcon, QImage, QKeySequence, QLinearGradient,
    QPainter, QPalette, QPixmap, QRadialGradient,
    QTransform)
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox, QFormLayout,
    QFrame, QGridLayout, QHBoxLayout, QHeaderView,
    QLabel, QLineEdit, QMainWindow, QMenu,
    QMenuBar, QPlainTextEdit, QProgressBar, QPushButton,
    QSizePolicy, QSpacerItem, QSplitter, QStackedWidget,
    QStatusBar, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget)
    
import os, sys
import numpy as np
from PySide6.QtWidgets import (QFileDialog, QMessageBox)
from PySide6.QtCore import (QTimer, QProcess)
import matplotlib
matplotlib.use('qt5agg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
import csv
import re
import shutil
import paraview.simple as pvsimple
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PySide6.QtCore import QUrl
from PySide6.QtGui import QDesktopServices
#import fierro_voxelizer
import fierro_mesh_builder
import tempfile
import time
import subprocess

from Explicit_SGH import *
from EVPFFT_Lattice import *
from Mesh_Builder_WInput import *
from FIERRO_Setup import *
from DeveloperInputs import *

class FIERRO_GUI(Ui_MainWindow):
    def setupUi(self, MainWindow):
        super().setupUi(MainWindow)

        # Paraview imports
        self.render_view = pvsimple.CreateRenderView()
        self.Paraview = QVTKRenderWindowInteractor(rw=self.render_view.GetRenderWindow(),iren=self.render_view.GetInteractor())
        self.verticalLayout_20.addWidget(self.Paraview)
        
        # Open Url in Help Menu
        def openUrl():
            url = QUrl('https://lanl.github.io/Fierro/')
            if not QDesktopServices.openUrl(url):
                QMessageBox.warning(self, 'Open Url', 'Could not open url')
        
        # Help menu
        self.actionManual.triggered.connect(openUrl)
        
        # Working Directory Menu
        self.actionChange_Working_Directory.triggered.connect(self.open_fierro_setup_dialog)
        
        # Upload Geometry
        def geometry_upload_click():
            if not self.INPartName.text():
                self.warning_message("Please name the part")
            else:
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
                    self.STLVoxelization.setEnabled(True)
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
                    self.STLVoxelization.setEnabled(False)
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
                else:
                    self.warning_message('ERROR: Incorrect file type')
        self.BUploadGeometryFile.clicked.connect(geometry_upload_click)
        
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
                # Run voxelization executable
                executable_path = fierro_voxelizer_exe
                vtk_location = self.voxelizer_dir + '/VTK_Geometry_' + str(self.INPartName.text()) + '.vtk'
                arguments = [b3_filename[0], vtk_location, self.INNumberOfVoxelsX.text(), self.INNumberOfVoxelsY.text(), self.INNumberOfVoxelsZ.text(), self.INOriginX.text(), self.INOriginY.text(), self.INOriginZ.text(), self.INLengthX.text(), self.INLengthY.text(), self.INLengthZ.text()]
                command = [executable_path] + arguments
                process = subprocess.Popen(command)
                process.wait()
                    
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
                
        self.BVoxelizeGeometry.clicked.connect(voxelize_geometry_click)
        
        # Delete any previously loaded geometries from table
        def delete_part():
            current_row = self.TParts.currentRow()
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
                self.newvar = "part_" + self.TParts.item(current_row,0).text()
                pvsimple.Delete(getattr(self, self.newvar))
                self.render_view.ResetCamera()
                self.render_view.StillRender()
                self.TParts.removeRow(current_row)
                
                # delete from material assignment options
                self.INPartMaterial.clear()
                self.INPartMaterial.addItem("global")
                for i in range(self.TParts.rowCount()):
                    self.INPartMaterial.addItem(self.TParts.item(i,0).text())
                for i in range(self.TBasicGeometries.rowCount()):
                    self.INPartMaterial.addItem(self.TBasicGeometries.item(i,0).text())
        self.BDeleteGeometry.clicked.connect(delete_part)
            
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
            
        # ======= EVPFFT SOLVER LATTICE PIPELINE =======
        EVPFFT_Lattice(self)
            
        # ======= EXPLICIT SOLVER SGH PIPELINE =======
        Explicit_SGH(self)
        
    # ========== WARNING MESSAGE ============
    def warning_message(self, msg):
        message = QMessageBox()
        message.setText(msg)
        message.exec()
        
    # ========== FIERRO SETUP ============
    def open_fierro_setup_dialog(self, gself):
        dialog = FierroSetup(gself)
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
    
    # ========== CHECK FOR BUILD PACKAGES ============
    def check_build_packages(self, executable_path, arguments):
        try:
            process = subprocess.Popen([executable_path] + arguments)
            process.wait()
            return True
        except Exception as e:
            print("Error occurred while running the process:", e)
            return False

