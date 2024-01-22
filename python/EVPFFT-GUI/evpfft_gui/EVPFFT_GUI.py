# To import images using resource file you must conver the .rc file to .py using (in command line):
# pyside6-rcc IconResourceFile.qrc -o IconResourceFile.py

from evpfft_gui.ui_EVPFFT_GUI import Ui_MainWindow

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
import fierro_voxelizer
import tempfile
import time

class LocalResource:
    FILE_PATH = os.path.abspath(
        os.path.join(*(os.path.split(os.path.expanduser(__file__))[:-1]))
    )

    @staticmethod
    def get_resource_name(relpath: str) -> str:
        return os.path.join(LocalResource.FILE_PATH, relpath)

VTK_OUTPUT = os.path.join(tempfile.gettempdir(), 'VTK_Geometry.vtk')
ELASTIC_PARAMETERS_0 = 'elastic_parameters_0.txt'
ELASTIC_PARAMETERS_1 = 'elastic_parameters_1.txt'
PLASTIC_PARAMETERS = LocalResource.get_resource_name('plastic_parameters.txt')
EVPFFT_INPUT = os.path.join(tempfile.gettempdir(), 'evpfft_lattice_input.txt')

class EVPFFT_GUI(Ui_MainWindow):
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
                
        
        # BUTTON SETUP
        # Connect tab buttons to settings windows
        self.BImportPart.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(1))
        self.BDefineMaterial.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(2))
        self.BApplyBC.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(3))
        self.BSolverSettings.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(4))
        self.BViewResults.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(5))
        self.VoxelResolution = (1., 1., 1.)
        
        # Help menu
        self.actionEVPFFT_Manual.triggered.connect(openUrl)
        
        # Upload Geometry
        def geometry_upload_click():
            try:
                self.stl
            except:
                print('')
            else:
                pvsimple.Delete(self.stl)
                
            try:
                self.voxel_reader
            except:
                print('')
            else:
                pvsimple.Delete(self.threshold)
                
            global b3_filename
            b3_filename = QFileDialog.getOpenFileName(
                filter="Geometry File (*.stl *.vtk)",
            )
            # Paraview window
            self.file_type = b3_filename[0][-4:-1]
            if self.file_type == '.st':
                self.stl = pvsimple.STLReader(FileNames = b3_filename)
                self.STLVoxelization.setEnabled(True)
                self.LNumberOfVoxelsX.setEnabled(True)
                self.INNumberOfVoxelsX.setEnabled(True)
                self.LNumberOfVoxelsY.setEnabled(True)
                self.INNumberOfVoxelsY.setEnabled(True)
                self.LNumberOfVoxelsZ.setEnabled(True)
                self.INNumberOfVoxelsZ.setEnabled(True)
                self.BVoxelizeGeometry.setEnabled(True)
            elif self.file_type == '.vt':
                self.stl = pvsimple.LegacyVTKReader(FileNames = b3_filename)
                pvsimple.SetDisplayProperties(Representation = "Surface")
                shutil.copy(b3_filename[0], VTK_OUTPUT)
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
            else:
                warning_message('ERROR: Incorrect file type')
            pvsimple.Show(self.stl, self.render_view)
            pvsimple.ResetCamera(view=None)
        self.BUploadGeometryFile.clicked.connect(geometry_upload_click)
        
        # Voxelize Geometry
        def voxelize_geometry_click():
            if not self.INNumberOfVoxelsX.text() or not self.INNumberOfVoxelsY.text() or not self.INNumberOfVoxelsZ.text():
                warning_message('ERROR: Number of voxels NOT defined')
            else:
                try:
                    self.voxel_reader
                except:
                    print('')
                else:
                    pvsimple.Delete(self.threshold)
                    
                self.VoxelResolution = fierro_voxelizer.create_voxel_vtk(
                    b3_filename[0],
                    VTK_OUTPUT,
                    int(self.INNumberOfVoxelsX.text()),
                    int(self.INNumberOfVoxelsY.text()),
                    int(self.INNumberOfVoxelsZ.text()),
                )
                # Paraview window
                pvsimple.Delete(self.stl)
                self.voxel_reader = pvsimple.LegacyVTKReader(FileNames = VTK_OUTPUT)
                pvsimple.SetDisplayProperties(Representation = "Surface")
                self.threshold = pvsimple.Threshold(Input = self.voxel_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
                pvsimple.Show(self.threshold, self.render_view)
                pvsimple.Hide(self.voxel_reader)
                self.render_view.ResetCamera()
                self.render_view.StillRender()
        self.BVoxelizeGeometry.clicked.connect(voxelize_geometry_click)
        
        # Apply Material
        def material_type():
            if str(self.INSolidGas.currentText()) == 'Solid':
                self.INMaterialType.clear()
                self.INMaterialType.addItem(QCoreApplication.translate("MainWindow", u"Isotropic", None))
                self.INMaterialType.addItem(QCoreApplication.translate("MainWindow", u"Transversely Isotropic", None))
                self.INMaterialType.addItem(QCoreApplication.translate("MainWindow", u"Orthotropic", None))
                self.INMaterialType.addItem(QCoreApplication.translate("MainWindow", u"Anisotropic", None))
                self.MaterialTypeTool.setEnabled(True)
            if str(self.INSolidGas.currentText()) == 'Gas':
                self.INMaterialType.clear()
                self.INMaterialType.addItem(QCoreApplication.translate("MainWindow", u"Ideal Gas", None))
                self.MaterialTypeTool.setEnabled(False)
        self.INSolidGas.currentIndexChanged.connect(material_type)
        
        def material_region():
            if str(self.INRegion.currentText()) == 'Void':
                pvsimple.Hide(self.threshold)
                self.threshold = pvsimple.Threshold(Input = self.voxel_reader, Scalars = "density", ThresholdMethod = "Below Lower Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
                pvsimple.Show(self.threshold, self.render_view)
                self.render_view.ResetCamera()
                self.render_view.StillRender()
            else:
                pvsimple.Hide(self.threshold)
                self.threshold = pvsimple.Threshold(Input = self.voxel_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
                pvsimple.Show(self.threshold, self.render_view)
                self.render_view.ResetCamera()
                self.render_view.StillRender()
        self.INRegion.currentIndexChanged.connect(material_region)
        
        def material_class():
            if str(self.INMaterialType.currentText()) == 'Isotropic':
                self.MaterialTypeTool.setCurrentIndex(0)
            if str(self.INMaterialType.currentText()) == 'Transversely Isotropic':
                self.MaterialTypeTool.setCurrentIndex(1)
            if str(self.INMaterialType.currentText()) == 'Orthotropic':
                self.MaterialTypeTool.setCurrentIndex(3)
            if str(self.INMaterialType.currentText()) == 'Anisotropic':
                self.MaterialTypeTool.setCurrentIndex(2)
        self.INMaterialType.currentIndexChanged.connect(material_class)
        
        def add_material():
            warning_flag = 0
            for i in range(self.TMaterials.rowCount()):
                if str(self.INRegion.currentText()) == self.TMaterials.item(i,1).text():
                    warning_message('ERROR: There is already a material assigned to this region')
                    warning_flag = 1
            
            if str(self.INSolidGas.currentText()) == 'Gas':
                if not self.INMaterialName.text():
                    warning_message('ERROR: Material definition incomplete')
                    warning_flag = 1
                else:
                    # Fill out material definition table
                    row = self.TMaterials.rowCount()
                    self.TMaterials.insertRow(row)
                    self.TMaterials.setItem(row, 0, QTableWidgetItem(
                        self.INMaterialName.text().strip())
                    )
                    self.TMaterials.setItem(row, 1, QTableWidgetItem(
                        str(self.INRegion.currentText()))
                    )
                    self.TMaterials.setItem(row, 2, QTableWidgetItem(
                        str(self.INMaterialType.currentText()))
                    )
                    self.INMaterialName.clear()
            else:
                if str(self.INMaterialType.currentText()) == 'Isotropic':
                    if not self.INYoungsModulus.text() or not self.INPoissonsRatio.text() or not self.INMaterialName.text():
                        warning_message('ERROR: Material definition incomplete')
                        warning_flag = 1
                    else:
                        # Calculate Stiffness Matrix
                        E = float(self.INYoungsModulus.text())
                        NU = float(self.INPoissonsRatio.text())
                        INCalcG = float(self.INYoungsModulus.text())/(2*(1+float(self.INPoissonsRatio.text())))
                        S11 = 1/E
                        S12 = -NU/E
                        S13 = -NU/E
                        S22 = 1/E
                        S23 = -NU/E
                        S33 = 1/E
                        S44 = 1/INCalcG
                        S55 = 1/INCalcG
                        S66 = 1/INCalcG
                        S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                        C11 = (1/S)*(S22*S33-S23*S23)
                        C12 = (1/S)*(S13*S23-S12*S33)
                        C13 = (1/S)*(S12*S23-S13*S22)
                        C22 = (1/S)*(S33*S11-S13*S13)
                        C23 = (1/S)*(S12*S13-S23*S11)
                        C33 = (1/S)*(S11*S22-S12*S12)
                        C44 = (1/S44)
                        C55 = (1/S55)
                        C66 = (1/S66)
                    
                        self.INYoungsModulus.clear()
                        self.INPoissonsRatio.clear()
                if str(self.INMaterialType.currentText()) == 'Transversely Isotropic':
                    if not self.INEip.text() or not self.INNUip.text() or not self.INEop.text() or not self.INNUop.text() or not self.INGop.text() or not self.INMaterialName.text():
                        warning_message('ERROR: Material definition incomplete')
                        warning_flag = 1
                    else:
                        if str(self.INIsotropicPlane.currentText()) == 'x-y plane':
                            # Calculate Stiffness Matrix
                            NUip = float(self.INNUip.text())
                            NUop = float(self.INNUop.text())
                            Eip = float(self.INEip.text())
                            Eop = float(self.INEop.text())
                            Gop = float(self.INGop.text())
                            S11 = 1/Eip
                            S12 = -NUip/Eip
                            S13 = -NUop/Eop
                            S22 = 1/Eip
                            S23 = -NUop/Eop
                            S33 = 1/Eop
                            S44 = 1/Gop
                            S55 = 1/Gop
                            S66 = 1/(Eip/(2*(1+NUip)))
                            S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                            C11 = (1/S)*(S22*S33-S23*S23)
                            C12 = (1/S)*(S13*S23-S12*S33)
                            C22 = (1/S)*(S33*S11-S13*S13)
                            C13 = (1/S)*(S12*S23-S13*S22)
                            C33 = (1/S)*(S11*S22-S12*S12)
                            C23 = (1/S)*(S12*S13-S23*S11)
                            C44 = (1/S44)
                            C55 = (1/S55)
                            C66 = (1/S66)
                            
                            self.INEip.clear()
                            self.INNUip.clear()
                            self.INEop.clear()
                            self.INNUop.clear()
                            self.INGop.clear()
                        if str(self.INIsotropicPlane.currentText()) == 'x-z plane':
                            # Calculate Stiffness Matrix
                            NUip = float(self.INNUip.text())
                            NUop = float(self.INNUop.text())
                            Eip = float(self.INEip.text())
                            Eop = float(self.INEop.text())
                            Gop = float(self.INGop.text())
                            S11 = 1/Eip
                            S12 = -NUop/Eop
                            S13 = -NUip/Eip
                            S22 = 1/Eop
                            S23 = -NUop/Eop
                            S33 = 1/Eip
                            S44 = 1/Gop
                            S55 = 1/(Eip/(2*(1+NUip)))
                            S66 = 1/Gop
                            S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                            C11 = (1/S)*(S22*S33-S23*S23)
                            C12 = (1/S)*(S13*S23-S12*S33)
                            C22 = (1/S)*(S33*S11-S13*S13)
                            C13 = (1/S)*(S12*S23-S13*S22)
                            C33 = (1/S)*(S11*S22-S12*S12)
                            C23 = (1/S)*(S12*S13-S23*S11)
                            C44 = (1/S44)
                            C55 = (1/S55)
                            C66 = (1/S66)
                            
                            self.INEip.clear()
                            self.INNUip.clear()
                            self.INEop.clear()
                            self.INNUop.clear()
                            self.INGop.clear()
                        if str(self.INIsotropicPlane.currentText()) == 'y-z plane':
                            # Calculate Stiffness Matrix
                            NUip = float(self.INNUip.text())
                            NUop = float(self.INNUop.text())
                            Eip = float(self.INEip.text())
                            Eop = float(self.INEop.text())
                            Gop = float(self.INGop.text())
                            S11 = 1/Eop
                            S12 = -NUop/Eop
                            S13 = -NUop/Eop
                            S22 = 1/Eip
                            S23 = -NUip/Eip
                            S33 = 1/Eip
                            S44 = 1/(Eip/(2*(1+NUip)))
                            S55 = 1/Gop
                            S66 = 1/Gop
                            S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                            C11 = (1/S)*(S22*S33-S23*S23)
                            C12 = (1/S)*(S13*S23-S12*S33)
                            C22 = (1/S)*(S33*S11-S13*S13)
                            C13 = (1/S)*(S12*S23-S13*S22)
                            C33 = (1/S)*(S11*S22-S12*S12)
                            C23 = (1/S)*(S12*S13-S23*S11)
                            C44 = (1/S44)
                            C55 = (1/S55)
                            C66 = (1/S66)
                            
                            self.INEip.clear()
                            self.INNUip.clear()
                            self.INEop.clear()
                            self.INNUop.clear()
                            self.INGop.clear()
        
                if str(self.INMaterialType.currentText()) == 'Orthotropic':
                    if not self.INEx.text() or not self.INEy.text() or not self.INEz.text() or not self.INNUxy.text() or not self.INNUxz.text() or not self.INNUyz.text() or not self.INGxy.text() or not self.INGxz.text() or not self.INGyz.text() or not self.INMaterialName.text():
                        warning_message('ERROR: Material definition incomplete')
                        warning_flag = 1
                    else:
                        # Calculate Stiffness Matrix
                        S11 = 1/float(self.INEx.text())
                        S12 = -float(self.INNUxy.text())/float(self.INEx.text())
                        S13 = -float(self.INNUxz.text())/float(self.INEx.text())
                        S22 = 1/float(self.INEy.text())
                        S23 = -float(self.INNUyz.text())/float(self.INEy.text())
                        S33 = 1/float(self.INEz.text())
                        S44 = 1/float(self.INGyz.text())
                        S55 = 1/float(self.INGxz.text())
                        S66 = 1/float(self.INGxy.text())
                        S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                        C11 = (1/S)*(S22*S33-S23*S23)
                        C12 = (1/S)*(S13*S23-S12*S33)
                        C22 = (1/S)*(S33*S11-S13*S13)
                        C13 = (1/S)*(S12*S23-S13*S22)
                        C33 = (1/S)*(S11*S22-S12*S12)
                        C23 = (1/S)*(S12*S13-S23*S11)
                        C44 = (1/S44)
                        C55 = (1/S55)
                        C66 = (1/S66)
                        
                        self.INEx.clear()
                        self.INEy.clear()
                        self.INEz.clear()
                        self.INNUxy.clear()
                        self.INNUxz.clear()
                        self.INNUyz.clear()
                        self.INGxy.clear()
                        self.INGxz.clear()
                        self.INGyz.clear()
                        
                if str(self.INMaterialType.currentText()) == 'Anisotropic':
                    if not self.TAnisotropic.item(0,0).text() or not self.TAnisotropic.item(0,1).text() or not self.TAnisotropic.item(0,2).text() or not self.TAnisotropic.item(0,3).text() or not self.TAnisotropic.item(0,4).text() or not self.TAnisotropic.item(0,5).text() or not self.TAnisotropic.item(1,1).text() or not self.TAnisotropic.item(1,2).text() or not self.TAnisotropic.item(1,3).text() or not self.TAnisotropic.item(1,4).text() or not self.TAnisotropic.item(1,5).text() or not self.TAnisotropic.item(2,2).text() or not self.TAnisotropic.item(2,3).text() or not self.TAnisotropic.item(2,4).text() or not self.TAnisotropic.item(2,5).text() or not self.TAnisotropic.item(3,3).text() or not self.TAnisotropic.item(3,4).text() or not self.TAnisotropic.item(3,5).text() or not self.TAnisotropic.item(4,4).text() or not self.TAnisotropic.item(4,5).text() or not self.TAnisotropic.item(5,5).text() or not self.INMaterialName.text():
                        warning_message('ERROR: Material definition incomplete')
                        warning_flag = 1
                    else:
                        # Fill out material definition table
                        row = self.TMaterials.rowCount()
                        self.TMaterials.insertRow(row)
                        self.TMaterials.setItem(row, 0, QTableWidgetItem(
                            self.INMaterialName.text().strip())
                        )
                        self.TMaterials.setItem(row, 1, QTableWidgetItem(
                            str(self.INMaterialType.currentText()))
                        )
                        k = 2
                        for i in [0,1,2,3,4,5,6]:
                            for j in range(i,6):
                                self.TMaterials.setItem(
                                    row, k, QTableWidgetItem(self.TAnisotropic.item(i,j).text())
                                )
                                self.TAnisotropic.item(i,j).setText('')
                                k += 1
                        self.INMaterialName.clear()
                        warning_flag = 1
                if warning_flag == 0:
                    # Fill out material definition table
                    row = self.TMaterials.rowCount()
                    self.TMaterials.insertRow(row)
                    self.TMaterials.setItem(row, 0, QTableWidgetItem(
                        self.INMaterialName.text().strip())
                    )
                    self.TMaterials.setItem(row, 1, QTableWidgetItem(
                        str(self.INRegion.currentText()))
                    )
                    if str(self.INMaterialType.currentText()) == 'Transversely Isotropic':
                        self.TMaterials.setItem(row, 2, QTableWidgetItem(
                            str(self.INMaterialType.currentText() + ' ' + self.INIsotropicPlane.currentText()))
                        )
                    else:
                        self.TMaterials.setItem(row, 2, QTableWidgetItem(
                            str(self.INMaterialType.currentText()))
                        )
                    self.TMaterials.setItem(
                        row, 3, QTableWidgetItem(str(C11))
                    )
                    self.TMaterials.setItem(
                        row, 4, QTableWidgetItem(str(C12))
                    )
                    self.TMaterials.setItem(
                        row, 5, QTableWidgetItem(str(C13))
                    )
                    self.TMaterials.setItem(
                        row, 9, QTableWidgetItem(str(C22))
                    )
                    self.TMaterials.setItem(
                        row, 10, QTableWidgetItem(str(C23))
                    )
                    self.TMaterials.setItem(
                        row, 14, QTableWidgetItem(str(C33))
                    )
                    self.TMaterials.setItem(
                        row, 18, QTableWidgetItem(str(C44))
                    )
                    self.TMaterials.setItem(
                        row, 21, QTableWidgetItem(str(C55))
                    )
                    self.TMaterials.setItem(
                        row, 23, QTableWidgetItem(str(C66))
                    )
                    for i in [6,7,8,11,12,13,15,16,17,19,20,22]:
                       self.TMaterials.setItem(row, i, QTableWidgetItem('0'))
                    self.INMaterialName.clear()
                else:
                    warning_flag = 0
                
        def delete_material():
            current_row = self.TMaterials.currentRow()
            if current_row < 0:
                return QMessageBox.warning(QMessageBox(),"Warning","Please select a record to delete")

            button = QMessageBox.question(
                QMessageBox(),
                'Confirmation',
                'Are you sure that you want to delete the selected row?',
                QMessageBox.Yes |
                QMessageBox.No
            )
            if button == QMessageBox.StandardButton.Yes:
                self.TMaterials.removeRow(current_row)
                
        def regenerate_elastic_constants():
            current_row = self.TMaterials.currentRow()
            if current_row < 0:
                return QMessageBox.warning(QMessageBox(),"Warning","Please select a material from the table")
                
            # Define Stiffness Matrix
            Mstiffness = [[float(self.TMaterials.item(current_row,3).text()), float(self.TMaterials.item(current_row,4).text()), float(self.TMaterials.item(current_row,5).text()),  float(self.TMaterials.item(current_row,6).text()), float(self.TMaterials.item(current_row,7).text()), float(self.TMaterials.item(current_row,8).text())], [float(self.TMaterials.item(current_row,4).text()), float(self.TMaterials.item(current_row,9).text()), float(self.TMaterials.item(current_row,10).text()),  float(self.TMaterials.item(current_row,11).text()), float(self.TMaterials.item(current_row,12).text()), float(self.TMaterials.item(current_row,13).text())], [float(self.TMaterials.item(current_row,5).text()), float(self.TMaterials.item(current_row,10).text()), float(self.TMaterials.item(current_row,14).text()), float(self.TMaterials.item(current_row,15).text()), float(self.TMaterials.item(current_row,16).text()), float(self.TMaterials.item(current_row,17).text())], [float(self.TMaterials.item(current_row,6).text()), float(self.TMaterials.item(current_row,11).text()), float(self.TMaterials.item(current_row,15).text()), float(self.TMaterials.item(current_row,18).text()), float(self.TMaterials.item(current_row,19).text()), float(self.TMaterials.item(current_row,20).text())], [float(self.TMaterials.item(current_row,7).text()), float(self.TMaterials.item(current_row,12).text()), float(self.TMaterials.item(current_row,16).text()), float(self.TMaterials.item(current_row,19).text()), float(self.TMaterials.item(current_row,21).text()), float(self.TMaterials.item(current_row,22).text())], [float(self.TMaterials.item(current_row,8).text()), float(self.TMaterials.item(current_row,13).text()), float(self.TMaterials.item(current_row,17).text()), float(self.TMaterials.item(current_row,20).text()), float(self.TMaterials.item(current_row,22).text()), float(self.TMaterials.item(current_row,23).text())]]
            Mcompliance = np.linalg.inv(Mstiffness)
            if self.TMaterials.item(current_row,2).text() == 'Isotropic':
                self.MaterialTypeTool.setCurrentIndex(0)
                self.INMaterialType.setCurrentIndex(0)
                self.INMaterialName.clear()
                self.INYoungsModulus.clear()
                self.INPoissonsRatio.clear()
                E = 1/Mcompliance[0][0]
                nu = -Mcompliance[0][1]*E
                self.INMaterialName.insert(self.TMaterials.item(current_row,0).text())
                self.INYoungsModulus.insert(str(E))
                self.INPoissonsRatio.insert(str(nu))
            elif 'Transversely Isotropic' in self.TMaterials.item(current_row,2).text():
                self.MaterialTypeTool.setCurrentIndex(1)
                self.INMaterialType.setCurrentIndex(1)
                self.INMaterialName.clear()
                self.INEip.clear()
                self.INNUip.clear()
                self.INEop.clear()
                self.INNUop.clear()
                self.INGop.clear()
                if 'x-y plane' in self.TMaterials.item(current_row,2).text():
                    Eip = 1/Mcompliance[0][0]
                    nuip = -Mcompliance[0][1]*Eip
                    Eop = 1/Mcompliance[2][2]
                    nuop = -Mcompliance[0][2]*Eop
                    Gop  = 1/Mcompliance[3][3]
                    self.INMaterialName.insert(self.TMaterials.item(current_row,0).text())
                    self.INEip.insert(str(Eip))
                    self.INNUip.insert(str(nuip))
                    self.INEop.insert(str(Eop))
                    self.INNUop.insert(str(nuop))
                    self.INGop.insert(str(Gop))
                    self.INIsotropicPlane.setCurrentIndex(0)
                elif 'x-z plane' in self.TMaterials.item(current_row,2).text():
                    Eip = 1/Mcompliance[0][0]
                    nuip = -Mcompliance[0][2]*Eip
                    Eop = 1/Mcompliance[1][1]
                    nuop = -Mcompliance[0][1]*Eop
                    Gop  = 1/Mcompliance[3][3]
                    self.INMaterialName.insert(self.TMaterials.item(current_row,0).text())
                    self.INEip.insert(str(Eip))
                    self.INNUip.insert(str(nuip))
                    self.INEop.insert(str(Eop))
                    self.INNUop.insert(str(nuop))
                    self.INGop.insert(str(Gop))
                    self.INIsotropicPlane.setCurrentIndex(1)
                elif 'y-z plane' in self.TMaterials.item(current_row,2).text():
                    Eip = 1/Mcompliance[1][1]
                    nuip = -Mcompliance[1][2]*Eip
                    Eop = 1/Mcompliance[0][0]
                    nuop = -Mcompliance[0][1]*Eop
                    Gop  = 1/Mcompliance[4][4]
                    self.INMaterialName.insert(self.TMaterials.item(current_row,0).text())
                    self.INEip.insert(str(Eip))
                    self.INNUip.insert(str(nuip))
                    self.INEop.insert(str(Eop))
                    self.INNUop.insert(str(nuop))
                    self.INGop.insert(str(Gop))
                    self.INIsotropicPlane.setCurrentIndex(2)
            elif self.TMaterials.item(current_row,2).text() == 'Orthotropic':
                self.MaterialTypeTool.setCurrentIndex(3)
                self.INMaterialType.setCurrentIndex(2)
                self.INMaterialName.clear()
                self.INEx.clear()
                self.INEy.clear()
                self.INEz.clear()
                self.INNUxy.clear()
                self.INNUxz.clear()
                self.INNUyz.clear()
                self.INGxy.clear()
                self.INGxz.clear()
                self.INGyz.clear()
                Ex = 1/Mcompliance[0][0]
                Ey = 1/Mcompliance[1][1]
                Ez = 1/Mcompliance[2][2]
                NUxy = -Mcompliance[0][1]*Ex
                NUxz = -Mcompliance[0][2]*Ex
                NUyz = -Mcompliance[1][2]*Ey
                Gxy = 1/Mcompliance[5][5]
                Gxz = 1/Mcompliance[4][4]
                Gyz = 1/Mcompliance[3][3]
                self.INMaterialName.insert(self.TMaterials.item(current_row,0).text())
                self.INEx.insert(str(Ex))
                self.INEy.insert(str(Ey))
                self.INEz.insert(str(Ez))
                self.INNUxy.insert(str(NUxy))
                self.INNUxz.insert(str(NUxz))
                self.INNUyz.insert(str(NUyz))
                self.INGxy.insert(str(Gxy))
                self.INGxz.insert(str(Gxz))
                self.INGyz.insert(str(Gyz))
            else:
                self.MaterialTypeTool.setCurrentIndex(2)
                self.INMaterialType.setCurrentIndex(3)
                self.INMaterialName.clear()
                self.INMaterialName.insert(self.TMaterials.item(current_row,0).text())
                k = 2
                for i in [0,1,2,3,4,5,6]:
                    for j in range(i,6):
                        self.TAnisotropic.item(i,j).setText('')
                        self.TAnisotropic.setItem(
                            i, j, QTableWidgetItem(self.TMaterials.item(current_row,k).text())
                        )
                        k += 1
                        
        self.BAddMaterial.clicked.connect(add_material)
        self.BDeleteMaterial.clicked.connect(delete_material)
        self.BRegenElasticConstants.clicked.connect(regenerate_elastic_constants)
        
        # Boundary Conditions
        def BC_direction():
            if self.INBoundaryCondition.currentText() == "Tension" or self.INBoundaryCondition.currentText() == "Compression":
                self.INBCDirection.clear()
                self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"x-direction", None))
                self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"y-direction", None))
                self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"z-direction", None))
            elif self.INBoundaryCondition.currentText() == "Shear":
                self.INBCDirection.clear()
                self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"xy-direction", None))
                self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"xz-direction", None))
                self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"yz-direction", None))
            elif self.INBoundaryCondition.currentText() == "Homogenization":
                self.INBCDirection.clear()
                self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"6 RVE BCs", None))
        self.INBoundaryCondition.currentTextChanged.connect(BC_direction)
        
        def add_bcs():
            row = self.TBCs.rowCount()
            if self.INBoundaryCondition.currentText() == 'Homogenization':
                # Tension x-direction
                self.TBCs.insertRow(row)
                self.TBCs.setItem(row, 0, QTableWidgetItem(str(self.INBoundaryCondition.currentText() + " Tension")))
                self.TBCs.setItem(row, 1, QTableWidgetItem(str("x-direction")))
                # Tension y-direction
                self.TBCs.insertRow(row+1)
                self.TBCs.setItem(row+1, 0, QTableWidgetItem(str(self.INBoundaryCondition.currentText() + " Tension")))
                self.TBCs.setItem(row+1, 1, QTableWidgetItem(str("y-direction")))
                # Tension z-direction
                self.TBCs.insertRow(row+2)
                self.TBCs.setItem(row+2, 0, QTableWidgetItem(str(self.INBoundaryCondition.currentText() + " Tension")))
                self.TBCs.setItem(row+2, 1, QTableWidgetItem(str("z-direction")))
                # Shear xy-direction
                self.TBCs.insertRow(row+3)
                self.TBCs.setItem(row+3, 0, QTableWidgetItem(str(self.INBoundaryCondition.currentText() + " Shear")))
                self.TBCs.setItem(row+3, 1, QTableWidgetItem(str("xy-direction")))
                # Shear xz-direction
                self.TBCs.insertRow(row+4)
                self.TBCs.setItem(row+4, 0, QTableWidgetItem(str(self.INBoundaryCondition.currentText() + " Shear")))
                self.TBCs.setItem(row+4, 1, QTableWidgetItem(str("xz-direction")))
                # Shear yz-direction
                self.TBCs.insertRow(row+5)
                self.TBCs.setItem(row+5, 0, QTableWidgetItem(str(self.INBoundaryCondition.currentText() + " Shear")))
                self.TBCs.setItem(row+5, 1, QTableWidgetItem(str("yz-direction")))
            else:
                self.TBCs.insertRow(row)
                self.TBCs.setItem(row, 0, QTableWidgetItem(str(
                    self.INBoundaryCondition.currentText()))
                )
                self.TBCs.setItem(
                    row, 1, QTableWidgetItem(str(self.INBCDirection.currentText()))
                )

        def delete_bcs():
            current_row = self.TBCs.currentRow()
            if current_row < 0:
                return QMessageBox.warning(self, 'Warning','Please select a record to delete')

            button = QMessageBox.question(
                QMessageBox(),
                'Confirmation',
                'Are you sure that you want to delete the selected row?',
                QMessageBox.StandardButton.Yes |
                QMessageBox.StandardButton.No
            )
            if button == QMessageBox.StandardButton.Yes:
                if 'Homogenization' in self.TBCs.item(current_row,0).text():
                    Hrmv = []
                    for i in range(self.TBCs.rowCount()):
                        if 'Homogenization' in self.TBCs.item(i,0).text():
                            Hrmv.append(i)
                    Hcnt = 0
                    for i in Hrmv:
                        self.TBCs.removeRow(i-Hcnt)
                        Hcnt += 1
                else:
                    self.TBCs.removeRow(current_row)
                
        self.BAddBC.clicked.connect(add_bcs)
        self.BDeleteBC.clicked.connect(delete_bcs)
        
        # Write input file
        def write_input_file(BC_index):
            # Elastic Input File
            for i in range(self.TMaterials.rowCount()):
                if self.TMaterials.item(i,2).text() == 'Isotropic' or 'Transversely Isotropic' in self.TMaterials.item(i,2).text() or self.TMaterials.item(i,2).text() == 'Orthotropic':
                    if i == 0:
                        print("generating EP0")
                        elastic_parameters = open(ELASTIC_PARAMETERS_0,"w")
                    else:
                        print("generating EP1")
                        elastic_parameters = open(ELASTIC_PARAMETERS_1,"w")
                    iso = '0\n'
                    elastic_parameters.write(iso)
                    stiffness = '  ' + self.TMaterials.item(i,3).text() + '  ' + self.TMaterials.item(i,4).text() + '  ' + self.TMaterials.item(i,5).text() + '  0  0  0     Cu (MPa)\n' + '  ' + self.TMaterials.item(i,4).text() + '  ' + self.TMaterials.item(i,9).text() + '  ' + self.TMaterials.item(i,10).text() + '  0  0  0\n' + '  ' + self.TMaterials.item(i,5).text() + '  ' + self.TMaterials.item(i,10).text() + '  ' + self.TMaterials.item(i,14).text() + '  0  0  0\n' + '  0  0  0   ' + self.TMaterials.item(i,18).text() + '  0  0\n' + '  0  0  0  0  ' + self.TMaterials.item(i,21).text() + '  0\n' + '  0  0  0  0  0  ' + self.TMaterials.item(i,23).text()
                    elastic_parameters.write(stiffness)
                    elastic_parameters.close()
                elif self.TMaterials.item(i,2).text() == 'Anisotropic':
                    if i == 0:
                        elastic_parameters = open(ELASTIC_PARAMETERS_0,"w")
                    else:
                        elastic_parameters = open(ELASTIC_PARAMETERS_1,"w")
                    iso = '0\n'
                    elastic_parameters.write(iso)
                    stiffness = '  ' + self.TMaterials.item(i,3).text() + '  ' + self.TMaterials.item(i,4).text() + '  ' + self.TMaterials.item(i,5).text() + '  ' + self.TMaterials.item(i,6).text() + '  ' + self.TMaterials.item(i,7).text() + '  ' + self.TMaterials.item(i,8).text() + '     Cu (MPa)\n' + '  ' + self.TMaterials.item(i,4).text() + '  ' + self.TMaterials.item(i,9).text() + '  ' + self.TMaterials.item(i,10).text() + '  ' + self.TMaterials.item(i,11).text() + '  ' + self.TMaterials.item(i,12).text() + '  ' + self.TMaterials.item(i,13).text() + '\n' + '  ' + self.TMaterials.item(i,5).text() + '  ' + self.TMaterials.item(i,10).text() + '  ' + self.TMaterials.item(i,14).text() + '  ' + self.TMaterials.item(i,15).text() + '  ' + self.TMaterials.item(i,16).text() + '  ' + self.TMaterials.item(i,17).text() + '\n' + '  ' + self.TMaterials.item(i,6).text() + '  ' + self.TMaterials.item(i,11).text() + '  ' + self.TMaterials.item(i,15).text() + '  ' + self.TMaterials.item(i,18).text() + '  ' + self.TMaterials.item(i,19).text() + '  ' + self.TMaterials.item(i,20).text() + '\n' + '  ' + self.TMaterials.item(i,7).text() + '  ' + self.TMaterials.item(i,12).text() + '  ' + self.TMaterials.item(i,16).text() + '  ' + self.TMaterials.item(i,19).text() + '  ' + self.TMaterials.item(i,21).text() + '  ' + self.TMaterials.item(i,22).text() + '\n' + '  ' + self.TMaterials.item(i,8).text() + '  ' + self.TMaterials.item(i,13).text() + '  ' + self.TMaterials.item(i,17).text() + '  ' + self.TMaterials.item(i,20).text() + '  ' + self.TMaterials.item(i,22).text() + '  ' +  self.TMaterials.item(i,23).text()
                    elastic_parameters.write(stiffness)
                    elastic_parameters.close()
            
            # EVPFFT input parameters file
            evpfft_lattice_input = open(EVPFFT_INPUT,"w")
            modes = '2 0 0 0               NPHMX, NMODMX, NTWMMX, NSYSMX\n'
            evpfft_lattice_input.write(modes)
            dimensions = str(int(self.INNumberOfVoxelsX.text())) + ' ' + str(int(self.INNumberOfVoxelsY.text())) + ' ' + str(int(self.INNumberOfVoxelsZ.text())) + '               x-dim, y-dim, z-dim\n'
            evpfft_lattice_input.write(dimensions)
            dx, dy, dz = self.VoxelResolution
            nph_delt = '2                      number of phases (nph)\n' + f'{dx:.4f} {dy:.4f} {dz:.4f}             RVE dimensions (delt)\n' + '* name and path of microstructure file (filetext)\n'
            evpfft_lattice_input.write(nph_delt)
            vtkfile = f'{VTK_OUTPUT}\n'
            evpfft_lattice_input.write(vtkfile)
            for i in range(2):
                if not self.TMaterials.item(i,2) or self.TMaterials.item(i,2).text() == 'Ideal Gas':
                    if not self.TMaterials.item(i,2) and i == 1 or self.TMaterials.item(i,1).text() == 'Void':
                        phase1 = '*INFORMATION ABOUT PHASE #1\n' + '1                          igas(iph)\n' + '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + 'dummy\n' + 'dummy\n'
                    else:
                        phase2 = '*INFORMATION ABOUT PHASE #2\n' + '1                          igas(iph)\n' + '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + 'dummy\n' + 'dummy\n'
                else:
                    if i == 0:
                        efile = f'{ELASTIC_PARAMETERS_0}'
                    else:
                        efile = f'{ELASTIC_PARAMETERS_1}'
                        
                    if self.TMaterials.item(i,1).text() == 'Void':
                        phase1 = '*INFORMATION ABOUT PHASE #1\n' + '0                          igas(iph)\n' + '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' +  f'{PLASTIC_PARAMETERS}\n' + efile + '\n'
                    else:
                        phase2 = '*INFORMATION ABOUT PHASE #1\n' + '0                          igas(iph)\n' + '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' +  f'{PLASTIC_PARAMETERS}\n' + efile + '\n'
            evpfft_lattice_input.write(phase1)
            evpfft_lattice_input.write(phase2)
            if self.TBCs.item(BC_index,1).text() == "x-direction":
                if "Tension" in self.TBCs.item(BC_index,0).text():
                    test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    1       1       1           iudot     |    flag for vel.grad.\n' + '    1       0       1                     |    (0:unknown-1:known)\n' + '    1       1       0                     |\n' + '                                          |\n' + '   1.0     0.        0.          udot     |    vel.grad\n' + '    0.      0.      0.                  |\n' + '    0.       0.         0.                |\n' + '                                          |\n' + '    0       0        0           iscau    |    flag for Cauchy\n' + '            1        0                    |\n' + '                     1                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
                elif "Compression" in self.TBCs.item(BC_index,0).text():
                    test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    1       1       1           iudot     |    flag for vel.grad.\n' + '    1       0       1                     |    (0:unknown-1:known)\n' + '    1       1       0                     |\n' + '                                          |\n' + '   -1.0     0.        0.          udot    |    vel.grad\n' + '    0.      0.      0.                  |\n' + '    0.       0.         0.                |\n' + '                                          |\n' + '    0       0        0           iscau    |    flag for Cauchy\n' + '            1        0                    |\n' + '                     1                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
                else:
                    print("INVALID BOUNDARY CONDITION")
            elif self.TBCs.item(BC_index,1).text() == "y-direction":
                if "Tension" in self.TBCs.item(BC_index,0).text():
                    test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    0       1       1           iudot     |    flag for vel.grad.\n' + '    1       1       1                     |    (0:unknown-1:known)\n' + '    1       1       0                     |\n' + '                                          |\n' + '   0.     0.        0.          udot    |    vel.grad\n' + '    0.      1.0      0.                  |\n' + '    0.       0.         0.                |\n' + '                                          |\n' + '    1       0        0           iscau    |    flag for Cauchy\n' + '            0        0                    |\n' + '                     1                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
                elif "Compression" in self.TBCs.item(BC_index,0).text():
                    test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    0       1       1           iudot     |    flag for vel.grad.\n' + '    1       1       1                     |    (0:unknown-1:known)\n' + '    1       1       0                     |\n' + '                                          |\n' + '   0.     0.        0.          udot    |    vel.grad\n' + '    0.      -1.0      0.                  |\n' + '    0.       0.         0.                |\n' + '                                          |\n' + '    1       0        0           iscau    |    flag for Cauchy\n' + '            0        0                    |\n' + '                     1                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
                else:
                    print("INVALID BOUNDARY CONDITION")
            elif self.TBCs.item(BC_index,1).text() == "z-direction":
                if "Tension" in self.TBCs.item(BC_index,0).text():
                    test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    0       1       1           iudot     |    flag for vel.grad.\n' + '    1       0       1                     |    (0:unknown-1:known)\n' + '    1       1       1                     |\n' + '                                          |\n' + '   0.     0.        0.          udot    |    vel.grad\n' + '    0.      0.      0.                  |\n' + '    0.       0.         1.0                |\n' + '                                          |\n' + '    1       0        0           iscau    |    flag for Cauchy\n' + '            1        0                    |\n' + '                     0                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
                elif "Compression" in self.TBCs.item(BC_index,0).text():
                    test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    0       1       1           iudot     |    flag for vel.grad.\n' + '    1       0       1                     |    (0:unknown-1:known)\n' + '    1       1       1                     |\n' + '                                          |\n' + '   0.0     0.        0.          udot    |    vel.grad\n' + '    0.      0.0      0.                  |\n' + '    0.       0.         -1.0                |\n' + '                                          |\n' + '    1       0        0           iscau    |    flag for Cauchy\n' + '            1        0                    |\n' + '                     0                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
                else:
                    print("INVALID BOUNDARY CONDITION")
            elif "Shear" in self.TBCs.item(BC_index,0).text():
                if "xy-direction" in self.TBCs.item(BC_index,1).text():
                    test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    1       1       1           iudot     |    flag for vel.grad.\n' + '    1       1       1                     |    (0:unknown-1:known)\n' + '    1       1       1                     |\n' + '                                          |\n' + '   0.     1.0        0.          udot    |    vel.grad\n' + '    1.0      0.      0.                  |\n' + '    0.       0.         0.                |\n' + '                                          |\n' + '    0       0        0           iscau    |    flag for Cauchy\n' + '            0        0                    |\n' + '                     0                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
                elif "xz-direction" in self.TBCs.item(BC_index,1).text():
                    test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    1       1       1           iudot     |    flag for vel.grad.\n' + '    1       1       1                     |    (0:unknown-1:known)\n' + '    1       1       1                     |\n' + '                                          |\n' + '   0.     0.        1.0          udot    |    vel.grad\n' + '    0.      0.      0.                  |\n' + '    1.0       0.         0.                |\n' + '                                          |\n' + '    0       0        0           iscau    |    flag for Cauchy\n' + '            0        0                    |\n' + '                     0                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
                elif "yz-direction" in self.TBCs.item(BC_index,1).text():
                    test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    1       1       1           iudot     |    flag for vel.grad.\n' + '    1       1       1                     |    (0:unknown-1:known)\n' + '    1       1       1                     |\n' + '                                          |\n' + '   0.     0.        0.          udot    |    vel.grad\n' + '    0.      0.      1.0                  |\n' + '    0.       1.0         0.                |\n' + '                                          |\n' + '    0       0        0           iscau    |    flag for Cauchy\n' + '            0        0                    |\n' + '                     0                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
            else:
                print("INVALID BOUNDARY CONDITION")
            evpfft_lattice_input.write(test_conditions)
            other = '* other\n' + '0.0001         eqincr (if ictrl>=0) or tdot (if ictrl=-1)\n' + '-1              ictrl (1-6: strain comp, 0: VM eq, -1: tdot)\n'
            evpfft_lattice_input.write(other)
            run_conditions = '*INFORMATION ABOUT RUN CONDITIONS\n' + self.INNumberOfSteps.text() + '              nsteps\n' + '0.00001         err\n' + '50              itmax\n' + '0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?\n' + '0               ISAVE write grain states in STRESS.OUT (1) or not (0)?\n' + '1               IUPDATE update tex & RVE dim (1) or not (0)?\n' + '0               IUPHARD\n' + '1               IWTEX\n' + '1 10            IWFIELDS,IWSTEP\n' + '0               ITHERMO (if ithermo=1, next line is filethermo)\n' + 'dummy\n'
            evpfft_lattice_input.write(run_conditions)
            evpfft_lattice_input.close()
            
            
        # Single Run of EVPFFT
        self.run_cnt = 0
        def single_EVPFFT(BC_index):
            if self.p == None:
                write_input_file(BC_index)
                self.p = QProcess()
                self.p.readyReadStandardOutput.connect(handle_stdout)
                self.p.readyReadStandardError.connect(handle_stderr)
                self.p.stateChanged.connect(handle_state)
                self.p.finished.connect(process_finished)
                self.p.start("evpfft",["-f", EVPFFT_INPUT, "-m", "2"])
                self.progress_re = re.compile("       Current  Time  STEP = (\d+)")
                self.run_cnt += 1
        def simple_percent_parser(output):
            m = self.progress_re.search(output)
            if m:
                pc_complete = m.group(1)
                return int(pc_complete)
        def process_finished():
            self.RunOutputProgress.setValue(100)
            self.p.close()
            self.p = None
        def handle_stdout():
            data = self.p.readAllStandardOutput()
            stdout = bytes(data).decode("utf8")
            progress = simple_percent_parser(stdout)
            if progress:
                self.RunOutputProgress.setValue((progress/int(self.INNumberOfSteps.text()))*100)
            self.RunOutputWindow.appendPlainText(stdout)
        def handle_stderr():
            data = self.p.readAllStandardError()
            stderr = bytes(data).decode("utf8")
            self.RunOutputWindow.appendPlainText(stderr)
        def handle_state(state):
            states = {
                QProcess.NotRunning: 'Finished',
                QProcess.Starting: 'Starting EVPFFT',
                QProcess.Running: 'Running EVPFFT',
            }
            self.state_name = states[state]
            self.RunOutputWindow.appendPlainText(f"{self.state_name}")
        
        # Batch Run of EVPFFT
        def batch_EVPFFT():
            for BC_index in range(self.TBCs.rowCount()):
                self.BRunEVPFFT.clicked.connect(single_EVPFFT(BC_index))
                self.p.waitForStarted()
                while self.p != None:
                    QApplication.processEvents()
                    
                # Save Output Files
                
                    
                # Generate Homogenized Elastic Constants
                if "Homogenization" in self.TBCs.item(BC_index,0).text():
                    self.BHomogenization.setEnabled(True)
                    self.THomogenization.setEnabled(True)
                    with open("str_str.out", newline='') as f:
                        reader = csv.reader(f)
                        self.ss_data = list(reader)
                    s11 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    s22 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    s33 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    s12 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    s13 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    s23 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    e11 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    e22 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    e33 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    e12 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    e13 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    e23 = [0 for i in range(int(self.INNumberOfSteps.text()))]
                    for i in range(int(self.INNumberOfSteps.text())):
                        s11[i] = float(self.ss_data[i+1][6])
                        s22[i] = float(self.ss_data[i+1][7])
                        s33[i] = float(self.ss_data[i+1][8])
                        s12[i] = float(self.ss_data[i+1][11])
                        s13[i] = float(self.ss_data[i+1][10])
                        s23[i] = float(self.ss_data[i+1][9])
                        e11[i] = float(self.ss_data[i+1][0])
                        e22[i] = float(self.ss_data[i+1][1])
                        e33[i] = float(self.ss_data[i+1][2])
                        e12[i] = float(self.ss_data[i+1][5])
                        e13[i] = float(self.ss_data[i+1][4])
                        e23[i] = float(self.ss_data[i+1][3])
                    if self.TBCs.item(BC_index,1).text() == "x-direction":
                        self.HE11 = np.polyfit(e11,s11,1)
                        self.HNU12 = np.polyfit(e11,e22,1)
                        self.HNU13 = np.polyfit(e11,e33,1)
                    if self.TBCs.item(BC_index,1).text() == "y-direction":
                        self.HE22 = np.polyfit(e22,s22,1)
                        self.HNU21 = np.polyfit(e22,e11,1)
                        self.HNU23 = np.polyfit(e22,e33,1)
                    if self.TBCs.item(BC_index,1).text() == "z-direction":
                        self.HE33 = np.polyfit(e33,s33,1)
                        self.HNU31 = np.polyfit(e33,e11,1)
                        self.HNU32 = np.polyfit(e33,e22,1)
                    if self.TBCs.item(BC_index,1).text() == "xy-direction":
                        self.HG12 = np.polyfit(np.multiply(e12,2.),s12,1)
                    if self.TBCs.item(BC_index,1).text() == "xz-direction":
                        self.HG13 = np.polyfit(np.multiply(e13,2.),s13,1)
                    if self.TBCs.item(BC_index,1).text() == "yz-direction":
                        self.HG23 = np.polyfit(np.multiply(e23,2.),s23,1)
        
        # Connect run button to indiviual or batch run
        self.p = None
        def run_click():
            batch_EVPFFT()
        self.BRunEVPFFT.clicked.connect(run_click)
        
        # Preview Results
        def preview_results_click():
            # Delete previous views
            try:
                self.stl
            except:
                print('')
            else:
                pvsimple.Delete(self.stl)
            try:
                self.threshold
            except:
                print('')
            else:
                pvsimple.Delete(self.threshold)
            try:
                self.threshold2
            except:
                print('')
            else:
                pvsimple.Delete(self.threshold2)
            # Render new view
            self.results_reader = pvsimple.XDMFReader(FileNames = "micro_state_timestep_10.xdmf")
            pvsimple.SetDisplayProperties(Representation = "Surface")
            self.threshold2 = pvsimple.Threshold(registrationName='results_threshold', Input = self.results_reader, Scalars = "phase_id", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 2, LowerThreshold = 1, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
            display = pvsimple.Show(self.threshold2, self.render_view)
            # Select what variable you want to display
            pvsimple.GetAnimationScene().GoToLast()
            pvsimple.ColorBy(display,('CELLS',str(self.INPreviewResults.currentText())))
            vmstressLUT = pvsimple.GetColorTransferFunction(str(self.INPreviewResults.currentText()))
            r = self.results_reader.CellData.GetArray(str(self.INPreviewResults.currentText())).GetRange()
            vmstressLUT.RescaleTransferFunction(r[0], r[1]/2)
            display.SetScalarBarVisibility(self.render_view, True)
            pvsimple.HideUnusedScalarBars(self.render_view)
            # Add time filter
            threshold1 = pvsimple.FindSource('results_threshold')
            annotateTimeFilter1 = pvsimple.AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=threshold1)
            annotateTimeFilter1Display = pvsimple.Show(annotateTimeFilter1, self.render_view, 'TextSourceRepresentation')
            # Remove old view / reset cameras
            pvsimple.Hide(self.results_reader)
            self.render_view.ResetCamera()
            self.render_view.StillRender()
        self.BPreviewResults.clicked.connect(preview_results_click)
        self.BPreviewResults.clicked.connect(lambda: self.OutputWindows.setCurrentIndex(0))
        
        # Stress vs Strain Plot
        def plot_ss_click():
            # Get the stress-strain data
            try:
                self.Plot.figure
            except:
                with open("str_str.out", newline='') as f:
                    reader = csv.reader(f)
                    self.ss_data = list(reader)
                x = [0 for i in range(int(self.INNumberOfSteps.text()))]
                y = [0 for i in range(int(self.INNumberOfSteps.text()))]
                for i in range(int(self.INNumberOfSteps.text())):
                    if str(self.INPlotSS.currentText()) == 'S11 vs E11':
                        xcol = 0
                        ycol = 6
                    elif str(self.INPlotSS.currentText()) == 'S22 vs E22':
                        xcol = 1
                        ycol = 7
                    elif str(self.INPlotSS.currentText()) == 'S33 vs E33':
                        xcol = 2
                        ycol = 8
                    x[i] = float(self.ss_data[i+1][xcol])
                    y[i] = float(self.ss_data[i+1][ycol])
                # Plot data
                self.Plot.figure = Figure()
                self.Plot.ax = self.Plot.figure.add_subplot()
                self.Plot.ax.plot(x,y)
                self.Plot.ax.set_xlabel('STRAIN')
                self.Plot.ax.set_ylabel('STRESS')
                self.Plot.figure.tight_layout()
                # Display plot and toolbar
                layout = QVBoxLayout()
                self.Plot.setLayout(layout)
                self.Plot.canvas = FigureCanvasQTAgg(self.Plot.figure)
                layout.addWidget(self.Plot.canvas)
                self.toolbar = NavigationToolbar2QT(self.Plot.canvas,self.Plot)
                layout.addWidget(self.toolbar)
            else:
                self.timer = QTimer()
                self.timer.setInterval(100)
                self.timer.timeout.connect(update_plot)
                self.timer.start()
        def update_plot():
            if self.run_cnt > 1:
                with open("str_str.out", newline='') as f:
                    reader = csv.reader(f)
                    self.ss_data = list(reader)
            x = [0 for i in range(int(self.INNumberOfSteps.text()))]
            y = [0 for i in range(int(self.INNumberOfSteps.text()))]
            for i in range(int(self.INNumberOfSteps.text())):
                if str(self.INPlotSS.currentText()) == 'S11 vs E11':
                    xcol = 0
                    ycol = 6
                elif str(self.INPlotSS.currentText()) == 'S22 vs E22':
                    xcol = 1
                    ycol = 7
                elif str(self.INPlotSS.currentText()) == 'S33 vs E33':
                    xcol = 2
                    ycol = 8
                x[i] = float(self.ss_data[i+1][xcol])
                y[i] = float(self.ss_data[i+1][ycol])
            self.Plot.ax.cla()
            self.Plot.ax.plot(x,y)
            self.Plot.ax.set_xlabel('STRAIN')
            self.Plot.ax.set_ylabel('STRESS')
            self.Plot.figure.tight_layout()
            self.Plot.canvas.draw()
            self.timer.stop()
        self.BPlotSS.clicked.connect(plot_ss_click)
        self.BPlotSS.clicked.connect(lambda: self.OutputWindows.setCurrentIndex(1))
        
        # Generate Homogenized Elastic Constants
        def homogenization_click():
            self.THomogenization.setItem(0,0,QTableWidgetItem(str(self.HE11[0])))
            self.THomogenization.setItem(1,0,QTableWidgetItem(str(self.HE22[0])))
            self.THomogenization.setItem(2,0,QTableWidgetItem(str(self.HE33[0])))
            self.THomogenization.setItem(3,0,QTableWidgetItem(str(-self.HNU12[0])))
            self.THomogenization.setItem(4,0,QTableWidgetItem(str(-self.HNU21[0])))
            self.THomogenization.setItem(5,0,QTableWidgetItem(str(-self.HNU13[0])))
            self.THomogenization.setItem(6,0,QTableWidgetItem(str(-self.HNU31[0])))
            self.THomogenization.setItem(7,0,QTableWidgetItem(str(-self.HNU23[0])))
            self.THomogenization.setItem(8,0,QTableWidgetItem(str(-self.HNU32[0])))
            self.THomogenization.setItem(9,0,QTableWidgetItem(str(self.HG12[0])))
            self.THomogenization.setItem(10,0,QTableWidgetItem(str(self.HG13[0])))
            self.THomogenization.setItem(11,0,QTableWidgetItem(str(self.HG23[0])))
        self.BHomogenization.clicked.connect(homogenization_click)
        
        # Open Paraview
        def open_paraview_click():
            os.system("open " + "micro_state_timestep_10.xdmf")
        self.BOpenParaview.clicked.connect(open_paraview_click)
        
        # Warning Message Popup
        def warning_message(msg):
            message = QMessageBox()
            message.setText(msg)
            message.exec()
