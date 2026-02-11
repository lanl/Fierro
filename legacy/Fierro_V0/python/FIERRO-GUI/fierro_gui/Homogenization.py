from PySide6.QtWidgets import (QTableWidgetItem, QMessageBox, QApplication, QFileDialog)
from PySide6.QtCore import (QCoreApplication, QProcess, Qt, QProcessEnvironment)
import re
import csv
import numpy as np
import subprocess
import os
import vtk
#import paraview.simple as paraview.simple
from paraview.simple import *
from Homogenization_WInput import *
import DeveloperInputs
from importlib import reload
from ImageToVTK import *
from Reload_Geometry import *
#from ReadHDF5 import *

# ==============================================
# ======= EVPFFT SOLVER LATTICE PIPELINE =======
# ==============================================

# Warning Message Popup
def warning_message(msg):
    message = QMessageBox()
    message.setText(msg)
    message.exec()

def Homogenization(self):
    # Connect tab buttons to settings windows
    #self.BImportPart.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(2))
    #self.BImportHDF5Part.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(5))
    #self.BDefineMaterial.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(4))
    #self.BViewResults.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(7))
    #self.BGlobalMesh.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(1))
    
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
        try:
            self.vtk_reader
            if str(self.INRegion_2.currentText()) == 'global':
                # Remove all objects from window view
                SetActiveView(self.render_view)
                renderer = self.render_view.GetRenderer()
                renderer.RemoveAllViewProps()
                self.render_view.Update()
                self.render_view.StillRender()
                
                # Show void region only
                setattr(self, self.variable_name, paraview.simple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Below Lower Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
                self.display = paraview.simple.Show(getattr(self, self.variable_name), self.render_view)
                self.render_view.ResetCamera()
                self.render_view.StillRender()
            else:
                # Remove all objects from window view
                SetActiveView(self.render_view)
                renderer = self.render_view.GetRenderer()
                renderer.RemoveAllViewProps()
                self.render_view.Update()
                self.render_view.StillRender()
                
                # Show material region only
                setattr(self, self.variable_name, paraview.simple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0))
                self.display = paraview.simple.Show(getattr(self, self.variable_name), self.render_view)
                self.render_view.ResetCamera()
                self.render_view.StillRender()
        except:
            print("Polycrystalline model was imported")
    self.INRegion_2.currentIndexChanged.connect(material_region)
    
    def material_class():
        if str(self.INMaterialType.currentText()) == 'Isotropic':
            self.MaterialTypeTool.setCurrentIndex(0)
        if str(self.INMaterialType.currentText()) == 'Transversely Isotropic':
            self.MaterialTypeTool.setCurrentIndex(1)
        if str(self.INMaterialType.currentText()) == 'Orthotropic':
            self.MaterialTypeTool.setCurrentIndex(3)
        if str(self.INMaterialType.currentText()) == 'Anisotropic':
            self.MaterialTypeTool.setCurrentIndex(2)
        if str(self.INMaterialType.currentText()) == 'Ideal Gas':
            self.MaterialTypeTool.setCurrentIndex(4)
    self.INMaterialType.currentIndexChanged.connect(material_class)
    
    def add_material():
        warning_flag = 0
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
                    str(self.INMaterialType.currentText()))
                )
                
                # Add material as an option for material assignment
                self.INMaterial_2.clear()
                for i in range(self.TMaterials.rowCount()):
                    self.INMaterial_2.addItem(self.TMaterials.item(i,0).text())
                
                # Clear fields
                self.INMaterialName.clear()
        else:
            if self.INMaterialType.currentText() == 'Isotropic':
                row = self.TParts.rowCount()
                if row > 0:
                    in_file_path = self.TParts.item(0,10).text()
                    file_type = os.path.splitext(in_file_path)[1].lower()
                    if ".txt" in file_type:
                        warning_message('ERROR: Crystal properties cannot be isotropic')
                        warning_flag = 1
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
                    # Add material as an option for material assignment
                    self.INMaterial_2.clear()
                    for i in range(self.TMaterials.rowCount()):
                        self.INMaterial_2.addItem(self.TMaterials.item(i,0).text())
            if warning_flag == 0:
                # Fill out material definition table
                row = self.TMaterials.rowCount()
                self.TMaterials.insertRow(row)
                self.TMaterials.setItem(row, 0, QTableWidgetItem(
                    self.INMaterialName.text().strip())
                )
                if str(self.INMaterialType.currentText()) == 'Transversely Isotropic':
                    self.TMaterials.setItem(row, 1, QTableWidgetItem(
                        str(self.INMaterialType.currentText() + ' ' + self.INIsotropicPlane.currentText()))
                    )
                else:
                    self.TMaterials.setItem(row, 1, QTableWidgetItem(
                        str(self.INMaterialType.currentText()))
                    )
                self.TMaterials.setItem(
                    row, 2, QTableWidgetItem(str(C11))
                )
                self.TMaterials.setItem(
                    row, 3, QTableWidgetItem(str(C12))
                )
                self.TMaterials.setItem(
                    row, 4, QTableWidgetItem(str(C13))
                )
                self.TMaterials.setItem(
                    row, 8, QTableWidgetItem(str(C22))
                )
                self.TMaterials.setItem(
                    row, 9, QTableWidgetItem(str(C23))
                )
                self.TMaterials.setItem(
                    row, 13, QTableWidgetItem(str(C33))
                )
                self.TMaterials.setItem(
                    row, 17, QTableWidgetItem(str(C44))
                )
                self.TMaterials.setItem(
                    row, 20, QTableWidgetItem(str(C55))
                )
                self.TMaterials.setItem(
                    row, 22, QTableWidgetItem(str(C66))
                )
                for i in [5,6,7,10,11,12,14,15,16,18,19,21]:
                   self.TMaterials.setItem(row, i, QTableWidgetItem('0'))
                
                # Add material as an option for material assignment
                self.INMaterial_2.clear()
                for i in range(self.TMaterials.rowCount()):
                    self.INMaterial_2.addItem(self.TMaterials.item(i,0).text())
                
                # Clear fields
                self.INMaterialName.clear()
            else:
                warning_flag = 0
    self.BAddMaterial.clicked.connect(add_material)
            
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
            # delete from material assignment options
            self.INMaterial_2.clear()
            for i in range(self.TMaterials.rowCount()):
                self.INMaterial_2.addItem(self.TMaterials.item(i,0).text())
    self.BDeleteMaterial.clicked.connect(delete_material)
            
    def regenerate_elastic_constants():
        current_row = self.TMaterials.currentRow()
        if current_row < 0:
            return QMessageBox.warning(QMessageBox(),"Warning","Please select a material from the table")
            
        # Define Stiffness Matrix
        Mstiffness = [[float(self.TMaterials.item(current_row,2).text()), float(self.TMaterials.item(current_row,3).text()), float(self.TMaterials.item(current_row,4).text()),  float(self.TMaterials.item(current_row,5).text()), float(self.TMaterials.item(current_row,6).text()), float(self.TMaterials.item(current_row,7).text())], [float(self.TMaterials.item(current_row,3).text()), float(self.TMaterials.item(current_row,8).text()), float(self.TMaterials.item(current_row,9).text()),  float(self.TMaterials.item(current_row,10).text()), float(self.TMaterials.item(current_row,11).text()), float(self.TMaterials.item(current_row,12).text())], [float(self.TMaterials.item(current_row,4).text()), float(self.TMaterials.item(current_row,9).text()), float(self.TMaterials.item(current_row,13).text()), float(self.TMaterials.item(current_row,14).text()), float(self.TMaterials.item(current_row,15).text()), float(self.TMaterials.item(current_row,16).text())], [float(self.TMaterials.item(current_row,5).text()), float(self.TMaterials.item(current_row,10).text()), float(self.TMaterials.item(current_row,14).text()), float(self.TMaterials.item(current_row,17).text()), float(self.TMaterials.item(current_row,18).text()), float(self.TMaterials.item(current_row,19).text())], [float(self.TMaterials.item(current_row,6).text()), float(self.TMaterials.item(current_row,11).text()), float(self.TMaterials.item(current_row,15).text()), float(self.TMaterials.item(current_row,18).text()), float(self.TMaterials.item(current_row,20).text()), float(self.TMaterials.item(current_row,21).text())], [float(self.TMaterials.item(current_row,7).text()), float(self.TMaterials.item(current_row,12).text()), float(self.TMaterials.item(current_row,16).text()), float(self.TMaterials.item(current_row,19).text()), float(self.TMaterials.item(current_row,21).text()), float(self.TMaterials.item(current_row,22).text())]]
        if self.TMaterials.item(current_row,1).text() == 'Isotropic':
            Mcompliance = np.linalg.inv(Mstiffness)
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
        elif 'Transversely Isotropic' in self.TMaterials.item(current_row,1).text():
            Mcompliance = np.linalg.inv(Mstiffness)
            self.MaterialTypeTool.setCurrentIndex(1)
            self.INMaterialType.setCurrentIndex(1)
            self.INMaterialName.clear()
            self.INEip.clear()
            self.INNUip.clear()
            self.INEop.clear()
            self.INNUop.clear()
            self.INGop.clear()
            if 'x-y plane' in self.TMaterials.item(current_row,1).text():
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
            elif 'x-z plane' in self.TMaterials.item(current_row,1).text():
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
            elif 'y-z plane' in self.TMaterials.item(current_row,1).text():
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
        elif self.TMaterials.item(current_row,1).text() == 'Orthotropic':
            Mcompliance = np.linalg.inv(Mstiffness)
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
    self.BRegenElasticConstants.clicked.connect(regenerate_elastic_constants)
    
    # Add material assignment to region
    def add_material_assignment():
        warning_flag = 0
        for i in range(self.TMaterialAssignment.rowCount()):
            if str(self.INRegion_2.currentText()) == self.TMaterialAssignment.item(i,0).text():
                warning_message('ERROR: There is already a material assigned to this region')
                warning_flag = 1
        if warning_flag == 0:
            row = self.TMaterialAssignment.rowCount()
            self.TMaterialAssignment.insertRow(row)
            self.TMaterialAssignment.setItem(row, 0, QTableWidgetItem(self.INRegion_2.currentText()))
            self.TMaterialAssignment.setItem(row, 1, QTableWidgetItem(self.INMaterial_2.currentText()))
        else:
            warning_flag = 0
    self.BAddMaterialAssignment.clicked.connect(add_material_assignment)
    
    # Delete regional material assignment
    def delete_material_assignemnt():
        current_row = self.TMaterialAssignment.currentRow()
        if current_row < 0:
            return QMessageBox.warning(QMessageBox(),"Warning","Please select a record to delete")

        button = QMessageBox.question(
            QMessageBox(),
            'Confirmation',
            'Are you sure that you want to delete the materail assignemnt?',
            QMessageBox.Yes |
            QMessageBox.No
        )
        if button == QMessageBox.StandardButton.Yes:
            self.TMaterialAssignment.removeRow(current_row)
    self.BDeleteMaterialAssignment.clicked.connect(delete_material_assignemnt)
    
    # Warn user if no material assignment was made
    def warning_no_material():
        index = self.NavigationMenu.currentIndex()
        if index > 4 and self.TMaterialAssignment.rowCount() == 0 and self.TMaterials.rowCount() > 0:
            warning_message('WARNING: No materials were assigned')
    self.NavigationMenu.currentChanged.connect(warning_no_material)
    
#    # Show results immediately when postprocessing tab is pressed
#    def show_results():
#        index = self.NavigationMenu.currentIndex()
#        if index == 8 and self.run == 1:
#            preview_results_click()
#    self.NavigationMenu.currentChanged.connect(show_results)
    
#    # Boundary Conditions
#    def BC_direction():
#        if self.INBoundaryCondition.currentText() == "Tension" or self.INBoundaryCondition.currentText() == "Compression":
#            self.INBCDirection.clear()
#            self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"x", None))
#            self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"y", None))
#            self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"z", None))
#        elif self.INBoundaryCondition.currentText() == "Shear":
#            self.INBCDirection.clear()
#            self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"xy", None))
#            self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"xz", None))
#            self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"yz", None))
#        elif self.INBoundaryCondition.currentText() == "Homogenization":
#            self.INBCDirection.clear()
#            self.INBCDirection.addItem(QCoreApplication.translate("MainWindow", u"6 RVE BCs", None))
#    self.INBoundaryCondition.currentTextChanged.connect(BC_direction)
    
    # Solver Settings
    def settings():
        if self.INHomogenizationSettingC.isChecked():
            self.INNumberOfSteps.setEnabled(True)
            self.INErrorTolerance.setEnabled(True)
            self.INMaxIterations.setEnabled(True)
        elif self.INHomogenizationSettingD.isChecked():
            self.INNumberOfSteps.setText('1')
            self.INErrorTolerance.setText('0.00001')
            self.INMaxIterations.setText('50')
            self.INNumberOfSteps.setEnabled(False)
            self.INErrorTolerance.setEnabled(False)
            self.INMaxIterations.setEnabled(False)
    self.INHomogenizationSettingC.toggled.connect(settings)
    self.INHomogenizationSettingD.toggled.connect(settings)
    
    # Single Run of EVPFFT
    self.run_cnt = 0
    self.IterationError = False
    self.ConvergenceError = False
    self.min_iterations_old = 0
    def single_EVPFFT(BC_index):
        # Restart progress bar
        if BC_index == 0:
            self.RunOutputProgress.setValue(0)
        
        # Create location to save files
        self.working_directory = os.path.join(self.evpfft_dir, f'{self.job_name}')
        if not os.path.exists(self.working_directory):
            os.makedirs(self.working_directory)
            
        # Create input files
        self.ELASTIC_PARAMETERS_0 = os.path.join(self.working_directory, 'elastic_parameters_0.txt')
        self.ELASTIC_PARAMETERS_1 = os.path.join(self.working_directory, 'elastic_parameters_1.txt')
        self.PLASTIC_PARAMETERS = os.path.join(self.working_directory, 'plastic_parameters.txt')
        self.EVPFFT_INPUT = os.path.join(self.working_directory, 'evpfft_lattice_input.txt')
        Homogenization_WInput(self,BC_index)
        
        # Define executable path
        in_file_path = self.TParts.item(0,10).text()
        file_type = os.path.splitext(in_file_path)[1].lower()
        reload(DeveloperInputs)
        if self.UserConfig == "Developer":
            if self.INHParallel.isChecked():
                executable_path = "mpirun"
            else:
                executable_path = DeveloperInputs.fierro_evpfft_exe
        elif self.UserConfig == "User":
            if self.INHParallel.isChecked():
                executable_path = "mpirun"
            else:
                executable_path = "evpfft"
        if ".vtk" in file_type:
            if self.INHParallel.isChecked():
                arguments = ["-np", self.INmpiRanks.value(), "evpfft", "-f", self.EVPFFT_INPUT, "-m", "2"]
            else:
                arguments = ["-f", self.EVPFFT_INPUT, "-m", "2"]
        else:
            if self.INHParallel.isChecked():
                arguments = ["-np", self.INmpiRanks.value(), "evpfft", "-f", self.EVPFFT_INPUT]
            else:
                arguments = ["-f", self.EVPFFT_INPUT]
        
        # Make a new directory to store all of the outputs
        if BC_index == 0:
            folder1 = "Tension"
            folder2 = "x"
        if BC_index == 1:
            folder1 = "Tension"
            folder2 = "y"
        if BC_index == 2:
            folder1 = "Tension"
            folder2 = "z"
        if BC_index == 3:
            folder1 = "Shear"
            folder2 = "xy"
        if BC_index == 4:
            folder1 = "Shear"
            folder2 = "xz"
        if BC_index == 5:
            folder1 = "Shear"
            folder2 = "yz"
        self.simulation_directory = os.path.join(self.working_directory, folder1, folder2)
        if not os.path.exists(self.simulation_directory):
            os.makedirs(self.simulation_directory)
        
        # Set up the environment
        self.p = QProcess()
        env = QProcessEnvironment.systemEnvironment()
        if not env.contains("OMP_PROC_BIND"):
            env.insert("OMP_PROC_BIND", "spread")
        if not env.contains("OMP_NUM_THREADS"):
            env.insert("OMP_NUM_THREADS", "1")
        self.p.setProcessEnvironment(env)
        # Set up the states
        self.p.setWorkingDirectory(self.simulation_directory)
#        self.p.readyReadStandardOutput.connect(handle_stdout)
#        self.p.readyReadStandardError.connect(handle_stderr)
        self.p.stateChanged.connect(handle_state)
        self.p.finished.connect(lambda: process_finished(BC_index))
        try:
            self.p.start(executable_path, arguments)
        except Exception as e:
            self.warning_message("ERROR: evpfft executable")
            return
        self.progress_re = re.compile("       Current  Time  STEP = (\\d+)")
        self.run_cnt += 1
        # Count how many iterations were taken towards the solution
        self.iterations_re = re.compile(r" ITER = (\d+)")
            
    def convergence_check(output):
        iterations = self.iterations_re.findall(output)
        if iterations:
            return int(iterations[-1])
    def process_finished(num):
        handle_stdout()
        handle_stderr()
        self.RunOutputProgress.setValue(((num+1)/6)*100)
        # If the number of iterations didn't exceed 10, write out an error
        if self.min_iterations_old < 10:
            self.IterationError = True
        # If the iterations reached the maximum number, write out an error
        if self.min_iterations_old >= int(self.INMaxIterations.text()):
            self.ConvergenceError = True
        self.min_iterations_old = 0
        self.p.close()
        self.p = None
        self.RunOutputWindow.appendPlainText("Finished")
    def handle_stdout():
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        # check for convergence
        self.min_iterations = convergence_check(stdout)
        if self.min_iterations is not None and self.min_iterations > self.min_iterations_old:
            self.min_iterations_old  = self.min_iterations
        self.RunOutputWindow.appendPlainText(stdout)
    def handle_stderr():
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.RunOutputWindow.appendPlainText(stderr)
    def handle_state(state):
        states = {
            QProcess.NotRunning: 'Finished',
            QProcess.Starting: 'Starting Fierro',
            QProcess.Running: 'Running Fierro',
        }
        self.state_name = states[state]
        self.RunOutputWindow.appendPlainText(f"{self.state_name}")
        
    # Terminate the solver
    self.terminate = False
    def kill_homogenization():
        self.terminate = True
        self.p.terminate()
        self.RunOutputWindow.appendPlainText("TERMINATED")
    self.BKillEVPFFT2.clicked.connect(kill_homogenization)
    
    # Run homogenization simulations (6 in total for generating orthotropic properties)
    def run_homogenization():
        for BC_index in range(6):
            if self.terminate:
                print("TERMINATED")
                self.terminate = False
                break
                
            single_EVPFFT(BC_index)
            self.p.waitForStarted()
            while self.p != None:
                QApplication.processEvents()
                
            # Check if tolerance is met
            if self.INHAutomatic.isChecked() and (self.IterationError == True or self.ConvergenceError == True):
                break

            # Generate Homogenized Elastic Constants
            self.THomogenization.setEnabled(True)
            with open(os.path.join(self.simulation_directory, 'str_str.out'), newline='') as f:
                reader = csv.reader(f)
                self.ss_data = list(reader)
            s11 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            s22 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            s33 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            s12 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            s13 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            s23 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            e11 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            e22 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            e33 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            e12 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            e13 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            e23 = [0 for i in range(int(self.INNumberOfSteps.text())+1)]
            for i in range(int(self.INNumberOfSteps.text())):
                s11[i+1] = float(self.ss_data[i+1][6])
                s22[i+1] = float(self.ss_data[i+1][7])
                s33[i+1] = float(self.ss_data[i+1][8])
                s12[i+1] = float(self.ss_data[i+1][11])
                s13[i+1] = float(self.ss_data[i+1][10])
                s23[i+1] = float(self.ss_data[i+1][9])
                e11[i+1] = float(self.ss_data[i+1][0])
                e22[i+1] = float(self.ss_data[i+1][1])
                e33[i+1] = float(self.ss_data[i+1][2])
                e12[i+1] = float(self.ss_data[i+1][5])
                e13[i+1] = float(self.ss_data[i+1][4])
                e23[i+1] = float(self.ss_data[i+1][3])
            # Tension in the x-direction
            if BC_index == 0:
                self.HE11 = np.polyfit(e11,s11,1)
                self.HNU12 = np.polyfit(e11,e22,1)
                self.HNU13 = np.polyfit(e11,e33,1)
                self.THomogenization.setItem(0,0,QTableWidgetItem(str(self.HE11[0])))
                self.THomogenization.setItem(3,0,QTableWidgetItem(str(-self.HNU12[0])))
#                self.THomogenization.setItem(5,0,QTableWidgetItem(str(-self.HNU13[0])))
            # Tension in the y-direction
            if BC_index == 1:
                self.HE22 = np.polyfit(e22,s22,1)
                self.HNU21 = np.polyfit(e22,e11,1)
                self.HNU23 = np.polyfit(e22,e33,1)
                self.THomogenization.setItem(1,0,QTableWidgetItem(str(self.HE22[0])))
#                self.THomogenization.setItem(4,0,QTableWidgetItem(str(-self.HNU21[0])))
                self.THomogenization.setItem(5,0,QTableWidgetItem(str(-self.HNU23[0])))
            # Tension in the z-direction
            if BC_index == 2:
                self.HE33 = np.polyfit(e33,s33,1)
                self.HNU31 = np.polyfit(e33,e11,1)
                self.HNU32 = np.polyfit(e33,e22,1)
                self.THomogenization.setItem(2,0,QTableWidgetItem(str(self.HE33[0])))
                self.THomogenization.setItem(4,0,QTableWidgetItem(str(-self.HNU31[0])))
#                self.THomogenization.setItem(8,0,QTableWidgetItem(str(-self.HNU32[0])))
            # Shear in the xy-direction
            if BC_index == 3:
                self.HG12 = np.polyfit(np.multiply(e12,2.),s12,1)
                self.THomogenization.setItem(6,0,QTableWidgetItem(str(self.HG12[0])))
            # Shear in the xz-direction
            if BC_index == 4:
                self.HG13 = np.polyfit(np.multiply(e13,2.),s13,1)
                self.THomogenization.setItem(7,0,QTableWidgetItem(str(self.HG13[0])))
            # Shear in the yz-direction
            if BC_index == 5:
                self.HG23 = np.polyfit(np.multiply(e23,2.),s23,1)
                self.THomogenization.setItem(8,0,QTableWidgetItem(str(self.HG23[0])))
                
                # Write output file of all homogenized constants
                self.homogenized_constants = os.path.join(self.evpfft_dir, 'MaterialConstants', "MaterialConstants_" + f'{self.mat_name}' + ".txt")
                if not os.path.exists(os.path.join(self.evpfft_dir, 'MaterialConstants')):
                        os.makedirs(os.path.join(self.evpfft_dir, 'MaterialConstants'))
                with open(self.homogenized_constants, "w") as file:
                    file.write(f'Exx, {self.THomogenization.item(0,0).text()}, {self.THomogenization.item(0,1).text()}\n')
                    file.write(f'Eyy, {self.THomogenization.item(1,0).text()}, {self.THomogenization.item(1,1).text()}\n')
                    file.write(f'Ezz, {self.THomogenization.item(2,0).text()}, {self.THomogenization.item(2,1).text()}\n')
                    file.write(f'NUxy, {self.THomogenization.item(3,0).text()}, \n')
                    file.write(f'NUyz, {self.THomogenization.item(5,0).text()}, \n')
                    file.write(f'NUzx, {self.THomogenization.item(4,0).text()}, \n')
                    file.write(f'Gxy, {self.THomogenization.item(6,0).text()}, {self.THomogenization.item(6,1).text()}\n')
                    file.write(f'Gyz, {self.THomogenization.item(8,0).text()}, {self.THomogenization.item(8,1).text()}\n')
                    file.write(f'Gzx, {self.THomogenization.item(7,0).text()}, {self.THomogenization.item(7,1).text()}\n')
        
        # update convergence parameters if convergence isn't met
        if self.IterationError == True:
            self.IterationError = False
            self.min_iterations_old = 0
            if self.INHAutomatic.isChecked():
                new_tol = float(self.INErrorTolerance.text())/10
                if new_tol < 1e-15:
                    warning_message("ERROR: Convergence could not be achieved.")
                else:
                    self.INErrorTolerance.setText(str(new_tol))
                    run_homogenization()
            elif self.INHManual.isChecked():
                warning_message("WARNING: It is recomended that you decrease your error tolerance. Solution convergence was NOT achieved.")
        if self.ConvergenceError == True:
            self.ConvergenceError = False
            self.min_iterations_old = 0
            if self.INHAutomatic.isChecked():
                new_iterations = int(self.INMaxIterations.text()) + 100
                if new_iterations > 1000:
                    warning_message("ERROR: Convergence could not be achieved.")
                else:
                    self.INMaxIterations.setText(str(new_iterations))
                    run_homogenization()
            elif self.INHManual.isChecked():
                warning_message("WARNING: It is recomended that you increase your maximum number of iterations. Solution convergence was NOT achieved.")
    
    # Connect homogenization run button
    self.p = None
    self.run = 0
    def run_click():
        if self.INSingleJob.isChecked():
            self.mat_name = self.TParts.item(0,0).text()
            self.job_name = "simulation_files"
            run_homogenization()
            self.run = 1
        elif self.INBatchJob.isChecked():
            # loop over all files
            for i in range(len(self.file_paths)):
                # Output information
                self.mat_name = f'{self.file_names[i]}'
                self.INFileNumber.setText(f'{i+1}/{self.number_of_files}')
                self.INHomogenizationBatchFile.setText(f'{self.file_paths[i]}')
                if self.INSaveMaterialFiles.isChecked():
                    self.INPartName.setText("Homogenization")
                    self.job_name = "simulation_files"
                elif self.INSaveAllFiles.isChecked():
                    self.INPartName.setText(f'{self.file_names[i]}')
                    self.job_name = f'{self.file_names[i]}'
                    
                # Get file type for image stacks
                if self.folders:
                    self.file_type = os.path.splitext(os.listdir(self.file_paths[i])[1])[1]
                    self.is_file_names = sorted([os.path.join(self.file_paths[i], f) for f in os.listdir(self.file_paths[i]) if f.endswith(self.file_type)])
            
                # Delete everything in paraview window
                sources = paraview.simple.GetSources()
                for source in sources.values():
                    paraview.simple.Delete(source)
                    
                # Reload geometry
                self.in_file_path = self.file_paths[i]
                Upload_Batch_Geometry(self)
                
                # Reset solver settings...
                if i == 0:
                    steps = self.INNumberOfSteps.text()
                    tol = self.INErrorTolerance.text()
                    iter = self.INMaxIterations.text()
                if i > 0:
                    self.INNumberOfSteps.setText(steps)
                    self.INErrorTolerance.setText(tol)
                    self.INMaxIterations.setText(iter)
                            
                # Run homogenization
                run_homogenization()
                
        # Save job directory
        self.INHomogenizationJobDir.setText(f'{self.working_directory}')
    self.BRunEVPFFT2.clicked.connect(run_click)
    
    # Select geometry files for batch run
    def batch_geometry():
        # Ask user to select a folder
        selected_folder = QFileDialog.getExistingDirectory(None, "Select Folder")

        # Get list of all items in the folder
        items = os.listdir(selected_folder)
        
        # Initialize lists to store folders and file types
        self.folders = []
        file_types = set()
        
        # Iterate through items
        for item in items:
            item_path = os.path.join(selected_folder, item)
            if os.path.isdir(item_path):
                self.folders.append(item)
            else:
                # Get file extension
                _, extension = os.path.splitext(item)
                if extension:
                    file_types.add(extension.lower())
        
        # Check if there are folders inside the selected folder
        if self.folders:
            self.file_paths = sorted([os.path.join(selected_folder, f) for f in items if not f.startswith('.')])
            self.file_names = sorted([os.path.splitext(f)[0] for f in items if not f.startswith('.')])
        else:
            # File information
            if not file_types:
                warning_message("WARNING: the file type could not be determined for the batch job")
            self.file_type = list(file_types)[0]
            self.file_paths = sorted([os.path.join(selected_folder, f) for f in items if f.endswith(self.file_type)])
            self.file_names = sorted([os.path.splitext(f)[0] for f in items if f.endswith(self.file_type)])

        # Output Information
        self.number_of_files = len(self.file_paths)
        self.INFileNumber.setText(f'1/{self.number_of_files}')
        self.INHomogenizationBatchFile.setText(f'{self.file_paths[0]}')
    self.BSelectGeometryFiles.clicked.connect(batch_geometry)
    
    # Write input files
    def Write_Input_Files():
        # Ask user to select a folder
        selected_directory = QFileDialog.getExistingDirectory(None, "Select Folder")
        # Create input files
        self.ELASTIC_PARAMETERS_0 = os.path.join(selected_directory, 'elastic_parameters_0.txt')
        self.ELASTIC_PARAMETERS_1 = os.path.join(selected_directory, 'elastic_parameters_1.txt')
        self.PLASTIC_PARAMETERS = os.path.join(selected_directory, 'plastic_parameters.txt')
        self.EVPFFT_INPUT = os.path.join(selected_directory, 'evpfft_lattice_input.txt')
        for BC_index in range(6):
            self.EVPFFT_INPUT = os.path.join(selected_directory, f'evpfft_lattice_input_{BC_index}.txt')
            Homogenization_WInput(self,BC_index)
        # Write a readme file
        readme = os.path.join(selected_directory, 'README.txt')
        wreadme = open(readme,"w")
        about = 'ABOUT: \n' \
                 'elastic_parameters_#.txt -> these files contain the elastic material properties\n' \
                 'plastic_parameters.txt -> this file contains null plastic material properties\n' \
                 'input_#.txt -> this file is your input file for each of the 6 homogenization steps\n' \
                 '            -> ** Ensure you go into these files and properly specify file paths\n'
        wreadme.write(about)
        run = 'TO RUN: \n' \
               'flags -> Please set the flags: export OMP_PROC_BIND=spread and export OMP_NUM_THREADS=1\n' \
               'parallel run -> (.txt input): mpirun -np # evpfft -f input_#.txt\n' \
               '             -> (.vtk input): mpirun -np # evpfft -f input_#.txt -m 2\n' \
               'serial run -> (.txt input): evpfft -f input_#.txt\n' \
               '           -> (.vtk input): evpfft -f input_#.txt -m 2\n' \
               'help -> evpfft --help'
        wreadme.write(run)
        wreadme.close()
    self.BHWriteFiles.clicked.connect(Write_Input_Files)
    
    # Select to run locally or write the input files
    self.INHRunLocally.toggled.connect(lambda: self.HomogenizationRunWrite.setCurrentIndex(0))
    self.INHWriteFiles.toggled.connect(lambda: self.HomogenizationRunWrite.setCurrentIndex(1))
    
    # Select single job or batch job
    self.INSingleJob.toggled.connect(lambda: self.HomogenizationBatch.setCurrentIndex(0))
    self.INBatchJob.toggled.connect(lambda: self.HomogenizationBatch.setCurrentIndex(1))
    
    # Select serial or parallel run type
    self.INHSerial.toggled.connect(lambda: self.HomogenizationRunType.setCurrentIndex(0))
    self.INHParallel.toggled.connect(lambda: self.HomogenizationRunType.setCurrentIndex(1))
    
    # Set number of mpi ranks based on the number of cores that are avaliable on the system
    num_cores = os.cpu_count()
    self.INmpiRanks.setMaximum(num_cores)
    
#    # Upload HDF5
#    def upload_hdf5_click():
#        global hdf5_filename
#        hdf5_filename = QFileDialog.getOpenFileName(
#            filter="HDF5 File (*.h5, *.hdf5, *.dream3d)",
#        )
#        ReadHDF5(hdf5_filename, "data")
    #self.BUploadHDF5.clicked.connect(upload_hdf5_click)
    
    # Preview Results
    def preview_results_click():
        self.vis = 1
        # Remove all objects from window view
        SetActiveView(self.render_view)
        renderer = self.render_view.GetRenderer()
        renderer.RemoveAllViewProps()
        self.render_view.Update()
        self.render_view.StillRender()
        try:
            self.display
        except:
            False
        else:
            self.display.SetScalarBarVisibility(self.render_view, False)

        # Display .xdmf data
        output_name = str(self.INBCFile.currentText())
        output_parts = output_name.split()
#        file_name = "micro_state_timestep_" + str(self.INNumberOfSteps.text()) + ".xdmf"
        file_name = f"MicroState_{int(self.INNumberOfSteps.text()):05}.pvtu"
        try:
            self.working_directory
        except:
            self.working_directory = self.INHomogenizationJobDir.text()
        else:
            False
        self.output_directory = os.path.join(self.working_directory, output_parts[0], output_parts[1], "pvtu", file_name)
#        self.output_directory = os.path.join(self.working_directory, output_parts[0], output_parts[1], file_name)
#        self.results_reader = paraview.simple.XDMFReader(FileNames=self.output_directory)
        self.results_reader = paraview.simple.XMLPartitionedUnstructuredGridReader(FileName=self.output_directory)

        # Apply warp filter as long as result wasn't run using more than 1 mpi rank
        if self.INHSerial.isChecked():
            # Enable deformation scale factor
            self.INHDeform.setEnabled(True)
            # Calculate transform filter stuff
            self.results_reader.UpdatePipeline()
            # Access the output from the reader
            output_data = self.results_reader.GetClientSideObject().GetOutput()
            # Get point data
            points = output_data.GetPoints()
            num_points = points.GetNumberOfPoints()
            # Extract coordinates
            coords = []
            for i in range(num_points):
                coord = points.GetPoint(i)
                coords.append(coord)
            # Extract differences
            diffX = []
            diffY = []
            diffZ = []
            count = 0
            for k in range(int(self.TParts.item(0,9).text())+1):
                for j in range(int(self.TParts.item(0,8).text())+1):
                    for i in range(int(self.TParts.item(0,7).text())+1):
                        diffX.append(coords[count][0]-(float(i)+0.5))
                        diffY.append(coords[count][1]-(float(j)+0.5))
                        diffZ.append(coords[count][2]-(float(k)+0.5))
                        count += 1
            # Create a new array for displacements
            displacement_array = vtk.vtkFloatArray()
            displacement_array.SetName("Displacement")
            displacement_array.SetNumberOfComponents(3)
            displacement_array.SetNumberOfTuples(num_points)
            for i in range(num_points):
                displacement_array.SetTuple3(i, diffX[i], diffY[i], diffZ[i])
            # Add the displacement array to the output
            output_data.GetPointData().AddArray(displacement_array)
            # Create a new source with the updated data
            temp_source = paraview.simple.TrivialProducer()
            temp_source.GetClientSideObject().SetOutput(output_data)
            # Get initial scale factor
            if self.INHDeform.value() == 0:
                # get the largest undeformed dimension
                udims = [float(self.TParts.item(0,7).text()),float(self.TParts.item(0,8).text()),float(self.TParts.item(0,9).text())]
                udimax = max(udims)
                # get the max difference in each direction and the
                dmaxX = max(abs(num) for num in diffX)
                dmaxY = max(abs(num) for num in diffY)
                dmaxZ = max(abs(num) for num in diffZ)
                dmax = max([dmaxX,dmaxY,dmaxZ])
                # get initial scale factor so that it the deformation is 10% of the largest dimension
                initial_scale_factor = (0.1*udimax)/dmax
                self.INHDeform.setValue(initial_scale_factor)
            # Scale the displacements
            scale_factor = self.INHDeform.value()  # Adjust this value to change the scaling
            scale_filter = paraview.simple.Calculator(Input=temp_source)
            scale_filter.ResultArrayName = 'ScaledDisplacement'
            scale_filter.Function = f'Displacement * {scale_factor}'
            paraview.simple.UpdatePipeline()
        else:
            # Disable deformation scale factor
            self.INHDeform.setEnabled(False)
        
        # Apply threshold and transform filters to view certain phase id's
        in_file_path = self.TParts.item(0,10).text()
        file_type = os.path.splitext(in_file_path)[1].lower()
        if ".txt" in file_type:
            self.INResultRegion.setCurrentIndex(2)
            self.INResultRegion.setEnabled(False)
        if hasattr(self, 'threshold'):
            paraview.simple.Hide(self.threshold)
        if str(self.INResultRegion.currentText()) == "Part + Void":
#                paraview.simple.SetDisplayProperties(Representation="Surface")
            # Apply warp filter as long as result wasn't run using more than 1 mpi rank
            if self.INHSerial.isChecked():
                # Warp Filter
                self.threshold = paraview.simple.WarpByVector(Input=scale_filter)
                self.threshold.Vectors = ['POINTS', 'ScaledDisplacement']
                # Display
                self.display = Show(self.threshold, self.render_view)
            else:
                self.display = Show(self.results_reader, self.render_view)
            paraview.simple.UpdatePipeline()
        elif str(self.INResultRegion.currentText()) == "Part":
            # Apply warp filter as long as result wasn't run using more than 1 mpi rank
            if self.INHSerial.isChecked():
                # Warp Filter
                self.transform = paraview.simple.WarpByVector(Input=scale_filter)
                self.transform.Vectors = ['POINTS', 'ScaledDisplacement']
                # Threshold Filter
                self.threshold = paraview.simple.Threshold(registrationName='results_threshold', Input = self.transform, Scalars = "phase_id", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 2, LowerThreshold = 1, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
            else:
                # Threshold Filter
                self.threshold = paraview.simple.Threshold(registrationName='results_threshold', Input = self.results_reader, Scalars = "phase_id", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 2, LowerThreshold = 1, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
            # Display
            self.display = Show(self.threshold, self.render_view)
            paraview.simple.UpdatePipeline()
        elif str(self.INResultRegion.currentText()) == "Void":
            # Apply warp filter as long as result wasn't run using more than 1 mpi rank
            if self.INHSerial.isChecked():
                # Warp Filter
                self.transform = paraview.simple.WarpByVector(Input=scale_filter)
                self.transform.Vectors = ['POINTS', 'ScaledDisplacement']
                # Threshold Filter
                self.threshold = paraview.simple.Threshold(registrationName='results_threshold', Input = self.transform, Scalars = "phase_id", ThresholdMethod = "Below Lower Threshold", UpperThreshold = 2, LowerThreshold = 1, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
            else:
                # Threshold Filter
                self.threshold = paraview.simple.Threshold(registrationName='results_threshold', Input = self.results_reader, Scalars = "phase_id", ThresholdMethod = "Below Lower Threshold", UpperThreshold = 2, LowerThreshold = 1, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
            # Display
            self.display = Show(self.threshold, self.render_view)
            paraview.simple.UpdatePipeline()
#        else:
#            self.INResultRegion.setEnabled(False)
#            paraview.simple.SetDisplayProperties(Representation="Surface")
#            self.display = Show(self.results_reader, self.render_view)
        
        # Color by the selected variable
        selected_variable = str(self.INPreviewResults.currentText())
        paraview.simple.ColorBy(self.display, ('CELLS', selected_variable))
        vmstressLUT = paraview.simple.GetColorTransferFunction(selected_variable)
        self.display.SetScalarBarVisibility(self.render_view, True)

        # View Results
        self.render_view.ResetCamera()
        self.render_view.StillRender()
        self.activated = 1
    
    # Restart deformation scale factor
    def restart_scale():
        self.INHDeform.setValue(0)
    self.INBCFile.currentIndexChanged.connect(restart_scale)
    
    # Show results immediately when postprocessing tab is pressed
    self.activated = 0
    def show_results():
        index = self.NavigationMenu.currentIndex()
        if index == 8 and self.activated == 0:
            self.INBCFile.currentIndexChanged.connect(preview_results_click)
            self.INPreviewResults.currentIndexChanged.connect(preview_results_click)
            self.INResultRegion.currentIndexChanged.connect(preview_results_click)
            self.INHDeform.valueChanged.connect(preview_results_click)
        if index == 8 and self.run == 1 and "Homogenization" in self.INSelectPostprocessing.currentText():
            preview_results_click()
    self.NavigationMenu.currentChanged.connect(show_results)
    
    # Go from results back to back to original input
    self.vis = 0
    def show_original():
        index = self.NavigationMenu.currentIndex()
        if index < 8 and self.vis == 1 and "Homogenization" in self.INSelectPostprocessing.currentText():
            # reset view
            self.INHDeform.setValue(.1)
            self.INResultRegion.setCurrentIndex(0)
            paraview.simple.ColorBy(self.display, ('CELLS', 'grain_id'))
            # View Results
            self.display.SetScalarBarVisibility(self.render_view, False)
            self.render_view.ResetCamera()
            self.render_view.StillRender()
            self.vis = 0
    self.NavigationMenu.currentChanged.connect(show_original)
    
    # Open Paraview
    def open_paraview_click():
        command = ["paraview", self.output_directory]
        subprocess.Popen(command)
    self.BOpenParaview.clicked.connect(open_paraview_click)
    
#    # Stress vs Strain Plot
#    def plot_ss_click():
#        # Get the stress-strain data
#        try:
#            self.Plot.figure
#        except:
#            with open("str_str.out", newline='') as f:
#                reader = csv.reader(f)
#                self.ss_data = list(reader)
#            x = [0 for i in range(int(self.INNumberOfSteps.text()))]
#            y = [0 for i in range(int(self.INNumberOfSteps.text()))]
#            for i in range(int(self.INNumberOfSteps.text())):
#                if str(self.INPlotSS.currentText()) == 'S11 vs E11':
#                    xcol = 0
#                    ycol = 6
#                elif str(self.INPlotSS.currentText()) == 'S22 vs E22':
#                    xcol = 1
#                    ycol = 7
#                elif str(self.INPlotSS.currentText()) == 'S33 vs E33':
#                    xcol = 2
#                    ycol = 8
#                x[i] = float(self.ss_data[i+1][xcol])
#                y[i] = float(self.ss_data[i+1][ycol])
#            # Plot data
#            self.Plot.figure = Figure()
#            self.Plot.ax = self.Plot.figure.add_subplot()
#            print("(",x[1],",",y[1],")")
#            self.Plot.ax.plot(x,y)
#            self.Plot.ax.set_xlabel('STRAIN')
#            self.Plot.ax.set_ylabel('STRESS')
#            self.Plot.figure.tight_layout()
#            # Display plot and toolbar
#            layout = QVBoxLayout()
#            self.Plot.setLayout(layout)
#            self.Plot.canvas = FigureCanvasQTAgg(self.Plot.figure)
#            layout.addWidget(self.Plot.canvas)
#            self.toolbar = NavigationToolbar2QT(self.Plot.canvas,self.Plot)
#            layout.addWidget(self.toolbar)
#        else:
#            self.timer = QTimer()
#            self.timer.setInterval(100)
#            self.timer.timeout.connect(update_plot)
#            self.timer.start()
#    def update_plot():
#        if self.run_cnt > 1:
#            with open("str_str.out", newline='') as f:
#                reader = csv.reader(f)
#                self.ss_data = list(reader)
#        x = [0 for i in range(int(self.INNumberOfSteps.text()))]
#        y = [0 for i in range(int(self.INNumberOfSteps.text()))]
#        for i in range(int(self.INNumberOfSteps.text())):
#            if str(self.INPlotSS.currentText()) == 'S11 vs E11':
#                xcol = 0
#                ycol = 6
#            elif str(self.INPlotSS.currentText()) == 'S22 vs E22':
#                xcol = 1
#                ycol = 7
#            elif str(self.INPlotSS.currentText()) == 'S33 vs E33':
#                xcol = 2
#                ycol = 8
#            x[i] = float(self.ss_data[i+1][xcol])
#            y[i] = float(self.ss_data[i+1][ycol])
#        self.Plot.ax.cla()
#        self.Plot.ax.plot(x,y)
#        self.Plot.ax.set_xlabel('STRAIN')
#        self.Plot.ax.set_ylabel('STRESS')
#        self.Plot.figure.tight_layout()
#        self.Plot.canvas.draw()
#        self.timer.stop()
#    self.BPlotSS.clicked.connect(plot_ss_click)
#    self.BPlotSS.clicked.connect(lambda: self.OutputWindows.setCurrentIndex(1))


