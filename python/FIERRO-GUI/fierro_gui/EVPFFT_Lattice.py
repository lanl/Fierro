from PySide6.QtWidgets import (QTableWidgetItem, QMessageBox, QApplication)
from PySide6.QtCore import (QCoreApplication, QProcess)
import re
import csv
import numpy as np
import subprocess
#import paraview.simple as paraview.simple
from paraview.simple import *
from EVPFFT_Lattice_WInput import *
from DeveloperInputs import *

# ==============================================
# ======= EVPFFT SOLVER LATTICE PIPELINE =======
# ==============================================

# Warning Message Popup
def warning_message(msg):
    message = QMessageBox()
    message.setText(msg)
    message.exec()

def EVPFFT_Lattice(self):
    # Connect tab buttons to settings windows
    self.BImportPart.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(2))
    self.BDefineMaterial.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(4))
    self.BApplyBC.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(5))
    self.BSolverSettings.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(6))
    self.BViewResults.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(7))
    self.BGlobalMesh.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(1))
    
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
            # Remove all objects from window view
            SetActiveView(self.render_view)
            renderer = self.render_view.GetRenderer()
            renderer.RemoveAllViewProps()
            self.render_view.Update()
            self.render_view.StillRender()
            
            # Show void region only
            self.threshold = paraview.simple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Below Lower Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
            paraview.simple.Show(self.threshold, self.render_view)
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
            self.threshold = paraview.simple.Threshold(Input = self.vtk_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
            paraview.simple.Show(self.threshold, self.render_view)
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
    
    # Single Run of EVPFFT
    self.run_cnt = 0
    def single_EVPFFT(BC_index):
        if self.p == None:
            EVPFFT_Lattice_WInput(self,BC_index)
            executable_path = fierro_evpfft_exe
            arguments = ["-f", self.EVPFFT_INPUT, "-m", "2"]
            self.p = QProcess()
            self.p.readyReadStandardOutput.connect(handle_stdout)
            self.p.readyReadStandardError.connect(handle_stderr)
            self.p.stateChanged.connect(handle_state)
            self.p.finished.connect(process_finished)
            self.p.start(executable_path, arguments)
            self.progress_re = re.compile("       Current  Time  STEP = (\\d+)")
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
            QProcess.Starting: 'Starting Fierro',
            QProcess.Running: 'Running Fierro',
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
    
    # Preview Results
    def preview_results_click():
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
        self.results_reader = paraview.simple.XDMFReader(FileNames="micro_state_timestep_10.xdmf")
        paraview.simple.SetDisplayProperties(Representation="Surface")
        self.display = Show(self.results_reader, self.render_view)
        
        # Color by the selected variable
        selected_variable = str(self.INPreviewResults.currentText())
        paraview.simple.ColorBy(self.display, ('CELLS', selected_variable))
        vmstressLUT = paraview.simple.GetColorTransferFunction(selected_variable)
        self.display.SetScalarBarVisibility(self.render_view, True)
        self.render_view.ResetCamera()
        self.render_view.StillRender()
    self.BPreviewResults.clicked.connect(preview_results_click)
    self.BPreviewResults.clicked.connect(lambda: self.OutputWindows.setCurrentIndex(0))
    
    # Open Paraview
    def open_paraview_click():
        command = ["paraview", "micro_state_timestep_10.xdmf"]
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


