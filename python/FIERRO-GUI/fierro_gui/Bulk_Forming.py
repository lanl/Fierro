from PySide6.QtWidgets import (QTableWidgetItem, QMessageBox, QApplication, QFileDialog)
from PySide6.QtCore import (QCoreApplication, QProcess, Qt)
from PySide6.QtGui import (QColor, QBrush)
import numpy as np
import re
import os
import subprocess
from paraview.simple import *
from importlib import reload
from Bulk_Forming_WInput import *
import DeveloperInputs

# Warning Message Popup
def warning_message(msg):
    message = QMessageBox()
    message.setText(msg)
    message.exec()

# Bulk Forming Functionalities
def Bulk_Forming(self):
    # Define Either Solid or Gas Material Input
    def material_type_2():
        if str(self.INSolidGas_2.currentText()) == 'Solid':
            self.INMaterialType_2.clear()
            self.INMaterialType_2.addItem(QCoreApplication.translate("MainWindow", u"Isotropic", None))
            self.INMaterialType_2.addItem(QCoreApplication.translate("MainWindow", u"Transversely Isotropic", None))
            self.INMaterialType_2.addItem(QCoreApplication.translate("MainWindow", u"Orthotropic", None))
            self.INMaterialType_2.addItem(QCoreApplication.translate("MainWindow", u"Anisotropic", None))
            self.MaterialTypeTool_2.setEnabled(True)
        if str(self.INSolidGas_2.currentText()) == 'Gas':
            self.INMaterialType_2.clear()
            self.INMaterialType_2.addItem(QCoreApplication.translate("MainWindow", u"Ideal Gas", None))
            self.MaterialTypeTool_2.setEnabled(False)
    self.INSolidGas_2.currentIndexChanged.connect(material_type_2)
    
    # Change inputs based on material input class
    def material_class_2():
        if str(self.INMaterialType_2.currentText()) == 'Isotropic':
            self.MaterialTypeTool_2.setCurrentIndex(0)
        if str(self.INMaterialType_2.currentText()) == 'Transversely Isotropic':
            self.MaterialTypeTool_2.setCurrentIndex(1)
        if str(self.INMaterialType_2.currentText()) == 'Orthotropic':
            self.MaterialTypeTool_2.setCurrentIndex(3)
        if str(self.INMaterialType_2.currentText()) == 'Anisotropic':
            self.MaterialTypeTool_2.setCurrentIndex(2)
        if str(self.INMaterialType_2.currentText()) == 'Ideal Gas':
            self.MaterialTypeTool_2.setCurrentIndex(4)
    self.INMaterialType_2.currentIndexChanged.connect(material_class_2)
    
    # Define material properties using predefined values
    def predefined_materials():
        if "Custom" in self.INMaterialDefinition.currentText():
            # Clear elastic parameters
            self.INYoungsModulus_2.clear()
            self.INPoissonsRatio_2.clear()
            self.INEip_2.clear()
            self.INNUip_2.clear()
            self.INEop_2.clear()
            self.INNUop_2.clear()
            self.INGop_2.clear()
            self.INEx_2.clear()
            self.INEy_2.clear()
            self.INEz_2.clear()
            self.INNUxy_2.clear()
            self.INNUxz_2.clear()
            self.INNUyz_2.clear()
            self.INGxy_2.clear()
            self.INGxz_2.clear()
            self.INGyz_2.clear()
            for i in [0,1,2,3,4,5,6]:
                for j in range(i,6):
                    if self.TAnisotropic_2.item(i,j):
                        self.TAnisotropic_2.item(i,j).setText('')
            # Clear plastic parameters
            self.INa.clear()
            self.INb.clear()
            self.INc.clear()
            self.INSlipSystems.clear()
            self.TSlipSystemParameters.setRowCount(0)
            self.INnrsx.clear()
            self.INgamd0x.clear()
            self.INtau0xb.clear()
            self.INtau0xf.clear()
            self.INtau1x.clear()
            self.INthet0.clear()
            self.INthet1.clear()
            self.INhselfx.clear()
            self.INhlatex.clear()
        elif "Single Crystal Copper" in self.INMaterialDefinition.currentText():
            # Define elastic properties
            for i in [0,1,2,3,4,5,6]:
                for j in range(i,6):
                    self.TAnisotropic_2.setItem(i,j,QTableWidgetItem('0.'))
            self.INSolidGas_2.setCurrentIndex(0)
            self.INMaterialType_2.setCurrentIndex(3)
            self.TAnisotropic_2.setItem(0,0,QTableWidgetItem('168400.'))
            self.TAnisotropic_2.setItem(1,1,QTableWidgetItem('168400.'))
            self.TAnisotropic_2.setItem(2,2,QTableWidgetItem('168400.'))
            self.TAnisotropic_2.setItem(0,1,QTableWidgetItem('121400.'))
            self.TAnisotropic_2.setItem(0,2,QTableWidgetItem('121400.'))
            self.TAnisotropic_2.setItem(1,2,QTableWidgetItem('121400.'))
            self.TAnisotropic_2.setItem(3,3,QTableWidgetItem('75400.'))
            self.TAnisotropic_2.setItem(4,4,QTableWidgetItem('75400.'))
            self.TAnisotropic_2.setItem(5,5,QTableWidgetItem('75400.'))
            # Define plastic properties
            self.INa.setText('1.')
            self.INb.setText('1.')
            self.INc.setText('1.')
            items = self.TSlipSystems.findItems('FCC', Qt.MatchExactly | Qt.MatchRecursive)
            if items:
                self.TSlipSystems.setCurrentItem(items[0])
                self.BAddSlipSystem.click()
                self.TSlipSystems.expandItem(items[0])
            self.INnrsx.setText('10')
            self.INgamd0x.setText('1.0')
            self.INtau0xf.setText('9.0')
            self.INtau0xb.setText('9.0')
            self.INtau1x.setText('5.')
            self.INthet0.setText('400.')
            self.INthet1.setText('250.')
            self.INhselfx.setText('1.0')
            self.INhlatex.setText('1.0')
    self.INMaterialDefinition.currentIndexChanged.connect(predefined_materials)
    
    # Add material to the table
    def add_material_2():
        # Elastic properties checks
        if 'Gas' in self.INSolidGas_2.currentText():
            if not self.INMaterialName_2.text():
                warning_message('ERROR: Elastic material definition incomplete')
                return
        if 'Solid' in self.INSolidGas_2.currentText():
            if self.INMaterialType_2.currentText() == 'Isotropic':
                if not self.INYoungsModulus_2.text().strip() or not self.INPoissonsRatio_2.text().strip() or not self.INMaterialName_2.text().strip():
                    warning_message('ERROR: Elastic material definition incomplete')
                    return
            if str(self.INMaterialType_2.currentText()) == 'Transversely Isotropic':
                if not self.INEip_2.text().strip() or not self.INNUip_2.text().strip() or not self.INEop_2.text().strip() or not self.INNUop_2.text().strip() or not self.INGop_2.text().strip() or not self.INMaterialName_2.text().strip():
                    warning_message('ERROR: Elastic material definition incomplete')
                    return
            if str(self.INMaterialType_2.currentText()) == 'Orthotropic':
                if not self.INEx_2.text().strip() or not self.INEy_2.text().strip() or not self.INEz_2.text().strip() or not self.INNUxy_2.text().strip() or not self.INNUxz_2.text().strip() or not self.INNUyz_2.text().strip() or not self.INGxy_2.text().strip() or not self.INGxz_2.text().strip() or not self.INGyz_2.text().strip() or not self.INMaterialName_2.text().strip():
                    warning_message('ERROR: Elastic material definition incomplete')
                    return
            if str(self.INMaterialType_2.currentText()) == 'Anisotropic':
                if not self.TAnisotropic_2.item(0,0).text().strip() or not self.TAnisotropic_2.item(0,1).text().strip() or not self.TAnisotropic_2.item(0,2).text().strip() or not self.TAnisotropic_2.item(0,3).text().strip() or not self.TAnisotropic_2.item(0,4).text().strip() or not self.TAnisotropic_2.item(0,5).text().strip() or not self.TAnisotropic_2.item(1,1).text().strip() or not self.TAnisotropic_2.item(1,2).text().strip() or not self.TAnisotropic_2.item(1,3).text().strip() or not self.TAnisotropic_2.item(1,4).text().strip() or not self.TAnisotropic_2.item(1,5).text().strip() or not self.TAnisotropic_2.item(2,2).text().strip() or not self.TAnisotropic_2.item(2,3).text().strip() or not self.TAnisotropic_2.item(2,4).text().strip() or not self.TAnisotropic_2.item(2,5).text().strip() or not self.TAnisotropic_2.item(3,3).text().strip() or not self.TAnisotropic_2.item(3,4).text().strip() or not self.TAnisotropic_2.item(3,5).text().strip() or not self.TAnisotropic_2.item(4,4).text().strip() or not self.TAnisotropic_2.item(4,5).text().strip() or not self.TAnisotropic_2.item(5,5).text().strip() or not self.INMaterialName_2.text().strip():
                    warning_message('ERROR: Elastic material definition incomplete')
                    return
        # Plastic parameters checks
        if not self.INa.text().strip() or not self.INb.text().strip() or not self.INc.text().strip():
            warning_message('ERROR: crystal axis definition incomplete')
            return
        if self.TSlipSystemParameters.rowCount() == 0:
            warning_message('ERROR: Voce parameters are incomplete')
            return
        for rowc in range(self.TSlipSystemParameters.rowCount()):
            for colc in range(self.TSlipSystemParameters.columnCount()):
                item = self.TSlipSystemParameters.item(rowc, colc)
                if item is None or not item.text().strip():
                    warning_message('ERROR: Voce parameters are incomplete')
                    return
        
        # Assign elastic parameters if all checks pass
        if 'Gas' in self.INSolidGas_2.currentText():
            # Fill out material definition table
            row = self.TMaterials_2.rowCount()
            self.TMaterials_2.insertRow(row)
            self.TMaterials_2.setItem(row, 0, QTableWidgetItem(
                self.INMaterialName_2.text().strip())
            )
            self.TMaterials_2.setItem(row, 1, QTableWidgetItem(
                str(self.INMaterialType_2.currentText()))
            )

            # Add material as an option for material assignment
            self.INMaterial_3.clear()
            for i in range(self.TMaterials_2.rowCount()):
                self.INMaterial_3.addItem(self.TMaterials_2.item(i,0).text())

            # Clear fields
            self.INMaterialName_2.clear()
            
        elif 'Solid' in self.INSolidGas_2.currentText():
            if self.INMaterialType_2.currentText() == 'Isotropic':
                # Calculate Stiffness Matrix
                E = float(self.INYoungsModulus_2.text())
                NU = float(self.INPoissonsRatio_2.text())
                INCalcG = float(self.INYoungsModulus_2.text())/(2*(1+float(self.INPoissonsRatio_2.text())))
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
                self.INYoungsModulus_2.clear()
                self.INPoissonsRatio_2.clear()
                
            if str(self.INMaterialType_2.currentText()) == 'Transversely Isotropic':
                if str(self.INIsotropicPlane_2.currentText()) == 'x-y plane':
                    # Calculate Stiffness Matrix
                    NUip = float(self.INNUip_2.text())
                    NUop = float(self.INNUop_2.text())
                    Eip = float(self.INEip_2.text())
                    Eop = float(self.INEop_2.text())
                    Gop = float(self.INGop_2.text())
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

                    self.INEip_2.clear()
                    self.INNUip_2.clear()
                    self.INEop_2.clear()
                    self.INNUop_2.clear()
                    self.INGop_2.clear()
                    
                if str(self.INIsotropicPlane_2.currentText()) == 'x-z plane':
                    # Calculate Stiffness Matrix
                    NUip = float(self.INNUip_2.text())
                    NUop = float(self.INNUop_2.text())
                    Eip = float(self.INEip_2.text())
                    Eop = float(self.INEop_2.text())
                    Gop = float(self.INGop_2.text())
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

                    self.INEip_2.clear()
                    self.INNUip_2.clear()
                    self.INEop_2.clear()
                    self.INNUop_2.clear()
                    self.INGop_2.clear()
                    
                if str(self.INIsotropicPlane_2.currentText()) == 'y-z plane':
                    # Calculate Stiffness Matrix
                    NUip = float(self.INNUip_2.text())
                    NUop = float(self.INNUop_2.text())
                    Eip = float(self.INEip_2.text())
                    Eop = float(self.INEop_2.text())
                    Gop = float(self.INGop_2.text())
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

                    self.INEip_2.clear()
                    self.INNUip_2.clear()
                    self.INEop_2.clear()
                    self.INNUop_2.clear()
                    self.INGop_2.clear()

            if str(self.INMaterialType_2.currentText()) == 'Orthotropic':
                # Calculate Stiffness Matrix
                S11 = 1/float(self.INEx_2.text())
                S12 = -float(self.INNUxy_2.text())/float(self.INEx_2.text())
                S13 = -float(self.INNUxz_2.text())/float(self.INEx_2.text())
                S22 = 1/float(self.INEy_2.text())
                S23 = -float(self.INNUyz_2.text())/float(self.INEy_2.text())
                S33 = 1/float(self.INEz_2.text())
                S44 = 1/float(self.INGyz_2.text())
                S55 = 1/float(self.INGxz_2.text())
                S66 = 1/float(self.INGxy_2.text())
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

                self.INEx_2.clear()
                self.INEy_2.clear()
                self.INEz_2.clear()
                self.INNUxy_2.clear()
                self.INNUxz_2.clear()
                self.INNUyz_2.clear()
                self.INGxy_2.clear()
                self.INGxz_2.clear()
                self.INGyz_2.clear()
                        
            # Fill out material definition table
            row = self.TMaterials_2.rowCount()
            self.TMaterials_2.insertRow(row)
            self.TMaterials_2.setItem(row, 0, QTableWidgetItem(
                self.INMaterialName_2.text())
            )
            if str(self.INMaterialType_2.currentText()) == 'Anisotropic':
                # Fill out material definition table
                self.TMaterials_2.setItem(row, 1, QTableWidgetItem(
                    str(self.INMaterialType_2.currentText()))
                )
                k = 2
                for i in range(7):
                    for j in range(i,6):
                        self.TMaterials_2.setItem(
                            row, k, QTableWidgetItem(self.TAnisotropic_2.item(i,j).text())
                        )
                        self.TAnisotropic_2.item(i,j).setText('')
                        k += 1
                self.INMaterialName_2.clear()
            else:
                if str(self.INMaterialType_2.currentText()) == 'Transversely Isotropic':
                    self.TMaterials_2.setItem(row, 1, QTableWidgetItem(
                        str(self.INMaterialType_2.currentText() + ' ' + self.INIsotropicPlane_2.currentText()))
                    )
                else:
                    self.TMaterials_2.setItem(row, 1, QTableWidgetItem(
                        str(self.INMaterialType_2.currentText()))
                    )
                self.TMaterials_2.setItem(
                    row, 2, QTableWidgetItem(str(C11))
                )
                self.TMaterials_2.setItem(
                    row, 3, QTableWidgetItem(str(C12))
                )
                self.TMaterials_2.setItem(
                    row, 4, QTableWidgetItem(str(C13))
                )
                self.TMaterials_2.setItem(
                    row, 8, QTableWidgetItem(str(C22))
                )
                self.TMaterials_2.setItem(
                    row, 9, QTableWidgetItem(str(C23))
                )
                self.TMaterials_2.setItem(
                    row, 13, QTableWidgetItem(str(C33))
                )
                self.TMaterials_2.setItem(
                    row, 17, QTableWidgetItem(str(C44))
                )
                self.TMaterials_2.setItem(
                    row, 20, QTableWidgetItem(str(C55))
                )
                self.TMaterials_2.setItem(
                    row, 22, QTableWidgetItem(str(C66))
                )
                for i in [5,6,7,10,11,12,14,15,16,18,19,21]:
                   self.TMaterials_2.setItem(row, i, QTableWidgetItem('0'))
                   
        # Add plastic parameters
        self.TMaterials_2.setItem(row,23,QTableWidgetItem(self.INa.text()))
        self.TMaterials_2.setItem(row,24,QTableWidgetItem(self.INb.text()))
        self.TMaterials_2.setItem(row,25,QTableWidgetItem(self.INc.text()))
        for colc in range(self.TSlipSystemParameters.columnCount()):
            for rowc in range(self.TSlipSystemParameters.rowCount()):
                item = self.TSlipSystemParameters.item(rowc, colc)
                if rowc == 0:
                    new_text = item.text()
                else:
                    new_text = new_text + ', ' + item.text()
            self.TMaterials_2.setItem(row,colc+26,QTableWidgetItem(new_text))
            
        # Clear plastic parameters
        self.INa.clear()
        self.INb.clear()
        self.INc.clear()
        self.INSlipSystems.clear()
        self.TSlipSystemParameters.setRowCount(0)
        self.INnrsx.clear()
        self.INgamd0x.clear()
        self.INtau0xb.clear()
        self.INtau0xf.clear()
        self.INtau1x.clear()
        self.INthet0.clear()
        self.INthet1.clear()
        self.INhselfx.clear()
        self.INhlatex.clear()

        # Add material as an option for material assignment
        self.INMaterial_3.clear()
        for i in range(self.TMaterials_2.rowCount()):
            self.INMaterial_3.addItem(self.TMaterials_2.item(i,0).text())

        self.INMaterialName_2.clear()
    self.BAddMaterial_2.clicked.connect(add_material_2)
    
    def delete_material_2():
        current_row = self.TMaterials_2.currentRow()
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
            self.TMaterials_2.removeRow(current_row)
            # delete from material assignment options
            self.INMaterial_3.clear()
            for i in range(self.TMaterials_2.rowCount()):
                self.INMaterial_3.addItem(self.TMaterials_2.item(i,0).text())
    self.BDeleteMaterial_2.clicked.connect(delete_material_2)
            
#    def regenerate_elastic_constants_2():
#        current_row = self.TMaterials_2.currentRow()
#        if current_row < 0:
#            return QMessageBox.warning(QMessageBox(),"Warning","Please select a material from the table")
#            
#        # Define Stiffness Matrix
#        Mstiffness = [[float(self.TMaterials_2.item(current_row,2).text()), float(self.TMaterials_2.item(current_row,3).text()), float(self.TMaterials_2.item(current_row,4).text()),  float(self.TMaterials_2.item(current_row,5).text()), float(self.TMaterials_2.item(current_row,6).text()), float(self.TMaterials_2.item(current_row,7).text())], [float(self.TMaterials_2.item(current_row,3).text()), float(self.TMaterials_2.item(current_row,8).text()), float(self.TMaterials_2.item(current_row,9).text()),  float(self.TMaterials_2.item(current_row,10).text()), float(self.TMaterials_2.item(current_row,11).text()), float(self.TMaterials_2.item(current_row,12).text())], [float(self.TMaterials_2.item(current_row,4).text()), float(self.TMaterials_2.item(current_row,9).text()), float(self.TMaterials_2.item(current_row,13).text()), float(self.TMaterials_2.item(current_row,14).text()), float(self.TMaterials_2.item(current_row,15).text()), float(self.TMaterials_2.item(current_row,16).text())], [float(self.TMaterials_2.item(current_row,5).text()), float(self.TMaterials_2.item(current_row,10).text()), float(self.TMaterials_2.item(current_row,14).text()), float(self.TMaterials_2.item(current_row,17).text()), float(self.TMaterials_2.item(current_row,18).text()), float(self.TMaterials_2.item(current_row,19).text())], [float(self.TMaterials_2.item(current_row,6).text()), float(self.TMaterials_2.item(current_row,11).text()), float(self.TMaterials_2.item(current_row,15).text()), float(self.TMaterials_2.item(current_row,18).text()), float(self.TMaterials_2.item(current_row,20).text()), float(self.TMaterials_2.item(current_row,21).text())], [float(self.TMaterials_2.item(current_row,7).text()), float(self.TMaterials_2.item(current_row,12).text()), float(self.TMaterials_2.item(current_row,16).text()), float(self.TMaterials_2.item(current_row,19).text()), float(self.TMaterials_2.item(current_row,21).text()), float(self.TMaterials_2.item(current_row,22).text())]]
#        if self.TMaterials_2.item(current_row,1).text() == 'Isotropic':
#            Mcompliance = np.linalg.inv(Mstiffness)
#            self.MaterialTypeTool_2.setCurrentIndex(0)
#            self.INMaterialType_2.setCurrentIndex(0)
#            self.INMaterialName_2.clear()
#            self.INYoungsModulus_2.clear()
#            self.INPoissonsRatio_2.clear()
#            E = 1/Mcompliance[0][0]
#            nu = -Mcompliance[0][1]*E
#            self.INMaterialName_2.insert(self.TMaterials_2.item(current_row,0).text())
#            self.INYoungsModulus_2.insert(str(E))
#            self.INPoissonsRatio_2.insert(str(nu))
#        elif 'Transversely Isotropic' in self.TMaterials_2.item(current_row,1).text():
#            Mcompliance = np.linalg.inv(Mstiffness)
#            self.MaterialTypeTool_2.setCurrentIndex(1)
#            self.INMaterialType_2.setCurrentIndex(1)
#            self.INMaterialName_2.clear()
#            self.INEip_2.clear()
#            self.INNUip_2.clear()
#            self.INEop_2.clear()
#            self.INNUop_2.clear()
#            self.INGop_2.clear()
#            if 'x-y plane' in self.TMaterials_2.item(current_row,1).text():
#                Eip = 1/Mcompliance[0][0]
#                nuip = -Mcompliance[0][1]*Eip
#                Eop = 1/Mcompliance[2][2]
#                nuop = -Mcompliance[0][2]*Eop
#                Gop  = 1/Mcompliance[3][3]
#                self.INMaterialName_2.insert(self.TMaterials_2.item(current_row,0).text())
#                self.INEip_2.insert(str(Eip))
#                self.INNUip_2.insert(str(nuip))
#                self.INEop_2.insert(str(Eop))
#                self.INNUop_2.insert(str(nuop))
#                self.INGop_2.insert(str(Gop))
#                self.INIsotropicPlane_2.setCurrentIndex(0)
#            elif 'x-z plane' in self.TMaterials_2.item(current_row,1).text():
#                Eip = 1/Mcompliance[0][0]
#                nuip = -Mcompliance[0][2]*Eip
#                Eop = 1/Mcompliance[1][1]
#                nuop = -Mcompliance[0][1]*Eop
#                Gop  = 1/Mcompliance[3][3]
#                self.INMaterialName_2.insert(self.TMaterials_2.item(current_row,0).text())
#                self.INEip_2.insert(str(Eip))
#                self.INNUip_2.insert(str(nuip))
#                self.INEop_2.insert(str(Eop))
#                self.INNUop_2.insert(str(nuop))
#                self.INGop_2.insert(str(Gop))
#                self.INIsotropicPlane_2.setCurrentIndex(1)
#            elif 'y-z plane' in self.TMaterials_2.item(current_row,1).text():
#                Eip = 1/Mcompliance[1][1]
#                nuip = -Mcompliance[1][2]*Eip
#                Eop = 1/Mcompliance[0][0]
#                nuop = -Mcompliance[0][1]*Eop
#                Gop  = 1/Mcompliance[4][4]
#                self.INMaterialName_2.insert(self.TMaterials_2.item(current_row,0).text())
#                self.INEip_2.insert(str(Eip))
#                self.INNUip_2.insert(str(nuip))
#                self.INEop_2.insert(str(Eop))
#                self.INNUop_2.insert(str(nuop))
#                self.INGop_2.insert(str(Gop))
#                self.INIsotropicPlane_2.setCurrentIndex(2)
#        elif self.TMaterials_2.item(current_row,1).text() == 'Orthotropic':
#            Mcompliance = np.linalg.inv(Mstiffness)
#            self.MaterialTypeTool_2.setCurrentIndex(3)
#            self.INMaterialType_2.setCurrentIndex(2)
#            self.INMaterialName_2.clear()
#            self.INEx_2.clear()
#            self.INEy_2.clear()
#            self.INEz_2.clear()
#            self.INNUxy_2.clear()
#            self.INNUxz_2.clear()
#            self.INNUyz_2.clear()
#            self.INGxy_2.clear()
#            self.INGxz_2.clear()
#            self.INGyz_2.clear()
#            Ex = 1/Mcompliance[0][0]
#            Ey = 1/Mcompliance[1][1]
#            Ez = 1/Mcompliance[2][2]
#            NUxy = -Mcompliance[0][1]*Ex
#            NUxz = -Mcompliance[0][2]*Ex
#            NUyz = -Mcompliance[1][2]*Ey
#            Gxy = 1/Mcompliance[5][5]
#            Gxz = 1/Mcompliance[4][4]
#            Gyz = 1/Mcompliance[3][3]
#            self.INMaterialName_2.insert(self.TMaterials_2.item(current_row,0).text())
#            self.INEx_2.insert(str(Ex))
#            self.INEy_2.insert(str(Ey))
#            self.INEz_2.insert(str(Ez))
#            self.INNUxy_2.insert(str(NUxy))
#            self.INNUxz_2.insert(str(NUxz))
#            self.INNUyz_2.insert(str(NUyz))
#            self.INGxy_2.insert(str(Gxy))
#            self.INGxz_2.insert(str(Gxz))
#            self.INGyz_2.insert(str(Gyz))
#        else:
#            self.MaterialTypeTool_2.setCurrentIndex(2)
#            self.INMaterialType_2.setCurrentIndex(3)
#            self.INMaterialName_2.clear()
#            self.INMaterialName_2.insert(self.TMaterials_2.item(current_row,0).text())
#            k = 2
#            for i in [0,1,2,3,4,5,6]:
#                for j in range(i,6):
#                    self.TAnisotropic_2.item(i,j).setText('')
#                    self.TAnisotropic_2.setItem(
#                        i, j, QTableWidgetItem(self.TMaterials_2.item(current_row,k).text())
#                    )
#                    k += 1
#    self.BRegenElasticConstants_2.clicked.connect(regenerate_elastic_constants_2)

    # View details about a slip system
    def slip_system_details():
        selected_items = self.TSlipSystems.selectedItems()
        selected_item = selected_items[0]
        if selected_item.parent() is not None:
            if '{111}<110>' in selected_item.text(0):
                self.SlipSystemInfo.setCurrentIndex(1)
            self.PlasticProperties.setCurrentIndex(2)
    self.BSlipSystemDetails.clicked.connect(slip_system_details)
    
    # Add slip system(s) to list
    def add_slip_system():
        # Get the currently selected items
        selected_items = self.TSlipSystems.selectedItems()

        # Check if there's a selected item
        if selected_items:
            selected_item = selected_items[0]  # Get the first selected item

            # Check if the selected item is a top-level item
            if selected_item.parent() is None:
                # Iterate through the children of the selected top-level item
                for i in range(selected_item.childCount()):
                    child = selected_item.child(i)
                    # Add the child's text to the list widget
                    self.INSlipSystems.addItem(child.text(0))
            else:
                # If it's not a top-level item, show a message or handle as needed
                self.INSlipSystems.addItem(selected_item.text(0))
            
            # Add slip systems as a selection to the Voce parameters combobox and table
            self.INSelectSlipSystem.clear()
            self.INSelectSlipSystem.addItem('ALL')
            self.TSlipSystemParameters.clearContents()
            self.TSlipSystemParameters.setRowCount(0)
            for i in range(self.INSlipSystems.count()):
                item = self.INSlipSystems.item(i)
                # Add the item's text to the QComboBox
                self.INSelectSlipSystem.addItem(item.text())
                # Add the item's text to the Table
                self.TSlipSystemParameters.insertRow(i)
                self.TSlipSystemParameters.setItem(i,0,QTableWidgetItem(item.text()))
            
    self.BAddSlipSystem.clicked.connect(add_slip_system)
    
    # Remove slip system from the list
    def remove_slip_system():
        for item in self.INSlipSystems.selectedItems():
            row = self.INSlipSystems.row(item)
            self.INSlipSystems.takeItem(row)
            del item
        # Remove slip systems as a selection to the Voce parameters combobox and table
        self.INSelectSlipSystem.clear()
        self.TSlipSystemParameters.clearContents()
        self.TSlipSystemParameters.setRowCount(0)
        if self.INSlipSystems.count() > 0:
            self.INSelectSlipSystem.addItem('ALL')
        for i in range(self.INSlipSystems.count()):
            item = self.INSlipSystems.item(i)
            # Add the item's text to the QComboBox
            self.INSelectSlipSystem.addItem(item.text())
            # Add the item's text to the Table
            self.TSlipSystemParameters.insertRow(i)
            self.TSlipSystemParameters.setItem(i,0,QTableWidgetItem(item.text()))
    self.BRemoveSlipSystem.clicked.connect(remove_slip_system)
    
    # Add Voce parameters to table based on assigned slip systems
    self.updating = False
    def voce_parameters():
        if self.updating:
            return
        self.updating = True
        slip_system = self.INSelectSlipSystem.currentText()
        rows = self.TSlipSystemParameters.rowCount()
        for i in range(rows):
            if "ALL" in slip_system:
                self.TSlipSystemParameters.setItem(i, 1, QTableWidgetItem(self.INnrsx.text()))
                self.TSlipSystemParameters.setItem(i, 2, QTableWidgetItem(self.INgamd0x.text()))
                self.TSlipSystemParameters.setItem(i, 3, QTableWidgetItem(self.INtau0xf.text()))
                self.TSlipSystemParameters.setItem(i, 4, QTableWidgetItem(self.INtau0xb.text()))
                self.TSlipSystemParameters.setItem(i, 5, QTableWidgetItem(self.INtau1x.text()))
                self.TSlipSystemParameters.setItem(i, 6, QTableWidgetItem(self.INthet0.text()))
                self.TSlipSystemParameters.setItem(i, 7, QTableWidgetItem(self.INthet1.text()))
                self.TSlipSystemParameters.setItem(i, 8, QTableWidgetItem(self.INhselfx.text()))
                self.TSlipSystemParameters.setItem(i, 9, QTableWidgetItem(self.INhlatex.text()))
            elif slip_system in self.TSlipSystemParameters.item(i, 0).text():
                self.TSlipSystemParameters.setItem(i, 1, QTableWidgetItem(self.INnrsx.text()))
                self.TSlipSystemParameters.setItem(i, 2, QTableWidgetItem(self.INgamd0x.text()))
                self.TSlipSystemParameters.setItem(i, 3, QTableWidgetItem(self.INtau0xf.text()))
                self.TSlipSystemParameters.setItem(i, 4, QTableWidgetItem(self.INtau0xb.text()))
                self.TSlipSystemParameters.setItem(i, 5, QTableWidgetItem(self.INtau1x.text()))
                self.TSlipSystemParameters.setItem(i, 6, QTableWidgetItem(self.INthet0.text()))
                self.TSlipSystemParameters.setItem(i, 7, QTableWidgetItem(self.INthet1.text()))
                self.TSlipSystemParameters.setItem(i, 8, QTableWidgetItem(self.INhselfx.text()))
                self.TSlipSystemParameters.setItem(i, 9, QTableWidgetItem(self.INhlatex.text()))
        self.updating = False
    self.INnrsx.textChanged.connect(voce_parameters)
    self.INgamd0x.textChanged.connect(voce_parameters)
    self.INtau0xf.textChanged.connect(voce_parameters)
    self.INtau0xb.textChanged.connect(voce_parameters)
    self.INtau1x.textChanged.connect(voce_parameters)
    self.INthet0.textChanged.connect(voce_parameters)
    self.INthet1.textChanged.connect(voce_parameters)
    self.INhselfx.textChanged.connect(voce_parameters)
    self.INhlatex.textChanged.connect(voce_parameters)
    
    # Add voce parameters to input items from table when qcombobox changes
    def voce_parameters_i():
        if self.updating:
            return
        self.updating = True
        slip_system = self.INSelectSlipSystem.currentText()
        rows = self.TSlipSystemParameters.rowCount()
        for i in range(rows):
            if "ALL" in slip_system:
                self.INnrsx.clear()
                self.INgamd0x.clear()
                self.INtau0xf.clear()
                self.INtau0xb.clear()
                self.INtau1x.clear()
                self.INthet0.clear()
                self.INthet1.clear()
                self.INhselfx.clear()
                self.INhlatex.clear()
            elif slip_system in self.TSlipSystemParameters.item(i, 0).text():
                self.INnrsx.setText(self.TSlipSystemParameters.item(i, 1).text() if self.TSlipSystemParameters.item(i, 1) is not None else "")
                self.INgamd0x.setText(self.TSlipSystemParameters.item(i, 2).text() if self.TSlipSystemParameters.item(i, 2) is not None else "")
                self.INtau0xf.setText(self.TSlipSystemParameters.item(i, 3).text() if self.TSlipSystemParameters.item(i, 3) is not None else "")
                self.INtau0xb.setText(self.TSlipSystemParameters.item(i, 4).text() if self.TSlipSystemParameters.item(i, 4) is not None else "")
                self.INtau1x.setText(self.TSlipSystemParameters.item(i, 5).text() if self.TSlipSystemParameters.item(i, 5) is not None else "")
                self.INthet0.setText(self.TSlipSystemParameters.item(i, 6).text() if self.TSlipSystemParameters.item(i, 6) is not None else "")
                self.INthet1.setText(self.TSlipSystemParameters.item(i, 7).text() if self.TSlipSystemParameters.item(i, 7) is not None else "")
                self.INhselfx.setText(self.TSlipSystemParameters.item(i, 8).text() if self.TSlipSystemParameters.item(i, 8) is not None else "")
                self.INhlatex.setText(self.TSlipSystemParameters.item(i, 9).text() if self.TSlipSystemParameters.item(i, 9) is not None else "")
        self.updating = False
    self.INSelectSlipSystem.currentIndexChanged.connect(voce_parameters_i)
    
    # Add material assignment to region
    def add_material_assignment_2():
        warning_flag = 0
        for i in range(self.TMaterialAssignment_2.rowCount()):
            if str(self.INRegion_3.currentText()) == self.TMaterialAssignment_2.item(i,0).text():
                warning_message('ERROR: There is already a material assigned to this region')
                warning_flag = 1
        if warning_flag == 0:
            row = self.TMaterialAssignment_2.rowCount()
            self.TMaterialAssignment_2.insertRow(row)
            self.TMaterialAssignment_2.setItem(row, 0, QTableWidgetItem(self.INRegion_3.currentText()))
            self.TMaterialAssignment_2.setItem(row, 1, QTableWidgetItem(self.INMaterial_3.currentText()))
        else:
            warning_flag = 0
    self.BAddMaterialAssignment_2.clicked.connect(add_material_assignment_2)
    
    # Delete regional material assignment
    def delete_material_assignemnt_2():
        current_row = self.TMaterialAssignment_2.currentRow()
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
            self.TMaterialAssignment_2.removeRow(current_row)
    self.BDeleteMaterialAssignment_2.clicked.connect(delete_material_assignemnt_2)
    
    # Warn user if no material assignment was made
    def warning_no_material_2():
        index = self.NavigationMenu.currentIndex()
        if index > 4 and self.TMaterialAssignment_2.rowCount() == 0 and self.TMaterials_2.rowCount() > 0:
            warning_message('WARNING: No materials were assigned')
    self.NavigationMenu.currentChanged.connect(warning_no_material_2)
    
    # Connect boundary conditions entered in initial condition table to solver starting estimate table
    def update_table(item):
        # get information from the first table
        row = item.row()
        col = item.column()
        value = item.text()
        
        # update the second and third tables
        table_item2 = QTableWidgetItem(value)
        table_item3 = QTableWidgetItem()
        if value.strip() == "":
            table_item2.setFlags(table_item2.flags() | Qt.ItemIsEditable)
            table_item2.setBackground(QBrush())
            table_item3.setFlags(table_item3.flags() | Qt.ItemIsEditable)
            table_item3.setBackground(QBrush())
        else:
            table_item2.setFlags(table_item2.flags() & ~Qt.ItemIsEditable)
            table_item2.setBackground(QBrush(QColor(214,214,214)))
            table_item3.setFlags(table_item3.flags() & ~Qt.ItemIsEditable)
            table_item3.setBackground(QBrush(QColor(214,214,214)))
        self.TVgradi.setItem(row, col, table_item2)
        self.TCstress.setItem(row, col, table_item3)
    self.TVgrad.itemChanged.connect(update_table)
    
    # Populate pre-defined boundary conditions
    def boundary_conditions():
        if "Custom" in self.INbulkBC.currentText():
            # Clear tables
            self.TVgrad.clearContents()
            self.TVgradi.clearContents()
            self.TCstress.clearContents()
        elif "Test" in self.INbulkBC.currentText():
            # Clear tables
            self.TVgrad.clearContents()
            self.TVgradi.clearContents()
            self.TCstress.clearContents()
            
            # Assign values
            self.TVgrad.setItem(0,1,QTableWidgetItem("0."))
            self.TVgrad.setItem(0,2,QTableWidgetItem("0."))
            self.TVgrad.setItem(1,0,QTableWidgetItem("0."))
            self.TVgrad.setItem(1,2,QTableWidgetItem("0."))
            self.TVgrad.setItem(2,0,QTableWidgetItem("0."))
            self.TVgrad.setItem(2,1,QTableWidgetItem("0."))
            self.TVgrad.setItem(2,2,QTableWidgetItem("1.0"))
            self.TVgradi.setItem(0,0,QTableWidgetItem("-0.35"))
            self.TVgradi.setItem(1,1,QTableWidgetItem("-0.35"))
            self.TCstress.setItem(0,0,QTableWidgetItem("0."))
            self.TCstress.setItem(1,1,QTableWidgetItem("0."))
    self.INbulkBC.currentIndexChanged.connect(boundary_conditions)
    
    # Run Bulk Formation
    self.run = 0
    def run_bulk_forming():
        self.run = 1
        # Create location to save files
        self.job_name = "simulation_files"
        self.working_directory = os.path.join(self.bulk_forming_dir, f'{self.job_name}')
        print(self.working_directory)
        if not os.path.exists(self.working_directory):
            os.makedirs(self.working_directory)
            
        # Create input files
        self.BULK_FORMING_ELASTIC_PARAMETERS = os.path.join(self.working_directory, 'elastic_parameters.txt')
        self.BULK_FORMING_PLASTIC_PARAMETERS = os.path.join(self.working_directory, 'plastic_parameters.txt')
        self.BULK_FORMING_INPUT = os.path.join(self.working_directory, 'bulk_forming_input.txt')
        Bulk_Forming_WInput(self)
        
        # Define executable path
        in_file_path = self.TParts.item(0,10).text()
        file_type = os.path.splitext(in_file_path)[1].lower()
        reload(DeveloperInputs)
        if self.UserConfig == "Developer":
            executable_path = DeveloperInputs.fierro_evpfft_exe
        elif self.UserConfig == "User":
            executable_path = "evpfft"
        arguments = ["-f", self.BULK_FORMING_INPUT]
        
        self.p = QProcess()
        self.p.setWorkingDirectory(self.working_directory)
        self.p.readyReadStandardOutput.connect(handle_stdout)
        self.p.readyReadStandardError.connect(handle_stderr)
        self.p.stateChanged.connect(handle_state)
        self.p.finished.connect(lambda: process_finished(1))
        try:
            self.p.start(executable_path, arguments)
        except Exception as e:
            self.warning_message("ERROR: evpfft executable")
            return
        self.progress_re = re.compile("       Current  Time  STEP = (\\d+)")
        
        # Save job directory
        self.INBFJobDir.setText(f'{self.working_directory}')
    self.BRunBulkForming.clicked.connect(run_bulk_forming)
            
    def simple_percent_parser(output):
        m = self.progress_re.search(output)
        if m:
            pc_complete = m.group(1)
            return int(pc_complete)
    def process_finished(num):
        self.RunOutputProgress.setValue(100)
        self.p.close()
        self.p = None
    def handle_stdout():
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        progress = simple_percent_parser(stdout)
        if progress:
            self.RunOutputProgress.setValue((progress/int(self.INBFloadsteps.text()))*100)
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
        
    # Preview Results
    def preview_results_click_2():
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
        file_name = "micro_state_timestep_" + str(self.INBFloadsteps.text()) + ".xdmf"
        self.output_directory = os.path.join(self.working_directory, file_name)
        self.results_reader = paraview.simple.XDMFReader(FileNames=self.output_directory)
        paraview.simple.SetDisplayProperties(Representation="Surface")
        self.display = Show(self.results_reader, self.render_view)
        
        # Color by the selected variable
        selected_variable = str(self.INBFResults.currentText())
        paraview.simple.ColorBy(self.display, ('CELLS', selected_variable))
        vmstressLUT = paraview.simple.GetColorTransferFunction(selected_variable)
        self.display.SetScalarBarVisibility(self.render_view, True)

        # View Results
        self.render_view.ResetCamera()
        self.render_view.StillRender()
    self.INBFResults.currentIndexChanged.connect(preview_results_click_2)
    
    # Show results immediately when postprocessing tab is pressed
    def show_results_2():
        index = self.NavigationMenu.currentIndex()
        if index == 8 and self.run == 1 and "Bulk Forming" in self.INSelectPostprocessing.currentText():
            preview_results_click_2()
    self.NavigationMenu.currentChanged.connect(show_results_2)
    
    # Open Paraview
    def open_paraview_click_2():
        command = ["paraview", self.output_directory]
        subprocess.Popen(command)
    self.BBFParaview.clicked.connect(open_paraview_click_2)
        
    