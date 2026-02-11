from PySide6.QtWidgets import (QTableWidgetItem, QMessageBox, QApplication, QFileDialog, QWidget, QHBoxLayout, QVBoxLayout, QLabel, QLineEdit, QTableWidget, QPushButton, QTreeWidgetItem, QAbstractItemView)
from PySide6.QtCore import (QCoreApplication, QProcess, Qt, QProcessEnvironment)
from PySide6.QtGui import (QColor, QBrush, QFont)
import numpy as np
import re
import os
import subprocess
import vtk
from paraview.simple import *
from importlib import reload
from Bulk_Forming_WInput import *
import DeveloperInputs
from Reload_Geometry import *
from LAFFT2VTK import *

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
    
    # Control addition of plasticity
    def plasticity():
        if self.BEnablePlasticity.isChecked():
            self.EnablePlasticity.setCurrentIndex(1)
        else:
            self.EnablePlasticity.setCurrentIndex(0)
    self.BEnablePlasticity.stateChanged.connect(plasticity)
    
    # Define material properties using predefined values
    self.clear_flag = 0
    def predefined_materials():
        if self.clear_flag == 1:
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
            self.clear_flag = 0
        if "Import Elastic Parameters File" in self.INMaterialDefinition.currentText():
            self.MaterialMenu_2.setCurrentIndex(0)
            try:
                self.elastic_filename
            except:
                self.elastic_filename, _ = QFileDialog.getOpenFileName(filter="Elastic Parameters File (*.txt)",)
            matrix = []
            with open(self.elastic_filename, 'r') as file:
                # Read the first line
                first_line = file.readline().strip()
                
                # Check if the first character of the first line is '0' - anisotropic
                if first_line[0] == '0':
                    # Set definition to anisotropic
                    self.INMaterialType_2.setCurrentIndex(3)
                    # Read next six lines for the matrix
                    for _ in range(6):
                        line = file.readline().strip()
                        if line:  # Ensure line is not empty
                            row = line.split()[:6]  # Convert to float
                            row = [float(val) for val in row if val.replace('.', '').isdigit()]
                            matrix.append(row)
                    # Add to anisotropic table
                    for i in range(6):
                        for j in range(i,6):
                            item = QTableWidgetItem(str(matrix[i][j]))
                            self.TAnisotropic_2.setItem(i, j, item)
                    # Turn page to anisotropic
                    self.MaterialTypeTool_2.setCurrentIndex(2)
                # Check if the first character of the first line is '1' - isotropic
                elif first_line[0] == '1':
                    # Set definition to isotropic
                    self.INMaterialType_2.setCurrentIndex(0)
                    line = file.readline().strip()
                    row = line.split()[:2]
                    row = [float(val) for val in row if val.replace('.', '').isdigit()]
                    matrix.append(row)
                    # Add to isotropic line edit definitions
                    self.INYoungsModulus_2.setText(str(matrix[0][0]))
                    self.INPoissonsRatio_2.setText(str(matrix[0][1]))
                    # Turn page to isotropic
                    self.MaterialTypeTool_2.setCurrentIndex(0)
                else:
                    warning_message("ERROR: The first character of the first line is not a '0' (anisotropic) or '1' (isotropic)")
                    return None
            del self.elastic_filename
        elif "Import Plastic Parameters File" in self.INMaterialDefinition.currentText():
            self.MaterialMenu_2.setCurrentIndex(1)
            self.BEnablePlasticity.setChecked(True)
            try:
                self.plastic_filename
            except:
                self.plastic_filename, _ = QFileDialog.getOpenFileName(filter="Plastic Parameters File (*.txt)",)
            nsmx = 0
            with open(self.plastic_filename, 'r') as file:
                iline = 1;
                for line in file:
                    # Find crystal axes
                    if iline == 3:
                        crystal_axes_line = line.strip().split()[:3]
                        if len(crystal_axes_line) < 3:
                            warning_message("ERROR: crystal axes was not found on line 3.")
                            return
                        else:
                            self.INa.setText(str(crystal_axes_line[0]))
                            self.INb.setText(str(crystal_axes_line[1]))
                            self.INc.setText(str(crystal_axes_line[2]))
                    # Find slip systems
                    if iline == 2:
                        slip_type = line.strip().split()[0]
                    if iline == 4:
                        nmodesx = line.strip().split()[:1]
                        nmodesx = int(nmodesx[0])
                    if iline == 5:
                        nmodes = line.strip().split()[:1]
                        nmodes = int(nmodes[0])
                    if iline == 6:
                        modei = line.strip().split()[:nmodes]
                        modei = [float(val) for val in modei if val.replace('.', '').isdigit()]
                    if iline > 6:
                        if iline == 7:
                            start = 7
                        if iline == start:
                            self.BCustomSlipSystem.click()
                            name = line.strip().split()
                            self.line_edit.setText(' '.join(name))
                        if iline == start+1:
                            param1 = line.strip().split()[:6]
                            modex = float(param1[0])
                            nsmx = float(param1[1])
                            self.table.setRowCount(nsmx)
                        if iline == start+2:
                            param2 = line.strip().split()[:5]
                        if iline == start+3:
                            param3 = line.strip().split()[:2]
                        if iline > start+3:
                            if "CUB" in slip_type or "cub" in slip_type or "ORT" in slip_type or "ort" in slip_type:
                                notation = 3
                            elif "HEX" in slip_type or "hex" in slip_type or "TRI" in slip_type or "tri" in slip_type:
                                notation = 4
                            else:
                                warning_message("ERROR: you must specify icryst (CUBIC, HEX, etc.)")
                                return
                            plane = line.strip().split()[:notation]
                            direction = line.strip().split()[notation:2*notation]
                            self.table.setItem(iline-start-4,0,QTableWidgetItem(','.join(str(x) for x in plane)))
                            self.table.setItem(iline-start-4,1,QTableWidgetItem(','.join(str(x) for x in direction)))
                            if iline == start+3+nsmx:
                                self.BSubmit.click()
                                nmodesx = nmodesx-1
                                custom_label = self.TSlipSystems.findItems("CUSTOM", Qt.MatchExactly, 0)
                                custom_item = next((item for item in custom_label if item.parent() is None), None)
                                self.TSlipSystems.setCurrentItem(custom_item.child(modex-1))
                                if modex in modei:
                                    self.BAddSlipSystem.click()
                                    row = self.TSlipSystemParameters.rowCount()
                                    self.TSlipSystemParameters.setItem(row-1,1,QTableWidgetItem(param1[2]))
                                    self.TSlipSystemParameters.setItem(row-1,2,QTableWidgetItem(param1[3]))
                                    self.TSlipSystemParameters.setItem(row-1,3,QTableWidgetItem(param2[0]))
                                    self.TSlipSystemParameters.setItem(row-1,4,QTableWidgetItem(param2[1]))
                                    self.TSlipSystemParameters.setItem(row-1,5,QTableWidgetItem(param2[2]))
                                    self.TSlipSystemParameters.setItem(row-1,6,QTableWidgetItem(param2[3]))
                                    self.TSlipSystemParameters.setItem(row-1,7,QTableWidgetItem(param2[4]))
                                    self.TSlipSystemParameters.setItem(row-1,8,QTableWidgetItem(param3[0]))
                                    self.TSlipSystemParameters.setItem(row-1,9,QTableWidgetItem(param3[1]))
                                if nmodesx == 0:
                                    return
                                else:
                                    start = iline+1
                    # Update line number
                    iline += 1
            # Turn page to plastic
            self.BEnablePlasticity.setChecked(True)
            self.MaterialMenu_2.setCurrentIndex(1)
        elif "Single Crystal FCC" in self.INMaterialDefinition.currentText():
            self.clear_flag = 1
            if 'MPa' in self.INUnits.currentText():
                m = 1
            elif 'Pa' in self.INUnits.currentText():
                m = 1000000
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
            self.BEnablePlasticity.setChecked(True)
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
        elif "Single Crystal BCC" in self.INMaterialDefinition.currentText():
            self.clear_flag = 1
            # Define elastic properties
            for i in [0,1,2,3,4,5,6]:
                for j in range(i,6):
                    self.TAnisotropic_2.setItem(i,j,QTableWidgetItem('0.'))
            self.INSolidGas_2.setCurrentIndex(0)
            self.INMaterialType_2.setCurrentIndex(3)
            self.TAnisotropic_2.setItem(0,0,QTableWidgetItem('266700.'))
            self.TAnisotropic_2.setItem(1,1,QTableWidgetItem('266700.'))
            self.TAnisotropic_2.setItem(2,2,QTableWidgetItem('266700.'))
            self.TAnisotropic_2.setItem(0,1,QTableWidgetItem('160800.'))
            self.TAnisotropic_2.setItem(0,2,QTableWidgetItem('160800.'))
            self.TAnisotropic_2.setItem(1,2,QTableWidgetItem('160800.'))
            self.TAnisotropic_2.setItem(3,3,QTableWidgetItem('82500.'))
            self.TAnisotropic_2.setItem(4,4,QTableWidgetItem('82500.'))
            self.TAnisotropic_2.setItem(5,5,QTableWidgetItem('82500.'))
            # Define plastic properties
            self.BEnablePlasticity.setChecked(True)
            self.INa.setText('1.')
            self.INb.setText('1.')
            self.INc.setText('1.')
            items = self.TSlipSystems.findItems('BCC', Qt.MatchExactly | Qt.MatchRecursive)
            if items:
                self.TSlipSystems.setCurrentItem(items[0])
                self.BAddSlipSystem.click()
                self.TSlipSystems.expandItem(items[0])
                self.INSlipSystems.setCurrentRow(2)
                self.BRemoveSlipSystem.click()
                self.INnrsx.setText('14')
                self.INgamd0x.setText('1.0')
                self.INtau0xf.setText('115.')
                self.INtau0xb.setText('115.')
                self.INtau1x.setText('30.')
                self.INthet0.setText('110.')
                self.INthet1.setText('5.')
                self.INhselfx.setText('1.0')
                self.INhlatex.setText('1.2')
        elif "Custom" in self.INMaterialDefinition.currentText():
            self.clear_flag = 1
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
        if self.BEnablePlasticity.isChecked():
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
        if self.BEnablePlasticity.isChecked():
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

    # View details about a slip system
    def slip_system_details():
        selected_items = self.TSlipSystems.selectedItems()
        if not selected_items:
            warning_message("ERROR: Please select a specific slip system to view it's details.")
        else:
            selected_item = selected_items[0]
            if selected_item.parent() is not None:
                page = selected_item.text(0).split('.', 1)[0]
                self.SlipSystemInfo.setCurrentIndex(int(page))
                self.PlasticProperties.setCurrentIndex(2)
            else:
                warning_message("ERROR: Please select a specific slip system to view it's details.")
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
            for i in range(self.INSlipSystems.count()):
                item = self.INSlipSystems.item(i)
                # Add the item's text to the QComboBox
                self.INSelectSlipSystem.addItem(item.text())
            row = self.TSlipSystemParameters.rowCount()
            self.TSlipSystemParameters.insertRow(row)
            self.TSlipSystemParameters.setItem(row,0,QTableWidgetItem(item.text()))
            
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
    
    # Define custom slip system
    def custom_slip_system():
        # Index of new page
        index = self.SlipSystemInfo.count()
        # Layout
        new_page = QWidget()
        page_layout = QVBoxLayout()
        page_layout.setContentsMargins(0, 0, 0, 0)
        page_layout.setSpacing(5)
        # Functionality
        # Slip System Name
        input_layout = QHBoxLayout()
        label = QLabel("Slip System Name:")
        input_layout.addWidget(label)
        self.line_edit = QLineEdit()
        input_layout.addWidget(self.line_edit)
        page_layout.addLayout(input_layout)
        # Add or remove colums from table
        arbuttons = QHBoxLayout()
        self.BAddRow = QPushButton("Add Row")
        arbuttons.addWidget(self.BAddRow)
        self.BRemoverow = QPushButton("Remove Row")
        arbuttons.addWidget(self.BRemoverow)
        page_layout.addLayout(arbuttons)
        # Functions to add or remove rows from table
        def add_row():
            current_row_count = self.table.rowCount()
            self.table.insertRow(current_row_count)
        self.BAddRow.clicked.connect(add_row)
        def remove_row():
            current_row = self.table.currentRow()  # Get the index of the selected row
            if current_row >= 0:  # Check if a row is selected
                self.table.removeRow(current_row)  # Remove the selected row
        self.BRemoverow.clicked.connect(remove_row)
        # Table
        setattr(self, f"T{index}", QTableWidget())
        self.table = getattr(self, f"T{index}")
        self.table.setRowCount(12)
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["(Slip Plane)", "<Slip Direction>"])
        self.table.horizontalHeader().setDefaultSectionSize(167)
        font = QFont()
        font.setPointSize(13)
        self.table.horizontalHeader().setFont(font)
        page_layout.addWidget(self.table)
        # Add custom slip system
        self.BSubmit = QPushButton("Add Custom Slip System")
        page_layout.addWidget(self.BSubmit)
        # Function to add the slip system to table
        def add_custom_system():
            # Check that everything is completely filled out
            warning_flag = 0
            if not self.line_edit.text().strip():
                warning_message("ERROR: Please assign a name to custom slip system.")
                return
            for row in range(self.table.rowCount()):
                for col in range(self.table.columnCount()):
                    item = self.table.item(row,col)
                    if item is None or not item.text().strip():
                        warning_flag = 1
            if warning_flag == 1:
                warning_message("ERROR: Ensure custom slip system is completely filled out, and delete any unused rows.")
                return
            # If there are no errors, add custom slip system to the table
            else:
                # count the number of slip systems already defined in the table
                num_children = 0
                for i in range(self.TSlipSystems.topLevelItemCount()):
                    top_level_item = self.TSlipSystems.topLevelItem(i)
                    num_children += top_level_item.childCount()
                # Find the custom heading and add the new system to it
                custom_label = self.TSlipSystems.findItems("CUSTOM", Qt.MatchExactly, 0)
                custom_item = next((item for item in custom_label if item.parent() is None), None)
                if custom_item:
                    # Add a subitem to the found "custom" item
                    custom_name = str(num_children+1) + '. ' + self.line_edit.text().strip()
                    new_subitem = QTreeWidgetItem(custom_item, [custom_name])
                    custom_item.addChild(new_subitem)
                    # Disable things so people can't change them after it has been added
                    self.BSubmit.setEnabled(False)
                    self.line_edit.setEnabled(False)
                    self.BAddRow.setEnabled(False)
                    self.BRemoverow.setEnabled(False)
                    self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
                    # Go back to main slip system screen
                    self.PlasticProperties.setCurrentIndex(1)
                    # Expand custom definition
                    self.TSlipSystems.expandItem(custom_label[0])
        self.BSubmit.clicked.connect(add_custom_system)
        # Setup New Page
        new_page.setLayout(page_layout)
        self.SlipSystemInfo.addWidget(new_page)
        # Display
        self.PlasticProperties.setCurrentIndex(2)
        self.SlipSystemInfo.setCurrentIndex(index)
    self.BCustomSlipSystem.clicked.connect(custom_slip_system)
    
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
#            if "ALL" in slip_system:
#                self.INnrsx.clear()
#                self.INgamd0x.clear()
#                self.INtau0xf.clear()
#                self.INtau0xb.clear()
#                self.INtau1x.clear()
#                self.INthet0.clear()
#                self.INthet1.clear()
#                self.INhselfx.clear()
#                self.INhlatex.clear()
            if slip_system in self.TSlipSystemParameters.item(i, 0).text():
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
        elif "Tension Z" in self.INbulkBC.currentText():
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
            self.TCstress.setItem(0,0,QTableWidgetItem("0."))
            self.TCstress.setItem(1,1,QTableWidgetItem("0."))
        elif "Shear XY" in self.INbulkBC.currentText():
            # Clear tables
            self.TVgrad.clearContents()
            self.TVgradi.clearContents()
            self.TCstress.clearContents()
            
            # Assign values
            self.TVgrad.setItem(0,0,QTableWidgetItem("0."))
            self.TVgrad.setItem(0,1,QTableWidgetItem("-1.0"))
            self.TVgrad.setItem(0,2,QTableWidgetItem("0."))
            self.TVgrad.setItem(1,0,QTableWidgetItem("0."))
            self.TVgrad.setItem(1,1,QTableWidgetItem("0."))
            self.TVgrad.setItem(1,2,QTableWidgetItem("0."))
            self.TVgrad.setItem(2,0,QTableWidgetItem("0."))
            self.TVgrad.setItem(2,1,QTableWidgetItem("0."))
            self.TVgrad.setItem(2,2,QTableWidgetItem("0."))
        elif "Compression Z" in self.INbulkBC.currentText():
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
            self.TVgrad.setItem(2,2,QTableWidgetItem("-1.0"))
            self.TCstress.setItem(0,0,QTableWidgetItem("0."))
            self.TCstress.setItem(1,1,QTableWidgetItem("0."))
    self.INbulkBC.currentIndexChanged.connect(boundary_conditions)
    
    # Select to run locally or write the input files
    self.INBFRunLocally.toggled.connect(lambda: self.BFRunWrite.setCurrentIndex(0))
    self.INBFWriteFiles.toggled.connect(lambda: self.BFRunWrite.setCurrentIndex(1))
    
    # Select serial or parallel run type
    self.INBFSerial.toggled.connect(lambda: self.BFRunType.setCurrentIndex(0))
    self.INBFParallel.toggled.connect(lambda: self.BFRunType.setCurrentIndex(1))
    
    # Set number of mpi ranks based on the number of cores that are avaliable on the system
    num_cores = os.cpu_count()
    self.INBFmpiRanks.setMaximum(num_cores)
    
    # Write input files
    def Write_Input_Files():
        # Ask user to select a folder
        selected_directory = QFileDialog.getExistingDirectory(None, "Select Folder")
        # Create input files
        self.BULK_FORMING_ELASTIC_PARAMETERS = os.path.join(selected_directory, 'elastic_parameters.txt')
        self.BULK_FORMING_PLASTIC_PARAMETERS = os.path.join(selected_directory, 'plastic_parameters.txt')
        self.BULK_FORMING_INPUT = os.path.join(selected_directory, 'bulk_forming_input.txt')
        Bulk_Forming_WInput(self)
        # Write a readme file
        readme = os.path.join(selected_directory, 'README.txt')
        wreadme = open(readme,"w")
        about = 'ABOUT: \n' \
                 'elastic_parameters.txt -> these files contain the elastic material properties\n' \
                 'plastic_parameters.txt -> this file contains the plastic material properties\n' \
                 'bulk_forming_input.txt -> this file is the input file to the evpfft solver\n' \
                 '            -> ** Ensure you go into this file and properly specify file paths\n'
        wreadme.write(about)
        run = 'TO RUN: \n' \
               'flags -> Please set the flags: export OMP_PROC_BIND=spread and export OMP_NUM_THREADS=1\n' \
               'parallel run -> (.txt input): mpirun -np # evpfft -f bulk_forming_input.txt\n' \
               '             -> (.vtk input): mpirun -np # evpfft -f bulk_forming_input.txt -m 2\n' \
               'serial run -> (.txt input): evpfft -f bulk_forming_input.txt\n' \
               '           -> (.vtk input): evpfft -f bulk_forming_input.txt -m 2\n' \
               'help -> evpfft --help'
        wreadme.write(run)
        wreadme.close()
    self.BBFWriteFiles.clicked.connect(Write_Input_Files)
    
    # Run Bulk Formation
    self.run = 0
    self.ConvergenceError = False
    self.min_iterations_old = 0
    def run_bulk_forming():
        self.run = 1
        # Restart progress bar
        self.RunOutputProgress.setValue(0)
        
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
            if self.INBFParallel.isChecked():
                executable_path = "mpirun"
            else:
                executable_path = DeveloperInputs.fierro_evpfft_exe
        elif self.UserConfig == "User":
            if self.INBFParallel.isChecked():
                executable_path = "mpirun"
            else:
                if self.INLargeStrain.isChecked():
                    executable_path = "ls-evpfft"
                else:
                    executable_path = "evpfft"
                
        # Define arguments
        if ".vtk" in file_type:
            if self.INBFParallel.isChecked():
                if self.INLargeStrain.isChecked():
                    arguments = ["-np", self.INBFmpiRanks.value(), "ls-evpfft", "-f", self.BULK_FORMING_INPUT, "-m", "2"]
                else:
                    arguments = ["-np", self.INBFmpiRanks.value(), "evpfft", "-f", self.BULK_FORMING_INPUT, "-m", "2"]
            else:
                arguments = ["-f", self.BULK_FORMING_INPUT, "-m", "2"]
        else:
            if self.INBFParallel.isChecked():
                if self.INLargeStrain.isChecked():
                    arguments = ["-np", self.INBFmpiRanks.value(), "ls-evpfft", "-f", self.BULK_FORMING_INPUT]
                else:
                    arguments = ["-np", self.INBFmpiRanks.value(), "evpfft", "-f", self.BULK_FORMING_INPUT]
            else:
                arguments = ["-f", self.BULK_FORMING_INPUT]
        
        # Set up the environment
        self.p = QProcess()
        env = QProcessEnvironment.systemEnvironment()
        if not env.contains("OMP_PROC_BIND"):
            env.insert("OMP_PROC_BIND", "spread")
        if not env.contains("OMP_NUM_THREADS"):
            env.insert("OMP_NUM_THREADS", "1")
        self.p.setProcessEnvironment(env)
        # Set up the states
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
        # Count how many iterations were taken towards the solution
        self.iterations_re = re.compile(r" ITER = (\d+)")
        
        # Save job directory
        self.INBFJobDir.setText(f'{self.working_directory}')
        
        # Wait for program to finish before checking for errors
        self.p.waitForStarted()
        while self.p != None:
            QApplication.processEvents()
        
        # Provide warnings for convergence
        if self.ConvergenceError == True:
            self.ConvergenceError = False
            self.min_iterations_old = 0
            warning_message("WARNING: It is recomended that you increase your maximum number of iterations. Solution convergence was NOT achieved.")
    self.BRunBulkForming.clicked.connect(run_bulk_forming)
            
    def convergence_check(output):
        iterations = self.iterations_re.findall(output)
        if iterations:
            return int(iterations[-1])
    def simple_percent_parser(output):
        m = self.progress_re.search(output)
        if m:
            pc_complete = m.group(1)
            return int(pc_complete)
    def process_finished(num):
        handle_stdout()
        handle_stderr()
        self.RunOutputProgress.setValue(100)
        # If the iterations reached the maximum number, write out an error
        if self.min_iterations_old >= int(self.INBFmaxiter.text()):
            self.ConvergenceError = True
        self.min_iterations_old = 0
        self.p.close()
        self.p = None
        # Restart deformation scale factor
        self.INBFDeform.setValue(0)
    def handle_stdout():
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        progress = simple_percent_parser(stdout)
        if progress:
            self.RunOutputProgress.setValue((progress/int(self.INBFloadsteps.text()))*100)
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
    def kill_bulk_forming():
        self.p.terminate()
        self.RunOutputWindow.appendPlainText("TERMINATED")
    self.BKillBulkForming.clicked.connect(kill_bulk_forming)
        
    # Preview Results
    def preview_results_click_2():
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
#        file_name = "micro_state_timestep_" + str(self.INBFloadsteps.text()) + ".xdmf"
        file_name = f"MicroState_{int(self.INBFloadsteps.text()):05}.pvtu"
#        self.output_directory = os.path.join(self.working_directory, file_name)
        try:
            self.working_directory
        except:
            self.working_directory = self.INBFJobDir.text()
        else:
            False
        self.output_directory = os.path.join(self.working_directory, "pvtu", file_name)
        self.results_reader = paraview.simple.XMLPartitionedUnstructuredGridReader(FileName=self.output_directory)
#        paraview.simple.SetDisplayProperties(Representation="Surface")
#        self.display = Show(self.results_reader, self.render_view)
        
        # Apply warp filter as long as result wasn't run using more than 1 mpi rank
        if self.INBFSerial.isChecked():
            # Enable deformation scale factor
            self.INBFDeform.setEnabled(True)
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
                        if self.INLargeStrain.isChecked():
                            diffX.append(coords[count][0]-(float(i)))
                            diffY.append(coords[count][1]-(float(j)))
                            diffZ.append(coords[count][2]-(float(k)))
                        else:
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
            if self.INBFDeform.value() == 0:
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
                self.INBFDeform.setValue(initial_scale_factor)
            # Scale the displacements
            scale_factor = self.INBFDeform.value()  # Adjust this value to change the scaling
            scale_filter = paraview.simple.Calculator(Input=temp_source)
            scale_filter.ResultArrayName = 'ScaledDisplacement'
            scale_filter.Function = f'Displacement * {scale_factor}'
            paraview.simple.UpdatePipeline()
            # Warp Filter
            self.threshold = paraview.simple.WarpByVector(Input=scale_filter)
            self.threshold.Vectors = ['POINTS', 'ScaledDisplacement']
            # Display
            self.display = Show(self.threshold, self.render_view)
        else:
            # Disable deformation scale factor
            self.INBFDeform.setEnabled(False)
            # Display
            self.display = Show(self.results_reader, self.render_view)
        
        # Color by the selected variable
        selected_variable = str(self.INBFResults.currentText())
        paraview.simple.ColorBy(self.display, ('CELLS', selected_variable))
        vmstressLUT = paraview.simple.GetColorTransferFunction(selected_variable)
        self.display.SetScalarBarVisibility(self.render_view, True)

        # View Results
        self.render_view.ResetCamera()
        self.render_view.StillRender()
    
    # Show results immediately when postprocessing tab is pressed
    def show_results_2():
        index = self.NavigationMenu.currentIndex()
        if index == 8:
            self.INBFResults.currentIndexChanged.connect(preview_results_click_2)
            self.INBFDeform.valueChanged.connect(preview_results_click_2)
        if index == 8 and self.run == 1 and "Bulk Forming" in self.INSelectPostprocessing.currentText():
            preview_results_click_2()
    self.NavigationMenu.currentChanged.connect(show_results_2)
    
    # Go from results back to back to original input
    self.vis = 0
    def show_original_2():
        index = self.NavigationMenu.currentIndex()
        if index < 8 and self.vis == 1 and "Bulk Forming" in self.INSelectPostprocessing.currentText():
            # reset view
            self.INBFDeform.setValue(.1)
            paraview.simple.ColorBy(self.display, ('CELLS', 'grain_id'))
            # View Results
            self.display.SetScalarBarVisibility(self.render_view, False)
            self.render_view.ResetCamera()
            self.render_view.StillRender()
            self.vis = 0
    self.NavigationMenu.currentChanged.connect(show_original_2)
    
    # Open Paraview
    def open_paraview_click_2():
        command = ["paraview", self.output_directory]
        subprocess.Popen(command)
    self.BBFParaview.clicked.connect(open_paraview_click_2)
    
    # Upload Legacy EVPFFT input file
    def legacyEVPFFT():
        # Get all necessary information from the input file and do stuff with it right away. The reason it has to be done right away is due to the fact that you might have multiple phases.
        input_filename, _ = QFileDialog.getOpenFileName(filter="Input File (*.txt)",)
        bc_flag = 0
        other_flag = 0
        settings_flag = 0
        bc_start = 1000000000
        other_start = 1000000000
        settings_start = 1000000000
        iudot = []
        udot = []
        iscau = []
        scauchy = []
        with open(input_filename, 'r') as file:
            iline = 1;
            phase_names = 1;
            for line in file:
                if iline == 1:
                    # Get the number of phases
                    phases = int(line.strip().split()[0])
                    phase_start = 10
                if iline == 2:
                    # Get the voxel numbers
                    dims = line.strip().split()[:3]
                    dimx = int(dims[0])
                    dimy = int(dims[1])
                    dimz = int(dims[2])
                if iline == 4:
                    # Get the RVE dimensions
                    delt = line.strip().split()[:3]
                    deltx = float(delt[0])
                    delty = float(delt[1])
                    deltz = float(delt[2])
                if iline == 6:
                    if '/' in line or '\\' in line:
                        # Get name of microstructure file
                        microstructure_filename = line.strip()
                        # Fill out geometry table
                        row = self.TParts.rowCount()
                        self.TParts.insertRow(row)
                        self.TParts.setItem(row,0,QTableWidgetItem(f"Phase{phase_names}"))
                        self.TParts.setItem(row,1,QTableWidgetItem("0"))
                        self.TParts.setItem(row,2,QTableWidgetItem("0"))
                        self.TParts.setItem(row,3,QTableWidgetItem("0"))
                        lengthx = dimx*deltx
                        lengthy = dimy*delty
                        lengthz = dimz*deltz
                        self.TParts.setItem(row,4,QTableWidgetItem(f"{lengthx}"))
                        self.TParts.setItem(row,5,QTableWidgetItem(f"{lengthy}"))
                        self.TParts.setItem(row,6,QTableWidgetItem(f"{lengthz}"))                        
                        self.TParts.setItem(row,7,QTableWidgetItem(f"{dimx}"))
                        self.TParts.setItem(row,8,QTableWidgetItem(f"{dimy}"))
                        self.TParts.setItem(row,9,QTableWidgetItem(f"{dimz}"))
                        self.TParts.setItem(row,10,QTableWidgetItem(f"{microstructure_filename}"))
                        # convert text file to vtk file for visualization
                        vtk_location = self.voxelizer_dir + '/VTK_Geometry_' + self.TParts.item(row,0).text() + '.vtk'
                        los_alamos_to_vtk(microstructure_filename, vtk_location)
                        self.in_file_path = self.TParts.item(row,10).text()
                        self.file_type = os.path.splitext(self.in_file_path)[1].lower()
                        Reload_Geometry(self)
                        # Add the geometry as an option for material assignment
                        self.INRegion_3.addItem(self.TParts.item(row,0).text())
                    else:
                        warning_message("ERROR: couldn't find the microstructure file path on line 6.")
                        return
                if iline == phase_start:
                    if '/' in line or '\\' in line: # or 'dummy' in line:
                        # Get plastic parameters file path
                        self.plastic_filename = line.strip()
                        # Upload the plastic parameters
                        self.INMaterialDefinition.setCurrentIndex(2)
                    else:
                        warning_message(f"ERROR: couldn't find the plastic file path on line {phase_start}.")
                        return
                if iline == phase_start+1:
                    if '/' in line or '\\' in line: # or 'dummy' in line:
                        # Get elastic parameters file path
                        self.elastic_filename = line.strip()
                        # Upload the elastic parameters
                        self.INMaterialDefinition.setCurrentIndex(1)
                        # Add material as an option for assignment
                        self.INMaterialName_2.setText(f"Phase{phase_names}Material")
                        self.BAddMaterial_2.click()
                        # Assign the material
                        self.INRegion_3.setCurrentIndex(phase_names-1)
                        self.INMaterial_3.setCurrentIndex(phase_names-1)
                        self.BAddMaterialAssignment_2.click()
                        # This update is for if there is more than one phase defined in the file
                        phase_names += 1
                        phases -= 1
                        if phases > 0:
                            phase_start = phase_start + 5
                    else:
                        warning_message(f"ERROR: couldn't find the elastic file path on line {phase_start+1}.")
                        return
                if "boundary conditions" in line:
                    bc_flag = 1
                    bc_start = iline
                # Velocity gradient BC information
                if iline > bc_start and iline < bc_start + 4:
                    iudot_ = line.strip().split()[:3]
                    iudot_ = [float(iudot_[0]), float(iudot_[1]), float(iudot_[2])]
                    iudot.append(iudot_)
                if iline > bc_start + 4 and iline < bc_start + 8:
                    udot_ = line.strip().split()[:3]
                    udot_ = [float(udot_[0]), float(udot_[1]), float(udot_[2])]
                    udot.append(udot_)
                # Cauchy stress BC information
                if iline == bc_start + 9:
                    iscau_ = line.strip().split()[:3]
                    iscau_ = [float(iscau_[0]), float(iscau_[1]), float(iscau_[2])]
                    iscau.append(iscau_)
                if iline == bc_start + 10:
                    iscau_ = line.strip().split()[:2]
                    iscau_ = [float(iscau_[0]), float(iscau_[1])]
                    iscau.append(iscau_)
                if iline == bc_start + 11:
                    iscau_ = line.strip().split()[:1]
                    iscau_ = [float(iscau_[0])]
                    iscau.append(iscau_)
                if iline == bc_start + 13:
                    scauchy_ = line.strip().split()[:3]
                    scauchy_ = [float(scauchy_[0]), float(scauchy_[1]), float(scauchy_[2])]
                    scauchy.append(scauchy_)
                if iline == bc_start + 14:
                    scauchy_ = line.strip().split()[:2]
                    scauchy_ = [float(scauchy_[0]), float(scauchy_[1])]
                    scauchy.append(scauchy_)
                if iline == bc_start + 15:
                    scauchy_ = line.strip().split()[:1]
                    scauchy_ = [float(scauchy_[0])]
                    scauchy.append(scauchy_)
                # Get solver information
                if "other" in line:
                    other_flag = 1
                    other_start = iline
                if iline == other_start + 1:
                    eqincr = line.strip().split()[:1]
                    self.INBFdt.setText(eqincr[0])
                if "INFORMATION ABOUT RUN CONDITIONS" in line:
                    settings_flag = 1
                    settings_start = iline
                if iline == settings_start + 1:
                    nsteps = line.strip().split()[:1]
                    self.INBFloadsteps.setText(nsteps[0])
                if iline == settings_start + 9:
                    iwstep = line.strip().split()[:2]
                    self.INBFoutputsteps.setText(iwstep[1])
                if iline == settings_start + 2:
                    errevp = line.strip().split()[:1]
                    self.INBFerrortol.setText(errevp[0])
                if iline == settings_start + 3:
                    itmax = line.strip().split()[:1]
                    self.INBFmaxiter.setText(itmax[0])
                iline += 1
        if bc_flag == 0:
            warning_message("ERROR: couldn't locate line \n'* boundary conditions'\n indicating the line before the start of the boundary conditions within the input file.")
        elif other_flag == 0:
            warning_message("ERROR: couldn't locate line \n'* other'\n indicating the line before the start of the step information within the input file")
        elif settings_flag == 0:
            warning_message("ERROR: couldn't locate line \n'* INFORMATION ABOUT RUN CONDITIONS'\n indicating the line before the starts of the run condition information within the input file.")
        else:
            # Fill out velocity gradient BCs
            for i in range(3):
                for j in range(3):
                    if iudot[i][j] == 1:
                        self.TVgrad.setItem(i,j,QTableWidgetItem(str(udot[i][j])))
                    else:
                        self.TVgradi.setItem(i,j,QTableWidgetItem(str(udot[i][j])))
            # Fill out Cauchy stress BCs
            for i in range(3):
                for j in range(3-i):
                    if iscau[i][j] == 1:
                        self.TCstress.setItem(i,j+i,QTableWidgetItem(str(scauchy[i][j])))
    self.BLegacyEVPFFT.clicked.connect(legacyEVPFFT)
        
    
