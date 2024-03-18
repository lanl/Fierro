from PySide6.QtWidgets import (QTableWidgetItem, QMessageBox)
from PySide6.QtCore import (QProcess)
from Explicit_SGH_WInput import *

# ============================================
# ======= EXPLICIT SOLVER SGH PIPELINE =======
# ============================================

# Warning Message Popup
def warning_message(msg):
    message = QMessageBox()
    message.setText(msg)
    message.exec()

def Explicit_SGH(self):
    # Connect tab buttons to settings windows
    self.BGlobalMesh.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(1))
    self.BImportPartSGH.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(2))
    self.BDefineMaterialSGH.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(7))
    self.BAssignMaterialSGH.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(8))
    self.BApplyBCSGH.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(9))
    self.BSolverSettingsSGH.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(10))

    # Material Definition
    def add_material_SGH():
        if self.INEOS.currentText() == "None" and not self.INMaterialNameSGH.text():
            warning_message("ERROR: Material definition incomplete")
        elif self.INEOS.currentText() != "None" and (not self.INMaterialNameSGH.text() or not self.INq1.text() or not self.INq2.text() or not self.INq1ex.text() or not self.INq2ex.text() or not self.INGamma.text() or not self.INMinSound.text() or not self.INSpecificHeat.text()):
            warning_message("ERROR: Material definition incomplete")
        else:
            row = self.TMaterialsSGH.rowCount()
            self.TMaterialsSGH.insertRow(row)
            self.TMaterialsSGH.setItem(row, 0, QTableWidgetItem(
                self.INMaterialNameSGH.text().strip())
            )
            self.TMaterialsSGH.setItem(
                row, 1, QTableWidgetItem(str(self.INEOS.currentText()))
            )
            self.TMaterialsSGH.setItem(
                row, 2, QTableWidgetItem(self.INq1.text())
            )
            self.TMaterialsSGH.setItem(
                row, 3, QTableWidgetItem(self.INq2.text())
            )
            self.TMaterialsSGH.setItem(
                row, 4, QTableWidgetItem(self.INq1ex.text())
            )
            self.TMaterialsSGH.setItem(
                row, 5, QTableWidgetItem(self.INq2ex.text())
            )
            self.TMaterialsSGH.setItem(
                row, 6, QTableWidgetItem(self.INGamma.text())
            )
            self.TMaterialsSGH.setItem(
                row, 7, QTableWidgetItem(self.INMinSound.text())
            )
            self.TMaterialsSGH.setItem(
                row, 8, QTableWidgetItem(self.INSpecificHeat.text())
            )
            self.TMaterialsSGH.setItem(
                row, 9, QTableWidgetItem(str(self.INStrengthModel.currentText()))
            )
            
            # Add material as an option for material assignment
            self.INMaterial.clear()
            for i in range(self.TMaterialsSGH.rowCount()):
                self.INMaterial.addItem(self.TMaterialsSGH.item(i,0).text())
            
            # Clear all inputs
            reset_material_SGH()
        
    def reset_material_SGH():
        self.INMaterialNameSGH.clear()
        self.INq1.clear()
        self.INq2.clear()
        self.INq1ex.clear()
        self.INq2ex.clear()
        self.INGamma.clear()
        self.INMinSound.clear()
        self.INSpecificHeat.clear()
            
    def delete_material_SGH():
        current_row = self.TMaterialsSGH.currentRow()
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
            # delete from table
            self.TMaterialsSGH.removeRow(current_row)
            
            # delete from material assignment options
            self.INMaterial.clear()
            for i in range(self.TMaterialsSGH.rowCount()):
                self.INMaterial.addItem(self.TMaterialsSGH.item(i,0).text())
            
    def EOS_Model():
        if str(self.INEOS.currentText()) == 'None':
            self.Lq1.setEnabled(False)
            self.INq1.setEnabled(False)
            self.Lq2.setEnabled(False)
            self.INq2.setEnabled(False)
            self.Lq1ex.setEnabled(False)
            self.INq1ex.setEnabled(False)
            self.Lq2ex.setEnabled(False)
            self.INq2ex.setEnabled(False)
            self.LGamma.setEnabled(False)
            self.INGamma.setEnabled(False)
            self.LMinSound.setEnabled(False)
            self.INMinSound.setEnabled(False)
            self.LSpecificHeat.setEnabled(False)
            self.INSpecificHeat.setEnabled(False)
        if str(self.INEOS.currentText()) == 'ideal_gas':
            self.Lq1.setEnabled(True)
            self.INq1.setEnabled(True)
            self.Lq2.setEnabled(True)
            self.INq2.setEnabled(True)
            self.Lq1ex.setEnabled(True)
            self.INq1ex.setEnabled(True)
            self.Lq2ex.setEnabled(True)
            self.INq2ex.setEnabled(True)
            self.LGamma.setEnabled(True)
            self.INGamma.setEnabled(True)
            self.LMinSound.setEnabled(True)
            self.INMinSound.setEnabled(True)
            self.LSpecificHeat.setEnabled(True)
            self.INSpecificHeat.setEnabled(True)
            
    self.INEOS.currentIndexChanged.connect(EOS_Model)
    self.BAddMaterialSGH.clicked.connect(add_material_SGH)
    self.BDeleteMaterialSGH.clicked.connect(delete_material_SGH)

    # Material assignment
    def assign_material_SGH():
        if not self.INDensity.text() or not self.INSIE.text() or not self.INVelocityX.text() or not self.INVelocityY.text() or not self.INVelocityZ.text():
            warning_message("ERROR: Material assignment incomplete")
        else:
            row = self.Tassignmat.rowCount()
            self.Tassignmat.insertRow(row)
            self.Tassignmat.setItem(row, 0, QTableWidgetItem(str(self.INPartMaterial.currentText()))
            )
            self.Tassignmat.setItem(
                row, 1, QTableWidgetItem(str(self.INMaterial.currentText()))
            )
            self.Tassignmat.setItem(
                row, 2, QTableWidgetItem(self.INDensity.text())
            )
            self.Tassignmat.setItem(
                row, 3, QTableWidgetItem(self.INSIE.text())
            )
            self.Tassignmat.setItem(
                row, 4, QTableWidgetItem(self.INVelocityX.text())
            )
            self.Tassignmat.setItem(
                row, 5, QTableWidgetItem(self.INVelocityY.text())
            )
            self.Tassignmat.setItem(
                row, 6, QTableWidgetItem(self.INVelocityZ.text())
            )
            
            reset_material_assign_SGH()

    def reset_material_assign_SGH():
        self.INDensity.clear()
        self.INSIE.clear()
        self.INVelocityX.clear()
        self.INVelocityY.clear()
        self.INVelocityZ.clear()

    def delete_assign_mat_SGH():
        current_row = self.Tassignmat.currentRow()
        if current_row < 0:
            return QMessageBox.warning(QMessageBox(), 'Warning','Please select a record to delete')

        button = QMessageBox.question(
            QMessageBox(),
            'Confirmation',
            'Are you sure that you want to delete the selected row?',
            QMessageBox.StandardButton.Yes |
            QMessageBox.StandardButton.No
        )
        if button == QMessageBox.StandardButton.Yes:
            self.Tassignmat.removeRow(current_row)
        
    self.Baddmaterialassignment.clicked.connect(assign_material_SGH)
    self.Bdeletematerialassignment.clicked.connect(delete_assign_mat_SGH)

    # Boundary Condition Assignment
    def selected_bcs(text):
        if text == "velocity":
            self.INVel0.setEnabled(True)
            self.INVel1.setEnabled(True)
            self.INVelstart.setEnabled(True)
            self.INVelend.setEnabled(True)
            self.LVel0.setEnabled(True)
            self.LVel1.setEnabled(True)
            self.LVelstart.setEnabled(True)
            self.LVelend.setEnabled(True)
        else:
            self.INVel0.setEnabled(False)
            self.INVel1.setEnabled(False)
            self.INVelstart.setEnabled(False)
            self.INVelend.setEnabled(False)
            self.LVel0.setEnabled(False)
            self.LVel1.setEnabled(False)
            self.LVelstart.setEnabled(False)
            self.LVelend.setEnabled(False)
    self.INType.currentTextChanged.connect(selected_bcs)
        
    def add_bcs_SGH():
        if not self.INPlanePosition.text():
            warning_message("ERROR: Boundary condition incomplete")
        elif self.INType.currentText() == "velocity" and (not self.INPlanePosition.text() or not self.INVel0.text() or not self.INVel1.text() or not self.INVelstart.text() or not self.INVelend.text()):
                                warning_message("ERROR: Boundary condition incomplete")
        else:
            row = self.TBoundaryConditions.rowCount()
            self.TBoundaryConditions.insertRow(row)
            self.TBoundaryConditions.setItem(row, 0, QTableWidgetItem(str(
                self.INBoundary.currentText()))
            )
            self.TBoundaryConditions.setItem(
                row, 1, QTableWidgetItem(self.INPlanePosition.text())
            )
            self.TBoundaryConditions.setItem(
                row, 2, QTableWidgetItem(str(self.INType.currentText()))
            )
            self.TBoundaryConditions.setItem(
                row, 3, QTableWidgetItem(self.INVel0.text())
            )
            self.TBoundaryConditions.setItem(
                row, 4, QTableWidgetItem(self.INVel1.text())
            )
            self.TBoundaryConditions.setItem(
                row, 5, QTableWidgetItem(self.INVelstart.text())
            )
            self.TBoundaryConditions.setItem(
                row, 6, QTableWidgetItem(self.INVelend.text())
            )

            reset_bcs_SGH()
        
    def reset_bcs_SGH():
        self.INPlanePosition.clear()
        self.INVel0.clear()
        self.INVel1.clear()
        self.INVelstart.clear()
        self.INVelend.clear()

    def delete_bcs_SGH():
        current_row = self.TBoundaryConditions.currentRow()
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
            self.TBoundaryConditions.removeRow(current_row)
            
    self.BaddBC.clicked.connect(add_bcs_SGH)
    self.BdeleteBC.clicked.connect(delete_bcs_SGH)
    
    def run_explicit_SGH():
        Explicit_SGH_WInput(self)
        import subprocess
        executable_path = "/Users/shankins/Documents/FY24/Github/XcodeFierro/Fierro/build-fierro-serial/bin/fierro-parallel-explicit"
        arguments = [self.EXPLICIT_SGH_INPUT]
        command = [executable_path] + arguments
        process = subprocess.Popen(command)
        process.wait()
        self.explicit_sgh = QProcess()
    self.BRunSGH.clicked.connect(run_explicit_SGH)
#        self.explicit_sgh.start("fierro-mesh-builder","/Users/shankins/Documents/FY24/Github/XcodeFierro/Fierro/python/FIERRO-GUI/global_mesh.yaml")
