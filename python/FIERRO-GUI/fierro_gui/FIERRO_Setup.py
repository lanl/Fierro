import sys
import os
import tempfile
from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton, QLabel, QVBoxLayout, QDialog, QFileDialog, QDialogButtonBox, QMessageBox, QRadioButton
from PySide6.QtCore import Qt, Signal
from ui_FIERRO_Setup import Ui_Dialog
from DeveloperInputs import *

# ===============================================
# ========== CREATE WORKING DIRECTORY ===========
# ===============================================

class FierroSetup(QDialog, Ui_Dialog):
    settingChanged = Signal(str)

    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.setWindowFlags(Qt.Window | Qt.WindowStaysOnTopHint)
        self.setupUi(self)
        self.MainWindow = main_window
        
        # Connect to user or development menu
        self.BUser.clicked.connect(lambda: self.Configuration.setCurrentIndex(0))
        self.BDeveloper.clicked.connect(lambda: self.Configuration.setCurrentIndex(1))
        
        # Set initial directory to temporary directory
        self.directory = self.create_temp_directory()
        self.BSelectWD.clicked.connect(self.select_directory)

        # Create confirm button
        self.BConfirm.clicked.connect(self.confirm)
        
        # Initially fill developer executables from input file
        self.INFierroMeshBuilder.setText(f"{fierro_mesh_builder_exe}")
        self.INFierroVoxelizer.setText(f"{fierro_voxelizer_exe}")
        self.INFierroParallelExplicit.setText(f"{fierro_parallel_explicit_exe}")
        self.INFierroEvpfft.setText(f"{fierro_evpfft_exe}")
#        self.INFierroMeshBuilder.setStyleSheet("color: blue")
#        self.INFierroVoxelizer.setStyleSheet("color: blue")
#        self.INFierroParallelExplicit.setStyleSheet("color: blue")
#        self.INFierroEvpfft.setStyleSheet("color: blue")
        
        # Select build executables
        self.BFierroMeshBuilder.clicked.connect(lambda: self.developer_executable("fierro-mesh-builder", self.INFierroMeshBuilder))
        self.BFierroVoxelizer.clicked.connect(lambda: self.developer_executable("fierro-voxelizer", self.INFierroVoxelizer))
        self.BFierroParallelExplicit.clicked.connect(lambda: self.developer_executable("fierro-parallel-explicit", self.INFierroParallelExplicit))
        self.BFierroEvpfft.clicked.connect(lambda: self.developer_executable("evpfft", self.INFierroEvpfft))
        
    def create_temp_directory(self):
        temp_dir = os.path.join(tempfile.gettempdir(), "fierro_temp_directory")
        os.makedirs(temp_dir, exist_ok=True)
        return temp_dir
    
    def select_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Working Directory", "")
        if directory:
            self.directory = directory
            self.INCurrentWD.setText(f"{self.directory}")
#            self.INCurrentWD.setStyleSheet("color: blue")
    
    def get_directory(self):
        return self.directory
        
    def closeEvent(self, event):
        message = QMessageBox()
        message.setText("ERROR: Please confirm working directory and profile configuration")
        message.exec()
        event.ignore()
        
    def confirm(self):
        if self.BDeveloper.isChecked():
            # Write executable locations to file for future use
            current_dir = os.path.dirname(os.path.abspath(__file__))
            executable_paths_file = os.path.join(current_dir, "DeveloperInputs.py")
            with open(executable_paths_file, "w") as file:
                file.write(f"fierro_mesh_builder_exe = '{self.INFierroMeshBuilder.text()}'\n")
                file.write(f"fierro_voxelizer_exe = '{self.INFierroVoxelizer.text()}'\n")
                file.write(f"fierro_parallel_explicit_exe = '{self.INFierroParallelExplicit.text()}'\n")
                file.write(f"fierro_evpfft_exe = '{self.INFierroEvpfft.text()}'\n")
                
        # Emit which type of user configuration is selected back to the main window
        if self.BUser.isChecked():
            self.settingChanged.emit("User")
        elif self.BDeveloper.isChecked():
            self.settingChanged.emit("Developer")
        
        # Re-enable the main window
        if hasattr(self, 'MainWindow') and hasattr(self.MainWindow, 'setEnabled'):
            self.MainWindow.setEnabled(True)
        self.accepted.emit()
        self.setResult(QDialog.Accepted)
        if self.isVisible():
            self.hide()
            
    def is_unix_executable(self, file_path):
        return os.access(file_path, os.X_OK)
        
    def is_valid_file_extension(self, file_path):
        return file_path.lower().endswith((".exe", ".exec"))
    
    def developer_executable(self, executable_name, input_var):
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFile)
        file_dialog.setNameFilter("Unix Executable Files (*)")
        file_dialog.setViewMode(QFileDialog.List)
        file_dialog.setOption(QFileDialog.ReadOnly, True)

        if file_dialog.exec_():
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                selected_file = selected_files[0]
                if not self.is_unix_executable(selected_file) or self.is_valid_file_extension(selected_file):
                    self.warning_message("ERROR: This is not an executable file")
                else:
                    file_name = os.path.basename(selected_file)
                    if file_name == executable_name:
                        input_var.setText(f"{selected_file}")
#                        input_var.setStyleSheet("color: blue")
                    else:
                        self.warning_message(f"WARNING: The executable file selected does not match the expected naming convention. \n\nExpected: {executable_name}   \n\nInput: {file_name}")
                        input_var.setText(f"{selected_file}")
#                        input_var.setStyleSheet("color: red")


    def warning_message(self, msg):
        message = QMessageBox()
        message.setText(msg)
        message.exec()

        
