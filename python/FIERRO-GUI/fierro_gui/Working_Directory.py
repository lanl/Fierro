import sys
import os
import tempfile
from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton, QLabel, QVBoxLayout, QDialog, QFileDialog, QDialogButtonBox, QMessageBox
from PySide6.QtCore import Qt

# ===============================================
# ========== CREATE WORKING DIRECTORY ===========
# ===============================================

class WorkingDirectoryDialog(QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.MainWindow = main_window
        
        self.setWindowTitle("Define Working Directory")
        
        self.directory = self.create_temp_directory()  # Set initial directory to temporary directory

        self.label = QLabel(f"Current Working Directory: <font color='blue'>/temp</font>")
        self.label.setWordWrap(True)
        self.select_button = QPushButton("Select Working Directory")
        self.select_button.clicked.connect(self.select_directory)
        
        self.confirm_button = QPushButton("Confirm")
        self.confirm_button.clicked.connect(self.confirm)

        self.button_box = QDialogButtonBox()
        self.button_box.addButton(self.confirm_button, QDialogButtonBox.ButtonRole.AcceptRole)

        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.select_button)
        layout.addWidget(self.button_box)
        self.setLayout(layout)
        
        self.setFixedSize(400,200)
        
    def create_temp_directory(self):
        temp_dir = os.path.join(tempfile.gettempdir(), "fierro_temp_directory")
        os.makedirs(temp_dir, exist_ok=True)
        return temp_dir
    
    def select_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Working Directory", "")
        if directory:
            self.directory = directory
            self.label.setText(f"Current Working Directory: <font color='blue'>{self.directory}</font>")
    
    def get_directory(self):
        return self.directory
        
    def closeEvent(self, event):
        message = QMessageBox()
        message.setText("You must confirm your working directory")
        message.exec()
        event.ignore()
        
    def confirm(self):
        if hasattr(self, 'MainWindow') and hasattr(self.MainWindow, 'setEnabled'):
            self.MainWindow.setEnabled(True)
        self.accepted.emit()
        self.setResult(QDialog.Accepted)
        if self.isVisible():
            self.hide()
        
