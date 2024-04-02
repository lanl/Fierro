import sys
import os

def main():
    sys.path.append(os.path.dirname(__file__))
    from fierro_gui.FIERRO_GUI import FIERRO_GUI
    from PySide6.QtWidgets import QMainWindow, QApplication
    from FIERRO_Setup import FierroSetup
    # MAIN WINDOW CLASS
    class MainWindow(QMainWindow):
        def __init__(self):
            QMainWindow.__init__(self)
            self.ui = FIERRO_GUI()
            self.ui.setupUi(self)
            self.setEnabled(False)
            
            # SHOW MAIN WINDOW
            self.show()
            
            # SHOW FIERRO SETUP WINDOW
            self.ui.open_fierro_setup_dialog(self)
#            self.dialog = WorkingDirectoryDialog(self)
#            self.dialog.exec()
            
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
