import sys
import os

def main():
    sys.path.append(os.path.dirname(__file__))
    from fierro_gui.FIERRO_GUI import FIERRO_GUI
    from PySide6.QtWidgets import QMainWindow, QApplication
    # MAIN WINDOW CLASS
    class MainWindow(QMainWindow):
        def __init__(self):
            QMainWindow.__init__(self)
            self.ui = FIERRO_GUI()
            self.ui.setupUi(self)
            
            # SHOW MAIN WINDOW
            self.show()
            
            # SHOW WORKING DIRECTORY WINDOW
            self.ui.open_working_directory_dialog()
            
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
