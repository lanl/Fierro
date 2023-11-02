import sys
import os

def main():
    from evpfft_gui.ui_EVPFFT_GUI import Ui_MainWindow
    from PySide6.QtWidgets import QMainWindow, QApplication
    # MAIN WINDOW CLASS
    class MainWindow(QMainWindow):
        def __init__(self):
            QMainWindow.__init__(self)
            self.ui = Ui_MainWindow()
            self.ui.setupUi(self)
            
            # SHOW WINDOW
            self.show()
            
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
