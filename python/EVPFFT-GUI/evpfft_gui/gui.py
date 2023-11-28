import sys
import os

def main():
    from evpfft_gui.EVPFFT_GUI import EVPFFT_GUI
    from PySide6.QtWidgets import QMainWindow, QApplication
    # MAIN WINDOW CLASS
    class MainWindow(QMainWindow):
        def __init__(self):
            QMainWindow.__init__(self)
            self.ui = EVPFFT_GUI()
            self.ui.setupUi(self)
            
            # SHOW WINDOW
            self.show()
            
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
