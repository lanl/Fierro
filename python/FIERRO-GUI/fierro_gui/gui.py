import sys
import os
import time

def main():
    t1 = time.perf_counter()
    sys.path.append(os.path.dirname(__file__))
    from fierro_gui.FIERRO_GUI import FIERRO_GUI
    from PySide6.QtWidgets import QMainWindow, QApplication, QSplashScreen
    from PySide6 import QtCore
    from PySide6.QtGui import QPixmap
    from fierro_gui.FIERRO_Setup import FierroSetup
    
    # MAIN WINDOW CLASS
    class MainWindow(QMainWindow):
        def __init__(self):
            super(MainWindow, self).__init__()
            self.setWindowTitle("Fierro")
            self.ui = FIERRO_GUI()
            self.ui.setupUi(self)
            self.setEnabled(False)
            
            # SHOW MAIN WINDOW
            self.show()
            
            # SHOW FIERRO SETUP WINDOW
            t2 = time.perf_counter()
            print(f"Started up in {t2 - t1:0.4f} seconds")
            self.ui.open_fierro_setup_dialog(self)
            
#            self.dialog = WorkingDirectoryDialog(self)
#            self.dialog.exec()
            
    app = QApplication(sys.argv)
    start = time.perf_counter()
    print ("Splash Screen")
    pixmap = QPixmap(":/Logos/Logos/FIERRO.png")
    pixmap = pixmap.scaledToWidth(500)
    print (pixmap)
    splash = QSplashScreen(pixmap)
    print (splash)
    splash.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint)
    #splash.showMessage("Loaded modules")
    splash.show()
    for i in range(100000):
        app.processEvents()
    window = MainWindow()
    window.show()
    splash.finish(window)
    
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
