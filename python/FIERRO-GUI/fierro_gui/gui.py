import sys
import os
import time
from PySide6.QtWidgets import QApplication, QMessageBox
from PySide6.QtGui import QPalette, QColor
from PySide6.QtCore import Qt, QEvent
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def main():
#    t1 = time.perf_counter()
    sys.path.append(os.path.dirname(__file__))
    from fierro_gui.FIERRO_GUI import FIERRO_GUI
    from PySide6.QtWidgets import QMainWindow, QApplication, QSplashScreen
    from PySide6 import QtCore
    from PySide6.QtGui import QPixmap, QIcon
    from fierro_gui.FIERRO_Setup import FierroSetup
    
    # MAIN WINDOW CLASS
    class MainWindow(QMainWindow):
        # Add restart capabilities
        def restart(self):
            custom_close_event = QEvent(QEvent.Close)
            self.closeEvent(custom_close_event)
            if self.close_accepted:
                self.is_restarting = True
                self.close()
                new_window = MainWindow()
                new_window.show()
        def setup_restart(self):
            self.ui.actionNew.triggered.connect(self.restart)
            
        # Add save changes pop up menu
        def closeEvent(self, event):
            if not self.is_restarting:
                reply = QMessageBox.question(self, 'Save Changes?', 'Do you want to save your changes before closing?', QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Save)
                if reply == QMessageBox.Save:
                    self.ui.open_save_dialog()
                    event.accept()
                    self.close_accepted = True
                elif reply == QMessageBox.Discard:
                    event.accept()
                    self.close_accepted = True
                else:
                    event.ignore()
                    self.close_accepted = False
            else:
                event.accept()
        
        # Run main GUI
        def __init__(self):
            super(MainWindow, self).__init__()
            self.setWindowTitle("Fierro")
            self.ui = FIERRO_GUI()
            self.ui.setupUi(self)
            self.setEnabled(False)
            
            # add restart
            self.setup_restart()
            self.close_accepted = False
            self.is_restarting = False
            
            # SHOW MAIN WINDOW
            self.show()
            
            # SHOW FIERRO SETUP WINDOW
#            t2 = time.perf_counter()
#            print(f"Started up in {t2 - t1:0.4f} seconds")
            self.ui.open_fierro_setup_dialog(self)
            
#            self.dialog = WorkingDirectoryDialog(self)
#            self.dialog.exec()

    app = QApplication(sys.argv)
    force_light_mode(app)
    app.setStyle('Fusion')
    app.setWindowIcon(QIcon(':/Logos/Logos/FierroAppIcon.png'))
    app.setApplicationDisplayName("Fierro")
#    start = time.perf_counter()
#    print ("Splash Screen")
#    pixmap = QPixmap(":/Logos/Logos/FIERRO.png")
#    pixmap = pixmap.scaledToWidth(500)
#    print (pixmap)
#    splash = QSplashScreen(pixmap)
#    print (splash)
#    splash.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint)
    #splash.showMessage("Loaded modules")
#    splash.show()
    #for i in range(100000):
        #app.processEvents()
    window = MainWindow()
    window.show()
#    splash.finish(window)
    
    sys.exit(app.exec())
    
def force_light_mode(app):
    palette = QPalette()
    
    # Set colors for active (enabled) state
    for state in [QPalette.Active, QPalette.Inactive, QPalette.Disabled]:
        palette.setColor(state, QPalette.Window, QColor(240, 240, 240))
        palette.setColor(state, QPalette.WindowText, Qt.black)
        palette.setColor(state, QPalette.Base, Qt.white)
        palette.setColor(state, QPalette.AlternateBase, QColor(233, 231, 227))
        palette.setColor(state, QPalette.ToolTipBase, Qt.white)
        palette.setColor(state, QPalette.ToolTipText, Qt.black)
        palette.setColor(state, QPalette.Text, Qt.black)
        palette.setColor(state, QPalette.Button, QColor(240, 240, 240))
        palette.setColor(state, QPalette.ButtonText, Qt.black)
        palette.setColor(state, QPalette.BrightText, Qt.red)
        palette.setColor(state, QPalette.Link, QColor(42, 130, 218))
        palette.setColor(state, QPalette.Highlight, QColor(42, 130, 218))
        palette.setColor(state, QPalette.HighlightedText, Qt.black)

    # Set colors for disabled state
    grey_color = QColor(128, 128, 128)  # Medium grey
    palette.setColor(QPalette.Disabled, QPalette.WindowText, grey_color)
    palette.setColor(QPalette.Disabled, QPalette.Text, grey_color)
    palette.setColor(QPalette.Disabled, QPalette.ButtonText, grey_color)
    palette.setColor(QPalette.Disabled, QPalette.Highlight, Qt.transparent)
    palette.setColor(QPalette.Disabled, QPalette.HighlightedText, grey_color)

    # You can adjust the grey color or use a slightly different color for disabled items
    disabled_background = QColor(240, 240, 240)  # Light grey background for disabled items
    palette.setColor(QPalette.Disabled, QPalette.Button, disabled_background)
    palette.setColor(QPalette.Disabled, QPalette.Base, disabled_background)
    palette.setColor(QPalette.Disabled, QPalette.Window, disabled_background)

    app.setPalette(palette)

if __name__ == "__main__":
    main()
