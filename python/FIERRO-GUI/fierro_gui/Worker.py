from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)


from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox, QFormLayout,
    QFrame, QGridLayout, QHBoxLayout, QHeaderView,
    QLabel, QLineEdit, QMainWindow, QMenu,
    QMenuBar, QPlainTextEdit, QProgressBar, QPushButton,
    QSizePolicy, QSpacerItem, QSplitter, QStackedWidget,
    QStatusBar, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget)

from fierro_gui.ui_FIERRO_GUI import Ui_MainWindow
from PySide6.QtWidgets import ( QTableWidgetItem, QSplashScreen)

import os, glob#, sys
import time
#import numpy as np
from PySide6.QtWidgets import (QFileDialog, QMessageBox)
#from PySide6.QtCore import (QTimer, QProcess)
import matplotlib
matplotlib.use('qt5agg')
#from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
#from matplotlib.figure import Figure
#from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
#import csv
#import re
import shutil
import paraview.simple as pvsimple
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PySide6.QtCore import *
import PySide6.QtConcurrent
from PySide6.QtGui import QDesktopServices, QMovie, QPixmap
from PySide6.QtCore import QRunnable, Slot, QThreadPool
#import fierro_voxelizer
#import fierro_mesh_builder
#import tempfile
#import time
import subprocess
from importlib import reload

from Geometry import *
from EVPFFT_Lattice import *
from Mesh_Builder_WInput import *
from FIERRO_Setup import *
from ImageToVTK import *
from TiffImageToVTK import *
import DeveloperInputs

class Worker(QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = None

    @Slot()  # QtCore.Slot
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        self.fn()