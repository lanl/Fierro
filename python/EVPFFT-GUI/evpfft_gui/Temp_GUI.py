# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'EVPFFT_GUINIjlCa.ui'
##
## Created by: Qt User Interface Compiler version 6.5.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QAction, QBrush, QColor, QConicalGradient,
    QCursor, QFont, QFontDatabase, QGradient,
    QIcon, QImage, QKeySequence, QLinearGradient,
    QPainter, QPalette, QPixmap, QRadialGradient,
    QTransform)
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox, QFormLayout,
    QFrame, QGridLayout, QHBoxLayout, QHeaderView,
    QLabel, QLineEdit, QMainWindow, QMenu,
    QMenuBar, QPlainTextEdit, QProgressBar, QPushButton,
    QSizePolicy, QSpacerItem, QSplitter, QStackedWidget,
    QStatusBar, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget)
    
# ================== New imports =====================
    
import fierro_voxelizer
import tempfile
import os
    
class LocalResource:
    FILE_PATH = os.path.abspath(
        os.path.join(*(os.path.split(os.path.expanduser(__file__))[:-1]))
    )

    @staticmethod
    def get_resource_name(relpath: str) -> str:
        return os.path.join(LocalResource.FILE_PATH, relpath)
    

VTK_OUTPUT = os.path.join(tempfile.gettempdir(), 'VTK_Geometry.vtk')
ELASTIC_PARAMETERS = LocalResource.get_resource_name('elastic_parameters.txt')
PLASTIC_PARAMETERS = LocalResource.get_resource_name('plastic_parameters.txt')
EVPFFT_INPUT = os.path.join(tempfile.gettempdir(), 'evpfft_lattice_input.txt')

# ======================================================

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1022, 883)
        icon = QIcon()
        icon.addFile(u"Icons/EVPFFT_logo_A2.png", QSize(), QIcon.Normal, QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setAutoFillBackground(False)
        MainWindow.setToolButtonStyle(Qt.ToolButtonIconOnly)
        MainWindow.setDockNestingEnabled(False)
        self.actionEVPFFT_Manual = QAction(MainWindow)
        self.actionEVPFFT_Manual.setObjectName(u"actionEVPFFT_Manual")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.centralwidget.setStyleSheet(u"#TitlePage, #GeometryInformationTool, #DefineMaterialTool, #BoundaryConditionsTool, #SolverSettingsTool, #ResultsTool, #Tools, #RunOutputs, #RunOutputWindow, #Main{\n"
"	background-color: rgb(235, 235, 235);\n"
"}\n"
"#ParaviewFrame{\n"
"	background-color: rgb(91, 97, 120);\n"
"}\n"
"#BImportPart:hover, #BDefineMaterial:hover, #BApplyBC:hover, #BSolverSettings:hover, #BRunEVPFFT:hover, #BViewResults:hover{\n"
"	background-color: rgb(192, 192, 192);\n"
"	border-radius: 15px;\n"
"}\n"
"#BImportPart, #BDefineMaterial, #BApplyBC, #BSolverSettings, #BRunEVPFFT, #BViewResults{\n"
"	border-style: flat;\n"
"}\n"
"#centralwidget{\n"
"	background-color: rgb(255, 255, 255);\n"
"}\n"
"\n"
"")
        self.verticalLayout = QVBoxLayout(self.centralwidget)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.HeaderMenu = QTabWidget(self.centralwidget)
        self.HeaderMenu.setObjectName(u"HeaderMenu")
        self.HeaderMenu.setMaximumSize(QSize(16777215, 125))
        font = QFont()
        font.setBold(True)
        self.HeaderMenu.setFont(font)
        self.Tools = QWidget()
        self.Tools.setObjectName(u"Tools")
        font1 = QFont()
        font1.setPointSize(13)
        font1.setBold(False)
        self.Tools.setFont(font1)
        self.horizontalLayout_2 = QHBoxLayout(self.Tools)
#ifndef Q_OS_MAC
        self.horizontalLayout_2.setSpacing(-1)
#endif
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.horizontalLayout_2.setContentsMargins(-1, 0, 5, 0)
        self.PartTools = QFrame(self.Tools)
        self.PartTools.setObjectName(u"PartTools")
        self.PartTools.setFrameShape(QFrame.NoFrame)
        self.PartTools.setFrameShadow(QFrame.Raised)
        self.verticalLayout_7 = QVBoxLayout(self.PartTools)
        self.verticalLayout_7.setSpacing(0)
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.verticalLayout_7.setContentsMargins(0, 0, 0, 0)
        self.ImportPart = QFrame(self.PartTools)
        self.ImportPart.setObjectName(u"ImportPart")
        self.ImportPart.setMinimumSize(QSize(70, 80))
        self.ImportPart.setMaximumSize(QSize(70, 80))
        self.ImportPart.setFrameShape(QFrame.NoFrame)
        self.ImportPart.setFrameShadow(QFrame.Raised)
        self.verticalLayout_3 = QVBoxLayout(self.ImportPart)
        self.verticalLayout_3.setSpacing(0)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 5)
        self.BImportPart = QPushButton(self.ImportPart)
        self.BImportPart.setObjectName(u"BImportPart")
        icon1 = QIcon()
        icon1.addFile(u"Icons/Blue/Cube.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BImportPart.setIcon(icon1)
        self.BImportPart.setIconSize(QSize(50, 50))
        self.BImportPart.setFlat(False)

        self.verticalLayout_3.addWidget(self.BImportPart)

        self.LImportPart = QLabel(self.ImportPart)
        self.LImportPart.setObjectName(u"LImportPart")
        self.LImportPart.setWordWrap(True)

        self.verticalLayout_3.addWidget(self.LImportPart, 0, Qt.AlignBottom)


        self.verticalLayout_7.addWidget(self.ImportPart)

        self.LinePartTools = QFrame(self.PartTools)
        self.LinePartTools.setObjectName(u"LinePartTools")
        self.LinePartTools.setFrameShape(QFrame.HLine)
        self.LinePartTools.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_7.addWidget(self.LinePartTools)

        self.LPartTools = QLabel(self.PartTools)
        self.LPartTools.setObjectName(u"LPartTools")

        self.verticalLayout_7.addWidget(self.LPartTools)


        self.horizontalLayout_2.addWidget(self.PartTools)

        self.line = QFrame(self.Tools)
        self.line.setObjectName(u"line")
        self.line.setFrameShape(QFrame.VLine)
        self.line.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout_2.addWidget(self.line)

        self.MaterialTools = QFrame(self.Tools)
        self.MaterialTools.setObjectName(u"MaterialTools")
        self.MaterialTools.setFrameShape(QFrame.NoFrame)
        self.MaterialTools.setFrameShadow(QFrame.Raised)
        self.verticalLayout_8 = QVBoxLayout(self.MaterialTools)
        self.verticalLayout_8.setSpacing(0)
        self.verticalLayout_8.setObjectName(u"verticalLayout_8")
        self.verticalLayout_8.setContentsMargins(0, 0, 0, 0)
        self.MaterialToolsBtns = QFrame(self.MaterialTools)
        self.MaterialToolsBtns.setObjectName(u"MaterialToolsBtns")
        self.MaterialToolsBtns.setFrameShape(QFrame.NoFrame)
        self.MaterialToolsBtns.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_4 = QHBoxLayout(self.MaterialToolsBtns)
        self.horizontalLayout_4.setSpacing(0)
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.DefineMaterial = QFrame(self.MaterialToolsBtns)
        self.DefineMaterial.setObjectName(u"DefineMaterial")
        self.DefineMaterial.setMinimumSize(QSize(80, 80))
        self.DefineMaterial.setMaximumSize(QSize(80, 80))
        self.DefineMaterial.setFrameShape(QFrame.NoFrame)
        self.DefineMaterial.setFrameShadow(QFrame.Raised)
        self.verticalLayout_4 = QVBoxLayout(self.DefineMaterial)
        self.verticalLayout_4.setSpacing(0)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 5)
        self.BDefineMaterial = QPushButton(self.DefineMaterial)
        self.BDefineMaterial.setObjectName(u"BDefineMaterial")
        icon2 = QIcon()
        icon2.addFile(u"Icons/Blue/mine.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BDefineMaterial.setIcon(icon2)
        self.BDefineMaterial.setIconSize(QSize(50, 50))
        self.BDefineMaterial.setFlat(True)

        self.verticalLayout_4.addWidget(self.BDefineMaterial)

        self.LDefineMaterial = QLabel(self.DefineMaterial)
        self.LDefineMaterial.setObjectName(u"LDefineMaterial")
        self.LDefineMaterial.setWordWrap(True)

        self.verticalLayout_4.addWidget(self.LDefineMaterial, 0, Qt.AlignBottom)


        self.horizontalLayout_4.addWidget(self.DefineMaterial, 0, Qt.AlignLeft)


        self.verticalLayout_8.addWidget(self.MaterialToolsBtns, 0, Qt.AlignLeft)

        self.LineMaterialTools = QFrame(self.MaterialTools)
        self.LineMaterialTools.setObjectName(u"LineMaterialTools")
        self.LineMaterialTools.setFrameShape(QFrame.HLine)
        self.LineMaterialTools.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_8.addWidget(self.LineMaterialTools)

        self.LMaterialTools = QLabel(self.MaterialTools)
        self.LMaterialTools.setObjectName(u"LMaterialTools")

        self.verticalLayout_8.addWidget(self.LMaterialTools)


        self.horizontalLayout_2.addWidget(self.MaterialTools)

        self.line_2 = QFrame(self.Tools)
        self.line_2.setObjectName(u"line_2")
        self.line_2.setFrameShape(QFrame.VLine)
        self.line_2.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout_2.addWidget(self.line_2)

        self.BCTools = QFrame(self.Tools)
        self.BCTools.setObjectName(u"BCTools")
        self.BCTools.setFrameShape(QFrame.NoFrame)
        self.BCTools.setFrameShadow(QFrame.Raised)
        self.verticalLayout_11 = QVBoxLayout(self.BCTools)
        self.verticalLayout_11.setSpacing(0)
        self.verticalLayout_11.setObjectName(u"verticalLayout_11")
        self.verticalLayout_11.setContentsMargins(0, 0, 0, 0)
        self.ApplyBC = QFrame(self.BCTools)
        self.ApplyBC.setObjectName(u"ApplyBC")
        self.ApplyBC.setMinimumSize(QSize(65, 80))
        self.ApplyBC.setMaximumSize(QSize(70, 80))
        self.ApplyBC.setFrameShape(QFrame.NoFrame)
        self.ApplyBC.setFrameShadow(QFrame.Raised)
        self.verticalLayout_5 = QVBoxLayout(self.ApplyBC)
        self.verticalLayout_5.setSpacing(0)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 5)
        self.BApplyBC = QPushButton(self.ApplyBC)
        self.BApplyBC.setObjectName(u"BApplyBC")
        icon3 = QIcon()
        icon3.addFile(u"Icons/Blue/brick.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BApplyBC.setIcon(icon3)
        self.BApplyBC.setIconSize(QSize(50, 50))
        self.BApplyBC.setFlat(True)

        self.verticalLayout_5.addWidget(self.BApplyBC)

        self.LApplyBC = QLabel(self.ApplyBC)
        self.LApplyBC.setObjectName(u"LApplyBC")
        self.LApplyBC.setWordWrap(True)

        self.verticalLayout_5.addWidget(self.LApplyBC)


        self.verticalLayout_11.addWidget(self.ApplyBC, 0, Qt.AlignLeft)

        self.LineBCTools = QFrame(self.BCTools)
        self.LineBCTools.setObjectName(u"LineBCTools")
        self.LineBCTools.setFrameShape(QFrame.HLine)
        self.LineBCTools.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_11.addWidget(self.LineBCTools)

        self.LBCTools = QLabel(self.BCTools)
        self.LBCTools.setObjectName(u"LBCTools")

        self.verticalLayout_11.addWidget(self.LBCTools)


        self.horizontalLayout_2.addWidget(self.BCTools)

        self.line_3 = QFrame(self.Tools)
        self.line_3.setObjectName(u"line_3")
        self.line_3.setFrameShape(QFrame.VLine)
        self.line_3.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout_2.addWidget(self.line_3)

        self.JobTools = QFrame(self.Tools)
        self.JobTools.setObjectName(u"JobTools")
        self.JobTools.setFrameShape(QFrame.NoFrame)
        self.JobTools.setFrameShadow(QFrame.Raised)
        self.verticalLayout_12 = QVBoxLayout(self.JobTools)
        self.verticalLayout_12.setSpacing(0)
        self.verticalLayout_12.setObjectName(u"verticalLayout_12")
        self.verticalLayout_12.setContentsMargins(0, 0, 0, 0)
        self.JobToolsBtns = QFrame(self.JobTools)
        self.JobToolsBtns.setObjectName(u"JobToolsBtns")
        self.JobToolsBtns.setMinimumSize(QSize(0, 0))
        self.JobToolsBtns.setMaximumSize(QSize(16777215, 16777215))
        self.JobToolsBtns.setFrameShape(QFrame.NoFrame)
        self.JobToolsBtns.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_3 = QHBoxLayout(self.JobToolsBtns)
        self.horizontalLayout_3.setSpacing(0)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.SolverSettings = QFrame(self.JobToolsBtns)
        self.SolverSettings.setObjectName(u"SolverSettings")
        self.SolverSettings.setMinimumSize(QSize(80, 80))
        self.SolverSettings.setMaximumSize(QSize(80, 80))
        self.SolverSettings.setFrameShape(QFrame.NoFrame)
        self.SolverSettings.setFrameShadow(QFrame.Raised)
        self.verticalLayout_6 = QVBoxLayout(self.SolverSettings)
        self.verticalLayout_6.setSpacing(0)
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.verticalLayout_6.setContentsMargins(0, 0, 0, 5)
        self.BSolverSettings = QPushButton(self.SolverSettings)
        self.BSolverSettings.setObjectName(u"BSolverSettings")
        icon4 = QIcon()
        icon4.addFile(u"Icons/Blue/gear.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BSolverSettings.setIcon(icon4)
        self.BSolverSettings.setIconSize(QSize(50, 50))
        self.BSolverSettings.setFlat(True)

        self.verticalLayout_6.addWidget(self.BSolverSettings)

        self.LSolverSettings = QLabel(self.SolverSettings)
        self.LSolverSettings.setObjectName(u"LSolverSettings")
        self.LSolverSettings.setWordWrap(True)

        self.verticalLayout_6.addWidget(self.LSolverSettings)


        self.horizontalLayout_3.addWidget(self.SolverSettings, 0, Qt.AlignLeft)

        self.RunEVPFFT = QFrame(self.JobToolsBtns)
        self.RunEVPFFT.setObjectName(u"RunEVPFFT")
        self.RunEVPFFT.setMinimumSize(QSize(75, 80))
        self.RunEVPFFT.setMaximumSize(QSize(65, 80))
        self.RunEVPFFT.setFrameShape(QFrame.NoFrame)
        self.RunEVPFFT.setFrameShadow(QFrame.Raised)
        self.verticalLayout_13 = QVBoxLayout(self.RunEVPFFT)
        self.verticalLayout_13.setSpacing(0)
        self.verticalLayout_13.setObjectName(u"verticalLayout_13")
        self.verticalLayout_13.setContentsMargins(0, 0, 0, 5)
        self.BRunEVPFFT = QPushButton(self.RunEVPFFT)
        self.BRunEVPFFT.setObjectName(u"BRunEVPFFT")
        icon5 = QIcon()
        icon5.addFile(u"Icons/Blue/Play.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BRunEVPFFT.setIcon(icon5)
        self.BRunEVPFFT.setIconSize(QSize(40, 40))
        self.BRunEVPFFT.setFlat(True)

        self.verticalLayout_13.addWidget(self.BRunEVPFFT)

        self.LRunEVPFFT = QLabel(self.RunEVPFFT)
        self.LRunEVPFFT.setObjectName(u"LRunEVPFFT")
        self.LRunEVPFFT.setWordWrap(True)

        self.verticalLayout_13.addWidget(self.LRunEVPFFT)


        self.horizontalLayout_3.addWidget(self.RunEVPFFT)

        self.ViewResults = QFrame(self.JobToolsBtns)
        self.ViewResults.setObjectName(u"ViewResults")
        self.ViewResults.setMinimumSize(QSize(80, 80))
        self.ViewResults.setMaximumSize(QSize(80, 80))
        self.ViewResults.setFrameShape(QFrame.NoFrame)
        self.ViewResults.setFrameShadow(QFrame.Raised)
        self.verticalLayout_14 = QVBoxLayout(self.ViewResults)
        self.verticalLayout_14.setSpacing(0)
        self.verticalLayout_14.setObjectName(u"verticalLayout_14")
        self.verticalLayout_14.setContentsMargins(0, 0, 0, 5)
        self.BViewResults = QPushButton(self.ViewResults)
        self.BViewResults.setObjectName(u"BViewResults")
        icon6 = QIcon()
        icon6.addFile(u"Icons/Blue/magnify.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BViewResults.setIcon(icon6)
        self.BViewResults.setIconSize(QSize(40, 40))
        self.BViewResults.setFlat(True)

        self.verticalLayout_14.addWidget(self.BViewResults)

        self.LViewResults = QLabel(self.ViewResults)
        self.LViewResults.setObjectName(u"LViewResults")
        self.LViewResults.setWordWrap(True)

        self.verticalLayout_14.addWidget(self.LViewResults)


        self.horizontalLayout_3.addWidget(self.ViewResults)


        self.verticalLayout_12.addWidget(self.JobToolsBtns, 0, Qt.AlignLeft)

        self.LineJobTools = QFrame(self.JobTools)
        self.LineJobTools.setObjectName(u"LineJobTools")
        self.LineJobTools.setFrameShape(QFrame.HLine)
        self.LineJobTools.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_12.addWidget(self.LineJobTools)

        self.LJobTools = QLabel(self.JobTools)
        self.LJobTools.setObjectName(u"LJobTools")

        self.verticalLayout_12.addWidget(self.LJobTools)


        self.horizontalLayout_2.addWidget(self.JobTools)

        self.HeaderMenu.addTab(self.Tools, "")

        self.verticalLayout.addWidget(self.HeaderMenu)

        self.Main = QFrame(self.centralwidget)
        self.Main.setObjectName(u"Main")
        self.Main.setMinimumSize(QSize(0, 0))
        self.Main.setFrameShape(QFrame.Box)
        self.Main.setFrameShadow(QFrame.Plain)
        self.Main.setLineWidth(1)
        self.horizontalLayout = QHBoxLayout(self.Main)
        self.horizontalLayout.setSpacing(0)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.splitter = QSplitter(self.Main)
        self.splitter.setObjectName(u"splitter")
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.splitter.sizePolicy().hasHeightForWidth())
        self.splitter.setSizePolicy(sizePolicy)
        self.splitter.setFrameShape(QFrame.NoFrame)
        self.splitter.setFrameShadow(QFrame.Plain)
        self.splitter.setOrientation(Qt.Horizontal)
        self.splitter.setOpaqueResize(True)
        self.splitter.setHandleWidth(7)
        self.ToolSettings = QStackedWidget(self.splitter)
        self.ToolSettings.setObjectName(u"ToolSettings")
        sizePolicy.setHeightForWidth(self.ToolSettings.sizePolicy().hasHeightForWidth())
        self.ToolSettings.setSizePolicy(sizePolicy)
        self.ToolSettings.setMinimumSize(QSize(300, 0))
        self.ToolSettings.setMaximumSize(QSize(360, 16777215))
        self.ToolSettings.setSizeIncrement(QSize(0, 0))
        self.ToolSettings.setBaseSize(QSize(300, 0))
        self.ToolSettings.setFrameShape(QFrame.NoFrame)
        self.TitlePage = QWidget()
        self.TitlePage.setObjectName(u"TitlePage")
        self.TitlePage.setMinimumSize(QSize(0, 0))
        self.TitlePage.setMaximumSize(QSize(16777215, 16777215))
        self.verticalLayout_2 = QVBoxLayout(self.TitlePage)
        self.verticalLayout_2.setSpacing(40)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.verticalLayout_2.setContentsMargins(12, 40, 12, 40)
        self.LosAlamosLogo = QLabel(self.TitlePage)
        self.LosAlamosLogo.setObjectName(u"LosAlamosLogo")
        self.LosAlamosLogo.setMaximumSize(QSize(16777215, 60))
        self.LosAlamosLogo.setPixmap(QPixmap(u"Icons/LANL Logo Ultramarine.png"))
        self.LosAlamosLogo.setScaledContents(True)

        self.verticalLayout_2.addWidget(self.LosAlamosLogo)

        self.verticalSpacer_7 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_2.addItem(self.verticalSpacer_7)

        self.EVPFFTLogo = QLabel(self.TitlePage)
        self.EVPFFTLogo.setObjectName(u"EVPFFTLogo")
        self.EVPFFTLogo.setMinimumSize(QSize(275, 175))
        self.EVPFFTLogo.setMaximumSize(QSize(275, 175))
        self.EVPFFTLogo.setPixmap(QPixmap(u"Icons/EVPFFT_logo_horse_ppt.png"))
        self.EVPFFTLogo.setScaledContents(True)
        self.EVPFFTLogo.setWordWrap(False)
        self.EVPFFTLogo.setIndent(-1)

        self.verticalLayout_2.addWidget(self.EVPFFTLogo, 0, Qt.AlignHCenter|Qt.AlignVCenter)

        self.verticalSpacer_6 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_2.addItem(self.verticalSpacer_6)

        self.AdditionalSoftware = QFrame(self.TitlePage)
        self.AdditionalSoftware.setObjectName(u"AdditionalSoftware")
        self.AdditionalSoftware.setFrameShape(QFrame.NoFrame)
        self.AdditionalSoftware.setFrameShadow(QFrame.Raised)
        self.verticalLayout_10 = QVBoxLayout(self.AdditionalSoftware)
        self.verticalLayout_10.setSpacing(0)
        self.verticalLayout_10.setObjectName(u"verticalLayout_10")
        self.verticalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.LAdditionalSoftware = QLabel(self.AdditionalSoftware)
        self.LAdditionalSoftware.setObjectName(u"LAdditionalSoftware")
        font2 = QFont()
        font2.setPointSize(16)
        self.LAdditionalSoftware.setFont(font2)

        self.verticalLayout_10.addWidget(self.LAdditionalSoftware, 0, Qt.AlignBottom)

        self.AdditionalSoftwareLogos = QFrame(self.AdditionalSoftware)
        self.AdditionalSoftwareLogos.setObjectName(u"AdditionalSoftwareLogos")
        self.AdditionalSoftwareLogos.setFrameShape(QFrame.NoFrame)
        self.AdditionalSoftwareLogos.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_5 = QHBoxLayout(self.AdditionalSoftwareLogos)
        self.horizontalLayout_5.setSpacing(10)
        self.horizontalLayout_5.setObjectName(u"horizontalLayout_5")
        self.horizontalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.MatarLogo = QLabel(self.AdditionalSoftwareLogos)
        self.MatarLogo.setObjectName(u"MatarLogo")
        self.MatarLogo.setMaximumSize(QSize(120, 60))
        self.MatarLogo.setPixmap(QPixmap(u"Icons/MATAR_logo2.png"))
        self.MatarLogo.setScaledContents(True)

        self.horizontalLayout_5.addWidget(self.MatarLogo)

        self.ParaviewLogo = QLabel(self.AdditionalSoftwareLogos)
        self.ParaviewLogo.setObjectName(u"ParaviewLogo")
        self.ParaviewLogo.setMaximumSize(QSize(130, 30))
        self.ParaviewLogo.setPixmap(QPixmap(u"Icons/ParaView_logo.png"))
        self.ParaviewLogo.setScaledContents(True)

        self.horizontalLayout_5.addWidget(self.ParaviewLogo)


        self.verticalLayout_10.addWidget(self.AdditionalSoftwareLogos)


        self.verticalLayout_2.addWidget(self.AdditionalSoftware, 0, Qt.AlignTop)

        self.ToolSettings.addWidget(self.TitlePage)
        self.GeometryInformationTool = QWidget()
        self.GeometryInformationTool.setObjectName(u"GeometryInformationTool")
        self.verticalLayout_15 = QVBoxLayout(self.GeometryInformationTool)
        self.verticalLayout_15.setObjectName(u"verticalLayout_15")
        self.LGeometryInformation = QLabel(self.GeometryInformationTool)
        self.LGeometryInformation.setObjectName(u"LGeometryInformation")

        self.verticalLayout_15.addWidget(self.LGeometryInformation)

        self.BUploadGeometryFile = QPushButton(self.GeometryInformationTool)
        self.BUploadGeometryFile.setObjectName(u"BUploadGeometryFile")

        self.verticalLayout_15.addWidget(self.BUploadGeometryFile)

        self.STLVoxelization = QLabel(self.GeometryInformationTool)
        self.STLVoxelization.setObjectName(u"STLVoxelization")
        self.STLVoxelization.setEnabled(False)

        self.verticalLayout_15.addWidget(self.STLVoxelization)

        self.GeometryInputs = QFrame(self.GeometryInformationTool)
        self.GeometryInputs.setObjectName(u"GeometryInputs")
        self.GeometryInputs.setFrameShape(QFrame.NoFrame)
        self.GeometryInputs.setFrameShadow(QFrame.Raised)
        self.formLayout = QFormLayout(self.GeometryInputs)
        self.formLayout.setObjectName(u"formLayout")
        self.formLayout.setContentsMargins(-1, 0, -1, 0)
        self.LNumberOfVoxelsX = QLabel(self.GeometryInputs)
        self.LNumberOfVoxelsX.setObjectName(u"LNumberOfVoxelsX")
        self.LNumberOfVoxelsX.setEnabled(False)

        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.LNumberOfVoxelsX)

        self.INNumberOfVoxelsX = QLineEdit(self.GeometryInputs)
        self.INNumberOfVoxelsX.setObjectName(u"INNumberOfVoxelsX")
        self.INNumberOfVoxelsX.setEnabled(False)

        self.formLayout.setWidget(0, QFormLayout.FieldRole, self.INNumberOfVoxelsX)

        self.LNumberOfVoxelsY = QLabel(self.GeometryInputs)
        self.LNumberOfVoxelsY.setObjectName(u"LNumberOfVoxelsY")
        self.LNumberOfVoxelsY.setEnabled(False)

        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.LNumberOfVoxelsY)

        self.INNumberOfVoxelsY = QLineEdit(self.GeometryInputs)
        self.INNumberOfVoxelsY.setObjectName(u"INNumberOfVoxelsY")
        self.INNumberOfVoxelsY.setEnabled(False)

        self.formLayout.setWidget(1, QFormLayout.FieldRole, self.INNumberOfVoxelsY)

        self.LNumberOfVoxelsZ = QLabel(self.GeometryInputs)
        self.LNumberOfVoxelsZ.setObjectName(u"LNumberOfVoxelsZ")
        self.LNumberOfVoxelsZ.setEnabled(False)

        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.LNumberOfVoxelsZ)

        self.INNumberOfVoxelsZ = QLineEdit(self.GeometryInputs)
        self.INNumberOfVoxelsZ.setObjectName(u"INNumberOfVoxelsZ")
        self.INNumberOfVoxelsZ.setEnabled(False)

        self.formLayout.setWidget(2, QFormLayout.FieldRole, self.INNumberOfVoxelsZ)


        self.verticalLayout_15.addWidget(self.GeometryInputs)

        self.BVoxelizeGeometry = QPushButton(self.GeometryInformationTool)
        self.BVoxelizeGeometry.setObjectName(u"BVoxelizeGeometry")
        self.BVoxelizeGeometry.setEnabled(False)

        self.verticalLayout_15.addWidget(self.BVoxelizeGeometry)

        self.verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_15.addItem(self.verticalSpacer)

        self.ToolSettings.addWidget(self.GeometryInformationTool)
        self.DefineMaterialTool = QWidget()
        self.DefineMaterialTool.setObjectName(u"DefineMaterialTool")
        self.verticalLayout_16 = QVBoxLayout(self.DefineMaterialTool)
        self.verticalLayout_16.setObjectName(u"verticalLayout_16")
        self.LDefineMaterials = QLabel(self.DefineMaterialTool)
        self.LDefineMaterials.setObjectName(u"LDefineMaterials")

        self.verticalLayout_16.addWidget(self.LDefineMaterials)

        self.MaterialInputs = QFrame(self.DefineMaterialTool)
        self.MaterialInputs.setObjectName(u"MaterialInputs")
        self.MaterialInputs.setFrameShape(QFrame.NoFrame)
        self.MaterialInputs.setFrameShadow(QFrame.Raised)
        self.gridLayout = QGridLayout(self.MaterialInputs)
        self.gridLayout.setObjectName(u"gridLayout")
        self.gridLayout.setContentsMargins(-1, 0, -1, 0)
        self.LMaterialName = QLabel(self.MaterialInputs)
        self.LMaterialName.setObjectName(u"LMaterialName")

        self.gridLayout.addWidget(self.LMaterialName, 0, 0, 1, 1)

        self.INMaterialName = QLineEdit(self.MaterialInputs)
        self.INMaterialName.setObjectName(u"INMaterialName")
        self.INMaterialName.setMinimumSize(QSize(93, 0))

        self.gridLayout.addWidget(self.INMaterialName, 0, 1, 1, 1)

        self.INMaterialType = QComboBox(self.MaterialInputs)
        self.INMaterialType.addItem("")
        self.INMaterialType.addItem("")
        self.INMaterialType.addItem("")
        self.INMaterialType.addItem("")
        self.INMaterialType.setObjectName(u"INMaterialType")

        self.gridLayout.addWidget(self.INMaterialType, 0, 2, 1, 1)


        self.verticalLayout_16.addWidget(self.MaterialInputs)

        self.MaterialTypeTool = QStackedWidget(self.DefineMaterialTool)
        self.MaterialTypeTool.setObjectName(u"MaterialTypeTool")
        self.Isotropic = QWidget()
        self.Isotropic.setObjectName(u"Isotropic")
        self.formLayout_6 = QFormLayout(self.Isotropic)
        self.formLayout_6.setObjectName(u"formLayout_6")
        self.formLayout_6.setContentsMargins(-1, -1, -1, 12)
        self.LYoungsModulus = QLabel(self.Isotropic)
        self.LYoungsModulus.setObjectName(u"LYoungsModulus")

        self.formLayout_6.setWidget(0, QFormLayout.LabelRole, self.LYoungsModulus)

        self.INYoungsModulus = QLineEdit(self.Isotropic)
        self.INYoungsModulus.setObjectName(u"INYoungsModulus")

        self.formLayout_6.setWidget(0, QFormLayout.FieldRole, self.INYoungsModulus)

        self.LPoissonsRatio = QLabel(self.Isotropic)
        self.LPoissonsRatio.setObjectName(u"LPoissonsRatio")

        self.formLayout_6.setWidget(1, QFormLayout.LabelRole, self.LPoissonsRatio)

        self.INPoissonsRatio = QLineEdit(self.Isotropic)
        self.INPoissonsRatio.setObjectName(u"INPoissonsRatio")

        self.formLayout_6.setWidget(1, QFormLayout.FieldRole, self.INPoissonsRatio)

        self.MaterialTypeTool.addWidget(self.Isotropic)
        self.TransverselyIsotropic = QWidget()
        self.TransverselyIsotropic.setObjectName(u"TransverselyIsotropic")
        self.verticalLayout_23 = QVBoxLayout(self.TransverselyIsotropic)
        self.verticalLayout_23.setSpacing(0)
        self.verticalLayout_23.setObjectName(u"verticalLayout_23")
        self.verticalLayout_23.setContentsMargins(0, 0, 0, 0)
        self.IsotropicPlane = QFrame(self.TransverselyIsotropic)
        self.IsotropicPlane.setObjectName(u"IsotropicPlane")
        self.IsotropicPlane.setFrameShape(QFrame.NoFrame)
        self.IsotropicPlane.setFrameShadow(QFrame.Raised)
        self.formLayout_2 = QFormLayout(self.IsotropicPlane)
        self.formLayout_2.setObjectName(u"formLayout_2")
        self.formLayout_2.setContentsMargins(0, 0, 0, 0)
        self.LIsotropicPlane = QLabel(self.IsotropicPlane)
        self.LIsotropicPlane.setObjectName(u"LIsotropicPlane")

        self.formLayout_2.setWidget(0, QFormLayout.LabelRole, self.LIsotropicPlane)

        self.INIsotropicPlane = QComboBox(self.IsotropicPlane)
        self.INIsotropicPlane.addItem("")
        self.INIsotropicPlane.addItem("")
        self.INIsotropicPlane.addItem("")
        self.INIsotropicPlane.setObjectName(u"INIsotropicPlane")

        self.formLayout_2.setWidget(0, QFormLayout.FieldRole, self.INIsotropicPlane)


        self.verticalLayout_23.addWidget(self.IsotropicPlane)

        self.TransverslyIsotropicMat = QFrame(self.TransverselyIsotropic)
        self.TransverslyIsotropicMat.setObjectName(u"TransverslyIsotropicMat")
        self.TransverslyIsotropicMat.setFrameShape(QFrame.NoFrame)
        self.TransverslyIsotropicMat.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_6 = QHBoxLayout(self.TransverslyIsotropicMat)
        self.horizontalLayout_6.setSpacing(0)
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.horizontalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.TransverseInPlane = QFrame(self.TransverslyIsotropicMat)
        self.TransverseInPlane.setObjectName(u"TransverseInPlane")
        self.TransverseInPlane.setFrameShape(QFrame.NoFrame)
        self.TransverseInPlane.setFrameShadow(QFrame.Raised)
        self.verticalLayout_25 = QVBoxLayout(self.TransverseInPlane)
        self.verticalLayout_25.setSpacing(0)
        self.verticalLayout_25.setObjectName(u"verticalLayout_25")
        self.verticalLayout_25.setContentsMargins(0, 12, 0, 0)
        self.LInPlane = QLabel(self.TransverseInPlane)
        self.LInPlane.setObjectName(u"LInPlane")

        self.verticalLayout_25.addWidget(self.LInPlane)

        self.TransverseInPlaneMat = QFrame(self.TransverseInPlane)
        self.TransverseInPlaneMat.setObjectName(u"TransverseInPlaneMat")
        self.TransverseInPlaneMat.setFrameShape(QFrame.NoFrame)
        self.TransverseInPlaneMat.setFrameShadow(QFrame.Raised)
        self.formLayout_7 = QFormLayout(self.TransverseInPlaneMat)
        self.formLayout_7.setObjectName(u"formLayout_7")
        self.LEip = QLabel(self.TransverseInPlaneMat)
        self.LEip.setObjectName(u"LEip")
        self.LEip.setMinimumSize(QSize(31, 0))

        self.formLayout_7.setWidget(0, QFormLayout.LabelRole, self.LEip)

        self.INEip = QLineEdit(self.TransverseInPlaneMat)
        self.INEip.setObjectName(u"INEip")

        self.formLayout_7.setWidget(0, QFormLayout.FieldRole, self.INEip)

        self.LNUip = QLabel(self.TransverseInPlaneMat)
        self.LNUip.setObjectName(u"LNUip")
        self.LNUip.setMinimumSize(QSize(31, 0))

        self.formLayout_7.setWidget(1, QFormLayout.LabelRole, self.LNUip)

        self.INNUip = QLineEdit(self.TransverseInPlaneMat)
        self.INNUip.setObjectName(u"INNUip")

        self.formLayout_7.setWidget(1, QFormLayout.FieldRole, self.INNUip)


        self.verticalLayout_25.addWidget(self.TransverseInPlaneMat)


        self.horizontalLayout_6.addWidget(self.TransverseInPlane, 0, Qt.AlignTop)

        self.TransverseOutOfPlane = QFrame(self.TransverslyIsotropicMat)
        self.TransverseOutOfPlane.setObjectName(u"TransverseOutOfPlane")
        self.TransverseOutOfPlane.setFrameShape(QFrame.NoFrame)
        self.TransverseOutOfPlane.setFrameShadow(QFrame.Raised)
        self.verticalLayout_24 = QVBoxLayout(self.TransverseOutOfPlane)
        self.verticalLayout_24.setSpacing(0)
        self.verticalLayout_24.setObjectName(u"verticalLayout_24")
        self.verticalLayout_24.setContentsMargins(0, 12, 0, 0)
        self.LOutOfPlane = QLabel(self.TransverseOutOfPlane)
        self.LOutOfPlane.setObjectName(u"LOutOfPlane")

        self.verticalLayout_24.addWidget(self.LOutOfPlane)

        self.TransverseOutOfPlaneMat = QFrame(self.TransverseOutOfPlane)
        self.TransverseOutOfPlaneMat.setObjectName(u"TransverseOutOfPlaneMat")
        self.TransverseOutOfPlaneMat.setFrameShape(QFrame.NoFrame)
        self.TransverseOutOfPlaneMat.setFrameShadow(QFrame.Raised)
        self.formLayout_8 = QFormLayout(self.TransverseOutOfPlaneMat)
        self.formLayout_8.setObjectName(u"formLayout_8")
        self.LEop = QLabel(self.TransverseOutOfPlaneMat)
        self.LEop.setObjectName(u"LEop")
        self.LEop.setMinimumSize(QSize(31, 0))
        self.LEop.setMaximumSize(QSize(31, 16777215))

        self.formLayout_8.setWidget(0, QFormLayout.LabelRole, self.LEop)

        self.INEop = QLineEdit(self.TransverseOutOfPlaneMat)
        self.INEop.setObjectName(u"INEop")
        self.INEop.setMinimumSize(QSize(0, 0))

        self.formLayout_8.setWidget(0, QFormLayout.FieldRole, self.INEop)

        self.LNUop = QLabel(self.TransverseOutOfPlaneMat)
        self.LNUop.setObjectName(u"LNUop")
        self.LNUop.setMinimumSize(QSize(31, 0))
        self.LNUop.setMaximumSize(QSize(31, 16777215))

        self.formLayout_8.setWidget(1, QFormLayout.LabelRole, self.LNUop)

        self.INNUop = QLineEdit(self.TransverseOutOfPlaneMat)
        self.INNUop.setObjectName(u"INNUop")

        self.formLayout_8.setWidget(1, QFormLayout.FieldRole, self.INNUop)

        self.LGop = QLabel(self.TransverseOutOfPlaneMat)
        self.LGop.setObjectName(u"LGop")
        self.LGop.setMinimumSize(QSize(31, 0))
        self.LGop.setMaximumSize(QSize(31, 16777215))

        self.formLayout_8.setWidget(2, QFormLayout.LabelRole, self.LGop)

        self.INGop = QLineEdit(self.TransverseOutOfPlaneMat)
        self.INGop.setObjectName(u"INGop")

        self.formLayout_8.setWidget(2, QFormLayout.FieldRole, self.INGop)


        self.verticalLayout_24.addWidget(self.TransverseOutOfPlaneMat)


        self.horizontalLayout_6.addWidget(self.TransverseOutOfPlane)


        self.verticalLayout_23.addWidget(self.TransverslyIsotropicMat)

        self.MaterialTypeTool.addWidget(self.TransverselyIsotropic)
        self.Anisotropic = QWidget()
        self.Anisotropic.setObjectName(u"Anisotropic")
        self.verticalLayout_26 = QVBoxLayout(self.Anisotropic)
        self.verticalLayout_26.setObjectName(u"verticalLayout_26")
        self.verticalLayout_26.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel(self.Anisotropic)
        self.label.setObjectName(u"label")

        self.verticalLayout_26.addWidget(self.label)

        self.TAnisotropic = QTableWidget(self.Anisotropic)
        if (self.TAnisotropic.columnCount() < 6):
            self.TAnisotropic.setColumnCount(6)
        if (self.TAnisotropic.rowCount() < 6):
            self.TAnisotropic.setRowCount(6)
        __qtablewidgetitem = QTableWidgetItem()
        self.TAnisotropic.setItem(0, 0, __qtablewidgetitem)
        brush = QBrush(QColor(235, 235, 235, 255))
        brush.setStyle(Qt.SolidPattern)
        __qtablewidgetitem1 = QTableWidgetItem()
        __qtablewidgetitem1.setBackground(brush);
        __qtablewidgetitem1.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(1, 0, __qtablewidgetitem1)
        __qtablewidgetitem2 = QTableWidgetItem()
        __qtablewidgetitem2.setBackground(brush);
        __qtablewidgetitem2.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(2, 0, __qtablewidgetitem2)
        __qtablewidgetitem3 = QTableWidgetItem()
        __qtablewidgetitem3.setBackground(brush);
        __qtablewidgetitem3.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(2, 1, __qtablewidgetitem3)
        __qtablewidgetitem4 = QTableWidgetItem()
        __qtablewidgetitem4.setBackground(brush);
        __qtablewidgetitem4.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 0, __qtablewidgetitem4)
        __qtablewidgetitem5 = QTableWidgetItem()
        __qtablewidgetitem5.setBackground(brush);
        __qtablewidgetitem5.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 1, __qtablewidgetitem5)
        __qtablewidgetitem6 = QTableWidgetItem()
        __qtablewidgetitem6.setBackground(brush);
        __qtablewidgetitem6.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 2, __qtablewidgetitem6)
        __qtablewidgetitem7 = QTableWidgetItem()
        __qtablewidgetitem7.setBackground(brush);
        __qtablewidgetitem7.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 0, __qtablewidgetitem7)
        __qtablewidgetitem8 = QTableWidgetItem()
        __qtablewidgetitem8.setBackground(brush);
        __qtablewidgetitem8.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 1, __qtablewidgetitem8)
        __qtablewidgetitem9 = QTableWidgetItem()
        __qtablewidgetitem9.setBackground(brush);
        __qtablewidgetitem9.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 2, __qtablewidgetitem9)
        __qtablewidgetitem10 = QTableWidgetItem()
        __qtablewidgetitem10.setBackground(brush);
        __qtablewidgetitem10.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 3, __qtablewidgetitem10)
        __qtablewidgetitem11 = QTableWidgetItem()
        __qtablewidgetitem11.setBackground(brush);
        __qtablewidgetitem11.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 0, __qtablewidgetitem11)
        __qtablewidgetitem12 = QTableWidgetItem()
        __qtablewidgetitem12.setBackground(brush);
        __qtablewidgetitem12.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 1, __qtablewidgetitem12)
        __qtablewidgetitem13 = QTableWidgetItem()
        __qtablewidgetitem13.setBackground(brush);
        __qtablewidgetitem13.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 2, __qtablewidgetitem13)
        __qtablewidgetitem14 = QTableWidgetItem()
        __qtablewidgetitem14.setBackground(brush);
        __qtablewidgetitem14.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 3, __qtablewidgetitem14)
        __qtablewidgetitem15 = QTableWidgetItem()
        __qtablewidgetitem15.setBackground(brush);
        __qtablewidgetitem15.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 4, __qtablewidgetitem15)
        self.TAnisotropic.setObjectName(u"TAnisotropic")
        self.TAnisotropic.setRowCount(6)
        self.TAnisotropic.setColumnCount(6)
        self.TAnisotropic.horizontalHeader().setMinimumSectionSize(21)
        self.TAnisotropic.horizontalHeader().setDefaultSectionSize(50)

        self.verticalLayout_26.addWidget(self.TAnisotropic)

        self.MaterialTypeTool.addWidget(self.Anisotropic)
        self.Orthotropic = QWidget()
        self.Orthotropic.setObjectName(u"Orthotropic")
        self.gridLayout_2 = QGridLayout(self.Orthotropic)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.gridLayout_2.setHorizontalSpacing(-1)
        self.gridLayout_2.setVerticalSpacing(0)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.INNUyz = QLineEdit(self.Orthotropic)
        self.INNUyz.setObjectName(u"INNUyz")

        self.gridLayout_2.addWidget(self.INNUyz, 2, 8, 1, 1)

        self.LGxy = QLabel(self.Orthotropic)
        self.LGxy.setObjectName(u"LGxy")

        self.gridLayout_2.addWidget(self.LGxy, 3, 1, 1, 1, Qt.AlignRight)

        self.INGxy = QLineEdit(self.Orthotropic)
        self.INGxy.setObjectName(u"INGxy")

        self.gridLayout_2.addWidget(self.INGxy, 3, 2, 1, 1)

        self.INEy = QLineEdit(self.Orthotropic)
        self.INEy.setObjectName(u"INEy")

        self.gridLayout_2.addWidget(self.INEy, 0, 5, 1, 1)

        self.LEx = QLabel(self.Orthotropic)
        self.LEx.setObjectName(u"LEx")

        self.gridLayout_2.addWidget(self.LEx, 0, 1, 1, 1, Qt.AlignRight)

        self.LNUyz = QLabel(self.Orthotropic)
        self.LNUyz.setObjectName(u"LNUyz")

        self.gridLayout_2.addWidget(self.LNUyz, 2, 7, 1, 1, Qt.AlignRight)

        self.LGyz = QLabel(self.Orthotropic)
        self.LGyz.setObjectName(u"LGyz")

        self.gridLayout_2.addWidget(self.LGyz, 3, 7, 1, 1, Qt.AlignRight)

        self.LEy = QLabel(self.Orthotropic)
        self.LEy.setObjectName(u"LEy")

        self.gridLayout_2.addWidget(self.LEy, 0, 4, 1, 1, Qt.AlignRight)

        self.LEz = QLabel(self.Orthotropic)
        self.LEz.setObjectName(u"LEz")

        self.gridLayout_2.addWidget(self.LEz, 0, 7, 1, 1, Qt.AlignRight)

        self.LGxz = QLabel(self.Orthotropic)
        self.LGxz.setObjectName(u"LGxz")

        self.gridLayout_2.addWidget(self.LGxz, 3, 4, 1, 1, Qt.AlignRight)

        self.LNUxz = QLabel(self.Orthotropic)
        self.LNUxz.setObjectName(u"LNUxz")

        self.gridLayout_2.addWidget(self.LNUxz, 2, 4, 1, 1, Qt.AlignRight)

        self.INNUxy = QLineEdit(self.Orthotropic)
        self.INNUxy.setObjectName(u"INNUxy")

        self.gridLayout_2.addWidget(self.INNUxy, 2, 2, 1, 1)

        self.INEz = QLineEdit(self.Orthotropic)
        self.INEz.setObjectName(u"INEz")

        self.gridLayout_2.addWidget(self.INEz, 0, 8, 1, 1)

        self.LNUxy = QLabel(self.Orthotropic)
        self.LNUxy.setObjectName(u"LNUxy")

        self.gridLayout_2.addWidget(self.LNUxy, 2, 1, 1, 1, Qt.AlignRight)

        self.INGyz = QLineEdit(self.Orthotropic)
        self.INGyz.setObjectName(u"INGyz")

        self.gridLayout_2.addWidget(self.INGyz, 3, 8, 1, 1)

        self.INEx = QLineEdit(self.Orthotropic)
        self.INEx.setObjectName(u"INEx")

        self.gridLayout_2.addWidget(self.INEx, 0, 2, 1, 1)

        self.INGxz = QLineEdit(self.Orthotropic)
        self.INGxz.setObjectName(u"INGxz")

        self.gridLayout_2.addWidget(self.INGxz, 3, 5, 1, 1)

        self.INNUxz = QLineEdit(self.Orthotropic)
        self.INNUxz.setObjectName(u"INNUxz")

        self.gridLayout_2.addWidget(self.INNUxz, 2, 5, 1, 1)

        self.MaterialTypeTool.addWidget(self.Orthotropic)

        self.verticalLayout_16.addWidget(self.MaterialTypeTool)

        self.BAddMaterial = QPushButton(self.DefineMaterialTool)
        self.BAddMaterial.setObjectName(u"BAddMaterial")

        self.verticalLayout_16.addWidget(self.BAddMaterial)

        self.TMaterials = QTableWidget(self.DefineMaterialTool)
        if (self.TMaterials.columnCount() < 11):
            self.TMaterials.setColumnCount(11)
        __qtablewidgetitem16 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(0, __qtablewidgetitem16)
        __qtablewidgetitem17 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(1, __qtablewidgetitem17)
        __qtablewidgetitem18 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(2, __qtablewidgetitem18)
        __qtablewidgetitem19 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(3, __qtablewidgetitem19)
        __qtablewidgetitem20 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(4, __qtablewidgetitem20)
        __qtablewidgetitem21 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(5, __qtablewidgetitem21)
        __qtablewidgetitem22 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(6, __qtablewidgetitem22)
        __qtablewidgetitem23 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(7, __qtablewidgetitem23)
        __qtablewidgetitem24 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(8, __qtablewidgetitem24)
        __qtablewidgetitem25 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(9, __qtablewidgetitem25)
        __qtablewidgetitem26 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(10, __qtablewidgetitem26)
        self.TMaterials.setObjectName(u"TMaterials")
        self.TMaterials.setEnabled(True)
        self.TMaterials.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.TMaterials.setRowCount(0)

        self.verticalLayout_16.addWidget(self.TMaterials)

        self.BDeleteMaterial = QPushButton(self.DefineMaterialTool)
        self.BDeleteMaterial.setObjectName(u"BDeleteMaterial")

        self.verticalLayout_16.addWidget(self.BDeleteMaterial)

        self.line_6 = QFrame(self.DefineMaterialTool)
        self.line_6.setObjectName(u"line_6")
        self.line_6.setFrameShadow(QFrame.Sunken)
        self.line_6.setLineWidth(1)
        self.line_6.setFrameShape(QFrame.HLine)

        self.verticalLayout_16.addWidget(self.line_6)

        self.line_7 = QFrame(self.DefineMaterialTool)
        self.line_7.setObjectName(u"line_7")
        self.line_7.setFrameShape(QFrame.HLine)
        self.line_7.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_16.addWidget(self.line_7)

        self.label_2 = QLabel(self.DefineMaterialTool)
        self.label_2.setObjectName(u"label_2")

        self.verticalLayout_16.addWidget(self.label_2)

        self.TStiffnessMatrix = QTableWidget(self.DefineMaterialTool)
        if (self.TStiffnessMatrix.columnCount() < 6):
            self.TStiffnessMatrix.setColumnCount(6)
        if (self.TStiffnessMatrix.rowCount() < 6):
            self.TStiffnessMatrix.setRowCount(6)
        __qtablewidgetitem27 = QTableWidgetItem()
        __qtablewidgetitem27.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(0, 0, __qtablewidgetitem27)
        __qtablewidgetitem28 = QTableWidgetItem()
        __qtablewidgetitem28.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(0, 1, __qtablewidgetitem28)
        __qtablewidgetitem29 = QTableWidgetItem()
        __qtablewidgetitem29.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(0, 2, __qtablewidgetitem29)
        __qtablewidgetitem30 = QTableWidgetItem()
        __qtablewidgetitem30.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(0, 3, __qtablewidgetitem30)
        __qtablewidgetitem31 = QTableWidgetItem()
        __qtablewidgetitem31.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(0, 4, __qtablewidgetitem31)
        __qtablewidgetitem32 = QTableWidgetItem()
        __qtablewidgetitem32.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(0, 5, __qtablewidgetitem32)
        __qtablewidgetitem33 = QTableWidgetItem()
        __qtablewidgetitem33.setBackground(brush);
        __qtablewidgetitem33.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(1, 0, __qtablewidgetitem33)
        __qtablewidgetitem34 = QTableWidgetItem()
        __qtablewidgetitem34.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(1, 1, __qtablewidgetitem34)
        __qtablewidgetitem35 = QTableWidgetItem()
        __qtablewidgetitem35.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(1, 2, __qtablewidgetitem35)
        __qtablewidgetitem36 = QTableWidgetItem()
        __qtablewidgetitem36.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(1, 3, __qtablewidgetitem36)
        __qtablewidgetitem37 = QTableWidgetItem()
        __qtablewidgetitem37.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(1, 4, __qtablewidgetitem37)
        __qtablewidgetitem38 = QTableWidgetItem()
        __qtablewidgetitem38.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(1, 5, __qtablewidgetitem38)
        __qtablewidgetitem39 = QTableWidgetItem()
        __qtablewidgetitem39.setBackground(brush);
        __qtablewidgetitem39.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(2, 0, __qtablewidgetitem39)
        __qtablewidgetitem40 = QTableWidgetItem()
        __qtablewidgetitem40.setBackground(brush);
        __qtablewidgetitem40.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(2, 1, __qtablewidgetitem40)
        __qtablewidgetitem41 = QTableWidgetItem()
        __qtablewidgetitem41.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(2, 2, __qtablewidgetitem41)
        __qtablewidgetitem42 = QTableWidgetItem()
        __qtablewidgetitem42.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(2, 3, __qtablewidgetitem42)
        __qtablewidgetitem43 = QTableWidgetItem()
        __qtablewidgetitem43.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(2, 4, __qtablewidgetitem43)
        __qtablewidgetitem44 = QTableWidgetItem()
        __qtablewidgetitem44.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(2, 5, __qtablewidgetitem44)
        __qtablewidgetitem45 = QTableWidgetItem()
        __qtablewidgetitem45.setBackground(brush);
        __qtablewidgetitem45.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(3, 0, __qtablewidgetitem45)
        __qtablewidgetitem46 = QTableWidgetItem()
        __qtablewidgetitem46.setBackground(brush);
        __qtablewidgetitem46.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(3, 1, __qtablewidgetitem46)
        __qtablewidgetitem47 = QTableWidgetItem()
        __qtablewidgetitem47.setBackground(brush);
        __qtablewidgetitem47.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(3, 2, __qtablewidgetitem47)
        __qtablewidgetitem48 = QTableWidgetItem()
        __qtablewidgetitem48.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(3, 3, __qtablewidgetitem48)
        __qtablewidgetitem49 = QTableWidgetItem()
        __qtablewidgetitem49.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(3, 4, __qtablewidgetitem49)
        __qtablewidgetitem50 = QTableWidgetItem()
        __qtablewidgetitem50.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(3, 5, __qtablewidgetitem50)
        __qtablewidgetitem51 = QTableWidgetItem()
        __qtablewidgetitem51.setBackground(brush);
        __qtablewidgetitem51.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(4, 0, __qtablewidgetitem51)
        __qtablewidgetitem52 = QTableWidgetItem()
        __qtablewidgetitem52.setBackground(brush);
        __qtablewidgetitem52.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(4, 1, __qtablewidgetitem52)
        __qtablewidgetitem53 = QTableWidgetItem()
        __qtablewidgetitem53.setBackground(brush);
        __qtablewidgetitem53.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(4, 2, __qtablewidgetitem53)
        __qtablewidgetitem54 = QTableWidgetItem()
        __qtablewidgetitem54.setBackground(brush);
        __qtablewidgetitem54.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(4, 3, __qtablewidgetitem54)
        __qtablewidgetitem55 = QTableWidgetItem()
        __qtablewidgetitem55.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(4, 4, __qtablewidgetitem55)
        __qtablewidgetitem56 = QTableWidgetItem()
        __qtablewidgetitem56.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(4, 5, __qtablewidgetitem56)
        __qtablewidgetitem57 = QTableWidgetItem()
        __qtablewidgetitem57.setBackground(brush);
        __qtablewidgetitem57.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(5, 0, __qtablewidgetitem57)
        __qtablewidgetitem58 = QTableWidgetItem()
        __qtablewidgetitem58.setBackground(brush);
        __qtablewidgetitem58.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(5, 1, __qtablewidgetitem58)
        __qtablewidgetitem59 = QTableWidgetItem()
        __qtablewidgetitem59.setBackground(brush);
        __qtablewidgetitem59.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(5, 2, __qtablewidgetitem59)
        __qtablewidgetitem60 = QTableWidgetItem()
        __qtablewidgetitem60.setBackground(brush);
        __qtablewidgetitem60.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(5, 3, __qtablewidgetitem60)
        __qtablewidgetitem61 = QTableWidgetItem()
        __qtablewidgetitem61.setBackground(brush);
        __qtablewidgetitem61.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TStiffnessMatrix.setItem(5, 4, __qtablewidgetitem61)
        __qtablewidgetitem62 = QTableWidgetItem()
        __qtablewidgetitem62.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        self.TStiffnessMatrix.setItem(5, 5, __qtablewidgetitem62)
        self.TStiffnessMatrix.setObjectName(u"TStiffnessMatrix")
        self.TStiffnessMatrix.setEnabled(True)
        self.TStiffnessMatrix.setRowCount(6)
        self.TStiffnessMatrix.setColumnCount(6)
        self.TStiffnessMatrix.horizontalHeader().setMinimumSectionSize(21)
        self.TStiffnessMatrix.horizontalHeader().setDefaultSectionSize(50)

        self.verticalLayout_16.addWidget(self.TStiffnessMatrix)

        self.verticalSpacer_2 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_16.addItem(self.verticalSpacer_2)

        self.ToolSettings.addWidget(self.DefineMaterialTool)
        self.BoundaryConditionsTool = QWidget()
        self.BoundaryConditionsTool.setObjectName(u"BoundaryConditionsTool")
        self.verticalLayout_17 = QVBoxLayout(self.BoundaryConditionsTool)
        self.verticalLayout_17.setObjectName(u"verticalLayout_17")
        self.LBoundaryConditions = QLabel(self.BoundaryConditionsTool)
        self.LBoundaryConditions.setObjectName(u"LBoundaryConditions")

        self.verticalLayout_17.addWidget(self.LBoundaryConditions)

        self.BoundaryConditionsInputs = QFrame(self.BoundaryConditionsTool)
        self.BoundaryConditionsInputs.setObjectName(u"BoundaryConditionsInputs")
        self.BoundaryConditionsInputs.setFrameShape(QFrame.NoFrame)
        self.BoundaryConditionsInputs.setFrameShadow(QFrame.Raised)
        self.formLayout_3 = QFormLayout(self.BoundaryConditionsInputs)
        self.formLayout_3.setObjectName(u"formLayout_3")
        self.LBoundaryCondition = QLabel(self.BoundaryConditionsInputs)
        self.LBoundaryCondition.setObjectName(u"LBoundaryCondition")

        self.formLayout_3.setWidget(0, QFormLayout.LabelRole, self.LBoundaryCondition)

        self.INBoundaryCondition = QComboBox(self.BoundaryConditionsInputs)
        self.INBoundaryCondition.addItem("")
        self.INBoundaryCondition.addItem("")
        self.INBoundaryCondition.setObjectName(u"INBoundaryCondition")

        self.formLayout_3.setWidget(0, QFormLayout.FieldRole, self.INBoundaryCondition)

        self.LBCDirection = QLabel(self.BoundaryConditionsInputs)
        self.LBCDirection.setObjectName(u"LBCDirection")

        self.formLayout_3.setWidget(1, QFormLayout.LabelRole, self.LBCDirection)

        self.INBCDirection = QComboBox(self.BoundaryConditionsInputs)
        self.INBCDirection.addItem("")
        self.INBCDirection.addItem("")
        self.INBCDirection.addItem("")
        self.INBCDirection.setObjectName(u"INBCDirection")

        self.formLayout_3.setWidget(1, QFormLayout.FieldRole, self.INBCDirection)


        self.verticalLayout_17.addWidget(self.BoundaryConditionsInputs)

        self.BAddBC = QPushButton(self.BoundaryConditionsTool)
        self.BAddBC.setObjectName(u"BAddBC")

        self.verticalLayout_17.addWidget(self.BAddBC)

        self.TBCs = QTableWidget(self.BoundaryConditionsTool)
        if (self.TBCs.columnCount() < 2):
            self.TBCs.setColumnCount(2)
        __qtablewidgetitem63 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(0, __qtablewidgetitem63)
        __qtablewidgetitem64 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(1, __qtablewidgetitem64)
        self.TBCs.setObjectName(u"TBCs")
        self.TBCs.setEditTriggers(QAbstractItemView.NoEditTriggers)

        self.verticalLayout_17.addWidget(self.TBCs)

        self.BDeleteBC = QPushButton(self.BoundaryConditionsTool)
        self.BDeleteBC.setObjectName(u"BDeleteBC")

        self.verticalLayout_17.addWidget(self.BDeleteBC)

        self.verticalSpacer_3 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_17.addItem(self.verticalSpacer_3)

        self.ToolSettings.addWidget(self.BoundaryConditionsTool)
        self.SolverSettingsTool = QWidget()
        self.SolverSettingsTool.setObjectName(u"SolverSettingsTool")
        self.verticalLayout_9 = QVBoxLayout(self.SolverSettingsTool)
        self.verticalLayout_9.setObjectName(u"verticalLayout_9")
        self.LSolverSettings_2 = QLabel(self.SolverSettingsTool)
        self.LSolverSettings_2.setObjectName(u"LSolverSettings_2")

        self.verticalLayout_9.addWidget(self.LSolverSettings_2)

        self.SolverSettingsInputs = QFrame(self.SolverSettingsTool)
        self.SolverSettingsInputs.setObjectName(u"SolverSettingsInputs")
        self.SolverSettingsInputs.setFrameShape(QFrame.NoFrame)
        self.SolverSettingsInputs.setFrameShadow(QFrame.Raised)
        self.formLayout_4 = QFormLayout(self.SolverSettingsInputs)
        self.formLayout_4.setObjectName(u"formLayout_4")
        self.LNumberOfSteps = QLabel(self.SolverSettingsInputs)
        self.LNumberOfSteps.setObjectName(u"LNumberOfSteps")

        self.formLayout_4.setWidget(0, QFormLayout.LabelRole, self.LNumberOfSteps)

        self.INNumberOfSteps = QLineEdit(self.SolverSettingsInputs)
        self.INNumberOfSteps.setObjectName(u"INNumberOfSteps")
        self.INNumberOfSteps.setReadOnly(False)

        self.formLayout_4.setWidget(0, QFormLayout.FieldRole, self.INNumberOfSteps)


        self.verticalLayout_9.addWidget(self.SolverSettingsInputs)

        self.verticalSpacer_4 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_9.addItem(self.verticalSpacer_4)

        self.ToolSettings.addWidget(self.SolverSettingsTool)
        self.ResultsTool = QWidget()
        self.ResultsTool.setObjectName(u"ResultsTool")
        self.verticalLayout_18 = QVBoxLayout(self.ResultsTool)
        self.verticalLayout_18.setObjectName(u"verticalLayout_18")
        self.LResults = QLabel(self.ResultsTool)
        self.LResults.setObjectName(u"LResults")

        self.verticalLayout_18.addWidget(self.LResults)

        self.PreviewResults = QFrame(self.ResultsTool)
        self.PreviewResults.setObjectName(u"PreviewResults")
        self.PreviewResults.setFrameShape(QFrame.NoFrame)
        self.PreviewResults.setFrameShadow(QFrame.Raised)
        self.formLayout_5 = QFormLayout(self.PreviewResults)
        self.formLayout_5.setObjectName(u"formLayout_5")
        self.formLayout_5.setContentsMargins(-1, 12, -1, 12)
        self.INPreviewResults = QComboBox(self.PreviewResults)
        self.INPreviewResults.addItem("")
        self.INPreviewResults.addItem("")
        self.INPreviewResults.setObjectName(u"INPreviewResults")
        self.INPreviewResults.setFrame(True)

        self.formLayout_5.setWidget(0, QFormLayout.LabelRole, self.INPreviewResults)

        self.BPreviewResults = QPushButton(self.PreviewResults)
        self.BPreviewResults.setObjectName(u"BPreviewResults")
        self.BPreviewResults.setMinimumSize(QSize(140, 0))

        self.formLayout_5.setWidget(0, QFormLayout.FieldRole, self.BPreviewResults)

        self.BPlotSS = QPushButton(self.PreviewResults)
        self.BPlotSS.setObjectName(u"BPlotSS")

        self.formLayout_5.setWidget(1, QFormLayout.FieldRole, self.BPlotSS)

        self.INPlotSS = QComboBox(self.PreviewResults)
        self.INPlotSS.addItem("")
        self.INPlotSS.addItem("")
        self.INPlotSS.addItem("")
        self.INPlotSS.setObjectName(u"INPlotSS")

        self.formLayout_5.setWidget(1, QFormLayout.LabelRole, self.INPlotSS)


        self.verticalLayout_18.addWidget(self.PreviewResults)

        self.BOpenParaview = QPushButton(self.ResultsTool)
        self.BOpenParaview.setObjectName(u"BOpenParaview")

        self.verticalLayout_18.addWidget(self.BOpenParaview)

        self.verticalSpacer_5 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_18.addItem(self.verticalSpacer_5)

        self.ToolSettings.addWidget(self.ResultsTool)
        self.splitter.addWidget(self.ToolSettings)
        self.ParaviewFrame = QFrame(self.splitter)
        self.ParaviewFrame.setObjectName(u"ParaviewFrame")
        self.ParaviewFrame.setFrameShape(QFrame.Box)
        self.ParaviewFrame.setFrameShadow(QFrame.Plain)
        self.ParaviewFrame.setLineWidth(1)
        self.verticalLayout_19 = QVBoxLayout(self.ParaviewFrame)
        self.verticalLayout_19.setSpacing(0)
        self.verticalLayout_19.setObjectName(u"verticalLayout_19")
        self.verticalLayout_19.setContentsMargins(0, 0, 0, 0)
        self.splitter_2 = QSplitter(self.ParaviewFrame)
        self.splitter_2.setObjectName(u"splitter_2")
        self.splitter_2.setOrientation(Qt.Vertical)
        self.OutputWindows = QStackedWidget(self.splitter_2)
        self.OutputWindows.setObjectName(u"OutputWindows")
        self.OutputWindows.setFrameShape(QFrame.NoFrame)
        self.OutputWindows.setLineWidth(6)
        self.OutputWindows.setMidLineWidth(0)
        self.ParaviewWindow = QWidget()
        self.ParaviewWindow.setObjectName(u"ParaviewWindow")
        self.verticalLayout_20 = QVBoxLayout(self.ParaviewWindow)
        self.verticalLayout_20.setObjectName(u"verticalLayout_20")
        self.verticalLayout_20.setContentsMargins(0, 0, 0, 0)
        self.Paraview = QWidget(self.ParaviewWindow)
        self.Paraview.setObjectName(u"Paraview")
        
        # ================== New imports =====================
        # Paraview imports
        import paraview.simple as pvsimple
        from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
        render_view = pvsimple.CreateRenderView()
        self.Paraview = QVTKRenderWindowInteractor(rw=render_view.GetRenderWindow(),iren=render_view.GetInteractor())
        
        # Open Url in Help Menu
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices
        def openUrl():
            url = QUrl('https://lanl.github.io/Fierro/')
            if not QDesktopServices.openUrl(url):
                QMessageBox.warning(self, 'Open Url', 'Could not open url')
        # ======================================================

        self.verticalLayout_20.addWidget(self.Paraview)

        self.OutputWindows.addWidget(self.ParaviewWindow)
        self.PlotWindow = QWidget()
        self.PlotWindow.setObjectName(u"PlotWindow")
        self.verticalLayout_21 = QVBoxLayout(self.PlotWindow)
        self.verticalLayout_21.setObjectName(u"verticalLayout_21")
        self.verticalLayout_21.setContentsMargins(0, 0, 0, 0)
        self.Plot = QWidget(self.PlotWindow)
        self.Plot.setObjectName(u"Plot")

        self.verticalLayout_21.addWidget(self.Plot)

        self.OutputWindows.addWidget(self.PlotWindow)
        self.splitter_2.addWidget(self.OutputWindows)
        self.RunOutputs = QFrame(self.splitter_2)
        self.RunOutputs.setObjectName(u"RunOutputs")
        self.RunOutputs.setMaximumSize(QSize(16777215, 175))
        self.RunOutputs.setFrameShape(QFrame.NoFrame)
        self.RunOutputs.setFrameShadow(QFrame.Raised)
        self.verticalLayout_22 = QVBoxLayout(self.RunOutputs)
        self.verticalLayout_22.setSpacing(0)
        self.verticalLayout_22.setObjectName(u"verticalLayout_22")
        self.verticalLayout_22.setContentsMargins(0, 0, 0, 0)
        self.RunOutputProgress = QProgressBar(self.RunOutputs)
        self.RunOutputProgress.setObjectName(u"RunOutputProgress")
        self.RunOutputProgress.setMinimumSize(QSize(500, 0))
        self.RunOutputProgress.setMaximumSize(QSize(16777215, 16777215))
        self.RunOutputProgress.setValue(0)

        self.verticalLayout_22.addWidget(self.RunOutputProgress, 0, Qt.AlignHCenter)

        self.RunOutputWindow = QPlainTextEdit(self.RunOutputs)
        self.RunOutputWindow.setObjectName(u"RunOutputWindow")
        self.RunOutputWindow.setFrameShape(QFrame.NoFrame)
        self.RunOutputWindow.setFrameShadow(QFrame.Plain)
        self.RunOutputWindow.setReadOnly(True)

        self.verticalLayout_22.addWidget(self.RunOutputWindow)

        self.splitter_2.addWidget(self.RunOutputs)

        self.verticalLayout_19.addWidget(self.splitter_2)

        self.splitter.addWidget(self.ParaviewFrame)

        self.horizontalLayout.addWidget(self.splitter)


        self.verticalLayout.addWidget(self.Main)

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1022, 24))
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setObjectName(u"menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuHelp.menuAction())
        self.menuHelp.addAction(self.actionEVPFFT_Manual)

        self.retranslateUi(MainWindow)

        self.HeaderMenu.setCurrentIndex(0)
        self.ToolSettings.setCurrentIndex(3)
        self.MaterialTypeTool.setCurrentIndex(0)
        self.OutputWindows.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
        
        # ===================== New imports =======================
        import os, sys
        import numpy as np
        from PySide6.QtWidgets import (QFileDialog, QMessageBox)
        from PySide6.QtCore import (QTimer, QProcess)
        import matplotlib
        matplotlib.use('qt5agg')
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        import csv
        import re
        import shutil
        
        # BUTTON SETUP
        # Connect tab buttons to settings windows
        self.BImportPart.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(1))
        self.BDefineMaterial.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(2))
        self.BApplyBC.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(3))
        self.BSolverSettings.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(4))
        self.BViewResults.clicked.connect(lambda: self.ToolSettings.setCurrentIndex(5))
        
        # Help menu
        self.actionEVPFFT_Manual.triggered.connect(openUrl)
        
        # Upload Geometry
        def geometry_upload_click():
            try:
                self.stl
            except:
                print('')
            else:
                pvsimple.Delete(self.stl)
                
            try:
                self.voxel_reader
            except:
                print('')
            else:
                pvsimple.Delete(self.threshold)
                
            global b3_filename
            b3_filename = QFileDialog.getOpenFileName(
                filter="Geometry File (*.stl *.vtk)",
            )
            # Paraview window
            self.file_type = b3_filename[0][-4:-1]
            if self.file_type == '.st':
                self.stl = pvsimple.STLReader(FileNames = b3_filename)
                self.STLVoxelization.setEnabled(True)
                self.LNumberOfVoxelsX.setEnabled(True)
                self.INNumberOfVoxelsX.setEnabled(True)
                self.LNumberOfVoxelsY.setEnabled(True)
                self.INNumberOfVoxelsY.setEnabled(True)
                self.LNumberOfVoxelsZ.setEnabled(True)
                self.INNumberOfVoxelsZ.setEnabled(True)
                self.BVoxelizeGeometry.setEnabled(True)
            elif self.file_type == '.vt':
                self.stl = pvsimple.LegacyVTKReader(FileNames = b3_filename)
                pvsimple.SetDisplayProperties(Representation = "Surface")
                shutil.copy(b3_filename[0],os.getcwd()+'/VTK_geometry.vtk')
                self.INNumberOfVoxelsX.setText(QCoreApplication.translate("MainWindow", u"32", None))
                self.INNumberOfVoxelsY.setText(QCoreApplication.translate("MainWindow", u"32", None))
                self.INNumberOfVoxelsZ.setText(QCoreApplication.translate("MainWindow", u"32", None))
                self.STLVoxelization.setEnabled(False)
                self.LNumberOfVoxelsX.setEnabled(False)
                self.INNumberOfVoxelsX.setEnabled(False)
                self.LNumberOfVoxelsY.setEnabled(False)
                self.INNumberOfVoxelsY.setEnabled(False)
                self.LNumberOfVoxelsZ.setEnabled(False)
                self.INNumberOfVoxelsZ.setEnabled(False)
                self.BVoxelizeGeometry.setEnabled(False)
            else:
                warning_message('ERROR: Incorrect file type')
            pvsimple.Show(self.stl, render_view)
            pvsimple.ResetCamera(view=None)
        self.BUploadGeometryFile.clicked.connect(geometry_upload_click)
        
        # Voxelize Geometry
        def voxelize_geometry_click():
            if not self.INNumberOfVoxelsX.text() or not self.INNumberOfVoxelsY.text() or not self.INNumberOfVoxelsZ.text():
                warning_message('ERROR: Number of voxels NOT defined')
            else:
                try:
                    self.voxel_reader
                except:
                    print('')
                else:
                    pvsimple.Delete(self.threshold)
                    
                os.system("./voxelvtk "+b3_filename[0]+" VTK_geometry.vtk "+self.INNumberOfVoxelsX.text()+" "+self.INNumberOfVoxelsY.text()+" "+self.INNumberOfVoxelsZ.text())
                # Paraview window
                pvsimple.Delete(self.stl)
                self.voxel_reader = pvsimple.LegacyVTKReader(FileNames = "VTK_geometry.vtk")
                pvsimple.SetDisplayProperties(Representation = "Surface")
                self.threshold = pvsimple.Threshold(Input = self.voxel_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
                pvsimple.Show(self.threshold, render_view)
                pvsimple.Hide(self.voxel_reader)
                render_view.ResetCamera()
                render_view.StillRender()
        self.BVoxelizeGeometry.clicked.connect(voxelize_geometry_click)
        
        # Apply Material
        def material_type():
            if str(self.INMaterialType.currentText()) == 'Isotropic':
                self.MaterialTypeTool.setCurrentIndex(0)
            if str(self.INMaterialType.currentText()) == 'Transversely Isotropic':
                self.MaterialTypeTool.setCurrentIndex(1)
            if str(self.INMaterialType.currentText()) == 'Orthotropic':
                self.MaterialTypeTool.setCurrentIndex(3)
            if str(self.INMaterialType.currentText()) == 'Anisotropic':
                self.MaterialTypeTool.setCurrentIndex(2)
        self.INMaterialType.currentIndexChanged.connect(material_type)
        
        def add_material():
            if str(self.INMaterialType.currentText()) == 'Isotropic':
                if not self.INYoungsModulus.text() or not self.INPoissonsRatio.text() or not self.INMaterialName.text():
                    warning_message('ERROR: Material definition incomplete')
                else:
                    # Fill out material definition table
                    row = self.TMaterials.rowCount()
                    self.TMaterials.insertRow(row)
                    self.TMaterials.setItem(row, 0, QTableWidgetItem(
                        self.INMaterialName.text().strip())
                    )
                    self.TMaterials.setItem(row, 1, QTableWidgetItem(
                        str(self.INMaterialType.currentText()))
                    )
                    self.TMaterials.setItem(
                        row, 2, QTableWidgetItem(self.INYoungsModulus.text())
                    )
                    self.TMaterials.setItem(
                        row, 3, QTableWidgetItem(self.INYoungsModulus.text())
                    )
                    self.TMaterials.setItem(
                        row, 4, QTableWidgetItem(self.INYoungsModulus.text())
                    )
                    self.TMaterials.setItem(
                        row, 5, QTableWidgetItem(self.INPoissonsRatio.text())
                    )
                    self.TMaterials.setItem(
                        row, 6, QTableWidgetItem(self.INPoissonsRatio.text())
                    )
                    self.TMaterials.setItem(
                        row, 7, QTableWidgetItem(self.INPoissonsRatio.text())
                    )
                    self.INCalcG = float(self.INYoungsModulus.text())/(2*(1+float(self.INPoissonsRatio.text())))
                    self.TMaterials.setItem(
                        row, 8, QTableWidgetItem(str(self.INCalcG))
                    )
                    self.TMaterials.setItem(
                        row, 9, QTableWidgetItem(str(self.INCalcG))
                    )
                    self.TMaterials.setItem(
                        row, 10, QTableWidgetItem(str(self.INCalcG))
                    )
                    
                    # Calculate Stiffness Matrix
                    E = float(self.INYoungsModulus.text())
                    NU = float(self.INPoissonsRatio.text())
                    S11 = 1/E
                    S12 = -NU/E
                    S13 = -NU/E
                    S22 = 1/E
                    S23 = -NU/E
                    S33 = 1/E
                    S44 = 1/self.INCalcG
                    S55 = 1/self.INCalcG
                    S66 = 1/self.INCalcG
                    S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                    self.C11 = (1/S)*(S22*S33-S23*S23)
                    self.C12 = (1/S)*(S13*S23-S12*S33)
                    self.C22 = (1/S)*(S33*S11-S13*S13)
                    self.C13 = (1/S)*(S12*S23-S13*S22)
                    self.C33 = (1/S)*(S11*S22-S12*S12)
                    self.C23 = (1/S)*(S12*S13-S23*S11)
                    self.C44 = (1/S44)
                    self.C55 = (1/S55)
                    self.C66 = (1/S66)
                    
                    # Fill out stiffness matrix table
                    self.TStiffnessMatrix.setItem(0, 0, QTableWidgetItem(str(self.C11)))
                    self.TStiffnessMatrix.setItem(0, 1, QTableWidgetItem(str(self.C12)))
                    self.TStiffnessMatrix.setItem(0, 2, QTableWidgetItem(str(self.C13)))
                    self.TStiffnessMatrix.setItem(0, 3, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(0, 4, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(0, 5, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(1, 1, QTableWidgetItem(str(self.C22)))
                    self.TStiffnessMatrix.setItem(1, 2, QTableWidgetItem(str(self.C23)))
                    self.TStiffnessMatrix.setItem(1, 3, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(1, 4, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(1, 5, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(2, 2, QTableWidgetItem(str(self.C33)))
                    self.TStiffnessMatrix.setItem(2, 3, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(2, 4, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(2, 5, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(3, 3, QTableWidgetItem(str(self.C44)))
                    self.TStiffnessMatrix.setItem(3, 4, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(3, 5, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(4, 4, QTableWidgetItem(str(self.C55)))
                    self.TStiffnessMatrix.setItem(4, 5, QTableWidgetItem('0'))
                    self.TStiffnessMatrix.setItem(5, 5, QTableWidgetItem(str(self.C66)))
                    
                    self.INMaterialName.clear()
                    self.INYoungsModulus.clear()
                    self.INPoissonsRatio.clear()
            if str(self.INMaterialType.currentText()) == 'Transversely Isotropic':
                if not self.INEip.text() or not self.INNUip.text() or not self.INEop.text() or not self.INNUop.text() or not self.INGop.text() or not self.INMaterialName.text():
                    warning_message('ERROR: Material definition incomplete')
                else:
                    row = self.TMaterials.rowCount()
                    self.TMaterials.insertRow(row)
                    if str(self.INIsotropicPlane.currentText()) == 'x-y plane':
                        self.TMaterials.setItem(row, 0, QTableWidgetItem(
                            self.INMaterialName.text().strip())
                        )
                        self.TMaterials.setItem(row, 1, QTableWidgetItem(
                            str(self.INMaterialType.currentText()) + ': ' + self.INIsotropicPlane.currentText())
                        )
                        self.TMaterials.setItem(
                            row, 2, QTableWidgetItem(self.INEip.text())
                        )
                        self.TMaterials.setItem(
                            row, 3, QTableWidgetItem(self.INEip.text())
                        )
                        self.TMaterials.setItem(
                            row, 4, QTableWidgetItem(self.INEop.text())
                        )
                        self.TMaterials.setItem(
                            row, 5, QTableWidgetItem(self.INNUip.text())
                        )
                        self.TMaterials.setItem(
                            row, 6, QTableWidgetItem(self.INNUop.text())
                        )
                        self.TMaterials.setItem(
                            row, 7, QTableWidgetItem(self.INNUop.text())
                        )
                        self.INCalcG = float(self.INEip.text())/(2*(1+float(self.INNUip.text())))
                        self.TMaterials.setItem(
                            row, 8, QTableWidgetItem(str(self.INCalcG))
                        )
                        self.TMaterials.setItem(
                            row, 9, QTableWidgetItem(self.INGop.text())
                        )
                        self.TMaterials.setItem(
                            row, 10, QTableWidgetItem(self.INGop.text())
                        )
                        self.INEip.clear()
                        self.INNUip.clear()
                        self.INEop.clear()
                        self.INNUop.clear()
                        self.INGop.clear()
                        self.INMaterialName.clear()
                    if str(self.INIsotropicPlane.currentText()) == 'x-z plane':
                        self.TMaterials.setItem(row, 0, QTableWidgetItem(
                            self.INMaterialName.text().strip())
                        )
                        self.TMaterials.setItem(row, 1, QTableWidgetItem(
                            str(self.INMaterialType.currentText()) + ': ' + self.INIsotropicPlane.currentText())
                        )
                        self.TMaterials.setItem(
                            row, 2, QTableWidgetItem(self.INEip.text())
                        )
                        self.TMaterials.setItem(
                            row, 3, QTableWidgetItem(self.INEop.text())
                        )
                        self.TMaterials.setItem(
                            row, 4, QTableWidgetItem(self.INEip.text())
                        )
                        self.TMaterials.setItem(
                            row, 5, QTableWidgetItem(self.INNUop.text())
                        )
                        self.TMaterials.setItem(
                            row, 6, QTableWidgetItem(self.INNUip.text())
                        )
                        self.TMaterials.setItem(
                            row, 7, QTableWidgetItem(self.INNUop.text())
                        )
                        self.INCalcG = float(self.INEip.text())/(2*(1+float(self.INNUip.text())))
                        self.TMaterials.setItem(
                            row, 8, QTableWidgetItem(self.INGop.text())
                        )
                        self.TMaterials.setItem(
                            row, 9, QTableWidgetItem(str(self.INCalcG))
                        )
                        self.TMaterials.setItem(
                            row, 10, QTableWidgetItem(self.INGop.text())
                        )
                        self.INEip.clear()
                        self.INNUip.clear()
                        self.INEop.clear()
                        self.INNUop.clear()
                        self.INGop.clear()
                        self.INMaterialName.clear()
                    if str(self.INIsotropicPlane.currentText()) == 'y-z plane':
                        self.TMaterials.setItem(row, 0, QTableWidgetItem(
                            self.INMaterialName.text().strip())
                        )
                        self.TMaterials.setItem(row, 1, QTableWidgetItem(
                            str(self.INMaterialType.currentText()) + ': ' + self.INIsotropicPlane.currentText())
                        )
                        self.TMaterials.setItem(
                            row, 2, QTableWidgetItem(self.INEop.text())
                        )
                        self.TMaterials.setItem(
                            row, 3, QTableWidgetItem(self.INEip.text())
                        )
                        self.TMaterials.setItem(
                            row, 4, QTableWidgetItem(self.INEip.text())
                        )
                        self.TMaterials.setItem(
                            row, 5, QTableWidgetItem(self.INNUop.text())
                        )
                        self.TMaterials.setItem(
                            row, 6, QTableWidgetItem(self.INNUop.text())
                        )
                        self.TMaterials.setItem(
                            row, 7, QTableWidgetItem(self.INNUip.text())
                        )
                        self.INCalcG = float(self.INEip.text())/(2*(1+float(self.INNUip.text())))
                        self.TMaterials.setItem(
                            row, 8, QTableWidgetItem(self.INGop.text())
                        )
                        self.TMaterials.setItem(
                            row, 9, QTableWidgetItem(self.INGop.text())
                        )
                        self.TMaterials.setItem(
                            row, 10, QTableWidgetItem(str(self.INCalcG))
                        )
                        self.INEip.clear()
                        self.INNUip.clear()
                        self.INEop.clear()
                        self.INNUop.clear()
                        self.INGop.clear()
                        self.INMaterialName.clear()
            if str(self.INMaterialType.currentText()) == 'Orthotropic':
                if not self.INEx.text() or not self.INEy.text() or not self.INEz.text() or not self.INNUxy.text() or not self.INNUxz.text() or not self.INNUyz.text() or not self.INGxy.text() or not self.INGxz.text() or not self.INGyz.text() or not self.INMaterialName.text():
                    warning_message('ERROR: Material definition incomplete')
                else:
                    row = self.TMaterials.rowCount()
                    self.TMaterials.insertRow(row)
                    self.TMaterials.setItem(row, 0, QTableWidgetItem(
                        self.INMaterialName.text().strip())
                    )
                    self.TMaterials.setItem(row, 1, QTableWidgetItem(
                        str(self.INMaterialType.currentText()))
                    )
                    self.INMaterialName.clear()
                
        def delete_material():
            current_row = self.TMaterials.currentRow()
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
                self.TMaterials.removeRow(current_row)
        
        self.BAddMaterial.clicked.connect(add_material)
        self.BDeleteMaterial.clicked.connect(delete_material)
        
        # Boundary Conditions
        def add_bcs():
#            if not self.INBCStrainRate.text():
#                on_clicked("ERROR: Boundary condition incomplete")
#            else:
            row = self.TBCs.rowCount()
            self.TBCs.insertRow(row)
            self.TBCs.setItem(row, 0, QTableWidgetItem(str(
                self.INBoundaryCondition.currentText()))
            )
            self.TBCs.setItem(
                row, 1, QTableWidgetItem(str(self.INBCDirection.currentText()))
            )
#            self.TBCs.setItem(
#                row, 2, QTableWidgetItem(self.INBCStrainRate.text())
#            )
#            reset_bcs()
#
#        def reset_bcs():
#            self.INBCStrainRate.clear()

        def delete_bcs():
            current_row = self.TBCs.currentRow()
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
                self.TBCs.removeRow(current_row)
                
        self.BAddBC.clicked.connect(add_bcs)
        self.BDeleteBC.clicked.connect(delete_bcs)
        
        # Write input file
        def write_input_file():
            # Get current file location
            current_file_loc = os.getcwd()
            
            # Elastic parameters file
            if self.TMaterials.item(0,1).text() == 'Isotropic':
                elastic_parameters = open("../example_input_files/lattice_input_files/elastic_parameters.txt","w")
                iso = '1                         ISO\n'
                elastic_parameters.write(iso)
                Young_Nu = self.TMaterials.item(0,2).text() + '   ' + self.TMaterials.item(0,5).text() + '                  YOUNG(MPa),NU (V+R/2)\n'
                elastic_parameters.write(Young_Nu)
                elastic_parameters.close()
            elif 'Transversely Isotropic' in self.TMaterials.item(0,1).text():
                if str(self.INIsotropicPlane.currentText()) == 'x-y plane':
                    NUip = float(self.TMaterials.item(0,5).text())
                    NUop = float(self.TMaterials.item(0,6).text())
                    Eip = float(self.TMaterials.item(0,2).text())
                    Eop = float(self.TMaterials.item(0,4).text())
                    Gop = float(self.TMaterials.item(0,10).text())
                    S11 = 1/Eip
                    S12 = -NUip/Eip
                    S13 = -NUop/Eop
                    S22 = 1/Eip
                    S23 = -NUop/Eop
                    S33 = 1/Eop
                    S44 = 1/Gop
                    S55 = 1/Gop
                    S66 = 1/(Eip/(2*(1+NUip)))
                    S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                    C11 = (1/S)*(S22*S33-S23*S23)
                    C12 = (1/S)*(S13*S23-S12*S33)
                    C22 = (1/S)*(S33*S11-S13*S13)
                    C13 = (1/S)*(S12*S23-S13*S22)
                    C33 = (1/S)*(S11*S22-S12*S12)
                    C23 = (1/S)*(S12*S13-S23*S11)
                    C44 = (1/S44)
                    C55 = (1/S55)
                    C66 = (1/S66)
                    elastic_parameters = open("../example_input_files/lattice_input_files/elastic_parameters.txt","w")
                    iso = '0\n'
                    elastic_parameters.write(iso)
                    stiffness = '  ' + str(C11) + '  ' + str(C12) + '  ' + str(C13) + '  0  0  0     Cu (MPa)\n' + '  ' + str(C12) + '  ' + str(C22) + '  ' + str(C23) + '  0  0  0\n' + '  ' + str(C13) + '  ' + str(C23) + '  ' + str(C33) + '  0  0  0\n' + '  0  0  0   ' + str(C44) + '  0  0\n' + '  0  0  0  0  ' + str(C55) + '  0\n' + '  0  0  0  0  0  ' + str(C66)
                    elastic_parameters.write(stiffness)
                elif str(self.INIsotropicPlane.currentText()) == 'x-z plane':
                    NUip = float(self.TMaterials.item(0,6).text())
                    NUop = float(self.TMaterials.item(0,5).text())
                    Eip = float(self.TMaterials.item(0,2).text())
                    Eop = float(self.TMaterials.item(0,3).text())
                    Gop = float(self.TMaterials.item(0,8).text())
                    S11 = 1/Eip
                    S12 = -NUop/Eop
                    S13 = -NUip/Eip
                    S22 = 1/Eop
                    S23 = -NUop/Eop
                    S33 = 1/Eip
                    S44 = 1/Gop
                    S55 = 1/(Eip/(2*(1+NUip)))
                    S66 = 1/Gop
                    S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                    C11 = (1/S)*(S22*S33-S23*S23)
                    C12 = (1/S)*(S13*S23-S12*S33)
                    C22 = (1/S)*(S33*S11-S13*S13)
                    C13 = (1/S)*(S12*S23-S13*S22)
                    C33 = (1/S)*(S11*S22-S12*S12)
                    C23 = (1/S)*(S12*S13-S23*S11)
                    C44 = (1/S44)
                    C55 = (1/S55)
                    C66 = (1/S66)
                    elastic_parameters = open("../example_input_files/lattice_input_files/elastic_parameters.txt","w")
                    iso = '0\n'
                    elastic_parameters.write(iso)
                    stiffness = '  ' + str(C11) + '  ' + str(C12) + '  ' + str(C13) + '  0  0  0     Cu (MPa)\n' + '  ' + str(C12) + '  ' + str(C22) + '  ' + str(C23) + '  0  0  0\n' + '  ' + str(C13) + '  ' + str(C23) + '  ' + str(C33) + '  0  0  0\n' + '  0  0  0   ' + str(C44) + '  0  0\n' + '  0  0  0  0  ' + str(C55) + '  0\n' + '  0  0  0  0  0  ' + str(C66)
                    elastic_parameters.write(stiffness)
                else:
                    NUip = float(self.TMaterials.item(0,7).text())
                    NUop = float(self.TMaterials.item(0,6).text())
                    Eip = float(self.TMaterials.item(0,3).text())
                    Eop = float(self.TMaterials.item(0,2).text())
                    Gop = float(self.TMaterials.item(0,8).text())
                    S11 = 1/Eop
                    S12 = -NUop/Eop
                    S13 = -NUop/Eop
                    S22 = 1/Eip
                    S23 = -NUip/Eip
                    S33 = 1/Eip
                    S44 = 1/(Eip/(2*(1+NUip)))
                    S55 = 1/Gop
                    S66 = 1/Gop
                    S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                    C11 = (1/S)*(S22*S33-S23*S23)
                    C12 = (1/S)*(S13*S23-S12*S33)
                    C22 = (1/S)*(S33*S11-S13*S13)
                    C13 = (1/S)*(S12*S23-S13*S22)
                    C33 = (1/S)*(S11*S22-S12*S12)
                    C23 = (1/S)*(S12*S13-S23*S11)
                    C44 = (1/S44)
                    C55 = (1/S55)
                    C66 = (1/S66)
                    elastic_parameters = open("../example_input_files/lattice_input_files/elastic_parameters.txt","w")
                    iso = '0\n'
                    elastic_parameters.write(iso)
                    stiffness = '  ' + str(C11) + '  ' + str(C12) + '  ' + str(C13) + '  0  0  0     Cu (MPa)\n' + '  ' + str(C12) + '  ' + str(C22) + '  ' + str(C23) + '  0  0  0\n' + '  ' + str(C13) + '  ' + str(C23) + '  ' + str(C33) + '  0  0  0\n' + '  0  0  0   ' + str(C44) + '  0  0\n' + '  0  0  0  0  ' + str(C55) + '  0\n' + '  0  0  0  0  0  ' + str(C66)
                    elastic_parameters.write(stiffness)
            elif self.TMaterials.item(0,1).text() == 'Orthotropic':
                    S11 = 1/float(self.TMaterials.item(0,2).text())
                    S12 = -float(self.TMaterials.item(0,5).text())/float(self.TMaterials.item(0,2).text())
                    S13 = -float(self.TMaterials.item(0,6).text())/float(self.TMaterials.item(0,2).text())
                    S22 = 1/float(self.TMaterials.item(0,3).text())
                    S23 = -float(self.TMaterials.item(0,7).text())/float(self.TMaterials.item(0,3).text())
                    S33 = 1/float(self.TMaterials.item(0,4).text())
                    S44 = 1/float(self.TMaterials.item(0,10).text())
                    S55 = 1/float(self.TMaterials.item(0,9).text())
                    S66 = 1/float(self.TMaterials.item(0,8).text())
                    S = S11*S22*S33 - S11*S23*S23 - S22*S13*S13 - S33*S12*S12 + 2*S12*S23*S13
                    C11 = (1/S)*(S22*S33-S23*S23)
                    C12 = (1/S)*(S13*S23-S12*S33)
                    C22 = (1/S)*(S33*S11-S13*S13)
                    C13 = (1/S)*(S12*S23-S13*S22)
                    C33 = (1/S)*(S11*S22-S12*S12)
                    C23 = (1/S)*(S12*S13-S23*S11)
                    C44 = (1/S44)
                    C55 = (1/S55)
                    C66 = (1/S66)
                    elastic_parameters = open("../example_input_files/lattice_input_files/elastic_parameters.txt","w")
                    iso = '0\n'
                    elastic_parameters.write(iso)
                    stiffness = '  ' + str(C11) + '  ' + str(C12) + '  ' + str(C13) + '  0  0  0     Cu (MPa)\n' + '  ' + str(C12) + '  ' + str(C22) + '  ' + str(C23) + '  0  0  0\n' + '  ' + str(C13) + '  ' + str(C23) + '  ' + str(C33) + '  0  0  0\n' + '  0  0  0   ' + str(C44) + '  0  0\n' + '  0  0  0  0  ' + str(C55) + '  0\n' + '  0  0  0  0  0  ' + str(C66)
                    elastic_parameters.write(stiffness)
            
            # EVPFFT lattice input parameters file
            evpfft_lattice_input = open("../example_input_files/lattice_input_files/evpfft_lattice_input.txt","w")
            modes = '2 0 0 0               NPHMX, NMODMX, NTWMMX, NSYSMX\n'
            evpfft_lattice_input.write(modes)
            dimensions = str(int(self.INNumberOfVoxelsX.text())) + ' ' + str(int(self.INNumberOfVoxelsY.text())) + ' ' + str(int(self.INNumberOfVoxelsZ.text())) + '               x-dim, y-dim, z-dim\n'
            evpfft_lattice_input.write(dimensions)
            nph_delt = '2                      number of phases (nph)\n' + '1.  1.  1.             RVE dimensions (delt)\n' + '* name and path of microstructure file (filetext)\n'
            evpfft_lattice_input.write(nph_delt)
            vtkfile = current_file_loc + '/VTK_geometry.vtk\n'
            evpfft_lattice_input.write(vtkfile)
            phases = '*INFORMATION ABOUT PHASE #1\n' + '1                          igas(iph)\n' + '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + 'dummy\n' + 'dummy\n' + '*INFORMATION ABOUT PHASE #2\n' + '0                          igas(iph)\n' + '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + os.path.normpath(os.getcwd() + os.sep + os.pardir) +  '/example_input_files/lattice_input_files/plastic_parameters.txt\n' + os.path.normpath(os.getcwd() + os.sep + os.pardir) + '/example_input_files/lattice_input_files/elastic_parameters.txt\n'
            evpfft_lattice_input.write(phases)
            if str(self.INBCDirection.currentText()) == "x-direction":
                test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    1       1       1           iudot     |    flag for vel.grad.\n' + '    1       0       1                     |    (0:unknown-1:known)\n' + '    1       1       0                     |\n' + '                                          |\n' + '   1.0     0.        0.          udot    |    vel.grad\n' + '    0.      -0.35      0.                  |\n' + '    0.       0.         -0.35                |\n' + '                                          |\n' + '    0       0        0           iscau    |    flag for Cauchy\n' + '            1        0                    |\n' + '                     1                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
            elif str(self.INBCDirection.currentText()) == "y-direction":
                test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    0       1       1           iudot     |    flag for vel.grad.\n' + '    1       1       1                     |    (0:unknown-1:known)\n' + '    1       1       0                     |\n' + '                                          |\n' + '   -0.35     0.        0.          udot    |    vel.grad\n' + '    0.      1.0      0.                  |\n' + '    0.       0.         -0.35                |\n' + '                                          |\n' + '    1       0        0           iscau    |    flag for Cauchy\n' + '            0        0                    |\n' + '                     1                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
            else:
                test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + '* boundary conditions\n' + '    0       1       1           iudot     |    flag for vel.grad.\n' + '    1       0       1                     |    (0:unknown-1:known)\n' + '    1       1       1                     |\n' + '                                          |\n' + '   -0.35     0.        0.          udot    |    vel.grad\n' + '    0.      -0.35      0.                  |\n' + '    0.       0.         1.0                |\n' + '                                          |\n' + '    1       0        0           iscau    |    flag for Cauchy\n' + '            1        0                    |\n' + '                     0                    |\n' + '                                          |\n' + '    0.      0.       0.          scauchy  |    Cauchy stress\n' + '            0.       0.                   |\n' + '                     0.                   @\n'
            evpfft_lattice_input.write(test_conditions)
            other = '* other\n' + '0.0001         eqincr (if ictrl>=0) or tdot (if ictrl=-1)\n' + '-1              ictrl (1-6: strain comp, 0: VM eq, -1: tdot)\n'
            evpfft_lattice_input.write(other)
            run_conditions = '*INFORMATION ABOUT RUN CONDITIONS\n' + self.INNumberOfSteps.text() + '              nsteps\n' + '0.00001         err\n' + '50              itmax\n' + '0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?\n' + '0               ISAVE write grain states in STRESS.OUT (1) or not (0)?\n' + '1               IUPDATE update tex & RVE dim (1) or not (0)?\n' + '0               IUPHARD\n' + '1               IWTEX\n' + '1 10            IWFIELDS,IWSTEP\n' + '0               ITHERMO (if ithermo=1, next line is filethermo)\n' + 'dummy\n'
            evpfft_lattice_input.write(run_conditions)
            evpfft_lattice_input.close()
            
        # Run EVPFFT
        def run_click():
            write_input_file()
            self.p = QProcess()
            self.p.readyReadStandardOutput.connect(handle_stdout)
            self.p.readyReadStandardError.connect(handle_stderr)
            self.p.stateChanged.connect(handle_state)
            self.p.finished.connect(process_finished)
            self.p.start("./evpfft",["-f", "../example_input_files/lattice_input_files/evpfft_lattice_input.txt", "-m", "2"])
            self.progress_re = re.compile("       Current  Time  STEP = (\d+)")
        def simple_percent_parser(output):
            m = self.progress_re.search(output)
            if m:
                pc_complete = m.group(1)
                return int(pc_complete)
        def process_finished():
            self.RunOutputProgress.setValue(100)
            self.p = None
        def handle_stdout():
            data = self.p.readAllStandardOutput()
            stdout = bytes(data).decode("utf8")
            progress = simple_percent_parser(stdout)
            if progress:
                self.RunOutputProgress.setValue((progress/int(self.INNumberOfSteps.text()))*100)
            self.RunOutputWindow.appendPlainText(stdout)
        def handle_stderr():
            data = self.p.readAllStandardError()
            stderr = bytes(data).decode("utf8")
            self.RunOutputWindow.appendPlainText(stderr)
        def handle_state(state):
            states = {
                QProcess.NotRunning: 'Finished',
                QProcess.Starting: 'Starting EVPFFT',
                QProcess.Running: 'Running EVPFFT',
            }
            state_name = states[state]
            self.RunOutputWindow.appendPlainText(f"{state_name}")
            
        self.BRunEVPFFT.clicked.connect(run_click)
        
        
        # Preview Results
        def preview_results_click():
            # Delete previous views
            try:
                self.stl
            except:
                print('')
            else:
                pvsimple.Delete(self.stl)
            try:
                self.threshold
            except:
                print('')
            else:
                pvsimple.Delete(self.threshold)
            try:
                self.threshold2
            except:
                print('')
            else:
                pvsimple.Delete(self.threshold2)
            # Render new view
            self.results_reader = pvsimple.XDMFReader(FileNames = "micro_state_timestep_10.xdmf")
            pvsimple.SetDisplayProperties(Representation = "Surface")
            self.threshold2 = pvsimple.Threshold(registrationName='results_threshold', Input = self.results_reader, Scalars = "phase_id", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 2, LowerThreshold = 1, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
            display = pvsimple.Show(self.threshold2, render_view)
            # Select what variable you want to display
            pvsimple.GetAnimationScene().GoToLast()
            pvsimple.ColorBy(display,('CELLS',str(self.INPreviewResults.currentText())))
            vmstressLUT = pvsimple.GetColorTransferFunction(str(self.INPreviewResults.currentText()))
            r = self.results_reader.CellData.GetArray(str(self.INPreviewResults.currentText())).GetRange()
            vmstressLUT.RescaleTransferFunction(r[0], r[1]/2)
            display.SetScalarBarVisibility(render_view, True)
            pvsimple.HideUnusedScalarBars(render_view)
            # Add time filter
            threshold1 = pvsimple.FindSource('results_threshold')
            annotateTimeFilter1 = pvsimple.AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=threshold1)
            annotateTimeFilter1Display = pvsimple.Show(annotateTimeFilter1, render_view, 'TextSourceRepresentation')
            # Remove old view / reset cameras
            pvsimple.Hide(self.results_reader)
            render_view.ResetCamera()
            render_view.StillRender()
        self.BPreviewResults.clicked.connect(preview_results_click)
        self.BPreviewResults.clicked.connect(lambda: self.OutputWindows.setCurrentIndex(0))
        
        # Stress vs Strain Plot
        def plot_ss_click():
            # Get the stress-strain data
            try:
                self.Plot.figure
            except:
                with open("str_str.out", newline='') as f:
                    reader = csv.reader(f)
                    self.ss_data = list(reader)
                x = [0 for i in range(int(self.INNumberOfSteps.text()))]
                y = [0 for i in range(int(self.INNumberOfSteps.text()))]
                for i in range(int(self.INNumberOfSteps.text())):
                    if str(self.INPlotSS.currentText()) == 'S11 vs E11':
                        xcol = 0
                        ycol = 6
                    elif str(self.INPlotSS.currentText()) == 'S22 vs E22':
                        xcol = 1
                        ycol = 7
                    elif str(self.INPlotSS.currentText()) == 'S33 vs E33':
                        xcol = 2
                        ycol = 8
                    x[i] = float(self.ss_data[i+1][xcol])
                    y[i] = float(self.ss_data[i+1][ycol])
                # Plot data
                self.Plot.figure = Figure()
                self.Plot.ax = self.Plot.figure.add_subplot()
                self.Plot.ax.plot(x,y)
                self.Plot.ax.set_xlabel('STRAIN')
                self.Plot.ax.set_ylabel('STRESS')
                self.Plot.figure.tight_layout()
                # Display plot and toolbar
                layout = QVBoxLayout()
                self.Plot.setLayout(layout)
                self.Plot.canvas = FigureCanvasQTAgg(self.Plot.figure)
                layout.addWidget(self.Plot.canvas)
                self.toolbar = NavigationToolbar2QT(self.Plot.canvas,self.Plot)
                layout.addWidget(self.toolbar)
            else:
                self.timer = QTimer()
                self.timer.setInterval(100)
                self.timer.timeout.connect(update_plot)
                self.timer.start()
        def update_plot():
            x = [0 for i in range(int(self.INNumberOfSteps.text()))]
            y = [0 for i in range(int(self.INNumberOfSteps.text()))]
            for i in range(int(self.INNumberOfSteps.text())):
                if str(self.INPlotSS.currentText()) == 'S11 vs E11':
                    xcol = 0
                    ycol = 6
                elif str(self.INPlotSS.currentText()) == 'S22 vs E22':
                    xcol = 1
                    ycol = 7
                elif str(self.INPlotSS.currentText()) == 'S33 vs E33':
                    xcol = 2
                    ycol = 8
                x[i] = float(self.ss_data[i+1][xcol])
                y[i] = float(self.ss_data[i+1][ycol])
            self.Plot.ax.cla()
            self.Plot.ax.plot(x,y)
            self.Plot.ax.set_xlabel('STRAIN')
            self.Plot.ax.set_ylabel('STRESS')
            self.Plot.figure.tight_layout()
            self.Plot.canvas.draw()
            self.timer.stop()
        self.BPlotSS.clicked.connect(plot_ss_click)
        self.BPlotSS.clicked.connect(lambda: self.OutputWindows.setCurrentIndex(1))
        
        # Open Paraview
        def open_paraview_click():
            os.system("open " + "micro_state_timestep_10.xdmf")
        self.BOpenParaview.clicked.connect(open_paraview_click)
        
        # Warning Message Popup
        def warning_message(msg):
            message = QMessageBox()
            message.setText(msg)
            message.exec()
        
    # =================================================================
        
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"EVPFFT", None))
        self.actionEVPFFT_Manual.setText(QCoreApplication.translate("MainWindow", u"EVPFFT Manual", None))
        self.BImportPart.setText("")
        self.LImportPart.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Import Part</span></p></body></html>", None))
        self.LPartTools.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Part Tools</span></p></body></html>", None))
        self.BDefineMaterial.setText("")
        self.LDefineMaterial.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Define Material</span></p></body></html>", None))
        self.LMaterialTools.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Material Tools</span></p></body></html>", None))
        self.BApplyBC.setText("")
        self.LApplyBC.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Apply BCs</span></p></body></html>", None))
        self.LBCTools.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Boundary Condition Tools</span></p></body></html>", None))
        self.BSolverSettings.setText("")
        self.LSolverSettings.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Solver Settings</span></p></body></html>", None))
        self.BRunEVPFFT.setText("")
        self.LRunEVPFFT.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Run EVPFFT</span></p></body></html>", None))
        self.BViewResults.setText("")
        self.LViewResults.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">View Results</span></p></body></html>", None))
        self.LJobTools.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Job Tools</span></p></body></html>", None))
        self.HeaderMenu.setTabText(self.HeaderMenu.indexOf(self.Tools), QCoreApplication.translate("MainWindow", u"Tools", None))
        self.LosAlamosLogo.setText("")
        self.EVPFFTLogo.setText("")
        self.LAdditionalSoftware.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#000000;\">Additional Software Packages:</span></p></body></html>", None))
        self.MatarLogo.setText("")
        self.ParaviewLogo.setText("")
        self.LGeometryInformation.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">IMPORT PART</span></p></body></html>", None))
        self.BUploadGeometryFile.setText(QCoreApplication.translate("MainWindow", u"Upload Geometry File", None))
        self.STLVoxelization.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">.STL VOXELIZATION</span></p></body></html>", None))
        self.LNumberOfVoxelsX.setText(QCoreApplication.translate("MainWindow", u"Number of voxels x: ", None))
        self.LNumberOfVoxelsY.setText(QCoreApplication.translate("MainWindow", u"Number of voxels y: ", None))
        self.LNumberOfVoxelsZ.setText(QCoreApplication.translate("MainWindow", u"Number of voxels z: ", None))
        self.BVoxelizeGeometry.setText(QCoreApplication.translate("MainWindow", u"Voxelize Geometry", None))
        self.LDefineMaterials.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">DEFINE MATERIALS</span></p></body></html>", None))
        self.LMaterialName.setText(QCoreApplication.translate("MainWindow", u"Name: ", None))
        self.INMaterialType.setItemText(0, QCoreApplication.translate("MainWindow", u"Isotropic", None))
        self.INMaterialType.setItemText(1, QCoreApplication.translate("MainWindow", u"Transversely Isotropic", None))
        self.INMaterialType.setItemText(2, QCoreApplication.translate("MainWindow", u"Orthotropic", None))
        self.INMaterialType.setItemText(3, QCoreApplication.translate("MainWindow", u"Anisotropic", None))

        self.LYoungsModulus.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E: </p></body></html>", None))
        self.LPoissonsRatio.setText(QCoreApplication.translate("MainWindow", u"nu: ", None))
        self.LIsotropicPlane.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Isotropic Plane: </p></body></html>", None))
        self.INIsotropicPlane.setItemText(0, QCoreApplication.translate("MainWindow", u"x-y plane", None))
        self.INIsotropicPlane.setItemText(1, QCoreApplication.translate("MainWindow", u"x-z plane", None))
        self.INIsotropicPlane.setItemText(2, QCoreApplication.translate("MainWindow", u"y-z plane", None))

        self.LInPlane.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" text-decoration: underline;\">In-Plane</span></p></body></html>", None))
        self.LEip.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">E<span style=\" vertical-align:sub;\">IP</span>: </p></body></html>", None))
        self.LNUip.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">nu<span style=\" vertical-align:sub;\">IP</span>: </p></body></html>", None))
        self.LOutOfPlane.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" text-decoration: underline;\">Out-Of-Plane</span></p></body></html>", None))
        self.LEop.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">E<span style=\" vertical-align:sub;\">OP</span>: </p></body></html>", None))
        self.LNUop.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">nu<span style=\" vertical-align:sub;\">OP</span>: </p></body></html>", None))
        self.LGop.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">G<span style=\" vertical-align:sub;\">OP</span>: </p></body></html>", None))
        self.label.setText(QCoreApplication.translate("MainWindow", u"Define The Stiffness Matrix [C]:", None))

        __sortingEnabled = self.TAnisotropic.isSortingEnabled()
        self.TAnisotropic.setSortingEnabled(False)
        self.TAnisotropic.setSortingEnabled(__sortingEnabled)

        self.LGxy.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>G<span style=\" vertical-align:sub;\">xy</span>: </p></body></html>", None))
        self.LEx.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E<span style=\" vertical-align:sub;\">x</span>: </p></body></html>", None))
        self.LNUyz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>nu<span style=\" vertical-align:sub;\">yz</span>: </p></body></html>", None))
        self.LGyz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>G<span style=\" vertical-align:sub;\">yz</span>: </p></body></html>", None))
        self.LEy.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E<span style=\" vertical-align:sub;\">y</span>: </p></body></html>", None))
        self.LEz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E<span style=\" vertical-align:sub;\">z</span>: </p></body></html>", None))
        self.LGxz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>G<span style=\" vertical-align:sub;\">xz</span>: </p></body></html>", None))
        self.LNUxz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>nu<span style=\" vertical-align:sub;\">xz</span>: </p></body></html>", None))
        self.LNUxy.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>nu<span style=\" vertical-align:sub;\">xy</span>: </p></body></html>", None))
        self.BAddMaterial.setText(QCoreApplication.translate("MainWindow", u"Add", None))
        ___qtablewidgetitem = self.TMaterials.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem1 = self.TMaterials.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("MainWindow", u"Definition", None));
        ___qtablewidgetitem2 = self.TMaterials.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("MainWindow", u"Ex", None));
        ___qtablewidgetitem3 = self.TMaterials.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("MainWindow", u"Ey", None));
        ___qtablewidgetitem4 = self.TMaterials.horizontalHeaderItem(4)
        ___qtablewidgetitem4.setText(QCoreApplication.translate("MainWindow", u"Ez", None));
        ___qtablewidgetitem5 = self.TMaterials.horizontalHeaderItem(5)
        ___qtablewidgetitem5.setText(QCoreApplication.translate("MainWindow", u"NUxy", None));
        ___qtablewidgetitem6 = self.TMaterials.horizontalHeaderItem(6)
        ___qtablewidgetitem6.setText(QCoreApplication.translate("MainWindow", u"NUxz", None));
        ___qtablewidgetitem7 = self.TMaterials.horizontalHeaderItem(7)
        ___qtablewidgetitem7.setText(QCoreApplication.translate("MainWindow", u"NUyz", None));
        ___qtablewidgetitem8 = self.TMaterials.horizontalHeaderItem(8)
        ___qtablewidgetitem8.setText(QCoreApplication.translate("MainWindow", u"Gxy", None));
        ___qtablewidgetitem9 = self.TMaterials.horizontalHeaderItem(9)
        ___qtablewidgetitem9.setText(QCoreApplication.translate("MainWindow", u"Gxz", None));
        ___qtablewidgetitem10 = self.TMaterials.horizontalHeaderItem(10)
        ___qtablewidgetitem10.setText(QCoreApplication.translate("MainWindow", u"Gyz", None));
        self.BDeleteMaterial.setText(QCoreApplication.translate("MainWindow", u"Delete", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"Stiffness Matrix [C]:", None))

        __sortingEnabled1 = self.TStiffnessMatrix.isSortingEnabled()
        self.TStiffnessMatrix.setSortingEnabled(False)
        self.TStiffnessMatrix.setSortingEnabled(__sortingEnabled1)

        self.LBoundaryConditions.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">BOUNDARY CONDITIONS</span></p></body></html>", None))
        self.LBoundaryCondition.setText(QCoreApplication.translate("MainWindow", u"Boundary Condition: ", None))
        self.INBoundaryCondition.setItemText(0, QCoreApplication.translate("MainWindow", u"Tension", None))
        self.INBoundaryCondition.setItemText(1, QCoreApplication.translate("MainWindow", u"Compression", None))

        self.LBCDirection.setText(QCoreApplication.translate("MainWindow", u"Direction: ", None))
        self.INBCDirection.setItemText(0, QCoreApplication.translate("MainWindow", u"x-direction", None))
        self.INBCDirection.setItemText(1, QCoreApplication.translate("MainWindow", u"y-direction", None))
        self.INBCDirection.setItemText(2, QCoreApplication.translate("MainWindow", u"z-direction", None))

        self.BAddBC.setText(QCoreApplication.translate("MainWindow", u"Add", None))
        ___qtablewidgetitem11 = self.TBCs.horizontalHeaderItem(0)
        ___qtablewidgetitem11.setText(QCoreApplication.translate("MainWindow", u"Boundary Condition", None));
        ___qtablewidgetitem12 = self.TBCs.horizontalHeaderItem(1)
        ___qtablewidgetitem12.setText(QCoreApplication.translate("MainWindow", u"Direction", None));
        self.BDeleteBC.setText(QCoreApplication.translate("MainWindow", u"Delete", None))
        self.LSolverSettings_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">SOLVER SETTINGS</span></p></body></html>", None))
        self.LNumberOfSteps.setText(QCoreApplication.translate("MainWindow", u"Number of steps: ", None))
        self.INNumberOfSteps.setInputMask("")
        self.INNumberOfSteps.setText(QCoreApplication.translate("MainWindow", u"10", None))
        self.INNumberOfSteps.setPlaceholderText("")
        self.LResults.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">RESULTS</span></p></body></html>", None))
        self.INPreviewResults.setItemText(0, QCoreApplication.translate("MainWindow", u"vm-stress", None))
        self.INPreviewResults.setItemText(1, QCoreApplication.translate("MainWindow", u"vm-strain", None))

        self.BPreviewResults.setText(QCoreApplication.translate("MainWindow", u"Preview Results", None))
        self.BPlotSS.setText(QCoreApplication.translate("MainWindow", u"Plot Stress vs Strain", None))
        self.INPlotSS.setItemText(0, QCoreApplication.translate("MainWindow", u"S11 vs E11", None))
        self.INPlotSS.setItemText(1, QCoreApplication.translate("MainWindow", u"S22 vs E22", None))
        self.INPlotSS.setItemText(2, QCoreApplication.translate("MainWindow", u"S33 vs E33", None))

        self.BOpenParaview.setText(QCoreApplication.translate("MainWindow", u"Open Paraview", None))
        self.menuHelp.setTitle(QCoreApplication.translate("MainWindow", u"Help", None))
    # retranslateUi

