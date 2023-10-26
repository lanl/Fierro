# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'EVPFFT_GUI_V10VZaXDX.ui'
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
from PySide6.QtWidgets import (QApplication, QComboBox, QFormLayout, QFrame,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit,
    QMainWindow, QMenu, QMenuBar, QPlainTextEdit,
    QProgressBar, QPushButton, QSizePolicy, QSpacerItem,
    QStackedWidget, QStatusBar, QTabWidget, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget)

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

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1009, 859)
        icon = QIcon()
        icon.addFile(LocalResource.get_resource_name("Icons/EVPFFT_logo_A2.png"), QSize(), QIcon.Normal, QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setAutoFillBackground(False)
        MainWindow.setToolButtonStyle(Qt.ToolButtonIconOnly)
        MainWindow.setDockNestingEnabled(False)
        self.actionEVPFFT_Manual = QAction(MainWindow)
        self.actionEVPFFT_Manual.setObjectName(u"actionEVPFFT_Manual")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.centralwidget.setStyleSheet(u"#TitlePage, #GeometryInformationTool, #DefineMaterialTool, #BoundaryConditionsTool, #SolverSettingsTool, #ResultsTool, #Tools, #RunOutputs, #RunOutputWindow{\n"
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
        icon1.addFile(LocalResource.get_resource_name("Icons/Blue/Cube.svg"), QSize(), QIcon.Normal, QIcon.Off)
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
        icon2.addFile(LocalResource.get_resource_name("Icons/Blue/mine.svg"), QSize(), QIcon.Normal, QIcon.Off)
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
        icon3.addFile(LocalResource.get_resource_name("Icons/Blue/brick.svg"), QSize(), QIcon.Normal, QIcon.Off)
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
        icon4.addFile(LocalResource.get_resource_name("Icons/Blue/gear.svg"), QSize(), QIcon.Normal, QIcon.Off)
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
        icon5.addFile(LocalResource.get_resource_name("Icons/Blue/Play.svg"), QSize(), QIcon.Normal, QIcon.Off)
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
        icon6.addFile(LocalResource.get_resource_name("Icons/Blue/magnify.svg"), QSize(), QIcon.Normal, QIcon.Off)
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

        self.line_4 = QFrame(self.centralwidget)
        self.line_4.setObjectName(u"line_4")
        self.line_4.setFrameShadow(QFrame.Plain)
        self.line_4.setLineWidth(50)
        self.line_4.setFrameShape(QFrame.HLine)

        self.verticalLayout.addWidget(self.line_4)

        self.Main = QFrame(self.centralwidget)
        self.Main.setObjectName(u"Main")
        self.Main.setMinimumSize(QSize(0, 0))
        self.Main.setFrameShape(QFrame.NoFrame)
        self.Main.setFrameShadow(QFrame.Raised)
        self.horizontalLayout = QHBoxLayout(self.Main)
        self.horizontalLayout.setSpacing(0)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.ToolSettings = QStackedWidget(self.Main)
        self.ToolSettings.setObjectName(u"ToolSettings")
        self.ToolSettings.setMinimumSize(QSize(0, 0))
        self.ToolSettings.setMaximumSize(QSize(300, 16777215))
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
        self.LosAlamosLogo.setPixmap(QPixmap(LocalResource.get_resource_name("Icons/LANL Logo Ultramarine.png")))
        self.LosAlamosLogo.setScaledContents(True)

        self.verticalLayout_2.addWidget(self.LosAlamosLogo)

        self.verticalSpacer_7 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_2.addItem(self.verticalSpacer_7)

        self.EVPFFTLogo = QLabel(self.TitlePage)
        self.EVPFFTLogo.setObjectName(u"EVPFFTLogo")
        self.EVPFFTLogo.setMinimumSize(QSize(275, 175))
        self.EVPFFTLogo.setMaximumSize(QSize(275, 175))
        self.EVPFFTLogo.setPixmap(QPixmap(LocalResource.get_resource_name("Icons/EVPFFT_logo_horse_ppt.png")))
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
        self.MatarLogo.setPixmap(QPixmap(LocalResource.get_resource_name("Icons/MATAR_logo2.png")))
        self.MatarLogo.setScaledContents(True)

        self.horizontalLayout_5.addWidget(self.MatarLogo)

        self.ParaviewLogo = QLabel(self.AdditionalSoftwareLogos)
        self.ParaviewLogo.setObjectName(u"ParaviewLogo")
        self.ParaviewLogo.setMaximumSize(QSize(130, 30))
        self.ParaviewLogo.setPixmap(QPixmap(LocalResource.get_resource_name("Icons/ParaView_logo.png")))
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
        self.formLayout_2 = QFormLayout(self.MaterialInputs)
        self.formLayout_2.setObjectName(u"formLayout_2")
        self.LMaterialName = QLabel(self.MaterialInputs)
        self.LMaterialName.setObjectName(u"LMaterialName")

        self.formLayout_2.setWidget(0, QFormLayout.LabelRole, self.LMaterialName)

        self.INMaterialName = QLineEdit(self.MaterialInputs)
        self.INMaterialName.setObjectName(u"INMaterialName")

        self.formLayout_2.setWidget(0, QFormLayout.FieldRole, self.INMaterialName)

        self.LYoungsModulus = QLabel(self.MaterialInputs)
        self.LYoungsModulus.setObjectName(u"LYoungsModulus")

        self.formLayout_2.setWidget(1, QFormLayout.LabelRole, self.LYoungsModulus)

        self.INYoungsModulus = QLineEdit(self.MaterialInputs)
        self.INYoungsModulus.setObjectName(u"INYoungsModulus")

        self.formLayout_2.setWidget(1, QFormLayout.FieldRole, self.INYoungsModulus)

        self.LPoissonsRatio = QLabel(self.MaterialInputs)
        self.LPoissonsRatio.setObjectName(u"LPoissonsRatio")

        self.formLayout_2.setWidget(2, QFormLayout.LabelRole, self.LPoissonsRatio)

        self.INPoissonsRatio = QLineEdit(self.MaterialInputs)
        self.INPoissonsRatio.setObjectName(u"INPoissonsRatio")

        self.formLayout_2.setWidget(2, QFormLayout.FieldRole, self.INPoissonsRatio)

        self.LApplyToRegion = QLabel(self.MaterialInputs)
        self.LApplyToRegion.setObjectName(u"LApplyToRegion")

        self.formLayout_2.setWidget(3, QFormLayout.LabelRole, self.LApplyToRegion)

        self.INApplyToRegion = QComboBox(self.MaterialInputs)
        self.INApplyToRegion.addItem("")
        self.INApplyToRegion.setObjectName(u"INApplyToRegion")

        self.formLayout_2.setWidget(3, QFormLayout.FieldRole, self.INApplyToRegion)


        self.verticalLayout_16.addWidget(self.MaterialInputs)

        self.BAddMaterial = QPushButton(self.DefineMaterialTool)
        self.BAddMaterial.setObjectName(u"BAddMaterial")

        self.verticalLayout_16.addWidget(self.BAddMaterial)

        self.TMaterials = QTableWidget(self.DefineMaterialTool)
        if (self.TMaterials.columnCount() < 4):
            self.TMaterials.setColumnCount(4)
        __qtablewidgetitem = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        __qtablewidgetitem2 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(2, __qtablewidgetitem2)
        __qtablewidgetitem3 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(3, __qtablewidgetitem3)
        self.TMaterials.setObjectName(u"TMaterials")
        self.TMaterials.setRowCount(0)

        self.verticalLayout_16.addWidget(self.TMaterials)

        self.BDeleteMaterial = QPushButton(self.DefineMaterialTool)
        self.BDeleteMaterial.setObjectName(u"BDeleteMaterial")

        self.verticalLayout_16.addWidget(self.BDeleteMaterial)

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
        if (self.TBCs.columnCount() < 3):
            self.TBCs.setColumnCount(3)
        __qtablewidgetitem4 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(0, __qtablewidgetitem4)
        __qtablewidgetitem5 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(1, __qtablewidgetitem5)
        __qtablewidgetitem6 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(2, __qtablewidgetitem6)
        self.TBCs.setObjectName(u"TBCs")

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

        self.horizontalLayout.addWidget(self.ToolSettings)

        self.line_5 = QFrame(self.Main)
        self.line_5.setObjectName(u"line_5")
        self.line_5.setFrameShadow(QFrame.Plain)
        self.line_5.setLineWidth(20)
        self.line_5.setFrameShape(QFrame.VLine)

        self.horizontalLayout.addWidget(self.line_5)

        self.ParaviewFrame = QFrame(self.Main)
        self.ParaviewFrame.setObjectName(u"ParaviewFrame")
        self.ParaviewFrame.setFrameShape(QFrame.NoFrame)
        self.ParaviewFrame.setFrameShadow(QFrame.Raised)
        self.verticalLayout_19 = QVBoxLayout(self.ParaviewFrame)
        self.verticalLayout_19.setSpacing(0)
        self.verticalLayout_19.setObjectName(u"verticalLayout_19")
        self.verticalLayout_19.setContentsMargins(0, 0, 0, 0)
        self.OutputWindows = QStackedWidget(self.ParaviewFrame)
        self.OutputWindows.setObjectName(u"OutputWindows")
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

        self.verticalLayout_19.addWidget(self.OutputWindows)

        self.RunOutputs = QFrame(self.ParaviewFrame)
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


        self.verticalLayout_19.addWidget(self.RunOutputs)


        self.horizontalLayout.addWidget(self.ParaviewFrame)


        self.verticalLayout.addWidget(self.Main)

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1009, 24))
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
        self.ToolSettings.setCurrentIndex(0)
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
                shutil.copy(b3_filename[0], VTK_OUTPUT)
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
                
                fierro_voxelizer.create_voxel_vtk(
                    b3_filename[0],
                    VTK_OUTPUT,
                    int(self.INNumberOfVoxelsX.text()),
                    int(self.INNumberOfVoxelsY.text()),
                    int(self.INNumberOfVoxelsZ.text()),
                    True
                )
                # Paraview window
                pvsimple.Delete(self.stl)
                self.voxel_reader = pvsimple.LegacyVTKReader(FileNames = VTK_OUTPUT)
                pvsimple.SetDisplayProperties(Representation = "Surface")
                self.threshold = pvsimple.Threshold(Input = self.voxel_reader, Scalars = "density", ThresholdMethod = "Above Upper Threshold", UpperThreshold = 1, LowerThreshold = 0, AllScalars = 1, UseContinuousCellRange = 0, Invert = 0)
                pvsimple.Show(self.threshold, render_view)
                pvsimple.Hide(self.voxel_reader)
                render_view.ResetCamera()
                render_view.StillRender()
        self.BVoxelizeGeometry.clicked.connect(voxelize_geometry_click)
        
        # Apply Material
        def add_material():
            if not self.INMaterialName.text() or not self.INYoungsModulus.text() or not self.INPoissonsRatio.text():
                warning_message('ERROR: Material definition incomplete')
            else:
                try:
                    self.stl
                except:
                    warning_message('ERROR: No part defined')
                else:
                    row = self.TMaterials.rowCount()
                    self.TMaterials.insertRow(row)
                    self.TMaterials.setItem(row, 0, QTableWidgetItem(
                        self.INMaterialName.text().strip())
                    )
                    self.TMaterials.setItem(
                        row, 1, QTableWidgetItem(self.INYoungsModulus.text())
                    )
                    self.TMaterials.setItem(
                        row, 2, QTableWidgetItem(self.INPoissonsRatio.text())
                    )
                    self.TMaterials.setItem(
                        row, 3, QTableWidgetItem(str(self.INApplyToRegion.currentText()))
                    )
                    reset_material()
            
        def reset_material():
            self.INMaterialName.clear()
            self.INYoungsModulus.clear()
            self.INPoissonsRatio.clear()
                
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
            
            # Elastic parameters file
            elastic_parameters = open(ELASTIC_PARAMETERS,"w")
            iso = '1                         ISO\n'
            elastic_parameters.write(iso)
            Young_Nu = self.TMaterials.item(0,1).text() + '   ' + self.TMaterials.item(0,2).text() + '                  YOUNG(MPa),NU (V+R/2)\n'
            elastic_parameters.write(Young_Nu)
            elastic_parameters.close()
            
            # EVPFFT lattice input parameters file
            evpfft_lattice_input = open(EVPFFT_INPUT,"w")
            modes = '2 0 0 0               NPHMX, NMODMX, NTWMMX, NSYSMX\n'
            evpfft_lattice_input.write(modes)
            dimensions = str(int(self.INNumberOfVoxelsX.text())) + ' ' + str(int(self.INNumberOfVoxelsY.text())) + ' ' + str(int(self.INNumberOfVoxelsZ.text())) + '               x-dim, y-dim, z-dim\n'
            evpfft_lattice_input.write(dimensions)
            nph_delt = '2                      number of phases (nph)\n' + '1.  1.  1.             RVE dimensions (delt)\n' + '* name and path of microstructure file (filetext)\n'
            evpfft_lattice_input.write(nph_delt)
            evpfft_lattice_input.write(f'{VTK_OUTPUT}\n')
            phases = '*INFORMATION ABOUT PHASE #1\n' + '1                          igas(iph)\n' + '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + 'dummy\n' + 'dummy\n' + '*INFORMATION ABOUT PHASE #2\n' + '0                          igas(iph)\n' + '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + f'{PLASTIC_PARAMETERS}\n' + f'{ELASTIC_PARAMETERS}\n'
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
            self.p.start("evpfft",["-f", EVPFFT_INPUT, "-m", "2"])
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
        self.LYoungsModulus.setText(QCoreApplication.translate("MainWindow", u"Young's Modulus: ", None))
        self.LPoissonsRatio.setText(QCoreApplication.translate("MainWindow", u"Poisson's Ratio: ", None))
        self.LApplyToRegion.setText(QCoreApplication.translate("MainWindow", u"Apply to Region: ", None))
        self.INApplyToRegion.setItemText(0, QCoreApplication.translate("MainWindow", u"Whole Geometry", None))

        self.BAddMaterial.setText(QCoreApplication.translate("MainWindow", u"Add", None))
        ___qtablewidgetitem = self.TMaterials.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem1 = self.TMaterials.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("MainWindow", u"Young's Modulus", None));
        ___qtablewidgetitem2 = self.TMaterials.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("MainWindow", u"Poisson's Ratio", None));
        ___qtablewidgetitem3 = self.TMaterials.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("MainWindow", u"Region", None));
        self.BDeleteMaterial.setText(QCoreApplication.translate("MainWindow", u"Delete", None))
        self.LBoundaryConditions.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">BOUNDARY CONDITIONS</span></p></body></html>", None))
        self.LBoundaryCondition.setText(QCoreApplication.translate("MainWindow", u"Boundary Condition: ", None))
        self.INBoundaryCondition.setItemText(0, QCoreApplication.translate("MainWindow", u"Tension", None))

        self.LBCDirection.setText(QCoreApplication.translate("MainWindow", u"Direction: ", None))
        self.INBCDirection.setItemText(0, QCoreApplication.translate("MainWindow", u"x-direction", None))
        self.INBCDirection.setItemText(1, QCoreApplication.translate("MainWindow", u"y-direction", None))
        self.INBCDirection.setItemText(2, QCoreApplication.translate("MainWindow", u"z-direction", None))

        self.BAddBC.setText(QCoreApplication.translate("MainWindow", u"Add", None))
        ___qtablewidgetitem4 = self.TBCs.horizontalHeaderItem(0)
        ___qtablewidgetitem4.setText(QCoreApplication.translate("MainWindow", u"Boundary Condition", None));
        ___qtablewidgetitem5 = self.TBCs.horizontalHeaderItem(1)
        ___qtablewidgetitem5.setText(QCoreApplication.translate("MainWindow", u"Direction", None));
        ___qtablewidgetitem6 = self.TBCs.horizontalHeaderItem(2)
        ___qtablewidgetitem6.setText(QCoreApplication.translate("MainWindow", u"Strain Rate", None));
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

