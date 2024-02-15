# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'FIERRO_GUIlBwFcJ.ui'
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
    QRadioButton, QSizePolicy, QSpacerItem, QSplitter,
    QStackedWidget, QStatusBar, QTabWidget, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget)
import IconResourceFile_rc

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1074, 957)
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
"    background-color: rgb(235, 235, 235);\n"
"}\n"
"#ParaviewFrame{\n"
"    background-color: rgb(91, 97, 120);\n"
"}\n"
"#BImportPart:hover, #BDefineMaterial:hover, #BApplyBC:hover, #BSolverSettings:hover, #BRunEVPFFT:hover, #BViewResults:hover, #BGlobalMesh:hover{\n"
"    background-color: rgb(192, 192, 192);\n"
"    border-radius: 15px;\n"
"}\n"
"#BImportPart, #BDefineMaterial, #BApplyBC, #BSolverSettings, #BRunEVPFFT, #BViewResults, #BGlobalMesh{\n"
"    border-style: flat;\n"
"}\n"
"#centralwidget{\n"
"    background-color: rgb(255, 255, 255);\n"
"}\n"
"\n"
"")
        self.verticalLayout = QVBoxLayout(self.centralwidget)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.verticalLayout.setContentsMargins(-1, 0, -1, -1)
        self.frame_6 = QFrame(self.centralwidget)
        self.frame_6.setObjectName(u"frame_6")
        self.frame_6.setFrameShape(QFrame.NoFrame)
        self.frame_6.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_10 = QHBoxLayout(self.frame_6)
        self.horizontalLayout_10.setSpacing(0)
        self.horizontalLayout_10.setObjectName(u"horizontalLayout_10")
        self.horizontalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.radioButton_2 = QRadioButton(self.frame_6)
        self.radioButton_2.setObjectName(u"radioButton_2")
        icon1 = QIcon()
        icon1.addFile(u":/Locks/Locks/GreenUnlock.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.radioButton_2.setIcon(icon1)
        self.radioButton_2.setIconSize(QSize(24, 24))
        self.radioButton_2.setChecked(True)

        self.horizontalLayout_10.addWidget(self.radioButton_2, 0, Qt.AlignBottom)

        self.radioButton = QRadioButton(self.frame_6)
        self.radioButton.setObjectName(u"radioButton")
        icon2 = QIcon()
        icon2.addFile(u":/Locks/Locks/RedLock.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.radioButton.setIcon(icon2)
        self.radioButton.setIconSize(QSize(24, 24))
        self.radioButton.setChecked(False)

        self.horizontalLayout_10.addWidget(self.radioButton)


        self.verticalLayout.addWidget(self.frame_6, 0, Qt.AlignLeft|Qt.AlignBottom)

        self.SolverTypeMenu = QTabWidget(self.centralwidget)
        self.SolverTypeMenu.setObjectName(u"SolverTypeMenu")
        self.SolverTypeMenu.setEnabled(True)
        font = QFont()
        font.setBold(True)
        self.SolverTypeMenu.setFont(font)
        self.SolverTypeMenu.setContextMenuPolicy(Qt.DefaultContextMenu)
        self.SolverTypeMenu.setTabShape(QTabWidget.Rounded)
        self.ChooseSolver = QWidget()
        self.ChooseSolver.setObjectName(u"ChooseSolver")
        self.ChooseSolver.setEnabled(True)
        self.verticalLayout_3 = QVBoxLayout(self.ChooseSolver)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 5)
        self.label_2 = QLabel(self.ChooseSolver)
        self.label_2.setObjectName(u"label_2")

        self.verticalLayout_3.addWidget(self.label_2)

        self.SolverTypeMenu.addTab(self.ChooseSolver, "")
        self.ExplicitSolver = QWidget()
        self.ExplicitSolver.setObjectName(u"ExplicitSolver")
        self.horizontalLayout_31 = QHBoxLayout(self.ExplicitSolver)
        self.horizontalLayout_31.setObjectName(u"horizontalLayout_31")
        self.horizontalLayout_31.setContentsMargins(0, 0, 0, 0)
        self.ExplicitPipelines = QTabWidget(self.ExplicitSolver)
        self.ExplicitPipelines.setObjectName(u"ExplicitPipelines")
        self.ExplicitPipelines.setMaximumSize(QSize(16777215, 125))
        font1 = QFont()
        font1.setBold(False)
        self.ExplicitPipelines.setFont(font1)
        self.Mesh_3 = QWidget()
        self.Mesh_3.setObjectName(u"Mesh_3")
        self.horizontalLayout_25 = QHBoxLayout(self.Mesh_3)
        self.horizontalLayout_25.setObjectName(u"horizontalLayout_25")
        self.horizontalLayout_25.setContentsMargins(5, 0, 5, 0)
        self.PartTools_4 = QFrame(self.Mesh_3)
        self.PartTools_4.setObjectName(u"PartTools_4")
        self.PartTools_4.setFrameShape(QFrame.NoFrame)
        self.PartTools_4.setFrameShadow(QFrame.Raised)
        self.verticalLayout_55 = QVBoxLayout(self.PartTools_4)
        self.verticalLayout_55.setSpacing(0)
        self.verticalLayout_55.setObjectName(u"verticalLayout_55")
        self.verticalLayout_55.setContentsMargins(0, 0, 0, 0)
        self.frame_5 = QFrame(self.PartTools_4)
        self.frame_5.setObjectName(u"frame_5")
        self.frame_5.setFrameShape(QFrame.NoFrame)
        self.frame_5.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_26 = QHBoxLayout(self.frame_5)
        self.horizontalLayout_26.setObjectName(u"horizontalLayout_26")
        self.horizontalLayout_26.setContentsMargins(0, 0, 0, 0)
        self.GlobaMesh = QFrame(self.frame_5)
        self.GlobaMesh.setObjectName(u"GlobaMesh")
        self.GlobaMesh.setMinimumSize(QSize(70, 80))
        self.GlobaMesh.setMaximumSize(QSize(70, 80))
        self.GlobaMesh.setFrameShape(QFrame.NoFrame)
        self.GlobaMesh.setFrameShadow(QFrame.Raised)
        self.verticalLayout_34 = QVBoxLayout(self.GlobaMesh)
        self.verticalLayout_34.setSpacing(0)
        self.verticalLayout_34.setObjectName(u"verticalLayout_34")
        self.verticalLayout_34.setContentsMargins(0, 0, 0, 5)
        self.BGlobalMesh = QPushButton(self.GlobaMesh)
        self.BGlobalMesh.setObjectName(u"BGlobalMesh")
        icon3 = QIcon()
        icon3.addFile(u":/Blue Icons/Blue Icons/mesh.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BGlobalMesh.setIcon(icon3)
        self.BGlobalMesh.setIconSize(QSize(40, 40))
        self.BGlobalMesh.setFlat(True)

        self.verticalLayout_34.addWidget(self.BGlobalMesh)

        self.LGlobalMesh = QLabel(self.GlobaMesh)
        self.LGlobalMesh.setObjectName(u"LGlobalMesh")
        self.LGlobalMesh.setWordWrap(True)

        self.verticalLayout_34.addWidget(self.LGlobalMesh)


        self.horizontalLayout_26.addWidget(self.GlobaMesh)

        self.ImportPart_4 = QFrame(self.frame_5)
        self.ImportPart_4.setObjectName(u"ImportPart_4")
        self.ImportPart_4.setMinimumSize(QSize(70, 80))
        self.ImportPart_4.setMaximumSize(QSize(70, 80))
        self.ImportPart_4.setFrameShape(QFrame.NoFrame)
        self.ImportPart_4.setFrameShadow(QFrame.Raised)
        self.verticalLayout_56 = QVBoxLayout(self.ImportPart_4)
        self.verticalLayout_56.setSpacing(0)
        self.verticalLayout_56.setObjectName(u"verticalLayout_56")
        self.verticalLayout_56.setContentsMargins(0, 0, 0, 5)
        self.BImportPart_4 = QPushButton(self.ImportPart_4)
        self.BImportPart_4.setObjectName(u"BImportPart_4")
        icon4 = QIcon()
        icon4.addFile(u":/Blue Icons/Blue Icons/Cube.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BImportPart_4.setIcon(icon4)
        self.BImportPart_4.setIconSize(QSize(50, 50))
        self.BImportPart_4.setFlat(True)

        self.verticalLayout_56.addWidget(self.BImportPart_4)

        self.LImportPart_4 = QLabel(self.ImportPart_4)
        self.LImportPart_4.setObjectName(u"LImportPart_4")
        self.LImportPart_4.setWordWrap(True)

        self.verticalLayout_56.addWidget(self.LImportPart_4)


        self.horizontalLayout_26.addWidget(self.ImportPart_4)


        self.verticalLayout_55.addWidget(self.frame_5, 0, Qt.AlignLeft)

        self.LinePartTools_4 = QFrame(self.PartTools_4)
        self.LinePartTools_4.setObjectName(u"LinePartTools_4")
        self.LinePartTools_4.setFrameShape(QFrame.HLine)
        self.LinePartTools_4.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_55.addWidget(self.LinePartTools_4)

        self.LPartTools_4 = QLabel(self.PartTools_4)
        self.LPartTools_4.setObjectName(u"LPartTools_4")

        self.verticalLayout_55.addWidget(self.LPartTools_4)


        self.horizontalLayout_25.addWidget(self.PartTools_4)

        self.MaterialTools_4 = QFrame(self.Mesh_3)
        self.MaterialTools_4.setObjectName(u"MaterialTools_4")
        self.MaterialTools_4.setFrameShape(QFrame.NoFrame)
        self.MaterialTools_4.setFrameShadow(QFrame.Raised)
        self.verticalLayout_57 = QVBoxLayout(self.MaterialTools_4)
        self.verticalLayout_57.setSpacing(0)
        self.verticalLayout_57.setObjectName(u"verticalLayout_57")
        self.verticalLayout_57.setContentsMargins(0, 0, 0, 0)
        self.MaterialToolsBtns_4 = QFrame(self.MaterialTools_4)
        self.MaterialToolsBtns_4.setObjectName(u"MaterialToolsBtns_4")
        self.MaterialToolsBtns_4.setFrameShape(QFrame.NoFrame)
        self.MaterialToolsBtns_4.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_27 = QHBoxLayout(self.MaterialToolsBtns_4)
        self.horizontalLayout_27.setSpacing(0)
        self.horizontalLayout_27.setObjectName(u"horizontalLayout_27")
        self.horizontalLayout_27.setContentsMargins(0, 0, 0, 0)
        self.DefineMaterial_4 = QFrame(self.MaterialToolsBtns_4)
        self.DefineMaterial_4.setObjectName(u"DefineMaterial_4")
        self.DefineMaterial_4.setMinimumSize(QSize(80, 80))
        self.DefineMaterial_4.setMaximumSize(QSize(80, 80))
        self.DefineMaterial_4.setFrameShape(QFrame.NoFrame)
        self.DefineMaterial_4.setFrameShadow(QFrame.Raised)
        self.verticalLayout_58 = QVBoxLayout(self.DefineMaterial_4)
        self.verticalLayout_58.setSpacing(0)
        self.verticalLayout_58.setObjectName(u"verticalLayout_58")
        self.verticalLayout_58.setContentsMargins(0, 0, 0, 5)
        self.BDefineMaterial_4 = QPushButton(self.DefineMaterial_4)
        self.BDefineMaterial_4.setObjectName(u"BDefineMaterial_4")
        icon5 = QIcon()
        icon5.addFile(u":/Blue Icons/Blue Icons/mine.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BDefineMaterial_4.setIcon(icon5)
        self.BDefineMaterial_4.setIconSize(QSize(50, 50))
        self.BDefineMaterial_4.setFlat(True)

        self.verticalLayout_58.addWidget(self.BDefineMaterial_4)

        self.LDefineMaterial_4 = QLabel(self.DefineMaterial_4)
        self.LDefineMaterial_4.setObjectName(u"LDefineMaterial_4")
        self.LDefineMaterial_4.setWordWrap(True)

        self.verticalLayout_58.addWidget(self.LDefineMaterial_4, 0, Qt.AlignBottom)


        self.horizontalLayout_27.addWidget(self.DefineMaterial_4, 0, Qt.AlignLeft)


        self.verticalLayout_57.addWidget(self.MaterialToolsBtns_4, 0, Qt.AlignLeft)

        self.LineMaterialTools_4 = QFrame(self.MaterialTools_4)
        self.LineMaterialTools_4.setObjectName(u"LineMaterialTools_4")
        self.LineMaterialTools_4.setFrameShape(QFrame.HLine)
        self.LineMaterialTools_4.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_57.addWidget(self.LineMaterialTools_4)

        self.LMaterialTools_4 = QLabel(self.MaterialTools_4)
        self.LMaterialTools_4.setObjectName(u"LMaterialTools_4")

        self.verticalLayout_57.addWidget(self.LMaterialTools_4)


        self.horizontalLayout_25.addWidget(self.MaterialTools_4)

        self.BCTools_4 = QFrame(self.Mesh_3)
        self.BCTools_4.setObjectName(u"BCTools_4")
        self.BCTools_4.setFrameShape(QFrame.NoFrame)
        self.BCTools_4.setFrameShadow(QFrame.Raised)
        self.verticalLayout_59 = QVBoxLayout(self.BCTools_4)
        self.verticalLayout_59.setSpacing(0)
        self.verticalLayout_59.setObjectName(u"verticalLayout_59")
        self.verticalLayout_59.setContentsMargins(0, 0, 0, 0)
        self.ApplyBC_4 = QFrame(self.BCTools_4)
        self.ApplyBC_4.setObjectName(u"ApplyBC_4")
        self.ApplyBC_4.setMinimumSize(QSize(65, 80))
        self.ApplyBC_4.setMaximumSize(QSize(70, 80))
        self.ApplyBC_4.setFrameShape(QFrame.NoFrame)
        self.ApplyBC_4.setFrameShadow(QFrame.Raised)
        self.verticalLayout_60 = QVBoxLayout(self.ApplyBC_4)
        self.verticalLayout_60.setSpacing(0)
        self.verticalLayout_60.setObjectName(u"verticalLayout_60")
        self.verticalLayout_60.setContentsMargins(0, 0, 0, 5)
        self.BApplyBC_4 = QPushButton(self.ApplyBC_4)
        self.BApplyBC_4.setObjectName(u"BApplyBC_4")
        icon6 = QIcon()
        icon6.addFile(u":/Blue Icons/Blue Icons/brick.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BApplyBC_4.setIcon(icon6)
        self.BApplyBC_4.setIconSize(QSize(50, 50))
        self.BApplyBC_4.setFlat(True)

        self.verticalLayout_60.addWidget(self.BApplyBC_4)

        self.LApplyBC_4 = QLabel(self.ApplyBC_4)
        self.LApplyBC_4.setObjectName(u"LApplyBC_4")
        self.LApplyBC_4.setWordWrap(True)

        self.verticalLayout_60.addWidget(self.LApplyBC_4)


        self.verticalLayout_59.addWidget(self.ApplyBC_4, 0, Qt.AlignLeft)

        self.LineBCTools_4 = QFrame(self.BCTools_4)
        self.LineBCTools_4.setObjectName(u"LineBCTools_4")
        self.LineBCTools_4.setFrameShape(QFrame.HLine)
        self.LineBCTools_4.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_59.addWidget(self.LineBCTools_4)

        self.LBCTools_4 = QLabel(self.BCTools_4)
        self.LBCTools_4.setObjectName(u"LBCTools_4")

        self.verticalLayout_59.addWidget(self.LBCTools_4)


        self.horizontalLayout_25.addWidget(self.BCTools_4)

        self.JobTools_6 = QFrame(self.Mesh_3)
        self.JobTools_6.setObjectName(u"JobTools_6")
        self.JobTools_6.setFrameShape(QFrame.NoFrame)
        self.JobTools_6.setFrameShadow(QFrame.Raised)
        self.verticalLayout_63 = QVBoxLayout(self.JobTools_6)
        self.verticalLayout_63.setSpacing(0)
        self.verticalLayout_63.setObjectName(u"verticalLayout_63")
        self.verticalLayout_63.setContentsMargins(0, 0, 0, 0)
        self.JobToolsBtns_6 = QFrame(self.JobTools_6)
        self.JobToolsBtns_6.setObjectName(u"JobToolsBtns_6")
        self.JobToolsBtns_6.setMinimumSize(QSize(0, 0))
        self.JobToolsBtns_6.setMaximumSize(QSize(16777215, 16777215))
        self.JobToolsBtns_6.setFrameShape(QFrame.NoFrame)
        self.JobToolsBtns_6.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_30 = QHBoxLayout(self.JobToolsBtns_6)
        self.horizontalLayout_30.setSpacing(0)
        self.horizontalLayout_30.setObjectName(u"horizontalLayout_30")
        self.horizontalLayout_30.setContentsMargins(0, 0, 0, 0)
        self.SolverSettings_6 = QFrame(self.JobToolsBtns_6)
        self.SolverSettings_6.setObjectName(u"SolverSettings_6")
        self.SolverSettings_6.setMinimumSize(QSize(80, 80))
        self.SolverSettings_6.setMaximumSize(QSize(80, 80))
        self.SolverSettings_6.setFrameShape(QFrame.NoFrame)
        self.SolverSettings_6.setFrameShadow(QFrame.Raised)
        self.verticalLayout_64 = QVBoxLayout(self.SolverSettings_6)
        self.verticalLayout_64.setSpacing(0)
        self.verticalLayout_64.setObjectName(u"verticalLayout_64")
        self.verticalLayout_64.setContentsMargins(0, 0, 0, 5)
        self.BSolverSettings_6 = QPushButton(self.SolverSettings_6)
        self.BSolverSettings_6.setObjectName(u"BSolverSettings_6")
        icon7 = QIcon()
        icon7.addFile(u":/Blue Icons/Blue Icons/gear.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BSolverSettings_6.setIcon(icon7)
        self.BSolverSettings_6.setIconSize(QSize(50, 50))
        self.BSolverSettings_6.setFlat(True)

        self.verticalLayout_64.addWidget(self.BSolverSettings_6)

        self.LSolverSettings_7 = QLabel(self.SolverSettings_6)
        self.LSolverSettings_7.setObjectName(u"LSolverSettings_7")
        self.LSolverSettings_7.setWordWrap(True)

        self.verticalLayout_64.addWidget(self.LSolverSettings_7)


        self.horizontalLayout_30.addWidget(self.SolverSettings_6, 0, Qt.AlignLeft)

        self.RunEVPFFT_6 = QFrame(self.JobToolsBtns_6)
        self.RunEVPFFT_6.setObjectName(u"RunEVPFFT_6")
        self.RunEVPFFT_6.setMinimumSize(QSize(75, 80))
        self.RunEVPFFT_6.setMaximumSize(QSize(65, 80))
        self.RunEVPFFT_6.setFrameShape(QFrame.NoFrame)
        self.RunEVPFFT_6.setFrameShadow(QFrame.Raised)
        self.verticalLayout_65 = QVBoxLayout(self.RunEVPFFT_6)
        self.verticalLayout_65.setSpacing(0)
        self.verticalLayout_65.setObjectName(u"verticalLayout_65")
        self.verticalLayout_65.setContentsMargins(0, 0, 0, 5)
        self.BRunEVPFFT_6 = QPushButton(self.RunEVPFFT_6)
        self.BRunEVPFFT_6.setObjectName(u"BRunEVPFFT_6")
        icon8 = QIcon()
        icon8.addFile(u":/Blue Icons/Blue Icons/Play.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BRunEVPFFT_6.setIcon(icon8)
        self.BRunEVPFFT_6.setIconSize(QSize(40, 40))
        self.BRunEVPFFT_6.setFlat(True)

        self.verticalLayout_65.addWidget(self.BRunEVPFFT_6)

        self.LRunEVPFFT_6 = QLabel(self.RunEVPFFT_6)
        self.LRunEVPFFT_6.setObjectName(u"LRunEVPFFT_6")
        self.LRunEVPFFT_6.setWordWrap(True)

        self.verticalLayout_65.addWidget(self.LRunEVPFFT_6)


        self.horizontalLayout_30.addWidget(self.RunEVPFFT_6)

        self.ViewResults_3 = QFrame(self.JobToolsBtns_6)
        self.ViewResults_3.setObjectName(u"ViewResults_3")
        self.ViewResults_3.setMinimumSize(QSize(80, 80))
        self.ViewResults_3.setMaximumSize(QSize(80, 80))
        self.ViewResults_3.setFrameShape(QFrame.NoFrame)
        self.ViewResults_3.setFrameShadow(QFrame.Raised)
        self.verticalLayout_66 = QVBoxLayout(self.ViewResults_3)
        self.verticalLayout_66.setSpacing(0)
        self.verticalLayout_66.setObjectName(u"verticalLayout_66")
        self.verticalLayout_66.setContentsMargins(0, 0, 0, 5)
        self.BViewResults_3 = QPushButton(self.ViewResults_3)
        self.BViewResults_3.setObjectName(u"BViewResults_3")
        icon9 = QIcon()
        icon9.addFile(u":/Blue Icons/Blue Icons/magnify.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BViewResults_3.setIcon(icon9)
        self.BViewResults_3.setIconSize(QSize(40, 40))
        self.BViewResults_3.setFlat(True)

        self.verticalLayout_66.addWidget(self.BViewResults_3)

        self.LViewResults_3 = QLabel(self.ViewResults_3)
        self.LViewResults_3.setObjectName(u"LViewResults_3")
        self.LViewResults_3.setWordWrap(True)

        self.verticalLayout_66.addWidget(self.LViewResults_3)


        self.horizontalLayout_30.addWidget(self.ViewResults_3)


        self.verticalLayout_63.addWidget(self.JobToolsBtns_6, 0, Qt.AlignLeft)

        self.LineJobTools_6 = QFrame(self.JobTools_6)
        self.LineJobTools_6.setObjectName(u"LineJobTools_6")
        self.LineJobTools_6.setFrameShape(QFrame.HLine)
        self.LineJobTools_6.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_63.addWidget(self.LineJobTools_6)

        self.LJobTools_6 = QLabel(self.JobTools_6)
        self.LJobTools_6.setObjectName(u"LJobTools_6")

        self.verticalLayout_63.addWidget(self.LJobTools_6)


        self.horizontalLayout_25.addWidget(self.JobTools_6)

        self.ExplicitPipelines.addTab(self.Mesh_3, "")

        self.horizontalLayout_31.addWidget(self.ExplicitPipelines)

        self.SolverTypeMenu.addTab(self.ExplicitSolver, "")
        self.ImplicitSolver = QWidget()
        self.ImplicitSolver.setObjectName(u"ImplicitSolver")
        self.verticalLayout_7 = QVBoxLayout(self.ImplicitSolver)
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.label_4 = QLabel(self.ImplicitSolver)
        self.label_4.setObjectName(u"label_4")

        self.verticalLayout_7.addWidget(self.label_4)

        self.SolverTypeMenu.addTab(self.ImplicitSolver, "")
        self.EVPFFTSolver = QWidget()
        self.EVPFFTSolver.setObjectName(u"EVPFFTSolver")
        self.verticalLayout_5 = QVBoxLayout(self.EVPFFTSolver)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.EVPFFTPipelines = QTabWidget(self.EVPFFTSolver)
        self.EVPFFTPipelines.setObjectName(u"EVPFFTPipelines")
        self.EVPFFTPipelines.setMaximumSize(QSize(16777215, 150))
        self.EVPFFTPipelines.setFont(font1)
        self.ChoosePipeline = QWidget()
        self.ChoosePipeline.setObjectName(u"ChoosePipeline")
        self.verticalLayout_4 = QVBoxLayout(self.ChoosePipeline)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.label_3 = QLabel(self.ChoosePipeline)
        self.label_3.setObjectName(u"label_3")

        self.verticalLayout_4.addWidget(self.label_3)

        self.EVPFFTPipelines.addTab(self.ChoosePipeline, "")
        self.EVPFFTGeneral = QWidget()
        self.EVPFFTGeneral.setObjectName(u"EVPFFTGeneral")
        font2 = QFont()
        font2.setPointSize(13)
        font2.setBold(False)
        self.EVPFFTGeneral.setFont(font2)
        self.horizontalLayout_2 = QHBoxLayout(self.EVPFFTGeneral)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.PartTools_5 = QFrame(self.EVPFFTGeneral)
        self.PartTools_5.setObjectName(u"PartTools_5")
        self.PartTools_5.setFrameShape(QFrame.NoFrame)
        self.PartTools_5.setFrameShadow(QFrame.Raised)
        self.verticalLayout_68 = QVBoxLayout(self.PartTools_5)
        self.verticalLayout_68.setSpacing(0)
        self.verticalLayout_68.setObjectName(u"verticalLayout_68")
        self.verticalLayout_68.setContentsMargins(0, 0, 0, 0)
        self.ImportPart_5 = QFrame(self.PartTools_5)
        self.ImportPart_5.setObjectName(u"ImportPart_5")
        self.ImportPart_5.setMinimumSize(QSize(70, 80))
        self.ImportPart_5.setMaximumSize(QSize(70, 80))
        self.ImportPart_5.setFrameShape(QFrame.NoFrame)
        self.ImportPart_5.setFrameShadow(QFrame.Raised)
        self.verticalLayout_69 = QVBoxLayout(self.ImportPart_5)
        self.verticalLayout_69.setSpacing(0)
        self.verticalLayout_69.setObjectName(u"verticalLayout_69")
        self.verticalLayout_69.setContentsMargins(0, 0, 0, 5)
        self.BImportPart_5 = QPushButton(self.ImportPart_5)
        self.BImportPart_5.setObjectName(u"BImportPart_5")
        self.BImportPart_5.setIcon(icon4)
        self.BImportPart_5.setIconSize(QSize(50, 50))
        self.BImportPart_5.setFlat(True)

        self.verticalLayout_69.addWidget(self.BImportPart_5)

        self.LImportPart_5 = QLabel(self.ImportPart_5)
        self.LImportPart_5.setObjectName(u"LImportPart_5")
        self.LImportPart_5.setWordWrap(True)

        self.verticalLayout_69.addWidget(self.LImportPart_5)


        self.verticalLayout_68.addWidget(self.ImportPart_5)

        self.LinePartTools_5 = QFrame(self.PartTools_5)
        self.LinePartTools_5.setObjectName(u"LinePartTools_5")
        self.LinePartTools_5.setFrameShape(QFrame.HLine)
        self.LinePartTools_5.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_68.addWidget(self.LinePartTools_5)

        self.LPartTools_5 = QLabel(self.PartTools_5)
        self.LPartTools_5.setObjectName(u"LPartTools_5")

        self.verticalLayout_68.addWidget(self.LPartTools_5)


        self.horizontalLayout_2.addWidget(self.PartTools_5, 0, Qt.AlignLeft)

        self.label_7 = QLabel(self.EVPFFTGeneral)
        self.label_7.setObjectName(u"label_7")

        self.horizontalLayout_2.addWidget(self.label_7)

        self.EVPFFTPipelines.addTab(self.EVPFFTGeneral, "")
        self.EVPFFTHomogenization = QWidget()
        self.EVPFFTHomogenization.setObjectName(u"EVPFFTHomogenization")
        self.verticalLayout_6 = QVBoxLayout(self.EVPFFTHomogenization)
        self.verticalLayout_6.setSpacing(0)
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.verticalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.frame_8 = QFrame(self.EVPFFTHomogenization)
        self.frame_8.setObjectName(u"frame_8")
        self.frame_8.setFrameShape(QFrame.NoFrame)
        self.frame_8.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_11 = QHBoxLayout(self.frame_8)
        self.horizontalLayout_11.setObjectName(u"horizontalLayout_11")
        self.horizontalLayout_11.setContentsMargins(0, 0, 0, 0)
        self.PartTools = QFrame(self.frame_8)
        self.PartTools.setObjectName(u"PartTools")
        self.PartTools.setFrameShape(QFrame.NoFrame)
        self.PartTools.setFrameShadow(QFrame.Raised)
        self.verticalLayout_49 = QVBoxLayout(self.PartTools)
        self.verticalLayout_49.setSpacing(0)
        self.verticalLayout_49.setObjectName(u"verticalLayout_49")
        self.verticalLayout_49.setContentsMargins(0, 0, 0, 0)
        self.ImportPart = QFrame(self.PartTools)
        self.ImportPart.setObjectName(u"ImportPart")
        self.ImportPart.setMinimumSize(QSize(70, 80))
        self.ImportPart.setMaximumSize(QSize(70, 80))
        self.ImportPart.setFrameShape(QFrame.NoFrame)
        self.ImportPart.setFrameShadow(QFrame.Raised)
        self.verticalLayout_50 = QVBoxLayout(self.ImportPart)
        self.verticalLayout_50.setSpacing(0)
        self.verticalLayout_50.setObjectName(u"verticalLayout_50")
        self.verticalLayout_50.setContentsMargins(0, 0, 0, 5)
        self.BImportPart = QPushButton(self.ImportPart)
        self.BImportPart.setObjectName(u"BImportPart")
        self.BImportPart.setIcon(icon4)
        self.BImportPart.setIconSize(QSize(50, 50))
        self.BImportPart.setFlat(True)

        self.verticalLayout_50.addWidget(self.BImportPart)

        self.LImportPart = QLabel(self.ImportPart)
        self.LImportPart.setObjectName(u"LImportPart")
        self.LImportPart.setWordWrap(True)

        self.verticalLayout_50.addWidget(self.LImportPart)


        self.verticalLayout_49.addWidget(self.ImportPart)

        self.LinePartTools = QFrame(self.PartTools)
        self.LinePartTools.setObjectName(u"LinePartTools")
        self.LinePartTools.setFrameShape(QFrame.HLine)
        self.LinePartTools.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_49.addWidget(self.LinePartTools)

        self.LPartTools = QLabel(self.PartTools)
        self.LPartTools.setObjectName(u"LPartTools")

        self.verticalLayout_49.addWidget(self.LPartTools)


        self.horizontalLayout_11.addWidget(self.PartTools)

        self.MaterialTools = QFrame(self.frame_8)
        self.MaterialTools.setObjectName(u"MaterialTools")
        self.MaterialTools.setFrameShape(QFrame.NoFrame)
        self.MaterialTools.setFrameShadow(QFrame.Raised)
        self.verticalLayout_53 = QVBoxLayout(self.MaterialTools)
        self.verticalLayout_53.setSpacing(0)
        self.verticalLayout_53.setObjectName(u"verticalLayout_53")
        self.verticalLayout_53.setContentsMargins(0, 0, 0, 0)
        self.MaterialToolsBtns = QFrame(self.MaterialTools)
        self.MaterialToolsBtns.setObjectName(u"MaterialToolsBtns")
        self.MaterialToolsBtns.setFrameShape(QFrame.NoFrame)
        self.MaterialToolsBtns.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_22 = QHBoxLayout(self.MaterialToolsBtns)
        self.horizontalLayout_22.setSpacing(0)
        self.horizontalLayout_22.setObjectName(u"horizontalLayout_22")
        self.horizontalLayout_22.setContentsMargins(0, 0, 0, 0)
        self.DefineMaterial = QFrame(self.MaterialToolsBtns)
        self.DefineMaterial.setObjectName(u"DefineMaterial")
        self.DefineMaterial.setMinimumSize(QSize(80, 80))
        self.DefineMaterial.setMaximumSize(QSize(80, 80))
        self.DefineMaterial.setFrameShape(QFrame.NoFrame)
        self.DefineMaterial.setFrameShadow(QFrame.Raised)
        self.verticalLayout_54 = QVBoxLayout(self.DefineMaterial)
        self.verticalLayout_54.setSpacing(0)
        self.verticalLayout_54.setObjectName(u"verticalLayout_54")
        self.verticalLayout_54.setContentsMargins(0, 0, 0, 5)
        self.BDefineMaterial = QPushButton(self.DefineMaterial)
        self.BDefineMaterial.setObjectName(u"BDefineMaterial")
        self.BDefineMaterial.setEnabled(True)
        self.BDefineMaterial.setIcon(icon5)
        self.BDefineMaterial.setIconSize(QSize(50, 50))
        self.BDefineMaterial.setFlat(True)

        self.verticalLayout_54.addWidget(self.BDefineMaterial)

        self.LDefineMaterial = QLabel(self.DefineMaterial)
        self.LDefineMaterial.setObjectName(u"LDefineMaterial")
        self.LDefineMaterial.setWordWrap(True)

        self.verticalLayout_54.addWidget(self.LDefineMaterial, 0, Qt.AlignBottom)


        self.horizontalLayout_22.addWidget(self.DefineMaterial, 0, Qt.AlignLeft)


        self.verticalLayout_53.addWidget(self.MaterialToolsBtns, 0, Qt.AlignLeft)

        self.LineMaterialTools = QFrame(self.MaterialTools)
        self.LineMaterialTools.setObjectName(u"LineMaterialTools")
        self.LineMaterialTools.setFrameShape(QFrame.HLine)
        self.LineMaterialTools.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_53.addWidget(self.LineMaterialTools)

        self.LMaterialTools = QLabel(self.MaterialTools)
        self.LMaterialTools.setObjectName(u"LMaterialTools")

        self.verticalLayout_53.addWidget(self.LMaterialTools)


        self.horizontalLayout_11.addWidget(self.MaterialTools)

        self.BCTools = QFrame(self.frame_8)
        self.BCTools.setObjectName(u"BCTools")
        self.BCTools.setFrameShape(QFrame.NoFrame)
        self.BCTools.setFrameShadow(QFrame.Raised)
        self.verticalLayout_51 = QVBoxLayout(self.BCTools)
        self.verticalLayout_51.setSpacing(0)
        self.verticalLayout_51.setObjectName(u"verticalLayout_51")
        self.verticalLayout_51.setContentsMargins(0, 0, 0, 0)
        self.ApplyBC = QFrame(self.BCTools)
        self.ApplyBC.setObjectName(u"ApplyBC")
        self.ApplyBC.setMinimumSize(QSize(65, 80))
        self.ApplyBC.setMaximumSize(QSize(70, 80))
        self.ApplyBC.setFrameShape(QFrame.NoFrame)
        self.ApplyBC.setFrameShadow(QFrame.Raised)
        self.verticalLayout_52 = QVBoxLayout(self.ApplyBC)
        self.verticalLayout_52.setSpacing(0)
        self.verticalLayout_52.setObjectName(u"verticalLayout_52")
        self.verticalLayout_52.setContentsMargins(0, 0, 0, 5)
        self.BApplyBC = QPushButton(self.ApplyBC)
        self.BApplyBC.setObjectName(u"BApplyBC")
        self.BApplyBC.setIcon(icon6)
        self.BApplyBC.setIconSize(QSize(50, 50))
        self.BApplyBC.setFlat(True)

        self.verticalLayout_52.addWidget(self.BApplyBC)

        self.LApplyBC = QLabel(self.ApplyBC)
        self.LApplyBC.setObjectName(u"LApplyBC")
        self.LApplyBC.setWordWrap(True)

        self.verticalLayout_52.addWidget(self.LApplyBC)


        self.verticalLayout_51.addWidget(self.ApplyBC, 0, Qt.AlignLeft)

        self.LineBCTools = QFrame(self.BCTools)
        self.LineBCTools.setObjectName(u"LineBCTools")
        self.LineBCTools.setFrameShape(QFrame.HLine)
        self.LineBCTools.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_51.addWidget(self.LineBCTools)

        self.LBCTools = QLabel(self.BCTools)
        self.LBCTools.setObjectName(u"LBCTools")

        self.verticalLayout_51.addWidget(self.LBCTools)


        self.horizontalLayout_11.addWidget(self.BCTools)

        self.JobTools = QFrame(self.frame_8)
        self.JobTools.setObjectName(u"JobTools")
        self.JobTools.setFrameShape(QFrame.NoFrame)
        self.JobTools.setFrameShadow(QFrame.Raised)
        self.verticalLayout_30 = QVBoxLayout(self.JobTools)
        self.verticalLayout_30.setSpacing(0)
        self.verticalLayout_30.setObjectName(u"verticalLayout_30")
        self.verticalLayout_30.setContentsMargins(0, 0, 0, 0)
        self.JobToolsBtns = QFrame(self.JobTools)
        self.JobToolsBtns.setObjectName(u"JobToolsBtns")
        self.JobToolsBtns.setMinimumSize(QSize(0, 0))
        self.JobToolsBtns.setMaximumSize(QSize(16777215, 16777215))
        self.JobToolsBtns.setFrameShape(QFrame.NoFrame)
        self.JobToolsBtns.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_4 = QHBoxLayout(self.JobToolsBtns)
        self.horizontalLayout_4.setSpacing(0)
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.SolverSettings = QFrame(self.JobToolsBtns)
        self.SolverSettings.setObjectName(u"SolverSettings")
        self.SolverSettings.setMinimumSize(QSize(80, 80))
        self.SolverSettings.setMaximumSize(QSize(80, 80))
        self.SolverSettings.setFrameShape(QFrame.NoFrame)
        self.SolverSettings.setFrameShadow(QFrame.Raised)
        self.verticalLayout_8 = QVBoxLayout(self.SolverSettings)
        self.verticalLayout_8.setSpacing(0)
        self.verticalLayout_8.setObjectName(u"verticalLayout_8")
        self.verticalLayout_8.setContentsMargins(0, 0, 0, 5)
        self.BSolverSettings = QPushButton(self.SolverSettings)
        self.BSolverSettings.setObjectName(u"BSolverSettings")
        self.BSolverSettings.setIcon(icon7)
        self.BSolverSettings.setIconSize(QSize(50, 50))
        self.BSolverSettings.setFlat(True)

        self.verticalLayout_8.addWidget(self.BSolverSettings)

        self.LSolverSettings = QLabel(self.SolverSettings)
        self.LSolverSettings.setObjectName(u"LSolverSettings")
        self.LSolverSettings.setWordWrap(True)

        self.verticalLayout_8.addWidget(self.LSolverSettings)


        self.horizontalLayout_4.addWidget(self.SolverSettings, 0, Qt.AlignLeft)

        self.RunEVPFFT = QFrame(self.JobToolsBtns)
        self.RunEVPFFT.setObjectName(u"RunEVPFFT")
        self.RunEVPFFT.setMinimumSize(QSize(75, 80))
        self.RunEVPFFT.setMaximumSize(QSize(65, 80))
        self.RunEVPFFT.setFrameShape(QFrame.NoFrame)
        self.RunEVPFFT.setFrameShadow(QFrame.Raised)
        self.verticalLayout_67 = QVBoxLayout(self.RunEVPFFT)
        self.verticalLayout_67.setSpacing(0)
        self.verticalLayout_67.setObjectName(u"verticalLayout_67")
        self.verticalLayout_67.setContentsMargins(0, 0, 0, 5)
        self.BRunEVPFFT = QPushButton(self.RunEVPFFT)
        self.BRunEVPFFT.setObjectName(u"BRunEVPFFT")
        self.BRunEVPFFT.setIcon(icon8)
        self.BRunEVPFFT.setIconSize(QSize(40, 40))
        self.BRunEVPFFT.setFlat(True)

        self.verticalLayout_67.addWidget(self.BRunEVPFFT)

        self.LRunEVPFFT = QLabel(self.RunEVPFFT)
        self.LRunEVPFFT.setObjectName(u"LRunEVPFFT")
        self.LRunEVPFFT.setWordWrap(True)

        self.verticalLayout_67.addWidget(self.LRunEVPFFT)


        self.horizontalLayout_4.addWidget(self.RunEVPFFT)

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
        self.BViewResults.setIcon(icon9)
        self.BViewResults.setIconSize(QSize(40, 40))
        self.BViewResults.setFlat(True)

        self.verticalLayout_14.addWidget(self.BViewResults)

        self.LViewResults = QLabel(self.ViewResults)
        self.LViewResults.setObjectName(u"LViewResults")
        self.LViewResults.setWordWrap(True)

        self.verticalLayout_14.addWidget(self.LViewResults)


        self.horizontalLayout_4.addWidget(self.ViewResults)


        self.verticalLayout_30.addWidget(self.JobToolsBtns)

        self.LineJobTools = QFrame(self.JobTools)
        self.LineJobTools.setObjectName(u"LineJobTools")
        self.LineJobTools.setFrameShape(QFrame.HLine)
        self.LineJobTools.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_30.addWidget(self.LineJobTools)

        self.LJobTools = QLabel(self.JobTools)
        self.LJobTools.setObjectName(u"LJobTools")

        self.verticalLayout_30.addWidget(self.LJobTools)


        self.horizontalLayout_11.addWidget(self.JobTools)


        self.verticalLayout_6.addWidget(self.frame_8)

        self.EVPFFTPipelines.addTab(self.EVPFFTHomogenization, "")

        self.verticalLayout_5.addWidget(self.EVPFFTPipelines)

        self.SolverTypeMenu.addTab(self.EVPFFTSolver, "")
        self.EVPFFTLSSolver = QWidget()
        self.EVPFFTLSSolver.setObjectName(u"EVPFFTLSSolver")
        self.verticalLayout_11 = QVBoxLayout(self.EVPFFTLSSolver)
        self.verticalLayout_11.setObjectName(u"verticalLayout_11")
        self.label_5 = QLabel(self.EVPFFTLSSolver)
        self.label_5.setObjectName(u"label_5")

        self.verticalLayout_11.addWidget(self.label_5)

        self.SolverTypeMenu.addTab(self.EVPFFTLSSolver, "")

        self.verticalLayout.addWidget(self.SolverTypeMenu, 0, Qt.AlignVCenter)

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
        sizePolicy1 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.ToolSettings.sizePolicy().hasHeightForWidth())
        self.ToolSettings.setSizePolicy(sizePolicy1)
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
        self.LosAlamosLogo.setPixmap(QPixmap(u":/Logos/Logos/LANL Logo Ultramarine.png"))
        self.LosAlamosLogo.setScaledContents(True)

        self.verticalLayout_2.addWidget(self.LosAlamosLogo)

        self.verticalSpacer_7 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_2.addItem(self.verticalSpacer_7)

        self.EVPFFTLogo = QLabel(self.TitlePage)
        self.EVPFFTLogo.setObjectName(u"EVPFFTLogo")
        self.EVPFFTLogo.setMinimumSize(QSize(225, 175))
        self.EVPFFTLogo.setMaximumSize(QSize(225, 175))
        self.EVPFFTLogo.setPixmap(QPixmap(u":/Logos/Logos/FIERRO.png"))
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
        font3 = QFont()
        font3.setPointSize(16)
        self.LAdditionalSoftware.setFont(font3)

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
        self.MatarLogo.setPixmap(QPixmap(u":/Logos/Logos/MATAR_logo2.png"))
        self.MatarLogo.setScaledContents(True)

        self.horizontalLayout_5.addWidget(self.MatarLogo)

        self.ParaviewLogo = QLabel(self.AdditionalSoftwareLogos)
        self.ParaviewLogo.setObjectName(u"ParaviewLogo")
        self.ParaviewLogo.setMaximumSize(QSize(130, 30))
        self.ParaviewLogo.setPixmap(QPixmap(u":/Logos/Logos/ParaView_logo.png"))
        self.ParaviewLogo.setScaledContents(True)

        self.horizontalLayout_5.addWidget(self.ParaviewLogo)


        self.verticalLayout_10.addWidget(self.AdditionalSoftwareLogos)


        self.verticalLayout_2.addWidget(self.AdditionalSoftware, 0, Qt.AlignTop)

        self.ToolSettings.addWidget(self.TitlePage)
        self.GlobalMesh = QWidget()
        self.GlobalMesh.setObjectName(u"GlobalMesh")
        self.verticalLayout_31 = QVBoxLayout(self.GlobalMesh)
        self.verticalLayout_31.setObjectName(u"verticalLayout_31")
        self.LDefineGlobalMesh = QLabel(self.GlobalMesh)
        self.LDefineGlobalMesh.setObjectName(u"LDefineGlobalMesh")

        self.verticalLayout_31.addWidget(self.LDefineGlobalMesh)

        self.MeshInputs = QFrame(self.GlobalMesh)
        self.MeshInputs.setObjectName(u"MeshInputs")
        self.MeshInputs.setFrameShape(QFrame.NoFrame)
        self.MeshInputs.setFrameShadow(QFrame.Raised)
        self.formLayout_9 = QFormLayout(self.MeshInputs)
        self.formLayout_9.setObjectName(u"formLayout_9")
        self.formLayout_9.setContentsMargins(-1, 0, -1, 0)
        self.LElementType = QLabel(self.MeshInputs)
        self.LElementType.setObjectName(u"LElementType")

        self.formLayout_9.setWidget(0, QFormLayout.LabelRole, self.LElementType)

        self.INElementType = QComboBox(self.MeshInputs)
        self.INElementType.addItem("")
        self.INElementType.addItem("")
        self.INElementType.addItem("")
        self.INElementType.setObjectName(u"INElementType")

        self.formLayout_9.setWidget(0, QFormLayout.FieldRole, self.INElementType)

        self.LCoordinateSystem = QLabel(self.MeshInputs)
        self.LCoordinateSystem.setObjectName(u"LCoordinateSystem")

        self.formLayout_9.setWidget(1, QFormLayout.LabelRole, self.LCoordinateSystem)

        self.INCoordinateSystem = QComboBox(self.MeshInputs)
        self.INCoordinateSystem.addItem("")
        self.INCoordinateSystem.addItem("")
        self.INCoordinateSystem.setObjectName(u"INCoordinateSystem")

        self.formLayout_9.setWidget(1, QFormLayout.FieldRole, self.INCoordinateSystem)

        self.LDimension = QLabel(self.MeshInputs)
        self.LDimension.setObjectName(u"LDimension")

        self.formLayout_9.setWidget(2, QFormLayout.LabelRole, self.LDimension)

        self.INDimension = QComboBox(self.MeshInputs)
        self.INDimension.addItem("")
        self.INDimension.addItem("")
        self.INDimension.setObjectName(u"INDimension")

        self.formLayout_9.setWidget(2, QFormLayout.FieldRole, self.INDimension)


        self.verticalLayout_31.addWidget(self.MeshInputs)

        self.MeshInputs2 = QStackedWidget(self.GlobalMesh)
        self.MeshInputs2.setObjectName(u"MeshInputs2")
        self.Rectangular2D = QWidget()
        self.Rectangular2D.setObjectName(u"Rectangular2D")
        self.gridLayout_3 = QGridLayout(self.Rectangular2D)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.INLengthXR2D = QLineEdit(self.Rectangular2D)
        self.INLengthXR2D.setObjectName(u"INLengthXR2D")

        self.gridLayout_3.addWidget(self.INLengthXR2D, 1, 2, 1, 1)

        self.LLengthR2D = QLabel(self.Rectangular2D)
        self.LLengthR2D.setObjectName(u"LLengthR2D")

        self.gridLayout_3.addWidget(self.LLengthR2D, 1, 0, 1, 1)

        self.INElementsYR2D = QLineEdit(self.Rectangular2D)
        self.INElementsYR2D.setObjectName(u"INElementsYR2D")

        self.gridLayout_3.addWidget(self.INElementsYR2D, 2, 4, 1, 1)

        self.label_14 = QLabel(self.Rectangular2D)
        self.label_14.setObjectName(u"label_14")

        self.gridLayout_3.addWidget(self.label_14, 2, 5, 1, 1)

        self.label_11 = QLabel(self.Rectangular2D)
        self.label_11.setObjectName(u"label_11")

        self.gridLayout_3.addWidget(self.label_11, 1, 5, 1, 1)

        self.label_10 = QLabel(self.Rectangular2D)
        self.label_10.setObjectName(u"label_10")

        self.gridLayout_3.addWidget(self.label_10, 1, 3, 1, 1)

        self.INLengthYR2D = QLineEdit(self.Rectangular2D)
        self.INLengthYR2D.setObjectName(u"INLengthYR2D")

        self.gridLayout_3.addWidget(self.INLengthYR2D, 1, 4, 1, 1)

        self.LElementsR2D = QLabel(self.Rectangular2D)
        self.LElementsR2D.setObjectName(u"LElementsR2D")

        self.gridLayout_3.addWidget(self.LElementsR2D, 2, 0, 1, 1)

        self.label_13 = QLabel(self.Rectangular2D)
        self.label_13.setObjectName(u"label_13")

        self.gridLayout_3.addWidget(self.label_13, 2, 3, 1, 1)

        self.INElementsXR2D = QLineEdit(self.Rectangular2D)
        self.INElementsXR2D.setObjectName(u"INElementsXR2D")

        self.gridLayout_3.addWidget(self.INElementsXR2D, 2, 2, 1, 1)

        self.INOriginXR2D = QLineEdit(self.Rectangular2D)
        self.INOriginXR2D.setObjectName(u"INOriginXR2D")

        self.gridLayout_3.addWidget(self.INOriginXR2D, 0, 2, 1, 1)

        self.INOriginYR2D = QLineEdit(self.Rectangular2D)
        self.INOriginYR2D.setObjectName(u"INOriginYR2D")

        self.gridLayout_3.addWidget(self.INOriginYR2D, 0, 4, 1, 1)

        self.label_6 = QLabel(self.Rectangular2D)
        self.label_6.setObjectName(u"label_6")

        self.gridLayout_3.addWidget(self.label_6, 0, 3, 1, 1)

        self.label_8 = QLabel(self.Rectangular2D)
        self.label_8.setObjectName(u"label_8")

        self.gridLayout_3.addWidget(self.label_8, 0, 5, 1, 1)

        self.LOriginR2D = QLabel(self.Rectangular2D)
        self.LOriginR2D.setObjectName(u"LOriginR2D")

        self.gridLayout_3.addWidget(self.LOriginR2D, 0, 0, 1, 1)

        self.label_15 = QLabel(self.Rectangular2D)
        self.label_15.setObjectName(u"label_15")

        self.gridLayout_3.addWidget(self.label_15, 0, 1, 1, 1)

        self.label_16 = QLabel(self.Rectangular2D)
        self.label_16.setObjectName(u"label_16")

        self.gridLayout_3.addWidget(self.label_16, 1, 1, 1, 1)

        self.label_17 = QLabel(self.Rectangular2D)
        self.label_17.setObjectName(u"label_17")

        self.gridLayout_3.addWidget(self.label_17, 2, 1, 1, 1)

        self.MeshInputs2.addWidget(self.Rectangular2D)
        self.Rectangular3D = QWidget()
        self.Rectangular3D.setObjectName(u"Rectangular3D")
        self.gridLayout_4 = QGridLayout(self.Rectangular3D)
        self.gridLayout_4.setObjectName(u"gridLayout_4")
        self.LLengthR3D = QLabel(self.Rectangular3D)
        self.LLengthR3D.setObjectName(u"LLengthR3D")

        self.gridLayout_4.addWidget(self.LLengthR3D, 1, 0, 1, 1)

        self.label_25 = QLabel(self.Rectangular3D)
        self.label_25.setObjectName(u"label_25")

        self.gridLayout_4.addWidget(self.label_25, 2, 3, 1, 1)

        self.LElementsR3D = QLabel(self.Rectangular3D)
        self.LElementsR3D.setObjectName(u"LElementsR3D")

        self.gridLayout_4.addWidget(self.LElementsR3D, 2, 0, 1, 1)

        self.label_23 = QLabel(self.Rectangular3D)
        self.label_23.setObjectName(u"label_23")

        self.gridLayout_4.addWidget(self.label_23, 2, 1, 1, 1)

        self.label_30 = QLabel(self.Rectangular3D)
        self.label_30.setObjectName(u"label_30")

        self.gridLayout_4.addWidget(self.label_30, 0, 5, 1, 1)

        self.label_20 = QLabel(self.Rectangular3D)
        self.label_20.setObjectName(u"label_20")

        self.gridLayout_4.addWidget(self.label_20, 1, 7, 1, 1)

        self.label_22 = QLabel(self.Rectangular3D)
        self.label_22.setObjectName(u"label_22")

        self.gridLayout_4.addWidget(self.label_22, 0, 1, 1, 1)

        self.INLengthYR3D = QLineEdit(self.Rectangular3D)
        self.INLengthYR3D.setObjectName(u"INLengthYR3D")

        self.gridLayout_4.addWidget(self.INLengthYR3D, 1, 4, 1, 1)

        self.INOriginYR3D = QLineEdit(self.Rectangular3D)
        self.INOriginYR3D.setObjectName(u"INOriginYR3D")

        self.gridLayout_4.addWidget(self.INOriginYR3D, 0, 4, 1, 1)

        self.LOriginR3D = QLabel(self.Rectangular3D)
        self.LOriginR3D.setObjectName(u"LOriginR3D")

        self.gridLayout_4.addWidget(self.LOriginR3D, 0, 0, 1, 1)

        self.INElementsYR3D = QLineEdit(self.Rectangular3D)
        self.INElementsYR3D.setObjectName(u"INElementsYR3D")

        self.gridLayout_4.addWidget(self.INElementsYR3D, 2, 4, 1, 1)

        self.label_26 = QLabel(self.Rectangular3D)
        self.label_26.setObjectName(u"label_26")

        self.gridLayout_4.addWidget(self.label_26, 1, 1, 1, 1)

        self.INLengthXR3D = QLineEdit(self.Rectangular3D)
        self.INLengthXR3D.setObjectName(u"INLengthXR3D")

        self.gridLayout_4.addWidget(self.INLengthXR3D, 1, 2, 1, 1)

        self.label_28 = QLabel(self.Rectangular3D)
        self.label_28.setObjectName(u"label_28")

        self.gridLayout_4.addWidget(self.label_28, 1, 3, 1, 1)

        self.label_18 = QLabel(self.Rectangular3D)
        self.label_18.setObjectName(u"label_18")

        self.gridLayout_4.addWidget(self.label_18, 2, 7, 1, 1)

        self.INOriginXR3D = QLineEdit(self.Rectangular3D)
        self.INOriginXR3D.setObjectName(u"INOriginXR3D")

        self.gridLayout_4.addWidget(self.INOriginXR3D, 0, 2, 1, 1)

        self.label_21 = QLabel(self.Rectangular3D)
        self.label_21.setObjectName(u"label_21")

        self.gridLayout_4.addWidget(self.label_21, 0, 3, 1, 1)

        self.INElementsXR3D = QLineEdit(self.Rectangular3D)
        self.INElementsXR3D.setObjectName(u"INElementsXR3D")

        self.gridLayout_4.addWidget(self.INElementsXR3D, 2, 2, 1, 1)

        self.label_19 = QLabel(self.Rectangular3D)
        self.label_19.setObjectName(u"label_19")

        self.gridLayout_4.addWidget(self.label_19, 0, 7, 1, 1)

        self.INOriginZR3D = QLineEdit(self.Rectangular3D)
        self.INOriginZR3D.setObjectName(u"INOriginZR3D")

        self.gridLayout_4.addWidget(self.INOriginZR3D, 0, 6, 1, 1)

        self.label_31 = QLabel(self.Rectangular3D)
        self.label_31.setObjectName(u"label_31")

        self.gridLayout_4.addWidget(self.label_31, 1, 5, 1, 1)

        self.INLengthZR3D = QLineEdit(self.Rectangular3D)
        self.INLengthZR3D.setObjectName(u"INLengthZR3D")

        self.gridLayout_4.addWidget(self.INLengthZR3D, 1, 6, 1, 1)

        self.label_32 = QLabel(self.Rectangular3D)
        self.label_32.setObjectName(u"label_32")

        self.gridLayout_4.addWidget(self.label_32, 2, 5, 1, 1)

        self.INElementsZR3D = QLineEdit(self.Rectangular3D)
        self.INElementsZR3D.setObjectName(u"INElementsZR3D")

        self.gridLayout_4.addWidget(self.INElementsZR3D, 2, 6, 1, 1)

        self.MeshInputs2.addWidget(self.Rectangular3D)
        self.Cylindrical2D = QWidget()
        self.Cylindrical2D.setObjectName(u"Cylindrical2D")
        self.verticalLayout_32 = QVBoxLayout(self.Cylindrical2D)
        self.verticalLayout_32.setObjectName(u"verticalLayout_32")
        self.Cylindrical2DInputs2 = QFrame(self.Cylindrical2D)
        self.Cylindrical2DInputs2.setObjectName(u"Cylindrical2DInputs2")
        self.Cylindrical2DInputs2.setFrameShape(QFrame.NoFrame)
        self.Cylindrical2DInputs2.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_13 = QHBoxLayout(self.Cylindrical2DInputs2)
        self.horizontalLayout_13.setObjectName(u"horizontalLayout_13")
        self.horizontalLayout_13.setContentsMargins(-1, 0, -1, 0)
        self.LInnerRadiusC2D = QLabel(self.Cylindrical2DInputs2)
        self.LInnerRadiusC2D.setObjectName(u"LInnerRadiusC2D")

        self.horizontalLayout_13.addWidget(self.LInnerRadiusC2D)

        self.INInnerRadiusC2D = QLineEdit(self.Cylindrical2DInputs2)
        self.INInnerRadiusC2D.setObjectName(u"INInnerRadiusC2D")

        self.horizontalLayout_13.addWidget(self.INInnerRadiusC2D)


        self.verticalLayout_32.addWidget(self.Cylindrical2DInputs2)

        self.Cylindrical2DInputs = QFrame(self.Cylindrical2D)
        self.Cylindrical2DInputs.setObjectName(u"Cylindrical2DInputs")
        self.Cylindrical2DInputs.setFrameShape(QFrame.NoFrame)
        self.Cylindrical2DInputs.setFrameShadow(QFrame.Raised)
        self.gridLayout_5 = QGridLayout(self.Cylindrical2DInputs)
        self.gridLayout_5.setObjectName(u"gridLayout_5")
        self.gridLayout_5.setContentsMargins(-1, 0, -1, 0)
        self.label_37 = QLabel(self.Cylindrical2DInputs)
        self.label_37.setObjectName(u"label_37")

        self.gridLayout_5.addWidget(self.label_37, 0, 1, 1, 1)

        self.LOriginC2D = QLabel(self.Cylindrical2DInputs)
        self.LOriginC2D.setObjectName(u"LOriginC2D")

        self.gridLayout_5.addWidget(self.LOriginC2D, 0, 0, 1, 1)

        self.INLengthThetaC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INLengthThetaC2D.setObjectName(u"INLengthThetaC2D")

        self.gridLayout_5.addWidget(self.INLengthThetaC2D, 1, 4, 1, 1)

        self.label_38 = QLabel(self.Cylindrical2DInputs)
        self.label_38.setObjectName(u"label_38")

        self.gridLayout_5.addWidget(self.label_38, 2, 1, 1, 1)

        self.INLengthOutRadC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INLengthOutRadC2D.setObjectName(u"INLengthOutRadC2D")

        self.gridLayout_5.addWidget(self.INLengthOutRadC2D, 1, 2, 1, 1)

        self.label_40 = QLabel(self.Cylindrical2DInputs)
        self.label_40.setObjectName(u"label_40")

        self.gridLayout_5.addWidget(self.label_40, 2, 3, 1, 1)

        self.INOriginYC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INOriginYC2D.setObjectName(u"INOriginYC2D")

        self.gridLayout_5.addWidget(self.INOriginYC2D, 0, 4, 1, 1)

        self.label_36 = QLabel(self.Cylindrical2DInputs)
        self.label_36.setObjectName(u"label_36")

        self.gridLayout_5.addWidget(self.label_36, 0, 3, 1, 1)

        self.label_41 = QLabel(self.Cylindrical2DInputs)
        self.label_41.setObjectName(u"label_41")

        self.gridLayout_5.addWidget(self.label_41, 1, 1, 1, 1)

        self.label_35 = QLabel(self.Cylindrical2DInputs)
        self.label_35.setObjectName(u"label_35")

        self.gridLayout_5.addWidget(self.label_35, 1, 5, 1, 1)

        self.LLengthC2D = QLabel(self.Cylindrical2DInputs)
        self.LLengthC2D.setObjectName(u"LLengthC2D")

        self.gridLayout_5.addWidget(self.LLengthC2D, 1, 0, 1, 1)

        self.INElementsArcC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INElementsArcC2D.setObjectName(u"INElementsArcC2D")

        self.gridLayout_5.addWidget(self.INElementsArcC2D, 2, 4, 1, 1)

        self.INOriginXC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INOriginXC2D.setObjectName(u"INOriginXC2D")

        self.gridLayout_5.addWidget(self.INOriginXC2D, 0, 2, 1, 1)

        self.INElementsRadialC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INElementsRadialC2D.setObjectName(u"INElementsRadialC2D")

        self.gridLayout_5.addWidget(self.INElementsRadialC2D, 2, 2, 1, 1)

        self.LElementsC2D = QLabel(self.Cylindrical2DInputs)
        self.LElementsC2D.setObjectName(u"LElementsC2D")

        self.gridLayout_5.addWidget(self.LElementsC2D, 2, 0, 1, 1)

        self.label_34 = QLabel(self.Cylindrical2DInputs)
        self.label_34.setObjectName(u"label_34")

        self.gridLayout_5.addWidget(self.label_34, 0, 5, 1, 1)

        self.label_33 = QLabel(self.Cylindrical2DInputs)
        self.label_33.setObjectName(u"label_33")

        self.gridLayout_5.addWidget(self.label_33, 2, 5, 1, 1)

        self.label_43 = QLabel(self.Cylindrical2DInputs)
        self.label_43.setObjectName(u"label_43")

        self.gridLayout_5.addWidget(self.label_43, 1, 3, 1, 1)


        self.verticalLayout_32.addWidget(self.Cylindrical2DInputs)

        self.MeshInputs2.addWidget(self.Cylindrical2D)
        self.Cylindrical3D = QWidget()
        self.Cylindrical3D.setObjectName(u"Cylindrical3D")
        self.verticalLayout_33 = QVBoxLayout(self.Cylindrical3D)
        self.verticalLayout_33.setObjectName(u"verticalLayout_33")
        self.verticalLayout_33.setContentsMargins(0, 0, 0, 0)
        self.Cylindrical3DInputs2 = QFrame(self.Cylindrical3D)
        self.Cylindrical3DInputs2.setObjectName(u"Cylindrical3DInputs2")
        self.Cylindrical3DInputs2.setFrameShape(QFrame.NoFrame)
        self.Cylindrical3DInputs2.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_14 = QHBoxLayout(self.Cylindrical3DInputs2)
        self.horizontalLayout_14.setObjectName(u"horizontalLayout_14")
        self.horizontalLayout_14.setContentsMargins(-1, 0, -1, 0)
        self.LInnerRadiusC3D = QLabel(self.Cylindrical3DInputs2)
        self.LInnerRadiusC3D.setObjectName(u"LInnerRadiusC3D")

        self.horizontalLayout_14.addWidget(self.LInnerRadiusC3D)

        self.INInnerRadiusC3D = QLineEdit(self.Cylindrical3DInputs2)
        self.INInnerRadiusC3D.setObjectName(u"INInnerRadiusC3D")

        self.horizontalLayout_14.addWidget(self.INInnerRadiusC3D)


        self.verticalLayout_33.addWidget(self.Cylindrical3DInputs2)

        self.Cylindrical3DInputs = QFrame(self.Cylindrical3D)
        self.Cylindrical3DInputs.setObjectName(u"Cylindrical3DInputs")
        self.Cylindrical3DInputs.setFrameShape(QFrame.NoFrame)
        self.Cylindrical3DInputs.setFrameShadow(QFrame.Raised)
        self.gridLayout_6 = QGridLayout(self.Cylindrical3DInputs)
        self.gridLayout_6.setObjectName(u"gridLayout_6")
        self.gridLayout_6.setContentsMargins(-1, 0, -1, 0)
        self.INElementsRadC3D = QLineEdit(self.Cylindrical3DInputs)
        self.INElementsRadC3D.setObjectName(u"INElementsRadC3D")

        self.gridLayout_6.addWidget(self.INElementsRadC3D, 2, 2, 1, 1)

        self.label_46 = QLabel(self.Cylindrical3DInputs)
        self.label_46.setObjectName(u"label_46")

        self.gridLayout_6.addWidget(self.label_46, 0, 1, 1, 1)

        self.LOriginC3D = QLabel(self.Cylindrical3DInputs)
        self.LOriginC3D.setObjectName(u"LOriginC3D")

        self.gridLayout_6.addWidget(self.LOriginC3D, 0, 0, 1, 1)

        self.label_49 = QLabel(self.Cylindrical3DInputs)
        self.label_49.setObjectName(u"label_49")

        self.gridLayout_6.addWidget(self.label_49, 2, 3, 1, 1)

        self.LElementsC3D = QLabel(self.Cylindrical3DInputs)
        self.LElementsC3D.setObjectName(u"LElementsC3D")

        self.gridLayout_6.addWidget(self.LElementsC3D, 2, 0, 1, 1)

        self.INElementsArcC3D = QLineEdit(self.Cylindrical3DInputs)
        self.INElementsArcC3D.setObjectName(u"INElementsArcC3D")

        self.gridLayout_6.addWidget(self.INElementsArcC3D, 2, 4, 1, 1)

        self.label_51 = QLabel(self.Cylindrical3DInputs)
        self.label_51.setObjectName(u"label_51")

        self.gridLayout_6.addWidget(self.label_51, 1, 1, 1, 1)

        self.label_50 = QLabel(self.Cylindrical3DInputs)
        self.label_50.setObjectName(u"label_50")

        self.gridLayout_6.addWidget(self.label_50, 0, 3, 1, 1)

        self.label_56 = QLabel(self.Cylindrical3DInputs)
        self.label_56.setObjectName(u"label_56")

        self.gridLayout_6.addWidget(self.label_56, 2, 7, 1, 1)

        self.label_48 = QLabel(self.Cylindrical3DInputs)
        self.label_48.setObjectName(u"label_48")

        self.gridLayout_6.addWidget(self.label_48, 2, 1, 1, 1)

        self.INLengthThetaC3D = QLineEdit(self.Cylindrical3DInputs)
        self.INLengthThetaC3D.setObjectName(u"INLengthThetaC3D")

        self.gridLayout_6.addWidget(self.INLengthThetaC3D, 1, 4, 1, 1)

        self.INOriginYC3D = QLineEdit(self.Cylindrical3DInputs)
        self.INOriginYC3D.setObjectName(u"INOriginYC3D")

        self.gridLayout_6.addWidget(self.INOriginYC3D, 0, 4, 1, 1)

        self.LLengthC3D = QLabel(self.Cylindrical3DInputs)
        self.LLengthC3D.setObjectName(u"LLengthC3D")

        self.gridLayout_6.addWidget(self.LLengthC3D, 1, 0, 1, 1)

        self.INLengthOutRadC3D = QLineEdit(self.Cylindrical3DInputs)
        self.INLengthOutRadC3D.setObjectName(u"INLengthOutRadC3D")

        self.gridLayout_6.addWidget(self.INLengthOutRadC3D, 1, 2, 1, 1)

        self.INOriginXC3D = QLineEdit(self.Cylindrical3DInputs)
        self.INOriginXC3D.setObjectName(u"INOriginXC3D")

        self.gridLayout_6.addWidget(self.INOriginXC3D, 0, 2, 1, 1)

        self.label_59 = QLabel(self.Cylindrical3DInputs)
        self.label_59.setObjectName(u"label_59")

        self.gridLayout_6.addWidget(self.label_59, 0, 5, 1, 1)

        self.label_52 = QLabel(self.Cylindrical3DInputs)
        self.label_52.setObjectName(u"label_52")

        self.gridLayout_6.addWidget(self.label_52, 1, 7, 1, 1)

        self.label_57 = QLabel(self.Cylindrical3DInputs)
        self.label_57.setObjectName(u"label_57")

        self.gridLayout_6.addWidget(self.label_57, 1, 3, 1, 1)

        self.label_55 = QLabel(self.Cylindrical3DInputs)
        self.label_55.setObjectName(u"label_55")

        self.gridLayout_6.addWidget(self.label_55, 0, 7, 1, 1)

        self.INOriginZC3D = QLineEdit(self.Cylindrical3DInputs)
        self.INOriginZC3D.setObjectName(u"INOriginZC3D")

        self.gridLayout_6.addWidget(self.INOriginZC3D, 0, 6, 1, 1)

        self.label_60 = QLabel(self.Cylindrical3DInputs)
        self.label_60.setObjectName(u"label_60")

        self.gridLayout_6.addWidget(self.label_60, 1, 5, 1, 1)

        self.INLengthZC3D = QLineEdit(self.Cylindrical3DInputs)
        self.INLengthZC3D.setObjectName(u"INLengthZC3D")

        self.gridLayout_6.addWidget(self.INLengthZC3D, 1, 6, 1, 1)

        self.label_61 = QLabel(self.Cylindrical3DInputs)
        self.label_61.setObjectName(u"label_61")

        self.gridLayout_6.addWidget(self.label_61, 2, 5, 1, 1)

        self.INElementsZC3D = QLineEdit(self.Cylindrical3DInputs)
        self.INElementsZC3D.setObjectName(u"INElementsZC3D")

        self.gridLayout_6.addWidget(self.INElementsZC3D, 2, 6, 1, 1)


        self.verticalLayout_33.addWidget(self.Cylindrical3DInputs)

        self.MeshInputs2.addWidget(self.Cylindrical3D)

        self.verticalLayout_31.addWidget(self.MeshInputs2, 0, Qt.AlignTop)

        self.BGenerateGlobalMesh = QPushButton(self.GlobalMesh)
        self.BGenerateGlobalMesh.setObjectName(u"BGenerateGlobalMesh")

        self.verticalLayout_31.addWidget(self.BGenerateGlobalMesh)

        self.verticalSpacer_8 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_31.addItem(self.verticalSpacer_8)

        self.ToolSettings.addWidget(self.GlobalMesh)
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

        self.LOriginX = QLabel(self.GeometryInputs)
        self.LOriginX.setObjectName(u"LOriginX")
        self.LOriginX.setEnabled(False)

        self.formLayout.setWidget(3, QFormLayout.LabelRole, self.LOriginX)

        self.INOriginX = QLineEdit(self.GeometryInputs)
        self.INOriginX.setObjectName(u"INOriginX")
        self.INOriginX.setEnabled(False)

        self.formLayout.setWidget(3, QFormLayout.FieldRole, self.INOriginX)

        self.LOriginY = QLabel(self.GeometryInputs)
        self.LOriginY.setObjectName(u"LOriginY")
        self.LOriginY.setEnabled(False)

        self.formLayout.setWidget(4, QFormLayout.LabelRole, self.LOriginY)

        self.INOriginY = QLineEdit(self.GeometryInputs)
        self.INOriginY.setObjectName(u"INOriginY")
        self.INOriginY.setEnabled(False)

        self.formLayout.setWidget(4, QFormLayout.FieldRole, self.INOriginY)

        self.LOriginZ = QLabel(self.GeometryInputs)
        self.LOriginZ.setObjectName(u"LOriginZ")
        self.LOriginZ.setEnabled(False)

        self.formLayout.setWidget(5, QFormLayout.LabelRole, self.LOriginZ)

        self.INOriginZ = QLineEdit(self.GeometryInputs)
        self.INOriginZ.setObjectName(u"INOriginZ")
        self.INOriginZ.setEnabled(False)

        self.formLayout.setWidget(5, QFormLayout.FieldRole, self.INOriginZ)


        self.verticalLayout_15.addWidget(self.GeometryInputs)

        self.Dimensions = QFrame(self.GeometryInformationTool)
        self.Dimensions.setObjectName(u"Dimensions")
        self.Dimensions.setFrameShape(QFrame.NoFrame)
        self.Dimensions.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_15 = QHBoxLayout(self.Dimensions)
        self.horizontalLayout_15.setObjectName(u"horizontalLayout_15")
        self.horizontalLayout_15.setContentsMargins(-1, 0, -1, 0)
        self.BStlDimensions = QRadioButton(self.Dimensions)
        self.BStlDimensions.setObjectName(u"BStlDimensions")
        self.BStlDimensions.setEnabled(False)
        self.BStlDimensions.setChecked(True)

        self.horizontalLayout_15.addWidget(self.BStlDimensions)

        self.BCustomDimensions = QRadioButton(self.Dimensions)
        self.BCustomDimensions.setObjectName(u"BCustomDimensions")
        self.BCustomDimensions.setEnabled(False)

        self.horizontalLayout_15.addWidget(self.BCustomDimensions)


        self.verticalLayout_15.addWidget(self.Dimensions, 0, Qt.AlignHCenter)

        self.Lengths = QFrame(self.GeometryInformationTool)
        self.Lengths.setObjectName(u"Lengths")
        self.Lengths.setFrameShape(QFrame.NoFrame)
        self.Lengths.setFrameShadow(QFrame.Raised)
        self.formLayout_10 = QFormLayout(self.Lengths)
        self.formLayout_10.setObjectName(u"formLayout_10")
        self.formLayout_10.setContentsMargins(-1, 0, -1, 0)
        self.LLengthX = QLabel(self.Lengths)
        self.LLengthX.setObjectName(u"LLengthX")
        self.LLengthX.setEnabled(False)

        self.formLayout_10.setWidget(0, QFormLayout.LabelRole, self.LLengthX)

        self.INLengthX = QLineEdit(self.Lengths)
        self.INLengthX.setObjectName(u"INLengthX")
        self.INLengthX.setEnabled(False)

        self.formLayout_10.setWidget(0, QFormLayout.FieldRole, self.INLengthX)

        self.LLengthY = QLabel(self.Lengths)
        self.LLengthY.setObjectName(u"LLengthY")
        self.LLengthY.setEnabled(False)

        self.formLayout_10.setWidget(1, QFormLayout.LabelRole, self.LLengthY)

        self.INLengthY = QLineEdit(self.Lengths)
        self.INLengthY.setObjectName(u"INLengthY")
        self.INLengthY.setEnabled(False)

        self.formLayout_10.setWidget(1, QFormLayout.FieldRole, self.INLengthY)

        self.LLengthZ = QLabel(self.Lengths)
        self.LLengthZ.setObjectName(u"LLengthZ")
        self.LLengthZ.setEnabled(False)

        self.formLayout_10.setWidget(2, QFormLayout.LabelRole, self.LLengthZ)

        self.INLengthZ = QLineEdit(self.Lengths)
        self.INLengthZ.setObjectName(u"INLengthZ")
        self.INLengthZ.setEnabled(False)

        self.formLayout_10.setWidget(2, QFormLayout.FieldRole, self.INLengthZ)


        self.verticalLayout_15.addWidget(self.Lengths)

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

        self.LType = QLabel(self.MaterialInputs)
        self.LType.setObjectName(u"LType")

        self.gridLayout.addWidget(self.LType, 2, 0, 1, 1)

        self.frame_2 = QFrame(self.MaterialInputs)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setFrameShape(QFrame.NoFrame)
        self.frame_2.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_8 = QHBoxLayout(self.frame_2)
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.horizontalLayout_8.setContentsMargins(0, 0, 0, 0)

        self.gridLayout.addWidget(self.frame_2, 1, 4, 1, 1, Qt.AlignLeft)

        self.LRegion = QLabel(self.MaterialInputs)
        self.LRegion.setObjectName(u"LRegion")

        self.gridLayout.addWidget(self.LRegion, 1, 0, 1, 1)

        self.frame = QFrame(self.MaterialInputs)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.NoFrame)
        self.frame.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_7 = QHBoxLayout(self.frame)
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.horizontalLayout_7.setContentsMargins(0, 0, 0, 0)

        self.gridLayout.addWidget(self.frame, 2, 4, 1, 1, Qt.AlignLeft)

        self.frame_3 = QFrame(self.MaterialInputs)
        self.frame_3.setObjectName(u"frame_3")
        self.frame_3.setFrameShape(QFrame.NoFrame)
        self.frame_3.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_9 = QHBoxLayout(self.frame_3)
        self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
        self.horizontalLayout_9.setContentsMargins(0, 0, 0, 0)
        self.INSolidGas = QComboBox(self.frame_3)
        self.INSolidGas.addItem("")
        self.INSolidGas.addItem("")
        self.INSolidGas.setObjectName(u"INSolidGas")
        self.INSolidGas.setMinimumSize(QSize(82, 0))

        self.horizontalLayout_9.addWidget(self.INSolidGas)

        self.INMaterialType = QComboBox(self.frame_3)
        self.INMaterialType.addItem("")
        self.INMaterialType.addItem("")
        self.INMaterialType.addItem("")
        self.INMaterialType.addItem("")
        self.INMaterialType.setObjectName(u"INMaterialType")

        self.horizontalLayout_9.addWidget(self.INMaterialType)


        self.gridLayout.addWidget(self.frame_3, 2, 1, 1, 1)

        self.INRegion = QComboBox(self.MaterialInputs)
        self.INRegion.addItem("")
        self.INRegion.addItem("")
        self.INRegion.setObjectName(u"INRegion")

        self.gridLayout.addWidget(self.INRegion, 1, 1, 1, 1, Qt.AlignLeft)

        self.INMaterialName = QLineEdit(self.MaterialInputs)
        self.INMaterialName.setObjectName(u"INMaterialName")
        self.INMaterialName.setMinimumSize(QSize(93, 0))

        self.gridLayout.addWidget(self.INMaterialName, 0, 1, 1, 1)


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

        self.formLayout_2.setWidget(1, QFormLayout.LabelRole, self.LIsotropicPlane)

        self.INIsotropicPlane = QComboBox(self.IsotropicPlane)
        self.INIsotropicPlane.addItem("")
        self.INIsotropicPlane.addItem("")
        self.INIsotropicPlane.addItem("")
        self.INIsotropicPlane.setObjectName(u"INIsotropicPlane")

        self.formLayout_2.setWidget(1, QFormLayout.FieldRole, self.INIsotropicPlane)


        self.verticalLayout_23.addWidget(self.IsotropicPlane, 0, Qt.AlignTop)

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
        if (self.TMaterials.columnCount() < 24):
            self.TMaterials.setColumnCount(24)
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
        __qtablewidgetitem27 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(11, __qtablewidgetitem27)
        __qtablewidgetitem28 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(12, __qtablewidgetitem28)
        __qtablewidgetitem29 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(13, __qtablewidgetitem29)
        __qtablewidgetitem30 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(14, __qtablewidgetitem30)
        __qtablewidgetitem31 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(15, __qtablewidgetitem31)
        __qtablewidgetitem32 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(16, __qtablewidgetitem32)
        __qtablewidgetitem33 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(17, __qtablewidgetitem33)
        __qtablewidgetitem34 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(18, __qtablewidgetitem34)
        __qtablewidgetitem35 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(19, __qtablewidgetitem35)
        __qtablewidgetitem36 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(20, __qtablewidgetitem36)
        __qtablewidgetitem37 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(21, __qtablewidgetitem37)
        __qtablewidgetitem38 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(22, __qtablewidgetitem38)
        __qtablewidgetitem39 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(23, __qtablewidgetitem39)
        self.TMaterials.setObjectName(u"TMaterials")
        self.TMaterials.setEnabled(True)
        self.TMaterials.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.TMaterials.setRowCount(0)

        self.verticalLayout_16.addWidget(self.TMaterials)

        self.BDeleteMaterial = QPushButton(self.DefineMaterialTool)
        self.BDeleteMaterial.setObjectName(u"BDeleteMaterial")

        self.verticalLayout_16.addWidget(self.BDeleteMaterial)

        self.BRegenElasticConstants = QPushButton(self.DefineMaterialTool)
        self.BRegenElasticConstants.setObjectName(u"BRegenElasticConstants")

        self.verticalLayout_16.addWidget(self.BRegenElasticConstants)

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
        __qtablewidgetitem40 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(0, __qtablewidgetitem40)
        __qtablewidgetitem41 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(1, __qtablewidgetitem41)
        self.TBCs.setObjectName(u"TBCs")
        self.TBCs.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.TBCs.horizontalHeader().setDefaultSectionSize(175)
        self.TBCs.horizontalHeader().setStretchLastSection(True)

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
        self.INNumberOfSteps.setEnabled(False)
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

        self.THomogenization = QTableWidget(self.ResultsTool)
        if (self.THomogenization.columnCount() < 1):
            self.THomogenization.setColumnCount(1)
        __qtablewidgetitem42 = QTableWidgetItem()
        self.THomogenization.setHorizontalHeaderItem(0, __qtablewidgetitem42)
        if (self.THomogenization.rowCount() < 12):
            self.THomogenization.setRowCount(12)
        __qtablewidgetitem43 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(0, __qtablewidgetitem43)
        __qtablewidgetitem44 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(1, __qtablewidgetitem44)
        __qtablewidgetitem45 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(2, __qtablewidgetitem45)
        __qtablewidgetitem46 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(3, __qtablewidgetitem46)
        __qtablewidgetitem47 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(4, __qtablewidgetitem47)
        __qtablewidgetitem48 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(5, __qtablewidgetitem48)
        __qtablewidgetitem49 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(6, __qtablewidgetitem49)
        __qtablewidgetitem50 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(7, __qtablewidgetitem50)
        __qtablewidgetitem51 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(8, __qtablewidgetitem51)
        __qtablewidgetitem52 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(9, __qtablewidgetitem52)
        __qtablewidgetitem53 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(10, __qtablewidgetitem53)
        __qtablewidgetitem54 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(11, __qtablewidgetitem54)
        self.THomogenization.setObjectName(u"THomogenization")
        self.THomogenization.setEnabled(False)
        self.THomogenization.setWordWrap(False)
        self.THomogenization.horizontalHeader().setMinimumSectionSize(40)
        self.THomogenization.horizontalHeader().setDefaultSectionSize(250)
        self.THomogenization.horizontalHeader().setStretchLastSection(True)
        self.THomogenization.verticalHeader().setDefaultSectionSize(30)

        self.verticalLayout_18.addWidget(self.THomogenization)

        self.BHomogenization = QPushButton(self.ResultsTool)
        self.BHomogenization.setObjectName(u"BHomogenization")
        self.BHomogenization.setEnabled(False)

        self.verticalLayout_18.addWidget(self.BHomogenization)

        self.verticalSpacer_5 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_18.addItem(self.verticalSpacer_5)

        self.ToolSettings.addWidget(self.ResultsTool)
        self.splitter.addWidget(self.ToolSettings)
        self.ParaviewFrame = QFrame(self.splitter)
        self.ParaviewFrame.setObjectName(u"ParaviewFrame")
        self.ParaviewFrame.setFocusPolicy(Qt.TabFocus)
        self.ParaviewFrame.setContextMenuPolicy(Qt.DefaultContextMenu)
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

        self.verticalSpacer_10 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout.addItem(self.verticalSpacer_10)

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1074, 24))
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setObjectName(u"menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuHelp.menuAction())
        self.menuHelp.addAction(self.actionEVPFFT_Manual)

        self.retranslateUi(MainWindow)

        self.SolverTypeMenu.setCurrentIndex(0)
        self.ExplicitPipelines.setCurrentIndex(0)
        self.EVPFFTPipelines.setCurrentIndex(0)
        self.ToolSettings.setCurrentIndex(0)
        self.MeshInputs2.setCurrentIndex(0)
        self.MaterialTypeTool.setCurrentIndex(0)
        self.OutputWindows.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"EVPFFT", None))
        self.actionEVPFFT_Manual.setText(QCoreApplication.translate("MainWindow", u"EVPFFT Manual", None))
        self.radioButton_2.setText("")
        self.radioButton.setText("")
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:700; color:#ff2600;\">\u2b06</span><span style=\" font-size:18pt; font-weight:700;\"> SELECT A TYPE OF SOLVER </span><span style=\" font-size:18pt; font-weight:700; color:#ff2600;\">\u2b06</span></p></body></html>", None))
        self.SolverTypeMenu.setTabText(self.SolverTypeMenu.indexOf(self.ChooseSolver), "")
        self.BGlobalMesh.setText("")
        self.LGlobalMesh.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Global        Mesh</span></p></body></html>", None))
        self.BImportPart_4.setText("")
        self.LImportPart_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Import Part</span></p></body></html>", None))
        self.LPartTools_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Geometry</span></p></body></html>", None))
        self.BDefineMaterial_4.setText("")
        self.LDefineMaterial_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Define Material</span></p></body></html>", None))
        self.LMaterialTools_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Material</span></p></body></html>", None))
        self.BApplyBC_4.setText("")
        self.LApplyBC_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Apply BCs</span></p></body></html>", None))
        self.LBCTools_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Boundary Conditions</span></p></body></html>", None))
        self.BSolverSettings_6.setText("")
        self.LSolverSettings_7.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Explicit Settings</span></p></body></html>", None))
        self.BRunEVPFFT_6.setText("")
        self.LRunEVPFFT_6.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Run Explicit</span></p></body></html>", None))
        self.BViewResults_3.setText("")
        self.LViewResults_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">View Results</span></p></body></html>", None))
        self.LJobTools_6.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Explicit Solver</span></p></body></html>", None))
        self.ExplicitPipelines.setTabText(self.ExplicitPipelines.indexOf(self.Mesh_3), QCoreApplication.translate("MainWindow", u"General", None))
        self.SolverTypeMenu.setTabText(self.SolverTypeMenu.indexOf(self.ExplicitSolver), QCoreApplication.translate("MainWindow", u"Explicit", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:700;\">COMING SOON </span><span style=\" font-size:18pt;\">\u203c\ufe0f</span></p></body></html>", None))
        self.SolverTypeMenu.setTabText(self.SolverTypeMenu.indexOf(self.ImplicitSolver), QCoreApplication.translate("MainWindow", u"Implicit", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:700; color:#ff2600;\">\u2b06</span><span style=\" font-size:18pt; font-weight:700;\"> SELECT A PIPELINE </span><span style=\" font-size:18pt; font-weight:700; color:#ff2600;\">\u2b06</span></p></body></html>", None))
        self.EVPFFTPipelines.setTabText(self.EVPFFTPipelines.indexOf(self.ChoosePipeline), "")
        self.BImportPart_5.setText("")
        self.LImportPart_5.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Import Part</span></p></body></html>", None))
        self.LPartTools_5.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Geometry</span></p></body></html>", None))
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-size:18pt; font-weight:700;\">UNDER CONSTRUCTION \u26a0\ufe0f</span></p></body></html>", None))
        self.EVPFFTPipelines.setTabText(self.EVPFFTPipelines.indexOf(self.EVPFFTGeneral), QCoreApplication.translate("MainWindow", u"General", None))
        self.BImportPart.setText("")
        self.LImportPart.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Import Part</span></p></body></html>", None))
        self.LPartTools.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Geometry</span></p></body></html>", None))
        self.BDefineMaterial.setText("")
        self.LDefineMaterial.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Define Material</span></p></body></html>", None))
        self.LMaterialTools.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Material</span></p></body></html>", None))
        self.BApplyBC.setText("")
        self.LApplyBC.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Apply BCs</span></p></body></html>", None))
        self.LBCTools.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Boundary Conditions</span></p></body></html>", None))
        self.BSolverSettings.setText("")
        self.LSolverSettings.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Solver Settings</span></p></body></html>", None))
        self.BRunEVPFFT.setText("")
        self.LRunEVPFFT.setText(QCoreApplication.translate("MainWindow", u"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:'.AppleSystemUIFont'; font-size:13pt; font-weight:400; font-style:normal;\">\n"
"<p align=\"center\" style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:700;\">Run Job</span></p></body></html>", None))
        self.BViewResults.setText("")
        self.LViewResults.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">View Results</span></p></body></html>", None))
        self.LJobTools.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Simulation</span></p></body></html>", None))
        self.EVPFFTPipelines.setTabText(self.EVPFFTPipelines.indexOf(self.EVPFFTHomogenization), QCoreApplication.translate("MainWindow", u"Homogenization", None))
        self.SolverTypeMenu.setTabText(self.SolverTypeMenu.indexOf(self.EVPFFTSolver), QCoreApplication.translate("MainWindow", u"EVPFFT", None))
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:700;\">COMING SOON </span><span style=\" font-size:18pt;\">\u203c\ufe0f</span></p></body></html>", None))
        self.SolverTypeMenu.setTabText(self.SolverTypeMenu.indexOf(self.EVPFFTLSSolver), QCoreApplication.translate("MainWindow", u"EVPFFT-LS", None))
        self.LosAlamosLogo.setText("")
        self.EVPFFTLogo.setText("")
        self.LAdditionalSoftware.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#000000;\">Additional Software Packages:</span></p></body></html>", None))
        self.MatarLogo.setText("")
        self.ParaviewLogo.setText("")
        self.LDefineGlobalMesh.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">DEFINE GLOBAL MESH</span></p></body></html>", None))
        self.LElementType.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Element Type:</p></body></html>", None))
        self.INElementType.setItemText(0, QCoreApplication.translate("MainWindow", u"Linear", None))
        self.INElementType.setItemText(1, QCoreApplication.translate("MainWindow", u"Quadratic", None))
        self.INElementType.setItemText(2, QCoreApplication.translate("MainWindow", u"Cubic", None))

        self.LCoordinateSystem.setText(QCoreApplication.translate("MainWindow", u"Coordinate System:", None))
        self.INCoordinateSystem.setItemText(0, QCoreApplication.translate("MainWindow", u"Rectangular", None))
        self.INCoordinateSystem.setItemText(1, QCoreApplication.translate("MainWindow", u"Cylindrical", None))

        self.LDimension.setText(QCoreApplication.translate("MainWindow", u"Dimension:", None))
        self.INDimension.setItemText(0, QCoreApplication.translate("MainWindow", u"2D", None))
        self.INDimension.setItemText(1, QCoreApplication.translate("MainWindow", u"3D", None))

        self.INLengthXR2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.LLengthR2D.setText(QCoreApplication.translate("MainWindow", u"Length:", None))
        self.INElementsYR2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.label_14.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.label_11.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.label_10.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INLengthYR2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.LElementsR2D.setText(QCoreApplication.translate("MainWindow", u"Elements:", None))
        self.label_13.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INElementsXR2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.INOriginXR2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.INOriginYR2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.label_6.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.label_8.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.LOriginR2D.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Origin:</p></body></html>", None))
        self.label_15.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.label_16.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.label_17.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.LLengthR3D.setText(QCoreApplication.translate("MainWindow", u"Length:", None))
        self.label_25.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.LElementsR3D.setText(QCoreApplication.translate("MainWindow", u"Elements:", None))
        self.label_23.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.label_30.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.label_20.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.label_22.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.INLengthYR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.INOriginYR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.LOriginR3D.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Origin:</p></body></html>", None))
        self.INElementsYR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.label_26.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.INLengthXR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.label_28.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.label_18.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.INOriginXR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.label_21.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INElementsXR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.label_19.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.INOriginZR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.label_31.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INLengthZR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.label_32.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INElementsZR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.LInnerRadiusC2D.setText(QCoreApplication.translate("MainWindow", u"Inner Radius:", None))
        self.label_37.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.LOriginC2D.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Origin:</p></body></html>", None))
        self.INLengthThetaC2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"theta (deg)", None))
        self.label_38.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.INLengthOutRadC2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"outer radius", None))
        self.label_40.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INOriginYC2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.label_36.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.label_41.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.label_35.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.LLengthC2D.setText(QCoreApplication.translate("MainWindow", u"Length:", None))
        self.INElementsArcC2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"arc ", None))
        self.INOriginXC2D.setText("")
        self.INOriginXC2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.INElementsRadialC2D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"radial axis", None))
        self.LElementsC2D.setText(QCoreApplication.translate("MainWindow", u"Elements:", None))
        self.label_34.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.label_33.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.label_43.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.LInnerRadiusC3D.setText(QCoreApplication.translate("MainWindow", u"Inner Radius:", None))
        self.INElementsRadC3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"radial axis", None))
        self.label_46.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.LOriginC3D.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Origin:</p></body></html>", None))
        self.label_49.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.LElementsC3D.setText(QCoreApplication.translate("MainWindow", u"Elements:", None))
        self.INElementsArcC3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"arc ", None))
        self.label_51.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.label_50.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.label_56.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.label_48.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.INLengthThetaC3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"theta (deg)", None))
        self.INOriginYC3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.LLengthC3D.setText(QCoreApplication.translate("MainWindow", u"Length:", None))
        self.INLengthOutRadC3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"outer radius", None))
        self.INOriginXC3D.setText("")
        self.INOriginXC3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.label_59.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.label_52.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.label_57.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.label_55.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.INOriginZC3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.label_60.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INLengthZC3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.label_61.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INElementsZC3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.BGenerateGlobalMesh.setText(QCoreApplication.translate("MainWindow", u"Generate Global Mesh", None))
        self.LGeometryInformation.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">IMPORT PART</span></p></body></html>", None))
        self.BUploadGeometryFile.setText(QCoreApplication.translate("MainWindow", u"Upload Geometry File", None))
        self.STLVoxelization.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">.STL VOXELIZATION</span></p></body></html>", None))
        self.LNumberOfVoxelsX.setText(QCoreApplication.translate("MainWindow", u"Number of voxels x: ", None))
        self.LNumberOfVoxelsY.setText(QCoreApplication.translate("MainWindow", u"Number of voxels y: ", None))
        self.LNumberOfVoxelsZ.setText(QCoreApplication.translate("MainWindow", u"Number of voxels z: ", None))
        self.LOriginX.setText(QCoreApplication.translate("MainWindow", u"Origin x:", None))
        self.LOriginY.setText(QCoreApplication.translate("MainWindow", u"Origin y:", None))
        self.LOriginZ.setText(QCoreApplication.translate("MainWindow", u"Origin z:", None))
        self.BStlDimensions.setText(QCoreApplication.translate("MainWindow", u"stl dimensions", None))
        self.BCustomDimensions.setText(QCoreApplication.translate("MainWindow", u"custom dimensions", None))
        self.LLengthX.setText(QCoreApplication.translate("MainWindow", u"Length x:", None))
        self.INLengthX.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.LLengthY.setText(QCoreApplication.translate("MainWindow", u"Length y:", None))
        self.INLengthY.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.LLengthZ.setText(QCoreApplication.translate("MainWindow", u"Length z:", None))
        self.INLengthZ.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.BVoxelizeGeometry.setText(QCoreApplication.translate("MainWindow", u"Voxelize Geometry", None))
        self.LDefineMaterials.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">DEFINE MATERIALS</span></p></body></html>", None))
        self.LMaterialName.setText(QCoreApplication.translate("MainWindow", u"Name:", None))
        self.LType.setText(QCoreApplication.translate("MainWindow", u"Type:", None))
        self.LRegion.setText(QCoreApplication.translate("MainWindow", u"Region:", None))
        self.INSolidGas.setItemText(0, QCoreApplication.translate("MainWindow", u"Solid", None))
        self.INSolidGas.setItemText(1, QCoreApplication.translate("MainWindow", u"Gas", None))

        self.INMaterialType.setItemText(0, QCoreApplication.translate("MainWindow", u"Isotropic", None))
        self.INMaterialType.setItemText(1, QCoreApplication.translate("MainWindow", u"Transversely Isotropic", None))
        self.INMaterialType.setItemText(2, QCoreApplication.translate("MainWindow", u"Orthotropic", None))
        self.INMaterialType.setItemText(3, QCoreApplication.translate("MainWindow", u"Anisotropic", None))

        self.INRegion.setItemText(0, QCoreApplication.translate("MainWindow", u"Imported Part", None))
        self.INRegion.setItemText(1, QCoreApplication.translate("MainWindow", u"Void", None))

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
        ___qtablewidgetitem1.setText(QCoreApplication.translate("MainWindow", u"Region", None));
        ___qtablewidgetitem2 = self.TMaterials.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("MainWindow", u"Type", None));
        ___qtablewidgetitem3 = self.TMaterials.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("MainWindow", u"C11", None));
        ___qtablewidgetitem4 = self.TMaterials.horizontalHeaderItem(4)
        ___qtablewidgetitem4.setText(QCoreApplication.translate("MainWindow", u"C12", None));
        ___qtablewidgetitem5 = self.TMaterials.horizontalHeaderItem(5)
        ___qtablewidgetitem5.setText(QCoreApplication.translate("MainWindow", u"C13", None));
        ___qtablewidgetitem6 = self.TMaterials.horizontalHeaderItem(6)
        ___qtablewidgetitem6.setText(QCoreApplication.translate("MainWindow", u"C14", None));
        ___qtablewidgetitem7 = self.TMaterials.horizontalHeaderItem(7)
        ___qtablewidgetitem7.setText(QCoreApplication.translate("MainWindow", u"C15", None));
        ___qtablewidgetitem8 = self.TMaterials.horizontalHeaderItem(8)
        ___qtablewidgetitem8.setText(QCoreApplication.translate("MainWindow", u"C16", None));
        ___qtablewidgetitem9 = self.TMaterials.horizontalHeaderItem(9)
        ___qtablewidgetitem9.setText(QCoreApplication.translate("MainWindow", u"C22", None));
        ___qtablewidgetitem10 = self.TMaterials.horizontalHeaderItem(10)
        ___qtablewidgetitem10.setText(QCoreApplication.translate("MainWindow", u"C23", None));
        ___qtablewidgetitem11 = self.TMaterials.horizontalHeaderItem(11)
        ___qtablewidgetitem11.setText(QCoreApplication.translate("MainWindow", u"C24", None));
        ___qtablewidgetitem12 = self.TMaterials.horizontalHeaderItem(12)
        ___qtablewidgetitem12.setText(QCoreApplication.translate("MainWindow", u"C25", None));
        ___qtablewidgetitem13 = self.TMaterials.horizontalHeaderItem(13)
        ___qtablewidgetitem13.setText(QCoreApplication.translate("MainWindow", u"C26", None));
        ___qtablewidgetitem14 = self.TMaterials.horizontalHeaderItem(14)
        ___qtablewidgetitem14.setText(QCoreApplication.translate("MainWindow", u"C33", None));
        ___qtablewidgetitem15 = self.TMaterials.horizontalHeaderItem(15)
        ___qtablewidgetitem15.setText(QCoreApplication.translate("MainWindow", u"C34", None));
        ___qtablewidgetitem16 = self.TMaterials.horizontalHeaderItem(16)
        ___qtablewidgetitem16.setText(QCoreApplication.translate("MainWindow", u"C35", None));
        ___qtablewidgetitem17 = self.TMaterials.horizontalHeaderItem(17)
        ___qtablewidgetitem17.setText(QCoreApplication.translate("MainWindow", u"C36", None));
        ___qtablewidgetitem18 = self.TMaterials.horizontalHeaderItem(18)
        ___qtablewidgetitem18.setText(QCoreApplication.translate("MainWindow", u"C44", None));
        ___qtablewidgetitem19 = self.TMaterials.horizontalHeaderItem(19)
        ___qtablewidgetitem19.setText(QCoreApplication.translate("MainWindow", u"C45", None));
        ___qtablewidgetitem20 = self.TMaterials.horizontalHeaderItem(20)
        ___qtablewidgetitem20.setText(QCoreApplication.translate("MainWindow", u"C46", None));
        ___qtablewidgetitem21 = self.TMaterials.horizontalHeaderItem(21)
        ___qtablewidgetitem21.setText(QCoreApplication.translate("MainWindow", u"C55", None));
        ___qtablewidgetitem22 = self.TMaterials.horizontalHeaderItem(22)
        ___qtablewidgetitem22.setText(QCoreApplication.translate("MainWindow", u"C56", None));
        ___qtablewidgetitem23 = self.TMaterials.horizontalHeaderItem(23)
        ___qtablewidgetitem23.setText(QCoreApplication.translate("MainWindow", u"C66", None));
        self.BDeleteMaterial.setText(QCoreApplication.translate("MainWindow", u"Delete", None))
        self.BRegenElasticConstants.setText(QCoreApplication.translate("MainWindow", u"Regenerate Elastic Constants", None))
        self.LBoundaryConditions.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">BOUNDARY CONDITIONS</span></p></body></html>", None))
        self.LBoundaryCondition.setText(QCoreApplication.translate("MainWindow", u"Boundary Condition: ", None))
        self.INBoundaryCondition.setItemText(0, QCoreApplication.translate("MainWindow", u"Tension", None))
        self.INBoundaryCondition.setItemText(1, QCoreApplication.translate("MainWindow", u"Compression", None))
        self.INBoundaryCondition.setItemText(2, QCoreApplication.translate("MainWindow", u"Shear", None))
        self.INBoundaryCondition.setItemText(3, QCoreApplication.translate("MainWindow", u"Homogenization", None))

        self.LBCDirection.setText(QCoreApplication.translate("MainWindow", u"Direction: ", None))
        self.INBCDirection.setItemText(0, QCoreApplication.translate("MainWindow", u"x-direction", None))
        self.INBCDirection.setItemText(1, QCoreApplication.translate("MainWindow", u"y-direction", None))
        self.INBCDirection.setItemText(2, QCoreApplication.translate("MainWindow", u"z-direction", None))

        self.BAddBC.setText(QCoreApplication.translate("MainWindow", u"Add", None))
        ___qtablewidgetitem24 = self.TBCs.horizontalHeaderItem(0)
        ___qtablewidgetitem24.setText(QCoreApplication.translate("MainWindow", u"Boundary Condition", None));
        ___qtablewidgetitem25 = self.TBCs.horizontalHeaderItem(1)
        ___qtablewidgetitem25.setText(QCoreApplication.translate("MainWindow", u"Direction", None));
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
        ___qtablewidgetitem26 = self.THomogenization.horizontalHeaderItem(0)
        ___qtablewidgetitem26.setText(QCoreApplication.translate("MainWindow", u"Homogenized Elastic Constants", None));
        ___qtablewidgetitem27 = self.THomogenization.verticalHeaderItem(0)
        ___qtablewidgetitem27.setText(QCoreApplication.translate("MainWindow", u"E11", None));
        ___qtablewidgetitem28 = self.THomogenization.verticalHeaderItem(1)
        ___qtablewidgetitem28.setText(QCoreApplication.translate("MainWindow", u"E22", None));
        ___qtablewidgetitem29 = self.THomogenization.verticalHeaderItem(2)
        ___qtablewidgetitem29.setText(QCoreApplication.translate("MainWindow", u"E33", None));
        ___qtablewidgetitem30 = self.THomogenization.verticalHeaderItem(3)
        ___qtablewidgetitem30.setText(QCoreApplication.translate("MainWindow", u"NU12", None));
        ___qtablewidgetitem31 = self.THomogenization.verticalHeaderItem(4)
        ___qtablewidgetitem31.setText(QCoreApplication.translate("MainWindow", u"NU21", None));
        ___qtablewidgetitem32 = self.THomogenization.verticalHeaderItem(5)
        ___qtablewidgetitem32.setText(QCoreApplication.translate("MainWindow", u"NU13", None));
        ___qtablewidgetitem33 = self.THomogenization.verticalHeaderItem(6)
        ___qtablewidgetitem33.setText(QCoreApplication.translate("MainWindow", u"NU31", None));
        ___qtablewidgetitem34 = self.THomogenization.verticalHeaderItem(7)
        ___qtablewidgetitem34.setText(QCoreApplication.translate("MainWindow", u"NU23", None));
        ___qtablewidgetitem35 = self.THomogenization.verticalHeaderItem(8)
        ___qtablewidgetitem35.setText(QCoreApplication.translate("MainWindow", u"NU32", None));
        ___qtablewidgetitem36 = self.THomogenization.verticalHeaderItem(9)
        ___qtablewidgetitem36.setText(QCoreApplication.translate("MainWindow", u"G12", None));
        ___qtablewidgetitem37 = self.THomogenization.verticalHeaderItem(10)
        ___qtablewidgetitem37.setText(QCoreApplication.translate("MainWindow", u"G13", None));
        ___qtablewidgetitem38 = self.THomogenization.verticalHeaderItem(11)
        ___qtablewidgetitem38.setText(QCoreApplication.translate("MainWindow", u"G23", None));
        self.BHomogenization.setText(QCoreApplication.translate("MainWindow", u"Generate Homogenized Elastic Constants", None))
        self.menuHelp.setTitle(QCoreApplication.translate("MainWindow", u"Help", None))
    # retranslateUi

