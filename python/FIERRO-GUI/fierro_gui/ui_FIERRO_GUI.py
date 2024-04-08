# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'FIERRO_GUIHQzYlF.ui'
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
    QTableWidgetItem, QToolButton, QVBoxLayout, QWidget)
import IconResourceFile_rc
import IconResourceFile_rc

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1156, 1012)
        icon = QIcon()
        icon.addFile(u":/Logos/Logos/FIERRO.png", QSize(), QIcon.Normal, QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setAutoFillBackground(False)
        MainWindow.setToolButtonStyle(Qt.ToolButtonIconOnly)
        MainWindow.setDockNestingEnabled(False)
        self.actionManual = QAction(MainWindow)
        self.actionManual.setObjectName(u"actionManual")
        self.actionChange_Working_Directory = QAction(MainWindow)
        self.actionChange_Working_Directory.setObjectName(u"actionChange_Working_Directory")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.centralwidget.setStyleSheet(u"#TitlePage, #GeometryInformationTool, #DefineMaterialTool, #BoundaryConditionsTool, #SolverSettingsTool, #ResultsTool, #Tools, #RunOutputs, #RunOutputWindow, #Main{\n"
"    background-color: rgb(235, 235, 235);\n"
"}\n"
"#ParaviewFrame{\n"
"    background-color: rgb(91, 97, 120);\n"
"}\n"
"#BImportPart:hover, #BDefineMaterial:hover, #BApplyBC:hover, #BSolverSettings:hover, #BRunEVPFFT:hover, #BViewResults:hover, #BGlobalMesh:hover, #BImportPartSGH:hover, #BDefineMaterialSGH:hover, #BAssignMaterialSGH:hover, #BApplyBCSGH:hover, #BSolverSettingsSGH:hover, #BViewResultsSGH:hover, #BRunSGH:hover, #BCreateBasicPart:hover{\n"
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
        self.horizontalLayout_10 = QHBoxLayout(self.frame_5)
        self.horizontalLayout_10.setObjectName(u"horizontalLayout_10")
        self.horizontalLayout_10.setContentsMargins(0, 0, 0, 0)
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
        icon1 = QIcon()
        icon1.addFile(u":/Blue Icons/Blue Icons/mesh.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BGlobalMesh.setIcon(icon1)
        self.BGlobalMesh.setIconSize(QSize(40, 40))
        self.BGlobalMesh.setFlat(True)

        self.verticalLayout_34.addWidget(self.BGlobalMesh)

        self.LGlobalMesh = QLabel(self.GlobaMesh)
        self.LGlobalMesh.setObjectName(u"LGlobalMesh")
        self.LGlobalMesh.setWordWrap(True)

        self.verticalLayout_34.addWidget(self.LGlobalMesh)


        self.horizontalLayout_10.addWidget(self.GlobaMesh)

        self.ImportPart_6 = QFrame(self.frame_5)
        self.ImportPart_6.setObjectName(u"ImportPart_6")
        self.ImportPart_6.setMinimumSize(QSize(70, 80))
        self.ImportPart_6.setMaximumSize(QSize(70, 80))
        self.ImportPart_6.setFrameShape(QFrame.NoFrame)
        self.ImportPart_6.setFrameShadow(QFrame.Raised)
        self.verticalLayout_70 = QVBoxLayout(self.ImportPart_6)
        self.verticalLayout_70.setSpacing(0)
        self.verticalLayout_70.setObjectName(u"verticalLayout_70")
        self.verticalLayout_70.setContentsMargins(0, 0, 0, 5)
        self.BImportPartSGH = QPushButton(self.ImportPart_6)
        self.BImportPartSGH.setObjectName(u"BImportPartSGH")
        icon2 = QIcon()
        icon2.addFile(u":/Blue Icons/Blue Icons/DownloadCube.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BImportPartSGH.setIcon(icon2)
        self.BImportPartSGH.setIconSize(QSize(50, 50))
        self.BImportPartSGH.setFlat(True)

        self.verticalLayout_70.addWidget(self.BImportPartSGH)

        self.LImportPart_7 = QLabel(self.ImportPart_6)
        self.LImportPart_7.setObjectName(u"LImportPart_7")
        self.LImportPart_7.setWordWrap(True)

        self.verticalLayout_70.addWidget(self.LImportPart_7)


        self.horizontalLayout_10.addWidget(self.ImportPart_6)

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
        self.BCreateBasicPart = QPushButton(self.ImportPart_4)
        self.BCreateBasicPart.setObjectName(u"BCreateBasicPart")
        icon3 = QIcon()
        icon3.addFile(u":/Blue Icons/Blue Icons/Shapes.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BCreateBasicPart.setIcon(icon3)
        self.BCreateBasicPart.setIconSize(QSize(50, 50))
        self.BCreateBasicPart.setFlat(True)

        self.verticalLayout_56.addWidget(self.BCreateBasicPart)

        self.LImportPart_4 = QLabel(self.ImportPart_4)
        self.LImportPart_4.setObjectName(u"LImportPart_4")
        self.LImportPart_4.setWordWrap(True)

        self.verticalLayout_56.addWidget(self.LImportPart_4)


        self.horizontalLayout_10.addWidget(self.ImportPart_4)


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
        self.horizontalLayout_3 = QHBoxLayout(self.MaterialToolsBtns_4)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
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
        self.BDefineMaterialSGH = QPushButton(self.DefineMaterial_4)
        self.BDefineMaterialSGH.setObjectName(u"BDefineMaterialSGH")
        icon4 = QIcon()
        icon4.addFile(u":/Blue Icons/Blue Icons/mine.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BDefineMaterialSGH.setIcon(icon4)
        self.BDefineMaterialSGH.setIconSize(QSize(50, 50))
        self.BDefineMaterialSGH.setFlat(True)

        self.verticalLayout_58.addWidget(self.BDefineMaterialSGH)

        self.LDefineMaterial_4 = QLabel(self.DefineMaterial_4)
        self.LDefineMaterial_4.setObjectName(u"LDefineMaterial_4")
        self.LDefineMaterial_4.setWordWrap(True)

        self.verticalLayout_58.addWidget(self.LDefineMaterial_4, 0, Qt.AlignBottom)


        self.horizontalLayout_3.addWidget(self.DefineMaterial_4)

        self.DefineMaterial_5 = QFrame(self.MaterialToolsBtns_4)
        self.DefineMaterial_5.setObjectName(u"DefineMaterial_5")
        self.DefineMaterial_5.setMinimumSize(QSize(80, 80))
        self.DefineMaterial_5.setMaximumSize(QSize(80, 80))
        self.DefineMaterial_5.setFrameShape(QFrame.NoFrame)
        self.DefineMaterial_5.setFrameShadow(QFrame.Raised)
        self.verticalLayout_61 = QVBoxLayout(self.DefineMaterial_5)
        self.verticalLayout_61.setSpacing(0)
        self.verticalLayout_61.setObjectName(u"verticalLayout_61")
        self.verticalLayout_61.setContentsMargins(0, 0, 0, 5)
        self.BAssignMaterialSGH = QPushButton(self.DefineMaterial_5)
        self.BAssignMaterialSGH.setObjectName(u"BAssignMaterialSGH")
        icon5 = QIcon()
        icon5.addFile(u":/Blue Icons/Blue Icons/Clipboard.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BAssignMaterialSGH.setIcon(icon5)
        self.BAssignMaterialSGH.setIconSize(QSize(50, 50))
        self.BAssignMaterialSGH.setFlat(True)

        self.verticalLayout_61.addWidget(self.BAssignMaterialSGH)

        self.LDefineMaterial_5 = QLabel(self.DefineMaterial_5)
        self.LDefineMaterial_5.setObjectName(u"LDefineMaterial_5")
        self.LDefineMaterial_5.setWordWrap(True)

        self.verticalLayout_61.addWidget(self.LDefineMaterial_5, 0, Qt.AlignBottom)


        self.horizontalLayout_3.addWidget(self.DefineMaterial_5)


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
        self.BApplyBCSGH = QPushButton(self.ApplyBC_4)
        self.BApplyBCSGH.setObjectName(u"BApplyBCSGH")
        icon6 = QIcon()
        icon6.addFile(u":/Blue Icons/Blue Icons/brick.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BApplyBCSGH.setIcon(icon6)
        self.BApplyBCSGH.setIconSize(QSize(50, 50))
        self.BApplyBCSGH.setFlat(True)

        self.verticalLayout_60.addWidget(self.BApplyBCSGH)

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
        self.BSolverSettingsSGH = QPushButton(self.SolverSettings_6)
        self.BSolverSettingsSGH.setObjectName(u"BSolverSettingsSGH")
        icon7 = QIcon()
        icon7.addFile(u":/Blue Icons/Blue Icons/gear.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BSolverSettingsSGH.setIcon(icon7)
        self.BSolverSettingsSGH.setIconSize(QSize(50, 50))
        self.BSolverSettingsSGH.setFlat(True)

        self.verticalLayout_64.addWidget(self.BSolverSettingsSGH)

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
        self.BRunSGH = QPushButton(self.RunEVPFFT_6)
        self.BRunSGH.setObjectName(u"BRunSGH")
        icon8 = QIcon()
        icon8.addFile(u":/Blue Icons/Blue Icons/Play.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BRunSGH.setIcon(icon8)
        self.BRunSGH.setIconSize(QSize(40, 40))
        self.BRunSGH.setFlat(True)

        self.verticalLayout_65.addWidget(self.BRunSGH)

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
        self.BViewResultsSGH = QPushButton(self.ViewResults_3)
        self.BViewResultsSGH.setObjectName(u"BViewResultsSGH")
        icon9 = QIcon()
        icon9.addFile(u":/Blue Icons/Blue Icons/magnify.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BViewResultsSGH.setIcon(icon9)
        self.BViewResultsSGH.setIconSize(QSize(40, 40))
        self.BViewResultsSGH.setFlat(True)

        self.verticalLayout_66.addWidget(self.BViewResultsSGH)

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
        icon10 = QIcon()
        icon10.addFile(u":/Blue Icons/Blue Icons/Cube.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BImportPart.setIcon(icon10)
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
        self.BDefineMaterial.setIcon(icon4)
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

        self.frame_4 = QFrame(self.GeometryInformationTool)
        self.frame_4.setObjectName(u"frame_4")
        self.frame_4.setFrameShape(QFrame.NoFrame)
        self.frame_4.setFrameShadow(QFrame.Raised)
        self.formLayout_11 = QFormLayout(self.frame_4)
        self.formLayout_11.setObjectName(u"formLayout_11")
        self.formLayout_11.setContentsMargins(0, 0, 0, 0)
        self.LPartName = QLabel(self.frame_4)
        self.LPartName.setObjectName(u"LPartName")

        self.formLayout_11.setWidget(0, QFormLayout.LabelRole, self.LPartName)

        self.INPartName = QLineEdit(self.frame_4)
        self.INPartName.setObjectName(u"INPartName")

        self.formLayout_11.setWidget(0, QFormLayout.FieldRole, self.INPartName)


        self.verticalLayout_15.addWidget(self.frame_4)

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


        self.verticalLayout_15.addWidget(self.Dimensions)

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

        self.TParts = QTableWidget(self.GeometryInformationTool)
        if (self.TParts.columnCount() < 10):
            self.TParts.setColumnCount(10)
        __qtablewidgetitem = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        __qtablewidgetitem2 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(2, __qtablewidgetitem2)
        __qtablewidgetitem3 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(3, __qtablewidgetitem3)
        __qtablewidgetitem4 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(4, __qtablewidgetitem4)
        __qtablewidgetitem5 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(5, __qtablewidgetitem5)
        __qtablewidgetitem6 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(6, __qtablewidgetitem6)
        __qtablewidgetitem7 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(7, __qtablewidgetitem7)
        __qtablewidgetitem8 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(8, __qtablewidgetitem8)
        __qtablewidgetitem9 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(9, __qtablewidgetitem9)
        self.TParts.setObjectName(u"TParts")
        self.TParts.setEnabled(True)
        self.TParts.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.TParts.setRowCount(0)

        self.verticalLayout_15.addWidget(self.TParts)

        self.BDeleteGeometry = QPushButton(self.GeometryInformationTool)
        self.BDeleteGeometry.setObjectName(u"BDeleteGeometry")

        self.verticalLayout_15.addWidget(self.BDeleteGeometry)

        self.verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_15.addItem(self.verticalSpacer)

        self.ToolSettings.addWidget(self.GeometryInformationTool)
        self.page_5 = QWidget()
        self.page_5.setObjectName(u"page_5")
        self.verticalLayout_29 = QVBoxLayout(self.page_5)
        self.verticalLayout_29.setObjectName(u"verticalLayout_29")
        self.LBasicGeometry = QLabel(self.page_5)
        self.LBasicGeometry.setObjectName(u"LBasicGeometry")

        self.verticalLayout_29.addWidget(self.LBasicGeometry)

        self.frame_6 = QFrame(self.page_5)
        self.frame_6.setObjectName(u"frame_6")
        self.frame_6.setFrameShape(QFrame.NoFrame)
        self.frame_6.setFrameShadow(QFrame.Raised)
        self.formLayout_14 = QFormLayout(self.frame_6)
        self.formLayout_14.setObjectName(u"formLayout_14")
        self.formLayout_14.setContentsMargins(0, 0, 0, 0)
        self.INSelectBasicGeometry = QComboBox(self.frame_6)
        self.INSelectBasicGeometry.addItem("")
        self.INSelectBasicGeometry.addItem("")
        self.INSelectBasicGeometry.addItem(u"cylinder")
        self.INSelectBasicGeometry.setObjectName(u"INSelectBasicGeometry")
        self.INSelectBasicGeometry.setEnabled(True)

        self.formLayout_14.setWidget(0, QFormLayout.FieldRole, self.INSelectBasicGeometry)

        self.LSelectGeometry = QLabel(self.frame_6)
        self.LSelectGeometry.setObjectName(u"LSelectGeometry")

        self.formLayout_14.setWidget(0, QFormLayout.LabelRole, self.LSelectGeometry)

        self.INBasicGeometryName = QLineEdit(self.frame_6)
        self.INBasicGeometryName.setObjectName(u"INBasicGeometryName")

        self.formLayout_14.setWidget(1, QFormLayout.FieldRole, self.INBasicGeometryName)

        self.LBasicGName = QLabel(self.frame_6)
        self.LBasicGName.setObjectName(u"LBasicGName")

        self.formLayout_14.setWidget(1, QFormLayout.LabelRole, self.LBasicGName)


        self.verticalLayout_29.addWidget(self.frame_6)

        self.BasicGeometries = QStackedWidget(self.page_5)
        self.BasicGeometries.setObjectName(u"BasicGeometries")
        self.BoxProperties = QWidget()
        self.BoxProperties.setObjectName(u"BoxProperties")
        self.verticalLayout_35 = QVBoxLayout(self.BoxProperties)
        self.verticalLayout_35.setObjectName(u"verticalLayout_35")
        self.verticalLayout_35.setContentsMargins(-1, 0, -1, -1)
        self.LBoxProperties = QLabel(self.BoxProperties)
        self.LBoxProperties.setObjectName(u"LBoxProperties")

        self.verticalLayout_35.addWidget(self.LBoxProperties)

        self.frame_10 = QFrame(self.BoxProperties)
        self.frame_10.setObjectName(u"frame_10")
        self.frame_10.setFrameShape(QFrame.NoFrame)
        self.frame_10.setFrameShadow(QFrame.Raised)
        self.formLayout_15 = QFormLayout(self.frame_10)
        self.formLayout_15.setObjectName(u"formLayout_15")
        self.formLayout_15.setContentsMargins(0, 0, 0, 0)
        self.LBoxx1 = QLabel(self.frame_10)
        self.LBoxx1.setObjectName(u"LBoxx1")

        self.formLayout_15.setWidget(0, QFormLayout.LabelRole, self.LBoxx1)

        self.INBoxx1 = QLineEdit(self.frame_10)
        self.INBoxx1.setObjectName(u"INBoxx1")

        self.formLayout_15.setWidget(0, QFormLayout.FieldRole, self.INBoxx1)

        self.LBoxx2 = QLabel(self.frame_10)
        self.LBoxx2.setObjectName(u"LBoxx2")

        self.formLayout_15.setWidget(1, QFormLayout.LabelRole, self.LBoxx2)

        self.INBoxx2 = QLineEdit(self.frame_10)
        self.INBoxx2.setObjectName(u"INBoxx2")

        self.formLayout_15.setWidget(1, QFormLayout.FieldRole, self.INBoxx2)

        self.LBoxy1 = QLabel(self.frame_10)
        self.LBoxy1.setObjectName(u"LBoxy1")

        self.formLayout_15.setWidget(2, QFormLayout.LabelRole, self.LBoxy1)

        self.INBoxy1 = QLineEdit(self.frame_10)
        self.INBoxy1.setObjectName(u"INBoxy1")

        self.formLayout_15.setWidget(2, QFormLayout.FieldRole, self.INBoxy1)

        self.LBoxy2 = QLabel(self.frame_10)
        self.LBoxy2.setObjectName(u"LBoxy2")

        self.formLayout_15.setWidget(3, QFormLayout.LabelRole, self.LBoxy2)

        self.INBoxy2 = QLineEdit(self.frame_10)
        self.INBoxy2.setObjectName(u"INBoxy2")

        self.formLayout_15.setWidget(3, QFormLayout.FieldRole, self.INBoxy2)

        self.LBoxz1 = QLabel(self.frame_10)
        self.LBoxz1.setObjectName(u"LBoxz1")

        self.formLayout_15.setWidget(4, QFormLayout.LabelRole, self.LBoxz1)

        self.INBoxz1 = QLineEdit(self.frame_10)
        self.INBoxz1.setObjectName(u"INBoxz1")

        self.formLayout_15.setWidget(4, QFormLayout.FieldRole, self.INBoxz1)

        self.LBoxz2 = QLabel(self.frame_10)
        self.LBoxz2.setObjectName(u"LBoxz2")

        self.formLayout_15.setWidget(5, QFormLayout.LabelRole, self.LBoxz2)

        self.INBoxz2 = QLineEdit(self.frame_10)
        self.INBoxz2.setObjectName(u"INBoxz2")

        self.formLayout_15.setWidget(5, QFormLayout.FieldRole, self.INBoxz2)


        self.verticalLayout_35.addWidget(self.frame_10)

        self.BasicGeometries.addWidget(self.BoxProperties)
        self.SphereProperties = QWidget()
        self.SphereProperties.setObjectName(u"SphereProperties")
        self.verticalLayout_36 = QVBoxLayout(self.SphereProperties)
        self.verticalLayout_36.setObjectName(u"verticalLayout_36")
        self.LSphereProperties = QLabel(self.SphereProperties)
        self.LSphereProperties.setObjectName(u"LSphereProperties")

        self.verticalLayout_36.addWidget(self.LSphereProperties)

        self.frame_11 = QFrame(self.SphereProperties)
        self.frame_11.setObjectName(u"frame_11")
        self.frame_11.setFrameShape(QFrame.NoFrame)
        self.frame_11.setFrameShadow(QFrame.Raised)
        self.formLayout_16 = QFormLayout(self.frame_11)
        self.formLayout_16.setObjectName(u"formLayout_16")
        self.formLayout_16.setContentsMargins(0, 0, 0, 0)
        self.LSphereri = QLabel(self.frame_11)
        self.LSphereri.setObjectName(u"LSphereri")

        self.formLayout_16.setWidget(0, QFormLayout.LabelRole, self.LSphereri)

        self.INSphereri = QLineEdit(self.frame_11)
        self.INSphereri.setObjectName(u"INSphereri")

        self.formLayout_16.setWidget(0, QFormLayout.FieldRole, self.INSphereri)

        self.LSpherero = QLabel(self.frame_11)
        self.LSpherero.setObjectName(u"LSpherero")

        self.formLayout_16.setWidget(1, QFormLayout.LabelRole, self.LSpherero)

        self.INSpherero = QLineEdit(self.frame_11)
        self.INSpherero.setObjectName(u"INSpherero")

        self.formLayout_16.setWidget(1, QFormLayout.FieldRole, self.INSpherero)

        self.LSphereox = QLabel(self.frame_11)
        self.LSphereox.setObjectName(u"LSphereox")

        self.formLayout_16.setWidget(2, QFormLayout.LabelRole, self.LSphereox)

        self.INSphereox = QLineEdit(self.frame_11)
        self.INSphereox.setObjectName(u"INSphereox")

        self.formLayout_16.setWidget(2, QFormLayout.FieldRole, self.INSphereox)

        self.INSphereoy = QLineEdit(self.frame_11)
        self.INSphereoy.setObjectName(u"INSphereoy")

        self.formLayout_16.setWidget(3, QFormLayout.FieldRole, self.INSphereoy)

        self.LSphereoy = QLabel(self.frame_11)
        self.LSphereoy.setObjectName(u"LSphereoy")

        self.formLayout_16.setWidget(3, QFormLayout.LabelRole, self.LSphereoy)

        self.LSphereoz = QLabel(self.frame_11)
        self.LSphereoz.setObjectName(u"LSphereoz")

        self.formLayout_16.setWidget(4, QFormLayout.LabelRole, self.LSphereoz)

        self.INSphereoz = QLineEdit(self.frame_11)
        self.INSphereoz.setObjectName(u"INSphereoz")

        self.formLayout_16.setWidget(4, QFormLayout.FieldRole, self.INSphereoz)


        self.verticalLayout_36.addWidget(self.frame_11)

        self.verticalSpacer_15 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_36.addItem(self.verticalSpacer_15)

        self.BasicGeometries.addWidget(self.SphereProperties)
        self.CylinderProperties = QWidget()
        self.CylinderProperties.setObjectName(u"CylinderProperties")
        self.verticalLayout_37 = QVBoxLayout(self.CylinderProperties)
        self.verticalLayout_37.setObjectName(u"verticalLayout_37")
        self.LCylinderProperties = QLabel(self.CylinderProperties)
        self.LCylinderProperties.setObjectName(u"LCylinderProperties")

        self.verticalLayout_37.addWidget(self.LCylinderProperties)

        self.frame_12 = QFrame(self.CylinderProperties)
        self.frame_12.setObjectName(u"frame_12")
        self.frame_12.setFrameShape(QFrame.NoFrame)
        self.frame_12.setFrameShadow(QFrame.Raised)
        self.formLayout_17 = QFormLayout(self.frame_12)
        self.formLayout_17.setObjectName(u"formLayout_17")
        self.formLayout_17.setContentsMargins(0, 0, 0, 0)
        self.LCylinderri = QLabel(self.frame_12)
        self.LCylinderri.setObjectName(u"LCylinderri")

        self.formLayout_17.setWidget(0, QFormLayout.LabelRole, self.LCylinderri)

        self.INCylinderri = QLineEdit(self.frame_12)
        self.INCylinderri.setObjectName(u"INCylinderri")

        self.formLayout_17.setWidget(0, QFormLayout.FieldRole, self.INCylinderri)

        self.LCylinderro = QLabel(self.frame_12)
        self.LCylinderro.setObjectName(u"LCylinderro")

        self.formLayout_17.setWidget(1, QFormLayout.LabelRole, self.LCylinderro)

        self.INCylinderro = QLineEdit(self.frame_12)
        self.INCylinderro.setObjectName(u"INCylinderro")

        self.formLayout_17.setWidget(1, QFormLayout.FieldRole, self.INCylinderro)


        self.verticalLayout_37.addWidget(self.frame_12)

        self.verticalSpacer_16 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_37.addItem(self.verticalSpacer_16)

        self.BasicGeometries.addWidget(self.CylinderProperties)

        self.verticalLayout_29.addWidget(self.BasicGeometries, 0, Qt.AlignTop)

        self.BGenerateBasicGeometry = QPushButton(self.page_5)
        self.BGenerateBasicGeometry.setObjectName(u"BGenerateBasicGeometry")

        self.verticalLayout_29.addWidget(self.BGenerateBasicGeometry)

        self.TBasicGeometries = QTableWidget(self.page_5)
        if (self.TBasicGeometries.columnCount() < 13):
            self.TBasicGeometries.setColumnCount(13)
        __qtablewidgetitem10 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(0, __qtablewidgetitem10)
        __qtablewidgetitem11 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(1, __qtablewidgetitem11)
        __qtablewidgetitem12 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(2, __qtablewidgetitem12)
        __qtablewidgetitem13 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(3, __qtablewidgetitem13)
        __qtablewidgetitem14 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(4, __qtablewidgetitem14)
        __qtablewidgetitem15 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(5, __qtablewidgetitem15)
        __qtablewidgetitem16 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(6, __qtablewidgetitem16)
        __qtablewidgetitem17 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(7, __qtablewidgetitem17)
        __qtablewidgetitem18 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(8, __qtablewidgetitem18)
        __qtablewidgetitem19 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(9, __qtablewidgetitem19)
        __qtablewidgetitem20 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(10, __qtablewidgetitem20)
        __qtablewidgetitem21 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(11, __qtablewidgetitem21)
        __qtablewidgetitem22 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(12, __qtablewidgetitem22)
        self.TBasicGeometries.setObjectName(u"TBasicGeometries")

        self.verticalLayout_29.addWidget(self.TBasicGeometries)

        self.BDeleteBasicGeometry = QPushButton(self.page_5)
        self.BDeleteBasicGeometry.setObjectName(u"BDeleteBasicGeometry")

        self.verticalLayout_29.addWidget(self.BDeleteBasicGeometry)

        self.verticalSpacer_14 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_29.addItem(self.verticalSpacer_14)

        self.ToolSettings.addWidget(self.page_5)
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
        __qtablewidgetitem23 = QTableWidgetItem()
        self.TAnisotropic.setItem(0, 0, __qtablewidgetitem23)
        brush = QBrush(QColor(235, 235, 235, 255))
        brush.setStyle(Qt.SolidPattern)
        __qtablewidgetitem24 = QTableWidgetItem()
        __qtablewidgetitem24.setBackground(brush);
        __qtablewidgetitem24.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(1, 0, __qtablewidgetitem24)
        __qtablewidgetitem25 = QTableWidgetItem()
        __qtablewidgetitem25.setBackground(brush);
        __qtablewidgetitem25.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(2, 0, __qtablewidgetitem25)
        __qtablewidgetitem26 = QTableWidgetItem()
        __qtablewidgetitem26.setBackground(brush);
        __qtablewidgetitem26.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(2, 1, __qtablewidgetitem26)
        __qtablewidgetitem27 = QTableWidgetItem()
        __qtablewidgetitem27.setBackground(brush);
        __qtablewidgetitem27.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 0, __qtablewidgetitem27)
        __qtablewidgetitem28 = QTableWidgetItem()
        __qtablewidgetitem28.setBackground(brush);
        __qtablewidgetitem28.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 1, __qtablewidgetitem28)
        __qtablewidgetitem29 = QTableWidgetItem()
        __qtablewidgetitem29.setBackground(brush);
        __qtablewidgetitem29.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 2, __qtablewidgetitem29)
        __qtablewidgetitem30 = QTableWidgetItem()
        __qtablewidgetitem30.setBackground(brush);
        __qtablewidgetitem30.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 0, __qtablewidgetitem30)
        __qtablewidgetitem31 = QTableWidgetItem()
        __qtablewidgetitem31.setBackground(brush);
        __qtablewidgetitem31.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 1, __qtablewidgetitem31)
        __qtablewidgetitem32 = QTableWidgetItem()
        __qtablewidgetitem32.setBackground(brush);
        __qtablewidgetitem32.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 2, __qtablewidgetitem32)
        __qtablewidgetitem33 = QTableWidgetItem()
        __qtablewidgetitem33.setBackground(brush);
        __qtablewidgetitem33.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 3, __qtablewidgetitem33)
        __qtablewidgetitem34 = QTableWidgetItem()
        __qtablewidgetitem34.setBackground(brush);
        __qtablewidgetitem34.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 0, __qtablewidgetitem34)
        __qtablewidgetitem35 = QTableWidgetItem()
        __qtablewidgetitem35.setBackground(brush);
        __qtablewidgetitem35.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 1, __qtablewidgetitem35)
        __qtablewidgetitem36 = QTableWidgetItem()
        __qtablewidgetitem36.setBackground(brush);
        __qtablewidgetitem36.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 2, __qtablewidgetitem36)
        __qtablewidgetitem37 = QTableWidgetItem()
        __qtablewidgetitem37.setBackground(brush);
        __qtablewidgetitem37.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 3, __qtablewidgetitem37)
        __qtablewidgetitem38 = QTableWidgetItem()
        __qtablewidgetitem38.setBackground(brush);
        __qtablewidgetitem38.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 4, __qtablewidgetitem38)
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
        __qtablewidgetitem39 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(0, __qtablewidgetitem39)
        __qtablewidgetitem40 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(1, __qtablewidgetitem40)
        __qtablewidgetitem41 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(2, __qtablewidgetitem41)
        __qtablewidgetitem42 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(3, __qtablewidgetitem42)
        __qtablewidgetitem43 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(4, __qtablewidgetitem43)
        __qtablewidgetitem44 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(5, __qtablewidgetitem44)
        __qtablewidgetitem45 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(6, __qtablewidgetitem45)
        __qtablewidgetitem46 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(7, __qtablewidgetitem46)
        __qtablewidgetitem47 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(8, __qtablewidgetitem47)
        __qtablewidgetitem48 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(9, __qtablewidgetitem48)
        __qtablewidgetitem49 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(10, __qtablewidgetitem49)
        __qtablewidgetitem50 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(11, __qtablewidgetitem50)
        __qtablewidgetitem51 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(12, __qtablewidgetitem51)
        __qtablewidgetitem52 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(13, __qtablewidgetitem52)
        __qtablewidgetitem53 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(14, __qtablewidgetitem53)
        __qtablewidgetitem54 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(15, __qtablewidgetitem54)
        __qtablewidgetitem55 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(16, __qtablewidgetitem55)
        __qtablewidgetitem56 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(17, __qtablewidgetitem56)
        __qtablewidgetitem57 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(18, __qtablewidgetitem57)
        __qtablewidgetitem58 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(19, __qtablewidgetitem58)
        __qtablewidgetitem59 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(20, __qtablewidgetitem59)
        __qtablewidgetitem60 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(21, __qtablewidgetitem60)
        __qtablewidgetitem61 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(22, __qtablewidgetitem61)
        __qtablewidgetitem62 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(23, __qtablewidgetitem62)
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
        __qtablewidgetitem63 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(0, __qtablewidgetitem63)
        __qtablewidgetitem64 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(1, __qtablewidgetitem64)
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


        self.verticalLayout_18.addWidget(self.PreviewResults)

        self.BOpenParaview = QPushButton(self.ResultsTool)
        self.BOpenParaview.setObjectName(u"BOpenParaview")

        self.verticalLayout_18.addWidget(self.BOpenParaview)

        self.THomogenization = QTableWidget(self.ResultsTool)
        if (self.THomogenization.columnCount() < 1):
            self.THomogenization.setColumnCount(1)
        __qtablewidgetitem65 = QTableWidgetItem()
        self.THomogenization.setHorizontalHeaderItem(0, __qtablewidgetitem65)
        if (self.THomogenization.rowCount() < 12):
            self.THomogenization.setRowCount(12)
        __qtablewidgetitem66 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(0, __qtablewidgetitem66)
        __qtablewidgetitem67 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(1, __qtablewidgetitem67)
        __qtablewidgetitem68 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(2, __qtablewidgetitem68)
        __qtablewidgetitem69 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(3, __qtablewidgetitem69)
        __qtablewidgetitem70 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(4, __qtablewidgetitem70)
        __qtablewidgetitem71 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(5, __qtablewidgetitem71)
        __qtablewidgetitem72 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(6, __qtablewidgetitem72)
        __qtablewidgetitem73 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(7, __qtablewidgetitem73)
        __qtablewidgetitem74 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(8, __qtablewidgetitem74)
        __qtablewidgetitem75 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(9, __qtablewidgetitem75)
        __qtablewidgetitem76 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(10, __qtablewidgetitem76)
        __qtablewidgetitem77 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(11, __qtablewidgetitem77)
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
        self.page = QWidget()
        self.page.setObjectName(u"page")
        self.verticalLayout_12 = QVBoxLayout(self.page)
        self.verticalLayout_12.setObjectName(u"verticalLayout_12")
        self.LDefineMaterials_2 = QLabel(self.page)
        self.LDefineMaterials_2.setObjectName(u"LDefineMaterials_2")

        self.verticalLayout_12.addWidget(self.LDefineMaterials_2)

        self.frame_7 = QFrame(self.page)
        self.frame_7.setObjectName(u"frame_7")
        self.frame_7.setFrameShape(QFrame.NoFrame)
        self.frame_7.setFrameShadow(QFrame.Raised)
        self.gridLayout_8 = QGridLayout(self.frame_7)
        self.gridLayout_8.setObjectName(u"gridLayout_8")
        self.gridLayout_8.setContentsMargins(0, 0, 0, 0)
        self.LEOS = QLabel(self.frame_7)
        self.LEOS.setObjectName(u"LEOS")

        self.gridLayout_8.addWidget(self.LEOS, 2, 0, 1, 1)

        self.LMaterialNameSGH = QLabel(self.frame_7)
        self.LMaterialNameSGH.setObjectName(u"LMaterialNameSGH")

        self.gridLayout_8.addWidget(self.LMaterialNameSGH, 1, 0, 1, 1)

        self.INEOS = QComboBox(self.frame_7)
        self.INEOS.addItem("")
        self.INEOS.setObjectName(u"INEOS")

        self.gridLayout_8.addWidget(self.INEOS, 2, 1, 1, 1)

        self.INMaterialNameSGH = QLineEdit(self.frame_7)
        self.INMaterialNameSGH.setObjectName(u"INMaterialNameSGH")
        self.INMaterialNameSGH.setMinimumSize(QSize(93, 0))

        self.gridLayout_8.addWidget(self.INMaterialNameSGH, 1, 1, 1, 1)

        self.LArtificialViscosity = QLabel(self.frame_7)
        self.LArtificialViscosity.setObjectName(u"LArtificialViscosity")

        self.gridLayout_8.addWidget(self.LArtificialViscosity, 3, 0, 1, 1, Qt.AlignRight)

        self.INArtificialViscosity = QComboBox(self.frame_7)
        self.INArtificialViscosity.addItem("")
        self.INArtificialViscosity.addItem("")
        self.INArtificialViscosity.setObjectName(u"INArtificialViscosity")

        self.gridLayout_8.addWidget(self.INArtificialViscosity, 3, 1, 1, 1)


        self.verticalLayout_12.addWidget(self.frame_7, 0, Qt.AlignHCenter)

        self.BAddMaterialSGH = QPushButton(self.page)
        self.BAddMaterialSGH.setObjectName(u"BAddMaterialSGH")

        self.verticalLayout_12.addWidget(self.BAddMaterialSGH)

        self.TMaterialsSGH = QTableWidget(self.page)
        if (self.TMaterialsSGH.columnCount() < 9):
            self.TMaterialsSGH.setColumnCount(9)
        __qtablewidgetitem78 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(0, __qtablewidgetitem78)
        __qtablewidgetitem79 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(1, __qtablewidgetitem79)
        __qtablewidgetitem80 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(2, __qtablewidgetitem80)
        __qtablewidgetitem81 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(3, __qtablewidgetitem81)
        __qtablewidgetitem82 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(4, __qtablewidgetitem82)
        __qtablewidgetitem83 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(5, __qtablewidgetitem83)
        __qtablewidgetitem84 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(6, __qtablewidgetitem84)
        __qtablewidgetitem85 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(7, __qtablewidgetitem85)
        __qtablewidgetitem86 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(8, __qtablewidgetitem86)
        self.TMaterialsSGH.setObjectName(u"TMaterialsSGH")
        self.TMaterialsSGH.setEnabled(True)
        self.TMaterialsSGH.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.TMaterialsSGH.setRowCount(0)

        self.verticalLayout_12.addWidget(self.TMaterialsSGH)

        self.BDeleteMaterialSGH = QPushButton(self.page)
        self.BDeleteMaterialSGH.setObjectName(u"BDeleteMaterialSGH")

        self.verticalLayout_12.addWidget(self.BDeleteMaterialSGH)

        self.ArtificialViscosity = QStackedWidget(self.page)
        self.ArtificialViscosity.setObjectName(u"ArtificialViscosity")
        sizePolicy.setHeightForWidth(self.ArtificialViscosity.sizePolicy().hasHeightForWidth())
        self.ArtificialViscosity.setSizePolicy(sizePolicy)
        self.page_6 = QWidget()
        self.page_6.setObjectName(u"page_6")
        sizePolicy.setHeightForWidth(self.page_6.sizePolicy().hasHeightForWidth())
        self.page_6.setSizePolicy(sizePolicy)
        self.page_6.setMaximumSize(QSize(0, 0))
        self.ArtificialViscosity.addWidget(self.page_6)
        self.page_7 = QWidget()
        self.page_7.setObjectName(u"page_7")
        sizePolicy.setHeightForWidth(self.page_7.sizePolicy().hasHeightForWidth())
        self.page_7.setSizePolicy(sizePolicy)
        self.verticalLayout_39 = QVBoxLayout(self.page_7)
        self.verticalLayout_39.setObjectName(u"verticalLayout_39")
        self.label_12 = QLabel(self.page_7)
        self.label_12.setObjectName(u"label_12")

        self.verticalLayout_39.addWidget(self.label_12)

        self.frame_16 = QFrame(self.page_7)
        self.frame_16.setObjectName(u"frame_16")
        self.frame_16.setFrameShape(QFrame.NoFrame)
        self.frame_16.setFrameShadow(QFrame.Raised)
        self.gridLayout_7 = QGridLayout(self.frame_16)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.gridLayout_7.setContentsMargins(0, 0, 0, 0)
        self.Lq1 = QLabel(self.frame_16)
        self.Lq1.setObjectName(u"Lq1")
        self.Lq1.setEnabled(True)
        self.Lq1.setWordWrap(True)

        self.gridLayout_7.addWidget(self.Lq1, 0, 0, 1, 1)

        self.INq1 = QLineEdit(self.frame_16)
        self.INq1.setObjectName(u"INq1")
        self.INq1.setEnabled(True)

        self.gridLayout_7.addWidget(self.INq1, 0, 1, 1, 1)

        self.Lq2 = QLabel(self.frame_16)
        self.Lq2.setObjectName(u"Lq2")
        self.Lq2.setEnabled(True)
        self.Lq2.setWordWrap(True)

        self.gridLayout_7.addWidget(self.Lq2, 1, 0, 1, 1)

        self.INq2 = QLineEdit(self.frame_16)
        self.INq2.setObjectName(u"INq2")
        self.INq2.setEnabled(True)

        self.gridLayout_7.addWidget(self.INq2, 1, 1, 1, 1)

        self.Lq1ex = QLabel(self.frame_16)
        self.Lq1ex.setObjectName(u"Lq1ex")
        self.Lq1ex.setEnabled(True)
        self.Lq1ex.setWordWrap(True)

        self.gridLayout_7.addWidget(self.Lq1ex, 2, 0, 1, 1)

        self.INq1ex = QLineEdit(self.frame_16)
        self.INq1ex.setObjectName(u"INq1ex")
        self.INq1ex.setEnabled(True)

        self.gridLayout_7.addWidget(self.INq1ex, 2, 1, 1, 1)

        self.Lq2ex = QLabel(self.frame_16)
        self.Lq2ex.setObjectName(u"Lq2ex")
        self.Lq2ex.setEnabled(True)
        self.Lq2ex.setWordWrap(True)

        self.gridLayout_7.addWidget(self.Lq2ex, 3, 0, 1, 1)

        self.INq2ex = QLineEdit(self.frame_16)
        self.INq2ex.setObjectName(u"INq2ex")
        self.INq2ex.setEnabled(True)

        self.gridLayout_7.addWidget(self.INq2ex, 3, 1, 1, 1)

        self.LGamma = QLabel(self.frame_16)
        self.LGamma.setObjectName(u"LGamma")
        self.LGamma.setEnabled(True)

        self.gridLayout_7.addWidget(self.LGamma, 4, 0, 1, 1)

        self.INGamma = QLineEdit(self.frame_16)
        self.INGamma.setObjectName(u"INGamma")
        self.INGamma.setEnabled(True)

        self.gridLayout_7.addWidget(self.INGamma, 4, 1, 1, 1)

        self.LMinSound = QLabel(self.frame_16)
        self.LMinSound.setObjectName(u"LMinSound")
        self.LMinSound.setEnabled(True)

        self.gridLayout_7.addWidget(self.LMinSound, 5, 0, 1, 1)

        self.INMinSound = QLineEdit(self.frame_16)
        self.INMinSound.setObjectName(u"INMinSound")
        self.INMinSound.setEnabled(True)

        self.gridLayout_7.addWidget(self.INMinSound, 5, 1, 1, 1)

        self.LSpecificHeat = QLabel(self.frame_16)
        self.LSpecificHeat.setObjectName(u"LSpecificHeat")
        self.LSpecificHeat.setEnabled(True)

        self.gridLayout_7.addWidget(self.LSpecificHeat, 6, 0, 1, 1)

        self.INSpecificHeat = QLineEdit(self.frame_16)
        self.INSpecificHeat.setObjectName(u"INSpecificHeat")
        self.INSpecificHeat.setEnabled(True)

        self.gridLayout_7.addWidget(self.INSpecificHeat, 6, 1, 1, 1)


        self.verticalLayout_39.addWidget(self.frame_16)

        self.ArtificialViscosity.addWidget(self.page_7)

        self.verticalLayout_12.addWidget(self.ArtificialViscosity)

        self.verticalSpacer_9 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_12.addItem(self.verticalSpacer_9)

        self.ToolSettings.addWidget(self.page)
        self.page_2 = QWidget()
        self.page_2.setObjectName(u"page_2")
        self.verticalLayout_13 = QVBoxLayout(self.page_2)
        self.verticalLayout_13.setObjectName(u"verticalLayout_13")
        self.LDefineMaterials_3 = QLabel(self.page_2)
        self.LDefineMaterials_3.setObjectName(u"LDefineMaterials_3")

        self.verticalLayout_13.addWidget(self.LDefineMaterials_3)

        self.frame_9 = QFrame(self.page_2)
        self.frame_9.setObjectName(u"frame_9")
        self.frame_9.setFrameShape(QFrame.NoFrame)
        self.frame_9.setFrameShadow(QFrame.Raised)
        self.gridLayout_9 = QGridLayout(self.frame_9)
        self.gridLayout_9.setObjectName(u"gridLayout_9")
        self.gridLayout_9.setContentsMargins(0, 0, 0, 0)
        self.LVelx = QLabel(self.frame_9)
        self.LVelx.setObjectName(u"LVelx")

        self.gridLayout_9.addWidget(self.LVelx, 5, 0, 1, 1)

        self.LRegion_3 = QLabel(self.frame_9)
        self.LRegion_3.setObjectName(u"LRegion_3")

        self.gridLayout_9.addWidget(self.LRegion_3, 2, 0, 1, 1)

        self.INVelocityZ = QLineEdit(self.frame_9)
        self.INVelocityZ.setObjectName(u"INVelocityZ")
        self.INVelocityZ.setEnabled(True)

        self.gridLayout_9.addWidget(self.INVelocityZ, 7, 1, 1, 1)

        self.LMaterialName_3 = QLabel(self.frame_9)
        self.LMaterialName_3.setObjectName(u"LMaterialName_3")

        self.gridLayout_9.addWidget(self.LMaterialName_3, 1, 0, 1, 1)

        self.LVelz = QLabel(self.frame_9)
        self.LVelz.setObjectName(u"LVelz")

        self.gridLayout_9.addWidget(self.LVelz, 7, 0, 1, 1)

        self.INMaterial = QComboBox(self.frame_9)
        self.INMaterial.setObjectName(u"INMaterial")

        self.gridLayout_9.addWidget(self.INMaterial, 2, 1, 1, 1)

        self.LSIE = QLabel(self.frame_9)
        self.LSIE.setObjectName(u"LSIE")
        self.LSIE.setWordWrap(False)

        self.gridLayout_9.addWidget(self.LSIE, 4, 0, 1, 1)

        self.LVely = QLabel(self.frame_9)
        self.LVely.setObjectName(u"LVely")

        self.gridLayout_9.addWidget(self.LVely, 6, 0, 1, 1)

        self.INPartMaterial = QComboBox(self.frame_9)
        self.INPartMaterial.setObjectName(u"INPartMaterial")

        self.gridLayout_9.addWidget(self.INPartMaterial, 1, 1, 1, 1)

        self.INDensity = QLineEdit(self.frame_9)
        self.INDensity.setObjectName(u"INDensity")
        self.INDensity.setEnabled(True)

        self.gridLayout_9.addWidget(self.INDensity, 3, 1, 1, 1)

        self.INSIE = QLineEdit(self.frame_9)
        self.INSIE.setObjectName(u"INSIE")
        self.INSIE.setEnabled(True)

        self.gridLayout_9.addWidget(self.INSIE, 4, 1, 1, 1)

        self.LDensity = QLabel(self.frame_9)
        self.LDensity.setObjectName(u"LDensity")

        self.gridLayout_9.addWidget(self.LDensity, 3, 0, 1, 1)

        self.INVelocityX = QLineEdit(self.frame_9)
        self.INVelocityX.setObjectName(u"INVelocityX")
        self.INVelocityX.setEnabled(True)

        self.gridLayout_9.addWidget(self.INVelocityX, 5, 1, 1, 1)

        self.INVelocityY = QLineEdit(self.frame_9)
        self.INVelocityY.setObjectName(u"INVelocityY")
        self.INVelocityY.setEnabled(True)

        self.gridLayout_9.addWidget(self.INVelocityY, 6, 1, 1, 1)


        self.verticalLayout_13.addWidget(self.frame_9)

        self.Baddmaterialassignment = QPushButton(self.page_2)
        self.Baddmaterialassignment.setObjectName(u"Baddmaterialassignment")

        self.verticalLayout_13.addWidget(self.Baddmaterialassignment)

        self.Tassignmat = QTableWidget(self.page_2)
        if (self.Tassignmat.columnCount() < 7):
            self.Tassignmat.setColumnCount(7)
        __qtablewidgetitem87 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(0, __qtablewidgetitem87)
        __qtablewidgetitem88 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(1, __qtablewidgetitem88)
        __qtablewidgetitem89 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(2, __qtablewidgetitem89)
        __qtablewidgetitem90 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(3, __qtablewidgetitem90)
        __qtablewidgetitem91 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(4, __qtablewidgetitem91)
        __qtablewidgetitem92 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(5, __qtablewidgetitem92)
        __qtablewidgetitem93 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(6, __qtablewidgetitem93)
        self.Tassignmat.setObjectName(u"Tassignmat")

        self.verticalLayout_13.addWidget(self.Tassignmat)

        self.frame_17 = QFrame(self.page_2)
        self.frame_17.setObjectName(u"frame_17")
        self.frame_17.setFrameShape(QFrame.NoFrame)
        self.frame_17.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_18 = QHBoxLayout(self.frame_17)
        self.horizontalLayout_18.setObjectName(u"horizontalLayout_18")
        self.horizontalLayout_18.setContentsMargins(0, 0, 0, 0)
        self.BUpMaterial = QToolButton(self.frame_17)
        self.BUpMaterial.setObjectName(u"BUpMaterial")
        self.BUpMaterial.setIconSize(QSize(32, 32))
        self.BUpMaterial.setAutoRaise(False)
        self.BUpMaterial.setArrowType(Qt.UpArrow)

        self.horizontalLayout_18.addWidget(self.BUpMaterial)

        self.BDownMaterial = QToolButton(self.frame_17)
        self.BDownMaterial.setObjectName(u"BDownMaterial")
        self.BDownMaterial.setIconSize(QSize(32, 32))
        self.BDownMaterial.setArrowType(Qt.DownArrow)

        self.horizontalLayout_18.addWidget(self.BDownMaterial)


        self.verticalLayout_13.addWidget(self.frame_17, 0, Qt.AlignHCenter)

        self.Bdeletematerialassignment = QPushButton(self.page_2)
        self.Bdeletematerialassignment.setObjectName(u"Bdeletematerialassignment")

        self.verticalLayout_13.addWidget(self.Bdeletematerialassignment)

        self.verticalSpacer_11 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_13.addItem(self.verticalSpacer_11)

        self.ToolSettings.addWidget(self.page_2)
        self.page_3 = QWidget()
        self.page_3.setObjectName(u"page_3")
        self.verticalLayout_27 = QVBoxLayout(self.page_3)
        self.verticalLayout_27.setObjectName(u"verticalLayout_27")
        self.Lbcs = QLabel(self.page_3)
        self.Lbcs.setObjectName(u"Lbcs")

        self.verticalLayout_27.addWidget(self.Lbcs)

        self.bcsettings = QFrame(self.page_3)
        self.bcsettings.setObjectName(u"bcsettings")
        self.bcsettings.setFrameShape(QFrame.NoFrame)
        self.bcsettings.setFrameShadow(QFrame.Raised)
        self.formLayout_12 = QFormLayout(self.bcsettings)
        self.formLayout_12.setObjectName(u"formLayout_12")
        self.formLayout_12.setContentsMargins(0, 0, 0, 0)
        self.Lbndry = QLabel(self.bcsettings)
        self.Lbndry.setObjectName(u"Lbndry")

        self.formLayout_12.setWidget(0, QFormLayout.LabelRole, self.Lbndry)

        self.INBoundary = QComboBox(self.bcsettings)
        self.INBoundary.addItem("")
        self.INBoundary.addItem("")
        self.INBoundary.addItem("")
        self.INBoundary.setObjectName(u"INBoundary")

        self.formLayout_12.setWidget(0, QFormLayout.FieldRole, self.INBoundary)

        self.Lvalue = QLabel(self.bcsettings)
        self.Lvalue.setObjectName(u"Lvalue")

        self.formLayout_12.setWidget(1, QFormLayout.LabelRole, self.Lvalue)

        self.INPlanePosition = QLineEdit(self.bcsettings)
        self.INPlanePosition.setObjectName(u"INPlanePosition")

        self.formLayout_12.setWidget(1, QFormLayout.FieldRole, self.INPlanePosition)

        self.Ltype = QLabel(self.bcsettings)
        self.Ltype.setObjectName(u"Ltype")

        self.formLayout_12.setWidget(2, QFormLayout.LabelRole, self.Ltype)

        self.INType = QComboBox(self.bcsettings)
        self.INType.addItem("")
        self.INType.addItem("")
        self.INType.addItem("")
        self.INType.setObjectName(u"INType")

        self.formLayout_12.setWidget(2, QFormLayout.FieldRole, self.INType)

        self.INVel0 = QLineEdit(self.bcsettings)
        self.INVel0.setObjectName(u"INVel0")
        self.INVel0.setEnabled(False)

        self.formLayout_12.setWidget(3, QFormLayout.FieldRole, self.INVel0)

        self.INVel1 = QLineEdit(self.bcsettings)
        self.INVel1.setObjectName(u"INVel1")
        self.INVel1.setEnabled(False)

        self.formLayout_12.setWidget(4, QFormLayout.FieldRole, self.INVel1)

        self.INVelstart = QLineEdit(self.bcsettings)
        self.INVelstart.setObjectName(u"INVelstart")
        self.INVelstart.setEnabled(False)

        self.formLayout_12.setWidget(5, QFormLayout.FieldRole, self.INVelstart)

        self.INVelend = QLineEdit(self.bcsettings)
        self.INVelend.setObjectName(u"INVelend")
        self.INVelend.setEnabled(False)

        self.formLayout_12.setWidget(6, QFormLayout.FieldRole, self.INVelend)

        self.LVel0 = QLabel(self.bcsettings)
        self.LVel0.setObjectName(u"LVel0")
        self.LVel0.setEnabled(False)

        self.formLayout_12.setWidget(3, QFormLayout.LabelRole, self.LVel0)

        self.LVel1 = QLabel(self.bcsettings)
        self.LVel1.setObjectName(u"LVel1")
        self.LVel1.setEnabled(False)

        self.formLayout_12.setWidget(4, QFormLayout.LabelRole, self.LVel1)

        self.LVelstart = QLabel(self.bcsettings)
        self.LVelstart.setObjectName(u"LVelstart")
        self.LVelstart.setEnabled(False)

        self.formLayout_12.setWidget(5, QFormLayout.LabelRole, self.LVelstart)

        self.LVelend = QLabel(self.bcsettings)
        self.LVelend.setObjectName(u"LVelend")
        self.LVelend.setEnabled(False)

        self.formLayout_12.setWidget(6, QFormLayout.LabelRole, self.LVelend)


        self.verticalLayout_27.addWidget(self.bcsettings)

        self.BaddBC = QPushButton(self.page_3)
        self.BaddBC.setObjectName(u"BaddBC")

        self.verticalLayout_27.addWidget(self.BaddBC)

        self.TBoundaryConditions = QTableWidget(self.page_3)
        if (self.TBoundaryConditions.columnCount() < 7):
            self.TBoundaryConditions.setColumnCount(7)
        __qtablewidgetitem94 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(0, __qtablewidgetitem94)
        __qtablewidgetitem95 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(1, __qtablewidgetitem95)
        __qtablewidgetitem96 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(2, __qtablewidgetitem96)
        __qtablewidgetitem97 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(3, __qtablewidgetitem97)
        __qtablewidgetitem98 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(4, __qtablewidgetitem98)
        __qtablewidgetitem99 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(5, __qtablewidgetitem99)
        __qtablewidgetitem100 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(6, __qtablewidgetitem100)
        self.TBoundaryConditions.setObjectName(u"TBoundaryConditions")

        self.verticalLayout_27.addWidget(self.TBoundaryConditions)

        self.BdeleteBC = QPushButton(self.page_3)
        self.BdeleteBC.setObjectName(u"BdeleteBC")

        self.verticalLayout_27.addWidget(self.BdeleteBC)

        self.verticalSpacer_12 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_27.addItem(self.verticalSpacer_12)

        self.ToolSettings.addWidget(self.page_3)
        self.page_4 = QWidget()
        self.page_4.setObjectName(u"page_4")
        self.verticalLayout_28 = QVBoxLayout(self.page_4)
        self.verticalLayout_28.setObjectName(u"verticalLayout_28")
        self.label_45 = QLabel(self.page_4)
        self.label_45.setObjectName(u"label_45")

        self.verticalLayout_28.addWidget(self.label_45)

        self.solversettings = QFrame(self.page_4)
        self.solversettings.setObjectName(u"solversettings")
        self.solversettings.setFrameShape(QFrame.NoFrame)
        self.solversettings.setFrameShadow(QFrame.Raised)
        self.formLayout_13 = QFormLayout(self.solversettings)
        self.formLayout_13.setObjectName(u"formLayout_13")
        self.Ltime = QLabel(self.solversettings)
        self.Ltime.setObjectName(u"Ltime")

        self.formLayout_13.setWidget(0, QFormLayout.LabelRole, self.Ltime)

        self.INTime = QLineEdit(self.solversettings)
        self.INTime.setObjectName(u"INTime")

        self.formLayout_13.setWidget(0, QFormLayout.FieldRole, self.INTime)

        self.Lmindt = QLabel(self.solversettings)
        self.Lmindt.setObjectName(u"Lmindt")

        self.formLayout_13.setWidget(1, QFormLayout.LabelRole, self.Lmindt)

        self.INMindt = QLineEdit(self.solversettings)
        self.INMindt.setObjectName(u"INMindt")

        self.formLayout_13.setWidget(1, QFormLayout.FieldRole, self.INMindt)

        self.Lmaxdt = QLabel(self.solversettings)
        self.Lmaxdt.setObjectName(u"Lmaxdt")

        self.formLayout_13.setWidget(2, QFormLayout.LabelRole, self.Lmaxdt)

        self.INMaxdt = QLineEdit(self.solversettings)
        self.INMaxdt.setObjectName(u"INMaxdt")

        self.formLayout_13.setWidget(2, QFormLayout.FieldRole, self.INMaxdt)

        self.Linitdt = QLabel(self.solversettings)
        self.Linitdt.setObjectName(u"Linitdt")

        self.formLayout_13.setWidget(3, QFormLayout.LabelRole, self.Linitdt)

        self.INInitialdt = QLineEdit(self.solversettings)
        self.INInitialdt.setObjectName(u"INInitialdt")

        self.formLayout_13.setWidget(3, QFormLayout.FieldRole, self.INInitialdt)

        self.Lmaxcycle = QLabel(self.solversettings)
        self.Lmaxcycle.setObjectName(u"Lmaxcycle")

        self.formLayout_13.setWidget(4, QFormLayout.LabelRole, self.Lmaxcycle)

        self.INmaxcycles = QLineEdit(self.solversettings)
        self.INmaxcycles.setObjectName(u"INmaxcycles")

        self.formLayout_13.setWidget(4, QFormLayout.FieldRole, self.INmaxcycles)

        self.verticalSpacer_13 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.formLayout_13.setItem(6, QFormLayout.FieldRole, self.verticalSpacer_13)

        self.INGraphicsOutput = QLineEdit(self.solversettings)
        self.INGraphicsOutput.setObjectName(u"INGraphicsOutput")

        self.formLayout_13.setWidget(5, QFormLayout.FieldRole, self.INGraphicsOutput)

        self.LGraphicsOutput = QLabel(self.solversettings)
        self.LGraphicsOutput.setObjectName(u"LGraphicsOutput")
        self.LGraphicsOutput.setWordWrap(True)

        self.formLayout_13.setWidget(5, QFormLayout.LabelRole, self.LGraphicsOutput)


        self.verticalLayout_28.addWidget(self.solversettings)

        self.ToolSettings.addWidget(self.page_4)
        self.ResultsSGH = QWidget()
        self.ResultsSGH.setObjectName(u"ResultsSGH")
        self.verticalLayout_38 = QVBoxLayout(self.ResultsSGH)
        self.verticalLayout_38.setObjectName(u"verticalLayout_38")
        self.LResultsSGH = QLabel(self.ResultsSGH)
        self.LResultsSGH.setObjectName(u"LResultsSGH")

        self.verticalLayout_38.addWidget(self.LResultsSGH)

        self.frame_15 = QFrame(self.ResultsSGH)
        self.frame_15.setObjectName(u"frame_15")
        self.frame_15.setFrameShape(QFrame.NoFrame)
        self.frame_15.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_17 = QHBoxLayout(self.frame_15)
        self.horizontalLayout_17.setObjectName(u"horizontalLayout_17")
        self.horizontalLayout_17.setContentsMargins(0, 0, 0, 0)
        self.INOuputVarSGH = QComboBox(self.frame_15)
        self.INOuputVarSGH.addItem("")
        self.INOuputVarSGH.setObjectName(u"INOuputVarSGH")

        self.horizontalLayout_17.addWidget(self.INOuputVarSGH)

        self.BPreviewResultsSGH = QPushButton(self.frame_15)
        self.BPreviewResultsSGH.setObjectName(u"BPreviewResultsSGH")

        self.horizontalLayout_17.addWidget(self.BPreviewResultsSGH)


        self.verticalLayout_38.addWidget(self.frame_15, 0, Qt.AlignHCenter)

        self.frame_13 = QFrame(self.ResultsSGH)
        self.frame_13.setObjectName(u"frame_13")
        self.frame_13.setFrameShape(QFrame.NoFrame)
        self.frame_13.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_12 = QHBoxLayout(self.frame_13)
        self.horizontalLayout_12.setObjectName(u"horizontalLayout_12")
        self.horizontalLayout_12.setContentsMargins(0, 0, 0, 0)
        self.BFirstFrame = QToolButton(self.frame_13)
        self.BFirstFrame.setObjectName(u"BFirstFrame")
        icon11 = QIcon()
        icon11.addFile(u":/Blue Icons/Blue Icons/FirstFrame.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BFirstFrame.setIcon(icon11)
        self.BFirstFrame.setIconSize(QSize(32, 32))
        self.BFirstFrame.setPopupMode(QToolButton.DelayedPopup)
        self.BFirstFrame.setToolButtonStyle(Qt.ToolButtonIconOnly)
        self.BFirstFrame.setAutoRaise(False)
        self.BFirstFrame.setArrowType(Qt.NoArrow)

        self.horizontalLayout_12.addWidget(self.BFirstFrame)

        self.BPreviousFrame = QToolButton(self.frame_13)
        self.BPreviousFrame.setObjectName(u"BPreviousFrame")
        icon12 = QIcon()
        icon12.addFile(u":/Blue Icons/Blue Icons/PreviousFrame.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BPreviousFrame.setIcon(icon12)
        self.BPreviousFrame.setIconSize(QSize(32, 32))

        self.horizontalLayout_12.addWidget(self.BPreviousFrame)

        self.BNextFrame = QToolButton(self.frame_13)
        self.BNextFrame.setObjectName(u"BNextFrame")
        icon13 = QIcon()
        icon13.addFile(u":/Blue Icons/Blue Icons/NextFrame.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BNextFrame.setIcon(icon13)
        self.BNextFrame.setIconSize(QSize(32, 32))

        self.horizontalLayout_12.addWidget(self.BNextFrame)

        self.BLastFrame = QToolButton(self.frame_13)
        self.BLastFrame.setObjectName(u"BLastFrame")
        icon14 = QIcon()
        icon14.addFile(u":/Blue Icons/Blue Icons/LastFrame.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BLastFrame.setIcon(icon14)
        self.BLastFrame.setIconSize(QSize(32, 32))

        self.horizontalLayout_12.addWidget(self.BLastFrame)


        self.verticalLayout_38.addWidget(self.frame_13, 0, Qt.AlignHCenter)

        self.frame_14 = QFrame(self.ResultsSGH)
        self.frame_14.setObjectName(u"frame_14")
        self.frame_14.setFrameShape(QFrame.NoFrame)
        self.frame_14.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_16 = QHBoxLayout(self.frame_14)
        self.horizontalLayout_16.setObjectName(u"horizontalLayout_16")
        self.horizontalLayout_16.setContentsMargins(0, 0, 0, 0)
        self.LThreshold = QLabel(self.frame_14)
        self.LThreshold.setObjectName(u"LThreshold")

        self.horizontalLayout_16.addWidget(self.LThreshold)

        self.INThreshold = QLineEdit(self.frame_14)
        self.INThreshold.setObjectName(u"INThreshold")

        self.horizontalLayout_16.addWidget(self.INThreshold)

        self.BThreshold = QToolButton(self.frame_14)
        self.BThreshold.setObjectName(u"BThreshold")
        icon15 = QIcon()
        icon15.addFile(u":/Blue Icons/Blue Icons/crop.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BThreshold.setIcon(icon15)
        self.BThreshold.setIconSize(QSize(32, 32))
        self.BThreshold.setToolButtonStyle(Qt.ToolButtonIconOnly)

        self.horizontalLayout_16.addWidget(self.BThreshold)


        self.verticalLayout_38.addWidget(self.frame_14, 0, Qt.AlignHCenter)

        self.BOpenParaviewSGH = QPushButton(self.ResultsSGH)
        self.BOpenParaviewSGH.setObjectName(u"BOpenParaviewSGH")

        self.verticalLayout_38.addWidget(self.BOpenParaviewSGH)

        self.verticalSpacer_17 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.verticalLayout_38.addItem(self.verticalSpacer_17)

        self.ToolSettings.addWidget(self.ResultsSGH)
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
        self.menubar.setGeometry(QRect(0, 0, 1156, 24))
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setObjectName(u"menuHelp")
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.menuHelp.addAction(self.actionManual)
        self.menuFile.addAction(self.actionChange_Working_Directory)

        self.retranslateUi(MainWindow)

        self.SolverTypeMenu.setCurrentIndex(0)
        self.ExplicitPipelines.setCurrentIndex(0)
        self.EVPFFTPipelines.setCurrentIndex(0)
        self.ToolSettings.setCurrentIndex(0)
        self.MeshInputs2.setCurrentIndex(0)
        self.BasicGeometries.setCurrentIndex(0)
        self.MaterialTypeTool.setCurrentIndex(0)
        self.ArtificialViscosity.setCurrentIndex(0)
        self.OutputWindows.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"Fierro", None))
        self.actionManual.setText(QCoreApplication.translate("MainWindow", u"Manual", None))
        self.actionChange_Working_Directory.setText(QCoreApplication.translate("MainWindow", u"Change Fierro Setup", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><br/></p></body></html>", None))
        self.SolverTypeMenu.setTabText(self.SolverTypeMenu.indexOf(self.ChooseSolver), "")
        self.BGlobalMesh.setText("")
        self.LGlobalMesh.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Global        Mesh</span></p></body></html>", None))
        self.BImportPartSGH.setText("")
        self.LImportPart_7.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Import Part</span></p></body></html>", None))
        self.BCreateBasicPart.setText("")
        self.LImportPart_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Basic Fill Region</span></p></body></html>", None))
        self.LPartTools_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Geometry</span></p></body></html>", None))
        self.BDefineMaterialSGH.setText("")
        self.LDefineMaterial_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Define Materials</span></p></body></html>", None))
        self.BAssignMaterialSGH.setText("")
        self.LDefineMaterial_5.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Assign Materials</span></p></body></html>", None))
        self.LMaterialTools_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Material</span></p></body></html>", None))
        self.BApplyBCSGH.setText("")
        self.LApplyBC_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Apply BCs</span></p></body></html>", None))
        self.LBCTools_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Boundary Conditions</span></p></body></html>", None))
        self.BSolverSettingsSGH.setText("")
        self.LSolverSettings_7.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Solver Settings</span></p></body></html>", None))
        self.BRunSGH.setText("")
        self.LRunEVPFFT_6.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Run Job</span></p></body></html>", None))
        self.BViewResultsSGH.setText("")
        self.LViewResults_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">View Results</span></p></body></html>", None))
        self.LJobTools_6.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#a9a9a9;\">Simulation</span></p></body></html>", None))
        self.ExplicitPipelines.setTabText(self.ExplicitPipelines.indexOf(self.Mesh_3), QCoreApplication.translate("MainWindow", u"SGH", None))
        self.SolverTypeMenu.setTabText(self.SolverTypeMenu.indexOf(self.ExplicitSolver), QCoreApplication.translate("MainWindow", u"Dynamic", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:700;\">COMING SOON </span><span style=\" font-size:18pt;\">\u203c\ufe0f</span></p></body></html>", None))
        self.SolverTypeMenu.setTabText(self.SolverTypeMenu.indexOf(self.ImplicitSolver), QCoreApplication.translate("MainWindow", u"Quasi-Static", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><br/></p></body></html>", None))
        self.EVPFFTPipelines.setTabText(self.EVPFFTPipelines.indexOf(self.ChoosePipeline), "")
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:700;\">COMING SOON </span><span style=\" font-size:18pt;\">\u203c\ufe0f</span></p></body></html>", None))
        self.EVPFFTPipelines.setTabText(self.EVPFFTPipelines.indexOf(self.EVPFFTGeneral), QCoreApplication.translate("MainWindow", u"Polycrystal", None))
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
        self.EVPFFTPipelines.setTabText(self.EVPFFTPipelines.indexOf(self.EVPFFTHomogenization), QCoreApplication.translate("MainWindow", u"Bulk Response", None))
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
        self.INDimension.setItemText(0, QCoreApplication.translate("MainWindow", u"3D", None))
        self.INDimension.setItemText(1, QCoreApplication.translate("MainWindow", u"2D", None))

        self.LLengthR3D.setText(QCoreApplication.translate("MainWindow", u"Length:", None))
        self.label_25.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.LElementsR3D.setText(QCoreApplication.translate("MainWindow", u"Elements:", None))
        self.label_23.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.label_30.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.label_20.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.label_22.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.INLengthYR3D.setText("")
        self.INLengthYR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.INOriginYR3D.setText("")
        self.INOriginYR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.LOriginR3D.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Origin:</p></body></html>", None))
        self.INElementsYR3D.setText("")
        self.INElementsYR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.label_26.setText(QCoreApplication.translate("MainWindow", u"(", None))
        self.INLengthXR3D.setText("")
        self.INLengthXR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.label_28.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.label_18.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.INOriginXR3D.setText("")
        self.INOriginXR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.label_21.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INElementsXR3D.setText("")
        self.INElementsXR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"X", None))
        self.label_19.setText(QCoreApplication.translate("MainWindow", u")", None))
        self.INOriginZR3D.setText("")
        self.INOriginZR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.label_31.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INLengthZR3D.setText("")
        self.INLengthZR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.label_32.setText(QCoreApplication.translate("MainWindow", u",", None))
        self.INElementsZR3D.setText("")
        self.INElementsZR3D.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Z", None))
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
        self.LPartName.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>        Part Name: </p></body></html>", None))
        self.INPartName.setText("")
        self.BUploadGeometryFile.setText(QCoreApplication.translate("MainWindow", u"Upload Geometry File", None))
        self.STLVoxelization.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">.STL VOXELIZATION</span></p></body></html>", None))
        self.LNumberOfVoxelsX.setText(QCoreApplication.translate("MainWindow", u"Number of voxels x: ", None))
        self.INNumberOfVoxelsX.setText("")
        self.LNumberOfVoxelsY.setText(QCoreApplication.translate("MainWindow", u"Number of voxels y: ", None))
        self.INNumberOfVoxelsY.setText("")
        self.LNumberOfVoxelsZ.setText(QCoreApplication.translate("MainWindow", u"Number of voxels z: ", None))
        self.INNumberOfVoxelsZ.setText("")
        self.LOriginX.setText(QCoreApplication.translate("MainWindow", u"Origin x:", None))
        self.INOriginX.setText("")
        self.LOriginY.setText(QCoreApplication.translate("MainWindow", u"Origin y:", None))
        self.INOriginY.setText("")
        self.LOriginZ.setText(QCoreApplication.translate("MainWindow", u"Origin z:", None))
        self.INOriginZ.setText("")
        self.BStlDimensions.setText(QCoreApplication.translate("MainWindow", u"stl dimensions", None))
        self.BCustomDimensions.setText(QCoreApplication.translate("MainWindow", u"custom dimensions", None))
        self.LLengthX.setText(QCoreApplication.translate("MainWindow", u"Length x:", None))
        self.INLengthX.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.LLengthY.setText(QCoreApplication.translate("MainWindow", u"Length y:", None))
        self.INLengthY.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.LLengthZ.setText(QCoreApplication.translate("MainWindow", u"Length z:", None))
        self.INLengthZ.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.BVoxelizeGeometry.setText(QCoreApplication.translate("MainWindow", u"Voxelize Geometry", None))
        ___qtablewidgetitem = self.TParts.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem1 = self.TParts.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("MainWindow", u"Origin x", None));
        ___qtablewidgetitem2 = self.TParts.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("MainWindow", u"Origin y", None));
        ___qtablewidgetitem3 = self.TParts.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("MainWindow", u"Origin z", None));
        ___qtablewidgetitem4 = self.TParts.horizontalHeaderItem(4)
        ___qtablewidgetitem4.setText(QCoreApplication.translate("MainWindow", u"Length x", None));
        ___qtablewidgetitem5 = self.TParts.horizontalHeaderItem(5)
        ___qtablewidgetitem5.setText(QCoreApplication.translate("MainWindow", u"Length y", None));
        ___qtablewidgetitem6 = self.TParts.horizontalHeaderItem(6)
        ___qtablewidgetitem6.setText(QCoreApplication.translate("MainWindow", u"Length z", None));
        ___qtablewidgetitem7 = self.TParts.horizontalHeaderItem(7)
        ___qtablewidgetitem7.setText(QCoreApplication.translate("MainWindow", u"Vox x", None));
        ___qtablewidgetitem8 = self.TParts.horizontalHeaderItem(8)
        ___qtablewidgetitem8.setText(QCoreApplication.translate("MainWindow", u"Vox y", None));
        ___qtablewidgetitem9 = self.TParts.horizontalHeaderItem(9)
        ___qtablewidgetitem9.setText(QCoreApplication.translate("MainWindow", u"Vox z", None));
        self.BDeleteGeometry.setText(QCoreApplication.translate("MainWindow", u"Delete Geometry", None))
        self.LBasicGeometry.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">CREATE BASIC PART</span></p></body></html>", None))
        self.INSelectBasicGeometry.setItemText(0, QCoreApplication.translate("MainWindow", u"box", None))
        self.INSelectBasicGeometry.setItemText(1, QCoreApplication.translate("MainWindow", u"sphere", None))

        self.LSelectGeometry.setText(QCoreApplication.translate("MainWindow", u"Select Geometry:", None))
        self.LBasicGName.setText(QCoreApplication.translate("MainWindow", u"Part Name:", None))
        self.LBoxProperties.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" text-decoration: underline;\">Box Properties</span></p></body></html>", None))
        self.LBoxx1.setText(QCoreApplication.translate("MainWindow", u"x1:", None))
        self.LBoxx2.setText(QCoreApplication.translate("MainWindow", u"x2:", None))
        self.LBoxy1.setText(QCoreApplication.translate("MainWindow", u"y1:", None))
        self.LBoxy2.setText(QCoreApplication.translate("MainWindow", u"y2:", None))
        self.LBoxz1.setText(QCoreApplication.translate("MainWindow", u"z1:", None))
        self.LBoxz2.setText(QCoreApplication.translate("MainWindow", u"z2:", None))
        self.LSphereProperties.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" text-decoration: underline;\">3D Sphere Properties</span></p></body></html>", None))
        self.LSphereri.setText(QCoreApplication.translate("MainWindow", u"Inner Radius:", None))
        self.LSpherero.setText(QCoreApplication.translate("MainWindow", u"Outer Radius:", None))
        self.LSphereox.setText(QCoreApplication.translate("MainWindow", u"Origin x:", None))
        self.LSphereoy.setText(QCoreApplication.translate("MainWindow", u"Origin y:", None))
        self.LSphereoz.setText(QCoreApplication.translate("MainWindow", u"Origin z:", None))
        self.LCylinderProperties.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" text-decoration: underline;\">2D Cylinder Properties</span></p></body></html>", None))
        self.LCylinderri.setText(QCoreApplication.translate("MainWindow", u"Inner Radius:", None))
        self.LCylinderro.setText(QCoreApplication.translate("MainWindow", u"Outer Radius:", None))
        self.BGenerateBasicGeometry.setText(QCoreApplication.translate("MainWindow", u"Generate Geometry", None))
        ___qtablewidgetitem10 = self.TBasicGeometries.horizontalHeaderItem(0)
        ___qtablewidgetitem10.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem11 = self.TBasicGeometries.horizontalHeaderItem(1)
        ___qtablewidgetitem11.setText(QCoreApplication.translate("MainWindow", u"Type", None));
        ___qtablewidgetitem12 = self.TBasicGeometries.horizontalHeaderItem(2)
        ___qtablewidgetitem12.setText(QCoreApplication.translate("MainWindow", u"x1", None));
        ___qtablewidgetitem13 = self.TBasicGeometries.horizontalHeaderItem(3)
        ___qtablewidgetitem13.setText(QCoreApplication.translate("MainWindow", u"x2", None));
        ___qtablewidgetitem14 = self.TBasicGeometries.horizontalHeaderItem(4)
        ___qtablewidgetitem14.setText(QCoreApplication.translate("MainWindow", u"y1", None));
        ___qtablewidgetitem15 = self.TBasicGeometries.horizontalHeaderItem(5)
        ___qtablewidgetitem15.setText(QCoreApplication.translate("MainWindow", u"y2", None));
        ___qtablewidgetitem16 = self.TBasicGeometries.horizontalHeaderItem(6)
        ___qtablewidgetitem16.setText(QCoreApplication.translate("MainWindow", u"z1", None));
        ___qtablewidgetitem17 = self.TBasicGeometries.horizontalHeaderItem(7)
        ___qtablewidgetitem17.setText(QCoreApplication.translate("MainWindow", u"z2", None));
        ___qtablewidgetitem18 = self.TBasicGeometries.horizontalHeaderItem(8)
        ___qtablewidgetitem18.setText(QCoreApplication.translate("MainWindow", u"Inner Radius", None));
        ___qtablewidgetitem19 = self.TBasicGeometries.horizontalHeaderItem(9)
        ___qtablewidgetitem19.setText(QCoreApplication.translate("MainWindow", u"Outer Radius", None));
        ___qtablewidgetitem20 = self.TBasicGeometries.horizontalHeaderItem(10)
        ___qtablewidgetitem20.setText(QCoreApplication.translate("MainWindow", u"Origin x", None));
        ___qtablewidgetitem21 = self.TBasicGeometries.horizontalHeaderItem(11)
        ___qtablewidgetitem21.setText(QCoreApplication.translate("MainWindow", u"Origin y", None));
        ___qtablewidgetitem22 = self.TBasicGeometries.horizontalHeaderItem(12)
        ___qtablewidgetitem22.setText(QCoreApplication.translate("MainWindow", u"Origin z", None));
        self.BDeleteBasicGeometry.setText(QCoreApplication.translate("MainWindow", u"Delete Geometry", None))
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
        ___qtablewidgetitem23 = self.TMaterials.horizontalHeaderItem(0)
        ___qtablewidgetitem23.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem24 = self.TMaterials.horizontalHeaderItem(1)
        ___qtablewidgetitem24.setText(QCoreApplication.translate("MainWindow", u"Region", None));
        ___qtablewidgetitem25 = self.TMaterials.horizontalHeaderItem(2)
        ___qtablewidgetitem25.setText(QCoreApplication.translate("MainWindow", u"Type", None));
        ___qtablewidgetitem26 = self.TMaterials.horizontalHeaderItem(3)
        ___qtablewidgetitem26.setText(QCoreApplication.translate("MainWindow", u"C11", None));
        ___qtablewidgetitem27 = self.TMaterials.horizontalHeaderItem(4)
        ___qtablewidgetitem27.setText(QCoreApplication.translate("MainWindow", u"C12", None));
        ___qtablewidgetitem28 = self.TMaterials.horizontalHeaderItem(5)
        ___qtablewidgetitem28.setText(QCoreApplication.translate("MainWindow", u"C13", None));
        ___qtablewidgetitem29 = self.TMaterials.horizontalHeaderItem(6)
        ___qtablewidgetitem29.setText(QCoreApplication.translate("MainWindow", u"C14", None));
        ___qtablewidgetitem30 = self.TMaterials.horizontalHeaderItem(7)
        ___qtablewidgetitem30.setText(QCoreApplication.translate("MainWindow", u"C15", None));
        ___qtablewidgetitem31 = self.TMaterials.horizontalHeaderItem(8)
        ___qtablewidgetitem31.setText(QCoreApplication.translate("MainWindow", u"C16", None));
        ___qtablewidgetitem32 = self.TMaterials.horizontalHeaderItem(9)
        ___qtablewidgetitem32.setText(QCoreApplication.translate("MainWindow", u"C22", None));
        ___qtablewidgetitem33 = self.TMaterials.horizontalHeaderItem(10)
        ___qtablewidgetitem33.setText(QCoreApplication.translate("MainWindow", u"C23", None));
        ___qtablewidgetitem34 = self.TMaterials.horizontalHeaderItem(11)
        ___qtablewidgetitem34.setText(QCoreApplication.translate("MainWindow", u"C24", None));
        ___qtablewidgetitem35 = self.TMaterials.horizontalHeaderItem(12)
        ___qtablewidgetitem35.setText(QCoreApplication.translate("MainWindow", u"C25", None));
        ___qtablewidgetitem36 = self.TMaterials.horizontalHeaderItem(13)
        ___qtablewidgetitem36.setText(QCoreApplication.translate("MainWindow", u"C26", None));
        ___qtablewidgetitem37 = self.TMaterials.horizontalHeaderItem(14)
        ___qtablewidgetitem37.setText(QCoreApplication.translate("MainWindow", u"C33", None));
        ___qtablewidgetitem38 = self.TMaterials.horizontalHeaderItem(15)
        ___qtablewidgetitem38.setText(QCoreApplication.translate("MainWindow", u"C34", None));
        ___qtablewidgetitem39 = self.TMaterials.horizontalHeaderItem(16)
        ___qtablewidgetitem39.setText(QCoreApplication.translate("MainWindow", u"C35", None));
        ___qtablewidgetitem40 = self.TMaterials.horizontalHeaderItem(17)
        ___qtablewidgetitem40.setText(QCoreApplication.translate("MainWindow", u"C36", None));
        ___qtablewidgetitem41 = self.TMaterials.horizontalHeaderItem(18)
        ___qtablewidgetitem41.setText(QCoreApplication.translate("MainWindow", u"C44", None));
        ___qtablewidgetitem42 = self.TMaterials.horizontalHeaderItem(19)
        ___qtablewidgetitem42.setText(QCoreApplication.translate("MainWindow", u"C45", None));
        ___qtablewidgetitem43 = self.TMaterials.horizontalHeaderItem(20)
        ___qtablewidgetitem43.setText(QCoreApplication.translate("MainWindow", u"C46", None));
        ___qtablewidgetitem44 = self.TMaterials.horizontalHeaderItem(21)
        ___qtablewidgetitem44.setText(QCoreApplication.translate("MainWindow", u"C55", None));
        ___qtablewidgetitem45 = self.TMaterials.horizontalHeaderItem(22)
        ___qtablewidgetitem45.setText(QCoreApplication.translate("MainWindow", u"C56", None));
        ___qtablewidgetitem46 = self.TMaterials.horizontalHeaderItem(23)
        ___qtablewidgetitem46.setText(QCoreApplication.translate("MainWindow", u"C66", None));
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
        ___qtablewidgetitem47 = self.TBCs.horizontalHeaderItem(0)
        ___qtablewidgetitem47.setText(QCoreApplication.translate("MainWindow", u"Boundary Condition", None));
        ___qtablewidgetitem48 = self.TBCs.horizontalHeaderItem(1)
        ___qtablewidgetitem48.setText(QCoreApplication.translate("MainWindow", u"Direction", None));
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
        self.BOpenParaview.setText(QCoreApplication.translate("MainWindow", u"Open Paraview", None))
        ___qtablewidgetitem49 = self.THomogenization.horizontalHeaderItem(0)
        ___qtablewidgetitem49.setText(QCoreApplication.translate("MainWindow", u"Homogenized Elastic Constants", None));
        ___qtablewidgetitem50 = self.THomogenization.verticalHeaderItem(0)
        ___qtablewidgetitem50.setText(QCoreApplication.translate("MainWindow", u"E11", None));
        ___qtablewidgetitem51 = self.THomogenization.verticalHeaderItem(1)
        ___qtablewidgetitem51.setText(QCoreApplication.translate("MainWindow", u"E22", None));
        ___qtablewidgetitem52 = self.THomogenization.verticalHeaderItem(2)
        ___qtablewidgetitem52.setText(QCoreApplication.translate("MainWindow", u"E33", None));
        ___qtablewidgetitem53 = self.THomogenization.verticalHeaderItem(3)
        ___qtablewidgetitem53.setText(QCoreApplication.translate("MainWindow", u"NU12", None));
        ___qtablewidgetitem54 = self.THomogenization.verticalHeaderItem(4)
        ___qtablewidgetitem54.setText(QCoreApplication.translate("MainWindow", u"NU21", None));
        ___qtablewidgetitem55 = self.THomogenization.verticalHeaderItem(5)
        ___qtablewidgetitem55.setText(QCoreApplication.translate("MainWindow", u"NU13", None));
        ___qtablewidgetitem56 = self.THomogenization.verticalHeaderItem(6)
        ___qtablewidgetitem56.setText(QCoreApplication.translate("MainWindow", u"NU31", None));
        ___qtablewidgetitem57 = self.THomogenization.verticalHeaderItem(7)
        ___qtablewidgetitem57.setText(QCoreApplication.translate("MainWindow", u"NU23", None));
        ___qtablewidgetitem58 = self.THomogenization.verticalHeaderItem(8)
        ___qtablewidgetitem58.setText(QCoreApplication.translate("MainWindow", u"NU32", None));
        ___qtablewidgetitem59 = self.THomogenization.verticalHeaderItem(9)
        ___qtablewidgetitem59.setText(QCoreApplication.translate("MainWindow", u"G12", None));
        ___qtablewidgetitem60 = self.THomogenization.verticalHeaderItem(10)
        ___qtablewidgetitem60.setText(QCoreApplication.translate("MainWindow", u"G13", None));
        ___qtablewidgetitem61 = self.THomogenization.verticalHeaderItem(11)
        ___qtablewidgetitem61.setText(QCoreApplication.translate("MainWindow", u"G23", None));
        self.BHomogenization.setText(QCoreApplication.translate("MainWindow", u"Generate Homogenized Elastic Constants", None))
        self.LDefineMaterials_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">DEFINE MATERIALS</span></p></body></html>", None))
        self.LEOS.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Equation Of State Model:</p></body></html>", None))
        self.LMaterialNameSGH.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Name:</p></body></html>", None))
        self.INEOS.setItemText(0, QCoreApplication.translate("MainWindow", u"ideal_gas", None))

        self.INMaterialNameSGH.setText("")
        self.LArtificialViscosity.setText(QCoreApplication.translate("MainWindow", u"Artificial Viscosity:", None))
        self.INArtificialViscosity.setItemText(0, QCoreApplication.translate("MainWindow", u"default", None))
        self.INArtificialViscosity.setItemText(1, QCoreApplication.translate("MainWindow", u"custom", None))

        self.BAddMaterialSGH.setText(QCoreApplication.translate("MainWindow", u"Add Material", None))
        ___qtablewidgetitem62 = self.TMaterialsSGH.horizontalHeaderItem(0)
        ___qtablewidgetitem62.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem63 = self.TMaterialsSGH.horizontalHeaderItem(1)
        ___qtablewidgetitem63.setText(QCoreApplication.translate("MainWindow", u"EOS", None));
        ___qtablewidgetitem64 = self.TMaterialsSGH.horizontalHeaderItem(2)
        ___qtablewidgetitem64.setText(QCoreApplication.translate("MainWindow", u"q1", None));
        ___qtablewidgetitem65 = self.TMaterialsSGH.horizontalHeaderItem(3)
        ___qtablewidgetitem65.setText(QCoreApplication.translate("MainWindow", u"q2", None));
        ___qtablewidgetitem66 = self.TMaterialsSGH.horizontalHeaderItem(4)
        ___qtablewidgetitem66.setText(QCoreApplication.translate("MainWindow", u"q1ex", None));
        ___qtablewidgetitem67 = self.TMaterialsSGH.horizontalHeaderItem(5)
        ___qtablewidgetitem67.setText(QCoreApplication.translate("MainWindow", u"q2ex", None));
        ___qtablewidgetitem68 = self.TMaterialsSGH.horizontalHeaderItem(6)
        ___qtablewidgetitem68.setText(QCoreApplication.translate("MainWindow", u"gamma", None));
        ___qtablewidgetitem69 = self.TMaterialsSGH.horizontalHeaderItem(7)
        ___qtablewidgetitem69.setText(QCoreApplication.translate("MainWindow", u"min sound speed", None));
        ___qtablewidgetitem70 = self.TMaterialsSGH.horizontalHeaderItem(8)
        ___qtablewidgetitem70.setText(QCoreApplication.translate("MainWindow", u"specific heat", None));
        self.BDeleteMaterialSGH.setText(QCoreApplication.translate("MainWindow", u"Delete Material", None))
        self.label_12.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">Custom Artificial Viscosity Variables</span></p></body></html>", None))
        self.Lq1.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Accoustic Coefficient in Compression (q1):</p></body></html>", None))
        self.INq1.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.Lq2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Linear Slope in Compression (q2):</p></body></html>", None))
        self.INq2.setText(QCoreApplication.translate("MainWindow", u"1.33333", None))
        self.Lq1ex.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Accoustic Coefficient in Expansion (q1ex):</p></body></html>", None))
        self.INq1ex.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.Lq2ex.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Linear Slope in Expansion (q2ex):</p></body></html>", None))
        self.INq2ex.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.LGamma.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Gamma:</p></body></html>", None))
        self.INGamma.setText(QCoreApplication.translate("MainWindow", u"1.66667", None))
        self.LMinSound.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Minimum Sound Speed:</p></body></html>", None))
        self.INMinSound.setText(QCoreApplication.translate("MainWindow", u"1E-14", None))
        self.LSpecificHeat.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Specific Heat:</p></body></html>", None))
        self.INSpecificHeat.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.LDefineMaterials_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">ASSIGN MATERIALS</span></p></body></html>", None))
        self.LVelx.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Velocity x:</p></body></html>", None))
        self.LRegion_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Material:</p></body></html>", None))
        self.INVelocityZ.setText("")
        self.LMaterialName_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Region</p></body></html>", None))
        self.LVelz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Velocity z:</p></body></html>", None))
        self.LSIE.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Specific Internal Energy:</p></body></html>", None))
        self.LVely.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Velocity y:</p></body></html>", None))
        self.INDensity.setText("")
        self.INSIE.setText("")
        self.LDensity.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Density:</p></body></html>", None))
        self.INVelocityX.setText("")
        self.INVelocityY.setText("")
        self.Baddmaterialassignment.setText(QCoreApplication.translate("MainWindow", u"Add Material Assignment", None))
        ___qtablewidgetitem71 = self.Tassignmat.horizontalHeaderItem(0)
        ___qtablewidgetitem71.setText(QCoreApplication.translate("MainWindow", u"Part", None));
        ___qtablewidgetitem72 = self.Tassignmat.horizontalHeaderItem(1)
        ___qtablewidgetitem72.setText(QCoreApplication.translate("MainWindow", u"Material", None));
        ___qtablewidgetitem73 = self.Tassignmat.horizontalHeaderItem(2)
        ___qtablewidgetitem73.setText(QCoreApplication.translate("MainWindow", u"Density", None));
        ___qtablewidgetitem74 = self.Tassignmat.horizontalHeaderItem(3)
        ___qtablewidgetitem74.setText(QCoreApplication.translate("MainWindow", u"Specific Internal Energy", None));
        ___qtablewidgetitem75 = self.Tassignmat.horizontalHeaderItem(4)
        ___qtablewidgetitem75.setText(QCoreApplication.translate("MainWindow", u"Velocity X", None));
        ___qtablewidgetitem76 = self.Tassignmat.horizontalHeaderItem(5)
        ___qtablewidgetitem76.setText(QCoreApplication.translate("MainWindow", u"Velocity Y", None));
        ___qtablewidgetitem77 = self.Tassignmat.horizontalHeaderItem(6)
        ___qtablewidgetitem77.setText(QCoreApplication.translate("MainWindow", u"Velocity Z", None));
        self.BUpMaterial.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.BDownMaterial.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.Bdeletematerialassignment.setText(QCoreApplication.translate("MainWindow", u"Delete Material Assignment", None))
        self.Lbcs.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline; color:#000000;\">BOUNDARY CONDITIONS</span></p></body></html>", None))
        self.Lbndry.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" color:#000000;\">Boundary: </span></p></body></html>", None))
        self.INBoundary.setItemText(0, QCoreApplication.translate("MainWindow", u"x_plane", None))
        self.INBoundary.setItemText(1, QCoreApplication.translate("MainWindow", u"y_plane", None))
        self.INBoundary.setItemText(2, QCoreApplication.translate("MainWindow", u"z_plane", None))

        self.Lvalue.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" color:#000000;\">Plane Position:</span></p></body></html>", None))
        self.INPlanePosition.setText("")
        self.Ltype.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" color:#000000;\">Type: </span></p></body></html>", None))
        self.INType.setItemText(0, QCoreApplication.translate("MainWindow", u"reflected", None))
        self.INType.setItemText(1, QCoreApplication.translate("MainWindow", u"fixed", None))
        self.INType.setItemText(2, QCoreApplication.translate("MainWindow", u"velocity", None))

        self.LVel0.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>vel_0: </p></body></html>", None))
        self.LVel1.setText(QCoreApplication.translate("MainWindow", u"vel_1: ", None))
        self.LVelstart.setText(QCoreApplication.translate("MainWindow", u"vel_t_start: ", None))
        self.LVelend.setText(QCoreApplication.translate("MainWindow", u"vel_t_end: ", None))
        self.BaddBC.setText(QCoreApplication.translate("MainWindow", u"Add Boundary Condition", None))
        ___qtablewidgetitem78 = self.TBoundaryConditions.horizontalHeaderItem(0)
        ___qtablewidgetitem78.setText(QCoreApplication.translate("MainWindow", u"Boundary", None));
        ___qtablewidgetitem79 = self.TBoundaryConditions.horizontalHeaderItem(1)
        ___qtablewidgetitem79.setText(QCoreApplication.translate("MainWindow", u"Plane Position", None));
        ___qtablewidgetitem80 = self.TBoundaryConditions.horizontalHeaderItem(2)
        ___qtablewidgetitem80.setText(QCoreApplication.translate("MainWindow", u"Type", None));
        ___qtablewidgetitem81 = self.TBoundaryConditions.horizontalHeaderItem(3)
        ___qtablewidgetitem81.setText(QCoreApplication.translate("MainWindow", u"vel_0", None));
        ___qtablewidgetitem82 = self.TBoundaryConditions.horizontalHeaderItem(4)
        ___qtablewidgetitem82.setText(QCoreApplication.translate("MainWindow", u"vel_1", None));
        ___qtablewidgetitem83 = self.TBoundaryConditions.horizontalHeaderItem(5)
        ___qtablewidgetitem83.setText(QCoreApplication.translate("MainWindow", u"vel_t_start", None));
        ___qtablewidgetitem84 = self.TBoundaryConditions.horizontalHeaderItem(6)
        ___qtablewidgetitem84.setText(QCoreApplication.translate("MainWindow", u"vel_t_end", None));
        self.BdeleteBC.setText(QCoreApplication.translate("MainWindow", u"Delete Boundary Condition", None))
        self.label_45.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline; color:#000000;\">SOLVER SETTINGS</span></p></body></html>", None))
        self.Ltime.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" color:#000000;\">Time T: </span></p></body></html>", None))
        self.INTime.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.Lmindt.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" color:#000000;\">Minimum dt: </span></p></body></html>", None))
        self.INMindt.setText(QCoreApplication.translate("MainWindow", u"1E-8", None))
        self.Lmaxdt.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" color:#000000;\">Maximum dt: </span></p></body></html>", None))
        self.INMaxdt.setText(QCoreApplication.translate("MainWindow", u"1E-2", None))
        self.Linitdt.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" color:#000000;\">Initial dt: </span></p></body></html>", None))
        self.INInitialdt.setText(QCoreApplication.translate("MainWindow", u"1E-5", None))
        self.Lmaxcycle.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" color:#000000;\">Maximum # of cycles: </span></p></body></html>", None))
        self.INmaxcycles.setText(QCoreApplication.translate("MainWindow", u"2000000", None))
        self.INGraphicsOutput.setText(QCoreApplication.translate("MainWindow", u"0.25", None))
        self.LGraphicsOutput.setText(QCoreApplication.translate("MainWindow", u"Graphics output step:", None))
        self.LResultsSGH.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline; color:#000000;\">RESULTS</span></p></body></html>", None))
        self.INOuputVarSGH.setItemText(0, QCoreApplication.translate("MainWindow", u"SIE", None))

        self.BPreviewResultsSGH.setText(QCoreApplication.translate("MainWindow", u"Preview Results", None))
        self.BFirstFrame.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.BPreviousFrame.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.BNextFrame.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.BLastFrame.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.LThreshold.setText(QCoreApplication.translate("MainWindow", u"Threshold Value:", None))
        self.INThreshold.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.BThreshold.setText(QCoreApplication.translate("MainWindow", u"threshold", None))
        self.BOpenParaviewSGH.setText(QCoreApplication.translate("MainWindow", u"Open Paraview", None))
        self.menuHelp.setTitle(QCoreApplication.translate("MainWindow", u"Help", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
    # retranslateUi

