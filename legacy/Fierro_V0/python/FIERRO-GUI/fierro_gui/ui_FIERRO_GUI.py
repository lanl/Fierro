# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'FIERRO_GUIHQOjTO.ui'
##
## Created by: Qt User Interface Compiler version 6.7.0
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
from PySide6.QtWidgets import (QAbstractItemView, QAbstractScrollArea, QAbstractSpinBox, QApplication,
    QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout,
    QFrame, QGridLayout, QHBoxLayout, QHeaderView,
    QLabel, QLayout, QLineEdit, QListWidget,
    QListWidgetItem, QMainWindow, QMenu, QMenuBar,
    QPlainTextEdit, QProgressBar, QPushButton, QRadioButton,
    QScrollArea, QSizePolicy, QSpacerItem, QSpinBox,
    QSplitter, QStackedWidget, QStatusBar, QTabWidget,
    QTableWidget, QTableWidgetItem, QToolButton, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout, QWidget)
import IconResourceFile_rc
import IconResourceFile_rc

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.setEnabled(True)
        MainWindow.resize(1326, 1168)
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setMaximumSize(QSize(16777215, 16777215))
        icon = QIcon()
        icon.addFile(u":/Logos/Logos/FIERRO.png", QSize(), QIcon.Normal, QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setAutoFillBackground(True)
        MainWindow.setStyleSheet(u"")
        MainWindow.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        MainWindow.setDockNestingEnabled(False)
        self.actionManual = QAction(MainWindow)
        self.actionManual.setObjectName(u"actionManual")
        self.actionChange_Working_Directory = QAction(MainWindow)
        self.actionChange_Working_Directory.setObjectName(u"actionChange_Working_Directory")
        self.actionSaveAs = QAction(MainWindow)
        self.actionSaveAs.setObjectName(u"actionSaveAs")
        self.actionHelp = QAction(MainWindow)
        self.actionHelp.setObjectName(u"actionHelp")
        self.actionOpen = QAction(MainWindow)
        self.actionOpen.setObjectName(u"actionOpen")
        self.actionNew = QAction(MainWindow)
        self.actionNew.setObjectName(u"actionNew")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.centralwidget.setEnabled(True)
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy1)
        self.centralwidget.setStyleSheet(u"#ParaviewWindow{\n"
"    background-color: rgb(91, 97, 120);\n"
"}\n"
"#BImportPart:hover, #BDefineMaterial:hover, #BApplyBC:hover, #BSolverSettings:hover, #BRunEVPFFT:hover, #BViewResults:hover, #BGlobalMesh:hover, #BImportPartSGH:hover, #BDefineMaterialSGH:hover, #BAssignMaterialSGH:hover, #BApplyBCSGH:hover, #BSolverSettingsSGH:hover, #BViewResultsSGH:hover, #BRunSGH:hover, #BCreateBasicPart:hover{\n"
"    background-color: rgb(192, 192, 192);\n"
"    border-radius: 15px;\n"
"}\n"
"#BImportPart, #BDefineMaterial, #BApplyBC, #BSolverSettings, #BRunEVPFFT, #BViewResults, #BGlobalMesh{\n"
"    border-style: flat;\n"
"}")
        self.verticalLayout_3 = QVBoxLayout(self.centralwidget)
        self.verticalLayout_3.setSpacing(0)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(-1, 9, -1, 9)
        self.NavigationMenu = QTabWidget(self.centralwidget)
        self.NavigationMenu.setObjectName(u"NavigationMenu")
        self.NavigationMenu.setEnabled(True)
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Preferred)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.NavigationMenu.sizePolicy().hasHeightForWidth())
        self.NavigationMenu.setSizePolicy(sizePolicy2)
        self.NavigationMenu.setMaximumSize(QSize(16777215, 16777215))
        font = QFont()
        font.setBold(True)
        self.NavigationMenu.setFont(font)
        self.NavigationMenu.setMouseTracking(False)
        self.NavigationMenu.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.NavigationMenu.setContextMenuPolicy(Qt.ContextMenuPolicy.NoContextMenu)
        self.NavigationMenu.setAcceptDrops(True)
        self.NavigationMenu.setToolTipDuration(-1)
        self.NavigationMenu.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.NavigationMenu.setAutoFillBackground(False)
        self.NavigationMenu.setStyleSheet(u"")
        self.NavigationMenu.setTabPosition(QTabWidget.TabPosition.North)
        self.NavigationMenu.setTabShape(QTabWidget.TabShape.Rounded)
        self.NavigationMenu.setIconSize(QSize(28, 28))
        self.NavigationMenu.setElideMode(Qt.TextElideMode.ElideRight)
        self.NavigationMenu.setUsesScrollButtons(False)
        self.NavigationMenu.setDocumentMode(True)
        self.NavigationMenu.setTabsClosable(False)
        self.NavigationMenu.setMovable(False)
        self.NavigationMenu.setTabBarAutoHide(True)
        self.Title = QWidget()
        self.Title.setObjectName(u"Title")
        icon1 = QIcon()
        icon1.addFile(u":/Logos/Logos/FIERRO_NoText.png", QSize(), QIcon.Normal, QIcon.Off)
        self.NavigationMenu.addTab(self.Title, icon1, "")
        self.Pipeline = QWidget()
        self.Pipeline.setObjectName(u"Pipeline")
        icon2 = QIcon()
        icon2.addFile(u":/Blue Icons/Blue Icons/Clipboard.svg", QSize(), QIcon.Selected, QIcon.On)
        self.NavigationMenu.addTab(self.Pipeline, icon2, "")
        self.Geometry = QWidget()
        self.Geometry.setObjectName(u"Geometry")
        self.Geometry.setEnabled(True)
        sizePolicy3 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.Geometry.sizePolicy().hasHeightForWidth())
        self.Geometry.setSizePolicy(sizePolicy3)
        icon3 = QIcon()
        icon3.addFile(u":/Blue Icons/Blue Icons/Cube.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.NavigationMenu.addTab(self.Geometry, icon3, "")
        self.Mesh = QWidget()
        self.Mesh.setObjectName(u"Mesh")
        icon4 = QIcon()
        icon4.addFile(u":/Blue Icons/Blue Icons/mesh.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.NavigationMenu.addTab(self.Mesh, icon4, "")
        self.Materials = QWidget()
        self.Materials.setObjectName(u"Materials")
        icon5 = QIcon()
        icon5.addFile(u":/Blue Icons/Blue Icons/mine.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.NavigationMenu.addTab(self.Materials, icon5, "")
        self.BoundaryConditions = QWidget()
        self.BoundaryConditions.setObjectName(u"BoundaryConditions")
        icon6 = QIcon()
        icon6.addFile(u":/Blue Icons/Blue Icons/brick.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.NavigationMenu.addTab(self.BoundaryConditions, icon6, "")
        self.Solver = QWidget()
        self.Solver.setObjectName(u"Solver")
        icon7 = QIcon()
        icon7.addFile(u":/Blue Icons/Blue Icons/gear.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.NavigationMenu.addTab(self.Solver, icon7, "")
        self.Run = QWidget()
        self.Run.setObjectName(u"Run")
        icon8 = QIcon()
        icon8.addFile(u":/Blue Icons/Blue Icons/Play.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.NavigationMenu.addTab(self.Run, icon8, "")
        self.Postprocessing = QWidget()
        self.Postprocessing.setObjectName(u"Postprocessing")
        icon9 = QIcon()
        icon9.addFile(u":/Blue Icons/Blue Icons/lens.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.NavigationMenu.addTab(self.Postprocessing, icon9, "")

        self.verticalLayout_3.addWidget(self.NavigationMenu)

        self.frame_27 = QFrame(self.centralwidget)
        self.frame_27.setObjectName(u"frame_27")
        sizePolicy3.setHeightForWidth(self.frame_27.sizePolicy().hasHeightForWidth())
        self.frame_27.setSizePolicy(sizePolicy3)
        self.frame_27.setMinimumSize(QSize(0, 200))
        self.frame_27.setFrameShape(QFrame.Shape.Box)
        self.frame_27.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_10 = QHBoxLayout(self.frame_27)
        self.horizontalLayout_10.setObjectName(u"horizontalLayout_10")
        self.horizontalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.ToolWindow = QStackedWidget(self.frame_27)
        self.ToolWindow.setObjectName(u"ToolWindow")
        self.ToolWindow.setEnabled(True)
        sizePolicy4 = QSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Expanding)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy4.setHeightForWidth(self.ToolWindow.sizePolicy().hasHeightForWidth())
        self.ToolWindow.setSizePolicy(sizePolicy4)
        self.ToolWindow.setMinimumSize(QSize(0, 0))
        self.ToolWindow.setMaximumSize(QSize(400, 16777215))
        self.ToolWindow.setSizeIncrement(QSize(0, 0))
        self.ToolWindow.setBaseSize(QSize(0, 0))
        font1 = QFont()
        font1.setFamilies([u"Hiragino Sans"])
        self.ToolWindow.setFont(font1)
        self.ToolWindow.setAutoFillBackground(False)
        self.ToolWindow.setFrameShape(QFrame.Shape.NoFrame)
        self.ToolWindow.setLineWidth(0)
        self.ToolWindow.setMidLineWidth(1)
        self.TitleTool = QWidget()
        self.TitleTool.setObjectName(u"TitleTool")
        self.verticalLayout_2 = QVBoxLayout(self.TitleTool)
        self.verticalLayout_2.setSpacing(0)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.verticalLayout_2.setContentsMargins(0, 20, 0, 20)
        self.LosAlamosLogo = QLabel(self.TitleTool)
        self.LosAlamosLogo.setObjectName(u"LosAlamosLogo")
        sizePolicy5 = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.MinimumExpanding)
        sizePolicy5.setHorizontalStretch(0)
        sizePolicy5.setVerticalStretch(0)
        sizePolicy5.setHeightForWidth(self.LosAlamosLogo.sizePolicy().hasHeightForWidth())
        self.LosAlamosLogo.setSizePolicy(sizePolicy5)
        self.LosAlamosLogo.setMaximumSize(QSize(376, 74))
        self.LosAlamosLogo.setPixmap(QPixmap(u":/Logos/Logos/LANL Logo Ultramarine.png"))
        self.LosAlamosLogo.setScaledContents(True)

        self.verticalLayout_2.addWidget(self.LosAlamosLogo, 0, Qt.AlignmentFlag.AlignHCenter)

        self.verticalSpacer_7 = QSpacerItem(20, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_2.addItem(self.verticalSpacer_7)

        self.EVPFFTLogo = QLabel(self.TitleTool)
        self.EVPFFTLogo.setObjectName(u"EVPFFTLogo")
        sizePolicy6 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.MinimumExpanding)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(0)
        sizePolicy6.setHeightForWidth(self.EVPFFTLogo.sizePolicy().hasHeightForWidth())
        self.EVPFFTLogo.setSizePolicy(sizePolicy6)
        self.EVPFFTLogo.setMinimumSize(QSize(0, 175))
        self.EVPFFTLogo.setMaximumSize(QSize(225, 175))
        self.EVPFFTLogo.setPixmap(QPixmap(u":/Logos/Logos/FIERRO.png"))
        self.EVPFFTLogo.setScaledContents(True)
        self.EVPFFTLogo.setWordWrap(False)
        self.EVPFFTLogo.setIndent(-1)

        self.verticalLayout_2.addWidget(self.EVPFFTLogo, 0, Qt.AlignmentFlag.AlignHCenter)

        self.verticalSpacer_6 = QSpacerItem(20, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_2.addItem(self.verticalSpacer_6)

        self.AdditionalSoftware = QFrame(self.TitleTool)
        self.AdditionalSoftware.setObjectName(u"AdditionalSoftware")
        self.AdditionalSoftware.setFrameShape(QFrame.Shape.NoFrame)
        self.verticalLayout_10 = QVBoxLayout(self.AdditionalSoftware)
        self.verticalLayout_10.setSpacing(8)
        self.verticalLayout_10.setObjectName(u"verticalLayout_10")
        self.verticalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.LAdditionalSoftware = QLabel(self.AdditionalSoftware)
        self.LAdditionalSoftware.setObjectName(u"LAdditionalSoftware")
        font2 = QFont()
        font2.setPointSize(16)
        self.LAdditionalSoftware.setFont(font2)

        self.verticalLayout_10.addWidget(self.LAdditionalSoftware)

        self.AdditionalSoftwareLogos = QFrame(self.AdditionalSoftware)
        self.AdditionalSoftwareLogos.setObjectName(u"AdditionalSoftwareLogos")
        sizePolicy7 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        sizePolicy7.setHorizontalStretch(0)
        sizePolicy7.setVerticalStretch(0)
        sizePolicy7.setHeightForWidth(self.AdditionalSoftwareLogos.sizePolicy().hasHeightForWidth())
        self.AdditionalSoftwareLogos.setSizePolicy(sizePolicy7)
        self.AdditionalSoftwareLogos.setFrameShape(QFrame.Shape.NoFrame)
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


        self.verticalLayout_2.addWidget(self.AdditionalSoftware)

        self.ToolWindow.addWidget(self.TitleTool)
        self.PipelineTool = QWidget()
        self.PipelineTool.setObjectName(u"PipelineTool")
        self.verticalLayout_20 = QVBoxLayout(self.PipelineTool)
        self.verticalLayout_20.setObjectName(u"verticalLayout_20")
        self.LPipeline = QLabel(self.PipelineTool)
        self.LPipeline.setObjectName(u"LPipeline")
        sizePolicy8 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        sizePolicy8.setHorizontalStretch(0)
        sizePolicy8.setVerticalStretch(0)
        sizePolicy8.setHeightForWidth(self.LPipeline.sizePolicy().hasHeightForWidth())
        self.LPipeline.setSizePolicy(sizePolicy8)
        font3 = QFont()
        font3.setFamilies([u"Sans Serif"])
        font3.setPointSize(16)
        self.LPipeline.setFont(font3)
        self.LPipeline.setMargin(0)

        self.verticalLayout_20.addWidget(self.LPipeline)

        self.frame_30 = QFrame(self.PipelineTool)
        self.frame_30.setObjectName(u"frame_30")
        self.frame_30.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_30.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_23 = QFormLayout(self.frame_30)
        self.formLayout_23.setObjectName(u"formLayout_23")
        self.formLayout_23.setContentsMargins(0, 0, 0, 0)
        self.LPipelineSelection = QLabel(self.frame_30)
        self.LPipelineSelection.setObjectName(u"LPipelineSelection")
        font4 = QFont()
        font4.setFamilies([u".AppleSystemUIFont"])
        self.LPipelineSelection.setFont(font4)
        self.LPipelineSelection.setWordWrap(False)

        self.formLayout_23.setWidget(0, QFormLayout.LabelRole, self.LPipelineSelection)

        self.INPipelineSelection = QComboBox(self.frame_30)
        self.INPipelineSelection.addItem("")
        self.INPipelineSelection.addItem("")
        self.INPipelineSelection.addItem("")
        self.INPipelineSelection.addItem("")
        self.INPipelineSelection.setObjectName(u"INPipelineSelection")
        sizePolicy8.setHeightForWidth(self.INPipelineSelection.sizePolicy().hasHeightForWidth())
        self.INPipelineSelection.setSizePolicy(sizePolicy8)
        font5 = QFont()
        font5.setFamilies([u".AppleSystemUIFont"])
        font5.setPointSize(13)
        self.INPipelineSelection.setFont(font5)
        self.INPipelineSelection.setIconSize(QSize(18, 18))

        self.formLayout_23.setWidget(0, QFormLayout.FieldRole, self.INPipelineSelection)


        self.verticalLayout_20.addWidget(self.frame_30)

        self.PipelineExtras = QStackedWidget(self.PipelineTool)
        self.PipelineExtras.setObjectName(u"PipelineExtras")
        self.page_27 = QWidget()
        self.page_27.setObjectName(u"page_27")
        self.PipelineExtras.addWidget(self.page_27)
        self.page_29 = QWidget()
        self.page_29.setObjectName(u"page_29")
        self.PipelineExtras.addWidget(self.page_29)
        self.page_30 = QWidget()
        self.page_30.setObjectName(u"page_30")
        self.PipelineExtras.addWidget(self.page_30)
        self.page_28 = QWidget()
        self.page_28.setObjectName(u"page_28")
        self.horizontalLayout_24 = QHBoxLayout(self.page_28)
        self.horizontalLayout_24.setObjectName(u"horizontalLayout_24")
        self.BLegacyEVPFFT = QPushButton(self.page_28)
        self.BLegacyEVPFFT.setObjectName(u"BLegacyEVPFFT")

        self.horizontalLayout_24.addWidget(self.BLegacyEVPFFT)

        self.pushButton_11 = QPushButton(self.page_28)
        self.pushButton_11.setObjectName(u"pushButton_11")
        self.pushButton_11.setEnabled(False)
        sizePolicy9 = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        sizePolicy9.setHorizontalStretch(0)
        sizePolicy9.setVerticalStretch(0)
        sizePolicy9.setHeightForWidth(self.pushButton_11.sizePolicy().hasHeightForWidth())
        self.pushButton_11.setSizePolicy(sizePolicy9)
        self.pushButton_11.setMaximumSize(QSize(17, 17))
        font6 = QFont()
        font6.setPointSize(18)
        self.pushButton_11.setFont(font6)
        icon10 = QIcon()
        icon10.addFile(u":/Blue Icons/Blue Icons/information.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.pushButton_11.setIcon(icon10)
        self.pushButton_11.setIconSize(QSize(20, 20))
        self.pushButton_11.setAutoDefault(False)
        self.pushButton_11.setFlat(True)

        self.horizontalLayout_24.addWidget(self.pushButton_11)

        self.PipelineExtras.addWidget(self.page_28)

        self.verticalLayout_20.addWidget(self.PipelineExtras, 0, Qt.AlignmentFlag.AlignTop)

        self.verticalSpacer_2 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_20.addItem(self.verticalSpacer_2)

        self.ToolWindow.addWidget(self.PipelineTool)
        self.ImportGeometryTool = QWidget()
        self.ImportGeometryTool.setObjectName(u"ImportGeometryTool")
        sizePolicy7.setHeightForWidth(self.ImportGeometryTool.sizePolicy().hasHeightForWidth())
        self.ImportGeometryTool.setSizePolicy(sizePolicy7)
        self.verticalLayout_15 = QVBoxLayout(self.ImportGeometryTool)
        self.verticalLayout_15.setObjectName(u"verticalLayout_15")
        self.verticalLayout_15.setSizeConstraint(QLayout.SizeConstraint.SetDefaultConstraint)
        self.verticalLayout_15.setContentsMargins(-1, 0, -1, 0)
        self.Import = QFrame(self.ImportGeometryTool)
        self.Import.setObjectName(u"Import")
        self.Import.setEnabled(True)
        sizePolicy8.setHeightForWidth(self.Import.sizePolicy().hasHeightForWidth())
        self.Import.setSizePolicy(sizePolicy8)
        self.Import.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.Import.setAutoFillBackground(False)
        self.Import.setFrameShape(QFrame.Shape.NoFrame)
        self._2 = QFormLayout(self.Import)
        self._2.setObjectName(u"_2")
        self._2.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        self._2.setRowWrapPolicy(QFormLayout.RowWrapPolicy.DontWrapRows)
        self._2.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        self._2.setHorizontalSpacing(0)
        self._2.setVerticalSpacing(0)
        self._2.setContentsMargins(0, 0, 0, 0)
        self.LGeometryInformation = QLabel(self.Import)
        self.LGeometryInformation.setObjectName(u"LGeometryInformation")
        sizePolicy8.setHeightForWidth(self.LGeometryInformation.sizePolicy().hasHeightForWidth())
        self.LGeometryInformation.setSizePolicy(sizePolicy8)
        self.LGeometryInformation.setFont(font3)
        self.LGeometryInformation.setMargin(0)

        self._2.setWidget(0, QFormLayout.SpanningRole, self.LGeometryInformation)

        self.line_3 = QFrame(self.Import)
        self.line_3.setObjectName(u"line_3")
        self.line_3.setFrameShape(QFrame.Shape.VLine)
        self.line_3.setFrameShadow(QFrame.Shadow.Sunken)

        self._2.setWidget(3, QFormLayout.SpanningRole, self.line_3)


        self.verticalLayout_15.addWidget(self.Import)

        self.INSelectGeometry = QComboBox(self.ImportGeometryTool)
        self.INSelectGeometry.addItem("")
        self.INSelectGeometry.addItem("")
        self.INSelectGeometry.addItem("")
        self.INSelectGeometry.setObjectName(u"INSelectGeometry")

        self.verticalLayout_15.addWidget(self.INSelectGeometry)

        self.frame_31 = QFrame(self.ImportGeometryTool)
        self.frame_31.setObjectName(u"frame_31")
        self.frame_31.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_31.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_24 = QFormLayout(self.frame_31)
        self.formLayout_24.setObjectName(u"formLayout_24")
        self.formLayout_24.setContentsMargins(0, 0, 0, 0)
        self.LPartSelection = QLabel(self.frame_31)
        self.LPartSelection.setObjectName(u"LPartSelection")
        self.LPartSelection.setFont(font4)
        self.LPartSelection.setWordWrap(False)

        self.formLayout_24.setWidget(0, QFormLayout.LabelRole, self.LPartSelection)

        self.INSelectGeometryImport = QComboBox(self.frame_31)
        self.INSelectGeometryImport.addItem("")
        self.INSelectGeometryImport.addItem("")
        self.INSelectGeometryImport.addItem("")
        self.INSelectGeometryImport.addItem("")
        self.INSelectGeometryImport.setObjectName(u"INSelectGeometryImport")
        sizePolicy8.setHeightForWidth(self.INSelectGeometryImport.sizePolicy().hasHeightForWidth())
        self.INSelectGeometryImport.setSizePolicy(sizePolicy8)
        self.INSelectGeometryImport.setFont(font5)
        self.INSelectGeometryImport.setIconSize(QSize(18, 18))

        self.formLayout_24.setWidget(0, QFormLayout.FieldRole, self.INSelectGeometryImport)

        self.LPartName = QLabel(self.frame_31)
        self.LPartName.setObjectName(u"LPartName")
        self.LPartName.setFont(font4)

        self.formLayout_24.setWidget(1, QFormLayout.LabelRole, self.LPartName)

        self.INPartName = QLineEdit(self.frame_31)
        self.INPartName.setObjectName(u"INPartName")
        self.INPartName.setMinimumSize(QSize(0, 0))
        self.INPartName.setFont(font1)

        self.formLayout_24.setWidget(1, QFormLayout.FieldRole, self.INPartName)


        self.verticalLayout_15.addWidget(self.frame_31, 0, Qt.AlignmentFlag.AlignTop)

        self.BUploadGeometryFile = QPushButton(self.ImportGeometryTool)
        self.BUploadGeometryFile.setObjectName(u"BUploadGeometryFile")
        self.BUploadGeometryFile.setFont(font4)

        self.verticalLayout_15.addWidget(self.BUploadGeometryFile)

        self.SAGeometryScrollArea = QScrollArea(self.ImportGeometryTool)
        self.SAGeometryScrollArea.setObjectName(u"SAGeometryScrollArea")
        self.SAGeometryScrollArea.setEnabled(True)
        sizePolicy3.setHeightForWidth(self.SAGeometryScrollArea.sizePolicy().hasHeightForWidth())
        self.SAGeometryScrollArea.setSizePolicy(sizePolicy3)
        self.SAGeometryScrollArea.setMinimumSize(QSize(0, 0))
        self.SAGeometryScrollArea.setMaximumSize(QSize(16777215, 16777215))
        self.SAGeometryScrollArea.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.SAGeometryScrollArea.setFrameShape(QFrame.Shape.NoFrame)
        self.SAGeometryScrollArea.setLineWidth(1)
        self.SAGeometryScrollArea.setMidLineWidth(1)
        self.SAGeometryScrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.SAGeometryScrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.SAGeometryScrollArea.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustToContents)
        self.SAGeometryScrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 376, 548))
        sizePolicy3.setHeightForWidth(self.scrollAreaWidgetContents.sizePolicy().hasHeightForWidth())
        self.scrollAreaWidgetContents.setSizePolicy(sizePolicy3)
        self.scrollAreaWidgetContents.setMinimumSize(QSize(0, 0))
        self.scrollAreaWidgetContents.setAcceptDrops(True)
        self.verticalLayout_40 = QVBoxLayout(self.scrollAreaWidgetContents)
        self.verticalLayout_40.setSpacing(0)
        self.verticalLayout_40.setObjectName(u"verticalLayout_40")
        self.verticalLayout_40.setSizeConstraint(QLayout.SizeConstraint.SetDefaultConstraint)
        self.verticalLayout_40.setContentsMargins(0, 0, 0, 0)
        self.GeometryOptions = QStackedWidget(self.scrollAreaWidgetContents)
        self.GeometryOptions.setObjectName(u"GeometryOptions")
        sizePolicy7.setHeightForWidth(self.GeometryOptions.sizePolicy().hasHeightForWidth())
        self.GeometryOptions.setSizePolicy(sizePolicy7)
        self.GeometryOptions.setMinimumSize(QSize(0, 0))
        self.GeometryOptions.setMouseTracking(False)
        self.GeometryOptions.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.GeometryOptions.setAutoFillBackground(False)
        self.page_14 = QWidget()
        self.page_14.setObjectName(u"page_14")
        self.GeometryOptions.addWidget(self.page_14)
        self.ImportPartTool = QWidget()
        self.ImportPartTool.setObjectName(u"ImportPartTool")
        self.verticalLayout_8 = QVBoxLayout(self.ImportPartTool)
        self.verticalLayout_8.setObjectName(u"verticalLayout_8")
        self.LSTLVoxelization = QLabel(self.ImportPartTool)
        self.LSTLVoxelization.setObjectName(u"LSTLVoxelization")
        font7 = QFont()
        font7.setPointSize(15)
        font7.setWeight(QFont.DemiBold)
        font7.setUnderline(False)
        font7.setStrikeOut(False)
        self.LSTLVoxelization.setFont(font7)

        self.verticalLayout_8.addWidget(self.LSTLVoxelization)

        self.LVoxelCount = QLabel(self.ImportPartTool)
        self.LVoxelCount.setObjectName(u"LVoxelCount")
        font8 = QFont()
        font8.setWeight(QFont.Medium)
        self.LVoxelCount.setFont(font8)
        self.LVoxelCount.setMargin(0)

        self.verticalLayout_8.addWidget(self.LVoxelCount)

        self.frame_8 = QFrame(self.ImportPartTool)
        self.frame_8.setObjectName(u"frame_8")
        self.frame_8.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_8.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_17 = QGridLayout(self.frame_8)
        self.gridLayout_17.setObjectName(u"gridLayout_17")
        self.gridLayout_17.setHorizontalSpacing(-1)
        self.gridLayout_17.setContentsMargins(0, 0, 0, 0)
        self.LNumberOfVoxelsX = QLabel(self.frame_8)
        self.LNumberOfVoxelsX.setObjectName(u"LNumberOfVoxelsX")
        self.LNumberOfVoxelsX.setEnabled(True)
        font9 = QFont()
        font9.setFamilies([u".AppleSystemUIFont"])
        font9.setPointSize(13)
        font9.setBold(False)
        font9.setItalic(False)
        font9.setKerning(True)
        self.LNumberOfVoxelsX.setFont(font9)
        self.LNumberOfVoxelsX.setTextFormat(Qt.TextFormat.PlainText)

        self.gridLayout_17.addWidget(self.LNumberOfVoxelsX, 0, 0, 1, 1)

        self.INNumberOfVoxelsX = QLineEdit(self.frame_8)
        self.INNumberOfVoxelsX.setObjectName(u"INNumberOfVoxelsX")
        self.INNumberOfVoxelsX.setEnabled(True)

        self.gridLayout_17.addWidget(self.INNumberOfVoxelsX, 0, 1, 1, 1)

        self.LNumberOfVoxelsY = QLabel(self.frame_8)
        self.LNumberOfVoxelsY.setObjectName(u"LNumberOfVoxelsY")
        self.LNumberOfVoxelsY.setEnabled(True)
        self.LNumberOfVoxelsY.setFont(font5)

        self.gridLayout_17.addWidget(self.LNumberOfVoxelsY, 1, 0, 1, 1)

        self.INNumberOfVoxelsY = QLineEdit(self.frame_8)
        self.INNumberOfVoxelsY.setObjectName(u"INNumberOfVoxelsY")
        self.INNumberOfVoxelsY.setEnabled(True)

        self.gridLayout_17.addWidget(self.INNumberOfVoxelsY, 1, 1, 1, 1)

        self.LNumberOfVoxelsZ = QLabel(self.frame_8)
        self.LNumberOfVoxelsZ.setObjectName(u"LNumberOfVoxelsZ")
        self.LNumberOfVoxelsZ.setEnabled(True)
        self.LNumberOfVoxelsZ.setFont(font5)

        self.gridLayout_17.addWidget(self.LNumberOfVoxelsZ, 2, 0, 1, 1)

        self.INNumberOfVoxelsZ = QLineEdit(self.frame_8)
        self.INNumberOfVoxelsZ.setObjectName(u"INNumberOfVoxelsZ")
        self.INNumberOfVoxelsZ.setEnabled(True)

        self.gridLayout_17.addWidget(self.INNumberOfVoxelsZ, 2, 1, 1, 1)


        self.verticalLayout_8.addWidget(self.frame_8)

        self.BStlDimensions = QRadioButton(self.ImportPartTool)
        self.BStlDimensions.setObjectName(u"BStlDimensions")
        self.BStlDimensions.setEnabled(True)
        self.BStlDimensions.setFont(font5)
        self.BStlDimensions.setChecked(True)

        self.verticalLayout_8.addWidget(self.BStlDimensions)

        self.BCustomDimensions = QRadioButton(self.ImportPartTool)
        self.BCustomDimensions.setObjectName(u"BCustomDimensions")
        self.BCustomDimensions.setEnabled(True)
        sizePolicy10 = QSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        sizePolicy10.setHorizontalStretch(0)
        sizePolicy10.setVerticalStretch(0)
        sizePolicy10.setHeightForWidth(self.BCustomDimensions.sizePolicy().hasHeightForWidth())
        self.BCustomDimensions.setSizePolicy(sizePolicy10)
        self.BCustomDimensions.setFont(font5)

        self.verticalLayout_8.addWidget(self.BCustomDimensions)

        self.LOriginPoint = QLabel(self.ImportPartTool)
        self.LOriginPoint.setObjectName(u"LOriginPoint")
        self.LOriginPoint.setFont(font8)
        self.LOriginPoint.setMargin(0)

        self.verticalLayout_8.addWidget(self.LOriginPoint)

        self.frame_9 = QFrame(self.ImportPartTool)
        self.frame_9.setObjectName(u"frame_9")
        self.frame_9.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_9.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_10 = QGridLayout(self.frame_9)
        self.gridLayout_10.setObjectName(u"gridLayout_10")
        self.gridLayout_10.setContentsMargins(0, 0, 0, 0)
        self.LOriginX = QLabel(self.frame_9)
        self.LOriginX.setObjectName(u"LOriginX")
        self.LOriginX.setEnabled(False)
        self.LOriginX.setFont(font5)

        self.gridLayout_10.addWidget(self.LOriginX, 0, 0, 1, 1)

        self.INOriginY = QLineEdit(self.frame_9)
        self.INOriginY.setObjectName(u"INOriginY")
        self.INOriginY.setEnabled(False)

        self.gridLayout_10.addWidget(self.INOriginY, 1, 1, 1, 1)

        self.INOriginZ = QLineEdit(self.frame_9)
        self.INOriginZ.setObjectName(u"INOriginZ")
        self.INOriginZ.setEnabled(False)

        self.gridLayout_10.addWidget(self.INOriginZ, 2, 1, 1, 1)

        self.INOriginX = QLineEdit(self.frame_9)
        self.INOriginX.setObjectName(u"INOriginX")
        self.INOriginX.setEnabled(False)

        self.gridLayout_10.addWidget(self.INOriginX, 0, 1, 1, 1)

        self.LOriginZ = QLabel(self.frame_9)
        self.LOriginZ.setObjectName(u"LOriginZ")
        self.LOriginZ.setEnabled(False)
        self.LOriginZ.setFont(font5)

        self.gridLayout_10.addWidget(self.LOriginZ, 2, 0, 1, 1)

        self.LOriginY = QLabel(self.frame_9)
        self.LOriginY.setObjectName(u"LOriginY")
        self.LOriginY.setEnabled(False)
        self.LOriginY.setFont(font5)

        self.gridLayout_10.addWidget(self.LOriginY, 1, 0, 1, 1)

        self.Ulength7 = QLabel(self.frame_9)
        self.Ulength7.setObjectName(u"Ulength7")

        self.gridLayout_10.addWidget(self.Ulength7, 0, 2, 1, 1)

        self.Ulength8 = QLabel(self.frame_9)
        self.Ulength8.setObjectName(u"Ulength8")

        self.gridLayout_10.addWidget(self.Ulength8, 1, 2, 1, 1)

        self.Ulength9 = QLabel(self.frame_9)
        self.Ulength9.setObjectName(u"Ulength9")

        self.gridLayout_10.addWidget(self.Ulength9, 2, 2, 1, 1)


        self.verticalLayout_8.addWidget(self.frame_9)

        self.LDimensions = QLabel(self.ImportPartTool)
        self.LDimensions.setObjectName(u"LDimensions")
        self.LDimensions.setFont(font8)
        self.LDimensions.setMargin(2)

        self.verticalLayout_8.addWidget(self.LDimensions)

        self.frame_18 = QFrame(self.ImportPartTool)
        self.frame_18.setObjectName(u"frame_18")
        self.frame_18.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_18.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_12 = QGridLayout(self.frame_18)
        self.gridLayout_12.setObjectName(u"gridLayout_12")
        self.gridLayout_12.setContentsMargins(0, 0, 0, 0)
        self.INLengthX = QLineEdit(self.frame_18)
        self.INLengthX.setObjectName(u"INLengthX")
        self.INLengthX.setEnabled(False)
        self.INLengthX.setMinimumSize(QSize(208, 0))

        self.gridLayout_12.addWidget(self.INLengthX, 0, 1, 1, 1)

        self.INLengthZ = QLineEdit(self.frame_18)
        self.INLengthZ.setObjectName(u"INLengthZ")
        self.INLengthZ.setEnabled(False)
        self.INLengthZ.setMinimumSize(QSize(208, 0))

        self.gridLayout_12.addWidget(self.INLengthZ, 2, 1, 1, 1)

        self.LLengthY = QLabel(self.frame_18)
        self.LLengthY.setObjectName(u"LLengthY")
        self.LLengthY.setEnabled(False)
        self.LLengthY.setFont(font5)

        self.gridLayout_12.addWidget(self.LLengthY, 1, 0, 1, 1)

        self.LLengthZ = QLabel(self.frame_18)
        self.LLengthZ.setObjectName(u"LLengthZ")
        self.LLengthZ.setEnabled(False)
        self.LLengthZ.setFont(font5)

        self.gridLayout_12.addWidget(self.LLengthZ, 2, 0, 1, 1)

        self.LLengthX = QLabel(self.frame_18)
        self.LLengthX.setObjectName(u"LLengthX")
        self.LLengthX.setEnabled(False)
        self.LLengthX.setFont(font5)

        self.gridLayout_12.addWidget(self.LLengthX, 0, 0, 1, 1)

        self.INLengthY = QLineEdit(self.frame_18)
        self.INLengthY.setObjectName(u"INLengthY")
        self.INLengthY.setEnabled(False)
        self.INLengthY.setMinimumSize(QSize(208, 0))

        self.gridLayout_12.addWidget(self.INLengthY, 1, 1, 1, 1)

        self.Ulength10 = QLabel(self.frame_18)
        self.Ulength10.setObjectName(u"Ulength10")

        self.gridLayout_12.addWidget(self.Ulength10, 0, 2, 1, 1)

        self.Ulength11 = QLabel(self.frame_18)
        self.Ulength11.setObjectName(u"Ulength11")

        self.gridLayout_12.addWidget(self.Ulength11, 1, 2, 1, 1)

        self.Ulength12 = QLabel(self.frame_18)
        self.Ulength12.setObjectName(u"Ulength12")

        self.gridLayout_12.addWidget(self.Ulength12, 2, 2, 1, 1)


        self.verticalLayout_8.addWidget(self.frame_18)

        self.BVoxelizeGeometry = QPushButton(self.ImportPartTool)
        self.BVoxelizeGeometry.setObjectName(u"BVoxelizeGeometry")
        self.BVoxelizeGeometry.setEnabled(True)
        self.BVoxelizeGeometry.setFont(font5)

        self.verticalLayout_8.addWidget(self.BVoxelizeGeometry)

        self.verticalSpacer_17 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_8.addItem(self.verticalSpacer_17)

        self.GeometryOptions.addWidget(self.ImportPartTool)
        self.ImportImageStackTool = QWidget()
        self.ImportImageStackTool.setObjectName(u"ImportImageStackTool")
        self.verticalLayout_16 = QVBoxLayout(self.ImportImageStackTool)
        self.verticalLayout_16.setObjectName(u"verticalLayout_16")
        self.LImageStack = QLabel(self.ImportImageStackTool)
        self.LImageStack.setObjectName(u"LImageStack")
        font10 = QFont()
        font10.setPointSize(15)
        font10.setBold(True)
        font10.setUnderline(False)
        self.LImageStack.setFont(font10)
        self.LImageStack.setLayoutDirection(Qt.LayoutDirection.LeftToRight)

        self.verticalLayout_16.addWidget(self.LImageStack, 0, Qt.AlignmentFlag.AlignTop)

        self.frame_29 = QFrame(self.ImportImageStackTool)
        self.frame_29.setObjectName(u"frame_29")
        self.frame_29.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_29.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_22 = QFormLayout(self.frame_29)
        self.formLayout_22.setObjectName(u"formLayout_22")
        self.formLayout_22.setContentsMargins(0, 0, 0, 0)
        self.LUploadedDirectory = QLabel(self.frame_29)
        self.LUploadedDirectory.setObjectName(u"LUploadedDirectory")
        font11 = QFont()
        font11.setPointSize(12)
        self.LUploadedDirectory.setFont(font11)
        self.LUploadedDirectory.setTextFormat(Qt.TextFormat.PlainText)

        self.formLayout_22.setWidget(0, QFormLayout.LabelRole, self.LUploadedDirectory)

        self.INDirectory = QLineEdit(self.frame_29)
        self.INDirectory.setObjectName(u"INDirectory")
        sizePolicy11 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        sizePolicy11.setHorizontalStretch(0)
        sizePolicy11.setVerticalStretch(0)
        sizePolicy11.setHeightForWidth(self.INDirectory.sizePolicy().hasHeightForWidth())
        self.INDirectory.setSizePolicy(sizePolicy11)

        self.formLayout_22.setWidget(0, QFormLayout.FieldRole, self.INDirectory)

        self.LImageFileFormat = QLabel(self.frame_29)
        self.LImageFileFormat.setObjectName(u"LImageFileFormat")
        self.LImageFileFormat.setFont(font11)
        self.LImageFileFormat.setTextFormat(Qt.TextFormat.PlainText)

        self.formLayout_22.setWidget(1, QFormLayout.LabelRole, self.LImageFileFormat)

        self.INFileFormat = QLineEdit(self.frame_29)
        self.INFileFormat.setObjectName(u"INFileFormat")
        sizePolicy11.setHeightForWidth(self.INFileFormat.sizePolicy().hasHeightForWidth())
        self.INFileFormat.setSizePolicy(sizePolicy11)

        self.formLayout_22.setWidget(1, QFormLayout.FieldRole, self.INFileFormat)


        self.verticalLayout_16.addWidget(self.frame_29)

        self.frame_4 = QFrame(self.ImportImageStackTool)
        self.frame_4.setObjectName(u"frame_4")
        self.frame_4.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_4.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_3 = QHBoxLayout(self.frame_4)
#ifndef Q_OS_MAC
        self.horizontalLayout_3.setSpacing(-1)
#endif
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.BImageToVTK = QPushButton(self.frame_4)
        self.BImageToVTK.setObjectName(u"BImageToVTK")
        self.BImageToVTK.setEnabled(False)
        sizePolicy7.setHeightForWidth(self.BImageToVTK.sizePolicy().hasHeightForWidth())
        self.BImageToVTK.setSizePolicy(sizePolicy7)

        self.horizontalLayout_3.addWidget(self.BImageToVTK)

        self.BTiffToStl = QPushButton(self.frame_4)
        self.BTiffToStl.setObjectName(u"BTiffToStl")
        self.BTiffToStl.setEnabled(False)
        sizePolicy7.setHeightForWidth(self.BTiffToStl.sizePolicy().hasHeightForWidth())
        self.BTiffToStl.setSizePolicy(sizePolicy7)

        self.horizontalLayout_3.addWidget(self.BTiffToStl)


        self.verticalLayout_16.addWidget(self.frame_4)

        self.verticalSpacer_3 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_16.addItem(self.verticalSpacer_3)

        self.GeometryOptions.addWidget(self.ImportImageStackTool)
        self.page_13 = QWidget()
        self.page_13.setObjectName(u"page_13")
        self.verticalLayout_27 = QVBoxLayout(self.page_13)
        self.verticalLayout_27.setObjectName(u"verticalLayout_27")
        self.LImageStack_2 = QLabel(self.page_13)
        self.LImageStack_2.setObjectName(u"LImageStack_2")
        self.LImageStack_2.setFont(font10)
        self.LImageStack_2.setLayoutDirection(Qt.LayoutDirection.LeftToRight)

        self.verticalLayout_27.addWidget(self.LImageStack_2, 0, Qt.AlignmentFlag.AlignTop)

        self.frame_20 = QFrame(self.page_13)
        self.frame_20.setObjectName(u"frame_20")
        self.frame_20.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_20.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_10 = QFormLayout(self.frame_20)
        self.formLayout_10.setObjectName(u"formLayout_10")
        self.formLayout_10.setContentsMargins(0, 0, 0, 0)
        self.LColorBy = QLabel(self.frame_20)
        self.LColorBy.setObjectName(u"LColorBy")
        sizePolicy12 = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
        sizePolicy12.setHorizontalStretch(0)
        sizePolicy12.setVerticalStretch(0)
        sizePolicy12.setHeightForWidth(self.LColorBy.sizePolicy().hasHeightForWidth())
        self.LColorBy.setSizePolicy(sizePolicy12)

        self.formLayout_10.setWidget(0, QFormLayout.LabelRole, self.LColorBy)

        self.INSelectColorBy = QComboBox(self.frame_20)
        self.INSelectColorBy.setObjectName(u"INSelectColorBy")
        self.INSelectColorBy.setEnabled(True)

        self.formLayout_10.setWidget(0, QFormLayout.FieldRole, self.INSelectColorBy)


        self.verticalLayout_27.addWidget(self.frame_20, 0, Qt.AlignmentFlag.AlignTop)

        self.verticalSpacer_4 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_27.addItem(self.verticalSpacer_4)

        self.GeometryOptions.addWidget(self.page_13)
        self.page = QWidget()
        self.page.setObjectName(u"page")
        self.verticalLayout_29 = QVBoxLayout(self.page)
        self.verticalLayout_29.setObjectName(u"verticalLayout_29")
        self.LBasicGeometry = QLabel(self.page)
        self.LBasicGeometry.setObjectName(u"LBasicGeometry")

        self.verticalLayout_29.addWidget(self.LBasicGeometry)

        self.frame_6 = QFrame(self.page)
        self.frame_6.setObjectName(u"frame_6")
        self.frame_6.setFrameShape(QFrame.Shape.NoFrame)
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

        self.BasicGeometries = QStackedWidget(self.page)
        self.BasicGeometries.setObjectName(u"BasicGeometries")
        self.BoxProperties = QWidget()
        self.BoxProperties.setObjectName(u"BoxProperties")
        self.verticalLayout_35 = QVBoxLayout(self.BoxProperties)
        self.verticalLayout_35.setObjectName(u"verticalLayout_35")
        self.verticalLayout_35.setContentsMargins(-1, 12, -1, -1)
        self.LBoxProperties = QLabel(self.BoxProperties)
        self.LBoxProperties.setObjectName(u"LBoxProperties")

        self.verticalLayout_35.addWidget(self.LBoxProperties)

        self.frame_10 = QFrame(self.BoxProperties)
        self.frame_10.setObjectName(u"frame_10")
        self.frame_10.setFrameShape(QFrame.Shape.NoFrame)
        self.formLayout_15 = QFormLayout(self.frame_10)
        self.formLayout_15.setObjectName(u"formLayout_15")
        self.formLayout_15.setContentsMargins(0, 0, 0, 0)
        self.LBoxx1 = QLabel(self.frame_10)
        self.LBoxx1.setObjectName(u"LBoxx1")

        self.formLayout_15.setWidget(1, QFormLayout.LabelRole, self.LBoxx1)

        self.INBoxx1 = QLineEdit(self.frame_10)
        self.INBoxx1.setObjectName(u"INBoxx1")

        self.formLayout_15.setWidget(1, QFormLayout.FieldRole, self.INBoxx1)

        self.LBoxx2 = QLabel(self.frame_10)
        self.LBoxx2.setObjectName(u"LBoxx2")

        self.formLayout_15.setWidget(2, QFormLayout.LabelRole, self.LBoxx2)

        self.INBoxx2 = QLineEdit(self.frame_10)
        self.INBoxx2.setObjectName(u"INBoxx2")

        self.formLayout_15.setWidget(2, QFormLayout.FieldRole, self.INBoxx2)

        self.LBoxy1 = QLabel(self.frame_10)
        self.LBoxy1.setObjectName(u"LBoxy1")

        self.formLayout_15.setWidget(3, QFormLayout.LabelRole, self.LBoxy1)

        self.INBoxy1 = QLineEdit(self.frame_10)
        self.INBoxy1.setObjectName(u"INBoxy1")

        self.formLayout_15.setWidget(3, QFormLayout.FieldRole, self.INBoxy1)

        self.LBoxy2 = QLabel(self.frame_10)
        self.LBoxy2.setObjectName(u"LBoxy2")

        self.formLayout_15.setWidget(4, QFormLayout.LabelRole, self.LBoxy2)

        self.INBoxy2 = QLineEdit(self.frame_10)
        self.INBoxy2.setObjectName(u"INBoxy2")

        self.formLayout_15.setWidget(4, QFormLayout.FieldRole, self.INBoxy2)

        self.LBoxz1 = QLabel(self.frame_10)
        self.LBoxz1.setObjectName(u"LBoxz1")

        self.formLayout_15.setWidget(5, QFormLayout.LabelRole, self.LBoxz1)

        self.INBoxz1 = QLineEdit(self.frame_10)
        self.INBoxz1.setObjectName(u"INBoxz1")

        self.formLayout_15.setWidget(5, QFormLayout.FieldRole, self.INBoxz1)

        self.LBoxz2 = QLabel(self.frame_10)
        self.LBoxz2.setObjectName(u"LBoxz2")

        self.formLayout_15.setWidget(6, QFormLayout.LabelRole, self.LBoxz2)

        self.INBoxz2 = QLineEdit(self.frame_10)
        self.INBoxz2.setObjectName(u"INBoxz2")

        self.formLayout_15.setWidget(6, QFormLayout.FieldRole, self.INBoxz2)


        self.verticalLayout_35.addWidget(self.frame_10, 0, Qt.AlignmentFlag.AlignTop)

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
        self.frame_11.setFrameShape(QFrame.Shape.NoFrame)
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

        self.verticalSpacer_15 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

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
        self.frame_12.setFrameShape(QFrame.Shape.NoFrame)
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

        self.verticalSpacer_16 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_37.addItem(self.verticalSpacer_16)

        self.BasicGeometries.addWidget(self.CylinderProperties)

        self.verticalLayout_29.addWidget(self.BasicGeometries)

        self.BGenerateBasicGeometry = QPushButton(self.page)
        self.BGenerateBasicGeometry.setObjectName(u"BGenerateBasicGeometry")

        self.verticalLayout_29.addWidget(self.BGenerateBasicGeometry)

        self.TBasicGeometries = QTableWidget(self.page)
        if (self.TBasicGeometries.columnCount() < 13):
            self.TBasicGeometries.setColumnCount(13)
        __qtablewidgetitem = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        __qtablewidgetitem2 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(2, __qtablewidgetitem2)
        __qtablewidgetitem3 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(3, __qtablewidgetitem3)
        __qtablewidgetitem4 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(4, __qtablewidgetitem4)
        __qtablewidgetitem5 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(5, __qtablewidgetitem5)
        __qtablewidgetitem6 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(6, __qtablewidgetitem6)
        __qtablewidgetitem7 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(7, __qtablewidgetitem7)
        __qtablewidgetitem8 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(8, __qtablewidgetitem8)
        __qtablewidgetitem9 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(9, __qtablewidgetitem9)
        __qtablewidgetitem10 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(10, __qtablewidgetitem10)
        __qtablewidgetitem11 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(11, __qtablewidgetitem11)
        __qtablewidgetitem12 = QTableWidgetItem()
        self.TBasicGeometries.setHorizontalHeaderItem(12, __qtablewidgetitem12)
        self.TBasicGeometries.setObjectName(u"TBasicGeometries")

        self.verticalLayout_29.addWidget(self.TBasicGeometries)

        self.BDeleteBasicGeometry = QPushButton(self.page)
        self.BDeleteBasicGeometry.setObjectName(u"BDeleteBasicGeometry")

        self.verticalLayout_29.addWidget(self.BDeleteBasicGeometry)

        self.verticalSpacer_5 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_29.addItem(self.verticalSpacer_5)

        self.GeometryOptions.addWidget(self.page)
        self.page_8 = QWidget()
        self.page_8.setObjectName(u"page_8")
        sizePolicy3.setHeightForWidth(self.page_8.sizePolicy().hasHeightForWidth())
        self.page_8.setSizePolicy(sizePolicy3)
        self.verticalLayout_14 = QVBoxLayout(self.page_8)
        self.verticalLayout_14.setObjectName(u"verticalLayout_14")
        self.label = QLabel(self.page_8)
        self.label.setObjectName(u"label")

        self.verticalLayout_14.addWidget(self.label)

        self.BVTKFileProperties = QRadioButton(self.page_8)
        self.BVTKFileProperties.setObjectName(u"BVTKFileProperties")
        self.BVTKFileProperties.setChecked(True)

        self.verticalLayout_14.addWidget(self.BVTKFileProperties)

        self.BVTKCustomProperties = QRadioButton(self.page_8)
        self.BVTKCustomProperties.setObjectName(u"BVTKCustomProperties")

        self.verticalLayout_14.addWidget(self.BVTKCustomProperties)

        self.label_7 = QLabel(self.page_8)
        self.label_7.setObjectName(u"label_7")
        self.label_7.setEnabled(False)
        self.label_7.setFont(font)

        self.verticalLayout_14.addWidget(self.label_7)

        self.frame_22 = QFrame(self.page_8)
        self.frame_22.setObjectName(u"frame_22")
        self.frame_22.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_22.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_16 = QGridLayout(self.frame_22)
        self.gridLayout_16.setObjectName(u"gridLayout_16")
        self.gridLayout_16.setContentsMargins(0, 0, 0, 0)
        self.INvtkvy = QLineEdit(self.frame_22)
        self.INvtkvy.setObjectName(u"INvtkvy")
        self.INvtkvy.setEnabled(False)

        self.gridLayout_16.addWidget(self.INvtkvy, 1, 1, 1, 1)

        self.Lvtkvy = QLabel(self.frame_22)
        self.Lvtkvy.setObjectName(u"Lvtkvy")
        self.Lvtkvy.setEnabled(False)

        self.gridLayout_16.addWidget(self.Lvtkvy, 1, 0, 1, 1)

        self.Lvtkvx = QLabel(self.frame_22)
        self.Lvtkvx.setObjectName(u"Lvtkvx")
        self.Lvtkvx.setEnabled(False)

        self.gridLayout_16.addWidget(self.Lvtkvx, 0, 0, 1, 1)

        self.INvtkvx = QLineEdit(self.frame_22)
        self.INvtkvx.setObjectName(u"INvtkvx")
        self.INvtkvx.setEnabled(False)

        self.gridLayout_16.addWidget(self.INvtkvx, 0, 1, 1, 1)

        self.Lvtkvz = QLabel(self.frame_22)
        self.Lvtkvz.setObjectName(u"Lvtkvz")
        self.Lvtkvz.setEnabled(False)

        self.gridLayout_16.addWidget(self.Lvtkvz, 2, 0, 1, 1)

        self.INvtkvz = QLineEdit(self.frame_22)
        self.INvtkvz.setObjectName(u"INvtkvz")
        self.INvtkvz.setEnabled(False)

        self.gridLayout_16.addWidget(self.INvtkvz, 2, 1, 1, 1)


        self.verticalLayout_14.addWidget(self.frame_22)

        self.label_4 = QLabel(self.page_8)
        self.label_4.setObjectName(u"label_4")

        self.verticalLayout_14.addWidget(self.label_4)

        self.frame_3 = QFrame(self.page_8)
        self.frame_3.setObjectName(u"frame_3")
        self.frame_3.setMinimumSize(QSize(0, 0))
        self.frame_3.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_3.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_2 = QGridLayout(self.frame_3)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.Lvoz = QLabel(self.frame_3)
        self.Lvoz.setObjectName(u"Lvoz")
        self.Lvoz.setEnabled(False)

        self.gridLayout_2.addWidget(self.Lvoz, 2, 0, 1, 1)

        self.Lvoy = QLabel(self.frame_3)
        self.Lvoy.setObjectName(u"Lvoy")
        self.Lvoy.setEnabled(False)

        self.gridLayout_2.addWidget(self.Lvoy, 1, 0, 1, 1)

        self.INvox = QLineEdit(self.frame_3)
        self.INvox.setObjectName(u"INvox")
        self.INvox.setEnabled(False)

        self.gridLayout_2.addWidget(self.INvox, 0, 1, 1, 1)

        self.INvoy = QLineEdit(self.frame_3)
        self.INvoy.setObjectName(u"INvoy")
        self.INvoy.setEnabled(False)

        self.gridLayout_2.addWidget(self.INvoy, 1, 1, 1, 1)

        self.INvoz = QLineEdit(self.frame_3)
        self.INvoz.setObjectName(u"INvoz")
        self.INvoz.setEnabled(False)

        self.gridLayout_2.addWidget(self.INvoz, 2, 1, 1, 1)

        self.Lvox = QLabel(self.frame_3)
        self.Lvox.setObjectName(u"Lvox")
        self.Lvox.setEnabled(False)

        self.gridLayout_2.addWidget(self.Lvox, 0, 0, 1, 1)

        self.Ulength1 = QLabel(self.frame_3)
        self.Ulength1.setObjectName(u"Ulength1")

        self.gridLayout_2.addWidget(self.Ulength1, 0, 2, 1, 1)

        self.Ulength2 = QLabel(self.frame_3)
        self.Ulength2.setObjectName(u"Ulength2")

        self.gridLayout_2.addWidget(self.Ulength2, 1, 2, 1, 1)

        self.Ulength3 = QLabel(self.frame_3)
        self.Ulength3.setObjectName(u"Ulength3")

        self.gridLayout_2.addWidget(self.Ulength3, 2, 2, 1, 1)


        self.verticalLayout_14.addWidget(self.frame_3)

        self.label_27 = QLabel(self.page_8)
        self.label_27.setObjectName(u"label_27")

        self.verticalLayout_14.addWidget(self.label_27)

        self.frame_7 = QFrame(self.page_8)
        self.frame_7.setObjectName(u"frame_7")
        self.frame_7.setMinimumSize(QSize(0, 0))
        self.frame_7.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_7.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_9 = QGridLayout(self.frame_7)
        self.gridLayout_9.setObjectName(u"gridLayout_9")
        self.gridLayout_9.setContentsMargins(0, 0, 0, 0)
        self.Lvlx = QLabel(self.frame_7)
        self.Lvlx.setObjectName(u"Lvlx")
        self.Lvlx.setEnabled(False)

        self.gridLayout_9.addWidget(self.Lvlx, 0, 0, 1, 1)

        self.INvly = QLineEdit(self.frame_7)
        self.INvly.setObjectName(u"INvly")
        self.INvly.setEnabled(False)

        self.gridLayout_9.addWidget(self.INvly, 1, 1, 1, 1)

        self.Lvlz = QLabel(self.frame_7)
        self.Lvlz.setObjectName(u"Lvlz")
        self.Lvlz.setEnabled(False)

        self.gridLayout_9.addWidget(self.Lvlz, 2, 0, 1, 1)

        self.Lvly = QLabel(self.frame_7)
        self.Lvly.setObjectName(u"Lvly")
        self.Lvly.setEnabled(False)

        self.gridLayout_9.addWidget(self.Lvly, 1, 0, 1, 1)

        self.INvlz = QLineEdit(self.frame_7)
        self.INvlz.setObjectName(u"INvlz")
        self.INvlz.setEnabled(False)

        self.gridLayout_9.addWidget(self.INvlz, 2, 1, 1, 1)

        self.INvlx = QLineEdit(self.frame_7)
        self.INvlx.setObjectName(u"INvlx")
        self.INvlx.setEnabled(False)

        self.gridLayout_9.addWidget(self.INvlx, 0, 1, 1, 1)

        self.Ulength4 = QLabel(self.frame_7)
        self.Ulength4.setObjectName(u"Ulength4")

        self.gridLayout_9.addWidget(self.Ulength4, 0, 2, 1, 1)

        self.Ulength5 = QLabel(self.frame_7)
        self.Ulength5.setObjectName(u"Ulength5")

        self.gridLayout_9.addWidget(self.Ulength5, 1, 2, 1, 1)

        self.Ulength6 = QLabel(self.frame_7)
        self.Ulength6.setObjectName(u"Ulength6")

        self.gridLayout_9.addWidget(self.Ulength6, 2, 2, 1, 1)


        self.verticalLayout_14.addWidget(self.frame_7)

        self.BAddVTKGeometry = QPushButton(self.page_8)
        self.BAddVTKGeometry.setObjectName(u"BAddVTKGeometry")

        self.verticalLayout_14.addWidget(self.BAddVTKGeometry)

        self.verticalSpacer_13 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_14.addItem(self.verticalSpacer_13)

        self.GeometryOptions.addWidget(self.page_8)

        self.verticalLayout_40.addWidget(self.GeometryOptions)

        self.SAGeometryScrollArea.setWidget(self.scrollAreaWidgetContents)

        self.verticalLayout_15.addWidget(self.SAGeometryScrollArea)

        self.TParts = QTableWidget(self.ImportGeometryTool)
        if (self.TParts.columnCount() < 11):
            self.TParts.setColumnCount(11)
        __qtablewidgetitem13 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(0, __qtablewidgetitem13)
        __qtablewidgetitem14 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(1, __qtablewidgetitem14)
        __qtablewidgetitem15 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(2, __qtablewidgetitem15)
        __qtablewidgetitem16 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(3, __qtablewidgetitem16)
        __qtablewidgetitem17 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(4, __qtablewidgetitem17)
        __qtablewidgetitem18 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(5, __qtablewidgetitem18)
        __qtablewidgetitem19 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(6, __qtablewidgetitem19)
        __qtablewidgetitem20 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(7, __qtablewidgetitem20)
        __qtablewidgetitem21 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(8, __qtablewidgetitem21)
        __qtablewidgetitem22 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(9, __qtablewidgetitem22)
        __qtablewidgetitem23 = QTableWidgetItem()
        self.TParts.setHorizontalHeaderItem(10, __qtablewidgetitem23)
        self.TParts.setObjectName(u"TParts")
        self.TParts.setEnabled(True)
        sizePolicy3.setHeightForWidth(self.TParts.sizePolicy().hasHeightForWidth())
        self.TParts.setSizePolicy(sizePolicy3)
        self.TParts.setMaximumSize(QSize(10000, 180))
        self.TParts.setFont(font1)
        self.TParts.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TParts.setRowCount(0)

        self.verticalLayout_15.addWidget(self.TParts)

        self.BDeleteGeometry = QPushButton(self.ImportGeometryTool)
        self.BDeleteGeometry.setObjectName(u"BDeleteGeometry")

        self.verticalLayout_15.addWidget(self.BDeleteGeometry)

        self.ToolWindow.addWidget(self.ImportGeometryTool)
        self.GenerateMeshTool = QWidget()
        self.GenerateMeshTool.setObjectName(u"GenerateMeshTool")
        self.verticalLayout = QVBoxLayout(self.GenerateMeshTool)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.LDefineGlobalMesh = QLabel(self.GenerateMeshTool)
        self.LDefineGlobalMesh.setObjectName(u"LDefineGlobalMesh")
        sizePolicy8.setHeightForWidth(self.LDefineGlobalMesh.sizePolicy().hasHeightForWidth())
        self.LDefineGlobalMesh.setSizePolicy(sizePolicy8)

        self.verticalLayout.addWidget(self.LDefineGlobalMesh)

        self.frame_5 = QFrame(self.GenerateMeshTool)
        self.frame_5.setObjectName(u"frame_5")
        self.frame_5.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_5.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_5 = QFormLayout(self.frame_5)
        self.formLayout_5.setObjectName(u"formLayout_5")
        self.LElementType = QLabel(self.frame_5)
        self.LElementType.setObjectName(u"LElementType")
        sizePolicy13 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        sizePolicy13.setHorizontalStretch(0)
        sizePolicy13.setVerticalStretch(0)
        sizePolicy13.setHeightForWidth(self.LElementType.sizePolicy().hasHeightForWidth())
        self.LElementType.setSizePolicy(sizePolicy13)

        self.formLayout_5.setWidget(0, QFormLayout.LabelRole, self.LElementType)

        self.INElementType = QComboBox(self.frame_5)
        self.INElementType.addItem("")
        self.INElementType.addItem("")
        self.INElementType.addItem("")
        self.INElementType.setObjectName(u"INElementType")
        sizePolicy11.setHeightForWidth(self.INElementType.sizePolicy().hasHeightForWidth())
        self.INElementType.setSizePolicy(sizePolicy11)

        self.formLayout_5.setWidget(0, QFormLayout.FieldRole, self.INElementType)

        self.LCoordinateSystem = QLabel(self.frame_5)
        self.LCoordinateSystem.setObjectName(u"LCoordinateSystem")
        sizePolicy13.setHeightForWidth(self.LCoordinateSystem.sizePolicy().hasHeightForWidth())
        self.LCoordinateSystem.setSizePolicy(sizePolicy13)

        self.formLayout_5.setWidget(1, QFormLayout.LabelRole, self.LCoordinateSystem)

        self.INCoordinateSystem = QComboBox(self.frame_5)
        self.INCoordinateSystem.addItem("")
        self.INCoordinateSystem.addItem("")
        self.INCoordinateSystem.setObjectName(u"INCoordinateSystem")
        sizePolicy11.setHeightForWidth(self.INCoordinateSystem.sizePolicy().hasHeightForWidth())
        self.INCoordinateSystem.setSizePolicy(sizePolicy11)

        self.formLayout_5.setWidget(1, QFormLayout.FieldRole, self.INCoordinateSystem)

        self.LDimension = QLabel(self.frame_5)
        self.LDimension.setObjectName(u"LDimension")
        sizePolicy13.setHeightForWidth(self.LDimension.sizePolicy().hasHeightForWidth())
        self.LDimension.setSizePolicy(sizePolicy13)

        self.formLayout_5.setWidget(2, QFormLayout.LabelRole, self.LDimension)

        self.INDimension = QComboBox(self.frame_5)
        self.INDimension.addItem("")
        self.INDimension.addItem("")
        self.INDimension.setObjectName(u"INDimension")
        sizePolicy11.setHeightForWidth(self.INDimension.sizePolicy().hasHeightForWidth())
        self.INDimension.setSizePolicy(sizePolicy11)

        self.formLayout_5.setWidget(2, QFormLayout.FieldRole, self.INDimension)


        self.verticalLayout.addWidget(self.frame_5, 0, Qt.AlignmentFlag.AlignTop)

        self.MeshInputs2 = QStackedWidget(self.GenerateMeshTool)
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
        self.Cylindrical2DInputs2.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_13 = QHBoxLayout(self.Cylindrical2DInputs2)
        self.horizontalLayout_13.setObjectName(u"horizontalLayout_13")
        self.horizontalLayout_13.setContentsMargins(-1, 0, -1, 0)
        self.LInnerRadiusC2D = QLabel(self.Cylindrical2DInputs2)
        self.LInnerRadiusC2D.setObjectName(u"LInnerRadiusC2D")

        self.horizontalLayout_13.addWidget(self.LInnerRadiusC2D)

        self.INInnerRadiusC2D = QLineEdit(self.Cylindrical2DInputs2)
        self.INInnerRadiusC2D.setObjectName(u"INInnerRadiusC2D")
        self.INInnerRadiusC2D.setInputMethodHints(Qt.InputMethodHint.ImhNone)

        self.horizontalLayout_13.addWidget(self.INInnerRadiusC2D)


        self.verticalLayout_32.addWidget(self.Cylindrical2DInputs2)

        self.Cylindrical2DInputs = QFrame(self.Cylindrical2D)
        self.Cylindrical2DInputs.setObjectName(u"Cylindrical2DInputs")
        self.Cylindrical2DInputs.setFrameShape(QFrame.Shape.NoFrame)
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
        self.INLengthThetaC2D.setInputMethodHints(Qt.InputMethodHint.ImhNone)

        self.gridLayout_5.addWidget(self.INLengthThetaC2D, 1, 4, 1, 1)

        self.label_38 = QLabel(self.Cylindrical2DInputs)
        self.label_38.setObjectName(u"label_38")

        self.gridLayout_5.addWidget(self.label_38, 2, 1, 1, 1)

        self.INLengthOutRadC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INLengthOutRadC2D.setObjectName(u"INLengthOutRadC2D")
        self.INLengthOutRadC2D.setInputMethodHints(Qt.InputMethodHint.ImhNone)

        self.gridLayout_5.addWidget(self.INLengthOutRadC2D, 1, 2, 1, 1)

        self.label_40 = QLabel(self.Cylindrical2DInputs)
        self.label_40.setObjectName(u"label_40")

        self.gridLayout_5.addWidget(self.label_40, 2, 3, 1, 1)

        self.INOriginYC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INOriginYC2D.setObjectName(u"INOriginYC2D")
        self.INOriginYC2D.setInputMethodHints(Qt.InputMethodHint.ImhNone)

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
        self.INElementsArcC2D.setInputMethodHints(Qt.InputMethodHint.ImhNone)

        self.gridLayout_5.addWidget(self.INElementsArcC2D, 2, 4, 1, 1)

        self.INOriginXC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INOriginXC2D.setObjectName(u"INOriginXC2D")
        self.INOriginXC2D.setInputMethodHints(Qt.InputMethodHint.ImhNone)

        self.gridLayout_5.addWidget(self.INOriginXC2D, 0, 2, 1, 1)

        self.INElementsRadialC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INElementsRadialC2D.setObjectName(u"INElementsRadialC2D")
        self.INElementsRadialC2D.setInputMethodHints(Qt.InputMethodHint.ImhNone)

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
        self.Cylindrical3DInputs2.setFrameShape(QFrame.Shape.NoFrame)
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
        self.Cylindrical3DInputs.setFrameShape(QFrame.Shape.NoFrame)
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

        self.verticalLayout.addWidget(self.MeshInputs2, 0, Qt.AlignmentFlag.AlignTop)

        self.BGenerateGlobalMesh = QPushButton(self.GenerateMeshTool)
        self.BGenerateGlobalMesh.setObjectName(u"BGenerateGlobalMesh")

        self.verticalLayout.addWidget(self.BGenerateGlobalMesh)

        self.verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout.addItem(self.verticalSpacer)

        self.ToolWindow.addWidget(self.GenerateMeshTool)
        self.MaterialsTab = QWidget()
        self.MaterialsTab.setObjectName(u"MaterialsTab")
        self.verticalLayout_5 = QVBoxLayout(self.MaterialsTab)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.verticalLayout_5.setContentsMargins(12, 12, 12, 12)
        self.LDefineMaterials = QLabel(self.MaterialsTab)
        self.LDefineMaterials.setObjectName(u"LDefineMaterials")
        sizePolicy8.setHeightForWidth(self.LDefineMaterials.sizePolicy().hasHeightForWidth())
        self.LDefineMaterials.setSizePolicy(sizePolicy8)

        self.verticalLayout_5.addWidget(self.LDefineMaterials)

        self.DefineAssignMats = QTabWidget(self.MaterialsTab)
        self.DefineAssignMats.setObjectName(u"DefineAssignMats")
        sizePolicy7.setHeightForWidth(self.DefineAssignMats.sizePolicy().hasHeightForWidth())
        self.DefineAssignMats.setSizePolicy(sizePolicy7)
        self.DefineAssignMats.setAutoFillBackground(False)
        self.DefineAssignMats.setIconSize(QSize(20, 20))
        self.DefineAssignMats.setElideMode(Qt.TextElideMode.ElideLeft)
        self.DefineMaterials = QWidget()
        self.DefineMaterials.setObjectName(u"DefineMaterials")
        self.verticalLayout_49 = QVBoxLayout(self.DefineMaterials)
        self.verticalLayout_49.setObjectName(u"verticalLayout_49")
        self.verticalLayout_49.setContentsMargins(0, -1, 0, 0)
        self.INSelectDefineMaterials = QComboBox(self.DefineMaterials)
        self.INSelectDefineMaterials.addItem("")
        self.INSelectDefineMaterials.addItem("")
        self.INSelectDefineMaterials.addItem("")
        self.INSelectDefineMaterials.setObjectName(u"INSelectDefineMaterials")

        self.verticalLayout_49.addWidget(self.INSelectDefineMaterials)

        self.DefineMaterialsOptions = QStackedWidget(self.DefineMaterials)
        self.DefineMaterialsOptions.setObjectName(u"DefineMaterialsOptions")
        self.DefineMaterialsOptions.setMinimumSize(QSize(0, 0))
        self.DefineMaterialsOptions.setMaximumSize(QSize(16777215, 16777215))
        self.DefineMaterialsSGH = QWidget()
        self.DefineMaterialsSGH.setObjectName(u"DefineMaterialsSGH")
        self.verticalLayout_50 = QVBoxLayout(self.DefineMaterialsSGH)
        self.verticalLayout_50.setObjectName(u"verticalLayout_50")
        self.verticalLayout_50.setContentsMargins(0, 0, 0, 0)
        self.frame_32 = QFrame(self.DefineMaterialsSGH)
        self.frame_32.setObjectName(u"frame_32")
        self.frame_32.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_32.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_25 = QFormLayout(self.frame_32)
        self.formLayout_25.setObjectName(u"formLayout_25")
        self.LMaterialNameSGH = QLabel(self.frame_32)
        self.LMaterialNameSGH.setObjectName(u"LMaterialNameSGH")
        sizePolicy12.setHeightForWidth(self.LMaterialNameSGH.sizePolicy().hasHeightForWidth())
        self.LMaterialNameSGH.setSizePolicy(sizePolicy12)

        self.formLayout_25.setWidget(0, QFormLayout.LabelRole, self.LMaterialNameSGH)

        self.INMaterialNameSGH = QLineEdit(self.frame_32)
        self.INMaterialNameSGH.setObjectName(u"INMaterialNameSGH")
        sizePolicy9.setHeightForWidth(self.INMaterialNameSGH.sizePolicy().hasHeightForWidth())
        self.INMaterialNameSGH.setSizePolicy(sizePolicy9)
        self.INMaterialNameSGH.setMinimumSize(QSize(150, 0))

        self.formLayout_25.setWidget(0, QFormLayout.FieldRole, self.INMaterialNameSGH)

        self.LEOS = QLabel(self.frame_32)
        self.LEOS.setObjectName(u"LEOS")

        self.formLayout_25.setWidget(1, QFormLayout.LabelRole, self.LEOS)

        self.INEOS = QComboBox(self.frame_32)
        self.INEOS.addItem("")
        self.INEOS.setObjectName(u"INEOS")
        self.INEOS.setMinimumSize(QSize(150, 0))

        self.formLayout_25.setWidget(1, QFormLayout.FieldRole, self.INEOS)

        self.LArtificialViscosity = QLabel(self.frame_32)
        self.LArtificialViscosity.setObjectName(u"LArtificialViscosity")
        self.LArtificialViscosity.setTextFormat(Qt.TextFormat.PlainText)

        self.formLayout_25.setWidget(2, QFormLayout.LabelRole, self.LArtificialViscosity)

        self.INArtificialViscosity = QComboBox(self.frame_32)
        self.INArtificialViscosity.addItem("")
        self.INArtificialViscosity.addItem("")
        self.INArtificialViscosity.setObjectName(u"INArtificialViscosity")
        self.INArtificialViscosity.setMinimumSize(QSize(150, 0))

        self.formLayout_25.setWidget(2, QFormLayout.FieldRole, self.INArtificialViscosity)


        self.verticalLayout_50.addWidget(self.frame_32, 0, Qt.AlignmentFlag.AlignTop)

        self.BAddMaterialSGH = QPushButton(self.DefineMaterialsSGH)
        self.BAddMaterialSGH.setObjectName(u"BAddMaterialSGH")

        self.verticalLayout_50.addWidget(self.BAddMaterialSGH)

        self.TMaterialsSGH = QTableWidget(self.DefineMaterialsSGH)
        if (self.TMaterialsSGH.columnCount() < 9):
            self.TMaterialsSGH.setColumnCount(9)
        __qtablewidgetitem24 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(0, __qtablewidgetitem24)
        __qtablewidgetitem25 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(1, __qtablewidgetitem25)
        __qtablewidgetitem26 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(2, __qtablewidgetitem26)
        __qtablewidgetitem27 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(3, __qtablewidgetitem27)
        __qtablewidgetitem28 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(4, __qtablewidgetitem28)
        __qtablewidgetitem29 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(5, __qtablewidgetitem29)
        __qtablewidgetitem30 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(6, __qtablewidgetitem30)
        __qtablewidgetitem31 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(7, __qtablewidgetitem31)
        __qtablewidgetitem32 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(8, __qtablewidgetitem32)
        self.TMaterialsSGH.setObjectName(u"TMaterialsSGH")
        self.TMaterialsSGH.setEnabled(True)
        sizePolicy7.setHeightForWidth(self.TMaterialsSGH.sizePolicy().hasHeightForWidth())
        self.TMaterialsSGH.setSizePolicy(sizePolicy7)
        self.TMaterialsSGH.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TMaterialsSGH.setRowCount(0)

        self.verticalLayout_50.addWidget(self.TMaterialsSGH)

        self.BDeleteMaterialSGH = QPushButton(self.DefineMaterialsSGH)
        self.BDeleteMaterialSGH.setObjectName(u"BDeleteMaterialSGH")

        self.verticalLayout_50.addWidget(self.BDeleteMaterialSGH)

        self.label_24 = QLabel(self.DefineMaterialsSGH)
        self.label_24.setObjectName(u"label_24")
        font12 = QFont()
        font12.setPointSize(11)
        self.label_24.setFont(font12)

        self.verticalLayout_50.addWidget(self.label_24)

        self.frame_33 = QFrame(self.DefineMaterialsSGH)
        self.frame_33.setObjectName(u"frame_33")
        self.frame_33.setFrameShape(QFrame.Shape.NoFrame)
        self.gridLayout_11 = QGridLayout(self.frame_33)
        self.gridLayout_11.setObjectName(u"gridLayout_11")
        self.gridLayout_11.setVerticalSpacing(6)
        self.gridLayout_11.setContentsMargins(0, 0, 0, 0)
        self.INq1 = QLineEdit(self.frame_33)
        self.INq1.setObjectName(u"INq1")
        self.INq1.setEnabled(True)
        sizePolicy1.setHeightForWidth(self.INq1.sizePolicy().hasHeightForWidth())
        self.INq1.setSizePolicy(sizePolicy1)

        self.gridLayout_11.addWidget(self.INq1, 0, 1, 1, 1)

        self.INSpecificHeat = QLineEdit(self.frame_33)
        self.INSpecificHeat.setObjectName(u"INSpecificHeat")
        self.INSpecificHeat.setEnabled(True)
        sizePolicy1.setHeightForWidth(self.INSpecificHeat.sizePolicy().hasHeightForWidth())
        self.INSpecificHeat.setSizePolicy(sizePolicy1)

        self.gridLayout_11.addWidget(self.INSpecificHeat, 6, 1, 1, 1)

        self.Lq1ex = QLabel(self.frame_33)
        self.Lq1ex.setObjectName(u"Lq1ex")
        self.Lq1ex.setEnabled(True)
        sizePolicy8.setHeightForWidth(self.Lq1ex.sizePolicy().hasHeightForWidth())
        self.Lq1ex.setSizePolicy(sizePolicy8)
        self.Lq1ex.setWordWrap(True)

        self.gridLayout_11.addWidget(self.Lq1ex, 2, 0, 1, 1)

        self.LGamma = QLabel(self.frame_33)
        self.LGamma.setObjectName(u"LGamma")
        self.LGamma.setEnabled(True)

        self.gridLayout_11.addWidget(self.LGamma, 4, 0, 1, 1)

        self.INq1ex = QLineEdit(self.frame_33)
        self.INq1ex.setObjectName(u"INq1ex")
        self.INq1ex.setEnabled(True)
        sizePolicy1.setHeightForWidth(self.INq1ex.sizePolicy().hasHeightForWidth())
        self.INq1ex.setSizePolicy(sizePolicy1)

        self.gridLayout_11.addWidget(self.INq1ex, 2, 1, 1, 1)

        self.INq2 = QLineEdit(self.frame_33)
        self.INq2.setObjectName(u"INq2")
        self.INq2.setEnabled(True)
        sizePolicy1.setHeightForWidth(self.INq2.sizePolicy().hasHeightForWidth())
        self.INq2.setSizePolicy(sizePolicy1)

        self.gridLayout_11.addWidget(self.INq2, 1, 1, 1, 1)

        self.INq2ex = QLineEdit(self.frame_33)
        self.INq2ex.setObjectName(u"INq2ex")
        self.INq2ex.setEnabled(True)
        sizePolicy1.setHeightForWidth(self.INq2ex.sizePolicy().hasHeightForWidth())
        self.INq2ex.setSizePolicy(sizePolicy1)

        self.gridLayout_11.addWidget(self.INq2ex, 3, 1, 1, 1)

        self.LMinSound = QLabel(self.frame_33)
        self.LMinSound.setObjectName(u"LMinSound")
        self.LMinSound.setEnabled(True)

        self.gridLayout_11.addWidget(self.LMinSound, 5, 0, 1, 1)

        self.Lq2ex = QLabel(self.frame_33)
        self.Lq2ex.setObjectName(u"Lq2ex")
        self.Lq2ex.setEnabled(True)
        self.Lq2ex.setWordWrap(True)

        self.gridLayout_11.addWidget(self.Lq2ex, 3, 0, 1, 1)

        self.Lq1 = QLabel(self.frame_33)
        self.Lq1.setObjectName(u"Lq1")
        self.Lq1.setEnabled(True)
        self.Lq1.setWordWrap(True)

        self.gridLayout_11.addWidget(self.Lq1, 0, 0, 1, 1)

        self.LSpecificHeat = QLabel(self.frame_33)
        self.LSpecificHeat.setObjectName(u"LSpecificHeat")
        self.LSpecificHeat.setEnabled(True)

        self.gridLayout_11.addWidget(self.LSpecificHeat, 6, 0, 1, 1)

        self.INGamma = QLineEdit(self.frame_33)
        self.INGamma.setObjectName(u"INGamma")
        self.INGamma.setEnabled(True)
        sizePolicy1.setHeightForWidth(self.INGamma.sizePolicy().hasHeightForWidth())
        self.INGamma.setSizePolicy(sizePolicy1)

        self.gridLayout_11.addWidget(self.INGamma, 4, 1, 1, 1)

        self.Lq2 = QLabel(self.frame_33)
        self.Lq2.setObjectName(u"Lq2")
        self.Lq2.setEnabled(True)
        self.Lq2.setWordWrap(True)

        self.gridLayout_11.addWidget(self.Lq2, 1, 0, 1, 1)

        self.INMinSound = QLineEdit(self.frame_33)
        self.INMinSound.setObjectName(u"INMinSound")
        self.INMinSound.setEnabled(True)
        sizePolicy1.setHeightForWidth(self.INMinSound.sizePolicy().hasHeightForWidth())
        self.INMinSound.setSizePolicy(sizePolicy1)

        self.gridLayout_11.addWidget(self.INMinSound, 5, 1, 1, 1)


        self.verticalLayout_50.addWidget(self.frame_33)

        self.DefineMaterialsOptions.addWidget(self.DefineMaterialsSGH)
        self.DefineMaterialsEVPFFT = QWidget()
        self.DefineMaterialsEVPFFT.setObjectName(u"DefineMaterialsEVPFFT")
        self.verticalLayout_51 = QVBoxLayout(self.DefineMaterialsEVPFFT)
        self.verticalLayout_51.setObjectName(u"verticalLayout_51")
        self.verticalLayout_51.setContentsMargins(0, -1, 0, -1)
        self.MaterialInputs = QFrame(self.DefineMaterialsEVPFFT)
        self.MaterialInputs.setObjectName(u"MaterialInputs")
        self.MaterialInputs.setFrameShape(QFrame.Shape.NoFrame)
        self.gridLayout_7 = QGridLayout(self.MaterialInputs)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.gridLayout_7.setContentsMargins(0, 0, 0, 0)
        self.LMaterialName_2 = QLabel(self.MaterialInputs)
        self.LMaterialName_2.setObjectName(u"LMaterialName_2")

        self.gridLayout_7.addWidget(self.LMaterialName_2, 0, 0, 1, 1)

        self.INMaterialName = QLineEdit(self.MaterialInputs)
        self.INMaterialName.setObjectName(u"INMaterialName")
        self.INMaterialName.setMinimumSize(QSize(93, 0))

        self.gridLayout_7.addWidget(self.INMaterialName, 0, 1, 1, 1)

        self.frame_35 = QFrame(self.MaterialInputs)
        self.frame_35.setObjectName(u"frame_35")
        self.frame_35.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_15 = QHBoxLayout(self.frame_35)
        self.horizontalLayout_15.setObjectName(u"horizontalLayout_15")
        self.horizontalLayout_15.setContentsMargins(0, 0, 0, 0)

        self.gridLayout_7.addWidget(self.frame_35, 0, 2, 1, 1)

        self.frame_34 = QFrame(self.MaterialInputs)
        self.frame_34.setObjectName(u"frame_34")
        self.frame_34.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_11 = QHBoxLayout(self.frame_34)
        self.horizontalLayout_11.setObjectName(u"horizontalLayout_11")
        self.horizontalLayout_11.setContentsMargins(0, 0, 0, 0)

        self.gridLayout_7.addWidget(self.frame_34, 0, 3, 1, 1)


        self.verticalLayout_51.addWidget(self.MaterialInputs)

        self.MaterialMenu = QTabWidget(self.DefineMaterialsEVPFFT)
        self.MaterialMenu.setObjectName(u"MaterialMenu")
        self.MaterialMenu.setEnabled(True)
        self.Elastic = QWidget()
        self.Elastic.setObjectName(u"Elastic")
        self.verticalLayout_23 = QVBoxLayout(self.Elastic)
        self.verticalLayout_23.setObjectName(u"verticalLayout_23")
        self.verticalLayout_23.setContentsMargins(6, 12, 6, 0)
        self.frame_24 = QFrame(self.Elastic)
        self.frame_24.setObjectName(u"frame_24")
        self.frame_24.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_24.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_19 = QGridLayout(self.frame_24)
        self.gridLayout_19.setObjectName(u"gridLayout_19")
        self.gridLayout_19.setContentsMargins(0, 0, 0, 0)
        self.LType = QLabel(self.frame_24)
        self.LType.setObjectName(u"LType")

        self.gridLayout_19.addWidget(self.LType, 0, 0, 1, 1)

        self.frame_36 = QFrame(self.frame_24)
        self.frame_36.setObjectName(u"frame_36")
        self.frame_36.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_18 = QHBoxLayout(self.frame_36)
        self.horizontalLayout_18.setObjectName(u"horizontalLayout_18")
        self.horizontalLayout_18.setContentsMargins(0, 0, 0, 0)
        self.INSolidGas = QComboBox(self.frame_36)
        self.INSolidGas.addItem("")
        self.INSolidGas.addItem("")
        self.INSolidGas.setObjectName(u"INSolidGas")
        self.INSolidGas.setMinimumSize(QSize(82, 0))

        self.horizontalLayout_18.addWidget(self.INSolidGas)

        self.INMaterialType = QComboBox(self.frame_36)
        self.INMaterialType.addItem("")
        self.INMaterialType.addItem("")
        self.INMaterialType.addItem("")
        self.INMaterialType.addItem("")
        self.INMaterialType.setObjectName(u"INMaterialType")

        self.horizontalLayout_18.addWidget(self.INMaterialType)


        self.gridLayout_19.addWidget(self.frame_36, 0, 1, 1, 1)


        self.verticalLayout_23.addWidget(self.frame_24)

        self.MaterialTypeTool = QStackedWidget(self.Elastic)
        self.MaterialTypeTool.setObjectName(u"MaterialTypeTool")
        self.Isotropic = QWidget()
        self.Isotropic.setObjectName(u"Isotropic")
        self.Isotropic.setEnabled(True)
        self.gridLayout_13 = QGridLayout(self.Isotropic)
        self.gridLayout_13.setObjectName(u"gridLayout_13")
        self.gridLayout_13.setContentsMargins(-1, -1, -1, 12)
        self.LPoissonsRatio = QLabel(self.Isotropic)
        self.LPoissonsRatio.setObjectName(u"LPoissonsRatio")

        self.gridLayout_13.addWidget(self.LPoissonsRatio, 1, 0, 1, 1)

        self.INPoissonsRatio = QLineEdit(self.Isotropic)
        self.INPoissonsRatio.setObjectName(u"INPoissonsRatio")

        self.gridLayout_13.addWidget(self.INPoissonsRatio, 1, 1, 1, 1)

        self.INYoungsModulus = QLineEdit(self.Isotropic)
        self.INYoungsModulus.setObjectName(u"INYoungsModulus")

        self.gridLayout_13.addWidget(self.INYoungsModulus, 0, 1, 1, 1)

        self.verticalSpacer_14 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.gridLayout_13.addItem(self.verticalSpacer_14, 2, 1, 1, 1)

        self.UPressure1 = QLabel(self.Isotropic)
        self.UPressure1.setObjectName(u"UPressure1")

        self.gridLayout_13.addWidget(self.UPressure1, 0, 2, 1, 1)

        self.LYoungsModulus = QLabel(self.Isotropic)
        self.LYoungsModulus.setObjectName(u"LYoungsModulus")

        self.gridLayout_13.addWidget(self.LYoungsModulus, 0, 0, 1, 1)

        self.MaterialTypeTool.addWidget(self.Isotropic)
        self.TransverselyIsotropic = QWidget()
        self.TransverselyIsotropic.setObjectName(u"TransverselyIsotropic")
        self.verticalLayout_52 = QVBoxLayout(self.TransverselyIsotropic)
        self.verticalLayout_52.setSpacing(0)
        self.verticalLayout_52.setObjectName(u"verticalLayout_52")
        self.verticalLayout_52.setContentsMargins(0, 0, 0, 0)
        self.frame_37 = QFrame(self.TransverselyIsotropic)
        self.frame_37.setObjectName(u"frame_37")
        self.frame_37.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_37.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_27 = QFormLayout(self.frame_37)
        self.formLayout_27.setObjectName(u"formLayout_27")
        self.formLayout_27.setContentsMargins(0, 0, -1, 0)
        self.LIsotropicPlane = QLabel(self.frame_37)
        self.LIsotropicPlane.setObjectName(u"LIsotropicPlane")

        self.formLayout_27.setWidget(0, QFormLayout.LabelRole, self.LIsotropicPlane)

        self.INIsotropicPlane = QComboBox(self.frame_37)
        self.INIsotropicPlane.addItem("")
        self.INIsotropicPlane.addItem("")
        self.INIsotropicPlane.addItem("")
        self.INIsotropicPlane.setObjectName(u"INIsotropicPlane")

        self.formLayout_27.setWidget(0, QFormLayout.FieldRole, self.INIsotropicPlane)


        self.verticalLayout_52.addWidget(self.frame_37, 0, Qt.AlignmentFlag.AlignTop)

        self.TransverslyIsotropicMat = QFrame(self.TransverselyIsotropic)
        self.TransverslyIsotropicMat.setObjectName(u"TransverslyIsotropicMat")
        self.TransverslyIsotropicMat.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_19 = QHBoxLayout(self.TransverslyIsotropicMat)
        self.horizontalLayout_19.setSpacing(0)
        self.horizontalLayout_19.setObjectName(u"horizontalLayout_19")
        self.horizontalLayout_19.setContentsMargins(0, 0, 0, 0)
        self.TransverseInPlane = QFrame(self.TransverslyIsotropicMat)
        self.TransverseInPlane.setObjectName(u"TransverseInPlane")
        self.TransverseInPlane.setFrameShape(QFrame.Shape.NoFrame)
        self.verticalLayout_53 = QVBoxLayout(self.TransverseInPlane)
        self.verticalLayout_53.setSpacing(0)
        self.verticalLayout_53.setObjectName(u"verticalLayout_53")
        self.verticalLayout_53.setContentsMargins(0, 12, 0, 0)
        self.LInPlane = QLabel(self.TransverseInPlane)
        self.LInPlane.setObjectName(u"LInPlane")

        self.verticalLayout_53.addWidget(self.LInPlane)

        self.TransverseInPlaneMat = QFrame(self.TransverseInPlane)
        self.TransverseInPlaneMat.setObjectName(u"TransverseInPlaneMat")
        self.TransverseInPlaneMat.setFrameShape(QFrame.Shape.NoFrame)
        self.gridLayout_14 = QGridLayout(self.TransverseInPlaneMat)
        self.gridLayout_14.setObjectName(u"gridLayout_14")
        self.INNUip = QLineEdit(self.TransverseInPlaneMat)
        self.INNUip.setObjectName(u"INNUip")

        self.gridLayout_14.addWidget(self.INNUip, 1, 1, 1, 1)

        self.LNUip = QLabel(self.TransverseInPlaneMat)
        self.LNUip.setObjectName(u"LNUip")
        self.LNUip.setMinimumSize(QSize(31, 0))

        self.gridLayout_14.addWidget(self.LNUip, 1, 0, 1, 1)

        self.LEip = QLabel(self.TransverseInPlaneMat)
        self.LEip.setObjectName(u"LEip")
        self.LEip.setMinimumSize(QSize(31, 0))

        self.gridLayout_14.addWidget(self.LEip, 0, 0, 1, 1)

        self.INEip = QLineEdit(self.TransverseInPlaneMat)
        self.INEip.setObjectName(u"INEip")

        self.gridLayout_14.addWidget(self.INEip, 0, 1, 1, 1)

        self.UPressure2 = QLabel(self.TransverseInPlaneMat)
        self.UPressure2.setObjectName(u"UPressure2")

        self.gridLayout_14.addWidget(self.UPressure2, 0, 2, 1, 1)


        self.verticalLayout_53.addWidget(self.TransverseInPlaneMat)


        self.horizontalLayout_19.addWidget(self.TransverseInPlane)

        self.TransverseOutOfPlane = QFrame(self.TransverslyIsotropicMat)
        self.TransverseOutOfPlane.setObjectName(u"TransverseOutOfPlane")
        self.TransverseOutOfPlane.setFrameShape(QFrame.Shape.NoFrame)
        self.verticalLayout_54 = QVBoxLayout(self.TransverseOutOfPlane)
        self.verticalLayout_54.setSpacing(0)
        self.verticalLayout_54.setObjectName(u"verticalLayout_54")
        self.verticalLayout_54.setContentsMargins(0, 12, 0, 0)
        self.LOutOfPlane = QLabel(self.TransverseOutOfPlane)
        self.LOutOfPlane.setObjectName(u"LOutOfPlane")

        self.verticalLayout_54.addWidget(self.LOutOfPlane)

        self.TransverseOutOfPlaneMat = QFrame(self.TransverseOutOfPlane)
        self.TransverseOutOfPlaneMat.setObjectName(u"TransverseOutOfPlaneMat")
        self.TransverseOutOfPlaneMat.setFrameShape(QFrame.Shape.NoFrame)
        self.gridLayout_15 = QGridLayout(self.TransverseOutOfPlaneMat)
        self.gridLayout_15.setObjectName(u"gridLayout_15")
        self.INNUop = QLineEdit(self.TransverseOutOfPlaneMat)
        self.INNUop.setObjectName(u"INNUop")

        self.gridLayout_15.addWidget(self.INNUop, 1, 1, 1, 1)

        self.INGop = QLineEdit(self.TransverseOutOfPlaneMat)
        self.INGop.setObjectName(u"INGop")

        self.gridLayout_15.addWidget(self.INGop, 2, 1, 1, 1)

        self.LEop = QLabel(self.TransverseOutOfPlaneMat)
        self.LEop.setObjectName(u"LEop")
        self.LEop.setMinimumSize(QSize(31, 0))
        self.LEop.setMaximumSize(QSize(31, 16777215))

        self.gridLayout_15.addWidget(self.LEop, 0, 0, 1, 1)

        self.LNUop = QLabel(self.TransverseOutOfPlaneMat)
        self.LNUop.setObjectName(u"LNUop")
        self.LNUop.setMinimumSize(QSize(31, 0))
        self.LNUop.setMaximumSize(QSize(31, 16777215))

        self.gridLayout_15.addWidget(self.LNUop, 1, 0, 1, 1)

        self.LGop = QLabel(self.TransverseOutOfPlaneMat)
        self.LGop.setObjectName(u"LGop")
        self.LGop.setMinimumSize(QSize(31, 0))
        self.LGop.setMaximumSize(QSize(31, 16777215))

        self.gridLayout_15.addWidget(self.LGop, 2, 0, 1, 1)

        self.INEop = QLineEdit(self.TransverseOutOfPlaneMat)
        self.INEop.setObjectName(u"INEop")
        self.INEop.setMinimumSize(QSize(0, 0))

        self.gridLayout_15.addWidget(self.INEop, 0, 1, 1, 1)

        self.UPressure3 = QLabel(self.TransverseOutOfPlaneMat)
        self.UPressure3.setObjectName(u"UPressure3")

        self.gridLayout_15.addWidget(self.UPressure3, 0, 2, 1, 1)

        self.UPressure4 = QLabel(self.TransverseOutOfPlaneMat)
        self.UPressure4.setObjectName(u"UPressure4")

        self.gridLayout_15.addWidget(self.UPressure4, 2, 2, 1, 1)


        self.verticalLayout_54.addWidget(self.TransverseOutOfPlaneMat)


        self.horizontalLayout_19.addWidget(self.TransverseOutOfPlane)


        self.verticalLayout_52.addWidget(self.TransverslyIsotropicMat)

        self.verticalSpacer_22 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_52.addItem(self.verticalSpacer_22)

        self.MaterialTypeTool.addWidget(self.TransverselyIsotropic)
        self.Anisotropic = QWidget()
        self.Anisotropic.setObjectName(u"Anisotropic")
        self.verticalLayout_55 = QVBoxLayout(self.Anisotropic)
        self.verticalLayout_55.setObjectName(u"verticalLayout_55")
        self.verticalLayout_55.setContentsMargins(0, 0, 0, 0)
        self.label_2 = QLabel(self.Anisotropic)
        self.label_2.setObjectName(u"label_2")

        self.verticalLayout_55.addWidget(self.label_2)

        self.TAnisotropic = QTableWidget(self.Anisotropic)
        if (self.TAnisotropic.columnCount() < 6):
            self.TAnisotropic.setColumnCount(6)
        if (self.TAnisotropic.rowCount() < 6):
            self.TAnisotropic.setRowCount(6)
        __qtablewidgetitem33 = QTableWidgetItem()
        self.TAnisotropic.setItem(0, 0, __qtablewidgetitem33)
        brush = QBrush(QColor(235, 235, 235, 255))
        brush.setStyle(Qt.SolidPattern)
        __qtablewidgetitem34 = QTableWidgetItem()
        __qtablewidgetitem34.setBackground(brush);
        __qtablewidgetitem34.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(1, 0, __qtablewidgetitem34)
        __qtablewidgetitem35 = QTableWidgetItem()
        __qtablewidgetitem35.setBackground(brush);
        __qtablewidgetitem35.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(2, 0, __qtablewidgetitem35)
        __qtablewidgetitem36 = QTableWidgetItem()
        __qtablewidgetitem36.setBackground(brush);
        __qtablewidgetitem36.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(2, 1, __qtablewidgetitem36)
        __qtablewidgetitem37 = QTableWidgetItem()
        __qtablewidgetitem37.setBackground(brush);
        __qtablewidgetitem37.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 0, __qtablewidgetitem37)
        __qtablewidgetitem38 = QTableWidgetItem()
        __qtablewidgetitem38.setBackground(brush);
        __qtablewidgetitem38.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 1, __qtablewidgetitem38)
        __qtablewidgetitem39 = QTableWidgetItem()
        __qtablewidgetitem39.setBackground(brush);
        __qtablewidgetitem39.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 2, __qtablewidgetitem39)
        __qtablewidgetitem40 = QTableWidgetItem()
        __qtablewidgetitem40.setBackground(brush);
        __qtablewidgetitem40.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 0, __qtablewidgetitem40)
        __qtablewidgetitem41 = QTableWidgetItem()
        __qtablewidgetitem41.setBackground(brush);
        __qtablewidgetitem41.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 1, __qtablewidgetitem41)
        __qtablewidgetitem42 = QTableWidgetItem()
        __qtablewidgetitem42.setBackground(brush);
        __qtablewidgetitem42.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 2, __qtablewidgetitem42)
        __qtablewidgetitem43 = QTableWidgetItem()
        __qtablewidgetitem43.setBackground(brush);
        __qtablewidgetitem43.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 3, __qtablewidgetitem43)
        __qtablewidgetitem44 = QTableWidgetItem()
        __qtablewidgetitem44.setBackground(brush);
        __qtablewidgetitem44.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 0, __qtablewidgetitem44)
        __qtablewidgetitem45 = QTableWidgetItem()
        __qtablewidgetitem45.setBackground(brush);
        __qtablewidgetitem45.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 1, __qtablewidgetitem45)
        __qtablewidgetitem46 = QTableWidgetItem()
        __qtablewidgetitem46.setBackground(brush);
        __qtablewidgetitem46.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 2, __qtablewidgetitem46)
        __qtablewidgetitem47 = QTableWidgetItem()
        __qtablewidgetitem47.setBackground(brush);
        __qtablewidgetitem47.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 3, __qtablewidgetitem47)
        __qtablewidgetitem48 = QTableWidgetItem()
        __qtablewidgetitem48.setBackground(brush);
        __qtablewidgetitem48.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 4, __qtablewidgetitem48)
        self.TAnisotropic.setObjectName(u"TAnisotropic")
        self.TAnisotropic.setRowCount(6)
        self.TAnisotropic.setColumnCount(6)
        self.TAnisotropic.horizontalHeader().setMinimumSectionSize(21)
        self.TAnisotropic.horizontalHeader().setDefaultSectionSize(50)

        self.verticalLayout_55.addWidget(self.TAnisotropic)

        self.MaterialTypeTool.addWidget(self.Anisotropic)
        self.Orthotropic = QWidget()
        self.Orthotropic.setObjectName(u"Orthotropic")
        self.gridLayout_8 = QGridLayout(self.Orthotropic)
        self.gridLayout_8.setObjectName(u"gridLayout_8")
        self.gridLayout_8.setHorizontalSpacing(6)
        self.gridLayout_8.setVerticalSpacing(0)
        self.gridLayout_8.setContentsMargins(0, 0, 0, 0)
        self.INGxy = QLineEdit(self.Orthotropic)
        self.INGxy.setObjectName(u"INGxy")

        self.gridLayout_8.addWidget(self.INGxy, 3, 2, 1, 1)

        self.LGxy = QLabel(self.Orthotropic)
        self.LGxy.setObjectName(u"LGxy")

        self.gridLayout_8.addWidget(self.LGxy, 3, 1, 1, 1)

        self.UPressure5 = QLabel(self.Orthotropic)
        self.UPressure5.setObjectName(u"UPressure5")

        self.gridLayout_8.addWidget(self.UPressure5, 0, 7, 1, 1)

        self.LEx = QLabel(self.Orthotropic)
        self.LEx.setObjectName(u"LEx")

        self.gridLayout_8.addWidget(self.LEx, 0, 1, 1, 1)

        self.INNUxz = QLineEdit(self.Orthotropic)
        self.INNUxz.setObjectName(u"INNUxz")

        self.gridLayout_8.addWidget(self.INNUxz, 2, 4, 1, 1)

        self.LNUxy = QLabel(self.Orthotropic)
        self.LNUxy.setObjectName(u"LNUxy")

        self.gridLayout_8.addWidget(self.LNUxy, 2, 1, 1, 1)

        self.LEy = QLabel(self.Orthotropic)
        self.LEy.setObjectName(u"LEy")

        self.gridLayout_8.addWidget(self.LEy, 0, 3, 1, 1)

        self.INGyz = QLineEdit(self.Orthotropic)
        self.INGyz.setObjectName(u"INGyz")

        self.gridLayout_8.addWidget(self.INGyz, 3, 6, 1, 1)

        self.LNUxz = QLabel(self.Orthotropic)
        self.LNUxz.setObjectName(u"LNUxz")

        self.gridLayout_8.addWidget(self.LNUxz, 2, 3, 1, 1)

        self.INEy = QLineEdit(self.Orthotropic)
        self.INEy.setObjectName(u"INEy")

        self.gridLayout_8.addWidget(self.INEy, 0, 4, 1, 1)

        self.LGyz = QLabel(self.Orthotropic)
        self.LGyz.setObjectName(u"LGyz")

        self.gridLayout_8.addWidget(self.LGyz, 3, 5, 1, 1)

        self.INEz = QLineEdit(self.Orthotropic)
        self.INEz.setObjectName(u"INEz")

        self.gridLayout_8.addWidget(self.INEz, 0, 6, 1, 1)

        self.INNUyz = QLineEdit(self.Orthotropic)
        self.INNUyz.setObjectName(u"INNUyz")

        self.gridLayout_8.addWidget(self.INNUyz, 2, 6, 1, 1)

        self.INGxz = QLineEdit(self.Orthotropic)
        self.INGxz.setObjectName(u"INGxz")

        self.gridLayout_8.addWidget(self.INGxz, 3, 4, 1, 1)

        self.LNUyz = QLabel(self.Orthotropic)
        self.LNUyz.setObjectName(u"LNUyz")

        self.gridLayout_8.addWidget(self.LNUyz, 2, 5, 1, 1)

        self.LEz = QLabel(self.Orthotropic)
        self.LEz.setObjectName(u"LEz")

        self.gridLayout_8.addWidget(self.LEz, 0, 5, 1, 1)

        self.INNUxy = QLineEdit(self.Orthotropic)
        self.INNUxy.setObjectName(u"INNUxy")

        self.gridLayout_8.addWidget(self.INNUxy, 2, 2, 1, 1)

        self.INEx = QLineEdit(self.Orthotropic)
        self.INEx.setObjectName(u"INEx")

        self.gridLayout_8.addWidget(self.INEx, 0, 2, 1, 1)

        self.LGxz = QLabel(self.Orthotropic)
        self.LGxz.setObjectName(u"LGxz")

        self.gridLayout_8.addWidget(self.LGxz, 3, 3, 1, 1)

        self.UPressure6 = QLabel(self.Orthotropic)
        self.UPressure6.setObjectName(u"UPressure6")

        self.gridLayout_8.addWidget(self.UPressure6, 3, 7, 1, 1)

        self.MaterialTypeTool.addWidget(self.Orthotropic)
        self.page_4 = QWidget()
        self.page_4.setObjectName(u"page_4")
        self.MaterialTypeTool.addWidget(self.page_4)

        self.verticalLayout_23.addWidget(self.MaterialTypeTool)

        icon11 = QIcon()
        icon11.addFile(u":/Blue Icons/Blue Icons/Elastic.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.MaterialMenu.addTab(self.Elastic, icon11, "")

        self.verticalLayout_51.addWidget(self.MaterialMenu, 0, Qt.AlignmentFlag.AlignTop)

        self.BAddMaterial = QPushButton(self.DefineMaterialsEVPFFT)
        self.BAddMaterial.setObjectName(u"BAddMaterial")

        self.verticalLayout_51.addWidget(self.BAddMaterial)

        self.TMaterials = QTableWidget(self.DefineMaterialsEVPFFT)
        if (self.TMaterials.columnCount() < 23):
            self.TMaterials.setColumnCount(23)
        __qtablewidgetitem49 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(0, __qtablewidgetitem49)
        __qtablewidgetitem50 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(1, __qtablewidgetitem50)
        __qtablewidgetitem51 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(2, __qtablewidgetitem51)
        __qtablewidgetitem52 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(3, __qtablewidgetitem52)
        __qtablewidgetitem53 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(4, __qtablewidgetitem53)
        __qtablewidgetitem54 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(5, __qtablewidgetitem54)
        __qtablewidgetitem55 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(6, __qtablewidgetitem55)
        __qtablewidgetitem56 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(7, __qtablewidgetitem56)
        __qtablewidgetitem57 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(8, __qtablewidgetitem57)
        __qtablewidgetitem58 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(9, __qtablewidgetitem58)
        __qtablewidgetitem59 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(10, __qtablewidgetitem59)
        __qtablewidgetitem60 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(11, __qtablewidgetitem60)
        __qtablewidgetitem61 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(12, __qtablewidgetitem61)
        __qtablewidgetitem62 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(13, __qtablewidgetitem62)
        __qtablewidgetitem63 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(14, __qtablewidgetitem63)
        __qtablewidgetitem64 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(15, __qtablewidgetitem64)
        __qtablewidgetitem65 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(16, __qtablewidgetitem65)
        __qtablewidgetitem66 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(17, __qtablewidgetitem66)
        __qtablewidgetitem67 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(18, __qtablewidgetitem67)
        __qtablewidgetitem68 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(19, __qtablewidgetitem68)
        __qtablewidgetitem69 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(20, __qtablewidgetitem69)
        __qtablewidgetitem70 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(21, __qtablewidgetitem70)
        __qtablewidgetitem71 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(22, __qtablewidgetitem71)
        self.TMaterials.setObjectName(u"TMaterials")
        self.TMaterials.setEnabled(True)
        self.TMaterials.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustToContents)
        self.TMaterials.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TMaterials.setRowCount(0)

        self.verticalLayout_51.addWidget(self.TMaterials)

        self.BDeleteMaterial = QPushButton(self.DefineMaterialsEVPFFT)
        self.BDeleteMaterial.setObjectName(u"BDeleteMaterial")

        self.verticalLayout_51.addWidget(self.BDeleteMaterial)

        self.BRegenElasticConstants = QPushButton(self.DefineMaterialsEVPFFT)
        self.BRegenElasticConstants.setObjectName(u"BRegenElasticConstants")

        self.verticalLayout_51.addWidget(self.BRegenElasticConstants)

        self.DefineMaterialsOptions.addWidget(self.DefineMaterialsEVPFFT)
        self.DefineMaterialsBulkForm = QWidget()
        self.DefineMaterialsBulkForm.setObjectName(u"DefineMaterialsBulkForm")
        self.verticalLayout_30 = QVBoxLayout(self.DefineMaterialsBulkForm)
        self.verticalLayout_30.setObjectName(u"verticalLayout_30")
        self.verticalLayout_30.setContentsMargins(0, -1, 0, 0)
        self.MaterialInputs_2 = QFrame(self.DefineMaterialsBulkForm)
        self.MaterialInputs_2.setObjectName(u"MaterialInputs_2")
        self.MaterialInputs_2.setFrameShape(QFrame.Shape.NoFrame)
        self.gridLayout_25 = QGridLayout(self.MaterialInputs_2)
        self.gridLayout_25.setObjectName(u"gridLayout_25")
        self.gridLayout_25.setContentsMargins(0, 0, 0, 0)
        self.frame_46 = QFrame(self.MaterialInputs_2)
        self.frame_46.setObjectName(u"frame_46")
        self.frame_46.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_22 = QHBoxLayout(self.frame_46)
        self.horizontalLayout_22.setObjectName(u"horizontalLayout_22")
        self.horizontalLayout_22.setContentsMargins(0, 0, 0, 0)

        self.gridLayout_25.addWidget(self.frame_46, 0, 2, 1, 1)

        self.LMaterialName_3 = QLabel(self.MaterialInputs_2)
        self.LMaterialName_3.setObjectName(u"LMaterialName_3")

        self.gridLayout_25.addWidget(self.LMaterialName_3, 0, 0, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.INMaterialName_2 = QLineEdit(self.MaterialInputs_2)
        self.INMaterialName_2.setObjectName(u"INMaterialName_2")
        self.INMaterialName_2.setMinimumSize(QSize(93, 0))

        self.gridLayout_25.addWidget(self.INMaterialName_2, 0, 1, 1, 1)

        self.LMaterialDefinition = QLabel(self.MaterialInputs_2)
        self.LMaterialDefinition.setObjectName(u"LMaterialDefinition")

        self.gridLayout_25.addWidget(self.LMaterialDefinition, 1, 0, 1, 1)

        self.INMaterialDefinition = QComboBox(self.MaterialInputs_2)
        self.INMaterialDefinition.addItem("")
        self.INMaterialDefinition.addItem("")
        self.INMaterialDefinition.addItem("")
        self.INMaterialDefinition.addItem("")
        self.INMaterialDefinition.addItem("")
        self.INMaterialDefinition.setObjectName(u"INMaterialDefinition")

        self.gridLayout_25.addWidget(self.INMaterialDefinition, 1, 1, 1, 1)

        self.frame_47 = QFrame(self.MaterialInputs_2)
        self.frame_47.setObjectName(u"frame_47")
        self.frame_47.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_23 = QHBoxLayout(self.frame_47)
        self.horizontalLayout_23.setObjectName(u"horizontalLayout_23")
        self.horizontalLayout_23.setContentsMargins(0, 0, 0, 0)

        self.gridLayout_25.addWidget(self.frame_47, 0, 3, 1, 1)


        self.verticalLayout_30.addWidget(self.MaterialInputs_2)

        self.MaterialMenu_2 = QTabWidget(self.DefineMaterialsBulkForm)
        self.MaterialMenu_2.setObjectName(u"MaterialMenu_2")
        self.MaterialMenu_2.setEnabled(True)
        self.Elastic_2 = QWidget()
        self.Elastic_2.setObjectName(u"Elastic_2")
        self.verticalLayout_25 = QVBoxLayout(self.Elastic_2)
        self.verticalLayout_25.setObjectName(u"verticalLayout_25")
        self.verticalLayout_25.setContentsMargins(6, 12, 6, 20)
        self.frame_25 = QFrame(self.Elastic_2)
        self.frame_25.setObjectName(u"frame_25")
        self.frame_25.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_25.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_20 = QGridLayout(self.frame_25)
        self.gridLayout_20.setObjectName(u"gridLayout_20")
        self.gridLayout_20.setContentsMargins(0, 0, 0, 0)
        self.LType_2 = QLabel(self.frame_25)
        self.LType_2.setObjectName(u"LType_2")

        self.gridLayout_20.addWidget(self.LType_2, 0, 0, 1, 1)

        self.frame_44 = QFrame(self.frame_25)
        self.frame_44.setObjectName(u"frame_44")
        self.frame_44.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_20 = QHBoxLayout(self.frame_44)
        self.horizontalLayout_20.setObjectName(u"horizontalLayout_20")
        self.horizontalLayout_20.setContentsMargins(0, 0, 0, 0)
        self.INSolidGas_2 = QComboBox(self.frame_44)
        self.INSolidGas_2.addItem("")
        self.INSolidGas_2.addItem("")
        self.INSolidGas_2.setObjectName(u"INSolidGas_2")
        self.INSolidGas_2.setMinimumSize(QSize(82, 0))

        self.horizontalLayout_20.addWidget(self.INSolidGas_2)

        self.INMaterialType_2 = QComboBox(self.frame_44)
        self.INMaterialType_2.addItem("")
        self.INMaterialType_2.addItem("")
        self.INMaterialType_2.addItem("")
        self.INMaterialType_2.addItem("")
        self.INMaterialType_2.setObjectName(u"INMaterialType_2")

        self.horizontalLayout_20.addWidget(self.INMaterialType_2)


        self.gridLayout_20.addWidget(self.frame_44, 0, 1, 1, 1)


        self.verticalLayout_25.addWidget(self.frame_25, 0, Qt.AlignmentFlag.AlignHCenter)

        self.MaterialTypeTool_2 = QStackedWidget(self.Elastic_2)
        self.MaterialTypeTool_2.setObjectName(u"MaterialTypeTool_2")
        self.Isotropic_2 = QWidget()
        self.Isotropic_2.setObjectName(u"Isotropic_2")
        self.Isotropic_2.setEnabled(True)
        self.gridLayout_21 = QGridLayout(self.Isotropic_2)
        self.gridLayout_21.setObjectName(u"gridLayout_21")
        self.gridLayout_21.setContentsMargins(-1, -1, -1, 12)
        self.LYoungsModulus_2 = QLabel(self.Isotropic_2)
        self.LYoungsModulus_2.setObjectName(u"LYoungsModulus_2")

        self.gridLayout_21.addWidget(self.LYoungsModulus_2, 0, 0, 1, 1)

        self.INYoungsModulus_2 = QLineEdit(self.Isotropic_2)
        self.INYoungsModulus_2.setObjectName(u"INYoungsModulus_2")

        self.gridLayout_21.addWidget(self.INYoungsModulus_2, 0, 1, 1, 1)

        self.UPressure7 = QLabel(self.Isotropic_2)
        self.UPressure7.setObjectName(u"UPressure7")

        self.gridLayout_21.addWidget(self.UPressure7, 0, 2, 1, 1)

        self.LPoissonsRatio_2 = QLabel(self.Isotropic_2)
        self.LPoissonsRatio_2.setObjectName(u"LPoissonsRatio_2")

        self.gridLayout_21.addWidget(self.LPoissonsRatio_2, 1, 0, 1, 1)

        self.INPoissonsRatio_2 = QLineEdit(self.Isotropic_2)
        self.INPoissonsRatio_2.setObjectName(u"INPoissonsRatio_2")

        self.gridLayout_21.addWidget(self.INPoissonsRatio_2, 1, 1, 1, 1)

        self.verticalSpacer_25 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.gridLayout_21.addItem(self.verticalSpacer_25, 2, 1, 1, 1)

        self.MaterialTypeTool_2.addWidget(self.Isotropic_2)
        self.TransverselyIsotropic_2 = QWidget()
        self.TransverselyIsotropic_2.setObjectName(u"TransverselyIsotropic_2")
        self.verticalLayout_60 = QVBoxLayout(self.TransverselyIsotropic_2)
        self.verticalLayout_60.setSpacing(0)
        self.verticalLayout_60.setObjectName(u"verticalLayout_60")
        self.verticalLayout_60.setContentsMargins(0, 0, 0, 0)
        self.frame_45 = QFrame(self.TransverselyIsotropic_2)
        self.frame_45.setObjectName(u"frame_45")
        self.frame_45.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_45.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_28 = QFormLayout(self.frame_45)
        self.formLayout_28.setObjectName(u"formLayout_28")
        self.formLayout_28.setContentsMargins(0, 0, -1, 0)
        self.LIsotropicPlane_2 = QLabel(self.frame_45)
        self.LIsotropicPlane_2.setObjectName(u"LIsotropicPlane_2")

        self.formLayout_28.setWidget(0, QFormLayout.LabelRole, self.LIsotropicPlane_2)

        self.INIsotropicPlane_2 = QComboBox(self.frame_45)
        self.INIsotropicPlane_2.addItem("")
        self.INIsotropicPlane_2.addItem("")
        self.INIsotropicPlane_2.addItem("")
        self.INIsotropicPlane_2.setObjectName(u"INIsotropicPlane_2")

        self.formLayout_28.setWidget(0, QFormLayout.FieldRole, self.INIsotropicPlane_2)


        self.verticalLayout_60.addWidget(self.frame_45, 0, Qt.AlignmentFlag.AlignTop)

        self.TransverslyIsotropicMat_2 = QFrame(self.TransverselyIsotropic_2)
        self.TransverslyIsotropicMat_2.setObjectName(u"TransverslyIsotropicMat_2")
        self.TransverslyIsotropicMat_2.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_21 = QHBoxLayout(self.TransverslyIsotropicMat_2)
        self.horizontalLayout_21.setSpacing(0)
        self.horizontalLayout_21.setObjectName(u"horizontalLayout_21")
        self.horizontalLayout_21.setContentsMargins(0, 0, 0, 0)
        self.TransverseInPlane_2 = QFrame(self.TransverslyIsotropicMat_2)
        self.TransverseInPlane_2.setObjectName(u"TransverseInPlane_2")
        self.TransverseInPlane_2.setFrameShape(QFrame.Shape.NoFrame)
        self.verticalLayout_61 = QVBoxLayout(self.TransverseInPlane_2)
        self.verticalLayout_61.setSpacing(0)
        self.verticalLayout_61.setObjectName(u"verticalLayout_61")
        self.verticalLayout_61.setContentsMargins(0, 12, 0, 0)
        self.LInPlane_2 = QLabel(self.TransverseInPlane_2)
        self.LInPlane_2.setObjectName(u"LInPlane_2")

        self.verticalLayout_61.addWidget(self.LInPlane_2)

        self.TransverseInPlaneMat_2 = QFrame(self.TransverseInPlane_2)
        self.TransverseInPlaneMat_2.setObjectName(u"TransverseInPlaneMat_2")
        self.TransverseInPlaneMat_2.setFrameShape(QFrame.Shape.NoFrame)
        self.gridLayout_22 = QGridLayout(self.TransverseInPlaneMat_2)
        self.gridLayout_22.setObjectName(u"gridLayout_22")
        self.INNUip_2 = QLineEdit(self.TransverseInPlaneMat_2)
        self.INNUip_2.setObjectName(u"INNUip_2")

        self.gridLayout_22.addWidget(self.INNUip_2, 1, 1, 1, 1)

        self.LNUip_2 = QLabel(self.TransverseInPlaneMat_2)
        self.LNUip_2.setObjectName(u"LNUip_2")
        self.LNUip_2.setMinimumSize(QSize(31, 0))

        self.gridLayout_22.addWidget(self.LNUip_2, 1, 0, 1, 1)

        self.LEip_2 = QLabel(self.TransverseInPlaneMat_2)
        self.LEip_2.setObjectName(u"LEip_2")
        self.LEip_2.setMinimumSize(QSize(31, 0))

        self.gridLayout_22.addWidget(self.LEip_2, 0, 0, 1, 1)

        self.INEip_2 = QLineEdit(self.TransverseInPlaneMat_2)
        self.INEip_2.setObjectName(u"INEip_2")

        self.gridLayout_22.addWidget(self.INEip_2, 0, 1, 1, 1)

        self.UPressure10 = QLabel(self.TransverseInPlaneMat_2)
        self.UPressure10.setObjectName(u"UPressure10")

        self.gridLayout_22.addWidget(self.UPressure10, 0, 2, 1, 1)


        self.verticalLayout_61.addWidget(self.TransverseInPlaneMat_2)


        self.horizontalLayout_21.addWidget(self.TransverseInPlane_2)

        self.TransverseOutOfPlane_2 = QFrame(self.TransverslyIsotropicMat_2)
        self.TransverseOutOfPlane_2.setObjectName(u"TransverseOutOfPlane_2")
        self.TransverseOutOfPlane_2.setFrameShape(QFrame.Shape.NoFrame)
        self.verticalLayout_62 = QVBoxLayout(self.TransverseOutOfPlane_2)
        self.verticalLayout_62.setSpacing(0)
        self.verticalLayout_62.setObjectName(u"verticalLayout_62")
        self.verticalLayout_62.setContentsMargins(0, 12, 0, 0)
        self.LOutOfPlane_2 = QLabel(self.TransverseOutOfPlane_2)
        self.LOutOfPlane_2.setObjectName(u"LOutOfPlane_2")

        self.verticalLayout_62.addWidget(self.LOutOfPlane_2)

        self.TransverseOutOfPlaneMat_2 = QFrame(self.TransverseOutOfPlane_2)
        self.TransverseOutOfPlaneMat_2.setObjectName(u"TransverseOutOfPlaneMat_2")
        self.TransverseOutOfPlaneMat_2.setFrameShape(QFrame.Shape.NoFrame)
        self.gridLayout_23 = QGridLayout(self.TransverseOutOfPlaneMat_2)
        self.gridLayout_23.setObjectName(u"gridLayout_23")
        self.INNUop_2 = QLineEdit(self.TransverseOutOfPlaneMat_2)
        self.INNUop_2.setObjectName(u"INNUop_2")

        self.gridLayout_23.addWidget(self.INNUop_2, 1, 1, 1, 1)

        self.INGop_2 = QLineEdit(self.TransverseOutOfPlaneMat_2)
        self.INGop_2.setObjectName(u"INGop_2")

        self.gridLayout_23.addWidget(self.INGop_2, 2, 1, 1, 1)

        self.LEop_2 = QLabel(self.TransverseOutOfPlaneMat_2)
        self.LEop_2.setObjectName(u"LEop_2")
        self.LEop_2.setMinimumSize(QSize(31, 0))
        self.LEop_2.setMaximumSize(QSize(31, 16777215))

        self.gridLayout_23.addWidget(self.LEop_2, 0, 0, 1, 1)

        self.LNUop_2 = QLabel(self.TransverseOutOfPlaneMat_2)
        self.LNUop_2.setObjectName(u"LNUop_2")
        self.LNUop_2.setMinimumSize(QSize(31, 0))
        self.LNUop_2.setMaximumSize(QSize(31, 16777215))

        self.gridLayout_23.addWidget(self.LNUop_2, 1, 0, 1, 1)

        self.LGop_2 = QLabel(self.TransverseOutOfPlaneMat_2)
        self.LGop_2.setObjectName(u"LGop_2")
        self.LGop_2.setMinimumSize(QSize(31, 0))
        self.LGop_2.setMaximumSize(QSize(31, 16777215))

        self.gridLayout_23.addWidget(self.LGop_2, 2, 0, 1, 1)

        self.INEop_2 = QLineEdit(self.TransverseOutOfPlaneMat_2)
        self.INEop_2.setObjectName(u"INEop_2")
        self.INEop_2.setMinimumSize(QSize(0, 0))

        self.gridLayout_23.addWidget(self.INEop_2, 0, 1, 1, 1)

        self.UPressure11 = QLabel(self.TransverseOutOfPlaneMat_2)
        self.UPressure11.setObjectName(u"UPressure11")

        self.gridLayout_23.addWidget(self.UPressure11, 0, 2, 1, 1)

        self.UPressure12 = QLabel(self.TransverseOutOfPlaneMat_2)
        self.UPressure12.setObjectName(u"UPressure12")

        self.gridLayout_23.addWidget(self.UPressure12, 2, 2, 1, 1)


        self.verticalLayout_62.addWidget(self.TransverseOutOfPlaneMat_2)


        self.horizontalLayout_21.addWidget(self.TransverseOutOfPlane_2)


        self.verticalLayout_60.addWidget(self.TransverslyIsotropicMat_2)

        self.verticalSpacer_26 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_60.addItem(self.verticalSpacer_26)

        self.MaterialTypeTool_2.addWidget(self.TransverselyIsotropic_2)
        self.Anisotropic_2 = QWidget()
        self.Anisotropic_2.setObjectName(u"Anisotropic_2")
        self.verticalLayout_63 = QVBoxLayout(self.Anisotropic_2)
        self.verticalLayout_63.setObjectName(u"verticalLayout_63")
        self.verticalLayout_63.setContentsMargins(0, 0, 0, 0)
        self.label_9 = QLabel(self.Anisotropic_2)
        self.label_9.setObjectName(u"label_9")

        self.verticalLayout_63.addWidget(self.label_9)

        self.TAnisotropic_2 = QTableWidget(self.Anisotropic_2)
        if (self.TAnisotropic_2.columnCount() < 6):
            self.TAnisotropic_2.setColumnCount(6)
        if (self.TAnisotropic_2.rowCount() < 6):
            self.TAnisotropic_2.setRowCount(6)
        __qtablewidgetitem72 = QTableWidgetItem()
        self.TAnisotropic_2.setItem(0, 0, __qtablewidgetitem72)
        __qtablewidgetitem73 = QTableWidgetItem()
        __qtablewidgetitem73.setBackground(brush);
        __qtablewidgetitem73.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(1, 0, __qtablewidgetitem73)
        __qtablewidgetitem74 = QTableWidgetItem()
        __qtablewidgetitem74.setBackground(brush);
        __qtablewidgetitem74.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(2, 0, __qtablewidgetitem74)
        __qtablewidgetitem75 = QTableWidgetItem()
        __qtablewidgetitem75.setBackground(brush);
        __qtablewidgetitem75.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(2, 1, __qtablewidgetitem75)
        __qtablewidgetitem76 = QTableWidgetItem()
        __qtablewidgetitem76.setBackground(brush);
        __qtablewidgetitem76.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(3, 0, __qtablewidgetitem76)
        __qtablewidgetitem77 = QTableWidgetItem()
        __qtablewidgetitem77.setBackground(brush);
        __qtablewidgetitem77.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(3, 1, __qtablewidgetitem77)
        __qtablewidgetitem78 = QTableWidgetItem()
        __qtablewidgetitem78.setBackground(brush);
        __qtablewidgetitem78.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(3, 2, __qtablewidgetitem78)
        __qtablewidgetitem79 = QTableWidgetItem()
        __qtablewidgetitem79.setBackground(brush);
        __qtablewidgetitem79.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(4, 0, __qtablewidgetitem79)
        __qtablewidgetitem80 = QTableWidgetItem()
        __qtablewidgetitem80.setBackground(brush);
        __qtablewidgetitem80.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(4, 1, __qtablewidgetitem80)
        __qtablewidgetitem81 = QTableWidgetItem()
        __qtablewidgetitem81.setBackground(brush);
        __qtablewidgetitem81.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(4, 2, __qtablewidgetitem81)
        __qtablewidgetitem82 = QTableWidgetItem()
        __qtablewidgetitem82.setBackground(brush);
        __qtablewidgetitem82.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(4, 3, __qtablewidgetitem82)
        __qtablewidgetitem83 = QTableWidgetItem()
        __qtablewidgetitem83.setBackground(brush);
        __qtablewidgetitem83.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(5, 0, __qtablewidgetitem83)
        __qtablewidgetitem84 = QTableWidgetItem()
        __qtablewidgetitem84.setBackground(brush);
        __qtablewidgetitem84.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(5, 1, __qtablewidgetitem84)
        __qtablewidgetitem85 = QTableWidgetItem()
        __qtablewidgetitem85.setBackground(brush);
        __qtablewidgetitem85.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(5, 2, __qtablewidgetitem85)
        __qtablewidgetitem86 = QTableWidgetItem()
        __qtablewidgetitem86.setBackground(brush);
        __qtablewidgetitem86.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(5, 3, __qtablewidgetitem86)
        __qtablewidgetitem87 = QTableWidgetItem()
        __qtablewidgetitem87.setBackground(brush);
        __qtablewidgetitem87.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic_2.setItem(5, 4, __qtablewidgetitem87)
        self.TAnisotropic_2.setObjectName(u"TAnisotropic_2")
        self.TAnisotropic_2.setRowCount(6)
        self.TAnisotropic_2.setColumnCount(6)
        self.TAnisotropic_2.horizontalHeader().setMinimumSectionSize(21)
        self.TAnisotropic_2.horizontalHeader().setDefaultSectionSize(50)

        self.verticalLayout_63.addWidget(self.TAnisotropic_2)

        self.verticalSpacer_32 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_63.addItem(self.verticalSpacer_32)

        self.MaterialTypeTool_2.addWidget(self.Anisotropic_2)
        self.Orthotropic_2 = QWidget()
        self.Orthotropic_2.setObjectName(u"Orthotropic_2")
        self.gridLayout_24 = QGridLayout(self.Orthotropic_2)
        self.gridLayout_24.setObjectName(u"gridLayout_24")
        self.gridLayout_24.setHorizontalSpacing(6)
        self.gridLayout_24.setVerticalSpacing(0)
        self.gridLayout_24.setContentsMargins(0, 0, 0, 0)
        self.INNUyz_2 = QLineEdit(self.Orthotropic_2)
        self.INNUyz_2.setObjectName(u"INNUyz_2")

        self.gridLayout_24.addWidget(self.INNUyz_2, 2, 6, 1, 1)

        self.LNUxz_2 = QLabel(self.Orthotropic_2)
        self.LNUxz_2.setObjectName(u"LNUxz_2")

        self.gridLayout_24.addWidget(self.LNUxz_2, 2, 3, 1, 1)

        self.LGyz_2 = QLabel(self.Orthotropic_2)
        self.LGyz_2.setObjectName(u"LGyz_2")

        self.gridLayout_24.addWidget(self.LGyz_2, 3, 5, 1, 1)

        self.INEx_2 = QLineEdit(self.Orthotropic_2)
        self.INEx_2.setObjectName(u"INEx_2")

        self.gridLayout_24.addWidget(self.INEx_2, 0, 2, 1, 1)

        self.LGxz_2 = QLabel(self.Orthotropic_2)
        self.LGxz_2.setObjectName(u"LGxz_2")

        self.gridLayout_24.addWidget(self.LGxz_2, 3, 3, 1, 1)

        self.LGxy_2 = QLabel(self.Orthotropic_2)
        self.LGxy_2.setObjectName(u"LGxy_2")

        self.gridLayout_24.addWidget(self.LGxy_2, 3, 1, 1, 1)

        self.INGyz_2 = QLineEdit(self.Orthotropic_2)
        self.INGyz_2.setObjectName(u"INGyz_2")

        self.gridLayout_24.addWidget(self.INGyz_2, 3, 6, 1, 1)

        self.LNUyz_2 = QLabel(self.Orthotropic_2)
        self.LNUyz_2.setObjectName(u"LNUyz_2")

        self.gridLayout_24.addWidget(self.LNUyz_2, 2, 5, 1, 1)

        self.INNUxy_2 = QLineEdit(self.Orthotropic_2)
        self.INNUxy_2.setObjectName(u"INNUxy_2")

        self.gridLayout_24.addWidget(self.INNUxy_2, 2, 2, 1, 1)

        self.INNUxz_2 = QLineEdit(self.Orthotropic_2)
        self.INNUxz_2.setObjectName(u"INNUxz_2")

        self.gridLayout_24.addWidget(self.INNUxz_2, 2, 4, 1, 1)

        self.INEz_2 = QLineEdit(self.Orthotropic_2)
        self.INEz_2.setObjectName(u"INEz_2")

        self.gridLayout_24.addWidget(self.INEz_2, 0, 6, 1, 1)

        self.INEy_2 = QLineEdit(self.Orthotropic_2)
        self.INEy_2.setObjectName(u"INEy_2")

        self.gridLayout_24.addWidget(self.INEy_2, 0, 4, 1, 1)

        self.LEx_2 = QLabel(self.Orthotropic_2)
        self.LEx_2.setObjectName(u"LEx_2")

        self.gridLayout_24.addWidget(self.LEx_2, 0, 1, 1, 1)

        self.INGxy_2 = QLineEdit(self.Orthotropic_2)
        self.INGxy_2.setObjectName(u"INGxy_2")

        self.gridLayout_24.addWidget(self.INGxy_2, 3, 2, 1, 1)

        self.UPressure8 = QLabel(self.Orthotropic_2)
        self.UPressure8.setObjectName(u"UPressure8")

        self.gridLayout_24.addWidget(self.UPressure8, 0, 7, 1, 1)

        self.UPressure9 = QLabel(self.Orthotropic_2)
        self.UPressure9.setObjectName(u"UPressure9")

        self.gridLayout_24.addWidget(self.UPressure9, 3, 7, 1, 1)

        self.INGxz_2 = QLineEdit(self.Orthotropic_2)
        self.INGxz_2.setObjectName(u"INGxz_2")

        self.gridLayout_24.addWidget(self.INGxz_2, 3, 4, 1, 1)

        self.LEy_2 = QLabel(self.Orthotropic_2)
        self.LEy_2.setObjectName(u"LEy_2")

        self.gridLayout_24.addWidget(self.LEy_2, 0, 3, 1, 1)

        self.LNUxy_2 = QLabel(self.Orthotropic_2)
        self.LNUxy_2.setObjectName(u"LNUxy_2")

        self.gridLayout_24.addWidget(self.LNUxy_2, 2, 1, 1, 1)

        self.LEz_2 = QLabel(self.Orthotropic_2)
        self.LEz_2.setObjectName(u"LEz_2")

        self.gridLayout_24.addWidget(self.LEz_2, 0, 5, 1, 1)

        self.verticalSpacer_33 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.gridLayout_24.addItem(self.verticalSpacer_33, 4, 4, 1, 1)

        self.MaterialTypeTool_2.addWidget(self.Orthotropic_2)
        self.page_12 = QWidget()
        self.page_12.setObjectName(u"page_12")
        self.MaterialTypeTool_2.addWidget(self.page_12)

        self.verticalLayout_25.addWidget(self.MaterialTypeTool_2)

        self.MaterialMenu_2.addTab(self.Elastic_2, icon11, "")
        self.Plastic = QWidget()
        self.Plastic.setObjectName(u"Plastic")
        self.Plastic.setStyleSheet(u"")
        self.verticalLayout_26 = QVBoxLayout(self.Plastic)
#ifndef Q_OS_MAC
        self.verticalLayout_26.setSpacing(-1)
#endif
        self.verticalLayout_26.setObjectName(u"verticalLayout_26")
        self.verticalLayout_26.setContentsMargins(0, 0, 0, 0)
        self.BEnablePlasticity = QCheckBox(self.Plastic)
        self.BEnablePlasticity.setObjectName(u"BEnablePlasticity")

        self.verticalLayout_26.addWidget(self.BEnablePlasticity, 0, Qt.AlignmentFlag.AlignHCenter)

        self.EnablePlasticity = QStackedWidget(self.Plastic)
        self.EnablePlasticity.setObjectName(u"EnablePlasticity")
        sizePolicy13.setHeightForWidth(self.EnablePlasticity.sizePolicy().hasHeightForWidth())
        self.EnablePlasticity.setSizePolicy(sizePolicy13)
        self.page_10 = QWidget()
        self.page_10.setObjectName(u"page_10")
        self.verticalLayout_71 = QVBoxLayout(self.page_10)
        self.verticalLayout_71.setObjectName(u"verticalLayout_71")
        self.EnablePlasticity.addWidget(self.page_10)
        self.page_26 = QWidget()
        self.page_26.setObjectName(u"page_26")
        self.verticalLayout_70 = QVBoxLayout(self.page_26)
        self.verticalLayout_70.setObjectName(u"verticalLayout_70")
        self.verticalLayout_70.setContentsMargins(0, 0, 0, 0)
        self.PlasticProperties = QTabWidget(self.page_26)
        self.PlasticProperties.setObjectName(u"PlasticProperties")
        self.PlasticProperties.setMinimumSize(QSize(0, 0))
        font13 = QFont()
        font13.setItalic(False)
        self.PlasticProperties.setFont(font13)
        self.PlasticProperties.setIconSize(QSize(16, 16))
        self.PlasticProperties.setMovable(False)
        self.PlasticProperties.setTabBarAutoHide(False)
        self.CrystalAxis = QWidget()
        self.CrystalAxis.setObjectName(u"CrystalAxis")
        self.CrystalAxis.setStyleSheet(u"")
        self.verticalLayout_65 = QVBoxLayout(self.CrystalAxis)
        self.verticalLayout_65.setObjectName(u"verticalLayout_65")
        self.frame_58 = QFrame(self.CrystalAxis)
        self.frame_58.setObjectName(u"frame_58")
        self.frame_58.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_58.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_33 = QGridLayout(self.frame_58)
        self.gridLayout_33.setObjectName(u"gridLayout_33")
        self.gridLayout_33.setContentsMargins(0, 0, 0, 0)
        self.Lc = QLabel(self.frame_58)
        self.Lc.setObjectName(u"Lc")

        self.gridLayout_33.addWidget(self.Lc, 0, 5, 1, 1)

        self.La = QLabel(self.frame_58)
        self.La.setObjectName(u"La")

        self.gridLayout_33.addWidget(self.La, 0, 0, 1, 1)

        self.Lb = QLabel(self.frame_58)
        self.Lb.setObjectName(u"Lb")

        self.gridLayout_33.addWidget(self.Lb, 0, 2, 1, 1)

        self.INc = QLineEdit(self.frame_58)
        self.INc.setObjectName(u"INc")

        self.gridLayout_33.addWidget(self.INc, 0, 6, 1, 1)

        self.INb = QLineEdit(self.frame_58)
        self.INb.setObjectName(u"INb")

        self.gridLayout_33.addWidget(self.INb, 0, 3, 1, 1)

        self.INa = QLineEdit(self.frame_58)
        self.INa.setObjectName(u"INa")

        self.gridLayout_33.addWidget(self.INa, 0, 1, 1, 1)

        self.pushButton_13 = QPushButton(self.frame_58)
        self.pushButton_13.setObjectName(u"pushButton_13")
        self.pushButton_13.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton_13.sizePolicy().hasHeightForWidth())
        self.pushButton_13.setSizePolicy(sizePolicy9)
        self.pushButton_13.setMaximumSize(QSize(17, 17))
        self.pushButton_13.setFont(font6)
        self.pushButton_13.setIcon(icon10)
        self.pushButton_13.setIconSize(QSize(20, 20))
        self.pushButton_13.setAutoDefault(False)
        self.pushButton_13.setFlat(True)

        self.gridLayout_33.addWidget(self.pushButton_13, 0, 7, 1, 1)


        self.verticalLayout_65.addWidget(self.frame_58, 0, Qt.AlignmentFlag.AlignTop)

        self.PlasticProperties.addTab(self.CrystalAxis, "")
        self.SlipSystems = QWidget()
        self.SlipSystems.setObjectName(u"SlipSystems")
        self.verticalLayout_66 = QVBoxLayout(self.SlipSystems)
        self.verticalLayout_66.setSpacing(0)
        self.verticalLayout_66.setObjectName(u"verticalLayout_66")
        self.verticalLayout_66.setContentsMargins(0, 0, 0, 0)
        self.TSlipSystems = QTreeWidget(self.SlipSystems)
        self.TSlipSystems.headerItem().setText(0, "")
        __qtreewidgetitem = QTreeWidgetItem(self.TSlipSystems)
        QTreeWidgetItem(__qtreewidgetitem)
        __qtreewidgetitem1 = QTreeWidgetItem(self.TSlipSystems)
        QTreeWidgetItem(__qtreewidgetitem1)
        QTreeWidgetItem(__qtreewidgetitem1)
        QTreeWidgetItem(__qtreewidgetitem1)
        QTreeWidgetItem(self.TSlipSystems)
        self.TSlipSystems.setObjectName(u"TSlipSystems")

        self.verticalLayout_66.addWidget(self.TSlipSystems)

        self.frame_56 = QFrame(self.SlipSystems)
        self.frame_56.setObjectName(u"frame_56")
        self.frame_56.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_56.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_7 = QHBoxLayout(self.frame_56)
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.horizontalLayout_7.setContentsMargins(0, 0, 0, 0)
        self.BSlipSystemDetails = QPushButton(self.frame_56)
        self.BSlipSystemDetails.setObjectName(u"BSlipSystemDetails")

        self.horizontalLayout_7.addWidget(self.BSlipSystemDetails)

        self.BCustomSlipSystem = QPushButton(self.frame_56)
        self.BCustomSlipSystem.setObjectName(u"BCustomSlipSystem")

        self.horizontalLayout_7.addWidget(self.BCustomSlipSystem)


        self.verticalLayout_66.addWidget(self.frame_56)

        self.INSlipSystems = QListWidget(self.SlipSystems)
        self.INSlipSystems.setObjectName(u"INSlipSystems")
        sizePolicy14 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        sizePolicy14.setHorizontalStretch(0)
        sizePolicy14.setVerticalStretch(0)
        sizePolicy14.setHeightForWidth(self.INSlipSystems.sizePolicy().hasHeightForWidth())
        self.INSlipSystems.setSizePolicy(sizePolicy14)
        self.INSlipSystems.setAutoFillBackground(False)
        self.INSlipSystems.setStyleSheet(u"background-color: rgb(194, 194, 194)")
        self.INSlipSystems.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustToContentsOnFirstShow)
        self.INSlipSystems.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)

        self.verticalLayout_66.addWidget(self.INSlipSystems)

        self.frame_50 = QFrame(self.SlipSystems)
        self.frame_50.setObjectName(u"frame_50")
        self.frame_50.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_50.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_8 = QHBoxLayout(self.frame_50)
#ifndef Q_OS_MAC
        self.horizontalLayout_8.setSpacing(-1)
#endif
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.horizontalLayout_8.setContentsMargins(0, 0, 0, 0)
        self.BAddSlipSystem = QPushButton(self.frame_50)
        self.BAddSlipSystem.setObjectName(u"BAddSlipSystem")

        self.horizontalLayout_8.addWidget(self.BAddSlipSystem)

        self.BRemoveSlipSystem = QPushButton(self.frame_50)
        self.BRemoveSlipSystem.setObjectName(u"BRemoveSlipSystem")

        self.horizontalLayout_8.addWidget(self.BRemoveSlipSystem)


        self.verticalLayout_66.addWidget(self.frame_50)

        self.PlasticProperties.addTab(self.SlipSystems, "")
        self.Details = QWidget()
        self.Details.setObjectName(u"Details")
        self.verticalLayout_42 = QVBoxLayout(self.Details)
        self.verticalLayout_42.setObjectName(u"verticalLayout_42")
        self.verticalLayout_42.setContentsMargins(0, 0, 0, 0)
        self.SlipSystemInfo = QStackedWidget(self.Details)
        self.SlipSystemInfo.setObjectName(u"SlipSystemInfo")
        self.SlipSystemInfo.setStyleSheet(u"")
        self.SlipSystemInfo.setFrameShape(QFrame.Shape.NoFrame)
        self.SlipSystemInfo.setLineWidth(0)
        self.BlankPage = QWidget()
        self.BlankPage.setObjectName(u"BlankPage")
        self.SlipSystemInfo.addWidget(self.BlankPage)
        self.FCC_111_110 = QWidget()
        self.FCC_111_110.setObjectName(u"FCC_111_110")
        self.verticalLayout_43 = QVBoxLayout(self.FCC_111_110)
        self.verticalLayout_43.setSpacing(5)
        self.verticalLayout_43.setObjectName(u"verticalLayout_43")
        self.verticalLayout_43.setContentsMargins(0, 0, 0, 0)
        self.label_42 = QLabel(self.FCC_111_110)
        self.label_42.setObjectName(u"label_42")
        font14 = QFont()
        font14.setBold(True)
        font14.setUnderline(True)
        self.label_42.setFont(font14)

        self.verticalLayout_43.addWidget(self.label_42, 0, Qt.AlignmentFlag.AlignHCenter)

        self.T1 = QTableWidget(self.FCC_111_110)
        if (self.T1.columnCount() < 2):
            self.T1.setColumnCount(2)
        font15 = QFont()
        font15.setPointSize(13)
        __qtablewidgetitem88 = QTableWidgetItem()
        __qtablewidgetitem88.setFont(font15);
        self.T1.setHorizontalHeaderItem(0, __qtablewidgetitem88)
        __qtablewidgetitem89 = QTableWidgetItem()
        __qtablewidgetitem89.setFont(font15);
        self.T1.setHorizontalHeaderItem(1, __qtablewidgetitem89)
        if (self.T1.rowCount() < 12):
            self.T1.setRowCount(12)
        __qtablewidgetitem90 = QTableWidgetItem()
        self.T1.setItem(0, 0, __qtablewidgetitem90)
        __qtablewidgetitem91 = QTableWidgetItem()
        self.T1.setItem(0, 1, __qtablewidgetitem91)
        __qtablewidgetitem92 = QTableWidgetItem()
        self.T1.setItem(1, 0, __qtablewidgetitem92)
        __qtablewidgetitem93 = QTableWidgetItem()
        self.T1.setItem(1, 1, __qtablewidgetitem93)
        __qtablewidgetitem94 = QTableWidgetItem()
        self.T1.setItem(2, 0, __qtablewidgetitem94)
        __qtablewidgetitem95 = QTableWidgetItem()
        self.T1.setItem(2, 1, __qtablewidgetitem95)
        __qtablewidgetitem96 = QTableWidgetItem()
        self.T1.setItem(3, 0, __qtablewidgetitem96)
        __qtablewidgetitem97 = QTableWidgetItem()
        self.T1.setItem(3, 1, __qtablewidgetitem97)
        __qtablewidgetitem98 = QTableWidgetItem()
        self.T1.setItem(4, 0, __qtablewidgetitem98)
        __qtablewidgetitem99 = QTableWidgetItem()
        self.T1.setItem(4, 1, __qtablewidgetitem99)
        __qtablewidgetitem100 = QTableWidgetItem()
        self.T1.setItem(5, 0, __qtablewidgetitem100)
        __qtablewidgetitem101 = QTableWidgetItem()
        self.T1.setItem(5, 1, __qtablewidgetitem101)
        __qtablewidgetitem102 = QTableWidgetItem()
        self.T1.setItem(6, 0, __qtablewidgetitem102)
        __qtablewidgetitem103 = QTableWidgetItem()
        self.T1.setItem(6, 1, __qtablewidgetitem103)
        __qtablewidgetitem104 = QTableWidgetItem()
        self.T1.setItem(7, 0, __qtablewidgetitem104)
        __qtablewidgetitem105 = QTableWidgetItem()
        self.T1.setItem(7, 1, __qtablewidgetitem105)
        __qtablewidgetitem106 = QTableWidgetItem()
        self.T1.setItem(8, 0, __qtablewidgetitem106)
        __qtablewidgetitem107 = QTableWidgetItem()
        self.T1.setItem(8, 1, __qtablewidgetitem107)
        __qtablewidgetitem108 = QTableWidgetItem()
        self.T1.setItem(9, 0, __qtablewidgetitem108)
        __qtablewidgetitem109 = QTableWidgetItem()
        self.T1.setItem(9, 1, __qtablewidgetitem109)
        __qtablewidgetitem110 = QTableWidgetItem()
        self.T1.setItem(10, 0, __qtablewidgetitem110)
        __qtablewidgetitem111 = QTableWidgetItem()
        self.T1.setItem(10, 1, __qtablewidgetitem111)
        __qtablewidgetitem112 = QTableWidgetItem()
        self.T1.setItem(11, 0, __qtablewidgetitem112)
        __qtablewidgetitem113 = QTableWidgetItem()
        self.T1.setItem(11, 1, __qtablewidgetitem113)
        self.T1.setObjectName(u"T1")
        self.T1.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.T1.horizontalHeader().setDefaultSectionSize(167)

        self.verticalLayout_43.addWidget(self.T1)

        self.SlipSystemInfo.addWidget(self.FCC_111_110)
        self.BCC_110_111 = QWidget()
        self.BCC_110_111.setObjectName(u"BCC_110_111")
        self.verticalLayout_44 = QVBoxLayout(self.BCC_110_111)
        self.verticalLayout_44.setSpacing(5)
        self.verticalLayout_44.setObjectName(u"verticalLayout_44")
        self.verticalLayout_44.setContentsMargins(0, 0, 0, 0)
        self.label_44 = QLabel(self.BCC_110_111)
        self.label_44.setObjectName(u"label_44")
        self.label_44.setFont(font14)

        self.verticalLayout_44.addWidget(self.label_44, 0, Qt.AlignmentFlag.AlignHCenter)

        self.T2 = QTableWidget(self.BCC_110_111)
        if (self.T2.columnCount() < 2):
            self.T2.setColumnCount(2)
        __qtablewidgetitem114 = QTableWidgetItem()
        __qtablewidgetitem114.setFont(font15);
        self.T2.setHorizontalHeaderItem(0, __qtablewidgetitem114)
        __qtablewidgetitem115 = QTableWidgetItem()
        __qtablewidgetitem115.setFont(font15);
        self.T2.setHorizontalHeaderItem(1, __qtablewidgetitem115)
        if (self.T2.rowCount() < 12):
            self.T2.setRowCount(12)
        __qtablewidgetitem116 = QTableWidgetItem()
        self.T2.setItem(0, 0, __qtablewidgetitem116)
        __qtablewidgetitem117 = QTableWidgetItem()
        self.T2.setItem(0, 1, __qtablewidgetitem117)
        __qtablewidgetitem118 = QTableWidgetItem()
        self.T2.setItem(1, 0, __qtablewidgetitem118)
        __qtablewidgetitem119 = QTableWidgetItem()
        self.T2.setItem(1, 1, __qtablewidgetitem119)
        __qtablewidgetitem120 = QTableWidgetItem()
        self.T2.setItem(2, 0, __qtablewidgetitem120)
        __qtablewidgetitem121 = QTableWidgetItem()
        self.T2.setItem(2, 1, __qtablewidgetitem121)
        __qtablewidgetitem122 = QTableWidgetItem()
        self.T2.setItem(3, 0, __qtablewidgetitem122)
        __qtablewidgetitem123 = QTableWidgetItem()
        self.T2.setItem(3, 1, __qtablewidgetitem123)
        __qtablewidgetitem124 = QTableWidgetItem()
        self.T2.setItem(4, 0, __qtablewidgetitem124)
        __qtablewidgetitem125 = QTableWidgetItem()
        self.T2.setItem(4, 1, __qtablewidgetitem125)
        __qtablewidgetitem126 = QTableWidgetItem()
        self.T2.setItem(5, 0, __qtablewidgetitem126)
        __qtablewidgetitem127 = QTableWidgetItem()
        self.T2.setItem(5, 1, __qtablewidgetitem127)
        __qtablewidgetitem128 = QTableWidgetItem()
        self.T2.setItem(6, 0, __qtablewidgetitem128)
        __qtablewidgetitem129 = QTableWidgetItem()
        self.T2.setItem(6, 1, __qtablewidgetitem129)
        __qtablewidgetitem130 = QTableWidgetItem()
        self.T2.setItem(7, 0, __qtablewidgetitem130)
        __qtablewidgetitem131 = QTableWidgetItem()
        self.T2.setItem(7, 1, __qtablewidgetitem131)
        __qtablewidgetitem132 = QTableWidgetItem()
        self.T2.setItem(8, 0, __qtablewidgetitem132)
        __qtablewidgetitem133 = QTableWidgetItem()
        self.T2.setItem(8, 1, __qtablewidgetitem133)
        __qtablewidgetitem134 = QTableWidgetItem()
        self.T2.setItem(9, 0, __qtablewidgetitem134)
        __qtablewidgetitem135 = QTableWidgetItem()
        self.T2.setItem(9, 1, __qtablewidgetitem135)
        __qtablewidgetitem136 = QTableWidgetItem()
        self.T2.setItem(10, 0, __qtablewidgetitem136)
        __qtablewidgetitem137 = QTableWidgetItem()
        self.T2.setItem(10, 1, __qtablewidgetitem137)
        __qtablewidgetitem138 = QTableWidgetItem()
        self.T2.setItem(11, 0, __qtablewidgetitem138)
        __qtablewidgetitem139 = QTableWidgetItem()
        self.T2.setItem(11, 1, __qtablewidgetitem139)
        self.T2.setObjectName(u"T2")
        self.T2.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.T2.horizontalHeader().setDefaultSectionSize(167)

        self.verticalLayout_44.addWidget(self.T2)

        self.SlipSystemInfo.addWidget(self.BCC_110_111)
        self.BCC_112_111 = QWidget()
        self.BCC_112_111.setObjectName(u"BCC_112_111")
        self.verticalLayout_68 = QVBoxLayout(self.BCC_112_111)
        self.verticalLayout_68.setSpacing(5)
        self.verticalLayout_68.setObjectName(u"verticalLayout_68")
        self.verticalLayout_68.setContentsMargins(0, 0, 0, 0)
        self.label_45 = QLabel(self.BCC_112_111)
        self.label_45.setObjectName(u"label_45")
        self.label_45.setFont(font14)

        self.verticalLayout_68.addWidget(self.label_45, 0, Qt.AlignmentFlag.AlignHCenter)

        self.T3 = QTableWidget(self.BCC_112_111)
        if (self.T3.columnCount() < 2):
            self.T3.setColumnCount(2)
        __qtablewidgetitem140 = QTableWidgetItem()
        __qtablewidgetitem140.setFont(font15);
        self.T3.setHorizontalHeaderItem(0, __qtablewidgetitem140)
        __qtablewidgetitem141 = QTableWidgetItem()
        __qtablewidgetitem141.setFont(font15);
        self.T3.setHorizontalHeaderItem(1, __qtablewidgetitem141)
        if (self.T3.rowCount() < 12):
            self.T3.setRowCount(12)
        __qtablewidgetitem142 = QTableWidgetItem()
        self.T3.setItem(0, 0, __qtablewidgetitem142)
        __qtablewidgetitem143 = QTableWidgetItem()
        self.T3.setItem(0, 1, __qtablewidgetitem143)
        __qtablewidgetitem144 = QTableWidgetItem()
        self.T3.setItem(1, 0, __qtablewidgetitem144)
        __qtablewidgetitem145 = QTableWidgetItem()
        self.T3.setItem(1, 1, __qtablewidgetitem145)
        __qtablewidgetitem146 = QTableWidgetItem()
        self.T3.setItem(2, 0, __qtablewidgetitem146)
        __qtablewidgetitem147 = QTableWidgetItem()
        self.T3.setItem(2, 1, __qtablewidgetitem147)
        __qtablewidgetitem148 = QTableWidgetItem()
        self.T3.setItem(3, 0, __qtablewidgetitem148)
        __qtablewidgetitem149 = QTableWidgetItem()
        self.T3.setItem(3, 1, __qtablewidgetitem149)
        __qtablewidgetitem150 = QTableWidgetItem()
        self.T3.setItem(4, 0, __qtablewidgetitem150)
        __qtablewidgetitem151 = QTableWidgetItem()
        self.T3.setItem(4, 1, __qtablewidgetitem151)
        __qtablewidgetitem152 = QTableWidgetItem()
        self.T3.setItem(5, 0, __qtablewidgetitem152)
        __qtablewidgetitem153 = QTableWidgetItem()
        self.T3.setItem(5, 1, __qtablewidgetitem153)
        __qtablewidgetitem154 = QTableWidgetItem()
        self.T3.setItem(6, 0, __qtablewidgetitem154)
        __qtablewidgetitem155 = QTableWidgetItem()
        self.T3.setItem(6, 1, __qtablewidgetitem155)
        __qtablewidgetitem156 = QTableWidgetItem()
        self.T3.setItem(7, 0, __qtablewidgetitem156)
        __qtablewidgetitem157 = QTableWidgetItem()
        self.T3.setItem(7, 1, __qtablewidgetitem157)
        __qtablewidgetitem158 = QTableWidgetItem()
        self.T3.setItem(8, 0, __qtablewidgetitem158)
        __qtablewidgetitem159 = QTableWidgetItem()
        self.T3.setItem(8, 1, __qtablewidgetitem159)
        __qtablewidgetitem160 = QTableWidgetItem()
        self.T3.setItem(9, 0, __qtablewidgetitem160)
        __qtablewidgetitem161 = QTableWidgetItem()
        self.T3.setItem(9, 1, __qtablewidgetitem161)
        __qtablewidgetitem162 = QTableWidgetItem()
        self.T3.setItem(10, 0, __qtablewidgetitem162)
        __qtablewidgetitem163 = QTableWidgetItem()
        self.T3.setItem(10, 1, __qtablewidgetitem163)
        __qtablewidgetitem164 = QTableWidgetItem()
        self.T3.setItem(11, 0, __qtablewidgetitem164)
        __qtablewidgetitem165 = QTableWidgetItem()
        self.T3.setItem(11, 1, __qtablewidgetitem165)
        self.T3.setObjectName(u"T3")
        self.T3.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.T3.horizontalHeader().setDefaultSectionSize(167)

        self.verticalLayout_68.addWidget(self.T3)

        self.SlipSystemInfo.addWidget(self.BCC_112_111)
        self.BCC_123_111 = QWidget()
        self.BCC_123_111.setObjectName(u"BCC_123_111")
        self.verticalLayout_69 = QVBoxLayout(self.BCC_123_111)
        self.verticalLayout_69.setSpacing(5)
        self.verticalLayout_69.setObjectName(u"verticalLayout_69")
        self.verticalLayout_69.setContentsMargins(0, 0, 0, 0)
        self.label_47 = QLabel(self.BCC_123_111)
        self.label_47.setObjectName(u"label_47")
        self.label_47.setFont(font14)

        self.verticalLayout_69.addWidget(self.label_47, 0, Qt.AlignmentFlag.AlignHCenter)

        self.T4 = QTableWidget(self.BCC_123_111)
        if (self.T4.columnCount() < 2):
            self.T4.setColumnCount(2)
        __qtablewidgetitem166 = QTableWidgetItem()
        __qtablewidgetitem166.setFont(font15);
        self.T4.setHorizontalHeaderItem(0, __qtablewidgetitem166)
        __qtablewidgetitem167 = QTableWidgetItem()
        __qtablewidgetitem167.setFont(font15);
        self.T4.setHorizontalHeaderItem(1, __qtablewidgetitem167)
        if (self.T4.rowCount() < 24):
            self.T4.setRowCount(24)
        __qtablewidgetitem168 = QTableWidgetItem()
        self.T4.setItem(0, 0, __qtablewidgetitem168)
        __qtablewidgetitem169 = QTableWidgetItem()
        self.T4.setItem(0, 1, __qtablewidgetitem169)
        __qtablewidgetitem170 = QTableWidgetItem()
        self.T4.setItem(1, 0, __qtablewidgetitem170)
        __qtablewidgetitem171 = QTableWidgetItem()
        self.T4.setItem(1, 1, __qtablewidgetitem171)
        __qtablewidgetitem172 = QTableWidgetItem()
        self.T4.setItem(2, 0, __qtablewidgetitem172)
        __qtablewidgetitem173 = QTableWidgetItem()
        self.T4.setItem(2, 1, __qtablewidgetitem173)
        __qtablewidgetitem174 = QTableWidgetItem()
        self.T4.setItem(3, 0, __qtablewidgetitem174)
        __qtablewidgetitem175 = QTableWidgetItem()
        self.T4.setItem(3, 1, __qtablewidgetitem175)
        __qtablewidgetitem176 = QTableWidgetItem()
        self.T4.setItem(4, 0, __qtablewidgetitem176)
        __qtablewidgetitem177 = QTableWidgetItem()
        self.T4.setItem(4, 1, __qtablewidgetitem177)
        __qtablewidgetitem178 = QTableWidgetItem()
        self.T4.setItem(5, 0, __qtablewidgetitem178)
        __qtablewidgetitem179 = QTableWidgetItem()
        self.T4.setItem(5, 1, __qtablewidgetitem179)
        __qtablewidgetitem180 = QTableWidgetItem()
        self.T4.setItem(6, 0, __qtablewidgetitem180)
        __qtablewidgetitem181 = QTableWidgetItem()
        self.T4.setItem(6, 1, __qtablewidgetitem181)
        __qtablewidgetitem182 = QTableWidgetItem()
        self.T4.setItem(7, 0, __qtablewidgetitem182)
        __qtablewidgetitem183 = QTableWidgetItem()
        self.T4.setItem(7, 1, __qtablewidgetitem183)
        __qtablewidgetitem184 = QTableWidgetItem()
        self.T4.setItem(8, 0, __qtablewidgetitem184)
        __qtablewidgetitem185 = QTableWidgetItem()
        self.T4.setItem(8, 1, __qtablewidgetitem185)
        __qtablewidgetitem186 = QTableWidgetItem()
        self.T4.setItem(9, 0, __qtablewidgetitem186)
        __qtablewidgetitem187 = QTableWidgetItem()
        self.T4.setItem(9, 1, __qtablewidgetitem187)
        __qtablewidgetitem188 = QTableWidgetItem()
        self.T4.setItem(10, 0, __qtablewidgetitem188)
        __qtablewidgetitem189 = QTableWidgetItem()
        self.T4.setItem(10, 1, __qtablewidgetitem189)
        __qtablewidgetitem190 = QTableWidgetItem()
        self.T4.setItem(11, 0, __qtablewidgetitem190)
        __qtablewidgetitem191 = QTableWidgetItem()
        self.T4.setItem(11, 1, __qtablewidgetitem191)
        __qtablewidgetitem192 = QTableWidgetItem()
        self.T4.setItem(12, 0, __qtablewidgetitem192)
        __qtablewidgetitem193 = QTableWidgetItem()
        self.T4.setItem(12, 1, __qtablewidgetitem193)
        __qtablewidgetitem194 = QTableWidgetItem()
        self.T4.setItem(13, 0, __qtablewidgetitem194)
        __qtablewidgetitem195 = QTableWidgetItem()
        self.T4.setItem(13, 1, __qtablewidgetitem195)
        __qtablewidgetitem196 = QTableWidgetItem()
        self.T4.setItem(14, 0, __qtablewidgetitem196)
        __qtablewidgetitem197 = QTableWidgetItem()
        self.T4.setItem(14, 1, __qtablewidgetitem197)
        __qtablewidgetitem198 = QTableWidgetItem()
        self.T4.setItem(15, 0, __qtablewidgetitem198)
        __qtablewidgetitem199 = QTableWidgetItem()
        self.T4.setItem(15, 1, __qtablewidgetitem199)
        __qtablewidgetitem200 = QTableWidgetItem()
        self.T4.setItem(16, 0, __qtablewidgetitem200)
        __qtablewidgetitem201 = QTableWidgetItem()
        self.T4.setItem(16, 1, __qtablewidgetitem201)
        __qtablewidgetitem202 = QTableWidgetItem()
        self.T4.setItem(17, 0, __qtablewidgetitem202)
        __qtablewidgetitem203 = QTableWidgetItem()
        self.T4.setItem(17, 1, __qtablewidgetitem203)
        __qtablewidgetitem204 = QTableWidgetItem()
        self.T4.setItem(18, 0, __qtablewidgetitem204)
        __qtablewidgetitem205 = QTableWidgetItem()
        self.T4.setItem(18, 1, __qtablewidgetitem205)
        __qtablewidgetitem206 = QTableWidgetItem()
        self.T4.setItem(19, 0, __qtablewidgetitem206)
        __qtablewidgetitem207 = QTableWidgetItem()
        self.T4.setItem(19, 1, __qtablewidgetitem207)
        __qtablewidgetitem208 = QTableWidgetItem()
        self.T4.setItem(20, 0, __qtablewidgetitem208)
        __qtablewidgetitem209 = QTableWidgetItem()
        self.T4.setItem(20, 1, __qtablewidgetitem209)
        __qtablewidgetitem210 = QTableWidgetItem()
        self.T4.setItem(21, 0, __qtablewidgetitem210)
        __qtablewidgetitem211 = QTableWidgetItem()
        self.T4.setItem(21, 1, __qtablewidgetitem211)
        __qtablewidgetitem212 = QTableWidgetItem()
        self.T4.setItem(22, 0, __qtablewidgetitem212)
        __qtablewidgetitem213 = QTableWidgetItem()
        self.T4.setItem(22, 1, __qtablewidgetitem213)
        __qtablewidgetitem214 = QTableWidgetItem()
        self.T4.setItem(23, 0, __qtablewidgetitem214)
        __qtablewidgetitem215 = QTableWidgetItem()
        self.T4.setItem(23, 1, __qtablewidgetitem215)
        self.T4.setObjectName(u"T4")
        self.T4.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.T4.setRowCount(24)
        self.T4.horizontalHeader().setDefaultSectionSize(167)

        self.verticalLayout_69.addWidget(self.T4)

        self.SlipSystemInfo.addWidget(self.BCC_123_111)

        self.verticalLayout_42.addWidget(self.SlipSystemInfo)

        self.PlasticProperties.addTab(self.Details, "")
        self.VoceParameters = QWidget()
        self.VoceParameters.setObjectName(u"VoceParameters")
        self.verticalLayout_67 = QVBoxLayout(self.VoceParameters)
#ifndef Q_OS_MAC
        self.verticalLayout_67.setSpacing(-1)
#endif
        self.verticalLayout_67.setObjectName(u"verticalLayout_67")
        self.verticalLayout_67.setContentsMargins(0, 0, 0, 0)
        self.TSlipSystemParameters = QTableWidget(self.VoceParameters)
        if (self.TSlipSystemParameters.columnCount() < 10):
            self.TSlipSystemParameters.setColumnCount(10)
        __qtablewidgetitem216 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(0, __qtablewidgetitem216)
        __qtablewidgetitem217 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(1, __qtablewidgetitem217)
        __qtablewidgetitem218 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(2, __qtablewidgetitem218)
        __qtablewidgetitem219 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(3, __qtablewidgetitem219)
        __qtablewidgetitem220 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(4, __qtablewidgetitem220)
        __qtablewidgetitem221 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(5, __qtablewidgetitem221)
        __qtablewidgetitem222 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(6, __qtablewidgetitem222)
        __qtablewidgetitem223 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(7, __qtablewidgetitem223)
        __qtablewidgetitem224 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(8, __qtablewidgetitem224)
        __qtablewidgetitem225 = QTableWidgetItem()
        self.TSlipSystemParameters.setHorizontalHeaderItem(9, __qtablewidgetitem225)
        self.TSlipSystemParameters.setObjectName(u"TSlipSystemParameters")
        sizePolicy14.setHeightForWidth(self.TSlipSystemParameters.sizePolicy().hasHeightForWidth())
        self.TSlipSystemParameters.setSizePolicy(sizePolicy14)
        self.TSlipSystemParameters.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustToContentsOnFirstShow)

        self.verticalLayout_67.addWidget(self.TSlipSystemParameters)

        self.frame_57 = QFrame(self.VoceParameters)
        self.frame_57.setObjectName(u"frame_57")
        self.frame_57.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_57.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_64 = QVBoxLayout(self.frame_57)
        self.verticalLayout_64.setSpacing(0)
        self.verticalLayout_64.setObjectName(u"verticalLayout_64")
        self.verticalLayout_64.setContentsMargins(0, 0, 0, 0)
        self.frame_62 = QFrame(self.frame_57)
        self.frame_62.setObjectName(u"frame_62")
        self.frame_62.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_62.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_32 = QGridLayout(self.frame_62)
        self.gridLayout_32.setObjectName(u"gridLayout_32")
        self.gridLayout_32.setContentsMargins(0, 0, 0, 0)
        self.LSelectSlipSystem = QLabel(self.frame_62)
        self.LSelectSlipSystem.setObjectName(u"LSelectSlipSystem")

        self.gridLayout_32.addWidget(self.LSelectSlipSystem, 0, 0, 1, 1)

        self.INSelectSlipSystem = QComboBox(self.frame_62)
        self.INSelectSlipSystem.setObjectName(u"INSelectSlipSystem")

        self.gridLayout_32.addWidget(self.INSelectSlipSystem, 0, 1, 1, 1)


        self.verticalLayout_64.addWidget(self.frame_62, 0, Qt.AlignmentFlag.AlignHCenter)


        self.verticalLayout_67.addWidget(self.frame_57)

        self.frame_55 = QFrame(self.VoceParameters)
        self.frame_55.setObjectName(u"frame_55")
        self.frame_55.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_55.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_29 = QGridLayout(self.frame_55)
        self.gridLayout_29.setObjectName(u"gridLayout_29")
        self.gridLayout_29.setContentsMargins(12, 0, 12, 10)
        self.pushButton_9 = QPushButton(self.frame_55)
        self.pushButton_9.setObjectName(u"pushButton_9")
        self.pushButton_9.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton_9.sizePolicy().hasHeightForWidth())
        self.pushButton_9.setSizePolicy(sizePolicy9)
        self.pushButton_9.setMaximumSize(QSize(17, 17))
        self.pushButton_9.setFont(font6)
        self.pushButton_9.setIcon(icon10)
        self.pushButton_9.setIconSize(QSize(20, 20))
        self.pushButton_9.setAutoDefault(False)
        self.pushButton_9.setFlat(True)

        self.gridLayout_29.addWidget(self.pushButton_9, 0, 2, 1, 1)

        self.Lgamd0x = QLabel(self.frame_55)
        self.Lgamd0x.setObjectName(u"Lgamd0x")

        self.gridLayout_29.addWidget(self.Lgamd0x, 0, 3, 1, 1)

        self.Lhlatex = QLabel(self.frame_55)
        self.Lhlatex.setObjectName(u"Lhlatex")

        self.gridLayout_29.addWidget(self.Lhlatex, 10, 3, 1, 1)

        self.INtau0xb = QLineEdit(self.frame_55)
        self.INtau0xb.setObjectName(u"INtau0xb")

        self.gridLayout_29.addWidget(self.INtau0xb, 3, 4, 1, 1)

        self.INgamd0x = QLineEdit(self.frame_55)
        self.INgamd0x.setObjectName(u"INgamd0x")

        self.gridLayout_29.addWidget(self.INgamd0x, 0, 4, 1, 1)

        self.pushButton_6 = QPushButton(self.frame_55)
        self.pushButton_6.setObjectName(u"pushButton_6")
        self.pushButton_6.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton_6.sizePolicy().hasHeightForWidth())
        self.pushButton_6.setSizePolicy(sizePolicy9)
        self.pushButton_6.setMaximumSize(QSize(17, 17))
        self.pushButton_6.setIcon(icon10)
        self.pushButton_6.setIconSize(QSize(20, 20))
        self.pushButton_6.setAutoDefault(False)
        self.pushButton_6.setFlat(True)

        self.gridLayout_29.addWidget(self.pushButton_6, 10, 2, 1, 1)

        self.Ltau0xb = QLabel(self.frame_55)
        self.Ltau0xb.setObjectName(u"Ltau0xb")

        self.gridLayout_29.addWidget(self.Ltau0xb, 3, 3, 1, 1)

        self.INthet1 = QLineEdit(self.frame_55)
        self.INthet1.setObjectName(u"INthet1")

        self.gridLayout_29.addWidget(self.INthet1, 5, 4, 1, 1)

        self.Ltau0xf = QLabel(self.frame_55)
        self.Ltau0xf.setObjectName(u"Ltau0xf")

        self.gridLayout_29.addWidget(self.Ltau0xf, 3, 0, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.Lthet0 = QLabel(self.frame_55)
        self.Lthet0.setObjectName(u"Lthet0")

        self.gridLayout_29.addWidget(self.Lthet0, 5, 0, 1, 1)

        self.Lhselfx = QLabel(self.frame_55)
        self.Lhselfx.setObjectName(u"Lhselfx")

        self.gridLayout_29.addWidget(self.Lhselfx, 10, 0, 1, 1)

        self.Ltau1x = QLabel(self.frame_55)
        self.Ltau1x.setObjectName(u"Ltau1x")

        self.gridLayout_29.addWidget(self.Ltau1x, 4, 0, 1, 1)

        self.pushButton_5 = QPushButton(self.frame_55)
        self.pushButton_5.setObjectName(u"pushButton_5")
        self.pushButton_5.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton_5.sizePolicy().hasHeightForWidth())
        self.pushButton_5.setSizePolicy(sizePolicy9)
        self.pushButton_5.setMaximumSize(QSize(17, 17))
        self.pushButton_5.setIcon(icon10)
        self.pushButton_5.setIconSize(QSize(20, 20))
        self.pushButton_5.setAutoDefault(False)
        self.pushButton_5.setFlat(True)

        self.gridLayout_29.addWidget(self.pushButton_5, 5, 5, 1, 1)

        self.INtau1x = QLineEdit(self.frame_55)
        self.INtau1x.setObjectName(u"INtau1x")

        self.gridLayout_29.addWidget(self.INtau1x, 4, 1, 1, 1)

        self.pushButton_7 = QPushButton(self.frame_55)
        self.pushButton_7.setObjectName(u"pushButton_7")
        self.pushButton_7.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton_7.sizePolicy().hasHeightForWidth())
        self.pushButton_7.setSizePolicy(sizePolicy9)
        self.pushButton_7.setMaximumSize(QSize(17, 17))
        self.pushButton_7.setIcon(icon10)
        self.pushButton_7.setIconSize(QSize(20, 20))
        self.pushButton_7.setAutoDefault(False)
        self.pushButton_7.setFlat(True)

        self.gridLayout_29.addWidget(self.pushButton_7, 10, 5, 1, 1)

        self.pushButton_8 = QPushButton(self.frame_55)
        self.pushButton_8.setObjectName(u"pushButton_8")
        self.pushButton_8.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton_8.sizePolicy().hasHeightForWidth())
        self.pushButton_8.setSizePolicy(sizePolicy9)
        self.pushButton_8.setMaximumSize(QSize(17, 17))
        self.pushButton_8.setFont(font6)
        self.pushButton_8.setIcon(icon10)
        self.pushButton_8.setIconSize(QSize(20, 20))
        self.pushButton_8.setAutoDefault(False)
        self.pushButton_8.setFlat(True)

        self.gridLayout_29.addWidget(self.pushButton_8, 0, 5, 1, 1)

        self.INhselfx = QLineEdit(self.frame_55)
        self.INhselfx.setObjectName(u"INhselfx")

        self.gridLayout_29.addWidget(self.INhselfx, 10, 1, 1, 1)

        self.INnrsx = QLineEdit(self.frame_55)
        self.INnrsx.setObjectName(u"INnrsx")

        self.gridLayout_29.addWidget(self.INnrsx, 0, 1, 1, 1)

        self.pushButton_2 = QPushButton(self.frame_55)
        self.pushButton_2.setObjectName(u"pushButton_2")
        self.pushButton_2.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton_2.sizePolicy().hasHeightForWidth())
        self.pushButton_2.setSizePolicy(sizePolicy9)
        self.pushButton_2.setMaximumSize(QSize(17, 17))
        self.pushButton_2.setIcon(icon10)
        self.pushButton_2.setIconSize(QSize(20, 20))
        self.pushButton_2.setAutoDefault(False)
        self.pushButton_2.setFlat(True)

        self.gridLayout_29.addWidget(self.pushButton_2, 3, 5, 1, 1)

        self.Lthet1 = QLabel(self.frame_55)
        self.Lthet1.setObjectName(u"Lthet1")

        self.gridLayout_29.addWidget(self.Lthet1, 5, 3, 1, 1)

        self.INthet0 = QLineEdit(self.frame_55)
        self.INthet0.setObjectName(u"INthet0")

        self.gridLayout_29.addWidget(self.INthet0, 5, 1, 1, 1)

        self.INhlatex = QLineEdit(self.frame_55)
        self.INhlatex.setObjectName(u"INhlatex")

        self.gridLayout_29.addWidget(self.INhlatex, 10, 4, 1, 1)

        self.pushButton = QPushButton(self.frame_55)
        self.pushButton.setObjectName(u"pushButton")
        self.pushButton.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton.sizePolicy().hasHeightForWidth())
        self.pushButton.setSizePolicy(sizePolicy9)
        self.pushButton.setMaximumSize(QSize(17, 17))
        self.pushButton.setFont(font6)
        self.pushButton.setIcon(icon10)
        self.pushButton.setIconSize(QSize(20, 20))
        self.pushButton.setAutoDefault(False)
        self.pushButton.setFlat(True)

        self.gridLayout_29.addWidget(self.pushButton, 3, 2, 1, 1)

        self.pushButton_3 = QPushButton(self.frame_55)
        self.pushButton_3.setObjectName(u"pushButton_3")
        self.pushButton_3.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton_3.sizePolicy().hasHeightForWidth())
        self.pushButton_3.setSizePolicy(sizePolicy9)
        self.pushButton_3.setMaximumSize(QSize(17, 17))
        self.pushButton_3.setIcon(icon10)
        self.pushButton_3.setIconSize(QSize(20, 20))
        self.pushButton_3.setAutoDefault(False)
        self.pushButton_3.setFlat(True)

        self.gridLayout_29.addWidget(self.pushButton_3, 4, 2, 1, 1)

        self.INtau0xf = QLineEdit(self.frame_55)
        self.INtau0xf.setObjectName(u"INtau0xf")

        self.gridLayout_29.addWidget(self.INtau0xf, 3, 1, 1, 1)

        self.pushButton_4 = QPushButton(self.frame_55)
        self.pushButton_4.setObjectName(u"pushButton_4")
        self.pushButton_4.setEnabled(False)
        sizePolicy9.setHeightForWidth(self.pushButton_4.sizePolicy().hasHeightForWidth())
        self.pushButton_4.setSizePolicy(sizePolicy9)
        self.pushButton_4.setMaximumSize(QSize(17, 17))
        self.pushButton_4.setIcon(icon10)
        self.pushButton_4.setIconSize(QSize(20, 20))
        self.pushButton_4.setAutoDefault(False)
        self.pushButton_4.setFlat(True)

        self.gridLayout_29.addWidget(self.pushButton_4, 5, 2, 1, 1)

        self.Lnrsx = QLabel(self.frame_55)
        self.Lnrsx.setObjectName(u"Lnrsx")

        self.gridLayout_29.addWidget(self.Lnrsx, 0, 0, 1, 1)


        self.verticalLayout_67.addWidget(self.frame_55)

        self.PlasticProperties.addTab(self.VoceParameters, "")

        self.verticalLayout_70.addWidget(self.PlasticProperties)

        self.EnablePlasticity.addWidget(self.page_26)

        self.verticalLayout_26.addWidget(self.EnablePlasticity)

        icon12 = QIcon()
        icon12.addFile(u":/Blue Icons/Blue Icons/Plastic.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.MaterialMenu_2.addTab(self.Plastic, icon12, "")

        self.verticalLayout_30.addWidget(self.MaterialMenu_2)

        self.frame_61 = QFrame(self.DefineMaterialsBulkForm)
        self.frame_61.setObjectName(u"frame_61")
        self.frame_61.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_61.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_9 = QHBoxLayout(self.frame_61)
        self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
        self.horizontalLayout_9.setContentsMargins(0, 0, 0, 0)
        self.BAddMaterial_2 = QPushButton(self.frame_61)
        self.BAddMaterial_2.setObjectName(u"BAddMaterial_2")

        self.horizontalLayout_9.addWidget(self.BAddMaterial_2)

        self.BDeleteMaterial_2 = QPushButton(self.frame_61)
        self.BDeleteMaterial_2.setObjectName(u"BDeleteMaterial_2")

        self.horizontalLayout_9.addWidget(self.BDeleteMaterial_2)


        self.verticalLayout_30.addWidget(self.frame_61)

        self.TMaterials_2 = QTableWidget(self.DefineMaterialsBulkForm)
        if (self.TMaterials_2.columnCount() < 36):
            self.TMaterials_2.setColumnCount(36)
        __qtablewidgetitem226 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(0, __qtablewidgetitem226)
        __qtablewidgetitem227 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(1, __qtablewidgetitem227)
        __qtablewidgetitem228 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(2, __qtablewidgetitem228)
        __qtablewidgetitem229 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(3, __qtablewidgetitem229)
        __qtablewidgetitem230 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(4, __qtablewidgetitem230)
        __qtablewidgetitem231 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(5, __qtablewidgetitem231)
        __qtablewidgetitem232 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(6, __qtablewidgetitem232)
        __qtablewidgetitem233 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(7, __qtablewidgetitem233)
        __qtablewidgetitem234 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(8, __qtablewidgetitem234)
        __qtablewidgetitem235 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(9, __qtablewidgetitem235)
        __qtablewidgetitem236 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(10, __qtablewidgetitem236)
        __qtablewidgetitem237 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(11, __qtablewidgetitem237)
        __qtablewidgetitem238 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(12, __qtablewidgetitem238)
        __qtablewidgetitem239 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(13, __qtablewidgetitem239)
        __qtablewidgetitem240 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(14, __qtablewidgetitem240)
        __qtablewidgetitem241 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(15, __qtablewidgetitem241)
        __qtablewidgetitem242 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(16, __qtablewidgetitem242)
        __qtablewidgetitem243 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(17, __qtablewidgetitem243)
        __qtablewidgetitem244 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(18, __qtablewidgetitem244)
        __qtablewidgetitem245 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(19, __qtablewidgetitem245)
        __qtablewidgetitem246 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(20, __qtablewidgetitem246)
        __qtablewidgetitem247 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(21, __qtablewidgetitem247)
        __qtablewidgetitem248 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(22, __qtablewidgetitem248)
        __qtablewidgetitem249 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(23, __qtablewidgetitem249)
        __qtablewidgetitem250 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(24, __qtablewidgetitem250)
        __qtablewidgetitem251 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(25, __qtablewidgetitem251)
        __qtablewidgetitem252 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(26, __qtablewidgetitem252)
        __qtablewidgetitem253 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(27, __qtablewidgetitem253)
        __qtablewidgetitem254 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(28, __qtablewidgetitem254)
        __qtablewidgetitem255 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(29, __qtablewidgetitem255)
        __qtablewidgetitem256 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(30, __qtablewidgetitem256)
        __qtablewidgetitem257 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(31, __qtablewidgetitem257)
        __qtablewidgetitem258 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(32, __qtablewidgetitem258)
        __qtablewidgetitem259 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(33, __qtablewidgetitem259)
        __qtablewidgetitem260 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(34, __qtablewidgetitem260)
        __qtablewidgetitem261 = QTableWidgetItem()
        self.TMaterials_2.setHorizontalHeaderItem(35, __qtablewidgetitem261)
        self.TMaterials_2.setObjectName(u"TMaterials_2")
        self.TMaterials_2.setEnabled(True)
        self.TMaterials_2.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustToContents)
        self.TMaterials_2.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TMaterials_2.setRowCount(0)

        self.verticalLayout_30.addWidget(self.TMaterials_2)

        self.frame_60 = QFrame(self.DefineMaterialsBulkForm)
        self.frame_60.setObjectName(u"frame_60")
        self.frame_60.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_60.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_6 = QHBoxLayout(self.frame_60)
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.horizontalLayout_6.setContentsMargins(0, 0, 0, 0)

        self.verticalLayout_30.addWidget(self.frame_60)

        self.DefineMaterialsOptions.addWidget(self.DefineMaterialsBulkForm)

        self.verticalLayout_49.addWidget(self.DefineMaterialsOptions)

        self.DefineAssignMats.addTab(self.DefineMaterials, icon5, "")
        self.AssignMaterials = QWidget()
        self.AssignMaterials.setObjectName(u"AssignMaterials")
        self.verticalLayout_56 = QVBoxLayout(self.AssignMaterials)
        self.verticalLayout_56.setObjectName(u"verticalLayout_56")
        self.verticalLayout_56.setContentsMargins(0, 12, 0, 0)
        self.INSelectAssignMaterials = QComboBox(self.AssignMaterials)
        self.INSelectAssignMaterials.addItem("")
        self.INSelectAssignMaterials.addItem("")
        self.INSelectAssignMaterials.addItem("")
        self.INSelectAssignMaterials.setObjectName(u"INSelectAssignMaterials")

        self.verticalLayout_56.addWidget(self.INSelectAssignMaterials)

        self.AssignMaterialsOptions = QStackedWidget(self.AssignMaterials)
        self.AssignMaterialsOptions.setObjectName(u"AssignMaterialsOptions")
        self.page_16 = QWidget()
        self.page_16.setObjectName(u"page_16")
        self.verticalLayout_57 = QVBoxLayout(self.page_16)
        self.verticalLayout_57.setObjectName(u"verticalLayout_57")
        self.frame_38 = QFrame(self.page_16)
        self.frame_38.setObjectName(u"frame_38")
        self.frame_38.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_38.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_30 = QFormLayout(self.frame_38)
        self.formLayout_30.setObjectName(u"formLayout_30")
        self.formLayout_30.setContentsMargins(0, 0, 0, 0)
        self.LMaterialName = QLabel(self.frame_38)
        self.LMaterialName.setObjectName(u"LMaterialName")

        self.formLayout_30.setWidget(0, QFormLayout.LabelRole, self.LMaterialName)

        self.INRegion = QComboBox(self.frame_38)
        self.INRegion.setObjectName(u"INRegion")

        self.formLayout_30.setWidget(0, QFormLayout.FieldRole, self.INRegion)

        self.LRegion_5 = QLabel(self.frame_38)
        self.LRegion_5.setObjectName(u"LRegion_5")

        self.formLayout_30.setWidget(1, QFormLayout.LabelRole, self.LRegion_5)

        self.INMaterial = QComboBox(self.frame_38)
        self.INMaterial.setObjectName(u"INMaterial")

        self.formLayout_30.setWidget(1, QFormLayout.FieldRole, self.INMaterial)


        self.verticalLayout_57.addWidget(self.frame_38)

        self.frame_39 = QFrame(self.page_16)
        self.frame_39.setObjectName(u"frame_39")
        self.frame_39.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_39.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_2 = QHBoxLayout(self.frame_39)
        self.horizontalLayout_2.setSpacing(0)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.frame_40 = QFrame(self.frame_39)
        self.frame_40.setObjectName(u"frame_40")
        self.frame_40.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_40.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_58 = QVBoxLayout(self.frame_40)
        self.verticalLayout_58.setObjectName(u"verticalLayout_58")
        self.verticalLayout_58.setContentsMargins(10, 0, 10, 0)
        self.LDensity = QLabel(self.frame_40)
        self.LDensity.setObjectName(u"LDensity")

        self.verticalLayout_58.addWidget(self.LDensity, 0, Qt.AlignmentFlag.AlignTop)

        self.INDensity = QLineEdit(self.frame_40)
        self.INDensity.setObjectName(u"INDensity")
        self.INDensity.setEnabled(True)

        self.verticalLayout_58.addWidget(self.INDensity)

        self.LSIE = QLabel(self.frame_40)
        self.LSIE.setObjectName(u"LSIE")
        self.LSIE.setWordWrap(False)

        self.verticalLayout_58.addWidget(self.LSIE)

        self.INSIE = QLineEdit(self.frame_40)
        self.INSIE.setObjectName(u"INSIE")
        self.INSIE.setEnabled(True)

        self.verticalLayout_58.addWidget(self.INSIE)


        self.horizontalLayout_2.addWidget(self.frame_40)

        self.frame_41 = QFrame(self.frame_39)
        self.frame_41.setObjectName(u"frame_41")
        self.frame_41.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_41.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_59 = QVBoxLayout(self.frame_41)
        self.verticalLayout_59.setObjectName(u"verticalLayout_59")
        self.verticalLayout_59.setContentsMargins(10, 0, 0, 0)
        self.label_5 = QLabel(self.frame_41)
        self.label_5.setObjectName(u"label_5")

        self.verticalLayout_59.addWidget(self.label_5)

        self.frame_42 = QFrame(self.frame_41)
        self.frame_42.setObjectName(u"frame_42")
        self.frame_42.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_42.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_31 = QFormLayout(self.frame_42)
        self.formLayout_31.setObjectName(u"formLayout_31")
        self.formLayout_31.setContentsMargins(0, 0, 10, 0)
        self.LVely = QLabel(self.frame_42)
        self.LVely.setObjectName(u"LVely")

        self.formLayout_31.setWidget(1, QFormLayout.LabelRole, self.LVely)

        self.LVelx = QLabel(self.frame_42)
        self.LVelx.setObjectName(u"LVelx")

        self.formLayout_31.setWidget(0, QFormLayout.LabelRole, self.LVelx)

        self.INVelocityX = QLineEdit(self.frame_42)
        self.INVelocityX.setObjectName(u"INVelocityX")
        self.INVelocityX.setEnabled(True)

        self.formLayout_31.setWidget(0, QFormLayout.FieldRole, self.INVelocityX)

        self.INVelocityY = QLineEdit(self.frame_42)
        self.INVelocityY.setObjectName(u"INVelocityY")
        self.INVelocityY.setEnabled(True)

        self.formLayout_31.setWidget(1, QFormLayout.FieldRole, self.INVelocityY)

        self.LVelz = QLabel(self.frame_42)
        self.LVelz.setObjectName(u"LVelz")

        self.formLayout_31.setWidget(2, QFormLayout.LabelRole, self.LVelz)

        self.INVelocityZ = QLineEdit(self.frame_42)
        self.INVelocityZ.setObjectName(u"INVelocityZ")
        self.INVelocityZ.setEnabled(True)

        self.formLayout_31.setWidget(2, QFormLayout.FieldRole, self.INVelocityZ)


        self.verticalLayout_59.addWidget(self.frame_42)


        self.horizontalLayout_2.addWidget(self.frame_41)


        self.verticalLayout_57.addWidget(self.frame_39)

        self.Baddmaterialassignment = QPushButton(self.page_16)
        self.Baddmaterialassignment.setObjectName(u"Baddmaterialassignment")

        self.verticalLayout_57.addWidget(self.Baddmaterialassignment)

        self.Tassignmat = QTableWidget(self.page_16)
        if (self.Tassignmat.columnCount() < 7):
            self.Tassignmat.setColumnCount(7)
        __qtablewidgetitem262 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(0, __qtablewidgetitem262)
        __qtablewidgetitem263 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(1, __qtablewidgetitem263)
        __qtablewidgetitem264 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(2, __qtablewidgetitem264)
        __qtablewidgetitem265 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(3, __qtablewidgetitem265)
        __qtablewidgetitem266 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(4, __qtablewidgetitem266)
        __qtablewidgetitem267 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(5, __qtablewidgetitem267)
        __qtablewidgetitem268 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(6, __qtablewidgetitem268)
        self.Tassignmat.setObjectName(u"Tassignmat")
        sizePolicy7.setHeightForWidth(self.Tassignmat.sizePolicy().hasHeightForWidth())
        self.Tassignmat.setSizePolicy(sizePolicy7)

        self.verticalLayout_57.addWidget(self.Tassignmat)

        self.frame_43 = QFrame(self.page_16)
        self.frame_43.setObjectName(u"frame_43")
        self.frame_43.setFrameShape(QFrame.Shape.NoFrame)
        self.formLayout_32 = QFormLayout(self.frame_43)
        self.formLayout_32.setObjectName(u"formLayout_32")
        self.formLayout_32.setContentsMargins(0, 0, 0, 0)
        self.BUpMaterial = QToolButton(self.frame_43)
        self.BUpMaterial.setObjectName(u"BUpMaterial")
        self.BUpMaterial.setIconSize(QSize(32, 32))
        self.BUpMaterial.setAutoRaise(False)
        self.BUpMaterial.setArrowType(Qt.ArrowType.UpArrow)

        self.formLayout_32.setWidget(0, QFormLayout.LabelRole, self.BUpMaterial)

        self.BDownMaterial = QToolButton(self.frame_43)
        self.BDownMaterial.setObjectName(u"BDownMaterial")
        self.BDownMaterial.setIconSize(QSize(32, 32))
        self.BDownMaterial.setArrowType(Qt.ArrowType.DownArrow)

        self.formLayout_32.setWidget(0, QFormLayout.FieldRole, self.BDownMaterial)


        self.verticalLayout_57.addWidget(self.frame_43)

        self.Bdeletematerialassignment = QPushButton(self.page_16)
        self.Bdeletematerialassignment.setObjectName(u"Bdeletematerialassignment")

        self.verticalLayout_57.addWidget(self.Bdeletematerialassignment)

        self.verticalSpacer_24 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_57.addItem(self.verticalSpacer_24)

        self.AssignMaterialsOptions.addWidget(self.page_16)
        self.page_17 = QWidget()
        self.page_17.setObjectName(u"page_17")
        self.verticalLayout_4 = QVBoxLayout(self.page_17)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.frame_48 = QFrame(self.page_17)
        self.frame_48.setObjectName(u"frame_48")
        self.frame_48.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_48.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_33 = QFormLayout(self.frame_48)
        self.formLayout_33.setObjectName(u"formLayout_33")
        self.formLayout_33.setContentsMargins(0, 0, 0, 0)
        self.LMaterialName_4 = QLabel(self.frame_48)
        self.LMaterialName_4.setObjectName(u"LMaterialName_4")

        self.formLayout_33.setWidget(0, QFormLayout.LabelRole, self.LMaterialName_4)

        self.INRegion_2 = QComboBox(self.frame_48)
        self.INRegion_2.setObjectName(u"INRegion_2")

        self.formLayout_33.setWidget(0, QFormLayout.FieldRole, self.INRegion_2)

        self.LRegion_6 = QLabel(self.frame_48)
        self.LRegion_6.setObjectName(u"LRegion_6")

        self.formLayout_33.setWidget(1, QFormLayout.LabelRole, self.LRegion_6)

        self.INMaterial_2 = QComboBox(self.frame_48)
        self.INMaterial_2.setObjectName(u"INMaterial_2")

        self.formLayout_33.setWidget(1, QFormLayout.FieldRole, self.INMaterial_2)


        self.verticalLayout_4.addWidget(self.frame_48)

        self.BAddMaterialAssignment = QPushButton(self.page_17)
        self.BAddMaterialAssignment.setObjectName(u"BAddMaterialAssignment")

        self.verticalLayout_4.addWidget(self.BAddMaterialAssignment)

        self.TMaterialAssignment = QTableWidget(self.page_17)
        if (self.TMaterialAssignment.columnCount() < 2):
            self.TMaterialAssignment.setColumnCount(2)
        __qtablewidgetitem269 = QTableWidgetItem()
        self.TMaterialAssignment.setHorizontalHeaderItem(0, __qtablewidgetitem269)
        __qtablewidgetitem270 = QTableWidgetItem()
        self.TMaterialAssignment.setHorizontalHeaderItem(1, __qtablewidgetitem270)
        self.TMaterialAssignment.setObjectName(u"TMaterialAssignment")
        sizePolicy14.setHeightForWidth(self.TMaterialAssignment.sizePolicy().hasHeightForWidth())
        self.TMaterialAssignment.setSizePolicy(sizePolicy14)
        self.TMaterialAssignment.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustIgnored)
        self.TMaterialAssignment.setEditTriggers(QAbstractItemView.EditTrigger.AnyKeyPressed|QAbstractItemView.EditTrigger.EditKeyPressed)
        self.TMaterialAssignment.horizontalHeader().setCascadingSectionResizes(False)
        self.TMaterialAssignment.horizontalHeader().setDefaultSectionSize(168)
        self.TMaterialAssignment.verticalHeader().setCascadingSectionResizes(False)

        self.verticalLayout_4.addWidget(self.TMaterialAssignment)

        self.BDeleteMaterialAssignment = QPushButton(self.page_17)
        self.BDeleteMaterialAssignment.setObjectName(u"BDeleteMaterialAssignment")

        self.verticalLayout_4.addWidget(self.BDeleteMaterialAssignment)

        self.verticalSpacer_8 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_4.addItem(self.verticalSpacer_8)

        self.AssignMaterialsOptions.addWidget(self.page_17)
        self.page_15 = QWidget()
        self.page_15.setObjectName(u"page_15")
        self.verticalLayout_24 = QVBoxLayout(self.page_15)
        self.verticalLayout_24.setObjectName(u"verticalLayout_24")
        self.frame_49 = QFrame(self.page_15)
        self.frame_49.setObjectName(u"frame_49")
        self.frame_49.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_49.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_34 = QFormLayout(self.frame_49)
        self.formLayout_34.setObjectName(u"formLayout_34")
        self.formLayout_34.setContentsMargins(0, 0, 0, 0)
        self.LMaterialName_5 = QLabel(self.frame_49)
        self.LMaterialName_5.setObjectName(u"LMaterialName_5")

        self.formLayout_34.setWidget(0, QFormLayout.LabelRole, self.LMaterialName_5)

        self.INRegion_3 = QComboBox(self.frame_49)
        self.INRegion_3.setObjectName(u"INRegion_3")

        self.formLayout_34.setWidget(0, QFormLayout.FieldRole, self.INRegion_3)

        self.LRegion_7 = QLabel(self.frame_49)
        self.LRegion_7.setObjectName(u"LRegion_7")

        self.formLayout_34.setWidget(1, QFormLayout.LabelRole, self.LRegion_7)

        self.INMaterial_3 = QComboBox(self.frame_49)
        self.INMaterial_3.setObjectName(u"INMaterial_3")

        self.formLayout_34.setWidget(1, QFormLayout.FieldRole, self.INMaterial_3)


        self.verticalLayout_24.addWidget(self.frame_49)

        self.BAddMaterialAssignment_2 = QPushButton(self.page_15)
        self.BAddMaterialAssignment_2.setObjectName(u"BAddMaterialAssignment_2")

        self.verticalLayout_24.addWidget(self.BAddMaterialAssignment_2)

        self.TMaterialAssignment_2 = QTableWidget(self.page_15)
        if (self.TMaterialAssignment_2.columnCount() < 2):
            self.TMaterialAssignment_2.setColumnCount(2)
        __qtablewidgetitem271 = QTableWidgetItem()
        self.TMaterialAssignment_2.setHorizontalHeaderItem(0, __qtablewidgetitem271)
        __qtablewidgetitem272 = QTableWidgetItem()
        self.TMaterialAssignment_2.setHorizontalHeaderItem(1, __qtablewidgetitem272)
        self.TMaterialAssignment_2.setObjectName(u"TMaterialAssignment_2")
        sizePolicy14.setHeightForWidth(self.TMaterialAssignment_2.sizePolicy().hasHeightForWidth())
        self.TMaterialAssignment_2.setSizePolicy(sizePolicy14)
        self.TMaterialAssignment_2.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustIgnored)
        self.TMaterialAssignment_2.setEditTriggers(QAbstractItemView.EditTrigger.AnyKeyPressed|QAbstractItemView.EditTrigger.EditKeyPressed)
        self.TMaterialAssignment_2.horizontalHeader().setCascadingSectionResizes(False)
        self.TMaterialAssignment_2.horizontalHeader().setDefaultSectionSize(168)
        self.TMaterialAssignment_2.verticalHeader().setCascadingSectionResizes(False)

        self.verticalLayout_24.addWidget(self.TMaterialAssignment_2)

        self.BDeleteMaterialAssignment_2 = QPushButton(self.page_15)
        self.BDeleteMaterialAssignment_2.setObjectName(u"BDeleteMaterialAssignment_2")

        self.verticalLayout_24.addWidget(self.BDeleteMaterialAssignment_2)

        self.verticalSpacer_27 = QSpacerItem(20, 379, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_24.addItem(self.verticalSpacer_27)

        self.AssignMaterialsOptions.addWidget(self.page_15)

        self.verticalLayout_56.addWidget(self.AssignMaterialsOptions)

        icon13 = QIcon()
        icon13.addFile(u":/Blue Icons/Blue Icons/Clipboard.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.DefineAssignMats.addTab(self.AssignMaterials, icon13, "")

        self.verticalLayout_5.addWidget(self.DefineAssignMats)

        self.verticalSpacer_23 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_5.addItem(self.verticalSpacer_23)

        self.ToolWindow.addWidget(self.MaterialsTab)
        self.BoundaryConditionsTab = QWidget()
        self.BoundaryConditionsTab.setObjectName(u"BoundaryConditionsTab")
        self.verticalLayout_7 = QVBoxLayout(self.BoundaryConditionsTab)
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.LBoundaryConditions = QLabel(self.BoundaryConditionsTab)
        self.LBoundaryConditions.setObjectName(u"LBoundaryConditions")
        sizePolicy8.setHeightForWidth(self.LBoundaryConditions.sizePolicy().hasHeightForWidth())
        self.LBoundaryConditions.setSizePolicy(sizePolicy8)

        self.verticalLayout_7.addWidget(self.LBoundaryConditions)

        self.INSelectBoundaryConditions = QComboBox(self.BoundaryConditionsTab)
        self.INSelectBoundaryConditions.addItem("")
        self.INSelectBoundaryConditions.addItem("")
        self.INSelectBoundaryConditions.addItem("")
        self.INSelectBoundaryConditions.setObjectName(u"INSelectBoundaryConditions")

        self.verticalLayout_7.addWidget(self.INSelectBoundaryConditions)

        self.BoundaryConditionsOptions = QStackedWidget(self.BoundaryConditionsTab)
        self.BoundaryConditionsOptions.setObjectName(u"BoundaryConditionsOptions")
        self.page_18 = QWidget()
        self.page_18.setObjectName(u"page_18")
        self.verticalLayout_31 = QVBoxLayout(self.page_18)
        self.verticalLayout_31.setObjectName(u"verticalLayout_31")
        self.bcsettings = QFrame(self.page_18)
        self.bcsettings.setObjectName(u"bcsettings")
        self.bcsettings.setFrameShape(QFrame.Shape.NoFrame)
        self.formLayout_18 = QFormLayout(self.bcsettings)
        self.formLayout_18.setObjectName(u"formLayout_18")
        self.formLayout_18.setContentsMargins(0, 0, 0, 0)
        self.Lbndry = QLabel(self.bcsettings)
        self.Lbndry.setObjectName(u"Lbndry")

        self.formLayout_18.setWidget(0, QFormLayout.LabelRole, self.Lbndry)

        self.INBoundary = QComboBox(self.bcsettings)
        self.INBoundary.addItem("")
        self.INBoundary.addItem("")
        self.INBoundary.addItem("")
        self.INBoundary.setObjectName(u"INBoundary")

        self.formLayout_18.setWidget(0, QFormLayout.FieldRole, self.INBoundary)

        self.Lvalue = QLabel(self.bcsettings)
        self.Lvalue.setObjectName(u"Lvalue")

        self.formLayout_18.setWidget(1, QFormLayout.LabelRole, self.Lvalue)

        self.INPlanePosition = QLineEdit(self.bcsettings)
        self.INPlanePosition.setObjectName(u"INPlanePosition")

        self.formLayout_18.setWidget(1, QFormLayout.FieldRole, self.INPlanePosition)

        self.Ltype = QLabel(self.bcsettings)
        self.Ltype.setObjectName(u"Ltype")

        self.formLayout_18.setWidget(2, QFormLayout.LabelRole, self.Ltype)

        self.INType = QComboBox(self.bcsettings)
        self.INType.addItem("")
        self.INType.addItem("")
        self.INType.addItem("")
        self.INType.setObjectName(u"INType")

        self.formLayout_18.setWidget(2, QFormLayout.FieldRole, self.INType)

        self.INVel0 = QLineEdit(self.bcsettings)
        self.INVel0.setObjectName(u"INVel0")
        self.INVel0.setEnabled(False)

        self.formLayout_18.setWidget(3, QFormLayout.FieldRole, self.INVel0)

        self.INVel1 = QLineEdit(self.bcsettings)
        self.INVel1.setObjectName(u"INVel1")
        self.INVel1.setEnabled(False)

        self.formLayout_18.setWidget(4, QFormLayout.FieldRole, self.INVel1)

        self.INVelstart = QLineEdit(self.bcsettings)
        self.INVelstart.setObjectName(u"INVelstart")
        self.INVelstart.setEnabled(False)

        self.formLayout_18.setWidget(5, QFormLayout.FieldRole, self.INVelstart)

        self.INVelend = QLineEdit(self.bcsettings)
        self.INVelend.setObjectName(u"INVelend")
        self.INVelend.setEnabled(False)

        self.formLayout_18.setWidget(6, QFormLayout.FieldRole, self.INVelend)

        self.LVel0 = QLabel(self.bcsettings)
        self.LVel0.setObjectName(u"LVel0")
        self.LVel0.setEnabled(False)

        self.formLayout_18.setWidget(3, QFormLayout.LabelRole, self.LVel0)

        self.LVel1 = QLabel(self.bcsettings)
        self.LVel1.setObjectName(u"LVel1")
        self.LVel1.setEnabled(False)

        self.formLayout_18.setWidget(4, QFormLayout.LabelRole, self.LVel1)

        self.LVelstart = QLabel(self.bcsettings)
        self.LVelstart.setObjectName(u"LVelstart")
        self.LVelstart.setEnabled(False)

        self.formLayout_18.setWidget(5, QFormLayout.LabelRole, self.LVelstart)

        self.LVelend = QLabel(self.bcsettings)
        self.LVelend.setObjectName(u"LVelend")
        self.LVelend.setEnabled(False)

        self.formLayout_18.setWidget(6, QFormLayout.LabelRole, self.LVelend)


        self.verticalLayout_31.addWidget(self.bcsettings)

        self.BaddBC = QPushButton(self.page_18)
        self.BaddBC.setObjectName(u"BaddBC")

        self.verticalLayout_31.addWidget(self.BaddBC)

        self.TBoundaryConditions = QTableWidget(self.page_18)
        if (self.TBoundaryConditions.columnCount() < 7):
            self.TBoundaryConditions.setColumnCount(7)
        __qtablewidgetitem273 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(0, __qtablewidgetitem273)
        __qtablewidgetitem274 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(1, __qtablewidgetitem274)
        __qtablewidgetitem275 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(2, __qtablewidgetitem275)
        __qtablewidgetitem276 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(3, __qtablewidgetitem276)
        __qtablewidgetitem277 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(4, __qtablewidgetitem277)
        __qtablewidgetitem278 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(5, __qtablewidgetitem278)
        __qtablewidgetitem279 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(6, __qtablewidgetitem279)
        self.TBoundaryConditions.setObjectName(u"TBoundaryConditions")

        self.verticalLayout_31.addWidget(self.TBoundaryConditions)

        self.BdeleteBC = QPushButton(self.page_18)
        self.BdeleteBC.setObjectName(u"BdeleteBC")

        self.verticalLayout_31.addWidget(self.BdeleteBC)

        self.verticalSpacer_10 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_31.addItem(self.verticalSpacer_10)

        self.BoundaryConditionsOptions.addWidget(self.page_18)
        self.page_22 = QWidget()
        self.page_22.setObjectName(u"page_22")
        self.BoundaryConditionsOptions.addWidget(self.page_22)
        self.page_19 = QWidget()
        self.page_19.setObjectName(u"page_19")
        self.verticalLayout_6 = QVBoxLayout(self.page_19)
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.frame_52 = QFrame(self.page_19)
        self.frame_52.setObjectName(u"frame_52")
        self.frame_52.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_52.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_28 = QGridLayout(self.frame_52)
        self.gridLayout_28.setObjectName(u"gridLayout_28")
        self.gridLayout_28.setContentsMargins(0, 0, 0, 0)
        self.LbulkBC = QLabel(self.frame_52)
        self.LbulkBC.setObjectName(u"LbulkBC")

        self.gridLayout_28.addWidget(self.LbulkBC, 0, 0, 1, 1)

        self.INbulkBC = QComboBox(self.frame_52)
        self.INbulkBC.addItem("")
        self.INbulkBC.addItem("")
        self.INbulkBC.addItem("")
        self.INbulkBC.addItem("")
        self.INbulkBC.setObjectName(u"INbulkBC")

        self.gridLayout_28.addWidget(self.INbulkBC, 0, 1, 1, 1)


        self.verticalLayout_6.addWidget(self.frame_52, 0, Qt.AlignmentFlag.AlignHCenter)

        self.line = QFrame(self.page_19)
        self.line.setObjectName(u"line")
        self.line.setFrameShape(QFrame.Shape.HLine)
        self.line.setFrameShadow(QFrame.Shadow.Sunken)

        self.verticalLayout_6.addWidget(self.line)

        self.LVgrad = QLabel(self.page_19)
        self.LVgrad.setObjectName(u"LVgrad")
        self.LVgrad.setFont(font14)

        self.verticalLayout_6.addWidget(self.LVgrad, 0, Qt.AlignmentFlag.AlignHCenter)

        self.label_12 = QLabel(self.page_19)
        self.label_12.setObjectName(u"label_12")

        self.verticalLayout_6.addWidget(self.label_12, 0, Qt.AlignmentFlag.AlignHCenter)

        self.TVgrad = QTableWidget(self.page_19)
        if (self.TVgrad.columnCount() < 3):
            self.TVgrad.setColumnCount(3)
        if (self.TVgrad.rowCount() < 3):
            self.TVgrad.setRowCount(3)
        self.TVgrad.setObjectName(u"TVgrad")
        sizePolicy15 = QSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)
        sizePolicy15.setHorizontalStretch(0)
        sizePolicy15.setVerticalStretch(0)
        sizePolicy15.setHeightForWidth(self.TVgrad.sizePolicy().hasHeightForWidth())
        self.TVgrad.setSizePolicy(sizePolicy15)
        self.TVgrad.setMinimumSize(QSize(0, 0))
        self.TVgrad.setMaximumSize(QSize(16777215, 115))
        self.TVgrad.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.TVgrad.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustToContents)
        self.TVgrad.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.TVgrad.setRowCount(3)
        self.TVgrad.setColumnCount(3)
        self.TVgrad.horizontalHeader().setDefaultSectionSize(115)

        self.verticalLayout_6.addWidget(self.TVgrad)

        self.LVgradi = QLabel(self.page_19)
        self.LVgradi.setObjectName(u"LVgradi")
        font16 = QFont()
        font16.setBold(False)
        self.LVgradi.setFont(font16)

        self.verticalLayout_6.addWidget(self.LVgradi, 0, Qt.AlignmentFlag.AlignHCenter)

        self.TVgradi = QTableWidget(self.page_19)
        if (self.TVgradi.columnCount() < 3):
            self.TVgradi.setColumnCount(3)
        if (self.TVgradi.rowCount() < 3):
            self.TVgradi.setRowCount(3)
        self.TVgradi.setObjectName(u"TVgradi")
        sizePolicy15.setHeightForWidth(self.TVgradi.sizePolicy().hasHeightForWidth())
        self.TVgradi.setSizePolicy(sizePolicy15)
        self.TVgradi.setMinimumSize(QSize(0, 0))
        self.TVgradi.setMaximumSize(QSize(16777215, 115))
        self.TVgradi.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustToContents)
        self.TVgradi.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.TVgradi.setRowCount(3)
        self.TVgradi.setColumnCount(3)
        self.TVgradi.horizontalHeader().setDefaultSectionSize(115)

        self.verticalLayout_6.addWidget(self.TVgradi)

        self.LCstress = QLabel(self.page_19)
        self.LCstress.setObjectName(u"LCstress")
        self.LCstress.setFont(font14)

        self.verticalLayout_6.addWidget(self.LCstress, 0, Qt.AlignmentFlag.AlignHCenter)

        self.label_29 = QLabel(self.page_19)
        self.label_29.setObjectName(u"label_29")

        self.verticalLayout_6.addWidget(self.label_29, 0, Qt.AlignmentFlag.AlignHCenter)

        self.TCstress = QTableWidget(self.page_19)
        if (self.TCstress.columnCount() < 3):
            self.TCstress.setColumnCount(3)
        if (self.TCstress.rowCount() < 3):
            self.TCstress.setRowCount(3)
        self.TCstress.setObjectName(u"TCstress")
        sizePolicy15.setHeightForWidth(self.TCstress.sizePolicy().hasHeightForWidth())
        self.TCstress.setSizePolicy(sizePolicy15)
        self.TCstress.setMinimumSize(QSize(0, 0))
        self.TCstress.setMaximumSize(QSize(16777215, 115))
        self.TCstress.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustToContents)
        self.TCstress.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.TCstress.setRowCount(3)
        self.TCstress.setColumnCount(3)
        self.TCstress.horizontalHeader().setDefaultSectionSize(115)

        self.verticalLayout_6.addWidget(self.TCstress)

        self.verticalSpacer_11 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_6.addItem(self.verticalSpacer_11)

        self.BoundaryConditionsOptions.addWidget(self.page_19)

        self.verticalLayout_7.addWidget(self.BoundaryConditionsOptions)

        self.ToolWindow.addWidget(self.BoundaryConditionsTab)
        self.SolverSettingsSGHTool = QWidget()
        self.SolverSettingsSGHTool.setObjectName(u"SolverSettingsSGHTool")
        self.verticalLayout_28 = QVBoxLayout(self.SolverSettingsSGHTool)
        self.verticalLayout_28.setObjectName(u"verticalLayout_28")
        self.verticalLayout_28.setContentsMargins(10, 7, -1, -1)
        self.LSolverSettings = QLabel(self.SolverSettingsSGHTool)
        self.LSolverSettings.setObjectName(u"LSolverSettings")
        sizePolicy8.setHeightForWidth(self.LSolverSettings.sizePolicy().hasHeightForWidth())
        self.LSolverSettings.setSizePolicy(sizePolicy8)

        self.verticalLayout_28.addWidget(self.LSolverSettings)

        self.INSelectSolverSettings = QComboBox(self.SolverSettingsSGHTool)
        self.INSelectSolverSettings.addItem("")
        self.INSelectSolverSettings.addItem("")
        self.INSelectSolverSettings.addItem("")
        self.INSelectSolverSettings.setObjectName(u"INSelectSolverSettings")

        self.verticalLayout_28.addWidget(self.INSelectSolverSettings)

        self.SolverSettingsOptions = QStackedWidget(self.SolverSettingsSGHTool)
        self.SolverSettingsOptions.setObjectName(u"SolverSettingsOptions")
        self.page_2 = QWidget()
        self.page_2.setObjectName(u"page_2")
        self.solversettings = QFrame(self.page_2)
        self.solversettings.setObjectName(u"solversettings")
        self.solversettings.setGeometry(QRect(5, 11, 373, 189))
        self.solversettings.setFrameShape(QFrame.Shape.NoFrame)
        self.formLayout_13 = QFormLayout(self.solversettings)
        self.formLayout_13.setObjectName(u"formLayout_13")
        self.formLayout_13.setVerticalSpacing(6)
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

        self.LGraphicsOutput = QLabel(self.solversettings)
        self.LGraphicsOutput.setObjectName(u"LGraphicsOutput")
        self.LGraphicsOutput.setWordWrap(True)

        self.formLayout_13.setWidget(5, QFormLayout.LabelRole, self.LGraphicsOutput)

        self.INGraphicsOutput = QLineEdit(self.solversettings)
        self.INGraphicsOutput.setObjectName(u"INGraphicsOutput")

        self.formLayout_13.setWidget(5, QFormLayout.FieldRole, self.INGraphicsOutput)

        self.SolverSettingsOptions.addWidget(self.page_2)
        self.page_3 = QWidget()
        self.page_3.setObjectName(u"page_3")
        self.verticalLayout_11 = QVBoxLayout(self.page_3)
        self.verticalLayout_11.setObjectName(u"verticalLayout_11")
        self.frame_2 = QFrame(self.page_3)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_2.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout = QGridLayout(self.frame_2)
        self.gridLayout.setObjectName(u"gridLayout")
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.INHomogenizationSettingD = QRadioButton(self.frame_2)
        self.INHomogenizationSettingD.setObjectName(u"INHomogenizationSettingD")
        self.INHomogenizationSettingD.setChecked(True)

        self.gridLayout.addWidget(self.INHomogenizationSettingD, 0, 0, 1, 1)

        self.INHomogenizationSettingC = QRadioButton(self.frame_2)
        self.INHomogenizationSettingC.setObjectName(u"INHomogenizationSettingC")

        self.gridLayout.addWidget(self.INHomogenizationSettingC, 1, 0, 1, 1)


        self.verticalLayout_11.addWidget(self.frame_2)

        self.frame = QFrame(self.page_3)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.Shape.NoFrame)
        self.frame.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_3 = QFormLayout(self.frame)
        self.formLayout_3.setObjectName(u"formLayout_3")
        self.formLayout_3.setContentsMargins(0, 0, 0, 0)
        self.LNumberOfSteps = QLabel(self.frame)
        self.LNumberOfSteps.setObjectName(u"LNumberOfSteps")

        self.formLayout_3.setWidget(0, QFormLayout.LabelRole, self.LNumberOfSteps)

        self.INNumberOfSteps = QLineEdit(self.frame)
        self.INNumberOfSteps.setObjectName(u"INNumberOfSteps")
        self.INNumberOfSteps.setEnabled(False)
        self.INNumberOfSteps.setReadOnly(False)

        self.formLayout_3.setWidget(0, QFormLayout.FieldRole, self.INNumberOfSteps)

        self.LErrorTolerance = QLabel(self.frame)
        self.LErrorTolerance.setObjectName(u"LErrorTolerance")

        self.formLayout_3.setWidget(1, QFormLayout.LabelRole, self.LErrorTolerance)

        self.INErrorTolerance = QLineEdit(self.frame)
        self.INErrorTolerance.setObjectName(u"INErrorTolerance")
        self.INErrorTolerance.setEnabled(False)

        self.formLayout_3.setWidget(1, QFormLayout.FieldRole, self.INErrorTolerance)

        self.LMaxIterations = QLabel(self.frame)
        self.LMaxIterations.setObjectName(u"LMaxIterations")

        self.formLayout_3.setWidget(2, QFormLayout.LabelRole, self.LMaxIterations)

        self.INMaxIterations = QLineEdit(self.frame)
        self.INMaxIterations.setObjectName(u"INMaxIterations")
        self.INMaxIterations.setEnabled(False)

        self.formLayout_3.setWidget(2, QFormLayout.FieldRole, self.INMaxIterations)


        self.verticalLayout_11.addWidget(self.frame)

        self.verticalSpacer_12 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_11.addItem(self.verticalSpacer_12)

        self.SolverSettingsOptions.addWidget(self.page_3)
        self.page_23 = QWidget()
        self.page_23.setObjectName(u"page_23")
        self.verticalLayout_34 = QVBoxLayout(self.page_23)
        self.verticalLayout_34.setObjectName(u"verticalLayout_34")
        self.label_69 = QLabel(self.page_23)
        self.label_69.setObjectName(u"label_69")
        self.label_69.setFont(font14)

        self.verticalLayout_34.addWidget(self.label_69, 0, Qt.AlignmentFlag.AlignHCenter)

        self.frame_73 = QFrame(self.page_23)
        self.frame_73.setObjectName(u"frame_73")
        self.frame_73.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_73.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_27 = QHBoxLayout(self.frame_73)
        self.horizontalLayout_27.setObjectName(u"horizontalLayout_27")
        self.horizontalLayout_27.setContentsMargins(0, 0, 0, 0)
        self.INLargeStrain = QRadioButton(self.frame_73)
        self.INLargeStrain.setObjectName(u"INLargeStrain")
        self.INLargeStrain.setChecked(True)

        self.horizontalLayout_27.addWidget(self.INLargeStrain)

        self.INSmallStrain = QRadioButton(self.frame_73)
        self.INSmallStrain.setObjectName(u"INSmallStrain")

        self.horizontalLayout_27.addWidget(self.INSmallStrain)


        self.verticalLayout_34.addWidget(self.frame_73, 0, Qt.AlignmentFlag.AlignHCenter)

        self.label_68 = QLabel(self.page_23)
        self.label_68.setObjectName(u"label_68")
        self.label_68.setFont(font14)

        self.verticalLayout_34.addWidget(self.label_68, 0, Qt.AlignmentFlag.AlignHCenter)

        self.frame_54 = QFrame(self.page_23)
        self.frame_54.setObjectName(u"frame_54")
        self.frame_54.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_54.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_8 = QFormLayout(self.frame_54)
        self.formLayout_8.setObjectName(u"formLayout_8")
        self.formLayout_8.setContentsMargins(0, 0, 0, 0)
        self.LNumberOfSteps_3 = QLabel(self.frame_54)
        self.LNumberOfSteps_3.setObjectName(u"LNumberOfSteps_3")

        self.formLayout_8.setWidget(1, QFormLayout.LabelRole, self.LNumberOfSteps_3)

        self.INBFloadsteps = QLineEdit(self.frame_54)
        self.INBFloadsteps.setObjectName(u"INBFloadsteps")
        self.INBFloadsteps.setEnabled(True)
        self.INBFloadsteps.setReadOnly(False)

        self.formLayout_8.setWidget(1, QFormLayout.FieldRole, self.INBFloadsteps)

        self.LErrorTolerance_3 = QLabel(self.frame_54)
        self.LErrorTolerance_3.setObjectName(u"LErrorTolerance_3")

        self.formLayout_8.setWidget(3, QFormLayout.LabelRole, self.LErrorTolerance_3)

        self.INBFerrortol = QLineEdit(self.frame_54)
        self.INBFerrortol.setObjectName(u"INBFerrortol")
        self.INBFerrortol.setEnabled(True)

        self.formLayout_8.setWidget(3, QFormLayout.FieldRole, self.INBFerrortol)

        self.LMaxIterations_3 = QLabel(self.frame_54)
        self.LMaxIterations_3.setObjectName(u"LMaxIterations_3")

        self.formLayout_8.setWidget(4, QFormLayout.LabelRole, self.LMaxIterations_3)

        self.INBFmaxiter = QLineEdit(self.frame_54)
        self.INBFmaxiter.setObjectName(u"INBFmaxiter")
        self.INBFmaxiter.setEnabled(True)

        self.formLayout_8.setWidget(4, QFormLayout.FieldRole, self.INBFmaxiter)

        self.LBFdt = QLabel(self.frame_54)
        self.LBFdt.setObjectName(u"LBFdt")

        self.formLayout_8.setWidget(0, QFormLayout.LabelRole, self.LBFdt)

        self.INBFdt = QLineEdit(self.frame_54)
        self.INBFdt.setObjectName(u"INBFdt")
        self.INBFdt.setEnabled(True)

        self.formLayout_8.setWidget(0, QFormLayout.FieldRole, self.INBFdt)

        self.LBFoutputsteps = QLabel(self.frame_54)
        self.LBFoutputsteps.setObjectName(u"LBFoutputsteps")

        self.formLayout_8.setWidget(2, QFormLayout.LabelRole, self.LBFoutputsteps)

        self.INBFoutputsteps = QLineEdit(self.frame_54)
        self.INBFoutputsteps.setObjectName(u"INBFoutputsteps")
        self.INBFoutputsteps.setEnabled(True)

        self.formLayout_8.setWidget(2, QFormLayout.FieldRole, self.INBFoutputsteps)


        self.verticalLayout_34.addWidget(self.frame_54, 0, Qt.AlignmentFlag.AlignTop)

        self.verticalSpacer_28 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_34.addItem(self.verticalSpacer_28)

        self.SolverSettingsOptions.addWidget(self.page_23)

        self.verticalLayout_28.addWidget(self.SolverSettingsOptions)

        self.ToolWindow.addWidget(self.SolverSettingsSGHTool)
        self.page_7 = QWidget()
        self.page_7.setObjectName(u"page_7")
        self.verticalLayout_9 = QVBoxLayout(self.page_7)
        self.verticalLayout_9.setObjectName(u"verticalLayout_9")
        self.LRun = QLabel(self.page_7)
        self.LRun.setObjectName(u"LRun")

        self.verticalLayout_9.addWidget(self.LRun)

        self.frame_26 = QFrame(self.page_7)
        self.frame_26.setObjectName(u"frame_26")
        self.frame_26.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_26.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_4 = QFormLayout(self.frame_26)
        self.formLayout_4.setObjectName(u"formLayout_4")
        self.formLayout_4.setContentsMargins(0, 0, 0, 0)
        self.label_3 = QLabel(self.frame_26)
        self.label_3.setObjectName(u"label_3")

        self.formLayout_4.setWidget(0, QFormLayout.LabelRole, self.label_3)

        self.INRunSelection = QComboBox(self.frame_26)
        self.INRunSelection.addItem("")
        self.INRunSelection.addItem("")
        self.INRunSelection.addItem("")
        self.INRunSelection.setObjectName(u"INRunSelection")

        self.formLayout_4.setWidget(0, QFormLayout.FieldRole, self.INRunSelection)


        self.verticalLayout_9.addWidget(self.frame_26, 0, Qt.AlignmentFlag.AlignVCenter)

        self.RunOptions = QStackedWidget(self.page_7)
        self.RunOptions.setObjectName(u"RunOptions")
        self.page_20 = QWidget()
        self.page_20.setObjectName(u"page_20")
        self.verticalLayout_12 = QVBoxLayout(self.page_20)
        self.verticalLayout_12.setObjectName(u"verticalLayout_12")
        self.BRunSGH = QPushButton(self.page_20)
        self.BRunSGH.setObjectName(u"BRunSGH")

        self.verticalLayout_12.addWidget(self.BRunSGH)

        self.verticalSpacer_19 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_12.addItem(self.verticalSpacer_19)

        self.RunOptions.addWidget(self.page_20)
        self.page_21 = QWidget()
        self.page_21.setObjectName(u"page_21")
        self.verticalLayout_13 = QVBoxLayout(self.page_21)
        self.verticalLayout_13.setObjectName(u"verticalLayout_13")
        self.frame_71 = QFrame(self.page_21)
        self.frame_71.setObjectName(u"frame_71")
        self.frame_71.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_71.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_26 = QFormLayout(self.frame_71)
        self.formLayout_26.setObjectName(u"formLayout_26")
        self.formLayout_26.setContentsMargins(0, 0, 0, 5)
        self.INHRunLocally = QRadioButton(self.frame_71)
        self.INHRunLocally.setObjectName(u"INHRunLocally")
        self.INHRunLocally.setChecked(True)

        self.formLayout_26.setWidget(0, QFormLayout.LabelRole, self.INHRunLocally)

        self.INHWriteFiles = QRadioButton(self.frame_71)
        self.INHWriteFiles.setObjectName(u"INHWriteFiles")

        self.formLayout_26.setWidget(0, QFormLayout.FieldRole, self.INHWriteFiles)


        self.verticalLayout_13.addWidget(self.frame_71, 0, Qt.AlignmentFlag.AlignHCenter)

        self.HomogenizationRunWrite = QStackedWidget(self.page_21)
        self.HomogenizationRunWrite.setObjectName(u"HomogenizationRunWrite")
        self.page_37 = QWidget()
        self.page_37.setObjectName(u"page_37")
        self.verticalLayout_73 = QVBoxLayout(self.page_37)
        self.verticalLayout_73.setObjectName(u"verticalLayout_73")
        self.label_53 = QLabel(self.page_37)
        self.label_53.setObjectName(u"label_53")
        self.label_53.setFont(font14)

        self.verticalLayout_73.addWidget(self.label_53, 0, Qt.AlignmentFlag.AlignHCenter)

        self.frame_16 = QFrame(self.page_37)
        self.frame_16.setObjectName(u"frame_16")
        self.frame_16.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_16.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout = QFormLayout(self.frame_16)
        self.formLayout.setObjectName(u"formLayout")
        self.formLayout.setContentsMargins(0, 5, 0, 5)
        self.INSingleJob = QRadioButton(self.frame_16)
        self.INSingleJob.setObjectName(u"INSingleJob")
        self.INSingleJob.setChecked(True)

        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.INSingleJob)

        self.INBatchJob = QRadioButton(self.frame_16)
        self.INBatchJob.setObjectName(u"INBatchJob")

        self.formLayout.setWidget(0, QFormLayout.FieldRole, self.INBatchJob)


        self.verticalLayout_73.addWidget(self.frame_16, 0, Qt.AlignmentFlag.AlignHCenter)

        self.HomogenizationBatch = QStackedWidget(self.page_37)
        self.HomogenizationBatch.setObjectName(u"HomogenizationBatch")
        self.page_5 = QWidget()
        self.page_5.setObjectName(u"page_5")
        self.HomogenizationBatch.addWidget(self.page_5)
        self.page_6 = QWidget()
        self.page_6.setObjectName(u"page_6")
        self.verticalLayout_17 = QVBoxLayout(self.page_6)
        self.verticalLayout_17.setObjectName(u"verticalLayout_17")
        self.verticalLayout_17.setContentsMargins(0, 0, 0, 0)
        self.frame_19 = QFrame(self.page_6)
        self.frame_19.setObjectName(u"frame_19")
        self.frame_19.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_19.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_7 = QFormLayout(self.frame_19)
        self.formLayout_7.setObjectName(u"formLayout_7")
        self.formLayout_7.setContentsMargins(0, 0, 0, 0)
        self.LBatchType = QLabel(self.frame_19)
        self.LBatchType.setObjectName(u"LBatchType")

        self.formLayout_7.setWidget(0, QFormLayout.LabelRole, self.LBatchType)

        self.INBatchType = QComboBox(self.frame_19)
        self.INBatchType.addItem("")
        self.INBatchType.setObjectName(u"INBatchType")

        self.formLayout_7.setWidget(0, QFormLayout.FieldRole, self.INBatchType)


        self.verticalLayout_17.addWidget(self.frame_19)

        self.frame_17 = QFrame(self.page_6)
        self.frame_17.setObjectName(u"frame_17")
        sizePolicy.setHeightForWidth(self.frame_17.sizePolicy().hasHeightForWidth())
        self.frame_17.setSizePolicy(sizePolicy)
        self.frame_17.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_17.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_2 = QFormLayout(self.frame_17)
        self.formLayout_2.setObjectName(u"formLayout_2")
        self.formLayout_2.setContentsMargins(0, 0, 0, 0)
        self.INSaveMaterialFiles = QRadioButton(self.frame_17)
        self.INSaveMaterialFiles.setObjectName(u"INSaveMaterialFiles")
        self.INSaveMaterialFiles.setChecked(True)

        self.formLayout_2.setWidget(0, QFormLayout.LabelRole, self.INSaveMaterialFiles)

        self.INSaveAllFiles = QRadioButton(self.frame_17)
        self.INSaveAllFiles.setObjectName(u"INSaveAllFiles")

        self.formLayout_2.setWidget(0, QFormLayout.FieldRole, self.INSaveAllFiles)


        self.verticalLayout_17.addWidget(self.frame_17, 0, Qt.AlignmentFlag.AlignHCenter)

        self.BSelectGeometryFiles = QPushButton(self.page_6)
        self.BSelectGeometryFiles.setObjectName(u"BSelectGeometryFiles")

        self.verticalLayout_17.addWidget(self.BSelectGeometryFiles)

        self.frame_21 = QFrame(self.page_6)
        self.frame_21.setObjectName(u"frame_21")
        self.frame_21.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_21.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout = QHBoxLayout(self.frame_21)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.LHomogenizingFile = QLabel(self.frame_21)
        self.LHomogenizingFile.setObjectName(u"LHomogenizingFile")
        self.LHomogenizingFile.setEnabled(False)
        font17 = QFont()
        font17.setItalic(True)
        self.LHomogenizingFile.setFont(font17)

        self.horizontalLayout.addWidget(self.LHomogenizingFile)

        self.INFileNumber = QLabel(self.frame_21)
        self.INFileNumber.setObjectName(u"INFileNumber")
        self.INFileNumber.setEnabled(False)

        self.horizontalLayout.addWidget(self.INFileNumber)


        self.verticalLayout_17.addWidget(self.frame_21, 0, Qt.AlignmentFlag.AlignHCenter)

        self.INHomogenizationBatchFile = QLineEdit(self.page_6)
        self.INHomogenizationBatchFile.setObjectName(u"INHomogenizationBatchFile")
        self.INHomogenizationBatchFile.setEnabled(False)

        self.verticalLayout_17.addWidget(self.INHomogenizationBatchFile)

        self.HomogenizationBatch.addWidget(self.page_6)

        self.verticalLayout_73.addWidget(self.HomogenizationBatch, 0, Qt.AlignmentFlag.AlignTop)

        self.label_54 = QLabel(self.page_37)
        self.label_54.setObjectName(u"label_54")
        self.label_54.setFont(font14)

        self.verticalLayout_73.addWidget(self.label_54, 0, Qt.AlignmentFlag.AlignHCenter)

        self.frame_59 = QFrame(self.page_37)
        self.frame_59.setObjectName(u"frame_59")
        self.frame_59.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_59.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_6 = QFormLayout(self.frame_59)
        self.formLayout_6.setObjectName(u"formLayout_6")
        self.formLayout_6.setContentsMargins(0, 0, 0, 5)
        self.INHSerial = QRadioButton(self.frame_59)
        self.INHSerial.setObjectName(u"INHSerial")
        self.INHSerial.setChecked(True)

        self.formLayout_6.setWidget(0, QFormLayout.LabelRole, self.INHSerial)

        self.INHParallel = QRadioButton(self.frame_59)
        self.INHParallel.setObjectName(u"INHParallel")

        self.formLayout_6.setWidget(0, QFormLayout.FieldRole, self.INHParallel)


        self.verticalLayout_73.addWidget(self.frame_59, 0, Qt.AlignmentFlag.AlignHCenter)

        self.HomogenizationRunType = QStackedWidget(self.page_37)
        self.HomogenizationRunType.setObjectName(u"HomogenizationRunType")
        sizePolicy8.setHeightForWidth(self.HomogenizationRunType.sizePolicy().hasHeightForWidth())
        self.HomogenizationRunType.setSizePolicy(sizePolicy8)
        self.HomogenizationRunType.setMinimumSize(QSize(0, 75))
        self.HomogenizationRunType.setMaximumSize(QSize(16777215, 75))
        self.page_31 = QWidget()
        self.page_31.setObjectName(u"page_31")
        self.HomogenizationRunType.addWidget(self.page_31)
        self.page_32 = QWidget()
        self.page_32.setObjectName(u"page_32")
        self.verticalLayout_72 = QVBoxLayout(self.page_32)
        self.verticalLayout_72.setObjectName(u"verticalLayout_72")
        self.verticalLayout_72.setContentsMargins(0, 0, 0, 0)
        self.frame_65 = QFrame(self.page_32)
        self.frame_65.setObjectName(u"frame_65")
        self.frame_65.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_65.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_74 = QVBoxLayout(self.frame_65)
        self.verticalLayout_74.setObjectName(u"verticalLayout_74")
        self.verticalLayout_74.setContentsMargins(0, 0, 0, 0)
        self.frame_70 = QFrame(self.frame_65)
        self.frame_70.setObjectName(u"frame_70")
        sizePolicy8.setHeightForWidth(self.frame_70.sizePolicy().hasHeightForWidth())
        self.frame_70.setSizePolicy(sizePolicy8)
        self.frame_70.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_70.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_9 = QFormLayout(self.frame_70)
        self.formLayout_9.setObjectName(u"formLayout_9")
        self.formLayout_9.setContentsMargins(0, 0, 0, 0)
        self.label_58 = QLabel(self.frame_70)
        self.label_58.setObjectName(u"label_58")

        self.formLayout_9.setWidget(0, QFormLayout.LabelRole, self.label_58)

        self.INmpiRanks = QSpinBox(self.frame_70)
        self.INmpiRanks.setObjectName(u"INmpiRanks")
        self.INmpiRanks.setReadOnly(False)
        self.INmpiRanks.setButtonSymbols(QAbstractSpinBox.ButtonSymbols.UpDownArrows)
        self.INmpiRanks.setMinimum(1)
        self.INmpiRanks.setMaximum(1)

        self.formLayout_9.setWidget(0, QFormLayout.FieldRole, self.INmpiRanks)


        self.verticalLayout_74.addWidget(self.frame_70)

        self.label_65 = QLabel(self.frame_65)
        self.label_65.setObjectName(u"label_65")
        self.label_65.setWordWrap(True)

        self.verticalLayout_74.addWidget(self.label_65, 0, Qt.AlignmentFlag.AlignTop)


        self.verticalLayout_72.addWidget(self.frame_65)

        self.HomogenizationRunType.addWidget(self.page_32)

        self.verticalLayout_73.addWidget(self.HomogenizationRunType, 0, Qt.AlignmentFlag.AlignTop)

        self.label_64 = QLabel(self.page_37)
        self.label_64.setObjectName(u"label_64")
        self.label_64.setFont(font14)

        self.verticalLayout_73.addWidget(self.label_64, 0, Qt.AlignmentFlag.AlignHCenter)

        self.frame_64 = QFrame(self.page_37)
        self.frame_64.setObjectName(u"frame_64")
        self.frame_64.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_64.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_19 = QFormLayout(self.frame_64)
        self.formLayout_19.setObjectName(u"formLayout_19")
        self.formLayout_19.setContentsMargins(0, 0, 0, 5)
        self.INHAutomatic = QRadioButton(self.frame_64)
        self.INHAutomatic.setObjectName(u"INHAutomatic")
        self.INHAutomatic.setChecked(True)

        self.formLayout_19.setWidget(0, QFormLayout.LabelRole, self.INHAutomatic)

        self.INHManual = QRadioButton(self.frame_64)
        self.INHManual.setObjectName(u"INHManual")

        self.formLayout_19.setWidget(0, QFormLayout.FieldRole, self.INHManual)


        self.verticalLayout_73.addWidget(self.frame_64, 0, Qt.AlignmentFlag.AlignHCenter)

        self.frame_68 = QFrame(self.page_37)
        self.frame_68.setObjectName(u"frame_68")
        self.frame_68.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_68.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_25 = QHBoxLayout(self.frame_68)
        self.horizontalLayout_25.setObjectName(u"horizontalLayout_25")
        self.horizontalLayout_25.setContentsMargins(0, 0, 0, 0)
        self.BRunEVPFFT2 = QPushButton(self.frame_68)
        self.BRunEVPFFT2.setObjectName(u"BRunEVPFFT2")
        self.BRunEVPFFT2.setFont(font)

        self.horizontalLayout_25.addWidget(self.BRunEVPFFT2)

        self.BKillEVPFFT2 = QPushButton(self.frame_68)
        self.BKillEVPFFT2.setObjectName(u"BKillEVPFFT2")
        self.BKillEVPFFT2.setFlat(False)

        self.horizontalLayout_25.addWidget(self.BKillEVPFFT2)


        self.verticalLayout_73.addWidget(self.frame_68)

        self.verticalSpacer_9 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_73.addItem(self.verticalSpacer_9)

        self.HomogenizationRunWrite.addWidget(self.page_37)
        self.page_38 = QWidget()
        self.page_38.setObjectName(u"page_38")
        self.verticalLayout_75 = QVBoxLayout(self.page_38)
        self.verticalLayout_75.setObjectName(u"verticalLayout_75")
        self.BHWriteFiles = QPushButton(self.page_38)
        self.BHWriteFiles.setObjectName(u"BHWriteFiles")

        self.verticalLayout_75.addWidget(self.BHWriteFiles)

        self.verticalSpacer_21 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_75.addItem(self.verticalSpacer_21)

        self.HomogenizationRunWrite.addWidget(self.page_38)

        self.verticalLayout_13.addWidget(self.HomogenizationRunWrite)

        self.RunOptions.addWidget(self.page_21)
        self.page_24 = QWidget()
        self.page_24.setObjectName(u"page_24")
        self.verticalLayout_39 = QVBoxLayout(self.page_24)
        self.verticalLayout_39.setObjectName(u"verticalLayout_39")
        self.frame_72 = QFrame(self.page_24)
        self.frame_72.setObjectName(u"frame_72")
        self.frame_72.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_72.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_29 = QFormLayout(self.frame_72)
        self.formLayout_29.setObjectName(u"formLayout_29")
        self.formLayout_29.setContentsMargins(0, 0, 0, 5)
        self.INBFRunLocally = QRadioButton(self.frame_72)
        self.INBFRunLocally.setObjectName(u"INBFRunLocally")
        self.INBFRunLocally.setChecked(True)

        self.formLayout_29.setWidget(0, QFormLayout.LabelRole, self.INBFRunLocally)

        self.INBFWriteFiles = QRadioButton(self.frame_72)
        self.INBFWriteFiles.setObjectName(u"INBFWriteFiles")

        self.formLayout_29.setWidget(0, QFormLayout.FieldRole, self.INBFWriteFiles)


        self.verticalLayout_39.addWidget(self.frame_72, 0, Qt.AlignmentFlag.AlignHCenter)

        self.BFRunWrite = QStackedWidget(self.page_24)
        self.BFRunWrite.setObjectName(u"BFRunWrite")
        self.page_39 = QWidget()
        self.page_39.setObjectName(u"page_39")
        self.verticalLayout_76 = QVBoxLayout(self.page_39)
        self.verticalLayout_76.setObjectName(u"verticalLayout_76")
        self.label_62 = QLabel(self.page_39)
        self.label_62.setObjectName(u"label_62")
        self.label_62.setFont(font14)

        self.verticalLayout_76.addWidget(self.label_62, 0, Qt.AlignmentFlag.AlignHCenter)

        self.frame_63 = QFrame(self.page_39)
        self.frame_63.setObjectName(u"frame_63")
        self.frame_63.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_63.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_11 = QFormLayout(self.frame_63)
        self.formLayout_11.setObjectName(u"formLayout_11")
        self.formLayout_11.setContentsMargins(0, 0, 0, 5)
        self.INBFSerial = QRadioButton(self.frame_63)
        self.INBFSerial.setObjectName(u"INBFSerial")
        self.INBFSerial.setChecked(True)

        self.formLayout_11.setWidget(0, QFormLayout.LabelRole, self.INBFSerial)

        self.INBFParallel = QRadioButton(self.frame_63)
        self.INBFParallel.setObjectName(u"INBFParallel")

        self.formLayout_11.setWidget(0, QFormLayout.FieldRole, self.INBFParallel)


        self.verticalLayout_76.addWidget(self.frame_63, 0, Qt.AlignmentFlag.AlignHCenter)

        self.BFRunType = QStackedWidget(self.page_39)
        self.BFRunType.setObjectName(u"BFRunType")
        sizePolicy6.setHeightForWidth(self.BFRunType.sizePolicy().hasHeightForWidth())
        self.BFRunType.setSizePolicy(sizePolicy6)
        self.page_33 = QWidget()
        self.page_33.setObjectName(u"page_33")
        self.BFRunType.addWidget(self.page_33)
        self.page_34 = QWidget()
        self.page_34.setObjectName(u"page_34")
        self.formLayout_12 = QFormLayout(self.page_34)
        self.formLayout_12.setObjectName(u"formLayout_12")
        self.formLayout_12.setContentsMargins(0, 0, 0, 0)
        self.INBFmpiRanks = QSpinBox(self.page_34)
        self.INBFmpiRanks.setObjectName(u"INBFmpiRanks")
        self.INBFmpiRanks.setReadOnly(False)
        self.INBFmpiRanks.setButtonSymbols(QAbstractSpinBox.ButtonSymbols.UpDownArrows)
        self.INBFmpiRanks.setMinimum(1)
        self.INBFmpiRanks.setMaximum(1)

        self.formLayout_12.setWidget(0, QFormLayout.FieldRole, self.INBFmpiRanks)

        self.label_63 = QLabel(self.page_34)
        self.label_63.setObjectName(u"label_63")

        self.formLayout_12.setWidget(0, QFormLayout.LabelRole, self.label_63)

        self.BFRunType.addWidget(self.page_34)

        self.verticalLayout_76.addWidget(self.BFRunType)

        self.frame_69 = QFrame(self.page_39)
        self.frame_69.setObjectName(u"frame_69")
        self.frame_69.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_69.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_26 = QHBoxLayout(self.frame_69)
        self.horizontalLayout_26.setObjectName(u"horizontalLayout_26")
        self.horizontalLayout_26.setContentsMargins(0, 0, 0, 0)
        self.BRunBulkForming = QPushButton(self.frame_69)
        self.BRunBulkForming.setObjectName(u"BRunBulkForming")
        self.BRunBulkForming.setFont(font)

        self.horizontalLayout_26.addWidget(self.BRunBulkForming)

        self.BKillBulkForming = QPushButton(self.frame_69)
        self.BKillBulkForming.setObjectName(u"BKillBulkForming")

        self.horizontalLayout_26.addWidget(self.BKillBulkForming)


        self.verticalLayout_76.addWidget(self.frame_69)

        self.BFRunWrite.addWidget(self.page_39)
        self.page_40 = QWidget()
        self.page_40.setObjectName(u"page_40")
        self.verticalLayout_77 = QVBoxLayout(self.page_40)
        self.verticalLayout_77.setObjectName(u"verticalLayout_77")
        self.BBFWriteFiles = QPushButton(self.page_40)
        self.BBFWriteFiles.setObjectName(u"BBFWriteFiles")

        self.verticalLayout_77.addWidget(self.BBFWriteFiles, 0, Qt.AlignmentFlag.AlignTop)

        self.verticalSpacer_31 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_77.addItem(self.verticalSpacer_31)

        self.BFRunWrite.addWidget(self.page_40)

        self.verticalLayout_39.addWidget(self.BFRunWrite, 0, Qt.AlignmentFlag.AlignTop)

        self.verticalSpacer_29 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_39.addItem(self.verticalSpacer_29)

        self.RunOptions.addWidget(self.page_24)

        self.verticalLayout_9.addWidget(self.RunOptions)

        self.ToolWindow.addWidget(self.page_7)
        self.ResultsEVPFFTTool = QWidget()
        self.ResultsEVPFFTTool.setObjectName(u"ResultsEVPFFTTool")
        self.verticalLayout_18 = QVBoxLayout(self.ResultsEVPFFTTool)
        self.verticalLayout_18.setObjectName(u"verticalLayout_18")
        self.LPostprocessing = QLabel(self.ResultsEVPFFTTool)
        self.LPostprocessing.setObjectName(u"LPostprocessing")

        self.verticalLayout_18.addWidget(self.LPostprocessing)

        self.tabWidget_2 = QTabWidget(self.ResultsEVPFFTTool)
        self.tabWidget_2.setObjectName(u"tabWidget_2")
        sizePolicy3.setHeightForWidth(self.tabWidget_2.sizePolicy().hasHeightForWidth())
        self.tabWidget_2.setSizePolicy(sizePolicy3)
        self.tabWidget_2.setIconSize(QSize(20, 20))
        self.tabWidget_2.setElideMode(Qt.TextElideMode.ElideLeft)
        self.Results = QWidget()
        self.Results.setObjectName(u"Results")
        self.verticalLayout_46 = QVBoxLayout(self.Results)
        self.verticalLayout_46.setObjectName(u"verticalLayout_46")
        self.INSelectPostprocessing = QComboBox(self.Results)
        self.INSelectPostprocessing.addItem("")
        self.INSelectPostprocessing.addItem("")
        self.INSelectPostprocessing.addItem("")
        self.INSelectPostprocessing.setObjectName(u"INSelectPostprocessing")

        self.verticalLayout_46.addWidget(self.INSelectPostprocessing)

        self.PostprocessingOptions = QStackedWidget(self.Results)
        self.PostprocessingOptions.setObjectName(u"PostprocessingOptions")
        self.page_11 = QWidget()
        self.page_11.setObjectName(u"page_11")
        self.verticalLayout_45 = QVBoxLayout(self.page_11)
        self.verticalLayout_45.setObjectName(u"verticalLayout_45")
        self.LResultsSGH = QLabel(self.page_11)
        self.LResultsSGH.setObjectName(u"LResultsSGH")

        self.verticalLayout_45.addWidget(self.LResultsSGH)

        self.frame_15 = QFrame(self.page_11)
        self.frame_15.setObjectName(u"frame_15")
        self.frame_15.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_17 = QHBoxLayout(self.frame_15)
        self.horizontalLayout_17.setSpacing(0)
        self.horizontalLayout_17.setObjectName(u"horizontalLayout_17")
        self.horizontalLayout_17.setContentsMargins(0, 0, 0, 0)
        self.INOuputVarSGH = QComboBox(self.frame_15)
        self.INOuputVarSGH.addItem("")
        self.INOuputVarSGH.setObjectName(u"INOuputVarSGH")

        self.horizontalLayout_17.addWidget(self.INOuputVarSGH)

        self.BPreviewResultsSGH = QPushButton(self.frame_15)
        self.BPreviewResultsSGH.setObjectName(u"BPreviewResultsSGH")

        self.horizontalLayout_17.addWidget(self.BPreviewResultsSGH)


        self.verticalLayout_45.addWidget(self.frame_15)

        self.frame_13 = QFrame(self.page_11)
        self.frame_13.setObjectName(u"frame_13")
        self.frame_13.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_12 = QHBoxLayout(self.frame_13)
        self.horizontalLayout_12.setObjectName(u"horizontalLayout_12")
        self.horizontalLayout_12.setContentsMargins(0, 0, 0, 0)
        self.BFirstFrame = QToolButton(self.frame_13)
        self.BFirstFrame.setObjectName(u"BFirstFrame")
        icon14 = QIcon()
        icon14.addFile(u":/Blue Icons/Blue Icons/FirstFrame.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BFirstFrame.setIcon(icon14)
        self.BFirstFrame.setIconSize(QSize(32, 32))
        self.BFirstFrame.setPopupMode(QToolButton.ToolButtonPopupMode.DelayedPopup)
        self.BFirstFrame.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        self.BFirstFrame.setAutoRaise(False)
        self.BFirstFrame.setArrowType(Qt.ArrowType.NoArrow)

        self.horizontalLayout_12.addWidget(self.BFirstFrame)

        self.BPreviousFrame = QToolButton(self.frame_13)
        self.BPreviousFrame.setObjectName(u"BPreviousFrame")
        icon15 = QIcon()
        icon15.addFile(u":/Blue Icons/Blue Icons/PreviousFrame.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BPreviousFrame.setIcon(icon15)
        self.BPreviousFrame.setIconSize(QSize(32, 32))

        self.horizontalLayout_12.addWidget(self.BPreviousFrame)

        self.BNextFrame = QToolButton(self.frame_13)
        self.BNextFrame.setObjectName(u"BNextFrame")
        icon16 = QIcon()
        icon16.addFile(u":/Blue Icons/Blue Icons/NextFrame.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BNextFrame.setIcon(icon16)
        self.BNextFrame.setIconSize(QSize(32, 32))

        self.horizontalLayout_12.addWidget(self.BNextFrame)

        self.BLastFrame = QToolButton(self.frame_13)
        self.BLastFrame.setObjectName(u"BLastFrame")
        icon17 = QIcon()
        icon17.addFile(u":/Blue Icons/Blue Icons/LastFrame.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BLastFrame.setIcon(icon17)
        self.BLastFrame.setIconSize(QSize(32, 32))

        self.horizontalLayout_12.addWidget(self.BLastFrame)


        self.verticalLayout_45.addWidget(self.frame_13)

        self.frame_14 = QFrame(self.page_11)
        self.frame_14.setObjectName(u"frame_14")
        self.frame_14.setFrameShape(QFrame.Shape.NoFrame)
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
        icon18 = QIcon()
        icon18.addFile(u":/Blue Icons/Blue Icons/crop.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BThreshold.setIcon(icon18)
        self.BThreshold.setIconSize(QSize(32, 32))
        self.BThreshold.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)

        self.horizontalLayout_16.addWidget(self.BThreshold)


        self.verticalLayout_45.addWidget(self.frame_14)

        self.BOpenParaviewSGH = QPushButton(self.page_11)
        self.BOpenParaviewSGH.setObjectName(u"BOpenParaviewSGH")

        self.verticalLayout_45.addWidget(self.BOpenParaviewSGH)

        self.verticalSpacer_18 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_45.addItem(self.verticalSpacer_18)

        self.PostprocessingOptions.addWidget(self.page_11)
        self.page_9 = QWidget()
        self.page_9.setObjectName(u"page_9")
        self.verticalLayout_47 = QVBoxLayout(self.page_9)
        self.verticalLayout_47.setObjectName(u"verticalLayout_47")
        self.PreviewResults = QFrame(self.page_9)
        self.PreviewResults.setObjectName(u"PreviewResults")
        self.PreviewResults.setFrameShape(QFrame.Shape.NoFrame)
        self.horizontalLayout_4 = QHBoxLayout(self.PreviewResults)
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.INBCFile = QComboBox(self.PreviewResults)
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.setObjectName(u"INBCFile")

        self.horizontalLayout_4.addWidget(self.INBCFile)

        self.INPreviewResults = QComboBox(self.PreviewResults)
        self.INPreviewResults.addItem("")
        self.INPreviewResults.addItem("")
        self.INPreviewResults.setObjectName(u"INPreviewResults")
        self.INPreviewResults.setFrame(True)

        self.horizontalLayout_4.addWidget(self.INPreviewResults)

        self.INResultRegion = QComboBox(self.PreviewResults)
        self.INResultRegion.addItem("")
        self.INResultRegion.addItem("")
        self.INResultRegion.addItem("")
        self.INResultRegion.setObjectName(u"INResultRegion")

        self.horizontalLayout_4.addWidget(self.INResultRegion)


        self.verticalLayout_47.addWidget(self.PreviewResults, 0, Qt.AlignmentFlag.AlignTop)

        self.frame_66 = QFrame(self.page_9)
        self.frame_66.setObjectName(u"frame_66")
        self.frame_66.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_66.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_20 = QFormLayout(self.frame_66)
        self.formLayout_20.setObjectName(u"formLayout_20")
        self.formLayout_20.setContentsMargins(0, 0, 0, 0)
        self.label_66 = QLabel(self.frame_66)
        self.label_66.setObjectName(u"label_66")

        self.formLayout_20.setWidget(0, QFormLayout.LabelRole, self.label_66)

        self.INHDeform = QDoubleSpinBox(self.frame_66)
        self.INHDeform.setObjectName(u"INHDeform")
        self.INHDeform.setButtonSymbols(QAbstractSpinBox.ButtonSymbols.UpDownArrows)
        self.INHDeform.setCorrectionMode(QAbstractSpinBox.CorrectionMode.CorrectToNearestValue)
        self.INHDeform.setDecimals(1)
        self.INHDeform.setMinimum(0.000000000000000)
        self.INHDeform.setMaximum(10000.000000000000000)
        self.INHDeform.setSingleStep(1.000000000000000)
        self.INHDeform.setValue(0.000000000000000)

        self.formLayout_20.setWidget(0, QFormLayout.FieldRole, self.INHDeform)


        self.verticalLayout_47.addWidget(self.frame_66)

        self.BOpenParaview = QPushButton(self.page_9)
        self.BOpenParaview.setObjectName(u"BOpenParaview")

        self.verticalLayout_47.addWidget(self.BOpenParaview)

        self.THomogenization = QTableWidget(self.page_9)
        if (self.THomogenization.columnCount() < 2):
            self.THomogenization.setColumnCount(2)
        __qtablewidgetitem280 = QTableWidgetItem()
        self.THomogenization.setHorizontalHeaderItem(0, __qtablewidgetitem280)
        __qtablewidgetitem281 = QTableWidgetItem()
        self.THomogenization.setHorizontalHeaderItem(1, __qtablewidgetitem281)
        if (self.THomogenization.rowCount() < 9):
            self.THomogenization.setRowCount(9)
        __qtablewidgetitem282 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(0, __qtablewidgetitem282)
        __qtablewidgetitem283 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(1, __qtablewidgetitem283)
        __qtablewidgetitem284 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(2, __qtablewidgetitem284)
        __qtablewidgetitem285 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(3, __qtablewidgetitem285)
        __qtablewidgetitem286 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(4, __qtablewidgetitem286)
        __qtablewidgetitem287 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(5, __qtablewidgetitem287)
        __qtablewidgetitem288 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(6, __qtablewidgetitem288)
        __qtablewidgetitem289 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(7, __qtablewidgetitem289)
        __qtablewidgetitem290 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(8, __qtablewidgetitem290)
        __qtablewidgetitem291 = QTableWidgetItem()
        self.THomogenization.setItem(0, 1, __qtablewidgetitem291)
        __qtablewidgetitem292 = QTableWidgetItem()
        self.THomogenization.setItem(1, 1, __qtablewidgetitem292)
        __qtablewidgetitem293 = QTableWidgetItem()
        self.THomogenization.setItem(2, 1, __qtablewidgetitem293)
        __qtablewidgetitem294 = QTableWidgetItem()
        self.THomogenization.setItem(6, 1, __qtablewidgetitem294)
        __qtablewidgetitem295 = QTableWidgetItem()
        self.THomogenization.setItem(7, 1, __qtablewidgetitem295)
        __qtablewidgetitem296 = QTableWidgetItem()
        self.THomogenization.setItem(8, 1, __qtablewidgetitem296)
        self.THomogenization.setObjectName(u"THomogenization")
        self.THomogenization.setEnabled(False)
        sizePolicy7.setHeightForWidth(self.THomogenization.sizePolicy().hasHeightForWidth())
        self.THomogenization.setSizePolicy(sizePolicy7)
        self.THomogenization.setMinimumSize(QSize(0, 300))
        self.THomogenization.setWordWrap(False)
        self.THomogenization.horizontalHeader().setMinimumSectionSize(40)
        self.THomogenization.horizontalHeader().setDefaultSectionSize(250)
        self.THomogenization.horizontalHeader().setStretchLastSection(True)
        self.THomogenization.verticalHeader().setDefaultSectionSize(30)

        self.verticalLayout_47.addWidget(self.THomogenization)

        self.frame_23 = QFrame(self.page_9)
        self.frame_23.setObjectName(u"frame_23")
        self.frame_23.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_23.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_18 = QGridLayout(self.frame_23)
        self.gridLayout_18.setObjectName(u"gridLayout_18")
        self.gridLayout_18.setContentsMargins(0, 0, 0, 0)
        self.LHomogenizationJobDir = QLabel(self.frame_23)
        self.LHomogenizationJobDir.setObjectName(u"LHomogenizationJobDir")
        self.LHomogenizationJobDir.setEnabled(False)

        self.gridLayout_18.addWidget(self.LHomogenizationJobDir, 0, 0, 1, 1)

        self.INHomogenizationJobDir = QLineEdit(self.frame_23)
        self.INHomogenizationJobDir.setObjectName(u"INHomogenizationJobDir")
        self.INHomogenizationJobDir.setEnabled(False)

        self.gridLayout_18.addWidget(self.INHomogenizationJobDir, 0, 1, 1, 1)


        self.verticalLayout_47.addWidget(self.frame_23)

        self.verticalSpacer_20 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_47.addItem(self.verticalSpacer_20)

        self.PostprocessingOptions.addWidget(self.page_9)
        self.page_25 = QWidget()
        self.page_25.setObjectName(u"page_25")
        self.verticalLayout_41 = QVBoxLayout(self.page_25)
        self.verticalLayout_41.setObjectName(u"verticalLayout_41")
        self.INBFResults = QComboBox(self.page_25)
        self.INBFResults.addItem("")
        self.INBFResults.addItem("")
        self.INBFResults.setObjectName(u"INBFResults")

        self.verticalLayout_41.addWidget(self.INBFResults)

        self.frame_67 = QFrame(self.page_25)
        self.frame_67.setObjectName(u"frame_67")
        self.frame_67.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_67.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_21 = QFormLayout(self.frame_67)
        self.formLayout_21.setObjectName(u"formLayout_21")
        self.formLayout_21.setContentsMargins(0, 0, 0, 0)
        self.label_67 = QLabel(self.frame_67)
        self.label_67.setObjectName(u"label_67")

        self.formLayout_21.setWidget(0, QFormLayout.LabelRole, self.label_67)

        self.INBFDeform = QDoubleSpinBox(self.frame_67)
        self.INBFDeform.setObjectName(u"INBFDeform")
        self.INBFDeform.setButtonSymbols(QAbstractSpinBox.ButtonSymbols.UpDownArrows)
        self.INBFDeform.setCorrectionMode(QAbstractSpinBox.CorrectionMode.CorrectToNearestValue)
        self.INBFDeform.setDecimals(1)
        self.INBFDeform.setMinimum(0.000000000000000)
        self.INBFDeform.setMaximum(10000.000000000000000)
        self.INBFDeform.setSingleStep(1.000000000000000)
        self.INBFDeform.setValue(0.000000000000000)

        self.formLayout_21.setWidget(0, QFormLayout.FieldRole, self.INBFDeform)


        self.verticalLayout_41.addWidget(self.frame_67)

        self.BBFParaview = QPushButton(self.page_25)
        self.BBFParaview.setObjectName(u"BBFParaview")

        self.verticalLayout_41.addWidget(self.BBFParaview)

        self.frame_53 = QFrame(self.page_25)
        self.frame_53.setObjectName(u"frame_53")
        self.frame_53.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_53.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_26 = QGridLayout(self.frame_53)
        self.gridLayout_26.setObjectName(u"gridLayout_26")
        self.gridLayout_26.setContentsMargins(0, 0, 0, 0)
        self.LBFJobDir = QLabel(self.frame_53)
        self.LBFJobDir.setObjectName(u"LBFJobDir")
        self.LBFJobDir.setEnabled(False)

        self.gridLayout_26.addWidget(self.LBFJobDir, 0, 0, 1, 1)

        self.INBFJobDir = QLineEdit(self.frame_53)
        self.INBFJobDir.setObjectName(u"INBFJobDir")
        self.INBFJobDir.setEnabled(False)

        self.gridLayout_26.addWidget(self.INBFJobDir, 0, 1, 1, 1)


        self.verticalLayout_41.addWidget(self.frame_53)

        self.verticalSpacer_30 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_41.addItem(self.verticalSpacer_30)

        self.PostprocessingOptions.addWidget(self.page_25)

        self.verticalLayout_46.addWidget(self.PostprocessingOptions)

        icon19 = QIcon()
        icon19.addFile(u":/Blue Icons/Blue Icons/magnify.svg", QSize(), QIcon.Normal, QIcon.Off)
        icon19.addFile(u":/Blue Icons/Blue Icons/magnify.svg", QSize(), QIcon.Selected, QIcon.On)
        self.tabWidget_2.addTab(self.Results, icon19, "")

        self.verticalLayout_18.addWidget(self.tabWidget_2)

        self.ToolWindow.addWidget(self.ResultsEVPFFTTool)
        self.ResultsSGHTool = QWidget()
        self.ResultsSGHTool.setObjectName(u"ResultsSGHTool")
        self.verticalLayout_38 = QVBoxLayout(self.ResultsSGHTool)
        self.verticalLayout_38.setObjectName(u"verticalLayout_38")
        self.ToolWindow.addWidget(self.ResultsSGHTool)

        self.horizontalLayout_10.addWidget(self.ToolWindow, 0, Qt.AlignmentFlag.AlignTop)

        self.frame_28 = QFrame(self.frame_27)
        self.frame_28.setObjectName(u"frame_28")
        self.frame_28.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_28.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_48 = QVBoxLayout(self.frame_28)
        self.verticalLayout_48.setObjectName(u"verticalLayout_48")
        self.verticalLayout_48.setContentsMargins(0, 0, 0, 0)
        self.ParaviewFrame = QFrame(self.frame_28)
        self.ParaviewFrame.setObjectName(u"ParaviewFrame")
        sizePolicy7.setHeightForWidth(self.ParaviewFrame.sizePolicy().hasHeightForWidth())
        self.ParaviewFrame.setSizePolicy(sizePolicy7)
        self.ParaviewFrame.setMinimumSize(QSize(0, 0))
        self.ParaviewFrame.setMaximumSize(QSize(16777215, 16777215))
        self.ParaviewFrame.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.ParaviewFrame.setContextMenuPolicy(Qt.ContextMenuPolicy.NoContextMenu)
        self.ParaviewFrame.setFrameShape(QFrame.Shape.NoFrame)
        self.ParaviewFrame.setLineWidth(1)
        self.verticalLayout_19 = QVBoxLayout(self.ParaviewFrame)
        self.verticalLayout_19.setSpacing(0)
        self.verticalLayout_19.setObjectName(u"verticalLayout_19")
        self.verticalLayout_19.setContentsMargins(0, 0, 0, 0)
        self.frame_75 = QFrame(self.ParaviewFrame)
        self.frame_75.setObjectName(u"frame_75")
        sizePolicy8.setHeightForWidth(self.frame_75.sizePolicy().hasHeightForWidth())
        self.frame_75.setSizePolicy(sizePolicy8)
        self.frame_75.setStyleSheet(u"background-color: rgb(82, 87, 110);")
        self.frame_75.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_75.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_29 = QHBoxLayout(self.frame_75)
        self.horizontalLayout_29.setObjectName(u"horizontalLayout_29")
        self.horizontalLayout_29.setContentsMargins(0, 0, 0, 0)
        self.BResetCamera = QPushButton(self.frame_75)
        self.BResetCamera.setObjectName(u"BResetCamera")
        icon20 = QIcon()
        icon20.addFile(u":/Blue Icons/Blue Icons/Zoom.svg", QSize(), QIcon.Normal, QIcon.Off)
        self.BResetCamera.setIcon(icon20)
        self.BResetCamera.setIconSize(QSize(32, 32))
        self.BResetCamera.setAutoDefault(False)
        self.BResetCamera.setFlat(True)

        self.horizontalLayout_29.addWidget(self.BResetCamera, 0, Qt.AlignmentFlag.AlignRight)


        self.verticalLayout_19.addWidget(self.frame_75)

        self.splitter_2 = QSplitter(self.ParaviewFrame)
        self.splitter_2.setObjectName(u"splitter_2")
        sizePolicy7.setHeightForWidth(self.splitter_2.sizePolicy().hasHeightForWidth())
        self.splitter_2.setSizePolicy(sizePolicy7)
        self.splitter_2.setMinimumSize(QSize(0, 0))
        self.splitter_2.setFrameShape(QFrame.Shape.NoFrame)
        self.splitter_2.setOrientation(Qt.Orientation.Vertical)
        self.splitter_2.setHandleWidth(6)
        self.OutputWindows = QStackedWidget(self.splitter_2)
        self.OutputWindows.setObjectName(u"OutputWindows")
        sizePolicy3.setHeightForWidth(self.OutputWindows.sizePolicy().hasHeightForWidth())
        self.OutputWindows.setSizePolicy(sizePolicy3)
        self.OutputWindows.setMinimumSize(QSize(0, 0))
        self.OutputWindows.setCursor(QCursor(Qt.OpenHandCursor))
        self.OutputWindows.setFrameShape(QFrame.Shape.NoFrame)
        self.OutputWindows.setFrameShadow(QFrame.Shadow.Plain)
        self.OutputWindows.setLineWidth(0)
        self.OutputWindows.setMidLineWidth(0)
        self.ParaviewWindow = QWidget()
        self.ParaviewWindow.setObjectName(u"ParaviewWindow")
        sizePolicy7.setHeightForWidth(self.ParaviewWindow.sizePolicy().hasHeightForWidth())
        self.ParaviewWindow.setSizePolicy(sizePolicy7)
        self.ParaviewWindow.setMinimumSize(QSize(0, 0))
        self.paraviewLayout = QVBoxLayout(self.ParaviewWindow)
        self.paraviewLayout.setObjectName(u"paraviewLayout")
        self.paraviewLayout.setContentsMargins(0, 0, 0, 0)
        self.OutputWindows.addWidget(self.ParaviewWindow)
        self.PlotWindow = QWidget()
        self.PlotWindow.setObjectName(u"PlotWindow")
        sizePolicy7.setHeightForWidth(self.PlotWindow.sizePolicy().hasHeightForWidth())
        self.PlotWindow.setSizePolicy(sizePolicy7)
        self.verticalLayout_21 = QVBoxLayout(self.PlotWindow)
        self.verticalLayout_21.setObjectName(u"verticalLayout_21")
        self.verticalLayout_21.setContentsMargins(0, 0, 0, 0)
        self.Plot = QWidget(self.PlotWindow)
        self.Plot.setObjectName(u"Plot")
        sizePolicy7.setHeightForWidth(self.Plot.sizePolicy().hasHeightForWidth())
        self.Plot.setSizePolicy(sizePolicy7)

        self.verticalLayout_21.addWidget(self.Plot)

        self.OutputWindows.addWidget(self.PlotWindow)
        self.splitter_2.addWidget(self.OutputWindows)
        self.RunOutputs = QFrame(self.splitter_2)
        self.RunOutputs.setObjectName(u"RunOutputs")
        sizePolicy3.setHeightForWidth(self.RunOutputs.sizePolicy().hasHeightForWidth())
        self.RunOutputs.setSizePolicy(sizePolicy3)
        self.RunOutputs.setMinimumSize(QSize(0, 0))
        self.RunOutputs.setMaximumSize(QSize(16777215, 250))
        self.RunOutputs.setFrameShape(QFrame.Shape.NoFrame)
        self.verticalLayout_22 = QVBoxLayout(self.RunOutputs)
        self.verticalLayout_22.setSpacing(0)
        self.verticalLayout_22.setObjectName(u"verticalLayout_22")
        self.verticalLayout_22.setContentsMargins(0, 0, 0, 0)
        self.RunOutputProgress = QProgressBar(self.RunOutputs)
        self.RunOutputProgress.setObjectName(u"RunOutputProgress")
        sizePolicy7.setHeightForWidth(self.RunOutputProgress.sizePolicy().hasHeightForWidth())
        self.RunOutputProgress.setSizePolicy(sizePolicy7)
        self.RunOutputProgress.setMinimumSize(QSize(0, 0))
        self.RunOutputProgress.setMaximumSize(QSize(16777215, 16777215))
        self.RunOutputProgress.setValue(0)

        self.verticalLayout_22.addWidget(self.RunOutputProgress)

        self.RunOutputWindow = QPlainTextEdit(self.RunOutputs)
        self.RunOutputWindow.setObjectName(u"RunOutputWindow")
        sizePolicy3.setHeightForWidth(self.RunOutputWindow.sizePolicy().hasHeightForWidth())
        self.RunOutputWindow.setSizePolicy(sizePolicy3)
        self.RunOutputWindow.setFrameShape(QFrame.Shape.NoFrame)
        self.RunOutputWindow.setTabChangesFocus(False)
        self.RunOutputWindow.setReadOnly(True)

        self.verticalLayout_22.addWidget(self.RunOutputWindow)

        self.splitter_2.addWidget(self.RunOutputs)

        self.verticalLayout_19.addWidget(self.splitter_2)


        self.verticalLayout_48.addWidget(self.ParaviewFrame)


        self.horizontalLayout_10.addWidget(self.frame_28)


        self.verticalLayout_3.addWidget(self.frame_27)

        self.frame_51 = QFrame(self.centralwidget)
        self.frame_51.setObjectName(u"frame_51")
        self.frame_51.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_51.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_27 = QGridLayout(self.frame_51)
        self.gridLayout_27.setObjectName(u"gridLayout_27")
        self.gridLayout_27.setContentsMargins(0, 0, 0, 0)
        self.INUnits = QComboBox(self.frame_51)
        self.INUnits.addItem("")
        self.INUnits.addItem("")
        self.INUnits.setObjectName(u"INUnits")
        self.INUnits.setMaximumSize(QSize(16777215, 20))
        self.INUnits.setFont(font11)
        self.INUnits.setFrame(True)

        self.gridLayout_27.addWidget(self.INUnits, 0, 1, 1, 1)

        self.label_39 = QLabel(self.frame_51)
        self.label_39.setObjectName(u"label_39")

        self.gridLayout_27.addWidget(self.label_39, 0, 0, 1, 1)


        self.verticalLayout_3.addWidget(self.frame_51, 0, Qt.AlignmentFlag.AlignRight)

        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1326, 24))
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setObjectName(u"menuHelp")
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        MainWindow.setMenuBar(self.menubar)
        QWidget.setTabOrder(self.INSphereri, self.INLengthYR2D)
        QWidget.setTabOrder(self.INLengthYR2D, self.TAnisotropic)
        QWidget.setTabOrder(self.TAnisotropic, self.INSpherero)
        QWidget.setTabOrder(self.INSpherero, self.INSphereox)
        QWidget.setTabOrder(self.INSphereox, self.INSphereoy)
        QWidget.setTabOrder(self.INSphereoy, self.INSphereoz)
        QWidget.setTabOrder(self.INSphereoz, self.INCylinderri)
        QWidget.setTabOrder(self.INCylinderri, self.INCylinderro)
        QWidget.setTabOrder(self.INCylinderro, self.BGenerateBasicGeometry)
        QWidget.setTabOrder(self.BGenerateBasicGeometry, self.TBasicGeometries)
        QWidget.setTabOrder(self.TBasicGeometries, self.BDeleteBasicGeometry)
        QWidget.setTabOrder(self.BDeleteBasicGeometry, self.INOriginXR3D)
        QWidget.setTabOrder(self.INOriginXR3D, self.INOriginYR3D)
        QWidget.setTabOrder(self.INOriginYR3D, self.INOriginZR3D)
        QWidget.setTabOrder(self.INOriginZR3D, self.INLengthXR3D)
        QWidget.setTabOrder(self.INLengthXR3D, self.INLengthYR3D)
        QWidget.setTabOrder(self.INLengthYR3D, self.INLengthZR3D)
        QWidget.setTabOrder(self.INLengthZR3D, self.INElementsXR3D)
        QWidget.setTabOrder(self.INElementsXR3D, self.INElementsYR3D)
        QWidget.setTabOrder(self.INElementsYR3D, self.INElementsZR3D)
        QWidget.setTabOrder(self.INElementsZR3D, self.INOriginXR2D)
        QWidget.setTabOrder(self.INOriginXR2D, self.INOriginYR2D)
        QWidget.setTabOrder(self.INOriginYR2D, self.INLengthXR2D)
        QWidget.setTabOrder(self.INLengthXR2D, self.INElementsXR2D)
        QWidget.setTabOrder(self.INElementsXR2D, self.INMaterialName)
        QWidget.setTabOrder(self.INMaterialName, self.INElementsYR2D)
        QWidget.setTabOrder(self.INElementsYR2D, self.INInnerRadiusC2D)
        QWidget.setTabOrder(self.INInnerRadiusC2D, self.INOriginXC2D)
        QWidget.setTabOrder(self.INOriginXC2D, self.INOriginYC2D)
        QWidget.setTabOrder(self.INOriginYC2D, self.INLengthOutRadC2D)
        QWidget.setTabOrder(self.INLengthOutRadC2D, self.INLengthThetaC2D)
        QWidget.setTabOrder(self.INLengthThetaC2D, self.INElementsRadialC2D)
        QWidget.setTabOrder(self.INElementsRadialC2D, self.INElementsArcC2D)
        QWidget.setTabOrder(self.INElementsArcC2D, self.INInnerRadiusC3D)
        QWidget.setTabOrder(self.INInnerRadiusC3D, self.INOriginXC3D)
        QWidget.setTabOrder(self.INOriginXC3D, self.INOriginYC3D)
        QWidget.setTabOrder(self.INOriginYC3D, self.INOriginZC3D)
        QWidget.setTabOrder(self.INOriginZC3D, self.INLengthOutRadC3D)
        QWidget.setTabOrder(self.INLengthOutRadC3D, self.INLengthThetaC3D)
        QWidget.setTabOrder(self.INLengthThetaC3D, self.INLengthZC3D)
        QWidget.setTabOrder(self.INLengthZC3D, self.INElementsRadC3D)
        QWidget.setTabOrder(self.INElementsRadC3D, self.INElementsArcC3D)
        QWidget.setTabOrder(self.INElementsArcC3D, self.INElementsZC3D)
        QWidget.setTabOrder(self.INElementsZC3D, self.INSelectSolverSettings)
        QWidget.setTabOrder(self.INSelectSolverSettings, self.INTime)
        QWidget.setTabOrder(self.INTime, self.INMindt)
        QWidget.setTabOrder(self.INMindt, self.INMaxdt)
        QWidget.setTabOrder(self.INMaxdt, self.INInitialdt)
        QWidget.setTabOrder(self.INInitialdt, self.INmaxcycles)
        QWidget.setTabOrder(self.INmaxcycles, self.INGraphicsOutput)
        QWidget.setTabOrder(self.INGraphicsOutput, self.RunOutputWindow)
        QWidget.setTabOrder(self.RunOutputWindow, self.tabWidget_2)
        QWidget.setTabOrder(self.tabWidget_2, self.INSelectPostprocessing)
        QWidget.setTabOrder(self.INSelectPostprocessing, self.INOuputVarSGH)
        QWidget.setTabOrder(self.INOuputVarSGH, self.BPreviewResultsSGH)
        QWidget.setTabOrder(self.BPreviewResultsSGH, self.BFirstFrame)
        QWidget.setTabOrder(self.BFirstFrame, self.BPreviousFrame)
        QWidget.setTabOrder(self.BPreviousFrame, self.BNextFrame)
        QWidget.setTabOrder(self.BNextFrame, self.BLastFrame)
        QWidget.setTabOrder(self.BLastFrame, self.INThreshold)
        QWidget.setTabOrder(self.INThreshold, self.BThreshold)
        QWidget.setTabOrder(self.BThreshold, self.BOpenParaviewSGH)
        QWidget.setTabOrder(self.BOpenParaviewSGH, self.INPreviewResults)
        QWidget.setTabOrder(self.INPreviewResults, self.INResultRegion)
        QWidget.setTabOrder(self.INResultRegion, self.BOpenParaview)
        QWidget.setTabOrder(self.BOpenParaview, self.THomogenization)
        QWidget.setTabOrder(self.THomogenization, self.INPipelineSelection)
        QWidget.setTabOrder(self.INPipelineSelection, self.INSelectGeometry)
        QWidget.setTabOrder(self.INSelectGeometry, self.INSelectGeometryImport)
        QWidget.setTabOrder(self.INSelectGeometryImport, self.INPartName)
        QWidget.setTabOrder(self.INPartName, self.INNumberOfVoxelsX)
        QWidget.setTabOrder(self.INNumberOfVoxelsX, self.INNumberOfVoxelsY)
        QWidget.setTabOrder(self.INNumberOfVoxelsY, self.INNumberOfVoxelsZ)
        QWidget.setTabOrder(self.INNumberOfVoxelsZ, self.BStlDimensions)
        QWidget.setTabOrder(self.BStlDimensions, self.BCustomDimensions)
        QWidget.setTabOrder(self.BCustomDimensions, self.INOriginX)
        QWidget.setTabOrder(self.INOriginX, self.INOriginY)
        QWidget.setTabOrder(self.INOriginY, self.INOriginZ)
        QWidget.setTabOrder(self.INOriginZ, self.INLengthX)
        QWidget.setTabOrder(self.INLengthX, self.INLengthY)
        QWidget.setTabOrder(self.INLengthY, self.INLengthZ)
        QWidget.setTabOrder(self.INLengthZ, self.INDirectory)
        QWidget.setTabOrder(self.INDirectory, self.INFileFormat)
        QWidget.setTabOrder(self.INFileFormat, self.BImageToVTK)
        QWidget.setTabOrder(self.BImageToVTK, self.BTiffToStl)
        QWidget.setTabOrder(self.BTiffToStl, self.INSelectColorBy)
        QWidget.setTabOrder(self.INSelectColorBy, self.INBoxx1)
        QWidget.setTabOrder(self.INBoxx1, self.INBoxx2)
        QWidget.setTabOrder(self.INBoxx2, self.INBoxy1)
        QWidget.setTabOrder(self.INBoxy1, self.INBoxy2)
        QWidget.setTabOrder(self.INBoxy2, self.INBoxz1)
        QWidget.setTabOrder(self.INBoxz1, self.INBoxz2)
        QWidget.setTabOrder(self.INBoxz2, self.BDeleteGeometry)
        QWidget.setTabOrder(self.BDeleteGeometry, self.INElementType)
        QWidget.setTabOrder(self.INElementType, self.INCoordinateSystem)
        QWidget.setTabOrder(self.INCoordinateSystem, self.INDimension)
        QWidget.setTabOrder(self.INDimension, self.BGenerateGlobalMesh)
        QWidget.setTabOrder(self.BGenerateGlobalMesh, self.INRunSelection)
        QWidget.setTabOrder(self.INRunSelection, self.INBCFile)
        QWidget.setTabOrder(self.INBCFile, self.DefineAssignMats)
        QWidget.setTabOrder(self.DefineAssignMats, self.INSelectDefineMaterials)
        QWidget.setTabOrder(self.INSelectDefineMaterials, self.INMaterialNameSGH)
        QWidget.setTabOrder(self.INMaterialNameSGH, self.INEOS)
        QWidget.setTabOrder(self.INEOS, self.INArtificialViscosity)
        QWidget.setTabOrder(self.INArtificialViscosity, self.BAddMaterialSGH)
        QWidget.setTabOrder(self.BAddMaterialSGH, self.TMaterialsSGH)
        QWidget.setTabOrder(self.TMaterialsSGH, self.BDeleteMaterialSGH)
        QWidget.setTabOrder(self.BDeleteMaterialSGH, self.INq1)
        QWidget.setTabOrder(self.INq1, self.INSpecificHeat)
        QWidget.setTabOrder(self.INSpecificHeat, self.INq1ex)
        QWidget.setTabOrder(self.INq1ex, self.INq2)
        QWidget.setTabOrder(self.INq2, self.INq2ex)
        QWidget.setTabOrder(self.INq2ex, self.INGamma)
        QWidget.setTabOrder(self.INGamma, self.INMinSound)
        QWidget.setTabOrder(self.INMinSound, self.INSolidGas)
        QWidget.setTabOrder(self.INSolidGas, self.INMaterialType)
        QWidget.setTabOrder(self.INMaterialType, self.INEx)
        QWidget.setTabOrder(self.INEx, self.INEy)
        QWidget.setTabOrder(self.INEy, self.INEz)
        QWidget.setTabOrder(self.INEz, self.INNUxy)
        QWidget.setTabOrder(self.INNUxy, self.INNUxz)
        QWidget.setTabOrder(self.INNUxz, self.INNUyz)
        QWidget.setTabOrder(self.INNUyz, self.INGxy)
        QWidget.setTabOrder(self.INGxy, self.INGxz)
        QWidget.setTabOrder(self.INGxz, self.INGyz)
        QWidget.setTabOrder(self.INGyz, self.INEip)
        QWidget.setTabOrder(self.INEip, self.INYoungsModulus)
        QWidget.setTabOrder(self.INYoungsModulus, self.INPoissonsRatio)
        QWidget.setTabOrder(self.INPoissonsRatio, self.INNUip)
        QWidget.setTabOrder(self.INNUip, self.INSelectBasicGeometry)
        QWidget.setTabOrder(self.INSelectBasicGeometry, self.INGop)
        QWidget.setTabOrder(self.INGop, self.INBasicGeometryName)
        QWidget.setTabOrder(self.INBasicGeometryName, self.INNUop)
        QWidget.setTabOrder(self.INNUop, self.BVoxelizeGeometry)
        QWidget.setTabOrder(self.BVoxelizeGeometry, self.INEop)
        QWidget.setTabOrder(self.INEop, self.INIsotropicPlane)
        QWidget.setTabOrder(self.INIsotropicPlane, self.BAddMaterial)
        QWidget.setTabOrder(self.BAddMaterial, self.TMaterials)
        QWidget.setTabOrder(self.TMaterials, self.BDeleteMaterial)
        QWidget.setTabOrder(self.BDeleteMaterial, self.BRegenElasticConstants)
        QWidget.setTabOrder(self.BRegenElasticConstants, self.INSelectAssignMaterials)
        QWidget.setTabOrder(self.INSelectAssignMaterials, self.INRegion)
        QWidget.setTabOrder(self.INRegion, self.INMaterial)
        QWidget.setTabOrder(self.INMaterial, self.INDensity)
        QWidget.setTabOrder(self.INDensity, self.INSIE)
        QWidget.setTabOrder(self.INSIE, self.INVelocityX)
        QWidget.setTabOrder(self.INVelocityX, self.INVelocityY)
        QWidget.setTabOrder(self.INVelocityY, self.INVelocityZ)
        QWidget.setTabOrder(self.INVelocityZ, self.Baddmaterialassignment)
        QWidget.setTabOrder(self.Baddmaterialassignment, self.Tassignmat)
        QWidget.setTabOrder(self.Tassignmat, self.BUpMaterial)
        QWidget.setTabOrder(self.BUpMaterial, self.BDownMaterial)
        QWidget.setTabOrder(self.BDownMaterial, self.Bdeletematerialassignment)
        QWidget.setTabOrder(self.Bdeletematerialassignment, self.INSelectBoundaryConditions)
        QWidget.setTabOrder(self.INSelectBoundaryConditions, self.INBoundary)
        QWidget.setTabOrder(self.INBoundary, self.INPlanePosition)
        QWidget.setTabOrder(self.INPlanePosition, self.INType)
        QWidget.setTabOrder(self.INType, self.INVel0)
        QWidget.setTabOrder(self.INVel0, self.INVel1)
        QWidget.setTabOrder(self.INVel1, self.INVelstart)
        QWidget.setTabOrder(self.INVelstart, self.INVelend)
        QWidget.setTabOrder(self.INVelend, self.BaddBC)
        QWidget.setTabOrder(self.BaddBC, self.TBoundaryConditions)
        QWidget.setTabOrder(self.TBoundaryConditions, self.BdeleteBC)
        QWidget.setTabOrder(self.BdeleteBC, self.TMaterialAssignment)
        QWidget.setTabOrder(self.TMaterialAssignment, self.BAddMaterialAssignment)
        QWidget.setTabOrder(self.BAddMaterialAssignment, self.BDeleteMaterialAssignment)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.menuHelp.addAction(self.actionManual)
        self.menuFile.addAction(self.actionNew)
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionSaveAs)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionChange_Working_Directory)

        self.retranslateUi(MainWindow)

        self.NavigationMenu.setCurrentIndex(0)
        self.ToolWindow.setCurrentIndex(0)
        self.INPipelineSelection.setCurrentIndex(0)
        self.PipelineExtras.setCurrentIndex(0)
        self.pushButton_11.setDefault(False)
        self.INSelectGeometryImport.setCurrentIndex(0)
        self.GeometryOptions.setCurrentIndex(0)
        self.BasicGeometries.setCurrentIndex(0)
        self.MeshInputs2.setCurrentIndex(0)
        self.DefineAssignMats.setCurrentIndex(0)
        self.DefineMaterialsOptions.setCurrentIndex(0)
        self.MaterialMenu.setCurrentIndex(0)
        self.MaterialTypeTool.setCurrentIndex(0)
        self.MaterialMenu_2.setCurrentIndex(0)
        self.MaterialTypeTool_2.setCurrentIndex(0)
        self.EnablePlasticity.setCurrentIndex(0)
        self.PlasticProperties.setCurrentIndex(0)
        self.pushButton_13.setDefault(False)
        self.SlipSystemInfo.setCurrentIndex(0)
        self.pushButton_9.setDefault(False)
        self.pushButton_6.setDefault(False)
        self.pushButton_5.setDefault(False)
        self.pushButton_7.setDefault(False)
        self.pushButton_8.setDefault(False)
        self.pushButton_2.setDefault(False)
        self.pushButton.setDefault(False)
        self.pushButton_3.setDefault(False)
        self.pushButton_4.setDefault(False)
        self.AssignMaterialsOptions.setCurrentIndex(0)
        self.BoundaryConditionsOptions.setCurrentIndex(0)
        self.SolverSettingsOptions.setCurrentIndex(0)
        self.RunOptions.setCurrentIndex(0)
        self.HomogenizationRunWrite.setCurrentIndex(0)
        self.HomogenizationBatch.setCurrentIndex(0)
        self.HomogenizationRunType.setCurrentIndex(0)
        self.BRunEVPFFT2.setDefault(True)
        self.BKillEVPFFT2.setDefault(False)
        self.BFRunWrite.setCurrentIndex(0)
        self.BFRunType.setCurrentIndex(0)
        self.BRunBulkForming.setDefault(True)
        self.tabWidget_2.setCurrentIndex(0)
        self.PostprocessingOptions.setCurrentIndex(0)
        self.BResetCamera.setDefault(False)
        self.OutputWindows.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"Fierro", None))
        self.actionManual.setText(QCoreApplication.translate("MainWindow", u"Manual", None))
        self.actionChange_Working_Directory.setText(QCoreApplication.translate("MainWindow", u"Change Working Directory", None))
        self.actionSaveAs.setText(QCoreApplication.translate("MainWindow", u"Save as...", None))
        self.actionHelp.setText(QCoreApplication.translate("MainWindow", u"Help", None))
        self.actionOpen.setText(QCoreApplication.translate("MainWindow", u"Open", None))
        self.actionNew.setText(QCoreApplication.translate("MainWindow", u"New...", None))
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Title), "")
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Pipeline), QCoreApplication.translate("MainWindow", u"Pipeline", None))
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Geometry), QCoreApplication.translate("MainWindow", u"Geometry", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.Geometry), QCoreApplication.translate("MainWindow", u"Geometry", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(whatsthis)
        self.NavigationMenu.setTabWhatsThis(self.NavigationMenu.indexOf(self.Geometry), QCoreApplication.translate("MainWindow", u"Geometry", None))
#endif // QT_CONFIG(whatsthis)
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Mesh), QCoreApplication.translate("MainWindow", u"Mesh", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.Mesh), QCoreApplication.translate("MainWindow", u"Mesh", None))
#endif // QT_CONFIG(tooltip)
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Materials), QCoreApplication.translate("MainWindow", u"Materials", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.Materials), QCoreApplication.translate("MainWindow", u"Materials", None))
#endif // QT_CONFIG(tooltip)
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.BoundaryConditions), QCoreApplication.translate("MainWindow", u"Boundary Conditions", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.BoundaryConditions), QCoreApplication.translate("MainWindow", u"Boundary Conditions", None))
#endif // QT_CONFIG(tooltip)
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Solver), QCoreApplication.translate("MainWindow", u"Solver Settings", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.Solver), QCoreApplication.translate("MainWindow", u"Solver", None))
#endif // QT_CONFIG(tooltip)
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Run), QCoreApplication.translate("MainWindow", u"Run", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.Run), QCoreApplication.translate("MainWindow", u"Run", None))
#endif // QT_CONFIG(tooltip)
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Postprocessing), QCoreApplication.translate("MainWindow", u"Postprocessing", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.Postprocessing), QCoreApplication.translate("MainWindow", u"Postprocessing", None))
#endif // QT_CONFIG(tooltip)
        self.LosAlamosLogo.setText("")
        self.EVPFFTLogo.setText("")
        self.LAdditionalSoftware.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; color:#000000;\">Additional Software Packages</span></p></body></html>", None))
        self.MatarLogo.setText("")
        self.ParaviewLogo.setText("")
        self.LPipeline.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/Clipboard.svg\" width=\"40\" height=\"40\"/></p><p align=\"center\"><span style=\" font-weight:700;\">Select Solver Pipeline</span></p></body></html>", None))
        self.LPipelineSelection.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Selection:</p></body></html>", None))
        self.INPipelineSelection.setItemText(0, "")
        self.INPipelineSelection.setItemText(1, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INPipelineSelection.setItemText(2, QCoreApplication.translate("MainWindow", u"Homogenization", None))
        self.INPipelineSelection.setItemText(3, QCoreApplication.translate("MainWindow", u"Bulk Deformation", None))

        self.BLegacyEVPFFT.setText(QCoreApplication.translate("MainWindow", u"Upload Legacy Input File (.txt)", None))
#if QT_CONFIG(tooltip)
        self.pushButton_11.setToolTip(QCoreApplication.translate("MainWindow", u"Legacy EVPFFT input file. Currently limited to a single phase.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_11.setText("")
        self.LGeometryInformation.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/Shapes.svg\" width=\"50\" height=\"40\"/></p><p align=\"center\"><span style=\" font-weight:700;\">Define Geometry</span></p></body></html>", None))
        self.INSelectGeometry.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectGeometry.setItemText(1, QCoreApplication.translate("MainWindow", u"Homogenization", None))
        self.INSelectGeometry.setItemText(2, QCoreApplication.translate("MainWindow", u"Bulk Forming", None))

        self.LPartSelection.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Selection:</p></body></html>", None))
        self.INSelectGeometryImport.setItemText(0, QCoreApplication.translate("MainWindow", u"Import Geometry (.stl, .vtk)", None))
        self.INSelectGeometryImport.setItemText(1, QCoreApplication.translate("MainWindow", u"Import Image Stack (.png, .jpg, .tif)", None))
        self.INSelectGeometryImport.setItemText(2, QCoreApplication.translate("MainWindow", u"Import Data Set (.txt)", None))
        self.INSelectGeometryImport.setItemText(3, QCoreApplication.translate("MainWindow", u"Create Basic Part", None))

        self.LPartName.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Part Name:</p></body></html>", None))
        self.INPartName.setText("")
        self.BUploadGeometryFile.setText(QCoreApplication.translate("MainWindow", u"Upload Geometry File", None))
        self.LSTLVoxelization.setText(QCoreApplication.translate("MainWindow", u"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:'.AppleSystemUIFont'; font-size:15pt; font-weight:600; font-style:normal;\">\n"
"<p align=\"center\" style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Input Geometry Settings</p></body></html>", None))
        self.LVoxelCount.setText(QCoreApplication.translate("MainWindow", u"Voxel Count", None))
        self.LNumberOfVoxelsX.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.LNumberOfVoxelsY.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.LNumberOfVoxelsZ.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.BStlDimensions.setText(QCoreApplication.translate("MainWindow", u"File Properties", None))
        self.BCustomDimensions.setText(QCoreApplication.translate("MainWindow", u"Custom Properties", None))
        self.LOriginPoint.setText(QCoreApplication.translate("MainWindow", u"Origin Point", None))
        self.LOriginX.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.LOriginZ.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.LOriginY.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.Ulength7.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.Ulength8.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.Ulength9.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.LDimensions.setText(QCoreApplication.translate("MainWindow", u"Dimensions", None))
        self.INLengthX.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.INLengthZ.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.LLengthY.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.LLengthZ.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.LLengthX.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.INLengthY.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.Ulength10.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.Ulength11.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.Ulength12.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.BVoxelizeGeometry.setText(QCoreApplication.translate("MainWindow", u"Add Geometry", None))
        self.LImageStack.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\">Image Stack Conversion</p></body></html>", None))
        self.LUploadedDirectory.setText(QCoreApplication.translate("MainWindow", u"Directory:", None))
        self.LImageFileFormat.setText(QCoreApplication.translate("MainWindow", u"Image File Format:", None))
        self.BImageToVTK.setText(QCoreApplication.translate("MainWindow", u"Voxelize", None))
        self.BTiffToStl.setText(QCoreApplication.translate("MainWindow", u"Convert to STL", None))
        self.LImageStack_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\">Data Viewer</p></body></html>", None))
        self.LColorBy.setText(QCoreApplication.translate("MainWindow", u"Color By:", None))
        self.LBasicGeometry.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Create Basic Part</span></p></body></html>", None))
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
        ___qtablewidgetitem = self.TBasicGeometries.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem1 = self.TBasicGeometries.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("MainWindow", u"Type", None));
        ___qtablewidgetitem2 = self.TBasicGeometries.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("MainWindow", u"x1", None));
        ___qtablewidgetitem3 = self.TBasicGeometries.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("MainWindow", u"x2", None));
        ___qtablewidgetitem4 = self.TBasicGeometries.horizontalHeaderItem(4)
        ___qtablewidgetitem4.setText(QCoreApplication.translate("MainWindow", u"y1", None));
        ___qtablewidgetitem5 = self.TBasicGeometries.horizontalHeaderItem(5)
        ___qtablewidgetitem5.setText(QCoreApplication.translate("MainWindow", u"y2", None));
        ___qtablewidgetitem6 = self.TBasicGeometries.horizontalHeaderItem(6)
        ___qtablewidgetitem6.setText(QCoreApplication.translate("MainWindow", u"z1", None));
        ___qtablewidgetitem7 = self.TBasicGeometries.horizontalHeaderItem(7)
        ___qtablewidgetitem7.setText(QCoreApplication.translate("MainWindow", u"z2", None));
        ___qtablewidgetitem8 = self.TBasicGeometries.horizontalHeaderItem(8)
        ___qtablewidgetitem8.setText(QCoreApplication.translate("MainWindow", u"Inner Radius", None));
        ___qtablewidgetitem9 = self.TBasicGeometries.horizontalHeaderItem(9)
        ___qtablewidgetitem9.setText(QCoreApplication.translate("MainWindow", u"Outer Radius", None));
        ___qtablewidgetitem10 = self.TBasicGeometries.horizontalHeaderItem(10)
        ___qtablewidgetitem10.setText(QCoreApplication.translate("MainWindow", u"Origin x", None));
        ___qtablewidgetitem11 = self.TBasicGeometries.horizontalHeaderItem(11)
        ___qtablewidgetitem11.setText(QCoreApplication.translate("MainWindow", u"Origin y", None));
        ___qtablewidgetitem12 = self.TBasicGeometries.horizontalHeaderItem(12)
        ___qtablewidgetitem12.setText(QCoreApplication.translate("MainWindow", u"Origin z", None));
        self.BDeleteBasicGeometry.setText(QCoreApplication.translate("MainWindow", u"Delete Geometry", None))
        self.label.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:14pt; font-weight:700;\">Geometry Settings</span></p></body></html>", None))
        self.BVTKFileProperties.setText(QCoreApplication.translate("MainWindow", u"File Properties", None))
        self.BVTKCustomProperties.setText(QCoreApplication.translate("MainWindow", u"Custom Properties", None))
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"Voxel Count", None))
        self.Lvtkvy.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.Lvtkvx.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.Lvtkvz.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-weight:700;\">Origin Point</span></p></body></html>", None))
        self.Lvoz.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.Lvoy.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.Lvox.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.Ulength1.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.Ulength2.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.Ulength3.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.label_27.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-weight:700;\">Dimensions</span></p></body></html>", None))
        self.Lvlx.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.Lvlz.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.Lvly.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.Ulength4.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.Ulength5.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.Ulength6.setText(QCoreApplication.translate("MainWindow", u"[m]", None))
        self.BAddVTKGeometry.setText(QCoreApplication.translate("MainWindow", u"Add Geometry", None))
        ___qtablewidgetitem13 = self.TParts.horizontalHeaderItem(0)
        ___qtablewidgetitem13.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem14 = self.TParts.horizontalHeaderItem(1)
        ___qtablewidgetitem14.setText(QCoreApplication.translate("MainWindow", u"Origin x", None));
        ___qtablewidgetitem15 = self.TParts.horizontalHeaderItem(2)
        ___qtablewidgetitem15.setText(QCoreApplication.translate("MainWindow", u"Origin y", None));
        ___qtablewidgetitem16 = self.TParts.horizontalHeaderItem(3)
        ___qtablewidgetitem16.setText(QCoreApplication.translate("MainWindow", u"Origin z", None));
        ___qtablewidgetitem17 = self.TParts.horizontalHeaderItem(4)
        ___qtablewidgetitem17.setText(QCoreApplication.translate("MainWindow", u"Length x", None));
        ___qtablewidgetitem18 = self.TParts.horizontalHeaderItem(5)
        ___qtablewidgetitem18.setText(QCoreApplication.translate("MainWindow", u"Length y", None));
        ___qtablewidgetitem19 = self.TParts.horizontalHeaderItem(6)
        ___qtablewidgetitem19.setText(QCoreApplication.translate("MainWindow", u"Length z", None));
        ___qtablewidgetitem20 = self.TParts.horizontalHeaderItem(7)
        ___qtablewidgetitem20.setText(QCoreApplication.translate("MainWindow", u"Vox x", None));
        ___qtablewidgetitem21 = self.TParts.horizontalHeaderItem(8)
        ___qtablewidgetitem21.setText(QCoreApplication.translate("MainWindow", u"Vox y", None));
        ___qtablewidgetitem22 = self.TParts.horizontalHeaderItem(9)
        ___qtablewidgetitem22.setText(QCoreApplication.translate("MainWindow", u"Vox z", None));
        ___qtablewidgetitem23 = self.TParts.horizontalHeaderItem(10)
        ___qtablewidgetitem23.setText(QCoreApplication.translate("MainWindow", u"File Path", None));
        self.BDeleteGeometry.setText(QCoreApplication.translate("MainWindow", u"Delete Geometry", None))
        self.LDefineGlobalMesh.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/mesh.svg\" width=\"40\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Generate Mesh</span></p></body></html>", None))
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
        self.LDefineMaterials.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/mine.svg\" width=\"45\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Materials</span></p></body></html>", None))
        self.INSelectDefineMaterials.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectDefineMaterials.setItemText(1, QCoreApplication.translate("MainWindow", u"Homogenization", None))
        self.INSelectDefineMaterials.setItemText(2, QCoreApplication.translate("MainWindow", u"Bulk Forming", None))

        self.LMaterialNameSGH.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Name:</p></body></html>", None))
        self.INMaterialNameSGH.setText("")
        self.LEOS.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Equation Of State Model:</p></body></html>", None))
        self.INEOS.setItemText(0, QCoreApplication.translate("MainWindow", u"ideal_gas", None))

        self.LArtificialViscosity.setText(QCoreApplication.translate("MainWindow", u"Artificial Viscosity:", None))
        self.INArtificialViscosity.setItemText(0, QCoreApplication.translate("MainWindow", u"default", None))
        self.INArtificialViscosity.setItemText(1, QCoreApplication.translate("MainWindow", u"custom", None))

        self.BAddMaterialSGH.setText(QCoreApplication.translate("MainWindow", u"Add Material", None))
        ___qtablewidgetitem24 = self.TMaterialsSGH.horizontalHeaderItem(0)
        ___qtablewidgetitem24.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem25 = self.TMaterialsSGH.horizontalHeaderItem(1)
        ___qtablewidgetitem25.setText(QCoreApplication.translate("MainWindow", u"EOS", None));
        ___qtablewidgetitem26 = self.TMaterialsSGH.horizontalHeaderItem(2)
        ___qtablewidgetitem26.setText(QCoreApplication.translate("MainWindow", u"q1", None));
        ___qtablewidgetitem27 = self.TMaterialsSGH.horizontalHeaderItem(3)
        ___qtablewidgetitem27.setText(QCoreApplication.translate("MainWindow", u"q2", None));
        ___qtablewidgetitem28 = self.TMaterialsSGH.horizontalHeaderItem(4)
        ___qtablewidgetitem28.setText(QCoreApplication.translate("MainWindow", u"q1ex", None));
        ___qtablewidgetitem29 = self.TMaterialsSGH.horizontalHeaderItem(5)
        ___qtablewidgetitem29.setText(QCoreApplication.translate("MainWindow", u"q2ex", None));
        ___qtablewidgetitem30 = self.TMaterialsSGH.horizontalHeaderItem(6)
        ___qtablewidgetitem30.setText(QCoreApplication.translate("MainWindow", u"gamma", None));
        ___qtablewidgetitem31 = self.TMaterialsSGH.horizontalHeaderItem(7)
        ___qtablewidgetitem31.setText(QCoreApplication.translate("MainWindow", u"min sound speed", None));
        ___qtablewidgetitem32 = self.TMaterialsSGH.horizontalHeaderItem(8)
        ___qtablewidgetitem32.setText(QCoreApplication.translate("MainWindow", u"specific heat", None));
        self.BDeleteMaterialSGH.setText(QCoreApplication.translate("MainWindow", u"Delete Material", None))
        self.label_24.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Custom Artificial Viscosity Variables</span></p></body></html>", None))
        self.INq1.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.INSpecificHeat.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.Lq1ex.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Accoustic Coefficient in Expansion (q1ex):</p></body></html>", None))
        self.LGamma.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Gamma:</p></body></html>", None))
        self.INq1ex.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.INq2.setText(QCoreApplication.translate("MainWindow", u"1.33333", None))
        self.INq2ex.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.LMinSound.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Minimum Sound Speed:</p></body></html>", None))
        self.Lq2ex.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Linear Slope in Expansion (q2ex):</p></body></html>", None))
        self.Lq1.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Accoustic Coefficient in Compression (q1):</p></body></html>", None))
        self.LSpecificHeat.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Specific Heat:</p></body></html>", None))
        self.INGamma.setText(QCoreApplication.translate("MainWindow", u"1.66667", None))
        self.Lq2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Linear Slope in Compression (q2):</p></body></html>", None))
        self.INMinSound.setText(QCoreApplication.translate("MainWindow", u"1E-14", None))
        self.LMaterialName_2.setText(QCoreApplication.translate("MainWindow", u"Name:", None))
        self.LType.setText(QCoreApplication.translate("MainWindow", u"Type:", None))
        self.INSolidGas.setItemText(0, QCoreApplication.translate("MainWindow", u"Solid", None))
        self.INSolidGas.setItemText(1, QCoreApplication.translate("MainWindow", u"Gas", None))

        self.INMaterialType.setItemText(0, QCoreApplication.translate("MainWindow", u"Isotropic", None))
        self.INMaterialType.setItemText(1, QCoreApplication.translate("MainWindow", u"Transversely Isotropic", None))
        self.INMaterialType.setItemText(2, QCoreApplication.translate("MainWindow", u"Orthotropic", None))
        self.INMaterialType.setItemText(3, QCoreApplication.translate("MainWindow", u"Anisotropic", None))

        self.LPoissonsRatio.setText(QCoreApplication.translate("MainWindow", u"nu: ", None))
        self.UPressure1.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.LYoungsModulus.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E: </p></body></html>", None))
        self.LIsotropicPlane.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Isotropic Plane: </p></body></html>", None))
        self.INIsotropicPlane.setItemText(0, QCoreApplication.translate("MainWindow", u"x-y plane", None))
        self.INIsotropicPlane.setItemText(1, QCoreApplication.translate("MainWindow", u"x-z plane", None))
        self.INIsotropicPlane.setItemText(2, QCoreApplication.translate("MainWindow", u"y-z plane", None))

        self.LInPlane.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" text-decoration: underline;\">In-Plane</span></p></body></html>", None))
        self.LNUip.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">nu<span style=\" vertical-align:sub;\">IP</span>: </p></body></html>", None))
        self.LEip.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">E<span style=\" vertical-align:sub;\">IP</span>: </p></body></html>", None))
        self.UPressure2.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.LOutOfPlane.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" text-decoration: underline;\">Out-Of-Plane</span></p></body></html>", None))
        self.LEop.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">E<span style=\" vertical-align:sub;\">OP</span>: </p></body></html>", None))
        self.LNUop.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">nu<span style=\" vertical-align:sub;\">OP</span>: </p></body></html>", None))
        self.LGop.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">G<span style=\" vertical-align:sub;\">OP</span>: </p></body></html>", None))
        self.UPressure3.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.UPressure4.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"Define The Stiffness Matrix [C]:", None))

        __sortingEnabled = self.TAnisotropic.isSortingEnabled()
        self.TAnisotropic.setSortingEnabled(False)
        self.TAnisotropic.setSortingEnabled(__sortingEnabled)

        self.LGxy.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>G<span style=\" vertical-align:sub;\">xy</span>: </p></body></html>", None))
        self.UPressure5.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.LEx.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E<span style=\" vertical-align:sub;\">x</span>: </p></body></html>", None))
        self.LNUxy.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>nu<span style=\" vertical-align:sub;\">xy</span>: </p></body></html>", None))
        self.LEy.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E<span style=\" vertical-align:sub;\">y</span>: </p></body></html>", None))
        self.LNUxz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>nu<span style=\" vertical-align:sub;\">xz</span>: </p></body></html>", None))
        self.LGyz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>G<span style=\" vertical-align:sub;\">yz</span>: </p></body></html>", None))
        self.LNUyz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>nu<span style=\" vertical-align:sub;\">yz</span>: </p></body></html>", None))
        self.LEz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E<span style=\" vertical-align:sub;\">z</span>: </p></body></html>", None))
        self.LGxz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>G<span style=\" vertical-align:sub;\">xz</span>: </p></body></html>", None))
        self.UPressure6.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.MaterialMenu.setTabText(self.MaterialMenu.indexOf(self.Elastic), QCoreApplication.translate("MainWindow", u"Elastic", None))
        self.BAddMaterial.setText(QCoreApplication.translate("MainWindow", u"Add", None))
        ___qtablewidgetitem33 = self.TMaterials.horizontalHeaderItem(0)
        ___qtablewidgetitem33.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem34 = self.TMaterials.horizontalHeaderItem(1)
        ___qtablewidgetitem34.setText(QCoreApplication.translate("MainWindow", u"Type", None));
        ___qtablewidgetitem35 = self.TMaterials.horizontalHeaderItem(2)
        ___qtablewidgetitem35.setText(QCoreApplication.translate("MainWindow", u"C11", None));
        ___qtablewidgetitem36 = self.TMaterials.horizontalHeaderItem(3)
        ___qtablewidgetitem36.setText(QCoreApplication.translate("MainWindow", u"C12", None));
        ___qtablewidgetitem37 = self.TMaterials.horizontalHeaderItem(4)
        ___qtablewidgetitem37.setText(QCoreApplication.translate("MainWindow", u"C13", None));
        ___qtablewidgetitem38 = self.TMaterials.horizontalHeaderItem(5)
        ___qtablewidgetitem38.setText(QCoreApplication.translate("MainWindow", u"C14", None));
        ___qtablewidgetitem39 = self.TMaterials.horizontalHeaderItem(6)
        ___qtablewidgetitem39.setText(QCoreApplication.translate("MainWindow", u"C15", None));
        ___qtablewidgetitem40 = self.TMaterials.horizontalHeaderItem(7)
        ___qtablewidgetitem40.setText(QCoreApplication.translate("MainWindow", u"C16", None));
        ___qtablewidgetitem41 = self.TMaterials.horizontalHeaderItem(8)
        ___qtablewidgetitem41.setText(QCoreApplication.translate("MainWindow", u"C22", None));
        ___qtablewidgetitem42 = self.TMaterials.horizontalHeaderItem(9)
        ___qtablewidgetitem42.setText(QCoreApplication.translate("MainWindow", u"C23", None));
        ___qtablewidgetitem43 = self.TMaterials.horizontalHeaderItem(10)
        ___qtablewidgetitem43.setText(QCoreApplication.translate("MainWindow", u"C24", None));
        ___qtablewidgetitem44 = self.TMaterials.horizontalHeaderItem(11)
        ___qtablewidgetitem44.setText(QCoreApplication.translate("MainWindow", u"C25", None));
        ___qtablewidgetitem45 = self.TMaterials.horizontalHeaderItem(12)
        ___qtablewidgetitem45.setText(QCoreApplication.translate("MainWindow", u"C26", None));
        ___qtablewidgetitem46 = self.TMaterials.horizontalHeaderItem(13)
        ___qtablewidgetitem46.setText(QCoreApplication.translate("MainWindow", u"C33", None));
        ___qtablewidgetitem47 = self.TMaterials.horizontalHeaderItem(14)
        ___qtablewidgetitem47.setText(QCoreApplication.translate("MainWindow", u"C34", None));
        ___qtablewidgetitem48 = self.TMaterials.horizontalHeaderItem(15)
        ___qtablewidgetitem48.setText(QCoreApplication.translate("MainWindow", u"C35", None));
        ___qtablewidgetitem49 = self.TMaterials.horizontalHeaderItem(16)
        ___qtablewidgetitem49.setText(QCoreApplication.translate("MainWindow", u"C36", None));
        ___qtablewidgetitem50 = self.TMaterials.horizontalHeaderItem(17)
        ___qtablewidgetitem50.setText(QCoreApplication.translate("MainWindow", u"C44", None));
        ___qtablewidgetitem51 = self.TMaterials.horizontalHeaderItem(18)
        ___qtablewidgetitem51.setText(QCoreApplication.translate("MainWindow", u"C45", None));
        ___qtablewidgetitem52 = self.TMaterials.horizontalHeaderItem(19)
        ___qtablewidgetitem52.setText(QCoreApplication.translate("MainWindow", u"C46", None));
        ___qtablewidgetitem53 = self.TMaterials.horizontalHeaderItem(20)
        ___qtablewidgetitem53.setText(QCoreApplication.translate("MainWindow", u"C55", None));
        ___qtablewidgetitem54 = self.TMaterials.horizontalHeaderItem(21)
        ___qtablewidgetitem54.setText(QCoreApplication.translate("MainWindow", u"C56", None));
        ___qtablewidgetitem55 = self.TMaterials.horizontalHeaderItem(22)
        ___qtablewidgetitem55.setText(QCoreApplication.translate("MainWindow", u"C66", None));
        self.BDeleteMaterial.setText(QCoreApplication.translate("MainWindow", u"Delete", None))
        self.BRegenElasticConstants.setText(QCoreApplication.translate("MainWindow", u"Regenerate Elastic Constants", None))
        self.LMaterialName_3.setText(QCoreApplication.translate("MainWindow", u"Name:", None))
        self.LMaterialDefinition.setText(QCoreApplication.translate("MainWindow", u"Material Definition:", None))
        self.INMaterialDefinition.setItemText(0, QCoreApplication.translate("MainWindow", u"Custom", None))
        self.INMaterialDefinition.setItemText(1, QCoreApplication.translate("MainWindow", u"Import Elastic Parameters File", None))
        self.INMaterialDefinition.setItemText(2, QCoreApplication.translate("MainWindow", u"Import Plastic Parameters File", None))
        self.INMaterialDefinition.setItemText(3, QCoreApplication.translate("MainWindow", u"Single Crystal FCC", None))
        self.INMaterialDefinition.setItemText(4, QCoreApplication.translate("MainWindow", u"Single Crystal BCC", None))

        self.LType_2.setText(QCoreApplication.translate("MainWindow", u"Elasticity:", None))
        self.INSolidGas_2.setItemText(0, QCoreApplication.translate("MainWindow", u"Solid", None))
        self.INSolidGas_2.setItemText(1, QCoreApplication.translate("MainWindow", u"Gas", None))

        self.INMaterialType_2.setItemText(0, QCoreApplication.translate("MainWindow", u"Isotropic", None))
        self.INMaterialType_2.setItemText(1, QCoreApplication.translate("MainWindow", u"Transversely Isotropic", None))
        self.INMaterialType_2.setItemText(2, QCoreApplication.translate("MainWindow", u"Orthotropic", None))
        self.INMaterialType_2.setItemText(3, QCoreApplication.translate("MainWindow", u"Anisotropic", None))

        self.LYoungsModulus_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E: </p></body></html>", None))
        self.UPressure7.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.LPoissonsRatio_2.setText(QCoreApplication.translate("MainWindow", u"nu: ", None))
        self.LIsotropicPlane_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Isotropic Plane: </p></body></html>", None))
        self.INIsotropicPlane_2.setItemText(0, QCoreApplication.translate("MainWindow", u"x-y plane", None))
        self.INIsotropicPlane_2.setItemText(1, QCoreApplication.translate("MainWindow", u"x-z plane", None))
        self.INIsotropicPlane_2.setItemText(2, QCoreApplication.translate("MainWindow", u"y-z plane", None))

        self.LInPlane_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" text-decoration: underline;\">In-Plane</span></p></body></html>", None))
        self.LNUip_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">nu<span style=\" vertical-align:sub;\">IP</span>: </p></body></html>", None))
        self.LEip_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">E<span style=\" vertical-align:sub;\">IP</span>: </p></body></html>", None))
        self.UPressure10.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.LOutOfPlane_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" text-decoration: underline;\">Out-Of-Plane</span></p></body></html>", None))
        self.LEop_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">E<span style=\" vertical-align:sub;\">OP</span>: </p></body></html>", None))
        self.LNUop_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">nu<span style=\" vertical-align:sub;\">OP</span>: </p></body></html>", None))
        self.LGop_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">G<span style=\" vertical-align:sub;\">OP</span>: </p></body></html>", None))
        self.UPressure11.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.UPressure12.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.label_9.setText(QCoreApplication.translate("MainWindow", u"Define The Stiffness Matrix [C]:", None))

        __sortingEnabled1 = self.TAnisotropic_2.isSortingEnabled()
        self.TAnisotropic_2.setSortingEnabled(False)
        self.TAnisotropic_2.setSortingEnabled(__sortingEnabled1)

        self.LNUxz_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>nu<span style=\" vertical-align:sub;\">xz</span>: </p></body></html>", None))
        self.LGyz_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>G<span style=\" vertical-align:sub;\">yz</span>: </p></body></html>", None))
        self.LGxz_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>G<span style=\" vertical-align:sub;\">xz</span>: </p></body></html>", None))
        self.LGxy_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>G<span style=\" vertical-align:sub;\">xy</span>: </p></body></html>", None))
        self.LNUyz_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>nu<span style=\" vertical-align:sub;\">yz</span>: </p></body></html>", None))
        self.LEx_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E<span style=\" vertical-align:sub;\">x</span>: </p></body></html>", None))
        self.UPressure8.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.UPressure9.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None))
        self.LEy_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E<span style=\" vertical-align:sub;\">y</span>: </p></body></html>", None))
        self.LNUxy_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>nu<span style=\" vertical-align:sub;\">xy</span>: </p></body></html>", None))
        self.LEz_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>E<span style=\" vertical-align:sub;\">z</span>: </p></body></html>", None))
        self.MaterialMenu_2.setTabText(self.MaterialMenu_2.indexOf(self.Elastic_2), QCoreApplication.translate("MainWindow", u"Elastic", None))
        self.BEnablePlasticity.setText(QCoreApplication.translate("MainWindow", u"Enable Plasticity", None))
        self.Lc.setText(QCoreApplication.translate("MainWindow", u"c", None))
        self.La.setText(QCoreApplication.translate("MainWindow", u"a", None))
        self.Lb.setText(QCoreApplication.translate("MainWindow", u"b", None))
#if QT_CONFIG(tooltip)
        self.pushButton_13.setToolTip(QCoreApplication.translate("MainWindow", u"Lattice parameters or unit cell dimensions of a crystal structure", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_13.setText("")
        self.PlasticProperties.setTabText(self.PlasticProperties.indexOf(self.CrystalAxis), QCoreApplication.translate("MainWindow", u"Crystal Axis", None))

        __sortingEnabled2 = self.TSlipSystems.isSortingEnabled()
        self.TSlipSystems.setSortingEnabled(False)
        ___qtreewidgetitem = self.TSlipSystems.topLevelItem(0)
        ___qtreewidgetitem.setText(0, QCoreApplication.translate("MainWindow", u"FCC", None));
        ___qtreewidgetitem1 = ___qtreewidgetitem.child(0)
        ___qtreewidgetitem1.setText(0, QCoreApplication.translate("MainWindow", u"1. (111)<110>", None));
        ___qtreewidgetitem2 = self.TSlipSystems.topLevelItem(1)
        ___qtreewidgetitem2.setText(0, QCoreApplication.translate("MainWindow", u"BCC", None));
        ___qtreewidgetitem3 = ___qtreewidgetitem2.child(0)
        ___qtreewidgetitem3.setText(0, QCoreApplication.translate("MainWindow", u"2. (110)<111>", None));
        ___qtreewidgetitem4 = ___qtreewidgetitem2.child(1)
        ___qtreewidgetitem4.setText(0, QCoreApplication.translate("MainWindow", u"3. (112)<111>", None));
        ___qtreewidgetitem5 = ___qtreewidgetitem2.child(2)
        ___qtreewidgetitem5.setText(0, QCoreApplication.translate("MainWindow", u"4. (123)<111>", None));
        ___qtreewidgetitem6 = self.TSlipSystems.topLevelItem(2)
        ___qtreewidgetitem6.setText(0, QCoreApplication.translate("MainWindow", u"CUSTOM", None));
        self.TSlipSystems.setSortingEnabled(__sortingEnabled2)

        self.BSlipSystemDetails.setText(QCoreApplication.translate("MainWindow", u"Details", None))
        self.BCustomSlipSystem.setText(QCoreApplication.translate("MainWindow", u"Custom", None))
        self.BAddSlipSystem.setText(QCoreApplication.translate("MainWindow", u"Add", None))
        self.BRemoveSlipSystem.setText(QCoreApplication.translate("MainWindow", u"Remove", None))
        self.PlasticProperties.setTabText(self.PlasticProperties.indexOf(self.SlipSystems), QCoreApplication.translate("MainWindow", u"Slip Systems", None))
        self.label_42.setText(QCoreApplication.translate("MainWindow", u"FCC (111)<110>", None))
        ___qtablewidgetitem56 = self.T1.horizontalHeaderItem(0)
        ___qtablewidgetitem56.setText(QCoreApplication.translate("MainWindow", u"(Slip Plane)", None));
        ___qtablewidgetitem57 = self.T1.horizontalHeaderItem(1)
        ___qtablewidgetitem57.setText(QCoreApplication.translate("MainWindow", u"<Slip Direction>", None));

        __sortingEnabled3 = self.T1.isSortingEnabled()
        self.T1.setSortingEnabled(False)
        ___qtablewidgetitem58 = self.T1.item(0, 0)
        ___qtablewidgetitem58.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem59 = self.T1.item(0, 1)
        ___qtablewidgetitem59.setText(QCoreApplication.translate("MainWindow", u"0, 1, 1", None));
        ___qtablewidgetitem60 = self.T1.item(1, 0)
        ___qtablewidgetitem60.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem61 = self.T1.item(1, 1)
        ___qtablewidgetitem61.setText(QCoreApplication.translate("MainWindow", u"1, 0, 1", None));
        ___qtablewidgetitem62 = self.T1.item(2, 0)
        ___qtablewidgetitem62.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem63 = self.T1.item(2, 1)
        ___qtablewidgetitem63.setText(QCoreApplication.translate("MainWindow", u"1, -1, 0", None));
        ___qtablewidgetitem64 = self.T1.item(3, 0)
        ___qtablewidgetitem64.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem65 = self.T1.item(3, 1)
        ___qtablewidgetitem65.setText(QCoreApplication.translate("MainWindow", u"0, 1, -1", None));
        ___qtablewidgetitem66 = self.T1.item(4, 0)
        ___qtablewidgetitem66.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem67 = self.T1.item(4, 1)
        ___qtablewidgetitem67.setText(QCoreApplication.translate("MainWindow", u"1, 0, 1", None));
        ___qtablewidgetitem68 = self.T1.item(5, 0)
        ___qtablewidgetitem68.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem69 = self.T1.item(5, 1)
        ___qtablewidgetitem69.setText(QCoreApplication.translate("MainWindow", u"1, 1, 0", None));
        ___qtablewidgetitem70 = self.T1.item(6, 0)
        ___qtablewidgetitem70.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem71 = self.T1.item(6, 1)
        ___qtablewidgetitem71.setText(QCoreApplication.translate("MainWindow", u"0, 1, 1", None));
        ___qtablewidgetitem72 = self.T1.item(7, 0)
        ___qtablewidgetitem72.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem73 = self.T1.item(7, 1)
        ___qtablewidgetitem73.setText(QCoreApplication.translate("MainWindow", u"1, 0, -1", None));
        ___qtablewidgetitem74 = self.T1.item(8, 0)
        ___qtablewidgetitem74.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem75 = self.T1.item(8, 1)
        ___qtablewidgetitem75.setText(QCoreApplication.translate("MainWindow", u"1, 1, 0", None));
        ___qtablewidgetitem76 = self.T1.item(9, 0)
        ___qtablewidgetitem76.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem77 = self.T1.item(9, 1)
        ___qtablewidgetitem77.setText(QCoreApplication.translate("MainWindow", u"0, 1, -1", None));
        ___qtablewidgetitem78 = self.T1.item(10, 0)
        ___qtablewidgetitem78.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem79 = self.T1.item(10, 1)
        ___qtablewidgetitem79.setText(QCoreApplication.translate("MainWindow", u"1, 0, -1", None));
        ___qtablewidgetitem80 = self.T1.item(11, 0)
        ___qtablewidgetitem80.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem81 = self.T1.item(11, 1)
        ___qtablewidgetitem81.setText(QCoreApplication.translate("MainWindow", u"1, -1, 0", None));
        self.T1.setSortingEnabled(__sortingEnabled3)

        self.label_44.setText(QCoreApplication.translate("MainWindow", u"BCC (110)<111>", None))
        ___qtablewidgetitem82 = self.T2.horizontalHeaderItem(0)
        ___qtablewidgetitem82.setText(QCoreApplication.translate("MainWindow", u"(Slip Plane)", None));
        ___qtablewidgetitem83 = self.T2.horizontalHeaderItem(1)
        ___qtablewidgetitem83.setText(QCoreApplication.translate("MainWindow", u"<Slip Direction>", None));

        __sortingEnabled4 = self.T2.isSortingEnabled()
        self.T2.setSortingEnabled(False)
        ___qtablewidgetitem84 = self.T2.item(0, 0)
        ___qtablewidgetitem84.setText(QCoreApplication.translate("MainWindow", u"0, 1, 1", None));
        ___qtablewidgetitem85 = self.T2.item(0, 1)
        ___qtablewidgetitem85.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem86 = self.T2.item(1, 0)
        ___qtablewidgetitem86.setText(QCoreApplication.translate("MainWindow", u"1, 0, 1", None));
        ___qtablewidgetitem87 = self.T2.item(1, 1)
        ___qtablewidgetitem87.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem88 = self.T2.item(2, 0)
        ___qtablewidgetitem88.setText(QCoreApplication.translate("MainWindow", u"1, -1, 0", None));
        ___qtablewidgetitem89 = self.T2.item(2, 1)
        ___qtablewidgetitem89.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem90 = self.T2.item(3, 0)
        ___qtablewidgetitem90.setText(QCoreApplication.translate("MainWindow", u"0, 1, -1", None));
        ___qtablewidgetitem91 = self.T2.item(3, 1)
        ___qtablewidgetitem91.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem92 = self.T2.item(4, 0)
        ___qtablewidgetitem92.setText(QCoreApplication.translate("MainWindow", u"1, 0, 1", None));
        ___qtablewidgetitem93 = self.T2.item(4, 1)
        ___qtablewidgetitem93.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem94 = self.T2.item(5, 0)
        ___qtablewidgetitem94.setText(QCoreApplication.translate("MainWindow", u"1, 1, 0", None));
        ___qtablewidgetitem95 = self.T2.item(5, 1)
        ___qtablewidgetitem95.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem96 = self.T2.item(6, 0)
        ___qtablewidgetitem96.setText(QCoreApplication.translate("MainWindow", u"0, 1, 1", None));
        ___qtablewidgetitem97 = self.T2.item(6, 1)
        ___qtablewidgetitem97.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem98 = self.T2.item(7, 0)
        ___qtablewidgetitem98.setText(QCoreApplication.translate("MainWindow", u"1, 0, -1", None));
        ___qtablewidgetitem99 = self.T2.item(7, 1)
        ___qtablewidgetitem99.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem100 = self.T2.item(8, 0)
        ___qtablewidgetitem100.setText(QCoreApplication.translate("MainWindow", u"1, 1, 0", None));
        ___qtablewidgetitem101 = self.T2.item(8, 1)
        ___qtablewidgetitem101.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem102 = self.T2.item(9, 0)
        ___qtablewidgetitem102.setText(QCoreApplication.translate("MainWindow", u"0, 1, -1", None));
        ___qtablewidgetitem103 = self.T2.item(9, 1)
        ___qtablewidgetitem103.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem104 = self.T2.item(10, 0)
        ___qtablewidgetitem104.setText(QCoreApplication.translate("MainWindow", u"1, 0, -1", None));
        ___qtablewidgetitem105 = self.T2.item(10, 1)
        ___qtablewidgetitem105.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem106 = self.T2.item(11, 0)
        ___qtablewidgetitem106.setText(QCoreApplication.translate("MainWindow", u"1, -1, 0", None));
        ___qtablewidgetitem107 = self.T2.item(11, 1)
        ___qtablewidgetitem107.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        self.T2.setSortingEnabled(__sortingEnabled4)

        self.label_45.setText(QCoreApplication.translate("MainWindow", u"BCC (112)<111>", None))
        ___qtablewidgetitem108 = self.T3.horizontalHeaderItem(0)
        ___qtablewidgetitem108.setText(QCoreApplication.translate("MainWindow", u"(Slip Plane)", None));
        ___qtablewidgetitem109 = self.T3.horizontalHeaderItem(1)
        ___qtablewidgetitem109.setText(QCoreApplication.translate("MainWindow", u"<Slip Direction>", None));

        __sortingEnabled5 = self.T3.isSortingEnabled()
        self.T3.setSortingEnabled(False)
        ___qtablewidgetitem110 = self.T3.item(0, 0)
        ___qtablewidgetitem110.setText(QCoreApplication.translate("MainWindow", u"-2, 1, -1", None));
        ___qtablewidgetitem111 = self.T3.item(0, 1)
        ___qtablewidgetitem111.setText(QCoreApplication.translate("MainWindow", u"-1, -1, 1", None));
        ___qtablewidgetitem112 = self.T3.item(1, 0)
        ___qtablewidgetitem112.setText(QCoreApplication.translate("MainWindow", u"1, -2, -1", None));
        ___qtablewidgetitem113 = self.T3.item(1, 1)
        ___qtablewidgetitem113.setText(QCoreApplication.translate("MainWindow", u"-1, -1, 1", None));
        ___qtablewidgetitem114 = self.T3.item(2, 0)
        ___qtablewidgetitem114.setText(QCoreApplication.translate("MainWindow", u"1, 1, 2", None));
        ___qtablewidgetitem115 = self.T3.item(2, 1)
        ___qtablewidgetitem115.setText(QCoreApplication.translate("MainWindow", u"-1, -1, 1", None));
        ___qtablewidgetitem116 = self.T3.item(3, 0)
        ___qtablewidgetitem116.setText(QCoreApplication.translate("MainWindow", u"-2, -1, -1", None));
        ___qtablewidgetitem117 = self.T3.item(3, 1)
        ___qtablewidgetitem117.setText(QCoreApplication.translate("MainWindow", u"-1, 1, 1", None));
        ___qtablewidgetitem118 = self.T3.item(4, 0)
        ___qtablewidgetitem118.setText(QCoreApplication.translate("MainWindow", u"1, 2, -1", None));
        ___qtablewidgetitem119 = self.T3.item(4, 1)
        ___qtablewidgetitem119.setText(QCoreApplication.translate("MainWindow", u"-1, 1, 1", None));
        ___qtablewidgetitem120 = self.T3.item(5, 0)
        ___qtablewidgetitem120.setText(QCoreApplication.translate("MainWindow", u"1, -1, 2", None));
        ___qtablewidgetitem121 = self.T3.item(5, 1)
        ___qtablewidgetitem121.setText(QCoreApplication.translate("MainWindow", u"-1, 1, 1", None));
        ___qtablewidgetitem122 = self.T3.item(6, 0)
        ___qtablewidgetitem122.setText(QCoreApplication.translate("MainWindow", u"2, 1, -1", None));
        ___qtablewidgetitem123 = self.T3.item(6, 1)
        ___qtablewidgetitem123.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem124 = self.T3.item(7, 0)
        ___qtablewidgetitem124.setText(QCoreApplication.translate("MainWindow", u"-1, -2, -1", None));
        ___qtablewidgetitem125 = self.T3.item(7, 1)
        ___qtablewidgetitem125.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem126 = self.T3.item(8, 0)
        ___qtablewidgetitem126.setText(QCoreApplication.translate("MainWindow", u"-1, 1, 2", None));
        ___qtablewidgetitem127 = self.T3.item(8, 1)
        ___qtablewidgetitem127.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem128 = self.T3.item(9, 0)
        ___qtablewidgetitem128.setText(QCoreApplication.translate("MainWindow", u"2, -1, -1", None));
        ___qtablewidgetitem129 = self.T3.item(9, 1)
        ___qtablewidgetitem129.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem130 = self.T3.item(10, 0)
        ___qtablewidgetitem130.setText(QCoreApplication.translate("MainWindow", u"-1, 2, -1", None));
        ___qtablewidgetitem131 = self.T3.item(10, 1)
        ___qtablewidgetitem131.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem132 = self.T3.item(11, 0)
        ___qtablewidgetitem132.setText(QCoreApplication.translate("MainWindow", u"-1, -1, 2", None));
        ___qtablewidgetitem133 = self.T3.item(11, 1)
        ___qtablewidgetitem133.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        self.T3.setSortingEnabled(__sortingEnabled5)

        self.label_47.setText(QCoreApplication.translate("MainWindow", u"BCC (123)<111>", None))
        ___qtablewidgetitem134 = self.T4.horizontalHeaderItem(0)
        ___qtablewidgetitem134.setText(QCoreApplication.translate("MainWindow", u"(Slip Plane)", None));
        ___qtablewidgetitem135 = self.T4.horizontalHeaderItem(1)
        ___qtablewidgetitem135.setText(QCoreApplication.translate("MainWindow", u"<Slip Direction>", None));

        __sortingEnabled6 = self.T4.isSortingEnabled()
        self.T4.setSortingEnabled(False)
        ___qtablewidgetitem136 = self.T4.item(0, 0)
        ___qtablewidgetitem136.setText(QCoreApplication.translate("MainWindow", u"1, 2, 3", None));
        ___qtablewidgetitem137 = self.T4.item(0, 1)
        ___qtablewidgetitem137.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem138 = self.T4.item(1, 0)
        ___qtablewidgetitem138.setText(QCoreApplication.translate("MainWindow", u"-1, 3, 2", None));
        ___qtablewidgetitem139 = self.T4.item(1, 1)
        ___qtablewidgetitem139.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem140 = self.T4.item(2, 0)
        ___qtablewidgetitem140.setText(QCoreApplication.translate("MainWindow", u"2, 1, 3", None));
        ___qtablewidgetitem141 = self.T4.item(2, 1)
        ___qtablewidgetitem141.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem142 = self.T4.item(3, 0)
        ___qtablewidgetitem142.setText(QCoreApplication.translate("MainWindow", u"-2, 3, 1", None));
        ___qtablewidgetitem143 = self.T4.item(3, 1)
        ___qtablewidgetitem143.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem144 = self.T4.item(4, 0)
        ___qtablewidgetitem144.setText(QCoreApplication.translate("MainWindow", u"3, -1, 2", None));
        ___qtablewidgetitem145 = self.T4.item(4, 1)
        ___qtablewidgetitem145.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem146 = self.T4.item(5, 0)
        ___qtablewidgetitem146.setText(QCoreApplication.translate("MainWindow", u"3, -2, 1", None));
        ___qtablewidgetitem147 = self.T4.item(5, 1)
        ___qtablewidgetitem147.setText(QCoreApplication.translate("MainWindow", u"1, 1, -1", None));
        ___qtablewidgetitem148 = self.T4.item(6, 0)
        ___qtablewidgetitem148.setText(QCoreApplication.translate("MainWindow", u"-1, 2, -3", None));
        ___qtablewidgetitem149 = self.T4.item(6, 1)
        ___qtablewidgetitem149.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem150 = self.T4.item(7, 0)
        ___qtablewidgetitem150.setText(QCoreApplication.translate("MainWindow", u"1, 3, -2", None));
        ___qtablewidgetitem151 = self.T4.item(7, 1)
        ___qtablewidgetitem151.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem152 = self.T4.item(8, 0)
        ___qtablewidgetitem152.setText(QCoreApplication.translate("MainWindow", u"2, -1, 3", None));
        ___qtablewidgetitem153 = self.T4.item(8, 1)
        ___qtablewidgetitem153.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem154 = self.T4.item(9, 0)
        ___qtablewidgetitem154.setText(QCoreApplication.translate("MainWindow", u"2, 3, -1", None));
        ___qtablewidgetitem155 = self.T4.item(9, 1)
        ___qtablewidgetitem155.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem156 = self.T4.item(10, 0)
        ___qtablewidgetitem156.setText(QCoreApplication.translate("MainWindow", u"3, 1, 2", None));
        ___qtablewidgetitem157 = self.T4.item(10, 1)
        ___qtablewidgetitem157.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem158 = self.T4.item(11, 0)
        ___qtablewidgetitem158.setText(QCoreApplication.translate("MainWindow", u"3, 2, 1", None));
        ___qtablewidgetitem159 = self.T4.item(11, 1)
        ___qtablewidgetitem159.setText(QCoreApplication.translate("MainWindow", u"1, -1, -1", None));
        ___qtablewidgetitem160 = self.T4.item(12, 0)
        ___qtablewidgetitem160.setText(QCoreApplication.translate("MainWindow", u"1, -2, -3", None));
        ___qtablewidgetitem161 = self.T4.item(12, 1)
        ___qtablewidgetitem161.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem162 = self.T4.item(13, 0)
        ___qtablewidgetitem162.setText(QCoreApplication.translate("MainWindow", u"1, 3, 2", None));
        ___qtablewidgetitem163 = self.T4.item(13, 1)
        ___qtablewidgetitem163.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem164 = self.T4.item(14, 0)
        ___qtablewidgetitem164.setText(QCoreApplication.translate("MainWindow", u"2, -1, -3", None));
        ___qtablewidgetitem165 = self.T4.item(14, 1)
        ___qtablewidgetitem165.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem166 = self.T4.item(15, 0)
        ___qtablewidgetitem166.setText(QCoreApplication.translate("MainWindow", u"2, 3, 1", None));
        ___qtablewidgetitem167 = self.T4.item(15, 1)
        ___qtablewidgetitem167.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem168 = self.T4.item(16, 0)
        ___qtablewidgetitem168.setText(QCoreApplication.translate("MainWindow", u"3, 1, -2", None));
        ___qtablewidgetitem169 = self.T4.item(16, 1)
        ___qtablewidgetitem169.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem170 = self.T4.item(17, 0)
        ___qtablewidgetitem170.setText(QCoreApplication.translate("MainWindow", u"3, 2, -1", None));
        ___qtablewidgetitem171 = self.T4.item(17, 1)
        ___qtablewidgetitem171.setText(QCoreApplication.translate("MainWindow", u"1, -1, 1", None));
        ___qtablewidgetitem172 = self.T4.item(18, 0)
        ___qtablewidgetitem172.setText(QCoreApplication.translate("MainWindow", u"1, 2, -3", None));
        ___qtablewidgetitem173 = self.T4.item(18, 1)
        ___qtablewidgetitem173.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem174 = self.T4.item(19, 0)
        ___qtablewidgetitem174.setText(QCoreApplication.translate("MainWindow", u"1, -3, 2", None));
        ___qtablewidgetitem175 = self.T4.item(19, 1)
        ___qtablewidgetitem175.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem176 = self.T4.item(20, 0)
        ___qtablewidgetitem176.setText(QCoreApplication.translate("MainWindow", u"2, 1, -3", None));
        ___qtablewidgetitem177 = self.T4.item(20, 1)
        ___qtablewidgetitem177.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem178 = self.T4.item(21, 0)
        ___qtablewidgetitem178.setText(QCoreApplication.translate("MainWindow", u"2, -3, 1", None));
        ___qtablewidgetitem179 = self.T4.item(21, 1)
        ___qtablewidgetitem179.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem180 = self.T4.item(22, 0)
        ___qtablewidgetitem180.setText(QCoreApplication.translate("MainWindow", u"-3, 1, 2", None));
        ___qtablewidgetitem181 = self.T4.item(22, 1)
        ___qtablewidgetitem181.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        ___qtablewidgetitem182 = self.T4.item(23, 0)
        ___qtablewidgetitem182.setText(QCoreApplication.translate("MainWindow", u"-3, 2, 1", None));
        ___qtablewidgetitem183 = self.T4.item(23, 1)
        ___qtablewidgetitem183.setText(QCoreApplication.translate("MainWindow", u"1, 1, 1", None));
        self.T4.setSortingEnabled(__sortingEnabled6)

        self.PlasticProperties.setTabText(self.PlasticProperties.indexOf(self.Details), QCoreApplication.translate("MainWindow", u"i", None))
        ___qtablewidgetitem184 = self.TSlipSystemParameters.horizontalHeaderItem(0)
        ___qtablewidgetitem184.setText(QCoreApplication.translate("MainWindow", u"Slip System", None));
        ___qtablewidgetitem185 = self.TSlipSystemParameters.horizontalHeaderItem(1)
        ___qtablewidgetitem185.setText(QCoreApplication.translate("MainWindow", u"nrsx", None));
        ___qtablewidgetitem186 = self.TSlipSystemParameters.horizontalHeaderItem(2)
        ___qtablewidgetitem186.setText(QCoreApplication.translate("MainWindow", u"gamd0x", None));
        ___qtablewidgetitem187 = self.TSlipSystemParameters.horizontalHeaderItem(3)
        ___qtablewidgetitem187.setText(QCoreApplication.translate("MainWindow", u"tau0xf", None));
        ___qtablewidgetitem188 = self.TSlipSystemParameters.horizontalHeaderItem(4)
        ___qtablewidgetitem188.setText(QCoreApplication.translate("MainWindow", u"tau0xb", None));
        ___qtablewidgetitem189 = self.TSlipSystemParameters.horizontalHeaderItem(5)
        ___qtablewidgetitem189.setText(QCoreApplication.translate("MainWindow", u"tau1x", None));
        ___qtablewidgetitem190 = self.TSlipSystemParameters.horizontalHeaderItem(6)
        ___qtablewidgetitem190.setText(QCoreApplication.translate("MainWindow", u"thet0", None));
        ___qtablewidgetitem191 = self.TSlipSystemParameters.horizontalHeaderItem(7)
        ___qtablewidgetitem191.setText(QCoreApplication.translate("MainWindow", u"thet1", None));
        ___qtablewidgetitem192 = self.TSlipSystemParameters.horizontalHeaderItem(8)
        ___qtablewidgetitem192.setText(QCoreApplication.translate("MainWindow", u"hselfx", None));
        ___qtablewidgetitem193 = self.TSlipSystemParameters.horizontalHeaderItem(9)
        ___qtablewidgetitem193.setText(QCoreApplication.translate("MainWindow", u"hlatex", None));
        self.LSelectSlipSystem.setText(QCoreApplication.translate("MainWindow", u"Slip System:", None))
#if QT_CONFIG(tooltip)
        self.pushButton_9.setToolTip(QCoreApplication.translate("MainWindow", u"Rate sensitivity exponent. Relates shear rate on a slip system to the resolved shear stress.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_9.setText("")
        self.Lgamd0x.setText(QCoreApplication.translate("MainWindow", u"gamd0x:", None))
        self.Lhlatex.setText(QCoreApplication.translate("MainWindow", u"hlatex:", None))
#if QT_CONFIG(tooltip)
        self.pushButton_6.setToolTip(QCoreApplication.translate("MainWindow", u"Self-hardening coefficient. How the current slip system affects the hardening of the same system.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_6.setText("")
        self.Ltau0xb.setText(QCoreApplication.translate("MainWindow", u"tau0xb:", None))
        self.Ltau0xf.setText(QCoreApplication.translate("MainWindow", u"tau0xf:", None))
        self.Lthet0.setText(QCoreApplication.translate("MainWindow", u"thet0:", None))
        self.Lhselfx.setText(QCoreApplication.translate("MainWindow", u"hselfx:", None))
        self.Ltau1x.setText(QCoreApplication.translate("MainWindow", u"tau1x:", None))
#if QT_CONFIG(tooltip)
        self.pushButton_5.setToolTip(QCoreApplication.translate("MainWindow", u"Asymptotic hardening rate. The rate of hardening with increasing strain.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_5.setText("")
#if QT_CONFIG(tooltip)
        self.pushButton_7.setToolTip(QCoreApplication.translate("MainWindow", u"Latent hardening coefficient. How the current slip system affects the hardening of other systems.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_7.setText("")
#if QT_CONFIG(tooltip)
        self.pushButton_8.setToolTip(QCoreApplication.translate("MainWindow", u"Reference shar rate. Sets the scale for the plastic strain rates in the model.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_8.setText("")
#if QT_CONFIG(tooltip)
        self.pushButton_2.setToolTip(QCoreApplication.translate("MainWindow", u"Initial critical resolved shear stress for backward slip. The stress required to initiate slip on a given slip system.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_2.setText("")
        self.Lthet1.setText(QCoreApplication.translate("MainWindow", u"thet1:", None))
#if QT_CONFIG(tooltip)
        self.pushButton.setToolTip(QCoreApplication.translate("MainWindow", u"Initial critical resolved shear stress for forward slip. The stress required to initiate slip on a given slip system.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton.setText("")
#if QT_CONFIG(tooltip)
        self.pushButton_3.setToolTip(QCoreApplication.translate("MainWindow", u"Saturation stress. The maximum value that the critical resolved shear stress can reach during hardening.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_3.setText("")
#if QT_CONFIG(tooltip)
        self.pushButton_4.setToolTip(QCoreApplication.translate("MainWindow", u"Initial hardening rate. The hardening rate at the beginning of plastic deformation.", None))
#endif // QT_CONFIG(tooltip)
        self.pushButton_4.setText("")
        self.Lnrsx.setText(QCoreApplication.translate("MainWindow", u"nrsx:", None))
        self.PlasticProperties.setTabText(self.PlasticProperties.indexOf(self.VoceParameters), QCoreApplication.translate("MainWindow", u"Voce Parameters", None))
        self.MaterialMenu_2.setTabText(self.MaterialMenu_2.indexOf(self.Plastic), QCoreApplication.translate("MainWindow", u"Plastic", None))
        self.BAddMaterial_2.setText(QCoreApplication.translate("MainWindow", u"Add Material", None))
        self.BDeleteMaterial_2.setText(QCoreApplication.translate("MainWindow", u"Delete Material", None))
        ___qtablewidgetitem194 = self.TMaterials_2.horizontalHeaderItem(0)
        ___qtablewidgetitem194.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem195 = self.TMaterials_2.horizontalHeaderItem(1)
        ___qtablewidgetitem195.setText(QCoreApplication.translate("MainWindow", u"Elasticity", None));
        ___qtablewidgetitem196 = self.TMaterials_2.horizontalHeaderItem(2)
        ___qtablewidgetitem196.setText(QCoreApplication.translate("MainWindow", u"C11", None));
        ___qtablewidgetitem197 = self.TMaterials_2.horizontalHeaderItem(3)
        ___qtablewidgetitem197.setText(QCoreApplication.translate("MainWindow", u"C12", None));
        ___qtablewidgetitem198 = self.TMaterials_2.horizontalHeaderItem(4)
        ___qtablewidgetitem198.setText(QCoreApplication.translate("MainWindow", u"C13", None));
        ___qtablewidgetitem199 = self.TMaterials_2.horizontalHeaderItem(5)
        ___qtablewidgetitem199.setText(QCoreApplication.translate("MainWindow", u"C14", None));
        ___qtablewidgetitem200 = self.TMaterials_2.horizontalHeaderItem(6)
        ___qtablewidgetitem200.setText(QCoreApplication.translate("MainWindow", u"C15", None));
        ___qtablewidgetitem201 = self.TMaterials_2.horizontalHeaderItem(7)
        ___qtablewidgetitem201.setText(QCoreApplication.translate("MainWindow", u"C16", None));
        ___qtablewidgetitem202 = self.TMaterials_2.horizontalHeaderItem(8)
        ___qtablewidgetitem202.setText(QCoreApplication.translate("MainWindow", u"C22", None));
        ___qtablewidgetitem203 = self.TMaterials_2.horizontalHeaderItem(9)
        ___qtablewidgetitem203.setText(QCoreApplication.translate("MainWindow", u"C23", None));
        ___qtablewidgetitem204 = self.TMaterials_2.horizontalHeaderItem(10)
        ___qtablewidgetitem204.setText(QCoreApplication.translate("MainWindow", u"C24", None));
        ___qtablewidgetitem205 = self.TMaterials_2.horizontalHeaderItem(11)
        ___qtablewidgetitem205.setText(QCoreApplication.translate("MainWindow", u"C25", None));
        ___qtablewidgetitem206 = self.TMaterials_2.horizontalHeaderItem(12)
        ___qtablewidgetitem206.setText(QCoreApplication.translate("MainWindow", u"C26", None));
        ___qtablewidgetitem207 = self.TMaterials_2.horizontalHeaderItem(13)
        ___qtablewidgetitem207.setText(QCoreApplication.translate("MainWindow", u"C33", None));
        ___qtablewidgetitem208 = self.TMaterials_2.horizontalHeaderItem(14)
        ___qtablewidgetitem208.setText(QCoreApplication.translate("MainWindow", u"C34", None));
        ___qtablewidgetitem209 = self.TMaterials_2.horizontalHeaderItem(15)
        ___qtablewidgetitem209.setText(QCoreApplication.translate("MainWindow", u"C35", None));
        ___qtablewidgetitem210 = self.TMaterials_2.horizontalHeaderItem(16)
        ___qtablewidgetitem210.setText(QCoreApplication.translate("MainWindow", u"C36", None));
        ___qtablewidgetitem211 = self.TMaterials_2.horizontalHeaderItem(17)
        ___qtablewidgetitem211.setText(QCoreApplication.translate("MainWindow", u"C44", None));
        ___qtablewidgetitem212 = self.TMaterials_2.horizontalHeaderItem(18)
        ___qtablewidgetitem212.setText(QCoreApplication.translate("MainWindow", u"C45", None));
        ___qtablewidgetitem213 = self.TMaterials_2.horizontalHeaderItem(19)
        ___qtablewidgetitem213.setText(QCoreApplication.translate("MainWindow", u"C46", None));
        ___qtablewidgetitem214 = self.TMaterials_2.horizontalHeaderItem(20)
        ___qtablewidgetitem214.setText(QCoreApplication.translate("MainWindow", u"C55", None));
        ___qtablewidgetitem215 = self.TMaterials_2.horizontalHeaderItem(21)
        ___qtablewidgetitem215.setText(QCoreApplication.translate("MainWindow", u"C56", None));
        ___qtablewidgetitem216 = self.TMaterials_2.horizontalHeaderItem(22)
        ___qtablewidgetitem216.setText(QCoreApplication.translate("MainWindow", u"C66", None));
        ___qtablewidgetitem217 = self.TMaterials_2.horizontalHeaderItem(23)
        ___qtablewidgetitem217.setText(QCoreApplication.translate("MainWindow", u"a", None));
        ___qtablewidgetitem218 = self.TMaterials_2.horizontalHeaderItem(24)
        ___qtablewidgetitem218.setText(QCoreApplication.translate("MainWindow", u"b", None));
        ___qtablewidgetitem219 = self.TMaterials_2.horizontalHeaderItem(25)
        ___qtablewidgetitem219.setText(QCoreApplication.translate("MainWindow", u"c", None));
        ___qtablewidgetitem220 = self.TMaterials_2.horizontalHeaderItem(26)
        ___qtablewidgetitem220.setText(QCoreApplication.translate("MainWindow", u"Slip Systems", None));
        ___qtablewidgetitem221 = self.TMaterials_2.horizontalHeaderItem(27)
        ___qtablewidgetitem221.setText(QCoreApplication.translate("MainWindow", u"nrsx", None));
        ___qtablewidgetitem222 = self.TMaterials_2.horizontalHeaderItem(28)
        ___qtablewidgetitem222.setText(QCoreApplication.translate("MainWindow", u"gamd0x", None));
        ___qtablewidgetitem223 = self.TMaterials_2.horizontalHeaderItem(29)
        ___qtablewidgetitem223.setText(QCoreApplication.translate("MainWindow", u"tau0xf", None));
        ___qtablewidgetitem224 = self.TMaterials_2.horizontalHeaderItem(30)
        ___qtablewidgetitem224.setText(QCoreApplication.translate("MainWindow", u"tau0xb", None));
        ___qtablewidgetitem225 = self.TMaterials_2.horizontalHeaderItem(31)
        ___qtablewidgetitem225.setText(QCoreApplication.translate("MainWindow", u"tau1x", None));
        ___qtablewidgetitem226 = self.TMaterials_2.horizontalHeaderItem(32)
        ___qtablewidgetitem226.setText(QCoreApplication.translate("MainWindow", u"thet0", None));
        ___qtablewidgetitem227 = self.TMaterials_2.horizontalHeaderItem(33)
        ___qtablewidgetitem227.setText(QCoreApplication.translate("MainWindow", u"thet1", None));
        ___qtablewidgetitem228 = self.TMaterials_2.horizontalHeaderItem(34)
        ___qtablewidgetitem228.setText(QCoreApplication.translate("MainWindow", u"hselfx", None));
        ___qtablewidgetitem229 = self.TMaterials_2.horizontalHeaderItem(35)
        ___qtablewidgetitem229.setText(QCoreApplication.translate("MainWindow", u"hlatex", None));
        self.DefineAssignMats.setTabText(self.DefineAssignMats.indexOf(self.DefineMaterials), QCoreApplication.translate("MainWindow", u"Define Materials", None))
        self.INSelectAssignMaterials.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectAssignMaterials.setItemText(1, QCoreApplication.translate("MainWindow", u"Homogenization", None))
        self.INSelectAssignMaterials.setItemText(2, QCoreApplication.translate("MainWindow", u"Bulk Forming", None))

        self.LMaterialName.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Region:</p></body></html>", None))
        self.LRegion_5.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Material:</p></body></html>", None))
        self.LDensity.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Density:</p></body></html>", None))
        self.INDensity.setText("")
        self.LSIE.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Specific Internal Energy:</p></body></html>", None))
        self.INSIE.setText("")
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"Velocity:", None))
        self.LVely.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Y</p></body></html>", None))
        self.LVelx.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">X</p></body></html>", None))
        self.INVelocityX.setText("")
        self.INVelocityY.setText("")
        self.LVelz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Z</p></body></html>", None))
        self.INVelocityZ.setText("")
        self.Baddmaterialassignment.setText(QCoreApplication.translate("MainWindow", u"Add Material Assignment", None))
        ___qtablewidgetitem230 = self.Tassignmat.horizontalHeaderItem(0)
        ___qtablewidgetitem230.setText(QCoreApplication.translate("MainWindow", u"Region", None));
        ___qtablewidgetitem231 = self.Tassignmat.horizontalHeaderItem(1)
        ___qtablewidgetitem231.setText(QCoreApplication.translate("MainWindow", u"Material", None));
        ___qtablewidgetitem232 = self.Tassignmat.horizontalHeaderItem(2)
        ___qtablewidgetitem232.setText(QCoreApplication.translate("MainWindow", u"Density", None));
        ___qtablewidgetitem233 = self.Tassignmat.horizontalHeaderItem(3)
        ___qtablewidgetitem233.setText(QCoreApplication.translate("MainWindow", u"Specific Internal Energy", None));
        ___qtablewidgetitem234 = self.Tassignmat.horizontalHeaderItem(4)
        ___qtablewidgetitem234.setText(QCoreApplication.translate("MainWindow", u"Velocity X", None));
        ___qtablewidgetitem235 = self.Tassignmat.horizontalHeaderItem(5)
        ___qtablewidgetitem235.setText(QCoreApplication.translate("MainWindow", u"Velocity Y", None));
        ___qtablewidgetitem236 = self.Tassignmat.horizontalHeaderItem(6)
        ___qtablewidgetitem236.setText(QCoreApplication.translate("MainWindow", u"Velocity Z", None));
        self.BUpMaterial.setText("")
        self.BDownMaterial.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.Bdeletematerialassignment.setText(QCoreApplication.translate("MainWindow", u"Delete Material Assignment", None))
        self.LMaterialName_4.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Region:</p></body></html>", None))
        self.LRegion_6.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Material:</p></body></html>", None))
        self.BAddMaterialAssignment.setText(QCoreApplication.translate("MainWindow", u"Add Material Assignment", None))
        ___qtablewidgetitem237 = self.TMaterialAssignment.horizontalHeaderItem(0)
        ___qtablewidgetitem237.setText(QCoreApplication.translate("MainWindow", u"Region", None));
        ___qtablewidgetitem238 = self.TMaterialAssignment.horizontalHeaderItem(1)
        ___qtablewidgetitem238.setText(QCoreApplication.translate("MainWindow", u"Material", None));
        self.BDeleteMaterialAssignment.setText(QCoreApplication.translate("MainWindow", u"Delete Material Assignment", None))
        self.LMaterialName_5.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Region:</p></body></html>", None))
        self.LRegion_7.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Material:</p></body></html>", None))
        self.BAddMaterialAssignment_2.setText(QCoreApplication.translate("MainWindow", u"Add Material Assignment", None))
        ___qtablewidgetitem239 = self.TMaterialAssignment_2.horizontalHeaderItem(0)
        ___qtablewidgetitem239.setText(QCoreApplication.translate("MainWindow", u"Region", None));
        ___qtablewidgetitem240 = self.TMaterialAssignment_2.horizontalHeaderItem(1)
        ___qtablewidgetitem240.setText(QCoreApplication.translate("MainWindow", u"Material", None));
        self.BDeleteMaterialAssignment_2.setText(QCoreApplication.translate("MainWindow", u"Delete Material Assignment", None))
        self.DefineAssignMats.setTabText(self.DefineAssignMats.indexOf(self.AssignMaterials), QCoreApplication.translate("MainWindow", u"Assign Materials", None))
        self.LBoundaryConditions.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/brick.svg\" width=\"40\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Boundary Conditions</span></p></body></html>", None))
        self.INSelectBoundaryConditions.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectBoundaryConditions.setItemText(1, QCoreApplication.translate("MainWindow", u"Homogenization", None))
        self.INSelectBoundaryConditions.setItemText(2, QCoreApplication.translate("MainWindow", u"Bulk Forming", None))

        self.Lbndry.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Boundary:</p></body></html>", None))
        self.INBoundary.setItemText(0, QCoreApplication.translate("MainWindow", u"x_plane", None))
        self.INBoundary.setItemText(1, QCoreApplication.translate("MainWindow", u"y_plane", None))
        self.INBoundary.setItemText(2, QCoreApplication.translate("MainWindow", u"z_plane", None))

        self.Lvalue.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Plane Position:</p></body></html>", None))
        self.INPlanePosition.setText("")
        self.Ltype.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Type:</p></body></html>", None))
        self.INType.setItemText(0, QCoreApplication.translate("MainWindow", u"reflected", None))
        self.INType.setItemText(1, QCoreApplication.translate("MainWindow", u"fixed", None))
        self.INType.setItemText(2, QCoreApplication.translate("MainWindow", u"velocity", None))

        self.LVel0.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>vel_0: </p></body></html>", None))
        self.LVel1.setText(QCoreApplication.translate("MainWindow", u"vel_1: ", None))
        self.LVelstart.setText(QCoreApplication.translate("MainWindow", u"vel_t_start: ", None))
        self.LVelend.setText(QCoreApplication.translate("MainWindow", u"vel_t_end: ", None))
        self.BaddBC.setText(QCoreApplication.translate("MainWindow", u"Add Boundary Condition", None))
        ___qtablewidgetitem241 = self.TBoundaryConditions.horizontalHeaderItem(0)
        ___qtablewidgetitem241.setText(QCoreApplication.translate("MainWindow", u"Boundary", None));
        ___qtablewidgetitem242 = self.TBoundaryConditions.horizontalHeaderItem(1)
        ___qtablewidgetitem242.setText(QCoreApplication.translate("MainWindow", u"Plane Position", None));
        ___qtablewidgetitem243 = self.TBoundaryConditions.horizontalHeaderItem(2)
        ___qtablewidgetitem243.setText(QCoreApplication.translate("MainWindow", u"Type", None));
        ___qtablewidgetitem244 = self.TBoundaryConditions.horizontalHeaderItem(3)
        ___qtablewidgetitem244.setText(QCoreApplication.translate("MainWindow", u"vel_0", None));
        ___qtablewidgetitem245 = self.TBoundaryConditions.horizontalHeaderItem(4)
        ___qtablewidgetitem245.setText(QCoreApplication.translate("MainWindow", u"vel_1", None));
        ___qtablewidgetitem246 = self.TBoundaryConditions.horizontalHeaderItem(5)
        ___qtablewidgetitem246.setText(QCoreApplication.translate("MainWindow", u"vel_t_start", None));
        ___qtablewidgetitem247 = self.TBoundaryConditions.horizontalHeaderItem(6)
        ___qtablewidgetitem247.setText(QCoreApplication.translate("MainWindow", u"vel_t_end", None));
        self.BdeleteBC.setText(QCoreApplication.translate("MainWindow", u"Delete Boundary Condition", None))
        self.LbulkBC.setText(QCoreApplication.translate("MainWindow", u"Type:", None))
        self.INbulkBC.setItemText(0, QCoreApplication.translate("MainWindow", u"Custom", None))
        self.INbulkBC.setItemText(1, QCoreApplication.translate("MainWindow", u"Tension Z", None))
        self.INbulkBC.setItemText(2, QCoreApplication.translate("MainWindow", u"Compression Z", None))
        self.INbulkBC.setItemText(3, QCoreApplication.translate("MainWindow", u"Shear XY", None))

        self.LVgrad.setText(QCoreApplication.translate("MainWindow", u"VELOCITY GRADIENT", None))
        self.label_12.setText(QCoreApplication.translate("MainWindow", u"--enforced condition--", None))
        self.LVgradi.setText(QCoreApplication.translate("MainWindow", u"--solver starting estimate--", None))
        self.LCstress.setText(QCoreApplication.translate("MainWindow", u"CAUCHY STRESS", None))
        self.label_29.setText(QCoreApplication.translate("MainWindow", u"--enforced condition--", None))
        self.LSolverSettings.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/gear.svg\" width=\"40\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Solver Settings</span></p></body></html>", None))
        self.INSelectSolverSettings.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectSolverSettings.setItemText(1, QCoreApplication.translate("MainWindow", u"Homogenization", None))
        self.INSelectSolverSettings.setItemText(2, QCoreApplication.translate("MainWindow", u"Bulk Forming", None))

        self.Ltime.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Time T:</p></body></html>", None))
        self.INTime.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.Lmindt.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Minimum dt:</p></body></html>", None))
        self.INMindt.setText(QCoreApplication.translate("MainWindow", u"1E-8", None))
        self.Lmaxdt.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Maximum dt:</p></body></html>", None))
        self.INMaxdt.setText(QCoreApplication.translate("MainWindow", u"1E-2", None))
        self.Linitdt.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Initial dt:</p></body></html>", None))
        self.INInitialdt.setText(QCoreApplication.translate("MainWindow", u"1E-5", None))
        self.Lmaxcycle.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Maximum # of cycles:</p></body></html>", None))
        self.INmaxcycles.setText(QCoreApplication.translate("MainWindow", u"2000000", None))
        self.LGraphicsOutput.setText(QCoreApplication.translate("MainWindow", u"Graphics output step:", None))
        self.INGraphicsOutput.setText(QCoreApplication.translate("MainWindow", u"0.25", None))
        self.INHomogenizationSettingD.setText(QCoreApplication.translate("MainWindow", u"Default", None))
        self.INHomogenizationSettingC.setText(QCoreApplication.translate("MainWindow", u"Custom ", None))
        self.LNumberOfSteps.setText(QCoreApplication.translate("MainWindow", u"Number of steps:", None))
        self.INNumberOfSteps.setInputMask("")
        self.INNumberOfSteps.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.INNumberOfSteps.setPlaceholderText("")
        self.LErrorTolerance.setText(QCoreApplication.translate("MainWindow", u"Error Tolerance:", None))
        self.INErrorTolerance.setText(QCoreApplication.translate("MainWindow", u"10", None))
        self.LMaxIterations.setText(QCoreApplication.translate("MainWindow", u"Maximum Iterations:", None))
        self.INMaxIterations.setText(QCoreApplication.translate("MainWindow", u"500", None))
        self.label_69.setText(QCoreApplication.translate("MainWindow", u"Solver Type", None))
        self.INLargeStrain.setText(QCoreApplication.translate("MainWindow", u"Large Strain (Default)", None))
        self.INSmallStrain.setText(QCoreApplication.translate("MainWindow", u"Small Strain", None))
        self.label_68.setText(QCoreApplication.translate("MainWindow", u"Solver Parameters", None))
        self.LNumberOfSteps_3.setText(QCoreApplication.translate("MainWindow", u"Load steps:", None))
        self.INBFloadsteps.setInputMask("")
        self.INBFloadsteps.setText(QCoreApplication.translate("MainWindow", u"30", None))
        self.INBFloadsteps.setPlaceholderText("")
        self.LErrorTolerance_3.setText(QCoreApplication.translate("MainWindow", u"Error Tolerance:", None))
        self.INBFerrortol.setText(QCoreApplication.translate("MainWindow", u"0.0000001", None))
        self.LMaxIterations_3.setText(QCoreApplication.translate("MainWindow", u"Maximum Iterations:", None))
        self.INBFmaxiter.setText(QCoreApplication.translate("MainWindow", u"500", None))
        self.LBFdt.setText(QCoreApplication.translate("MainWindow", u"dt:", None))
        self.INBFdt.setText(QCoreApplication.translate("MainWindow", u"0.00005", None))
        self.LBFoutputsteps.setText(QCoreApplication.translate("MainWindow", u"Output every X steps:", None))
        self.INBFoutputsteps.setText(QCoreApplication.translate("MainWindow", u"30", None))
        self.LRun.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/Play.svg\" width=\"45\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Run</span></p></body></html>", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"Select Solver", None))
        self.INRunSelection.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INRunSelection.setItemText(1, QCoreApplication.translate("MainWindow", u"Homogenization", None))
        self.INRunSelection.setItemText(2, QCoreApplication.translate("MainWindow", u"Bulk Forming", None))

        self.BRunSGH.setText(QCoreApplication.translate("MainWindow", u"Run SGH Solver", None))
        self.INHRunLocally.setText(QCoreApplication.translate("MainWindow", u"Run Locally", None))
        self.INHWriteFiles.setText(QCoreApplication.translate("MainWindow", u"Write Input Files", None))
        self.label_53.setText(QCoreApplication.translate("MainWindow", u"Job Type", None))
        self.INSingleJob.setText(QCoreApplication.translate("MainWindow", u"Single Job", None))
        self.INBatchJob.setText(QCoreApplication.translate("MainWindow", u"Batch Job", None))
        self.LBatchType.setText(QCoreApplication.translate("MainWindow", u"Batch Type:", None))
        self.INBatchType.setItemText(0, QCoreApplication.translate("MainWindow", u"Geometry", None))

        self.INSaveMaterialFiles.setText(QCoreApplication.translate("MainWindow", u"Save Material Constants Files", None))
        self.INSaveAllFiles.setText(QCoreApplication.translate("MainWindow", u"Save All Files", None))
        self.BSelectGeometryFiles.setText(QCoreApplication.translate("MainWindow", u"Select Geometry Files", None))
        self.LHomogenizingFile.setText(QCoreApplication.translate("MainWindow", u"Homogenizing File:", None))
        self.INFileNumber.setText(QCoreApplication.translate("MainWindow", u"1/1", None))
        self.label_54.setText(QCoreApplication.translate("MainWindow", u"Run Type", None))
        self.INHSerial.setText(QCoreApplication.translate("MainWindow", u"Serial", None))
        self.INHParallel.setText(QCoreApplication.translate("MainWindow", u"Parallel", None))
        self.label_58.setText(QCoreApplication.translate("MainWindow", u"MPI Ranks:", None))
        self.label_65.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" color:#797979;\">Warning: mpi could be slower than serial for small models due to initialization and overhead.</span></p></body></html>", None))
        self.label_64.setText(QCoreApplication.translate("MainWindow", u"Convergence Tuning", None))
#if QT_CONFIG(tooltip)
        self.INHAutomatic.setToolTip(QCoreApplication.translate("MainWindow", u"Automatically adjust model parameters if solution convergence is not met.", None))
#endif // QT_CONFIG(tooltip)
        self.INHAutomatic.setText(QCoreApplication.translate("MainWindow", u"Automatic", None))
#if QT_CONFIG(tooltip)
        self.INHManual.setToolTip(QCoreApplication.translate("MainWindow", u"Outputs a convergence warning. Requires manual parameter adjustment and manual re-run.", None))
#endif // QT_CONFIG(tooltip)
        self.INHManual.setText(QCoreApplication.translate("MainWindow", u"Manual", None))
        self.BRunEVPFFT2.setText(QCoreApplication.translate("MainWindow", u"Run Homogenization", None))
        self.BKillEVPFFT2.setText(QCoreApplication.translate("MainWindow", u"Terminate", None))
        self.BHWriteFiles.setText(QCoreApplication.translate("MainWindow", u"Write Input Files", None))
        self.INBFRunLocally.setText(QCoreApplication.translate("MainWindow", u"Run Locally", None))
        self.INBFWriteFiles.setText(QCoreApplication.translate("MainWindow", u"Write Input Files", None))
        self.label_62.setText(QCoreApplication.translate("MainWindow", u"Run Type", None))
        self.INBFSerial.setText(QCoreApplication.translate("MainWindow", u"Serial", None))
        self.INBFParallel.setText(QCoreApplication.translate("MainWindow", u"Parallel", None))
        self.label_63.setText(QCoreApplication.translate("MainWindow", u"MPI Ranks:", None))
        self.BRunBulkForming.setText(QCoreApplication.translate("MainWindow", u"Run Bulk Forming Simulation", None))
        self.BKillBulkForming.setText(QCoreApplication.translate("MainWindow", u"Terminate", None))
        self.BBFWriteFiles.setText(QCoreApplication.translate("MainWindow", u"Write Input Files", None))
        self.LPostprocessing.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/lens.svg\" width=\"40\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Postprocessing</span></p></body></html>", None))
        self.INSelectPostprocessing.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectPostprocessing.setItemText(1, QCoreApplication.translate("MainWindow", u"Homogenization", None))
        self.INSelectPostprocessing.setItemText(2, QCoreApplication.translate("MainWindow", u"Bulk Forming", None))

        self.LResultsSGH.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">RESULTS</span></p></body></html>", None))
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
        self.INBCFile.setItemText(0, QCoreApplication.translate("MainWindow", u"Tension x", None))
        self.INBCFile.setItemText(1, QCoreApplication.translate("MainWindow", u"Tension y", None))
        self.INBCFile.setItemText(2, QCoreApplication.translate("MainWindow", u"Tension z", None))
        self.INBCFile.setItemText(3, QCoreApplication.translate("MainWindow", u"Shear xy", None))
        self.INBCFile.setItemText(4, QCoreApplication.translate("MainWindow", u"Shear xz", None))
        self.INBCFile.setItemText(5, QCoreApplication.translate("MainWindow", u"Shear yz", None))

        self.INPreviewResults.setItemText(0, QCoreApplication.translate("MainWindow", u"vm-stress", None))
        self.INPreviewResults.setItemText(1, QCoreApplication.translate("MainWindow", u"vm-strain", None))

        self.INResultRegion.setItemText(0, QCoreApplication.translate("MainWindow", u"Part", None))
        self.INResultRegion.setItemText(1, QCoreApplication.translate("MainWindow", u"Void", None))
        self.INResultRegion.setItemText(2, QCoreApplication.translate("MainWindow", u"Part + Void", None))

        self.label_66.setText(QCoreApplication.translate("MainWindow", u"Deformation Scale Factor:", None))
        self.BOpenParaview.setText(QCoreApplication.translate("MainWindow", u"Open Paraview", None))
        ___qtablewidgetitem248 = self.THomogenization.horizontalHeaderItem(0)
        ___qtablewidgetitem248.setText(QCoreApplication.translate("MainWindow", u"Homogenized Elastic Constants", None));
        ___qtablewidgetitem249 = self.THomogenization.verticalHeaderItem(0)
        ___qtablewidgetitem249.setText(QCoreApplication.translate("MainWindow", u"Exx", None));
        ___qtablewidgetitem250 = self.THomogenization.verticalHeaderItem(1)
        ___qtablewidgetitem250.setText(QCoreApplication.translate("MainWindow", u"Eyy", None));
        ___qtablewidgetitem251 = self.THomogenization.verticalHeaderItem(2)
        ___qtablewidgetitem251.setText(QCoreApplication.translate("MainWindow", u"Ezz", None));
        ___qtablewidgetitem252 = self.THomogenization.verticalHeaderItem(3)
        ___qtablewidgetitem252.setText(QCoreApplication.translate("MainWindow", u"NUxy", None));
        ___qtablewidgetitem253 = self.THomogenization.verticalHeaderItem(4)
        ___qtablewidgetitem253.setText(QCoreApplication.translate("MainWindow", u"NUzx", None));
        ___qtablewidgetitem254 = self.THomogenization.verticalHeaderItem(5)
        ___qtablewidgetitem254.setText(QCoreApplication.translate("MainWindow", u"NUyz", None));
        ___qtablewidgetitem255 = self.THomogenization.verticalHeaderItem(6)
        ___qtablewidgetitem255.setText(QCoreApplication.translate("MainWindow", u"Gxy", None));
        ___qtablewidgetitem256 = self.THomogenization.verticalHeaderItem(7)
        ___qtablewidgetitem256.setText(QCoreApplication.translate("MainWindow", u"Gxz", None));
        ___qtablewidgetitem257 = self.THomogenization.verticalHeaderItem(8)
        ___qtablewidgetitem257.setText(QCoreApplication.translate("MainWindow", u"Gyz", None));

        __sortingEnabled7 = self.THomogenization.isSortingEnabled()
        self.THomogenization.setSortingEnabled(False)
        ___qtablewidgetitem258 = self.THomogenization.item(0, 1)
        ___qtablewidgetitem258.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None));
        ___qtablewidgetitem259 = self.THomogenization.item(1, 1)
        ___qtablewidgetitem259.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None));
        ___qtablewidgetitem260 = self.THomogenization.item(2, 1)
        ___qtablewidgetitem260.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None));
        ___qtablewidgetitem261 = self.THomogenization.item(6, 1)
        ___qtablewidgetitem261.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None));
        ___qtablewidgetitem262 = self.THomogenization.item(7, 1)
        ___qtablewidgetitem262.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None));
        ___qtablewidgetitem263 = self.THomogenization.item(8, 1)
        ___qtablewidgetitem263.setText(QCoreApplication.translate("MainWindow", u"[Pa]", None));
        self.THomogenization.setSortingEnabled(__sortingEnabled7)

        self.LHomogenizationJobDir.setText(QCoreApplication.translate("MainWindow", u"Job Directory:", None))
        self.INBFResults.setItemText(0, QCoreApplication.translate("MainWindow", u"vm-stress", None))
        self.INBFResults.setItemText(1, QCoreApplication.translate("MainWindow", u"vm-strain", None))

        self.label_67.setText(QCoreApplication.translate("MainWindow", u"Deformation Scale Factor:", None))
        self.BBFParaview.setText(QCoreApplication.translate("MainWindow", u"Open Paraview", None))
        self.LBFJobDir.setText(QCoreApplication.translate("MainWindow", u"Job Directory:", None))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.Results), QCoreApplication.translate("MainWindow", u"Results", None))
#if QT_CONFIG(tooltip)
        self.BResetCamera.setToolTip(QCoreApplication.translate("MainWindow", u"Reset camera.", None))
#endif // QT_CONFIG(tooltip)
        self.BResetCamera.setText("")
        self.INUnits.setItemText(0, QCoreApplication.translate("MainWindow", u"SI (m, Pa, kg)", None))
        self.INUnits.setItemText(1, QCoreApplication.translate("MainWindow", u"SI (mm, MPa, kg)", None))

        self.label_39.setText(QCoreApplication.translate("MainWindow", u"Unit System:", None))
        self.menuHelp.setTitle(QCoreApplication.translate("MainWindow", u"Help", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
    # retranslateUi

