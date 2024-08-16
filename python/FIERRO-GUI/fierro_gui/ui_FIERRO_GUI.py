# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'FIERRO_GUIOuunyT.ui'
##
## Created by: Qt User Interface Compiler version 6.7.2
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
from PySide6.QtWidgets import (QAbstractItemView, QAbstractScrollArea, QApplication, QComboBox,
    QFormLayout, QFrame, QGraphicsView, QGridLayout,
    QHBoxLayout, QHeaderView, QLabel, QLayout,
    QLineEdit, QMainWindow, QMenu, QMenuBar,
    QPlainTextEdit, QProgressBar, QPushButton, QRadioButton,
    QScrollArea, QSizePolicy, QSpacerItem, QSplitter,
    QStackedWidget, QStatusBar, QTabWidget, QTableWidget,
    QTableWidgetItem, QToolButton, QVBoxLayout, QWidget)
import IconResourceFile_rc
import IconResourceFile_rc
import IconResourceFile_rc
import IconResourceFile_rc

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.setEnabled(True)
        MainWindow.resize(1200, 901)
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setMaximumSize(QSize(1980, 1080))
        icon = QIcon()
        icon.addFile(u":/Logos/Logos/FIERRO.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setAutoFillBackground(True)
        MainWindow.setStyleSheet(u"")
        MainWindow.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        MainWindow.setDockNestingEnabled(True)
        MainWindow.setDockOptions(QMainWindow.DockOption.AllowNestedDocks|QMainWindow.DockOption.AllowTabbedDocks|QMainWindow.DockOption.AnimatedDocks|QMainWindow.DockOption.ForceTabbedDocks|QMainWindow.DockOption.GroupedDragging|QMainWindow.DockOption.VerticalTabs)
        self.actionManual = QAction(MainWindow)
        self.actionManual.setObjectName(u"actionManual")
        self.actionChange_Working_Directory = QAction(MainWindow)
        self.actionChange_Working_Directory.setObjectName(u"actionChange_Working_Directory")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.centralwidget.setEnabled(True)
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.centralwidget.setStyleSheet(u"#ParaviewWindow{\n"
"	background-color: rgb(91, 97, 120);\n"
"}\n"
"#BImportPart:hover, #BDefineMaterial:hover, #BApplyBC:hover, #BSolverSettings:hover, #BRunEVPFFT:hover, #BViewResults:hover, #BGlobalMesh:hover, #BImportPartSGH:hover, #BDefineMaterialSGH:hover, #BAssignMaterialSGH:hover, #BApplyBCSGH:hover, #BSolverSettingsSGH:hover, #BViewResultsSGH:hover, #BRunSGH:hover, #BCreateBasicPart:hover{\n"
"	background-color: rgb(192, 192, 192);\n"
"	border-radius: 15px;\n"
"}\n"
"#BImportPart, #BDefineMaterial, #BApplyBC, #BSolverSettings, #BRunEVPFFT, #BViewResults, #BGlobalMesh{\n"
"	border-style: flat;\n"
"}")
        self.verticalLayout_3 = QVBoxLayout(self.centralwidget)
        self.verticalLayout_3.setSpacing(0)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(-1, 9, -1, 9)
        self.NavigationMenu = QTabWidget(self.centralwidget)
        self.NavigationMenu.setObjectName(u"NavigationMenu")
        self.NavigationMenu.setEnabled(True)
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.NavigationMenu.sizePolicy().hasHeightForWidth())
        self.NavigationMenu.setSizePolicy(sizePolicy1)
        self.NavigationMenu.setMaximumSize(QSize(16777215, 1000000))
        font = QFont()
        font.setBold(True)
        self.NavigationMenu.setFont(font)
        self.NavigationMenu.setMouseTracking(False)
        self.NavigationMenu.setFocusPolicy(Qt.FocusPolicy.TabFocus)
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
        self.NavigationMenu.setMovable(True)
        self.NavigationMenu.setTabBarAutoHide(True)
        self.Title = QWidget()
        self.Title.setObjectName(u"Title")
        icon1 = QIcon()
        icon1.addFile(u":/Logos/Logos/FIERRO_NoText.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.NavigationMenu.addTab(self.Title, icon1, "")
        self.Pipeline = QWidget()
        self.Pipeline.setObjectName(u"Pipeline")
        icon2 = QIcon()
        icon2.addFile(u":/Blue Icons/Blue Icons/Clipboard.svg", QSize(), QIcon.Mode.Selected, QIcon.State.On)
        self.NavigationMenu.addTab(self.Pipeline, icon2, "")
        self.Geometry = QWidget()
        self.Geometry.setObjectName(u"Geometry")
        self.Geometry.setEnabled(True)
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.Geometry.sizePolicy().hasHeightForWidth())
        self.Geometry.setSizePolicy(sizePolicy2)
        self.verticalLayout_6 = QVBoxLayout(self.Geometry)
        self.verticalLayout_6.setSpacing(0)
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.verticalLayout_6.setContentsMargins(0, 0, 0, 0)
        icon3 = QIcon()
        icon3.addFile(u":/Blue Icons/Blue Icons/Cube.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.NavigationMenu.addTab(self.Geometry, icon3, "")
        self.Mesh = QWidget()
        self.Mesh.setObjectName(u"Mesh")
        icon4 = QIcon()
        icon4.addFile(u":/Blue Icons/Blue Icons/mesh.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.NavigationMenu.addTab(self.Mesh, icon4, "")
        self.Solver = QWidget()
        self.Solver.setObjectName(u"Solver")
        self.verticalLayout_7 = QVBoxLayout(self.Solver)
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        icon5 = QIcon()
        icon5.addFile(u":/Blue Icons/Blue Icons/gear.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.NavigationMenu.addTab(self.Solver, icon5, "")
        self.BoundaryConditions = QWidget()
        self.BoundaryConditions.setObjectName(u"BoundaryConditions")
        self.verticalLayout_5 = QVBoxLayout(self.BoundaryConditions)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        icon6 = QIcon()
        icon6.addFile(u":/Blue Icons/Blue Icons/brick.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.NavigationMenu.addTab(self.BoundaryConditions, icon6, "")
        self.Materials = QWidget()
        self.Materials.setObjectName(u"Materials")
        self.verticalLayout_11 = QVBoxLayout(self.Materials)
        self.verticalLayout_11.setObjectName(u"verticalLayout_11")
        icon7 = QIcon()
        icon7.addFile(u":/Blue Icons/Blue Icons/mine.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.NavigationMenu.addTab(self.Materials, icon7, "")
        self.Run = QWidget()
        self.Run.setObjectName(u"Run")
        icon8 = QIcon()
        icon8.addFile(u":/Blue Icons/Blue Icons/Play.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.NavigationMenu.addTab(self.Run, icon8, "")
        self.Postprocessing = QWidget()
        self.Postprocessing.setObjectName(u"Postprocessing")
        icon9 = QIcon()
        icon9.addFile(u":/Blue Icons/Blue Icons/lens.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.NavigationMenu.addTab(self.Postprocessing, icon9, "")

        self.verticalLayout_3.addWidget(self.NavigationMenu)

        self.Windows = QFormLayout()
        self.Windows.setObjectName(u"Windows")
        self.Windows.setHorizontalSpacing(0)
        self.Windows.setVerticalSpacing(0)
        self.ParaviewFrame = QFrame(self.centralwidget)
        self.ParaviewFrame.setObjectName(u"ParaviewFrame")
        sizePolicy.setHeightForWidth(self.ParaviewFrame.sizePolicy().hasHeightForWidth())
        self.ParaviewFrame.setSizePolicy(sizePolicy)
        self.ParaviewFrame.setMinimumSize(QSize(0, 0))
        self.ParaviewFrame.setMaximumSize(QSize(1000000, 1000000))
        self.ParaviewFrame.setFocusPolicy(Qt.FocusPolicy.TabFocus)
        self.ParaviewFrame.setContextMenuPolicy(Qt.ContextMenuPolicy.DefaultContextMenu)
        self.ParaviewFrame.setFrameShape(QFrame.Shape.NoFrame)
        self.ParaviewFrame.setFrameShadow(QFrame.Shadow.Plain)
        self.ParaviewFrame.setLineWidth(1)
        self.verticalLayout_19 = QVBoxLayout(self.ParaviewFrame)
        self.verticalLayout_19.setSpacing(6)
        self.verticalLayout_19.setObjectName(u"verticalLayout_19")
        self.verticalLayout_19.setContentsMargins(0, 0, 0, 0)
        self.splitter_2 = QSplitter(self.ParaviewFrame)
        self.splitter_2.setObjectName(u"splitter_2")
        self.splitter_2.setFrameShape(QFrame.Shape.Box)
        self.splitter_2.setOrientation(Qt.Orientation.Vertical)
        self.splitter_2.setHandleWidth(6)
        self.OutputWindows = QStackedWidget(self.splitter_2)
        self.OutputWindows.setObjectName(u"OutputWindows")
        sizePolicy3 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.OutputWindows.sizePolicy().hasHeightForWidth())
        self.OutputWindows.setSizePolicy(sizePolicy3)
        self.OutputWindows.setCursor(QCursor(Qt.CursorShape.OpenHandCursor))
        self.OutputWindows.setFrameShape(QFrame.Shape.NoFrame)
        self.OutputWindows.setLineWidth(6)
        self.OutputWindows.setMidLineWidth(0)
        self.ParaviewWindow = QWidget()
        self.ParaviewWindow.setObjectName(u"ParaviewWindow")
        self.paraviewLayout = QVBoxLayout(self.ParaviewWindow)
        self.paraviewLayout.setObjectName(u"paraviewLayout")
        self.paraviewLayout.setContentsMargins(0, 0, 0, 0)
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
        sizePolicy2.setHeightForWidth(self.RunOutputs.sizePolicy().hasHeightForWidth())
        self.RunOutputs.setSizePolicy(sizePolicy2)
        self.RunOutputs.setMaximumSize(QSize(16777215, 250))
        self.RunOutputs.setFrameShape(QFrame.Shape.NoFrame)
        self.RunOutputs.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_22 = QVBoxLayout(self.RunOutputs)
        self.verticalLayout_22.setSpacing(0)
        self.verticalLayout_22.setObjectName(u"verticalLayout_22")
        self.verticalLayout_22.setContentsMargins(0, 0, 0, 0)
        self.RunOutputProgress = QProgressBar(self.RunOutputs)
        self.RunOutputProgress.setObjectName(u"RunOutputProgress")
        self.RunOutputProgress.setMinimumSize(QSize(500, 0))
        self.RunOutputProgress.setMaximumSize(QSize(16777215, 16777215))
        self.RunOutputProgress.setValue(0)

        self.verticalLayout_22.addWidget(self.RunOutputProgress, 0, Qt.AlignmentFlag.AlignHCenter)

        self.graphicsView = QGraphicsView(self.RunOutputs)
        self.graphicsView.setObjectName(u"graphicsView")

        self.verticalLayout_22.addWidget(self.graphicsView)

        self.RunOutputWindow = QPlainTextEdit(self.RunOutputs)
        self.RunOutputWindow.setObjectName(u"RunOutputWindow")
        self.RunOutputWindow.setFrameShape(QFrame.Shape.NoFrame)
        self.RunOutputWindow.setFrameShadow(QFrame.Shadow.Plain)
        self.RunOutputWindow.setTabChangesFocus(False)
        self.RunOutputWindow.setReadOnly(True)

        self.verticalLayout_22.addWidget(self.RunOutputWindow)

        self.splitter_2.addWidget(self.RunOutputs)

        self.verticalLayout_19.addWidget(self.splitter_2)


        self.Windows.setWidget(1, QFormLayout.FieldRole, self.ParaviewFrame)

        self.ToolWindow = QStackedWidget(self.centralwidget)
        self.ToolWindow.setObjectName(u"ToolWindow")
        self.ToolWindow.setEnabled(True)
        sizePolicy4 = QSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Expanding)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(1)
        sizePolicy4.setHeightForWidth(self.ToolWindow.sizePolicy().hasHeightForWidth())
        self.ToolWindow.setSizePolicy(sizePolicy4)
        self.ToolWindow.setMinimumSize(QSize(0, 700))
        self.ToolWindow.setMaximumSize(QSize(400, 1000000))
        self.ToolWindow.setSizeIncrement(QSize(0, 0))
        self.ToolWindow.setBaseSize(QSize(300, 300))
        font1 = QFont()
        font1.setFamilies([u"Hiragino Sans"])
        self.ToolWindow.setFont(font1)
        self.ToolWindow.setAutoFillBackground(False)
        self.ToolWindow.setFrameShape(QFrame.Shape.Box)
        self.ToolWindow.setFrameShadow(QFrame.Shadow.Plain)
        self.ToolWindow.setLineWidth(0)
        self.ToolWindow.setMidLineWidth(1)
        self.TitleTool = QWidget()
        self.TitleTool.setObjectName(u"TitleTool")
        self.verticalLayout_2 = QVBoxLayout(self.TitleTool)
        self.verticalLayout_2.setSpacing(10)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.verticalLayout_2.setContentsMargins(12, 40, 20, 40)
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

        self.verticalLayout_2.addWidget(self.LosAlamosLogo)

        self.verticalSpacer_7 = QSpacerItem(20, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.MinimumExpanding)

        self.verticalLayout_2.addItem(self.verticalSpacer_7)

        self.EVPFFTLogo = QLabel(self.TitleTool)
        self.EVPFFTLogo.setObjectName(u"EVPFFTLogo")
        sizePolicy2.setHeightForWidth(self.EVPFFTLogo.sizePolicy().hasHeightForWidth())
        self.EVPFFTLogo.setSizePolicy(sizePolicy2)
        self.EVPFFTLogo.setMinimumSize(QSize(0, 0))
        self.EVPFFTLogo.setMaximumSize(QSize(225, 175))
        self.EVPFFTLogo.setPixmap(QPixmap(u":/Logos/Logos/FIERRO.png"))
        self.EVPFFTLogo.setScaledContents(True)
        self.EVPFFTLogo.setWordWrap(False)
        self.EVPFFTLogo.setIndent(-1)

        self.verticalLayout_2.addWidget(self.EVPFFTLogo, 0, Qt.AlignmentFlag.AlignHCenter|Qt.AlignmentFlag.AlignVCenter)

        self.verticalSpacer_6 = QSpacerItem(20, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.MinimumExpanding)

        self.verticalLayout_2.addItem(self.verticalSpacer_6)

        self.AdditionalSoftware = QFrame(self.TitleTool)
        self.AdditionalSoftware.setObjectName(u"AdditionalSoftware")
        self.AdditionalSoftware.setFrameShape(QFrame.Shape.NoFrame)
        self.AdditionalSoftware.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_10 = QVBoxLayout(self.AdditionalSoftware)
        self.verticalLayout_10.setSpacing(8)
        self.verticalLayout_10.setObjectName(u"verticalLayout_10")
        self.verticalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.LAdditionalSoftware = QLabel(self.AdditionalSoftware)
        self.LAdditionalSoftware.setObjectName(u"LAdditionalSoftware")
        font2 = QFont()
        font2.setPointSize(16)
        self.LAdditionalSoftware.setFont(font2)

        self.verticalLayout_10.addWidget(self.LAdditionalSoftware, 0, Qt.AlignmentFlag.AlignBottom)

        self.AdditionalSoftwareLogos = QFrame(self.AdditionalSoftware)
        self.AdditionalSoftwareLogos.setObjectName(u"AdditionalSoftwareLogos")
        sizePolicy2.setHeightForWidth(self.AdditionalSoftwareLogos.sizePolicy().hasHeightForWidth())
        self.AdditionalSoftwareLogos.setSizePolicy(sizePolicy2)
        self.AdditionalSoftwareLogos.setFrameShape(QFrame.Shape.NoFrame)
        self.AdditionalSoftwareLogos.setFrameShadow(QFrame.Shadow.Raised)
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


        self.verticalLayout_2.addWidget(self.AdditionalSoftware, 0, Qt.AlignmentFlag.AlignTop)

        self.ToolWindow.addWidget(self.TitleTool)
        self.PipelineTool = QWidget()
        self.PipelineTool.setObjectName(u"PipelineTool")
        self.Import_2 = QFrame(self.PipelineTool)
        self.Import_2.setObjectName(u"Import_2")
        self.Import_2.setEnabled(True)
        self.Import_2.setGeometry(QRect(9, 0, 382, 189))
        sizePolicy1.setHeightForWidth(self.Import_2.sizePolicy().hasHeightForWidth())
        self.Import_2.setSizePolicy(sizePolicy1)
        self.Import_2.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.Import_2.setAutoFillBackground(False)
        self.Import_2.setFrameShape(QFrame.Shape.NoFrame)
        self.Import_2.setFrameShadow(QFrame.Shadow.Raised)
        self._3 = QFormLayout(self.Import_2)
        self._3.setObjectName(u"_3")
        self._3.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        self._3.setRowWrapPolicy(QFormLayout.RowWrapPolicy.DontWrapRows)
        self._3.setLabelAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignTop)
        self._3.setContentsMargins(-1, 6, -1, -1)
        self.LPipeline = QLabel(self.Import_2)
        self.LPipeline.setObjectName(u"LPipeline")
        sizePolicy1.setHeightForWidth(self.LPipeline.sizePolicy().hasHeightForWidth())
        self.LPipeline.setSizePolicy(sizePolicy1)
        font3 = QFont()
        font3.setFamilies([u"Sans Serif"])
        font3.setPointSize(16)
        self.LPipeline.setFont(font3)
        self.LPipeline.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft)
        self.LPipeline.setMargin(0)

        self._3.setWidget(0, QFormLayout.SpanningRole, self.LPipeline)

        self.INPipelineSelection = QComboBox(self.Import_2)
        self.INPipelineSelection.addItem("")
        self.INPipelineSelection.addItem("")
        self.INPipelineSelection.addItem("")
        self.INPipelineSelection.setObjectName(u"INPipelineSelection")
        sizePolicy1.setHeightForWidth(self.INPipelineSelection.sizePolicy().hasHeightForWidth())
        self.INPipelineSelection.setSizePolicy(sizePolicy1)
        font4 = QFont()
        font4.setFamilies([u"Hiragino Sans"])
        font4.setPointSize(10)
        self.INPipelineSelection.setFont(font4)
        self.INPipelineSelection.setIconSize(QSize(18, 18))

        self._3.setWidget(1, QFormLayout.FieldRole, self.INPipelineSelection)

        self.line_4 = QFrame(self.Import_2)
        self.line_4.setObjectName(u"line_4")
        self.line_4.setFrameShape(QFrame.Shape.HLine)
        self.line_4.setFrameShadow(QFrame.Shadow.Sunken)

        self._3.setWidget(2, QFormLayout.SpanningRole, self.line_4)

        self.LPipelineSelection = QLabel(self.Import_2)
        self.LPipelineSelection.setObjectName(u"LPipelineSelection")
        self.LPipelineSelection.setFont(font1)
        self.LPipelineSelection.setWordWrap(False)

        self._3.setWidget(1, QFormLayout.LabelRole, self.LPipelineSelection)

        self.ToolWindow.addWidget(self.PipelineTool)
        self.ImportGeometryTool = QWidget()
        self.ImportGeometryTool.setObjectName(u"ImportGeometryTool")
        sizePolicy3.setHeightForWidth(self.ImportGeometryTool.sizePolicy().hasHeightForWidth())
        self.ImportGeometryTool.setSizePolicy(sizePolicy3)
        self.verticalLayout_15 = QVBoxLayout(self.ImportGeometryTool)
        self.verticalLayout_15.setObjectName(u"verticalLayout_15")
        self.verticalLayout_15.setSizeConstraint(QLayout.SizeConstraint.SetDefaultConstraint)
        self.verticalLayout_15.setContentsMargins(-1, 0, -1, 0)
        self.Import = QFrame(self.ImportGeometryTool)
        self.Import.setObjectName(u"Import")
        self.Import.setEnabled(True)
        sizePolicy1.setHeightForWidth(self.Import.sizePolicy().hasHeightForWidth())
        self.Import.setSizePolicy(sizePolicy1)
        self.Import.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.Import.setAutoFillBackground(False)
        self.Import.setFrameShape(QFrame.Shape.NoFrame)
        self.Import.setFrameShadow(QFrame.Shadow.Raised)
        self._2 = QFormLayout(self.Import)
        self._2.setObjectName(u"_2")
        self._2.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        self._2.setRowWrapPolicy(QFormLayout.RowWrapPolicy.DontWrapRows)
        self._2.setLabelAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignTop)
        self._2.setContentsMargins(-1, 6, -1, -1)
        self.LGeometryInformation = QLabel(self.Import)
        self.LGeometryInformation.setObjectName(u"LGeometryInformation")
        sizePolicy1.setHeightForWidth(self.LGeometryInformation.sizePolicy().hasHeightForWidth())
        self.LGeometryInformation.setSizePolicy(sizePolicy1)
        self.LGeometryInformation.setFont(font3)
        self.LGeometryInformation.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft)
        self.LGeometryInformation.setMargin(0)

        self._2.setWidget(0, QFormLayout.SpanningRole, self.LGeometryInformation)

        self.LPartName_3 = QLabel(self.Import)
        self.LPartName_3.setObjectName(u"LPartName_3")
        self.LPartName_3.setFont(font1)
        self.LPartName_3.setWordWrap(False)

        self._2.setWidget(1, QFormLayout.LabelRole, self.LPartName_3)

        self.INSelectGeometryImport = QComboBox(self.Import)
        icon10 = QIcon()
        icon10.addFile(u":/Blue Icons/Blue Icons/DownloadCube.svg", QSize(), QIcon.Mode.Selected, QIcon.State.On)
        self.INSelectGeometryImport.addItem(icon10, "")
        self.INSelectGeometryImport.addItem("")
        self.INSelectGeometryImport.addItem("")
        self.INSelectGeometryImport.addItem("")
        self.INSelectGeometryImport.setObjectName(u"INSelectGeometryImport")
        sizePolicy1.setHeightForWidth(self.INSelectGeometryImport.sizePolicy().hasHeightForWidth())
        self.INSelectGeometryImport.setSizePolicy(sizePolicy1)
        self.INSelectGeometryImport.setFont(font4)
        self.INSelectGeometryImport.setIconSize(QSize(18, 18))

        self._2.setWidget(1, QFormLayout.FieldRole, self.INSelectGeometryImport)

        self.LPartName = QLabel(self.Import)
        self.LPartName.setObjectName(u"LPartName")
        self.LPartName.setFont(font1)

        self._2.setWidget(2, QFormLayout.LabelRole, self.LPartName)

        self.INPartName = QLineEdit(self.Import)
        self.INPartName.setObjectName(u"INPartName")
        self.INPartName.setFont(font1)

        self._2.setWidget(2, QFormLayout.FieldRole, self.INPartName)

        self.BUploadGeometryFile = QPushButton(self.Import)
        self.BUploadGeometryFile.setObjectName(u"BUploadGeometryFile")
        self.BUploadGeometryFile.setFont(font1)

        self._2.setWidget(3, QFormLayout.SpanningRole, self.BUploadGeometryFile)

        self.line_3 = QFrame(self.Import)
        self.line_3.setObjectName(u"line_3")
        self.line_3.setFrameShape(QFrame.Shape.HLine)
        self.line_3.setFrameShadow(QFrame.Shadow.Sunken)

        self._2.setWidget(4, QFormLayout.SpanningRole, self.line_3)


        self.verticalLayout_15.addWidget(self.Import)

        self.SAGeometryScrollArea = QScrollArea(self.ImportGeometryTool)
        self.SAGeometryScrollArea.setObjectName(u"SAGeometryScrollArea")
        self.SAGeometryScrollArea.setEnabled(True)
        sizePolicy3.setHeightForWidth(self.SAGeometryScrollArea.sizePolicy().hasHeightForWidth())
        self.SAGeometryScrollArea.setSizePolicy(sizePolicy3)
        self.SAGeometryScrollArea.setMinimumSize(QSize(0, 0))
        self.SAGeometryScrollArea.setMaximumSize(QSize(16777215, 1000))
        self.SAGeometryScrollArea.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.SAGeometryScrollArea.setFrameShape(QFrame.Shape.NoFrame)
        self.SAGeometryScrollArea.setLineWidth(1)
        self.SAGeometryScrollArea.setMidLineWidth(1)
        self.SAGeometryScrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.SAGeometryScrollArea.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustToContentsOnFirstShow)
        self.SAGeometryScrollArea.setWidgetResizable(False)
        self.SAGeometryScrollArea.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignTop)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 368, 700))
        sizePolicy6 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(200)
        sizePolicy6.setHeightForWidth(self.scrollAreaWidgetContents.sizePolicy().hasHeightForWidth())
        self.scrollAreaWidgetContents.setSizePolicy(sizePolicy6)
        self.scrollAreaWidgetContents.setMinimumSize(QSize(0, 700))
        self.scrollAreaWidgetContents.setAcceptDrops(True)
        self.verticalLayout_40 = QVBoxLayout(self.scrollAreaWidgetContents)
        self.verticalLayout_40.setObjectName(u"verticalLayout_40")
        self.verticalLayout_40.setSizeConstraint(QLayout.SizeConstraint.SetDefaultConstraint)
        self.GeometryOptions = QStackedWidget(self.scrollAreaWidgetContents)
        self.GeometryOptions.setObjectName(u"GeometryOptions")
        self.GeometryOptions.setMinimumSize(QSize(0, 1400))
        self.GeometryOptions.setMouseTracking(False)
        self.GeometryOptions.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.GeometryOptions.setAutoFillBackground(False)
        self.ImportPartTool = QWidget()
        self.ImportPartTool.setObjectName(u"ImportPartTool")
        self.BVoxelizeGeometry = QPushButton(self.ImportPartTool)
        self.BVoxelizeGeometry.setObjectName(u"BVoxelizeGeometry")
        self.BVoxelizeGeometry.setEnabled(False)
        self.BVoxelizeGeometry.setGeometry(QRect(4, 420, 340, 21))
        self.BVoxelizeGeometry.setFont(font1)
        self.LDimensions = QLabel(self.ImportPartTool)
        self.LDimensions.setObjectName(u"LDimensions")
        self.LDimensions.setGeometry(QRect(62, 263, 77, 20))
        font5 = QFont()
        font5.setWeight(QFont.Medium)
        self.LDimensions.setFont(font5)
        self.LDimensions.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignTop)
        self.LDimensions.setMargin(2)
        self.LSTLVoxelization = QLabel(self.ImportPartTool)
        self.LSTLVoxelization.setObjectName(u"LSTLVoxelization")
        self.LSTLVoxelization.setGeometry(QRect(76, 7, 200, 16))
        font6 = QFont()
        font6.setPointSize(15)
        font6.setWeight(QFont.DemiBold)
        font6.setUnderline(False)
        font6.setStrikeOut(False)
        self.LSTLVoxelization.setFont(font6)
        self.LSTLVoxelization.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.BCustomDimensions = QRadioButton(self.ImportPartTool)
        self.BCustomDimensions.setObjectName(u"BCustomDimensions")
        self.BCustomDimensions.setEnabled(False)
        self.BCustomDimensions.setGeometry(QRect(64, 284, 200, 18))
        sizePolicy7 = QSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        sizePolicy7.setHorizontalStretch(0)
        sizePolicy7.setVerticalStretch(0)
        sizePolicy7.setHeightForWidth(self.BCustomDimensions.sizePolicy().hasHeightForWidth())
        self.BCustomDimensions.setSizePolicy(sizePolicy7)
        font7 = QFont()
        font7.setFamilies([u"Hiragino Sans"])
        font7.setPointSize(11)
        self.BCustomDimensions.setFont(font7)
        self.BStlDimensions = QRadioButton(self.ImportPartTool)
        self.BStlDimensions.setObjectName(u"BStlDimensions")
        self.BStlDimensions.setEnabled(False)
        self.BStlDimensions.setGeometry(QRect(64, 304, 200, 18))
        self.BStlDimensions.setFont(font7)
        self.BStlDimensions.setChecked(True)
        self.layoutWidget = QWidget(self.ImportPartTool)
        self.layoutWidget.setObjectName(u"layoutWidget")
        self.layoutWidget.setGeometry(QRect(65, 152, 233, 106))
        self.formLayout = QFormLayout(self.layoutWidget)
        self.formLayout.setObjectName(u"formLayout")
        self.formLayout.setHorizontalSpacing(5)
        self.formLayout.setVerticalSpacing(5)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.LOriginX = QLabel(self.layoutWidget)
        self.LOriginX.setObjectName(u"LOriginX")
        self.LOriginX.setEnabled(False)
        self.LOriginX.setFont(font1)

        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.LOriginX)

        self.INOriginX = QLineEdit(self.layoutWidget)
        self.INOriginX.setObjectName(u"INOriginX")
        self.INOriginX.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.INOriginX.sizePolicy().hasHeightForWidth())
        self.INOriginX.setSizePolicy(sizePolicy3)
        self.INOriginX.setFont(font1)
        self.INOriginX.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.INOriginX.setInputMethodHints(Qt.InputMethodHint.ImhHiddenText|Qt.InputMethodHint.ImhLatinOnly|Qt.InputMethodHint.ImhPreferNumbers|Qt.InputMethodHint.ImhUrlCharactersOnly)

        self.formLayout.setWidget(2, QFormLayout.FieldRole, self.INOriginX)

        self.LOriginY = QLabel(self.layoutWidget)
        self.LOriginY.setObjectName(u"LOriginY")
        self.LOriginY.setEnabled(False)
        self.LOriginY.setFont(font1)

        self.formLayout.setWidget(3, QFormLayout.LabelRole, self.LOriginY)

        self.INOriginY = QLineEdit(self.layoutWidget)
        self.INOriginY.setObjectName(u"INOriginY")
        self.INOriginY.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.INOriginY.sizePolicy().hasHeightForWidth())
        self.INOriginY.setSizePolicy(sizePolicy3)
        self.INOriginY.setFont(font1)
        self.INOriginY.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.INOriginY.setInputMethodHints(Qt.InputMethodHint.ImhHiddenText|Qt.InputMethodHint.ImhLatinOnly|Qt.InputMethodHint.ImhPreferNumbers|Qt.InputMethodHint.ImhUrlCharactersOnly)

        self.formLayout.setWidget(3, QFormLayout.FieldRole, self.INOriginY)

        self.LOriginZ = QLabel(self.layoutWidget)
        self.LOriginZ.setObjectName(u"LOriginZ")
        self.LOriginZ.setEnabled(False)
        self.LOriginZ.setFont(font1)

        self.formLayout.setWidget(4, QFormLayout.LabelRole, self.LOriginZ)

        self.INOriginZ = QLineEdit(self.layoutWidget)
        self.INOriginZ.setObjectName(u"INOriginZ")
        self.INOriginZ.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.INOriginZ.sizePolicy().hasHeightForWidth())
        self.INOriginZ.setSizePolicy(sizePolicy3)
        self.INOriginZ.setFont(font1)
        self.INOriginZ.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.INOriginZ.setInputMethodHints(Qt.InputMethodHint.ImhHiddenText|Qt.InputMethodHint.ImhLatinOnly|Qt.InputMethodHint.ImhPreferNumbers|Qt.InputMethodHint.ImhUrlCharactersOnly)

        self.formLayout.setWidget(4, QFormLayout.FieldRole, self.INOriginZ)

        self.LOriginPoint = QLabel(self.layoutWidget)
        self.LOriginPoint.setObjectName(u"LOriginPoint")
        self.LOriginPoint.setFont(font5)
        self.LOriginPoint.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignTop)
        self.LOriginPoint.setMargin(0)

        self.formLayout.setWidget(1, QFormLayout.SpanningRole, self.LOriginPoint)

        self.layoutWidget1 = QWidget(self.ImportPartTool)
        self.layoutWidget1.setObjectName(u"layoutWidget1")
        self.layoutWidget1.setGeometry(QRect(65, 38, 233, 106))
        self.formLayout_23 = QFormLayout(self.layoutWidget1)
        self.formLayout_23.setObjectName(u"formLayout_23")
        self.formLayout_23.setHorizontalSpacing(5)
        self.formLayout_23.setVerticalSpacing(5)
        self.formLayout_23.setContentsMargins(0, 0, 0, 0)
        self.LNumberOfVoxelsX = QLabel(self.layoutWidget1)
        self.LNumberOfVoxelsX.setObjectName(u"LNumberOfVoxelsX")
        self.LNumberOfVoxelsX.setEnabled(False)
        font8 = QFont()
        font8.setFamilies([u"Hiragino Sans"])
        font8.setPointSize(9)
        font8.setBold(False)
        font8.setItalic(False)
        font8.setKerning(True)
        self.LNumberOfVoxelsX.setFont(font8)
        self.LNumberOfVoxelsX.setTextFormat(Qt.TextFormat.AutoText)

        self.formLayout_23.setWidget(1, QFormLayout.LabelRole, self.LNumberOfVoxelsX)

        self.INNumberOfVoxelsX = QLineEdit(self.layoutWidget1)
        self.INNumberOfVoxelsX.setObjectName(u"INNumberOfVoxelsX")
        self.INNumberOfVoxelsX.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.INNumberOfVoxelsX.sizePolicy().hasHeightForWidth())
        self.INNumberOfVoxelsX.setSizePolicy(sizePolicy3)
        font9 = QFont()
        font9.setFamilies([u"Hiragino Sans"])
        font9.setPointSize(9)
        self.INNumberOfVoxelsX.setFont(font9)
        self.INNumberOfVoxelsX.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.INNumberOfVoxelsX.setInputMethodHints(Qt.InputMethodHint.ImhHiddenText|Qt.InputMethodHint.ImhLatinOnly|Qt.InputMethodHint.ImhPreferNumbers|Qt.InputMethodHint.ImhUrlCharactersOnly)
        self.INNumberOfVoxelsX.setDragEnabled(False)

        self.formLayout_23.setWidget(1, QFormLayout.FieldRole, self.INNumberOfVoxelsX)

        self.LNumberOfVoxelsY = QLabel(self.layoutWidget1)
        self.LNumberOfVoxelsY.setObjectName(u"LNumberOfVoxelsY")
        self.LNumberOfVoxelsY.setEnabled(False)
        self.LNumberOfVoxelsY.setFont(font9)

        self.formLayout_23.setWidget(2, QFormLayout.LabelRole, self.LNumberOfVoxelsY)

        self.INNumberOfVoxelsY = QLineEdit(self.layoutWidget1)
        self.INNumberOfVoxelsY.setObjectName(u"INNumberOfVoxelsY")
        self.INNumberOfVoxelsY.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.INNumberOfVoxelsY.sizePolicy().hasHeightForWidth())
        self.INNumberOfVoxelsY.setSizePolicy(sizePolicy3)
        self.INNumberOfVoxelsY.setFont(font9)
        self.INNumberOfVoxelsY.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.INNumberOfVoxelsY.setInputMethodHints(Qt.InputMethodHint.ImhHiddenText|Qt.InputMethodHint.ImhLatinOnly|Qt.InputMethodHint.ImhPreferNumbers|Qt.InputMethodHint.ImhUrlCharactersOnly)

        self.formLayout_23.setWidget(2, QFormLayout.FieldRole, self.INNumberOfVoxelsY)

        self.LNumberOfVoxelsZ = QLabel(self.layoutWidget1)
        self.LNumberOfVoxelsZ.setObjectName(u"LNumberOfVoxelsZ")
        self.LNumberOfVoxelsZ.setEnabled(False)
        self.LNumberOfVoxelsZ.setFont(font9)

        self.formLayout_23.setWidget(3, QFormLayout.LabelRole, self.LNumberOfVoxelsZ)

        self.INNumberOfVoxelsZ = QLineEdit(self.layoutWidget1)
        self.INNumberOfVoxelsZ.setObjectName(u"INNumberOfVoxelsZ")
        self.INNumberOfVoxelsZ.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.INNumberOfVoxelsZ.sizePolicy().hasHeightForWidth())
        self.INNumberOfVoxelsZ.setSizePolicy(sizePolicy3)
        self.INNumberOfVoxelsZ.setFont(font9)
        self.INNumberOfVoxelsZ.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.INNumberOfVoxelsZ.setInputMethodHints(Qt.InputMethodHint.ImhHiddenText|Qt.InputMethodHint.ImhLatinOnly|Qt.InputMethodHint.ImhPreferNumbers|Qt.InputMethodHint.ImhUrlCharactersOnly)

        self.formLayout_23.setWidget(3, QFormLayout.FieldRole, self.INNumberOfVoxelsZ)

        self.LVoxelCount = QLabel(self.layoutWidget1)
        self.LVoxelCount.setObjectName(u"LVoxelCount")
        self.LVoxelCount.setFont(font5)
        self.LVoxelCount.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignTop)
        self.LVoxelCount.setMargin(0)

        self.formLayout_23.setWidget(0, QFormLayout.SpanningRole, self.LVoxelCount)

        self.layoutWidget2 = QWidget(self.ImportPartTool)
        self.layoutWidget2.setObjectName(u"layoutWidget2")
        self.layoutWidget2.setGeometry(QRect(65, 328, 233, 80))
        self.formLayout_39 = QFormLayout(self.layoutWidget2)
        self.formLayout_39.setObjectName(u"formLayout_39")
        self.formLayout_39.setContentsMargins(0, 0, 0, 0)
        self.LLengthX = QLabel(self.layoutWidget2)
        self.LLengthX.setObjectName(u"LLengthX")
        self.LLengthX.setEnabled(False)
        self.LLengthX.setFont(font1)

        self.formLayout_39.setWidget(0, QFormLayout.LabelRole, self.LLengthX)

        self.LLengthY = QLabel(self.layoutWidget2)
        self.LLengthY.setObjectName(u"LLengthY")
        self.LLengthY.setEnabled(False)
        self.LLengthY.setFont(font1)

        self.formLayout_39.setWidget(1, QFormLayout.LabelRole, self.LLengthY)

        self.INLengthY = QLineEdit(self.layoutWidget2)
        self.INLengthY.setObjectName(u"INLengthY")
        self.INLengthY.setEnabled(False)
        self.INLengthY.setMinimumSize(QSize(208, 0))
        self.INLengthY.setFont(font1)
        self.INLengthY.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.INLengthY.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

        self.formLayout_39.setWidget(1, QFormLayout.FieldRole, self.INLengthY)

        self.LLengthZ = QLabel(self.layoutWidget2)
        self.LLengthZ.setObjectName(u"LLengthZ")
        self.LLengthZ.setEnabled(False)
        self.LLengthZ.setFont(font1)

        self.formLayout_39.setWidget(2, QFormLayout.LabelRole, self.LLengthZ)

        self.INLengthZ = QLineEdit(self.layoutWidget2)
        self.INLengthZ.setObjectName(u"INLengthZ")
        self.INLengthZ.setEnabled(False)
        sizePolicy8 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        sizePolicy8.setHorizontalStretch(0)
        sizePolicy8.setVerticalStretch(0)
        sizePolicy8.setHeightForWidth(self.INLengthZ.sizePolicy().hasHeightForWidth())
        self.INLengthZ.setSizePolicy(sizePolicy8)
        self.INLengthZ.setMinimumSize(QSize(210, 0))
        self.INLengthZ.setFont(font1)
        self.INLengthZ.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.INLengthZ.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

        self.formLayout_39.setWidget(2, QFormLayout.FieldRole, self.INLengthZ)

        self.INLengthX = QLineEdit(self.layoutWidget2)
        self.INLengthX.setObjectName(u"INLengthX")
        self.INLengthX.setEnabled(False)
        sizePolicy8.setHeightForWidth(self.INLengthX.sizePolicy().hasHeightForWidth())
        self.INLengthX.setSizePolicy(sizePolicy8)
        self.INLengthX.setMinimumSize(QSize(208, 0))
        self.INLengthX.setFont(font1)
        self.INLengthX.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.INLengthX.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

        self.formLayout_39.setWidget(0, QFormLayout.FieldRole, self.INLengthX)

        self.TParts = QTableWidget(self.ImportPartTool)
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
        self.TParts.setGeometry(QRect(4, 446, 341, 180))
        sizePolicy2.setHeightForWidth(self.TParts.sizePolicy().hasHeightForWidth())
        self.TParts.setSizePolicy(sizePolicy2)
        self.TParts.setMaximumSize(QSize(10000, 180))
        self.TParts.setFont(font1)
        self.TParts.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TParts.setRowCount(0)
        self.BDeleteGeometry = QPushButton(self.ImportPartTool)
        self.BDeleteGeometry.setObjectName(u"BDeleteGeometry")
        self.BDeleteGeometry.setGeometry(QRect(4, 632, 341, 22))
        sizePolicy8.setHeightForWidth(self.BDeleteGeometry.sizePolicy().hasHeightForWidth())
        self.BDeleteGeometry.setSizePolicy(sizePolicy8)
        self.BDeleteGeometry.setMaximumSize(QSize(450, 22))
        self.GeometryOptions.addWidget(self.ImportPartTool)
        self.ImportImageStackTool = QWidget()
        self.ImportImageStackTool.setObjectName(u"ImportImageStackTool")
        self.LImageFileFormat = QLabel(self.ImportImageStackTool)
        self.LImageFileFormat.setObjectName(u"LImageFileFormat")
        self.LImageFileFormat.setGeometry(QRect(36, 60, 245, 25))
        font10 = QFont()
        font10.setPointSize(12)
        self.LImageFileFormat.setFont(font10)
        self.LImageFileFormat.setTextFormat(Qt.TextFormat.RichText)
        self.LImageStack = QLabel(self.ImportImageStackTool)
        self.LImageStack.setObjectName(u"LImageStack")
        self.LImageStack.setGeometry(QRect(28, 5, 300, 25))
        font11 = QFont()
        font11.setPointSize(15)
        font11.setBold(True)
        font11.setUnderline(False)
        self.LImageStack.setFont(font11)
        self.LImageStack.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.LImageStack.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.layoutWidget3 = QWidget(self.ImportImageStackTool)
        self.layoutWidget3.setObjectName(u"layoutWidget3")
        self.layoutWidget3.setGeometry(QRect(36, 92, 277, 24))
        self.horizontalLayout = QHBoxLayout(self.layoutWidget3)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.BImageToVTK = QPushButton(self.layoutWidget3)
        self.BImageToVTK.setObjectName(u"BImageToVTK")
        self.BImageToVTK.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.BImageToVTK.sizePolicy().hasHeightForWidth())
        self.BImageToVTK.setSizePolicy(sizePolicy3)

        self.horizontalLayout.addWidget(self.BImageToVTK)

        self.BTiffToStl = QPushButton(self.layoutWidget3)
        self.BTiffToStl.setObjectName(u"BTiffToStl")
        self.BTiffToStl.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.BTiffToStl.sizePolicy().hasHeightForWidth())
        self.BTiffToStl.setSizePolicy(sizePolicy3)

        self.horizontalLayout.addWidget(self.BTiffToStl)

        self.LUploadedDirectory = QLabel(self.ImportImageStackTool)
        self.LUploadedDirectory.setObjectName(u"LUploadedDirectory")
        self.LUploadedDirectory.setGeometry(QRect(36, 33, 245, 25))
        self.LUploadedDirectory.setFont(font10)
        self.LUploadedDirectory.setTextFormat(Qt.TextFormat.RichText)
        self.GeometryOptions.addWidget(self.ImportImageStackTool)
        self.page_13 = QWidget()
        self.page_13.setObjectName(u"page_13")
        self.LImageStack_2 = QLabel(self.page_13)
        self.LImageStack_2.setObjectName(u"LImageStack_2")
        self.LImageStack_2.setGeometry(QRect(28, 13, 300, 25))
        self.LImageStack_2.setFont(font11)
        self.LImageStack_2.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.LImageStack_2.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.GeometryOptions.addWidget(self.page_13)
        self.page = QWidget()
        self.page.setObjectName(u"page")
        self.BDeleteBasicGeometry = QPushButton(self.page)
        self.BDeleteBasicGeometry.setObjectName(u"BDeleteBasicGeometry")
        self.BDeleteBasicGeometry.setGeometry(QRect(12, 550, 336, 22))
        self.BasicGeometries = QStackedWidget(self.page)
        self.BasicGeometries.setObjectName(u"BasicGeometries")
        self.BasicGeometries.setGeometry(QRect(12, 93, 336, 218))
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
        self.frame_10.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_10.setFrameShadow(QFrame.Shadow.Raised)
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
        self.frame_11.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_11.setFrameShadow(QFrame.Shadow.Raised)
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
        self.frame_12.setFrameShadow(QFrame.Shadow.Raised)
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
        self.frame_6 = QFrame(self.page)
        self.frame_6.setObjectName(u"frame_6")
        self.frame_6.setGeometry(QRect(12, 37, 336, 50))
        self.frame_6.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_6.setFrameShadow(QFrame.Shadow.Raised)
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

        self.TBasicGeometries = QTableWidget(self.page)
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
        self.TBasicGeometries.setGeometry(QRect(12, 345, 336, 199))
        self.LBasicGeometry = QLabel(self.page)
        self.LBasicGeometry.setObjectName(u"LBasicGeometry")
        self.LBasicGeometry.setGeometry(QRect(28, 5, 300, 25))
        self.BGenerateBasicGeometry = QPushButton(self.page)
        self.BGenerateBasicGeometry.setObjectName(u"BGenerateBasicGeometry")
        self.BGenerateBasicGeometry.setGeometry(QRect(12, 317, 336, 22))
        self.GeometryOptions.addWidget(self.page)

        self.verticalLayout_40.addWidget(self.GeometryOptions)

        self.SAGeometryScrollArea.setWidget(self.scrollAreaWidgetContents)

        self.verticalLayout_15.addWidget(self.SAGeometryScrollArea)

        self.ToolWindow.addWidget(self.ImportGeometryTool)
        self.GenerateMeshTool = QWidget()
        self.GenerateMeshTool.setObjectName(u"GenerateMeshTool")
        self.layoutWidget4 = QWidget(self.GenerateMeshTool)
        self.layoutWidget4.setObjectName(u"layoutWidget4")
        self.layoutWidget4.setGeometry(QRect(16, 6, 369, 1346))
        self.verticalLayout = QVBoxLayout(self.layoutWidget4)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.LDefineGlobalMesh = QLabel(self.layoutWidget4)
        self.LDefineGlobalMesh.setObjectName(u"LDefineGlobalMesh")
        sizePolicy1.setHeightForWidth(self.LDefineGlobalMesh.sizePolicy().hasHeightForWidth())
        self.LDefineGlobalMesh.setSizePolicy(sizePolicy1)

        self.verticalLayout.addWidget(self.LDefineGlobalMesh)

        self.verticalLayout_31 = QVBoxLayout()
        self.verticalLayout_31.setObjectName(u"verticalLayout_31")
        self.verticalLayout_31.setContentsMargins(6, -1, 6, -1)
        self.MeshInputs = QFrame(self.layoutWidget4)
        self.MeshInputs.setObjectName(u"MeshInputs")
        self.MeshInputs.setFrameShape(QFrame.Shape.NoFrame)
        self.MeshInputs.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_9 = QFormLayout(self.MeshInputs)
        self.formLayout_9.setObjectName(u"formLayout_9")
        self.formLayout_9.setHorizontalSpacing(0)
        self.formLayout_9.setContentsMargins(0, 0, 0, 0)
        self.formLayout_11 = QFormLayout()
        self.formLayout_11.setObjectName(u"formLayout_11")
        self.formLayout_11.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        self.formLayout_11.setLabelAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignTop)
        self.formLayout_11.setHorizontalSpacing(6)
        self.formLayout_11.setVerticalSpacing(6)
        self.formLayout_11.setContentsMargins(0, -1, 0, 6)
        self.LElementType = QLabel(self.MeshInputs)
        self.LElementType.setObjectName(u"LElementType")
        sizePolicy9 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        sizePolicy9.setHorizontalStretch(0)
        sizePolicy9.setVerticalStretch(0)
        sizePolicy9.setHeightForWidth(self.LElementType.sizePolicy().hasHeightForWidth())
        self.LElementType.setSizePolicy(sizePolicy9)

        self.formLayout_11.setWidget(1, QFormLayout.LabelRole, self.LElementType)

        self.INElementType = QComboBox(self.MeshInputs)
        self.INElementType.addItem("")
        self.INElementType.addItem("")
        self.INElementType.addItem("")
        self.INElementType.setObjectName(u"INElementType")
        sizePolicy8.setHeightForWidth(self.INElementType.sizePolicy().hasHeightForWidth())
        self.INElementType.setSizePolicy(sizePolicy8)

        self.formLayout_11.setWidget(1, QFormLayout.FieldRole, self.INElementType)

        self.LCoordinateSystem = QLabel(self.MeshInputs)
        self.LCoordinateSystem.setObjectName(u"LCoordinateSystem")
        sizePolicy9.setHeightForWidth(self.LCoordinateSystem.sizePolicy().hasHeightForWidth())
        self.LCoordinateSystem.setSizePolicy(sizePolicy9)

        self.formLayout_11.setWidget(2, QFormLayout.LabelRole, self.LCoordinateSystem)

        self.INCoordinateSystem = QComboBox(self.MeshInputs)
        self.INCoordinateSystem.addItem("")
        self.INCoordinateSystem.addItem("")
        self.INCoordinateSystem.setObjectName(u"INCoordinateSystem")
        sizePolicy8.setHeightForWidth(self.INCoordinateSystem.sizePolicy().hasHeightForWidth())
        self.INCoordinateSystem.setSizePolicy(sizePolicy8)

        self.formLayout_11.setWidget(2, QFormLayout.FieldRole, self.INCoordinateSystem)

        self.LDimension = QLabel(self.MeshInputs)
        self.LDimension.setObjectName(u"LDimension")
        sizePolicy9.setHeightForWidth(self.LDimension.sizePolicy().hasHeightForWidth())
        self.LDimension.setSizePolicy(sizePolicy9)

        self.formLayout_11.setWidget(3, QFormLayout.LabelRole, self.LDimension)

        self.INDimension = QComboBox(self.MeshInputs)
        self.INDimension.addItem("")
        self.INDimension.addItem("")
        self.INDimension.setObjectName(u"INDimension")
        sizePolicy8.setHeightForWidth(self.INDimension.sizePolicy().hasHeightForWidth())
        self.INDimension.setSizePolicy(sizePolicy8)

        self.formLayout_11.setWidget(3, QFormLayout.FieldRole, self.INDimension)

        self.line = QFrame(self.MeshInputs)
        self.line.setObjectName(u"line")
        self.line.setFrameShape(QFrame.Shape.HLine)
        self.line.setFrameShadow(QFrame.Shadow.Sunken)

        self.formLayout_11.setWidget(4, QFormLayout.SpanningRole, self.line)


        self.formLayout_9.setLayout(0, QFormLayout.FieldRole, self.formLayout_11)


        self.verticalLayout_31.addWidget(self.MeshInputs)

        self.MeshInputs2 = QStackedWidget(self.layoutWidget4)
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
        self.Cylindrical2DInputs2.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_13 = QHBoxLayout(self.Cylindrical2DInputs2)
        self.horizontalLayout_13.setObjectName(u"horizontalLayout_13")
        self.horizontalLayout_13.setContentsMargins(-1, 0, -1, 0)
        self.LInnerRadiusC2D = QLabel(self.Cylindrical2DInputs2)
        self.LInnerRadiusC2D.setObjectName(u"LInnerRadiusC2D")

        self.horizontalLayout_13.addWidget(self.LInnerRadiusC2D)

        self.INInnerRadiusC2D = QLineEdit(self.Cylindrical2DInputs2)
        self.INInnerRadiusC2D.setObjectName(u"INInnerRadiusC2D")
        self.INInnerRadiusC2D.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

        self.horizontalLayout_13.addWidget(self.INInnerRadiusC2D)


        self.verticalLayout_32.addWidget(self.Cylindrical2DInputs2)

        self.Cylindrical2DInputs = QFrame(self.Cylindrical2D)
        self.Cylindrical2DInputs.setObjectName(u"Cylindrical2DInputs")
        self.Cylindrical2DInputs.setFrameShape(QFrame.Shape.NoFrame)
        self.Cylindrical2DInputs.setFrameShadow(QFrame.Shadow.Raised)
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
        self.INLengthThetaC2D.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

        self.gridLayout_5.addWidget(self.INLengthThetaC2D, 1, 4, 1, 1)

        self.label_38 = QLabel(self.Cylindrical2DInputs)
        self.label_38.setObjectName(u"label_38")

        self.gridLayout_5.addWidget(self.label_38, 2, 1, 1, 1)

        self.INLengthOutRadC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INLengthOutRadC2D.setObjectName(u"INLengthOutRadC2D")
        self.INLengthOutRadC2D.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

        self.gridLayout_5.addWidget(self.INLengthOutRadC2D, 1, 2, 1, 1)

        self.label_40 = QLabel(self.Cylindrical2DInputs)
        self.label_40.setObjectName(u"label_40")

        self.gridLayout_5.addWidget(self.label_40, 2, 3, 1, 1)

        self.INOriginYC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INOriginYC2D.setObjectName(u"INOriginYC2D")
        self.INOriginYC2D.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

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
        self.INElementsArcC2D.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

        self.gridLayout_5.addWidget(self.INElementsArcC2D, 2, 4, 1, 1)

        self.INOriginXC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INOriginXC2D.setObjectName(u"INOriginXC2D")
        self.INOriginXC2D.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

        self.gridLayout_5.addWidget(self.INOriginXC2D, 0, 2, 1, 1)

        self.INElementsRadialC2D = QLineEdit(self.Cylindrical2DInputs)
        self.INElementsRadialC2D.setObjectName(u"INElementsRadialC2D")
        self.INElementsRadialC2D.setInputMethodHints(Qt.InputMethodHint.ImhPreferNumbers)

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
        self.Cylindrical3DInputs2.setFrameShadow(QFrame.Shadow.Raised)
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
        self.Cylindrical3DInputs.setFrameShadow(QFrame.Shadow.Raised)
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

        self.verticalLayout_31.addWidget(self.MeshInputs2, 0, Qt.AlignmentFlag.AlignTop)

        self.BGenerateGlobalMesh = QPushButton(self.layoutWidget4)
        self.BGenerateGlobalMesh.setObjectName(u"BGenerateGlobalMesh")

        self.verticalLayout_31.addWidget(self.BGenerateGlobalMesh)

        self.verticalSpacer_8 = QSpacerItem(20, 1000, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_31.addItem(self.verticalSpacer_8)


        self.verticalLayout.addLayout(self.verticalLayout_31)

        self.ToolWindow.addWidget(self.GenerateMeshTool)
        self.SolverSettingsSGHTool = QWidget()
        self.SolverSettingsSGHTool.setObjectName(u"SolverSettingsSGHTool")
        self.verticalLayout_28 = QVBoxLayout(self.SolverSettingsSGHTool)
        self.verticalLayout_28.setObjectName(u"verticalLayout_28")
        self.verticalLayout_28.setContentsMargins(10, 7, -1, -1)
        self.LSolverSettings = QLabel(self.SolverSettingsSGHTool)
        self.LSolverSettings.setObjectName(u"LSolverSettings")
        sizePolicy1.setHeightForWidth(self.LSolverSettings.sizePolicy().hasHeightForWidth())
        self.LSolverSettings.setSizePolicy(sizePolicy1)

        self.verticalLayout_28.addWidget(self.LSolverSettings)

        self.INSelectSolverSettings = QComboBox(self.SolverSettingsSGHTool)
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
        self.solversettings.setFrameShadow(QFrame.Shadow.Raised)
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
        self.LNumberOfSteps = QLabel(self.page_3)
        self.LNumberOfSteps.setObjectName(u"LNumberOfSteps")
        self.LNumberOfSteps.setGeometry(QRect(13, 20, 108, 22))
        self.INNumberOfSteps = QLineEdit(self.page_3)
        self.INNumberOfSteps.setObjectName(u"INNumberOfSteps")
        self.INNumberOfSteps.setEnabled(False)
        self.INNumberOfSteps.setGeometry(QRect(124, 20, 250, 22))
        self.INNumberOfSteps.setReadOnly(False)
        self.SolverSettingsOptions.addWidget(self.page_3)

        self.verticalLayout_28.addWidget(self.SolverSettingsOptions)

        self.ToolWindow.addWidget(self.SolverSettingsSGHTool)
        self.BoundaryConditionsTool = QWidget()
        self.BoundaryConditionsTool.setObjectName(u"BoundaryConditionsTool")
        self.verticalLayout_17 = QVBoxLayout(self.BoundaryConditionsTool)
        self.verticalLayout_17.setObjectName(u"verticalLayout_17")
        self.verticalLayout_17.setContentsMargins(-1, 6, -1, -1)
        self.LBoundaryConditions = QLabel(self.BoundaryConditionsTool)
        self.LBoundaryConditions.setObjectName(u"LBoundaryConditions")
        sizePolicy1.setHeightForWidth(self.LBoundaryConditions.sizePolicy().hasHeightForWidth())
        self.LBoundaryConditions.setSizePolicy(sizePolicy1)

        self.verticalLayout_17.addWidget(self.LBoundaryConditions)

        self.INSelectBoundaryConditions = QComboBox(self.BoundaryConditionsTool)
        self.INSelectBoundaryConditions.addItem("")
        self.INSelectBoundaryConditions.addItem("")
        self.INSelectBoundaryConditions.setObjectName(u"INSelectBoundaryConditions")

        self.verticalLayout_17.addWidget(self.INSelectBoundaryConditions)

        self.BoundaryConditionsOptions = QStackedWidget(self.BoundaryConditionsTool)
        self.BoundaryConditionsOptions.setObjectName(u"BoundaryConditionsOptions")
        self.page_8 = QWidget()
        self.page_8.setObjectName(u"page_8")
        self.bcsettings = QFrame(self.page_8)
        self.bcsettings.setObjectName(u"bcsettings")
        self.bcsettings.setGeometry(QRect(4, 20, 376, 190))
        self.bcsettings.setFrameShape(QFrame.Shape.NoFrame)
        self.bcsettings.setFrameShadow(QFrame.Shadow.Raised)
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

        self.BaddBC = QPushButton(self.page_8)
        self.BaddBC.setObjectName(u"BaddBC")
        self.BaddBC.setGeometry(QRect(4, 216, 376, 22))
        self.BdeleteBC = QPushButton(self.page_8)
        self.BdeleteBC.setObjectName(u"BdeleteBC")
        self.BdeleteBC.setGeometry(QRect(4, 501, 376, 22))
        self.TBoundaryConditions = QTableWidget(self.page_8)
        if (self.TBoundaryConditions.columnCount() < 7):
            self.TBoundaryConditions.setColumnCount(7)
        __qtablewidgetitem23 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(0, __qtablewidgetitem23)
        __qtablewidgetitem24 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(1, __qtablewidgetitem24)
        __qtablewidgetitem25 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(2, __qtablewidgetitem25)
        __qtablewidgetitem26 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(3, __qtablewidgetitem26)
        __qtablewidgetitem27 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(4, __qtablewidgetitem27)
        __qtablewidgetitem28 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(5, __qtablewidgetitem28)
        __qtablewidgetitem29 = QTableWidgetItem()
        self.TBoundaryConditions.setHorizontalHeaderItem(6, __qtablewidgetitem29)
        self.TBoundaryConditions.setObjectName(u"TBoundaryConditions")
        self.TBoundaryConditions.setGeometry(QRect(4, 244, 376, 251))
        self.BoundaryConditionsOptions.addWidget(self.page_8)
        self.page_6 = QWidget()
        self.page_6.setObjectName(u"page_6")
        self.layoutWidget5 = QWidget(self.page_6)
        self.layoutWidget5.setObjectName(u"layoutWidget5")
        self.layoutWidget5.setGeometry(QRect(0, 4, 380, 341))
        self.verticalLayout_4 = QVBoxLayout(self.layoutWidget5)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.BoundaryConditionsInputs = QFrame(self.layoutWidget5)
        self.BoundaryConditionsInputs.setObjectName(u"BoundaryConditionsInputs")
        self.BoundaryConditionsInputs.setFrameShape(QFrame.Shape.NoFrame)
        self.BoundaryConditionsInputs.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout_3 = QFormLayout(self.BoundaryConditionsInputs)
        self.formLayout_3.setObjectName(u"formLayout_3")
        self.formLayout_3.setContentsMargins(0, -1, 0, -1)
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


        self.verticalLayout_4.addWidget(self.BoundaryConditionsInputs)

        self.BAddBC = QPushButton(self.layoutWidget5)
        self.BAddBC.setObjectName(u"BAddBC")

        self.verticalLayout_4.addWidget(self.BAddBC)

        self.TBCs = QTableWidget(self.layoutWidget5)
        if (self.TBCs.columnCount() < 2):
            self.TBCs.setColumnCount(2)
        __qtablewidgetitem30 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(0, __qtablewidgetitem30)
        __qtablewidgetitem31 = QTableWidgetItem()
        self.TBCs.setHorizontalHeaderItem(1, __qtablewidgetitem31)
        self.TBCs.setObjectName(u"TBCs")
        self.TBCs.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TBCs.horizontalHeader().setDefaultSectionSize(175)
        self.TBCs.horizontalHeader().setStretchLastSection(True)

        self.verticalLayout_4.addWidget(self.TBCs)

        self.BDeleteBC = QPushButton(self.layoutWidget5)
        self.BDeleteBC.setObjectName(u"BDeleteBC")

        self.verticalLayout_4.addWidget(self.BDeleteBC)

        self.BoundaryConditionsOptions.addWidget(self.page_6)

        self.verticalLayout_17.addWidget(self.BoundaryConditionsOptions)

        self.ToolWindow.addWidget(self.BoundaryConditionsTool)
        self.MaterialsTool = QWidget()
        self.MaterialsTool.setObjectName(u"MaterialsTool")
        self.verticalLayout_13 = QVBoxLayout(self.MaterialsTool)
        self.verticalLayout_13.setObjectName(u"verticalLayout_13")
        self.LDefineMaterials = QLabel(self.MaterialsTool)
        self.LDefineMaterials.setObjectName(u"LDefineMaterials")
        sizePolicy1.setHeightForWidth(self.LDefineMaterials.sizePolicy().hasHeightForWidth())
        self.LDefineMaterials.setSizePolicy(sizePolicy1)

        self.verticalLayout_13.addWidget(self.LDefineMaterials)

        self.tabWidget = QTabWidget(self.MaterialsTool)
        self.tabWidget.setObjectName(u"tabWidget")
        sizePolicy3.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy3)
        self.tabWidget.setIconSize(QSize(20, 20))
        self.tabWidget.setElideMode(Qt.TextElideMode.ElideNone)
        self.DefineMaterials = QWidget()
        self.DefineMaterials.setObjectName(u"DefineMaterials")
        self.verticalLayoutWidget = QWidget(self.DefineMaterials)
        self.verticalLayoutWidget.setObjectName(u"verticalLayoutWidget")
        self.verticalLayoutWidget.setGeometry(QRect(0, 0, 369, 692))
        self.verticalLayout_12 = QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout_12.setSpacing(0)
        self.verticalLayout_12.setObjectName(u"verticalLayout_12")
        self.verticalLayout_12.setSizeConstraint(QLayout.SizeConstraint.SetDefaultConstraint)
        self.verticalLayout_12.setContentsMargins(6, 6, 6, 0)
        self.INSelectDefineMaterials = QComboBox(self.verticalLayoutWidget)
        self.INSelectDefineMaterials.addItem("")
        self.INSelectDefineMaterials.addItem("")
        self.INSelectDefineMaterials.setObjectName(u"INSelectDefineMaterials")

        self.verticalLayout_12.addWidget(self.INSelectDefineMaterials)

        self.DefineMaterialsOptions = QStackedWidget(self.verticalLayoutWidget)
        self.DefineMaterialsOptions.setObjectName(u"DefineMaterialsOptions")
        self.DefineMaterialsOptions.setMinimumSize(QSize(0, 600))
        self.DefineMaterialsOptions.setMaximumSize(QSize(16777215, 16777215))
        self.page_4 = QWidget()
        self.page_4.setObjectName(u"page_4")
        self.BAddMaterialSGH = QPushButton(self.page_4)
        self.BAddMaterialSGH.setObjectName(u"BAddMaterialSGH")
        self.BAddMaterialSGH.setGeometry(QRect(3, 96, 351, 22))
        self.frame_16 = QFrame(self.page_4)
        self.frame_16.setObjectName(u"frame_16")
        self.frame_16.setGeometry(QRect(3, 348, 351, 241))
        self.frame_16.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_16.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_10 = QGridLayout(self.frame_16)
        self.gridLayout_10.setObjectName(u"gridLayout_10")
        self.gridLayout_10.setVerticalSpacing(12)
        self.gridLayout_10.setContentsMargins(0, 0, 0, 0)
        self.INq1 = QLineEdit(self.frame_16)
        self.INq1.setObjectName(u"INq1")
        self.INq1.setEnabled(True)
        sizePolicy.setHeightForWidth(self.INq1.sizePolicy().hasHeightForWidth())
        self.INq1.setSizePolicy(sizePolicy)

        self.gridLayout_10.addWidget(self.INq1, 0, 1, 1, 1)

        self.INSpecificHeat = QLineEdit(self.frame_16)
        self.INSpecificHeat.setObjectName(u"INSpecificHeat")
        self.INSpecificHeat.setEnabled(True)
        sizePolicy.setHeightForWidth(self.INSpecificHeat.sizePolicy().hasHeightForWidth())
        self.INSpecificHeat.setSizePolicy(sizePolicy)

        self.gridLayout_10.addWidget(self.INSpecificHeat, 6, 1, 1, 1)

        self.Lq1ex = QLabel(self.frame_16)
        self.Lq1ex.setObjectName(u"Lq1ex")
        self.Lq1ex.setEnabled(True)
        self.Lq1ex.setWordWrap(True)

        self.gridLayout_10.addWidget(self.Lq1ex, 2, 0, 1, 1)

        self.LGamma = QLabel(self.frame_16)
        self.LGamma.setObjectName(u"LGamma")
        self.LGamma.setEnabled(True)

        self.gridLayout_10.addWidget(self.LGamma, 4, 0, 1, 1)

        self.INq1ex = QLineEdit(self.frame_16)
        self.INq1ex.setObjectName(u"INq1ex")
        self.INq1ex.setEnabled(True)
        sizePolicy.setHeightForWidth(self.INq1ex.sizePolicy().hasHeightForWidth())
        self.INq1ex.setSizePolicy(sizePolicy)

        self.gridLayout_10.addWidget(self.INq1ex, 2, 1, 1, 1)

        self.INq2 = QLineEdit(self.frame_16)
        self.INq2.setObjectName(u"INq2")
        self.INq2.setEnabled(True)
        sizePolicy.setHeightForWidth(self.INq2.sizePolicy().hasHeightForWidth())
        self.INq2.setSizePolicy(sizePolicy)

        self.gridLayout_10.addWidget(self.INq2, 1, 1, 1, 1)

        self.INq2ex = QLineEdit(self.frame_16)
        self.INq2ex.setObjectName(u"INq2ex")
        self.INq2ex.setEnabled(True)
        sizePolicy.setHeightForWidth(self.INq2ex.sizePolicy().hasHeightForWidth())
        self.INq2ex.setSizePolicy(sizePolicy)

        self.gridLayout_10.addWidget(self.INq2ex, 3, 1, 1, 1)

        self.LMinSound = QLabel(self.frame_16)
        self.LMinSound.setObjectName(u"LMinSound")
        self.LMinSound.setEnabled(True)

        self.gridLayout_10.addWidget(self.LMinSound, 5, 0, 1, 1)

        self.Lq2ex = QLabel(self.frame_16)
        self.Lq2ex.setObjectName(u"Lq2ex")
        self.Lq2ex.setEnabled(True)
        self.Lq2ex.setWordWrap(True)

        self.gridLayout_10.addWidget(self.Lq2ex, 3, 0, 1, 1)

        self.Lq1 = QLabel(self.frame_16)
        self.Lq1.setObjectName(u"Lq1")
        self.Lq1.setEnabled(True)
        self.Lq1.setWordWrap(True)

        self.gridLayout_10.addWidget(self.Lq1, 0, 0, 1, 1)

        self.LSpecificHeat = QLabel(self.frame_16)
        self.LSpecificHeat.setObjectName(u"LSpecificHeat")
        self.LSpecificHeat.setEnabled(True)

        self.gridLayout_10.addWidget(self.LSpecificHeat, 6, 0, 1, 1)

        self.INGamma = QLineEdit(self.frame_16)
        self.INGamma.setObjectName(u"INGamma")
        self.INGamma.setEnabled(True)
        sizePolicy.setHeightForWidth(self.INGamma.sizePolicy().hasHeightForWidth())
        self.INGamma.setSizePolicy(sizePolicy)

        self.gridLayout_10.addWidget(self.INGamma, 4, 1, 1, 1)

        self.Lq2 = QLabel(self.frame_16)
        self.Lq2.setObjectName(u"Lq2")
        self.Lq2.setEnabled(True)
        self.Lq2.setWordWrap(True)

        self.gridLayout_10.addWidget(self.Lq2, 1, 0, 1, 1)

        self.INMinSound = QLineEdit(self.frame_16)
        self.INMinSound.setObjectName(u"INMinSound")
        self.INMinSound.setEnabled(True)
        sizePolicy.setHeightForWidth(self.INMinSound.sizePolicy().hasHeightForWidth())
        self.INMinSound.setSizePolicy(sizePolicy)

        self.gridLayout_10.addWidget(self.INMinSound, 5, 1, 1, 1)

        self.BDeleteMaterialSGH = QPushButton(self.page_4)
        self.BDeleteMaterialSGH.setObjectName(u"BDeleteMaterialSGH")
        self.BDeleteMaterialSGH.setGeometry(QRect(3, 292, 351, 22))
        self.TMaterialsSGH = QTableWidget(self.page_4)
        if (self.TMaterialsSGH.columnCount() < 9):
            self.TMaterialsSGH.setColumnCount(9)
        __qtablewidgetitem32 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(0, __qtablewidgetitem32)
        __qtablewidgetitem33 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(1, __qtablewidgetitem33)
        __qtablewidgetitem34 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(2, __qtablewidgetitem34)
        __qtablewidgetitem35 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(3, __qtablewidgetitem35)
        __qtablewidgetitem36 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(4, __qtablewidgetitem36)
        __qtablewidgetitem37 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(5, __qtablewidgetitem37)
        __qtablewidgetitem38 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(6, __qtablewidgetitem38)
        __qtablewidgetitem39 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(7, __qtablewidgetitem39)
        __qtablewidgetitem40 = QTableWidgetItem()
        self.TMaterialsSGH.setHorizontalHeaderItem(8, __qtablewidgetitem40)
        self.TMaterialsSGH.setObjectName(u"TMaterialsSGH")
        self.TMaterialsSGH.setEnabled(True)
        self.TMaterialsSGH.setGeometry(QRect(3, 124, 351, 161))
        self.TMaterialsSGH.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TMaterialsSGH.setRowCount(0)
        self.label_12 = QLabel(self.page_4)
        self.label_12.setObjectName(u"label_12")
        self.label_12.setGeometry(QRect(4, 320, 353, 21))
        font12 = QFont()
        font12.setPointSize(11)
        self.label_12.setFont(font12)
        self.frame_7 = QFrame(self.page_4)
        self.frame_7.setObjectName(u"frame_7")
        self.frame_7.setGeometry(QRect(3, 8, 351, 78))
        self.frame_7.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_7.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_9 = QGridLayout(self.frame_7)
        self.gridLayout_9.setObjectName(u"gridLayout_9")
        self.gridLayout_9.setContentsMargins(0, 0, 0, 0)
        self.LMaterialNameSGH = QLabel(self.frame_7)
        self.LMaterialNameSGH.setObjectName(u"LMaterialNameSGH")

        self.gridLayout_9.addWidget(self.LMaterialNameSGH, 1, 0, 1, 1)

        self.INMaterialNameSGH = QLineEdit(self.frame_7)
        self.INMaterialNameSGH.setObjectName(u"INMaterialNameSGH")
        self.INMaterialNameSGH.setMinimumSize(QSize(93, 0))

        self.gridLayout_9.addWidget(self.INMaterialNameSGH, 1, 1, 1, 1)

        self.INArtificialViscosity = QComboBox(self.frame_7)
        self.INArtificialViscosity.addItem("")
        self.INArtificialViscosity.addItem("")
        self.INArtificialViscosity.setObjectName(u"INArtificialViscosity")

        self.gridLayout_9.addWidget(self.INArtificialViscosity, 3, 1, 1, 1)

        self.INEOS = QComboBox(self.frame_7)
        self.INEOS.addItem("")
        self.INEOS.setObjectName(u"INEOS")

        self.gridLayout_9.addWidget(self.INEOS, 2, 1, 1, 1)

        self.LEOS = QLabel(self.frame_7)
        self.LEOS.setObjectName(u"LEOS")

        self.gridLayout_9.addWidget(self.LEOS, 2, 0, 1, 1)

        self.LArtificialViscosity = QLabel(self.frame_7)
        self.LArtificialViscosity.setObjectName(u"LArtificialViscosity")
        self.LArtificialViscosity.setTextFormat(Qt.TextFormat.RichText)

        self.gridLayout_9.addWidget(self.LArtificialViscosity, 3, 0, 1, 1)

        self.DefineMaterialsOptions.addWidget(self.page_4)
        self.page_5 = QWidget()
        self.page_5.setObjectName(u"page_5")
        self.BAddMaterial = QPushButton(self.page_5)
        self.BAddMaterial.setObjectName(u"BAddMaterial")
        self.BAddMaterial.setGeometry(QRect(0, 302, 353, 22))
        self.MaterialTypeTool = QStackedWidget(self.page_5)
        self.MaterialTypeTool.setObjectName(u"MaterialTypeTool")
        self.MaterialTypeTool.setGeometry(QRect(4, 90, 350, 212))
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
        self.IsotropicPlane.setFrameShape(QFrame.Shape.NoFrame)
        self.IsotropicPlane.setFrameShadow(QFrame.Shadow.Raised)
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


        self.verticalLayout_23.addWidget(self.IsotropicPlane, 0, Qt.AlignmentFlag.AlignTop)

        self.TransverslyIsotropicMat = QFrame(self.TransverselyIsotropic)
        self.TransverslyIsotropicMat.setObjectName(u"TransverslyIsotropicMat")
        self.TransverslyIsotropicMat.setFrameShape(QFrame.Shape.NoFrame)
        self.TransverslyIsotropicMat.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_6 = QHBoxLayout(self.TransverslyIsotropicMat)
        self.horizontalLayout_6.setSpacing(0)
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.horizontalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.TransverseInPlane = QFrame(self.TransverslyIsotropicMat)
        self.TransverseInPlane.setObjectName(u"TransverseInPlane")
        self.TransverseInPlane.setFrameShape(QFrame.Shape.NoFrame)
        self.TransverseInPlane.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_25 = QVBoxLayout(self.TransverseInPlane)
        self.verticalLayout_25.setSpacing(0)
        self.verticalLayout_25.setObjectName(u"verticalLayout_25")
        self.verticalLayout_25.setContentsMargins(0, 12, 0, 0)
        self.LInPlane = QLabel(self.TransverseInPlane)
        self.LInPlane.setObjectName(u"LInPlane")

        self.verticalLayout_25.addWidget(self.LInPlane)

        self.TransverseInPlaneMat = QFrame(self.TransverseInPlane)
        self.TransverseInPlaneMat.setObjectName(u"TransverseInPlaneMat")
        self.TransverseInPlaneMat.setFrameShape(QFrame.Shape.NoFrame)
        self.TransverseInPlaneMat.setFrameShadow(QFrame.Shadow.Raised)
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


        self.horizontalLayout_6.addWidget(self.TransverseInPlane, 0, Qt.AlignmentFlag.AlignTop)

        self.TransverseOutOfPlane = QFrame(self.TransverslyIsotropicMat)
        self.TransverseOutOfPlane.setObjectName(u"TransverseOutOfPlane")
        self.TransverseOutOfPlane.setFrameShape(QFrame.Shape.NoFrame)
        self.TransverseOutOfPlane.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_24 = QVBoxLayout(self.TransverseOutOfPlane)
        self.verticalLayout_24.setSpacing(0)
        self.verticalLayout_24.setObjectName(u"verticalLayout_24")
        self.verticalLayout_24.setContentsMargins(0, 12, 0, 0)
        self.LOutOfPlane = QLabel(self.TransverseOutOfPlane)
        self.LOutOfPlane.setObjectName(u"LOutOfPlane")

        self.verticalLayout_24.addWidget(self.LOutOfPlane)

        self.TransverseOutOfPlaneMat = QFrame(self.TransverseOutOfPlane)
        self.TransverseOutOfPlaneMat.setObjectName(u"TransverseOutOfPlaneMat")
        self.TransverseOutOfPlaneMat.setFrameShape(QFrame.Shape.NoFrame)
        self.TransverseOutOfPlaneMat.setFrameShadow(QFrame.Shadow.Raised)
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
        __qtablewidgetitem41 = QTableWidgetItem()
        self.TAnisotropic.setItem(0, 0, __qtablewidgetitem41)
        brush = QBrush(QColor(235, 235, 235, 255))
        brush.setStyle(Qt.SolidPattern)
        __qtablewidgetitem42 = QTableWidgetItem()
        __qtablewidgetitem42.setBackground(brush);
        __qtablewidgetitem42.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(1, 0, __qtablewidgetitem42)
        __qtablewidgetitem43 = QTableWidgetItem()
        __qtablewidgetitem43.setBackground(brush);
        __qtablewidgetitem43.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(2, 0, __qtablewidgetitem43)
        __qtablewidgetitem44 = QTableWidgetItem()
        __qtablewidgetitem44.setBackground(brush);
        __qtablewidgetitem44.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(2, 1, __qtablewidgetitem44)
        __qtablewidgetitem45 = QTableWidgetItem()
        __qtablewidgetitem45.setBackground(brush);
        __qtablewidgetitem45.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 0, __qtablewidgetitem45)
        __qtablewidgetitem46 = QTableWidgetItem()
        __qtablewidgetitem46.setBackground(brush);
        __qtablewidgetitem46.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 1, __qtablewidgetitem46)
        __qtablewidgetitem47 = QTableWidgetItem()
        __qtablewidgetitem47.setBackground(brush);
        __qtablewidgetitem47.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(3, 2, __qtablewidgetitem47)
        __qtablewidgetitem48 = QTableWidgetItem()
        __qtablewidgetitem48.setBackground(brush);
        __qtablewidgetitem48.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 0, __qtablewidgetitem48)
        __qtablewidgetitem49 = QTableWidgetItem()
        __qtablewidgetitem49.setBackground(brush);
        __qtablewidgetitem49.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 1, __qtablewidgetitem49)
        __qtablewidgetitem50 = QTableWidgetItem()
        __qtablewidgetitem50.setBackground(brush);
        __qtablewidgetitem50.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 2, __qtablewidgetitem50)
        __qtablewidgetitem51 = QTableWidgetItem()
        __qtablewidgetitem51.setBackground(brush);
        __qtablewidgetitem51.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(4, 3, __qtablewidgetitem51)
        __qtablewidgetitem52 = QTableWidgetItem()
        __qtablewidgetitem52.setBackground(brush);
        __qtablewidgetitem52.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 0, __qtablewidgetitem52)
        __qtablewidgetitem53 = QTableWidgetItem()
        __qtablewidgetitem53.setBackground(brush);
        __qtablewidgetitem53.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 1, __qtablewidgetitem53)
        __qtablewidgetitem54 = QTableWidgetItem()
        __qtablewidgetitem54.setBackground(brush);
        __qtablewidgetitem54.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 2, __qtablewidgetitem54)
        __qtablewidgetitem55 = QTableWidgetItem()
        __qtablewidgetitem55.setBackground(brush);
        __qtablewidgetitem55.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 3, __qtablewidgetitem55)
        __qtablewidgetitem56 = QTableWidgetItem()
        __qtablewidgetitem56.setBackground(brush);
        __qtablewidgetitem56.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEditable|Qt.ItemIsDragEnabled|Qt.ItemIsDropEnabled|Qt.ItemIsUserCheckable);
        self.TAnisotropic.setItem(5, 4, __qtablewidgetitem56)
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
        self.gridLayout_2.setHorizontalSpacing(6)
        self.gridLayout_2.setVerticalSpacing(0)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.INNUyz = QLineEdit(self.Orthotropic)
        self.INNUyz.setObjectName(u"INNUyz")

        self.gridLayout_2.addWidget(self.INNUyz, 2, 8, 1, 1)

        self.LGxy = QLabel(self.Orthotropic)
        self.LGxy.setObjectName(u"LGxy")

        self.gridLayout_2.addWidget(self.LGxy, 3, 1, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.INGxy = QLineEdit(self.Orthotropic)
        self.INGxy.setObjectName(u"INGxy")

        self.gridLayout_2.addWidget(self.INGxy, 3, 2, 1, 1)

        self.INEy = QLineEdit(self.Orthotropic)
        self.INEy.setObjectName(u"INEy")

        self.gridLayout_2.addWidget(self.INEy, 0, 5, 1, 1)

        self.LEx = QLabel(self.Orthotropic)
        self.LEx.setObjectName(u"LEx")

        self.gridLayout_2.addWidget(self.LEx, 0, 1, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.LNUyz = QLabel(self.Orthotropic)
        self.LNUyz.setObjectName(u"LNUyz")

        self.gridLayout_2.addWidget(self.LNUyz, 2, 7, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.LGyz = QLabel(self.Orthotropic)
        self.LGyz.setObjectName(u"LGyz")

        self.gridLayout_2.addWidget(self.LGyz, 3, 7, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.LEy = QLabel(self.Orthotropic)
        self.LEy.setObjectName(u"LEy")

        self.gridLayout_2.addWidget(self.LEy, 0, 4, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.LEz = QLabel(self.Orthotropic)
        self.LEz.setObjectName(u"LEz")

        self.gridLayout_2.addWidget(self.LEz, 0, 7, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.LGxz = QLabel(self.Orthotropic)
        self.LGxz.setObjectName(u"LGxz")

        self.gridLayout_2.addWidget(self.LGxz, 3, 4, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.LNUxz = QLabel(self.Orthotropic)
        self.LNUxz.setObjectName(u"LNUxz")

        self.gridLayout_2.addWidget(self.LNUxz, 2, 4, 1, 1, Qt.AlignmentFlag.AlignRight)

        self.INNUxy = QLineEdit(self.Orthotropic)
        self.INNUxy.setObjectName(u"INNUxy")

        self.gridLayout_2.addWidget(self.INNUxy, 2, 2, 1, 1)

        self.INEz = QLineEdit(self.Orthotropic)
        self.INEz.setObjectName(u"INEz")

        self.gridLayout_2.addWidget(self.INEz, 0, 8, 1, 1)

        self.LNUxy = QLabel(self.Orthotropic)
        self.LNUxy.setObjectName(u"LNUxy")

        self.gridLayout_2.addWidget(self.LNUxy, 2, 1, 1, 1, Qt.AlignmentFlag.AlignRight)

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
        self.TMaterials = QTableWidget(self.page_5)
        if (self.TMaterials.columnCount() < 24):
            self.TMaterials.setColumnCount(24)
        __qtablewidgetitem57 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(0, __qtablewidgetitem57)
        __qtablewidgetitem58 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(1, __qtablewidgetitem58)
        __qtablewidgetitem59 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(2, __qtablewidgetitem59)
        __qtablewidgetitem60 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(3, __qtablewidgetitem60)
        __qtablewidgetitem61 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(4, __qtablewidgetitem61)
        __qtablewidgetitem62 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(5, __qtablewidgetitem62)
        __qtablewidgetitem63 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(6, __qtablewidgetitem63)
        __qtablewidgetitem64 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(7, __qtablewidgetitem64)
        __qtablewidgetitem65 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(8, __qtablewidgetitem65)
        __qtablewidgetitem66 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(9, __qtablewidgetitem66)
        __qtablewidgetitem67 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(10, __qtablewidgetitem67)
        __qtablewidgetitem68 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(11, __qtablewidgetitem68)
        __qtablewidgetitem69 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(12, __qtablewidgetitem69)
        __qtablewidgetitem70 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(13, __qtablewidgetitem70)
        __qtablewidgetitem71 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(14, __qtablewidgetitem71)
        __qtablewidgetitem72 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(15, __qtablewidgetitem72)
        __qtablewidgetitem73 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(16, __qtablewidgetitem73)
        __qtablewidgetitem74 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(17, __qtablewidgetitem74)
        __qtablewidgetitem75 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(18, __qtablewidgetitem75)
        __qtablewidgetitem76 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(19, __qtablewidgetitem76)
        __qtablewidgetitem77 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(20, __qtablewidgetitem77)
        __qtablewidgetitem78 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(21, __qtablewidgetitem78)
        __qtablewidgetitem79 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(22, __qtablewidgetitem79)
        __qtablewidgetitem80 = QTableWidgetItem()
        self.TMaterials.setHorizontalHeaderItem(23, __qtablewidgetitem80)
        self.TMaterials.setObjectName(u"TMaterials")
        self.TMaterials.setEnabled(True)
        self.TMaterials.setGeometry(QRect(0, 330, 353, 192))
        self.TMaterials.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TMaterials.setRowCount(0)
        self.BDeleteMaterial = QPushButton(self.page_5)
        self.BDeleteMaterial.setObjectName(u"BDeleteMaterial")
        self.BDeleteMaterial.setGeometry(QRect(0, 528, 353, 22))
        self.MaterialInputs = QFrame(self.page_5)
        self.MaterialInputs.setObjectName(u"MaterialInputs")
        self.MaterialInputs.setGeometry(QRect(3, 8, 357, 78))
        self.MaterialInputs.setFrameShape(QFrame.Shape.NoFrame)
        self.MaterialInputs.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout = QGridLayout(self.MaterialInputs)
        self.gridLayout.setObjectName(u"gridLayout")
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.LMaterialName = QLabel(self.MaterialInputs)
        self.LMaterialName.setObjectName(u"LMaterialName")

        self.gridLayout.addWidget(self.LMaterialName, 0, 0, 1, 1)

        self.LType = QLabel(self.MaterialInputs)
        self.LType.setObjectName(u"LType")

        self.gridLayout.addWidget(self.LType, 2, 0, 1, 1)

        self.frame_2 = QFrame(self.MaterialInputs)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_2.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_8 = QHBoxLayout(self.frame_2)
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.horizontalLayout_8.setContentsMargins(0, 0, 0, 0)

        self.gridLayout.addWidget(self.frame_2, 1, 4, 1, 1, Qt.AlignmentFlag.AlignLeft)

        self.LRegion = QLabel(self.MaterialInputs)
        self.LRegion.setObjectName(u"LRegion")

        self.gridLayout.addWidget(self.LRegion, 1, 0, 1, 1)

        self.frame = QFrame(self.MaterialInputs)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.Shape.NoFrame)
        self.frame.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_7 = QHBoxLayout(self.frame)
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.horizontalLayout_7.setContentsMargins(0, 0, 0, 0)

        self.gridLayout.addWidget(self.frame, 2, 4, 1, 1, Qt.AlignmentFlag.AlignLeft)

        self.frame_3 = QFrame(self.MaterialInputs)
        self.frame_3.setObjectName(u"frame_3")
        self.frame_3.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_3.setFrameShadow(QFrame.Shadow.Raised)
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

        self.gridLayout.addWidget(self.INRegion, 1, 1, 1, 1, Qt.AlignmentFlag.AlignLeft)

        self.INMaterialName = QLineEdit(self.MaterialInputs)
        self.INMaterialName.setObjectName(u"INMaterialName")
        self.INMaterialName.setMinimumSize(QSize(93, 0))

        self.gridLayout.addWidget(self.INMaterialName, 0, 1, 1, 1)

        self.BRegenElasticConstants = QPushButton(self.page_5)
        self.BRegenElasticConstants.setObjectName(u"BRegenElasticConstants")
        self.BRegenElasticConstants.setGeometry(QRect(0, 556, 353, 22))
        self.DefineMaterialsOptions.addWidget(self.page_5)

        self.verticalLayout_12.addWidget(self.DefineMaterialsOptions)

        self.tabWidget.addTab(self.DefineMaterials, icon7, "")
        self.AssignMaterials = QWidget()
        self.AssignMaterials.setObjectName(u"AssignMaterials")
        self.Tassignmat = QTableWidget(self.AssignMaterials)
        if (self.Tassignmat.columnCount() < 7):
            self.Tassignmat.setColumnCount(7)
        __qtablewidgetitem81 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(0, __qtablewidgetitem81)
        __qtablewidgetitem82 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(1, __qtablewidgetitem82)
        __qtablewidgetitem83 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(2, __qtablewidgetitem83)
        __qtablewidgetitem84 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(3, __qtablewidgetitem84)
        __qtablewidgetitem85 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(4, __qtablewidgetitem85)
        __qtablewidgetitem86 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(5, __qtablewidgetitem86)
        __qtablewidgetitem87 = QTableWidgetItem()
        self.Tassignmat.setHorizontalHeaderItem(6, __qtablewidgetitem87)
        self.Tassignmat.setObjectName(u"Tassignmat")
        self.Tassignmat.setGeometry(QRect(4, 206, 360, 192))
        self.frame_17 = QFrame(self.AssignMaterials)
        self.frame_17.setObjectName(u"frame_17")
        self.frame_17.setGeometry(QRect(4, 394, 360, 69))
        self.frame_17.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_17.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_18 = QHBoxLayout(self.frame_17)
        self.horizontalLayout_18.setObjectName(u"horizontalLayout_18")
        self.horizontalLayout_18.setContentsMargins(0, 0, 0, 0)
        self.BUpMaterial = QToolButton(self.frame_17)
        self.BUpMaterial.setObjectName(u"BUpMaterial")
        self.BUpMaterial.setIconSize(QSize(32, 32))
        self.BUpMaterial.setAutoRaise(False)
        self.BUpMaterial.setArrowType(Qt.ArrowType.UpArrow)

        self.horizontalLayout_18.addWidget(self.BUpMaterial)

        self.BDownMaterial = QToolButton(self.frame_17)
        self.BDownMaterial.setObjectName(u"BDownMaterial")
        self.BDownMaterial.setIconSize(QSize(32, 32))
        self.BDownMaterial.setArrowType(Qt.ArrowType.DownArrow)

        self.horizontalLayout_18.addWidget(self.BDownMaterial)

        self.Baddmaterialassignment = QPushButton(self.AssignMaterials)
        self.Baddmaterialassignment.setObjectName(u"Baddmaterialassignment")
        self.Baddmaterialassignment.setGeometry(QRect(4, 178, 360, 22))
        self.Bdeletematerialassignment = QPushButton(self.AssignMaterials)
        self.Bdeletematerialassignment.setObjectName(u"Bdeletematerialassignment")
        self.Bdeletematerialassignment.setGeometry(QRect(4, 459, 360, 22))
        self.frame_9 = QFrame(self.AssignMaterials)
        self.frame_9.setObjectName(u"frame_9")
        self.frame_9.setGeometry(QRect(4, 8, 361, 289))
        self.frame_9.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.frame_9.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_9.setFrameShadow(QFrame.Shadow.Raised)
        self.layoutWidget_2 = QWidget(self.frame_9)
        self.layoutWidget_2.setObjectName(u"layoutWidget_2")
        self.layoutWidget_2.setGeometry(QRect(235, 59, 125, 103))
        self.formLayout_5 = QFormLayout(self.layoutWidget_2)
        self.formLayout_5.setObjectName(u"formLayout_5")
        self.formLayout_5.setVerticalSpacing(7)
        self.formLayout_5.setContentsMargins(0, 0, 0, 0)
        self.LVelx = QLabel(self.layoutWidget_2)
        self.LVelx.setObjectName(u"LVelx")
        self.LVelx.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.formLayout_5.setWidget(1, QFormLayout.LabelRole, self.LVelx)

        self.INVelocityX = QLineEdit(self.layoutWidget_2)
        self.INVelocityX.setObjectName(u"INVelocityX")
        self.INVelocityX.setEnabled(True)

        self.formLayout_5.setWidget(1, QFormLayout.FieldRole, self.INVelocityX)

        self.LVely = QLabel(self.layoutWidget_2)
        self.LVely.setObjectName(u"LVely")
        self.LVely.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.formLayout_5.setWidget(2, QFormLayout.LabelRole, self.LVely)

        self.INVelocityY = QLineEdit(self.layoutWidget_2)
        self.INVelocityY.setObjectName(u"INVelocityY")
        self.INVelocityY.setEnabled(True)

        self.formLayout_5.setWidget(2, QFormLayout.FieldRole, self.INVelocityY)

        self.LVelz = QLabel(self.layoutWidget_2)
        self.LVelz.setObjectName(u"LVelz")
        self.LVelz.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.formLayout_5.setWidget(3, QFormLayout.LabelRole, self.LVelz)

        self.INVelocityZ = QLineEdit(self.layoutWidget_2)
        self.INVelocityZ.setObjectName(u"INVelocityZ")
        self.INVelocityZ.setEnabled(True)

        self.formLayout_5.setWidget(3, QFormLayout.FieldRole, self.INVelocityZ)

        self.label_2 = QLabel(self.layoutWidget_2)
        self.label_2.setObjectName(u"label_2")

        self.formLayout_5.setWidget(0, QFormLayout.SpanningRole, self.label_2)

        self.layoutWidget_3 = QWidget(self.frame_9)
        self.layoutWidget_3.setObjectName(u"layoutWidget_3")
        self.layoutWidget_3.setGeometry(QRect(0, 0, 361, 53))
        self.formLayout_10 = QFormLayout(self.layoutWidget_3)
        self.formLayout_10.setObjectName(u"formLayout_10")
        self.formLayout_10.setContentsMargins(0, 0, 0, 0)
        self.LMaterialName_3 = QLabel(self.layoutWidget_3)
        self.LMaterialName_3.setObjectName(u"LMaterialName_3")
        self.LMaterialName_3.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.formLayout_10.setWidget(0, QFormLayout.LabelRole, self.LMaterialName_3)

        self.INPartMaterial = QComboBox(self.layoutWidget_3)
        self.INPartMaterial.setObjectName(u"INPartMaterial")

        self.formLayout_10.setWidget(0, QFormLayout.FieldRole, self.INPartMaterial)

        self.LRegion_3 = QLabel(self.layoutWidget_3)
        self.LRegion_3.setObjectName(u"LRegion_3")
        self.LRegion_3.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.formLayout_10.setWidget(1, QFormLayout.LabelRole, self.LRegion_3)

        self.INMaterial = QComboBox(self.layoutWidget_3)
        self.INMaterial.setObjectName(u"INMaterial")

        self.formLayout_10.setWidget(1, QFormLayout.FieldRole, self.INMaterial)

        self.layoutWidget_4 = QWidget(self.frame_9)
        self.layoutWidget_4.setObjectName(u"layoutWidget_4")
        self.layoutWidget_4.setGeometry(QRect(0, 59, 205, 101))
        self.verticalLayout_8 = QVBoxLayout(self.layoutWidget_4)
        self.verticalLayout_8.setObjectName(u"verticalLayout_8")
        self.verticalLayout_8.setContentsMargins(0, 0, 0, 0)
        self.LDensity = QLabel(self.layoutWidget_4)
        self.LDensity.setObjectName(u"LDensity")
        self.LDensity.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.verticalLayout_8.addWidget(self.LDensity)

        self.INDensity = QLineEdit(self.layoutWidget_4)
        self.INDensity.setObjectName(u"INDensity")
        self.INDensity.setEnabled(True)

        self.verticalLayout_8.addWidget(self.INDensity)

        self.LSIE = QLabel(self.layoutWidget_4)
        self.LSIE.setObjectName(u"LSIE")
        self.LSIE.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)
        self.LSIE.setWordWrap(False)

        self.verticalLayout_8.addWidget(self.LSIE)

        self.INSIE = QLineEdit(self.layoutWidget_4)
        self.INSIE.setObjectName(u"INSIE")
        self.INSIE.setEnabled(True)

        self.verticalLayout_8.addWidget(self.INSIE)

        icon11 = QIcon()
        icon11.addFile(u":/Blue Icons/Blue Icons/Clipboard.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.tabWidget.addTab(self.AssignMaterials, icon11, "")

        self.verticalLayout_13.addWidget(self.tabWidget)

        self.ToolWindow.addWidget(self.MaterialsTool)
        self.page_7 = QWidget()
        self.page_7.setObjectName(u"page_7")
        self.RunOptions = QStackedWidget(self.page_7)
        self.RunOptions.setObjectName(u"RunOptions")
        self.RunOptions.setGeometry(QRect(8, 116, 385, 97))
        self.page_10 = QWidget()
        self.page_10.setObjectName(u"page_10")
        self.BRunSGH = QPushButton(self.page_10)
        self.BRunSGH.setObjectName(u"BRunSGH")
        self.BRunSGH.setGeometry(QRect(0, 8, 384, 30))
        sizePolicy10 = QSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)
        sizePolicy10.setHorizontalStretch(0)
        sizePolicy10.setVerticalStretch(0)
        sizePolicy10.setHeightForWidth(self.BRunSGH.sizePolicy().hasHeightForWidth())
        self.BRunSGH.setSizePolicy(sizePolicy10)
        font13 = QFont()
        font13.setPointSize(12)
        font13.setBold(True)
        self.BRunSGH.setFont(font13)
        self.BRunSGH.setIcon(icon8)
        self.BRunSGH.setIconSize(QSize(24, 24))
        self.RunOptions.addWidget(self.page_10)
        self.page_12 = QWidget()
        self.page_12.setObjectName(u"page_12")
        self.BRunEVPFFT = QPushButton(self.page_12)
        self.BRunEVPFFT.setObjectName(u"BRunEVPFFT")
        self.BRunEVPFFT.setGeometry(QRect(0, 8, 384, 30))
        sizePolicy10.setHeightForWidth(self.BRunEVPFFT.sizePolicy().hasHeightForWidth())
        self.BRunEVPFFT.setSizePolicy(sizePolicy10)
        self.BRunEVPFFT.setFont(font13)
        self.BRunEVPFFT.setIcon(icon8)
        self.BRunEVPFFT.setIconSize(QSize(24, 24))
        self.RunOptions.addWidget(self.page_12)
        self.layoutWidget6 = QWidget(self.page_7)
        self.layoutWidget6.setObjectName(u"layoutWidget6")
        self.layoutWidget6.setGeometry(QRect(8, 1, 385, 113))
        self.verticalLayout_9 = QVBoxLayout(self.layoutWidget6)
        self.verticalLayout_9.setObjectName(u"verticalLayout_9")
        self.verticalLayout_9.setContentsMargins(0, 0, 0, 0)
        self.LRun = QLabel(self.layoutWidget6)
        self.LRun.setObjectName(u"LRun")

        self.verticalLayout_9.addWidget(self.LRun)

        self.formLayout_4 = QFormLayout()
        self.formLayout_4.setObjectName(u"formLayout_4")
        self.label_3 = QLabel(self.layoutWidget6)
        self.label_3.setObjectName(u"label_3")

        self.formLayout_4.setWidget(0, QFormLayout.LabelRole, self.label_3)

        self.INRunSelection = QComboBox(self.layoutWidget6)
        self.INRunSelection.addItem("")
        self.INRunSelection.addItem("")
        self.INRunSelection.setObjectName(u"INRunSelection")

        self.formLayout_4.setWidget(0, QFormLayout.FieldRole, self.INRunSelection)


        self.verticalLayout_9.addLayout(self.formLayout_4)

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
        sizePolicy2.setHeightForWidth(self.tabWidget_2.sizePolicy().hasHeightForWidth())
        self.tabWidget_2.setSizePolicy(sizePolicy2)
        self.tabWidget_2.setIconSize(QSize(20, 20))
        self.tabWidget_2.setElideMode(Qt.TextElideMode.ElideNone)
        self.Results = QWidget()
        self.Results.setObjectName(u"Results")
        self.INSelectPostprocessing = QComboBox(self.Results)
        self.INSelectPostprocessing.addItem("")
        self.INSelectPostprocessing.addItem("")
        self.INSelectPostprocessing.setObjectName(u"INSelectPostprocessing")
        self.INSelectPostprocessing.setGeometry(QRect(6, 8, 364, 22))
        self.PostprocessingOptions = QStackedWidget(self.Results)
        self.PostprocessingOptions.setObjectName(u"PostprocessingOptions")
        self.PostprocessingOptions.setGeometry(QRect(6, 36, 364, 605))
        self.page_11 = QWidget()
        self.page_11.setObjectName(u"page_11")
        self.frame_14 = QFrame(self.page_11)
        self.frame_14.setObjectName(u"frame_14")
        self.frame_14.setGeometry(QRect(56, 117, 247, 37))
        self.frame_14.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_14.setFrameShadow(QFrame.Shadow.Raised)
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
        icon12 = QIcon()
        icon12.addFile(u":/Blue Icons/Blue Icons/crop.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.BThreshold.setIcon(icon12)
        self.BThreshold.setIconSize(QSize(32, 32))
        self.BThreshold.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)

        self.horizontalLayout_16.addWidget(self.BThreshold)

        self.LResultsSGH = QLabel(self.page_11)
        self.LResultsSGH.setObjectName(u"LResultsSGH")
        self.LResultsSGH.setGeometry(QRect(4, 24, 352, 12))
        self.frame_13 = QFrame(self.page_11)
        self.frame_13.setObjectName(u"frame_13")
        self.frame_13.setGeometry(QRect(89, 72, 180, 37))
        self.frame_13.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_13.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_12 = QHBoxLayout(self.frame_13)
        self.horizontalLayout_12.setObjectName(u"horizontalLayout_12")
        self.horizontalLayout_12.setContentsMargins(0, 0, 0, 0)
        self.BFirstFrame = QToolButton(self.frame_13)
        self.BFirstFrame.setObjectName(u"BFirstFrame")
        icon13 = QIcon()
        icon13.addFile(u":/Blue Icons/Blue Icons/FirstFrame.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.BFirstFrame.setIcon(icon13)
        self.BFirstFrame.setIconSize(QSize(32, 32))
        self.BFirstFrame.setPopupMode(QToolButton.ToolButtonPopupMode.DelayedPopup)
        self.BFirstFrame.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        self.BFirstFrame.setAutoRaise(False)
        self.BFirstFrame.setArrowType(Qt.ArrowType.NoArrow)

        self.horizontalLayout_12.addWidget(self.BFirstFrame)

        self.BPreviousFrame = QToolButton(self.frame_13)
        self.BPreviousFrame.setObjectName(u"BPreviousFrame")
        icon14 = QIcon()
        icon14.addFile(u":/Blue Icons/Blue Icons/PreviousFrame.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.BPreviousFrame.setIcon(icon14)
        self.BPreviousFrame.setIconSize(QSize(32, 32))

        self.horizontalLayout_12.addWidget(self.BPreviousFrame)

        self.BNextFrame = QToolButton(self.frame_13)
        self.BNextFrame.setObjectName(u"BNextFrame")
        icon15 = QIcon()
        icon15.addFile(u":/Blue Icons/Blue Icons/NextFrame.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.BNextFrame.setIcon(icon15)
        self.BNextFrame.setIconSize(QSize(32, 32))

        self.horizontalLayout_12.addWidget(self.BNextFrame)

        self.BLastFrame = QToolButton(self.frame_13)
        self.BLastFrame.setObjectName(u"BLastFrame")
        icon16 = QIcon()
        icon16.addFile(u":/Blue Icons/Blue Icons/LastFrame.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.BLastFrame.setIcon(icon16)
        self.BLastFrame.setIconSize(QSize(32, 32))

        self.horizontalLayout_12.addWidget(self.BLastFrame)

        self.BOpenParaviewSGH = QPushButton(self.page_11)
        self.BOpenParaviewSGH.setObjectName(u"BOpenParaviewSGH")
        self.BOpenParaviewSGH.setGeometry(QRect(4, 162, 352, 20))
        self.frame_15 = QFrame(self.page_11)
        self.frame_15.setObjectName(u"frame_15")
        self.frame_15.setGeometry(QRect(89, 44, 178, 20))
        self.frame_15.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_15.setFrameShadow(QFrame.Shadow.Raised)
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

        self.PostprocessingOptions.addWidget(self.page_11)
        self.page_9 = QWidget()
        self.page_9.setObjectName(u"page_9")
        self.PreviewResults = QFrame(self.page_9)
        self.PreviewResults.setObjectName(u"PreviewResults")
        self.PreviewResults.setGeometry(QRect(0, 20, 364, 22))
        self.PreviewResults.setFrameShape(QFrame.Shape.NoFrame)
        self.PreviewResults.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_19 = QHBoxLayout(self.PreviewResults)
        self.horizontalLayout_19.setObjectName(u"horizontalLayout_19")
        self.horizontalLayout_19.setContentsMargins(0, 0, 0, 0)
        self.INBCFile = QComboBox(self.PreviewResults)
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.addItem("")
        self.INBCFile.setObjectName(u"INBCFile")

        self.horizontalLayout_19.addWidget(self.INBCFile)

        self.INPreviewResults = QComboBox(self.PreviewResults)
        self.INPreviewResults.addItem("")
        self.INPreviewResults.addItem("")
        self.INPreviewResults.setObjectName(u"INPreviewResults")
        self.INPreviewResults.setFrame(True)

        self.horizontalLayout_19.addWidget(self.INPreviewResults)

        self.INResultRegion = QComboBox(self.PreviewResults)
        self.INResultRegion.addItem("")
        self.INResultRegion.addItem("")
        self.INResultRegion.addItem("")
        self.INResultRegion.setObjectName(u"INResultRegion")

        self.horizontalLayout_19.addWidget(self.INResultRegion)

        self.BOpenParaview = QPushButton(self.page_9)
        self.BOpenParaview.setObjectName(u"BOpenParaview")
        self.BOpenParaview.setGeometry(QRect(0, 76, 364, 22))
        self.THomogenization = QTableWidget(self.page_9)
        if (self.THomogenization.columnCount() < 1):
            self.THomogenization.setColumnCount(1)
        __qtablewidgetitem88 = QTableWidgetItem()
        self.THomogenization.setHorizontalHeaderItem(0, __qtablewidgetitem88)
        if (self.THomogenization.rowCount() < 12):
            self.THomogenization.setRowCount(12)
        __qtablewidgetitem89 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(0, __qtablewidgetitem89)
        __qtablewidgetitem90 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(1, __qtablewidgetitem90)
        __qtablewidgetitem91 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(2, __qtablewidgetitem91)
        __qtablewidgetitem92 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(3, __qtablewidgetitem92)
        __qtablewidgetitem93 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(4, __qtablewidgetitem93)
        __qtablewidgetitem94 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(5, __qtablewidgetitem94)
        __qtablewidgetitem95 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(6, __qtablewidgetitem95)
        __qtablewidgetitem96 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(7, __qtablewidgetitem96)
        __qtablewidgetitem97 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(8, __qtablewidgetitem97)
        __qtablewidgetitem98 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(9, __qtablewidgetitem98)
        __qtablewidgetitem99 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(10, __qtablewidgetitem99)
        __qtablewidgetitem100 = QTableWidgetItem()
        self.THomogenization.setVerticalHeaderItem(11, __qtablewidgetitem100)
        self.THomogenization.setObjectName(u"THomogenization")
        self.THomogenization.setEnabled(False)
        self.THomogenization.setGeometry(QRect(0, 104, 364, 369))
        self.THomogenization.setWordWrap(False)
        self.THomogenization.horizontalHeader().setMinimumSectionSize(40)
        self.THomogenization.horizontalHeader().setDefaultSectionSize(250)
        self.THomogenization.horizontalHeader().setStretchLastSection(True)
        self.THomogenization.verticalHeader().setDefaultSectionSize(30)
        self.BPreviewResults = QPushButton(self.page_9)
        self.BPreviewResults.setObjectName(u"BPreviewResults")
        self.BPreviewResults.setGeometry(QRect(0, 48, 364, 22))
        self.BPreviewResults.setMinimumSize(QSize(140, 0))
        self.BPreviewResults.setMaximumSize(QSize(16777215, 16777215))
        self.PostprocessingOptions.addWidget(self.page_9)
        icon17 = QIcon()
        icon17.addFile(u":/Blue Icons/Blue Icons/magnify.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        icon17.addFile(u":/Blue Icons/Blue Icons/magnify.svg", QSize(), QIcon.Mode.Selected, QIcon.State.On)
        self.tabWidget_2.addTab(self.Results, icon17, "")

        self.verticalLayout_18.addWidget(self.tabWidget_2)

        self.ToolWindow.addWidget(self.ResultsEVPFFTTool)
        self.ResultsSGHTool = QWidget()
        self.ResultsSGHTool.setObjectName(u"ResultsSGHTool")
        self.verticalLayout_38 = QVBoxLayout(self.ResultsSGHTool)
        self.verticalLayout_38.setObjectName(u"verticalLayout_38")
        self.verticalSpacer_17 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_38.addItem(self.verticalSpacer_17)

        self.ToolWindow.addWidget(self.ResultsSGHTool)

        self.Windows.setWidget(1, QFormLayout.LabelRole, self.ToolWindow)


        self.verticalLayout_3.addLayout(self.Windows)

        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1200, 19))
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setObjectName(u"menuHelp")
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        MainWindow.setMenuBar(self.menubar)
        QWidget.setTabOrder(self.INSelectGeometryImport, self.INPartName)
        QWidget.setTabOrder(self.INPartName, self.BUploadGeometryFile)
        QWidget.setTabOrder(self.BUploadGeometryFile, self.INNumberOfVoxelsX)
        QWidget.setTabOrder(self.INNumberOfVoxelsX, self.INNumberOfVoxelsY)
        QWidget.setTabOrder(self.INNumberOfVoxelsY, self.INNumberOfVoxelsZ)
        QWidget.setTabOrder(self.INNumberOfVoxelsZ, self.INOriginX)
        QWidget.setTabOrder(self.INOriginX, self.INOriginY)
        QWidget.setTabOrder(self.INOriginY, self.INOriginZ)
        QWidget.setTabOrder(self.INOriginZ, self.BCustomDimensions)
        QWidget.setTabOrder(self.BCustomDimensions, self.BStlDimensions)
        QWidget.setTabOrder(self.BStlDimensions, self.INLengthX)
        QWidget.setTabOrder(self.INLengthX, self.INLengthY)
        QWidget.setTabOrder(self.INLengthY, self.INLengthZ)
        QWidget.setTabOrder(self.INLengthZ, self.BVoxelizeGeometry)
        QWidget.setTabOrder(self.BVoxelizeGeometry, self.INSelectBasicGeometry)
        QWidget.setTabOrder(self.INSelectBasicGeometry, self.INBasicGeometryName)
        QWidget.setTabOrder(self.INBasicGeometryName, self.INBoxx1)
        QWidget.setTabOrder(self.INBoxx1, self.INBoxx2)
        QWidget.setTabOrder(self.INBoxx2, self.INBoxy1)
        QWidget.setTabOrder(self.INBoxy1, self.INBoxy2)
        QWidget.setTabOrder(self.INBoxy2, self.INBoxz1)
        QWidget.setTabOrder(self.INBoxz1, self.INBoxz2)
        QWidget.setTabOrder(self.INBoxz2, self.BGenerateBasicGeometry)
        QWidget.setTabOrder(self.BGenerateBasicGeometry, self.RunOutputWindow)
        QWidget.setTabOrder(self.RunOutputWindow, self.NavigationMenu)
        QWidget.setTabOrder(self.NavigationMenu, self.INGraphicsOutput)
        QWidget.setTabOrder(self.INGraphicsOutput, self.INInitialdt)
        QWidget.setTabOrder(self.INInitialdt, self.INmaxcycles)
        QWidget.setTabOrder(self.INmaxcycles, self.INBCFile)
        QWidget.setTabOrder(self.INBCFile, self.INBoundaryCondition)
        QWidget.setTabOrder(self.INBoundaryCondition, self.INPreviewResults)
        QWidget.setTabOrder(self.INPreviewResults, self.INBCDirection)
        QWidget.setTabOrder(self.INBCDirection, self.INResultRegion)
        QWidget.setTabOrder(self.INResultRegion, self.BImageToVTK)
        QWidget.setTabOrder(self.BImageToVTK, self.BTiffToStl)
        QWidget.setTabOrder(self.BTiffToStl, self.BDeleteBasicGeometry)
        QWidget.setTabOrder(self.BDeleteBasicGeometry, self.INMaxdt)
        QWidget.setTabOrder(self.INMaxdt, self.ParaviewFrame)
        QWidget.setTabOrder(self.ParaviewFrame, self.graphicsView)
        QWidget.setTabOrder(self.graphicsView, self.INTime)
        QWidget.setTabOrder(self.INTime, self.GeometryOptions)
        QWidget.setTabOrder(self.GeometryOptions, self.INSphereri)
        QWidget.setTabOrder(self.INSphereri, self.INSpherero)
        QWidget.setTabOrder(self.INSpherero, self.INSphereox)
        QWidget.setTabOrder(self.INSphereox, self.INSphereoy)
        QWidget.setTabOrder(self.INSphereoy, self.INSphereoz)
        QWidget.setTabOrder(self.INSphereoz, self.INCylinderri)
        QWidget.setTabOrder(self.INCylinderri, self.INCylinderro)
        QWidget.setTabOrder(self.INCylinderro, self.INMindt)
        QWidget.setTabOrder(self.INMindt, self.TBasicGeometries)
        QWidget.setTabOrder(self.TBasicGeometries, self.SAGeometryScrollArea)
        QWidget.setTabOrder(self.SAGeometryScrollArea, self.INElementType)
        QWidget.setTabOrder(self.INElementType, self.INCoordinateSystem)
        QWidget.setTabOrder(self.INCoordinateSystem, self.INDimension)
        QWidget.setTabOrder(self.INDimension, self.INLengthYR3D)
        QWidget.setTabOrder(self.INLengthYR3D, self.INOriginYR3D)
        QWidget.setTabOrder(self.INOriginYR3D, self.INElementsYR3D)
        QWidget.setTabOrder(self.INElementsYR3D, self.INLengthXR3D)
        QWidget.setTabOrder(self.INLengthXR3D, self.INOriginXR3D)
        QWidget.setTabOrder(self.INOriginXR3D, self.INElementsXR3D)
        QWidget.setTabOrder(self.INElementsXR3D, self.INOriginZR3D)
        QWidget.setTabOrder(self.INOriginZR3D, self.INLengthZR3D)
        QWidget.setTabOrder(self.INLengthZR3D, self.INElementsZR3D)
        QWidget.setTabOrder(self.INElementsZR3D, self.INLengthXR2D)
        QWidget.setTabOrder(self.INLengthXR2D, self.INElementsYR2D)
        QWidget.setTabOrder(self.INElementsYR2D, self.INLengthYR2D)
        QWidget.setTabOrder(self.INLengthYR2D, self.INElementsXR2D)
        QWidget.setTabOrder(self.INElementsXR2D, self.INOriginXR2D)
        QWidget.setTabOrder(self.INOriginXR2D, self.INOriginYR2D)
        QWidget.setTabOrder(self.INOriginYR2D, self.INInnerRadiusC2D)
        QWidget.setTabOrder(self.INInnerRadiusC2D, self.INLengthThetaC2D)
        QWidget.setTabOrder(self.INLengthThetaC2D, self.INLengthOutRadC2D)
        QWidget.setTabOrder(self.INLengthOutRadC2D, self.INOriginYC2D)
        QWidget.setTabOrder(self.INOriginYC2D, self.INElementsArcC2D)
        QWidget.setTabOrder(self.INElementsArcC2D, self.INOriginXC2D)
        QWidget.setTabOrder(self.INOriginXC2D, self.INElementsRadialC2D)
        QWidget.setTabOrder(self.INElementsRadialC2D, self.INInnerRadiusC3D)
        QWidget.setTabOrder(self.INInnerRadiusC3D, self.INElementsRadC3D)
        QWidget.setTabOrder(self.INElementsRadC3D, self.INElementsArcC3D)
        QWidget.setTabOrder(self.INElementsArcC3D, self.INLengthThetaC3D)
        QWidget.setTabOrder(self.INLengthThetaC3D, self.INOriginYC3D)
        QWidget.setTabOrder(self.INOriginYC3D, self.INLengthOutRadC3D)
        QWidget.setTabOrder(self.INLengthOutRadC3D, self.INOriginXC3D)
        QWidget.setTabOrder(self.INOriginXC3D, self.INOriginZC3D)
        QWidget.setTabOrder(self.INOriginZC3D, self.INLengthZC3D)
        QWidget.setTabOrder(self.INLengthZC3D, self.INElementsZC3D)
        QWidget.setTabOrder(self.INElementsZC3D, self.BGenerateGlobalMesh)
        QWidget.setTabOrder(self.BGenerateGlobalMesh, self.BAddBC)
        QWidget.setTabOrder(self.BAddBC, self.TBCs)
        QWidget.setTabOrder(self.TBCs, self.BDeleteBC)
        QWidget.setTabOrder(self.BDeleteBC, self.INSelectSolverSettings)
        QWidget.setTabOrder(self.INSelectSolverSettings, self.INNumberOfSteps)
        QWidget.setTabOrder(self.INNumberOfSteps, self.BAddMaterialSGH)
        QWidget.setTabOrder(self.BAddMaterialSGH, self.INq1)
        QWidget.setTabOrder(self.INq1, self.INq2)
        QWidget.setTabOrder(self.INq2, self.INq1ex)
        QWidget.setTabOrder(self.INq1ex, self.INq2ex)
        QWidget.setTabOrder(self.INq2ex, self.INGamma)
        QWidget.setTabOrder(self.INGamma, self.INMinSound)
        QWidget.setTabOrder(self.INMinSound, self.INSpecificHeat)
        QWidget.setTabOrder(self.INSpecificHeat, self.BDeleteMaterialSGH)
        QWidget.setTabOrder(self.BDeleteMaterialSGH, self.TMaterialsSGH)
        QWidget.setTabOrder(self.TMaterialsSGH, self.INMaterialNameSGH)
        QWidget.setTabOrder(self.INMaterialNameSGH, self.INArtificialViscosity)
        QWidget.setTabOrder(self.INArtificialViscosity, self.INEOS)
        QWidget.setTabOrder(self.INEOS, self.BAddMaterial)
        QWidget.setTabOrder(self.BAddMaterial, self.INYoungsModulus)
        QWidget.setTabOrder(self.INYoungsModulus, self.INPoissonsRatio)
        QWidget.setTabOrder(self.INPoissonsRatio, self.INIsotropicPlane)
        QWidget.setTabOrder(self.INIsotropicPlane, self.INEip)
        QWidget.setTabOrder(self.INEip, self.INNUip)
        QWidget.setTabOrder(self.INNUip, self.INEop)
        QWidget.setTabOrder(self.INEop, self.INNUop)
        QWidget.setTabOrder(self.INNUop, self.INGop)
        QWidget.setTabOrder(self.INGop, self.TAnisotropic)
        QWidget.setTabOrder(self.TAnisotropic, self.INNUyz)
        QWidget.setTabOrder(self.INNUyz, self.INGxy)
        QWidget.setTabOrder(self.INGxy, self.INEy)
        QWidget.setTabOrder(self.INEy, self.INNUxy)
        QWidget.setTabOrder(self.INNUxy, self.INEz)
        QWidget.setTabOrder(self.INEz, self.INGyz)
        QWidget.setTabOrder(self.INGyz, self.INEx)
        QWidget.setTabOrder(self.INEx, self.INGxz)
        QWidget.setTabOrder(self.INGxz, self.INNUxz)
        QWidget.setTabOrder(self.INNUxz, self.TMaterials)
        QWidget.setTabOrder(self.TMaterials, self.BDeleteMaterial)
        QWidget.setTabOrder(self.BDeleteMaterial, self.INSolidGas)
        QWidget.setTabOrder(self.INSolidGas, self.INMaterialType)
        QWidget.setTabOrder(self.INMaterialType, self.INRegion)
        QWidget.setTabOrder(self.INRegion, self.INMaterialName)
        QWidget.setTabOrder(self.INMaterialName, self.BRegenElasticConstants)
        QWidget.setTabOrder(self.BRegenElasticConstants, self.INSelectBoundaryConditions)
        QWidget.setTabOrder(self.INSelectBoundaryConditions, self.INBoundary)
        QWidget.setTabOrder(self.INBoundary, self.INPlanePosition)
        QWidget.setTabOrder(self.INPlanePosition, self.INType)
        QWidget.setTabOrder(self.INType, self.INVel0)
        QWidget.setTabOrder(self.INVel0, self.INVel1)
        QWidget.setTabOrder(self.INVel1, self.INVelstart)
        QWidget.setTabOrder(self.INVelstart, self.INVelend)
        QWidget.setTabOrder(self.INVelend, self.BaddBC)
        QWidget.setTabOrder(self.BaddBC, self.BdeleteBC)
        QWidget.setTabOrder(self.BdeleteBC, self.TBoundaryConditions)
        QWidget.setTabOrder(self.TBoundaryConditions, self.INSelectPostprocessing)
        QWidget.setTabOrder(self.INSelectPostprocessing, self.INThreshold)
        QWidget.setTabOrder(self.INThreshold, self.BThreshold)
        QWidget.setTabOrder(self.BThreshold, self.BFirstFrame)
        QWidget.setTabOrder(self.BFirstFrame, self.BPreviousFrame)
        QWidget.setTabOrder(self.BPreviousFrame, self.BNextFrame)
        QWidget.setTabOrder(self.BNextFrame, self.BLastFrame)
        QWidget.setTabOrder(self.BLastFrame, self.BOpenParaviewSGH)
        QWidget.setTabOrder(self.BOpenParaviewSGH, self.INOuputVarSGH)
        QWidget.setTabOrder(self.INOuputVarSGH, self.BPreviewResultsSGH)
        QWidget.setTabOrder(self.BPreviewResultsSGH, self.BOpenParaview)
        QWidget.setTabOrder(self.BOpenParaview, self.THomogenization)
        QWidget.setTabOrder(self.THomogenization, self.BPreviewResults)
        QWidget.setTabOrder(self.BPreviewResults, self.tabWidget_2)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.menuHelp.addAction(self.actionManual)
        self.menuFile.addAction(self.actionChange_Working_Directory)

        self.retranslateUi(MainWindow)

        self.NavigationMenu.setCurrentIndex(2)
        self.OutputWindows.setCurrentIndex(0)
        self.ToolWindow.setCurrentIndex(2)
        self.INPipelineSelection.setCurrentIndex(0)
        self.INSelectGeometryImport.setCurrentIndex(0)
        self.GeometryOptions.setCurrentIndex(2)
        self.BasicGeometries.setCurrentIndex(0)
        self.MeshInputs2.setCurrentIndex(0)
        self.SolverSettingsOptions.setCurrentIndex(0)
        self.BoundaryConditionsOptions.setCurrentIndex(0)
        self.tabWidget.setCurrentIndex(0)
        self.DefineMaterialsOptions.setCurrentIndex(1)
        self.MaterialTypeTool.setCurrentIndex(3)
        self.RunOptions.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(0)
        self.PostprocessingOptions.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"Fierro", None))
        self.actionManual.setText(QCoreApplication.translate("MainWindow", u"Manual", None))
        self.actionChange_Working_Directory.setText(QCoreApplication.translate("MainWindow", u"Change Fierro Setup", None))
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
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Solver), QCoreApplication.translate("MainWindow", u"Solver Settings", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.Solver), QCoreApplication.translate("MainWindow", u"Solver", None))
#endif // QT_CONFIG(tooltip)
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.BoundaryConditions), QCoreApplication.translate("MainWindow", u"Boundary Conditions", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.BoundaryConditions), QCoreApplication.translate("MainWindow", u"Boundary Conditions", None))
#endif // QT_CONFIG(tooltip)
        self.NavigationMenu.setTabText(self.NavigationMenu.indexOf(self.Materials), QCoreApplication.translate("MainWindow", u"Materials", None))
#if QT_CONFIG(tooltip)
        self.NavigationMenu.setTabToolTip(self.NavigationMenu.indexOf(self.Materials), QCoreApplication.translate("MainWindow", u"Materials", None))
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
        self.LPipeline.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/Clipboard.svg\" width=\"40\" height=\"40\"/></p><p align=\"center\"><span style=\" font-weight:700;\">Select Pipeline</span></p></body></html>", None))
        self.INPipelineSelection.setItemText(0, QCoreApplication.translate("MainWindow", u"No Pipeline", None))
        self.INPipelineSelection.setItemText(1, QCoreApplication.translate("MainWindow", u"Explicit SGH", None))
        self.INPipelineSelection.setItemText(2, QCoreApplication.translate("MainWindow", u"EVPFFT", None))

        self.LPipelineSelection.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Selection:</p></body></html>", None))
        self.LGeometryInformation.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/Shapes.svg\" width=\"50\" height=\"40\"/></p><p align=\"center\"><span style=\" font-weight:700;\">Define Geometry</span></p></body></html>", None))
        self.LPartName_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Selection:</p></body></html>", None))
        self.INSelectGeometryImport.setItemText(0, QCoreApplication.translate("MainWindow", u"Import Part", None))
        self.INSelectGeometryImport.setItemText(1, QCoreApplication.translate("MainWindow", u"Import Image Stack", None))
        self.INSelectGeometryImport.setItemText(2, QCoreApplication.translate("MainWindow", u"Import Polycrystal", None))
        self.INSelectGeometryImport.setItemText(3, QCoreApplication.translate("MainWindow", u"Create Basic Part", None))

        self.LPartName.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Part Name:</p></body></html>", None))
        self.INPartName.setText("")
        self.BUploadGeometryFile.setText(QCoreApplication.translate("MainWindow", u"Upload Geometry", None))
        self.BVoxelizeGeometry.setText(QCoreApplication.translate("MainWindow", u"Voxelize Geometry", None))
        self.LDimensions.setText(QCoreApplication.translate("MainWindow", u"Dimensions", None))
        self.LSTLVoxelization.setText(QCoreApplication.translate("MainWindow", u"STL Voxelization", None))
        self.BCustomDimensions.setText(QCoreApplication.translate("MainWindow", u"Custom Dimensions", None))
        self.BStlDimensions.setText(QCoreApplication.translate("MainWindow", u"STL Dimensions", None))
        self.LOriginX.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.INOriginX.setText("")
        self.LOriginY.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.INOriginY.setText("")
        self.LOriginZ.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.INOriginZ.setText("")
        self.LOriginPoint.setText(QCoreApplication.translate("MainWindow", u"Origin Point", None))
        self.LNumberOfVoxelsX.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.INNumberOfVoxelsX.setText("")
        self.LNumberOfVoxelsY.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.INNumberOfVoxelsY.setText("")
        self.LNumberOfVoxelsZ.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.INNumberOfVoxelsZ.setText("")
        self.LVoxelCount.setText(QCoreApplication.translate("MainWindow", u"Voxel Count", None))
        self.LLengthX.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.LLengthY.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.INLengthY.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.LLengthZ.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.INLengthZ.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.INLengthX.setText(QCoreApplication.translate("MainWindow", u"0", None))
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
        self.LImageFileFormat.setText(QCoreApplication.translate("MainWindow", u"Image File Format:", None))
        self.LImageStack.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:15pt; font-weight:700;\">STL/VTK Conversion</span></p><p align=\"center\"><br/></p></body></html>", None))
        self.BImageToVTK.setText(QCoreApplication.translate("MainWindow", u"Voxelize", None))
        self.BTiffToStl.setText(QCoreApplication.translate("MainWindow", u"Convert to STL", None))
        self.LUploadedDirectory.setText(QCoreApplication.translate("MainWindow", u"Directory:", None))
        self.LImageStack_2.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:15pt; font-weight:700;\">Dream3D</span></p><p align=\"center\"><br/></p></body></html>", None))
        self.BDeleteBasicGeometry.setText(QCoreApplication.translate("MainWindow", u"Delete Geometry", None))
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
        self.INSelectBasicGeometry.setItemText(0, QCoreApplication.translate("MainWindow", u"box", None))
        self.INSelectBasicGeometry.setItemText(1, QCoreApplication.translate("MainWindow", u"sphere", None))

        self.LSelectGeometry.setText(QCoreApplication.translate("MainWindow", u"Select Geometry:", None))
        self.LBasicGName.setText(QCoreApplication.translate("MainWindow", u"Part Name:", None))
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
        self.LBasicGeometry.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:15pt; font-weight:700;\">Create Basic Part</span></p><p align=\"center\"><br/></p></body></html>", None))
        self.BGenerateBasicGeometry.setText(QCoreApplication.translate("MainWindow", u"Generate Geometry", None))
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
        self.LSolverSettings.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/gear.svg\" width=\"40\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Solver Settings</span></p></body></html>", None))
        self.INSelectSolverSettings.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectSolverSettings.setItemText(1, QCoreApplication.translate("MainWindow", u"EVPFFT", None))

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
        self.LNumberOfSteps.setText(QCoreApplication.translate("MainWindow", u"Number of steps: ", None))
        self.INNumberOfSteps.setInputMask("")
        self.INNumberOfSteps.setText(QCoreApplication.translate("MainWindow", u"10", None))
        self.INNumberOfSteps.setPlaceholderText("")
        self.LBoundaryConditions.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/brick.svg\" width=\"40\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Boundary Conditions</span></p></body></html>", None))
        self.INSelectBoundaryConditions.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectBoundaryConditions.setItemText(1, QCoreApplication.translate("MainWindow", u"EVPFFT", None))

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
        self.BdeleteBC.setText(QCoreApplication.translate("MainWindow", u"Delete Boundary Condition", None))
        ___qtablewidgetitem23 = self.TBoundaryConditions.horizontalHeaderItem(0)
        ___qtablewidgetitem23.setText(QCoreApplication.translate("MainWindow", u"Boundary", None));
        ___qtablewidgetitem24 = self.TBoundaryConditions.horizontalHeaderItem(1)
        ___qtablewidgetitem24.setText(QCoreApplication.translate("MainWindow", u"Plane Position", None));
        ___qtablewidgetitem25 = self.TBoundaryConditions.horizontalHeaderItem(2)
        ___qtablewidgetitem25.setText(QCoreApplication.translate("MainWindow", u"Type", None));
        ___qtablewidgetitem26 = self.TBoundaryConditions.horizontalHeaderItem(3)
        ___qtablewidgetitem26.setText(QCoreApplication.translate("MainWindow", u"vel_0", None));
        ___qtablewidgetitem27 = self.TBoundaryConditions.horizontalHeaderItem(4)
        ___qtablewidgetitem27.setText(QCoreApplication.translate("MainWindow", u"vel_1", None));
        ___qtablewidgetitem28 = self.TBoundaryConditions.horizontalHeaderItem(5)
        ___qtablewidgetitem28.setText(QCoreApplication.translate("MainWindow", u"vel_t_start", None));
        ___qtablewidgetitem29 = self.TBoundaryConditions.horizontalHeaderItem(6)
        ___qtablewidgetitem29.setText(QCoreApplication.translate("MainWindow", u"vel_t_end", None));
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
        ___qtablewidgetitem30 = self.TBCs.horizontalHeaderItem(0)
        ___qtablewidgetitem30.setText(QCoreApplication.translate("MainWindow", u"Boundary Condition", None));
        ___qtablewidgetitem31 = self.TBCs.horizontalHeaderItem(1)
        ___qtablewidgetitem31.setText(QCoreApplication.translate("MainWindow", u"Direction", None));
        self.BDeleteBC.setText(QCoreApplication.translate("MainWindow", u"Delete", None))
        self.LDefineMaterials.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/mine.svg\" width=\"45\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Materials</span></p></body></html>", None))
        self.INSelectDefineMaterials.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectDefineMaterials.setItemText(1, QCoreApplication.translate("MainWindow", u"EVPFFT", None))

        self.BAddMaterialSGH.setText(QCoreApplication.translate("MainWindow", u"Add Material", None))
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
        self.BDeleteMaterialSGH.setText(QCoreApplication.translate("MainWindow", u"Delete Material", None))
        ___qtablewidgetitem32 = self.TMaterialsSGH.horizontalHeaderItem(0)
        ___qtablewidgetitem32.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem33 = self.TMaterialsSGH.horizontalHeaderItem(1)
        ___qtablewidgetitem33.setText(QCoreApplication.translate("MainWindow", u"EOS", None));
        ___qtablewidgetitem34 = self.TMaterialsSGH.horizontalHeaderItem(2)
        ___qtablewidgetitem34.setText(QCoreApplication.translate("MainWindow", u"q1", None));
        ___qtablewidgetitem35 = self.TMaterialsSGH.horizontalHeaderItem(3)
        ___qtablewidgetitem35.setText(QCoreApplication.translate("MainWindow", u"q2", None));
        ___qtablewidgetitem36 = self.TMaterialsSGH.horizontalHeaderItem(4)
        ___qtablewidgetitem36.setText(QCoreApplication.translate("MainWindow", u"q1ex", None));
        ___qtablewidgetitem37 = self.TMaterialsSGH.horizontalHeaderItem(5)
        ___qtablewidgetitem37.setText(QCoreApplication.translate("MainWindow", u"q2ex", None));
        ___qtablewidgetitem38 = self.TMaterialsSGH.horizontalHeaderItem(6)
        ___qtablewidgetitem38.setText(QCoreApplication.translate("MainWindow", u"gamma", None));
        ___qtablewidgetitem39 = self.TMaterialsSGH.horizontalHeaderItem(7)
        ___qtablewidgetitem39.setText(QCoreApplication.translate("MainWindow", u"min sound speed", None));
        ___qtablewidgetitem40 = self.TMaterialsSGH.horizontalHeaderItem(8)
        ___qtablewidgetitem40.setText(QCoreApplication.translate("MainWindow", u"specific heat", None));
        self.label_12.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700;\">Custom Artificial Viscosity Variables</span></p></body></html>", None))
        self.LMaterialNameSGH.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Name:</p></body></html>", None))
        self.INMaterialNameSGH.setText("")
        self.INArtificialViscosity.setItemText(0, QCoreApplication.translate("MainWindow", u"default", None))
        self.INArtificialViscosity.setItemText(1, QCoreApplication.translate("MainWindow", u"custom", None))

        self.INEOS.setItemText(0, QCoreApplication.translate("MainWindow", u"ideal_gas", None))

        self.LEOS.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Equation Of State Model:</p></body></html>", None))
        self.LArtificialViscosity.setText(QCoreApplication.translate("MainWindow", u"Artificial Viscosity:", None))
        self.BAddMaterial.setText(QCoreApplication.translate("MainWindow", u"Add", None))
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
        ___qtablewidgetitem41 = self.TMaterials.horizontalHeaderItem(0)
        ___qtablewidgetitem41.setText(QCoreApplication.translate("MainWindow", u"Name", None));
        ___qtablewidgetitem42 = self.TMaterials.horizontalHeaderItem(1)
        ___qtablewidgetitem42.setText(QCoreApplication.translate("MainWindow", u"Region", None));
        ___qtablewidgetitem43 = self.TMaterials.horizontalHeaderItem(2)
        ___qtablewidgetitem43.setText(QCoreApplication.translate("MainWindow", u"Type", None));
        ___qtablewidgetitem44 = self.TMaterials.horizontalHeaderItem(3)
        ___qtablewidgetitem44.setText(QCoreApplication.translate("MainWindow", u"C11", None));
        ___qtablewidgetitem45 = self.TMaterials.horizontalHeaderItem(4)
        ___qtablewidgetitem45.setText(QCoreApplication.translate("MainWindow", u"C12", None));
        ___qtablewidgetitem46 = self.TMaterials.horizontalHeaderItem(5)
        ___qtablewidgetitem46.setText(QCoreApplication.translate("MainWindow", u"C13", None));
        ___qtablewidgetitem47 = self.TMaterials.horizontalHeaderItem(6)
        ___qtablewidgetitem47.setText(QCoreApplication.translate("MainWindow", u"C14", None));
        ___qtablewidgetitem48 = self.TMaterials.horizontalHeaderItem(7)
        ___qtablewidgetitem48.setText(QCoreApplication.translate("MainWindow", u"C15", None));
        ___qtablewidgetitem49 = self.TMaterials.horizontalHeaderItem(8)
        ___qtablewidgetitem49.setText(QCoreApplication.translate("MainWindow", u"C16", None));
        ___qtablewidgetitem50 = self.TMaterials.horizontalHeaderItem(9)
        ___qtablewidgetitem50.setText(QCoreApplication.translate("MainWindow", u"C22", None));
        ___qtablewidgetitem51 = self.TMaterials.horizontalHeaderItem(10)
        ___qtablewidgetitem51.setText(QCoreApplication.translate("MainWindow", u"C23", None));
        ___qtablewidgetitem52 = self.TMaterials.horizontalHeaderItem(11)
        ___qtablewidgetitem52.setText(QCoreApplication.translate("MainWindow", u"C24", None));
        ___qtablewidgetitem53 = self.TMaterials.horizontalHeaderItem(12)
        ___qtablewidgetitem53.setText(QCoreApplication.translate("MainWindow", u"C25", None));
        ___qtablewidgetitem54 = self.TMaterials.horizontalHeaderItem(13)
        ___qtablewidgetitem54.setText(QCoreApplication.translate("MainWindow", u"C26", None));
        ___qtablewidgetitem55 = self.TMaterials.horizontalHeaderItem(14)
        ___qtablewidgetitem55.setText(QCoreApplication.translate("MainWindow", u"C33", None));
        ___qtablewidgetitem56 = self.TMaterials.horizontalHeaderItem(15)
        ___qtablewidgetitem56.setText(QCoreApplication.translate("MainWindow", u"C34", None));
        ___qtablewidgetitem57 = self.TMaterials.horizontalHeaderItem(16)
        ___qtablewidgetitem57.setText(QCoreApplication.translate("MainWindow", u"C35", None));
        ___qtablewidgetitem58 = self.TMaterials.horizontalHeaderItem(17)
        ___qtablewidgetitem58.setText(QCoreApplication.translate("MainWindow", u"C36", None));
        ___qtablewidgetitem59 = self.TMaterials.horizontalHeaderItem(18)
        ___qtablewidgetitem59.setText(QCoreApplication.translate("MainWindow", u"C44", None));
        ___qtablewidgetitem60 = self.TMaterials.horizontalHeaderItem(19)
        ___qtablewidgetitem60.setText(QCoreApplication.translate("MainWindow", u"C45", None));
        ___qtablewidgetitem61 = self.TMaterials.horizontalHeaderItem(20)
        ___qtablewidgetitem61.setText(QCoreApplication.translate("MainWindow", u"C46", None));
        ___qtablewidgetitem62 = self.TMaterials.horizontalHeaderItem(21)
        ___qtablewidgetitem62.setText(QCoreApplication.translate("MainWindow", u"C55", None));
        ___qtablewidgetitem63 = self.TMaterials.horizontalHeaderItem(22)
        ___qtablewidgetitem63.setText(QCoreApplication.translate("MainWindow", u"C56", None));
        ___qtablewidgetitem64 = self.TMaterials.horizontalHeaderItem(23)
        ___qtablewidgetitem64.setText(QCoreApplication.translate("MainWindow", u"C66", None));
        self.BDeleteMaterial.setText(QCoreApplication.translate("MainWindow", u"Delete", None))
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

        self.BRegenElasticConstants.setText(QCoreApplication.translate("MainWindow", u"Regenerate Elastic Constants", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.DefineMaterials), QCoreApplication.translate("MainWindow", u"Define Materials", None))
        ___qtablewidgetitem65 = self.Tassignmat.horizontalHeaderItem(0)
        ___qtablewidgetitem65.setText(QCoreApplication.translate("MainWindow", u"Part", None));
        ___qtablewidgetitem66 = self.Tassignmat.horizontalHeaderItem(1)
        ___qtablewidgetitem66.setText(QCoreApplication.translate("MainWindow", u"Material", None));
        ___qtablewidgetitem67 = self.Tassignmat.horizontalHeaderItem(2)
        ___qtablewidgetitem67.setText(QCoreApplication.translate("MainWindow", u"Density", None));
        ___qtablewidgetitem68 = self.Tassignmat.horizontalHeaderItem(3)
        ___qtablewidgetitem68.setText(QCoreApplication.translate("MainWindow", u"Specific Internal Energy", None));
        ___qtablewidgetitem69 = self.Tassignmat.horizontalHeaderItem(4)
        ___qtablewidgetitem69.setText(QCoreApplication.translate("MainWindow", u"Velocity X", None));
        ___qtablewidgetitem70 = self.Tassignmat.horizontalHeaderItem(5)
        ___qtablewidgetitem70.setText(QCoreApplication.translate("MainWindow", u"Velocity Y", None));
        ___qtablewidgetitem71 = self.Tassignmat.horizontalHeaderItem(6)
        ___qtablewidgetitem71.setText(QCoreApplication.translate("MainWindow", u"Velocity Z", None));
        self.BUpMaterial.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.BDownMaterial.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.Baddmaterialassignment.setText(QCoreApplication.translate("MainWindow", u"Add Material Assignment", None))
        self.Bdeletematerialassignment.setText(QCoreApplication.translate("MainWindow", u"Delete Material Assignment", None))
        self.LVelx.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">X</p></body></html>", None))
        self.INVelocityX.setText("")
        self.LVely.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Y</p></body></html>", None))
        self.INVelocityY.setText("")
        self.LVelz.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Z</p></body></html>", None))
        self.INVelocityZ.setText("")
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"Velocity", None))
        self.LMaterialName_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Region:</p></body></html>", None))
        self.LRegion_3.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"right\">Material:</p></body></html>", None))
        self.LDensity.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Density:</p></body></html>", None))
        self.INDensity.setText("")
        self.LSIE.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"left\">Specific Internal Energy:</p></body></html>", None))
        self.INSIE.setText("")
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.AssignMaterials), QCoreApplication.translate("MainWindow", u"Assign Materials", None))
        self.BRunSGH.setText(QCoreApplication.translate("MainWindow", u"Run Explicit SGH Solver", None))
        self.BRunEVPFFT.setText(QCoreApplication.translate("MainWindow", u"Run EVPFFT Solver", None))
        self.LRun.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/Play.svg\" width=\"45\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Run</span></p></body></html>", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"Select Solver", None))
        self.INRunSelection.setItemText(0, QCoreApplication.translate("MainWindow", u"Explicit SGH", None))
        self.INRunSelection.setItemText(1, QCoreApplication.translate("MainWindow", u"EVPFFT", None))

        self.LPostprocessing.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><img src=\":/Blue Icons/Blue Icons/lens.svg\" width=\"40\" height=\"40\"/></p><p align=\"center\"><span style=\" font-size:16pt; font-weight:700;\">Postprocessing</span></p></body></html>", None))
        self.INSelectPostprocessing.setItemText(0, QCoreApplication.translate("MainWindow", u"SGH", None))
        self.INSelectPostprocessing.setItemText(1, QCoreApplication.translate("MainWindow", u"EVPFFT", None))

        self.LThreshold.setText(QCoreApplication.translate("MainWindow", u"Threshold Value:", None))
        self.INThreshold.setText(QCoreApplication.translate("MainWindow", u"1", None))
        self.BThreshold.setText(QCoreApplication.translate("MainWindow", u"threshold", None))
        self.LResultsSGH.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-weight:700; text-decoration: underline;\">RESULTS</span></p></body></html>", None))
        self.BFirstFrame.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.BPreviousFrame.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.BNextFrame.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.BLastFrame.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.BOpenParaviewSGH.setText(QCoreApplication.translate("MainWindow", u"Open Paraview", None))
        self.INOuputVarSGH.setItemText(0, QCoreApplication.translate("MainWindow", u"SIE", None))

        self.BPreviewResultsSGH.setText(QCoreApplication.translate("MainWindow", u"Preview Results", None))
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

        self.BOpenParaview.setText(QCoreApplication.translate("MainWindow", u"Open Paraview", None))
        ___qtablewidgetitem72 = self.THomogenization.horizontalHeaderItem(0)
        ___qtablewidgetitem72.setText(QCoreApplication.translate("MainWindow", u"Homogenized Elastic Constants", None));
        ___qtablewidgetitem73 = self.THomogenization.verticalHeaderItem(0)
        ___qtablewidgetitem73.setText(QCoreApplication.translate("MainWindow", u"Exx", None));
        ___qtablewidgetitem74 = self.THomogenization.verticalHeaderItem(1)
        ___qtablewidgetitem74.setText(QCoreApplication.translate("MainWindow", u"Eyy", None));
        ___qtablewidgetitem75 = self.THomogenization.verticalHeaderItem(2)
        ___qtablewidgetitem75.setText(QCoreApplication.translate("MainWindow", u"Ezz", None));
        ___qtablewidgetitem76 = self.THomogenization.verticalHeaderItem(3)
        ___qtablewidgetitem76.setText(QCoreApplication.translate("MainWindow", u"NUxy", None));
        ___qtablewidgetitem77 = self.THomogenization.verticalHeaderItem(4)
        ___qtablewidgetitem77.setText(QCoreApplication.translate("MainWindow", u"NUyx", None));
        ___qtablewidgetitem78 = self.THomogenization.verticalHeaderItem(5)
        ___qtablewidgetitem78.setText(QCoreApplication.translate("MainWindow", u"NUxz", None));
        ___qtablewidgetitem79 = self.THomogenization.verticalHeaderItem(6)
        ___qtablewidgetitem79.setText(QCoreApplication.translate("MainWindow", u"NUzx", None));
        ___qtablewidgetitem80 = self.THomogenization.verticalHeaderItem(7)
        ___qtablewidgetitem80.setText(QCoreApplication.translate("MainWindow", u"NUyz", None));
        ___qtablewidgetitem81 = self.THomogenization.verticalHeaderItem(8)
        ___qtablewidgetitem81.setText(QCoreApplication.translate("MainWindow", u"NUzy", None));
        ___qtablewidgetitem82 = self.THomogenization.verticalHeaderItem(9)
        ___qtablewidgetitem82.setText(QCoreApplication.translate("MainWindow", u"Gxy", None));
        ___qtablewidgetitem83 = self.THomogenization.verticalHeaderItem(10)
        ___qtablewidgetitem83.setText(QCoreApplication.translate("MainWindow", u"Gxz", None));
        ___qtablewidgetitem84 = self.THomogenization.verticalHeaderItem(11)
        ___qtablewidgetitem84.setText(QCoreApplication.translate("MainWindow", u"Gyz", None));
        self.BPreviewResults.setText(QCoreApplication.translate("MainWindow", u"Preview Results", None))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.Results), QCoreApplication.translate("MainWindow", u"Results", None))
        self.menuHelp.setTitle(QCoreApplication.translate("MainWindow", u"Help", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
    # retranslateUi

