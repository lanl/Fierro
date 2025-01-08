# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'FIERRO_SetuppuYJnm.ui'
##
## Created by: Qt User Interface Compiler version 6.7.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QDialog, QFrame, QGridLayout,
    QHBoxLayout, QLabel, QPushButton, QRadioButton,
    QSizePolicy, QStackedWidget, QVBoxLayout, QWidget)
import IconResourceFile_rc

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        if not Dialog.objectName():
            Dialog.setObjectName(u"Dialog")
        Dialog.resize(627, 450)
        Dialog.setMinimumSize(QSize(0, 450))
        self.verticalLayout = QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.frame_5 = QFrame(Dialog)
        self.frame_5.setObjectName(u"frame_5")
        self.frame_5.setMaximumSize(QSize(16777215, 16777215))
        self.frame_5.setFrameShape(QFrame.Shape.Panel)
        self.frame_5.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_2 = QGridLayout(self.frame_5)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.label_3 = QLabel(self.frame_5)
        self.label_3.setObjectName(u"label_3")
        self.label_3.setMaximumSize(QSize(100, 100))
        self.label_3.setPixmap(QPixmap(u":/Blue Icons/Blue Icons/UserProfile.svg"))
        self.label_3.setScaledContents(True)

        self.gridLayout_2.addWidget(self.label_3, 0, 0, 1, 1, Qt.AlignmentFlag.AlignTop)

        self.frame_6 = QFrame(self.frame_5)
        self.frame_6.setObjectName(u"frame_6")
        self.frame_6.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_6.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_3 = QVBoxLayout(self.frame_6)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 5, 0, 5)
        self.label = QLabel(self.frame_6)
        self.label.setObjectName(u"label")

        self.verticalLayout_3.addWidget(self.label)

        self.frame_2 = QFrame(self.frame_6)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_2.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout = QHBoxLayout(self.frame_2)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.BUser = QRadioButton(self.frame_2)
        self.BUser.setObjectName(u"BUser")
        self.BUser.setEnabled(True)
        self.BUser.setChecked(True)

        self.horizontalLayout.addWidget(self.BUser)

        self.BDeveloper = QRadioButton(self.frame_2)
        self.BDeveloper.setObjectName(u"BDeveloper")
        self.BDeveloper.setChecked(False)

        self.horizontalLayout.addWidget(self.BDeveloper)


        self.verticalLayout_3.addWidget(self.frame_2, 0, Qt.AlignmentFlag.AlignHCenter)

        self.Configuration = QStackedWidget(self.frame_6)
        self.Configuration.setObjectName(u"Configuration")
        self.Configuration.setFrameShape(QFrame.Shape.NoFrame)
        self.page = QWidget()
        self.page.setObjectName(u"page")
        self.verticalLayout_5 = QVBoxLayout(self.page)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.Configuration.addWidget(self.page)
        self.page_2 = QWidget()
        self.page_2.setObjectName(u"page_2")
        self.verticalLayout_4 = QVBoxLayout(self.page_2)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.frame_7 = QFrame(self.page_2)
        self.frame_7.setObjectName(u"frame_7")
        self.frame_7.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_7.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_3 = QGridLayout(self.frame_7)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.gridLayout_3.setContentsMargins(0, 0, -1, 0)
        self.INFierroMeshBuilder = QLabel(self.frame_7)
        self.INFierroMeshBuilder.setObjectName(u"INFierroMeshBuilder")
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.INFierroMeshBuilder.sizePolicy().hasHeightForWidth())
        self.INFierroMeshBuilder.setSizePolicy(sizePolicy)
        self.INFierroMeshBuilder.setMinimumSize(QSize(252, 0))
        self.INFierroMeshBuilder.setWordWrap(True)

        self.gridLayout_3.addWidget(self.INFierroMeshBuilder, 1, 2, 1, 1)

        self.LFierroMeshBuilder = QLabel(self.frame_7)
        self.LFierroMeshBuilder.setObjectName(u"LFierroMeshBuilder")

        self.gridLayout_3.addWidget(self.LFierroMeshBuilder, 1, 0, 1, 1)

        self.BFierroParallelExplicit = QPushButton(self.frame_7)
        self.BFierroParallelExplicit.setObjectName(u"BFierroParallelExplicit")
        self.BFierroParallelExplicit.setAutoDefault(False)

        self.gridLayout_3.addWidget(self.BFierroParallelExplicit, 3, 3, 1, 1)

        self.LExecutables = QLabel(self.frame_7)
        self.LExecutables.setObjectName(u"LExecutables")

        self.gridLayout_3.addWidget(self.LExecutables, 0, 2, 1, 1)

        self.INFierroVoxelizer = QLabel(self.frame_7)
        self.INFierroVoxelizer.setObjectName(u"INFierroVoxelizer")
        sizePolicy.setHeightForWidth(self.INFierroVoxelizer.sizePolicy().hasHeightForWidth())
        self.INFierroVoxelizer.setSizePolicy(sizePolicy)
        self.INFierroVoxelizer.setMinimumSize(QSize(252, 0))
        self.INFierroVoxelizer.setWordWrap(True)

        self.gridLayout_3.addWidget(self.INFierroVoxelizer, 2, 2, 1, 1)

        self.BFierroMeshBuilder = QPushButton(self.frame_7)
        self.BFierroMeshBuilder.setObjectName(u"BFierroMeshBuilder")
        self.BFierroMeshBuilder.setAutoDefault(False)

        self.gridLayout_3.addWidget(self.BFierroMeshBuilder, 1, 3, 1, 1)

        self.INFierroParallelExplicit = QLabel(self.frame_7)
        self.INFierroParallelExplicit.setObjectName(u"INFierroParallelExplicit")
        sizePolicy.setHeightForWidth(self.INFierroParallelExplicit.sizePolicy().hasHeightForWidth())
        self.INFierroParallelExplicit.setSizePolicy(sizePolicy)
        self.INFierroParallelExplicit.setMinimumSize(QSize(252, 0))
        self.INFierroParallelExplicit.setWordWrap(True)

        self.gridLayout_3.addWidget(self.INFierroParallelExplicit, 3, 2, 1, 1)

        self.BFierroVoxelizer = QPushButton(self.frame_7)
        self.BFierroVoxelizer.setObjectName(u"BFierroVoxelizer")
        self.BFierroVoxelizer.setAutoDefault(False)

        self.gridLayout_3.addWidget(self.BFierroVoxelizer, 2, 3, 1, 1)

        self.LFierroParallelExplicit = QLabel(self.frame_7)
        self.LFierroParallelExplicit.setObjectName(u"LFierroParallelExplicit")

        self.gridLayout_3.addWidget(self.LFierroParallelExplicit, 3, 0, 1, 1)

        self.LFierroVoxelizer = QLabel(self.frame_7)
        self.LFierroVoxelizer.setObjectName(u"LFierroVoxelizer")

        self.gridLayout_3.addWidget(self.LFierroVoxelizer, 2, 0, 1, 1)

        self.Levpfft = QLabel(self.frame_7)
        self.Levpfft.setObjectName(u"Levpfft")

        self.gridLayout_3.addWidget(self.Levpfft, 4, 0, 1, 1)

        self.BFierroEvpfft = QPushButton(self.frame_7)
        self.BFierroEvpfft.setObjectName(u"BFierroEvpfft")

        self.gridLayout_3.addWidget(self.BFierroEvpfft, 4, 3, 1, 1)

        self.INFierroEvpfft = QLabel(self.frame_7)
        self.INFierroEvpfft.setObjectName(u"INFierroEvpfft")
        sizePolicy.setHeightForWidth(self.INFierroEvpfft.sizePolicy().hasHeightForWidth())
        self.INFierroEvpfft.setSizePolicy(sizePolicy)
        self.INFierroEvpfft.setMinimumSize(QSize(252, 0))
        self.INFierroEvpfft.setWordWrap(True)

        self.gridLayout_3.addWidget(self.INFierroEvpfft, 4, 2, 1, 1)


        self.verticalLayout_4.addWidget(self.frame_7)

        self.Configuration.addWidget(self.page_2)

        self.verticalLayout_3.addWidget(self.Configuration)


        self.gridLayout_2.addWidget(self.frame_6, 0, 1, 1, 1)


        self.verticalLayout.addWidget(self.frame_5)

        self.frame_3 = QFrame(Dialog)
        self.frame_3.setObjectName(u"frame_3")
        self.frame_3.setFrameShape(QFrame.Shape.Panel)
        self.frame_3.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout = QGridLayout(self.frame_3)
        self.gridLayout.setObjectName(u"gridLayout")
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.frame_4 = QFrame(self.frame_3)
        self.frame_4.setObjectName(u"frame_4")
        self.frame_4.setFrameShape(QFrame.Shape.NoFrame)
        self.frame_4.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_2 = QVBoxLayout(self.frame_4)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.LDirections = QLabel(self.frame_4)
        self.LDirections.setObjectName(u"LDirections")
        self.LDirections.setWordWrap(True)

        self.verticalLayout_2.addWidget(self.LDirections)

        self.BSelectWD = QPushButton(self.frame_4)
        self.BSelectWD.setObjectName(u"BSelectWD")
        self.BSelectWD.setAutoDefault(False)

        self.verticalLayout_2.addWidget(self.BSelectWD, 0, Qt.AlignmentFlag.AlignHCenter)

        self.frame = QFrame(self.frame_4)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.Shape.NoFrame)
        self.frame.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout_4 = QGridLayout(self.frame)
        self.gridLayout_4.setObjectName(u"gridLayout_4")
        self.gridLayout_4.setContentsMargins(0, 0, 0, 0)
        self.label_4 = QLabel(self.frame)
        self.label_4.setObjectName(u"label_4")

        self.gridLayout_4.addWidget(self.label_4, 0, 0, 1, 1)

        self.INCurrentWD = QLabel(self.frame)
        self.INCurrentWD.setObjectName(u"INCurrentWD")
        sizePolicy.setHeightForWidth(self.INCurrentWD.sizePolicy().hasHeightForWidth())
        self.INCurrentWD.setSizePolicy(sizePolicy)
        self.INCurrentWD.setMinimumSize(QSize(280, 0))
        self.INCurrentWD.setWordWrap(True)

        self.gridLayout_4.addWidget(self.INCurrentWD, 0, 1, 1, 1)


        self.verticalLayout_2.addWidget(self.frame, 0, Qt.AlignmentFlag.AlignTop)


        self.gridLayout.addWidget(self.frame_4, 0, 1, 1, 1)

        self.label_2 = QLabel(self.frame_3)
        self.label_2.setObjectName(u"label_2")
        self.label_2.setMaximumSize(QSize(100, 100))
        self.label_2.setPixmap(QPixmap(u":/Blue Icons/Blue Icons/File.svg"))
        self.label_2.setScaledContents(True)

        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)


        self.verticalLayout.addWidget(self.frame_3)

        self.BConfirm = QPushButton(Dialog)
        self.BConfirm.setObjectName(u"BConfirm")
        self.BConfirm.setAutoDefault(False)

        self.verticalLayout.addWidget(self.BConfirm, 0, Qt.AlignmentFlag.AlignRight)


        self.retranslateUi(Dialog)

        self.Configuration.setCurrentIndex(0)
        self.BConfirm.setDefault(True)


        QMetaObject.connectSlotsByName(Dialog)
    # setupUi

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QCoreApplication.translate("Dialog", u"Fierro Setup", None))
        self.label_3.setText("")
        self.label.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p align=\"center\"><span style=\" font-style:italic;\">Select profile configuration:</span></p></body></html>", None))
        self.BUser.setText(QCoreApplication.translate("Dialog", u"User", None))
        self.BDeveloper.setText(QCoreApplication.translate("Dialog", u"Developer", None))
        self.INFierroMeshBuilder.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p>/</p></body></html>", None))
        self.LFierroMeshBuilder.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p>fierro-mesh-builder:</p></body></html>", None))
        self.BFierroParallelExplicit.setText(QCoreApplication.translate("Dialog", u"Select", None))
        self.LExecutables.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p align=\"center\"><span style=\" font-style:italic;\">Specifiy the path to the build executables:</span></p></body></html>", None))
        self.INFierroVoxelizer.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p>/</p></body></html>", None))
        self.BFierroMeshBuilder.setText(QCoreApplication.translate("Dialog", u"Select", None))
        self.INFierroParallelExplicit.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p>/</body></html>", None))
        self.BFierroVoxelizer.setText(QCoreApplication.translate("Dialog", u"Select", None))
        self.LFierroParallelExplicit.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p>fierro-parallel-explicit:</p></body></html>", None))
        self.LFierroVoxelizer.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p>fierro-voxelizer:</p></body></html>", None))
        self.Levpfft.setText(QCoreApplication.translate("Dialog", u"evpfft:", None))
        self.BFierroEvpfft.setText(QCoreApplication.translate("Dialog", u"Select", None))
        self.INFierroEvpfft.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p>/</p></body></html>", None))
        self.LDirections.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p align=\"center\"><span style=\" font-style:italic;\">Specify the directory path for saving input and output files:</span></p></body></html>", None))
        self.BSelectWD.setText(QCoreApplication.translate("Dialog", u"Select Working Directory", None))
        self.label_4.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p><span style=\" font-weight:700;\">Current Working Directory:</span></p></body></html>", None))
        self.INCurrentWD.setText(QCoreApplication.translate("Dialog", u"<html><head/><body><p>/temporary</p></body></html>", None))
        self.label_2.setText("")
        self.BConfirm.setText(QCoreApplication.translate("Dialog", u"Confirm", None))
    # retranslateUi

