#-------------------------------------------------
#
# Project created by QtCreator 2019-01-04T00:53:30
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = VRM_Editor
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Linux is simple
#--------------------
LIBS += -lnetcdf 
#
# For My MAC
#-----------------
#INCLUDEPATH += /Users/patc/Projects/Qt/VRM_Editor/netcdf/include
#LIBS += -L/Users/patc/Projects/Qt/VRM_Editor/netcdf/lib -lnetcdf -lhdf5_hl -lhdf5 -lm /Users/patc/Projects/Qt/VRM_Editor/netcdf/lib/libz.a

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    lonbar.cpp \
    vrmview.cpp \
    latbar.cpp \
    vrmwin.cpp \
    CubeGrid.cpp \
    SQuadGen/GridElements.cpp \
    RefineCubeGrid.cpp \
    RefinementCube.cpp \
    RefinementMap.cpp \
    SQuadGen/SpringDynamics.cpp \
    SQuadGen/Tessellate.cpp \
    SQuadGen/IcosahedralFlagGrid.cpp \
    SQuadGen/CubedSphereGrid.cpp \
    SQuadGen/RefinementTemplateCUBIT.cpp \
    SQuadGen/RefinementTemplateLOWCONN.cpp \
    SQuadGen/RefinementTemplateLOWCONNOLD.cpp \
    vrmgrid.cpp \
    cubegriditem.cpp \
    polyedititem.cpp \
    polynodeitem.cpp \
    rectedititem.cpp \
    NeighborAdjustments.cpp \
    SQuadGen/netcdf.cpp \
    SQuadGen/ncvalues.cpp \
    ReadWrite.cpp \
    gridinfoitem.cpp

HEADERS += \
        mainwindow.h \
    lonbar.h \
    vrmview.h \
    latbar.h \
    vrmwin.h \
    CubeGrid.h \
    SQuadGen/GridElements.h \
    SQuadGen/Exception.h \
    RefineCubeGrid.h \
    RefinementCube.h \
    RefinementMap.h \
    SQuadGen/DataMatrix3D.h \
    SQuadGen/CubedSphereTrans.h \
    SQuadGen/MathHelper.h \
    SQuadGen/DataVector.h \
    SQuadGen/DataMatrix.h \
    SQuadGen/SpringDynamics.h \
    SQuadGen/IcosahedralFlagGrid.h \
    SQuadGen/CubedSphereGrid.h \
    SQuadGen/Tessellate.h \
    SQuadGen/RefinementTemplateCUBIT.h \
    SQuadGen/RefinementTemplateLOWCONN.h \
    SQuadGen/RefinementTemplateLOWCONNOLD.h \
    vrmgrid.h \
    cubegriditem.h \
    polyedititem.h \
    polynodeitem.h \
    rectedititem.h \
    NeighborAdjustments.h \
    SQuadGen/netcdfcpp.h \
    SQuadGen/ncvalues.h \
    ReadWrite.h \
    gridinfoitem.h

FORMS += \
        mainwindow.ui

RESOURCES += \
    vrm_resources.qrc
