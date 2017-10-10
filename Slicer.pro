QT += core
QT -= gui

CONFIG += c++11 precompile_header

TARGET = Slicer
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

HEADERS += \
    #cgalmachine.h \
    #cgalutility.h \
    #cgaltool.h \
    #pch_cgal.h \
    #geometry.h \
    #tiny_obj_loader.h \
    #iglmachine.h \
    #cgalnef.h \
    tool.h \
    pch.h \
    Matrices.h \
    Vectors.h \
    mesh.h \
    utility.h \
    topo.h \
    datastructure.h

SOURCES += \
    #cgalmachine.cpp \
    #cgalremesher.cpp \
    #cgalholefiller.cpp \
    #objloader.cpp \
    #iglmachine.cpp \
    #cgalnef.cpp \
    #cgaltool.cpp \
    main.cpp \
    Matrices.cpp \
    mesh.cpp \
    geometry.cpp \
    utility.cpp \
    utility_bkup.cpp \
    topo.cpp \
    datastructure.cpp

win32 {

    #INCLUDEPATH += $$PWD/externals/include/CGAL/
    #INCLUDEPATH += $$PWD/externals/include/boost/
    #INCLUDEPATH += $$PWD/externals/include/glm/
    #INCLUDEPATH += $$PWD/externals/include/eigen3/
    #INCLUDEPATH += $$PWD/externals/include/gmp/

    #LIBS += $$PWD/externals/lib/x64/CGAL-vc140-mt-gd-4.10.lib
    #LIBS += $$PWD/externals/lib/x64/CGAL_Core-vc140-mt-gd-4.10.lib
    #LIBS += $$PWD/externals/lib/x64/libboost_system-vc140-mt-gd-1_64.lib
    #LIBS += $$PWD/externals/lib/x64/libboost_thread-vc140-mt-gd-1_64.lib

    #-------------------------------------------------

#INCLUDEPATH += $$PWD/externals/include/glew/
#INCLUDEPATH += $$PWD/externals/include/
#INCLUDEPATH += $$PWD/externals/include/gmp/

#LIBS += $$PWD/externals/lib/x64/glew32s.lib
#LIBS += $$PWD/externals/lib/x64/CGAL-vc140-mt-4.10.lib
#LIBS += $$PWD/externals/lib/x64/CGAL_Core-vc140-mt-4.10.lib
#LIBS += $$PWD/externals/lib/x64/boost_system-vc140-mt-1_64.lib
#LIBS += $$PWD/externals/lib/x64/boost_thread-vc140-mt-1_64.lib
#LIBS += $$PWD/externals/lib/x64/libgmp-10.lib
#LIBS += $$PWD/externals/lib/x64/libmpfr-4.lib
#LIBS += $$PWD/externals/lib/x64/tinyobjloader.lib

    #-------------------------------------------------

    #QMAKE_LFLAGS_DEBUG = /NODEFAULTLIB:libcmt.lib
    #QMAKE_LFLAGS_RELEASE = /NODEFAULTLIB:libcmt.lib

    #-------------------------------------------------

#QMAKE_CXXFLAGS += -bigobj

}
macx {
    #INCLUDEPATH += /usr/local/include
    #INCLUDEPATH += /usr/local/include/eigen3
    #LIBS += -L/usr/local/lib -lCGAL -lboost_thread-mt -lboost_system -lgmp -lmpfr
    #LIBS += $$PWD/externals/lib/mac/libtinyobjloader.a
}
#INCLUDEPATH += $$PWD/externals/include/igl

PRECOMPILED_HEADER = pch.h

#QMAKE_CXXFLAGS_DEBUG -= -g
#QMAKE_CXXFLAGS_RELEASE -= -O2

