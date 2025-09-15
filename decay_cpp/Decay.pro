#-------------------------------------------------
#
# Project created by QtCreator 2025-07-23T18:18:07
#
#-------------------------------------------------

QT       -=core gui
QT +=core

TARGET = decaysum7
TEMPLATE = lib
CONFIG+=shared
CONFIG+=c++17

DEFINES += DECAY_LIBRARY

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.



SOURCES +=
SOURCES += \
    ccalculateelevation.cpp \
    ccalculateheight.cpp \
    cdegtorad.cpp \
    cedgerectangle1.cpp \
    cedgerectangle2.cpp \
    cgreatcirclesolver.cpp \
    cisplanebetween.cpp \
    cisvaluebetween.cpp \
    cisvaluebetween1.cpp \
    cloudatt.cpp \
    decaysum.cpp \
    fogatt.cpp \
    freespaceatt.cpp \
    gasatt.cpp \
    hazeatt.cpp \
    rainatt.cpp \
    reconfig.cpp \
    sciatt.cpp \
    snowatt.cpp
HEADERS += \
    ccalculateelevation.h \
    ccalculateheight.h \
    cdegtorad.h \
    cedgerectangle1.h \
    cedgerectangle2.h \
    cgreatcirclesolver.h \
    cisplanebetween.h \
    cisvaluebetween.h \
    cisvaluebetween1.h \
    cloudatt.h \
    decaysum.h \
    enum.h \
    fogatt.h \
    freespaceatt.h \
    gasatt.h \
    hazeatt.h \
    rainatt.h \
    reconfig.h \
    sciatt.h \
    snowatt.h


unix {
    target.path = /usr/lib
    INSTALLS += target
}
