TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        functions.cpp \
        main.cpp

#HEADERS +=


INCLUDEPATH += ../agrolib/mathFunctions 

SUBDIRS = ../agrolib/mathFunctions 

unix:{
    CONFIG(debug, debug|release) {
        TARGET = debug/marquardt_nDimension
    } else {
        TARGET = release/marquardt_nDimension
    }
}
win32:{
    TARGET = marquardt_nDimension
}

CONFIG(release, debug|release) {
    LIBS += -L../agrolib/mathFunctions/release -lmathFunctions
} else {
    LIBS += -L../agrolib/mathFunctions/debug -lmathFunctions
}

HEADERS += \
    functions.h
