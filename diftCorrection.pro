QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += $$PWD/include
INCLUDEPATH += $$PWD/include/opencv2

CONFIG += debug_and_release
CONFIG(debug,debug|release){
TARGET = driftCorrectiond
LIBS += -L$$PWD/lib -ltiffd  -llibfftw3-3 -llibfftw3f-3 -llibfftw3l-3 -lopencv_world480d
}
CONFIG(release,debug|release){
TARGET = driftCorrection
LIBS += -L$$PWD/lib -ltiff  -llibfftw3-3 -llibfftw3f-3 -llibfftw3l-3 -lopencv_world480
}

DESTDIR = $$PWD/bin

SOURCES += \
    accord.cpp \
    imageutils.cpp \
    main.cpp \
    mainwindow.cpp \
    myprocess.cpp \
    process.cpp

HEADERS += \
    accord.h \
    imageutils.h \
    mainwindow.h \
    myprocess.h \
    process.h

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
