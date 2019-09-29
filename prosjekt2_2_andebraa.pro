TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    functions.cpp \
    main.cpp


INCLUDEPATH += C:\Users\ander\fys3150\armadillo-9.200.5\include
DEPENDPATH += C:\Users\ander\fys3150\armadillo-9.200.5\include


LIBS += \
    -LC:\Users\ander\fys3150\armadillo-9.200.5\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT

HEADERS += \
    header.h
