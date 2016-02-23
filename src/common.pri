#CONFIG -= qt
#DEFINES -= QT_WEBKIT
QT += core network

CONFIG(debug, debug|release) {
} else {
    DEFINES += NDEBUG
}

DEFINES +=  PBRT_HAS_OPENEXR \
            PBRT_PROBES_NONE    # PBRT_PROBES_DTRACE

QMAKE_CXXFLAGS += -std=c++0x -g

DEPENDPATH += ../

INCLUDEPATH +=  ../ \
                ../core \
                /usr/include/OpenEXR
