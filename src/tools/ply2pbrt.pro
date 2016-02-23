TEMPLATE = app
DESTDIR = ../

include( ../common.pri )

HEADERS +=  ply.h

SOURCES +=  ply.c \
            ply2pbrt.c
