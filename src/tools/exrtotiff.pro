include( ../common.pri )

LIBS += -L../core -lcore -lIlmImf -lHalf -lpthread -ltiff

SOURCES += exrtotiff.cpp
