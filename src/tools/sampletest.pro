TEMPLATE = app
DESTDIR = ../

include( ../common.pri )

LIBS += -L../core -lcore -lIlmImf -lHalf -lpthread
unix:!macx: LIBS += -L/home/jonas/Documents/Doutorado/Documents/Benchmark/Code/build-BenchmarkSystem-Desktop-Debug/ -lRenderingServer -lCore

SOURCES += sampletest.cpp
