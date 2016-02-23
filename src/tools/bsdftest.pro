include( ../common.pri )

DESTDIR = ../

PRE_TARGETDEPS += ../core/libcore.a

LIBS += -L../core -lcore -lIlmImf -lHalf -lpthread -llapacke -llapack -lblas -lgfortran

SOURCES += bsdftest.cpp

unix:!macx: LIBS += -L/home/jonas/Documents/Doutorado/Documents/Benchmark/Code/build-BenchmarkSystem-Desktop-Debug/ -lRenderingServer -lCore
