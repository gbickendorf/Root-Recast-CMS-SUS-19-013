CC=g++
CFLAGS=-c -g -Wall `root-config --cflags` -llibPhysics -I/home/gerrit/Documents/PHD/Tools/MG5_aMC_v3_2_0/ExRootAnalysis -I/home/gerrit/Documents/PHD/Tools/Delphes-3.5.0
LDFLAGS=`root-config --glibs ` -L/home/gerrit/Documents/PHD/Tools/Delphes-3.5.0 -l Delphes -l Physics -L/home/gerrit/Documents/PHD/Tools/MG5_aMC_v3_2_0/ExRootAnalysis -l ExRootAnalysis
SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm ./*~ ./*.o ./main
