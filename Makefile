SHELL = /bin/bash

PROG1 = mls_disp main
SOURCES = $(PROG1:%=%.cpp)
OBJECTS = $(PROG1:%=%.o)
HEADERS = $(PROG1).h
TARGET = mls_disp

CXX_FLAGS = -O3 -Wall ${FLAGS_OPENMP} -fopenmp -DPROFILE

all: clean $(TARGET)

intel: CXX_FLAGS = ${FLAGS_OPENMP} -g -DPROFILE ${FLAGS_FAST} -ansi-alias -ipo
intel: clean $(TARGET)

noforce_fast: CXX_FLAGS = ${FLAGS_OPENMP} -g -DNOFORCE -DPROFILE ${FLAGS_FAST} -ansi-alias -ipo
noforce_fast: clean $(TARGET)

noforce: CXX_FLAGS = ${FLAGS_OPENMP} -g -DNOFORCE -DPROFILE 
noforce: clean $(TARGET)

debug: CXX_FLAGS = -O0 -Wall -g -DDEBUG -fopenmp
debug: clean $(TARGET)

clean-all: clean $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) ${CXX_FLAGS} $(OBJECTS) -o $(TARGET)

%.o: %.cpp
	${CXX} ${CXX_FLAGS} -c $<

run:
	./$(TARGET) ./input/Ref4000atoms.input ./config4000.cfg ./input/LJe0Ref4000.lammps

run500:
	./$(TARGET) ./input/Ref500atoms.input ./config500.cfg ./input/LJe0Ref500.lammps

run32:
	./$(TARGET) ./input/Ref32000atoms.input ./config32000.cfg ./input/LJe0Ref32000.lammps

run256:
	./$(TARGET) ./input/Ref256000atoms.input ./config256000.cfg ./input/LJe0Ref256000.lammps

run_sh:
	bash test.sh

clean:
	rm -f $(TARGET)
	rm -f $(OBJECTS)
	rm -f *~



