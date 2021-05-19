CXXFLAGS += -Wall -std=c++11 \
   -Werror \
   -Wextra \
   -Wconversion \
   -Wno-deprecated \
   -Winit-self \
   -Wsign-conversion \
   -Wredundant-decls \
   -Wvla -Wshadow -Wctor-dtor-privacy -Wnon-virtual-dtor -Woverloaded-virtual \
   -Winit-self \
   -Wpointer-arith \
   -Wcast-qual \
   -Wcast-align \
   -Wdouble-promotion \
   -Wold-style-cast -Wno-error=old-style-cast \
   -Wsign-promo \
   -Wswitch-enum \
   -Wswitch-default \
   -Wundef
EXE = pic2d
HEADERS=$(wildcard *.h)
SRCS=$(wildcard *.cpp)
OBJS=Cell.o MacroParticle.o main.o Mesh.o System.o Solve.o Matrix.o Node.o

Cell.o: Node.h Cell.h Constant.h
MacroParticle.o: MacroParticle.h Constant.h
main.o: $(HEADERS) $(SRCS)
Matrix.o: Matrix.h Constant.h
Mesh.o: Node.h Cell.h Cell.cpp Constant.h Mesh.h Mesh.cpp
Node.o: Node.h Constant.h
Solve.o: Solve.h Matrix.h Constant.h
System.o:$(filter-out Matrix.h,$(HEADERS)) $(filter-out main.cpp,$(SRCS))

.PHONY: all build clean distclean run

.DEFAULT_GOAL = build

all: build

build: $(EXE)

$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $^ $(OUTPUT_OPTION)

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) $(EXE)
	$(RM) *.csv *.txt
	$(RM) *.png *.fig

run:
	./$(EXE)
