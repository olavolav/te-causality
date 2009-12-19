TEMPDIR = maketemp

.PHONY: all 
all: granger

CXXFLAGS = -Wall -O3 -I/usr/local/include
# -g
# -O3

# debugger ddd
# valgrind --tool=memcheck ./prog

granger: $(TEMPDIR)/main.o
	g++ -lgsl -lm -lgslcblas -o granger $(TEMPDIR)/main.o $(TEMPDIR)/olav.o

$(TEMPDIR)/main.o: main.cpp $(TEMPDIR)/olav.o
	g++ -c $(CXXFLAGS) -o $(TEMPDIR)/main.o main.cpp

$(TEMPDIR)/olav.o: ../../Simulationen/olav.cpp
	g++ -c $(CXXFLAGS) -o $(TEMPDIR)/olav.o ../../Simulationen/olav.cpp

.PHONY: clean
clean:
	rm -f $(TEMPDIR)/*.o granger
