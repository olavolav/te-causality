TEMPDIR = maketemp

.PHONY: all 
all: transferentropy

CXXFLAGS = -Wall -O3
# -g
# -O3

# debugger ddd
# valgrind --tool=memcheck ./prog

transferentropy: $(TEMPDIR)/main.o
	g++ -o transferentropy $(TEMPDIR)/main.o $(TEMPDIR)/olav.o

te-global: $(TEMPDIR)/te-global.o
	g++ -o te-global $(TEMPDIR)/te-global.o $(TEMPDIR)/olav.o

te-detrended: $(TEMPDIR)/te-detrended.o
	g++ -o te-detrended $(TEMPDIR)/te-detrended.o $(TEMPDIR)/olav.o

$(TEMPDIR)/main.o: main.cpp $(TEMPDIR)/olav.o
	g++ -c $(CXXFLAGS) -o $(TEMPDIR)/main.o main.cpp

$(TEMPDIR)/te-global.o: te-global.cpp $(TEMPDIR)/olav.o
	g++ -c $(CXXFLAGS) -o $(TEMPDIR)/te-global.o te-global.cpp

$(TEMPDIR)/te-detrended.o: te-detrended.cpp $(TEMPDIR)/olav.o
	g++ -c $(CXXFLAGS) -o $(TEMPDIR)/te-detrended.o te-detrended.cpp

$(TEMPDIR)/olav.o: ../../Simulationen/olav.cpp
	g++ -c $(CXXFLAGS) -o $(TEMPDIR)/olav.o ../../Simulationen/olav.cpp

.PHONY: clean
clean:
	rm -f $(TEMPDIR)/*.o transferentropy te-global
