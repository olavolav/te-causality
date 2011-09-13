# ----------------------------------- #
# Makefile for Sim										#
# ----------------------------------- #

# path to Sim libraries
SIM_PATH = ../../../Sonstiges/SimKernel
# path to object files
MAKETEMP_PATH = ./maketemp

# Georgs debugging:
# CFLAGS plus -g

# debugger ddd (or gdb on Mac OSX)
# valgrind --tool=memcheck ./prog


# single or parallel
MPI = off

ifeq ($(MPI), on)
CC		 = mpic++
CFLAGS = -DSIM_MPI -O3
else
CC		 = g++
CFLAGS = -O3
# CFLAGS = -g
endif

LDFLAGS	 = -lgslcblas -lgsl -lm -L.
LDFLAGS_YAML = -lyaml-cpp
LDFLAGS_FLANN = -lflann_cpp
OBJBASE = $(MAKETEMP_PATH)/olav.o $(MAKETEMP_PATH)/te-datainit.o $(MAKETEMP_PATH)/expression.o $(MAKETEMP_PATH)/expression_basic.o $(MAKETEMP_PATH)/expression_sim.o $(MAKETEMP_PATH)/expression_extra.o $(MAKETEMP_PATH)/expression_parser.o $(MAKETEMP_PATH)/sim_io_manager.o $(MAKETEMP_PATH)/sim_io.o $(MAKETEMP_PATH)/sim_signals.o $(MAKETEMP_PATH)/sim.o

EXPRH = $(SIM_PATH)/expression.h $(SIM_PATH)/expression_basic.h $(SIM_PATH)/expression_sim.h $(SIM_PATH)/expression_extra.h $(SIM_PATH)/expression_parser.h
SIMHBASE =	$(SIM_PATH)/sim.h $(SIM_PATH)/sim_io.h $(SIM_PATH)/sim_signals.h $(SIM_PATH)/sim_control.h $(SIM_PATH)/sim_io_manager.h 

ifeq ($(MPI), on)
OBJ = $(MAKETEMP_PATH)/sim_mpi.o $(OBJBASE)
SIMH = $(SIM_PATH)/sim_mpi.h $(SIM_PATH)/sim_control_mpi.h $(SIMHBASE)
else
OBJ = $(OBJBASE)
SIMH = $(SIMHBASE)
endif


%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

$(MAKETEMP_PATH)/%.o : $(SIM_PATH)/%.cpp
	$(CC) -c $(CFLAGS) $< -o $@

$(MAKETEMP_PATH)/%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@
	
# .PHONY: all 
#		all: transferentropy-sim te-global-sim te-extended

# transferentropy : $(OBJ) $(MAKETEMP_PATH)/main.o
#   $(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/main.o -o transferentropy $(LDFLAGS)

# te-global : $(OBJ) $(MAKETEMP_PATH)/te-global.o
#   $(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/te-global.o -o te-global $(LDFLAGS)

# te-extended : $(OBJ) $(MAKETEMP_PATH)/te-extended.o
#   $(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/te-extended.o -o te-extended $(LDFLAGS) $(LDFLAGS_YAML)

te-extended : $(OBJ) $(MAKETEMP_PATH)/te-extended.o
	$(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/te-extended.o $(MAKETEMP_PATH)/multidimarray.o -o te-extended $(LDFLAGS) $(LDFLAGS_YAML)

mi : $(OBJ) $(MAKETEMP_PATH)/mi.o
	$(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/mi.o $(MAKETEMP_PATH)/multidimarray.o -o mi $(LDFLAGS) $(LDFLAGS_YAML)

te-binless-Kraskov : $(OBJ) $(MAKETEMP_PATH)/te-binless-Kraskov.o
	$(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/te-binless-Kraskov.o -o te-binless-Kraskov $(LDFLAGS) $(LDFLAGS_YAML) $(LDFLAGS_FLANN)

te-binless-Leonenko : $(OBJ) $(MAKETEMP_PATH)/te-binless-Leonenko.o
	$(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/te-binless-Leonenko.o -o te-binless-Leonenko $(LDFLAGS) $(LDFLAGS_YAML) $(LDFLAGS_FLANN)

# te-extended-localcond : $(OBJ) $(MAKETEMP_PATH)/te-extended-localcond.o
#   $(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/te-extended-localcond.o -o te-extended-localcond $(LDFLAGS)

# te-separated : $(OBJ) $(MAKETEMP_PATH)/te-separated.o
#   $(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/te-separated.o -o te-separated $(LDFLAGS)
	
# te-iterated : $(OBJ) $(MAKETEMP_PATH)/te-iterated.o
#   $(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/te-iterated.o -o te-iterated $(LDFLAGS)

xcorr : $(OBJ) $(MAKETEMP_PATH)/xcorr.o
	$(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/xcorr.o -o xcorr $(LDFLAGS) $(LDFLAGS_YAML)

granger : $(OBJ) $(MAKETEMP_PATH)/granger.o
	$(CC) $(CFLAGS) $(OBJ) $(MAKETEMP_PATH)/granger.o -o granger $(LDFLAGS) $(LDFLAGS_YAML)


# Dependencies

# $(MAKETEMP_PATH)/main.o : main.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
# $(MAKETEMP_PATH)/te-global.o : te-global.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
# $(MAKETEMP_PATH)/te-extended.o : te-extended.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
$(MAKETEMP_PATH)/te-extended.o : te-extended.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o $(MAKETEMP_PATH)/multidimarray.o
$(MAKETEMP_PATH)/mi.o : mi.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o $(MAKETEMP_PATH)/multidimarray.o
$(MAKETEMP_PATH)/te-binless-Kraskov.o : te-binless-Kraskov.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
$(MAKETEMP_PATH)/te-binless-Leonenko.o : te-binless-Leonenko.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
# $(MAKETEMP_PATH)/te-extended-localcond.o : te-extended-localcond.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
# $(MAKETEMP_PATH)/te-separated.o : te-extended.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
# $(MAKETEMP_PATH)/te-iterated.o : te-iterated.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
$(MAKETEMP_PATH)/xcorr.o : xcorr.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
$(MAKETEMP_PATH)/granger.o : granger.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o

$(MAKETEMP_PATH)/sim.o : $(SIM_PATH)/sim.cpp $(SIMH)
$(MAKETEMP_PATH)/sim_io.o : $(SIM_PATH)/sim_io.cpp $(SIM_PATH)/sim_io.h 
$(MAKETEMP_PATH)/sim_mpi.o : $(SIM_PATH)/sim_mpi.cpp $(SIM_PATH)/sim_mpi.h
$(MAKETEMP_PATH)/sim_signals.o : $(SIM_PATH)/sim_signals.cpp $(SIM_PATH)/sim_signals.h
$(MAKETEMP_PATH)/sim_io_manager.o : $(SIM_PATH)/sim_io_manager.cpp $(SIM_PATH)/sim_io_manager.h $(SIM_PATH)/sim_signals.h

$(MAKETEMP_PATH)/expression.o : $(SIM_PATH)/expression.cpp $(SIM_PATH)/expression.h
$(MAKETEMP_PATH)/expression_basic.o : $(SIM_PATH)/expression_basic.cpp $(SIM_PATH)/expression_basic.h $(SIM_PATH)/expression.h
$(MAKETEMP_PATH)/expression_sim.o : $(SIM_PATH)/expression_sim.cpp $(SIM_PATH)/expression_sim.h $(SIM_PATH)/expression.h $(SIM_PATH)/expression_basic.h 
$(MAKETEMP_PATH)/expression_extra.o : $(SIM_PATH)/expression_extra.cpp $(SIM_PATH)/expression_extra.h $(SIM_PATH)/expression.h $(SIM_PATH)/expression_basic.h 
$(MAKETEMP_PATH)/expression_parser.o : $(SIM_PATH)/expression_parser.cpp $(SIM_PATH)/expression_parser.h $(SIM_PATH)/expression.h $(SIM_PATH)/expression_basic.h $(SIM_PATH)/expression_sim.h $(SIM_PATH)/expression_extra.h

$(MAKETEMP_PATH)/olav.o: ../../Simulationen/olav.cpp
	$(CC) -c $(CFLAGS) -o $(MAKETEMP_PATH)/olav.o ../../Simulationen/olav.cpp
$(MAKETEMP_PATH)/te-datainit.o: ../te-datainit.cpp
	$(CC) -c $(CFLAGS) -o $(MAKETEMP_PATH)/te-datainit.o ../te-datainit.cpp
$(MAKETEMP_PATH)/multidimarray.o: ../multidimarray.cpp
	$(CC) -c $(CFLAGS) -o $(MAKETEMP_PATH)/multidimarray.o ../multidimarray.cpp


clean:
	rm -f transferentropy te-global te-separated te-extended te-extended-multidim mi te-binless-Kraskov te-binless-Leonenko te-extended-localcond te-iterated xcorr granger maketemp/*
	rm -rf *.dSYM ../*.dSYM
