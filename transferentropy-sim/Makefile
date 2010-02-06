# ----------------------------------- #
# Makefile for Sim                    #
# ----------------------------------- #

# path to Sim libraries
SIM_PATH = ../../Simulationen/SimKernel
# path to object files
MAKETEMP_PATH = ./maketemp

# single or parallel
MPI = off
EXE = transferentropy-sim

ifeq ($(MPI), on)
CC     = mpic++
CFLAGS = -DSIM_MPI -O3
else
CC     = g++
CFLAGS = -O3
# -O2 
endif

LDFLAGS  = -lgsl -lm -L.
OBJ  = $(MAKETEMP_PATH)/olav.o $(MAKETEMP_PATH)/main.o $(MAKETEMP_PATH)/expression.o $(MAKETEMP_PATH)/expression_basic.o $(MAKETEMP_PATH)/expression_sim.o $(MAKETEMP_PATH)/expression_extra.o $(MAKETEMP_PATH)/expression_parser.o $(MAKETEMP_PATH)/sim_io_manager.o $(MAKETEMP_PATH)/sim_io.o $(MAKETEMP_PATH)/sim_signals.o $(MAKETEMP_PATH)/sim.o

EXPRH = $(SIM_PATH)/expression.h $(SIM_PATH)/expression_basic.h $(SIM_PATH)/expression_sim.h $(SIM_PATH)/expression_extra.h $(SIM_PATH)/expression_parser.h
SIMH =  $(SIM_PATH)/sim.h $(SIM_PATH)/sim_io.h $(SIM_PATH)/sim_signals.h $(SIM_PATH)/sim_control.h $(SIM_PATH)/sim_io_manager.h 

ifeq ($(MPI), on)
OBJ = $(MAKETEMP_PATH)/sim_mpi.o $(OBJ)
SIMH = $(SIM_PATH)/sim_mpi.h $(SIM_PATH)/sim_control_mpi.h  $(SIMH)
endif


%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

$(MAKETEMP_PATH)/%.o : $(SIM_PATH)/%.cpp
	$(CC) -c $(CFLAGS) $< -o $@

$(MAKETEMP_PATH)/%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

all : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LDFLAGS)

# Dependencies

$(MAKETEMP_PATH)/main.o : main.cpp $(EXPRH) $(SIMH) $(MAKETEMP_PATH)/olav.o
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


clean:
	rm -f $(EXE) maketemp/* maketemp/*~

