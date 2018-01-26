CC = mpic++
MPI_FLAGS = 
CPP_FLAGS =
DEBUG_FLAGS = -W -Wall -g
OBJ = parallel_conduction.o util.o
MAIN = main.o 
PROG = run
HEADER_DIRS = header Eigen/Eigen
HEADER_FLAGS = $(foreach dir, $(HEADER_DIRS),-I$(dir))

vpath %.cpp test src
vpath %.hpp header

all:$(PROG)

debug : CPP_FLAGS += $(DEBUG_FLAGS)
debug : clean
debug : all

$(PROG):$(OBJ) $(MAIN)
	$(CC) $^ -o $@

%.o: %.cpp %.hpp
	$(CC) $(CPP_FLAGS) $(HEADER_FLAGS) -c $< -o $@

%.o: %.cpp
	$(CC) $(CPP_FLAGS) $(HEADER_FLAGS) -c $< -o $@

test: test_unitaire

test_unitaire : $(OBJ) test_unitaire.o
	$(CC) $(CPP_FLAGS) $(HEADER_FLAGS) -o $@ $^ $(LDFLAGS)

courbe :
	./script.sh
	python plot_csv.py
clean:
	rm -f *.o *~ $(PROG) output test_unitaire *.csv
