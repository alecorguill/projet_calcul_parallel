CC = g++
CPP_FLAGS =
DEBUG_FLAGS = -W -Wall -g
OBJ = parallel_conduction.o
MAIN = main.o 
PROG = run
HEADER_DIR = header

vpath %.cpp test src
vpath %.hpp header

all:$(PROG)

debug : CPP_FLAGS += $(DEBUG_FLAGS)
debug : clean
debug : all

$(PROG):$(OBJ) $(MAIN)
	$(CC) $^ -o $@

%.o: %.cpp
	$(CC) $(CPP_FLAGS) -I$(HEADER_DIR) -c $< -o $@

test: test_unitaire

test_unitaire : $(OBJ) test_unitaire.o
	$(CC) $(CPP_FLAGS) -I$(HEADER_DIR) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o *~ $(PROG)
