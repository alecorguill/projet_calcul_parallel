FC = mpif90
OBJ = conjugate.o fonction.o main.o

PROG = run

$(PROG):$(OBJ)
	$(FC) $^ -o $@

%.o:%.f90
	$(FC) -c $< -o $@

conjugate.o:fonction.o

main.o:conjugate.o fonction.o

all:$(PROG)

clean:
	rm -f *.o *.mod *~ $(PROG)
