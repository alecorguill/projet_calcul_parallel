#include <mpi.h>
#include "parallel_conduction.hpp"
#include "main.hpp"
#include "util.hpp"


#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char** argv){
  if(argc < 2){
    fprintf(stderr, "ERROR COMMAND LINE\n");
    fprintf(stderr, "USAGE : %s config_file\n", argv[0]);
    exit(0);
  }
  MPI_Init(NULL,NULL);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  /* Parsing config file */
  config_t c;
  parse_file(argv[1],c);
  c.np = size;
  c.dx = c.Lx/(c.Nx+1);
  c.dy = c.Ly/(c.Ny+1);
  c.dt = 0.1;
  
  
  MPI_Finalize();

  return 0;
}
