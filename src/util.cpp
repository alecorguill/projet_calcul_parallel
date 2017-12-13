#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <Eigen>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>


#include "util.hpp"

#define LINESIZE 82

using namespace std;

void parse_file(char* filename, config_t&c){
  ifstream file(filename);

  string line;
  while(getline(file, line)){
    stringstream iss(line);
    string key;
    if( getline(iss, key, '=') )
      {
	string value;
	getline(iss, value);
	const char * c_value = value.c_str();
	if(key == "Lx")
	  c.Lx = atof(c_value);
	else if(key == "Ly")
	  c.Ly = atof(c_value);
	else if(key == "D")
	  c.D = atof(c_value);
	else if(key == "Nx")
	  c.Nx = atoi(c_value);
	else if(key == "Ny")
	  c.Ny = atoi(c_value);
	else if(key == "choix")
	  c.choix = atoi(c_value);
      }
  }
}

void log_result(int output, Eigen::VectorXd& u){
  int rank;
  off_t offset;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //MPI_Comm_size(MPI_COMM_WORLD, &size);
  int size =  u.size();
  offset = LINESIZE * size * rank;
  lseek(output,offset,SEEK_SET);
  for(int i = 0; i<size; ++i){
    dprintf(output, "%lf\n", u(i));
  }
}
