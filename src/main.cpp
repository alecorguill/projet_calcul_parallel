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
  /* Parsing config file */
  config_t c;
  parse_file(argv[1],c);
  

  return 1;
}
