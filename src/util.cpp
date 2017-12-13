#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "util.hpp"

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
	else if(key == "Nx"){
	  c.Nx = atoi(c_value);
	}
	else if(key == "Ny")
	  c.Ny = atoi(c_value);
	else if(key == "choix")
	  c.choix = atoi(c_value);
      }
  }
}
