#ifndef UTIL_H
#define UTIL_H

typedef struct {
  double Lx;
  double Ly;
  int Nx;
  int Ny;
  double D;
  int choix;
  int np;
  double dx;
  double dy;
  double dt;
} config_t;

void parse_file(char* filename, config_t& c);

#endif
