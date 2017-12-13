#ifndef UTIL_H
#define UTIL_H

#include <Eigen>

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
  int output;
} config_t;

void parse_file(char* filename, config_t& c);
void log_result(int output, Eigen::VectorXd& u);
#endif
