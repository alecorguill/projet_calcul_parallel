/* parallel_conduction.cpp */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>

double f(double dx, double dy, double Lx, double Ly, int i, int j, double dt, int choix){
  double x = i*dx;
  double y = j*dy;
  if(choix == 0)
    return 2*(x-x*x+y-y*y);
  else if(choix == 1)
    return sin(x)+cos(y);
  else if(choix == 2)
    return exp(-pow((x-Lx*0.5),2))*exp(-pow((y-Ly*0.5),2))*cos((M_PI*0.5)*dt);
  else{
    fprintf(stderr, "CHOIX FONCTION INVALID");
    exit(EXIT_FAILURE);
  }    
}

double g(double dx, double dy, int i, int j, int choix){
  double x = i*dx;
  double y = j*dy;
  if(choix == 0 || choix == 2)
    return 0;
  else if(choix == 1)
    return sin(x)+cos(y);
  else{
    fprintf(stderr, "CHOIX FONCTION INVALID");
    exit(EXIT_FAILURE);
  }      
}

double h(double dx, double dy, int i, int j, int choix){
  double x = i*dx;
  double y = j*dy;
  if(choix == 0)
    return 0;
  else if(choix == 1)
    return sin(x)+cos(y);
  else if(choix == 2)
    return 1;
  else{
    fprintf(stderr, "CHOIX FONCTION INVALID");
    exit(EXIT_FAILURE);
  }      
}

void indice(int k, int nx, int ny, int& i, int& j){
  if(k > nx*ny){
    fprintf(stderr, "INDICE TROP GRAND : %d", k);
    exit(EXIT_FAILURE);
  }
  i = k/nx;
  j = k%nx;
}

void charge(int me, int ny, int np, int &i0, int &i1){
  int q = ny/np;
  int r = ny%np;
  if(r == 0){
    i0 = me*q;
    i1 = (me+1)*q-1;
  }
  else{
    if(me<r){
      i0 = me*q+me;
      i1 = q*(me+1)+me;
    }
    else{
      i0 = me*q+r;
      i1 = i0 + q-1;
    }
  }
}


