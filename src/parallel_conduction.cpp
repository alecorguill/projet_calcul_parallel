/* parallel_conduction.cpp */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen>

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

Eigen::VectorXd second_membre(int me, int nx, int ny, int np, double dx, double dy, double D, double Lx, double Ly, double dt, int choix)
{
  MPI_Status Status;

  int i0, i1, i, j, nyloc;
  double alpha, beta_y, beta_x;
  alpha = 1/dt + 2D/(dx*dx) + 2D/(dy*dy);
  beta_x = -D/(dx*dx);
  beta_y = -D/(dy*dy);

  charge(me, ny, np, i0, i1);
  nyloc = i1-i0+1;

  Eigen::VectorXd floc, g1, g2, env2;

  floc.resize(nyloc*nx);
  floc = 0.;
  g1.resize(nx);
  g2.resize(nx);

  if(me!=0)
  {
    Eigen::VectorXd env1;
    env1.resize(nx);
    for(int i=0; i<nx; i++)
    {
      env1(i) =
    }

    MPI_Send()
  }
  if(me!=np-1)
  {
    Eigen::VectorXd env2;
    env2.resize(nx);
  }

  // Envoie des informations des processeurs précédants le bord haut (env1) et le bas (env2)



  // Appel des vecteurs pour le bord haut (g1) et bas (g2) de la matrice initiale

  for(int i=0; i<nx; i++)
  {
    if (me == 0)
    {
      g1(i) = D*g(dx, dy, i, 0, choix)/(dy*dy);
      g2(i) = ;
    }
    else if (me == np-1)
    {
      g2(i) = D*g(dx, dy, i, np-1, choix)/(dy*dy);
      g1(i) = ;
    }
    else
    {
      g1(i) = ;
      g2(i) = ;
    }


  }



  for (int j=0; j<nyloc; j++) // Numérotation ligne par ligne donc j en première boucle
  {
    for (int i=0; i<nx; i++)
    {
      // Rajoute des termes pour la matrice initiale gauche et droite
      if ((i==0) || (i==nx-1))
      {
        floc(i+j*nx) = f(dx, dy, Lx, Ly, i, j, dt, choix) + D*h(dx, dy, i, j, choix)/(dx*dx);
      }
      else
      {
      floc(i+j*nx) = f(dx, dy, Lx, Ly, i, j, dt, choix);
      }
      // Rajoute des termes pour la matrice initiale haut et bas
      if (j==0)
      {
        floc(i+j*nx) =+ ;
      }
      else if (j==nyloc-1)
      {
        floc(i+j*nx) =+ ;
      }
    }
  }


  f.resize()
}
