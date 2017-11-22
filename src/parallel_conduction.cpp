/* parallel_conduction.cpp */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen>

double f(int i, int j, config_t& c){
  double x = i*c.dx;
  double y = j*c.dy;
  if(c.choix == 0)
    return 2*(x-x*x+y-y*y);
  else if(c.choix == 1)
    return sin(x)+cos(y);
  else if(c.choix == 2)
    return exp(-pow((x-c.Lx*0.5),2))*exp(-pow((y-c.Ly*0.5),2))*cos((M_PI*0.5)*c.dt);
  else{
    fprintf(stderr, "CHOIX FONCTION INVALID");
    exit(EXIT_FAILURE);
  }
}

double g(int i, int j, config_t& c){
  double x = i*c.dx;
  double y = j*c.dy;
  if(c.choix == 0 || c.choix == 2)
    return 0;
  else if(c.choix == 1)
    return sin(x)+cos(y);
  else{
    fprintf(stderr, "CHOIX FONCTION INVALID");
    exit(EXIT_FAILURE);
  }
}

double h(int i, int j, config_t& c){
  double x = i*c.dx;
  double y = j*c.dy;
  if(c.choix == 0)
    return 0;
  else if(c.choix == 1)
    return sin(x)+cos(y);
  else if(c.choix == 2)
    return 1;
  else{
    fprintf(stderr, "CHOIX FONCTION INVALID");
    exit(EXIT_FAILURE);
  }
}

void indice(int k, int& i, int& j, config_t& c){
  if(k > c.Nx*c.Ny){
    fprintf(stderr, "INDICE TROP GRAND : %d", k);
    exit(EXIT_FAILURE);
  }
  i = k/c.Nx;
  j = k%c.Nx;
}

void charge(int me, int &i0, int &i1, config_t& c){
  int q = c.Ny/c.np;
  int r = c.Ny%c.np;
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

Eigen::VectorXd second_membre(int me, Eigen::VectorXd u, config_t& c)
{
  MPI_Status Status;

  int i0, i1, i, j, Nyloc;
  double alpha, beta_y, beta_x;
  alpha = 1/c.dt + 2*c.D/(c.dx*c.dx) + 2*c.D/(c.dy*c.dy);
  beta_x = -c.D/(c.dx*c.dx);
  beta_y = -c.D/(c.dy*c.dy);

  charge(me, i0, i1, c);
  Nyloc = i1-i0+1;

  Eigen::VectorXd floc, g1, g2, env2;

  floc.resize(Nyloc*c.Nx);
  floc = 0.;
  g1.resize(c.Nx);
  g2.resize(c.Nx);

  if(me!=0)
  {
    Eigen::VectorXd env1;
    env1.resize(c.Nx);
    for(int i=0; i<c.Nx; i++)
    {
      env1(i) = ;
    }

    MPI_Send()
  }
  if(me!=c.np-1)
  {
    Eigen::VectorXd env2;
    env2.resize(c.Nx);
  }

  // Envoie des informations des processeurs précédants le bord haut (env1) et le bas (env2)



  // Appel des vecteurs pour le bord haut (g1) et bas (g2) de la matrice initiale

  for(int i=0; i<c.Nx; i++)
  {
    if (me == 0)
    {
      g1(i) = c.D*g(c.dx, c.dy, i, 0, c.choix)/(c.dy*c.dy);
      g2(i) = ;
    }
    else if (me == c.np-1)
    {
      g2(i) = c.D*g(c.dx, c.dy, i, c.np-1, c.choix)/(c.dy*c.dy);
      g1(i) = ;
    }
    else
    {
      g1(i) = ;
      g2(i) = ;
    }


  }



  for (int j=0; j<Nyloc; j++) // Numérotation ligne par ligne donc j en première boucle
  {
    for (int i=0; i<c.Nx; i++)
    {
      // Rajoute des termes pour la matrice initiale gauche et droite
      if ((i==0) || (i==c.Nx-1))
      {
        floc(i+j*c.Nx) = f(c.dx, c.dy, c.Lx, c.Ly, i, j, c.dt, c.choix) + c.D*h(c.dx, c.dy, i, j, c.choix)/(c.dx*c.dx);
      }
      else
      {
      floc(i+j*c.Nx) = f(c.dx, c.dy, c.Lx, c.Ly, i, j, c.dt, c.choix);
      }
      // Rajoute des termes pour la matrice initiale haut et bas
      if (j==0)
      {
        floc(i+j*c.Nx) =+ ;
      }
      else if (j==Nyloc-1)
      {
        floc(i+j*c.Nx) =+ ;
      }
    }
  }
}
