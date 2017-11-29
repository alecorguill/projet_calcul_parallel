/* parallel_conduction.cpp */

#include "util.hpp"

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen>
#include <mpi.h>

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

  int i0, i1, i, j, Nyloc, tag(10);
  double alpha, beta_y, beta_x;
  alpha = 1/c.dt + 2*c.D/(c.dx*c.dx) + 2*c.D/(c.dy*c.dy);
  beta_x = -c.D/(c.dx*c.dx);
  beta_y = -c.D/(c.dy*c.dy);

  charge(me, i0, i1, c);
  Nyloc = i1-i0+1;

  Eigen::VectorXd g1, g2, env2;
  Eigen::VectorXd floc = Eigen::VectorXd::Zero(Nyloc*c.Nx);
  g1.resize(c.Nx);
  g2.resize(c.Nx);

  if(me!=0)
  {
    double env1[c.Nx];
    for(int i=0; i<c.Nx; i++)
    {
      env1[i] = u(i);
    }
    MPI_Send(env1, c.Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
  }

  if(me!=c.np-1)
  {
    double env2[c.Nx];
    for(int i=Nyloc*(c.Nx-2); i<Nyloc*c.Nx-1; i++)
    {
      env2[i] = u(i);
    }
    MPI_Send(env2, c.Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);
  }

  // Envoie des informations des processeurs précédants le bord haut (env1) et le bas (env2)



  // Appel des vecteurs pour le bord haut (g1) et bas (g2) de la matrice initiale


  if (me == 0)
  {
    double rev1[c.Nx];
    MPI_Recv(rev1, c.Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &Status);
    for(int i=0; i<c.Nx; i++)
    {
      g1(i) = c.D*g(i, 0, c)/(c.dy*c.dy);
      g2(i) = c.D*rev1[i]/(c.dy*c.dy);
    }
  }
  else if (me == c.np-1)
  {
    double rev2[c.Nx];
    MPI_Recv(rev2, c.Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &Status);
    for(int i=0; i<c.Nx; i++)
    {
      g2(i) = c.D*g(i, c.np-1, c)/(c.dy*c.dy);
      g1(i) = c.D*rev2[i]/(c.dy*c.dy);
    }
  }
  else
  {
    double rev1[c.Nx], rev2[c.Nx];
    MPI_Recv(rev1, c.Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &Status);
    MPI_Recv(rev2, c.Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &Status);
    for(int i=0; i<c.Nx; i++)
    {
      g1(i) = c.D*rev2[i]/(c.dy*c.dy);
      g2(i) = c.D*rev1[i]/(c.dy*c.dy);
    }


  }



  for (int j=0; j<Nyloc; j++) // Numérotation ligne par ligne donc j en première boucle
  {
    for (int i=0; i<c.Nx; i++)
    {
      // Rajoute des termes pour la matrice initiale gauche et droite
      if ((i==0) || (i==c.Nx-1))
      {
        floc(i+j*c.Nx) = f(i,j,c) + c.D*h(i, j, c)/(c.dx*c.dx);
      }
      else
      {
      floc(i+j*c.Nx) = f(i,j,c);
      }
      // Rajoute des termes pour la matrice initiale haut et bas
      if (j==0)
      {
        floc(i+j*c.Nx) =+ g1(i);
      }
      else if (j==Nyloc-1)
      {
        floc(i+j*c.Nx) =+ g2(i);
      }
    }
  }
}
