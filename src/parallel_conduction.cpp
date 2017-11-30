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

  int i0, i1, Nyloc, tag(10);

  charge(me, i0, i1, c);
  Nyloc = i1-i0+1;
  std::cout<<"Voila la valeur de Nyloc " <<Nyloc<<"pour me = "<< me <<std::endl;

  Eigen::VectorXd g1, g2;
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
    MPI_Send(env1, c.Nx, MPI_DOUBLE, me-1, tag+2*me, MPI_COMM_WORLD);
  }

  if(me!=c.np-1)
  {
    double env2[c.Nx];
    for(int i=0; i<c.Nx; i++)
    {
      env2[i] = u(i+c.Nx*(Nyloc-1));
    }
    MPI_Send(env2, c.Nx, MPI_DOUBLE, me+1, tag+2*me+1, MPI_COMM_WORLD);
  }

  // Envoie des informations des processeurs précédants le bord haut (env1) et le bas (env2)



  // Appel des vecteurs pour le bord haut (g1) et bas (g2) de la matrice initiale

  if (me==0)
  {
    double rev1[c.Nx];
    MPI_Recv(rev1, c.Nx, MPI_DOUBLE, me+1, tag+2*(me+1), MPI_COMM_WORLD, &Status);
    for(int i=0; i<c.Nx; i++)
    {
      g1(i) = c.D*g(i, 0, c)/(c.dy*c.dy);
      g2(i) = c.D*rev1[i]/(c.dy*c.dy);
    }
  }
  else if (me == c.np-1)
  {
    double rev2[c.Nx];
    MPI_Recv(rev2, c.Nx, MPI_DOUBLE, me-1, tag+2*(me-1)+1, MPI_COMM_WORLD, &Status);
    for(int i=0; i<c.Nx; i++)
    {
      g2(i) = c.D*g(i, c.np-1, c)/(c.dy*c.dy);
      g1(i) = c.D*rev2[i]/(c.dy*c.dy);
    }
  }
  else
  {
    double rev1[c.Nx], rev2[c.Nx];
    MPI_Recv(rev1, c.Nx, MPI_DOUBLE, me+1, tag+2*(me+1), MPI_COMM_WORLD, &Status);
    MPI_Recv(rev2, c.Nx, MPI_DOUBLE, me-1, tag+2*(me-1)+1, MPI_COMM_WORLD, &Status);
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

  return floc;
}

void Gradientconjugue(Eigen::MatrixXd A, Eigen::VectorXd u, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax)
{
  int N(0);
  N=u.size();
  Eigen::VectorXd Z;
  Eigen::VectorXd rk;
  Eigen::VectorXd p;
  Eigen::VectorXd rk1; //rk1=rk+1

  int k(0);
  double a(0);
  double g(0);

  rk.resize(N);
  rk1.resize(N);
  p.resize(N);
  rk=b-A*x0;
  p=rk;

  Z.resize(N);

  while( (rk.squaredNorm() > tolerance) && (k <= kmax ))
  {
    Z = A*p;
    a = (rk.dot(rk))/(Z.dot(p));
    u = u + a*p;
    rk1 = rk -a*Z;
    g =(rk1.dot(rk1))/(rk.dot(rk));
    p = rk1 +g*p;
    rk=rk1;
    k++;
  }
  //cout << _x << endl;
  // il faut retourner x c'est donc pas un Void !!!!!!!!!!!!!!!!
  if(k > kmax)
  {
    std::cout << "tolérance non atteinte " << std::endl;
  }
}

void Remplissage(Eigen::MatrixXd& A, int Nx, int Ny, config_t& c)
{
  int N;
  N = Nx*Ny;
  for(int i=0; i<N; i++)
  {
      A(i,i) = 1/c.dt + 2*c.D/(c.dx*c.dx) + 2*c.D/(c.dy*c.dy);
      if((i%Ny != 0) && ( i!=N-1))
      {
        A(i+1,i) = -c.D/(c.dy*c.dy);
        A(i,i+1) = -c.D/(c.dy*c.dy);
      }
      if (i>Ny-1)
      {
        A(i,i-Ny) = -c.D/(c.dx*c.dx);
        A(i-Ny,i) = -c.D/(c.dx*c.dx);
      }
  }
}
