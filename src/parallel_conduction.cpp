/* parallel_conduction.cpp */

#include "util.hpp"

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include "float.h"
#include <fcntl.h>


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
  int k1;
  k1=k;
  i=0;
  j=0;
  while(k1>=c.Nx)
  {
    k1 = k1-c.Nx;
  }
  i = k1;
  j = (k-i)/c.Nx;
}

void charge(int me, int &i0, int &i1, config_t& c){
  int r = c.Ny%c.np;
  int q = (c.Ny-r)/c.np;

  if(r == 0){
    i0 = (me*q)*c.Nx;
    i1 = ((me+1)*q)*c.Nx-1;
  }
  else{
    if(me<r){
      i0 = (me*q+me)*c.Nx;
      i1 = (q+1)*(me+1)*c.Nx-1;
    }
    else{
      i0 = (me*q+r)*c.Nx;
      i1 = (me*q+r)*c.Nx + q*c.Nx-1;
    }
  }
}


Eigen::VectorXd second_membre(int me, Eigen::VectorXd u, config_t& c)
{
  MPI_Status Status;

  int i0, i1, Nyloc, tag(10), i, j;

  charge(me, i0, i1, c);

  Nyloc = (i1-i0+1)/c.Nx;

  Eigen::VectorXd g1, g2;
  Eigen::VectorXd floc = Eigen::VectorXd::Zero(Nyloc*c.Nx);


  if(c.np!=1)
  {
    g1.resize(c.Nx);
    g2.resize(c.Nx);

    double env1[c.Nx];
    double env2[c.Nx];

    if(me!=0)
    {
      for(int i=0; i<c.Nx; i++)
      {
        env1[i] = u(i);//0.2*u(i) + 0.5*(-u(i)+u(i+c.Nx))/c.dy;//u(i);//0.5*u(i) + 0.5*(u(i)-u(i+c.Nx))/c.dy;
        // alpha u(i) + beta u(i+c.Nx)
      }
      MPI_Send(env1, c.Nx, MPI_DOUBLE, me-1, tag+2*me, MPI_COMM_WORLD);
    }

    if(me!=c.np-1)
    {
      for(int i=0; i<c.Nx; i++)
      {
        env2[i] = u(i+c.Nx*(Nyloc-1));//0.2*u(i+c.Nx*(Nyloc-1)) + 0.5*(-u(i+c.Nx*(Nyloc-1))+ u(i+c.Nx*(Nyloc-1)-c.Nx))/c.dy;
        // alpha u(i+c.Nx*(Nyloc-1)) + beta u(i+c.Nx*(Nyloc-1)-c.Nx)

      }
      MPI_Send(env2, c.Nx, MPI_DOUBLE, me+1, tag+2*me+1, MPI_COMM_WORLD);
    }

    // Envoie des informations des processeurs précédants le bord haut (env1) et le bas (env2)



    // Appel des vecteurs pour le bord haut (g1) et bas (g2) de la matrice initiale
    double rev1[c.Nx];
    double rev2[c.Nx];

    if (me==0)
    {

      MPI_Recv(rev1, c.Nx, MPI_DOUBLE, me+1, tag+2*(me+1), MPI_COMM_WORLD, &Status);
      for(int i=0; i<c.Nx; i++)
      {
        g1(i) = c.D*g(i, 0, c)/(c.dy*c.dy);
        g2(i) = c.D*rev1[i]/(c.dy*c.dy);
      }
    }
    else if (me == c.np-1)
    {
      MPI_Recv(rev2, c.Nx, MPI_DOUBLE, me-1, tag+2*(me-1)+1, MPI_COMM_WORLD, &Status);
      for(int i=0; i<c.Nx; i++)
      {
        g2(i) = c.D*g(i, c.Ny-1, c)/(c.dy*c.dy);
        g1(i) = c.D*rev2[i]/(c.dy*c.dy);
      }
    }
    else
    {
      MPI_Recv(rev1, c.Nx, MPI_DOUBLE, me+1, tag+2*(me+1), MPI_COMM_WORLD, &Status);
      MPI_Recv(rev2, c.Nx, MPI_DOUBLE, me-1, tag+2*(me-1)+1, MPI_COMM_WORLD, &Status);
      for(int i=0; i<c.Nx; i++)
      {
        g1(i) = c.D*rev2[i]/(c.dy*c.dy);
        g2(i) = c.D*rev1[i]/(c.dy*c.dy);
      }
    }
  }

  if(c.np==1)
  {
    g1.resize(c.Nx);
    g2.resize(c.Nx);

    for(int i=0; i<c.Nx; i++)
    {
      g1(i) = c.D*g(i, 0, c)/(c.dy*c.dy);
      g2(i) = c.D*g(i, c.Ny-1, c)/(c.dy*c.dy);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  /* //pour avec 1 proc
  for (int k=i0; k<i1+1; k++)
  {
  indice(k, i, j, c);
  floc(k-i0) = u(k-i0)/c.dt + f(i,j,c) + c.D*h(i,j,c)/(c.dx*c.dx);
}
*/

for (int k=i0; k<i1+1; k++)
{
  i=0;
  j=0;
  indice(k, i, j, c);

  // Rajoute des termes pour la matrice initiale gauche et droite
  if ((i==0) || (i==c.Nx-1))
  {
    floc(k-i0) = u(k-i0)/c.dt + f(i,j,c) + c.D*h(i,j,c)/(c.dx*c.dx);
  }
  else
  {
    floc(k-i0) = u(k-i0)/c.dt + f(i,j,c);
  }
  // Rajoute des termes pour la matrice initiale haut et bas
  if (j==i0/c.Nx)
  {
    floc(k-i0) = floc(k-i0) + g1(i);
  }
  if (j==((i1+1)/c.Nx) - 1)
  {

    floc(k-i0) = floc(k-i0) + g2(i);
  }
}
return floc;
}

int Convergence(Eigen::VectorXd utemp, Eigen::VectorXd u, double e)
{
  int i, N;
  N = u.size();

    if (fabs(u.norm()-utemp.norm())/u.norm()>e)
    {
      return 1;
    }
  return 0;
}

void Gradientconjugue(Eigen::MatrixXd A, Eigen::VectorXd& u, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax, config_t& c, int me)
{
  int _k(0), itemax(250);
  Eigen::VectorXd rk(A.rows()), rk1(A.rows()), z(A.rows()), p(A.rows());
  double alpha(1.0);
  double gamma(1.0);
  rk = b;
  double beta = rk.norm();
  p = rk;
  beta = rk.norm();
  for(int i=0; i<A.rows(); i++)
  {
    rk1(i)=0.;
    u(i)=0.;
  }

  while ((beta > tolerance) && (_k < itemax))
  {
    z = A*p;
    alpha = (rk.dot(rk))/(z.dot(p)+DBL_EPSILON);
    u = u + alpha*p;
    rk1 = rk - alpha*z;
    gamma = (rk1.dot(rk1))/(rk.dot(rk)+DBL_EPSILON);
    p = rk1 + gamma*p;
    beta = rk.norm();
    rk=rk1;
    _k = _k + 1;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

  void BIGradientconjugue(Eigen::MatrixXd A, Eigen::VectorXd& u, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax, config_t& c, int me)
    {
      int _k(0), itemax(150);
      Eigen::MatrixXd K(A.rows(),A.cols());
      Eigen::VectorXd rk(A.rows()), rk1(A.rows()), z(A.rows()), p(A.rows()),s(A.rows()),z3(A.rows()),z2(A.rows()),r0(A.rows()),y(A.rows()),t(A.rows());
      double alpha,w,bet,rho,rho1;
      double beta = rk.norm();
      rk = b;
      r0 = b;

      alpha=1;
      rho=1;
      rho1=1;
      w=1;
      beta = rk.norm();
      for(int i=0; i<A.rows(); i++)
      {
          rk1(i)=0.;
          u(i)=0.;
          p(i)=0.;
          K(i,i)=1/A(i,i);
      }

      while ((beta > tolerance) && (_k < itemax))
      {
        rho1=rk.dot(r0);
        bet=alpha*rho1/(rho*w);
        p=rk+bet*(p-w*z);
        y= K*p;
        z = A*y;
        alpha = (rho1)/(z.dot(b));

        s=rk-alpha*z;
        t=K*s;
        z2= A*t;
        z3= K*z2;
        w=(z3.dot(t)/(z3.dot(z3)));
        u = u + alpha*y+w*t;
        rk1 = s - w*z2;

        beta = rk.norm();
        rk=rk1;
        rho=rho1;
        _k = _k + 1;
      }
      MPI_Barrier(MPI_COMM_WORLD);

    }


// Fonction vérifiée -> Elle fonctionne
void Remplissage(Eigen::MatrixXd& A, int Nx, int Ny, config_t& c)
{
  int N;
  N = Nx*Ny;
  for(int i=0; i<N; i++)
  {
    A(i,i) = 1/c.dt + 2*c.D/(c.dx*c.dx) + 2*c.D/(c.dy*c.dy);
    if(((i+1)%Nx != 0) && (i!=N-1))
    {
      A(i+1,i) = -c.D/(c.dy*c.dy);
      A(i,i+1) = -c.D/(c.dy*c.dy);
    }
    if (i>Nx-1)
    {
      A(i,i-Nx) = -c.D/(c.dx*c.dx);
      A(i-Nx,i) = -c.D/(c.dx*c.dx);
    }
  }
}
