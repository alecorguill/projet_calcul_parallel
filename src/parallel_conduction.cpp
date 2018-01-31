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


// Fonction f pouvant être modifiée selon 3 choix possibles
double f(int i, int j, config_t& c){
  double x = i*c.dx;
  double y = j*c.dy;
  // Cas numéro 1
  if(c.choix == 0)
  return 2*(x-x*x+y-y*y);
  // Cas numéro 2
  else if(c.choix == 1)
  return sin(x)+cos(y);
  // Cas numéro 3
  else if(c.choix == 2)
  return exp(-pow((x-c.Lx*0.5),2))*exp(-pow((y-c.Ly*0.5),2))*cos((M_PI*0.5)*c.dt);
  else{
    fprintf(stderr, "CHOIX FONCTION INVALID");
    exit(EXIT_FAILURE);
  }
}

// Fonction g qui permet de gérer les conditions limites du domaine
double g(int i, int j, config_t& c){
  double x = i*c.dx;
  double y = j*c.dy;
  // Cas numéro 1 et 3
  if(c.choix == 0 || c.choix == 2)
  return 0;
  // Cas numéro 2
  else if(c.choix == 1)
  return sin(x)+cos(y);
  else{
    fprintf(stderr, "CHOIX FONCTION INVALID");
    exit(EXIT_FAILURE);
  }
}

// Fonction g qui permet de gérer les conditions limites du domaine
double h(int i, int j, config_t& c){
  double x = i*c.dx;
  double y = j*c.dy;
  // Cas numéro 1
  if(c.choix == 0)
  return 0;
  // Cas numéro 2
  else if(c.choix == 1)
  return sin(x)+cos(y);
  // Cas numéro 3
  else if(c.choix == 2)
  return 1;
  else{
    fprintf(stderr, "CHOIX FONCTION INVALID");
    exit(EXIT_FAILURE);
  }
}

// Procédure permettant de rentrer un indice k qui est sous forme vectoriel et l'écrire sous forme matriciel
// On donne en entrée k et on resort en sortie i et j
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

// Procédure permettant en connaissant un processeur me de connaître les indices de début et de fin de ses coefficients
// Cette procédure prend donc comme paramètre me et donne l'indice de début i0 et celui de fin i1
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

// Fonction permettant la modification du second membre pour envoyer les informations de recouvrement
// Elle prend en paramètre le processeur et le vecteur u et renvoie un vecteur local du second membre
Eigen::VectorXd second_membre(int me, Eigen::VectorXd u, config_t& c)
{
  MPI_Status Status;

  int i0, i1, Nyloc, tag(10), i, j;

  charge(me, i0, i1, c);

  Nyloc = (i1-i0+1)/c.Nx;

  Eigen::VectorXd g1, g2;
  Eigen::VectorXd floc = Eigen::VectorXd::Zero(Nyloc*c.Nx);

  // Nécessaire si on veut faire fonctionner le programme en séquentiel
  if(c.np!=1)
  {
    g1.resize(c.Nx);
    g2.resize(c.Nx);

    double env1[c.Nx];
    double env2[c.Nx];
    // Envoie du recouvrement du haut de chaque processeur
    if(me!=0)
    {
      for(int i=0; i<c.Nx; i++)
      {
        // Pour la premiere gestion du recouvrement
        env1[i] = u(i);
        //Pour la seconde gestion du recouvrement
        //env1[i]=0.2*u(i) + 0.5*(-u(i)+u(i+c.Nx))/c.dy;//u(i);//0.5*u(i) + 0.5*(u(i)-u(i+c.Nx))/c.dy;
        // alpha u(i) + beta u(i+c.Nx)
      }
      MPI_Send(env1, c.Nx, MPI_DOUBLE, me-1, tag+2*me, MPI_COMM_WORLD);
    }
    // Envoie du recouvrement du bas de chaque processeur
    if(me!=c.np-1)
    {
      for(int i=0; i<c.Nx; i++)
      {
        //Pour la première gestion du recouvrement
        env2[i] = u(i+c.Nx*(Nyloc-1));
        //Pour la seconde gestion du recouvrement
        //envi2[i]=0.2*u(i+c.Nx*(Nyloc-1)) + 0.5*(-u(i+c.Nx*(Nyloc-1))+ u(i+c.Nx*(Nyloc-1)-c.Nx))/c.dy;
        // alpha u(i+c.Nx*(Nyloc-1)) + beta u(i+c.Nx*(Nyloc-1)-c.Nx)

      }
      MPI_Send(env2, c.Nx, MPI_DOUBLE, me+1, tag+2*me+1, MPI_COMM_WORLD);
    }

    // Envoie des informations des processeurs précédants le bord haut (env1) et le bas (env2)



    // Appel des vecteurs pour le bord haut (g1) et bas (g2) de la matrice initiale
    double rev1[c.Nx];
    double rev2[c.Nx];
    // Reception du vecteur gérant le recouvrement pour me = 0
    if (me==0)
    {

      MPI_Recv(rev1, c.Nx, MPI_DOUBLE, me+1, tag+2*(me+1), MPI_COMM_WORLD, &Status);
      for(int i=0; i<c.Nx; i++)
      {
        g1(i) = c.D*g(i, 0, c)/(c.dy*c.dy);
        g2(i) = c.D*rev1[i]/(c.dy*c.dy);
      }
    }
    // Reception du vecteur gérant le recouvrement pour me = np-1
    else if (me == c.np-1)
    {
      MPI_Recv(rev2, c.Nx, MPI_DOUBLE, me-1, tag+2*(me-1)+1, MPI_COMM_WORLD, &Status);
      for(int i=0; i<c.Nx; i++)
      {
        g2(i) = c.D*g(i, c.Ny-1, c)/(c.dy*c.dy);
        g1(i) = c.D*rev2[i]/(c.dy*c.dy);
      }
    }
    // Reception du vecteur gérant le recouvrement pour me différent de 0 et np-1
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
  // Conditions pour faire fonctionner le code en séquentiel
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

// Fonction permettant de vérifier la Convergence pour sortir des itérations de Schwarz
// Elle prend en paramètre 2 vecteurs et une tolérance et renvoie un nombre
// Si le processeur a une bonne convergence on renvoie 1 sinon c'est 0
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

// Procédure pour faire le gradient conjugué
// Prend en paramètre une matrice A, un vecteur u et un second membre b : on résoud Au = b
// Le vecteur x0 est nécessaire pour le gradient conjugué, dans notre cas il sera prit égale à 0
// Tolerance et kmax sont la pour gérer la sortie du gradient conjugué, on tend vers une tolérance et un garde fou kmax qui est un nombre d'itérations maximas
// La résolution se fait en locale, nécessitant ainsi la présence du me
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

  // Fonction Bi-Gradient Conjugué permettant de faire comme le gradient conjugué pour des matrices non symétriques
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


// Procédure permettant de remplir la matrice A connaissant son Nx et son Ny
// Prend en paramètre Nx et Ny et renvoie A 
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
