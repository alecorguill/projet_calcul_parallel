#include "Solveurs.h"

using namespace std;
using namespace Eigen;


//constructeur
Solveur::Solveur(Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0,double tolerance, int kmax)
:_A(A),_x(x),_b(b),_x0(x0),_tolerance(tolerance),_kmax(kmax)
{

}




//****************************************************************************************
//                                       Solveurs
//****************************************************************************************

//Gradientoptimal
void Solveur::Gradientoptimal(Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax)
{
  int N(0);
  N=_x.size();
  Eigen::VectorXd Z;
  Eigen::VectorXd rk;

  int k(0);
  double a(0);
  double m(0);
  //cout << b << endl;
  rk.resize(N);
  rk=_b-_A*_x0;
  //cout << rk << endl;
  m=rk.squaredNorm();
  //cout << m << endl;
  Z.resize(N);

  //cout << _tolerance << endl;
  //cout << k << endl;
  //cout << rk.squaredNorm() << endl;
  while( (m > _tolerance) && (k <= _kmax ) )
  {
    Z = _A*rk;
    //  cout << Z << endl;
    a = rk.dot(rk)/ (Z.dot(rk)); //normalement le produit scalaire, à tester

    //cout << a << endl;
    _x = _x + a*rk;
    //cout << _x << endl;
    rk = rk - a*Z;
    m=rk.squaredNorm() ;
    k++;
    //cout << k << endl;
  }
  //cout << _x << endl;
  // il faut retourner x c'est donc pas un Void !!!!!!!!!!!!!!!!
  if(k > _kmax)
  {
    cout << "tolérance non atteinte : " << rk.squaredNorm() << endl;
  }

}


//RésiduMinimum
void Solveur::Residumini(Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax)
{
  int N(0);
  N=_x.size();
  Eigen::VectorXd dk;
  Eigen::VectorXd rk;
  Eigen::VectorXd Z;
  Eigen::MatrixXd At; //At sera la transposée de A
  int k(0);
  double a(0);

  rk.resize(N);
  rk=_b-_A*_x0;

  dk.resize(N);

  Z.resize(N);
  At.resize(N,N);

  At=_A.transpose();

  dk=At*rk;

  while( (rk.squaredNorm() > _tolerance) && (k <= _kmax ))
  {
    Z=_A*dk;
    a= (rk.dot(Z))/(Z.dot(Z));
    _x = _x + a*dk;
    rk=rk -a*Z; //PENSER A SAUVEGARDER LE RESIDU POUR LE TRACER PLUS TARD!!!
    dk=At*rk;
    k++;

  }
  //cout << _x << endl;
  // il faut retourner x c'est donc pas un Void !!!!!!!!!!!!!!!!
  if(k > _kmax)
  {
    cout << "tolérance non atteinte " << endl;
  }

}




//GradientConjugue
void Solveur::Gradientconjugue(Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax)
{
  int N(0);
  N=_x.size();
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
  rk=_b-_A*_x0;
  p=rk;

  Z.resize(N);

  while( (rk.squaredNorm() > _tolerance) && (k <= _kmax ))
  {
    Z = _A*p;
    a = (rk.dot(rk))/(Z.dot(p));
    _x = _x + a*p;
    rk1 = rk -a*Z;
    g =(rk1.dot(rk1))/(rk.dot(rk));
    p = rk1 +g*p;
    rk=rk1;
    k++;
  }
  //cout << _x << endl;
  // il faut retourner x c'est donc pas un Void !!!!!!!!!!!!!!!!
  if(k > _kmax)
  {
    cout << "tolérance non atteinte " << endl;
  }

}
// résidu minimum préconditionné à gauche
void Solveur::Residuminiprecond(Eigen::MatrixXd M, Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax)
{
  //cout << M << endl;
  int N(0);
  int k(0);
  double a(0);
  N=_x.size();
  Eigen::VectorXd Z;

  Eigen::VectorXd q;
  Eigen::VectorXd r;
  Eigen::VectorXd w;



  r.resize(N);
  q.resize(N);
  w.resize(N);
  Z.resize(N);

  _x=_x0;
  r=_b-_A*_x0;

  Solveur resol(_M, q, r,_x0,_tolerance,_kmax);



  //cout << r << endl;

  //Résolution de Mq=r
  resol.Gradientconjugue( _M, q, r, _x0 , _tolerance, _kmax);
  //cout << resol.Getx() << endl;
  q=resol.Getx();

  while( (r.squaredNorm() > _tolerance) && (k<=_kmax) )
  {
    //cout << r << endl;
    w=_A*q;
    //cout << w << endl;
    Solveur resol2(_M,Z,w,_x0,_tolerance,_kmax);
    //Résolution de MZ=w
    resol2.Gradientconjugue( _M, Z, w, _x0 , _tolerance, _kmax);
    Z=resol2.Getx();
  //  cout << q << endl;
    //cout << Z << endl;
    a = (q.dot(Z))/(Z.dot(Z));
    //cout << a << endl;
    _x += a*q;
    //cout << "vite" << endl;
    //cout << _x << endl;
    r -= a*w;
    //cout << r.squaredNorm() << endl;
    q = q - a*Z;
    //cout << a*Z << endl;
    //cout << q << endl;
    k++;
  }
    if(k > _kmax)
  {
    cout << "tolérance non atteinte "<< r.squaredNorm() << endl;
  }
}

//****************************************************************************************
//                                       Préconditionneur
//****************************************************************************************
//préconditionneur diagonal
void Solveur::preconddiag(Eigen::MatrixXd A,Eigen::MatrixXd M)
{
  int N(0);
  N=_A.rows();
  _M.resize(N,N);
  _M.setIdentity(N,N);
  for (int i(0); i<N; i++)
  {
      _M(i,i)=_A(i,i);
  }
}
//préconditionneur Gauss-Seidel
void Solveur::precondgauss(Eigen::MatrixXd A,Eigen::MatrixXd M)
{
  int N(0);
  Eigen::MatrixXd D;
  Eigen::MatrixXd D1;
  Eigen::MatrixXd E;
  Eigen::MatrixXd F;
  N=_A.rows();
  _M.resize(N,N);
  D.resize(N,N);
  D1.resize(N,N);
  E.resize(N,N);
  F.resize(N,N);
  _M.setIdentity(N,N);
  for (int i(0); i<N; i++)
  {
      D(i,i)=A(i,i);
      D1(i,i)=1./(A(i,i));
      for (int j(0); j<N; j++)
      {
          if(j>i)
          {
            F(i,j)=-A(i,j);
          }
          else
          {
          E(i,j)=-A(i,j);
          }
      }

  }
  M=(D-E)*D1*(D-F);
}
//préconditionneur Cholesky incomplet


//permet de renvoyer _x
const VectorXd & Solveur::Getx() const
{
  return _x;
}
const MatrixXd & Solveur::GetM() const
{
  return _M;
}
