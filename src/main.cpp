#include <mpi.h>
#include "parallel_conduction.hpp"
#include "main.hpp"
#include "util.hpp"

#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char** argv){
  if(argc < 2){
    fprintf(stderr, "ERROR COMMAND LINE\n");
    fprintf(stderr, "USAGE : %s config_file\n", argv[0]);
    exit(0);
  }
  MPI_Init(NULL,NULL);
  int me, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* output file */
  int output = open("output", O_CREAT | O_WRONLY | O_TRUNC,0744);
  if (!output){
    perror("open : fichier output\n");
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  /* Parsing config file */
  config_t c;
  parse_file(argv[1],c);
  c.np = size;
  c.dx = c.Lx/(c.Nx-1);
  c.dy = c.Ly/(c.Ny-1);
  c.dt = 10;
  c.output = output;
  int i0,i1,kmax,m,n,i,Nyloc,convloc, conv;
  double tolerance(1.E-6);
  kmax =10000;
  charge(me,i0,i1,c);
  //cout << "Je suis me "<< me << " et i0 et i1 valent " << i0 << " " << i1 << endl;
  indice(i1, m, n, c);
  //cout << m << "   " << n << endl;
  Nyloc=(i1-i0+1)/c.Nx;
  conv =1;
  //

  Eigen::VectorXd u = Eigen::VectorXd::Zero(c.Nx*Nyloc);
  Eigen::VectorXd utemp = Eigen::VectorXd::Zero(c.Nx*Nyloc);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(c.Nx*Nyloc);
  Eigen::VectorXd x0 = Eigen::VectorXd::Zero(c.Nx*Nyloc);
  Eigen::MatrixXd Aloc = Eigen::MatrixXd::Zero(c.Nx*Nyloc,c.Nx*Nyloc);

  Remplissage(Aloc,c.Nx,Nyloc,c);
  //cout << "Je suis me : " << me << " et voilà Aloc : "<<Aloc << endl;
  //second_membre(me,u,c);
  i=0;
  while((i<10000) && (conv!=0))
  {
    w =second_membre(me, u, c);
    //cout<<"Je suis me "<<me << " voilà x0 : " << x0 << endl;
    //cout<< "Voila la conv : " << convloc;
    utemp = u;
    //Gradientconjugue(Aloc,u,w,x0 ,tolerance, kmax, c, me);
    BIGradientconjugue(Aloc,u,w,x0 ,tolerance, kmax, c, me);

    //cout << "Voila utemp :" <<utemp;
    //cout << "Je suis me 2 : " << me << " et voila u : " << u << endl;
    convloc = Convergence(utemp, u, tolerance);
    //cout<< "Voila la convloc : " << convloc;
    MPI_Allreduce(&convloc, &conv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    i++;
    // cout<<i<<endl;
    //fflush(stdout);
  }
  cout<<"Je suis sorti avec tant d'itérations : " <<i<<endl;

  log_result(c.output,u);

  if (c.choix == 0)
  {
    for(int i=0; i<i1-i0+1; i++)
    {
      indice(i0+i, m, n, c);
      //cout<< "test " << i0  << " " << i1 << " " << i << " " << " " << m << " "  << n <<  endl;
      //cout << 2*( m*c.dx*(1-m*c.dx)*n*c.dy*(1-n*c.dy)) << endl;
    }
  }
  MPI_Finalize();

  return 0;
}
