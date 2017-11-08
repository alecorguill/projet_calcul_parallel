#ifndef PARALLEL_CONDUCTION
#define PARALLEL_CONDUCTION

/* flux */
double f(double dx, double dy, double Lx, double Ly, double dt);

/* conditions de bord */
double g(double dx, double dy);
double h(double dx, double dy);

/* traduction des indices matricielles en indice linéaires */
int indice(int nx, int ny, int i, int j);

/* renvoie le morceau de domaine en fonction du proc  */
void charge(int me, int N, int np, int &i1, int &iN);

#endif