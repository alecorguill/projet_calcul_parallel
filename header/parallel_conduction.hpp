#ifndef PARALLEL_CONDUCTION
#define PARALLEL_CONDUCTION


/* flux */
double f(double dx, double dy, double Lx, double Ly, int i, int j, double dt, int choix);

/* conditions de bord */
double g(double dx, double dy, int i, int j, int choix);
double h(double dx, double dy, int i, int j, int choix);

/* traduction des indices matricielles en indice linéaires */
void indice(int k, int nx, int ny, int& i, int& j);

/* renvoie le morceau de domaine en fonction du proc  */
void charge(int me, int ny, int np, int &i0, int &i1);

#endif
