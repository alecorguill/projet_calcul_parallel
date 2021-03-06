MODULE fonction
  IMPLICIT NONE
  INCLUDE 'mpif.h'

CONTAINS

  ! On donne un k et il sort les i, j correspondantes 
  SUBROUTINE Indice(k,nx,ny,i,j)
    INTEGER, INTENT(in)  :: k, nx, ny 
    INTEGER, INTENT(out) :: i, j
    INTEGER              :: k1 

    IF (k<=nx*ny) THEN 
       k1 = k 
       DO WHILE (k1>nx)
          k1 = k1-nx
       END DO
       i = k1 
       j = (k-i)/nx + 1 
    ELSE 
       PRINT*, "Il est pas possible de faire l'indicage" 
    END IF
  END SUBROUTINE Indice


  SUBROUTINE charge(me, N, Np, i1, iN)
    INTEGER, INTENT(in) :: me, N, Np
    INTEGER, INTENT(out) :: i1, iN
    INTEGER :: l,r

    r=MOD(N,Np)
    l= N/Np ! division entière

    IF (r==0) THEN
       i1 = me*l +1
       iN = (me+1) *l

    ELSE
       IF (me <r) THEN
          i1 = me*(l+1) + 1
          iN = (me+1) *(l+1)
       ELSE
          i1 = r*(l+1) + (me-r)*l + 1
          iN = i1 + l - 1
       END IF

    END IF
  END SUBROUTINE charge


  !Subroutine pour le calcul du vecteur F en fonction des conditions au bord 
  SUBROUTINE fonct(nx, ny, Lx, Ly, dt, D, i1, in, F, choix)
    INTEGER,                INTENT(in) :: nx, ny, i1, in, choix
    REAL,                   INTENT(in) :: D, dt, Lx, Ly 
    REAL, DIMENSION(i1:in), INTENT(out):: F
    INTEGER                            :: i,j,k
    REAL                               :: dx,dy
    REAL, DIMENSION(i1:in)             :: G,H

    F=0.; G=0.; H=0.
    dx=1./(nx+1)
    dy=1./(ny+1)
    F = 0.; G = 0.; H = 0.
    DO k=i1,in
       CALL Indice(k,nx,ny,i,j) ! on determine les indices i,j en fct de k

       ! Choix si on est dans le cas stationnaire
       If (choix == 1) Then 
          F(k)= 2*(j*dy-j*dy*j*dy+i*dx-i*dx*i*dx)

       ! Choix si on est dans le cas stationnaire périodique
       Else if (choix == 2) Then 
          F(k) = sin(i*dx) + cos(j*dy) 
          IF(j==1) THEN
             G(k)= sin(i*dx) + cos((j-1)*dy)
          ELSE IF (j==ny) THEN
             G(k)= sin(i*dx) + cos((j+1)*dy)
          END IF
          IF (i==1) THEN
             H(k)= sin((i-1)*dx) + cos(j*dy)
          ELSE IF (i==nx) THEN
             H(k)= sin((i+1)*dx) + cos(j*dy)
          END IF

       ! Choix si on est dans le cas instationnaire
       Else if (choix == 3) Then 
          F(k) = exp(-(i*dx-0.5*Lx)*(i*dx-0.5*Lx))*exp(-(j*dy-0.5*Ly)*(j*dy-0.5*Ly))*cos(2.*atan(1.)*dt)
          IF (i==1) THEN
             H(k)= 1.
          ELSE IF (i==nx) THEN
             H(k)= 1.
          END IF

       ! Choix par défaut, le choix numéro 1 dans notre cas
       Else 
          Print*, "Ce choix n'est pas possible, le choix 1 est attribué par défaut"
          F(k)= 2*(j*dy-j*dy*j*dy+i*dx-i*dx*i*dx)
       End If
    END DO
    F(i1:in) = F(i1:in) + D*H(i1:in)/(dx*dx) + D*G(i1:in)/(dy*dy)
  END SUBROUTINE fonct


  ! Calcule le produit matrice vecteur x avec la matrice creuse et sort le vecteur solution, calcule sol = A*F 
  SUBROUTINE Prod(F, sol, nx, ny, Lx, Ly, dt, D, i1, in)
    INTEGER,                      INTENT(in)  :: nx, ny, i1, in
    REAL, DIMENSION(i1-ny:in+ny), INTENT(in)  :: F
    REAL,                         INTENT(in)  :: Lx, Ly, D, dt 
    REAL, DIMENSION(i1:in),       INTENT(out) :: sol
    INTEGER                                   :: N, i
    REAL                                      :: a, bx, by, dx, dy

    ! Initialisation du vecteur solution et définition de N 
    sol = 0.
    N=nx*ny

    ! Initialisation de a, bx, by, dx et dy  
     dx = Lx / (nx+1)
     dy = Ly / (ny+1)
     a = 1/dt + (2*D)/(dx*dx) + (2*D)/(dy*dy)
     bx = -D/(dx*dx)
     by = -D/(dy*dy)
  

    DO i=i1,iN

       ! Gère la première ligne de chaque bloc
       IF (MOD(i-1,nx)==0) THEN ! Il n'y a pas de i-1 -> conditions aux limites a gauche
          IF (i < ny+1) THEN
             sol(i) = a*F(i) + by*F(i+1) + bx*F(i+ny)
          ELSE IF (i>N-ny) THEN
             sol(i) = a*F(i) + by*F(i+1) + bx*F(i-ny) 
          ELSE 
             sol(i) = a*F(i) + by*F(i+1) + bx*F(i-ny) + bx*F(i+ny)
          END IF

       ! Gère la dernière ligne de chaque bloc 
       ELSE IF (MOD(i,nx)==0) THEN ! Il n'y a pas de i+1 -> conditions aux limites a droite
          IF (i<ny+1) THEN
             sol(i) = a*F(i) + by*F(i-1) + bx*F(i+ny) 
          ELSE IF (i>N-ny) THEN 
             sol(i) = a*F(i) + by*F(i-1) + bx*F(i-ny) 
          ELSE 
             sol(i) = a*F(i) + by*F(i-1) + bx*F(i-ny) + bx*F(i+ny)
          END IF
          
       ! Gère le milieu du bloc    
       ELSE
          IF (i<ny+1) THEN
             sol(i)= a*F(i) + by*F(i-1) + by*F(i+1) + bx*F(i+ny)
          ELSE IF (i>N-ny) THEN
             sol(i)= a*F(i) + by*F(i-1) + by*F(i+1) + bx*F(i-ny) 
          ELSE 
             sol(i) = a*F(i) + by*F(i-1) + by*F(i+1) + bx*F(i-ny) + bx*F(i+ny) 
          END IF
       END IF
    END DO
  END SUBROUTINE Prod

END MODULE fonction
