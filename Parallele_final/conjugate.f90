MODULE conjugate
  USE fonction
  IMPLICIT NONE

CONTAINS

  SUBROUTINE GC(U0, nx, ny, Lx, Ly, dt, D, i1, in, me, np, status, statinfo, choix)

    INTEGER,PARAMETER                              :: tag=100
    INTEGER,                            INTENT(in) :: nx, ny, i1, in, me, np, statinfo, choix 
    REAL,                               INTENT(in) :: D, dt, Lx, Ly 
    REAL, DIMENSION(i1:iN),             INTENT(OUT):: U0
    REAL, DIMENSION(i1:iN)                         :: U, diff, pr, rk, rk1, F
    REAL, DIMENSION(i1-ny:in+ny)                   :: x
    INTEGER, DIMENSION(MPI_STATUS_SIZE),INTENT(in) :: status
    INTEGER                                        :: k, i
    REAL                                           :: Alpha, Beta, Norme, alphah, alphab, betah, s1
    REAL, DIMENSION(ny)                            :: p1, p2


    ! On determine la fonction qui sert à exprimer les conditions du terme source et des conditions aux limites
    CALL Fonct(nx, ny, Lx, Ly, dt, D, i1, in, F, choix)  

    ! Initialisation des différents paramètres
    x(i1-ny:in+ny) = 0.; Norme = 1.; 
    rk(i1:iN)= F(i1:iN); U0(i1:in)=0.
    k = 0
    x(i1:in) = F(i1:in)

    ! Boucle du gradient conjugué, avec condition de sortie sur la norme ou le nombre d'itérations  
    DO WHILE (k<100000 .AND. Norme > 0.000001)
       x(i1:in)=x(i1:in)+U0(i1:in)/dt

       ! Envoie le haut du vecteur x au processeur me-1
       IF (me /= 0) THEN       
          CALL MPI_SEND(x(i1:i1+ny-1), ny, MPI_FLOAT, me-1, tag+2*me, MPI_COMM_WORLD, statinfo)
       ENDIF

       ! Envoie le bas du vecteur x au processeur me+1
       IF(me /= np-1)THEN
          CALL MPI_SEND(x(in-ny+1:in), ny, MPI_FLOAT, me+1, tag+2*me+1, MPI_COMM_WORLD, statinfo)
       ENDIF

       ! Reception par me du message de me-1
       IF (me/=0)THEN
          p1=0.
          CALL MPI_RECV(p1(1:ny), ny, MPI_FLOAT, me-1, tag+2*(me-1)+1, MPI_COMM_WORLD, status,statinfo)
          x(i1-ny:i1-1)=p1(1:ny)
       ENDIF

       ! Reception par me du message de me+1
       IF (me/=np-1)THEN
          p2=0.
          CALL MPI_RECV(p2(1:ny), ny, MPI_FLOAT, me+1, tag+2*(me+1), MPI_COMM_WORLD, status,statinfo)
          x(in+1:in+ny)=p2(1:ny)
       ENDIF

       ! Fin des communications et début d'une itération du gradient conjugué en parallèle

       ! Calcul du produit matrice vecteur A*x, le résultat sort avec pr 
       CALL Prod(x, pr, nx, ny, Lx, Ly, dt, D, i1, in)

       ! Numérateur du produit alpha en local, noté alphah qui est ensuite envoyé à tous les processeurs
       alphah=0.
       s1=0. 
       Do i = i1, in
          s1 =s1 + rk(i)*rk(i)
       End Do
       CALL mpi_allreduce(s1,alphah,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD,statinfo)


       ! Calcul du produit scalaire entre pr et x entre i1 et in 
       ! On utilise pas DOT_PRODUCT qui ne marche pas pour des vecteurs de taille différentes
       s1=0.
       DO i=i1,in
          s1=s1+pr(i)*x(i)
       ENDDO

       ! Dénominateur du produit alpha en local,noté alphab qui est ensuite envoyé à tous les processeurs
       alphab=0.
       CALL mpi_allreduce(s1,alphab,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD,statinfo)

       ! Calcul du alpha global pour tous les processeurs
       Alpha = alphah/alphab

       U(i1:in) = U0(i1:in) + Alpha*x(i1:in)
       rk1(i1:in) = rk(i1:in) - Alpha*pr(i1:in)

       ! Numérateur du produit beta en local, noté betah qui est ensuite envoyé à tous les processeurs
       s1 =0. 
       Do i =i1, in 
          s1 = s1 + rk1(i)*rk1(i)
       End Do
       CALL mpi_allreduce(s1,betah,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD,statinfo)

       ! Calcul du beta global pour tous les processeurs
       Beta = betah/alphah
       x(i1:iN) = rk1(i1:in) + Beta*x(i1:iN)
       rk(i1:in) = rk1(i1:in)
       k = k+1

       ! Pour le calcul de la norme locale
       DO i=i1,iN
          diff(i)=U(i)-U0(i)
       ENDDO

       Norme=0.
       s1=0.
       Do i = i1, in 
          s1 = s1 + diff(i)*diff(i)
       End Do

       ! Calcul de la norme globale  
       CALL mpi_allreduce(s1, Norme,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD,statinfo)

       U0(i1:in)=U(i1:in)

    END DO
  END SUBROUTINE GC

END MODULE conjugate
