PROGRAM projetpara
  USE conjugate
  USE fonction
  IMPLICIT NONE

  INTEGER,PARAMETER                   :: tag=100
  INTEGER                             :: me, np, statinfo, i1, in, N, nx, ny, i, j, k, choix 
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
  REAL                                :: D, Lx, Ly, dx, dy, dt
  REAL*8                              :: t1,t2
  REAL, DIMENSION(:), ALLOCATABLE     :: sol
  character(20)                       :: name

  CALL MPI_INIT(statinfo)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,statinfo)

  ! Définition des paramètres par l'utilisateur 
  open(unit = 9, file = 'data.txt')
  read(9,*) nx, ny
  read(9,*) Lx, Ly
  read(9,*) dt
  read(9,*) D
  read(9,*) choix
  close(9)

  ! Paramètres de discrétisation 
  N  = nx*ny
  dx = Lx / (nx+1)
  dy = Ly / (ny+1)

  ! Allocation du vecteur solution à la taille de i1, in 
  CALL charge(me, N, Np, i1, in)
  allocate (sol(i1:iN))
  sol = 0.
  choix = 1


  ! Calcul du temps pour le gradient conjugué pour chaque processeur
  t1 = MPI_WTIME()
  CALL GC(sol, nx, ny, Lx, Ly, dt, D, i1, in, me, np, status, statinfo, choix)
  t2 = MPI_WTIME()
  print*,"Pour ",Np,"processeurs,le temps de calcul est le suivant",t2-t1,"s"


  ! Ecriture dans un fichier solution
  call rename(me,name)
  open(5,file = name)

  IF (choix == 2) Then 
     Do  k=i1, iN
        Call Indice(k, nx, ny, i, j)
        write(5,*) sol(k), sin(i*dx) + cos(j*dy), sol(k) - sin(i*dx) - cos(j*dy)
     end do
  Else If (choix == 3) Then 
     Do  k=i1, iN
        Call Indice(k, nx, ny, i, j)
        write(5,*) sol(k)
     end do
  Else
     do k=i1, iN
        Call Indice(k, nx, ny, i, j)
        write(5,*) sol(k), i*dx*(1-i*dx)*j*dy*(1-j*dy), sol(k) -  i*dx*(1-i*dx)*j*dy*(1-j*dy)
     end do
  End if

  close(5)
  deallocate(sol)
  CALL MPI_FINALIZE(statinfo)

CONTAINS

  ! Subroutine pour nommer le fichier en fonction du processeur
  SUBROUTINE Rename(Me,name)
    IMPLICIT NONE
    INTEGER      :: Me
    CHARACTER*13 ::name
    CHARACTER*3  :: tn
    INTEGER      :: i1,i2,i3
    i1 = Me/100
    i2 =( Me - 100*i1)/10
    i3 = Me - 100*i1 -10*i2
    tn = CHAR(i1+48)//CHAR(i2+48)//CHAR(i3+48)
    name='sol'//tn//'.dat'
  END SUBROUTINE Rename

END PROGRAM projetpara
