PROGRAM finalProjDriver
  USE finalProjMod
  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: s, lVector, uVector, dVector, bVector, xInitial
  REAL(KIND=8) :: gamma, w
  INTEGER :: n, nMinus2, i
  INTEGER, DIMENSION(:), ALLOCATABLE :: iterations
  CHARACTER (LEN = 13) :: filename
  CHARACTER (LEN = 14) :: outputFile
  
  !OUTPUTS THE PROGRAM HEADER
  CALL title()
  
  filename = "FILTER_IN.DAT"
  outputFile = "FILTER_OUT.DAT"
  
  !READS IN THE FILTER DATA
  CALL readInFile(filename, n, s)
  
  !LETS THE USER KNOW THAT THE FILE WAS SUCCESSFULLY LOADED
  WRITE(*,*) "The file ", filename, " was successfully loaded with ", n, " data values"
  WRITE(*,*) ""
  
  !ALLOCATES THE MEMORY FOR ALL OF THE VECTORS THAT COMPRISE THE TRIDIAGONAL SYSTEM
  ALLOCATE(lVector(n-2), uVector(n-2), dVector(n-2), bVector(n-2), xInitial(n-2))
  
  !INITIALIZES ALL VECTORS TO 0.0
  lVector = 0.0D0
  uVector = 0.0D0
  dVector = 0.0D0
  bVector = 0.0D0
  xInitial = 0.0D0
  
  !PROMPTS THE USER FOR THE VALUE OF GAMMA
  WRITE(*,*) "Please input a value for gamma (-1, 1): "
  READ(*,*) gamma
  
  !OBTAINS THE RIGHT HAND SIDE OF THE EQUATION (IE E*sUnfiltered + b)
  CALL getRHS(gamma, n, s, bVector)
  
  !SETS UP THE TRIDIAGONAL MATRIX WITH THEIR APPROPRIATE VALUES
  nMinus2 = n-2
  CALL initializeTridiagMatrix(uVector, dVector, lVector, gamma, nMinus2)
  
  !w = 2.0D0 / (1.0D0 + SQRT(1.0D0 - (gamma / 2.0D0)**2))
  
  ALLOCATE(iterations(99))
  
  OPEN(UNIT=4, FILE="iterations1.dat", STATUS="REPLACE")
  
  DO i=1, 99
    w = 1.0D0 + i*0.01D0
    WRITE(*,*) "W= ", w
    CALL succOverRelaxation(lVector, dVector, uVector, bVector, xInitial, w, nMinus2, iterations(i))
    xInitial = 0.0D0
    
    write(*,*) iterations(i)
    
  END DO
  
  
END PROGRAM finalProjDriver