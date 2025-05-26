!STEVE WINWARD

MODULE finalProjMod
  CONTAINS
  
  !RETURNS THE APPROXIMATED CONDITION NUMBER FOR THE PADE FILTER MATRIX
  FUNCTION conditionNum(gamma)
    IMPLICIT NONE
    REAL(KIND=8) :: gamma, conditionNum
    
    conditionNum = (1 + ABS(gamma)) / (1 - ABS(gamma))
    
    RETURN
  END FUNCTION conditionNum
  
  !CALCULATES THE RIGHT HAND SIDE OF THE PADE FILTER EQUATION
  SUBROUTINE getRHS(gamma, n, sUnfiltered, rhs)
    REAL(KIND=8), INTENT(INOUT), DIMENSION(:) :: sUnfiltered, rhs
    REAL(KIND=8), INTENT(INOUT) :: gamma
    REAL(KIND=8) :: a, aOver2
    INTEGER, INTENT(INOUT) :: n
    
    a = (1.0D0+gamma) / 2.0D0
    aOver2 = a / 2.0D0
    
    rhs(1) = a*sUnfiltered(2) + aOver2*sUnfiltered(3) + (aOver2 - gamma/2.0D0) * sUnfiltered(1)
    
    DO i=2, n-3
      rhs(i) = aOver2*sUnfiltered(i) + a*sUnfiltered(i+1) + aOver2*sUnfiltered(i+2)
    END DO
    
    rhs(i) = aOver2*sUnfiltered(i) + a*sUnfiltered(i+1) + (aOver2 - gamma/2.0D0)*sUnfiltered(i+2)
    
  END SUBROUTINE getRHS
  
  !INITIALIZES ALL THREE VECTORS USED IN THE PADE FILTER TRIDIAGONAL MATRIX
  SUBROUTINE initializeTridiagMatrix(uVector, dVector, lVector, gamma, n)
    REAL(KIND=8), INTENT(INOUT), DIMENSION(:) :: uVector, dVector, lVector
    REAL(KIND=8), INTENT(IN) :: gamma
    INTEGER, INTENT(INOUT) :: n
    
    uVector(1) = gamma / 2.0D0
    dVector(1) = 1.0D0
    lVector(1) = 0.0 !SINCE THERE IS NO L VECTOR ELEMENT IN THE FIRST ROW
    DO i=2, n-1
      lVector(i) = gamma / 2.0D0
      uVector(i) = lVector(i)
      dVector(i) = 1.0D0
    END DO
    
    lVector(i) = gamma / 2.0D0
    dVector(i) = 1.0D0
    uVector(i) = 0.0 !SINCE THERE IS NO U-VECTOR ELEMENT IN THE LAST ROW
    
  END SUBROUTINE initializeTridiagMatrix
  
  !PERFORMS JACOBI ITERATION FOR A TRIDIAGONAL MATRIX
  SUBROUTINE jacobiIteration(lVector, dVector, uVector, bVector, xInitial, n, tolerance, maxIterations)
    REAL(KIND=8), INTENT(INOUT), DIMENSION(:) :: lVector, dVector, uVector, bVector, xInitial
    REAL(KIND=8), INTENT(INOUT) :: tolerance
    REAL(KIND=8), DIMENSION(n) :: xNew, residual, tempVector
    REAL(KIND=8) :: residualNorm, bVectorNorm, toleranceProduct
    INTEGER, INTENT(INOUT) :: n, maxIterations
    INTEGER :: i, k
    
    !INITIALIZE ALL ARRAYS TO 0
    residual = 0.0D0
    xNew = 0.0D0
    
    bVectorNorm = normInfinity(bVector, n)
    toleranceProduct = tolerance*bVectorNorm
    
    !UPDATES THE RESIDUAL MATRIX
    CALL triDiagMatrixTimesVector(lVector, dVector, uVector, xInitial, tempVector, n)
    residual = bVector - tempVector
    
    DO k=1, maxIterations
      !UPDATES NORM OF VECTORS
      residualNorm = normInfinity(residual, n)
      
      !CHECKS IF THE RELATIVE ERROR IS LESS THAN THE TOLERANCE
      IF(residualNorm < toleranceProduct) EXIT
      
      xNew(1) = bVector(1) - uVector(1)*xInitial(2)
      residual(1) = xNew(1) - xInitial(1)
      
      DO i=2, n-1
        xNew(i) = bVector(i) - lVector(i)*xInitial(i-1) - uVector(i)*xInitial(i+1)
        residual(i) = xNew(i) - xInitial(i)
      END DO
      
      !NOTE THAT I=N AFTER THE DO LOOP
      xNew(i) = bVector(i) - lVector(i)*xInitial(i-1)
      residual(i) = xNew(i) - xInitial(i)
      
      !UPDATES THE INITIAL GUESS
      xInitial = xNew
      
    END DO
    
    !OVERWRITES THE MAX # OF ITERATIONS TO BE THE NUMBER OF ITERATIONS TO CONVERGENCE
    maxIterations = k-1
    
    !TESTING CASE
    !WRITE(*,*) "The jacobi iteration took ", (k-1), " iterations to get within a tolerance of ", tolerance
    
  END SUBROUTINE jacobiIteration
  
  !RETURNS THE INFINITY NORM FOR ANY COLUMN VECTOR (IE THE MAX ELEMENT)
  FUNCTION normInfinity(x, n)
    IMPLICIT NONE
    REAL(KIND=8) :: normInfinity
    REAL(KIND=8), DIMENSION(:) :: X
    INTEGER :: n, i
    
    normInfinity = ABS(x(1))
    
    DO i=2, n
      IF(ABS(x(i)) > normInfinity) THEN
        normInfinity = ABS(x(i))
      END IF
    END DO
    
    RETURN
  END FUNCTION normInfinity
  
  !READS IN A FILE WITH THE SPECIFIED "FILENAME" AND STORES ALL THE DATA VALUES TO THE
  !ARRAY S.
  SUBROUTINE readInFile(filename, n, s)
    REAL(KIND=8), INTENT(INOUT), DIMENSION(:), ALLOCATABLE :: s
    INTEGER, INTENT(INOUT) :: n
    CHARACTER (LEN=*) :: filename
    
    OPEN(3, FILE=filename, STATUS = 'UNKNOWN')
    READ(3, *) n !NOTE THERE ARE N+1 ELEMENTS IN S
    
    n = n+1
    
    !ALLOWS THE NUMBER OF DATA POINTS TO BE DYNAMIC
    ALLOCATE(s(n))
    
    DO i=1, n
      READ(3, *) s(i)
      !TESTING CASE
      !WRITE(*,*) s(i)
    END DO
    
  END SUBROUTINE readInFile
  
  !SUCCESSIVE OVER RELAXATION
  SUBROUTINE succOverRelaxation(lVector, dVector, uVector, bVector, xInitial, w, n, iterations)
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: lVector, dVector, uVector, bVector, xInitial
    REAL(KIND=8), INTENT(INOUT) :: w
    REAL(KIND=8), DIMENSION(n) :: xNew, residual, tempVector
    REAL(KIND=8) :: residualNorm, toleranceProduct, bVectorNorm, tolerance
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(INOUT) :: iterations
    INTEGER :: i, j, maxIterations
    
    maxIterations = 1000
    
    tolerance = SQRT(EPSILON(tolerance))
    
    bVectorNorm = normInfinity(bVector,n)
    toleranceProduct = tolerance*bVectorNorm
    
    DO i=1, maxIterations
      CALL triDiagMatrixTimesVector(lVector, dVector, uVector, xInitial, tempVector, n)
      
      residual = bVector - tempVector
      
      residualNorm = normInfinity(residual, n)
      
      IF(residualNorm < toleranceProduct) EXIT
      
      !FIRST ELEMENT IN THE NEW VECTOR
      xNew(1) = w*(bVector(1)-uVector(1)*xInitial(2))
      xNew(1) = (1.0D0 - w)*xInitial(1) + xNew(1)
      
      DO j=2, n-1
        xNew(j) = w*(bVector(j)-lVector(j)*xNew(j-1)-uVector(j)*xInitial(j+1))
        xNew(j) = (1.0D0 - w)*xInitial(j) + xNew(j)
      END DO
      
      !LAST ELEMENT IN THE NEW VECTOR
      xNew(j) = w*(bVector(j)-lVector(j)*xNew(j-1))
      xNew(j) = (1.0D0 - w)*xInitial(j) + xNew(j)
      
      !UPDATES THE INITIAL GUESS
      xInitial = xNew
    END DO
    
    iterations = i
    
  END SUBROUTINE succOverRelaxation
  
  !PRINTS PROGRAM HEADER
  SUBROUTINE title()
    WRITE(*,*)"_________________________________________________________"
    WRITE(*,*)"" 
    WRITE(*,*)"This program was written and compiled by Steve Winward"
    WRITE(*,*)"James Madison University, Mathematics Department"
    WRITE(*,*)"_________________________________________________________"
    WRITE(*,*)"" 
  END SUBROUTINE title
  
  !PERFORMS GAUSSIAN ELIMIANTION AND BACK SUBSITUTION ON A TRIDIAGONAL MATRIX
  !THIS SUBROUTINE STORES THE FINAL VALUES INTO THE B-VECTOR
  SUBROUTINE triDiagGE(lVector, dVector, uVector, bVector, n)
    REAL(KIND=8), INTENT(INOUT), DIMENSION(:) :: lVector, dVector, uVector, bVector
    REAL(KIND=8) :: temp
    INTEGER, INTENT(IN) :: n
    INTEGER :: i
    
    !GAUSSIAN ELIMINATION FOR A TRIDIAGONAL MATRIX
    DO i=2, n
      temp = lVector(i) / dVector(i-1)
      dVector(i) = dVector(i) - temp*uVector(i-1)
      bVector(i) = bVector(i) - temp*bVector(i-1)
    END DO
    
    
    !BACK SUBSTITUTION, REPLACING THE B-VECTOR WITH THE FILTERED VALUES
    bVector(n) = bVector(n) / dVector(n)
    
    DO i=n-1, 1, -1
      bVector(i) = bVector(i) - uVector(i)*bVector(i+1)
      bVector(i) = bVector(i) / dVector(i)
    END DO
    
  END SUBROUTINE triDiagGE
  
  !CALCULATES THE RESULTING COLUMN VECTOR FOR A TRIDIAGONAL MATRIX MULTIPLIED BY A COLUMN VECTOR
  SUBROUTINE triDiagMatrixTimesVector(lVector, dVector, uVector, xVector, resultantVector, n)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(INOUT), DIMENSION(:) :: lVector, dVector, uVector, xVector, resultantVector
    INTEGER, INTENT(IN) :: n
    INTEGER :: i
    
    resultantVector(1) = dVector(1)*xVector(1) + uVector(1)*xVector(2)
    
    DO i=2, n-1
      resultantVector(i) = lVector(i)*xVector(i-1)+dVector(i)*xVector(i)+uVector(i)*xVector(i+1)
    END DO
    
    resultantVector(i) = lVector(i)*xVector(i-1)+dVector(i)*xVector(i)
    
  END SUBROUTINE triDiagMatrixTimesVector
  
  !WRITES THE MATRIX CONTAINING THE JACOBI ITERATION INFORMATION (IE CN #, ITERATIONS, GAMMA)
  SUBROUTINE writeJacobiInfo(matrix, n, filename)
    REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: matrix
    CHARACTER(LEN=*), INTENT(INOUT) :: filename
    INTEGER, INTENT(IN) :: n
    
    OPEN(UNIT=4, FILE=filename, STATUS="REPLACE")
    
    WRITE(4,*) "GAMMA        ITERATIONS    CONDITIION#"
    
    DO i=1, n
      WRITE(4, '(1X, F5.2, F17.0, F27.2)' ) matrix(i, 1), matrix(i,2), matrix(i, 3)
    END DO
    
  END SUBROUTINE writeJacobiInfo
  
  !WRITES THE MATRIX TO THE SPECIFIED FILENAME
  SUBROUTINE writeToFile(matrix, r, c, filename)
    CHARACTER(LEN=*), INTENT(INOUT) :: filename
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: matrix
    INTEGER, INTENT(IN) :: r, c
    INTEGER :: i, j    
    
    OPEN(UNIT=4, FILE=filename, STATUS="REPLACE")
    
    DO i=1, r
      DO j=1, c
        WRITE(4, *) matrix(i, j)
      END DO
      WRITE(*,*) " "
    END DO
    
  END SUBROUTINE writeToFile
  
END MODULE finalProjMod