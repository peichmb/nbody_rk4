PROGRAM I2S
  IMPLICIT NONE
   
  INTEGER :: a,n
  CHARACTER(LEN=20) :: str
  READ *, a
  n=4
    
  CALL int2str(a,str,n)
    
  PRINT *, str
    
CONTAINS

  SUBROUTINE int2str(a,str,n)
    
    INTEGER :: a,n,temp,i
    CHARACTER(LEN=n+5) :: str
    
    DO i=1,n
      temp=a-a/10*10
      print *,a,temp
      str(n-i+1:n-i+1)=CHAR(temp+48)
      a=a/10
    END DO
    
  END SUBROUTINE int2str

END PROGRAM I2S

