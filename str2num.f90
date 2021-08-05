MODULE STR2NUM

CONTAINS

  SUBROUTINE int2str(a,str,n)
    
    INTEGER :: a,n,temp,i,a_copy
    CHARACTER(LEN=n+5) :: str
   
    a_copy=a
    DO i=1,n
      temp=a_copy-a_copy/10*10
      str(n-i+1:n-i+1)=CHAR(temp+48)
      a_copy=a_copy/10
    END DO
    
  END SUBROUTINE int2str
  
  FUNCTION str2int(str,longstr)
  
    CHARACTER(LEN=longstr) :: str
    CHARACTER :: curcar
    INTEGER :: i,longint,longstr,start,str2int

    longint=LEN(TRIM(str))
    IF(str(1:1).EQ.'-') THEN
      start=2
    ELSE
      start=1
    END IF
    
    str2int=0
    
    DO i=start,longflt
      curcar=str(i:i)
      str2int=str2int*10
      str2int=str2int+IACHAR(curcar)-48
    END DO
    
    IF(start.EQ.2) THEN
      str2int=-str2int
    END IF
    
  END FUNCTION str2int

  FUNCTION str2flt(str,longstr)
    IMPLICIT NONE
     
    CHARACTER(LEN=longstr) :: str
    CHARACTER :: curcar
    INTEGER :: i,bruto,pointpos,longflt,longstr,start
    DOUBLE PRECISION :: str2flt
    
    longflt=LEN(TRIM(str))
    pointpos=0
    bruto=0
    IF(str(1:1).EQ.'-') THEN
      start=2
    ELSE
      start=1
    END IF
    
    DO i=start,longflt
      curcar=str(i:i)
      IF (curcar.EQ.'.') THEN
        pointpos=i
        CYCLE
      END IF
      bruto=bruto*10
      bruto=bruto+IACHAR(curcar)-48
    END DO
    IF (pointpos.NE.0) THEN
      str2flt=bruto/(10.d0**(longflt-pointpos))
    ELSE
      str2flt=bruto
    END IF
    IF(start.EQ.2) THEN
      str2flt=-str2flt
    END IF
        
  END FUNCTION str2flt
      
END MODULE STR2NUM
