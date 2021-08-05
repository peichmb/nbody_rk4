MODULE NBINTEGRATORS
  
CONTAINS

    SUBROUTINE rk4nb(n,m,r,v,dt,eps)
        INTEGER :: n
        DOUBLE PRECISION :: dt,eps
        DOUBLE PRECISION, DIMENSION(1:n) :: m
        DOUBLE PRECISION, DIMENSION(1:n,1:3) :: r,v,a
        DOUBLE PRECISION, DIMENSION(1:n,1:3) :: kr1,kr2,kr3,kr4
        DOUBLE PRECISION, DIMENSION(1:n,1:3) :: kv1,kv2,kv3,kv4
        
        CALL comp_acc(n,m,r,a,eps)
        kr1=dt*v
        kv1=dt*a
        
        CALL comp_acc(n,m,r+0.5d0*kr1,a,eps)
        kr2=dt*(v+0.5d0*kv1)
        kv2=dt*a
        
        CALL comp_acc(n,m,r+0.5d0*kr2,a,eps)
        kr3=dt*(v+0.5d0*kv2)
        kv3=dt*a
        
        CALL comp_acc(n,m,r+kr3,a,eps)
        kr4=dt*(v+kv3)
        kv4=dt*a
        
        r=r+(kr1+2.d0*(kr2+kr3)+kr4)/6.d0
        v=v+(kv1+2.d0*(kv2+kv3)+kv4)/6.d0               
    END SUBROUTINE rk4nb
    
    SUBROUTINE comp_acc(n,m,r,a,eps)
        INTEGER :: n,i,j
        DOUBLE PRECISION, DIMENSION(1:n) :: m
        DOUBLE PRECISION, DIMENSION(1:n,1:3) :: r,a
        DOUBLE PRECISION, DIMENSION(1:3) :: rij_
        DOUBLE PRECISION :: rij2,eps
        a=0.d0
        DO i=1,n
            DO j=i+1,n
                rij_=r(j,:)-r(i,:)
                rij2=DOT_PRODUCT(rij_,rij_)
                a(i,:)=a(i,:)+m(j)/( (rij2+eps**2)**1.5d0 )*rij_
                a(j,:)=a(j,:)-m(i)/( (rij2+eps**2)**1.5d0 )*rij_
            END DO
        END DO
    END SUBROUTINE comp_acc
    
END MODULE NBINTEGRATORS
