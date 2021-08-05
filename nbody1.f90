! Código N-BODY con algoritmo Runge-Kutta 4 con softening
! Versión modular, con barra de progreso.

PROGRAM RKUTTANB

    USE STR2NUM
    USE NBINTEGRATORS
    IMPLICIT NONE
    
    INTEGER :: n ! Número de partículas
    DOUBLE PRECISION :: dt, dt_pt, dt_dia, t_run ! Paso, paso para imprimir, paso diagnóst., tiempo de simulación.
    INTEGER :: i,j,k,n_steps,n_pt,cnt_pt,n_dia,n_progbar,cnt_dia,cnt_shot,cnt_progbar
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: m ! Masas de todas las partículas
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: r, v, a ! Posiciones, velocidades y aceleraciones
    DOUBLE PRECISION :: e_kin0, e_pot0, e_kin1, e_pot1, e0, e1 ! Energías cinética y potencial antes y después
    DOUBLE PRECISION :: t, rij2, eps
    CHARACTER(LEN=20) :: shot_name,infile,outfile
    
    READ *, dt
    READ *, dt_pt
    READ *, dt_dia
    READ *, t_run
    READ *, eps
    READ *, infile
    READ *, outfile
    
    OPEN(10,FILE=infile)
    READ(10,*) n
    ALLOCATE(m(1:n))
    ALLOCATE(r(1:n,1:3))
    ALLOCATE(v(1:n,1:3))
    ALLOCATE(a(1:n,1:3))
    DO i=1,n
        READ (10,*) m(i), r(i,1), r(i,2), r(i,3), v(i,1), v(i,2), v(i,3)
    END DO
    CLOSE(10)
    
    n_steps=t_run/dt
    n_pt=dt_pt/dt
    n_dia=dt_dia/dt
    n_progbar=n_steps/50
    cnt_pt=0
    cnt_dia=0
    cnt_shot=1
    cnt_progbar=0
    
    ! Energías iniciales
    
    e_kin0=0.d0
    e_pot0=0.d0
    DO j=1,n
        e_kin0=e_kin0+m(j)*DOT_PRODUCT(v(j,:),v(j,:))
        DO i=1,n
            IF (i /= j) THEN
                rij2=DOT_PRODUCT(r(i,:)-r(j,:),r(i,:)-r(j,:))
                e_pot0=e_pot0-m(i)*m(j)/SQRT(rij2+eps**2)
            END IF
        END DO
    END DO
    e0=e_kin0+e_pot0
    
    OPEN(10,FILE='simul/shot0000.dat')
    DO j=1,n
        WRITE(10,*) t,j,r(j,:),v(j,:)
    END DO
    CLOSE(10)
    
    IF(dt_dia > 0.d0) THEN
        OPEN(20,FILE='simul/rk4dia.dat')    
        WRITE(20,*) 'FICHERO DE DIAGNOSTICO'
        WRITE(20,*) '----------------------'
        WRITE(20,*) ''
        WRITE(20,*) 'PARAMETROS INICIALES'
        WRITE(20,*) '--------------------'
        WRITE(20,*) '  dt = ',dt
        WRITE(20,*) '  dt_pt = ',dt_pt
        WRITE(20,*) '  t_run = ',t_run
        WRITE(20,*) '  n = ',n
        WRITE(20,*) '  eps = ',eps
        WRITE(20,*) '  method = rk4'
        WRITE(20,*) ''
        WRITE(20,*) 'ENERGIA INICIAL'
        WRITE(20,*) '---------------'
        WRITE(20,*) 'At time t = 0.d0'
        WRITE(20,*) '  e_pot0 = ',e_pot0,'; e_kin0 = ',e_kin0
        WRITE(20,*) '     e_tot1 = ',e0,'; err = 0'
        WRITE(20,*) ''    
    END IF
    
    PRINT *
    PRINT *,'PARAMETROS INICIALES'
    PRINT *,'--------------------'
    PRINT *,"  dt = ",dt
    PRINT *,"  dt_pt = ",dt_pt
    PRINT *,"  t_run = ",t_run
    PRINT *,"  n = ",n
    PRINT *,"  eps = ",eps
    PRINT *,"  method = rk4"
    
    ! Preparo barra de progreso
    
    PRINT *
    PRINT *, 'Progreso:'
    PRINT *
    PRINT *, '  0%                                              100%'
    PRINT *, '  |--------------------------------------------------|'
    WRITE(*,'(1A4)',ADVANCE='NO') '   |'

    ! Main loop
    
    DO i=1,n_steps
        CALL rk4nb(n,m,r,v,dt,eps) ! Devuelve r(i+1) y v(i+1). No reciclaje aproximado.
        t=t+dt 
        cnt_pt=cnt_pt+1
        cnt_dia=cnt_dia+1
        cnt_progbar=cnt_progbar+1
        
        IF(cnt_pt == n_pt) THEN
            shot_name=''
            CALL int2str(cnt_shot,shot_name,4)
            shot_name='simul/shot'//TRIM(shot_name)//'.dat'
            OPEN(10,FILE=shot_name)
            DO j=1,n
                WRITE(10,*) t,j,r(j,:),v(j,:)
            END DO
            CLOSE(10)
            cnt_pt=0
            cnt_shot=cnt_shot+1
        END IF

        IF(cnt_dia == n_dia .AND. dt_dia > 0.d0) THEN
            e_kin1=0.d0
            e_pot1=0.d0
            DO j=1,n
                e_kin1=e_kin1+m(j)*DOT_PRODUCT(v(j,:),v(j,:))
                DO k=1,n
                    IF (k /= j) THEN
                        rij2=DOT_PRODUCT(r(k,:)-r(j,:),r(k,:)-r(j,:))
                        e_pot1=e_pot1-m(k)*m(j)/SQRT(rij2+eps**2)
                    END IF
                END DO
            END DO
            e1=e_kin1+e_pot1
            WRITE(20,*) 'At time t = ',t
            WRITE(20,*) '  e_pot1 = ',e_pot1,'; e_kin1 = ',e_kin1
            WRITE(20,*) '     e_tot1 = ',e1,'; err = ', ABS(e0-e1)/ABS(e0)
            WRITE(20,*) ''
            cnt_dia=0                
        END IF
        
        IF (cnt_progbar == n_progbar) THEN
            WRITE(*,'(1A1)',ADVANCE='NO') '='
            cnt_progbar=0
        END IF
              
    END DO
    
    WRITE(*,'(1A1)') '|'
    
    e_kin1=0.d0
    e_pot1=0.d0
    DO j=1,n
        e_kin1=e_kin1+m(j)*DOT_PRODUCT(v(j,:),v(j,:))
        DO i=1,n
            IF (i /= j) THEN
                rij2=DOT_PRODUCT(r(i,:)-r(j,:),r(i,:)-r(j,:))
                e_pot1=e_pot1-m(i)*m(j)/SQRT(rij2+eps**2)
            END IF
        END DO
    END DO    
    
    ! Diagnóstico
    
    e1=e_kin1+e_pot1
    
    IF(dt_dia > 0.d0) THEN
        WRITE(20,*)
        WRITE(20,*) 'ENERGIA FINAL'
        WRITE(20,*) '-------------'
        WRITE(20,*) 'At time t = ',t
        WRITE(20,*) '  e_pot1 = ',e_pot1,'; e_kin1 = ',e_kin1
        WRITE(20,*) '     e_tot1 = ',e1,'; err = ', ABS(e0-e1)/ABS(e0)
        CLOSE(20)
    END IF
    
    OPEN(30,FILE=outfile)
    WRITE(30,*) n
    DO j=1,n
        WRITE(30,*) m(j),r(j,:),v(j,:)
    END DO
    CLOSE(30)

    PRINT *
    PRINT *,'DIAGNOSTICO ENERGIAS'
    PRINT *,'--------------------'    
    PRINT *,"  e_pot0 = ",e_pot0,"; e_kin0 = ",e_kin0
    PRINT *,"     e_tot0 = ",e0
    PRINT *,"  e_pot1 = ",e_pot1,"; e_kin1 = ",e_kin1
    PRINT *,"     e_tot1 = ",e1    
    PRINT *,"  err = ", ABS(e0-e1)/ABS(e0)
    PRINT *
    PRINT *,"Pasos simulacion: ", n_steps
    PRINT *,"Ficheros escritos:  ", FLOOR(t_run/dt_pt)+1
    PRINT *

END PROGRAM RKUTTANB
