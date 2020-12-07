!***********************************************************************
! Código que calcula la evolución de trayectorias de electrones en 
! potenciales centrales. Problema 2-dimensional
! Usa las Rutinas RK4 para integrar las ecuaciones diferenciales
! y Ran1 y Gasdev para generar numeros aleatorios (Numerical Recipes)
!
! Sebastián D. López (Version Diciembre 2020).
! Esta versión incorpora potenciales de Yukawa entre los modelos de 
! interacción que tiene el electrón con el núcleo
!************************************************************************

!******************************************************************
!modulo CTMCparfix
!modulo de parametros globales N: numero de coordenadas (q,p)
!nstep: numero de trayectorias con tiempo de ionizacion diferentes
!ntraj: numero de trayectorias con el mismo tiempo de ionizacion
MODULE CTMCparfix
  implicit none
  integer, parameter :: N=4,nstep=30000,ntraj=1
end MODULE
!******************************************************************

!******************************************************************
!modulo models, variables globales que se actualizan en el runtime
!model: 0 for coulombic, 1 for Yukawa
!in Yukawa, V(r)=-beta*Exp(-alfa*r)
MODULE models
  implicit none
  integer model,alfa,beta
end MODULE
!****************************************************************** 


!******************************************************************
! modulo laseratom, variables globales que se actualizan en el runtime
! cantidades en unidades atomicas
! envelope: 'flat', 'sin2' - envolvente del pulso laser
! nc: numero de ciclos del pulso, w0: frecuencia del pulso
! F0: maximo valor del campo electrico
! tau: duracion del pulso
! phi: fase absoluta del pulso
! Zt: carga del core
! IP: energia de ionizacion
MODULE laseratom
  implicit none
  integer nc
  REAL(8) w0,F0,tau,phi,Zt,IP,pi
  character(10) envelope
END MODULE
!******************************************************************

!******************************************************************
!module random para variable dummy de rutinas ran1
!y gasdev (Numerical Recipes, Press)
module random
  implicit none
  integer idum
end module
!******************************************************************

!******************************************************************
!Cuerpo del programa
PROGRAM CTMC_varstep
!modulos a utilizar
  use CTMCparfix
  use laseratom
  use random
  use models
  implicit none
  logical :: stat
  INTEGER :: i,j,k,jtr,t_n,t_f,t_i
  real(8) :: ti,tf
  integer :: ctraj,nbound,nevolp,cont_resc
  real(8), ALLOCATABLE, DIMENSION (:) :: y1,y,dydx,yout 
  real(8) h,hevol,t0,RAN1,bpar,gpar 
  real(8) campo,tp
  real(8) fi3,rmod,En
  real(8) eps,proyAs, ev_lenght

  Allocate (y1(N),y(N),dydx(N),yout(N))
  pi=dacos(-1d0)
  
!************************************************************
! Con el clock del sistema genero una semilla aleatoria 
! entre (0,1000) y la imprimo en la salida estandar para 
! verificar que para dos conjuntos de simulaciones en 
! paralelo las semillas son diferentes.
  call SYSTEM_CLOCK(t_i,t_n,t_f)
  idum=-modulo(t_i,1000)
  print*,idum  

!
! Valores de los parametros del problema
!

  Zt=1d0
  IP=0.5d0
  envelope = 'flat'
  nc=1
  w0=0.05d0
  F0=0.075d0
  phi=0d0*pi/2d0
  tau=nc*2d0*pi/w0

!
! parametros de la evolucion temporal
!
  eps=1d-4 !para formulita Avrines y Percival

!
! modelo de potencial a simular
!
  model=1
  
  if(model==0) then
    alfa=0d0
    beta=Zt
  elseif(model==1) then
    alfa=1d0
    beta=1.90863d0
  end if

!
! archivos de salida
!

  open(100,file='ctmc_Test.dat')
  open(300,file='field_Test.dat')

! ************************************************
! cuerpo del calculo
!
!
! puesta a 0 de algunos contadores
! ctraj: trayectorias aceptadas en el sampleo por importancia
! nbound: trayectorias que finalizan ligadas
!
  ctraj=0
  nbound=0
!
! Bucle para cada una de las trayectorias 
! con diferente tiempo de ionizacion
!
  call cpu_time(ti)
  do i=1,nstep
!
! t0: tiempo de ionizacion
!
    t0=tau*Ran1(idum)
    write(300,555)t0,campo(t0)
!
! Bucle con tiempos de ionizacion iguales
!
    do jtr=1,ntraj
!
! condiciones iniciales y(1..4):x,z,vx,vz en tiempo t0
! stat: variable logica de aceptacion/rechazo
!
      call initial(t0,y,N,stat)
      if(stat) then
        ctraj=ctraj+1
!
! guardo las condiciones iniciales en y1
!
        y1 = y 
        yout = 0.d0       
!
! variables relacionadas con la distancia al nucleo
!
        rmod = dsqrt(y(1)*y(1)+y(2)*y(2))
!        
! parametros temporales iniciales
        tp = t0
        h = 0d0
!
! calculo las derivadas para tener el paso temporal
!
        call derivs(t0,y,dydx) 
!
! fi3 fases acumuladas
!
        fi3=0d0

!
! Hasta encontrar una expresion analitica para la evolucion clasica en un potencial de Yukawa
! evolucionaremos una cantidad de tiempo considerable luego de apagado el pulso hasta obtener
! convergencia en las distribuciones (la evolucion deberia hacerse hasta infinito).
!
        if (model == 0) then
          ev_lenght = 1d0
        else if(model == 1) then
          ev_lenght = 100d0
        end if
!
! bucle de evolucion temporal h y hevol corresponden al paso temporal anterior y posterior,
! que se calcula a partir de las fuerzas (Avrines y Percival). Hay dos para hacer integrales 
! por el metodo del trapecio
!
! comienzo de la evolucion de cada trayectoria
!    
        do while (tp<=ev_lenght*tau)
          hevol=dsqrt((eps**2/((y(3)**2+y(4)**2)/rmod**2)+eps*dsqrt(dydx(3)**2+dydx(4)**2)/2d0/rmod))        
!
! límites del paso temporal. Es importante setear de manera adecuada para manejar la relacion
! entre precision de las trayectorias y el tiempo computacional de simulacion
!
          if(hevol < 1d-20) hevol=1d-20
          if(hevol > 1d-1) hevol=1d-1
          tp = tp + hevol
!
! evolucion con el algoritmo RK4 
!
          call derivs(tp,y,dydx)
          call rk4(y,dydx,N,tp,hevol,yout,derivs)
!
! las variables el t+h salen en yout
!
          y = yout
          rmod = dsqrt(y(1)*y(1)+y(2)*y(2))
          fi3 = fi3-((y(3)*y(3)+y(4)*y(4))/2d0-(beta*dexp(-alfa*rmod)* (alfa + 2d0/rmod)))*(hevol+h)/2d0
          h = hevol
        end do ! fin de la evolucion de una trayectoria
!
! calculo h un paso mas para completar la integracion
!
        hevol=dsqrt((eps**2/((y(3)**2+y(4)**2)/rmod**2)+eps*dsqrt(dydx(3)**2+dydx(4)**2)/2d0/rmod))
        if(hevol < 1d-20) hevol=1d-20
        if(hevol > 1d-1) hevol=1d-1
        tp=tp+hevol
        call derivs(tp,y,dydx)
        call rk4(y,dydx,N,tp,hevol,yout,derivs)
        rmod=dsqrt(y(1)*y(1)+y(2)*y(2))
!        
! fase para SCTS
!
        fi3=fi3-((y(3)*y(3)+y(4)*y(4))/2d0-(beta*dexp(-alfa*rmod)* (alfa + 2d0/rmod)))*(hevol)/2d0  
!
! Calculos para la parte asintotica de la trayectoria (solo caso coulombiano)
! Energía de la trayectoria. Ver referencia: Eq(40) PRA94, 013415
!
        En=(y(3)*y(3)+y(4)*y(4))/2d0-(beta*dexp(-alfa*rmod)* (1d0/rmod))
        if(En<=0)nbound=nbound+1
        if (model==0) then       
!
! parametros para la proyección asintotica sin campo PARA LA FASE (solo caso coulombiano)
!
          bpar=dsqrt(1d0/2d0/En)
          gpar=dsqrt(1d0+(y(1)*y(4)-y(2)*y(3))**2/bpar)
          proyAs = Zt*dsqrt(bpar)*(dlog(gpar)+dasinh((y(1)*y(3)+y(2)*y(4))/gpar/dsqrt(bpar))) ! Eq(40) PRA94, 013415 
        end if
        fi3=fi3-(y1(3)*y1(1)+y1(4)*y1(2)-IP*t0)
!
!Escribiendo en archivos de salida
!
! trajectory output
! ** columnas 1:t0,2..5:condiciones iniciales,6..9:estado en t=ev_lenght*tau,
! **          10:E(ev_lenght*tau),11:Fi3 (ev_lenght*tau) 
!
! si el potencial es coulombiano: 12: asympt Coul phase 
!

        if (model==0) write(100,555)t0,y1,y,en,fi3,proyAs 
        if (model==1) write(100,555)t0,y1,y,en,fi3

!
! alguna info en runtime para saber el estado de ejecucion
!
        if (mod(i,1000) == 0) then
          print*,'varstep_grendel trajs:',ctraj,'ligadas',nbound
          call cpu_time(tf)
          print*,'tiempo transcurrido', tf-ti
          call cpu_time(ti)
        end if
      end if
    end do 
  end do 
  print*,'fin trayectorias. Aceptadas:',ctraj
  print*,'fin trayectorias. Ligadas:',nbound

555   Format(g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,& 
     g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6,g14.6)
END PROGRAM


!******************************************************************
! Subrutina que calcula una distribucion inicial mediante un 
! metodo de aceptacion-rechazo (Similar a Metropolis)
subroutine initial(t,y,n,stat)
  use random
  use laseratom
  implicit none
  integer N,m
  real*8 t,y(N),gasdev,kapp,campo,alpha,IPF,gam,eta
  Real*8 vrmax,beta,rvr,wdistn,wdistd,ran1,FF,randu,preexp
  logical stat

  kapp=dsqrt(2d0*IP)
  FF=campo(t)+1d-15
!
!preexp value is the modulation to distribution see N. Shvetsov-Shilovski article (SFA)
!
  preexp=dsqrt(dexp(-2d0*kapp**3d0/3d0/dabs(FF))/dexp(-2d0*kapp**3d0/3d0/(F0)))
  randu=ran1(idum)

  if(randu <= preexp)then
    stat=.true.
    y(3)=gasdev(idum)*dsqrt(dabs(FF/2d0/kapp)) !(the correct is this, perpendicular velocity)
!   y(3)=gasdev(idum)*dsqrt(dabs(FF/kapp))     !(this is what N. Shvetsov-Shilovski  used in paper)
!
!the simplest tunneling distribution
!
  y(1)=0d0    !x coordinate
  y(2)=-IP/FF !z coordinate
  y(4)=0d0    !parallel velocity
!the perpendicular velocity correspond to the Gaussian given before

  else
   stat=.false.
  end if
end

!*********************************************************************
!---------------------------------------------------------------------
! Subrutina de derivadas 
SUBROUTINE derivs(x,y,dydx)
  use CTMCparfix
  use laseratom
  use models
  implicit none
  real(8) :: x,y(N),dydx(N),campo,Aan,rmod,yuk
  rmod=dsqrt(y(1)*y(1)+y(2)*y(2))
  dydx(1)=y(3)
  dydx(2)=y(4) 

  if (model == 0) then
    dydx(3)=-Zt*y(1)/rmod**3
    dydx(4)=-campo(x)-Zt*y(2)/rmod**3
  elseif (model==1) then
    yuk = -beta*dexp(-alfa*rmod)/rmod**2*(alfa+1d0/rmod)
    dydx(3)=y(1)*yuk
    dydx(4)=-campo(x)-y(2)*yuk  
  end if
  return
END





!******************************************************
!Subroutine that provides electric field
function campo(t)
  use laseratom
  implicit none
  real*8 t,campo
  if((t >= 0).and.(t <= tau/1d0)) then
    if (envelope == 'sin2') then
      campo=F0*dSin((t)*pi/tau)**2*dSin(w0*t+phi) !sin squared envelope
    else if(envelope == 'flat') then
      campo=F0*dsin(w0*t+phi)
    end if
  else
    campo = 0d0
  end if
return
end

