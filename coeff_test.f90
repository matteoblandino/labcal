      program coeff_test
!
!     Test coeffcients:
!     A = int dg/dn dl
!     B = int g dl
!
!______________________________________________________________________
      implicit none
      integer ::i,n
      real(8) :: theta,dtheta
      real(8) :: R
      real(8) :: pi
      real(8) :: xa,xb,ya,yb
      real(8) :: xs,ys
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: A,B,sumA,sumB
      character(10) :: cmn_string
      character(43) :: cmn_string2
!______________________________________________________________________
      pi = acos(-1.d0)
!     input
      cmn_string = ' input n'
      write(*,1000) cmn_string
      read(*,*) n
      cmn_string = ' input R'
      write(*,1000) cmn_string
      read(*,*) R
      cmn_string = ' input xs'
      write(*,1000) cmn_string
      read(*,*) xs
      cmn_string = ' input ys'
      write(*,1000) cmn_string
      read(*,*) ys
!
      cmn_string = '    n = '
      write(*,1001) cmn_string,n
      cmn_string = '    R = '
      write(*,1002) cmn_string,R
      cmn_string = '   xs = '
      write(*,1002) cmn_string,xs
      cmn_string = '   ys = '
      write(*,1002) cmn_string,ys
!
1000  format(a)
1001  format(a,i4.4)
1002  format(a,g11.4)
!
      dtheta = 2.d0*pi/float(n)
      theta = 0.d0
      xa = R*cos(theta)
      ya = R*sin(theta)
      open(unit=11,file='./data/nodes_out.dat')
      open(unit=12,file='./data/tangents.dat')
      open(unit=13,file='./data/normales.dat')
      cmn_string2 = '  nodes n = '
      write(11,2000) cmn_string2,n
      cmn_string2 = '  i      theta          x            y      '
      write(11,2001) cmn_string2
      cmn_string2 = 'tangents n= '
      write(12,2000) cmn_string2,n
      cmn_string2 = '  i      theta         taux        tauy     '
      write(12,2001) cmn_string2
      cmn_string2 = ' normals n= '
      write(13,2000) cmn_string2,n
      cmn_string2 = '  i      theta          hnx         hny     '
      write(13,2001) cmn_string2
!
      write(11,2002) i,theta,xa,ya
      sumA = 0.d0
      sumB = 0.d0
      do i=1,n
         theta = float(i)*dtheta
         xb = R*cos(theta)
         yb = R*sin(theta)
!
         call panel(xa,ya,xb,yb,hl,taux,tauy,hnx,hny)
         write(11,2002) i,theta,xb,yb
         write(12,2002) i,theta-dtheta*0.5,taux,tauy
         write(13,2002) i,theta-dtheta*0.5,hnx,hny
!
         call coeffA(xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny,A)
         call coeffB(xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny,B)
         sumA = sumA + A
         sumB = sumB + B
!       
         xa = xb
         ya = yb
      end do
!
      cmn_string = ' sumA = '
      write(*,1002) cmn_string,sumA
      cmn_string = ' sumB = '
      write(*,1002) cmn_string,sumB
!
      close(11)
      close(12)
      close(13)
!
2000  format(a,i4.4)
2001  format(a)
2002  format(i4.4,3(1x,g11.4))
!______________________________________________________________________
      end
!______________________________________________________________________
!______________________________________________________________________
      subroutine panel(xa,ya,xb,yb,hl,taux,tauy,hnx,hny)
!______________________________________________________________________
!
!     defines panel geometry:
!     lenght hl, tangent, taux,tauy and normal hnx,hny
!
!______________________________________________________________________
      implicit none
      real(8) :: xa,xb,ya,yb
      real(8) :: hl,taux,tauy,hnx,hny
!______________________________________________________________________
      hl = sqrt((xb-xa)**2 + (yb-ya)**2)
      taux = (xb-xa)/hl
      tauy = (yb-ya)/hl
      hnx  = (yb-ya)/hl
      hny  = (xb-xa)/hl
      hny  = -1.d0*hny
!______________________________________________________________________
      return
      end 
!______________________________________________________________________
!______________________________________________________________________
      subroutine coeffA(xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny,A)
!______________________________________________________________________
!
!     A = int dg/dn dl
!
!______________________________________________________________________
      implicit none
      real(8) :: xa,xb,ya,yb
      real(8) :: xs,ys
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: A
      real(8) :: sn,cs,atg,atan_S
      real(8) :: xm,ym
      real(8) :: taus,hns
      real(8) :: tau,d
      real(8) :: atgA,atgB
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
      xm = (xa+xb)*0.5d0
      ym = (ya+yb)*0.5d0
      hns = (xs-xm)*hnx + (ys-ym)*hny
      taus = (xs-xm)*taux+(ys-ym)*tauy
!
!     B - integration limit
      tau = hl*0.5d0 -taus
      d = sqrt(tau**2 + hns**2)
      sn = hns/d
      cs = tau/d
      atgB = atan_S(sn,cs)
!
!     A - integration limit
      tau = -hl*0.5d0 -taus
      d = sqrt(tau**2 + hns**2)
      sn = hns/d
      cs = tau/d
      atgA = atan_S(sn,cs)
!
      A = -1.d0*(atgB-atgA)/(2.d0*pi)
!______________________________________________________________________
      return
      end 
!______________________________________________________________________
!______________________________________________________________________
      function atan_S(sn,cs)
!______________________________________________________________________
!
!     Given sine, sn, and cosine, cs produces atan(sn/cs) with
!     values in [0, 2.*pi]
!
!______________________________________________________________________
      implicit none
      real(8) :: atan_S,sn,cs
      real(8) :: atg
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
      atg = atan2(sn,cs)
      if(atg.lt.0.d0) atg = 2.d0*pi+atg
      atan_S = atg
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!______________________________________________________________________
      subroutine coeffB(xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny,B)
!______________________________________________________________________
!
!     B = int g dl
!
!______________________________________________________________________
      implicit none
      real(8) :: xa,xb,ya,yb
      real(8) :: xs,ys
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: B
      real(8) :: pi
      real(8) :: xm,ym,hns,taus
      real(8) :: tauA,tauB,diffA,diffB
      real(8) :: sn,cs,atg,atan_S
      real(8) :: dA,atanA
      real(8) :: dB,atanB
      real(8) :: int1a,int1b,int2a,int2b,int3a,int3b
!______________________________________________________________________
      pi   = acos(-1.d0)
      xm   = (xa+xb)*0.5d0
      ym   = (ya+yb)*0.5d0
      hns  = (xs-xm)*hnx + (ys-ym)*hny
      taus = (xs-xm)*taux + (ys-ym)*tauy

!!    definisco termini utili per integrale

      tauA  = -hl*0.5d0 - taus
      tauB  = hl*0.5d0 - taus
      diffA = (tauA - taus)/hns
      diffB = (tauB - taus)/hns
      
!     definisco arcotangenti
      dA = sqrt(tauA**2 + hns**2)
      sn = hns/dA
      cs = tauA/dA
      atanA = atan_S(sn,cs)
     
      dB = sqrt(tauB**2 + hns**2)
      sn = hns/dB
      cs = tauB/dB
      atanB = atan_S(sn,cs)

!     calcolo parti dell'integrale

      int1b = diffB*log(1 + diffB**2)
      int1a = diffA*log(1 + diffA**2)
      int2b = 2.d0*diffB
      int2a = 2.d0*diffA
      int3b = 2.d0*atanB
      int3a = 2.d0*atanB

!     calcolo dell'integrale
   
      B=(tauB - tauA)*1/(4.d0*pi)*log(hns)**2 +(hns/(pi*4.d0))*(int1b &
      - int2b + int3b - int1a  + int2a - int3a)       
!________________________________________________________________________
      return
      end
