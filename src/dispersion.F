!------------------------------------------------------------------------------------
!
!      FILE dispersion.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
! 
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!  
!-------------------------------------------------------------------------------------
!
!    CAL_DISPERSION is subroutine to calculation dispersion terms
!    so far V^4 and V^1
!    called by
!       MAIN
!    call DERIVATIVE_XX
!         DERIVATIVE_XY
!    
!    HISTORY: 
!      05/01/2010 Fengyan Shi
!      10/14/2012 Fengyan Shi, added coupling bc,
!                              change derivative_xx_high to second order
!                              according to Harris suggestion
!      08/06/2015 - 08/18/2015 Young-Kwang Choi, modified t-derivatives
!                   corrected U1p,V1p, V2, V3 and omega_1 terms
!
!-------------------------------------------------------------------------------------
SUBROUTINE CAL_DISPERSION
     USE GLOBAL
     IMPLICIT NONE

     REAL(SP),Dimension(Mloc,Nloc) :: DU,DV,DUt,DVt
     REAL(SP) :: UxxVxy,UxyVyy,HUxxHVxy,HUxyHVyy, &
                 UxxVxy_x,UxxVxy_y,UxyVyy_x,UxyVyy_y, &
                 HUxxHVxy_x,HUxxHVxy_y,HUxyHVyy_x,HUxyHVyy_y, &
                 rh,rhx,rhy,reta,ken1,ken2,ken3,ken4,ken5
# if defined (CARTESIAN)
     REAL(SP) :: omega_0,omega_1
     REAL(SP),Dimension(Mloc,Nloc) :: omega            
# endif
    
! uxx
    CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,U,Uxx)
! uxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,U,Uxy)
! vxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,V,Vxy)
! vyy
    CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,V,Vyy)

    IF(SHOW_BREAKING)THEN
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,Eta,ETAx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,Eta,ETAy)
    ENDIF

! DU DV
     DO J=1,Nloc-1
     DO I=1,Mloc-1
       DU(I,J)=Max(Depth(I,J),MinDepthFrc)*U(I,J)
       DV(I,J)=Max(Depth(I,J),MinDepthFrc)*V(I,J)

! ykchoi (15. 08. 18)
! Computation of Etat for ( U1p )_t in conservative form of Shi et al. (2012)
! and for viscosity of wave maker
! 05/25/2018, ykchoi pointed a bug: J=1,Nloc, I=1,Mloc

# if defined (CARTESIAN)
       ETAT(I,J)=-(P(I+1,J)-P(I,J))/DX-(Q(I,J+1)-Q(I,J))/DY
# else
       ETAT(I,J)=-(P(I+1,J)-P(I,J))/DX(I,J)-(Q(I,J+1)-Q(I,J))/DY(I,J)
# endif
     ENDDO
     ENDDO 

! DUxx
    CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DU,DUxx)
! DUxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DU,DUxy)
! DVxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DV,DVxy)
! DVyy
    CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DV,DVyy)
            
! this may affect parallel version
! I added coupling boundary 10/14/2012
!  add left_bc wavemaker 09/12/2017

!  boundary conditions
# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
    IF (WaveMaker(1:11)=='LEFT_BC_IRR')THEN
     ! do nothing
    ELSE
     DO J=1,Nloc
       Uxy(Ibeg,J)=ZERO
       DUxy(Ibeg,J)=ZERO
       Vxy(Ibeg,J)=ZERO
       DVxy(Ibeg,J)=ZERO
       Utxy(Ibeg,J)=ZERO
       DUtxy(Ibeg,J)=ZERO
       Vtxy(Ibeg,J)=ZERO
       DVtxy(Ibeg,J)=ZERO
     ENDDO
    ENDIF ! left_bc wavemaker
# if defined(COUPLING)  
   ENDIF 
# endif
  
# if defined (PARALLEL)
    endif  
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     DO J=1,Nloc
       Uxy(Iend,J)=ZERO
       DUxy(Iend,J)=ZERO
       Vxy(Iend,J)=ZERO
       DVxy(Iend,J)=ZERO
       Utxy(Iend,J)=ZERO
       DUtxy(Iend,J)=ZERO
       Vtxy(Iend,J)=ZERO
       DVtxy(Iend,J)=ZERO
     ENDDO 
# if defined(COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif  
# endif
  
# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     DO I=1,Mloc
       Uxy(I,Jbeg)=ZERO
       DUxy(I,Jbeg)=ZERO
       Vxy(I,Jbeg)=ZERO
       DVxy(I,Jbeg)=ZERO
       Utxy(I,Jbeg)=ZERO
       DUtxy(I,Jbeg)=ZERO
       Vtxy(I,Jbeg)=ZERO
       DVtxy(I,Jbeg)=ZERO
     ENDDO   
# if defined(COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif  
# endif

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     DO I=1,Mloc
       Uxy(I,Jend)=ZERO
       DUxy(I,Jend)=ZERO
       Vxy(I,Jend)=ZERO
       DVxy(I,Jend)=ZERO
       Utxy(I,Jend)=ZERO
       DUtxy(I,Jend)=ZERO
       Vtxy(I,Jend)=ZERO
       DVtxy(I,Jend)=ZERO
     ENDDO 
# if defined(COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif  
# endif
     
    CALL EXCHANGE_DISPERSION
     
! calculate V1p  without nonlinear dispersion
     DO J=1,Nloc
     DO I=1,Mloc
# if defined (CARTESIAN)
       U4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*DEPTH(I,J)*DEPTH(I,J)*(Uxx(I,J)+Vxy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))
       V4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*DEPTH(I,J)*DEPTH(I,J)*(Uxy(I,J)+Vyy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))      

# endif

# if defined (ZALPHA)
       U4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)  &
               *DEPTH(I,J)*DEPTH(I,J)*(Uxx(I,J)+Vxy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))
       V4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)  &
               *DEPTH(I,J)*DEPTH(I,J)*(Uxy(I,J)+Vyy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))
# endif


! ykchoi( 15. 08. 06.)
! U1p, V1p terms are modified to 0.5_SP*(1.0_SP-Beta_1) --> 0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)
               
# if defined (CARTESIAN)
       !U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
	 U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)  &
                *(Uxx(I,J)+Vxy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))
       
	 !V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
	 V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)  &
                *(Uxy(I,J)+Vyy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))
# else
       !U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
	 U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)*(Uxx(I,J)+Vxy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))

       !V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
       V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)  &
                *(Uxy(I,J)+Vyy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))
# if defined(SPH_EXTRA_DISP)
! extra terms for spherical za^2/2(-1/r0 tan theta vx +za(-1/r0 tan theta HVx), if this option if off, it will be included in a1 and b1 …
 
       U1p(I,J)=U1p(I,J)+0.5_SP*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)*(-1.0_SP/R_earth*  &
                TAN(Lat_theta(I,J))*0.5_SP*(VxL(I,J)+VxR(I,J))) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(-1.0_SP/R_earth* &
                TAN(Lat_theta(I,J))*0.5_SP*(HVxL(I,J)+HVxR(I,J)))

! extra terms in y: za^2/2(1/r0 tan theta ux -1/r0 tan theta vy -1/r0^2/cos^2theta v)
!                  +za (1/r0 tan theta HUx -1/r0 tan theta HVy -1/r0^2/cos^2 HV)
       V1p(I,J)=V1p(I,J)+0.5_SP*(1.0_SP-Beta_1) &
                 *DEPTH(I,J)*DEPTH(I,J) &
                /R_earth*( TAN(Lat_theta(I,J))*0.5_SP*(UxL(I,J)+UxR(I,J)) &
                          -TAN(Lat_theta(I,J))*0.5_SP*(VyL(I,J)+VyR(I,J)) &
                          -V(I,J)/R_earth/COS(Lat_theta(I,J))/COS(Lat_theta(I,J)) ) &
               +(Beta_1-1.0_SP)*DEPTH(I,J) &
                /R_earth*( 0.5_SP*(HUxL(I,J)+HUxR(I,J))*TAN(Lat_theta(I,J)) &
                          -0.5_SP*(HVyL(I,J)+HVyR(I,J))*TAN(Lat_theta(I,J)) &
                          -HV(I,J)/R_earth/COS(Lat_theta(I,J))/COS(Lat_theta(I,J)) &        
                 )
# endif
! end spherical extra dispersion terms which can be implemented in TIME_LEFT
# endif

     ENDDO
     ENDDO

END SUBROUTINE CAL_DISPERSION
