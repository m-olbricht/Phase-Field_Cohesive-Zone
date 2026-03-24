!DEC$ FREEFORM
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
! MODULE provides interface to Abaqus-implementation
! - especially the kinds of real and integer variables that Abaqus defines
!   implicitely in ABA_PARAM.INC. The corresponding variables AbqRK and AbqIK are
!   read from the implcitely defined variables i and r. Using AbqRK and AbqIK own
!   routines can use explcit type declaration
! - headers to the Abaqus utility routines
!
! Geralf Hütter, 2013/09/04
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

MODULE ABQINTERFACE
   INCLUDE 'ABA_PARAM.INC'
   PRIVATE
   INTEGER,PARAMETER,PUBLIC::AbqRK=KIND(r),AbqIK=KIND(i)
   PUBLIC::SPRINC,SINV,SPRIND,ROTSIG
   INTERFACE
     SUBROUTINE SPRINC(S,PS,LSTR,NDI,NSHR)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION::S(*),PS(*)
     END SUBROUTINE
     SUBROUTINE SINV(STRESS,SINV1,SINV2,NDI,NSHR)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION::STRESS(*)
     END SUBROUTINE
     SUBROUTINE SPRIND(S,PS,AN,LSTR,NDI,NSHR)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION::S(*),PS(3),AN(3,3)
     END SUBROUTINE
     SUBROUTINE ROTSIG(S,R,SPRIME,LSTR,NDI,NSHR)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION::S(*),R(*),SPRIME(*)
     END SUBROUTINE
   END INTERFACE
END MODULE
