!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! contact energy module
!
! Stephan Roth, TU Bergakademie Freiberg, 08.01.2024
!
! 17.01.2024: contact energy to prevent interpenetration (negative normal separation)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ContactEnergyModule

  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION ContactED(D,normalDirectionIndex,sep,npar,par)
    ! contact energy density

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex, npar
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), par(npar)
      REAL(kind=AbqRK) :: normalSep, prop_delta0, prop_t0, prop_compressionstiffness
      !
      prop_delta0               = par(1)
      prop_t0                   = par(2)
      prop_compressionstiffness = par(3)
      !
      normalSep = sep(normalDirectionIndex)
      !
      ContactED = zero
      IF (normalSep .LT. zero) THEN
        ContactED = prop_t0/prop_delta0*prop_compressionstiffness/two * normalSep**2
      END IF

    END FUNCTION ContactED

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_ContactED_d_normalSep(D,normalDirectionIndex,sep,npar,par)
    ! contact energy density: first derivative

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex, npar
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), par(npar)
      REAL(kind=AbqRK) :: normalSep, prop_delta0, prop_t0, prop_compressionstiffness
      !
      prop_delta0               = par(1)
      prop_t0                   = par(2)
      prop_compressionstiffness = par(3)
      !
      normalSep = sep(normalDirectionIndex)
      !
      d_ContactED_d_normalSep = zero
      IF (normalSep .LT. zero) THEN
        d_ContactED_d_normalSep = prop_t0/prop_delta0*prop_compressionstiffness * normalSep
      END IF

    END FUNCTION d_ContactED_d_normalSep

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_ContactED_d_normalSep_d_normalSep(D,normalDirectionIndex,sep,npar,par)
    ! contact energy density: second derivative

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex, npar
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), par(npar)
      REAL(kind=AbqRK) :: normalSep, prop_delta0, prop_t0, prop_compressionstiffness
      !
      prop_delta0               = par(1)
      prop_t0                   = par(2)
      prop_compressionstiffness = par(3)
      !
      normalSep = sep(normalDirectionIndex)
      !
      d_ContactED_d_normalSep_d_normalSep = zero
      IF (normalSep .LT. zero) THEN
        d_ContactED_d_normalSep_d_normalSep = prop_t0/prop_delta0*prop_compressionstiffness
      END IF

    END FUNCTION d_ContactED_d_normalSep_d_normalSep

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION tracWithContact(D,normalDirectionIndex,trac,sep,npar,par)
    ! consider negative normal separation / compression in traction vector

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex, npar
      REAL(kind=AbqRK), INTENT(IN) :: trac(D), sep(D), par(npar)
      REAL(kind=AbqRK) :: normalSep
      DIMENSION tracWithContact(D)
      !
      normalSep = sep(normalDirectionIndex)
      !
      tracWithContact = trac
      IF (normalSep .LT. zero) THEN
        tracWithContact(normalDirectionIndex) = d_ContactED_d_normalSep(D,normalDirectionIndex,sep,npar,par)
      END IF

    END FUNCTION tracWithContact

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION tangentWithContact(D,normalDirectionIndex,K,sep,npar,par)
    ! consider negative normal separation / compression in tangent

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex, npar
      REAL(kind=AbqRK), INTENT(IN) :: K(D,D), sep(D), par(npar)
      REAL(kind=AbqRK) :: normalSep
      DIMENSION tangentWithContact(D,D)
      !
      normalSep = sep(normalDirectionIndex)
      !
      tangentWithContact = K
      IF (normalSep .LT. zero) THEN
        tangentWithContact(normalDirectionIndex,:)                    = zero
        tangentWithContact(:,normalDirectionIndex)                    = zero
        tangentWithContact(normalDirectionIndex,normalDirectionIndex) = d_ContactED_d_normalSep_d_normalSep(D,normalDirectionIndex,sep,npar,par)
      END IF

    END FUNCTION tangentWithContact

!------------------------------------------------------------------------------------

END MODULE ContactEnergyModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

