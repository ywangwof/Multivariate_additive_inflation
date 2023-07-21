module ran
  public ran1
contains

 FUNCTION RAN1(IDUM)
!
!
  IMPLICIT NONE        ! Force explicit declarations
  INTEGER :: IDUM        ! The seed for random number generation
  REAL :: RAN1           ! The function to generate random number.
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  REAL,SAVE :: R(97)
  INTEGER :: IX1,IX2,IX3,J,IFF
  INTEGER :: M1,M2,M3,IA1,IA2,IA3,IC1,IC2,IC3
  REAL :: RM1,RM2
  SAVE IX1,IX2,IX3

  PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
  PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
  PARAMETER (M3=243000,IA3=4561,IC3=51349)
  DATA IFF /0/
!
  IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
    IFF=1
    IX1=MOD(IC1-IDUM,M1)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX2=MOD(IX1,M2)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX3=MOD(IX1,M3)
    DO J=1,97
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
    ENDDO
    IDUM=1
  ENDIF
  IX1=MOD(IA1*IX1+IC1,M1)
  IX2=MOD(IA2*IX2+IC2,M2)
  IX3=MOD(IA3*IX3+IC3,M3)
  J=1+(97*IX3)/M3
  IF(J.GT.97.OR.J.LT.1)THEN
    WRITE(*,*)'J is greater than 97 or less than 1','IDUM=',IDUM
    STOP
  ENDIF
  RAN1=R(J)
  R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1

  RETURN
  END FUNCTION RAN1

end module ran
