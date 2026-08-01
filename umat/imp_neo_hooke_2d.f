C======================================================================
C THIS VUMAT USER MATERIAL SUBROUTINE IS WRITTEN FOR THE IMPLEMENTATION
C OF THE CONSTITUTIVE LINEAR ELASTIC MODEL. THE STRESS HAVE TO BE GIVEN
C IN SPATIAL CO-ROTATED QUANTITIES AT THE END OF THE ROUTINE. KEEP IN
C MIND THAT FOR STRUCTURAL MATERIALS ABAQUS GIVES THE ROTATED (LOCAL)
C DEFORMATION GRADIENT (F_R=R^T*F).
C
C
C  WRITTEN BY JD LAUBRIE
C======================================================================
C
C PARAMETER AND SOME VARIABLES MODULE FOR UMAT
      MODULE VARIABLES
        IMPLICIT NONE
        REAL(KIND=8),PARAMETER :: PI=3.141592653589D0
        REAL(KIND=8),PARAMETER :: RAD2GRAD=180.D0/PI
        REAL(KIND=8),PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0,
     1                            THREE=3.D0, FOUR=4.D0, FIVE=5.D0,
     2                            SIX=6.D0, SEVEN=7.D0, EIGHT=8.D0,
     3                            NINE=9.D0, TEN=10.D0, HALF=0.5D0
C
      END MODULE VARIABLES
C
C====================================================================
      SUBROUTINE UMAT(
C VARIABLES TO BE DEFINED (MODIFIABLE)
     1                STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     2                RPL,DDSDDT,DRPLDE,DRPLDT,
C PASSED IN VARIABLES (UNMODIFIABLE)
     3                STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     4                CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     5                DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     6                KSPT,JSTEP,KINC)
C
      USE VARIABLES
      INCLUDE 'ABA_PARAM.INC'
C
C--------------------------------------------------------------------
C GLOBAL ARGUMENTS
C--------------------------------------------------------------------
C VARIABLES PASSED IN FOR INFORMATION
      DIMENSION
C STRESS TENSOR
     1  STRESS(NTENS),
C STATE VARIABLES
     2  STATEV(NSTATV),
C ELASTICITY TENSOR
     3  DDSDDE(NTENS,NTENS),
C VARIATION OF STRESS INCREMENTS RESPECT TO TEMPERATURE
     4  DDSDDT(NTENS),
C VARIATION OF RPL WITH RESPECT TO THE STRAIN INCREMENTS
     5  DRPLDE(NTENS),
C TOTAL STRAINS AT THE BEGINNING OF THE INCREMENT
     6  STRAN(NTENS),
C STRAINS INCREMENTS
     7  DSTRAN(NTENS),
C TIME STEP (1) AND TOTAL TIME (2) AT THE BEGINNING OF THE INCREMENT
     8  TIMEA(2),
C PREDEFINED FIELD VARIABLES AT THE START OF THE INCREMENT
     9  PREDEF(1),
C INCREMENTS OF PREDEFINED FIELD VARIABLES
     1  DPRED(1),
C USER DEFINED MATERIAL PROPERTIES
     2  PROPS(NPROPS),
C COORDINATES OF THIS POINT
     3  COORDS(3),
C ROTATION INCREMENT MATRIX
     4  DROT(3,3),
C DEFORMATION GRADIENT AT THE BEGINNING OF THE INCREMENT
     5  DFGRD0(3,3),
C DEFORMATION GRADIENT AT THE END OF THE INCREMENT
     6  DFGRD1(3,3),
C STEP NUMBER (1), PROCEDURE TYPE KEY (2), BOOLEAN NLGEOM (3), BOOLEAN LINEAR PERTURBATION (4)
     7  JSTEP(4)
C USER DEFINED MATERIAL NAME
      CHARACTER*8 CMNAME
C
C--------------------------------------------------------------------
C LOCAL ARGUMENTS
C--------------------------------------------------------------------
      DIMENSION DISTGR(3,3), BBAR(6)
C
C--------------------------------------------------------------------
C WARNING FOR THE NUMBER OF DIMENSIONS
C THIS UMAT IS FOR PLANE STRESS
C--------------------------------------------------------------------
C      WRITE(7,*) 'MATERIAL NAME = ', CMNAME              ! string with length = 8
C      WRITE(7,*) 'NUMBER OF STRESS COMPONENTS = ', NTENS
      IF ((NDIR.NE.3).AND.(NSHR.NE.1)) THEN
        WRITE(7,*) 'NUMBER OF DIRECT STRESS COMPONENTS = ', NDI
        WRITE(7,*) 'NUMBER OF SHEAR STRESS COMPONENTS = ', NSHR
        WRITE(7,*) 'THIS UMAT MAY ONLY BE USED FOR ELEMENTS'
        WRITE(7,*) 'WITH THREE DIRECT AND ONE SHEAR STRESS COMPONENTS'
        CALL XIT
      END IF
C
C--------------------------------------------------------------------
C ELASTIC PROPERTIES FOR THE MODEL
C--------------------------------------------------------------------
      EMODULO=PROPS(1)
      POISSON=PROPS(2)
      ASTRETCH=ONE
C
C--------------------------------------------------------------------
C KINEMATICS
C--------------------------------------------------------------------
C JACOBIAN OF THE DEFORMATION GRADIENT
      DET=DFGRD1(1,1)*DFGRD1(2,2)*DFGRD1(3,3)
     1   -DFGRD1(1,2)*DFGRD1(2,1)*DFGRD1(3,3)
      IF(NSHR.EQ.3) THEN
        DET=DET+DFGRD1(1,2)*DFGRD1(2,3)*DFGRD1(3,1)
     1         +DFGRD1(1,3)*DFGRD1(3,2)*DFGRD1(2,1)
     2         -DFGRD1(1,3)*DFGRD1(3,1)*DFGRD1(2,2)
     3         -DFGRD1(2,3)*DFGRD1(3,2)*DFGRD1(1,1)
      END IF
C
C DISTORTED DEFORMATION
      DISTGR=DET**(-ONE/THREE)*DFGRD1
C
C COMPUTATION OF MODIFIED B
      BBAR(1)=DISTGR(1,1)**2+DISTGR(1,2)**2+DISTGR(1,3)**2
      BBAR(2)=DISTGR(2,1)**2+DISTGR(2,2)**2+DISTGR(2,3)**2
      BBAR(3)=DISTGR(3,3)**2+DISTGR(3,1)**2+DISTGR(3,2)**2
      BBAR(4)=DISTGR(1,1)*DISTGR(2,1)+DISTGR(1,2)*DISTGR(2,2)
     1       +DISTGR(1,3)*DISTGR(2,3)
      IF (NSHR.EQ.3) THEN
        BBAR(5)=DISTGR(1,1)*DISTGR(3,1)+DISTGR(1,2)*DISTGR(3,2)
     1         +DISTGR(1,3)*DISTGR(3,3)
        BBAR(6)=DISTGR(2,1)*DISTGR(3,1)+DISTGR(2,2)*DISTGR(3,2)
     1         +DISTGR(2,3)*DISTGR(3,3)
      END IF
C
      BINV1=(BBAR(1)+BBAR(2)+BBAR(3))/THREE
C
C--------------------------------------------------------------------
C ELASTIC PROPERTIES FOR NEO-HOOKEAN MODEL
C--------------------------------------------------------------------
      SHEAR2=EMODULO/(ONE+POISSON)
      AKAPPA=EMODULO/((ONE-TWO*POISSON)*THREE)
C
      SHEAR=HALF*SHEAR2/DET
      PR=AKAPPA*(DET-ONE)
      SHEAR23=SHEAR*(TWO/THREE)
      EK=AKAPPA*(TWO*DET-ONE)
C
C--------------------------------------------------------------------
C CALCULATE CAUCHY STRESS
C--------------------------------------------------------------------
      DO K1=1,NDI
        STRESS(K1)=SHEAR*(BBAR(K1)-BINV1) + PR
      END DO
      DO K1=NDI+1,NDI+NSHR
        STRESS(K1)=SHEAR*BBAR(K1)
      END DO
C
C--------------------------------------------------------------------
C COMPUTATION OF ELASTICITY
C--------------------------------------------------------------------
      DDSDDE(1,1)= SHEAR23*(BBAR(1)+BINV1)+EK
      DDSDDE(2,2)= SHEAR23*(BBAR(2)+BINV1)+EK
      DDSDDE(3,3)= SHEAR23*(BBAR(3)+BINV1)+EK
      DDSDDE(1,2)=-SHEAR23*(BBAR(1)+BBAR(2)-BINV1)+EK
      DDSDDE(1,3)=-SHEAR23*(BBAR(1)+BBAR(3)-BINV1)+EK
      DDSDDE(2,3)=-SHEAR23*(BBAR(2)+BBAR(3)-BINV1)+EK
      DDSDDE(1,4)= SHEAR23*BBAR(4)/TWO
      DDSDDE(2,4)= SHEAR23*BBAR(4)/TWO
      DDSDDE(3,4)=-SHEAR23*BBAR(4)
      DDSDDE(4,4)= SHEAR*(BBAR(1)+BBAR(2))/TWO
      IF (NSHR.EQ.3) THEN
        DDSDDE(1,5)= SHEAR23*BBAR(5)/TWO
        DDSDDE(2,5)=-SHEAR23*BBAR(5)
        DDSDDE(3,5)= SHEAR23*BBAR(5)/TWO
        DDSDDE(4,5)= SHEAR*BBAR(6)/TWO
        DDSDDE(5,5)= SHEAR*(BBAR(1)+BBAR(3))/TWO
        DDSDDE(1,6)=-SHEAR23*BBAR(6)
        DDSDDE(2,6)= SHEAR23*BBAR(6)/TWO
        DDSDDE(3,6)= SHEAR23*BBAR(6)/TWO
        DDSDDE(4,6)= SHEAR*BBAR(5)/TWO
        DDSDDE(5,6)= SHEAR*BBAR(4)/TWO
        DDSDDE(6,6)= SHEAR*(BBAR(2)+BBAR(3))/TWO
      END IF
C
      DO K1=1,NTENS
        DO K2=1,K1-1
          DDSDDE(K1,K2)=DDSDDE(K2,K1)
        END DO
      END DO
C
C--------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE UMAT
C
C====================================================================
