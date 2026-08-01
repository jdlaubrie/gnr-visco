C======================================================================
C THIS UMAT USER MATERIAL SUBROUTINE IS WRITTEN FOR THE IMPLEMENTATION
C OF CONSTITUTIVES MODELS FOR A LAMELLAR UNIT. THE STRESS HAVE TO BE GIVEN
C IN SPATIAL CO-ROTATED QUANTITIES AT THE END OF THE ROUTINE. KEEP IN
C MIND THAT FOR CONTINUUM ELEMENTS ABAQUS GIVES THE CO-ROTATED (LOCAL)
C DEFORMATION GRADIENT (F_al=R^T*F*R=U*R).
C
C THE REMODELING HERE IS WITHOUT CORRECTION
C
C  WRITTEN BY JD LAUBRIE
C======================================================================
C
C PARAMETER AND SOME VARIABLES MODULE FOR UMAT
      MODULE VARIABLES
        IMPLICIT NONE
        LOGICAL, PARAMETER :: VECN_BOOL=.FALSE.                     ! radial growth (false) normal to cell(true). latter is unstable
        INTEGER, PARAMETER :: NCELL=128,NGAUS=2                      ! number of cell elements (truss) and number of integration points in it
        REAL(KIND=8),PARAMETER :: PI=3.141592653589793D0
        REAL(KIND=8),PARAMETER :: RAD2GRAD=180.D0/PI
        REAL(KIND=8),PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0,
     1                            THREE=3.D0, FOUR=4.D0, FIVE=5.D0,
     2                            SIX=6.D0, SEVEN=7.D0, EIGHT=8.D0,
     3                            NINE=9.D0, TEN=10.D0, HALF=0.5D0,
     4                            EPS=1.D-9
C INITIAL LAGRANGIAN DENSITY = EULERIAN DENSITY (CONSTANT)
        REAL(KIND=8), PARAMETER :: DNSTY_TOT_0=1.0D-6       ! density | [kg/m3] | total
C
      END MODULE VARIABLES
C
C====================================================================
C SOLUTION-DEPENDENT STATE VARIABLES. INITIALIZATION ROUTINE
C====================================================================
      SUBROUTINE SDVINI(
     1                STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT)
C
      USE VARIABLES
      INCLUDE 'ABA_PARAM.INC'
C
C--------------------------------------------------------------------
C GLOBAL ARGUMENTS
C--------------------------------------------------------------------
	    DIMENSION STATEV(NSTATV),COORDS(NCRDS)
C
C--------------------------------------------------------------------
      STATEV(1) = ONE        ! inverse inelastic**2, circumferential
      STATEV(2) = ONE        ! inverse inelastic**2, radial
      STATEV(3) = ONE        ! inverse inelastic**2, longitudinal
      STATEV(4) = ZERO       ! inverse inelastic**2, circ-radial
      STATEV(5) = ZERO       ! inverse inelastic**2, circ-long
      STATEV(6) = ZERO       ! inverse inelastic**2, radial-long
      STATEV(7) = ONE        ! inelastic, fibres
      STATEV(8) = ONE        ! determinant inelastic, error
      STATEV(9) = ONE        ! determinant inelastic, actual growth
      STATEV(10) = ONE       ! determinant elastic
      STATEV(11) = ONE       ! determinant total
      STATEV(12) = 1.0D-6    ! total density
      STATEV(13) = HALF      ! matrix mass fraction
      STATEV(14) = HALF      ! fibres mass fraction
      STATEV(15) = ZERO      ! reference total stress magnitude, step-1
      STATEV(16) = ZERO      ! reference volumetric stress magnitude, step-1
      STATEV(17) = ZERO      ! reference matrix stress magnitude, step-1
      STATEV(18) = ZERO      ! reference fibre stress magnitude
      STATEV(19) = ZERO      ! total stress magnitude
      STATEV(20) = ZERO      ! volumetric stress magnitude
      STATEV(21) = ZERO      ! matrix stress magnitude
      STATEV(22) = ZERO      ! fibre stress magnitude
      STATEV(23) = ZERO      ! cell stress 1-d
C
      RETURN
      END
C
C====================================================================
C USER-MATERIAL ROUTINE
C ***** MIXTURE *****
C PSI = 0.5*K*(J-1)**2 + 0.5*MU*(I1bar-3) + 0.5*(K1/K2)*(EXP(K2*(I4_F-1)**2)-1)
C I4_F = (H:C)/LAMBDA_P**2 = LAMBDA_F**2/LAMBDA_P**2
C H = A0_X_A0
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
     8  TIME(2),
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
C UMAT SELECTOR, AS THIS MODEL HAS TWO CONSTITUTIVE MODELS IN IT
C 1. CELLS MODELED AS TRUSS ELEMENTS
C 2. GROUND SUBSTANCE (INTERLAMELLA) OR ELASTIN LAMELLA MODELED AS
C    A MIXTURE, THERE SHOULD BE COLLAGEN FIBERS AND SOMETHING ELSE
C--------------------------------------------------------------------
      IF (CMNAME(1:4).EQ.'CELL') THEN
        CALL UMAT_CELL(STRESS,STATEV,DDSDDE,SSE,TIME(1),DTIME,NDI,
     1               NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DFGRD1,
     2               NOEL,NPT,JSTEP(1))
      ELSE
        CALL UMAT_GRND(STRESS,STATEV,DDSDDE,SSE,TIME(1),DTIME,NDI,
     1               NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DFGRD1,
     2               NOEL,NPT,JSTEP(1))
      END IF
C
C--------------------------------------------------------------------
      RETURN
      END SUBROUTINE UMAT
C
C====================================================================
C PSI_PAS:truss = 0.5*MU*(I_1-3) + P*(J-1); P: LAGRANGE MULTIPLIER
C PSI_ACT(L_ACT) = S_MAX*(L_ACT+(1/3)*(L_MAX-L_ACT)**3/(L_MAX-L_MIN)**2)
C L_ACT = STRETCH but it could also be L_ACT=STRETCH/Ĺ_ACT
C====================================================================
      SUBROUTINE UMAT_CELL(STRESS,STATEV,DDSDDE,SSE,TIME1,DTIME,NDI,
     1                  NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DFGRD1,
     2                  NOEL,NPT,JSTEP1)
C
      USE VARIABLES
      INCLUDE 'ABA_PARAM.INC'
C--------------------------------------------------------------------
C GLOBAL ARGUMENTS
C--------------------------------------------------------------------
C VARIABLES PASSED IN FOR INFORMATION
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),
     1          PROPS(NPROPS), COORDS(3), DFGRD1(3,3)
C
C--------------------------------------------------------------------
C WARNING WITH THE NUMBER OF CELLS
C--------------------------------------------------------------------
      IF (NOEL.GT.NCELL) THEN
        WRITE(7,*) 'IN THIS UMAT THE NUMBER OF CELLS DECLARED IN THE'
        WRITE(7,*) 'MODULE MUST BE EQUAL TO THE CELLS IN THE ASSEMBLY'
        CALL XIT
      END IF
C
C--------------------------------------------------------------------
C WARNING WITH THE NUMBER OF DIMENSIONS
C THIS UMAT IS FOR TRUSS ELEMENTS
C--------------------------------------------------------------------
      IF (NSHR.NE.0) THEN
        WRITE(7,*) 'NUMBER OF DIRECT STRESS COMPONENTS = ', NDI
        WRITE(7,*) 'NUMBER OF SHEAR STRESS COMPONENTS = ', NSHR
        WRITE(7,*) 'THIS UMAT MAY ONLY BE USED FOR ELEMENTS'
        WRITE(7,*) 'WITH ONE DIRECT AND ZERO SHEAR STRESS COMPONENTS'
        CALL XIT
      END IF
C
C--------------------------------------------------------------------
C DEFINITON OF CONSTANTS FOR THE MODEL
C--------------------------------------------------------------------
C CELL PROPERTIES
      EMOD_M=PROPS(1)
C
C APPLICATION OF ACTIVATION STRESS IN THE FIRST STEP
      TIME2 = TIME1 + DTIME
      IF (JSTEP1.EQ.1) THEN
        S_MAX=TIME2*PROPS(2)
        RSTRETCH0=ONE               !+(1.D0-ONE)*TIME2   !residual stretch
      ELSE
        S_MAX=PROPS(2)
        RSTRETCH0=STATEV(7)
      END IF
      ALAMDA_MIN=PROPS(3)
      ALAMDA_MAX=PROPS(4)
C
      SHEAR0=HALF*EMOD_M/(ONE+HALF)
      SHEAR=SHEAR0
C
C--------------------------------------------------------------------
C KINEMATICS
C--------------------------------------------------------------------
      RSTRETCH=RSTRETCH0
C
      TSTRETCH = DFGRD1(1,1)
      ESTRETCH = TSTRETCH/RSTRETCH
C
C--------------------------------------------------------------------
C CALCULATE PASSIVE AND ACTIVE STRESS
C--------------------------------------------------------------------
C MUSCLE PASSIVE STRESS
      STRESS_P = SHEAR*(ESTRETCH**TWO-ONE/ESTRETCH)
      DSDE_P = SHEAR*(TWO*ESTRETCH**TWO+ONE/ESTRETCH)

C MUSCLE ACTIVE STRESS
      KMODEL=2
      STRESS_A=ZERO
      DSDE_A=ZERO
      ACTMO3=ZERO
      IF (KMODEL.EQ.1) THEN
C ACTIVATION COMPUTED FROM TOTAL STRETCH OF THE CELL
        AL_DIFF = ALAMDA_MAX-TSTRETCH
        AL_RANG = ALAMDA_MAX-ALAMDA_MIN
        ACTMO1 = TWO*S_MAX*(AL_DIFF/AL_RANG**TWO)
        ACTMO2 = S_MAX*(ONE-(AL_DIFF/AL_RANG)**TWO)
        ACTMO3 = S_MAX*(TSTRETCH
     1    + (ONE/THREE)*AL_DIFF**THREE/AL_RANG**TWO)
C
        STRESS_A = ACTMO2*TSTRETCH
        DSDE_A = ACTMO1*TSTRETCH**TWO + ACTMO2*TSTRETCH
      ELSE IF (KMODEL.EQ.2) THEN
C ACTIVATION COMPUTED FROM ELASTIC STRETCH OF THE CELL
        AL_DIFF = ALAMDA_MAX-ESTRETCH
        AL_RANG = ALAMDA_MAX-ALAMDA_MIN
        ACTMO1 = TWO*S_MAX*(AL_DIFF/AL_RANG**TWO)
        ACTMO2 = S_MAX*(ONE-(AL_DIFF/AL_RANG)**TWO)
        ACTMO3 = S_MAX*(ESTRETCH
     1    + (ONE/THREE)*AL_DIFF**THREE/AL_RANG**TWO)
C
        STRESS_A = ACTMO2*ESTRETCH
        DSDE_A = ACTMO1*ESTRETCH**TWO + ACTMO2*ESTRETCH
      ELSE
C ACTIVATION COMPUTED FROM STRETCH_ACT=1, ASSUMING ACTIVATION HAPPENED IN SHORT-TIME SCALE
        AL_DIFF = ALAMDA_MAX-ONE
        AL_RANG = ALAMDA_MAX-ALAMDA_MIN
        ACTMO1 = TWO*S_MAX*(AL_DIFF/AL_RANG**TWO)
        ACTMO2 = S_MAX*(ONE-(AL_DIFF/AL_RANG)**TWO)
        ACTMO3 = S_MAX*(ONE+(ONE/THREE)*AL_DIFF**THREE/AL_RANG**TWO)
C
        STRESS_A = ACTMO2
        DSDE_A = ACTMO1
      END IF
C--------------------------------------------------------------------
C STRESS AND ELASTICITY OUTPUT
C--------------------------------------------------------------------
      STRESS_Q1 = STRESS_P + STRESS_A
      STIFF = DSDE_P + DSDE_A
C
      STRESS(1) = STRESS_Q1
      DDSDDE(1,1) = STIFF
      SSE=HALF*SHEAR*(ESTRETCH**TWO+TWO/ESTRETCH-THREE)+ACTMO3
C
C VOLUMETRIC AND DEVIATORIC STRESS CRITERIONS
      SPRESS=(STRESS_P+STRESS_A)/THREE
      SMISES=(STRESS_P+STRESS_A)*DSQRT(TWO/THREE)
C
C UPDATE OF SOLUTION-DEPENDENT STATE VARIABLES
      STATEV(7)=RSTRETCH0
      STATEV(12)=DNSTY_TOT_0
      STATEV(15)=SPRESS
      STATEV(17)=SMISES
      STATEV(21) = STRESS_P
      STATEV(22) = STRESS_A
      STATEV(23)=STRESS_Q1
C
C--------------------------------------------------------------------
      RETURN
      END SUBROUTINE UMAT_CELL
C
C====================================================================
C PSI = 0.5*K*(J-1)**2 + 0.5*MU*(I1bar-3) + 0.5*(k1/k2)*(EXP(k2*(I4-1)**2)-1)
C====================================================================
      SUBROUTINE UMAT_GRND(STRESS,STATEV,DDSDDE,SSE,TIME1,DTIME,NDI,
     1                  NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DFGRD1,
     2                  NOEL,NPT,JSTEP1)
C
      USE VARIABLES
      INCLUDE 'ABA_PARAM.INC'
C     WARNING - the aba_param.inc file declares: Implicit real*8(a-h,o-z)
C--------------------------------------------------------------------
C GLOBAL ARGUMENTS
C--------------------------------------------------------------------
C VARIABLES PASSED IN FOR INFORMATION
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),
     1          PROPS(NPROPS), COORDS(3), DFGRD1(3,3)
C
C--------------------------------------------------------------------
C LOCAL ARGUMENTS
C--------------------------------------------------------------------
C VARIABLES CREATED INTO THIS ROUTINE
      DIMENSION ROT(3,3), DFTOT(3,3), DFTOTT(3,3), DFTOT_INV(3,3),
     1  DFTOT_INVT(3,3), BEL_AUX(3,3), BEL(3,3), DISTGR(3,3),
     2  CPL_INV(3,3), CPL_AUX(3,3), BBAR(NTENS),BBAR_NEW(NTENS),
     3  STRESS_VOL(NTENS), STRESS_ICH(NTENS), STRESS_FIB(NTENS),
     4  DSDE_VOL(NTENS,NTENS), DSDE_ICH(NTENS,NTENS), FRAC_FF(4),
     5  DSDE_FIB(NTENS,NTENS), SNORMAL_M(NTENS)
C
C--------------------------------------------------------------------
C WARNING FOR THE NUMBER OF DIMENSIONS
C THIS UMAT IS FOR TRUSS ELEMENTS
C--------------------------------------------------------------------
      IF ((NDI.NE.3).AND.(NSHR.NE.1)) THEN
        WRITE(7,*) 'NUMBER OF DIRECT STRESS COMPONENTS = ', NDI
        WRITE(7,*) 'NUMBER OF SHEAR STRESS COMPONENTS = ', NSHR
        WRITE(7,*) 'THIS UMAT MAY ONLY BE USED FOR ELEMENTS'
        WRITE(7,*) 'WITH ONE DIRECT AND ZERO SHEAR STRESS COMPONENTS'
        CALL XIT
      END IF
C
C--------------------------------------------------------------------
C PARAMETERS FOR THE MODEL
C--------------------------------------------------------------------
C ELASTIC BEHAVIOR PARAMETERS
      EMOD_0=PROPS(1)                                         ! Young modulus  | [J/kg]     | matrix
      POIS=PROPS(2)                                           ! Poisson ratio  | [-]        | matrix
      CTE1_0=PROPS(3)                                         ! k1-HGO         | [J/kg]     | fibres
      CTE2=PROPS(4)                                           ! k2-HGO         | [-]        | fibres
      ANGLE=PROPS(5)*PI/180.D0                                ! angle-diagonal | [°]->[rad] | fibres
C
      TIME2 = TIME1 + DTIME
C SELECT PARAMETERS, EITHER REFERENCE-STEP (JSTEP1=1), EITHER REMODELING-STEP (JSTEP1/=1)
      IF (JSTEP1.EQ.1) THEN
        STATEV(1)=(ONE/(ONE+(PROPS(6)-ONE)*TIME2))**TWO       ! circumferential | inverse pre-stretch | matrix
        STATEV(3)=(ONE/(ONE+(PROPS(7)-ONE)*TIME2))**TWO       ! longitudinal    | inverse pre-stretch | matrix
        STATEV(2)=ONE/(STATEV(1)*STATEV(3))                   ! thickness       | inverse pre-stretch | matrix
        RSTRETCH_F0=ONE+(PROPS(8)-ONE)*TIME2                  ! fibre           |         pre-stretch | fibres
      ELSE
        RSTRETCH_F0=STATEV(7)                                 ! fibre           |         pre-stretch | fibres
      END IF
C MATRIX AND FIBRES MASS FRACTIONS
      FRAC_MT_0=PROPS(9)                                      ! mass-fraction   | matrix
      FRAC_FT_0=PROPS(10)                                     ! mass-fraction   | fibres
C COLLAGEN FIBRES FOUR FAMILIES FRACTIONS
      FRAC_FF(1) = 0.1D0               ! mass-fraction circumferential family fiber
      FRAC_FF(2) = 0.1D0               ! mass-fraction axial family fiber
      FRAC_FF(3) = 0.4D0               ! mass-fraction diagonal symmetric family fiber (+)
      FRAC_FF(4) = 0.4D0               ! mass-fraction diagonal symmetric family fiber (-)
C
C REDEFINITION OF MASS-DEPENDENT ELASTIC PARAMETERS
      EMOD=EMOD_0*FRAC_MT_0                                 ! Young modulus  | [Pa]     | matrix
      CTE1=CTE1_0*FRAC_FT_0                                 ! k1-HGO         | [Pa]     | fibres
C
C--------------------------------------------------------------------
C KINEMATICS
C--------------------------------------------------------------------
C RIGHT DEFORMATION (U) AND LEFT ROTATION (R) ACCORDING TO ABAQUS
C CO-ROTATED LOCAL CONVENTION. F=RU and F^al=UR. NOLAN ETAL 2020.
      CALL RIGHTROTATION(ROT,DFTOT,DFGRD1)
      IF (ISNAN(DFTOT(1,1))) THEN
        WRITE(7,*) 'POLAR DECOMPOSITION FAILED, F(1,1)= ', DFTOT(1,1)
        CALL XIT
      END IF
C
C JACOBIAN OF THE TOTAL DEFORMATION GRADIENT
      DET=(DFTOT(1,1)*DFTOT(2,2)-DFTOT(1,2)*DFTOT(2,1))*DFTOT(3,3)
      IF(NSHR.EQ.3) THEN
        DET=DET+DFTOT(1,2)*DFTOT(2,3)*DFTOT(3,1)
     1         +DFTOT(1,3)*DFTOT(3,2)*DFTOT(2,1)
     2         -DFTOT(1,3)*DFTOT(3,1)*DFTOT(2,2)
     3         -DFTOT(2,3)*DFTOT(3,2)*DFTOT(1,1)
      END IF
C
C--------------------------------------------------------------------
C DEFORMATION DECOMPOSITION, REMODELING AND ELASTIC PARTS
C GROWTH HERE IS CONSIDERED PERPENDICULAR TO THE FIBER, THEN IN DIRECTION 22
C OVERALL INELASTIC DEFORMATION, C_p=F_g.T*C_r*F_g
C
C MATRIX INVERSE RIGHT-CACHY REMODELING, OR INVERSE ISOCHORIC INELASTIC DEFORMATION
      CPL_INV(1,1)=STATEV(1)
      CPL_INV(2,2)=STATEV(2)
      CPL_INV(3,3)=STATEV(3)
      CPL_INV(1,2)=STATEV(4)
      CPL_INV(1,3)=STATEV(5)
      CPL_INV(2,3)=STATEV(6)
      CPL_INV(2,1)=STATEV(4)
      CPL_INV(3,1)=STATEV(5)
      CPL_INV(3,2)=STATEV(6)
C
C--------------------------------------------------------------------
C CALCULATE STRESS FROM CONSTITUTIVE LAW
C--------------------------------------------------------------------
      CALL ELASTIC_LAW(BBAR,STRESS_VOL,STRESS_ICH,STRESS_FIB,DSDE_VOL,
     1      DSDE_ICH,DSDE_FIB,ENERGY,DETEL,EKAPA,SHEAR,STIFF,NDI,NSHR,
     2      NTENS,DFTOT,CPL_INV,EMOD,POIS,FRAC_FF,CTE1,CTE2,
     3      ANGLE,RSTRETCH_F0)
C
C--------------------------------------------------------------------
C ADDITION OF STRESSES BY COMPONENT
C--------------------------------------------------------------------
      DO K1=1,NTENS
        STRESS(K1)=(STRESS_VOL(K1)+STRESS_ICH(K1)+STRESS_FIB(K1))/DET
      END DO
C
      DO K1=1,NTENS
        DO K2=K1,NTENS
          DDSDDE(K1,K2)=(DSDE_VOL(K1,K2)+DSDE_ICH(K1,K2)
     1                  +DSDE_FIB(K1,K2))/DET
        END DO
      END DO
      DO K1=1,NTENS
        DO K2=1,K1-1
          DDSDDE(K1,K2)=DDSDDE(K2,K1)
        END DO
      END DO
C
C--------------------------------------------------------------------
C NORM OF THE STRESS BY COMPONENT. TOTAL, VOLUME, MATRIX AND FIBRE
C--------------------------------------------------------------------
C QUADRATIC STRESSES
C      STRESS_TQ2=SNORM2_VOIGT(STRESS,NDI,NSHR,NTENS)
C      STRESS_VQ2=SNORM2_VOIGT(STRESS_VOL,NDI,NSHR,NTENS)   ! this may not have much sense, better use stress_vsph, spheric component
      STRESS_MQ2=SNORM2_VOIGT(STRESS_ICH,NDI,NSHR,NTENS)
      STRESS_FQ2=SNORM2_VOIGT(STRESS_FIB,NDI,NSHR,NTENS)
C
C SPHERICAL STRESS, FIRST INVARIANT, tr(stress_ich)=0
      STRESS_I1_MAT=STRESS_VOL(1)+STRESS_VOL(2)+STRESS_VOL(3)
      STRESS_I1_FIB=STRESS_FIB(1)+STRESS_FIB(2)+STRESS_FIB(3)
C
C DEVIATOR STRESS, SECOND INVARIANT
      STRESS_J2_MAT=STRESS_MQ2
      STRESS_J2_FIB=STRESS_FQ2-(STRESS_I1_FIB**TWO)/THREE
C
C VOLUMETRIC AND DEVIATORIC STRESS CRITERIONS
      SPRESS_MAT=STRESS_I1_MAT/THREE
      SPRESS_FIB=STRESS_I1_FIB/THREE
      SMISES_MAT=DSQRT(STRESS_J2_MAT)
      SMISES_FIB=DSQRT(STRESS_J2_FIB)
C
      STRESS_MQ1=DSQRT(STRESS_MQ2)
C--------------------------------------------------------------------
C STATE VARIABLES UPDATE. FIBRE STRESS
C--------------------------------------------------------------------
      SSE=ENERGY
      STATEV(1) = CPL_INV(1,1)
      STATEV(2) = CPL_INV(2,2)
      STATEV(3) = CPL_INV(3,3)
      STATEV(4) = HALF*(CPL_INV(1,2)+CPL_INV(2,1))
      STATEV(5) = HALF*(CPL_INV(1,3)+CPL_INV(3,1))
      STATEV(6) = HALF*(CPL_INV(2,3)+CPL_INV(3,2))
      STATEV(7) = RSTRETCH_F0
      STATEV(10) = DETEL
      STATEV(11) = DET
      STATEV(12) = DNSTY_TOT_0
      STATEV(13) = FRAC_MT_0
      STATEV(14) = FRAC_FT_0
      IF (JSTEP1.EQ.1) THEN
        STATEV(15) = SPRESS_MAT/DET
        STATEV(16) = SPRESS_FIB/DET
        STATEV(17) = SMISES_MAT/DET
        STATEV(18) = SMISES_FIB/DET
      END IF
      STATEV(19) = SPRESS_MAT/DET
      STATEV(20) = SPRESS_FIB/DET
      STATEV(21) = SMISES_MAT/DET
      STATEV(22) = SMISES_FIB/DET
      STATEV(23) = ZERO  !ESTRES_SM(1)
C
C--------------------------------------------------------------------
      RETURN
      END SUBROUTINE UMAT_GRND
C
C====================================================================
C ISOTROPIC PART OF THE STRESS AND ELASTICITY REPRESENTED BY A
C NEO-HOOKEAN MODEL.
C PSI = 0.5*K*(J-1)**2 + 0.5*MU*(I1bar-3) 0.5*(k1/k2)*{exp[k2*(I4-1)**2]-1}
C====================================================================
      SUBROUTINE ELASTIC_LAW(BBAR,STRESS_VOL,STRESS_ICH,STRESS_FIB,
     1    DSDE_VOL,DSDE_ICH,DSDE_FIB,SEF_TOT,DETEL,EK,SHEAR,STIFF,NDI,
     2    NSHR,NTENS,DFTOT,CPL_INV,EMOD,POIS,FRAC,CTE1,CTE2,
     3    ANGLE,RSTRETCH_F0)
C
      USE VARIABLES
      INCLUDE 'ABA_PARAM.INC'
C--------------------------------------------------------------------
C GLOBAL ARGUMENTS
C--------------------------------------------------------------------
C VARIABLES PASSED IN FOR INFORMATION
      DIMENSION DFTOT(3,3), CPL_INV(3,3), BBAR(NTENS),
     1  STRESS_VOL(NTENS), STRESS_ICH(NTENS), STRESS_FIB(NTENS),
     2  DSDE_VOL(NTENS,NTENS), DSDE_ICH(NTENS,NTENS), FRAC(4),
     3  DSDE_FIB(NTENS,NTENS)
C
C--------------------------------------------------------------------
C LOCAL ARGUMENTS
C--------------------------------------------------------------------
      DIMENSION DFTOTT(3,3), BEL_AUX(3,3), BEL(3,3), DISTGR(3,3),
     1  DFREM(3,3), DFGNR(3,3), DFREM_VEC(3), THETA(4), AM0(3), AM(3),
     1  GST(NTENS), GST0(3,3), AIDEN3(3,3)
C
C--------------------------------------------------------------------
C INITIALIZATION OF STRESS AND ELASTICITY IN THE FIBER
      DO K1=1,NTENS
        STRESS_VOL(K1)=ZERO
        STRESS_ICH(K1)=ZERO
        STRESS_FIB(K1)=ZERO
        DO K2=1,NTENS
          DSDE_VOL(K1,K2)=ZERO
          DSDE_ICH(K1,K2)=ZERO
          DSDE_FIB(K1,K2)=ZERO
        END DO
      END DO
C
C DEFINE THE 2-RANK IDENTITY TENSOR FOR 3-DIMENSIONS
      DO K1=1,3
        DO K2=1,3
          AIDEN3(K1,K2)=ZERO
        END DO
        AIDEN3(K1,K1)=ONE
      END DO
C--------------------------------------------------------------------
C CALCULATE ISOTROPIC PART FOR STRESS AND ELASTICITY
C--------------------------------------------------------------------
C TOTAL ELASTIC LEFT-CAUCHY GREEN DEFORMATION. NOT FORCING THIRD DEFORMATION TO ONE
      DFTOTT=TRANSPOSE(DFTOT)
      BEL_AUX=MATMUL(DFTOT,CPL_INV)
      BEL=MATMUL(BEL_AUX,DFTOTT)
C JACOBIAN OF THE ELASTIC PART (MATRIX) DEFORMATION GRADIENT
      DETEL2=(BEL(1,1)*BEL(2,2)-BEL(1,2)*BEL(2,1))*BEL(3,3)
      IF(NSHR.EQ.3) THEN
        DETEL2=DETEL2+BEL(1,2)*BEL(2,3)*BEL(3,1)
     1               +BEL(1,3)*BEL(3,2)*BEL(2,1)
     2               -BEL(1,3)*BEL(3,1)*BEL(2,2)
     3               -BEL(2,3)*BEL(3,2)*BEL(1,1)
      END IF
      DETEL=DSQRT(DETEL2)
C
C DISTORTED DEFORMATION
      DISTGR=DETEL2**(-ONE/THREE)*BEL
C
C ISOCHORIC ELASTIC LEFT-CAUCHY GREEN DEFORMATION, B-BAR
      BBAR(1)=DISTGR(1,1)
      BBAR(2)=DISTGR(2,2)
      BBAR(3)=DISTGR(3,3)
      BBAR(4)=HALF*(DISTGR(2,1)+DISTGR(1,2))
      IF (NSHR.EQ.3) THEN
        BBAR(5)=HALF*(DISTGR(3,1)+DISTGR(1,3))
        BBAR(6)=HALF*(DISTGR(3,2)+DISTGR(2,3))
      END IF
C
C STIFFNESS MODIFIED BY THE VOLUME FRACTION
      SHEAR2=EMOD/(ONE+POIS)
      AKAPPA=EMOD/(THREE-SIX*POIS)
C
C PARAMETERS FOR ISOTROPIC BEHAVIOR
      SHEAR=HALF*SHEAR2
      PR=AKAPPA*(DETEL-ONE)*DETEL
      SHEAR23=SHEAR*(TWO/THREE)
      EK=AKAPPA*(TWO*DETEL-ONE)*DETEL
C
C FIRST INVARIANT FOR B
      BINV1=(BBAR(1)+BBAR(2)+BBAR(3))/THREE
C
C STRESS VOIGT VECTOR
      DO K1=1,NDI
        STRESS_VOL(K1)=PR
        STRESS_ICH(K1)=SHEAR*(BBAR(K1)-BINV1)
      END DO
      DO K1=NDI+1,NDI+NSHR
        STRESS_ICH(K1)=SHEAR*BBAR(K1)
      END DO
C
C ELASTICITY VOIGT MATRIX
      DO K1=1,3
        DO K2=1,3
          DSDE_VOL(K1,K2)= EK
        END DO
      END DO
C
      DSDE_ICH(1,1)= SHEAR23*(BBAR(1)+BINV1)
      DSDE_ICH(2,2)= SHEAR23*(BBAR(2)+BINV1)
      DSDE_ICH(3,3)= SHEAR23*(BBAR(3)+BINV1)
      DSDE_ICH(1,2)=-SHEAR23*(BBAR(1)+BBAR(2)-BINV1)
      DSDE_ICH(1,3)=-SHEAR23*(BBAR(1)+BBAR(3)-BINV1)
      DSDE_ICH(2,3)=-SHEAR23*(BBAR(2)+BBAR(3)-BINV1)
      DSDE_ICH(1,4)= SHEAR23*BBAR(4)/TWO
      DSDE_ICH(2,4)= SHEAR23*BBAR(4)/TWO
      DSDE_ICH(3,4)=-SHEAR23*BBAR(4)
      DSDE_ICH(4,4)= SHEAR*(BBAR(1)+BBAR(2))/TWO
      IF (NSHR.EQ.3) THEN
        DSDE_ICH(1,5)= SHEAR23*BBAR(5)/TWO
        DSDE_ICH(2,5)=-SHEAR23*BBAR(5)
        DSDE_ICH(3,5)= SHEAR23*BBAR(5)/TWO
        DSDE_ICH(4,5)= SHEAR*BBAR(6)/TWO
        DSDE_ICH(5,5)= SHEAR*(BBAR(1)+BBAR(3))/TWO
        DSDE_ICH(1,6)=-SHEAR23*BBAR(6)
        DSDE_ICH(2,6)= SHEAR23*BBAR(6)/TWO
        DSDE_ICH(3,6)= SHEAR23*BBAR(6)/TWO
        DSDE_ICH(4,6)= SHEAR*BBAR(5)/TWO
        DSDE_ICH(5,6)= SHEAR*BBAR(4)/TWO
        DSDE_ICH(6,6)= SHEAR*(BBAR(2)+BBAR(3))/TWO
      END IF
C
      SEF_VOL=HALF*AKAPPA*(DETEL-ONE)**TWO
      SEF_ICH=HALF*SHEAR*THREE*(BINV1-ONE)
C
C--------------------------------------------------------------------
C FAMILY FIBERS, ANISOTROPIC PART FOR STRESS AND ELASTICITY
C--------------------------------------------------------------------
      THETA(1)=ZERO                  ! circumferential fiber
      THETA(2)=HALF*PI               ! axial fiber
      THETA(3)=ANGLE                 ! symmetric diagonal fiber (+)
      THETA(4)=-ANGLE                ! symmetric diagonal fiber (-)
C
      STIFF=ZERO
      DO K3=1,4
C FIBER DIRECTION IN REFERENCE CONFIGURATION
        AM0(1)=DCOS(THETA(K3))
        AM0(2)=ZERO
        AM0(3)=DSIN(THETA(K3))
C
        DO K1=1,3
          DO K2=1,3
            GST0(K1,K2)=AM0(K1)*AM0(K2)
          END DO
        END DO
C
C FIBRE REMODELING DEFORMATION GRADIENT, OR ISOCHORIC INELASTIC DEFORMATION
        RSTRETCH_F=RSTRETCH_F0        !no growth for fibres
C
C FIBRE DIRECTION VECTOR
        DO K1=1,3
          AM(K1) = ZERO
          DO K2=1,3
            AM(K1) = AM(K1) + DFTOT(K1,K2)*AM0(K2)
          END DO
        END DO
C
C GENERAL STRUCTURAL TENSOR. (C:m0xm0)/lambda_p^2. total/plastic
        GST(1) = AM(1)*AM(1)/RSTRETCH_F**TWO
        GST(2) = AM(2)*AM(2)/RSTRETCH_F**TWO
        GST(3) = AM(3)*AM(3)/RSTRETCH_F**TWO
        GST(4) = AM(1)*AM(2)/RSTRETCH_F**TWO
        IF (NSHR.EQ.3) THEN
          GST(5) = AM(1)*AM(3)/RSTRETCH_F**TWO
          GST(6) = AM(2)*AM(3)/RSTRETCH_F**TWO
        END IF
C
C FOURTH INVARIANT FOR THE STRETCH OF THE FIBER. ELASTIC
        AINV4 = GST(1) + GST(2) + GST(3)
C
C FIBER STRAIN
        AINV41=AINV4-ONE
        AINV412=AINV41/TWO
        AINV414=AINV41/FOUR
        EXPO=DEXP(CTE2*AINV41**TWO)
        CTE3=ONE+TWO*CTE2*AINV41**TWO
        IF (AINV41.LT.ZERO) THEN
          EFIB=ZERO
        ELSE
          EFIB=TWO*FRAC(K3)*CTE1*EXPO
        END IF
        EFIB2=EFIB*TWO
C
C STRESS IN THE FIBER
        DO K1=1,NTENS
          STRESS_FIB(K1)=STRESS_FIB(K1) + EFIB*AINV41*GST(K1)
        END DO
C
C ELASTICITY IN THE FIBER
        DSDE_FIB(1,1)=DSDE_FIB(1,1)+EFIB2*(CTE3*GST(1)+AINV41)*GST(1)
        DSDE_FIB(2,2)=DSDE_FIB(2,2)+EFIB2*(CTE3*GST(2)+AINV41)*GST(2)
        DSDE_FIB(3,3)=DSDE_FIB(3,3)+EFIB2*(CTE3*GST(3)+AINV41)*GST(3)
        DSDE_FIB(1,2)=DSDE_FIB(1,2)+EFIB2*CTE3*GST(1)*GST(2)
        DSDE_FIB(1,3)=DSDE_FIB(1,3)+EFIB2*CTE3*GST(1)*GST(3)
        DSDE_FIB(2,3)=DSDE_FIB(2,3)+EFIB2*CTE3*GST(2)*GST(3)
        DSDE_FIB(1,4)=DSDE_FIB(1,4)+EFIB2*(CTE3*GST(1)+AINV412)*GST(4)
        DSDE_FIB(2,4)=DSDE_FIB(2,4)+EFIB2*(CTE3*GST(2)+AINV412)*GST(4)
        DSDE_FIB(3,4)=DSDE_FIB(3,4)+EFIB2*CTE3*GST(3)*GST(4)
        DSDE_FIB(4,4)=DSDE_FIB(4,4)+EFIB2*(CTE3*GST(4)*GST(4)
     1                     +AINV414*(GST(1)+GST(2)))
        IF (NSHR.EQ.3) THEN
          DSDE_FIB(1,5)=DSDE_FIB(1,5)
     1                  +EFIB2*(CTE3*GST(1)+AINV412)*GST(5)
          DSDE_FIB(2,5)=DSDE_FIB(2,5)+EFIB2*CTE3*GST(2)*GST(5)
          DSDE_FIB(3,5)=DSDE_FIB(3,5)
     1                  +EFIB2*(CTE3*GST(3)+AINV412)*GST(5)
          DSDE_FIB(4,5)=DSDE_FIB(4,5)
     1              +EFIB2*(CTE3*GST(4)*GST(5)+AINV414*GST(6))
          DSDE_FIB(5,5)=DSDE_FIB(5,5)+EFIB2*(CTE3*GST(5)*GST(5)
     1              +AINV414*(GST(1)+GST(3)))
          DSDE_FIB(1,6)=DSDE_FIB(1,6)+EFIB2*CTE3*GST(1)*GST(6)
          DSDE_FIB(2,6)=DSDE_FIB(2,6)
     1                  +EFIB2*(CTE3*GST(2)+AINV412)*GST(6)
          DSDE_FIB(3,6)=DSDE_FIB(3,6)
     1                  +EFIB2*(CTE3*GST(3)+AINV412)*GST(6)
          DSDE_FIB(4,6)=DSDE_FIB(4,6)
     1              +EFIB2*(CTE3*GST(4)*GST(6)+AINV414*GST(5))
          DSDE_FIB(5,6)=DSDE_FIB(5,6)
     1              +EFIB2*(CTE3*GST(5)*GST(6)+AINV414*GST(4))
          DSDE_FIB(6,6)=DSDE_FIB(6,6)+EFIB2*(CTE3*GST(6)*GST(6)
     1              +AINV414*(GST(2)+GST(3)))
        END IF
C
C STRAIN ENERGY FUNCTION FOR FIBRES
      SEF_FIB=FRAC(K3)*HALF*(CTE1/CTE2)*(EXPO-ONE)
C STIFFNESS MAGNITUDE FOR REMODELING FLOW
        STIFF=STIFF+FOUR*FRAC(K3)*CTE1*(ONE-TWO*AINV4
     1      -TWO*CTE2*AINV4*AINV41**TWO)*EXPO*AINV4/RSTRETCH_F   !this it is with -> *lambda**2
C
      END DO
C
C--------------------------------------------------------------------
C TOTAL STRAIN ENERGY FUNCTION
      SEF_TOT=SEF_VOL+SEF_ICH+SEF_FIB
C
C--------------------------------------------------------------------
      RETURN
      END SUBROUTINE ELASTIC_LAW
C
C====================================================================
C COMPUTE THE INVERSE OF A MATRIX USING GAUSS-JORDAN ELIMINATION
C====================================================================
      SUBROUTINE XMATRIXINV(AMAT,BMAT,NUM)
C THE INVERSE OF MATRIX A(N,N) IS CALCULATED AND STORED IN THE MATRIX B(N,N)
      USE VARIABLES
      INCLUDE 'ABA_PARAM.INC'
C
C--------------------------------------------------------------------
C GLOBAL ARGUMENTS
C--------------------------------------------------------------------
C VARIABLES PASSED IN FOR INFORMATION
      DIMENSION AMAT(NUM,NUM), BMAT(NUM,NUM)
C
C--------------------------------------------------------------------
C LOCAL ARGUMENTS
C--------------------------------------------------------------------
C
C--------------------------------------------------------------------
C BUILD THE IDENTITY MATRIX
      DO I1=1,NUM
        DO J1=1,NUM
          BMAT(I1,J1) = ZERO
        END DO
  	    BMAT(I1,I1) = ONE
      END DO
C
C THIS IS THE BIG LOOP OVER ALL THE COLUMNS OF A(N,N)
      DO I1=1,NUM
C IN CASE THE ENTRY A(I,I) IS ZERO, WE NEED TO FIND A GOOD PIVOT; THIS PIVOT
C IS CHOSEN AS THE LARGEST VALUE ON THE COLUMN I FROM A(J,I) WITH J = 1,N
        BIG=AMAT(I1,I1)
        DO J1=I1,NUM
          IF (AMAT(J1,I1).GT.BIG) THEN
            BIG = AMAT(J1,I1)
            IROW=J1
          END IF
        END DO
C INTERCHANGE LINES I WITH IROW FOR BOTH A() AND B() MATRICES
        IF (BIG.GT.AMAT(I1,I1)) THEN
          DO K1=1,NUM
            DUM=AMAT(I1,K1) ! MATRIX A()
            AMAT(I1,K1)=AMAT(IROW,K1)
            AMAT(IROW,K1)=DUM
            DUM=BMAT(I1,K1) ! MATRIX B()
            BMAT(I1,K1)=BMAT(IROW,K1)
            BMAT(IROW,K1)=DUM
          ENDDO
        ENDIF
C DIVIDE ALL ENTRIES IN LINE I FROM A(I,J) BY THE VALUE A(I,I);
C SAME OPERATION FOR THE IDENTITY MATRIX
        DUM=AMAT(I1,I1)
        DO J1=1,NUM
          AMAT(I1,J1)=AMAT(I1,J1)/DUM
          BMAT(I1,J1)=BMAT(I1,J1)/DUM
        ENDDO
C MAKE ZERO ALL ENTRIES IN THE COLUMN A(J,I); SAME OPERATION FOR INDENT()
        DO J1=I1+1,NUM
          DUM=AMAT(J1,I1)
          DO K1=1,NUM
            AMAT(J1,K1)=AMAT(J1,K1) - DUM*AMAT(I1,K1)
            BMAT(J1,K1)=BMAT(J1,K1) - DUM*BMAT(I1,K1)
          END DO
        END DO
      END DO
C
C SUBSTRACT APPROPIATE MULTIPLE OF ROW J FROM ROW J-1
      DO I1=1,NUM-1
        DO J1=I1+1,NUM
          DUM=AMAT(I1,J1)
          DO L1=1,NUM
            AMAT(I1,L1)=AMAT(I1,L1)-DUM*AMAT(J1,L1)
            BMAT(I1,L1)=BMAT(I1,L1)-DUM*BMAT(J1,L1)
          END DO
        END DO
	    END DO
C
      RETURN
      END SUBROUTINE XMATRIXINV
C
C====================================================================
C POLAR DECOMPOSITION OF F=UR (LEFT DECOMPOSITION).
C THIS DECOMPOSTIOIN OF F IS UNCONVENTIONAL AS ABAQUS DOUBLE ROTATE
C F WHEN GIVING IT TO UMAT.
C BOX 7.1 SIMO AND HUGHES BOOK
C====================================================================
      SUBROUTINE RIGHTROTATION(ROT,STRU,DFGRD)
C
      USE VARIABLES
      INCLUDE 'ABA_PARAM.INC'
C--------------------------------------------------------------------
C GLOBAL ARGUMENTS
C--------------------------------------------------------------------
C VARIABLES PASSED IN FOR INFORMATION
      DIMENSION ROT(3,3), STRU(3,3), DFGRD(3,3)
C
C--------------------------------------------------------------------
C LOCAL ARGUMENTS
C--------------------------------------------------------------------
      DIMENSION ALAMDA(3), XA(3), BLCG(3,3), BLCG2(3,3),
     1  AID2(3,3), STRUINV(3,3)
C
      TOL = 1.0D-8
C--------------------------------------------------------------------
C
      DO K1=1,3
        DO K2=1,3
          AID2(K1,K2) = ZERO
        END DO
        AID2(K1,K1) = ONE
        XA(K1) = ZERO
      END DO
C
C COMPUTE THE LEFT CAUCHY-GREEN TENSOR, B=U^2
      BLCG = MATMUL(DFGRD,TRANSPOSE(DFGRD))
C
C COMPUTE THE INVARIANTS OF LEFT CAUCHY-GREEN TENSOR, B
      BLCG2 = MATMUL(BLCG,BLCG)
      BINV1 = BLCG(1,1) + BLCG(2,2) + BLCG(3,3)
      BINV2 = HALF*(BINV1**2-BLCG2(1,1)-BLCG2(2,2)-BLCG2(3,3))
      BINV3 = BLCG(1,1)*BLCG(2,2)*BLCG(3,3)
     1         -BLCG(1,2)*BLCG(2,1)*BLCG(3,3)
     2         +BLCG(1,2)*BLCG(2,3)*BLCG(3,1)
     3         +BLCG(1,3)*BLCG(3,2)*BLCG(2,1)
     4         -BLCG(1,3)*BLCG(3,1)*BLCG(2,2)
     5         -BLCG(2,3)*BLCG(3,2)*BLCG(1,1)
C
C--------------------------------------------------------------------
C
C COMPUTE THE SQUARES OF THE PRINCIPAL STRETCHES, ALAMDA**2
      CTB=(BINV1**TWO - THREE*BINV2)/9.0D0                              !Q
      CTC=(TWO*BINV1**THREE - 9.0D0*BINV1*BINV2 + 27.0D0*BINV3)/54.0D0  !R
C
      IF (DABS(CTB).LT.EPS) THEN
        CTB = ZERO
      END IF
      IF (DABS(CTC).LT.EPS) THEN
        CTC = ZERO
      END IF
C
C COMPUTE DISCRIMINANT. DELTA=CTB**3-CTC**2
      IF (CTB**THREE .GT. CTC**TWO) THEN
        ! Three real eigenvalues (Delta > 0)
        SQRTB = DSQRT(CTB)
        CTN = CTC / (SQRTB**THREE)
        CTT1 = DACOS(CTN) / THREE
c        CTT2 = DATAN2(DSQRT(ONE-CTN**TWO),CTN) / THREE
        DO K1=1,3
          XA(K1) = TWO*SQRTB*DCOS(CTT1+TWO*(REAL(K1,8)-ONE)*PI/THREE)
        END DO
      ELSE
        ! Two or more equal eigenvalues (Delta <= 0)
        DO K1=1,3
          XA(K1) = CUBIC_ROOT(CTC)
        END DO
      END IF
C
      DO K1=1,3
        ALAMDA(K1) = DSQRT(XA(K1) + BINV1/THREE)
      END DO
C
C COMPUTE THE STRETCH TENSOR, U
      UINV1 = ALAMDA(1) + ALAMDA(2) + ALAMDA(3)
      UINV2 = ALAMDA(1)*ALAMDA(2) + ALAMDA(1)*ALAMDA(3)
     1         + ALAMDA(2)*ALAMDA(3)
      UINV3 = ALAMDA(1)*ALAMDA(2)*ALAMDA(3)
C
      CTD = (ALAMDA(1)+ALAMDA(2))*(ALAMDA(1)+ALAMDA(3))
     1        *(ALAMDA(2)+ALAMDA(3))
C
      STRU = (-BLCG2+(UINV1**2-UINV2)*BLCG+UINV1*UINV3*AID2)/CTD
      STRUINV = (BLCG-UINV1*STRU+UINV2*AID2)/UINV3
C
C COMPUTE THE ROTATION TENSOR R FROM RIGHT STRETCH TENSOR U
      ROT = MATMUL(STRUINV,DFGRD)
C
      RETURN
      END SUBROUTINE
C
C====================================================================
C DEFINITION OF A CUBIC ROOT FUNCTION TO HANDLE NEGATIVE ARGUMENT
C====================================================================
      FUNCTION CUBIC_ROOT(X)
C
      REAL(KIND=8) X, CUBIC_ROOT
C
C HANDLES THE SIGN OF THE ARGUMENT SEPARATELY
      IF (X.GE.0.D0) THEN
        CUBIC_ROOT = X**(1.D0/3.D0)
      ELSE
        CUBIC_ROOT = -((-X)**(1.D0/3.D0))
      END IF
C
      RETURN
      END FUNCTION CUBIC_ROOT
C
C====================================================================
C NORM OF A SYMMETRIC SECOND-RANK TENSOR IN VOIGT NOTATION
C====================================================================
      FUNCTION SNORM2_VOIGT(TENSOR,NDI,NSHR,NTENS)
C
      INTEGER K1, NDI, NSHR, NTENS
      REAL(KIND=8) TENSOR(NTENS), NORM_VOIGT2, SNORM2_VOIGT
C
      NORM_VOIGT2 = 0.D0
      DO K1=1,NDI
        NORM_VOIGT2 = NORM_VOIGT2 + TENSOR(K1)*TENSOR(K1)
      END DO
      DO K1=NDI+1,NDI+NSHR
        NORM_VOIGT2 = NORM_VOIGT2 + TWO*TENSOR(K1)*TENSOR(K1)
      END DO
c      SNORM_VOIGT = DSQRT(NORM_VOIGT2)
      SNORM2_VOIGT = NORM_VOIGT2
C
      RETURN
      END FUNCTION SNORM2_VOIGT
C
C====================================================================