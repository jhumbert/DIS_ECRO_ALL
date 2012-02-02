      SUBROUTINE DINON11(NEQ,UL,DUL,UTL,NNO,
     &                  NBCOMP,NBVARI,VARIMO,RAIDE,NBPAR,PARAM,
     &                  VARIPL)
C ----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  NEQ,NBCOMP,NNO,NBPAR,NBVARI
      REAL*8   UL(NEQ),DUL(NEQ),UTL(NEQ)
      REAL*8   VARIMO(NBVARI),VARIPL(NBVARI)
      REAL*8   RAIDE(NBCOMP),PARAM(6,NBPAR)
C ----------------------------------------------------------------------
C            CONFIGURATION MANAGEMENT OF EDF VERSION
C MODIF ELEMENTS  DATE 05/02/2008   AUTEUR FLEJOU J-L.FLEJOU 
C ======================================================================
C COPYRIGHT (C) 1991 - 2006  EDF R&D                  WWW.CODE-ASTER.ORG
C THIS PROGRAM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C (AT YOUR OPTION) ANY LATER VERSION.
C
C THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT
C WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. SEE THE GNU
C GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C ALONG WITH THIS PROGRAM; IF NOT, WRITE TO EDF R&D CODE_ASTER,
C   1 AVENUE DU GENERAL DE GAULLE, 92141 CLAMART CEDEX, FRANCE.
C ======================================================================

C ======================================================================
C
C     DIS_ECRO_ISO3D
C
C     RELATION DE COMPORTEMENT ISOTROPE 3D EXACTE (RETURN MAPPING)
C
C     f = |F| - R - Fy
C     .      .         .   .
C     F = Ky.Uel = Ky.(U - Upl)
C
C     R = Ci.p
C
C        f    : surface de charge
C        Ky   : raideur elastique
C        Fy   : limite elastique
C        Ci   : coefficient isotrope
C        Uel  : deplacement relatif elastique
C        Upl  : deplacement relatif plastique
C        U    : deplacement relatif total U = Uel + Upl
C        R    : rayon de la surface de charge
C        p    : deformation plastique cumulee p = SUM ( Upl )
C
C======================================================================
C
C IN  :
C       NEQ    : NOMBRE DE DDL DE L'ELEMENT
C       UL     : DEPLACEMENT PRECEDENT REPERE LOCAL (DIM NEQ)
C       DUL    : INCREMENT DE DEPLACEMENT REPERE LOCAL (DIM NEQ)
C       UTL    : DEPLACEMENT COURANT REPERE LOCAL (DIM NEQ)
C       NNO    : NOMBRE DE NOEUDS
C       NBCOMP : NOMBRE DE COMPOSANTES
C       NBVARI : NOMBRE DE VARIABLES INTERNES
C       VARIMO : VARIABLES INTERNES A T- (3 PAR COMPOSANTES)
C       RAIDE  : RAIDEUR ELASTIQUE DES DISCRETS
C       NBPAR  : NOMBRE DE PARAMETRE MAXIMUM TOUTES LOIS CONFONDUES
C       PARAM  : PARAMETRES DE LA LOI (DIM = 6,NBPAR)
C
C OUT :
C       RAIDE  : RAIDEUR QUASI-TANGENTE AU COMPORTEMENT DES DISCRETS
C       VARIPL : VARIABLES INTERNES INTERNES A T+ (3 PAR COMPOSANTES)
C
C PARA:
C     1 FY     : LIMITE ELASTIQUE                            Fy
C     2 CI     : COEFFICIENT ISOTROPE                        Ci
C
C VARI:
C  1~ 6 UPL    : DEFORMATION PLASTIQUE                       Upl    .
C  7~12 VUPL   : VITESSE DE EPL                              dUpl = Upl
C    13 P      : DEFORMATION PLASTIQUE CUMULEE               p    .
C    14 VP     : VITESSE DE P                                dp = p
C
C***************** DECLARATION DES VARIABLES LOCALES *******************
C
      INTEGER II
      REAL*8  R8MIEM,DNRM2
      REAL*8  R8MIN
      REAL*8  ULEL(NBCOMP),DULEL(NBCOMP),UTLEL(NBCOMP)
      REAL*8  NORM,TMP,VTMP(NBCOMP)
      REAL*8  KY,FY,CI
      REAL*8  FM(NBCOMP),FP(NBCOMP),RM(NBCOMP),RP(NBCOMP)
C
C************ FIN DES DECLARATIONS DES VARIABLES LOCALES ***************

C ----------------------------------------------------------------------
      R8MIN  = R8MIEM()

C
C     1. VALEURS PAR DEFAUT (COMPORTEMENT ELASTIQUE)
C

C     PAR DEFAUT LES VARIABLES N'EVOLUENT PAS
      CALL DCOPY(NBVARI, VARIMO, 1, VARIPL, 1)
      CALL DSCAL(NBCOMP, 0.0D0, VARIPL(7), 1)
      VARIPL(14) = 0.0D0
      WRITE(6,*) "NBCOMP = ",NBCOMP
      WRITE(6,*) "RAIDE  = ",RAIDE
      WRITE(6,*) "Upl-   = ",VARIMO(1:2)
      WRITE(6,*) "dUpl-  = ",VARIMO(7:8)
      WRITE(6,*) "p-     = ",VARIMO(13)
      WRITE(6,*) "dp-    = ",VARIMO(14)


C
C     2. CALCUL DES DEPLACEMENTS
C

C     CONSTRUCTION DES VECTEURS DEPLACEMENT ELEMENTAIRE
      WRITE(6,*) "DUL   = ",DUL(7:8)
      WRITE(6,*) "UL    = ",UL(7:8)
      WRITE(6,*) "UTL   = ",UTL(7:8)
      IF ( NNO .EQ. 1 ) THEN
C       ULEL = UL
        CALL DCOPY(NBCOMP, DUL, 1, DULEL, 1)
        CALL DCOPY(NBCOMP, UL,  1, ULEL,  1)
        CALL DCOPY(NBCOMP, UTL, 1, UTLEL, 1)
      ELSE
C       ULEL = UL(n2)
        CALL DCOPY(NBCOMP, DUL(NBCOMP+1), 1, DULEL, 1)
        CALL DCOPY(NBCOMP, UL(NBCOMP+1),  1, ULEL,  1)
        CALL DCOPY(NBCOMP, UTL(NBCOMP+1), 1, UTLEL, 1)
C       ULEL = -1 * UL(n1) + ULEL = UL(n2) - UL(n1)
        CALL DAXPY(NBCOMP, -1.0D0, DUL, 1, DULEL, 1)
        CALL DAXPY(NBCOMP, -1.0D0, UL,  1, ULEL,  1)
        CALL DAXPY(NBCOMP, -1.0D0, UTL, 1, UTLEL, 1)
      ENDIF
      WRITE(6,*) "DULEL = ",DULEL(1:2)
      WRITE(6,*) "ULEL  = ",ULEL(1:2)
      WRITE(6,*) "UTLEL = ",UTLEL(1:2)

C     NORME DE DULEL
      NORM = DNRM2(NBCOMP, DULEL, 1)

C     TEST NO-OP ("NO OPERATION"; PAS DE DEPLACEMENT)
      IF ( NORM .LE. R8MIN ) THEN
C        NOTE: ON NE TESTE PAS RIGI_MECA_TANG
C        ICI ON SORT AVEC UNE MATRICE ELASTIQUE
C        MEME SI ON A PLASTIFIE AU PAS PRECEDENT
C        ON POURRAIT AMELIORER EN PASSANT DIRECTEMENT
C        LA MATRICE TANGENTE POUR GAGNER 1 ITERATION
         WRITE(6,*) "NO-OP -> END"
         RETURN
      ENDIF


C
C     3. RECUPERATION DES PARAMETRES DE LA LOI
C

C     RAIDEUR ELASTIQUE
      KY   = RAIDE(1)
C     LIMITE ELASTIQUE
      FY   = PARAM(1,1)
C     CONSTANTE ISOTROPE
      CI   = PARAM(1,2)
      WRITE(6,*) "KY    = ",KY
      WRITE(6,*) "FY    = ",FY
      WRITE(6,*) "CI    = ",CI


C
C     4. DETERMINATION DU COMPORTEMENT
C

C     VECTEUR FORCE F- = Ky * ( U- - Upl- )
      CALL DCOPY(NBCOMP, VARIMO(1), 1, FM, 1)
      CALL DSCAL(NBCOMP, -KY, FM, 1)
      CALL DAXPY(NBCOMP, KY, ULEL, 1, FM, 1)
      WRITE(6,*) "FM    = ",FM(1:2)

C     VECTEUR FORCE F+ = Ky * ( U+ - Upl+ )
C     PREDICTION ELASTIQUE (Upl+ = Upl-) => F+ = Ky * ( U+ - Upl- )
      CALL DCOPY(NBCOMP, VARIMO(1), 1, FP, 1)
      CALL DSCAL(NBCOMP, -KY, FP, 1)
      CALL DAXPY(NBCOMP, KY, UTLEL, 1, FP, 1)
      WRITE(6,*) "FP    = ",FP(1:2)," (prediction elastique)"

C     RAYON DE LA SURFACE DE CHARGE R-
!       CALL DCOPY(NBCOMP, VARIMO(13), 1, RM, 1)
!       CALL DSCAL(NBCOMP, CI, RM, 1)
      RM(1) = VARIMO(13) * CI
      WRITE(6,*) "RM    = ",RM(1)

C     VALEUR DE LA FONCTION DE SURFACE DE CHARGE
C     f = |F+| - R+ - Fy = |F+| - R- - Fy (approx)
      NORM = DNRM2(NBCOMP, FP, 1) - RM(1)
      WRITE(6,*) "NORM  = ",NORM
      WRITE(6,*) "FY    = ",FY

C     TEST DE COMPORTEMENT ELASTIQUE
      IF ( NORM .LE. FY ) THEN
         WRITE(6,*) "ELASTIC BEHAVIOR -> END"
         RETURN
      ENDIF
      WRITE(6,*) "||F|| - R > Fy : PLASTIC BEHAVIOR"


C
C     5. COMPORTEMENT PLASTIQUE
C

C     CALCUL DU MULTIPLICATEUR PLASTIQUE dL = dp
C     dLambda = (||Ky.dU + F-|| - Fy - Ci.p-) / (Ky + Ci)
      CALL DCOPY(NBCOMP, FM, 1, VTMP, 1)
      CALL DAXPY(NBCOMP, KY, DULEL, 1, VTMP, 1)
      NORM = DNRM2(NBCOMP, VTMP, 1)
      VARIPL(14) = (NORM - FY - CI * VARIMO(13)) / (KY + CI)
      WRITE(6,*) "dp    = ",VARIPL(14)

C     CALCUL DE L'INCREMENT DE DEPLACEMENT PLASTIQUE
C     dUpl = dLambda.df/dF = dLambda.F/||F|| (normale a 1 sphere)
      CALL DCOPY(NBCOMP, FM, 1, VARIPL(7), 1)
      NORM = DNRM2(NBCOMP, VARIPL(7), 1)
      CALL DSCAL(NBCOMP, VARIPL(14) / NORM, VARIPL(7), 1)
      WRITE(6,*) "dUpl  = ",VARIPL(7:8)

C     CALCUL DU DEPLACEMENT PLASTIQUE
C     Upl+ = Upl- + dUpl
      CALL DCOPY(NBCOMP, VARIPL(7), 1, VARIPL(1), 1)
      CALL DAXPY(NBCOMP, 1.0D0, VARIMO(1), 1, VARIPL(1), 1)
      WRITE(6,*) "Upl+  = ",VARIPL(1:2)

C     CALCUL DE LA DEFORMATION PLASTIQUE CUMULEE p
C     p+ = p- + dp = p- + dLambda
      CALL DAXPY(NBCOMP, 1.0D0, VARIPL(14), 1, VARIPL(13), 1)
      WRITE(6,*) "p+    = ",VARIPL(13)

C     CALCUL DE dF
C     dF = Ky.(dU - dUpl)

C     CALCUL DE F+
C     F+ = F- + dF = F- + Ky.dUel= F- + Ky.(dU - dUpl)
      CALL DCOPY(NBCOMP, FM, 1, FP, 1)
      CALL DCOPY(NBCOMP, DULEL, 1, VTMP, 1)
      CALL DAXPY(NBCOMP, -1.0D0, VARIPL(7), 1, VTMP, 1)
      CALL DAXPY(NBCOMP, KY, VTMP, 1, FP, 1)
      WRITE(6,*) "FP    = ",FP(1:2)

C     CALCUL DE LA RAIDEUR TANGENTE
C     Kt = dF / dU
      DO 20, II=1,NBCOMP
         IF ( ABS(DULEL(II)) .GT. R8MIN ) THEN
            RAIDE(II) = ABS((FP(II) - FM(II)) / DULEL(II))
         ENDIF
20    CONTINUE
      WRITE(6,*) "RAIDE = ",RAIDE(1:2)

      END
