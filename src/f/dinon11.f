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
C     R = Ci.p/(1+(Ci.p/FU)^n)^(1/n)
C
C        f    : surface de charge
C        Ky   : raideur elastique
C        Fy   : limite elastique
C        Ci   : coefficient isotrope
C        n    : exposant isotrope
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
C     3 N      : EXPOSANT ISOTROPE                           n
C     4 FU     : LIMITE ULTIME                               Fu
C
C VARI:
C  1~ 6 UPL    : DEFORMATION PLASTIQUE                       Upl    .
C  7~12 VUPL   : VITESSE DE EPL                              dUpl = Upl
C    13 P      : DEFORMATION PLASTIQUE CUMULEE               p    .
C    14 VP     : VITESSE DE P                                dp = p
C
C***************** DECLARATION DES VARIABLES LOCALES *******************
C
      INTEGER II,ITER,JJ,ITER2
      REAL*8  R8MIEM,DNRM2
      REAL*8  R8MIN,ERRE,ERRE2,OLDERRE2
      REAL*8  ULEL(NBCOMP),DULEL(NBCOMP),UTLEL(NBCOMP)
      REAL*8  NORM,TMP,VTMP(NBCOMP)
      REAL*8  KY,FY,CI,N,FU
      REAL*8  FM(NBCOMP),FP(NBCOMP),RM(NBCOMP),RP(NBCOMP)
      REAL*8  A0,A1,A2,DLN,DLNN,FDLN,FPDLN,DLNNEG,DLNPOS,FDLNNEG,FDLNPOS
      REAL*8  FK(NBCOMP),FKOLD(NBCOMP),FFK,FFKOLD,ETA,FKEST(NBCOMP)
      REAL*8  AK(NBCOMP),BK(NBCOMP),DEN,BETAJ,BETAJOLD
C
C************ FIN DES DECLARATIONS DES VARIABLES LOCALES ***************

      REAL*8  DDOT
      REAL*8  DINON11_FDL
      REAL*8  DINON11_FPDL

C ----------------------------------------------------------------------
      R8MIN  = R8MIEM()

C
C     1. VALEURS PAR DEFAUT (COMPORTEMENT ELASTIQUE)
C

C     PAR DEFAUT LES VARIABLES N'EVOLUENT PAS
      CALL DCOPY(NBVARI, VARIMO, 1, VARIPL, 1)
      CALL DSCAL(NBCOMP, 0.0D0, VARIPL(7), 1)
      VARIPL(14) = 0.0D0
!       WRITE(6,*) "NBCOMP = ",NBCOMP
!       WRITE(6,*) "RAIDE  = ",RAIDE
!       WRITE(6,*) "Upl-   = ",VARIMO(1:2)
!       WRITE(6,*) "dUpl-  = ",VARIMO(7:8)
!       WRITE(6,*) "p-     = ",VARIMO(13)
!       WRITE(6,*) "dp-    = ",VARIMO(14)


C
C     2. CALCUL DES DEPLACEMENTS
C

C     CONSTRUCTION DES VECTEURS DEPLACEMENT ELEMENTAIRE
!       WRITE(6,*) "DUL   = ",DUL(7:8)
!       WRITE(6,*) "UL    = ",UL(7:8)
!       WRITE(6,*) "UTL   = ",UTL(7:8)
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
!       WRITE(6,*) "DULEL = ",DULEL(1:2)
!       WRITE(6,*) "ULEL  = ",ULEL(1:2)
!       WRITE(6,*) "UTLEL = ",UTLEL(1:2)

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
C     EXPOSANT ISOTROPE
      N    = PARAM(1,3)
C     LIMITE ULTIME
      FU   = PARAM(1,4)
!       WRITE(6,*) "KY    = ",KY
!       WRITE(6,*) "FY    = ",FY
!       WRITE(6,*) "CI    = ",CI
!       WRITE(6,*) "N     = ",N
!       WRITE(6,*) "FU    = ",FU


C
C     4. RECALCUL DES QUANTITES PRECENTES
C

C     VECTEUR FORCE F- = Ky * ( U- - Upl- )
      CALL DCOPY(NBCOMP, VARIMO(1), 1, FM, 1)
      CALL DSCAL(NBCOMP, -KY, FM, 1)
      CALL DAXPY(NBCOMP, KY, ULEL, 1, FM, 1)
!       WRITE(6,*) "F-    = ",FM(1:2)

C     RAYON DE LA SURFACE DE CHARGE R-
C     R = Ci.p/(1+(Ci.p/FU)^n)^(1/n)
      RM(1) = CI * VARIMO(13)
      RM(1) = RM(1) / ((1 + ((RM(1) / FU) ** N)) ** (1.0D0 / N))
!       WRITE(6,*) "R-    = ",RM(1)


C
C     5. PREDICTION ELASTIQUE DE F+
C

C     R+ = R- (PREDICTION ELASTIQUE INITIALE)
      CALL DCOPY(NBCOMP, RM, 1, RP, 1)

! C     Upl+ = Upl- (PREDICTION ELASTIQUE INITIALE)
!       CALL DCOPY(NBCOMP, VARIMO(1), 1, VARIPL(1), 1)

C     VECTEUR FORCE F+ = F- + dF
C                   dF = Ky . dU               dSigma = E : dEpsilon <=> dF = E : dUel
C                   F- = Ky * ( U- - Upl- )    dU = dUel + dUpl
C                   Upl+ = Upl-                par hypothese de prediction elastique
C               =>  F+ = Ky * ( U+ - Upl- )
      CALL DCOPY(NBCOMP, VARIMO(1), 1, FP, 1)
      CALL DSCAL(NBCOMP, -KY, FP, 1)
      CALL DAXPY(NBCOMP, KY, UTLEL, 1, FP, 1)
!       WRITE(6,*) "F+    = ",FP(1:2)," (prediction elastique)"


C
C     6. TEST DE LA FONCTION SEUIL
C

C     VALEUR DE LA FONCTION DE SURFACE DE CHARGE
C     f = |F+| - R+ - Fy
      NORM = DNRM2(NBCOMP, FP, 1) - RP(1)
!       WRITE(6,*) "Prediction elastique (R+ = R-):"
!       WRITE(6,*) "||F+|| - (R-) = ",NORM
!       WRITE(6,*) "FY    = ",FY

C     TEST DE COMPORTEMENT ELASTIQUE
      IF ( NORM .LE. FY ) THEN
         WRITE(6,*) "ELASTIC BEHAVIOR -> END"
         RETURN
      ENDIF
!       WRITE(6,*) "||F+|| - (R-) > Fy : PLASTIC BEHAVIOR"


C
C     7. INITIALISATION DES VARIABLES DE BOUCLE
C

C     ERREUR
      ERRE  = 1.0D6

C     NUMERO D'ITERATION ITER=k
      ITER  = 0

C     FORCE FK=F(k) ET FKOLD=F(k-1)
      CALL DCOPY(NBCOMP, FP, 1, FK, 1)
      CALL DCOPY(NBCOMP, FM, 1, FKOLD, 1)

C     B0 = df/dF | (F=F-) = F- / ||F-||
      CALL DCOPY(NBCOMP, FM, 1, BK, 1)
      NORM = DNRM2(NBCOMP, FM, 1)
      CALL DSCAL(NBCOMP, 1.0D0 / NORM, BK, 1)
!       WRITE(6,*) "B0    = ",BK(1:2)

C     A0 = Ky . df/dF | (F=F-) = Ky . F- / ||F-|| = Ky . B0
      CALL DCOPY(NBCOMP, BK, 1, AK, 1)
      CALL DSCAL(NBCOMP, KY, AK, 1)
!       WRITE(6,*) "A0    = ",AK(1:2)


C
C     8. BOUCLE
C

10    IF ( ERRE .GT. 1.0D-3 ) THEN

         ITER = ITER + 1
         WRITE(6,*) "====== ITERATION (1) #",ITER
         IF ( ITER .GT. 100 ) THEN
!             CALL U2MESx
            WRITE(6,*) "@@@@@@@@@@ TOO MANY ITERATIONS (1) @@@@@@@@@@"
            RETURN
         ENDIF


         ITER2 = 0
         ERRE2 = 1.0D6
         OLDERRE2 = ERRE2
15       IF ( ERRE2 .GT. 1.0D-1 ) THEN

            ITER2 = ITER2 + 1
            WRITE(6,*) "------ ITERATION (2) #",ITER2
            IF ( ITER2 .GT. 30 ) THEN
!                CALL U2MESx
               WRITE(6,*) "@@@@@@@@@@ TOO MANY ITERATIONS (2) @@@@@@@@@"
               RETURN
            ENDIF


C
C           9. COMPUTE THE FRACTION OUTSIDE THE YIELD SURFACE
C

C           f(Fk-1) = ||Fk-1|| - R+ - Fy
            WRITE(6,*) "Fk<   = ",FKOLD(1:2)
            FFKOLD = DNRM2(NBCOMP, FKOLD, 1) - RP(1) - FY
            WRITE(6,*) "f(Fk<)= ",FFKOLD

C           f(Fk) = ||Fk|| - R+ - Fy
            WRITE(6,*) "Fk>   = ",FK(1:2)
            FFK = DNRM2(NBCOMP, FK, 1) - RP(1) - FY
            WRITE(6,*) "f(Fk>)= ",FFK
            
C           FRACTION eta = 1 / [ 1 - f(Fk-1) / f(Fk) ]
            ETA = 1.0D0 / (1.0D0 - FFKOLD / FFK)
            WRITE(6,*) "ETA   = ",ETA

C           Fk(estimation) = eta * Fk-1 + (1 - eta) * Fk
            CALL DCOPY(NBCOMP, FKOLD, 1, FKEST, 1)
            CALL DSCAL(NBCOMP, ETA, FKEST, 1)
            CALL DAXPY(NBCOMP, (1.0D0 - ETA), FK, 1, FKEST, 1)
            WRITE(6,*) "Fk*   = ",FKEST(1:2)


C
C           9. UPDATE INTERNAL STATE VARIABLES
C

C           CALCUL DU MULTIPLICATEUR PLASTIQUE dL
C           dL = (B : dF) / [A : B - G : M] | F=Ftrial
C             dF = Ky . dU
C             - G : M = -mk M : M = -mk = - df/dR = +1
C        => dL = Ky . (B : dU) / (A : B + 1)
C             B = F / ||F|| => Ky . (B : dU) = F : Ky.dU / ||F||
C             A : B = Ky . B : B = Ky car B unitaire dans notre cas
C        => dL = F : dU / ||F|| * Ky / (Ky + 1)
C
C           BRUTE FORCE WAY...
            CALL DCOPY(NBCOMP, FP, 1, VTMP, 1)
            CALL DAXPY(NBCOMP, -1.0D0, FKEST, 1, VTMP, 1)
            NORM = DDOT(NBCOMP, BK, 1, VTMP, 1)
            VARIPL(14) = NORM / (DDOT(NBCOMP, AK, 1, BK, 1) + 1.0D0)
            WRITE(6,*) "dL    = ",VARIPL(14)," (= dp)"

C           CALCUL DE LA DEFORMATION PLASTIQUE CUMULEE p
C           p+ = p- + dp = p- + dLambda
            VARIPL(13) = VARIMO(13) + VARIPL(14)
            WRITE(6,*) "p+    = ",VARIPL(13)

C           RAYON DE LA SURFACE DE CHARGE R+
C           R = Ci.p/(1+(Ci.p/FU)^n)^(1/n)
            RP(1) = CI * VARIPL(13)
            RP(1) = RP(1) / ((1 + ((RP(1) / FU) ** N)) ** (1.0D0 / N))
            WRITE(6,*) "R+    = ",RP(1)


C
C           11. UPDATE THE YIELD SURFACE
C

C           Bk = df/dF | (F=Fk) = Fest / ||Fest||
            CALL DCOPY(NBCOMP, FKEST, 1, BK, 1)
            NORM = DNRM2(NBCOMP, FKEST, 1)
            CALL DSCAL(NBCOMP, 1.0D0 / NORM, BK, 1)
            WRITE(6,*) "Bk    = ",BK(1:2)

C           Ak = Ky . df/dF | (F=Fk) = Ky . Fest / ||Fest|| = Ky . Bk
            CALL DCOPY(NBCOMP, BK, 1, AK, 1)
            CALL DSCAL(NBCOMP, KY, AK, 1)
            WRITE(6,*) "Ak    = ",AK(1:2)


C           COMPUTE THE ERROR
            OLDERRE2 = ERRE2
            ERRE2 = ABS(DNRM2(NBCOMP, FKEST, 1) - RP(1) - FY)
            WRITE(6,*) "FKest = ",FKEST
            WRITE(6,*) "RP    = ",RP(1)
            WRITE(6,*) "FY    = ",FY
            WRITE(6,*) "erre2 = ",ERRE2

C           CHECK CONVERGENCE
            IF ( ABS((OLDERRE2 - ERRE2) / OLDERRE2) .LT. 0.2D0 ) THEN

C              LINEAR SEARCH
               WRITE(6,*) "*** LINEAR SEARCH ***"
               DO 17, JJ=1,6

C                 ETA
                  ETA = 1.0D0 / 7.0D0 * JJ
                  WRITE(6,*) "ETA   = ",ETA

C                 Fk(estimation) = eta * Fk-1 + (1 - eta) * Fk
                  CALL DCOPY(NBCOMP, FKOLD, 1, FKEST, 1)
                  CALL DSCAL(NBCOMP, ETA, FKEST, 1)
                  CALL DAXPY(NBCOMP, (1.0D0 - ETA), FK, 1, FKEST, 1)
                  WRITE(6,*) "Fk*   = ",FKEST(1:2)

C                 BRUTE FORCE WAY...
                  CALL DCOPY(NBCOMP, FP, 1, VTMP, 1)
                  CALL DAXPY(NBCOMP, -1.0D0, FKEST, 1, VTMP, 1)
                  NORM = DDOT(NBCOMP, BK, 1, VTMP, 1)
                  VARIPL(14) = NORM / (DDOT(NBCOMP, AK, 1, BK, 1) + 1.0D0)
                  WRITE(6,*) "dL    = ",VARIPL(14)," (= dp)"

C                 CALCUL DE LA DEFORMATION PLASTIQUE CUMULEE p
C                 p+ = p- + dp = p- + dLambda
                  VARIPL(13) = VARIMO(13) + VARIPL(14)
                  WRITE(6,*) "p+    = ",VARIPL(13)

C                 RAYON DE LA SURFACE DE CHARGE R+
C                 R = Ci.p/(1+(Ci.p/FU)^n)^(1/n)
                  RP(1) = CI * VARIPL(13)
                  RP(1) = RP(1) / ((1 + ((RP(1) / FU) ** N)) ** (1.0D0 / N))
                  WRITE(6,*) "R+    = ",RP(1)

C                 Bk = df/dF | (F=Fk) = Fest / ||Fest||
                  CALL DCOPY(NBCOMP, FKEST, 1, BK, 1)
                  NORM = DNRM2(NBCOMP, FKEST, 1)
                  CALL DSCAL(NBCOMP, 1.0D0 / NORM, BK, 1)
                  WRITE(6,*) "Bk    = ",BK(1:2)

C                 Ak = Ky . df/dF | (F=Fk) = Ky . Fest / ||Fest|| = Ky . Bk
                  CALL DCOPY(NBCOMP, BK, 1, AK, 1)
                  CALL DSCAL(NBCOMP, KY, AK, 1)
                  WRITE(6,*) "Ak    = ",AK(1:2)

C                 COMPUTE THE ERROR
!                   OLDERRE2 = ERRE2
                  ERRE2 = ABS(DNRM2(NBCOMP, FKEST, 1) - RP(1) - FY)
                  WRITE(6,*) "FKest = ",FKEST
                  WRITE(6,*) "RP    = ",RP(1)
                  WRITE(6,*) "FY    = ",FY
                  WRITE(6,*) "erre2 = ",ERRE2

C              LOOP ON LINEAR SEARCH
17             CONTINUE
               WRITE(6,*) "*** END OF LINEAR SEARCH ***"


            ENDIF

C           UPDATE QUANTITIES
            IF ( DNRM2(NBCOMP, FKEST, 1) - RP(1) - FY .LE. 0.D0 ) THEN
               CALL DCOPY(NBCOMP, FKEST, 1, FKOLD, 1)
            ELSE
               CALL DCOPY(NBCOMP, FKEST, 1, FK, 1)
            ENDIF

            GOTO 15
         ENDIF


C
C        12. UPDATE THE STRESS ESTIMATE BY PROJECTING THE TRIAL
C            ONTO THE NEW YIELD SURFACE (AS COMPUTED IN STEP 10b)
C

C        Fk+1 = Fk - eta * [Ak * (Bk : (Fk - Fk-1))] / [Ak : Bk]

C        CALCUL DU DENOMINATEUR
         DEN = DDOT(NBCOMP, AK, 1, BK, 1)
         WRITE(6,*) "DEN   = ",DEN

C        CALCUL DU SCALAIRE DU NUMERATEUR [ Bk : (Fk - Fk-1) ] * eta
         CALL DCOPY(NBCOMP, FK, 1, VTMP, 1)
         CALL DAXPY(NBCOMP, -1.0D0, FKOLD, 1, VTMP, 1)
         NORM = DDOT(NBCOMP, BK, 1, VTMP, 1) * ETA
         WRITE(6,*) "NUM   = ",NORM
         
C        CALCUL DE Fk+1 = Fk - Ak . NORM / DEN
         CALL DCOPY(NBCOMP, FK, 1, FKOLD, 1)
         CALL DAXPY(NBCOMP, - NORM / DEN, AK, 1, FK, 1)
         WRITE(6,*) "Fk+1  = ",FK

C        CALCUL DE L'INCREMENT DE DEPLACEMENT PLASTIQUE
C        dUpl = dLambda.df/dF = dLambda.F+/||F+|| (normale a 1 sphere)
         CALL DCOPY(NBCOMP, FK, 1, VARIPL(7), 1)
         NORM = DNRM2(NBCOMP, VARIPL(7), 1)
         CALL DSCAL(NBCOMP, VARIPL(14) / NORM, VARIPL(7), 1)
         WRITE(6,*) "dUpl  = ",VARIPL(7:8)

C        CALCUL DU DEPLACEMENT PLASTIQUE
C        Upl+ = Upl- + dUpl
         CALL DCOPY(NBCOMP, VARIPL(7), 1, VARIPL(1), 1)
         CALL DAXPY(NBCOMP, 1.0D0, VARIMO(1), 1, VARIPL(1), 1)
         WRITE(6,*) "Upl+  = ",VARIPL(1:2)


C
C        13. TEST DE LA FONCTION SEUIL
C

C        VALEUR DE LA FONCTION DE SURFACE DE CHARGE
C        f = |F+| - R+ - Fy
         ERRE = DNRM2(NBCOMP, FK, 1) - RP(1) - FY
         WRITE(6,*) "f(Fk+1) = ",ERRE
         ERRE = ABS(ERRE)

C        TEST DE COMPORTEMENT ELASTIQUE
         IF ( ERRE .LE. R8MIN ) THEN
            ERRE = 0.0D0
            GOTO 10
         ENDIF

C        BOUCLE
         GOTO 10

      ENDIF


C     COPIE DE LA SOLUTION FINALE EN FORCE
      CALL DCOPY(NBCOMP, FK, 1, FP, 1)

C     CALCUL DE LA RAIDEUR TANGENTE
C     Kt = dF / dU
      DO 20, II=1,NBCOMP
         IF ( ABS(DULEL(II)) .GT. R8MIN ) THEN
            RAIDE(II) = ABS((FP(II) - FM(II)) / DULEL(II))
         ENDIF
20    CONTINUE
!       WRITE(6,*) "RAIDE = ",RAIDE(1:2)

      END








! C ----------------------------------------------------------------------
! C
! C     CALCUL DE f(dLambda) POUR LA RESOLUTION NUMERIQUE DE f(dLambda)=0
! C
! C     f(x) = Ci.dL + (Cst + Ky.dL).[1+(Ci.dL/Fu)^n]^(1/n)
! C
! C     avec Cst = Fy - || Ky.dU + F- ||
! C
! C ----------------------------------------------------------------------
!       FUNCTION DINON11_FDL(DL,CST,KY,FY,CI,N,FU,R8MIN)
!       IMPLICIT NONE
!       REAL*8 DINON11_FDL
!       REAL*8 DL,CST,KY,FY,CI,N,FU,R8MIN
!       REAL*8 TMP
! 
!       IF ( ABS(DL) .LE. R8MIN ) THEN
!          DINON11_FDL = CST
!       ELSE
!          TMP = 1.0D0 + ((CI * DL / FU) ** N)
!          DINON11_FDL = CI * DL + (CST + KY * DL) * (TMP ** (1.0D0 / N))
!       ENDIF
! 
!       END
! 
! C ----------------------------------------------------------------------
! C
! C     CALCUL DE f'(dLambda) POUR LA RESOLUTION NUMERIQUE DE f(dLambda)=0
! C
! C     f'(x) = Ci + Ky.[1+(Ci.dL/Fu)^n]^(1/n)
! C           + (Cst + Ky.dL).[1+(Ci.dL/Fu)^n]^((1/n)-1).(Ci/Fu)^n.x^(n-1)
! C
! C     avec Cst = Fy - || Ky.dU + F- ||
! C
! C ----------------------------------------------------------------------
!       FUNCTION DINON11_FPDL(DL,CST,KY,FY,CI,N,FU,R8MIN)
!       IMPLICIT NONE
!       REAL*8 DINON11_FPDL
!       REAL*8 DL,CST,KY,FY,CI,N,FU,R8MIN
!       REAL*8 TMP,TMP2,TMP3
! 
!       IF ( ABS(DL) .LE. R8MIN ) THEN
!          DINON11_FPDL = CI + KY
!       ELSE
!          TMP = 1.0D0 + ((CI / FU * DL) ** N)
!          TMP2 = CI + KY * (TMP ** (1.0D0 / N))
!          TMP3 = (TMP ** (1.0D0 / N - 1.0D0)) * (DL ** (N - 1.0D0))
!          DINON11_FPDL = TMP2 + (CST + KY * DL) * TMP3 * ((CI / FU) ** N)
!       ENDIF
! 
!       END
