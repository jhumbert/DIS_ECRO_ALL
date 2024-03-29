      SUBROUTINE DINONC3D(NOMTE,ICODRE,VALRE,KLV,RAIDE,NBPAR,PARAM,
     &                  NPLOI)
C ----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*16  NOMTE
      INTEGER   ICODRE(*)
      INTEGER       NBPAR, NPLOI
      REAL*8        VALRE(*),KLV(*),RAIDE(*),PARAM(6,NBPAR)

C ----------------------------------------------------------------------
C            CONFIGURATION MANAGEMENT OF EDF VERSION
C MODIF ELEMENTS  DATE 20/04/2011   AUTEUR COURTOIS M.COURTOIS 
C ======================================================================
C COPYRIGHT (C) 1991 - 2011  EDF R&D                  WWW.CODE-ASTER.ORG
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
C           AFFECTATION DES VALEURS ISSUES DU COMPORTEMENT
C
C   DANS CETTE VERSION 3D LES DIRECTIONS SONT COUPLEES DONC TOUTES
C   ACTIVES. PAS DE 'OKDIRE' NECESSAIRE DONC, NI DE RANGEMENT DES
C   PARAMETRES DANS UN ORDRE PARTICULIER DANS LE DATA 'NOMRE'.
C
C ======================================================================
C
C  IN
C     NOMTE : NOM DE L'ELEMENT
C     ICODRE: 0 SI LE COEFF EST PRESENT SINON 1
C     VALRE : VALEUR DES COEFFICIENTS
C     KLV   : RAIDEUR ELASTIQUE DU DISCRET
C     NBPAR : NOMBRE DE PARAMETRE MAXIMUM TOUTES LOIS CONFONDUES
C     NPLOI : NOMBRE DE PARAMETRE DE LA LOI PAR DDL
C  OUT
C     RAIDE  : RAIDEUR AU COMPORTEMENT
C     PARAM  : PARAMETRES DE LA LOI
C
C***************** DECLARATION DES VARIABLES LOCALES *******************
C
      INTEGER II,JJ
C
C************ FIN DES DECLARATIONS DES VARIABLES LOCALES ***************

      IF  ( (NOMTE .EQ. 'MECA_DIS_TR_N').OR.
     &      (NOMTE .EQ. 'MECA_DIS_TR_L') ) THEN
         DO 101 II = 0 , 5
            DO 102 JJ = 1 , NPLOI
               IF ( ICODRE(NPLOI*II+JJ).EQ.0 ) THEN
                  PARAM(II+1,JJ) = VALRE(NPLOI*II+JJ)
               ENDIF
102         CONTINUE
101      CONTINUE
         RAIDE(1)= KLV(1)
         RAIDE(2)= KLV(3)
         RAIDE(3)= KLV(6)
         RAIDE(4)= KLV(10)
         RAIDE(5)= KLV(15)
         RAIDE(6)= KLV(21)
      ENDIF
      IF  ( (NOMTE .EQ. 'MECA_DIS_T_N').OR.
     &      (NOMTE .EQ. 'MECA_DIS_T_L') ) THEN
         DO 105 II= 0 , 2
            DO 106 JJ = 1 , NPLOI
               IF ( ICODRE(NPLOI*II+JJ).EQ.0 ) THEN
                  PARAM(II+1,JJ) = VALRE(NPLOI*II+JJ)
               ENDIF
106         CONTINUE
105      CONTINUE
         DO 107 II= 3 , 5
            DO 108 JJ = 1 , NPLOI
               IF ( ICODRE(NPLOI*II+JJ).EQ.0 ) THEN
                  CALL U2MESK('F','DISCRETS_1',1,NOMTE)
               ENDIF
108         CONTINUE
107      CONTINUE
         RAIDE(1)= KLV(1)
         RAIDE(2)= KLV(3)
         RAIDE(3)= KLV(6)
      ENDIF
      IF  ( (NOMTE .EQ. 'MECA_2D_DIS_TR_N').OR.
     &      (NOMTE .EQ. 'MECA_2D_DIS_TR_L') ) THEN
         DO 110 II= 0 , 1
            DO 111 JJ = 1 , NPLOI
               IF ( ICODRE(NPLOI*II+JJ).EQ.0 ) THEN
                  PARAM(II+1,JJ) = VALRE(NPLOI*II+JJ)
               ENDIF
111         CONTINUE
110      CONTINUE
         II= 5
            DO 112 JJ = 1 , NPLOI
               IF ( ICODRE(NPLOI*II+JJ).EQ.0 ) THEN
                  PARAM(3,JJ) = VALRE(NPLOI*II+JJ)
               ENDIF
112         CONTINUE
         DO 113 II= 2 , 4
            DO 114 JJ = 1 , NPLOI
               IF ( ICODRE(NPLOI*II+JJ).EQ.0 ) THEN
                  CALL U2MESK('F','DISCRETS_2',1,NOMTE)
               ENDIF
114         CONTINUE
113      CONTINUE
         RAIDE(1)= KLV(1)
         RAIDE(2)= KLV(3)
         RAIDE(3)= KLV(6)
      ENDIF
      IF  ( (NOMTE .EQ. 'MECA_2D_DIS_T_N').OR.
     &      (NOMTE .EQ. 'MECA_2D_DIS_T_L') ) THEN
         DO 115 II= 0 , 1
            DO 116 JJ = 1 , NPLOI
               IF ( ICODRE(NPLOI*II+JJ).EQ.0 ) THEN
                  PARAM(II+1,JJ) = VALRE(NPLOI*II+JJ)
               ENDIF
116         CONTINUE
115      CONTINUE
         DO 117 II= 2 , 5
            DO 118 JJ = 1 , NPLOI
               IF ( ICODRE(NPLOI*II+JJ).EQ.0 ) THEN
                  CALL U2MESK('F','DISCRETS_3',1,NOMTE)
               ENDIF
118         CONTINUE
117      CONTINUE
         RAIDE(1)= KLV(1)
         RAIDE(2)= KLV(3)
      ENDIF

      END
