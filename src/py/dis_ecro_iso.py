#@ AJOUT dis_ecro_iso Comportement  DATE 07/12/2010   AUTEUR HUMBERT J.HUMBERT 
# -*- coding: iso-8859-1 -*-
#            CONFIGURATION MANAGEMENT OF EDF VERSION
# ======================================================================
# COPYRIGHT (C) 1991 - 2008  EDF R&D                  WWW.CODE-ASTER.ORG
# THIS PROGRAM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY  
# IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY  
# THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR     
# (AT YOUR OPTION) ANY LATER VERSION.                                                  
#                                                                       
# THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT   
# WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF            
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. SEE THE GNU      
# GENERAL PUBLIC LICENSE FOR MORE DETAILS.                              
#                                                                       
# YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE     
# ALONG WITH THIS PROGRAM; IF NOT, WRITE TO EDF R&D CODE_ASTER,         
#    1 AVENUE DU GENERAL DE GAULLE, 92141 CLAMART CEDEX, FRANCE.        
# ======================================================================
# RESPONSABLE VOLDOIRE F.VOLDOIRE

from cata_comportement import LoiComportement

loi = LoiComportement(
   nom            = 'DIS_ECRO_ISO',
   doc            = """Relation de comportement à écrouissage isotrope des elements discrets""",
   num_lc         = 9999,
   nb_vari        = 26,
   nom_vari       = (
		    'UPL_DX','UPL_DY','UPL_DZ','UPL_DRX','UPL_DRY','UPL_DRZ',
		    'VUPL_DX','VUPL_DY','VUPL_DZ','VUPL_DRX','VUPL_DRY','VUPL_DRZ',
		    'P','VP',
		    'X_DX','X_DY','X_DZ','X_DRX','X_DRY','X_DRZ',
		    'VX_DX','VX_DY','VX_DZ','VX_DRX','VX_DRY','VX_DRZ',
		    ),
   modelisation   = ('DIS_T','DIS_TR','2D_DIS_T','2D_DIS_TR'),
   deformation    = ('PETIT','PETIT_REAC', 'GROT_GDEP'),
   nom_varc       = None,
   algo_inte      = ('ANALYTIQUE'),
   type_matr_tang = None,
   proprietes     = None,
)

