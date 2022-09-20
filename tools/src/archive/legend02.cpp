/* ------------------------------------------------------------------------
  
  Programme      :  legend02
  Fichier        :  legend02.c
  
  Auteur         :  Florent Lyard
  Modif(ie par   :  Loren Carrere (passage en C) MARS 2000
 
  Description    :  Gestion des legendes 
  Liste des fonctions : 
 
     
     
-------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "legend.def"

#include "legend.h"

/*--------------------------------------------------------------------*/
/* Fonctions appelees depuis Fortran */
/*--------------------------------------------------------------------*/
#ifdef _add_

extern void lgd_nearestpoint_ ();
#pragma weak lgd_nearestpoint_ = lgd_nearestpoint
extern void lgd_findpoint_ ();
#pragma weak lgd_findpoint_ = lgd_findpoint
extern void lgd_createlegend_ ();
#pragma weak lgd_createlegend_ = lgd_createlegend
extern void lgd_editlegend_ ();
#pragma weak lgd_editlegend_ = lgd_editlegend
extern void lgd_deletelegend_ ();
#pragma weak lgd_deletelegend_ = lgd_deletelegend
extern void lgd_move_drag_ ();
#pragma weak lgd_move_drag_ = lgd_move_drag
extern void lgd_move_end_ ();
#pragma weak lgd_move_end_ = lgd_move_end
extern void lgd_move_init_ ();
#pragma weak lgd_move_init_ = lgd_move_init
extern void lgd_move_cancel_ ();
#pragma weak lgd_move_cancel_ = lgd_move_cancel

#endif
/*--------------------------------------------------------------------*/


