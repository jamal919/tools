
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief Decoding routines and NetCDF parsing definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "poc-netcdf.def"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "fe.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static int sfind(cdfvar_t variable, cdfgbl_t global, const char *c)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos,d,dim,n,verbose=0,status;
  for(pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,c)==0) {
      return(pos);
      break;
      }
    }
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_UGdecode_discretisation(cdfvar_t variable, cdfgbl_t global, UGdecoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int xdim,ydim,zdim,tdim,mdim,ndim,chk;
  int pos,d,dim,n,verbose=0,status;
  int length,shift=0;
  int axis_id;

  decoded->axis=NULL;

  decoded->hdiscretisation=-1;
  decoded->vdiscretisation=-1;

  chk=0;

/*------------------------------------------------------------------------
 */

/* *------------------------------------------------------------------------
  T-UGOm standard*/
  for(pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"T")==0) {
      shift++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"N")==0) {
      decoded->hdiscretisation=LGP1;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"E")==0) {
      decoded->hdiscretisation=NCP1;
      continue;
      }

    if (strcmp(variable.dim[pos].name,LGP2_standards[NNODES_DIMNAME])==0) {
      decoded->hdiscretisation=LGP2;
      continue;
      }
      
//     if (strcmp(variable.dim[pos].name,"CGP2")==0) {
//       decoded->hdiscretisation=LGP2;
//       continue;
//       }

    if (strcmp(variable.dim[pos].name,LGP0_standards[NNODES_DIMNAME])==0) {
      decoded->hdiscretisation=LGP0;
      continue;
      }
    if (strcmp(variable.dim[pos].name,DGP1_standards[NNODES_DIMNAME])==0) {
      decoded->hdiscretisation=D_LGP1;
      continue;
      }
    if (strcmp(variable.dim[pos].name,DNP1_standards[NNODES_DIMNAME])==0) {
      decoded->hdiscretisation=D_NCP1;
      continue;
      }
      
    if (strcmp(variable.dim[pos].name,"M")==0) {
      decoded->hdiscretisation=LGP0;
      if(sfind(variable,global,"P")!=-1) {
        decoded->hdiscretisation=D_LGP1;
        shift++;
        }
      continue;
      }
    if (strcmp(variable.dim[pos].name,"NQUADRANGLES")==0) {
      decoded->hdiscretisation=LGP0;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"L")==0) {
      decoded->vdiscretisation=LGP1;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"K")==0) {
      decoded->vdiscretisation=LGP0;
      continue;
      }
    }

  chk=((decoded->hdiscretisation!=-1)+(decoded->vdiscretisation!=-1));

  if (chk==variable.ndim-shift) goto finished;

/* *------------------------------------------------------------------------
  FVCOM standard*/
  for(pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0) {
      shift=1;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"node")==0) {
      decoded->hdiscretisation=LGP1;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"nele")==0) {
      decoded->hdiscretisation=LGP0;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"siglev")==0) {
      decoded->vdiscretisation=LGP1;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"siglay")==0) {
      decoded->vdiscretisation=LGP0;
      continue;
      }
    }

  chk=((decoded->hdiscretisation!=-1)+(decoded->vdiscretisation!=-1));

  if (chk==variable.ndim-shift) goto finished;


/*------------------------------------------------------------------------
  to be completed*/
  goto error;

  finished:

  status=0;
  return(status);

  error:
  status=-1;
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_UGdecode_axis(cdfvar_t variable, cdfgbl_t global, UGdecoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  {
  int xdim,ydim,zdim,tdim,mdim,ndim,chk;
  int pos,d,dim,n,verbose=0,status;
  int length,shift=0;
  int axis_id;

  decoded->axis=NULL;

  decoded->hdim=-1;
  decoded->vdim=-1;
  decoded->tdim=-1;
  decoded->fdim=-1;

  chk=0;

/*------------------------------------------------------------------------
  Find dimensions order: use "AXIS" attribute if available */

/* *------------------------------------------------------------------------
  T-UGOm standard*/
  axis_id=-1;
  for(n=0;n<variable.natt;n++) {
/* *------------------------------------------------------------------------
    attribute name changed to comply with CF standard*/
    if(strcmp(variable.att[n].name,"content")==0) {
      axis_id=n;
      break;
      }
    }

//   if(axis_id==-1) {
// /**------------------------------------------------------------------------
//     old stuff, quite buggy*/
//     for(n=0;n<variable.natt;n++) {
//       if(strcmp(variable.att[n].name,"axis")==0) {
//         axis_id=n;
//         break;
//         }
//       }
//     }

  if(axis_id!=-1) {
/* *------------------------------------------------------------------------
    T-UGOm standard*/
    decoded->axis=strdup(variable.att[axis_id].data);
    }
  else {
/* *------------------------------------------------------------------------
    FVCOM standard*/
    decoded->axis=(char *) malloc(variable.ndim+1);
    memset(decoded->axis,0,variable.ndim+1);
    pos=0;
    if (strcmp(variable.dim[0].name,"time")==0) {
      decoded->axis[pos]='T';
      pos++;
      }
    if (strcmp(variable.dim[0].name,"f")==0) {
      decoded->axis[pos]='F';
      pos++;
      }
/* *------------------------------------------------------------------------
    pour NoNo*/
    if (strcmp(variable.dim[0].name,"N")==0) {
      decoded->axis[pos]='N';
      pos++;
      }
    if (strcmp(variable.dim[0].name,"node")==0) {
      decoded->axis[pos]='N';
      pos++;
      }
    if (strcmp(variable.dim[0].name,"nele")==0) {
      decoded->axis[pos]='N';
      pos++;
      }
    if (strcmp(variable.dim[0].name,"siglev")==0) {
      decoded->axis[pos]='Z';
      pos++;
      }
    if (strcmp(variable.dim[0].name,"siglay")==0) {
      decoded->axis[pos]='Z';
      pos++;
      }
    if (strcmp(variable.dim[0].name,"three")==0) {
      decoded->axis[pos]='N';
      pos++;
      }
    if (strcmp(variable.dim[0].name,"four")==0) {
      decoded->axis[pos]='N';
      pos++;
      }
    if(pos==variable.ndim) goto parse;
    if(variable.ndim<2) goto parse;
    if (strcmp(variable.dim[1].name,"node")==0) {
      decoded->axis[pos]='N';
      pos++;
      }
    if (strcmp(variable.dim[1].name,"f")==0) {
      decoded->axis[pos]='F';
      pos++;
      }
    if (strcmp(variable.dim[1].name,"nele")==0) {
      decoded->axis[pos]='N';
      pos++;
      }
    if (strcmp(variable.dim[1].name,"siglev")==0) {
      decoded->axis[pos]='Z';
      pos++;
      }
    if (strcmp(variable.dim[1].name,"siglay")==0) {
      decoded->axis[pos]='Z';
      pos++;
      }
    if(pos==variable.ndim) goto parse;
    if(variable.ndim<3) goto parse;
    if (strcmp(variable.dim[2].name,"node")==0) {
      decoded->axis[pos]='N';
      pos++;
      }
    if (strcmp(variable.dim[2].name,"nele")==0) {
      decoded->axis[pos]='N';
      pos++;
      }
    if (strcmp(variable.dim[2].name,"siglev")==0) {
      decoded->axis[pos]='Z';
      pos++;
      }
    if (strcmp(variable.dim[2].name,"siglay")==0) {
      decoded->axis[pos]='Z';
      pos++;
      }
    decoded->axis[pos]=0;
    }

/* *------------------------------------------------------------------------
  T-UGOm standard*/

parse:

  length=strlen(decoded->axis);
  if(length==variable.ndim) {
    for(dim=0;dim<variable.ndim;dim++) {
      d=variable.dim[dim].id;
      if((decoded->axis[dim]=='M') || (decoded->axis[dim]=='m'))
        decoded->hdim=dim;
      if((decoded->axis[dim]=='N') || (decoded->axis[dim]=='n'))
        decoded->hdim=dim;
      if((decoded->axis[dim]=='E') || (decoded->axis[dim]=='e'))
        decoded->hdim=dim;
      if((decoded->axis[dim]=='Z') || (decoded->axis[dim]=='z'))
        decoded->vdim=dim;
      if((decoded->axis[dim]=='T') || (decoded->axis[dim]=='t'))
        decoded->tdim=dim;
      if((decoded->axis[dim]=='F') || (decoded->axis[dim]=='f'))
        decoded->fdim=dim;
      if((decoded->axis[dim]=='P') || (decoded->axis[dim]=='p'))
        shift=-1;
/* *------------------------------------------------------------------------
      WW3 */
//       if((decoded->axis[dim]=='X') || (decoded->axis[dim]=='x')) {
//         decoded->axis[dim]=='N';
//         decoded->hdim=dim;
//         }
//       if((decoded->axis[dim]=='Y') || (decoded->axis[dim]=='y')) {
//         decoded->axis[dim]=='N';
//         decoded->hdim=dim;
//         }
      }
    }

  chk=((decoded->hdim!=-1)+(decoded->vdim!=-1)+(decoded->tdim!=-1)+(decoded->fdim!=-1));

  if (chk==variable.ndim+shift) goto finished;


/*------------------------------------------------------------------------
  to be completed*/
  goto error;

  finished:

  if(decoded->hdim!=-1) decoded->hlen=variable.dim[decoded->hdim].length;
  if(decoded->vdim!=-1) decoded->vlen=variable.dim[decoded->vdim].length;
  if(decoded->tdim!=-1) decoded->tlen=variable.dim[decoded->tdim].length;
  if(decoded->fdim!=-1) decoded->flen=variable.dim[decoded->fdim].length;

/*------------------------------------------------------------------------
  to be completed*/

  status=0;
  return(status);

  error:
  status=-1;
  return(status);
  }


/*------------------------------------------------------------------------

  Nom         :  poc_UGdecodeinfo

--------------------------------------------------------------------------*/
extern "C" int poc_UGdecodeinfo(cdfvar_t variable, cdfgbl_t global, int* nx,int* ny,int* nz,int* nd,int* nt, int *nm, int *nn,char *name, variable_t *details)
  {
  int xdim,ydim,zdim,tdim,mdim,ndim,chk;
  int d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char *axis;
  UGdecoded_t decoded;

  strcpy(name, variable.name);

  *nx=1;
  *ny=1;
  *nz=1;
  *nt=1;
  *nm=1;
  *nn=1;

  xdim=-1;
  ydim=-1;
  zdim=-1;
  tdim=-1;
  mdim=-1;
  ndim=-1;

  chk=0;

  status= poc_UGdecode_axis(variable, global, &decoded);
  if (status!=0) goto error;

  ndim=decoded.hdim;
  zdim=decoded.vdim;
  tdim=decoded.tdim;

  status= poc_UGdecode_discretisation(variable, global, &decoded);
  if (status!=0) goto error;

  finished:

  if(xdim!=-1) *nx=variable.dim[xdim].length;
  if(ydim!=-1) *ny=variable.dim[ydim].length;
  if(zdim!=-1) *nz=variable.dim[zdim].length;
  if(tdim!=-1) *nt=variable.dim[tdim].length;
  if(mdim!=-1) *nm=variable.dim[mdim].length;
  if(ndim!=-1) *nn=variable.dim[ndim].length;

  details->time=0;
  details->name       =strdup(variable.name);
  details->origin     =strdup("none");
  details->units      =strdup("none");
  details->long_name  =strdup(variable.name);
  details->ndim       =variable.ndim;
  details->id         =variable.id;

  status=0;
  return(status);

  error:
  status=-1;
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool isT(const char *name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Checks if the name is for time (i.e. whether the name is one of timeNames)
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,ok=false;
  
  for(i=0; !ok and i<sizeof(timeNames)/sizeof(char *); i++){
    /** case insensitive */
    ok=(strcasecmp(name,timeNames[i])==0);
    }
  return ok;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool isT(const cdfdim_t &dim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Checks if a dimension is a time dimension

  dim : reference to the dimension

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  return isT(dim.name);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_axis_MarsPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_axis_MarsPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"latitude")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"longitude")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
        }
    if (strcmp(variable.dim[pos].name,"level")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis_RomsPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  const char *names[]={//dimnames->vnames
    "time","T",
    "eta_","Y",
    "xi_","X",
    "level","Z",
    NULL};
  int i,ok;//name index, boolean

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    for(i=0;names[i]!=NULL;i+=2){
      ok=(strncmp(names[i],variable.dim[pos].name)==0);
      if(ok) {
        decoded->axis[pos]=names[i+1][0];
        chk++;
        break;
        }
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis_WRFPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
//         float U(Time, bottom_top, south_north, west_east_stag) ;
//                 U:FieldType = 104 ;
//                 U:MemoryOrder = "XYZ" ;
//                 U:description = "x-wind component" ;
//                 U:units = "m s-1" ;
//                 U:stagger = "X" ;
//                 U:coordinates = "XLONG_U XLAT_U" ;
  axis_id=-1;
  for(n=0;n<variable.natt;n++) {
    if(strcmp(variable.att[n].name,"MemoryOrder")==0) {
      axis_id=n;
      break;
      }
    }

  if(axis_id!=-1) {
//    decoded->axis=strdup(variable.att[axis_id].data);
    status=0;
    }
  else{
    status=-1;
    return(status);
    }

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  length=strlen(variable.att[axis_id].data);
  while(variable.att[axis_id].data[length-1]==' ') length--;

  chk=0;
  if(length==variable.ndim-1) {
    decoded->axis[chk]='T';
    chk++;
    }

  for (pos=0;pos<length;pos++) {
    decoded->axis[chk]=variable.att[axis_id].data[length-pos-1];
    chk++;
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_axis_WW3Patch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_axis_WW3Patch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"f")==0) {
      decoded->axis[pos]='W';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"latitude")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"longitude")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_axis_XCONVPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_axis_XCONVPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"t")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"latitude")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"longitude")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
        }
    if (strcmp(variable.dim[pos].name,"entire_atomosphere")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_axis_ECMWFPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_axis_ECMWFPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"lat")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"lon")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_axis_ECMWFPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_axis_WOAPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"lat")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"lon")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"depth")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis_ECCOPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"TIME")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"LATITUDE_T")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"LONGITUDE_T")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"DEPTH_T")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_axis_LevitusPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_axis_LevitusPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"T")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"Z")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"Y")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"X")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_axis_LegosPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_axis_LegosPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"w")==0 or strcmp(variable.dim[pos].name,"W")==0) {
      decoded->axis[pos]='W';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"time")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"z")==0 or strcmp(variable.dim[pos].name,"Z")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"y")==0 or strcmp(variable.dim[pos].name,"Y")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"x")==0 or strcmp(variable.dim[pos].name,"X")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// int poc_decode_axis_MercatorPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   {
//   int chk;
//   int pos,d,dim,n,verbose=0,status;
//   int length;
//   int axis_id;
//
//   decoded->axis=new char[variable.ndim+1];
//   memset(decoded->axis,0,variable.ndim+1);
//
//   chk=0;
//   for (pos=0;pos<variable.ndim;pos++) {
//     if (strcmp(variable.dim[pos].name,"w")==0) {
//       decoded->axis[pos]='W';
//       chk++;
//       continue;
//       }
//     if (strcmp(variable.dim[pos].name,"z")==0) {
//       decoded->axis[pos]='Z';
//       chk++;
//       continue;
//       }
//     if (strcmp(variable.dim[pos].name,"y")==0) {
//       decoded->axis[pos]='Y';
//       chk++;
//       continue;
//       }
//     if (strcmp(variable.dim[pos].name,"x")==0) {
//       decoded->axis[pos]='X';
//       chk++;
//       continue;
//       }
//     }
//
//   if(chk==variable.ndim) return(0);
//   else {
//     delete[] decoded->axis;
//     decoded->axis=0;
//     return(-1);
//     }
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis_POMPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
   {
   int chk;
   int pos,d,dim,n,verbose=0,status;
   int length;
   int axis_id;

   decoded->axis=(char *) malloc(variable.ndim+1);
   for(pos=0;pos<variable.ndim+1;pos++) decoded->axis[pos]=0;

   chk=0;
   for (pos=0;pos<variable.ndim;pos++) {
     if (strcmp(variable.dim[pos].name,"time")==0) {
       decoded->axis[pos]='T';
       chk++;
       continue;
       }
     if (strcmp(variable.dim[pos].name,"zpos")==0) {
       decoded->axis[pos]='Z';
       chk++;
       continue;
       }
     if (strcmp(variable.dim[pos].name,"ypos")==0) {
       decoded->axis[pos]='Y';
       chk++;
       continue;
       }
     if (strcmp(variable.dim[pos].name,"xpos")==0) {
       decoded->axis[pos]='X';
       chk++;
       continue;
       }
     }

   if(chk==variable.ndim) return(0);
   else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis_EMODNETPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
   {
   int chk;
   int pos,d,dim,n,verbose=0,status;
   int length;
   int axis_id;

   decoded->axis=(char *) malloc(variable.ndim+1);
   for(pos=0;pos<variable.ndim+1;pos++) decoded->axis[pos]=0;

   chk=0;
   for (pos=0;pos<variable.ndim;pos++) {
     if (strcmp(variable.dim[pos].name,"position_lat")==0) {
       decoded->axis[pos]='Y';
       chk++;
       continue;
       }
     if (strcmp(variable.dim[pos].name,"position_long")==0) {
       decoded->axis[pos]='X';
       chk++;
       continue;
       }
     }

   if(chk==variable.ndim) return(0);
   else return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis_SymphoniePatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
   {
   int chk;
   int pos,d,dim,n,verbose=0,status;
   int length;
   int axis_id;

   decoded->axis=(char *) malloc(variable.ndim+1);
   for(pos=0;pos<variable.ndim+1;pos++) decoded->axis[pos]=0;

   chk=0;
   for (pos=0;pos<variable.ndim;pos++) {
     if (strcmp(variable.dim[pos].name,"ni_t")==0) {
       decoded->axis[pos]='X';
       chk++;
       continue;
       }
     if (strcmp(variable.dim[pos].name,"nj_t")==0) {
       decoded->axis[pos]='Y';
       chk++;
       continue;
       }
     }

   if(chk==variable.ndim) return(0);
   else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_axis_GenericPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_axis_GenericPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;

  decoded->axis=(char *) malloc(variable.ndim+1);
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"nT")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"nZ")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"nY")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"nX")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"y1")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"x1")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  
  for (pos=0;pos<variable.natt;pos++) {
    if (strcmp(variable.att[pos].name,"coordinates")==0) {
      char **coordinates;
      n=get_token(variable.att[pos].data,&coordinates);
      chk=0;
      for(d=0;d<n;d++){
        if(strncmp("lon",coordinates[d])==0) {
          decoded->axis[d]='X';
          chk++;
          continue;
          }
        if(strncmp("lat",coordinates[d])==0) {
          decoded->axis[d]='Y';
          chk++;
          continue;
          }
        }
      break;
      }
    }
  
  if(chk==variable.ndim) return(0);
  
  return NC_ENOTATT;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis_MercatorPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"w")==0) {
      decoded->axis[pos]='W';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"time_counter")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"t")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"z")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"deptht")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"y")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"x")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis_ESRIPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"LINES")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"COLUMNS")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis_HYCOMPatch2(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=new char[variable.ndim+1];
  memset(decoded->axis,0,variable.ndim+1);

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"MT")==0) {
      decoded->axis[pos]='T';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"X")==0) {
      decoded->axis[pos]='X';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"Y")==0) {
      decoded->axis[pos]='Y';
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"Z")==0) {
      decoded->axis[pos]='Z';
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) return(0);
  else {
    delete[] decoded->axis;
    decoded->axis=0;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_axis(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int xdim,ydim,zdim,tdim,mdim,ndim,chk;
  int dim,n,verbose=0,status;
  int length;
  int axis_id;

  decoded->axis=NULL;

  decoded->xdim=-1;
  decoded->ydim=-1;
  decoded->zdim=-1;
  decoded->tdim=-1;

  chk=0;

/*------------------------------------------------------------------------
  Find dimensions order: use "AXIS" attribute if available */
  axis_id=-1;
  ///\todo use variable.getattr("content") instead
  for(n=0;n<variable.natt;n++) {
/* *------------------------------------------------------------------------
    attribute name changed to comply with CF standard*/
    if(strcmp(variable.att[n].name, "content") == 0) {
      axis_id = n;
      break;
      }
    }
    
// /*------------------------------------------------------------------------
//   latest comodo standard */
//   if(axis_id==-1) {
//     for(dim=0;dim<variable.ndim;dim++) {
//       for(n=0;n<global.ndimsp;n++){
//         if(strcmp(global.dimension[n].name,variable.dim[dim].name)==0)
//           break;
//         }
//       cdfdim_t *dimvar=&global.dimension[n];
//       for(n=0;n<dimvar->natts;n++) {
//         if(strcmp(dimvar->att[n].name,"axis")==0) {
//           axis_id=n;
//           break;
//           }
//         }
//       }
//     }

  bool allow_obsolete=false;
  
  if(allow_obsolete) {
    if(axis_id==-1) {
    for(n=0;n<variable.natt;n++) {
/* *------------------------------------------------------------------------
      attribute name changed to comply with CF standard*/
    ///\todo get a better approach of old standard
      if(strcmp(variable.att[n].name,"axis")==0) {
        axis_id=n;
        break;
        }
      }
    }
  }
  
  
  if(axis_id!=-1) {
    decoded->axis=poc_strdup(variable.att[axis_id].data);
    status=0;
    }
  else{
    status=NC_ENOTATT; /* Attribute not found */
    }

  if(status!=0) {
    status=poc_decode_axis_WRFPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_RomsPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_MarsPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_WW3Patch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_ECMWFPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_LevitusPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_XCONVPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_LegosPatch( variable,  global, decoded);
    }

   if(status!=0) {
    status=poc_decode_axis_POMPatch( variable,  global, decoded);
    }

   if(status!=0) {
    status=poc_decode_axis_EMODNETPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_GenericPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_WOAPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_MercatorPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_ESRIPatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_SymphoniePatch( variable,  global, decoded);
    }

  if(status!=0) {
    status=poc_decode_axis_ECCOPatch( variable,  global, decoded);
    }
    
  if(status!=0) {
    status=poc_decode_axis_HYCOMPatch2( variable,  global, decoded);
    }
    
//   if(axis_id!=-1) {
//     decoded->axis=strdup(variable.att[axis_id].data);
  if(status==0) {
    length=strlen(decoded->axis);
    if(length==variable.ndim) {
      for(dim=0;dim<variable.ndim;dim++) {
        if(toupper(decoded->axis[dim])=='X')
          decoded->xdim=dim;
        if(toupper(decoded->axis[dim])=='Y')
          decoded->ydim=dim;
        if(toupper(decoded->axis[dim])=='Z')
          decoded->zdim=dim;
        if(toupper(decoded->axis[dim])=='T')
          decoded->tdim=dim;
/* *------------------------------------------------------------------------
        patch due to SYMPHONIE's NETCDF bug*/
        if(toupper(decoded->axis[dim])=='F')
          decoded->tdim=dim;
/* *------------------------------------------------------------------------
        patch due to old LEGOS NETCDF*/
        if(toupper(decoded->axis[dim])=='W')
          decoded->tdim=dim;
        }
      }
    }

  chk=((decoded->xdim!=-1)+(decoded->ydim!=-1)+(decoded->zdim!=-1)+(decoded->tdim!=-1));

  if (chk==variable.ndim) goto finished;

/*------------------------------------------------------------------------
  to be completed*/
  goto error;

  finished:

  if(decoded->xdim!=-1) decoded->xlen=variable.dim[decoded->xdim].length;
  if(decoded->ydim!=-1) decoded->ylen=variable.dim[decoded->ydim].length;
  if(decoded->zdim!=-1) decoded->zlen=variable.dim[decoded->zdim].length;
  if(decoded->tdim!=-1) decoded->tlen=variable.dim[decoded->tdim].length;
  if(decoded->fdim!=-1) decoded->tlen=variable.dim[decoded->fdim].length;

/*------------------------------------------------------------------------
  to be completed*/

  status=0;
  return(status);

error:
  status=NC_ENOTATT; /* Attribute not found */
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_units(cdfvar_t variable, decoded_t *decoded, double *factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int n,status;
  char *text=0;
  for(n=0;n<variable.natt;n++) {
    if(strcmp(variable.att[n].name,"units")==0){
      decoded->units=strdup(variable.att[n].data);
      text=strdup(variable.att[n].data);
      if(strcmp(text,"radian_east") ==0) *factor=180./M_PI;;
      if(strcmp(text,"degrees_east")==0) *factor=  1.;
      if(strcmp(text,"degrees_west")==0) *factor= -1.;
      if(strcmp(text,"degree_east") ==0) *factor=  1.;
      if(strcmp(text,"degrees") ==0)     *factor=1.;
      if(strcmp(text,"radian_north")==0) *factor=180./M_PI;
      if(strcmp(text,"degree_north")==0) *factor=1.;
      if(strcmp(text,"degree_north")==0) *factor=1.;
      status=0;
      break;
      }
    }
  if(text!=0) free(text);
  return(status);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// int poc_decode_mask(const cdfvar_t variable, decoded_t *decoded)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   {
//   int xdim,ydim,zdim,tdim,mdim,ndim,chk;
//   int d,dim,n,verbose=0,status;
//
//   decoded->scale=1;
//   decoded->offset=0;
//   decoded->spec=1./0.;
//
// /*------------------------------------------------------------------------
//   Find dimensions order: use "AXIS" attribute if available */
//   for(n=0;n<variable.natt;n++) {
//     if(strcmp(variable.att[n].name,"scale_factor")==0)  decoded->scale ==(float) *((float *) variable.att[n].data);
//     if(strcmp(variable.att[n].name,"add_offset")  ==0)  decoded->offset==(float) *((float *) variable.att[n].data);
//     if(strcmp(variable.att[n].name,"missing_value")==0) {
//       decoded->spec=(float) *((float *) variable.att[n].data);
//       }
//     if(strcmp(variable.att[n].name,"_FillValue")==0) {
//       decoded->spec=(float) *((float *) variable.att[n].data);
//       }
//     }
//
// /*------------------------------------------------------------------------
//   to be completed*/
//   status=0;
//   return(status);
//
//   error:
//   status=-1;
//   return(status);
//   }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int decoded_nvalues(decoded_t decoded, int *ndims)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Gives number of values and optionally number of dimensions
/**
\date 2011-11-21 Damien ALLAIN : creation

\param decoded
\param *ndims Pointer to number of dimensions. If NULL, without effect. Default: NULL.
\returns the number of values
*/
/*----------------------------------------------------------------------------*/
{
  int nvalues=1,ndims_=0;
  if(decoded.xlen) {
    nvalues*=decoded.xlen;
    ndims_++;
    }
  if(decoded.ylen) {
    nvalues*=decoded.ylen;
    ndims_++;
    }
  if(decoded.zlen) {
    nvalues*=decoded.zlen;
    ndims_++;
    }
  if(ndims_==0)
    nvalues=0;
  if(ndims!=NULL)
    *ndims=ndims_;
  return nvalues;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_decode_mask_template(const cdfvar_t &variable, T *scale, T *offset, T *spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/** \brief gets scale, offset and mask value

\param variable
\param *scale pointer to scale value. Defaults to 1 if not found in the attributes. Without effect if NULL
\param *offset pointer to mask value. Defaults to 0 if not found in the attributes. Without effect if NULL
\param *spec pointer to mask value. Defaults to 0 if not found in the attributes. Without effect if NULL
\returns NC_ENOTATT if one of the requested attributes is not found (which is not a big deal) or NC_NOERR otherwise.
*/
{
  int n,status=NC_NOERR;
  T scale_,offset_,spec_,valid_min,valid_max;
  int has_scale,has_offset,has_spec,has_min,has_max;
  
  has_scale=has_offset=has_spec=has_min=has_max=0;
  
  scale_ =1.;
  offset_=0.;
  
  ///\date 2011-09-27 Damien Allain : changed default value from Inf to 0.
  ///\date 2011-12-20 Florent Lyard : changed default value from 0 to 1.e+10.
  
  spec_=1.e+10;

  for(n=0;n<variable.natt;n++) {
    switch (variable.att[n].type) {

      case NC_BYTE:/* signed 1 byte integer */
        if(strcmp(variable.att[n].name,"scale_factor")==0) {
          scale_ = (T) *((signed char *) variable.att[n].data);
          has_scale=1;
          }
        if(strcmp(variable.att[n].name,"add_offset")  ==0) {
          offset_=(T) *((signed char *) variable.att[n].data);
          has_offset=1;
          }
        if(strcmp(variable.att[n].name,"missing_value")==0) {
          spec_=(T) *((signed char *) variable.att[n].data);
          has_spec=1;
          }
        if(strcmp(variable.att[n].name,"_FillValue")==0) {
          spec_=(T) *((signed char *) variable.att[n].data);
          has_spec=1;
          }
        break;

      case NC_SHORT:
        if(strcmp(variable.att[n].name,"scale_factor")==0) {
          scale_ = (T) *((short *) variable.att[n].data);
          has_scale=1;
          }
        if(strcmp(variable.att[n].name,"add_offset")  ==0) {
          offset_=(T) *((short *) variable.att[n].data);
          has_offset=1;
          }
        if(strcmp(variable.att[n].name,"missing_value")==0) {
          spec_=(T) *((short *) variable.att[n].data);
          has_spec=1;
          }
        if(strcmp(variable.att[n].name,"_FillValue")==0) {
          spec_=(T) *((short *) variable.att[n].data);
          has_spec=1;
          }
        break;

      case NC_INT:
        if(strcmp(variable.att[n].name,"scale_factor")==0) {
          scale_ = (T) *((int *) variable.att[n].data);
          has_scale=1;
          }
        if(strcmp(variable.att[n].name,"add_offset")  ==0) {
          offset_=(T) *((int *) variable.att[n].data);
          has_offset=1;
          }
        if(strcmp(variable.att[n].name,"missing_value")==0) {
          spec_=(T) *((int *) variable.att[n].data);
          has_spec=1;
          }
        if(strcmp(variable.att[n].name,"_FillValue")==0) {
          spec_=(T) *((int *) variable.att[n].data);
          has_spec=1;
          }
        break;

      case NC_FLOAT:
        if(strcmp(variable.att[n].name,"scale_factor")==0) {
          scale_ = (T) *((float *) variable.att[n].data);
          has_scale=1;
          }
        if(strcmp(variable.att[n].name,"add_offset")  ==0) {
          offset_=(T) *((float *) variable.att[n].data);
          has_offset=1;
          }
        if(strcmp(variable.att[n].name,"missing_value")==0) {
          spec_=(T) *((float *) variable.att[n].data);
          has_spec=1;
          }
        if(strcmp(variable.att[n].name,"mask")==0) {
//          mask_=strdup(variable.att[n].data);
          }
        if(strcmp(variable.att[n].name,"_FillValue")==0) {
          spec_=(T) *((float *) variable.att[n].data);
          has_spec=1;
          }
        if(strcmp(variable.att[n].name,"valid_min")==0) {
          valid_min=(T) *((float *) variable.att[n].data);
          has_min=1;
          }
        if(strcmp(variable.att[n].name,"valid_max")==0) {
          valid_max=(T) *((float *) variable.att[n].data);
          has_max=1;
          }
        break;

      case NC_DOUBLE:
        if(strcmp(variable.att[n].name,"scale_factor")==0) {
          scale_ = (T) *((double *) variable.att[n].data);
          has_scale=1;
          }
        if(strcmp(variable.att[n].name,"add_offset")  ==0) {
          offset_=(T) *((double *) variable.att[n].data);
          has_offset=1;
          }
        if(strcmp(variable.att[n].name,"missing_value")==0) {
          spec_=(T) *((double *) variable.att[n].data);
          has_spec=1;
          }
        if(strcmp(variable.att[n].name,"mask")==0) {
//          mask_=strdup(variable.att[n].data);
          }
        if(strcmp(variable.att[n].name,"_FillValue")==0) {
          spec_=(T) *((double *) variable.att[n].data);
          has_spec=1;
          }
        break;

      }
    }

  if(scale!=NULL){
    if(!has_scale) status=NC_ENOTATT;/* Attribute not found */
    *scale=scale_;
    }
  if(offset!=NULL){
    if(!has_offset) status=NC_ENOTATT;/* Attribute not found */
    *offset=offset_;
    }
    
  if(!has_spec) {
    if(has_min && has_max) {
      spec_=valid_max-valid_min+1;
      has_spec=1;
      }
    }
  if(spec!=NULL){
    if(!has_spec) status=NC_ENOTATT;
    else status=0;
    *spec=spec_;
    }

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_mask(const cdfvar_t &variable, char *scale, char *offset, char *spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(variable,scale,offset,spec);
  
  return result;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_mask(const cdfvar_t &variable, signed char *scale, signed char *offset, signed char *spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(variable,scale,offset,spec);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_mask(const cdfvar_t &variable, short int *scale, short int *offset, short int *spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(variable,scale,offset,spec);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_mask(const cdfvar_t &variable, int *scale, int *offset, int *spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(variable,scale,offset,spec);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_mask(const cdfvar_t &variable, float *scale, float *offset, float *spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(variable,scale,offset,spec);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_mask(const cdfvar_t &variable, double *scale, double *offset, double *spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(variable,scale,offset,spec);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_decode_mask(const cdfvar_t &variable, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///wrapper for poc_decode_mask()
/*----------------------------------------------------------------------------*/
{
  return poc_decode_mask(variable,&decoded->scale,&decoded->offset,&decoded->spec);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_names(cdfvar_t variable, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int xdim,ydim,zdim,tdim,mdim,ndim,chk;
  int d,dim,n,verbose=0,status;

  decoded->long_name    =0;
  decoded->standard_name=0;
  decoded->short_name   =0;
  decoded->units        =0;

/*------------------------------------------------------------------------
  Find dimensions order: use "AXIS" attribute if available */
  for(n=0;n<variable.natt;n++) {
    if(strcmp(variable.att[n].name,"long_name")==0)      decoded->long_name     = strdup(variable.att[n].data);
    if(strcmp(variable.att[n].name,"standard_name")==0)  decoded->standard_name = strdup(variable.att[n].data);
    if(strcmp(variable.att[n].name,"short_name")==0)     decoded->short_name    = strdup(variable.att[n].data);
    if(strcmp(variable.att[n].name,"units")==0)          decoded->units         = strdup(variable.att[n].data);
    }

/*------------------------------------------------------------------------
  to be completed*/
  status=0;
  return(status);

  error:
  status=-1;
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_names(cdfgbl_t variable, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int xdim,ydim,zdim,tdim,mdim,ndim,chk;
  int d,dim,n,verbose=0,status;

  decoded->production  = 0;
  
/*------------------------------------------------------------------------
  Find dimensions order: use "AXIS" attribute if available */
  for(n=0;n<variable.ngattsp;n++) {
    if(strcmp(variable.attribute[n].name,"production")==0) decoded->production=strdup(variable.attribute[n].data);
    }

/*------------------------------------------------------------------------
  to be completed*/
  status=0;
  return(status);

  error:
  status=-1;
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_associates_MarsPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_associates_MarsPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0) {
      strcat(vnames,"time ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"latitude")==0) {
      strcat(vnames,"latitude ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"longitude")==0) {
      strcat(vnames,"longitude ");
      chk++;
      continue;
        }
    if (strcmp(variable.dim[pos].name,"level")==0) {
      strcat(vnames,"SIG ");
      chk++;
      continue;
      }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_associates_RomsPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};
  const char *names[]={//dimnames->vnames
    "time","scrum_time",
    "eta_rho","lat_rho",
    "xi_rho","lon_rho",
    "eta_u","lat_u",
    "xi_u","lon_u",
    "eta_v","lat_v",
    "xi_v","lon_v",
    "eta_psi","lat_psi",
    "xi_psi","lon_psi",
    "level","SIG",
    NULL};
  int i,ok;//name index, boolean

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    for(i=0;names[i]!=NULL;i+=2){
      ok=(strcmp(variable.dim[pos].name,names[i])==0);
      if(ok) {
        strcat(vnames,names[i+1]);
        strcat(vnames," ");
        chk++;
        break;
        }
      }
    }

  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_associates_WRFPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int associate_id;
  char vnames[1024]={"\0"},xname[64],yname[64];
//         float U(Time, bottom_top, south_north, west_east_stag) ;
//                 U:FieldType = 104 ;
//                 U:MemoryOrder = "XYZ" ;
//                 U:description = "x-wind component" ;
//                 U:units = "m s-1" ;
//                 U:stagger = "X" ;
//                 U:coordinates = "XLONG_U XLAT_U" ;

// 	float XLONG(Time, south_north, west_east) ;
// 		XLONG:FieldType = 104 ;
// 		XLONG:MemoryOrder = "XY " ;
// 		XLONG:description = "LONGITUDE, WEST IS NEGATIVE" ;
// 		XLONG:units = "degree_east" ;
// 		XLONG:stagger = "" ;

  associate_id=-1;
  for(n=0;n<variable.natt;n++) {
    if(strcmp(variable.att[n].name,"coordinates")==0) {
      associate_id=n;
      break;
      }
    }

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/
  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"Time")==0) {
      strcat(vnames,"Times ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"bottom_top")==0) {
      strcat(vnames,"ZNU ");
      chk++;
      continue;
      }
    }
/* *----------------------------------------------------------------------------
  Patch based on (partial) attribute*/
  if(associate_id!=-1) {
    sscanf(variable.att[associate_id].data,"%s %s",xname,yname);
    strcat(vnames,yname);
    strcat(vnames," ");
    chk++;
    strcat(vnames,xname);
    strcat(vnames," ");
    chk++;
    }
  else{
    for (pos=0;pos<variable.ndim;pos++) {
      if (strcmp(variable.dim[pos].name,"south_north")==0) {
        strcat(vnames,"XLAT ");
        chk++;
        continue;
        }
      if (strcmp(variable.dim[pos].name,"west_east")==0) {
        strcat(vnames,"XLONG ");
        chk++;
        continue;
        }
      }
    }

  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_associates_WW3Patch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_associates_WW3Patch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0) {
      strcat(vnames,"time ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"f")==0) {
      strcat(vnames,"f ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"latitude")==0) {
      strcat(vnames,"latitude ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"longitude")==0) {
      strcat(vnames,"longitude ");
      chk++;
      continue;
        }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_associates_ECMWFPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_associates_ECMWFPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0) {
      strcat(vnames,"time ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"lat")==0) {
      strcat(vnames,"lat ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"lon")==0) {
      strcat(vnames,"lon ");
      chk++;
      continue;
        }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=poc_strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_associates_ECCOPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"TIME")==0) {
      strcat(vnames,"TIME ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"LATITUDE_T")==0) {
      strcat(vnames,"LATITUDE_T ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"LONGITUDE_T")==0) {
      strcat(vnames,"LONGITUDE_T ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"DEPTH_T")==0) {
      strcat(vnames,"DEPTH_T ");
      chk++;
      continue;
      }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=poc_strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_associates_EMODNETPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"position_lat")==0) {
      strcat(vnames,"lat ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"position_long")==0) {
      strcat(vnames,"lon ");
      chk++;
      continue;
        }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=poc_strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_associates_LevitusPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_associates_LevitusPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"T")==0) {
      strcat(vnames,"T ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"Y")==0) {
      strcat(vnames,"Y ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"X")==0) {
      strcat(vnames,"X ");
      chk++;
      continue;
        }
    if (strcmp(variable.dim[pos].name,"Z")==0) {
      strcat(vnames,"Z ");
      chk++;
      continue;
      }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=poc_strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_associates_LegosPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_associates_LegosPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"W")==0) {
      strcat(vnames,"spectrum ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"T")==0) {
      strcat(vnames,"time ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"Y")==0) {
      strcat(vnames,"lat ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"X")==0) {
      strcat(vnames,"lon ");
      chk++;
      continue;
        }
    if (strcmp(variable.dim[pos].name,"Z")==0) {
      strcat(vnames,"Z ");
      chk++;
      continue;
      }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_associates_MercatorPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status,id;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0 or strcmp(variable.dim[pos].name,"t")==0) {
      id=cdf_identify(global,"time");
      if(id!=-1) {
        strcat(vnames,"time ");
        chk++;
        continue;
        }
      id=cdf_identify(global,"time_counter");
      if(id!=-1) {
        strcat(vnames,"time_counter ");
        chk++;
        continue;
        }
      }
    if (strcmp(variable.dim[pos].name,"y")==0) {
      strcat(vnames,"nav_lat ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"x")==0) {
      strcat(vnames,"nav_lon ");
      chk++;
      continue;
        }
    if (strcmp(variable.dim[pos].name,"z")==0) {
      strcat(vnames,"nav_lev ");
      chk++;
      continue;
      }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_associates_ORCAPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id,id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time_counter")==0) {
      strcat(vnames,"time_counter ");
      id=cdf_identify(global,"time_counter");
      if(id==-1) break;
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"y")==0) {
      strcat(vnames,"nav_lat ");
      id=cdf_identify(global,"nav_lat");
      if(id==-1) break;
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"x")==0) {
      strcat(vnames,"nav_lon ");
      id=cdf_identify(global,"nav_lon");
      if(id==-1) break;
      chk++;
      continue;
        }
    if (strcmp(variable.dim[pos].name,"deptht")==0) {
      strcat(vnames,"deptht ");
      chk++;
      continue;
      }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  else return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_associates_ESRIPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id,id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"LINES")==0) {
      strcat(vnames,"LINES ");
      id=cdf_identify(global,"LINES");
      if(id==-1) break;
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"COLUMNS")==0) {
      strcat(vnames,"COLUMNS ");
      id=cdf_identify(global,"COLUMNS");
      if(id==-1) break;
      chk++;
      continue;
        }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_associates_SymphoniePatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id,id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"ni_t")==0) {
      strcat(vnames,"longitude_t ");
      id=cdf_identify(global,"longitude_t");
      if(id==-1) break;
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"nj_t")==0) {
      strcat(vnames,"latitude_t ");
      id=cdf_identify(global,"latitude_t");
      if(id==-1) break;
      chk++;
      continue;
        }
    }


  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//int poc_decode_associates_GenericPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded) __attribute__((warning("This does not use dimname SO WILL NOT WORK IF DATA AND GRID FILES DO NOT HAVE THEIR DIMENSIONS AT THE SAME INDICES !!!")));
int poc_decode_associates_GenericPatch(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int chk;
  int pos,d,dim,n,verbose=0,status;
  int length;
  int axis_id;
  char vnames[1024]={"\0"};

/* *----------------------------------------------------------------------------
  Patch based on dimension name identification with variables*/

  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"time")==0) {
      strcat(vnames,"time ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"lat")==0) {
      strcat(vnames,"lat ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"lon")==0) {
      strcat(vnames,"lon ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"depth")==0) {
      strcat(vnames,"depth ");
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
    
  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"nT")==0) {
      strcat(vnames,"T ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"nY")==0) {
      strcat(vnames,"Lat ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"nX")==0) {
      strcat(vnames,"Lon ");
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
    
  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"y1")==0) {
      strcat(vnames,"y1 ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"x1")==0) {
      strcat(vnames,"x1 ");
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  
  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"y")==0) {
      int id=cdf_identify(global,"y");
      if(id==-1) break;
      strcat(vnames,"y ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"x")==0) {
      int id=cdf_identify(global,"x");
      if(id==-1) break;
      strcat(vnames,"x ");
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  
  chk=0;
  for (pos=0;pos<variable.ndim;pos++) {
    if (strcmp(variable.dim[pos].name,"nj")==0) {
      int id=cdf_identify(global,"latitude");
      if(id==-1) break;
      strcat(vnames,"latitude ");
      chk++;
      continue;
      }
    if (strcmp(variable.dim[pos].name,"ni")==0) {
      int id=cdf_identify(global,"longitude");
      if(id==-1) break;
      strcat(vnames,"longitude ");
      chk++;
      continue;
      }
    }

  if(chk==variable.ndim) {
    vnames[strlen(vnames)-1]=0;
    decoded->associate=strdup(vnames);
    return(0);
    }
  
  return(-1);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_associates(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded, int initialize, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  {
  int xdim,ydim,zdim,tdim,mdim,ndim,chk;
  int d,dim,n,verbose=0,status;
  int length;
  int associate_id,v;
  char *cdum,*text;

  status= poc_decode_axis( variable,  global, decoded);
  if(status!=0) return(status);

  decoded->associate=NULL;

//   decoded->vx=-1;
//   decoded->vy=-1;
//   decoded->vz=-1;
//   decoded->vt=-1;
  if(initialize==1) {
    decoded->vx=-1;
    decoded->vy=-1;
    decoded->vz=-1;
    decoded->vt=-1;
    }


  chk=0;

/*------------------------------------------------------------------------
  Find dimensions order: use "AXIS" attribute if available */
  associate_id=-1;
  for(n=0;n<variable.natt;n++) {
    if(strcmp(variable.att[n].name,"associate")==0) {
      associate_id=n;
      break;
      }
    }
  if(associate_id==-1) {
    for(n=0;n<variable.natt;n++) {
      if(strcmp(variable.att[n].name,"coordinates")==0) {
        associate_id=n;
        break;
        }
      }
    }

//   if(associate_id!=-1) {
//     decoded->associate=strdup(variable.att[associate_id].data);
//     length=strlen(decoded->associate);
//     }
//   else {
//     goto error;
//     }
  if(associate_id!=-1) {
    if(strcmp(variable.att[associate_id].data,"lon lat")==0) {
      printf("untrusted coordinates attribute <lon lat> \n");
      strcpy(variable.att[associate_id].data,"lat lon");
      status=0;
      }
    if(strcmp(variable.att[associate_id].data,"Longitude Latitude Date")==0) {
      printf("untrusted coordinates attribute <Longitude Latitude Date> \n");
      strcpy(variable.att[associate_id].data,"MT Latitude Longitude");
      status=0;
      }
    decoded->associate=poc_strdup(variable.att[associate_id].data);
    length=strlen(decoded->associate);
    }
  else{
    status=-1;
    }

  /// \todo 2011-09-26 Damien Allain : re-write what follows with a do{}while(0) loop

  if(status!=0) {
    status=poc_decode_associates_GenericPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }
    
  if(status!=0) {
    status=poc_decode_associates_RomsPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_MarsPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_WW3Patch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_WRFPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_LevitusPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_ECMWFPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_EMODNETPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_MercatorPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_ORCAPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_ESRIPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_SymphoniePatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

  if(status!=0) {
    status=poc_decode_associates_ECCOPatch( variable,  global, decoded);
    if(status!=0) decoded->associate=0;
    }

/*------------------------------------------------------------------------
  find grid information from "associate" attribute */
  if(debug) 
    STDERR_BASE_LINE_FUNC("*** BECAUSE OF THE CODE BELOW, LON AND LAT WILL BE SWAPPED WITH ONE ANOTHER WHEN THE coordinates ATTRIBUTE IS NOT IN THE SAME ORDER AS THE DIMENSIONS! ***\n",__FUNCTION__);
  if(decoded->associate!=NULL) {
    cdum=strdup(decoded->associate);
    for(dim=0;dim<variable.ndim;dim++) {
      if( dim==0) text=strtok(cdum," ");
      else text=strtok(NULL," ");
      if(text !=NULL) {
        v=cdf_identify(global, text);
        if(v!=-1) {
          if(dim==decoded->xdim) decoded->vx=v;
          if(dim==decoded->ydim) decoded->vy=v;
          if(dim==decoded->zdim) decoded->vz=v;
          if(dim==decoded->tdim) decoded->vt=v;
//          printf("associate %d, name %s %d\n",dim,text,v);
          }
        }
      }
    }

/*------------------------------------------------------------------------
  to be completed*/

  status=0;
  return(status);

  error:
  status=-1;
  return(status);
  }

