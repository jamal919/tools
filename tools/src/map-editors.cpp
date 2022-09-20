
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

/******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief grid edition functions

Apart from operators, that are used e.g. within poc-formula-parse.ypp,
functions defined in this file are declared directly in poc-data-operators.cpp
before poc_formula_implement_func(), the only function through which they will be used.
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>

#include "poc-data-operators.hpp"
#include "polygones.h" //for plg_*
#include "geo.h" //for geo_mercator_*
#include "functions.h" //for poc_getr*


/* see LAPACK's doxygen documentation on
  http://www.netlib.org/lapack/explore-html/index.html */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_getrf(int n,double *M,int *pivot,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
#if LAPACKC
#warning using LAPACKC
  dgetrf(n, n, M, n, pivot, &status);
#elif LAPACKF_
#warning using LAPACKF_
  dgetrf_ (&n, &n, M, &n, pivot, &status);
#else
  /* NOTE: LINPACK and ATLAS are deprecated */
  #error no linear algebra library available
#endif
  
  if(verbose>0 and status>0)
    STDERR_BASE_LINE_FUNC("dgetrf*() warning: factor %d is exactly zero\n",status-1);
  else if(verbose>=0 and status<0)
    STDERR_BASE_LINE_FUNC("dgetrf*() error: argument %d had an illegal value\n",-status);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_getrf(int n,complex<double> *M,int *pivot,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
#if LAPACKC
  zgetrf(n, n, M, n, pivot, &status);
#elif LAPACKF_
  zgetrf_(&n, &n, M, &n, pivot, &status);
#else
  /* NOTE: LINPACK and ATLAS are deprecated */
  #error no linear algebra library available
#endif
  
  if(verbose>0 and status>0)
    STDERR_BASE_LINE_FUNC("dgetrf*() warning: factor %d is exactly zero\n",status-1);
  else if(verbose>=0 and status<0)
    STDERR_BASE_LINE_FUNC("dgetrf*() error: argument %d had an illegal value\n",-status);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_getrs(int n,int nrhs,const double *A,const int *pivot,double *b,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  /* N : No transpose */
#if LAPACKC
  dgetrs('N',n,nrhs,A,n,pivot,b,n,&status);
#elif LAPACKF_
  dgetrs_("N",&n,&nrhs,A,&n,pivot,b,&n,&status);
#else
  /* NOTE: LINPACK and ATLAS are deprecated */
  #error no linear algebra library available
#endif
  
  if(verbose>=0 and status!=0)
    STDERR_BASE_LINE_FUNC("dgetrs*() error: argument %d had an illegal value\n",-status);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_getrs(int n,int nrhs,const complex<double> *A,const int *pivot,complex<double> *b,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  /* N : No transpose */
#if LAPACKC
  zgetrs('N',n,nrhs,A,n,pivot,b,n,&status);
#elif LAPACKF_ == 1
  zgetrs_("N",&n,&nrhs,A,&n,pivot,b,&n,&status);
#else
  /* NOTE: LINPACK and ATLAS are deprecated */
  #error no linear algebra library available
#endif
  
  if(verbose>=0 and status!=0)
    STDERR_BASE_LINE_FUNC("zgetrs*() error: argument %d had an illegal value\n",-status);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg2grd(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP(
    "polynomial grid from a 4-points polygone\n"
    "Syntax: "+funcName+"(plg&k&n)\n"
    "  plg  4 points polygone\n"
    "  k    (real) aspect-ratio factor, default:1\n"
    "  n    (complex) the grid size, default:9+9j\n"
    "The order of k and n in the input is irrelevant.\n"
    "  If any part of n is negative, the absolute value is taken and no Mercator projection is done. "
    "This is for when the Mercator projection has already been done.\n"
    );
  
  const int n=3,n2=n*n;
  int i,j,k,status;
  const complex<double> *d;
  
  int nx=0x9,ny=nx;
  double l=1.;
  bool doMerc=true;
  
  if(input.length<4){
    STDERR_BASE_LINE_FUNC("%d points\n",input.length);
    return;
    }
  
  for(i=4;i<input.length && i<6;i++){
    d=&input.data[i];
    
    if(imag(*d)==0){
      l=real(*d);
      }
    else{
      nx=fabs(real(*d));
      ny=fabs(imag(*d));
      doMerc= real(*d)>0. && imag(*d)>0.;
      }
    
    }
  
/*------------------------------------------------------------------------------
  vector from projection to POC-made Mercator */
  const complex<double> d0=input.data[0];
  geo_t merc=geo_mercator_init(real(d0),imag(d0),r2d);
  complex<double> bx[n+1];
  
#define debug_plg2grd 0
#if debug_plg2grd
  STDERR_BASE_LINE_FUNC("bx=[");
#endif
  for(i=0;i<n;i++){
    j=i+1;
    d=&input.data[j];
    if(doMerc){
      double rbxi,ibxi;
      geo_mercator_directe(merc,real(*d),imag(*d),&rbxi,&ibxi);
      bx[i]=complex<double>(rbxi,ibxi);
      }
    else
      bx[i]=*d-d0;
#if debug_plg2grd
    cerr << bx[i];
#endif
    }
#if debug_plg2grd
  cerr << "]\n";
#endif
  
  double dx,dy;
  dx=abs(bx[0]+bx[1]-bx[2]);
  dy=abs(bx[0]-bx[1]-bx[2]);
  l*=dx/dy;
  STDERR_BASE_LINE("%s:l*=(%g/%g=%g);l=%g\n",__func__,dx*.5,dy*.5,dx/dy,l);
  
  double area=0;
  for(j=1;j<n;j++){
    i=j-1;;
    area+=real(bx[i])*imag(bx[j])-imag(bx[i])*real(bx[j]);
    }
  if(area<0){
    l*=-1;
    STDERR_BASE_LINE("area=%g;l*=-1;l=%g\n",__func__,area,l);
    }
  
/*------------------------------------------------------------------------------
  matrix */
  int m,pivot[n];
  complex<double> A[n2];
  
  i=0;
  A[i++]=complex<double>(l,0);
  A[i++]=complex<double>(l,1);
  A[i++]=complex<double>(0,1);
  
  for(j=0;j<n;j++){
#if debug_plg2grd
    STDERR_BASE_LINE_FUNC("");
    cerr << "A[" << j<<"]=" << A[j];
#endif
    for(i=1;i<n;i++){
      m=i*n+j;
      A[m]=pow(A[j],i+1);
#if debug_plg2grd
      cerr << " A[" << m<<"]=" << A[m];
#endif
      }
#if debug_plg2grd
    cerr << '\n';
#endif
    }
  
  status=poc_getrf(n,A,pivot,1);
  
  if(status!=0)TRAP_ERR_RETURN(,1,"poc_getrf() had an error and returned %d\n",status);
  
/*------------------------------------------------------------------------------
  solution */
  status=poc_getrs(n,1,A,pivot,bx,1);
  
#if debug_plg2grd
  STDERR_BASE_LINE_FUNC("");
  for(j=0;j<n;j++){
    cerr << " b[" << j<<"]=" << bx[j];
    }
  cerr << '\n';
#endif
  
  if(status!=0)TRAP_ERR_RETURN(,1,"poc_getrs() had an error and returned %d\n",status);
  
/*------------------------------------------------------------------------------
  grid */
  if(ny<=1)
    ny=nx;
  
  output->info.dimensions << poc_dim_t("nj",ny) << poc_dim_t("ni",nx);
  output->init();
  
  complex<double> z,w;
  double rz=0.,iz=0.;
  
  dx=l/(nx-1);
  dy=1./(ny-1);
  
#if debug_plg2grd
  STDERR_BASE_LINE_FUNC("");
#endif
  for(j=0;j<ny;j++,iz+=dy){
    rz=0;
    m=j*nx;
    for(i=0;i<nx;i++,rz+=dx,m++){
      complex<double> *d=&output->data[m];
      w=0;
      for(k=0;k<n;k++){
        z=complex<double>(rz,iz);
        w+=bx[k]*pow(z,k+1);
        }
#if debug_plg2grd
      if((j==0 && i==nx-1)||(j==nx-1 && (i==0 || i==nx-1)))
        cerr << " w[" << i<<',' << j<<"]=" << w;
#endif
      if(doMerc){
        double rd,id;
        geo_mercator_inverse(merc,&rd,&id,real(w),imag(w));
        *d=complex<double>(rd,id);
        }
      else
        *d=w+d0;
      }
    }
#if debug_plg2grd
  cerr << '\n';
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void var2plg(const poc_cdata_t & var,plg_t *plg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  
  plg->init(var.length,PLG_INIT_SHARED);
  
  for(j=0;j<plg->npt;j++){
    complex<double> *z=&var.data[j];
    plg_move_point(*plg,j,real(*z),imag(*z));
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void resample(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* See also plg_resample() */
{
  POC_FORMULA_FUNC_HELP("resample a polygone to a given number of equally-spaced points\n"
    "  "+funcName+"(p,n)\n"
    "  p: given polygone\n"
    "  n: number of points\n"
    );
  
/*------------------------------------------------------------------------------
  input polygon */
  plg_t inPlg,outPlg;
  
  var2plg(input,&inPlg);
  
/*------------------------------------------------------------------------------
  argument */
  string name;
  const poc_att_t *att;
  att=input.info.attributes.findP("");
  if(att!=0)
    name=att->as_string();
  
  if(att!=0){
    char *c=&name[name.size()-1];
    
    const poc_cdata_t
      *var=poc_formula_check_var(vars,name);
    const complex<double>
      *var0=&var->data[0];
    
    if(isfinite(*var0))
      outPlg.npt=real(*var0);
    
    if(*c>'a')
      STDERR_BASE_LINE_FUNC("%s():NOTE:all %d last arguments are ignored\n",__func__,*c-'a');
    
    if(outPlg.npt<=0)
      TRAP_ERR_RETURN(,1,"%s():number of output point is %d<=0\n",__func__,outPlg.npt);
    
    }
  else
    TRAP_ERR_RETURN(,1,"%s():number of output point must be given!\n",__func__);
  
/*----------------------------------------------------------------------------*/
  outPlg=plg_resample(inPlg,outPlg.npt);
  
  asprintf(name="","n%d",outPlg.npt);
  output->info << poc_dim_t(name,outPlg.npt);
  output->init();
  
  for(int i=0;i<outPlg.npt;i++){
    output->data[i]=complex<double>(outPlg.t[i],outPlg.p[i]);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator << (const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// distort a grid with a border polygone
/**
\param d1 grid
\param d2 border polygone
*/
/*----------------------------------------------------------------------------*/
{
  poc_cdata_t d0;
  
  int i,j,l,m,intersections;
  
  const poc_list_t<poc_dim_t> &d1ds=d1.info.dimensions;
  
  m=d1ds.size();
  if(m!=2) TRAP_ERR_RETURN(d0,1,"%s:%d dimensions when 2 are expected\n",__func__,m);
  m=d2.info.dimensions.size();
  if(m!=1) TRAP_ERR_RETURN(d0,1,"%s:%d dimensions when 1 is expected\n",__func__,m);
  
  plg_t border,journey;
  
  border.init(d2.length,PLG_INIT_SHARED);
  for(m=0;m<d2.length;m++)
    plg_move_point(border,m,real(d2.data[m]),imag(d2.data[m]));
  
  journey.init(2,PLG_INIT_SHARED);
  
/*------------------------------------------------------------------------------
  mask */
  bool *mask,mm=false;
#define USE_nmasked 0
#if USE_nmasked
  int nmasked=0;
#endif
  
  plg_move_point(journey,0,real(d1.data[0]),imag(d1.data[0]));
  
  mask=new bool[d1.length];
  mask[0]=mm;
  
  for(m=1;m<d1.length;m++){
    const complex<double> *z=&d1.data[m];
    plg_move_point(journey,1,real(*z),imag(*z));
    
    intersections=plg_secantpoints(journey,border);
    mm^=intersections & 1;
    
    mask[m]=mm;
    
#if USE_nmasked
    if(mm)
      nmasked++;
#endif
    
    plg_copy_point(&journey,0,journey,1);
    }
  
  journey.destroy();
  
#if USE_nmasked
  mm=nmasked>(d1.length/2);
  STDERR_BASE_LINE_FUNC("%d/%d on one side :%s swiching.\n",nmasked,d1.length,mm?"":" NOT");
#else
  mm=mask[0]==mask[d1.length-1];
  STDERR_BASE_LINE_FUNC("first and last corners on %s :%s swiching.\n",mm?"same side":"either sides",mm?"":" NOT");
#endif
  if(mm)
    for(m=0;m<d1.length;m++)
      mask[m]^=1;
  
  const int
    nx=d1ds[1].len,
    ny=d1ds[0].len;
  
  d0.info.dimensions=d1ds;
  d0.init();
  const string
    formula=getNameOrFormula(d1.info)+" << "+getNameOrFormula(d2.info);
  d0.info << poc_att_t("formula",formula);
  
  plg_t line;
  double *lineIs=0,dl,p;
  int fp,cp;
  
/*------------------------------------------------------------------------------
  vertical */
  line.init(nx,PLG_INIT_SHARED);
  
  for(j=0;j<ny;j++){
    l=j*nx;
    m=l+nx-1;
    
    m=l;
    
    if(!mask[l] && !mask[m]){
      for(i=0;i<nx;i++,m++){
        d0.data[m]=d1.data[m];
        }
      continue;
      }
    
    for(i=0;i<nx;i++,m++){
      const complex<double> *z=&d1.data[m];
      plg_move_point(line,i,real(*z),imag(*z));
      }
    
    deletep(&lineIs);
    intersections=plg_secantpoints(border,line,0,0,0,&lineIs);
    
    if(intersections==2);
    else if(intersections==1){
      double lineIs0=lineIs[0];
      delete[]lineIs;
      lineIs=new double[2];
      lineIs[0]=lineIs0;
      lineIs[1]=mask[l]?nx-1:0;
      }
    else continue;
    
    if(lineIs[1]<lineIs[0]) swapValues(&lineIs[1],&lineIs[0]);
    
    dl=(lineIs[1]-lineIs[0])/(nx-1);
    
    p=lineIs[0];
    for(i=0;i<nx;i++,p+=dl){
      fp=floor(p);
      cp=fp+1;
      
      m=l+i;
      complex<double> *z=&d0.data[m];
      double rz,iz;
      
      rz=line.x[fp]*(cp-p)+line.x[cp]*(p-fp);
      iz=line.y[fp]*(cp-p)+line.y[cp]*(p-fp);
      
      *z=complex<double>(rz,iz);
      }
    }
  
  deletep(&lineIs);
  delete[]mask;
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_flag(bool nn[8])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// compute flag from neighbour states
/**
\param nn array[8] of neighbour states
For more information, see
- TEST_compute_flag() (activated with #DO_TEST_compute_flag )
- and orthogonalise()
*/
/*----------------------------------------------------------------------------*/
{
  int j;
  
  for(j=0;j<8 && !nn[j];j++);
  if(j>=8)
    return 0;
  
  int i,l,borders,flag=-1;
  bool nnl[3],isCorner,isVertical;
  
  isCorner=false;
  borders=0;
  
  /* corners */
  for(i=0;i<8 && !isCorner;i+=2){
    for(j=0;j<3;j++){
      l=i+j-1;
      while(l<0)l+=8;
      while(l>8)l-=8;
      nnl[j]=nn[l];
      }
    isCorner=
      (nnl[1] && !nnl[0] && !nnl[2]) ||
      (nnl[0] && nnl[2]);
    }
  
  /* borders */
  if(!isCorner) for(i=1;i<8;i+=2){
    for(j=0;j<3;j++){
      l=i+2*(j-1);
      while(l<0)l+=8;
      while(l>8)l-=8;
      nnl[j]=nn[l];
      }
    if(!nnl[0] && nnl[1] && !nnl[2]){
      borders++;
      isVertical= i==3 || i==7;
      }
    }
  
  if(isCorner){
    flag=1;
    }
  else if(borders>0){
    flag=2*borders+isVertical;
    }
  else TRAP_ERR_EXIT(ENOEXEC,"coding error\n");
  
  return flag;
}


#define DO_TEST_compute_flag 0


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#if DO_TEST_compute_flag
  void TEST_compute_flag()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  bool nn[8];/* NAN neighbours */
  const int nflags=6;
  const char
    *nanClasses[2]={"false","true"},
    *flagClasses[nflags]={"inside","corner","h1","v1","h2","v2"},
    *flagColors[nflags]={"fff","333","ff0","0ff","0f0","00f"};
  int flag,flagDist[nflags];
  aset(flagDist,nflags,0);
  
  FILE *f=fopen("orthogonalise-test.html","w");
  fprintf(f,"<style><!-- table{display:inline} td{width:7px;height:9px}\n"
    ".false{background-color:#aaa}\n"
    ".true{background-color:#000}\n");
  for(flag=0;flag<nflags;flag++)
    fprintf(f,".%s{background-color:#%s}\n",flagClasses[flag],flagColors[flag]);
  fprintf(f,"--></style>\n");
  for(int ii=0;ii<256;ii++){
    int i,j,l;
    
    if(ii>0 && (ii%16==0))
      fprintf(f,"<br>\n");
    
    for(i=0;i<8;i++)
      nn[i]=(ii>>i)&1;
    
    flag=compute_flag(nn);
    
    if(flag<0 || nflags<=flag)
      fprintf(f,"%d:%d\n",ii,flag);
    else{
      fprintf(f,"<table title=%d:%d><tr><td class=%s><td class=%s><td class=%s>"
        "<tr><td class=%s><td class=%s><td class=%s>"
        "<tr><td class=%s><td class=%s><td class=%s></table>\n",ii,flag,
        nanClasses[nn[0]],nanClasses[nn[1]],nanClasses[nn[2]],
        nanClasses[nn[7]],flagClasses[flag],nanClasses[nn[3]],
        nanClasses[nn[6]],nanClasses[nn[5]],nanClasses[nn[4]]);
      flagDist[flag]++;
      }
    }
  for(flag=0;flag<nflags;flag++)
    fprintf(f,"<br>%d flag=%d\n",flagDist[flag],flag);
  fclose(f);
  TRAP_ERR_EXIT(ENOEXEC,"testing\n");
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void orthogonalise(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP(
    "make a grid more orthogonal\n"
    "  Set `iterations' variable to change the number of iterations, 256 by default.\n"
    );
  
  const poc_list_t<poc_dim_t> &dimensions=input.info.dimensions;
  if(dimensions.size()!=2) TRAP_ERR_RETURN(,1,"%s:%d dimensions when 2 are expected\n",__func__,dimensions.size());
  
  output->info.dimensions=dimensions;
  output->init();
  
  const int
    nx=dimensions[1].len,
    ny=dimensions[0].len;
  
  vector2D_t *v0,*v;
  v0=new vector2D_t[input.length];
  v =new vector2D_t[input.length];
  
/*------------------------------------------------------------------------------
  initialisation */
  
  /* conformal projection */
  const complex<double> *z0=&input.data[nx/2+nx*(ny/2)];
  geo_t merc=geo_mercator_init(real(*z0),imag(*z0),r2d);
  
  for(int m=0;m<input.length;m++){
    vector2D_t *v0m=&v0[m];
    const complex<double> *z=&input.data[m];
    geo_mercator_directe(merc,real(*z),imag(*z),&v0m->x,&v0m->y);
    }
  
#if DO_TEST_compute_flag
  TEST_compute_flag();
#endif
  
/*------------------------------------------------------------------------------
  iterations */
  poc_cdata_t *iterationsVar=&poc_formula_get_var(vars,"iterations");
  if(iterationsVar->length<=0 || real(iterationsVar->data[0])<1)
    *iterationsVar=0x100;
  const int iterations=real(iterationsVar->data[0]);
  
  for(int ii=0;ii<iterations;ii++){
    //#pragma omp parallel for
    for(int j=0;j<ny;j++){
      int i,m=j*nx,flag;
      
      bool nn[8];/* NAN neighbours */
      int di,dj,il,jl,l;
      
      vector2D_t v1,v2,vi;
      range_t<double> vpR(0.,+INFINITY);
      double vp,vp2, hh,hv,k;
      
      for(i=0;i<nx;i++,m++){
        int assert=0;
        //if(ii==0) assert=(6<j && j<11)|(i==4 && j==8) << 1;
        if(assert) STDERR_BASE_LINE_FUNC("[%d,%d];",i,j);
        
        const vector2D_t *v0m=&v0[m];
        vector2D_t *vm=&v[m];
        if(isnan(v0m->x) || isnan(v0m->y)){
          if(ii==0) *vm=*v0m;
          if(assert) STDERR_LINE("\n");
          continue;
          }
        
        aset(nn,8,true);
        for(dj=-1;dj<=1;dj++){
          jl=j+dj;
          if( jl<0 || ny<=jl ) continue;
          
          for(di=-1;di<=1;di++){
            il=i+di;
            if( il<0 || nx<=il ) continue;
            
            switch(dj){
            case -1:
              l=1+di;
              break;
            case 0:
              if(di==0) continue;
              l=5-2*di;
              break;
            case 1:
              l=5-di;
              break;
            default: TRAP_ERR_EXIT(ENOEXEC,"coding error\n");
              }
            
            const vector2D_t *v0l=&v0[m+dj*nx+di];
            nn[l]=isnan(v0l->x) || isnan(v0l->y);
            }
          }
        
        flag=compute_flag(nn);
        if(assert) {STDERR_LINE("%d(",flag);for(int i=0;i<8;i++)fputc(nn[i]?'1':'0',stderr);fprintf(stderr,");");}
        
        switch(flag){
        case 0:
          /* inside */
          *vm=0.5* *v0m;
#if 0
          /* simple algorithm */
          *vm+=0.125*v0[m+1];
          *vm+=0.125*v0[m-1];
          *vm+=0.125*v0[m+nx];
          *vm+=0.125*v0[m-nx];
#else
          /* anti-collapses algorithm */
          hh=hypot(v0[m+1]-v0[m-1]);
          hv=hypot(v0[m+nx]-v0[m-nx]);
          k=.25/(hh+hv);
          hh*=k;
          hv*=k;
          *vm+=hv*v0[m+1];
          *vm+=hv*v0[m-1];
          *vm+=hh*v0[m+nx];
          *vm+=hh*v0[m-nx];
#endif
          if(assert) STDERR_LINE("(%g,%g)->(%g,%g)\n",v0m->x,v0m->y,vm->x,vm->y);
          continue;
        
        case 4:
        case 5:
          /* double edges */
        case 1:
          /* corners */
          *vm=*v0m;/* do not move */
          if(assert) STDERR_LINE("(%g,%g)->(%g,%g)\n",v0m->x,v0m->y,vm->x,vm->y);
          continue;
        
        case 2:
        case 3:
          /* single edges */
          switch(flag){
          case 2:
            v1=v0[m-1];
            v2=v0[m+1];
            if(nn[1])
              vi=v0[m+nx];
            else
              vi=v0[m-nx];
            break;
          case 3:
            v1=v0[m-nx];
            v2=v0[m+nx];
            if(nn[7])
              vi=v0[m+1];
            else
              vi=v0[m-1];
            break;
          default:
            TRAP_ERR_EXIT(ENOEXEC,"coding error\n");
            }
          
          v1-=*v0m;
          v2-=*v0m;
          vi-=*v0m;
          
          vpR.max=hypot(v1);
          v1/=vpR.max;
          vp=v1%vi;
          vpR.constrain(&vp);
          
          vpR.max=hypot(v2);
          v2/=vpR.max;
          vp2=v2%vi;
          vpR.constrain(&vp2);
          
          vp-=vp2;
          vp*=.5;
          
          *vm=*v0m;
          if(vp<0)
            *vm-=vp*v2;
          if(vp>0)
            *vm+=vp*v1;
          if(assert) STDERR_LINE("(%g,%g)->(%g,%g)\n",v0m->x,v0m->y,vm->x,vm->y);
          continue;
          }
        
        TRAP_ERR_EXIT(ENOEXEC,"programming error in %s: flag=%d\n",__func__,flag);
        }
      }
    
    swapValues(&v,&v0);
    }
  
/*------------------------------------------------------------------------------
  conclusion */
  for(int m=0;m<input.length;m++){
    vector2D_t *v0m=&v0[m];
    complex<double> *z=&output->data[m];
    double rz,iz;
    geo_mercator_inverse(merc,&rz,&iz,v0m->x,v0m->y);
    *z=complex<double>(rz,iz);
    }
  
  delete[]v0;
  delete[]v;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void crange(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP(
    "make a matrix from a range of complex\n"
    "Syntax 1: "+funcName+"(endVal)\n"
    "Syntax 2: "+funcName+"(startVal&endVal)\n"
    "Syntax 3: "+funcName+"(startVal&endVal&n)\n"
    "The first value is startVal or 0. The last value is endVal.\n"
    "We take here dz = endVal-startVal\n"
    "If not given, n = dz+1+1j\n"
    "If real(n)<0, n is taken as the increment.\n"
    "Otherwise the size is imag(n)*real(n).\n"
    );
  
  int i,j,m;
  
/*------------------------------------------------------------------------------
  interpret input */
  int nx,ny;
  complex<double> z0(0.,0.),z1(1.,1.);
  if(input.length==1){
    z1=input.data[0];
    }
  else if(input.length>1){
    z0=input.data[0];
    z1=input.data[1];
    }
  complex<double> dz,n;
  
  dz=z1-z0;
  
  if(input.length>=3){
    n=input.data[2];
    }
  else{
    n=dz+complex<double>(1.,1.);
    }
  
  if(real(n)<0.){
    if(imag(n)==0.)
      n=complex<double>(real(n),real(n));
    
    nx=1+real(dz)/fabs(real(n));
    ny=1+imag(dz)/fabs(imag(n));
    dz=complex<double>(fabs(real(n)),fabs(imag(n)));
    }
  else{
    nx=max(1,(int)real(n));
    ny=max(1,(int)imag(n));
    dz=complex<double>( real(dz)/max(1,nx-1) , imag(dz)/max(1,ny-1) );
    }
  
/*------------------------------------------------------------------------------
  build matrix */
  if(nx>1 && ny>1)
    output->info.dimensions << poc_dim_t("nj",ny) << poc_dim_t("ni",nx);
  else
    output->info.dimensions << poc_dim_t(max(ny,nx));
  output->init();
  
  m=0;
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++,m++){
      output->data[m]=z0+complex<double>(real(dz)*i,imag(dz)*j);
      }
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fmerc(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("forward Mercator projection\n");
  
  geo_t merc=geo_mercator_init(0.,0.,r2d);
  
  output->info.dimensions=input.info.dimensions;
  output->init();
  
  for(int m=0;m<input.length;m++){
    const complex<double> *zi=&input.data[m];
    complex<double> *zo=&output->data[m];
    double rzo,izo;
    geo_mercator_directe(merc,real(*zi),imag(*zi),&rzo,&izo);
    *zo=complex<double>(rzo,izo);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void imerc(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("inverse Mercator projection\n");
  
  geo_t merc=geo_mercator_init(0.,0.,r2d);
  
  output->info.dimensions=input.info.dimensions;
  output->init();
  
  for(int m=0;m<input.length;m++){
    const complex<double> *zi=&input.data[m];
    complex<double> *zo=&output->data[m];
    double rzo,izo;
    geo_mercator_inverse(merc,&rzo,&izo,real(*zi),imag(*zi));
    *zo=complex<double>(rzo,izo);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void geokm(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("polygone length\n");
  
  const int dl=input.info.dimensions.size();
  if(dl!=1) TRAP_ERR_RETURN(,1,"%s:%d dimensions when 1 is expected\n",__func__,dl);
  
  double total=0.,distance;
  complex<double> *zi0=0,*zi1=0;
  int m;
  
  for(m=0;m<input.length;m++){
    zi0=zi1;
    zi1=&input.data[m];
    
    if(zi0==0) continue;
    if(isnan(*zi1)) TRAP_ERR_RETURN(,1,"%s:nan value found: "+input.info.name+"[%d]\n",__func__,m);
    
    distance=geo_km(zi0->real(),zi0->imag(),zi1->real(),zi1->imag());
    total+=distance;
    }
  
  *output=total;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fliph(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("flip horizontally or flip 1D data\n");
  
  const poc_list_t<poc_dim_t> &dimensions=input.info.dimensions;
  const int dl=dimensions.size();
  if(dl!=2 && dl!=1) TRAP_ERR_RETURN(,1,"%s:%d dimensions when 1 or 2 are expected\n",__func__,dl);
  
  int nx,ny;
  get_poc_cdata_dims(input,&nx,&ny);
  
  output->info.dimensions=dimensions;
  output->init();
  
  int i,j,k=0,l;
  
  for(j=0;j<ny;j++){
    l=j*nx+nx-1;
    for(i=0;i<nx;i++,k++,l--){
      output->data[l]=input.data[k];
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void flipv(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("flip vertically\n");
  
  const poc_list_t<poc_dim_t> &dimensions=input.info.dimensions;
  const int dl=dimensions.size();
  if(dl!=2) TRAP_ERR_RETURN(,1,"%s:%d dimensions when 2 are expected\n",__func__,dl);
  
  const int
    ny=dimensions[0].len,
    nx=dimensions[1].len;
  
  output->info.dimensions=dimensions;
  output->init();
  
  int i,j,k=0,l;
  
  for(j=0;j<ny;j++){
    l=(ny-1-j)*nx;
    for(i=0;i<nx;i++,k++,l++){
      output->data[l]=input.data[k];
      }
    }
}
