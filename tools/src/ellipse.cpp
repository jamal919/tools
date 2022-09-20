
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\date 2011-11-29 Damien ALLAIN : creation

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Calculates current ellipses

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>

#include "tools-structures.h"

#include "functions.h" //fct_echo
#include "poc-netcdf-data.hpp"
#include "map.h"
#include "geo.h"
#include "tides.h"
#include "poc-time.h"

#include "zapper.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : as[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s grid_file gridded_varname atlas_convention lon_speed_varname lat_speed_varname wave1 [wave2 ...]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Calculates current ellipses.\n"
    "  The first argument is the grid file.\n"
    "  The following argument is a variable name at whose points the speed components will be interpolated. So that can be anything (bathymetry, temperature, ...).\n"
    "  The 3rd argument is the atlas file name convention. See below.\n"
    "  The following 2 arguments are the variable names for the longitudinal and and latitudinal speed components.\n"
    "  The following arguments are the waves you want to analyse.\n"
    "  For each wave, it will calculate current ellipses.\n");
  print_decode_atlasname_help();
  printf("\n"
    "OPTION :\n"
    "  -h,--help  Show this help and exit.\n"
    "  -t  for testing\n"
    ); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int aC, char *as[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief main function

Document this to show how this works.
See the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/
{
  int aI,wI=0,status;//argument index, wave index, status

  char *a;//argument
  char *aPF=NULL,*aP,*bP,*bN;//atlas path convention,atlas path,bathymetry path and variable name
  const int vC=2;//number of variables
  char *vNs[vC]/*,*&uN=vNs[0],*&vN=vNs[1]*/;//names of input variables
  char **wNs;//list of wave names
  spectrum_t ws;//list of waves

  grid_t grid;
  int gridn;//<grid size and value indexes
  int vI;//variable index

  string cmd=fct_echo(aC,as);

  ws.n=aC-4-vC;
  if(ws.n<1){
    if(aC<2 ||
      strcmp(as[1],"--help")==0 ||
      strcmp(as[1],"-h")==0){
      print_help(as[0]);
      exit(0);
      }
    if(strcmp(as[1],"-t")==0){
      double M,a,d,m,p;
      int i,j,k,l;
      for(i=-1;i<=1;i++)
      for(j=-1;j<=1;j++)
      for(k=-1;k<=1;k++)
      for(l=-1;l<=1;l++){
        const complex<double> u(i,j),v(k,l);
        STDERR_BASE_LINE("[%d] u=(%g;%g) v=(%g;%g) :\n",i,real(u),imag(u),real(v),imag(v));
        ellipse_t ellipse=ellipse_parameter(u,v);
        STDERR_BASE_LINE("%g %g %g %g %g\n",ellipse.a,ellipse.b,ellipse.phase*r2d,ellipse.inclination*r2d,ellipse.polarisation);
        ellipse_Madmp(u,v,&M,&a,&d,&m,&p);
        STDERR_BASE_LINE("%g %g %g %g %g\n",M,m,a,d,p);
        M=-1.;
        a=d=NAN;
        for(double ia=-179.;ia<=180.;ia++){
          const complex<double> ca=polar(1.,ia*d2r),uv(real(v*ca),real(u*ca));
          const double m_=abs(uv);
          if(M<m_){
            a=ia;
            M=m_;
            d=arg(uv)*r2d;
            }
          }
        if(d<0.)d+=360.;
        STDERR_BASE_LINE("%g %g %g %g %g\n",M,m,a,d,p);
        STDERR_BASE_LINE("%g %g %g %g %g\n",M,m,degree_recale(a+180.,0.),degree_recale(d+180.,180.),p);
        }
      return ENOEXEC;
      }
    STDERR_BASE_LINE("Missing %d arguments\n",1-ws.n);
    if(ws.n<0)
      fprintf(stderr,"*** Please provide a grid path and a gridDED variable name, an atlas file name convention, 2 variable names and a list of waves. ***\n");
    else
      fprintf(stderr,"*** Please provide a list of waves. ***\n");
    print_help(as[0]);
    exit(-1);
    }
  wNs=new char*[ws.n+1];
  for(aI=1;aI<aC;aI++){
    a=strdup(as[aI]);//DO NOT FREE a !!!
    if(a[0]=='-'){
      if(strcmp(a,"--help")==0 ||
        strcmp(a,"-h")==0){
        print_help(as[0]);
        exit(0);
        }
      //else
      printf("unknown option %s\n",a);
      print_help(as[0]);
      exit(-1);
      }
    switch(aI){
    case 1:
      bP=a;
      break;
    case 2:
      bN=a;
      break;
    case 3:
      aPF=a;
      break;
    case 4:
    case 5:
      vNs[aI-4]=a;
      break;
    default:
      wNs[wI]=a;
      wI++;
      }
    }
  wNs[wI]=NULL;
  
  ///It initialises the list of waves
  ws.init(initialize_tide(),wNs);
  printf("%d valid wave names given.\n",ws.n);
  
/*---------------------------------------------------------------------*//**<h1>
  reads the grid </h1>*/
  printf("Grid from %s gridDED variable %s : reading...",bP,bN);fflush(stdout);
  status=poc_get_grid(bP,bN,&grid);
  NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\"%s\",\"%s\",) error",bP,bN);
  gridn=grid.nx*grid.ny;
#if 0
  printf(" done:\n");
  grid.print(cout,__FILE__ ":" TOSTRING(__LINE__) ":");
#else
  printf(" done.\n");
#endif
  
/*---------------------------------------------------------------------*//**<h1>
  carry out ellipse calculations for each wave </h1>*/
  for(wI=0;wI<ws.n;wI++){//for each wave
    int verbose=0;
    poc_var_t gridVar,var;
    const int nparts=2;
    const string
#if 1
      pNs[nparts]={AMPSUF,PHASUF},//<part names
#else
      pNs[nparts]={"a","g"},//<part names
#endif
      pSNs[nparts]={"tidal_amplitude","tidal_phase"};//<part standard names
    const string oNs[vC]={"U","V"},//<output variable names
      oSNs[vC]={"sea_water_x_velocity","sea_water_y_velocity"};//<output standard names. See #standard_names_URL
    complex<double> *vars[vC],*&u=vars[0],*&v=vars[1];
    double *c=0,*s=0;
    
    string oP=ws.waves[wI].name+string("-ellipses.nc");//output path
    double *M,*a,*d,*m,*p;
    ellipse_t ellipse;
    
    status=poc_save_grid(oP,&gridVar,grid,1,verbose);
    status=poc_def_att(oP,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
    status=poc_def_att(oP,poc_att_t("history",cmd));
    
    /// <h2>read speed components from the relevant atlases</h2>
    for(vI=0;vI<vC;vI++){//for each component
      grid_t vgrid;
      complex<double> *buf,*varsVi;
      string ampName=vNs[vI]+pNs[0];
      bool getRotation;
      
      aP=decode_atlasname(aPF,ws.waves[wI].name,vNs[vI]);
      printf("Wave %3s constants for %3s in %s: ",ws.waves[wI].name,vNs[vI],aP);fflush(stdout);
      
      printf("grid, ");fflush(stdout);
      status=poc_get_grid(aP,ampName,&vgrid,1);
      NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\"%s\",\""+ampName+"\",) error",aP);
      
      set_grid_list(&vgrid,1,0);
      
      printf("read, ");fflush(stdout);
      buf=new complex<double>[vgrid.nx*vgrid.ny];
      status=poc_get_cvara(aP,ampName,vNs[vI]+pNs[1],-1,buf,1);
      
      printf("interpolation, ");fflush(stdout);
      varsVi=new complex<double>[gridn];
      vars[vI]=varsVi;
      getRotation= vgrid.modeH==2 and vI==0 ;
      if(getRotation){
        c=new double[gridn];
        s=new double[gridn];
        }
#pragma omp parallel
      {
      int i,j,m, vi0,vi,vj;
      int64_t vm;
      complex<double> value;
      double lon,lat,clat,dLon,dLat,dih;
#pragma omp for
      for(j=0;j<grid.ny;j++){
        vm=-1;
        m=j*grid.nx;
        
        for(i=0;i<grid.nx;i++,m++){
          
          grid.xy(m,lon,lat);
          
          index_interpolation(vgrid,lon,lat,&vm,buf,NC_FILL_COMPLEX,&value,verbose);
          varsVi[m]=value;
          
          if(not getRotation)
            continue;
          if(vm<0){
            c[m]=s[m]=NC_FILL_DOUBLE;
            continue;
            }
          
          clat=cos(lat*d2r);
          vgrid.ij(vm,&vi0,&vj);
          
          vi=max(0,vi0-1);
          vgrid.xy(vi,vj,lon,lat);
          vi=min(vgrid.nx,vi0+1);
          vgrid.xy(vi,vj,dLon,dLat);
          
          dLon-=lon;
          dLon*=clat;
          dLat-=lat;
          
          dih=1./hypot(dLon,dLat);
          c[m]=dLon*dih;
          s[m]=dLat*dih;
          
          }
        
        }/* EO omp for */
      }/* EO omp parallel */
      delete[]buf;
      vgrid.free();
      
      printf("write, ");fflush(stdout);
      poc_var_t av,gv;
      av=gridVar;
      av.init(oNs[vI]+pNs[0],NC_DOUBLE,oSNs[vI]+"_"+pSNs[0],"m s-1",NAN,long_name_from_varname());
      gv=gridVar;
      gv.init(oNs[vI]+pNs[1],NC_DOUBLE,oSNs[vI]+"_"+pSNs[1],"degrees",NAN,long_name_from_varname());
      status=poc_put_cvara(oP,av,gv,0,vars[vI],1);
      
      delete[]aP;
      printf("done.\n");
      }
    
    if(c!=0 and s!=0){
      /* rotate */
      printf("Rotating, ");fflush(stdout);
      
      complex<double> ui0,vi0;
      
      for(int i=0;i<gridn;i++){
        complex<double> *ui=&u[i],*vi=&v[i];
        const double ci=c[i],si=s[i];
        
        if(*ui==NC_FILL_COMPLEX or *vi==NC_FILL_COMPLEX){
          *ui=*vi=NC_FILL_COMPLEX;
          continue;
          }
        
        ui0=*ui;
        vi0=*vi;
        
        *ui=ui0*ci-vi0*si;
        *vi=ui0*si+vi0*ci;
        }
      
      printf("done. ");fflush(stdout);
      
      printf("Saving rotation to "+oP+" : ");fflush(stdout);
      var=gridVar;
      
      var.init("c",NC_DOUBLE,"cosine_of_grid_rotation","",NAN,"cos(grid_rotation)");
      printf(var.name+", ");fflush(stdout);
      poc_put_vara(oP,var,0,c,1);
      
      var.init("s",NC_DOUBLE,"sine_of_grid_rotation","",NAN,"sin(grid_rotation)");
      printf(var.name+", ");fflush(stdout);
      poc_put_vara(oP,var,0,s,1);
     
      const string roNs[vC]={"E","N"},//<output variable names
        roSNs[vC]={"eastward_sea_water_velocity","northward_sea_water_velocity"};//<output standard names. See #standard_names_URL
      
      for(vI=0;vI<vC;vI++){
        poc_var_t av,gv;
        av=gridVar;
        av.init(roNs[vI]+pNs[0],NC_DOUBLE,roSNs[vI]+"_"+pSNs[0],"m s-1",NAN,long_name_from_varname());
        gv=gridVar;
        gv.init(roNs[vI]+pNs[1],NC_DOUBLE,roSNs[vI]+"_"+pSNs[1],"degrees",NAN,long_name_from_varname());
        printf("("+av.name+","+gv.name+"), ");fflush(stdout);
        status=poc_put_cvara(oP,av,gv,0,vars[vI],1);
        }
      
      printf("done.\n");
      }
    
    /// <h2>call ellipse_Madmp() for each node</h2>
    
    printf("Allocating buffers, ");fflush(stdout);
    M=new double[gridn];
    a=new double[gridn];
    d=new double[gridn];
    m=new double[gridn];
    p=new double[gridn];
    
    printf("computing ellipses, ");fflush(stdout);
    for(int i=0;i<gridn;i++){
      const complex<double> *ui=&u[i],*vi=&v[i];
      if(*ui==NC_FILL_COMPLEX || *vi==NC_FILL_COMPLEX){
        M[i]=a[i]=d[i]=m[i]=p[i]=NC_FILL_DOUBLE;
        continue;
        }
      ellipse_Madmp(*ui,*vi,&M[i],&a[i],&d[i],&m[i],&p[i]);
      }
    printf("done.\n");
    
    /// <h2>save</h2>
    
    printf("Saving to "+oP+" : ");fflush(stdout);
    var=gridVar;
    
    var.init("M",NC_DOUBLE,"sea_water_velocity_tidal_maximum","m s-1",NAN,"maximum_velocity");
    printf(var.name+", ");fflush(stdout);
    poc_put_vara(oP,var,0,M,1);
    
    var.init("a",NC_DOUBLE,"sea_water_velocity_tidal_phase_of_maximum","degrees",NAN,"maximum_velocity_phase");
    printf(var.name+", ");fflush(stdout);
    poc_put_vara(oP,var,0,a,1);
    
    var.init("d",NC_DOUBLE,"sea_water_velocity_tidal_direction_of_maximum","degrees",NAN,"maximum_velocity_direction");
    printf(var.name+", ");fflush(stdout);
    poc_put_vara(oP,var,0,d,1);
    
    var.init("m",NC_DOUBLE,"sea_water_velocity_tidal_minimum","m s-1",NAN,"minimum_velocity");
    printf(var.name+", ");fflush(stdout);
    poc_put_vara(oP,var,0,m,1);
    
    var.init("p",NC_DOUBLE,"sea_water_velocity_tidal_clockwise_polarisation","",NAN,"clockwise_polarisation");
    printf(var.name+", ");fflush(stdout);
    poc_put_vara(oP,var,0,p,1);
    
    printf("done.\n");
    
    delete[]M;
    delete[]a;
    delete[]d;
    delete[]m;
    }//EO for each wave
  
  return 0;
}
