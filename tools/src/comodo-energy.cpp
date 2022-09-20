
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
\brief Calculates the spectral energy budget from structured-grid tidal atlases

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
#include "maths.h" //operator %
#include "constants.h" //MeanEarthRadius
#include "poc-netcdf-data.hpp"
#include "map.h"
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
    "  %s bathymetry bathymetry_varname atlas_convention elevation_varname lon_speed_varname lat_speed_varname wave1 [wave2 ...]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Calculates the spectral energy budget from structured-grid tidal atlases.\n"
    "  The first argument is the bathymetry.\n"
    "  The following argument is the bathymetry variable name.\n"
    "  The 3rd argument is the atlas file name convention. See below.\n"
    "  The following 3 arguments are the variable names for the elevation and the longitudinal and and latitudinal speed components.\n"
    "  The following arguments are the waves you want to analyse.\n"
    "  For each wave, it will calculate Stokes'transport, energy flux and dissipation rate. The elevations and speed components will be interpolated at bathymetry points.\n");
  print_decode_atlasname_help();
  printf("\n"
    "OPTION :\n"
    "  -h,--help  Show this help and exit.\n"
    );
  print_OPENMP_help(prog_name); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void gradient(const grid_t &grid,T *z,T mask,T *dzdx,T *dzdy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///gradient of structured data
/*----------------------------------------------------------------------------*/
{
  int j,gridn=grid.nx*grid.ny;
  
  aset(dzdx,gridn,mask);
  aset(dzdy,gridn,mask);
  
  const int gridnx=grid.nx;/* because of icc pure madness */
  #pragma omp parallel for
  for(j=grid.nx;j<gridn-grid.nx;j+=gridnx){
    int i;
    
    for(i=1;i<grid.nx-1;i++){
      vector3_t d;
      
      if(z[j+i-1]==mask
       ||z[j+i+1]==mask){
        }
      else{
        d=math_polar2cartesian(grid.x[j+i-1]*d2r,grid.y[j+i-1]*d2r,MeanEarthRadius)-
          math_polar2cartesian(grid.x[j+i+1]*d2r,grid.y[j+i+1]*d2r,MeanEarthRadius);
        dzdx[j+i]=(z[j+i-1]-z[j+i+1])/hypot(d);
        }
      
      if(z[j+i-grid.nx]==mask
       ||z[j+i+grid.nx]==mask){
        }
      else{
        d=math_polar2cartesian(grid.x[j+i-grid.nx]*d2r,grid.y[j+i-grid.nx]*d2r,MeanEarthRadius)-
          math_polar2cartesian(grid.x[j+i+grid.nx]*d2r,grid.y[j+i+grid.nx]*d2r,MeanEarthRadius);
        dzdy[j+i]=(z[j+i-grid.nx]-z[j+i+grid.nx])/hypot(d);
        }
      
      }
    }
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
  int aI,wI=0,M2wI=-1,status;//argument index, wave index, M2 wave index, status

  char *a;//argument
  char *aPF=NULL,*aP,*hP=NULL,*hN=NULL;//atlas path convention,atlas path,bathymetry path and variable name
  const int vC=3;//number of variables
  char *vNs[vC]/*,*&eN=vNs[0],*&uN=vNs[1],*&vN=vNs[2]*/;//names of input variables
  char **wNs;//list of wave names
  spectrum_t ws;//list of waves

  grid_t grid;
  string oP;
  poc_var_t gridVar,var;
  int gridn,i,j;//<grid size and value indexex
  int vI;//variable index
  
  poc_data_t<double> h,dhdx,dhdy;
  zmatrix2x2_t *fric1=NULL,*fric2=NULL;
  double *hd;
  double *totalDissipation=NULL;

  struct timeval before;
  string cmd=fct_echo(aC,as);

  ws.n=aC-4-vC;
  if(ws.n<1){
    if(aC<2 ||
      strcmp(as[1],"--help")==0 ||
      strcmp(as[1],"-h")==0){
      print_help(as[0]);
      exit(0);
      }
    __ERR_BASE_LINE__("Missing %d arguments\n",1-ws.n);
    if(ws.n<0)
      fprintf(stderr,"*** Please provide a bathymetry path and variable name, an atlas file name convention, 3 variable names and a list of waves. ***\n");
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
      hP=a;
      break;
    case 2:
      hN=a;
      break;
    case 3:
      aPF=a;
      break;
    case 4:
    case 5:
    case 6:
      vNs[aI-4]=a;
      break;
    default:
      wNs[wI]=a;
      if(!strcmp(a,"M2"))M2wI=wI;
      wI++;
      }
    }
  wNs[wI]=NULL;

  ///It initialises the list of waves
  if(M2wI>=0){///making sure M2 is at the start
    char *w0=wNs[0];
    wNs[0]=wNs[M2wI];
    wNs[M2wI]=w0;
    }
  ws.init(initialize_tide(),wNs);
  printf("%d valid wave names given.\n",ws.n);
  
  gettimeofday(&before);
  
/*---------------------------------------------------------------------*//**<h1>
  reads the depth </h1>*/
  printf("Bathymetry from %s variable %s: allocating, ",hP,hN);fflush(stdout);
  status=h.init(hP,hN);
  if(status) __NC_TRAP_ERROR__(wexit,status,1,"poc_data_t::init(\"%s\",\"%s\") error",hP,hN);
  
  printf("grid, ");fflush(stdout);
  status=poc_get_grid(hP,h.info,&grid);
  if(status) __NC_TRAP_ERROR__(wexit,status,1,"poc_get_grid(,,\"%s\") error",hP);
  gridn=grid.nx*grid.ny;
  
  printf("data, ");fflush(stdout);
  status=h.read_data(hP,0,0,1);
  
  hd=h.data;
  printf("done in %gs.\n",difftime(&before));
  grid.print(cout,__FILE__ ":" TOSTRING(__LINE__) ":");
  
  h.info.type=NC_DOUBLE;
  h.info.attributes.erase("scale_factor").erase("add_offset").erase(_FillValue)
    .erase("valid_min").erase("valid_max");
  h.info.name="h";
  
/*---------------------------------------------------------------------*//**<h1>
  calculate depth gradient </h1>*/
  printf("bathymetry gradients: ");fflush(stdout);
  
  dhdx.info=h.info;
  dhdx.info.init("dhdx",NC_DOUBLE,"x_derivative_of_model_sea_floor_depth_below_geoid");
  dhdx.init();

  dhdy.info=h.info;
  dhdy.info.init("dhdy",NC_DOUBLE,"y_derivative_of_model_sea_floor_depth_below_geoid");
  dhdy.init();
  
  gradient(grid,hd,NC_FILL_DOUBLE,dhdx.data,dhdy.data);
  
  printf("done in %gs.\n",difftime(before));
  
  totalDissipation=new double[gridn];//this should set all values to 0.
  
/*---------------------------------------------------------------------*//**<h1>
  carry out the energy budget for each wave </h1>*/
  for(wI=0;wI<ws.n;wI++){//for each wave
    int verbose=0;
    const int nparts=2;
    const string pNs[nparts]={AMPSUF,PHASUF};//<part names
    const string waveName=ws.waves[wI].name,
      SNS="at_"+waveName+"_frequency_due_to_non_equilibrium_ocean_tide";
    const string oNs[vC]={"E","U","V"},//<output variable names
      vUs[vC]={"m","m s-1","m s-1"};//<variable units
    complex<double> *vars[vC],*&e=vars[0],*&u=vars[1],*&v=vars[2];
    complex<double> *dedx,*dedy;
    
    const double density=1025.,gravity=9.81,rho_g=density*gravity;
    double factor,total;
    double *fu,*fv;//<energy flux components
    double *dissipation,*a,*Ee,*Eu,*Q;//<dissipation rate
    
    double tt,ft,dft,wt,drt,pet,cet,qt;
    
    gettimeofday(&before);
    
    oP=waveName+"-energy.nc";//output path
    printf("preparing "+oP+": grid,");fflush(stdout);
    status=poc_save_grid(oP,&gridVar,grid,1,verbose);
    printf(" attributes,");fflush(stdout);
    status=poc_def_att(oP,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
    status=poc_def_att(oP,poc_att_t("history",cmd));
    printf(" bathymetry,");fflush(stdout);
    h.info.dimensions=gridVar.dimensions;
    dhdx.info.dimensions=gridVar.dimensions;
    dhdy.info.dimensions=gridVar.dimensions;
    
    status=h.write_data(oP);
    __NC_CHKERR_LINE_FILE__(status,"poc_data_t::write_data(\""+oP+"\",) error");
    status=dhdx.write_data(oP);
    __NC_CHKERR_LINE_FILE__(status,"poc_data_t::write_data(\""+oP+"\",) error");
    status=dhdy.write_data(oP);
    __NC_CHKERR_LINE_FILE__(status,"poc_data_t::write_data(\""+oP+"\",) error");
    
    printf(" done in %gs.\n",difftime(before));
    
/*---------------------------------------------------------------------*//**<h2>
    Read elevation and both speed components from the relevant atlases</h2>*/
    for(vI=0;vI<vC;vI++){//for each component
      grid_t vgrid;
      string ampName=vNs[vI]+pNs[0];
      
      gettimeofday(&before);
      
      ///decode the atlas name with decode_atlasname()
      aP=decode_atlasname(aPF,ws.waves[wI].name,vNs[vI]);
      printf("Wave %3s constants for %3s in %s: ",ws.waves[wI].name,vNs[vI],aP);fflush(stdout);
      
      ///get the component grid with poc_get_grid()
      printf("grid, ");fflush(stdout);
      status=poc_get_grid(aP,ampName,&vgrid,1);
      __NC_CHKERR_LINE_FILE__(status,"poc_get_grid(\"%s\",\""+ampName+"\",) error",aP);
      
      ///read the data with poc_get_vara()
      printf("read, ");fflush(stdout);
      complex<double> *buf=new complex<double>[vgrid.nx*vgrid.ny];
      status=poc_get_cvara(aP,ampName,vNs[vI]+pNs[1],-1,buf,1);
      
      ///interpolate with map_export
      printf("interpolation, ");fflush(stdout);
      vars[vI]=new complex<double>[gridn];
      map_export(vgrid,buf,NC_FILL_DOUBLE,grid,vars[vI],NC_FILL_DOUBLE);
      delete[]buf;
      
      ///write the data with poc_put_vara()
      printf("write, ");fflush(stdout);
      poc_var_t av,gv;
      av=gv=gridVar;
      av.init(oNs[vI]+pNs[0],NC_DOUBLE,comodo_standard_name("",waveName),vUs[vI],NAN,long_name_from_varname());
      gv.init(oNs[vI]+pNs[1],NC_DOUBLE,comodo_standard_name("",waveName),"degrees",NAN,long_name_from_varname());
      status=poc_put_cvara(oP,av,gv,0,vars[vI],1);
      
      delete[]aP;
      printf("done in %gs.\n",difftime(before));
      }
    
    /** We take :
    - the depth \f$ h \f$
    - the complex elevation \f$ \eta \f$
    - the complex speed components \f$ u \f$
    */
    
/*---------------------------------------------------------------------*//**<h2>
    Stokes'transport</h2>
    \f$ S \f$ with: \f[ \overrightarrow{S} = \eta \cdot \overrightarrow{u} \f] */
    printf("transport");fflush(stdout);
    fu=new double[gridn];
    fv=new double[gridn];
    gettimeofday(&before);
    for(i=0;i<gridn;i++){
      if(e[i]==NC_FILL_DOUBLE){
        fu[i]=NC_FILL_DOUBLE;
        fv[i]=NC_FILL_DOUBLE;
        continue;
        }
      if(u[i]!=NC_FILL_DOUBLE)
        fu[i]=e[i]%u[i];
      else
        fu[i]=NC_FILL_DOUBLE;
      if(v[i]!=NC_FILL_DOUBLE)
        fv[i]=e[i]%v[i];
      else
        fv[i]=NC_FILL_DOUBLE;
      }
    tt=difftime(before);
    printf(", ");fflush(stdout);
    var=gridVar;
    var.init("utransport",NC_DOUBLE,"sea_water_x_transport_"+SNS,"m2 s-1");
    poc_put_vara(oP,var,0,fu,1);
    var.init("vtransport",NC_DOUBLE,"sea_water_y_transport_"+SNS,"m2 s-1");
    poc_put_vara(oP,var,0,fv,1);
    
/*---------------------------------------------------------------------*//**<h2>
    energy flux</h2>
    \f$ F \f$ with: \f[ \overrightarrow{F} = 0.5 \rho g h \overrightarrow{S} \f] */
    printf("flux");fflush(stdout);
    gettimeofday(&before);
    for(i=0;i<gridn;i++){
      if(hd[i]==NC_FILL_DOUBLE){
        fu[i]=NC_FILL_DOUBLE;
        fv[i]=NC_FILL_DOUBLE;
        continue;
        }
      factor=rho_g*hd[i]*.5;
      if(fu[i]!=NC_FILL_DOUBLE)
        fu[i]*=factor;
      if(fv[i]!=NC_FILL_DOUBLE)
        fv[i]*=factor;
      }
    ft=difftime(before);
    printf(", ");fflush(stdout);
    var.init("uflux",NC_DOUBLE,"sea_water_x_energy_flux_"+SNS,"W m-1",NAN,"x_energy_flux");
    poc_put_vara(oP,var,0,fu,1);
    var.init("vflux",NC_DOUBLE,"sea_water_y_energy_flux_"+SNS,"W m-1",NAN,"y_energy_flux");
    poc_put_vara(oP,var,0,fv,1);
    
/*---------------------------------------------------------------------*//**<h2>
    divergence of flux</h2>
    \note Because it is obtained from a spatial derivation, is not calculated for values at the boundary of the domain */
    printf("divergence of flux");fflush(stdout);
    dissipation=new double[gridn];
    a=new double[gridn];
    aset(dissipation,gridn,NC_FILL_DOUBLE);
    gettimeofday(&before);
    #pragma omp parallel for private(i)
    for(j=grid.nx;j<gridn-grid.nx;j+=grid.nx){
      for(i=1;i<grid.nx-1;i++){
        vector3_t d;
        double dx,dy;
        d=math_polar2cartesian(grid.x[j+i-1]*d2r,grid.y[j+i-1]*d2r,MeanEarthRadius)-
          math_polar2cartesian(grid.x[j+i+1]*d2r,grid.y[j+i+1]*d2r,MeanEarthRadius);
        dx=hypot(d);
        d=math_polar2cartesian(grid.x[j+i-grid.nx]*d2r,grid.y[j+i-grid.nx]*d2r,MeanEarthRadius)-
          math_polar2cartesian(grid.x[j+i+grid.nx]*d2r,grid.y[j+i+grid.nx]*d2r,MeanEarthRadius);
        dy=hypot(d);
        a[j+i]=.25*dx*dy;
        if(fu[j+i-1]==NC_FILL_DOUBLE
         ||fu[j+i+1]==NC_FILL_DOUBLE
         ||fv[j+i-grid.nx]==NC_FILL_DOUBLE
         ||fv[j+i+grid.nx]==NC_FILL_DOUBLE
         ){
          totalDissipation[j+i]=NC_FILL_DOUBLE;
          continue;
          }
        dissipation[j+i]=0.;
        dissipation[j+i]-=(fu[j+i-1]-fu[j+i+1])/dx;
        dissipation[j+i]-=(fv[j+i-grid.nx]-fv[j+i+grid.nx])/dy;
        totalDissipation[j+i]+=dissipation[j+i];
        }
      }
    dft=difftime(before);
    printf(", ");fflush(stdout);
    var.init("divFlux",NC_DOUBLE,"horizontal_divergence_of_sea_water_energy_flux_"+SNS,"W m-2",NAN,"divergence_of_energy_flux");
    poc_put_vara(oP,var,0,dissipation,1);
    
/*---------------------------------------------------------------------*//**<h2>
    pressure work</h2>
    \f$ W_p \f$ with: \f[ W_p = 0.5 \rho g h (\nabla \eta) \cdot u \f] */
    printf("pressure work");fflush(stdout);
    dedx=new complex<double>[gridn];
    dedy=new complex<double>[gridn];
    gradient(grid,e,NC_FILL_COMPLEX,dedx,dedy);
    total=0.;
    gettimeofday(&before);
    for(i=0;i<gridn;i++){
      if(h.data[i]==NC_FILL_DOUBLE ||
         dedx[i]==NC_FILL_DOUBLE ||
         dedy[i]==NC_FILL_DOUBLE ||
         u[i]==NC_FILL_DOUBLE ||
         v[i]==NC_FILL_DOUBLE){
        dissipation[i]=NC_FILL_DOUBLE;
        continue;
        }
      factor=rho_g*hd[i]*.5;
      dissipation[i]=factor*(dedx[i]%u[i]+dedy[i]%v[i]);
      total+=dissipation[i]*a[i];
      }
    wt=difftime(before);
    printf(": %gW, ",total);fflush(stdout);
    var.init("pressureWork",NC_DOUBLE,"sea_water_pressure_work_"+SNS,"W m-2");
    poc_put_vara(oP,var,0,dissipation,1);
    
/*---------------------------------------------------------------------*//**<h2>
    dissipation</h2>
    \f$ W_d \f$ with: \f[ W_d = C_d \rho |u| u \cdot u \f]*/
    printf("dissipation rate");fflush(stdout);
    ///taking \f$ C_d = \simeq 2.5e-3 \f$
    factor=-density*2.5e-3*.5;
    if(!wI){
      fric1=new zmatrix2x2_t[gridn];
      fric2=new zmatrix2x2_t[gridn];
      }
    gettimeofday(&before);
    for(i=0;i<gridn;i++){
      dcomplex uu,vv;//force components
      if(u[i]==NC_FILL_DOUBLE ||
         v[i]==NC_FILL_DOUBLE){
        dissipation[i]=NC_FILL_DOUBLE;
        continue;
        }
      ///and \f$ |u| \f$ from:
      zmatrix2x2_t *fric;
      if(!wI){
        ///- spectral_friction01() for the dominant wave
        fric=&fric1[i];
        spectral_friction01( fric    ,u[i],v[i]);
        spectral_friction02(&fric2[i],u[i],v[i]);
        }
      else{
        ///- spectral_friction02() for all others
        fric=&fric2[i];
        }
      uu=fric->c[0][0]*u[i]+fric->c[1][0]*v[i];
      vv=fric->c[0][1]*u[i]+fric->c[1][1]*v[i];
      dissipation[i]=factor*(uu%u[i]+vv%v[i]);
      }
    drt=difftime(before);
    printf(", ");fflush(stdout);
    var.init("dissipation",NC_DOUBLE,"sea_water_dissipation_"+SNS,"W m-2");
    poc_put_vara(oP,var,0,dissipation,1);
    
/*---------------------------------------------------------------------*//**<h2>
    potential energy</h2>
    \f$ E_\eta \f$ with: \f[ E_\eta = .5 \rho g \eta^2 \f] */
    printf("potential energy");fflush(stdout);
    Ee=new double[gridn];
    factor=.5*rho_g;
    gettimeofday(&before);
    for(i=0;i<gridn;i++){
      if(e[i]==NC_FILL_DOUBLE){
        Ee[i]=NC_FILL_DOUBLE;
        continue;
        }
      Ee[i]=factor*(e[i]%e[i]);
      }
    pet=difftime(before);
    printf(", ");fflush(stdout);
    var.init("potE",NC_DOUBLE,"gravitational_potential_energy_"+SNS,"J m-2",NAN,"gravitational_potential_energy");
    poc_put_vara(oP,var,0,Ee,1);
    
/*---------------------------------------------------------------------*//**<h2>
    kinetic energy</h2>
    \f$ E_u \f$ with: \f[ E_u = .5 \rho h (|u|^2 + |v|^2) \f] */
    printf("kinetic energy");fflush(stdout);
    Eu=new double[gridn];
    factor=.5*density;
    gettimeofday(&before);
    for(i=0;i<gridn;i++){
      if(u[i]==NC_FILL_DOUBLE ||
         v[i]==NC_FILL_DOUBLE ||
         hd[i]==NC_FILL_DOUBLE){
        Eu[i]=NC_FILL_DOUBLE;
        continue;
        }
      Eu[i]=factor*hd[i]*(u[i]%u[i]+v[i]%v[i]);
      }
    cet=difftime(before);
    printf(", ");fflush(stdout);
    var.init("cinE",NC_DOUBLE,"kinetic_energy_content_"+SNS,"J m-2",NAN,"kinetic_energy");
    poc_put_vara(oP,var,0,Eu,1);
    
/*---------------------------------------------------------------------*//**<h2>
    Q factor</h2>
    \f$ Q \f$ with: \f[ Q = (E_\eta + E_u) / (W T) \f] */
    printf("Q factor");fflush(stdout);
    Q=new double[gridn];
    factor=ws.waves[wI].omega*dph2rps/M_PIx2;//frequency
    total=0.;
    gettimeofday(&before);
    for(i=0;i<gridn;i++){
      if(Ee[i]==NC_FILL_DOUBLE ||
         Eu[i]==NC_FILL_DOUBLE ||
         dissipation[i]==NC_FILL_DOUBLE ||
         /* as per IFREMER's request, because it messes-up ferret */
         dissipation[i]==0.){
        Q[i]=NC_FILL_DOUBLE;
        continue;
        }
      Q[i]=(Ee[i]+Eu[i])*factor/dissipation[i];
      total+=(Ee[i]+Eu[i])*a[i];
      }
    qt=difftime(before);
    printf(": kinetic+potential=%gJ, ",total);fflush(stdout);
    var.init("Q",NC_DOUBLE,"quality_factor_of_non_equilibrium_ocean_tide","",NAN,"quality_factor");
    poc_put_vara(oP,var,0,Q,1);
    
    printf("done.\n",tt,ft,dft,wt,drt,pet,cet,qt);
    delete[]fu;
    delete[]fv;
    delete[]dissipation;
    delete[]dedx;
    delete[]dedy;
    delete[]Ee;
    delete[]Eu;
    delete[]Q;
    }//EO for each wave
  
  oP="total-energy.nc";
  poc_save_grid(oP,grid,__FILE__,__LINE__);
  var.init("divFlux",NC_DOUBLE,"horizontal_divergence_of_sea_water_energy_flux_due_to_non_equilibrium_ocean_tide","W m-2",NAN,"divergence_of_total_tidal_energy_flux");
  poc_put_vara(oP,var,0,totalDissipation,1);
  return 0;
}
