
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Expands the spectrum of a set of atlases using the admittance method

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <sys/stat.h> /* for fstat() */
#include <unistd.h> // for access

#include "tools-structures.h"

#include "functions.h" //fct_echo
#include "poc-netcdf-data.hpp"
#include "fe-proto.h"
#include "map.h"
#include "tides.h"
#include "admittance-mts.h"
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
    "  %s [OPTIONS] wave1 wave2 wave3 [wave4 ...]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Expands the spectrum of a set of atlases using the admittance method.\n"
    "  The list of waves (at LEAST 3) are the waves you have and the waves you want.\n"
    "  Missing atlases will be created with their constants interpolated or extrapolated with the admittance method.\n"
    "  If only waves are given, prints the coefficients.\n"
    "\n"
    "TIP\n"
    "  You will have better atlases if you use comodo-detidor with a complete list of waves instead.\n"
    "\n"
    "IMPORTANT WARNINGS\n"
    "  When adding frequencies, we have among other : K1+O1=M2 , K1+P1=S2 and M1+O1=N2. See the output of `showarg --combinations' for more.\n"
    "SO THE MORE NON-LINEAR AND MIXED SEMIDIURNAL THE ZONE IS, THE LESS RELIABLE THE ADMITTANCE METHOD IS. USE AT YOUR OWN RISK.\n"
    "The input waves must be strongest in their category: only M2, K2, N2, K1, O1, Q1, Mf, Mm and Mtm should be allowed as input!\n"
    "  S2 has got a strong radiative component, so\n"
    "DO NOT PUT S2 IN THE LIST OF WAVES UNLESS YOU ARE AWARE OF WHAT RISKS YOU ARE TAKING.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help : Show this help and exit.\n"
    "  -a : followed by atlas file name convention. See below.\n"
    "  -v : followed by 2 variable names, repectively for amplitude and phase\n"
    "  -d : followed by the discretisation, e.g. LGP2. If unspecified, a structured grid will be assumed\n"
    "  -p : prefix for forcing admittance of last wave. BUG: NOT AVAILABLE FOR UNSTRUCTURED GRIDS YET!\n");
  print_tide_decode_atlasname_help();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LGP1_area(mesh_t *mesh,double **area, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,m,n;
  
  fe_vertex_element_tables(mesh);/* for vertex->elmts */
  fe_initaffine(mesh);/* for triangles[m].TrueArea */
  
  *area=new double[mesh->nvtxs];
  
  for(n=0;n<mesh->nvtxs;n++){
    const vertex_t *vertex=&mesh->vertices[n];
    
    if(vertex->nelmts<=0 and verbose>0)
      STDOUT_BASE_LINE("surprising : vertices[%d].nelmts=%d<=0\n",n,vertex->nelmts);
    
    double vArea=0.;
    
    for(i=0;i<vertex->nelmts;i++){
      m=vertex->elmts[i];
      vArea+=mesh->triangles[m].TrueArea;
      }
    if(not isfinite(vArea)) TRAP_ERR_EXIT(-1,"vertices[%d] has area=%g\n",n,vArea);
    vArea*=1./3;
    
    (*area)[n]=vArea;
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int remove_average_template(const double *area,int nvtxs,complex<T> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  
  double totalArea=0.;
  complex<double> average,zsum=0.;
  
  for(n=0;n<nvtxs;n++){
    double arean=area[n];
    complex<T> *zn=&z[n];
    
    totalArea+=arean;
    zsum+=((complex<double>)*zn)*arean;
    }
  
  average=zsum/totalArea;
  
  if(not isfinite(average)) TRAP_ERR_EXIT(-1,"average=%g/%g=%g\n",zsum/totalArea,average);
  
  for(n=0;n<nvtxs;n++){
    z[n]-=average;
    }
}


/*------------------------------------------------------------------------------
  option variables*/

string cmd;
spectrum_t aws;     //list of astronomical waves within those given by the user
char *aPF=NULL;     //atlas path convention
const int vC=2;     //number of variables
char *vNs[vC];      //names of variables
char *discN=NULL;
string prefix;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int admittance_UG()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Expands the spectrum of a set of unstructured atlases using the admittance method
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int *wA,missingC,availableC;//wave availability, number of missing, number of available
  int wI,status;//wave index, status
  char *aP=0;//atlas path
  int vI;//variable index
  
  int discretisation=discretisation_from_name(discN);
  discretisation_t descriptor;
  
  mesh_t mesh;
  int nC,nI;//number of nodes and index
  double *buf[vC];
  hconstant_t *constants=NULL;
  
  poc_var_t var;
  string units;
  
/*---------------------------------------------------------------------*//**<h1>
  checks available and missing atlases </h1>*/
  wA=new int[aws.n];
  missingC=0;
  availableC=0;
  bool doEquilibrium=false;
  
  for(wI=0;wI<aws.n;wI++){
    int *wAI=&wA[wI];
    tidal_wave *waveI=&aws.waves[wI];
    
    tide_decode_atlasname(NULL,aPF,waveI->name, 0,&aP);
    status=access(aP,R_OK);//NOTE: access rights may change during the run
    *wAI = (status==0);
    
    printf("Constants for %10s (%8.5g deg/h) : %s ",waveI->name,waveI->omega,aP);
    int thisDisc;
    
    if(*wAI){
      printf("exists : reading... ");fflush(stdout);
      ///reading available atlases
      
      if(constants==NULL){
        printf("(mesh:");fflush(stdout);
        
        status=poc_inq_var(aP,vNs[0],&var);
        
        poc_att_t *ua;
        ua=var.attributes.findP("units");
        if(ua!=0)
          units=ua->as_string();
        
        status=poc_get_mesh(aP,var,&mesh,0,-1,&thisDisc);
        NC_TRAP_ERROR(return,status,1,"poc_get_mesh error with %s",aP);
        
        if(thisDisc!=discretisation and discretisation==LGP1)
          doEquilibrium=true;
        
        discretisation_init(&mesh,discretisation,0);
        descriptor=get_descriptor(mesh,discretisation);
        nC=descriptor.nnodes;
        printf("%d %s points",nC,discN);fflush(stdout);
        if(nC==0){
          printf("\n");fflush(stdout);
          TRAP_ERR_EXIT(ENOEXEC,"Not coded for 0 point fields\n");
          }
        
        printf(",allocating buffers");fflush(stdout);
        constants=new hconstant_t[nC];
        for(nI=0;nI<nC;nI++){
          constants[nI].init_polar(aws.n);
          }
        for(vI=0;vI<vC;vI++){
          buf[vI]=new double[nC];
          }
        printf(") ");fflush(stdout);
        }
      
      if(thisDisc==discretisation){
        for(vI=0;vI<vC;vI++){
          status=poc_get_UG3D(aP, -1, vNs[vI], buf[vI]);
          if(status==NC_ENOTVAR){
            printf("%s not found.\n",vNs[vI]);
            *wAI=0;
            goto missingVar;
            }
          if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_get_UG3D error with %s",vNs[vI]);
          }
        
        for(nI=0;nI<nC;nI++){
          constants[nI].a[wI]=buf[0][nI];
          constants[nI].G[wI]=buf[1][nI];
          }
        }
      
      printf("done.\n");
      availableC++;
      }
    else{
      printf("does not exist.\n");
missingVar:
      missingC++;
      }
    
    deletep(&aP);
    }
  
  if(availableC==0)
    TRAP_ERR_EXIT(-1,"*** No constant available. ***\n");
  
  if(missingC==0)
    TRAP_ERR_EXIT(0,"No wave missing.\n");/* behave like /bin/true */
  
/*------------------------------------------------------------------------------
  calculates missing constants */
  admittance_t admittance;
  
  delete[] admittance_mts_check(&admittance,aws,wA);
  
  status=admittance_mts_verify(admittance,aws,wA,1);
  
  if(status!=0){
    if(missingC==aws.n-1 and discretisation==LGP1){
      doEquilibrium=true;
      STDERR_BASE_LINE("LGP1 discretisation and only %d wave available and %d of the %d other(s) missing : EQUILIBRIUM\n",
        availableC,missingC,aws.n-availableC,status);
      }
    else
      TRAP_ERR_EXIT(-1,"admittance_mts_verify() error : %d waves can not be resolved\n",status);
    }
  
  double *area=0;
  
  if(doEquilibrium){
    
    for(wI=0;wI<aws.n;wI++){
      if(wA[wI]!=0)
        continue;
      tidal_wave *waveI=&aws.waves[wI];
      
      for(nI=0;nI<nC;nI++){
        hconstant_t *ci=&constants[nI];
        const vertex_t *vertex=&mesh.vertices[nI];
        status=tidal_equilibrium(*waveI,vertex->lon,vertex->lat,&ci->a[wI],&ci->G[wI]);
        if(status!=0){
          wA[wI]=status;
          break;
          }
        }
      }
    
    status=LGP1_area(&mesh,&area,1);
    }
  else
   for(nI=0;nI<nC;nI++){
    admittance_mts_compute(admittance, aws, constants[nI], wA, NC_FILL_FLOAT, false, 0);
    }

/*------------------------------------------------------------------------------
  creates missing atlases */
  complex<float> *z=new complex<float>[nC];
  
  for(wI=0;wI<aws.n;wI++){
    if(wA[wI])continue;
    for(nI=0;nI<nC;nI++){
      hconstant_t *ci=&constants[nI];
      z[nI]=polar<float>(ci->a[wI],-(ci->G[wI])*d2r);
      }
    if(doEquilibrium)
      status=remove_average_template(area,mesh.nvtxs,z);
    tide_decode_atlasname(NULL,aPF,aws.waves[wI].name, 0,&aP);
    printf("writing to %s... ",aP);fflush(stdout);
    status=poc_put_mesh_cvara(aP,mesh,discretisation,vNs[0],vNs[1],-1,z,units,0);
    printf("done.\n");
    delete[]aP;
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int admittance_RG()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Expands the spectrum of a set of structured atlases using the admittance method
/*----------------------------------------------------------------------------*/
{
  int *wA,missingC,availableC;//wave availability, number of missing, number of available
  int wI,status;//wave index, status
  char *aP;//atlas path
  int vI;//variable index
  int verbose=0;
  
  poc_grid_data_t gd;
  int nC,nI;//number of nodes and index
  poc_data_t<float> data[vC];
  hconstant_t *constants=NULL;
  
/*---------------------------------------------------------------------*//**<h1>
  checks available and missing atlases </h1>*/
  wA=new int[aws.n];
  missingC=availableC=0;
  
  for(wI=0;wI<aws.n;wI++){
    poc_global_t global;
    tide_decode_atlasname(NULL,aPF,aws.waves[wI].name, 0,&aP);
    wA[wI]=!access(aP,R_OK);//NOTE: access rights may change during the run
    printf("Constants for %10s (%8.5g deg/h) : %s ",aws.waves[wI].name,aws.waves[wI].omega,aP);
    if(wA[wI]){
      printf("exists : reading... ");fflush(stdout);
      ///reading available atlases
      if(constants==NULL){
        printf("(allocating");fflush(stdout);
        poc_inq(aP,&global,0);
        for(vI=0;vI<vC;vI++){
          status=data[vI].init(global,vNs[vI]);
          if(status<0){
            NC_CHKERR_LINE_FILE(status,"poc_data_t::init(\"%s\",\"%s\") error",aP,vNs[vI]);
            wexit(-2);
            }
          }
        printf(" %dD buffers",data[0].info.dimensions.size());fflush(stdout);
        nC=data[0].length;
        constants=new hconstant_t[nC];
        for(nI=0;nI<nC;nI++){
          constants[nI].init_polar(aws.n);
          }
        printf(", grid data");fflush(stdout);
        status=poc_get_grid_data(aP,data[0].info,&gd);
        if(status!=NC_NOERR){
          NC_CHKERR_LINE_FILE(status,"poc_getgrid_data(\""+data[0].info.name+"\",,\""+aP+"\") error");
          wexit(-1);
          }
        printf(") ");fflush(stdout);
        }
      for(vI=0;vI<vC;vI++){
        data[vI].read_data(aP,0,0,1);
        }
      for(nI=0;nI<nC;nI++){
        constants[nI].a[wI]=data[0].data[nI];
        constants[nI].G[wI]=data[1].data[nI];
        }
      printf("done.\n");
      availableC++;
      }
    else{
      printf("does not exist.\n");
      missingC++;
      }
    delete[]aP;
    }
  
  if(availableC==0)
    TRAP_ERR_EXIT(-1,"*** No constant available. ***\n");
  if(missingC==0 and prefix=="")
    TRAP_ERR_EXIT(0,"No wave missing.\n");//behave like /bin/true

/*---------------------------------------------------------------------*//**<h1>
  calculates missing constants with admittance_mts_compute() </h1>*/
  admittance_t admittance;
  delete[] admittance_mts_check(&admittance,aws,wA);
  
  status=admittance_mts_verify(admittance,aws,wA,1);
  if(status!=0) TRAP_ERR_EXIT(-1,"admittance_mts_verify() error : %d waves can not be resolved\n",status);
  
  struct timeval before;
  gettimeofday(&before);
  
  int done=0;
  #pragma omp parallel for private(nI) if(nC>1000000)
  for(nI=0;nI<nC;nI++){
    if(data[0].data[nI]==data[0].mask)
      continue;
    if(verbose>0 and not done){
      for(wI=0;wI<aws.n;wI++)
        STDOUT_BASE_LINE("%d %g %g %g\n",wA[wI],aws.waves[wI].omega,constants[nI].a[wI]/aws.waves[wI].Ap,constants[nI].G[wI]);
      }
//     bool force=prefix!="";
    admittance_mts_compute(admittance,aws,constants[nI],wA, NC_FILL_FLOAT, false, 0);
    if(verbose>0 and not done){
      for(wI=0;wI<aws.n;wI++)
        STDOUT_BASE_LINE("%d %g %g %g\n",wA[wI],aws.waves[wI].omega,constants[nI].a[wI]/aws.waves[wI].Ap,constants[nI].G[wI]);
      done=1;
      }
    }
  STDERR_BASE_LINE("%d points took %gs.\n",nC,difftime(before));

/*---------------------------------------------------------------------*//**<h1>
  creates missing atlases. </h1>*/
  for(wI=0;wI<aws.n;wI++){
    if(wA[wI] and not(prefix!="" and wI==0))continue;
    poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
    global << poc_att_t("history",cmd);
    for(nI=0;nI<nC;nI++){
      if(data[0].data[nI]==data[0].mask){
        data[1].data[nI]=data[1].mask;
        continue;
        }
      data[0].data[nI]=constants[nI].a[wI];
      data[1].data[nI]=constants[nI].G[wI];
      }
    tide_decode_atlasname(NULL,aPF,aws.waves[wI].name, 0,&aP);
    if(wI==0){
      char *junk=aP;
      aP=strdup((prefix+aP).c_str());
      delete[]junk;
      }
    printf("writing to %s (",aP);fflush(stdout);
    status=poc_create(aP,global,1);
    if(status!=NC_NOERR){NC_CHKERR_LINE_FILE(status,"poc_create(\"%s\",,) error",aP);goto EOwFor;}
    for(vI=0;vI<gd.vc;vI++){
      if(gd[vI].data==NULL)continue;
      printf(gd[vI].info.name+" ");fflush(stdout);
      status=gd[vI].write_data(aP);
      if(status!=NC_NOERR){NC_CHKERR_LINE_FILE(status,"poc_data_t<>::write_data(\"%s\") error with "+gd[vI].info.name,aP);goto EOwFor;}
      }
    for(vI=0;vI<vC;vI++){
      printf(data[vI].info.name+" ");fflush(stdout);
      data[vI].info.type=NC_FLOAT;
      data[vI].info << poc_att_t(_FillValue,(float)data[vI].mask);//for ncview
      status=data[vI].write_data(aP,0,0,1);
      if(status!=NC_NOERR){NC_CHKERR_LINE_FILE(status,"poc_data_t<>::write_data(\"%s\") error with ",aP,vNs[vI]);goto EOwFor;}
      }
    printf(") done.\n");
EOwFor:
    delete[]aP;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void admittance_show()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int *wA=0,wI,nC=0,nI;
  
  wA=aset(aws.n,0);
  
  for(wI=0;wI<aws.n;wI++){
    tidal_wave *waveI=&aws.waves[wI];
    if( strcmp(waveI->name,"M2")!=0 and
        strcmp(waveI->name,"N2")!=0 and
        strcmp(waveI->name,"K2")!=0 and
        strcmp(waveI->name,"K1")!=0 and
        strcmp(waveI->name,"O1")!=0 and
        strcmp(waveI->name,"Q1")!=0 and
        strcmp(waveI->name,"Mf")!=0 and
        strcmp(waveI->name,"Mm")!=0 and
        strcmp(waveI->name,"Mtm")!=0 )
      continue;
    wA[wI]=1;
    nC++;
    }
  
  admittance_t admittance;
  
  delete[] admittance_mts_check(&admittance,aws,wA);
  
  status=admittance_mts_verify(admittance,aws,wA,1);
  
  string *names;
  hconstant_t *constants;
  
  names=new string[nC];
  constants=new hconstant_t[nC];
  
  nI=-1;
  for(wI=0;wI<aws.n;wI++){
    if(wA[wI]==0)
      continue;
    nI++;
    
    names[nI]=aws.waves[wI].name;
    
    hconstant_t *constantI=&constants[nI];
    
    constantI->init_polar(aws.n);
    aset(constantI->a,aws.n,0.f);
    aset(constantI->G,aws.n,0.f);
    constantI->a[wI]=1;
    
    admittance_mts_compute(admittance, aws, *constantI, wA, NC_FILL_FLOAT, false, 0);
    }
  
  for(wI=0;wI<aws.n;wI++){
    
    if(wA[wI]==1)
      continue;
    
    printf("%s=",aws.waves[wI].name);
    for(nI=0;nI<nC;nI++){
      float aII=real(polar(constants[nI].a[wI],(float)(-constants[nI].G[wI]*d2r)));
      if(aII==0.f)
        continue;
      printf(" %+g "+names[nI],aII);
      }
    printf("\n");
    
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
  int aI,wI=0,status;//argument index, wave index, status

  char *a;//argument
  char **wNs;//list of wave names
  spectrum_t ws;//list of waves given by the user

  int vI;//variable index
  
/*---------------------------------------------------------------------*//**<h1>
  parses the arguments </h1>*/
  cmd=fct_echo( aC,as);
  
  aI=1;
  wNs=new char*[aC];
  for(vI=0;vI<vC;vI++){
    vNs[vI]=NULL;
    }
  while (aI < aC) {
    a=strdup(as[aI]);
    switch (a[0]) {
      case '-':
        if(strcmp(a,"-h")==0 || strcmp(a,"--help")==0 ) {
          print_help(as[0]);
          return 0;
          }
        switch (a[1]) {

        case 'a' :
          aPF= strdup(as[aI+1]);
          aI++;
          aI++;
          break;

        case 'v' :
          aI++;
          for(vI=0;vI<vC && aI<aC;vI++){
            vNs[vI]=strdup(as[aI]);
            aI++;
            }
          break;

        case 'd' :
          discN=strdup(as[aI+1]);
          aI++;
          aI++;
          break;

        case 'p' :
          prefix=as[aI+1];
          aI++;
          aI++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",a);
          wexit(-1);
        }
        break;

      default:
        wNs[wI]=strdup(a);
        wI++;
        aI++;
        break;
      }
      free(a);
    }
  wNs[wI]=NULL;
  
/*------------------------------------------------------------------------------
  check everything compulsory is provided for */
  status=0;
  
  for(vI=0;vI<vC;vI++){
    if(vNs[vI]==NULL)
      status=-1;
    }
  
  ws.n=wI;
  
  if( ws.n<3 or
      ( (status!=0) != (aPF==NULL) )
      ){
    
    if(status!=0){
      fprintf(stderr,"*** Please provide 2 variable names. ***\n");
      }
    
    if(aPF==NULL){
      fprintf(stderr,"*** Please provide an atlas file name convention. ***\n");
      }
    
    if(ws.n<3){
      fprintf(stderr,"*** Please provide a list of at least 3 waves. ***\n");
      }
    
    print_help(as[0]);
    wexit(-1);
    }

/*------------------------------------------------------------------------------
  print warning */
  #define STAR_LINE "\n**************************************************************\n"
  #define ADMITTANCE_METHOD_WARNING STAR_LINE "MAKE SURE YOU HAVE READ THE HELP OF THIS PROGRAMME\n AND ARE AWARE OF THE RISKS OF USING THE ADMITTANCE METHOD." STAR_LINE
  ///-on stderr
  STDERR_BASE_LINE(ADMITTANCE_METHOD_WARNING);
  struct stat outstat,errstat;
  fstat(fileno(stderr),&errstat);
  fstat(fileno(stdout),&outstat);
  ///-and on stdout if there is a redirection
  if(!S_ISCHR(errstat.st_mode) || !S_ISCHR(outstat.st_mode) || errstat.st_rdev!=outstat.st_rdev)
    STDOUT_BASE_LINE(ADMITTANCE_METHOD_WARNING);
  #undef ADMITTANCE_METHOD_WARNING
  #undef STAR_LINE

/*---------------------------------------------------------------------*//**<h1>
  initialises the list of waves </h1>*/
  ws.init(initialize_tide(),wNs);
  printf("# %d valid wave names given.\n",ws.n);
  ///stripping the non-astronomical waves
  for(wI=0;wI<ws.n;wI++)
    if(ws.waves[wI].Ap!=0.)
      aws.n++;
  
  aws.init(aws.n);
  for(wI=0;wI<ws.n;wI++)
    if(ws.waves[wI].Ap!=0.)
      aws.add(ws.waves[wI]);
    else{
      printf("%s is not astronomical : removing.\n",ws.waves[wI].name);
      }
  printf("# %d astronomical waves given.\n",aws.n);
  
  if(status!=0 and aPF==NULL){
    
    admittance_show();
    
    TRAP_ERR_EXIT(0,"exiting\n");
    }
  
  if(discN==NULL){
    status=admittance_RG();
    }
  else{
    status=admittance_UG();
    }
  
  return status;
}
