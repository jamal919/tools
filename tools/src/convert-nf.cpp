
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Converts TUGOm outputs. USED BY CTOH.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <unistd.h> //for unlink
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "tides.def"
#include "tides.h"
#include "fe.h"
#include "map.h"
#include "poc-time.h"
#include "netcdf-proto.h"
#include "poc-netcdf-data.hpp"
#include "filter.h"
#include "statistic.h"
#include "ascii.h"
#include "archive.def"
#include "archive.h"
#include "functions.h"
#include "datastream.h"


#define ARCHIVE_FORMAT_ASCII  0
#define ARCHIVE_FORMAT_NETCDF 1
#define ARCHIVE_FORMAT_BINARY 2

#define GLOBAL      0
#define PER_YEAR    1
#define PER_MONTH   2
#define PER_WEEK    3
#define PER_DAY     4
#define PER_HOUR    5
#define PER_MINUTE  6
#define PER_SECOND  7
#define CUSTOM      8


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_name(date_t actual, const char *name_template, const char *varname, char **filename, int *packing)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  status;
  char dummy[256], *pointer;
  date_t cdf_reference;
  FILE *out;

/*------------------------------------------------------------------------------
  Reconstruct the input file name*/
  (*filename) = new char[strlen(name_template) + 256];
  sprintf((*filename), "%s", name_template);

  out = fopen(*filename, "r");
  status = (out == NULL);

  switch (status) {
    case 0:
/*------------------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      *packing=GLOBAL;
      break;

    default:
/*------------------------------------------------------------------------------
      use format convention information*/
/*------------------------------------------------------------------------------
      substitute YYYY with current year*/
      pointer = strstr((*filename), "YYYY");
      if(pointer != NULL) {
        sprintf(dummy, "%4d", actual.year);
        strncpy(pointer, dummy, 4);
        *packing=PER_YEAR;
        }
/*------------------------------------------------------------------------------
      substitute MM with current month*/
      pointer = strstr((*filename), "MM");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.month);
        strncpy(pointer, dummy, 2);
        *packing=PER_MONTH;
        }
/*------------------------------------------------------------------------------
      substitute DD with current day*/
      pointer = strstr((*filename), "DD");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.day);
        strncpy(pointer, dummy, 2);
        *packing=PER_DAY;
        }
/*------------------------------------------------------------------------------
      substitute DD with current day*/
      pointer = strstr((*filename), "NNN");
      if(pointer != NULL) {
        int day=day_in_year(actual);
        sprintf(dummy, "%3.3d", day);
        strncpy(pointer, dummy, 3);
        *packing=PER_DAY;
        }
/*------------------------------------------------------------------------------
      substitute HH with current hour*/
      pointer = strstr((*filename), "HH");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", (int) floor(actual.second/3600.));
        strncpy(pointer, dummy, 2);
        *packing=PER_HOUR;
        }
/*------------------------------------------------------------------------------
      substitute MN with current hour*/
      pointer = strstr((*filename), "MN");
      if(pointer != NULL) {
        double hour=floor(actual.second/3600.);
        sprintf(dummy, "%2.2d", (int) floor((actual.second-3600.*hour)/60.));
        strncpy(pointer, dummy, 2);
        *packing=PER_HOUR;
        }
// /*------------------------------------------------------------------------------
//       lookup file fitting ???? with minute and seconds index*/
//       pointer = strstr((*filename), "????");
//       if(pointer != NULL) {
//         struct stat buf;
//         int i,j;
//         for(j=0;j<60;j++) {
//           for(i=0;i<60;i++) {
//             sprintf(dummy, "%2.2d%2.2d", j,i);
//             strncpy(pointer, dummy, 4);
//             status= get_file_size(*filename, 0);
// /*------------------------------------------------------------------------------
//             status=0 if file found*/
//             if(status==0) break;
//             }
//           if(status==0) break;
//           }
//         }
/*------------------------------------------------------------------------------
      substitute VARNAME*/
      pointer = strstr((*filename), "VARNAME");
      if(pointer != NULL) {
        sprintf(dummy, "%s", varname);
        strncpy(pointer, dummy, strlen(varname));
        }
      break;

      break;
    }

  return (status);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int filemapr1 (char* input,char* output,grid_t grid,int minstep,int maxstep,meta_archive_t info)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  {
  FILE *out;
  double time;
  int status,t;
  int dmax=1,kmax=1, shift;
  size_t offset;
  bmgheader_t header;
  float *levels,mask=999.9,*buf;
  int *element;
  date_t date=info.reference;

  status=map_completegridaxis_2(&grid);

  element=fe_scan_elements(info.mesh,grid,0);
  if (element == NULL) return(-1);

  buf=new float[grid.Hsize()];
  if (buf == NULL) return(-1);

  out=fopen(output,"wb");

  offset=4;

  strcpy(header.comment[0],"mog2d simulation");
  strcpy(header.comment[1],"no comment");
  strcpy(header.comment[2],input);
  strcpy(header.comment[3],"no comment");
  header.ni=grid.nx;
  header.nj=grid.ny;
  header.xmin=grid.xmin;
  header.ymin=grid.ymin;
  header.dx=grid.dx;
  header.dy=grid.dy;
  header.nk=kmax;
  header.nd=dmax;
  header.nt=maxstep-minstep+1;
  header.spec=mask;
  levels=(float *) malloc(kmax*sizeof(float));
  header.levels=levels;

  status=ascii_writeheader (out, header);
  if(status!=0) goto error;

  shift=mjd(date.month,date.day,date.year)-mjd(1,1,1950);

  for(t=minstep;t<=maxstep;t++) {
    status=fe_framemapr1(input, info, grid, mask, element, t, buf, &time);
    time=time/d2s+shift;
    if(status!=0) goto error;
    status=ascii_save_r (out, t, grid,  buf, mask, time);
    if(status!=0) goto error;
    }

  fclose(out);
  free(element);
  free(buf);
  return(0);

error:
  fclose(out);
  free(element);
  free(buf);
  return(-1);
  }

/*----------------------------------------------------------------------------*/

int filemapr2 (char* input,char* output,grid_t grid,int minstep,int maxstep,meta_archive_t info)

/*----------------------------------------------------------------------------*/
  {
  FILE *out;
  double time;
  int status,t;
  int dmax=2,kmax=1, shift;
  size_t offset;
  bmgheader_t header;
  float *levels,mask=999.9,*bufx,*bufy;
  int *element;
  date_t date=info.reference;

  status=map_completegridaxis_2(&grid);

  element=fe_scan_elements(info.mesh,grid,0);
  if (element == NULL) return(-1);

  bufx=new float[grid.Hsize()];
  if (bufx == NULL) return(-1);

  bufy=new float[grid.Hsize()];
  if (bufy == NULL) return(-1);

  out=fopen(output,"wb");
  
  offset=4;

  strcpy(header.comment[0],"mog2d simulation");
  strcpy(header.comment[1],"no comment");
  strcpy(header.comment[2],input);
  strcpy(header.comment[3],"no comment");
  header.ni=grid.nx;
  header.nj=grid.ny;
  header.xmin=grid.xmin;
  header.ymin=grid.ymin;
  header.dx=grid.dx;
  header.dy=grid.dy;
  header.nk=kmax;
  header.nd=dmax;
  header.nt=maxstep-minstep+1;
  header.spec=mask;
  levels=(float *) malloc(kmax*sizeof(float));
  header.levels=levels;

  status=ascii_writeheader (out, header);
  if(status!=0) goto error;

  shift=mjd(date.month,date.day,date.year)-mjd(1,1,1950);

  for(t=minstep;t<=maxstep;t++) {
    status=fe_framemapr2(input, info, grid, mask, element, t, bufx, bufy, &time);
    if(status!=0) goto error;
    time=time/d2s+shift;
    status=ascii_save_r (out, t, grid,  bufx, mask, time);
    if(status!=0) goto error;
    status=ascii_save_r (out, t, grid,  bufy, mask, time);
    if(status!=0) goto error;
    }


  fclose(out);
  free(element);
  free(bufx);
  free(bufy);
  return(0);

error:
  fclose(out);
  free(element);
  free(bufx);
  free(bufy);
  return(-1);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
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
    "  %s OPTIONS\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Converts TUGOm outputs. USED BY CTOH.\n"
    "  If no grid is given with either -z or -g, produces unstructured outputs.\n"
    "Otherwise, interpolates to structured grid.\n"
    "  Input format is TUGOm's binary, or, if no binary file can be found, TUGOm's NetCDF.\n"
    "\n"
    "OPTIONS :\n"
    "  --help  Show this help and exit.\n"
    "  -pressure  include pressure in the output\n"
    "  -currents  include currents in the output\n"
    "  -ibd  include inverse barometer in the output\n"
    "  -s1  enable detiding of S1\n"
    "  -s2  enable detiding of S2\n"
    "  -sa  enable detiding of Sa\n"
    "  -cm  scale to centimeters\n"
    "  -p  followed by root name for forcing files. Default: forcing-\n"
    "  -h  followed by root name for analysis files. Default: analysis-\n"
    "      When used as last argument, show this help and exit.\n"
    "  -d  followed by input directory. Default: .\n"
    "  -z  followed by output grid zone name\n"
    "  -g  followed by output grid file\n"
    "  -s  followed by start date in mm/yyyy format\n"
    "  -e  followed by end date in mm/yyyy format\n"
    "  -o  followed by output directory. Default: .\n"
    "  -f  followed by output format: cdf(default) or ascii\n"
    "  -r  input files are reduced and their names end with .reduced-6\n"
    "  -i  followed by frame stride. Default: 1\n"
    ); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *buf=0,*bufx=0,*bufy=0,mask=1.e+35;

  double time;
  int month,year;
  int frame;
  int ncid;
  int *element=0,m;

  int i,k,minstep,maxstep,nndes,stride=1,n,status,count;
  char *h_root=NULL,*p_root=NULL,*zone=NULL,*fmt=NULL;
  const char *keyword,*next_arg;
  char *hfile=0,*pfile=0,*nfile=0,*pathname=NULL,*outpath=NULL;
  char *gridfile=NULL;
  char filename[1024]="\0",tidefile[1024]="\0",*s;
  grid_t grid;
  date_t start_date,end_date,actual;
  start_date.year=end_date.year=0;
  meta_archive_t hinfo,pinfo;
  int reduced=0;
  
  int options[create_ncfile2d_OPTIONS];
  options[0]=1;
  aset(&options[1],create_ncfile2d_OPTIONS-1,0);
  
  int detide_s2=0,detide_s1=0,detide_sa=0;
  fcomplex **constants=NULL;
  double  *tide=NULL,cnestime;
  double scale=1.0;
  float *buffer[2][3]={{NULL,NULL,NULL},{NULL,NULL,NULL}};
  spectrum_t spectrum;
  FILE *in1,*in2;
  statistic_t check;

  /*warning, not initialzed*/
  int time_id;
  int h_id;
  int u_id;
  int v_id;
  int p_id;
  int ibd_id;

  fct_echo( argc, argv);

  spectrum.n=0;
  n=1;

  while (n < argc) {
    keyword=strdup(argv[n]);
    
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      exit(0);
      }
    
    if(strcmp(keyword,"-pressure")==0) {
      options[PRESSURE]=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-currents")==0) {
      options[CURRENTS]=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-ibd")==0) {
      options[IBD]=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-s2")==0) {
      detide_s2=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-s1")==0) {
      detide_s1=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-sa")==0) {
      detide_sa=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-cm")==0) {
      scale=0.01;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
/*------------------------------------------------------------------------------
        rootname for  model simulations forcing files */
        case 'p' :
          p_root= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        rootname for  model simulations sea state files */
        case 'h' :
          next_arg=argv[n+1];
          if(next_arg==0){
            print_help(argv[0]);
            exit(-1);
            }
          h_root= strdup(next_arg);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        pathname for model simulations files */
        case 'd' :
          pathname= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        region of extraction [1] */
        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        grid of extraction   [2]*/
        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        starting date */
        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d",&start_date.month,&start_date.year);
          free(s);
          break;

/*------------------------------------------------------------------------------
        ending date */
        case 'e' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d",&end_date.month,&end_date.year);
          free(s);
          break;

/*------------------------------------------------------------------------------
        pathname for extraction files */
        case 'o' :
          outpath= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        extraction file format */
        case 'f' :
          fmt= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        input are reduced archives */
        case 'r' :
          reduced=1;
          n++;
          break;

/*------------------------------------------------------------------------------
        extraction interval (2 = 1frame over 2)*/
        case 'i' :
          sscanf(argv[n+1],"%d",&stride);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        STDOUT_BASE_LINE("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
    
    }
  
  if(isnad(start_date) or isnad(end_date)){
    printf("*** PLEASE SPECIFY BOTH START AND END DATES WITH -s AND -e ***\n");
    print_help(argv[0]);
    exit(-1);
    }

  if(fmt==NULL) {
    printf("output archive format default: netcdf\n");
    fmt= strdup("cdf");
    }
  
  if(zone==0 and gridfile==0){
    printf("*** NO GRID SPECIFIED WITH EITHER -z OR -g : NO INTERPOLATION ***\n");
    zone=strdup("mesh");
    grid.free();
    }
  else if(zone != NULL) {
    bool identified;
    grid=get_zonegrid(zone, & identified);
    if (not identified) TRAP_ERR_EXIT(-1,"get_zonegrid(\"%s\",) error\n");
    status=map_completegridaxis_2(&grid);
    }
  else
    zone=strdup("grid");

  if(gridfile != NULL) status=cdf_loadvargrid (gridfile,0,&grid);
  
  const int gsize=grid.Hsize();
  const char * varname;
  
  if(gsize>0){
    exitIfNull(
      buf=new float[gsize]
      );
    
    exitIfNull(
      bufx=new float[gsize]
      );
    
    exitIfNull(
      bufy=new float[gsize]
      );
    }

  if(h_root==NULL) {
//     h_root= strdup("analysis.YYYY-MM");
    h_root= strdup("analysis-");
    printf("use %s as root name for sea state archive files\n", h_root);
    }
  else {
    printf("use %s as root name for sea state archive files\n",h_root);
    }

  if(p_root==NULL) {
    p_root= strdup("forcing-");
    printf("use <forcing> as root name for forcing archive files\n");
    }
  else {
    printf("use %s as root name for forcing archive files\n",p_root);
    }

  if(pathname==NULL) {
    pathname= strdup(".");
    printf("use <.> as default path for archive files\n");
    }
  else {
    printf("use %s as default path for archive files\n",pathname );
    }

 if(outpath==NULL) {
    outpath = strdup(".");
    printf("use <.> as default path for output files\n");
    }

  if(detide_s2==1) {
    status=tide_addwave(&spectrum, wS2);
    }

  if(detide_s1==1) {
    status=tide_addwave(&spectrum, wS1);
    }

  if(detide_sa==1) {
    status=tide_addwave(&spectrum, wSa);
    }

  for (year=start_date.year;year<=end_date.year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start_date.year) && (month < start_date.month)) continue;
      if((year==end_date.year)   && (month > end_date.month))  break;
      
      deletep_void(&hfile,free);
      deletep_void(&pfile,free);
      deletep_void(&nfile,free);
      
      if(reduced==0) {
        asprintf(&hfile,"%s/%s%4.4d.%2.2d",pathname,h_root,year,month);
        asprintf(&pfile,"%s/%s%4.4d.%2.2d",pathname,p_root,year,month);
        }
      else {
        asprintf(&hfile,"%s/%s%4.4d.%2.2d.reduced-6",pathname,h_root,year,month);
        asprintf(&pfile,"%s/%s%4.4d.%2.2d.reduced-6",pathname,p_root,year,month);
        }
      
      int haccess,paccess,naccess=-1;
      haccess=access(hfile,R_OK);
      paccess=access(pfile,R_OK);
      
      if(haccess!=0 and paccess!=0) {
        asprintf(&nfile,"%s.nc",hfile);
        naccess=access(nfile,R_OK);
        }
      
      if(naccess!=0){
        status=clx_archiveinfo(hfile,&hinfo);
        if(status!=0) TRAP_ERR_EXIT(status,"clx_archiveinfo(\"%s\",) error (%d %s)\n",hfile,status,strerror(status));
        status=clx_archiveinfo(pfile,&pinfo);
        if(status!=0) TRAP_ERR_EXIT(status,"clx_archiveinfo(\"%s\",) error (%d %s)\n",pfile,status,strerror(status));
        }
      else{
        STDERR_BASE_LINE("Fail-safing to %s as input file\n",nfile);
        status=cdf_archiveinfo(nfile,&hinfo);
        NC_TRAP_ERROR(wexit,status,1,"cdf_archiveinfo(\"%s\",) error",nfile);
        }
      
      status=fe_list(&hinfo.mesh);
      if(gsize>0){
        element=fe_scan_elements(hinfo.mesh,grid,0,1);
        if (element == 0) TRAP_ERR_EXIT(-1,"fe_scan_elements((\"%s\"),,) error\n",hfile);
        }

      nndes=hinfo.mesh.nvtxs;

      minstep=1;
      maxstep=hinfo.nframe;
      
      if(naccess!=0){
        in1=fopen(hfile,"rb");
        in2=fopen(pfile,"rb");
        
        if(in2==NULL) {
          options[PRESSURE]=0;
          options[IBD]=0;
          }
        }
      else
        in1=in2=0;
      
/*------------------------------------------------------------------------------
      !!! not secure if sea state and forcing files are from different FE mesh */
      for(k=0;k<3;k++) {
        for(i=0;i<2;i++) {
          buffer[i][k]=new float[nndes];
          }
        }
/*    !!! not secure if sea state and forcing files are from different FE mesh
------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
      !!! not secure if consecutive archive files are from different FE mesh */
      if((constants==NULL) && (spectrum.n !=0)) {
        status=initialize_omega(&spectrum);
        constants=(fcomplex **) malloc(spectrum.n*sizeof(fcomplex *));
        tide=(double *) malloc(nndes*sizeof(double));
        for(k=0;k<spectrum.n;k++) {
          constants[k]=(fcomplex *) malloc(nndes*sizeof(fcomplex));
          sprintf(tidefile,"%s-elevation-atlas.nc",spectrum.waves[k].name);
          printf("Trying climatology from %s\n",tidefile);fflush(stdout);
          status=poc_get_cvara(tidefile,"elevation_a","elevation_G",0,constants[k]);
          if(status!=0){
            NC_CHKERR_BASE_LINE(status,"poc_get_cvara error");
            sprintf(tidefile,"%s.ele.s2c",spectrum.waves[k].name);
            printf("Trying climatology from %s\n",tidefile);fflush(stdout);
            status=quoddy_loadc1(tidefile, nndes,constants[k]);
            }
          if(status!=0) TRAP_ERR_EXIT(status,"quoddy_loadc1(\"%s\",,) error (%d %s)\n",tidefile,status,strerror(status));
          }
        }
/*    !!! not secure if consecutive archive files are from different FE mesh
------------------------------------------------------------------------*/

      if(strcmp(fmt,"cdf")==0) {
        
        sprintf(filename,"%s/%s-%04d.%02d.huv.nc",outpath,zone,year,month);
        fprintf(stderr,"creating %s\n",filename);
        status=unlink(filename);
        if(status!=0 && errno!=ENOENT) TRAP_ERR_EXIT(errno,"unlink(\"%s\") error (%d %s)\n",filename,errno,strerror(errno));

        if(gsize>0){
          status=create_ncfile2d(filename,(size_t) grid.nx, (size_t) grid.ny, grid,options);
          status=nc_open(filename,NC_NOWRITE,&ncid);
          NC_TRAP_ERROR(wexit,status,1,"nc_open(\"%s\",NC_NOWRITE,) error\n",filename);

/*------------------------------------------------------------------------------
          get variable's id */
          status=nc_inq_varid(ncid,"time",&time_id);
          if(status!=0) time_id=-1;
          status=nc_inq_varid(ncid,"h",&h_id);
          if(status!=0) h_id=-1;
          status=nc_inq_varid(ncid,"u",&u_id);
          if(status!=0) u_id=-1;
          status=nc_inq_varid(ncid,"v",&v_id);
          if(status!=0) v_id=-1;
          status=nc_inq_varid(ncid,"p",&p_id);
          if(status!=0) p_id=-1;
          status=nc_inq_varid(ncid,"ibd",&ibd_id);
          if(status!=0) ibd_id=-1;
          
          status=nc_close(ncid);
          }
        else{
          time_id=0;
          }
        
        count=-1;
        for(frame=minstep;frame<=maxstep;frame+=stride) {
          count++;
          if(naccess!=0){
            status=clx_archiveread(in1, hinfo, frame, buffer[0], &time);
            if(status!=0) TRAP_ERR_EXIT(status,"clx_archiveread((\"%s\"),,%d,,) error\n",hfile,frame);
            }
          else{
            hinfo.flag=FLAG_STATE;
            status=cdf_archiveread(nfile, hinfo, frame-1, buffer[0], &time);
            NC_TRAP_ERROR(wexit,status,1,"cdf_archiveread(\"%s\",,%d,,) error",nfile,frame-1);
            }
          printf("#treating frame %d for %s\r",frame-1,sgetnewdate(hinfo.reference,time));fflush(stdout);
          
/*------------------------------------------------------------------------------
          treat elevations */
          if(options[ELEVATION]==1) {
            if(spectrum.n !=0) {
              status=harmonic_prediction(tide, cnestime, nndes, spectrum, constants, 0);
              check= get_statistics(buffer[0][0], 1.e+10, nndes,1);
              for(m=0;m<nndes;m++) buffer[0][0][m]=buffer[0][0][m]-tide[m]/scale;
              check= get_statistics(buffer[0][0], 1.e+10, nndes,1);
              }
            if(gsize>0){
              status=fe_map(hinfo.mesh,buffer[0][0],grid,element,buf,mask);
              if(status!=0) TRAP_ERR_EXIT(status,"fe_map() error %d\n",status);
/*------------------------------------------------------------------------------
              convert elevation in meters */
              for(m=0;m<grid.nx*grid.ny;m++) if(buf[m]!=mask) buf[m]*=scale;
              status=write_ncfile2d_obsolete(filename,  grid, count, h_id, buf);
              NC_TRAP_ERROR(wexit,status,1,"write_ncfile2d_obsolete(\"%s\",,,,) error",filename);
              }
            else{
              varname="h";
              status=poc_put_mesh_vara(filename,hinfo.mesh,LGP1,varname,count,buffer[0][0],"m",0);
              NC_TRAP_ERROR(wexit,status,1,"poc_put_mesh_vara(\"%s\",,LGP1,\"%s\",%d,...) error",filename,varname,count);
              }
            }
          
          actual=poctime_getdate(time,hinfo.reference,'s', &cnestime);
          status=poc_writetime(filename, count, time_id,  cnestime*24.*3600.);
          if(status!=0) TRAP_ERR_EXIT(status,"poc_writetime(\"%s\",,,) error\n",filename);
          
/*------------------------------------------------------------------------------
          treat currents */
          if(options[CURRENTS]==1 and gsize>0) {
            status=fe_map(hinfo.mesh,buffer[0][1],grid,element,bufx,mask);
            if(status!=0) TRAP_ERR_EXIT(status,"fe_map() error %d\n",status);
            status=fe_map(hinfo.mesh,buffer[0][2],grid,element,bufy,mask);
            if(status!=0) TRAP_ERR_EXIT(status,"fe_map() error %d\n",status);
/*------------------------------------------------------------------------------
            convert currents in meters/s */
            for(m=0;m<grid.nx*grid.ny;m++) if(bufx[m]!=mask) bufx[m]*=scale;
            for(m=0;m<grid.nx*grid.ny;m++) if(bufy[m]!=mask) bufy[m]*=scale;
            status=write_ncfile2d_obsolete(filename,  grid, count, u_id, bufx);
            NC_TRAP_ERROR(wexit,status,1,"write_ncfile2d_obsolete(\"%s\",,,,) error",filename);
            status=write_ncfile2d_obsolete(filename,  grid, count, v_id, bufy);
            NC_TRAP_ERROR(wexit,status,1,"write_ncfile2d_obsolete(\"%s\",,,,) error",filename);
            }

          if((options[PRESSURE]==0)&&(options[IBD]==0)) continue;
          
          if(naccess!=0){
            status=clx_archiveread(in2, pinfo, frame, buffer[1], &time);
            if(status!=0) TRAP_ERR_EXIT(status,"clx_archiveread((\"%s\"),,,) error\n",pfile);
            }
          else{
            hinfo.flag=FLAG_FORCING;
            status=cdf_archiveread(nfile, hinfo, frame-1, buffer[1], &time);
            NC_TRAP_ERROR(wexit,status,1,"cdf_archiveread(\"%s\",,%d,,) error",nfile,frame-1);
            }

/*------------------------------------------------------------------------------
          treat surface pressure */
          if(options[PRESSURE]==1) {
            if(gsize>0){
              status=fe_map(hinfo.mesh,buffer[1][0],grid,element,buf,mask);
              if(status!=0) TRAP_ERR_EXIT(status,"fe_map() error %d\n",status);
/*------------------------------------------------------------------------------
              leave pressure in HPa */
              status=write_ncfile2d_obsolete(filename,  grid, count, p_id, buf);
              NC_TRAP_ERROR(wexit,status,1,"write_ncfile2d_obsolete(\"%s\",,,,) error",filename);
              }
            else{
              varname="p";
              status=poc_put_mesh_vara(filename,hinfo.mesh,LGP1,varname,count,buffer[1][0],"m",0);
              NC_TRAP_ERROR(wexit,status,1,"poc_put_mesh_vara(\"%s\",,LGP1,\"%s\",%d,...) error",filename,varname,count);
              }
            }

/*------------------------------------------------------------------------------
          treat inverted barometer departure */
          if(options[IBD]==1 and gsize>0) {
            for(m=0;m<nndes;m++) buffer[1][0][m]=buffer[1][0][m]+buffer[0][0][m];
            status=fe_map(hinfo.mesh,buffer[1][0],grid,element,buf,mask);
            if(status!=0) TRAP_ERR_EXIT(status,"fe_map() error %d\n",status);
/*------------------------------------------------------------------------------
            leave ibd in cm */
            status=write_ncfile2d_obsolete(filename,  grid, count, ibd_id, buf);
            NC_TRAP_ERROR(wexit,status,1,"write_ncfile2d_obsolete(\"%s\",,,,) error",filename);
            }
          
          }

        if(in1 != NULL) fclose(in1);
        if(in2 != NULL) fclose(in2);
        
        fprintf(stderr,"%s",el);
        }

/*------------------------------------------------------------------------------
      Format obsolete !!!
      if(strcmp(fmt,"bmg")==0) {
        fprintf(stderr,"create the elevation sequence gridded file: %2d%4d\n",month,year);
        sprintf(filename,"%s/%s-%4.4d.%2.2d.h.bmg",outpath,zone,year,month);
        status=fe_filemapr1(hfile,filename,grid,minstep,maxstep,info);

        fprintf(stderr,"create the velocity sequence gridded file: %2d%4d\n",month,year);
        sprintf(filename,"%s/%s-%4.4d.%2.2d.u.bmg",outpath,zone,year,month);
        status=fe_filemapr2 (hfile,filename,grid,minstep,maxstep,info);
  
        fprintf(stderr,"create the pressure sequence gridded file: %2d%4d\n",month,year);
        sprintf(filename,"%s/%s-%4.4d.%2.2d.p.bmg",outpath,zone,year,month);
        status=fe_filemapr1(pfile,filename,grid,minstep,maxstep,info);

        fprintf(stderr,"create the wind stress sequence gridded file\n");
        sprintf(filename,"%s/%s-%4.4d.%2.2d.w.bmg",outpath,zone,year,month);
        status=fe_filemapr2 (pfile,filename,grid,minstep,maxstep,info);
        }
*/

      if(strcmp(fmt,"ascii")==0) {
        fprintf(stderr,"create the elevation sequence gridded file: %2d%4d\n",month,year);
        sprintf(filename,"%s/%s-%4.4d.%2.2d.h.asc",outpath,zone,year,month);
        status=filemapr1(hfile,filename,grid,minstep,maxstep,hinfo);

        fprintf(stderr,"create the velocity sequence gridded file: %2d%4d\n",month,year);
        sprintf(filename,"%s/%s-%4.4d.%2.2d.u.asc",outpath,zone,year,month);
        status=filemapr2 (hfile,filename,grid,minstep,maxstep,hinfo);
/*

        fprintf(stderr,"create the pressure sequence gridded file: %2d%4d\n",month,year);
        sprintf(filename,"%s/%s-%4.4d.%2.2d.p.asc",outpath,zone,year,month);
        status=filemapr1(pfile,filename,grid,minstep,maxstep,info);

        fprintf(stderr,"create the wind stress sequence gridded file\n");
        sprintf(filename,"%s/%s-%4.4d.%2.2d.w.bmg",outpath,zone,year,month);
        status=fe_filemapr2 (pfile,filename,grid,minstep,maxstep,info);
*/
        }

        archive_freeinfo(&hinfo);
        archive_freeinfo(&pinfo);
        for(k=0;k<3;k++) {
          for(i=0;i<2;i++) deletep(&buffer[i][k]);
          }
      }
    }

  STDOUT_BASE_LINE("end of convert ... \n");
  exit(0);
}
