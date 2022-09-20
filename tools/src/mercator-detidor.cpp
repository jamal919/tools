
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

/* *----------------------------------------------------------------------

Detide model outputs

----------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "tides.h"
#include "functions.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "matrix.h"

#define nstat 9


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *buffer[10],mask,*a,*G,*hmean=0;
  double tau;
  double time,duration;
  int count;
  int fmt,i,j,k,l;
  int n,status;
  int day,month,year,nday,loop,elapsed,frame;
  char *filename[3],*hroot,*proot,*keyword,*pathname=NULL,*s;
  char *varname=NULL,*convention=NULL;
  char *onde[11],out[256],test[256];
  char rootname[256]="\0",tmp[256],output[256],vname[256];
  grid_t grid;
  int nonde;
  date_t start,final,reference,start_date;
  cdfgbl_t data_info,grid_info[10];
  cdfvar_t info[10],detided[10],variable;
  spectrum_t WaveList;
  date_t actual;
  const int nrhs=3;
  size_t nvalues[nrhs];
  int nodal_corrections=0;
  astro_angles_t astro_angles;
  harmonic_t harmonic;
  hconstant_t **constants;
  pocgrd_t z_ncgrid, u_ncgrid, v_ncgrid;
  float *rbufx,*rbufy,rmask;
  grid_t zgrid,ugrid,vgrid;
  double t;

  fct_echo( argc, argv);

  i=0;
  //onde[i++]=strdup("Z0");
  onde[i]=NULL;

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'c' :
          convention= strdup(argv[n+1]);  /* directory */
          n++;
          n++;
          break;

        case 'p' :
          pathname= strdup(argv[n+1]);  /* directory */
          n++;
          n++;
          break;

        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&start.day,&start.month,&start.year);
          break;

       case 'f' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&final.day,&final.month,&final.year);
          break;

        case 'v' :
          varname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          n++;
          i=pos((char*)NULL,onde,11);
          for(;n<argc;i++,n++) {
            onde[i]= strdup(argv[n]);
            }
          onde[i]=NULL;
          nonde=i;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        __OUT_BASE_LINE__("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
      free(keyword);
    }

  if(pathname==NULL) pathname=strdup(".");

/* --------------- init data for harmonic analysis ------------------ */
/*ATTENTION: l'egalite de structures ne marche pas bien sur pccls ... */

  WaveList.init(initialize_tide(),onde);
  printf("# nbre d'ondes a analyser: %d \n",WaveList.n);

  for (i=0; i<WaveList.n; i++) {
    if(WaveList.waves[i].omega == 0.) {
      WaveList.waves[i].init();
      }
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", WaveList.waves[i].name,WaveList.waves[i].omega);
    }

  for (i=0; i<WaveList.n; i++) {
    for (j=i+1; j<WaveList.n; j++) {
      tau=deltaOmega2separation(WaveList.waves[j].omega-WaveList.waves[i].omega);
      printf ("wave: %10s %10s, separation: %9.3f days \n",
        WaveList.waves[i].name,WaveList.waves[j].name,tau);
      }
    }

/*-----------------------------------------------------------------------------
  */

  start.second=0;
  final.second=0;
  reference=start;
  init_argument(&astro_angles,reference);

  count=0;

/* *----------------------------------------------------------------------------
  Get grid dimension once and for all (assumes T,U and V grids have same size !!!) */
  status=decode_name(start, "T", convention, &(filename[0]));
  status= cdf_globalinfo(filename[0],&grid_info[0],0);
  status= cdf_varinfo  (filename[0],"sossheig",&(info[0]),0);
  status= poc_getgrid2d (filename[0], grid_info[0], info[0], &grid);

  status= poc_getgrid2d (filename[0], grid_info[0], info[0], &zgrid);

  status=decode_name(start, "U", convention, &(filename[1]));
  status= cdf_globalinfo(filename[1],&grid_info[1],0);
  status= cdf_varinfo  (filename[1],"sozocrtx",&(info[1]),0);
  status= poc_getgrid2d (filename[1], grid_info[1], info[1], &ugrid);

  status=decode_name(start, "V", convention, &(filename[2]));
  status= cdf_globalinfo(filename[2],&grid_info[2],0);
  status= cdf_varinfo  (filename[2],"somecrty",&(info[2]),0);
  status= poc_getgrid2d (filename[2], grid_info[2], info[2], &vgrid);

  aset(nvalues,nrhs,(size_t) grid.nx*grid.ny);

/* *-------------------------------------------------------------------------------------
  Initialize harmonic analysis */
  printAndGetSpectrumDetails(WaveList);
  harmonic=harmonic_start(WaveList, nvalues, nrhs);

/*-----------------------------------------------------------------------------
  allocate memory for vectors */
  for(k=0;k<3;k++) {
    buffer[k]=new float[grid.nx*grid.ny];
    if(buffer[k] == NULL) {
      printf("#memory allocation error for buffer N= %d \n",nvalues);
      goto error;
      }
    }

/* *-------------------------------------------------------------------------------------
  Compute harmonic matrix and RHS vectors */
  for (year=start.year;year<=final.year;year++) {
    for (month=1;month<=12;month++) {
      if ((year==start.year) && (month<start.month)) continue;
      if ((year==final.year) && (month>final.month)) break ;
      for (day=1;day<=poctime_dpm(month,year);day++) {
        if ((year==start.year) && (month=start.month) && (day<start.day)) continue ;
        if ((year==final.year) && (month=final.month) && (day>final.day)) break ;
        actual.year=year;
        actual.month=month;
        actual.day=day;
        actual.second=0.0;
        status=decode_name(actual, "T", convention, &(filename[0]));
        status= cdf_varinfo  (filename[0],"sossheig",&(info[0]),0);
        status=decode_name(actual, "U", convention, &(filename[1]));
        status= cdf_varinfo  (filename[1],"sozocrtx",&(info[1]),0);
        status=decode_name(actual, "V", convention, &(filename[2]));
        status= cdf_varinfo  (filename[1],"sozocrtx",&(info[2]),0);
        grid.free();
        status= poc_getgrid2d (filename[0], grid_info[0], info[0], &grid);
/*-------------------------------------------------------------------------------------
        read the output files (MOG2D format)*/
        printf("read %s %s %s ...\n",filename[0],filename[1],filename[2]);
        for(frame=0;frame<24;frame++) {
          for(k=0;k<3;k++) {
//            printf("read frame=%d in %s ...\n",frame,filename[k]);
            status= poc_getvar2d (filename[k], info[k].id, frame,(float *) buffer[k], &mask ,info[k]);
            if(status !=0) goto read_error;
            }
          t=grid.time[frame]+cnes_time(grid.origine,'s')-cnes_time(start,'s');
          harmonic_storage(t, harmonic, nodal_corrections, buffer,astro_angles);
          }
        }
      }
    }

/* *-------------------------------------------------------------------------------------
  HARMONIC ANALYSIS */
  duration= ellapsed_time( start,  final, 's');

  constants=harmonic_analysis_core(harmonic, duration,1);
  for(k=0;k<3;k++) {
    for(n=0;n<grid.nx*grid.ny;n++) {
       constants[k][n].set_complex();
       }
    }

//  goto atlas;

/* *-------------------------------------------------------------------------------------
  Detiding */
  for (year=start.year;year<=final.year;year++) {
    for (month=1;month<=12;month++) {
      if ((year==start.year) && (month<start.month)) continue;
      if ((year==final.year) && (month>final.month)) break ;
      for (day=1;day<=poctime_dpm(month,year);day++) {
        if ((year==start.year) && (month=start.month) && (day<start.day)) continue ;
        if ((year==final.year) && (month=final.month) && (day>final.day)) break ;
        actual.year=year;
        actual.month=month;
        actual.day=day;
        actual.second=0.0;
        status=decode_name(actual, "T", convention, &(filename[0]));
        status= cdf_varinfo  (filename[0],"sossheig",&(info[0]),0);
        status=decode_name(actual, "U", convention, &(filename[1]));
        status= cdf_varinfo  (filename[1],"sozocrtx",&(info[1]),0);
        status=decode_name(actual, "V", convention, &(filename[2]));
        status= cdf_varinfo  (filename[1],"sozocrtx",&(info[2]),0);
        printf("read %s %s %s  ...\n",filename[0],filename[1],filename[2]);
        grid.free();
        status= poc_getgrid2d (filename[0], grid_info[0], info[0], &grid);
        for(k=0;k<3;k++) {
          sprintf(test,"%s-%s",info[k].name,"detided");
          status= cdf_varinfo  (filename[k],test,&(detided[k]),0);
          if(status==-1){
            detided[k]=info[k];
            free(detided[k].name);
            detided[k].name=new char[strlen(info[k].name)+9];
            sprintf(detided[k].name,"%s-%s",info[k].name,"detided");
            status=create_ncvariable(filename[k], &(detided[k]));
            }
          }
/*-------------------------------------------------------------------------------------
        read the output files (MOG2D format)*/
        for(frame=0;frame<24;frame++) {
          for(k=0;k<3;k++) {
            status= poc_getvar2d (filename[k], info[k].id, frame,(float *) buffer[k], &mask ,info[k]);
            }
          if(status !=0) goto read_error;
          t=grid.time[frame]+cnes_time(grid.origine,'s')-cnes_time(start,'s');
          harmonic_correction(t,WaveList,constants,buffer,nvalues,nrhs,nodal_corrections,astro_angles);
          for(k=0;k<3;k++) {
            size_t *start, *count;
            start=new size_t[info[k].ndim];
            count=new size_t[info[k].ndim];
            start[0]=frame;
            count[0]=1;
            start[1]=0;
            count[1]=grid.ny;
            start[2]=0;
            count[2]=grid.nx;
            status=poc_write(filename[k], detided[k], start, count, buffer[k]);
            }
          }
        }
      }
    }

atlas:
/* *-------------------------------------------------------------------------------------
  Archive harmonic onstants */
  for (k=0;k<WaveList.n;k++) {

    sprintf(output,"%s.nc",WaveList.waves[k].name);
    status= poc_createfile(output);

    status=poc_sphericalgrid_xy(output,"ZHL",zgrid,&z_ncgrid);
    status=poc_sphericalgrid_xy(output,"XHL",ugrid,&u_ncgrid);
    status=poc_sphericalgrid_xy(output,"YHL",vgrid,&v_ncgrid);

    rbufx =new float[zgrid.nx*zgrid.ny];
    rbufy =new float[zgrid.nx*zgrid.ny];

    for(n=0;n<grid.nx*grid.ny;n++) rbufx[n]=constants[0][n].a[k];
    for(n=0;n<grid.nx*grid.ny;n++) rbufy[n]=constants[0][n].G[k];

    sprintf(vname,"Ha");
    poc_standardvariable_xy(&variable,vname,rmask,"m",1., 0.,"tidal_elevation_amplitude","tidal elevation amplitude",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufx);
    variable.destroy();

    sprintf(vname,"Hg");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_elevation_phase_lag","tidal elevation phase lag",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufy);
    variable.destroy();

    for(n=0;n<grid.nx*grid.ny;n++) rbufx[n]=constants[1][n].a[k];
    for(n=0;n<grid.nx*grid.ny;n++) rbufy[n]=constants[1][n].G[k];

    sprintf(vname,"Ua");
    poc_standardvariable_xy(&variable,vname,rmask,"m/s",1., 0.,"tidal_eastward_current_amplitude","tidal eastward current amplitude",vname,u_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, ugrid, variable.id,rbufx);
    variable.destroy();

    sprintf(vname,"Ug");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_eastward_current_phase_lag","tidal eastward current phase lag",vname,u_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, ugrid, variable.id,rbufy);
    variable.destroy();

    for(n=0;n<grid.nx*grid.ny;n++) rbufx[n]=constants[2][n].a[k];
    for(n=0;n<grid.nx*grid.ny;n++) rbufy[n]=constants[2][n].G[k];

    sprintf(vname,"Va");
    poc_standardvariable_xy(&variable,vname,rmask,"m/s",1., 0.,"tidal_northward_current_amplitude","tidal north current amplitude",vname,v_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, vgrid, variable.id,rbufx);
    variable.destroy();

    sprintf(vname,"Vg");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_northward_current_phase_lag","tidal north current phase lag",vname,v_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, vgrid, variable.id,rbufy);
    variable.destroy();

    delete[] rbufx;
    delete[] rbufy;
    }
 
/*----free vectors------------------------------------------*/
  for(k=0;k<3;k++) {
    delete[] buffer[k];
    }
  if(hmean!=0) delete[] hmean;

  goto end;

error:
  printf(" error in allocating memory  ...\n");
  read_error:
  __OUT_BASE_LINE__(" error in reading %s ...\n",filename);
  exit(-1);

end:
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
