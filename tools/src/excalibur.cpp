#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "tools-structures.h"
    
#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "poc-time.h"
#include "fe.h"
#include "archive.h"
#include "datastream.h"
#include "imbrication.h"

typedef struct
   {
     char    name[60];
     double  lon,lat,z;
   } maregraphe_t;

extern int identify_BINARY(const char *file, const char *name, int *id, int verbose);
extern int stream_UGxBINARY(UGfield_t<double> & field, double t);

#define ARCHIVE_FORMAT_ASCII  0
#define ARCHIVE_FORMAT_NETCDF 1
#define ARCHIVE_FORMAT_BINARY 2

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_intpl_xyt2(char * datafile, meta_archive_t info, double t, double p, double time,float serie[3], memory_t *memory)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double beta[3];
  float  z,z1,z2;
  float  count,start,end;
  double time1,time2,sampling=0.0,t0;
  int origine, nodes[3];
  int i,j,k,node;
  int status;
  int nndes,elt,read;
  date_t reference;

  FILE *in;
  
  reference=info.reference;

  origine=julian_day(1,1,1950);
  t0=(julian_day(reference.month,reference.day,reference.year)-origine)*24.*3600.+reference.second;

  start=t0+info.start;

  status=fe_beta(info.mesh, t, p,&elt,nodes,beta);
  if(status!=0) return(-1);

  nndes=info.mesh.nvtxs;

  read=1;

  if(memory->file==NULL) {
/*     memory->h1=smatrix(0,2,0,nndes-1); */
/*     memory->h2=smatrix(0,2,0,nndes-1); */
    goto skip;
    }

  if(strcmp(memory->file,datafile)==0) {
    t=time*24.*3600.;
    i=(int)( (t-start)/info.sampling );
    if(i==memory->frame) read=0;
    else read=1;
    }
   
    
skip:
   if (read) {
     in=fopen(datafile,"r");
     if(memory->file!=NULL) free(memory->file);
     memory->file=strdup(datafile);
     }
     
    t=time*24.*3600.;
    i=(int)( (t-start)/info.sampling);
    if (read) {
      status=clx_archiveread(in,info,i+1,memory->h1,&memory->time1);
      if(status!=0) goto error;
      status=clx_archiveread(in,info,i+2,memory->h2,&memory->time2);
      if(status!=0) goto error;
      memory->frame=i;
      }
/*     printf("t= %6.1f hours status= %d %s\n",time1/3600.0,status, sgetnewdate(info.reference,time1)); */
/*     printf("t= %6.1f hours status= %d %s\n",time2/3600.0,status, sgetnewdate(info.reference,time2)); */
    t=t-t0;
    for (j=0;j<3;j++) {
      z1=0;
      z2=0;
      for (i=0;i<3;i++) {
        z1+=memory->h1[j][nodes[i]-1]*beta[i];
        z2+=memory->h2[j][nodes[i]-1]*beta[i];
        }
      serie[j] = ((memory->time2-t)*z1+(t-memory->time1)*z2)/info.sampling;
      }
 
  if (read) fclose(in);
  return(0);

  error:

  fclose(in);
  return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_datastream(double time, double x, double y, int format, string & path, string & convention, sequence_t< imbrication_t <double> > & sequence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
 
  provide simulation extraction for a moving probe
  
  As time discretisation is discontinuous, and geographical is continuous, time
  interpolation must come first 
  
  - get lon,lat,time
  
  - load simulation adequate frames if required
  
    -> build archive filename
    -> update physical field
  
  - perform time interplation 
  
  - perform geographical interpolation
 
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int nprobes=1;
  
  int status;
  mesh_t mesh;
  discretisation_t *descriptor=new discretisation_t;
  int discretisation;
  int nframes, count;
  nesting_t<double> *nesting=new nesting_t<double>;
  meta_archive_t info;
  int (*processor) (UGfield_t<double> & field, double t);
  int (*identify)(const char *, const char *, int *, int);
  
  switch(format) {
    case ARCHIVE_FORMAT_BINARY:
      processor=stream_UGxBINARY;
      identify=identify_BINARY;
      break;
    case ARCHIVE_FORMAT_NETCDF:
      processor=stream_UGxNETCDF;
      identify=0;
      break;
    default:
      TRAP_ERR_EXIT(-1, "illegal format\n\n");
      break;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialise imbrication file stream (use to build archive filename from a
  given convention and time, and select target variable)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  filestream_t<double> *filestream=new filestream_t<double> (path.c_str(),convention.c_str(),"elevation",identify);
  status=filestream->finalize(time);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  read imbrication mesh (should be done later in datastream-fashion way)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  switch(format) {
    case ARCHIVE_FORMAT_BINARY:
//       status=archive_readheader(filestream->filename.c_str(), &info);
      status=clx_archivereadheader(filestream->filename.c_str(), &info);
      if(status != 0) {
        return(-1);
        }
      status=fe_list(&(info.mesh));
      mesh=info.mesh;
      break;
    case ARCHIVE_FORMAT_NETCDF:
      status=fe_readmesh3d(filestream->filename.c_str(), &mesh, 0);
      if(status!=0) TRAP_ERR_EXIT(-1, "error by reading mesh in %s \n",filestream->filename.c_str());
      break;
    default:
      TRAP_ERR_EXIT(-1, "illegal format\n\n");
      break;
    }
    
  if(status != 0) {
    return(-1);
    }
    
/*-----------------------------------------------------------------------------
  initialise descriptor */
  discretisation=LGP1;
  
  status=discretisation_init(&mesh,discretisation,false,0);
  *descriptor=get_descriptor(mesh,discretisation);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialise nesting; probes can be occasionally out of the mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nesting->allocate(nprobes, descriptor, 0);
  
  count=0;
  for(int i = 0; i < nprobes; i++) {
//     double x,y;
    int hint=-1;
    status = fe_beta(mesh, x, y, hint, &(nesting->basics[i].element), nesting->basics[i].nodes, nesting->basics[i].beta);
    if(status != 0) {
      printf("probe position (%6d) at %lf %lf outside of UG domain\n", i, x, y);
      }
    else count++;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialise elevation UG stream (use to load a UG field at the appropriate
  archive time frame)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  fieldstream_t< UGfield_t<double> > *fieldstream;
  fieldstream=new fieldstream_t< UGfield_t<double> >;
  fieldstream->data.stream=filestream;
  fieldstream->data.descriptor=descriptor;
  
  fieldstream->data.allocate();
  fieldstream->data.processor=processor;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialise elevation sequence (use to interpolate the UG filed at a given
  <position, time> location along probe trajectory)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  nframes=2;
  status=sequence.allocate(nframes);
  
  for(int k=0;k<sequence.nframes;k++) {
/*------------------------------------------------------------------------------
    stream is a generic imbrication datastream*/
    imbrication_t <double> *stream= new imbrication_t <double> (nprobes, nesting, (datastream_t<double> *) fieldstream, stream_imbricationxUG);
    sequence.frames[k]=*stream;
    }
  status=sequence.init(time);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_imbricationxUG_03(imbrication_t<double> & imbrication, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  double  *data,mask;
  
  sequence_t <UGfield_t<double> > *stream;
  stream=( sequence_t <UGfield_t<double> > *) imbrication.stream;
  
/*------------------------------------------------------------------------------
  update stream, i.e. source unstructured field*/
  status=stream->check(t);
  if(status!=0) return(-1);
  
  size_t size=stream->frames[0].nvalues();
  
  data=new double[size];
  
  status=stream->interpolate(t, data);
  
/**-----------------------------------------------------------------------------
  then interpolate on model unstructured field*/
  for(n = 0; n < imbrication.nvalues; n++) {
    imbrication.x[n] = imbrication.nesting->basics[n].interpolate(data, mask, (double) 0);
    }
  
  delete[] data;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_datastream_new(double time, double x, double y, int format, string & path, string & convention,imbrication_t<double>* & imbrication)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
 
  provide simulation extraction for a moving probe
  
  As time discretisation is discontinuous, and geographical is continuous, time
  interpolation must come first 
  
  - get lon,lat,time
  
  - load simulation adequate frames if required
  
    -> build archive filename
    -> update physical field
  
  - perform time interplation 
  
  - perform geographical interpolation
 
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int nprobes=1;
  sequence_t< UGfield_t <double> > *sequence=new sequence_t< UGfield_t <double> >;
  
  int status;
  mesh_t *mesh=new mesh_t;
  discretisation_t *descriptor=new discretisation_t;
  int discretisation;
  int nframes, count;
  nesting_t<double> *nesting=new nesting_t<double>;
  meta_archive_t info;
  int (*processor) (UGfield_t<double> & field, double t);
  int (*identify)(const char *, const char *, int *, int);
  
  switch(format) {
    case ARCHIVE_FORMAT_BINARY:
      processor=stream_UGxBINARY;
      identify=identify_BINARY;
      break;
    case ARCHIVE_FORMAT_NETCDF:
      processor=stream_UGxNETCDF;
      identify=0;
      break;
    default:
      TRAP_ERR_EXIT(-1, "illegal format\n\n");
      break;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialise imbrication file stream (use to build archive filename from a
  given convention and time, and select target variable)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  filestream_t<double> *filestream=new filestream_t<double> (path.c_str(),convention.c_str(),"elevation",identify);
  status=filestream->finalize(time);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  read imbrication mesh (should be done later in datastream-fashion way)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  switch(format) {
    case ARCHIVE_FORMAT_BINARY:
//       status=archive_readheader(filestream->filename.c_str(), &info);
      status=clx_archivereadheader(filestream->filename.c_str(), &info);
      if(status != 0) {
        return(-1);
        }
      status=fe_list(&(info.mesh));
      *mesh=info.mesh;
      break;
    case ARCHIVE_FORMAT_NETCDF:
      status=fe_readmesh3d(filestream->filename.c_str(), mesh, 0);
      if(status!=0) TRAP_ERR_EXIT(-1, "error by reading mesh in %s \n",filestream->filename.c_str());
      break;
    default:
      TRAP_ERR_EXIT(-1, "illegal format\n\n");
      break;
    }
    
  if(status != 0) {
    return(-1);
    }
    
/*-----------------------------------------------------------------------------
  initialise descriptor */
  discretisation=LGP1;
  
  status=discretisation_init(mesh,discretisation,false,0);
  *descriptor=get_descriptor(*mesh,discretisation);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialise nesting; probes can be occasionally out of the mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nesting->allocate(nprobes, descriptor, mesh);
  
  count=0;
  for(int i = 0; i < nprobes; i++) {
//     double x,y;
    int hint=-1;
    status = fe_beta(*mesh, x, y, hint, &(nesting->basics[i].element), nesting->basics[i].nodes, nesting->basics[i].beta);
    nesting->basics[i].lon=x;
    nesting->basics[i].lat=y;
    if(status != 0) {
      printf("probe position (%6d) at %lf %lf outside of UG domain\n", i, x, y);
      }
    else count++;
    }

// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//  
//   initialise elevation UG stream (use to load a UG field at the appropriate
//   archive time frame)
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//  
//   fieldstream_t< UGfield_t<double> > *fieldstream;
//   fieldstream=new fieldstream_t< UGfield_t<double> >;
//   fieldstream->data.stream=filestream;
//   fieldstream->data.descriptor=descriptor;
//   
//   fieldstream->data.allocate();
//   fieldstream->data.processor=processor;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialise elevation sequence (use to interpolate the UG filed at a given
  time)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  nframes=2;
  status=sequence->allocate(nframes);
  
  for(int k=0;k<sequence->nframes;k++) {
/*------------------------------------------------------------------------------
    stream is a generic imbrication datastream*/
    UGfield_t <double> *frame= new UGfield_t <double> (path.c_str(), convention.c_str(), "elevation", descriptor, processor, identify);
    sequence->frames[k]=*frame;
    }
  status=sequence->init(time);
  
  imbrication= new imbrication_t <double> (nprobes, nesting, (datastream_t<double> *) sequence, stream_imbricationxUG_03);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int invoke(string & inputfile, string & outputfile, string & path, string & convention, int & format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *input, *output;
  int cgps=16,clon=15,clat=14,ctime=0;
  
  sequence_t< imbrication_t <double> > sequence; 
  
  imbrication_t <double> *imbrication;
  
  int nitems, status;
  double lon,lat,time,gps,lonp,latp, d;
  char *s, *line;
  vector<string> tokens,keys,values;
  string delimiter;
  size_t count;
   
  input =fopen(inputfile.c_str(),"r");
  output=fopen(outputfile.c_str(),"w");
 
  line=new char[1024];
  
  d=0;
  count=0;
  do {
    s=fgets(line,1024,input);
    if(s==0) break;
/*------------------------------------------------------------------------------
    warning: '\n' will be accounted as token */
    nitems=count_token(line);
    delimiter=" ";
    tokens=string_split(line, delimiter);
    nitems=tokens.size();
    nitems=sscanf(tokens[clon].c_str(), "%lf",&lon);
    nitems=sscanf(tokens[clat].c_str(), "%lf",&lat);
    nitems=sscanf(tokens[ctime].c_str(),"%lf",&time);
    nitems=sscanf(tokens[cgps].c_str(), "%lf",&gps);
    
//     if(count==0) status=initialise_datastream(time, lon, lat, format, path, convention,sequence);
//     else status=sequence.check(time);
    
    if(count==0) status=initialise_datastream_new(time, lon, lat, format, path, convention, imbrication);
    vector<double> x, y;
    x.push_back(lon);
    y.push_back(lat);
    
//     sequence.frames[0].locate(x,y);
//     double z=sequence.interpolate(time,0,status);
    imbrication->locate(x,y);
    imbrication->check(time);
    double z=imbrication->x[0];
//     printf("%lf %lf %lf %lf\n", time, lon, lat, z);
    if(count==0) {
      lonp=lon;
      latp=lat;
      }
    d+=geo_haversin(lon,lat,lonp,latp);
    fprintf(output, "%lf %lf %lf %lf %lf %lf %lf\n", time, lon, lat, z, gps, gps-z,d);
    count++;
    tokens.clear();
    lonp=lon;
    latp=lat;
    } while (true);

  delete[] line;
  
  fclose(input);
  fclose(output);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  
  string path, convention, trajectory, outputfile, formatname;
  int format;
  char *keyword;
  
  fct_echo(argc,argv);
  
  format=ARCHIVE_FORMAT_NETCDF;
  path="/home/models/seine/fleuve_trace_alti/simulation-tides+rivers+DAC/calibration-2017";
  convention="tugo.state-YYYY.MM.nc";
  
  trajectory="/home/data/sealevel/seine/nappe-GPS/cngc170622_all_air_fran_edit_ssv.track_fmt_300_lis";
  outputfile="/home/data/sealevel/seine/nappe-GPS/cngc170622.mod";
  
//   trajectory="/home/data/sealevel/seine/nappe-GPS/cngc170623_all_air_fran_edit_ssv.track_fmt_300_lis";
//   outputfile="/home/data/sealevel/seine/nappe-GPS/cngc170623.mod";
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'c' :
          convention=argv[n+1];
          n++;
          n++;
          break;

        case 'p' :
          path=argv[n+1];
          n++;
          n++;
          break;

        case 'f' :
          formatname=argv[n+1];
          n++;
          n++;
          break;

        case 'o' :
          outputfile=argv[n+1];
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        trajectory=argv[n];
        n++;
        break;
        }
      free(keyword);
    }

  status=invoke(trajectory, outputfile, path, convention, format);

 end:
  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  fprintf(stderr,"%s -computation aborted ^^^^^^^^^\n",argv[0]);
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}

