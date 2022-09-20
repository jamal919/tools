
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Convert CM93 database to XYZ file.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <glob.h>

#include "poc-netcdf-data.hpp"
#include "polygones.h"
#include "topo.h"
#include "xyz.h"

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

##########################
CODE INSPIRED FROM OpenCPN
##########################

https://github.com/OpenCPN/OpenCPN/blob/master/src/cm93.cpp

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

static unsigned char Table_0[] =
{
0x0CD,0x0EA,0x0DC,0x048,0x03E,0x06D,0x0CA,0x07B,0x052,0x0E1,0x0A4,0x08E,0x0AB,0x005,0x0A7,0x097,
0x0B9,0x060,0x039,0x085,0x07C,0x056,0x07A,0x0BA,0x068,0x06E,0x0F5,0x05D,0x002,0x04E,0x00F,0x0A1,
0x027,0x024,0x041,0x034,0x000,0x05A,0x0FE,0x0CB,0x0D0,0x0FA,0x0F8,0x06C,0x074,0x096,0x09E,0x00E,
0x0C2,0x049,0x0E3,0x0E5,0x0C0,0x03B,0x059,0x018,0x0A9,0x086,0x08F,0x030,0x0C3,0x0A8,0x022,0x00A,
0x014,0x01A,0x0B2,0x0C9,0x0C7,0x0ED,0x0AA,0x029,0x094,0x075,0x00D,0x0AC,0x00C,0x0F4,0x0BB,0x0C5,
0x03F,0x0FD,0x0D9,0x09C,0x04F,0x0D5,0x084,0x01E,0x0B1,0x081,0x069,0x0B4,0x009,0x0B8,0x03C,0x0AF,
0x0A3,0x008,0x0BF,0x0E0,0x09A,0x0D7,0x0F7,0x08C,0x067,0x066,0x0AE,0x0D4,0x04C,0x0A5,0x0EC,0x0F9,
0x0B6,0x064,0x078,0x006,0x05B,0x09B,0x0F2,0x099,0x0CE,0x0DB,0x053,0x055,0x065,0x08D,0x007,0x033,
0x004,0x037,0x092,0x026,0x023,0x0B5,0x058,0x0DA,0x02F,0x0B3,0x040,0x05E,0x07F,0x04B,0x062,0x080,
0x0E4,0x06F,0x073,0x01D,0x0DF,0x017,0x0CC,0x028,0x025,0x02D,0x0EE,0x03A,0x098,0x0E2,0x001,0x0EB,
0x0DD,0x0BC,0x090,0x0B0,0x0FC,0x095,0x076,0x093,0x046,0x057,0x02C,0x02B,0x050,0x011,0x00B,0x0C1,
0x0F0,0x0E7,0x0D6,0x021,0x031,0x0DE,0x0FF,0x0D8,0x012,0x0A6,0x04D,0x08A,0x013,0x043,0x045,0x038,
0x0D2,0x087,0x0A0,0x0EF,0x082,0x0F1,0x047,0x089,0x06A,0x0C8,0x054,0x01B,0x016,0x07E,0x079,0x0BD,
0x06B,0x091,0x0A2,0x071,0x036,0x0B7,0x003,0x03D,0x072,0x0C6,0x044,0x08B,0x0CF,0x015,0x09F,0x032,
0x0C4,0x077,0x083,0x063,0x020,0x088,0x0F6,0x0AD,0x0F3,0x0E8,0x04A,0x0E9,0x035,0x01C,0x05F,0x019,
0x01F,0x07D,0x070,0x0FB,0x0D1,0x051,0x010,0x0D3,0x02E,0x061,0x09D,0x05C,0x02A,0x042,0x0BE,0x0E6
};

static unsigned char Encode_table[256];
static unsigned char Decode_table[256];

void CreateDecodeTable ( void )
{
      int i;
      for ( i=0 ; i < 256 ; i++ )
      {
            Encode_table[i] = Table_0[i] ^ 8;
      }

      for ( i=0 ; i < 256 ; i++ )
      {
            unsigned char a = Encode_table[i];
            Decode_table[ ( int ) a] = ( unsigned char ) i;
      }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int read_and_decode(FILE *f, T *x,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  /* inspired from read_and_decode_*() */
  
  const size_t size=sizeof(T);
  int count,i;
  
  count=fread(x,size,1,f);
  if(count!=1) TRAP_ERR_RETURN(count,verbose,"read error: returned %d!=1\n",count);
  
  if(verbose>0)STDERR_BASE_LINE_FUNC("size:%u\n",size);
  
  unsigned char
    *xc=(unsigned char*)x,
    *xci;
  for(i=0;i<size;i++){
    xci=&xc[i];
    *xci=Decode_table[*xci];
    }
  
  return count;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int skip(FILE *f,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  because if something is not used, there is :
    - no need to store it
    - and even less need to decode it !
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  const size_t size=sizeof(T);
  int count;
  T dummy;
  
  count=fread(&dummy,size,1,f);
  if(verbose>0)STDERR_BASE_LINE_FUNC("size:%u\n",size);
  
  return count;
}


#define PI M_PI
#define DEGREE d2r
const double CM93_semimajor_axis_meters        = 6378388.0;
typedef struct{
  double x_origin,y_origin,x_rate,y_rate;
  } transform_t;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void Transform(const transform_t & t,uint16_t sx,uint16_t sy,double trans_x, double trans_y, double *lat, double *lon)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  /* reworked from cm93chart::Transform() */
  //    Simple linear transform
  double valx = ( sx * t.x_rate ) + t.x_origin;
  double valy = ( sy * t.y_rate ) + t.y_origin;
  
  //    Add in the WGS84 offset corrections
  valx -= trans_x;
  valy -= trans_y;
  
  //    Convert to lat/lon
  *lat = ( 2.0 * atan ( exp ( valy/CM93_semimajor_axis_meters ) ) - PI/2. ) / DEGREE;
  *lon = ( valx / ( DEGREE * CM93_semimajor_axis_meters ) );
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_hydro(const char *hydro, int signus, const char *XYZfile, int persistence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,m,n,status;
  float *buffer,mask,pbma;
  grid_t grid;
  bool debug=false;
  double *x,*y,*z,dmask,xx,yy;
  string output;
  int ndata;

  status=topo_loadfield(hydro, &grid, &buffer, &mask, debug);
  if(status!=0) return(status);
  
  for (int k=0;k<persistence;k++) status=map_persistence(grid, buffer, mask, 0.0);
      
  status=xyz_loadraw (XYZfile, "", 0, x, y, z, &dmask, ndata, debug);

/*------------------------------------------------------------------------------
  Interpolate minimum low tide level and add to topo*/
  for(j=0;j<ndata;j++) {
    if(z[j]==dmask) {
      continue;
      }
    xx=x[j];
    yy=y[j];
    status=map_interpolation(grid, buffer,mask,xx,yy,&pbma);
    if (pbma!=mask) {
/*------------------------------------------------------------------------------
      assume negative depths and negative lowest astronomical tide*/
      z[j]+=signus*pbma;
      }
    }

  delete[] buffer;

  output=(string)XYZfile+(string)".msl";
  status=xyz_save (output.c_str(), x, y, z, dmask, ndata,0);

  return(0);
}


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

##############
MAIN FUNCTIONS
##############

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
////prints help of this programme
/**
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
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Convert CM93 database to XYZ file.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help : show this help and exit\n"
    "  -i : followed by input directory\n"
    "  -o : followed by output XYZ file\n"
    "  --lon : followed by 2 non-space-separated longitude boundary values\n"
    "  --lat : followed by 2 non-space-separated latitude boundary values\n"
    "  --pbma : followed by lowest astronomical tide level.\n"
    "    The name of the additional output will be that of the output file with '.msl' appended.\n"
    "  -negative : switch from positive depth to negative depth\n"
    ); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
  int i,j,k;/*< general indexes */
  int n,status;/*< argument index, status */
  bool negative=false;
  
  const char *keyword=0,*outPath=0;
  string inDir,globPattern,PBMA;
  glob_t globbed;
  globbed.gl_pathv=0;
  
  frame_t frame(-INFINITY,INFINITY,-90,90);
  
/*------------------------------------------------------------------------------
  parse arguments */
  string cmd=fct_echo(argc,argv);
  
  n=1;
  while (n < argc) {
    keyword=argv[n];
    
    if( strcmp(keyword,"-h")==0 or
        strcmp(keyword,"--help")==0 ){
      print_help(argv[0]);
      return(0);
      }
    
    if( strcmp(keyword,"--pbma")==0 ){
      PBMA=argv[n+1];
      n++;
      n++;
      continue;
      }
    
    if( strcmp(keyword,"-negative")==0 ){
      negative=true;
      n++;
      continue;
      }
    
    if( strcmp(keyword,"--lon")==0 ){
      sscanf(argv[n+1],"%lg%*c%lg",&frame.xmin,&frame.xmax);
      n++;
      n++;
      continue;
      }
    
    if( strcmp(keyword,"--lat")==0 ){
      sscanf(argv[n+1],"%lg%*c%lg",&frame.ymin,&frame.ymax);
      n++;
      n++;
      continue;
      }
    
    switch (keyword[0]) {
      
      case '-':
        
        switch (keyword[1]) {
          
          case 'i':
            inDir=argv[n+1];
            n++;
            n++;
            break;
          
          case 'o':
            outPath=argv[n+1];
            n++;
            n++;
            break;
          
          default:
            STDERR_BASE_LINE("unknown option %s\n",keyword);
            print_help(argv[0]);
            return(-1);
          }
        
        break;
      
      default:
        STDERR_BASE_LINE("unknown option %s\n",keyword);
        print_help(argv[0]);
        return(-1);
      }
    
    }
  
/*------------------------------------------------------------------------------
  check arguments */
  bool argsOK=true;
  
  if(inDir==""){
    fprintf(stderr,"*** please spectify input dir with -i ***\n");
    argsOK=false;
    }
  
  if(outPath==0){
    fprintf(stderr,"*** please spectify output file with -o ***\n");
    argsOK=false;
    }
  
  if(argsOK==false){
    print_help(argv[0]);
    return(-1);
    }
  
  FILE *outFile,*inFile;
  outFile=fopen(outPath,"w");
  if(outFile==0) TRAP_ERR_EXIT(errno,"fopen(\"%s\",\"w\") error (%d %s)\n",outPath,errno,strerror(errno));
  
  const bool
    isPlg= plg_find_format(outPath)==PLG_FORMAT_SCAN;
  
  i=inDir.length()-1;
  if(inDir[i]=='/') inDir.erase(i);
  
/*------------------------------------------------------------------------------
  search files */
  globPattern=inDir+"/*/?/*.?";
  printf("searching "+globPattern+"\n");
  status=glob(globPattern.c_str(),0,0,&globbed);
  switch(status){
  case GLOB_NOSPACE:
    TRAP_ERR_EXIT(ENOMEM,"out of memory\n");
  case GLOB_ABORTED:
    TRAP_ERR_EXIT(EIO,"read error\n");
  case GLOB_NOMATCH:
    TRAP_ERR_EXIT(-1,""+globPattern+" did not match\n");
    }
  printf("found %u files:\n",globbed.gl_pathc);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  read files
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  CreateDecodeTable();
  
  const char *inPath;
  int verbose=0;
  
  long int fileLength,totalLength,position;
  uint16_t headerLength;
  int32_t table1Length,table2Length;
  
  frame_t lonlat,eastnorth;
  uint16_t vectorCount;
  int32_t vectorPointCount,
    vectorDescriptorCount1,vectorDescriptorCount2;
  uint16_t point3DDescriptorCount;
  int32_t point3DCount;
  
  uint16_t count,x,y,z;
  uint32_t total;
  
  double delta_x;
  transform_t t;
  double lon,lat,depth;
  
  for(k=0;;k++){
/*------------------------------------------------------------------------------
    open file */
    inPath=globbed.gl_pathv[k];
    if(inPath==0)
      break;
    inFile=fopen(inPath,"r");
    if(inFile==0){
      fprintf(stderr,"fopen(\"%s\",\"r\") error (%d %s)\n",inFile,errno,strerror(errno));
      continue;
      }
    
/*------------------------------------------------------------------------------
    read file */

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  inspired from Ingest_CM93_Cell(), that calls:
    - read_header_and_populate_cib()
    - read_vector_record_table()
    - read_3dpoint_table()
  also searching for ``cm93_point_3d``.
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    if(verbose>=0) printf("[%d/%d]reading %s%s%c",k+1,globbed.gl_pathc,inPath,el,verbose>0?'\n':'\r');fflush(stdout);
    
    /* lengths */
    headerLength=0u;
    table1Length=table2Length=-1;
    
    fseek(inFile,0L,SEEK_END);
    fileLength=ftell(inFile);
    rewind(inFile);
    
    status=read_and_decode(inFile,&headerLength,0);
    status=read_and_decode(inFile,&table1Length,0);
    status=read_and_decode(inFile,&table2Length,0);
    
    totalLength=headerLength+table1Length+table2Length;
    
    if(headerLength!=138u or fileLength!=totalLength){
      if(verbose>=0) STDERR_BASE_LINE("error with %s:\n"
        "(%u==138u) or (%ld != %ld = %u + %d + %d)\n",
        inPath,headerLength,fileLength,totalLength,headerLength,table1Length,table2Length);
      goto close;
      }
    if(verbose>3) STDERR_BASE_LINE("%ld == (%ld = %u(==138) + %d + %d)\n",fileLength,totalLength,headerLength,table1Length,table2Length);
    
    /* header */
    
    read_and_decode(inFile,&lonlat.xmin,0);
    read_and_decode(inFile,&lonlat.ymin,0);
    read_and_decode(inFile,&lonlat.xmax,0);
    read_and_decode(inFile,&lonlat.ymax,0);
    
    if(verbose>0){
      STDERR_BASE_LINE("");
      lonlat.print(0);
      }
    
    if( not lonlat.intersect(frame) )
      goto close;
    
    if(isPlg and verbose<4){
      fprintf(outFile,"1 5\n1 %g %g\n2 %g %g\n3 %g %g\n4 %g %g\n5 %g %g\n",
        lonlat.xmin,lonlat.ymin,lonlat.xmin,lonlat.ymax,
        lonlat.xmax,lonlat.ymax,lonlat.xmax,lonlat.ymin,
        lonlat.xmin,lonlat.ymin);
      goto close;
      }
    
    read_and_decode(inFile,&eastnorth.xmin,0);
    read_and_decode(inFile,&eastnorth.ymin,0);
    read_and_decode(inFile,&eastnorth.xmax,0);
    read_and_decode(inFile,&eastnorth.ymax,0);
    
    read_and_decode(inFile,&vectorCount,0);
    read_and_decode(inFile,&vectorPointCount,0);
    read_and_decode(inFile,&vectorDescriptorCount1,0);
    read_and_decode(inFile,&vectorDescriptorCount2,0);
    
    read_and_decode(inFile,&point3DDescriptorCount,0);
    read_and_decode(inFile,&point3DCount,0);
    
    degree_recale(&lonlat.xmax,lonlat.xmin+180.);
    lonlat.dilatation(0.1);
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  Replaced:
    - header.northing_ → eastnorth.y
    - header.easting_ → eastnorth.x
    - pCIB->transform_ -> t.
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    delta_x = eastnorth.xmax - eastnorth.xmin;
    if(delta_x<0)
      delta_x += CM93_semimajor_axis_meters * 2.0 * PI;
    
    t.x_rate = delta_x / 65535;
    t.y_rate = ( eastnorth.ymax - eastnorth.ymin ) / 65535;
    
    t.x_origin = eastnorth.xmin;
    t.y_origin = eastnorth.ymin;
    
    /* skip vectors */
    
#if 0 /* safe way */
    fseek(inFile,headerLength,SEEK_SET);
    
    total=0;
    
    for(j=0;j<vectorCount;j++){
      read_and_decode(inFile,&count,0);
      if(verbose>1) STDERR_BASE_LINE("%u 2D points\n",count);
      total+=count;
      for(i=0;i<count;i++){
        skip<uint16_t>(inFile,0);/* x */
        skip<uint16_t>(inFile,0);/* y */
        }
      }
    
    position=ftell(inFile);
    if(verbose>0) STDERR_BASE_LINE("reached byte %ld (%ld,%ld) "
      "by skipping %u 2D points (%d)\n",
      position,position-headerLength-table1Length,
        headerLength+vectorCount*2L+vectorPointCount*4L,
      total,vectorPointCount);
    
    if(total!=vectorPointCount){
      if(verbose>=0) STDERR_BASE_LINE("error with %s: %d != %d\n",
        inPath,total,vectorPointCount);
      if(verbose>0) TRAP_ERR_EXIT(ENOEXEC,"testing\n");
      goto close;
      }
    
    if(position!=headerLength+vectorCount*2L+vectorPointCount*4L){
      if(verbose>=0) STDERR_BASE_LINE("error with %s: %ld != %ld\n",
        inPath,position,headerLength+vectorCount*2L+vectorPointCount*4L);
      if(verbose>0) TRAP_ERR_EXIT(ENOEXEC,"testing\n");
      goto close;
      }
#else /* way that is tested above */
    fseek(inFile,headerLength+vectorCount*2L+vectorPointCount*4L,SEEK_SET);
#endif
    
/*------------------------------------------------------------------------------
    read bathy */
    
    total=0;
    
    for(j=0;j<point3DDescriptorCount;j++){
      read_and_decode(inFile,&count,0);
      if(verbose>1) STDERR_BASE_LINE("%u 3D points\n",count);
      total+=count;
      for(i=0;i<count;i++){
        read_and_decode(inFile,&x,0);
        read_and_decode(inFile,&y,0);
        read_and_decode(inFile,&z,0);
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  inspired from cm93chart::BuildGeom()
  also searching for ``cm93_point_3d``.
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
        if(z>=12000)
          depth=z-12000.;
        else
          depth=z/10.;
        
        if(negative)
          depth=-depth;
        
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  inspired from cm93chart::CreateS57Obj()
  also searching for ``OGRPoint``.
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
        Transform(t,x,y,0,0,&lat,&lon);
        if(not lonlat.inside(lon,lat)){
          if(verbose>=0) STDERR_BASE_LINE("error with %s: point %d(%g;%g) out of boundaries\n",inPath,i,lon,lat);
          if(verbose>0) TRAP_ERR_EXIT(ENOEXEC,"testing\n");
          }
        
        if(not frame.inside(lon,lat))
          continue;
        
        if(isPlg){
          fprintf(outFile,"1 1\n1 %g %g\n",lon,lat);
          continue;
          }
        
        fprintf(outFile,"%g %g %g\n",lon,lat,depth);
        }
      }
    
    position=ftell(inFile);
    if(verbose>0) STDERR_BASE_LINE("reached byte %ld "
      "by reading %u 3D points (%d)\n",
      position,total,point3DCount);
    
    if(total!=point3DCount){
      if(verbose>=0) STDERR_BASE_LINE("error with %s: %d != %d\n",
        inPath,total,point3DCount);
      if(verbose>0) TRAP_ERR_EXIT(ENOEXEC,"testing\n");
      goto close;
      }
    
close:
    if(verbose>2) TRAP_ERR_EXIT(ENOEXEC,"testing\n");
    fclose(inFile);
    }
  
  fclose(outFile);
  printf("Finished convertion of %d files. Output is in %s%s\n",k,outPath,el);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  convert "zero-hydro related" depths to "mean level related" depths
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  if(PBMA!="") {
    int signus;
    if(negative) signus=+1;
    else         signus=-1;
    status= xyz_hydro(PBMA.c_str(), signus, outPath, 5);
//     if(status !=0) {
//       goto error;
//       }
    }
  
/*------------------------------------------------------------------------------
  clean-up */
  if(globbed.gl_pathv!=0)
    globfree(&globbed);
  
  return(0);
}
