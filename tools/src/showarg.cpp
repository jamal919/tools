
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Shows informations about a list of waves : Astronomic potential, pulsation, period, critical latitude, Doodson number and separation

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "tides.h"
#include "tides.def"
#include "poc-time.h"
#include "functions.h"


int compWaveOmega(const void *a,const void *b){
  const tidal_wave *a_=(const tidal_wave *)a,*b_=(const tidal_wave *)b;
  return a_->omega>b_->omega;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool is_astronomic_or_radiational(const tidal_wave & wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const bool
    astronomic= (wave.Ap!=0.),
    radiational=
      strcmp(wave.name,"Sa")==0 or
      strcmp(wave.name,"S1")==0,
    result=astronomic or radiational;
  return result;
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
    " %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTIONS] [ wave1 ... ]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Shows informations about a list of waves : Astronomic potential, pulsation, period, critical latitude, Doodson number and separation.\n"
    "  If no list of waves is given, show informations about all coded waves, sorted by increasing pulsation, but do not show waves separation.\n"
    "  The wave names are taken from the Darwin convention, described in " WAVELISTPARTIALREF "\n"
    "" WAVELISTURL "\n"
    "This convention is followed throughout this software.\n"
    "  It shows the Doodson number. A list of waves with their doodson number is available from the IHO at\n"
    IHO_LIST_URL "\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  show this help\n"
    "  --alias  followed by aliasing frequency in days. Not compatible with --combinations\n"
    "  --combinations  Show non-linear combinations. Not compatible with --alias\n"
    "  -s  followed by the start date in dd/mm/yyyy format\n"
    "  -f  followed by the end date in dd/mm/yyyy format\n"
    "  -a  show arguments periodically, by default every day\n"
    "  -i  followed by the daily frequency at which the arguments are shown. Enforced minimum is daily. Implies -a.\n"
    "\n"
    "EXAMPLE\n"
    "From https://www.aviso.altimetry.fr/en/missions/past-missions/jason-1.html , we have:\n"
    "  showarg --alias 9.9156\n"
    ); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double tau;

  int nonde=0,i,j,k;//number of waves, indexes
  int n,N;//argument index,number of days
  int day,month,year;
  date_t start=NADate,final=NADate,reference;
  const char *argvn,*argvn1;/*< option and its argument */
  char *onde[100],*doodson;/* list of waves */
  spectrum_t WaveList,basic;
  tidal_wave *wp,*wpj,*wpk;/* pointers to wave */
  double V,u,f,increment=1.;
  bool showArguments=false;
  bool showCombinations=false;
  double period,alias=NAN;
  const char *aliased="",*unaliased="";
  const double two_Omega_deg_h=two_Omega/dph2rps;//in deg/h
  double critical;//critical latitude in degrees

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    argvn=argv[n];
    
    if( strcmp(argvn,"-h")==0 or
        strcmp(argvn,"--help")==0 ){
      print_help(argv[0]);
      exit(0);
      }
    
    switch (argvn[0]) {
      case '-':
        switch (argvn[1]) {

        case '-' :
          if( strcmp(argvn,"--alias")==0 ){
            argvn1=argv[n+1];
            n++;
            n++;
            alias=atof(argvn1);
            }
          else if( strncmp("--combi",argvn)==0 ){
            showCombinations=true;
            n++;
            }
          else{
            STDOUT_BASE_LINE("unknown option %s\n",argvn);
            print_help(argv[0]);
            exit(-1);
            }
          break;

        case 's' :
          argvn1=argv[n+1];
          n++;
          n++;
          sscanf(argvn,"%d/%d/%d",&start.day,&start.month,&start.year);
          start.second=0.;
          break;

       case 'f' :
          argvn1=argv[n+1];
          n++;
          n++;
          sscanf(argvn1,"%d/%d/%d",&final.day,&final.month,&final.year);
          final.second=0.;
          break;

       case 'i' :
          argvn1=argv[n+1];
          n++;
          n++;
          sscanf(argvn1,"%lf",&increment);
          showArguments=true;
          break;

       case 'a' :
          n++;
          showArguments=true;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",argvn);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        i=0;
        while(n+i< argc) {
          onde[i]= strdup(argv[n+i]);
          i++;
          }
        onde[i]=NULL;
        nonde=i;
        n+=i;
        break;
      }
    }


/* --------------- init data for harmonic analysis ------------------ */
/*ATTENTION: l'egalite de structures ne marche pas bien sur pccls ... */

//   double e=0.0549;
//   double U=0.5582e-7;
//   
//   double a=0.9154*0.1759*MeanEarthRadius;
//   printf(" N2: %lf %lf\n", 3./4.*U*a, 3./4.*U*0.1759*MeanEarthRadius);
// 
//   a=0.9154*0.9085*MeanEarthRadius;
//   printf(" M2: %lf %lf\n", 3./4.*U*a, 3./4.*U*0.9085*MeanEarthRadius);
// 
//   a=0.1565*0.0786*MeanEarthRadius;
//   printf(" K2: %lf %lf\n", 3./4.*U*a, 3./4.*U*0.9085*MeanEarthRadius);

/*------------------------------------------------------------------------------
  initialise wave list */
  
  reference=CNES0;

  astro_angles_t astro_angles;
  basic=initialize_tide(&astro_angles,reference);
  
  if(nonde<=0){
    WaveList=basic;
    }
  else{
    WaveList.init(basic,onde);
    if(nonde!=WaveList.n)
      nonde=0;
    }
  
  WaveList.waves_init();
  
  if(WaveList.n+122<=WaveList.nmax){
extern tidal_wave wave_from_doodson(int tau, int s, int h, int p, int N_, int p1, float coef, const char *name);
    
    /* see set_w2nd() in src/prediction.c in libfes */
    WaveList.add(wave_from_doodson( 0, 0,  0,  0,  1,  0,  0.02793f , "libfes0"));
    WaveList.add(wave_from_doodson( 0, 0,  0,  0,  2,  0, -0.00027f , "libfes1"));
    WaveList.add(wave_from_doodson( 0, 0,  0,  2,  1,  0,  0.00004f , "libfes2"));
    WaveList.add(wave_from_doodson( 0, 0,  1,  0, -1, -1, -0.00004f , "libfes3"));
    WaveList.add(wave_from_doodson( 0, 0,  1,  0,  0, -1, -0.00492f , "libfes4"));
    WaveList.add(wave_from_doodson( 0, 0,  1,  0,  0,  1,  0.00026f , "libfes5"));
    WaveList.add(wave_from_doodson( 0, 0,  1,  0,  1, -1,  0.00005f , "libfes6"));
    WaveList.add(wave_from_doodson( 0, 0,  2, -2, -1,  0,  0.00002f , "libfes7"));
    WaveList.add(wave_from_doodson( 0, 0,  2, -2,  0,  0, -0.00031f , "libfes8"));
    WaveList.add(wave_from_doodson( 0, 0,  2,  0,  0,  0, -0.03095f , "libfes9"));
    WaveList.add(wave_from_doodson( 0, 0,  2,  0,  0, -2, -0.00008f , "libfes10"));
    WaveList.add(wave_from_doodson( 0, 0,  2,  0,  1,  0,  0.00077f , "libfes11"));
    WaveList.add(wave_from_doodson( 0, 0,  2,  0,  2,  0,  0.00017f , "libfes12"));
    WaveList.add(wave_from_doodson( 0, 0,  3,  0,  0, -1, -0.00181f , "libfes13"));
    WaveList.add(wave_from_doodson( 0, 0,  3,  0,  1, -1,  0.00003f , "libfes14"));
    WaveList.add(wave_from_doodson( 0, 0,  4,  0,  0, -2, -0.00007f , "libfes15"));
    WaveList.add(wave_from_doodson( 0, 1, -3,  1, -1,  1,  0.00002f , "libfes16"));
    WaveList.add(wave_from_doodson( 0, 1, -3,  1,  0,  1, -0.00029f , "libfes17"));
    WaveList.add(wave_from_doodson( 0, 1, -3,  1,  1,  1,  0.00002f , "libfes18"));
    WaveList.add(wave_from_doodson( 0, 1, -2, -1, -2,  0,  0.00003f , "libfes19"));
    WaveList.add(wave_from_doodson( 0, 1, -2, -1, -1,  0,  0.00007f , "libfes20"));
    WaveList.add(wave_from_doodson( 0, 1, -2,  1, -1,  0,  0.00048f , "libfes21"));
    WaveList.add(wave_from_doodson( 0, 1, -2,  1,  0,  0, -0.00673f , "libfes22"));
    WaveList.add(wave_from_doodson( 0, 1, -2,  1,  1,  0,  0.00043f , "libfes23"));
    WaveList.add(wave_from_doodson( 0, 1, -1, -1, -1,  1,  0.00002f , "libfes24"));
    WaveList.add(wave_from_doodson( 0, 1, -1, -1,  0,  1, -0.00021f , "libfes25"));
    WaveList.add(wave_from_doodson( 0, 1, -1, -1,  1,  1,  0.00000f , "libfes26"));
    WaveList.add(wave_from_doodson( 0, 1, -1,  0,  0,  0,  0.00020f , "libfes27"));
    WaveList.add(wave_from_doodson( 0, 1, -1,  1,  0, -1,  0.00005f , "libfes28"));
    WaveList.add(wave_from_doodson( 0, 1,  0, -1, -2,  0, -0.00003f , "libfes29"));
    WaveList.add(wave_from_doodson( 0, 1,  0, -1, -1,  0,  0.00231f , "libfes30"));
    WaveList.add(wave_from_doodson( 0, 1,  0, -1,  0,  0, -0.03518f , "libfes31"));
    WaveList.add(wave_from_doodson( 0, 1,  0, -1,  1,  0,  0.00228f , "libfes32"));
    WaveList.add(wave_from_doodson( 0, 1,  0,  1,  0,  0,  0.00189f , "libfes33"));
    WaveList.add(wave_from_doodson( 0, 1,  0,  1,  1,  0,  0.00077f , "libfes34"));
    WaveList.add(wave_from_doodson( 0, 1,  0,  1,  2,  0,  0.00021f , "libfes35"));
    WaveList.add(wave_from_doodson( 0, 1,  1, -1,  0, -1,  0.00018f , "libfes36"));
    WaveList.add(wave_from_doodson( 0, 1,  2, -1,  0,  0,  0.00049f , "libfes37"));
    WaveList.add(wave_from_doodson( 0, 1,  2, -1,  1,  0,  0.00024f , "libfes38"));
    WaveList.add(wave_from_doodson( 0, 1,  2, -1,  2,  0,  0.00004f , "libfes39"));
    WaveList.add(wave_from_doodson( 0, 1,  3, -1,  0, -1,  0.00003f , "libfes40"));
    WaveList.add(wave_from_doodson( 0, 2, -4,  2,  0,  0, -0.00011f , "libfes41"));
    WaveList.add(wave_from_doodson( 0, 2, -3,  0,  0,  1, -0.00038f , "libfes42"));
    WaveList.add(wave_from_doodson( 0, 2, -3,  0,  1,  1,  0.00002f , "libfes43"));
    WaveList.add(wave_from_doodson( 0, 2, -2,  0, -1,  0, -0.00042f , "libfes44"));
    WaveList.add(wave_from_doodson( 0, 2, -2,  0,  0,  0, -0.00582f , "libfes45"));
    WaveList.add(wave_from_doodson( 0, 2, -2,  0,  1,  0,  0.00037f , "libfes46"));
    WaveList.add(wave_from_doodson( 0, 2, -2,  2,  0,  0,  0.00004f , "libfes47"));
    WaveList.add(wave_from_doodson( 0, 2, -1, -2,  0,  1, -0.00004f , "libfes48"));
    WaveList.add(wave_from_doodson( 0, 2, -1, -1,  0,  0,  0.00003f , "libfes49"));
    WaveList.add(wave_from_doodson( 0, 2, -1,  0,  0, -1,  0.00007f , "libfes50"));
    WaveList.add(wave_from_doodson( 0, 2, -1,  0,  0,  1, -0.00020f , "libfes51"));
    WaveList.add(wave_from_doodson( 0, 2, -1,  0,  1,  1, -0.00004f , "libfes52"));
    WaveList.add(wave_from_doodson( 0, 2,  0, -2, -1,  0,  0.00015f , "libfes53"));
    WaveList.add(wave_from_doodson( 0, 2,  0, -2,  0,  0, -0.00288f , "libfes54"));
    WaveList.add(wave_from_doodson( 0, 2,  0, -2,  1,  0,  0.00019f , "libfes55"));
    WaveList.add(wave_from_doodson( 0, 2,  0,  0,  0,  0, -0.06662f , "libfes56"));
    WaveList.add(wave_from_doodson( 0, 2,  0,  0,  1,  0, -0.02762f , "libfes57"));
    WaveList.add(wave_from_doodson( 0, 2,  0,  0,  2,  0, -0.00258f , "libfes58"));
    WaveList.add(wave_from_doodson( 0, 2,  0,  0,  3,  0,  0.00007f , "libfes59"));
    WaveList.add(wave_from_doodson( 0, 2,  1, -2,  0, -1,  0.00003f , "libfes60"));
    WaveList.add(wave_from_doodson( 0, 2,  1,  0,  0, -1,  0.00023f , "libfes61"));
    WaveList.add(wave_from_doodson( 0, 2,  1,  0,  1, -1,  0.00006f , "libfes62"));
    WaveList.add(wave_from_doodson( 0, 2,  2, -2,  0,  0,  0.00020f , "libfes63"));
    WaveList.add(wave_from_doodson( 0, 2,  2, -2,  1,  0,  0.00008f , "libfes64"));
    WaveList.add(wave_from_doodson( 0, 2,  2,  0,  2,  0,  0.00003f , "libfes65"));
    WaveList.add(wave_from_doodson( 0, 3, -5,  1,  0,  1, -0.00002f , "libfes66"));
    WaveList.add(wave_from_doodson( 0, 3, -4,  1,  0,  0, -0.00017f , "libfes67"));
    WaveList.add(wave_from_doodson( 0, 3, -3, -1,  0,  1, -0.00007f , "libfes68"));
    WaveList.add(wave_from_doodson( 0, 3, -3,  1,  0,  1, -0.00012f , "libfes69"));
    WaveList.add(wave_from_doodson( 0, 3, -3,  1,  1,  1, -0.00004f , "libfes70"));
    WaveList.add(wave_from_doodson( 0, 3, -2, -1, -1,  0, -0.00010f , "libfes71"));
    WaveList.add(wave_from_doodson( 0, 3, -2, -1,  0,  0, -0.00091f , "libfes72"));
    WaveList.add(wave_from_doodson( 0, 3, -2, -1,  1,  0,  0.00006f , "libfes73"));
    WaveList.add(wave_from_doodson( 0, 3, -2,  1,  0,  0, -0.00242f , "libfes74"));
    WaveList.add(wave_from_doodson( 0, 3, -2,  1,  1,  0, -0.00100f , "libfes75"));
    WaveList.add(wave_from_doodson( 0, 3, -2,  1,  2,  0, -0.00009f , "libfes76"));
    WaveList.add(wave_from_doodson( 0, 3, -1, -1,  0,  1, -0.00013f , "libfes77"));
    WaveList.add(wave_from_doodson( 0, 3, -1, -1,  1,  1, -0.00004f , "libfes78"));
    WaveList.add(wave_from_doodson( 0, 3, -1,  0,  0,  0,  0.00006f , "libfes79"));
    WaveList.add(wave_from_doodson( 0, 3, -1,  0,  1,  0,  0.00003f , "libfes80"));
    WaveList.add(wave_from_doodson( 0, 3, -1,  1,  0, -1,  0.00003f , "libfes81"));
    WaveList.add(wave_from_doodson( 0, 3,  0, -3,  0,  0, -0.00023f , "libfes82"));
    WaveList.add(wave_from_doodson( 0, 3,  0, -3,  1, -1,  0.00004f , "libfes83"));
    WaveList.add(wave_from_doodson( 0, 3,  0, -3,  1,  1,  0.00004f , "libfes84"));
    WaveList.add(wave_from_doodson( 0, 3,  0, -1,  0,  0, -0.01275f , "libfes85"));
    WaveList.add(wave_from_doodson( 0, 3,  0, -1,  1,  0, -0.00528f , "libfes86"));
    WaveList.add(wave_from_doodson( 0, 3,  0, -1,  2,  0, -0.00051f , "libfes87"));
    WaveList.add(wave_from_doodson( 0, 3,  0,  1,  2,  0,  0.00005f , "libfes88"));
    WaveList.add(wave_from_doodson( 0, 3,  0,  1,  3,  0,  0.00002f , "libfes89"));
    WaveList.add(wave_from_doodson( 0, 3,  1, -1,  0, -1,  0.00011f , "libfes90"));
    WaveList.add(wave_from_doodson( 0, 3,  1, -1,  1, -1,  0.00004f , "libfes91"));
    WaveList.add(wave_from_doodson( 0, 4, -4,  0,  0,  0, -0.00008f , "libfes92"));
    WaveList.add(wave_from_doodson( 0, 4, -4,  2,  0,  0, -0.00006f , "libfes93"));
    WaveList.add(wave_from_doodson( 0, 4, -4,  2,  1,  0, -0.00002f , "libfes94"));
    WaveList.add(wave_from_doodson( 0, 4, -3,  0,  0,  1, -0.00014f , "libfes95"));
    WaveList.add(wave_from_doodson( 0, 4, -3,  0,  1,  1, -0.00006f , "libfes96"));
    WaveList.add(wave_from_doodson( 0, 4, -2, -2,  0,  0, -0.00011f , "libfes97"));
    WaveList.add(wave_from_doodson( 0, 4, -2,  0,  0,  0, -0.00205f , "libfes98"));
    WaveList.add(wave_from_doodson( 0, 4, -2,  0,  1,  0, -0.00085f , "libfes99"));
    WaveList.add(wave_from_doodson( 0, 4, -2,  0,  2,  0, -0.00008f , "libfes100"));
    WaveList.add(wave_from_doodson( 0, 4, -1, -2,  0,  1, -0.00003f , "libfes101"));
    WaveList.add(wave_from_doodson( 0, 4, -1,  0,  0, -1,  0.00003f , "libfes102"));
    WaveList.add(wave_from_doodson( 0, 4,  0, -2,  0,  0, -0.00169f , "libfes103"));
    WaveList.add(wave_from_doodson( 0, 4,  0, -2,  1,  0, -0.00070f , "libfes104"));
    WaveList.add(wave_from_doodson( 0, 4,  0, -2,  2,  0, -0.00006f , "libfes105"));
    
    /* see lpe_minus_n_waves() in src/prediction.c in libfes */
    WaveList.add(wave_from_doodson( 0,0, 0, 1, 0, 0, -0.00021f , "libfes301"));
    WaveList.add(wave_from_doodson( 0,0, 2, -1, 0, 0, -0.00004f , "libfes302"));
    WaveList.add(wave_from_doodson( 0,1, -2, 0, 0, 0, 0.00004f , "libfes303"));
    WaveList.add(wave_from_doodson( 0,1, 0, 0, -1, 0, 0.00019f , "libfes304"));
    WaveList.add(wave_from_doodson( 0,1, 0, 0, 0, 0, -0.00375f , "libfes305"));
    WaveList.add(wave_from_doodson( 0,1, 0, 0, 1, 0, -0.00059f , "libfes306"));
    WaveList.add(wave_from_doodson( 0,1, 0, 0, 2, 0, 0.00005f , "libfes307"));
    WaveList.add(wave_from_doodson( 0,2, -2, 1, 0, 0, -0.00012f , "libfes308"));
    WaveList.add(wave_from_doodson( 0,2, 0, -1, 0, 0, -0.00061f , "libfes309"));
    WaveList.add(wave_from_doodson( 0,2, 0, -1, 1, 0, -0.00010f , "libfes310"));
    WaveList.add(wave_from_doodson( 0,3, -2, 0, 0, 0, -0.00010f , "libfes311"));
    WaveList.add(wave_from_doodson( 0,3, 0, -2, 0, 0, -0.00007f , "libfes312"));
    WaveList.add(wave_from_doodson( 0,3, 0, 0, 0, 0, -0.00030f , "libfes313"));
    WaveList.add(wave_from_doodson( 0,3, 0, 0, 1, 0, -0.00019f , "libfes314"));
    WaveList.add(wave_from_doodson( 0,3, 0, 0, 2, 0, -0.00004f , "libfes315"));
    WaveList.add(wave_from_doodson( 0,4, 0, -1, 0, 0, -0.00008f , "libfes316"));
    WaveList.add(wave_from_doodson( 0,4, 0, -1, 1, 0, -0.00005f , "libfes317"));
    }
  
  if(isnormal(alias)){
    
    if(showCombinations){
      printf("*** --combinations and --alias are not compatible ***\n");
      print_help(argv[0]);
      wexit(-1);
      }
    
    printf("# ALIASING: %g days = ",alias);
    alias=360./24./alias;/* convert days -> deg/h */
    printf(" = %g degrees/h\n",alias);
    
    for (i=0; i<WaveList.n; i++) {
      wp=&WaveList.waves[i];
      alias_frequency(alias,&wp->omega);
      }
    
    aliased="ALIASED ";
    unaliased="UNALIASED ";
    }
  
  if(nonde<=0){
    qsort(WaveList.waves,WaveList.n,sizeof(*WaveList.waves),compWaveOmega);
    }

  printf("# nbre d'ondes: %d\n",WaveList.n);
/*------------------------------------------------------------------------------
  print table */
  
  printf("--------astronomic potential amplitude in cm, %spulsation in degrees/h",aliased);
  if(not showCombinations)
    printf(" and in /h, %speriod in hours and in days, %scritical latitude in degrees, Doodson number",aliased,unaliased);
  printf("\n");
  
  printf("%10s %7s %15s","wave ","Ap","deg/h ","/h ","h ","days ","deg ","Doodson");
  if(showCombinations)
    printf(" non-linear combinations");
  else
    printf(" %15s %15s %15s %10s %8s","/h ","h ","days ","deg ","Doodson");
  printf("\n");
  
  //printf("----------|---------------|---------------|---------------\n");
  for (i=0; i<WaveList.n; i++) {
    wp=&WaveList.waves[i];
    if(showCombinations and not is_astronomic_or_radiational(*wp))
        continue;
    printf("%10s %7.4f %15.9f", wp->name,wp->Ap*1e2, wp->omega);
    if(showCombinations){
      for(j=0;j<i;j++){
        wpj=&WaveList.waves[j];
        if(not is_astronomic_or_radiational(*wpj))
          continue;
        for(k=0;k<j;k++){
          wpk=&WaveList.waves[k];
          if(not is_astronomic_or_radiational(*wpk))
            continue;
          if( fabs(wpj->omega+wpk->omega-wp->omega)>0.02 )
            continue;
          printf(" %s+%s",wpk->name,wpj->name);
          }
        }
      }
    else{
      period=tide_period(*wp,false);
      critical= asin(wp->omega / two_Omega_deg_h)*r2d;
      asprintf(&doodson,"%07d",doodson_from_wave(*wp));
      printf(" %15.9f %15.9f %15.9f %10.4f %8s", wp->omega/360., period,period/24, critical,doodson);
      free(doodson);
      }
    printf("\n");
    }

  if(nonde>0){
    double duration=cnes_time(final,'d')-cnes_time(start,'d');
    if(isnan(duration)){
      duration=2.;
      }
    printf("----------------------separations above %g days\n",duration/2);
/*---------------------------------------------------------------------*//**<h1>
    For every pair of waves </h1>*/
    for (i=0; i<WaveList.n; i++) {
      for (j=i+1; j<WaveList.n; j++) {
        tau=deltaOmega2separation(WaveList.waves[j].omega-WaveList.waves[i].omega);
        ///unless the separation is well below the duration,
        if(2*tau<duration)//will be false if isnan(duration)
          continue;
        ///it shows the separation
        printf ("wave: %10s %10s, separation: %9.3f days \n",
          WaveList.waves[i].name,WaveList.waves[j].name,tau);
        }
      }
    }

 /*-----------------------------------------------------------------------*/

  if(showArguments && isad(start) && isad(final)){
    printf("\n----------------------arguments (deg, deg, dimensionless) \n");
    for (year=start.year;year<=final.year;year++) {
      reference.year=year;
      for (month=1;month<=12;month++) {
        if ((year==start.year) && (month<start.month)) continue;
        if ((year==final.year) && (month>final.month)) break ;
        for (day=1;day<=poctime_dpm(month,year);day++) {
          if ((year==start.year) && (month==start.month) && (day<start.day)) continue;
          if ((year==final.year) && (month==final.month) && (day>final.day)) break ;
          reference.month=month;
          reference.day=day;
          for (reference.second=0;reference.second<24.*3600.;reference.second+=24.*3600./increment) {
            N=julian_day(reference)-julian_day(start);
            printf("\n----------------------> %s, %d \n",sgetdate(reference),N);
            init_argument(&astro_angles,reference);
            for (i=0; i<WaveList.n; i++) {
              V=greenwhich_argument(astro_angles,WaveList.waves[i]);
              if(V<0.) V+=pi2;
              u=nodal_phase(astro_angles,WaveList.waves[i]);
              if(u<0.) u+=pi2;
              f=nodal_factor(astro_angles,WaveList.waves[i].formula);
              printf("%s : V=%9.3f u= %9.3f f=%9.3f\n",WaveList.waves[i].name,V*r2d,u*r2d,f);
              }
            }
          }
        }
      }
    }

/*------------------------------------------------------------------------------*/

  printf(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
}
