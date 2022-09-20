/**************************************************************************

 T-UGO tools, 2006-2011

 Unstructured Ocean Grid initiative

 ***************************************************************************/
/** \file

 \author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
 \author  Yves Soufflet      LEGOS, Toulouse, France

 VERSION :
 \date May 29, 2012 : Yves Soufflet : First draft

 <!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
 \brief Read LGD files (ascii) or ADCP binaries and compute the vertical diffusivity
 associated with that profile

 \todo at the moment , it only reads a set of harmonic coefficient, needs to be extended
      to read currents time series

 */
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "version-macros.def" //for VERSION and REVISION
#include <iostream>
#include <string>
#include <fstream>


#include "tools-define.h"
#include "tools-structures.h"
#include "legend.h"
#include "spectrum.h"

/* Table of constant values */

static double c_b2 = 16.6;
static double c_b3 = 0.6666667;
static double c_b4 = 1.5;

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
                        " NAME AND VERSION\n"
                        "  %s version " VERSION " " REVISION "\n", prog_name);
        printf("\n USE");
        printf("\n   %s -f FILE  \n",prog_name);
        printf("\n DESCRIPTION");
        printf(
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
             "\n   This program compute the vertical diffusivity associated with a current profile");
        printf("\n OPTIONS :\n"
             "   -f       specify the file name of the lgd file\n"
             "   -o       output file name\n"
              "   -i       for the number of iterations to done in TKE (INT)\n"
              "   -d       specify the dt (time step) for the explicit time stepping\n"
              "   -t       specify which diagnostic you want to do: \n"
              "           test=1 (default) read a LGD file and use Kepsilon\n"
              "           test=2 read LGD file and use the gaspar scheme\n"
              "           test=5 read real velocity time series (not harmonic coefficients) and use kepsilon\n"
              "           test=0 use a specified/hard coded linear velocity profile \n"
            "   -h       print this help\n");

                        __ERR_BASE_LINE__("exiting\n");
                        exit(-1); /** \endcode */
                }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int trisolverd(int la, double *a, double *b, double *diag, double *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    /* Local variables */
    int k;
    double tmp, tmpdiag;

/* -------------------------------------------------------------- */
/* TRIAGONAL SYSTEM SOLVER */
/* -------------------------------------------------------------- */
/* first Gauss elimination */

/* At surface */

    /* Function Body */
    tmp = 1. / diag[0];
    a[0] *= tmp;
    rhs[0] *= tmp;

/* at mid depth */

    for (k = 1; k < la-1; k++) {
        tmpdiag = diag[k] - a[k - 1] * b[k];
        tmp = 1. / tmpdiag;
        a[k] *= tmp;
        rhs[k] = (rhs[k] - rhs[k - 1] * b[k]) * tmp;
    }

/* at bottom, solution in rhs */

    tmpdiag = diag[la-1] - a[la-2] * b[la-1];
    tmp = 1. / tmpdiag;
    rhs[la-1] = (rhs[la-1] - rhs[la-2] * b[la-1]) * tmp;
/* -------------------------------------------------------------- */
/* second and final Gauss elimination to surface */
/* solution in rhs */

    for (k = la-2; k >= 0; k--) {
        rhs[k] -= rhs[k + 1] * a[k];
    }
    return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int gaspar(int kb,     // number of levels
        double grav,   // gravity constant
        double dti,    // baroclinic time step
        double ft,
        double fb,
        double z0b,
        double *q2z,   // q²
        double *q2lz,  // q²l
        double *uzz,   // x-axis horizontal velocity
        double *vzz,   // y-axis horizontal velocity
        double *rhoz,  // levels' density
        double *zi,    // levels' immersion
        double *dzi,
        double *dzzi,
        double *kmz,   // vertical diffusion coefficient (momentum)
        double *khz,   // vertical diffusion coefficient (tracers)
        double *kqz,   // vertical diffusion coefficient (TKE)
        double q2min,
        double lenqmin,
        double kmmin)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* ----------------------------------------------------------------------- */
/* Gaspar's turbulence closer scheme */

/* ----------------------------------------------------------------------- */
/* MODELE SYMPHONIE : VERSION 2003 */
/* DERNIERE REVISION: 1 NOVEMBRE 2001 */
/*      Equipe d'Oc�nographie Coti�e du Laboratoire d'A�ologie */
/*                Laboratoire d'A�ologie */
/*    CNRS - Universit�Paul Sabatier - Observatoire Midi Pyr��s */
/*         14 Avenue Edouard Belin - 31400 Toulouse - FRANCE */
/* Contact:     Patrick.Marsaleix@aero.obs-mip.fr */
/*               Claude.Estournel@aero.obs-mip.fr */
/*                Francis.Auclair@aero.obs-mip.fr */
/* ----------------------------------------------------------------------- */
{
    /**----------------------------------------------------------------------------
  Initialized constants */
    static double ctke1 =  0.1;
    static double ctke2 =  0.7;
    static double b1    = 16.6;

    /**----------------------------------------------------------------------------
  Local variables */
    int kbm1;
    double ddzk, tmp;
    double lenq, suma, dudz, dvdz, sumb, buoy;
    int k;
    double rbase, shear, ldown;
    int itest, k2;
    double const1, dz, drhodz, qoverb1, rap, lup, fkh;
    double *a, *c, *diag, *rhsq2, *sq2;

    /*----------------------------------------------------------------------
  allocate memory for solving the triagonal matrix problem */
    a     = new double[kb];
    c     = new double[kb];
    diag  = new double[kb];
    rhsq2 = new double[kb];
    sq2   = new double[kb];

    /*----------------------------------------------------------------------
  initialize the right hand side */
    for (k = 0; k < kb; k++) rhsq2[k]=q2z[k];

    /*----------------------------------------------------------------------
  define new linit */
    const1 = 1. / sqrt(ctke1 * ctke2);
    kbm1 = kb - 1;

    /**----------------------------------------------------------------------
  loop on mesh levels*/
    for (k = 1; k < kbm1; k++) {
        itest = 0;
        lup = 0.;
        sumb = 0.;
        rap = (zi[k - 1] - zi[k]) * .5 / dzzi[k - 1];
        rbase = (1. - rap) * rhoz[k - 1] + rap * rhoz[k];
        /**----------------------------------------------------------------------
    CALCUL DE L HAUT: LUP
-----------------------------------------------------------------------*/
        /*  calcul de l'integrale: */
        for (k2 = k - 1; k2 >= 0; k2--) {
            /*----------------------------------------------------------------------
      par la suite l'integrale s'effectue sur des intervales entiers. */
            suma = sumb;
            sumb -= grav * (rhoz[k2] - rbase) * (zi[k2] - zi[k2 + 1]);
            if (sumb >= q2z[k]) {
                /*----------------------------------------------------------------------
        test: si l'integrale devient plus grande que l'energie,
        alors la longueur est comprise entre les niveaux K2 et K2-1.
        on sort de la boucle. L'avant derniere valeur de l'integrale
        est archivee dans SUMA.
----------------------------------------------------------------------*/
                itest = 1;
                goto L55;
            }
        }
        /*----------------------------------------------------------------------
   rdv itest.eq.1: */
        L55:
        if (itest == 1) {
            /*----------------------------------------------------------------------
      la longueur de mélange est entre zi(k2)/zi(k2+1) et zi(k)
      lup = zi(k2+1) - zi(k) + dz
      tke = suma - (rho(k2)-rbase) g dz => dz = -(tke - suma) /((rho(k2)-rbase)g)
----------------------------------------------------------------------*/
            dz = -(q2z[k] - suma) / ((rhoz[k2] - rbase) * grav);
            lup = zi[k2 + 1] - zi[k] + dz;
        }
        else {
            lup = zi[1] - zi[k] + z0b;
        }
        /*----------------------------------------------------------------------
    CALCUL DE L BAS  LDOWN */
        ldown = 0.;
        itest = 0;
        sumb = 0.;
        for (k2 = k; k2 < kbm1; k2++) {
            suma = sumb;
            sumb += grav * (rhoz[k2] - rbase) * (zi[k2] - zi[k2 + 1]);
            if (sumb >= q2z[k]) {
                itest = 1;
                goto L77;
            }
        }
        /*----------------------------------------------------------------------
     rdv itest.eq.1: */
        L77:
        if (itest == 1) {
            /*----------------------------------------------------------------------
      la longueur de mélange est entre zi(k2)/zi(k2+1) et zi(k)
      ldown = zi(k) - zi(k2) + dz (dz>0)
      tke = suma + (rho(k2)-rbase) g dz => dz = (tke - suma) /((rho(k2)-rbase)g)
----------------------------------------------------------------------*/
            dz = (q2z[k] - suma) / ((rhoz[k2] - rbase) * grav);
            ldown = zi[k] - zi[k2] + dz;
        }
        else {
            ldown = zi[k] - zi[kb-1];
        }
        lenq = MIN(ldown,lup);
        ddzk = abs(ldown * lup);
        tmp  = sqrt(abs(ddzk));
        /*----------------------------------------------------------------------
    ici q2lz=lq */
        q2lz[k] = MAX(lenqmin,tmp);

        /**----------------------------------------------------------------------------
    diffusion coeficients */
        kmz[k] = ctke1 * lenq * sqrt(q2z[k]);
        //  kmz(K)=0.5d0*(kmz(k)+CTKE1*lenq*SQRT(q2z(K)))
        kqz[k] = kmz[k];
        khz[k] = kmz[k];
        /*----------------------------------------------------------------------
    production terms due to shear (gm) and buoyancy (gh)
                   note gm not used in galperin
----------------------------------------------------------------------*/
        dudz = (uzz[k] - uzz[k - 1]) / dzzi[k - 1];
        dvdz = (vzz[k] - vzz[k - 1]) / dzzi[k - 1];
        shear  = dudz * dudz + dvdz * dvdz;
        drhodz = (rhoz[k - 1] - rhoz[k]) / dzzi[k - 1];
        buoy   = grav * drhodz;

        /*----------------------------------------------------------------------
    finally, compute node right-hand side terms */
        rhsq2[k] += dti * (kmz[k] * shear + khz[k] * buoy);

        /*----------------------------------------------------------------------
    compute node=based dissipation/decay terms
    as quasilinear first-order decay for q2, q2l
----------------------------------------------------------------------*/
        qoverb1 = sqrt(q2z[k]) / q2lz[k];
        sq2[k] = dti * ctke2 * qoverb1;
    }

    lenq = lenqmin;

    /**----------------------------------------------------------------------------
  vertical boundary conditions */
    k = 0;
    rhsq2[k] = fb * const1;
    qoverb1  = sqrt(q2z[k]) / b1;
    sq2[k]   = dti * 2. * qoverb1 / lenq;
    kmz[k]   = kmz[k + 1];
    q2lz[k]  = z0b;

    k = kb-1;
    rhsq2[k] = ft * const1;
    qoverb1  = sqrt(q2z[k]) / b1;
    sq2[k]   = dti * 2. * qoverb1 / lenq;
    kmz[k]   = kmz[k - 1];
    q2lz[k]  = z0b;

    rbase = (1. - rap) * rhoz[k - 1] + rap * rhoz[k];

    /**----------------------------------------------------------------------------
  construct the triagonal matrix and solve for q2 */
    a[0] = 0.;
    c[0] = 0.;
    a[kb-1]=0.;
    c[kb-1]=0.;

    for (k = 1; k < kbm1; k++) {
        fkh  =  0.5 * (kqz[k+1]+kqz[k]+2.0*kmmin);
        a[k] = -dti * fkh / (dzzi[k-1]*dzi[k]);
        fkh  =  0.5 * (kqz[k-1]+kqz[k]+2.0*kmmin);
        c[k] = -dti * fkh / (dzzi[k-1]*dzi[k-1]);
    }

    /**----------------------------------------------------------------------------
  matrix diagonal */
    diag[0] = 1.;
    for (k = 1; k < kbm1; k++) diag[k] = 1. - a[k] - c[k] + sq2[k];
    diag[kb-1]=1.;

    /**----------------------------------------------------------------------------
  solve turbulent kinetic energy system */
    trisolverd(kb,a,c,diag,rhsq2);
    for (k = 0; k < kb; k++) q2z[k]=MAX(q2min,rhsq2[k]);

    delete[] a;
    delete[] c;
    delete[] diag;
    delete[] rhsq2;
    delete[] sq2;

    return 0;

} /* gaspar_ */

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/* ---------------------------------------------------------------------- */
/* k-e turbulence closure scheme */
/* following Burchard et al. (1998), JGR 103, p. 10543-10554 */
/* and Warner et al. (2005) Ocean Modelling 8, 81-113 */
/* ---------------------------------------------------------------------- */

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */



int kepsilon(int kb, //number of levels

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double grav,     //gravity constant
        double dti,      //baroclinic time step
        double ft,       //top boundary condition
        double fb,       //bottom boundary condition
        double z0b,      //roughness length
        double *kz,      //Kinetic energy
        double *ez,      //dissipation rate of kinetic energy
        double *uzz,     //x-axis horizontal velocity
        double *vzz,     //y-axis horizontal velocity
        double *rhoz,    //level density
        double *dzi,     //mid-level depth
        double *dzzi,    //layer thickness
        double *kmz,     //vertical diffusion coefficient (momentum)
        double *khz,     //vertical diffusion coefficient (tracer)
        double *kqz,     //vertical diffusion coefficient (TKE)
        double *kez,     //vertical diffusion coefficient (dissipation rate)
        double kmin,     //minimum Kz
        double emin,     //minimum epsilon
        double kmmin         //minimun Kv
    )
{
    /* Initialized data */

  static double b1 = 16.6;
  static double sk = 1.;
  static double se = 1.3;
  static double e1 = 1.44;
  static double e2 = 1.92;
  static double e3p = 1.;
  static double e3m = -.52;
  static double c0mu = .5544;
  static double g5 = 6.1272;
  static double g6 = .49393;
  static double g8 = 1.83816;
  static double g9 = 15.2352;
  static double g11 = 30.192;

/**----------------------------------------------------------------------------
  System generated locals */
  int kbm1;
  double ddzk, tmp, c0mu3;

/**----------------------------------------------------------------------------
  Local variables */
  static double lenq, dudz, dvdz, buoy, b;
  static int k;
  static double p, shear, lenqq, e3, lenqe1, const1, gh, gm, sh, sm,
                lenq2q2, drhodz, bvf, ftr, fkh;
  double *a, *c, *diag, *rhsk, *rhse, *sq, *sd;


/**----------------------------------------------------------------------------
  allocate memory for solving the tri-diagonal matrix problem */
  a = new double[kb];
  c = new double[kb];
  diag  = new double[kb];
  rhsk = new double[kb];
  rhse = new double[kb];
  sq = new double[kb];
  sd = new double[kb];

/**----------------------------------------------------------------------------
  initialize the right hand side */

  for (k = 0; k < kb; k++) rhsk[k]=kz[k];
  for (k = 0; k < kb; k++) rhse[k]=ez[k];

/**----------------------------------------------------------------------------
  define new linit */

  kbm1 = kb - 1;

/**----------------------------------------------------------------------------
    galperin constants (DRL notes, 24 sept. 1993) */

/*         g0=(1.-6.*a1/b1) */
/*         g1=6.*a1+b2 */
/*         g2=a1*(g0-3.*c1) */
/*         g3=3.*a1*a2*( (b2-3.*a2)*g0 - 3.*c1*g1 ) */
/*         g4=3.*a2*g1 */
/*         g5=9.*a1*a2 */
/*         g6=a2*g0 */
/*         g7=g4 */
/*         g8=g5*(1.-c2) */
/*         g9=18.*a1*a1 */
/*         g10=6.*a1+b2*(1.-c3) */
/*         g11=3.*a2*g10 */

  const1 = exp(log(c_b2)*c_b3);
  c0mu3= c0mu * c0mu * c0mu;

  for (k = 1; k < kbm1; k++) {
/**----------------------------------------------------------------------------
    test if stratification unstable */
    drhodz = (rhoz[k - 1] - rhoz[k]) / dzzi[k - 1];
    if (drhodz > 0.) {
      e3 = e3p;
      }
    else {
      e3 = e3m;
      bvf = sqrt(-grav * drhodz) + 1e-8;
      ftr = kz[k] * 1.3363 * bvf * c0mu3;
      ez[k] = MAX(ez[k],ftr);
      }
/**----------------------------------------------------------------------------
    compute nodal values of mixing length */
    lenq = c0mu3 * kz[k] * sqrt(kz[k]) / ez[k];
    lenqq = lenq * sqrt(2*kz[k]); //see Marsaleix notes

/**----------------------------------------------------------------------------
    production terms due to shear (gm) and buoyancy (gh); note gm not used in galperin */
    dudz = (uzz[k] - uzz[k - 1]) / dzzi[k - 1];
    dvdz = (vzz[k] - vzz[k - 1]) / dzzi[k - 1];
    shear = dudz * dudz + dvdz * dvdz;
    buoy = grav * drhodz;

/**----------------------------------------------------------------------------
   finally, compute node right-hand side terms */
    p = kmz[k] * shear;
    b = khz[k] * buoy;
    rhsk[k] += dti * (p + b);
    sq[k] = dti * ez[k] / kz[k];
    lenqe1 = ez[k] / kz[k];
    rhse[k] += dti * lenqe1 * (e1 * p + e3 * b);
    sd[k] = dti * lenqe1 * e2;

/**----------------------------------------------------------------------------
    compute stability functions */
    lenq2q2 = lenq * lenq / kz[k] * .5;
    gm = lenq2q2 * shear;

/**----------------------------------------------------------------------------
    upper limit on gh (galperin '88, eqn. 29) */
    gh = lenq2q2 * buoy;
/*         gh = (gh-(gh-0.02d0)**2)/(gh+0.028d0-0.04d0) */
    tmp = MIN(gh,.028);
    gh = MAX(-.28,tmp);
/*         gh=min(gh,.028d0) */

/**----------------------------------------------------------------------------
   next compute sm and sh (galperin et al 1988) */
/*         sh=g6/(1.d0-g7*gh) */
/*         sm=(g2-g3*gh)/((1.d0-g4*gh)*(1.d0-g5*gh)) */
/**----------------------------------------------------------------------------
    Kantha and Clayson 1994 eq 29 */
    sh = g6 / (1. - g11 * gh);
    sm = (const1 / b1 + (g9 + g8) * sh * gh) / (1. - g5 * gh);

/**----------------------------------------------------------------------------
    now compute the mixing coefficients and impose minimum values */
    kmz[k] = lenqq * sm;
    khz[k] = lenqq * sh;
    kqz[k] = kmz[k] / sk;
    kez[k] = kmz[k] / se;
    }
/**----------------------------------------------------------------------------
  no decay l.h.s. terms */
/*         sq=0.d0 */
/*         sd=0.d0 */
  k = 0;
  rhsk[k] = fb/sqrt(sqrt(2)*c0mu3*sm);
  rhse[k] = 0.;
  sq[k] = 0.;
  sd[k] = 0.;
  kqz[k] = kqz[k + 1];
  kez[k] = kez[k + 1];
  kmz[k] = kmz[k + 1];
  khz[k] = khz[k + 1];
  k = kb-1;
  rhsk[k] = ft;
  rhse[k] = 0.;
  sq[k] = 0.;
  sd[k] = 0.;
  kqz[k] = kqz[k - 1];
  kez[k] = kez[k - 1];
  kmz[k] = kmz[k - 1];
  khz[k] = khz[k - 1];
  
/**----------------------------------------------------------------------------
  construct the triagonal matrix and solve for k */

  a[0] = 0.;
  c[0] = 0.;
  a[kb-1]=0.;
  c[kb-1]=0.;

  for (k = 1; k < kbm1; k++) {
    fkh  = 0.5 * (kqz[k+1]+kqz[k]+2.0*kmmin);
    a[k] = -dti * fkh / (dzzi[k-1]*dzi[k]);
    fkh  = 0.5 * (kqz[k-1]+kqz[k]+2.0*kmmin);
    c[k] = -dti * fkh / (dzzi[k-1]*dzi[k-1]);
    }

  diag[0] = 1.;
  for (k = 1; k < kbm1; k++) diag[k] = 1. - a[k] - c[k] + sq[k];
  diag[kb-1]=1.;

  trisolverd(kb,a,c,diag,rhsk);
  for (k = 0; k < kb; k++) {
    kz[k]=MAX(kmin,rhsk[k]);
    }

/**----------------------------------------------------------------------------
  construct the triagonal matrix and solve for epsilon  */

  a[0] = 0.;
  c[0] = 0.;
  a[kb-1]=0.;
  c[kb-1]=0.;

// get the BC conditions from kz computed above
  rhse[0]=pow(kz[1],1.5)*c0mu3/(0.4*(dzzi[0]+z0b));
//rhse[0]=kz[1]*kz[1]*sqrt(2)*c0mu3*sm/kmz[1];
//rhse[kb-1]=kz[kb-2]*kz[kb-2]*sqrt(2)*c0mu3*sm/kmz[kb-2];
  rhse[kb-1]=pow(kz[kb-2],1.5)*c0mu3/(0.4*(kb*dzzi[kb-2]+z0b));
//rhse[kb-1]=kz[kb-2]*c0mu3/(0.4*dzzi[kb-2]);

  for (k = 1; k < kbm1; k++) {
    fkh  = 0.5 * (kez[k+1]+kez[k]+2.0*kmmin);
    a[k] = -dti * fkh / (dzzi[k-1]*dzi[k]);
    fkh  = 0.5 * (kez[k-1]+kez[k]+2.0*kmmin);
    c[k] = -dti * fkh / (dzzi[k-1]*dzi[k-1]);
    }

  diag[0] = 1.;
  for (k = 1; k < kbm1; k++) diag[k] = 1. - a[k] - c[k] + sd[k];
  diag[kb-1]=1.;

  trisolverd(kb,a,c,diag,rhse);
  
  for (k = 0; k < kb; k++){
    ez[k]=MAX(emin,rhse[k]);
    kmz[k]=MAX(kmmin,kmz[k]);
    }

  delete[] a;
  delete[] c;
  delete[] diag;
  delete[] rhsk;
  delete[] rhse;
  delete[] sq;
  delete[] sd;
  return 0;
} /* kepsilon_ */

void lecture3(char * filename, double ***U ,double ***V ,double **time,double **Z, int& nsample,int& nlevel){
    char *filename1,*filename2;
    char *ext = "U.plt";
    int len   = strlen(filename) + strlen(ext) + 1;
    filename1=new char[len];

    strcpy(filename1,filename);
    strcat(filename1,ext);
    char *ext2 = "V.plt";
    filename2=new char[len];
    strcpy(filename2,filename);
    strcat(filename2,ext2);
    ifstream fichier1(filename1);
    ifstream fichier2(filename2);
    double lon,lat,datestart,dateend;
    double MeanPressureDepth,firstdepth,binsize;
    double *depth,*dates,**profilesU,**profilesV;
    int i,j;
    string line,l,keystring,stationname;


    if(fichier1) {
        //L'ouverture s'est bien passée, on peut donc lire

        //Une variable pour stocker les lignes lues

        getline(fichier1, line); //get first line
        stringstream iss( line );
        iss>>keystring>>stationname>>keystring>>lon>>keystring>>lat>>keystring>>datestart>>keystring>>dateend>>keystring>>nsample;
        printf(" coucou %d ",nsample);
        iss.flush ();
        iss.clear ();
        iss.str ("");
        getline(fichier1, line); //get second line
        iss<<line;
        iss>>keystring>>MeanPressureDepth>>keystring>>firstdepth>>keystring>>binsize>>keystring>>nlevel;

        printf(" %lf",MeanPressureDepth);

        printf(" coucou %d",nlevel);
        depth = new double[nlevel];
        dates = new double[nsample];
        profilesU = new double *[nsample];

        //iss.flush();
        getline(fichier1, line); //get third line
        iss.clear ();
        iss.str ("");
        iss<<line;
        for (i =0; i<nlevel;i++) iss>>depth[i];
        i=0;
        while(getline(fichier1, line)) //Tant qu'on n'est pas à la fin, on lit
        {     //iss.flush();
            iss.clear ();
            iss.str ("");
            profilesU[i] = new double [nlevel];
            iss<<line;
            iss>>dates[i];
            for (j =0; j<nlevel;j++) iss>>profilesU[i][j];
            i++;

        }
        printf(" et voilou %lf et jour final %lf \n",dates[nsample -1],dateend );
        printf(" et le dernier profile:\n");
        for (j =0; j<nlevel;j++) printf(" %lf",profilesU[nsample-1][j]);
        printf(" \n");
    }
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    if (fichier2) {
        getline(fichier2, line);
        stringstream iss( line );
        getline(fichier2, line);
        iss.clear ();
        iss.str ("");
        iss<<line;
        getline(fichier2, line);
        profilesV = new double *[nsample];
        double tmp;

        i=0;
        while(getline(fichier2, line)) //Tant qu'on n'est pas à la fin, on lit
        {     //iss.flush();
            iss.clear ();
            iss.str ("");
            profilesV[i] = new double [nlevel];
            iss<<line;
            iss>>tmp;
            for (j =0; j<nlevel;j++) iss>>profilesV[i][j];
            i++;

        }
        printf(" et le dernier profile:\n");
        for (j =0; j<nlevel;j++) printf(" %lf",profilesV[nsample-1][j]);
        printf(" \n");
    }
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    *Z=depth;
    *time=dates;
    *U=profilesU;
    *V=profilesV;
    //return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    int status, *rstatus,i,j,k;
    char *filename=NULL, *keyword, *output=NULL,*maxiteration=NULL,*dtime=NULL,*testnumber=NULL;
    int n = 1,nframe,maxiter;
    double *time,dt2;
    spectrum_t prediction;
    i=j=k=0;
    FILE   *fp=NULL;
    legend03_t *legend,*outlegend;
    legend_t *OutLegend;
    output=new char[500];
    maxiter=30;
    dt2=30.0;
    int test=1;
  
  fct_echo(argc,argv);
 
     if(argc==1) {
    __OUT_BASE_LINE__("turbulence -h (for help)\n");
    exit(-1);
    }
   
    while (n < argc) {
                keyword = strdup(argv[n]);
                switch (keyword[0]) {
                case '-':
                        switch (keyword[1]) {

                        case 'f':
                                filename = strdup(argv[n + 1]);
                                n++;
                                n++;
                                printf("FILENAME %s\n", filename);
                                break;
                        case 'o':
                                output = strdup(argv[n + 1]);
                                n++;
                                n++;
                                break;
                        case 'i':
                                maxiteration = strdup(argv[n + 1]);
                                maxiter=(int) NINT(strtod(maxiteration,NULL));
                                n++;
                                n++;
                                break;
                        case 'd':
                                dtime = strdup(argv[n + 1]);
                                dt2=(double) strtod(dtime,NULL);
                                n++;
                                n++;
                                break;
                       case 't':
                                testnumber = strdup(argv[n + 1]);
                                test=(int) NINT(strtod(testnumber,NULL));
                                n++;
                                n++;
                                break;
                        case 'h':
                                print_help(argv[0]);
                                break;
                        }
                        break;

                default:
                        printf("unknown option %s\n", keyword);
                        print_help(argv[0]);
                        exit(-1);
                }
                free (keyword);
        }

        printf("on est laaaaa! \n");
//fp = fopen(filename,"r");
/*         if ((fp = fopen(filename, "r")) == NULL) {
              printf("\n Cannot read the file: %s\n",filename);
              return(EXIT_FAILURE);
            }*/
        //read LGD ascii File
//        fscanf(fp,"%d",i);
//        printf("on est ici %d \n",i);
//        std::ifstream input(filename);
//        std::string line;
//        std::getline( input, line );
//        int nitems=sscanf(line.c_str(), "%d %d %d",&i,&j,&k);
//        printf("on est ici %d %d %d %d \n",nitems,i,j,k);
    int np,ndepths,ncolumns,nsample;
    const double ezmin=1e-12;
    const double kzmin=7.6e-6;
    const double kmmin=1e-4;
    double *meanez,*meankz,*meankmz;
    double Cd=0.001;
    double ft=0.0;
    double Z0=0.001;
    double ustar,fb;
    double *z,*zz,*dz,*dzz,*rhoz,*kmz,*khz,*kqz,*kz,*kez,*ez,*q2lz;
    char buf[1000];
    double **profileU=NULL,**profileV=NULL,*dates=NULL;
    /**-----------------------------------------------------------------------
      minimum turbulent macro-scale*/
    double lenqmin=1.;
    double *u,*v;

    if (test!=5){
        lgd_load(filename, rstatus);
        legend=(legend03_t*) legendstab[0].ptr;
        printf("number of records: %d , levels: %d , and data columns: %d \n",legend->np,legend->ndepths,legend->ncolumns);
        printf("position of profile: %f %f \n",legend->points[0].t,legend->points[0].p);
        printf("first lines of data: (depth , amp phase...) %f %f %f %f %f \n",legend->points[0].depths[0],legend->points[0].z[0][0],legend->points[0].z[0][1],legend->points[0].z[0][2],legend->points[0].z[0][3]);
        printf("number of iterations per frame %d  \n",maxiter);
        printf("corresponding to a delta t: %f  \n",dt2);

        np=legend->np;
        ndepths=legend->ndepths;
        ncolumns=legend->ncolumns;


        /*-----------------------------------------------------------------------
      immersion of levels at ncP1 nodes horizontal position*/

        if (test==0) ndepths=101;

        z  = new double[ndepths];
        if (test==0){
            for (k=0;k<ndepths;k++) {
                z[k]=-100.0+k;
            }
            //printf("Hola %f \n",z[0]);
        }
        else{
            for (k=0;k<ndepths;k++) {
                z[k]=legend->points[0].depths[k];
            }
        }

    }
    else{    //read real profile

        lecture3(filename,&profileU,&profileV,&dates, &z,nsample,ndepths);
        printf("apres lecture du fichier vitesse %d \n",ndepths);
        
    }

    u =new double[ndepths];
    v= new double[ndepths];
    double H=-z[0];

    /*layer thickness*/
    zz = new double[ndepths -1];
    dz = new double[ndepths-1];
    dzz= new double[ndepths-1];


    /*density profile*/
    rhoz= new double[ndepths];

    /**-----------------------------------------------------------------------
    diffusion coefficient for momentum*/
    kmz= new double[ndepths];

  /**-----------------------------------------------------------------------
    diffusion coefficient for tracer*/
    khz= new double[ndepths];
  /**-----------------------------------------------------------------------
    diffusion coefficient for turbulent energy*/
    kqz= new double[ndepths];
    kez= new double[ndepths];
    ez= new double[ndepths];
    kz= new double[ndepths];
    q2lz=new double[ndepths];

    meankz= new double[ndepths];
    meanez= new double[ndepths];
    meankmz= new double[ndepths];
  /**-----------------------------------------------------------------------
    dissipation rate for turbulent energy*/

  /**----------------------------------------------------------------------
    initialisation*/
    for (i=0; i<ndepths; i++) {
            rhoz[i]= 26.;//*(0.2*z[i]+1.1);
        kz[i]=kzmin;
        meankz[i]=kzmin;
        ez[i]= ezmin;
        meanez[i]= ezmin;
        kqz[i]=kmmin;
        kez[i]=kmmin;
        kmz[i]=kmmin;
        khz[i]=kmmin;
        meankmz[i]=kmmin;
        q2lz[i]=Z0;
    }

    for (i=0; i<ndepths-1; i++) {
           zz[i]=0.5*(z[i]+z[i+1]);
           dz[i]=z[i+1]-z[i];
           }
    for (i=0; i<ndepths-1; i++) {
        dzz[i]=zz[i+1]-zz[i];
        }
          double uamp,vamp,uphase,vphase;
          double cs[3],sn[3];
          double mkz=1.e10,mez=1.e10;


    /**----------------------------------------------------------------------
                initialisation of the prediction parameters */
    double dt=3600.0;
    double duration=24*3600;
    nframe=(int) NINT(duration/dt);
    if (test==5) {
        nframe=nsample;
        int secondes=int((dates[1]-dates[0])*24*3600);
        dt2=10.0;
        maxiter=secondes/10;
        printf("el dt esta= %d",secondes);
    }
    
    harmonic_t h2;
    date_t start;
    if ((test>0)&&(test!=5)){
        printf("Number of levels %i \n",ndepths);
    time = new double[nframe];
    for(k=0;k<nframe;k++) time[k]=(double)k*dt;

    start.year=1990;
    start.month=1;
    start.day=1;
    start.second=0;
    prediction.waves=new tidal_wave[100];
    prediction.nmax=100;
    bool Z0=false;
    prediction = spectrum_init_ref("COASTAL", Z0);
    status=harmonic_coefficients(nframe, start, time, prediction, &h2,0);
    }
    int nsamples;

    /**----------------------------------------------------------------------
      Loop over each frame of the prediction */
    for(int m=0;m<nframe;m++) {
    nsamples=0;
    /**----------------------------------------------------------------------
    reinitialisation of currents and turbulence variables
    as we compute a value of Kv for each "frame" , i.e.: each value of (u,v)*/
    for (k=0;k<ndepths;k++){
        u[k]=0.0;
        v[k]=0.0;
        kmz[k]=kqz[k]=kez[k]=khz[k]=kmmin;
        meankmz[k]=kmmin;
        kz[k]=meankz[k]=kzmin;
        ez[k]=meanez[k]=ezmin;
        }
    /**----------------------------------------------------------------------
          computing (u,v) for this frame from the harmonic coefficient*/
        if (test==0){
            for (k=0;k<ndepths;k++){
                u[k]=(z[k]+H)/H;
                v[k]=0.0;
            }
            //printf("Que t'al2 %lf \n",u[0]);
        }
        if (test==5){
            for (k=0;k<ndepths;k++){
                u[k]=profileU[m][k];
                v[k]=profileV[m][k];
            }
        }
        else {
          for(int w = 0; w < prediction.n; w++) {
          if (strcmp(prediction.waves[w].name,"M2")==0){
              for (k=0;k<ndepths-1;k++) {
               uamp=legend->points[0].z[k][0];
               uphase=legend->points[0].z[k][1];
               vamp=legend->points[0].z[k][2];
               vphase=legend->points[0].z[k][3];
               //printf("Read Data: %d %lf %lf %lf %lf \n",k,uamp,uphase,vamp,vphase);
          cs[0]=  uamp*cos(uphase*M_PI/180.) ;
          cs[1]=  vamp*cos(vphase*M_PI/180.) ;
          sn[0]= uamp*sin(-uphase*M_PI/180.) ;
          sn[1]= vamp*sin(-vphase*M_PI/180.) ;
          u[k] += cs[0] * h2.cs[w][m] + sn[0] * h2.sn[w][m];
          v[k] += cs[1] * h2.cs[w][m] + sn[1] * h2.sn[w][m];
          }
      u[ndepths-1]=u[ndepths-2];
      v[ndepths-1]=v[ndepths-2];
      }
      }
      //printf("Que t'al3 %lf \n",u[0]);
      }

      //compute TKE
      ustar=Cd*sqrt(u[0]*u[0]+v[0]*v[0]);
      fb=ustar*ustar/Cd;
      if (test!=5){
      for (int iter=0; iter<maxiter;iter++){
          if (test==2) gaspar(ndepths,9.81,dt2,ft,fb,Z0,kz,q2lz,u,v,rhoz, z,dz,dzz,kmz,khz,kqz,kzmin,lenqmin,kmmin);
          else if (test!=2) kepsilon(ndepths,9.81,dt2,ft,fb,Z0,kz,ez,u,v,rhoz,dz,dzz,kmz,khz,kqz,kez,kzmin,ezmin,kmmin);

          for (k = 0; k < ndepths; k++){
                    kmz[k]=MAX(kmmin,kmz[k]);
                }

          for (k=0; k<ndepths;k++) {
                    meankz[k]= (meankz[k]*nsamples+kz[k])/(nsamples+1);
                    meankmz[k]= (meankmz[k]*nsamples+kmz[k])/(nsamples+1);
                    meanez[k]= (meanez[k]*nsamples+ez[k])/(nsamples+1);
                }
          nsamples++;
            }
        }
        else {
            for (int iter=0; iter<maxiter;iter++){
                kepsilon(ndepths,9.81,dt2,ft,fb,Z0,kz,ez,u,v,rhoz,dz,dzz,kmz,khz,kqz,kez,kzmin,ezmin,kmmin);

            for (k = 0; k < ndepths; k++){
                kmz[k]=MAX(kmmin,kmz[k]);
                }
            if (iter>30){
                for (k=0; k<ndepths;k++) {
                        meankz[k]= (meankz[k]*nsamples+kz[k])/(nsamples+1);
                        meankmz[k]= (meankmz[k]*nsamples+kmz[k])/(nsamples+1);
                        meanez[k]= (meanez[k]*nsamples+ez[k])/(nsamples+1);
                        }
                nsamples++;
                }
            }
        }

        /**----------------------------------------------------------------------
                      initialisation of the legend format output */
        if (test==5){
            lgd_alloc(1,rstatus);
            legendstab[0].Type=20;//LGD_PROFILE
            legendstab[0].ptr=new char[sizeof(legend03_t)];
        }
        OutLegend=legendstab;
        outlegend=new legend03_t;
        outlegend->np=1;
        outlegend->ID=1;
        outlegend->ndepths=ndepths-1;
        outlegend->ncolumns=6;
        outlegend->points=new point3D_t[outlegend->np];
        outlegend->points[0].t=0.0;
        outlegend->points[0].p=0.0;
        outlegend->points[0].depths=new float[ndepths-1];
        outlegend->points[0].z=new float *[ndepths-1];
        for (int n=0; n < outlegend->ndepths;n++) outlegend->points[0].z[n]=new float[outlegend->ncolumns];
        char *fmt=NULL,*comments=NULL;

                for (int n=0; n < outlegend->ndepths;n++) {
                            outlegend->points[0].depths[n]=zz[n];
                            outlegend->points[0].z[n][0]=u[n];
                            outlegend->points[0].z[n][1]=v[n];
                            outlegend->points[0].z[n][2]=kqz[n];
                            outlegend->points[0].z[n][3]=kz[n];
                            outlegend->points[0].z[n][4]=ez[n];
                            outlegend->points[0].z[n][5]=kmz[n];
                          }

                OutLegend[0].ptr= (char *) outlegend;
                //char *outfilename="essai.lgd";
                snprintf(buf, sizeof buf, "%s.%d", output,m);
                lgd_save(buf,OutLegend,1,fmt,comments);
                //output k and epsilon
        free(outlegend);
    }

printf("un quatrieme essai %f %f \n", meankz[0],meankmz[0]);





printf("C'est fini!!!");
}

