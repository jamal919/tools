

/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

***************************************************************************/


#define FREE_SLIP 0
#define NO_SLIP   1
#define DISCHARGE 2


class river_t {
private :
public :
  vector<int> edges;                                     /* edge path                             */
  vector<int> nodes;                                     /* edge path                             */
  vector<int> elements;                                  /* elements                              */
  string DischargeFilename,     ElevationFilename;       /* input data                            */
  string DischargeTimeTemplate, ElevationTimeTemplate;   /* input data                            */
  int    DischargeID,           ElevationID;             /*                                       */
  int    RefreshMode;                                    /* time evolution mode                   */
  int    PrescriptionMode;                               /* model prescripyion mode               */
  double tau;                                            /* relaxation time                       */
  double debit, elevation;                               /* instant mass elevation                */
  double *Fx, *Fy;                                       /* instant mass flux                     */
  double T, S, R;                                        /* temperature, salinity,density         */
  double *nx, *ny;                                       /* to convert m³/s into directional m²/s */
  double *h, length, section;                            /*                                       */
  hconstant_t Fperiodic, Zperiodic;                      /*                                       */
  /*STS can not find sequence_t and does not use DischargeStream, ElevationStream*/
  double ElevationOffset;
  
  void init() {
//    edge=-1;
    RefreshMode=RefreshMode=-1;
    DischargeFilename=(string) "NONE";
    ElevationFilename=(string) "NONE";
    debit=elevation=0;
    Fx=Fy=nx=ny=h=0;
    length=section=0;
    ElevationOffset=0;
    T=S=R=0;
    }
    
  river_t() {
//    edge=-1;
    RefreshMode=RefreshMode=-1;
//     DischargeFilename=(string) "NONE";
//     ElevationFilename=(string) "NONE";
    debit=elevation=0;
    Fx=Fy=nx=ny=h=0;
    length=section=0;
    ElevationOffset=0;
    T=S=R=0;
    }
    
  void destroy() {
    deletep(&Fx);
    deletep(&Fy);
    deletep(&nx);
    deletep(&ny);
    deletep(&h);
    init();
    edges.clear();
    elements.clear();
//    ~edges();
    }
  ~river_t() {
//     deletep(&Fx);
//     deletep(&Fy);
//     deletep(&nx);
//     deletep(&ny);
//     deletep(&h);
//    ~edges();
    }
};

// class  riverspec_t {
// private :
// public :
//   int    edge,orientation; /*element edge, ???                      */
//   int    mode;             /* time evolution mode                   */
//   double debit,elevation;  /* instant mass elevation                */
//   double Fx,Fy;            /* instant mass flux                     */
//   double T,S,R;            /* temperature, salinity,density         */
//   double cx,cy;            /* to convert m³/s into directional m²/s */
//   hconstant_t Fperiodic,Zperiodic;
//   
//   riverspec_t() {
//     edge=-1;
//     mode=-1;
//     }
//   
// } ;

class specification_t {
private :
public :
  int *nodes, nnodes;           /* nodes index list                  */
//  discretisation_t *descriptor;
  int code;                     /* boundary code and type             */
//   tidespec_t *tides;
//   fluxspec_t *fluxes;


  specification_t() {
    nodes=0;
    code=-1;
    }
  void destroy() {
    if(nodes!=0) {
      delete[] nodes;
      nodes=0;
      };
//    deletep(&tides);
//    deletep(&fluxes);
    nnodes=0;
//    descriptor=0;
    code=-1;
    };
};

class dataspec_t {                      /*                                                       */
private :
public :
  int node;                             /*  Node index in the current mesh                       */
  int ancestor;                         /*  Node index in the current mesh                       */
  double nx,ny,L;                       /*  flux direction, boundary size                        */
};

class tidespec_t {                      /*                                                       */
private :
public :
  int node;                             /*  Node index in the current mesh                       */
  hconstant_t ztide,utide,vtide;        /*  Specified tidal constants                            */
  spectrum_t *spectrum;
//   tidespec_t (int n) {
//     
//     };
  void destroy() {
    ztide.destroy();
    utide.destroy();
    vtide.destroy();
    }
};

class statespec_t {                     /*                                                       */
private :
public :
  int node;                             /*  Node index in the current mesh                       */
  float z,u,v;
  float H,Hu,Hv;
  float h;
};

class fluxspec_t {                      /*                                                       */
private :
public :
  int node;                             /*  Node index in the current mesh                       */
  int option;                           /*  Node index in the current mesh                       */
  double value;                         /*  Specified flux value                                 */
};

class specification_obsolete_t {
private :
public :
  int *edges, nedges;   /* edges index list                   */
  int *nodes, nnodes;   /* vertex index list                  */
  int *elmts, nelmts;   /* element index list                 */
  int discretisation;
  int *PrescriptionMode;
  
  dataspec_t  *data;
  tidespec_t  *tides;
  fluxspec_t  *fluxes;
  statespec_t *external;
  
  int *r_edge, *l_edge; /* neighbour index list               */
  int *r_node, *l_node; /* vertex index list                  */
  int *table;           /* vertex index list                  */
  int code;             /* boundary code and type             */
  int hbw;              /* half bandwidth                     */

  specification_obsolete_t() {
    edges=0;
    nodes=0;
    elmts=0;
    nedges=nnodes=nelmts=0;
    discretisation=-1;
    PrescriptionMode=0;
    
    data=0;
    tides=0;
    fluxes=0;
    external=0;
    
    r_edge=l_edge=0;
    r_node=l_node=0;
    table=0;
    code=-1;
    hbw=-1;
    }
    
  void destroy() {
    int n;
    if(edges!=0) {
      delete[] edges;
      edges=0;
      }
    if(nodes!=0) {
      delete[] nodes;
      nodes=0;
      }
    if(elmts!=0) {
      delete[] elmts;
      elmts=0;
      }
    if(PrescriptionMode!=0) {
      delete[] PrescriptionMode;
      PrescriptionMode=0;
      }
      
    deletep(&r_edge);
    deletep(&l_edge);
    deletep(&r_node);
    deletep(&l_node);
    deletep(&table);
    
    if(tides!=0) {
      for(n=0;n<nnodes;n++) {
//        tides[n].destroy(); /// HERE !!!
        }
      }
    
    deletep(&data);
    deletep(&tides);
    deletep(&fluxes);
    deletep(&external);
    
    nedges=nnodes=nelmts=0;
    discretisation=-1;
    }
};

