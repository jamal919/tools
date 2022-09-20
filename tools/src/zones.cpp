
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief zone option parsing functions
*/
/*----------------------------------------------------------------------------*/

#include "map.h"
#include "functions.h" // for scan_angle


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_zonegrid(const char *zone, bool *identified)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  grid_t grid;
  int zone_initialised=0;
/*
To carry out the clean-up, I used this piece of javascript code :
.replace(/(?:([ \t]*)grid\.d?[xy](?:min|max|)[^;]*; *\n){4,6}/g,function(m,i){grid={};eval(m.replace(/= *([^;]+?) *;/g,"='$1';"));try{with(grid){return i+'map_set2Dgrid(&grid,'+[xmin,ymin,xmax,ymax,dx,dy].join(',')+');\n'}}catch(e){return m;}})
*/

  if(strcmp(zone,"FRANCE")==0) {
    map_set2Dgrid(&grid,-6.,+42.5,+9.,+51.5,1./60.,1./60.);
    zone_initialised=1;
    }

/* NEA grid */
  if(strcmp(zone,"NEA")==0) {
    map_set2Dgrid(&grid,-20.,+29.5,+15.,+65.,1./30.,1./30.);
    zone_initialised=1;
    }

/* NEA_hr grid haute resolution */
  else if(strcmp(zone,"NEA-LR")==0) {
    map_set2Dgrid(&grid,-20.,+22.5,+15.,+65.,5./60.,5./60.);
    zone_initialised=1;
    }

/* NEA_hr grid haute resolution */
  else if(strcmp(zone,"NEA-HR")==0) {
    map_set2Dgrid(&grid,-20.,+22.5,+15.,+65.,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"NEA-VHR")==0) {
    map_set2Dgrid(&grid,-20.,+22.5,+15.,+65.,1./120.,1./120.);
    zone_initialised=1;
    }

/* NEA_hr grid projet REF */
  else if(strcmp(zone,"NEA-SURFREF")==0) {
    map_set2Dgrid(&grid,-12.5,+42.5,+3.,+52.5,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"NEA-COMAPI")==0) {
    map_set2Dgrid(&grid,-22.5,+22.5,13.,+65.,1./60.,1./60.);
    zone_initialised=1;
    }

/* Manche_hr grid haute resolution */
  else if(strcmp(zone,"Manche-HR")==0) {
    map_set2Dgrid(&grid,-6.,+48.,+6.,+53.,1./240.,1./240.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"seine-estuary-LR")==0) {
    map_set2Dgrid(&grid,0.,+49.25,1.25,+49.75,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"seine-estuary")==0) {
    map_set2Dgrid(&grid,0.,+49.25,1.25,+49.75,1./3600.,1./3600.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"seine-estuary-HR")==0) {
    map_set2Dgrid(&grid,0.,+49.25,1.25,+49.75,1./3600./4.,1./3600./4.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"North-Sea")==0) {
    map_set2Dgrid(&grid,+5.,+51.,+13.,+60.,2./60.,2./60.);
    zone_initialised=1;
    }

 /* NEA_puertos grid distribution puertos  haute resolution */
  else if(strcmp(zone,"NEA_puertos")==0) {
    map_set2Dgrid(&grid,-18.,+30.,+0.,+50.,2./60.,2./60.);
    zone_initialised=1;
    }

 /* NEA_puertos grid distribution puertos  haute resolution */
  else if(strcmp(zone,"NEA-ESEO")==0) {
    map_set2Dgrid(&grid,-21.,+21.,+14.,+64.,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"NEA-ESEO-2011")==0) {
    map_set2Dgrid(&grid,-35.,+35.,+20.,+49.,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"NEA-ESEO-LR")==0) {
    map_set2Dgrid(&grid,-21.,+23.,-9.,+35.,1./20.,1./20.);
    zone_initialised=1;
    }

 /* NEA_shom grid distribution shom  haute resolution */
  else if(strcmp(zone,"NEA_shom")==0) {
    map_set2Dgrid(&grid,-20.,+40.,+15.,+55.,1./30.,1./30.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"NEA_shom2013")==0) {
    map_set2Dgrid(&grid,-30.,+30.,+20.,+65.,1./30.,1./30.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"norway")==0) {
    map_set2Dgrid(&grid,0.,30.,55.,75.,1./120.,1./120.);
    zone_initialised=1;
    }

 /* NEA_shom grid distribution shom  haute resolution */
  else if(strcmp(zone,"azores")==0) {
    map_set2Dgrid(&grid,-40.,+31.,-20.,+45.,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"gascogne")==0) {
    map_set2Dgrid(&grid,-20.,+40.,+5.,+55.,2./60.,2./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"gascogne-COMAPI")==0) {
    map_set2Dgrid(&grid,-19.,+39.,+2.5,+52.5,2./60.,2./60.);
    zone_initialised=1;
    }

/* Iroise grid */
  else if(strcmp(zone,"iroise")==0) {
    map_set2Dgrid(&grid,-6.,+47.25,-3.5,+49.,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"iroise-VHR")==0) {
    map_set2Dgrid(&grid,-6.,+47.25,-3.5,+49.,1./300.,1./300.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"gironde")==0) {
    /*  */
    map_set2Dgrid(&grid,-2.5,+44.30,0.,+46.30,1./120.,1./120.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"gironde-estuary")==0) {
    /* Gironde estuary, ~10 m resolution */
    map_set2Dgrid(&grid,-1.5,+44.30,0.1,+46.30,1./3600./3.,1./3600./3.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"anglet")==0) {
    /*l */
    map_set2Dgrid(&grid,-2.25,+43.25,-1.25,+44.25,0.025,0.025);
    zone_initialised=1;
    }

/* global grids : keep 0<=lon<=360. for backward consistency */
  else if(strcmp(zone,"global")==0 ||
      strcmp(zone,"global-0.25")==0) {
    map_set2Dgrid(&grid,0.,-90.,+360.,+90.,0.25,0.25);
    zone_initialised=1;
    }

  else if(strcmp(zone,"global-HR")==0) {
    map_set2Dgrid(&grid,0.,-90.,+360.,+90.,1./8.,1./8.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"global-VHR")==0) {
    map_set2Dgrid(&grid,0.,-90.,+360.,+90.,1./16.,1./16.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"global-VHR-180")==0) {
    map_set2Dgrid(&grid,-180.,-90.,+180.,+90.,1./16.,1./16.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"global-0.5")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,0.,-90.,+360.,+90.,0.5,0.5);
    zone_initialised=1;
    }

  else if(strcmp(zone,"global-1.0")==0) {
    map_set2Dgrid(&grid,0.,-90.,+360.,+90.,1.,1.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"global-1mn")==0) {
    map_set2Dgrid(&grid,0.,-90.,+360.,+90.,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"global-30s")==0) {
    map_set2Dgrid(&grid,0.,-90.,+360.,+90.,1./120.,1./120.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"caspian")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,46.,36.,55.,47.5,0.05,0.05);
    zone_initialised=1;
    }

/* madridge grid, normal resolution */
  else if(strcmp(zone,"madridge")==0) {
    map_set2Dgrid(&grid,40.,-36.4,50.0,-25.0,1./120.,1./120.);
    zone_initialised=1;
    }

/* medsea grid, normal resolution */
  else if(strcmp(zone,"medsea")==0) {
    map_set2Dgrid(&grid,-10.,27.5,40.,47.5,1./12.,1./12.);
    zone_initialised=1;
    }

/* medsea grid, high resolution */
  else if(strcmp(zone,"medsea-HR")==0) {
    map_set2Dgrid(&grid,-10.,27.5,40.,47.5,1./60.,1./60.);
    zone_initialised=1;
    }

/* medsea grid, high resolution */
  else if(strcmp(zone,"medsea-SURFREF")==0) {
    map_set2Dgrid(&grid,3.,41.,13.,44.5,1./60.,1./60.);
    zone_initialised=1;
    }

/* medsea grid, high resolution */
  else if(strcmp(zone,"medsea-COMAPI")==0) {
    map_set2Dgrid(&grid,3.,41.,13.,44.5,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"medsea-INGV")==0) {
    map_set2Dgrid(&grid,-20.,29.,40.,46,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"gibraltar")==0) {
    /* Zoom Gibraltar */
    map_set2Dgrid(&grid,-6.5,+35.5,-4.5,+36.5,.01,.01);
    zone_initialised=1;
    }
  else if(strcmp(zone,"gibraltar-strait")==0) {
    /* Zoom Gibraltar */
    map_set2Dgrid(&grid,-10.0,+32.5,0.0,+38.0,.01,.01);
    zone_initialised=1;
    }
/* alobran grid, normal resolution */
  else if(strcmp(zone,"alboran")==0) {
    map_set2Dgrid(&grid,-9.,32,-2+1./300,38+1./300,1./60.,1./60.);
    zone_initialised=1;
    }

/* alobran grid, normal resolution */
  else if(strcmp(zone,"strait-of-sicily")==0) {
    map_set2Dgrid(&grid,-9.,32,-2+1./300,38+1./300,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"sicily")==0) {
    map_set2Dgrid(&grid,8.,30,20+1./300,40+1./300,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"caledonie")==0) {
    map_set2Dgrid(&grid,140.,-35.,185.,-5.,1./8.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"pacific")==0) {
    map_set2Dgrid(&grid,140.,-35.,185.,-5.,1./8.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"western-medsea")==0) {
    map_set2Dgrid(&grid,-5.6,31.,16.3,45.,0.1,0.1);
    zone_initialised=1;
    }

  else if(strcmp(zone,"eastern-medsea")==0) {
    map_set2Dgrid(&grid,15.,30.,36.5,42.5,1./60.,1./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"albicocca")==0) {
    map_set2Dgrid(&grid,5.,40.,12.,45.,0.02,0.02);
    zone_initialised=1;
    }

  else if(strcmp(zone,"cataluna")==0) {
    map_set2Dgrid(&grid,-0.5,38.5,6.5,44.,0.02,0.02);
    zone_initialised=1;
    }

  else if(strcmp(zone,"patagonia")==0) {
    map_set2Dgrid(&grid,-70.,-65.,-10.,0.,1./120.,1./120.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"Indien_austral")==0) {
    map_set2Dgrid(&grid,+43.,-74.,+90.,-35.,1./30.,1./30.);
    grid.modeH=0;
    zone_initialised=1;
    }

 else if(strcmp(zone,"kerguelen")==0) {
     /* Frame kerguelen */
    map_set2Dgrid(&grid,66,-52,78,-46,.01,.01);
    zone_initialised=1;
    }
  else if(strcmp(zone,"crozet")==0) {
    map_set2Dgrid(&grid,+51.5,-47.25,+52.5,-45.75,1./120.,1./120.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"plateau_Ker")==0) {
    map_set2Dgrid(&grid,+60.,-64.,+88.,-43.,1./30.,1./30.);
    grid.modeH=0;
    zone_initialised=1;
    }

  else if(strcmp(zone,"madagascar")==0) {
    map_set2Dgrid(&grid,30.,-35.,70.,5.,1./8.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"morbihan")==0) {
     /* Frame golfe-kerguelen */
    map_set2Dgrid(&grid,69,-50,71,-49,.01,.01);
    zone_initialised=1;
     }
   else if(strcmp(zone,"AIS-front")==0) {
    map_set2Dgrid(&grid,60,-70,82.5,-65,.05,.05);
    zone_initialised=1;
     }

  else if(strcmp(zone,"neiwp")==0) {
    map_set2Dgrid(&grid,86.,-24.,149.,42.,2./60.,2./60.);
    zone_initialised=1;
    }

  else if(strcmp(zone,"extended")==0) {
    map_set2Dgrid(&grid,-46.,21.,+30.,71.,2./60.,2./60.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"amazone")==0) {
    map_set2Dgrid(&grid,-60.5,-3.5,-39.,11.,1./60.,1./60.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"amazone-estuary")==0) {
    map_set2Dgrid(&grid,-52.5,-2.5,-45.,3.5,1./120.,1./120.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"amazone-river")==0) {
    map_set2Dgrid(&grid,-52.786,-1.720,-49.157,1.048,1./60.,1./60.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"para-estuary")==0) {
    map_set2Dgrid(&grid,-49.631,-2.448,-47.500,-2.448,1./60.,1./60.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"amazone-extended")==0) {
    map_set2Dgrid(&grid,-60.,-5.,-30.,15.,2./60.,2./60.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"persian")==0) {
    map_set2Dgrid(&grid,45.,22.58,60.,32.5,1./60.,1./60.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"macquarie")==0) {
    map_set2Dgrid(&grid,153.,-60.,160.,-50,0.05,0.05);
    zone_initialised=1;
    }
  else if(strcmp(zone,"macquarie_island")==0) {
    map_set2Dgrid(&grid,158.,-55.,160.,-54,0.002,0.002);
    zone_initialised=1;
    }
  else if(strcmp(zone,"terre_adelie")==0) {
    map_set2Dgrid(&grid,138.5,-69.,150.,-63.5,0.005,0.005);
    zone_initialised=1;
    }
  else if(strcmp(zone,"mertz")==0) {
    map_set2Dgrid(&grid,142.5,-69.,147.5,-65.5,0.002,0.002);
    zone_initialised=1;
    }
  else if(strcmp(zone,"mertz_large")==0) {
    map_set2Dgrid(&grid,135.5,-70.,160.,-63.5,0.01,0.004);
    zone_initialised=1;
    }
  else if(strcmp(zone,"east_antarctic")==0) {
    map_set2Dgrid(&grid,90.,-60.,160,-70.,1./60.,1./60.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"antarctic")==0) {
    map_set2Dgrid(&grid,-180.,-85.,180.,-55.,1./120.,1./120.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"arctic")==0) {
    map_set2Dgrid(&grid,-180.,55.,180.,90.,1./120.,1./120.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"arctic_ESA")==0) {
    map_set2Dgrid(&grid,-180.,54.5,180.,90.,1./120.,1./120.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"hudson")==0) {
    map_set2Dgrid(&grid,-96.0,50.0,-58.0,72.0,1./240.,1./240.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"siberia")==0) {
    map_set2Dgrid(&grid,30.0,60.0,210,80.0,1./120.,1./120.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"japan")==0) {
    map_set2Dgrid(&grid,122.0,24.0,148.0,46.0,1./120.,1./120.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"salomon")==0) {
    map_set2Dgrid(&grid,145.,-13.,162,0.,1./12.,1./12.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"new-zealand")==0) {
    map_set2Dgrid(&grid,157.0,-57.5,193.0,-24.0,1./240.,1./240.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"indochina-sea")==0) {
    map_set2Dgrid(&grid,95.0,-20,150.0,25.0,1./120.,1./120.);
    zone_initialised=1;
    }
  else if(strcmp(zone,"indeso")==0) {
    map_set2Dgrid(&grid,85.0,-30,150.0,30.0,1./60.,1./60.);
    zone_initialised=1;
    }

  if(identified!=0)
    *identified=(zone_initialised==1);
    
  if(zone_initialised!=1) {
    double resolution=atof(zone);
    if(resolution==0) STDERR_BASE_LINE("%s is neither a zone nor a resolution\n",zone);
    map_set2Dgrid(&grid,0.,-90.,+360.,+90.,resolution);
    }

  return(grid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool read_zone_arg(const char *keyword,const char *arg,frame_t *prescribed,double *dx,double *dy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(strcmp(keyword,"-xmin")==0) {
    sscanf(arg,"%lf",&prescribed->xmin);
    return true;
    }
  if(strcmp(keyword,"-xmax")==0) {
    sscanf(arg,"%lf",&prescribed->xmax);
    return true;
    }
  if(strcmp(keyword,"-ymin")==0) {
    sscanf(arg,"%lf",&prescribed->ymin);
    return true;
    }
  if(strcmp(keyword,"-ymax")==0) {
    sscanf(arg,"%lf",&prescribed->ymax);
    return true;
    }
  if(strcmp(keyword,"-dx")==0) {
    *dx=scan_angle(arg);
    return true;
    }
  if(strcmp(keyword,"-dy")==0) {
    *dy=scan_angle(arg);
    return true;
    }
  
  return false;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_zone_arg_help()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints zone options help for a programme help function. */
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf(
    "  -xmin,-xmax,-ymin,-ymax : followed by grid span parameters in degrees. Default: input grid span.\n"
    "  -dx,-dy : followed by output resolution parameters. Units recognised : d,m,',s,\". Default unit is degrees. Default value is the lowest reasonable value from 15\" to 30' while keeping the number of points >=1440e3.\n"
    );
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int apply_zone_arg(grid_t *grid, const frame_t &prescribed,double dx,double dy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// interpret zone options
/**
\param *grid Note: \c grid->modeH will NOT be updated!
\sa print_zone_arg_help()
\return 0 if nothing has been done and a bit mask otherwise.
*/
/*----------------------------------------------------------------------------*/
{
/*----------------------------------------------------------------------------------------------------
                     15"       20"       30"       1'     2'     5'     6'     10'     15'     30'    */
  const double preset[10]={1./60./4.,1./60./3.,1./60./2.,1./60.,2./60.,5./60.,6./60.,10./60.,15./60.,30./60.};
  int i;
  int status=0;
  
/*-----------------------------------------------------------------------------
  apply prescribed limits */
  if(isfinite(prescribed.xmin)){
    grid->xmin=prescribed.xmin;
    status|=1;
    }
  
  if(isfinite(prescribed.xmax)){
    grid->xmax=prescribed.xmax;
    status|=2;
    }
  
  if(isfinite(prescribed.ymin)){
    grid->ymin=prescribed.ymin;
    status|=4;
    }
  
  if(isfinite(prescribed.ymax)){
    grid->ymax=prescribed.ymax;
    status|=8;
    }

/*-----------------------------------------------------------------------------
  apply prescribed resolution*/
  if(isnormal(dx)){
    grid->dx=dx;
    status|=0x10;
    }
  if(isnormal(dy)){
    grid->dy=dy;
    status|=0x10;
    }
  
  if(!isnormal(grid->dx) || !isnormal(grid->dy)){
/*-----------------------------------------------------------------------------
    default resolution*/
    for(i=0; i<sizeof(preset)/sizeof(*preset); i++) {
      grid->dy=preset[i];
      grid->ny=round((grid->ymax-grid->ymin)/grid->dy)+1;
      grid->dx=preset[i];
      grid->nx=round((grid->xmax-grid->xmin)/grid->dx)+1;
      if(grid->nx*grid->ny>1200*1200) break;
      }

// /*-----------------------------------------------------------------------------
//     apply prescribed resolution*/
//     if(isnormal(dx)) grid->dx=dx;
//     if(isnormal(dy)) grid->dy=dy;
    
    status|=0x10;
    }
  
/*-----------------------------------------------------------------------------
  a possible issue */
  grid->xmin=round(grid->xmin/grid->dx)*grid->dx;
  grid->xmax=round(grid->xmax/grid->dx)*grid->dx;
  grid->ymin=round(grid->ymin/grid->dy)*grid->dy;
  grid->ymax=round(grid->ymax/grid->dy)*grid->dy;
  
  grid->nx=round((grid->xmax-grid->xmin)/grid->dx)+1;
  grid->ny=round((grid->ymax-grid->ymin)/grid->dy)+1;
  
  return status;
}
