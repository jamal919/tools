
#include <stdio.h>
#include <string.h>

typedef double *point;

int edges, inpoints, nBounds, *Bound_length;
double **coords;

/******************************************************************************/
int inbounds( point ref, int index )
/* Check to see if a point is inside a polygon by calculating the angles between
line segments from the point to successive polygon vertices. */
/*  NO LONGER USED!!!!!! */
{
  double x1, x2, y1, y2, opp, hyp, angles;
  int i, strtnode;
  
  printf( "A array %d is %d long.\n", index, Bound_length[index] );
  printf( "Ref. pt is %f %f.\n", ref[0], ref[1] );
  angles = 0.0;
  strtnode = 0;
  
  for( i = 0; i < index; i++ ) strtnode += Bound_length[i];
  
  printf( "Strtnode is %d.\n", strtnode );
  for( i = 0; i < Bound_length[index]; i++ ) {
    x1 = coords[strtnode + i][0];
    y1 = coords[strtnode + i][1];
    x2 = coords[strtnode + i + 1][0];
    y2 = coords[strtnode + i + 1][1];
    opp = sqrt( (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
    hyp = sqrt( (x1 - ref[0])*(x1 - ref[0]) + (y1 - ref[1])*(y1 - ref[1]));
    angles += asin( opp / hyp );
  }
  
  x1 = coords[strtnode + Bound_length[index]][0];
  y1 = coords[strtnode + Bound_length[index]][1];
  x2 = coords[strtnode][0];
  y2 = coords[strtnode][1];
  opp = sqrt( (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
  hyp = sqrt( (x1 - ref[0])*(x1 - ref[0]) + (y1 - ref[1])*(y1 - ref[1]));
  angles += asin( opp / hyp );
  
  if(( angles - floor( angles )) < 5.0 )
    return -1;
  else
    return 0;
}

/******************************************************************************/
int raybound( point ref, int index )
/*  Subroutine to check wether or not a point is inside a polygon.
The process is as follows:
        Use an arbitrary ray (here, y = constant and x >= xref), starting from
the point and going off to infinity.
        Count the number of polygon boundaries it crosses.
        If an odd number, the point is inside the polygon, otherwise it is
outside.   */
{
  int i, bcross, strt;
  double b, m, x;
  point pt1, pt2;

  bcross = 0; /* NUmber of boundary crossings. */
  strt = 0;   /* reference # for node starting the island */

  for( i = 0; i < index; i++ ) strt += Bound_length[i];

  for( i = 0; i < (Bound_length[index] - 1); i++ ) {
/* for each line segment around the island */
/*    printf( " Point number %d of %d.\n", i, Bound_length[index] ); */
    x = 0;
    m = 0;
    b = 0;
/* If the line segment's y values bracket the point's y value... */
    if((( coords[strt + i][1] > ref[1] ) &&
                ( coords[strt + i + 1][1] <= ref[1] )) ||
                (( coords[strt + i][1] < ref[1] ) &&
                ( coords[strt + i + 1][1] >= ref[1] ))) {
      pt1 = coords[strt + i];
      pt2 = coords[strt + i + 1];
/*
      printf(" pts are %f %f, %f %f.\n", pt1[0], pt1[1], pt2[0], pt2[1] );
      printf(" Reference point is %f, %f\n", ref[0], ref[1] );
*/
      /* If the line is not y = constant... */
      if( pt1[1] != pt2[1] ) {
        /* If the line is NOT x = const... */
        if( pt1[0] != pt2[0] ) {
        /* Find the equation of the line segment */
          m = (double)(( pt2[1] - pt1[1] ) / ( pt2[0] - pt1[0] ));
/*          printf( "M = %f\n", m ); */
          b = (double)(pt1[1] - m * pt1[0]);
/*          printf( "B = %f\n", b ); */
          x = (double)(( ref[1] - b ) / m);
/*          printf( "X = %f\n", x );  */
          if( x > ref[0]+1.e-06 ) {
            bcross++;
          }
        } else {
          /* the line is x = const and x > xref */
          if( pt1[0] > ref[0] ) {
            bcross++;
          }
        }
      }
    }
  }
          
/*  Check the line segment from the last point in the island boundary
 to the first point */
  if((( coords[strt][1] >= ref[1] ) &&
        ( coords[strt+Bound_length[index]-1][1] < ref[1] )) ||
        (( coords[strt][1] <= ref[1] ) &&
        ( coords[strt+Bound_length[index]-1][1] > ref[1] ))) {
    x = 0;
    m = 0;
    b = 0;
    pt1 = coords[strt + Bound_length[index] - 1];
    pt2 = coords[strt];
    if( pt1[0] != pt2[0] ) {
      m = (double)(( pt2[1] - pt1[1] ) / ( pt2[0] - pt1[0] ));
      b = (double)(pt1[1] - m * pt1[0]);
      x = (double)(( ref[1] - b ) / m);
      if( x > ref[0]+1.e-06  )  bcross++;
    } else {
      if(( pt1[0] > ref[0] ) && ((( pt1[1] <= ref[1] ) &&
           ( pt2[1] >= ref[1] )) || (( pt1[1] >= ref[1] ) &&
           ( pt2[1] <= ref[1] )))) {
        bcross++;
      }
    }
  }

/*  printf(" RAYBOUND: bcross = %d.\n", bcross%2 ); */
/*  Return the evenness/oddness of the boundary crossings
        i.e. the remainder from division by two. */
  return( bcross%2 );
}

/******************************************************************************/
double triarea( point val1, point val2, point val3 )
/*Calculates the area of a triangle using the vertex coordinates.
     The order of the points is important in the calculation!!!
     This is used to find in what order the points are traversed:
                Clockwise or Counter-Clockwise.
     A negative answer denotes clockwise traversal, correct for islands.
*/
{
  double ans;
  
  printf( "TRIAREA\n" );
  printf( "val1 = %f %f\n val2 = %f %f\n val3 = %f %f\n", val1[0], val1[1],
  val2[0], val2[1], val3[0], val3[1] );
  ans = ( val1[0] * val2[1] - val2[0] * val1[1] + val2[0] * val3[1] - val3[0] *
          val2[1] + val3[0] * val1[1] - val1[0] * val3[1] ) * 0.5;
          
  printf( "ans = %f.\n", ans );
  return ans;
}

/******************************************************************************/
void read_dot_nod( char *nodfile )
/* Read the data from the .nod file:
        -> Coordinates, # nodes, start nodes and length of boundaries, etc */
{
  int i, curr_bound, nodsread, freenodes, j;
  FILE *infile;
  double x, y, z;
  double *vals;
  
  printf(" INFILE = %s\n", nodfile );
  infile = fopen( nodfile, "r" );
  if( infile == (FILE *)NULL ) {
    __OUT_BASE_LINE__("   ERROR:  Cannot access file %s.\n", nodfile );
    exit( 1 );
  }
  
  fscanf( infile, "%d", &j );
  printf( "Number of nodes: %d\n", j );
  
/* Allocate memory for the coordinates:
        an array of pointers to arrays for the points  */
  vals = (double *)malloc( 3 * j * sizeof( double ));
  if( vals == (double *)NULL ) {
    __OUT_BASE_LINE__( "Error in assigning memeory to vals.\n" );
    exit( 1 );
  }
  coords = (double **)malloc( j * sizeof( double * ));
  if( coords == (double **)NULL ) {
    __OUT_BASE_LINE__( "Error in assigning memeory.\n" );
    exit( 1 );
  }
  for( i = 0; i < j; i++ ) {
    coords[i] = vals;
    vals += 3;
  }
  inpoints = j;
  
  fscanf( infile, "%d", &nBounds );
  Bound_length = (int *)malloc( nBounds * sizeof( int ));
  nodsread = 0;
  printf( "num BOUNDS = %d.\n", nBounds );
  
  i = 0;
  nodsread = 0;
/*  READ in the boundary nodes  */
  for( curr_bound = 0; curr_bound < nBounds; curr_bound++ ) {
    fscanf( infile, "%d", &Bound_length[curr_bound] );
/*    printf( "#nodes in boundary #%d is %d.\n", curr_bound,
         Bound_length[curr_bound] ); */
    for( i = 0; i < Bound_length[curr_bound]; i++ ) {
      fscanf( infile, " %lf %lf %lf", &x, &y, &z );
      /*Values are increased for calculations to be done later */
      coords[nodsread + i][0] = x;
      coords[nodsread + i][1] = y;
      coords[nodsread + i][2] = z;
/*      printf( "Point set at %f %f %f.\n", coords[nodsread + i][0], coords[nodsread +
      i][1], coords[nodsread + i][2] );*/
    }
    nodsread += Bound_length[curr_bound];
  }
  fscanf( infile, "%d", &freenodes );
/*  READ in the free nodes  */
  for( i = nodsread; i < inpoints; i++ ) {
    fscanf( infile, " %lf %lf %lf", &x, &y, &z );
/*    printf( "pt read is %f %f %f\n", x, y, z ); */
    coords[i][0] = x;
    coords[i][1] = y;
    coords[i][2] = z;
  }
  
  fclose( infile );
  edges = nodsread;
}

/******************************************************************************/

int node2poly(char *input, char *output)
{
  FILE *outfile;
  int i, j, noded, segmt, strtnode,flag;
  float x, y;
  point isle;
  double *h1, **holes,dx[2],dy[2],dn[2],alpha,limit;

  read_dot_nod(input);

  outfile = fopen(output, "w" );
  if( outfile == (FILE *)NULL ) {
    __OUT_BASE_LINE__( "Error opening output file.  Exitting.\n" );
    exit( 1 );
    }

//  fprintf( outfile, "%d 2 1 0\n", inpoints );
/*  Header for poly file:  #nodes total, #dimensions, #attribs, ?  */

  /* Output all the node data */
  for( i = 0; i < inpoints; i++ ) {
    fprintf( outfile, "%d %f %f %f\n", i, coords[i][0],
    coords[i][1], coords[i][2]);
    }

  /*Output #boundary nodes and list node #s that are connected */
  fprintf( outfile, "%d 0\n", edges );
  noded = 0;
  for( segmt = 0; segmt < nBounds; segmt++ ) {
    for( i = 0; i < (Bound_length[segmt]-1); i++ ) {
      fprintf( outfile, "%d %d %d\n", noded + i, noded + i, noded + i + 1 );
      }
    fprintf( outfile, "%d %d %d\n", noded + i, noded + i, noded );
    noded += Bound_length[segmt];
    }
  printf( "Segments are finished.\n" );

/*  ADD HOLES  */
/* Holes are used to eliminate the mesh inside boundaries -- islands */
  strtnode = 0;
  isle = (double *)malloc( 2 * sizeof( double ));

  /*allocate memory for coordinates of the holes */
  h1 = (double *)malloc( 2 * nBounds * sizeof( double ));
  if( h1 == (double *)NULL ) {
    __OUT_BASE_LINE__( "Error in assigning memory to vals.\n" );
    exit( 1 );
    }
  holes = (double **)malloc( nBounds * sizeof( double * ));
  if( coords == (double **)NULL ) {
    __OUT_BASE_LINE__( "Error in assigning memory.\n" );
    exit( 1 );
  }
  for( i = 0; i < nBounds; i++ ) {
    holes[i] = h1;
    h1 += 2;
  }

/*  for( segmt = 1; segmt < nBounds; segmt++ ) {
    strtnode += Bound_length[segmt - 1];
    for( j = 0; j < Bound_length[segmt]; j++ ) {
      fprintf( outfile, "%f %f\n", coords[strtnode+j][0], coords[strtnode+j][1] );
    }
    fprintf( outfile, "\n" );
  }*/

  for( segmt = 1; segmt < nBounds; segmt++ ) {
    strtnode += Bound_length[segmt - 1];
    limit=M_PI/10.0;
  redo:
    if( Bound_length[segmt] > 3 ) {
/* *----------------------------------------------------------------------------
     Find a triangle (3 consecutive boundary points) that has clockwise
     orientation.  The initial hole point will be the midpoint of the inner
     (on the island) side of the triangle.
----------------------------------------------------------------------------*/
      j = 1;
/*
     printf( "Isle is about to be defined. %f %f %f %f\n", coords[strtnode][0],     coords[strtnode][1],
                                                           coords[strtnode + 2][0], coords[strtnode + 2][1]);
*/
/* *----------------------------------------------------------------------------
     initiate process*/
      dx[0] = coords[strtnode + 1][0] - coords[strtnode + 0][0];
      dy[0] = coords[strtnode + 1][1] - coords[strtnode + 0][1];
      dx[1] = coords[strtnode + 1][0] - coords[strtnode + 2][0];
      dy[1] = coords[strtnode + 1][1] - coords[strtnode + 2][1];
      dn[0] = sqrt(dx[0]*dx[0]+dy[0]*dy[0]);
      dn[1] = sqrt(dx[1]*dx[1]+dy[1]*dy[1]);
      dx[0] /= dn[0];
      dy[0] /= dn[0];
      dx[1] /= dn[1];
      dy[1] /= dn[1];
      alpha=asin(dx[0]*dy[1]-dy[0]*dx[1]);
/*
     isle[0] = ( coords[strtnode][0] + coords[strtnode + 2][0] ) * 0.5;
     isle[1] = ( coords[strtnode][1] + coords[strtnode + 2][1] ) * 0.5;
*/
/* *----------------------------------------------------------------------------
      mid position point*/
      isle[0] = ( coords[strtnode][0] +coords[strtnode+1][0] + coords[strtnode + 2][0] ) /3.0;
      isle[1] = ( coords[strtnode][1] +coords[strtnode+1][1] + coords[strtnode + 2][1] ) /3.0;
      flag=raybound( isle, segmt );
/*
     printf( "Isle #%d %d, is %f %f.\n", segmt, j, isle[0], isle[1] );
*/
      while( !( (raybound( isle, segmt ) == 1) && ( fabs(alpha) > limit) ) && ( j < Bound_length[segmt] )) {
        j++;
/* *----------------------------------------------------------------------------
       mid position point*/
        isle[0] = (coords[strtnode+j-1][0] + coords[strtnode+j+1][0]) * 0.5;
        isle[1] = (coords[strtnode+j-1][1] + coords[strtnode+j+1][1]) * 0.5;
        flag=raybound( isle, segmt );
/* *----------------------------------------------------------------------------
       */
        dx[0] = coords[strtnode + j][0] - coords[strtnode + j-1][0];
        dy[0] = coords[strtnode + j][1] - coords[strtnode + j-1][1];
        dx[1] = coords[strtnode + j][0] - coords[strtnode + j+1][0];
        dy[1] = coords[strtnode + j][1] - coords[strtnode + j+1][1];
        dn[0] = sqrt(dx[0]*dx[0]+dy[0]*dy[0]);
        dn[1] = sqrt(dx[1]*dx[1]+dy[1]*dy[1]);
        dx[0] /= dn[0];
        dy[0] /= dn[0];
        dx[1] /= dn[1];
        dy[1] /= dn[1];
/* *----------------------------------------------------------------------------
       */
        alpha=asin(dx[0]*dy[1]-dy[0]*dx[1]);
        }
      if(j== Bound_length[segmt]) {
        limit/=2.0;
        goto redo;
        }
      printf( "segmt #%d, pt #%d, %f %f %f %d %d\n", segmt, j, isle[0], isle[1],alpha,j,Bound_length[segmt] );
      holes[segmt][0] = isle[0];
      holes[segmt][1] = isle[1];
      }
    else {
      /* If island has 3 or fewer points, no hole point needed */
      holes[segmt][0] = 0.0;
      holes[segmt][1] = 0.0;
    }
  }

  i = 0;
  for( segmt = 1; segmt < nBounds; segmt++ ) {
    if( holes[segmt][0] != 0.0 )  i++;
    }
  j = 1;
  fprintf( outfile, "%d\n", i );
  for( segmt = 1; segmt < (i+1); segmt++ ) {
    if( holes[segmt][0] != 0.0 ) {
      fprintf( outfile, "%d %lf %lf\n", j, holes[segmt][0], holes[segmt][1]);
      j++;
      }
    else {
      i++;
      }
    }

  printf( "Finished.\n" );
  close( outfile );
  return(0);
}

