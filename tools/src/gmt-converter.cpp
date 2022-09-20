
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

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  .

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "netcdf-proto.h"
#include "grd.h"
#include "map.h"

//#include "gmt_grd.h"

#define XYZ 0
#define YXZ 1

#define BOOLEAN int

/*
MGD-2000 Formats: GRD98

U.S. DEPARTMENT OF COMMERCE
NATIONAL OCEANIC AND ATMOSPHERIC ADMINISTRATION
NATIONAL ENVIRONMENTAL SATELLITE, DATA, AND INFORMATION SERVICE

THE GEODAS GRIDDED DATA FORMAT - "GRD98"

Dan R. Metzger

National Geophysical Data Center
Boulder, Colorado
February 1998


INTRODUCTION                   I

GENERAL DESCRIPTION            II

THE HEADER                     III

THE DATA                       IV

NGDC CONTACTS                  APPENDIX A

I.   INTRODUCTION



     During 1998 the Marine Geology & Geophysics Division of the
National Geophysical Data Center undertook an innovative and
exciting project. The NGDC Coastal Relief Model Project involves
gridding the hydrographic survey data compiled by the National
Ocean Service for US waters and, after combining this with US
Geological Survey land topography grids, creating a series of 1
degree grids at a 3 arc-second resolution in US coastal areas.
Thus a high-quality, high-resolution set of grids covering the
U.S land/sea coastal zone have become available for the first
time, allowing detailed study of this important area. This data
is available on CD/DVD sets with NGDC's GEODAS (GEophysical DAta
System) software which contains functionality for combining,
translating and sub-sampling the grids, as well as for creating
and viewing screen plots of the grid data.
     During this project a format for storing the grids was
naturally developed. It was decided, for space, speed and
simplicity to save the grids in a very utilitarian, no-nonsense
format, which became known as GRD98. This format is part of the
MGD-2000 (Marine Geophysical Data 2000) formats, a series of Year
2000 Compliant formats which also includes MGD77 (Marine
Geophysical Data Exchange Format), HYD93 (Hydrographic Surveys
Data Exchange Format) and ARO88 (Aeromagnetic Survey Header
Format).


II.   GENERAL DESCRIPTION

     The GRD98 Format, is a digital format for the storage of
gridded data. Though developed for bathymetric/topographic data,
the format can handle virtually any type of gridded data. It is
very utilitarian format and contains no documentation about the
grids (such as information about references, methods and datums
used, etc.). Rather GRD98 formatted files only contain grid-
structure information followed by the grid cell data values.
     GRD98 formatted files consist of header information followed
by a series of grid cell data values. The files contain binary
data only. The Header values are 4-byte signed integers which
describe the structure, size and extent of the grid cell values
that follow. The grid cells themselves can be 1-byte signed
integer, 2-byte signed  integer, 4-byte signed integer or 4-byte
floating point values.
     The grids described by GRD98 are node based. See IV. THE
DATA below for more information.
     The GRD98 format can be used for the exchange of grid data,
using virtually any media type.  The National Geophysical Data
Center uses CD/DVD disks as its chief method of distribution of
these data.





III.   THE HEADER

     The purpose of the Header  is to enumerate the structure,
size and extent of the grid cell values which follow it.

     The GRD98 Header is 128 bytes in length and consists of 32
binary signed 4-byte integers.

     In descriptions below, "original data values" refers to the
data used to construct the grid.

     Grid-Radius is a method for qualifying the data in grids.
Simply put, if a grid-radius of n cells is applied to a grid,
this means that any cells for which the nearest original data
value is more than n cells distant will be filled with the Empty
Grid Cell value. If grid-radius is not applied to a grid, all the
cells in the grid will contain values, no matter how far away the
nearest original data was.

     The upper-left corner is the origin of the grid. The grid
progresses row by row, top to bottom until the last row, with no
special terminating data. Within each row the grid progresses
column by column, left to right.

     The following is a detailed description of the Header fields


Name of Field                  Description
_____________________________________________________


Version
                    1,000,000,001 = version 1

Length
                    Length of the Header in bytes (128)
                    
Data Type
                    Describes what the cell values represent.
                    1 = Data, e.g. interpolated depths
                    2 = Data Density - density values for each cell,
                        (number of original data values falling in
                        the cell, as centered on cell node)
                    3 - Grid-Radius - grid-radius values for each cell,
                        (distance, in units of cells, to closest
                        original data point)

Latitude Degrees
                    Degrees portion of uppermost cell's latitude

Latitude Minutes
                    Minutes portion of uppermost cell's latitude

Latitude Seconds
                    Seconds portion of uppermost cell's latitude

Latitude Cell Size
                    Latitudinal size (height) of each cell in
                    seconds (i.e. distance in seconds between
                    cells)

Latitude Number of Cells
                    Number of rows in grid

Longitude Degrees
                    Degrees portion of leftmost cell's longitude

Longitude Minutes
                    Minutes portion of leftmost cell's longitude

Longitude Seconds
                    Seconds portion of leftmost cell's longitude

Longitude Cell Size
                    Longitudinal size (width) of each cell in
                    seconds (i.e. distance in seconds between
                    cells)

Longitude Number of Cells
                    Number columns in grid

Minimum Value
                    Minimum value of all cells in grid, excluding
                    empty grid cells. Per precision, i.e. based
                    on actual numbers found in cells.  E.g.  if
                    the lowest cell value found is -123 and
                    precision is 10ths of meters (-12.3 meters)
                    the Minimum Value is -123

Maximum Value
                    Maximum value of all cells in grid, excluding
                    empty grid cells. Per precision (see above)

Grid Radius
                    The Grid-Radius which was applied to the grid.
                    When a Grid-Radius of n is applied this means
                    that only cells which are within n cells of
                    actual data will be filled with data values.
                    Cells for which real data is more than n cells
                    away will be given the Empty Grid Cell Value.
                    If the Grid-Radius equals -1 then Grid-Radius
                    was not applied.

Precision
                    Precision of the cell data values.
                    1 = whole units
                    10 = tenths of units

Empty Grid Cell Value
                    Value placed in a cell with no data (e.g.
                    grid-radius was applied to a Data Grid). For
                    Density Grids or Grid-Radius Grids, land (as
                    opposed to water) cells would contain this
                    value if density and grid-radius were not
                    calculated for land cells.

Number Type
                    Byte size of cell data values, positive =
                    integers, negative = floats
                    (e.g. +2 = 2-byte integers)

Water Datum
                    The vertical datum used for non-land cell depths.
                    Local datums could be used for inland lakes, with
                    land values tied to mean sea level.  Using local
                    water datums means that water shore values will be
                    zero. Using MSL for water means that water shore
                    values will match up with land shore values.
                    0 = Mean Sea Level
                    1 = Local Vertical Datum used for depths.


Data Value Limit
                    This is the maximum possible value for cell data.
                    This is used when there is a limit to the
                    calculated values, e.g. to keep values within
                    a the range of Number Type. Cell values containing
                    this number mean that the value is this large or
                    larger.  0 = not applied

Cell Registration
                   Gridline-registered = 0
                   Pixel-registered (cell centered) = 1
                   e.g. a 1 min cell size grid with upper-left
                   lat/lon = 60 deg 0 min 0 sec / 45 deg 0 min 0 deg
                   Gridline-registered: 1st value at 60 0 0 / 45 0 0
                   Pixel-registered: 1st value at 59 59 30 / 45 0 30


Unused (10 fields)
                    Unused field; set to zero. Repeat 10 times





IV.   THE DATA



The GRD98 Header is followed immediately by the grid cell data
values. The cells progress in a row by row manner, starting with
the topmost (northernmost) row and continuing downward to the
bottommost row. Within each row the cells progress column by
column from leftmost (westernmost) to rightmost. The exact
structure, size and extent of the data values themselves are
described by the Header.

What the data values actually represent depends on the Grid Type
as enumerated in the Header Version. Generally the GRD98 format
is used for Data Grids, but it can also be used to describe
Density Grids  and Grid-Radius Grids (see Version in III THE
HEADER above for details).

GRD98 formatted grids are node based. The lat/lon position for a
specific grid data value represents the center of a grid "cell"
which extends half the grid-cell-size in 4 directions. The data
value corresponds to this exact position (NOT to a position a
half cell horizontal and a half cell vertical away). These
positions are set up so that they will match up with exact
latitude (and longitude) whole-degree values as they progress
across the rows ( or columns). For example, a grid whose upper
left corner was at 0 degrees latitude, 0 degrees longitude would
have as it's first grid cell a lat/lon position at exactly 0,0.
If this grid had an extent of 1 degree by 1 degree, and latitude
and longitude cell sizes of 60 and 60 (seconds), the grid would
be 61 rows by 61 columns.

The data values are one of the following number types (as
described in the Header):
     1-byte signed integer
     2-byte signed integer
     4-byte signed integer
     4-byte floating point

The data values must be interpreted as per the Precision
enumerated in the Header. For example if a data cell value is an
integer number equal to 12345 and the precision in the Header
equals 10 (tenths of units), then 1234.5 is the actual real-world
value for that cell. For floating point grids the data values do
not depend on the Precision; the value found is the actual real-
world value as is.

Where the grid cell is meant to contain "no data", the data value
for that cell will be the Empty Grid Cell Value as enumerated in
the Header. This can occur when Grid-Radius was applied to the
grid in order to keep out data values which are too far away from
original (real) data. For Density Grids and Grid-Radius Grids,
(see III THE HEADER above) this Empty Grid Cell Value could
correspond to "land" positions, where density or grid-radius
values were not calculated.




APPENDIX A   NGDC CONTACTS
Dan R. Metzger: (303) 497-6542  Dan.R.Metzger@noaa.gov
  or
David L. Divins (303) 497-6505  David.Divins@noaa.gov
  
National Geophysical Data Center
NOAA, E/GC3
325 Broadway
Boulder, CO 80305-3328

TELEX 592811 NOAA MASC BDR
FAX (303) 497-6513

____________________________________________________________

U.S. DEPARTMENT OF COMMERCE
NATIONAL OCEANIC AND ATMOSPHERIC ADMINISTRATION
NATIONAL ENVIRONMENTAL SATELLITE, DATA, AND INFORMATION SERVICE

THE GEODAS 2-D VECTOR DATA FORMAT - "VCT00"

Compiled By:
Dan R. Metzger


National Geophysical Data Center
Boulder, Colorado
August 2000


INTRODUCTION                   I

GENERAL DESCRIPTION            II

HEADERS (Binary VCT00)         III

DATA RECORDS (Binary VCT00)    IV

DATA RECORDS (ASCII VCT00)     V

NGDC CONTACTS                  APPENDIX A



I.   INTRODUCTION

     During 2000 the Marine Geology & Geophysics Division of the
National Geophysical Data Center undertook the development of
software for coastline data.  In doing so a generalized format
for vector data was implemented. This format is geared toward
coastline data, but can be used for other vector applications.



II.   GENERAL DESCRIPTION

     The VCT00 Format, is a digital format for the storage of 2-D
vector data. Though developed for high resolution boundary
(coastline) data, the format can handle other types of vector
data. It is very utilitarian format and contains no documentation
about the data (such as information about references, methods and
datums used, etc.). Rather VCT00 formatted files only contain
vector data records and headers with minimal geographic and
file-indexing information.
     VCT00 formatted files can be either binary or ASCII.  VCT00
binary files consist of a series of Header Records followed by
Data Records.  VCT00 ASCII files contain only the Data Records.
The Binary VCT00 "physical records", both Header and Data, are 10
bytes each, and contain two signed 4-byte integers followed by
one signed 2 byte integer. The Headers contain four of these
physical records; the Data records consist of one record for each
"point".  Each Header represents a geographic box or "Block"
containing a series of line segments. Each Header contains the
"Address" (record number) where those data segments start in the
file, the "Type" of data the segments represent, an optional
"Value" for the data, the total Number of Points in the data
segments, and the latitude/longitude limits of the Block of data.
Data Records for both Binary and ASCII VCT00 files contain
resolution information, so a single file can be plotted at
several different resolutions, by sub-setting the data for
different resolutions.
     The VCT00 format can be used for the exchange of 2-D vector
data, using virtually any media type.  The National Geophysical
Data Center uses CD-ROM disks as its chief method of distribution
of these data.



III.   HEADERS (Binary VCT00)

     The purpose of each binary Header is to locate and describe
the vector data which it (geographically) encompasses.

     Each VCT00 Header is 40 bytes in length and consists of 4 10-
byte records, each of which contains two signed 4-byte integers
followed by one signed 2 byte integer.

     The Header records are all together at the start of the
file.

     The last Header Record is an "empty" Header. This indicates
that there are no more Headers. Data records immediately follow
this Header.

     The following is a detailed description of the Header fields


FIELD NAME         NUMBER TYPE     DESCRIPTION
________________________________________________________

Address of Data    4 Byte Integer  File record number of
                                   the first point of the vector data
                                   which the Header represents. Record
                                   numbers start at 1 and are10 bytes
                                   each.

Number of Records  4 Byte Integer  Total number of
                                   points (vertices) represented by
                                   this Header

Type of Data       2 Byte Integer  For delineating
                                   different types of vector data. e.g.
                                   for geographic boundary data,
                                   1=coastlines 2=primary international
                                   boundaries, 3= secondary boundaries
                                   such as US states, etc. A Type of -1
                                   is used for the last (no records)
                                   Block of the file.

Data Attribute    4 Byte Integer   A value can be
or Value                           attached to the data which the
                                   Header represents. E.g. for
                                   toptographic contours, this could
                                   represent a contour value for all
                                   the contours in this Header.

Unused           4 Byte Integer    Unused at this time.

Unused           2 Byte Integer    Unused at this time.

Upper Left       4 Byte Integer    Topmost Latitude of
Latitude                           this Header block, in millionths of
                                   degrees (decimal degrees times
                                   1,000,000).

Upper Left       4 Byte Integer    Leftmost Longitude of
Longitude                          this Header block, in millionths of
                                   degrees.

Unused           2 Byte Integer    Unused at this time.

Lower Right      4 Byte Integer    Bottommost Latitude
Latitude                           of this Header block, in millionths
                                   of degrees.

Lower Right      4 Byte Integer    Rightmost Longitude
Longitude                          of this Header block, in millionths
                                   of degrees.

Unused           2 Byte Integer    Unused at this time.



IV.   DATA (Binary VCT00)
         
The binary VCT00 Headers are followed immediately by the binary
vector data records. The data records for a specific Header can
only be delineated by the information in the Header record. The
data records for a specific Header are broken up into a number of
line segments. The "Pencode" field is used to determine extent of
these line segments. The Type, (optional) Value and number of
points of the line segment data for a particular Header are
described in the Header Record.

The data records are 10 bytes and each record contains two signed
4-byte integers followed by one signed 2 byte integer.

The plotting instructions (MoveTo or LineTo) are contained within
the Pencode Field. The Pencode also contains resolution
information. Each Non-Zero pencode value is evenly divisible by a
number of Resolution Values, meaning that the lat/lon point
"belongs" to those particular resolutions.  This Resolution Value
is one of the first 7 prime numbers, the value of which increases
with lower Resolution. To plot at a particular resolution you
would only plot those points which have a Pencode which is evenly
divisible by your particular Resolution Value. For full
resolution this Resolution Value would be 1. The remaining
Resolution Values increase through the Prime Numbers as the
desired resolution decreases (i.e. 2, 3, 5, 7, 11, 13). Thus a
simple Modulus Function can be used to determine whether to plot
specific point.  The usable resolutions can be thought of as
Full, High, Medium-High, Medium, Medium-Low, Low, Crude.

     The following is a detailed description of the Binary Data
Record fields


FIELD NAME     NUMBER TYPE     DESCRIPTION
________________________________________________________

Latitude       4 Byte Integer  Latitude of the
                               point, in millionths of degrees
                               (decimal degrees times 1,000,000).

Longitude      4 Byte Integer  Longitude of the
                               point, in millionths of degrees.

Pencode        2 Byte Integer  Delineates the flow of the segments:
                               0 - Start, MoveTo Pt (Pen up Move)
                               1 - Full Res, LineTo Pt (Pen Down Move)
                               Factor of 2 - Second Level Res LineTo Pt
                               Factor of 3 - Third Level Res LineTo Pt
                               Factor of 5 - Fourth Level Res LineTo Pt
                               etc.



                             
V.   DATA (ASCII VCT00)

ASCII VCT00 data files do not contain Headers. The ASCII VCT00
files are simply as series of records, broken up into line
segments according to the pencode.

The data records are 28 ASCII characters in fixed fields, plus
computer's end-of-record indicator. There will always be at least
one space between fields, allowing the data to be input as space-
delimited data.  Decimal points for latitude and longitude are
explicit.

The plotting instructions (MoveTo or LineTo) are contained within
the Pencode Field. The Pencode also contains resolution
information. Each Non-Zero pencode value is evenly divisible by a
number of Resolution Values, meaning that the lat/lon point
"belongs" to those particular resolutions.  This Resolution Value
is one of the first 7 prime numbers, the value of which increases
with lower Resolution. To plot at a particular resolution you
would only plot those points which have a Pencode which is evenly
divisible by your particular Resolution Value. For full
resolution this Resolution Value would be 1. The remaining
Resolution Values increase through the Prime Numbers as the
desired resolution decreases (i.e. 2, 3, 5, 7, 11, 13). Thus a
simple Modulus Function can be used to determine whether to plot
specific point.  The usable resolutions can be thought of as
Full, High, Medium-High, Medium, Medium-Low, Low, Crude.

     The following is a detailed description of the ASCII Data
Record fields


FIELD NAME     NUMBER OF BYTES     DESCRIPTION
________________________________________________________

Longitude         11          Longitude of the point,
                              with precision up to 6 decimals.

Latitude          11          Latitude of the point,
                              with precision up to 6 decimals.

Pencode            6          Delineates the flow of the segments
                              0 - Start, MoveTo Pt (Pen up Move)
                              1 - Full Res, LineTo Pt (Pen Down Move)
                              Factor of 2 - Second Level Res LineTo Pt
                              Factor of 3 - Third Level Res LineTo Pt
                              Factor of 5 - Fourth Level Res LineTo Pt
                              etc.


-----------------------------------------------------------------

APPENDIX A   NGDC CONTACTS

Dan R. Metzger: (303) 497-6542  Dan.R.Metzger@noaa.gov
  or
David L. Divins (303) 497-6505  David.Divins@noaa.gov


National Geophysical Data Center
NOAA, E/GC3
325 Broadway
Boulder, CO 80305

TELEX 592811 NOAA MASC BDR
FAX (303) 497-6513

____________________________________________________________
*/

typedef struct {
  frame_t  frame;
  short    type;
  size_t   pos;
  size_t   ndata,nlines;
} vct00_header_t;

typedef struct {
  short   pen;
  double  lon,lat;
} vct00_data_t;

typedef struct {
  vct00_header_t *headers;
  vct00_data_t   **data;
  size_t   nblocks;
} vct00_t;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int grd98_loadgrid (char *filename, grid_t *grid, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,option,status;
  int i,j,k,n;
  char *masked=0;
  char line[1024];
  int  val[32];

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

  for (k=0;k<32;k++) {
    fread(&(val[k]),sizeof(int),1,in);
    }

  grid->nx=val[12];
  grid->ny=val[7];

  grid->dx=(double) val[11] / (double) 3600.0;
  grid->dy=(double) val[6]  / (double) 3600.0;

  grid->xmin=(double) val[8] + (double) val[9] / (double) 60.0 + (double) val[10] / (double) 3600.0;
  grid->ymax=(double) val[3] + (double) val[4] / (double) 60.0 + (double) val[5]  / (double) 3600.0;

  grid->xmax=grid->xmin+grid->dx*(grid->nx-1);
  grid->ymin=grid->ymax-grid->dy*(grid->ny-1);

  grid->modeH=0;

  grid->nz=1;

  fclose(in);

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int grd98_load_r1(char *filename, grid_t grid, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,option,status;
  int i,j,k,n;
  char *masked=0;
  char line[1024];
  int  type,data,val[32];
  float precision;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

  for (k=0;k<32;k++) {
    fread(&(val[k]),sizeof(int),1,in);
    }

  precision=(float) val[16];

  *mask=(float) val[17];
  type=val[18];

  switch (type) {
    case -4:
      for (j=0;j<grid.ny;j++) {
        for (i=0;i<grid.nx;i++) {
           n=(grid.ny-j-1)*grid.nx+i;
           fread(&buf[n],sizeof(int),1,in);
           if(buf[n]!=*mask) buf[n]=buf[n]/precision;
           }
         }
       break;

    case 4:
      for (j=0;j<grid.ny;j++) {
        for (i=0;i<grid.nx;i++) {
           n=(grid.ny-j-1)*grid.nx+i;
           fread(&data,sizeof(int),1,in);
           buf[n]=data/precision;
           }
         }
       break;
     }

  fclose(in);

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int vct00_inquire(char *filename, vct00_t *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,option,status;
  int i,j,k,m,n;
  char *masked=0;
  char line[1024];
  int  type,data,val[32];
  short pen,sval[10];
  float precision;
  double x,y;
  grid_t grid;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

/* *------------------------------------------------------------------------------
  scan numbre of blocks */
  buffer->nblocks=0;
  type=0;
  while (type != -1) {
    for (k=0;k<2;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[0]),sizeof(short),1,in);

    for (k=2;k<4;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[1]),sizeof(short),1,in);

    for (k=4;k<6;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[2]),sizeof(short),1,in);

    for (k=6;k<8;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[3]),sizeof(short),1,in);

    type=sval[0];
    buffer->nblocks++;
    }

  buffer->nblocks--;

  buffer->headers=new vct00_header_t[buffer->nblocks];
  buffer->data=new vct00_data_t *[buffer->nblocks];

  rewind(in);

  for (n=0;n<buffer->nblocks;n++) {
    for (k=0;k<2;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[0]),sizeof(short),1,in);

    buffer->headers[n].pos  =val[0];
    buffer->headers[n].ndata=val[1];

    for (k=2;k<4;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[1]),sizeof(short),1,in);

    for (k=4;k<6;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[2]),sizeof(short),1,in);

    for (k=6;k<8;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[3]),sizeof(short),1,in);

    buffer->headers[n].frame.ymax=(double) val[4]/1000000.0;
    buffer->headers[n].frame.xmin=(double) val[5]/1000000.0;

    buffer->headers[n].frame.ymin=(double) val[6]/1000000.0;
    buffer->headers[n].frame.xmax=(double) val[7]/1000000.0;

    buffer->headers[n].type=sval[0];
    }

  for (n=0;n<buffer->nblocks;n++) {
    for (m=0;m<buffer->headers[n].ndata;m++) {
      for (k=8;k<10;k++) {
        fread(&(val[k]),sizeof(int),1,in);
        }
      fread(&(sval[3]),sizeof(short),1,in);
      pen=sval[3];
/* *----------------------------------------------------------------------------
Pencode            6          Delineates the flow of the segments
                              0 - Start, MoveTo Pt (Pen up Move)
                              1 - Full Res, LineTo Pt (Pen Down Move)
                              Factor of 2 - Second Level Res LineTo Pt
                              Factor of 3 - Third Level Res LineTo Pt
                              Factor of 5 - Fourth Level Res LineTo Pt
                              etc.
-----------------------------------------------------------------------------*/
      switch (pen) {
        case 0:
          break;

        default:
          break;
        }
      }
    }

  fclose(in);

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int gmt_load_cst(char *filename, vct00_t *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,option,status;
  int i,j,k,m,n;
  char *masked=0;
  char line[1024];
  int  type,data,val[32];
  short pen,sval[10];
  float precision;
  double x,y;
  grid_t grid;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

/* *------------------------------------------------------------------------------
  scan numbre of blocks */
  buffer->nblocks=0;
  type=0;
  while (type != -1) {
    for (k=0;k<2;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[0]),sizeof(short),1,in);

    for (k=2;k<4;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[1]),sizeof(short),1,in);

    for (k=4;k<6;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[2]),sizeof(short),1,in);

    for (k=6;k<8;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[3]),sizeof(short),1,in);

    type=sval[0];
    buffer->nblocks++;
    }

  buffer->nblocks--;

  buffer->headers=new vct00_header_t[buffer->nblocks];
  buffer->data=new vct00_data_t *[buffer->nblocks];

  rewind(in);

  for (n=0;n<buffer->nblocks;n++) {
    for (k=0;k<2;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[0]),sizeof(short),1,in);

    buffer->headers[n].pos  =val[0];
    buffer->headers[n].ndata=val[1];

    for (k=2;k<4;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[1]),sizeof(short),1,in);

    for (k=4;k<6;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[2]),sizeof(short),1,in);

    for (k=6;k<8;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[3]),sizeof(short),1,in);

    buffer->headers[n].frame.ymax=(double) val[4]/1000000.0;
    buffer->headers[n].frame.xmin=(double) val[5]/1000000.0;

    buffer->headers[n].frame.ymin=(double) val[6]/1000000.0;
    buffer->headers[n].frame.xmax=(double) val[7]/1000000.0;

    buffer->headers[n].type=sval[0];
    }

  for (n=0;n<buffer->nblocks;n++) {
    buffer->data[n]=new vct00_data_t[buffer->headers[n].ndata];
    for (m=0;m<buffer->headers[n].ndata;m++) {
      for (k=8;k<10;k++) {
        fread(&(val[k]),sizeof(int),1,in);
        }
      fread(&(sval[3]),sizeof(short),1,in);
      buffer->data[n][m].lon=(double) val[9]/1000000.0;
      buffer->data[n][m].lat=(double) val[8]/1000000.0;
      buffer->data[n][m].pen=sval[3];
      }
    }

  fclose(in);

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  geo_t projection;

  float  topo_mask=1.e+10,mask=1.e+10;
  double  x,y,t,p;
  double  dx,dy;

  int i,j,k,l,m,n,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
  int flag;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL;
  char *shorelines=NULL,*bathymetry=NULL;
  char *bathy=NULL,*notebook=NULL,*meshbook=NULL,*poly=NULL;

  grid_t topogrid;
  grid_t grid;
  mesh_t mesh;
  short *topo,smask=256*127+255;;
  float *ftopo,*tmp,*buffer;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  int dimlgth[4]={10,20,30,40};

  vct00_t shorelines_data;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'c' :
          shorelines= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
          n++;
          }
        else
          {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

/*-----------------------------------------------------------------------------
  read shoreline database */
  status=gmt_load_cst(shorelines, &shorelines_data);

/*-----------------------------------------------------------------------------
  read topo modifiable database */
  status=grd98_loadgrid(input,&topogrid,buffer,&mask);
//  status=gmt_read_grid(input,&topogrid,buffer,&mask);

  if(status !=0) {
    __OUT_BASE_LINE__("cannot load grid in bathymetry file=%s\n",input);
    exit(-1);
    }

  buffer=new float[topogrid.nx*topogrid.ny];
  status=grd98_load_r1(input,topogrid,buffer,&mask);
  if(status !=0) {
    __OUT_BASE_LINE__("cannot load data in bathymetry file=%f\n",input);
    exit(-1);
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data, to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  save rectified topo */
  if(rootname==0) rootname=strdup("out");
  output=new char[1024];

  ftopo=buffer;
   grid=topogrid;

  topo=new short[grid.nx*grid.ny];
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      n=(grid.ny-j-1)*grid.nx+i;
      if(ftopo[m]!=mask) {
        topo[n]=(short) floor(-ftopo[m]+0.5);
        }
      else {
        topo[n]=256*127+255;
        }
      }
    }

  sprintf(output,"%s-short.grd",rootname);
  status=grd_save(output,topogrid,topogrid.nx, topo,smask);

  sprintf(output,"%s-float.grd",rootname);
  status=grd_save(output,topogrid,topogrid.nx, ftopo, mask);

  free(topo);

  __OUT_BASE_LINE__("bathymetry sucessfully completed\n");

  exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);

}
