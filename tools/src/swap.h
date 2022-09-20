#include "tools-structures.h"
#include "poc-time-classes.h"
#include "map.h"

extern short swap(short v);
extern int swap(int v);
extern unsigned swap(unsigned v);
extern float swap(float v);
extern double swap(double v);

extern short shswap( short shrtint );
extern int lswap( int long_word );
extern float fswap( float long_word );
extern void f3swap(float *long_word);
extern double dswap2( double dword );
extern fcomplex cswap( fcomplex dword );
extern grid_t gswap( grid_t grid);
extern date_t date_swap( date_t date);
extern void swap_header( bmgheader_t *header );
