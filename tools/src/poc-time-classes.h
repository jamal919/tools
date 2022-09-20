#if POC_TIME_CLASSES_H == 0
#define POC_TIME_CLASSES_H 1


class date_t {
private :
public :
  int   day,month,year;   /* Gregorian date (define uniquely a given day)         */
  float second;           /* */
  
  void init(int Y=1950, int M=1, int D=1, float S=0.f){
    year=Y;
    month=M;
    day=D;
    second=S;
    }
    
  date_t(int Y=1950, int M=1, int D=1, float S=0.f){
    init(Y,M,D,S);
    }
    
  void set(int Y=1950, int M=1, int D=1, int H=0, int MN=0, double S=0.f) {
    S+=3600.*H+60.0*MN;
    init(Y,M,D,(float)S);
    }
    
  void add(double elapsed) {
    extern double cnes_time(date_t actual, const char units);
    extern  date_t poctime_getdatecnes(double t, char fmt);
    double time=cnes_time(*this,'d');
    time+=elapsed;
    *this=poctime_getdatecnes(time,'d');
    }
};

extern int isnad(const date_t &date);
extern int isad(const date_t &date);

  inline bool operator > (const date_t& __x, const date_t& __y) {
    if(isnad(__x) || isnad(__y))return 0;
    if(__x.year>__y.year) return(1);
    if(__x.year<__y.year) return(0);
    if(__x.month>__y.month) return(1);
    if(__x.month<__y.month) return(0);
    if(__x.day>__y.day) return(1);
    if(__x.day<__y.day) return(0);
    if(__x.second>__y.second) return(1);
    if(__x.second<__y.second) return(0);
    return(0);
    }

  inline bool operator < (const date_t& __x, const date_t& __y) {
    if(isnad(__x) || isnad(__y))return 0;
    if(__x.year<__y.year) return(1);
    if(__x.year>__y.year) return(0);
    if(__x.month<__y.month) return(1);
    if(__x.month>__y.month) return(0);
    if(__x.day<__y.day) return(1);
    if(__x.day>__y.day) return(0);
    if(__x.second<__y.second) return(1);
    if(__x.second>__y.second) return(0);
    return(0);
    }

  inline bool operator == (const date_t& __x, const date_t& __y) {
    if(isnad(__x) || isnad(__y))return 0;
    if(__x.year==__y.year && __x.month==__y.month &&
       __x.day==__y.day && __x.second==__y.second) return(1);
    return(0);
    }

  inline bool operator != (const date_t& __x, const date_t& __y) {
    if(isnad(__x) || isnad(__y))return 0;
    return !(__x==__y);
    }

#endif
