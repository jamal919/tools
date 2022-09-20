
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Thierry Letellier  LEGOS, Toulouse, France (PhD)
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include <list> //for std::list<>
#include <vector> //for std::vector<>

#include "tools-structures.h"

#include "tides.h"
#include "tides.def"

#include <gsl/gsl_matrix.h> //for gsl_matrix


/* ---------------------------------------------------------------------- */
#define SPECTRUM_LIST_DEFAULT {\
  "ESTUARINE",\
  "ESTUARINE-HF",\
  "COASTAL",\
  "COASTAL-HF",\
  "SHELF",\
  "DEEP"\
}
#define SPECTRUM_DEFAULT_NB 6

extern const char* spectrum_list_default[];

/* *----------------------------------------------------------------------------
  commentaires ? */
typedef std::pair<tidal_wave,double> AdjacencyRule;
typedef std::list<AdjacencyRule> AdjacencyRuleList;
typedef std::pair<AdjacencyRuleList,double>  EnvirWave;


extern spectrum_t spectrum_init_ref_deep(void);
extern spectrum_t spectrum_init_ref_shelf(void);
extern spectrum_t spectrum_init_coastal(void);
extern spectrum_t spectrum_init_ref(string refClass, bool Z0, bool exitOnError=true);
extern spectrum_t spectrum_init_from_file(const string & filename,bool exitOnError);
extern spectrum_t spectrum_match(spectrum_t reference);
extern void spectrum_order(spectrum_t *s, string orderCrit);
extern void spectrum_order(spectrum_t *s, string orderCrit, std::vector<size_t> &perm);
extern void spectrum_reduce(spectrum_t &s, int* keep, int* deduce, FILE *out=0);
extern void spectrum_print(spectrum_t s, FILE *out);

extern gsl_matrix *adjacencyMatrix_build(spectrum_t s, double duration);
extern std::vector<EnvirWave> adjacencyList_build(spectrum_t s, double duration);
extern std::vector<tidal_wave> adjacencyList_reduce(std::vector<EnvirWave> &source, int nmes);
extern void adjacencyMatrix_plot(gsl_matrix *adjacencyMatrix);
extern void adjacencyMatrix_weight(gsl_matrix *adjacencyMatrix, spectrum_t s, double duration);
extern void edge_setWeight_byRaileigh(spectrum_t s, double duration, gsl_matrix **weight);
extern void vertex_setWeight_byMagnitude(spectrum_t s, double* &weight);
extern void vertex_setWeight_gotPotential(spectrum_t s, double* &weight);
extern void vertex_setWeight_gotNoPotential(spectrum_t s, double* &weight);
extern void vertex_setWeight_isBasewave(spectrum_t s, double* &weight);
extern void vertex_setWeight_asIndividual(spectrum_t s, double* &weight);

extern spectrum_t convert2spectrum_t(std::vector<tidal_wave> list);
//extern spectrum_t convert2spectrum_t(std::vector<EnvirWave> list);
extern void spectrum_check_bugged(spectrum_t *optimal, double duration, int nmes, int **keep, int **deduce, FILE *out=0);
