#ifndef RESULT_SET_H
#define RESULT_SET_H
struct result_set {
	QString casename;
	int nsteps;
 	double *tnow;      // simulation time (mins)
    double *nDC;       // number of DCs
    double *act;       // total DC antigen activity level
    double *ntot_LN;   // total T cell population in the LN
    double *ncog_PER;   // total T cell population in the periphery
    double *ncogseed;  // number of naive cognate T cells that have arrived
    double *ncog_LN;      // current number of cognate T cells
    double *ndead;     // number of cognate T cells that have died
    double *teffgen;   // number of activated cognate T cells that have left the LN
    double *nbnd;     // number of cognate T cells that are bound to a DC
    double max_tnow;      // simulation time (mins)
    double max_nDC;       // number of DCs
    double max_act;       // total DC antigen activity level
    double max_ntot;      // total T cell population
    double max_ncogseed;  // number of naive cognate T cells that have arrived
    double max_ncog_LN;      // current number of cognate T cells
    double max_ncog_PER;      // current number of cognate T cells
    double max_ndead;     // number of cognate T cells that have died
    double max_teffgen;   // number of activated cognate T cells that have left the LN
    double max_nbnd;     // number of cognate T cells that are bound to a DC

    double *pData[32];
    double maxValue[32];
};
typedef result_set RESULT_SET;
#endif
