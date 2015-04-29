#ifndef GLOBAL_H
#define GLOBAL_H

#include <QtGui>
#include <QMutex.h>
#include "profile.h"

#define MAX_TC 100000

namespace Global
{
    extern QMutex mutex1, mutex2;

    extern int NX, NY, NZ;
    extern double DELTA_T;
    extern int ncpu;
    extern int nsteps;
    extern int nt_vtk;
    extern int istep;  // ???
    extern bool leftb;
    extern int delay;
    extern int summary_interval;

    extern int nvars_used;
    extern int GUI_to_DLL_index[32];
    extern int DLL_to_GUI_index[32];
    extern QString var_string[32];

    extern double *FACS_data;
    extern int nFACS_cells;
    extern int nFACS_dim;
    extern int nFACS_vars;

    extern double *histo_data;
    extern double *histo_data_log;
    extern int nhisto_bins;
    extern int nhisto_dim;
    extern double histo_vmin[3*32];
    extern double histo_vmax[3*32];
    extern double histo_vmin_log[3*32];
    extern double histo_vmax_log[3*32];
    extern int histo_celltype;

    extern double *profile_x[20];
    extern double *profile_y[20];
    extern int profile_n[20];

    extern int summaryData[100];

//    double concData[4000];
//    int conc_nc;
//    double conc_dx;

    extern bool showingVTK;
    extern bool recordingVTK;
    extern bool showingFACS;
    extern bool recordingFACS;

    extern int VTKbuffer[100];
    extern int *TC_list;
    extern int nTC_list;
    extern int tubeLength, tubeRadius, nLines;

    extern bool redimflag;
    extern int recordfrom;
    extern int recordto;

} // namespace Global

#endif // GLOBAL_H
