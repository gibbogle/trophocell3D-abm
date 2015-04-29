#include "global.h"

namespace Global
{

    QMutex mutex1, mutex2;

    int NX, NY, NZ;
    double DELTA_T;
    int ncpu;
    int nsteps;
    int nt_vtk;
    int istep;  // ???
    bool leftb;
    int delay=0;
    int summary_interval;

    int nvars_used;
    int GUI_to_DLL_index[32];
    int DLL_to_GUI_index[32];
    QString var_string[32];

    double *FACS_data;
    int nFACS_cells;
    int nFACS_dim;
    int nFACS_vars;

    double *histo_data;
    double *histo_data_log;
    int nhisto_bins;
    int nhisto_dim;
    double histo_vmin[3*32];
    double histo_vmax[3*32];
    double histo_vmin_log[3*32];
    double histo_vmax_log[3*32];
    int histo_celltype;

    double *profile_x[20];
    double *profile_y[20];
    int profile_n[20];

    int summaryData[100];

//    double concData[4000];
//    int conc_nc;
//    double conc_dx;

    bool showingVTK;
    bool recordingVTK;
    bool showingFACS;
    bool recordingFACS;

    int VTKbuffer[100];
    int *TC_list;
    int nTC_list;
    int tubeLength, tubeRadius, nLines;

    bool redimflag;
    int recordfrom;
    int recordto;

} // namespace Global
