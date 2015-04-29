#include <qstring.h>
#include "params.h"

Params::Params()
{
	PARAM_SET params[] = {

{"DELTA_X", 10, 0, 0,
"Grid spacing",
"Size of lattice grid in um."},

{"DELTA_T", 1.0, 0, 0,
"Time step",
"Time step in minutes"},

{"MOTILITY_BETA", 0.01, 0, 0,
"Motility speed parameter",
"Cell motility is described by speed and persistence parameters, each in the range 0 - 1. \n\
 Median cell speed is roughly proportional to MOTILITY_BETA."},

{"MOTILITY_RHO", 0.25, 0, 0,
"Motility persistence parameter",
"Cell motility is described by speed and persistence parameters, each in the range 0 - 1. \n\
 MOTILITY_RHO determines the extent to which motion is in the same direction from one time step to the next."},

{"NDAYS", 1.0, 0, 0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"SEED1", 12345, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 56789, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 1, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation (currently only one used)."},

{"NLENGTH", 100, 0, 0,
"Tube length",
"Length of the tube in grids."},

{"NRADIUS", 10, 0, 0,
"Tube radius",
"Tube radius in grids"},

{"NPLUG", 10, 0, 0,
"Plug length",
"Plug length in grids"},

{"NT_ANIMATION", 10, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"DELAY", 100, 0, 0,
 "Simulation delay",
 "Sleep interval after each simulation time step (ms)"},

{"CHEMO_USED_1", 1, 0, 0,
 "Chemokine used",
 "Chemokine used"},

{"CHEMO_GRAD_AMP_1", 1.0, 0, 0,
 "Gradient amplitude",
 "Amplitude of the chemokine gradient (constant)"},

//{"CHEMO_USED_2", 0, 0, 0,
// "Chemokine used",
// "Chemokine used"},

//{"CHEMO_GRAD_AMP_2", 0.1, 0, 0,
// "Gradient amplitude",
// "Amplitude of the chemokine gradient (constant)"},

{"BG_FLOW_AMP", 1, 0, 0,
 "Background flow amplitude",
 "Amplitude of the background flow (constant)"},

{"KDRAG", 1, 0, 0,
"Drag weight: Kdrag",
"Weighting for drag vs. chemotaxis"},

{"KADHESION", 0.1, 0, 0,
"Adhesion weight: Kadh",
"Weighting for site attractiveness based on proximity"},

{"KSTAY", 0.2, 0, 0,
"Stay weight: Kstay",
"Weighting for extra attractiveness of current site - i.e. immobility"},

{"SAVE_CELL_POSITIONS", 0, 0, 0,
 "Saved cell positions",
 "Number of cells to save positions for at each time step in the log file"},

// Time-series plots
    {"NTcells",                 1, 0,1,"",""},
// Profile plots
    {"CD69",                    0, 0,1,"",""}

};
    nParams = sizeof(params)/sizeof(PARAM_SET);
    workingParameterList = new PARAM_SET[nParams];
    for (int i=0; i<nParams; i++) {
        workingParameterList[i] = params[i];
    }
}


PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
	workingParameterList[k].value = v;
}

void Params::set_label(int k, QString str)
{
	workingParameterList[k].label = str;
}
