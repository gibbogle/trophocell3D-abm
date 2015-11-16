#include <qstring.h>
#include "params.h"

Params::Params()
{
	PARAM_SET params[] = {

{"DELTA_T", 10.0, 0, 0,
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

{"SETTLING_TIME", 1.0, 0, 0,
"Settling time",
"Settling time in hours."},

{"SEED1", 12345, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 56789, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 4, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation (currently only one used)."},

{"TUBE_LENGTH", 1000, 0, 0,
"Tube length",
"Length of the tube (um)."},

{"TUBE_RADIUS", 100, 0, 0,
"Tube radius",
"Tube radius (um)"},

{"PLUG_ZMIN", 20, 0, 0,
"Plug minimum z",
"Plug minimum z (um)"},

{"PLUG_ZMAX", 120, 0, 0,
"Plug maximum z",
"Plug maximum z (um)"},

{"PLUG_HMAX", 50, 0, 0,
"Plug maximum height",
"Plug maximum height (um)"},

{"RAVERAGE", 5, 0, 0,
"Average cell radius",
"Average cell radius (um)"},

{"NT_ANIMATION", 1, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).Global::"},

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

    {"A_SEPARATION", 1.0, 0, 0,
    "Separation force factor",
    "During mitosis the two capped spheres are effectively connected by a nonlinear spring. \n\
    The length of the spring s is determined by the mitosis level, and if the centre-centre distance is d \n\
    the contribution to the force of repulsion between the spheres is a_separation*(s-d)^3."},

    {"A_FORCE", 1., 0, 0,
    "Repulsion force factor 'a'",
    "The cell-cell force is a function of x = distance/(sum of radii): x = d/(R1+R2). \n\
    The force function is: F(x) = a/((x-x0)(x1-x)) + b, where x0 and x1 are the locations of the bounding asymptotes. \n\
    The parameter 'b' is calculated by setting the minimum value of F(x) (which occurs at x = (x0+x1)/2 ) equal to -c, this is the maximum attraction force. \n\
    The function F(x) is zero at two points, xc1 and xc2.  There is a hysteresis loop for x > xc1: for two cells not in contact, the force is zero for x > xc1. \n\
    After contact is made the force is non-zero until x > xc2 - this is the effect of cell-cell adhesion."},

    {"C_FORCE", 0.2, 0, 0,
    "Attraction force factor 'c'",
    "The cell-cell force is a function of x = distance/(sum of radii): x = d/(R1+R2). \n\
    The force function is: F(x) = a/((x-x0)(x1-x)) + b, where x0 and x1 are the locations of the bounding asymptotes. \n\
    The parameter 'b' is calculated by setting the minimum value of F(x) (which occurs at x = (x0+x1)/2 ) equal to -c, this is the maximum attraction force. \n\
    The function F(x) is zero at two points, xc1 and xc2.  There is a hysteresis loop for x > xc1: for two cells not in contact, the force is zero for x > xc1. \n\
    After contact is made the force is non-zero until x > xc2 - this is the effect of cell-cell adhesion."},

    {"X0_FORCE", 0.7, 0, 0,
    "Left asymptote 'x0'",
    "The cell-cell force is a function of x = distance/(sum of radii): x = d/(R1+R2). \n\
    The force function is: F(x) = a/((x-x0)(x1-x)) + b, where x0 and x1 are the locations of the bounding asymptotes. \n\
    The parameter 'b' is calculated by setting the minimum value of F(x) (which occurs at x = (x0+x1)/2 ) equal to -c, this is the maximum attraction force. \n\
    The function F(x) is zero at two points, xc1 and xc2.  There is a hysteresis loop for x > xc1: for two cells not in contact, the force is zero for x > xc1. \n\
    After contact is made the force is non-zero until x > xc2 - this is the effect of cell-cell adhesion."},

    {"X1_FORCE", 1.3, 0, 0,
    "Right asymptote 'x1'",
    "The cell-cell force is a function of x = distance/(sum of radii): x = d/(R1+R2). \n\
    The force function is: F(x) = a/((x-x0)(x1-x)) + b, where x0 and x1 are the locations of the bounding asymptotes. \n\
    The parameter 'b' is calculated by setting the minimum value of F(x) (which occurs at x = (x0+x1)/2 ) equal to -c, this is the maximum attraction force. \n\
    The function F(x) is zero at two points, xc1 and xc2.  There is a hysteresis loop for x > xc1: for two cells not in contact, the force is zero for x > xc1. \n\
    After contact is made the force is non-zero until x > xc2 - this is the effect of cell-cell adhesion."},

    {"KDRAG", 5, 0, 0,
     "Drag factor",
     "Displacement = dt*F/drag"},

    {"FRANDOM",0.1, 0, 0,
     "Random force factor",
     "Magnitude of random additive force"},

//{"KDRAG", 1, 0, 0,
//"Drag weight: Kdrag",
//"Weighting for drag vs. chemotaxis"},

//{"KADHESION", 0.1, 0, 0,
//"Adhesion weight: Kadh",
//"Weighting for site attractiveness based on proximity"},

//{"KSTAY", 0.2, 0, 0,
//"Stay weight: Kstay",
//"Weighting for extra attractiveness of current site - i.e. immobility"},

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
