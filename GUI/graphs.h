#ifndef GRAPHS_H
#define GRAPHS_H

#include "profile.h"

#define maxGraphs 12     // This is the number of graphs that can be displayed (tiled)

struct graph_set {
	QString tag;
	QString title;
	QString yAxisTitle;
	int dataIndex;		// this must be consistent with the ordering of summary_data[]
	bool active;		// false for dummy graphs
	double maxValue;
	double scaling;		// multiplying factor for scaling of summary_data
    double yscale;      // if != 0, this is the yscale for the plot
    bool ts;
};

typedef graph_set GRAPH_SET;

class Graphs
{
	GRAPH_SET *graphList;

public:

	Graphs();
	~Graphs();
    GRAPH_SET *tsGraphs;
    GRAPH_SET get_graph(int);
    int n_tsGraphs;
    int nGraphs;
	int get_dataIndex(int);
	QString get_tag(int);
	QString get_title(int);
	QString get_yAxisTitle(int);
	double get_maxValue(int);
	double get_scaling(int);
    double get_yscale(int);
    double get_xscale(double xmax);
    void set_activity(int, bool);
    bool isActive(int);
    bool isTimeseries(int);
    bool isProfile(int);
    void set_maxValue(int, double);
//    void makeGraphList(int);
    void makeGraphList();

};

#endif
