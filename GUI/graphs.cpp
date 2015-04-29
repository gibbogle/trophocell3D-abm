//ABM
#include <qstring.h>
#include "graphs.h"
#include "log.h"

LOG_USE();

Graphs::Graphs()
{
    GRAPH_SET tsGraphSet[] = {

{"NTcells",
"Number of Cells",
"",  //No. of cells ",
2, true, 0, 1, 0, true},

{"CD69",
"CD69 Profile",
"",  //Fraction",
PROFILE_CD69, false, 0, 1, 1.0, false}

};
    // Note: tsGraphs[] is constant = tsGraphSet[]
    // .active is either true or false, depending on checkBox_selectgraph
    // Note: the "ts" prefix is misleading because the list can include profile plots.

    n_tsGraphs = sizeof(tsGraphSet)/sizeof(GRAPH_SET);
    tsGraphs = new GRAPH_SET[n_tsGraphs];
    for (int i=0; i<n_tsGraphs; i++) {
        tsGraphs[i] = tsGraphSet[i];
    }
    graphList = new GRAPH_SET[maxGraphs];
    nGraphs = maxGraphs;
}


GRAPH_SET Graphs::get_graph(int k)
{
	return graphList[k];
}

int Graphs::get_dataIndex(int k)
{
	return graphList[k].dataIndex;
}

QString Graphs::get_tag(int k)
{
	return graphList[k].tag;
}

QString Graphs::get_title(int k)
{
//    QString title;
//    if (!graphList[k].active)
//        title = "";
//    else
//        title = graphList[k].title;
    return graphList[k].title;
}

QString Graphs::get_yAxisTitle(int k)
{
//    QString title;
//    if (!graphList[k].active)
//        title = "";
//    else
//        title = graphList[k].yAxisTitle;
    return graphList[k].yAxisTitle;
}

double Graphs::get_maxValue(int k) {
	return graphList[k].maxValue;
}

double Graphs::get_scaling(int k) {
	return graphList[k].scaling;
}

double Graphs::get_yscale(int k) {
    return graphList[k].yscale;
}

double Graphs::get_xscale(double xmax) {
    int n = 1;
    for (;;) {
        if (xmax <= n) break;
        n++;
    }
    return double(n);
}

void Graphs::set_activity(int k, bool activity)
{
    graphList[k].active = activity;
}

bool Graphs::isActive(int k)
{
	return graphList[k].active;
}

bool Graphs::isTimeseries(int k)
{
    return graphList[k].ts;
}

bool Graphs::isProfile(int k)
{
    return !graphList[k].ts;
}

void Graphs::set_maxValue(int k, double v)
{
	graphList[k].maxValue = v;
}

void Graphs::makeGraphList()
{
    char msg[128];
    int k = maxGraphs;
    int nts = 0;
    for (int i=0; i<n_tsGraphs; i++) {
        if (tsGraphs[i].active) {
            k--;
            graphList[k] = tsGraphs[i];
            nts++;
            if (nts == maxGraphs) break;
        }
    }
    int ndummy = maxGraphs - nts;
    for (k=0; k<ndummy; k++) {
        graphList[k].active = false;
        graphList[k].ts = true;
        graphList[k].tag = "dummy";
        graphList[k].scaling = 1;
        graphList[k].title = "";
        graphList[k].yAxisTitle = "";
    }
    nGraphs = maxGraphs;
//    sprintf(msg,"nGraphs: %d",nGraphs);
//    LOG_MSG(msg);
//    for (k=0; k<nGraphs; k++) {
//        LOG_QMSG(graphList[k].tag);
//        sprintf(msg,"k: %d  active: %d  ts: %d",k,graphList[k].active,graphList[k].ts);
//        LOG_MSG(msg);
//    }
}
