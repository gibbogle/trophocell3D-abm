#ifndef PARAMS_H
#define PARAMS_H

struct param_set {
	QString tag;
	double value;
	double minvalue;
	double maxvalue;
	QString label;
	QString text;
};

typedef param_set PARAM_SET;

class Params 
{
	PARAM_SET *workingParameterList;

public:

	Params();
	~Params();
	PARAM_SET get_param(int);
	int nParams;
	void set_value(int, double);
	void set_label(int, QString);
};

#endif
