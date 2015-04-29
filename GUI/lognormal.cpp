    //-----------------------------------------------------------------------------------------
    // from: http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html
    // Implements the Gauss error function.
    //   erf(z) = 2 / sqrt(pi) * integral(exp(-t*t), t = 0..z)
    //
    // fractional error in math formula less than 1.2 * 10 ^ -7.
    // although subject to catastrophic cancellation when z in very close to 0
    // from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
    //-----------------------------------------------------------------------------------------
#include <math.h>
       double erf(double z)
	{
       double t = 1.0 / (1.0 + 0.5 * fabs(z));
       // use Horner's method
       double ans = 1 - t * exp( -z*z -  1.26551223 +
             t * ( 1.00002368 +
             t * ( 0.37409196 +
             t * ( 0.09678418 +
             t * (-0.18628806 +
             t * ( 0.27886807 +
             t * (-1.13520398 +
             t * ( 1.48851587 +
             t * (-0.82215223 +
             t * ( 0.17087277))))))))));
       if (z >= 0.0)
           return ans;
       else
           return -ans;
	}

	
    //-----------------------------------------------------------------------------------------
    // When X is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
    //-----------------------------------------------------------------------------------------
    double pnorm(double x1, double x2, double mu, double sig)
	{
        double z1, z2, e1, e2;
			
		z1 = (x1-mu)/sig;
        e1 = erf(z1/sqrt(2.0))/2;
        z2 = (x2-mu)/sig;
        e2 = erf(z2/sqrt(2.0))/2;
        return e2 - e1;
	}

	
    //-----------------------------------------------------------------------------------------
    // When log(X) is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
    //-----------------------------------------------------------------------------------------
    double plognorm(double x1, double x2, double mu, double sig)
	{
        double z1, z2, e1, e2;

        z1 = 0;
        z2 = 0;
        if (x1 == 0)
            e1 = -0.5;
		else {
            z1 = (log(x1)-mu)/sig;
            e1 = erf(z1/sqrt(2.0))/2;
		}
        if (x2 == 0)
            e2 = -0.5;
		else {
            z2 = (log(x2)-mu)/sig;
            e2 = erf(z2/sqrt(2.0))/2;
		}
        return e2 - e1;
	}


    //-----------------------------------------------------------------------------------------
    // Create the lognormal distribution with median = p1, shape = p2
    // at n points stored in x[], probability values stored in prob[].
    // Note that x[0] = 0.
    // The range of x is currently just less than 4*median.  This should be
    // OK for values of shape < 2.
    // Convert probability into probability density
    //-----------------------------------------------------------------------------------------
    void create_lognorm_dist(double p1, double p2,int n, double *x, double *prob)
	{
		double xmax, dx, mu_l, sig_l, x1, x2;

        if (p1 >= 0.5)
            xmax = p1*4;
        else
            xmax = p1*8;
            
        dx = xmax/n;
        mu_l = log(p1);
        sig_l = log(p2);
		for (int ix=0; ix<n; ix++) {
            x1 = (ix - 0.5)*dx;
            x2 = x1 + dx;
            x[ix] = (x1+x2)/2;
            prob[ix] = plognorm(x1,x2,mu_l,sig_l)/(x2-x1);
		}
	}

	
