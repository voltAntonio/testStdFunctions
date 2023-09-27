#include <stdio.h>
#include <math.h>

// Function to calculate weighted variance
double weighted_variance(double *data, double *weights, int n)
{
    if (n <= 0) {
        // Handle the case when there is no data or weights
        return 0.0;
    }

    double sum_weights = 0.0;
    double sum_weighted_data = 0.0;
    double sum_squared_weighted_data = 0.0;

    // Calculate the necessary sums for weighted variance calculation
    for (int i = 0; i < n; i++)
    {
        sum_weights += weights[i];
        sum_weighted_data += data[i] * weights[i];
        sum_squared_weighted_data += data[i] * data[i] * weights[i];
    }

    // Calculate the weighted variance
    double weighted_mean = sum_weighted_data / sum_weights;
    double variance = (sum_squared_weighted_data / sum_weights) - (weighted_mean * weighted_mean);

    return variance;
}

// Function to calculate the weighted R-squared (coefficient of determination)
double calculate_weighted_r_squared(double *observed, double *predicted, double *weights, int n)
{
    double sum_weighted_squared_residuals = 0.0;
    double sum_weighted_squared_total = 0.0;
    double weighted_mean_observed = 0.0;

    // Calculate the weighted mean of the observed values
    double sum_weights = 0.0;
    for (int i = 0; i < n; i++)
    {
        weighted_mean_observed += observed[i] * weights[i];
        sum_weights += weights[i];
    }
    weighted_mean_observed /= sum_weights;

    // Calculate the sums needed for weighted R-squared calculation
    for (int i = 0; i < n; i++)
    {
        double weighted_residual = weights[i] * (observed[i] - predicted[i]);
        sum_weighted_squared_residuals += weighted_residual * weighted_residual;

        double weighted_total_deviation = weights[i] * (observed[i] - weighted_mean_observed);
        sum_weighted_squared_total += weighted_total_deviation * weighted_total_deviation;
    }

    // Calculate weighted R-squared
    double weighted_r_squared = 1.0 - (sum_weighted_squared_residuals / sum_weighted_squared_total);

    return weighted_r_squared;
}


double tempVsHeightSigmoidal(double* x, double* par, int xDim, int nrPar)
{
    double y;
    y = par[0] + par[1]*x[0] + par[2]*(1/(1+exp(-par[3]*(x[0] - par[4]))));
    return y;
}

double tempVsHeightFrei(double* x, double* par, int xDim, int nrPar)
{
    /*
    par[0] = T0;
    par[1] = gamma;
    par[2] = a;
    par[3] = h0;
    par[4] = h1;
    */
    double y;
    y = par[0] - par[1]*x[0];
    if (x[0] <= par[3])
    {
       return y - par[2];
    }
    else if (x[0] >= par[4])
    {
        return y;
    }
    double PI = 3.1415926536;
    return y - 0.5*par[2]*(1 + cos(PI*(x[0]-par[3])/(par[4]-par[3])));
}

double tempVsHeightPiecewise(double* x, double* par, int xDim, int nrPar)
{
    //double y,m,q;
    //double y,q;
    double xb;
    // par[2] means the delta between the two quotes. It must be positive.
    xb = par[0]+par[2];
    if (x[0]>xb)
    {
        //m = par[4];
        //q = par[3]-par[4]*xb;
        return par[4]*x[0] + par[3]-par[4]*xb;
    }
    else if (x[0] < par[0])
    {
        //m = par[4];
        //q = par[1]-par[4]*par[0];
        return par[4]*x[0] + par[1]-par[4]*par[0];
    }
    else
    {
        //m = (par[3]-par[1])/par[2];
        //q = par[1]-(par[3]-par[1])/par[2]*par[0];
        return ((par[3]-par[1])/par[2])*x[0]+ par[1]-(par[3]-par[1])/par[2]*par[0];
    }
    //y = m*x[0]+q;
    //return y;
}

double functionMualem(double* h, double* par, int xDim, int nrPar)
{
    double k;
    k = (par[0]*pow((1-pow(par[1]*h[0],par[2]*par[3])*pow(1+pow(par[1]*h[0],par[3]),-par[2])),2))/pow(1+pow(par[1]*h[0],par[3]),par[2]*par[4]);
    return k;
}

double modifiedVanGenuchtenRestricted_nDimensional(double* psi, double *parameters, int nrPsi, int nrParameters)
{
    *psi = fabs(*psi);
    double thetaS, thetaR, he,deltaTheta;
    double alpha, n, m;

    thetaS = parameters[0];         // water content at saturation [m^3 m^-3]
    deltaTheta = parameters[1];
    thetaR = thetaS - deltaTheta;         // water content residual [m^3 m^-3]
    he = parameters[2];             // air entry [kPa]

    if (*psi <= he) return thetaS;

    alpha = parameters[3];          // Van Genuchten curve parameter [kPa^-1]
    n = parameters[4];              // Van Genuchten curve parameter [-]

    m = 1 - 1/n;                // Van Genuchten curve parameter (restricted: 1-1/n) [-]

    // reduction factor for modified VG (Ippisch, 2006) [-]
    double sc = pow(1 + pow(alpha * he, n), -m);

    // degree of saturation [-]
    double Se = pow(1 + pow(alpha * (*psi), n), -m) / sc;

    // volumetric water content [m^3 m^-3]
    return Se * (deltaTheta) + thetaR;
}

double modifiedVanGenuchtenNotRestricted_nDimensional(double* psi, double *parameters, int nrPsi, int nrParameters)
{
    *psi = fabs(*psi);
    double thetaS, thetaR, he;
    double alpha, n, m;

    thetaS = parameters[0];         // water content at saturation [m^3 m^-3]
    thetaR = parameters[1];         // water content residual [m^3 m^-3]
    he = parameters[2];             // air entry [kPa]

    if (*psi <= he) return thetaS;

    alpha = parameters[3];          // Van Genuchten curve parameter [kPa^-1]
    n = parameters[4];              // Van Genuchten curve parameter [-]

    m = parameters[5];

    // reduction factor for modified VG (Ippisch, 2006) [-]
    double sc = pow(1 + pow(alpha * he, n), -m);

    // degree of saturation [-]
    double Se = pow(1 + pow(alpha * (*psi), n), -m) / sc;

    // volumetric water content [m^3 m^-3]
    return Se * (thetaS - thetaR) + thetaR;
}
