#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double functionMualem(double* h, double* par, int xDim, int nrPar);
double modifiedVanGenuchtenRestricted_nDimensional(double* psi, double *parameters, int nrPsi, int nrParameters);
double modifiedVanGenuchtenNotRestricted_nDimensional(double* psi, double *parameters, int nrPsi, int nrParameters);

double tempVsHeightSigmoidal(double* x, double* par, int xDim, int nrPar);
double tempVsHeightFrei(double* x, double* par, int xDim, int nrPar);
double tempVsHeightPiecewise(double* x, double* par, int xDim, int nrPar);

#endif // FUNCTIONS_H
