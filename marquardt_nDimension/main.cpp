#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <string.h>


#include "statistics.h"
#include "furtherMathFunctions.h"
#include "commonConstants.h"
#include "functions.h"

using namespace std;




#define MAX_LINE_LENGTH 1000



// Function to find the minimum and maximum using the bisection method
void find_min_max(double (*func)(double*,double*,int,int),double* parameters,int xDim, int nrParameters,double a, double b, double *min, double *max)
{
    double epsilon = 1e-6; // Tolerance for precision
    double aStart = a;
    double bStart = b;
    // Bisection to find the minimum
    while (fabs(b - a) > epsilon) {
        double c = (a + b) / 2.0;
        double left = (a + c) / 2.0;
        double right = (c + b) / 2.0;

        if (func(&left,parameters,xDim,nrParameters) < func(&right,parameters,xDim,nrParameters))
            b = c;
        else
            a = c;
    }
    *min = (a + b) / 2.0;

    // Bisection to find the maximum
    //a = *min - 1.0;
    //b = *min + 1.0;
    a = aStart;
    b = bStart;
    while (fabs(b - a) > epsilon) {
        double c = (a + b) / 2.0;
        double left = (a + c) / 2.0;
        double right = (c + b) / 2.0;

        if (func(&left,parameters,xDim,nrParameters) > func(&right,parameters,xDim,nrParameters))
            b = c;
        else
            a = c;
    }
    *max = (a + b) / 2.0;
}


typedef struct {
    double field1;
    double field2;
    double field13;
    double field14;
} Fields;

typedef struct {
    Fields *fields;
    int line_count;
} CSVData;

CSVData* readCSV(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file %s.\n", filename);
        return NULL;
    }

    char line[MAX_LINE_LENGTH];
    int line_count = 0;

    // Skip the header
    if (fgets(line, MAX_LINE_LENGTH, file) == NULL) {
        fprintf(stderr, "Error reading the header.\n");
        fclose(file);
        return NULL;
    }

    // Allocate memory for CSVData
    CSVData *csv_data = (CSVData *)malloc(sizeof(CSVData));
    if (csv_data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        fclose(file);
        return NULL;
    }

    // Count the number of lines in the file
    while (fgets(line, MAX_LINE_LENGTH, file) != NULL) {
        line_count++;
    }

    // Reset file pointer to read from the beginning
    fseek(file, 0, SEEK_SET);

    // Allocate memory for Fields
    csv_data->fields = (Fields *)malloc(line_count * sizeof(Fields));
    if (csv_data->fields == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        fclose(file);
        free(csv_data);
        return NULL;
    }

    int current_line = 0;

    // Skip the header
    fgets(line, MAX_LINE_LENGTH, file);

    // Read and save fields 13 and 14
    while (fgets(line, MAX_LINE_LENGTH, file) != NULL && current_line < line_count) {
        char *token = strtok(line, ",");
        int field_count = 1;

        while (token != NULL) {
            if (field_count == 13) {
                csv_data->fields[current_line].field13 = strtod(token, NULL);
            } else if (field_count == 14) {
                csv_data->fields[current_line].field14 = strtod(token, NULL);
            }
            token = strtok(NULL, ",");
            field_count++;
        }

        current_line++;
    }

    // Save the line count
    csv_data->line_count = line_count;

    fclose(file);
    return csv_data;
}

void freeCSVData(CSVData *csv_data) {
    if (csv_data != NULL) {
        free(csv_data->fields);
        free(csv_data);
    }
}

int main()
{
    //const char *filename = "../dataset/data2.csv";
    const char *filename = "../dataset/erg5_dailytmin_20230201.csv";
    //const char *filename = "../dataset/erg5_dailytmin_20230724.csv";
    //const char *filename = "../dataset/erg5_dailytmax_20230201.csv";
    //const char *filename = "../dataset/erg5_dailytmax_20230724.csv";
    CSVData *csv_data = readCSV(filename);

    int nrParameters0 = 5; // to be parameterized
    int nrParameters1 = 1; // to be parameterized
    int maxIterationsNr = 10000; // to be parameterized
    int nrMinima = 5; // to be parameterized
    int nrPredictors = 1;
    double myEpsilon = EPSILON;
    std::vector <int> nrParameters(nrPredictors);


    nrParameters[0] = nrParameters0;
    //nrParameters[1] = nrParameters1;
    std::vector <std::vector <double>> parametersMin(nrPredictors);
    std::vector <std::vector <double>> parametersMax(nrPredictors);
    std::vector <std::vector <double>> parameters(nrPredictors);
    std::vector <std::vector <double>> parametersDelta(nrPredictors);


    int nrSteps=0;
    int nrData;
    if (csv_data != NULL)
    {
        printf("Number of lines: %d\n", csv_data->line_count);

        // Access and print values of fields 13 and 14 for each line
        /*for (int i = 0; i < csv_data->line_count; i++)
        {
            printf("Line %d - Field 13: %.2lf, Field 14: %.2lf\n", i+1, csv_data->fields[i].field13, csv_data->fields[i].field14);
        }*/
        nrData = csv_data->line_count;
        //nrData = 4;
        std::vector<double> value;
        value.resize(nrData);
        std::vector<double> weights;
        weights.resize(nrData);
        std::vector <std::vector <double>> predictors(nrData);
        for (int i=0;i<nrData;i++)
        {
            predictors[i].resize(nrPredictors);
        }
        for (int i=0;i<nrData;i++)
        {
            weights[i] = 1.0;
        }

        // prova doppia lineare

        /*predictors[0][0] = 0;
        predictors[1][0] = 1;
        predictors[2][0] = 2;
        predictors[3][0] = 3;
        value[0] = 0;
        value[1] = 1.;
        value[2] = 2.;
        value[3] = 3.;
        for (int i=0;i<nrData;i++)
        {
            weights[i] = 1.0;
        }
        weights[3] = 0.0001;*/

        //predictors[0][1] = 0;
        //predictors[1][1] = -0.5;
        //predictors[2][1] = 1;
        //predictors[3][1] = 1.2;


        for (int i=0;i<nrData;i++)
        {
            //predictors[i][0] = csv_data->fields[i].field13;
            //value[i] = csv_data->fields[i].field14;
            predictors[i][0] = csv_data->fields[i].field13;
            value[i] = csv_data->fields[i].field14;
            //predictors[i][1] = 0.1*i;
        }

            // Free allocated memory
        freeCSVData(csv_data);
        parameters[0].resize(nrParameters0);
        //parameters[1].resize(nrParameters1);
        parametersMin[0].resize(nrParameters0);
        //parametersMin[1].resize(nrParameters1);
        parametersMax[0].resize(nrParameters0);
        //parametersMax[1].resize(nrParameters1);
        parametersDelta[0].resize(nrParameters0);
        //parametersDelta[1].resize(nrParameters1);*/



        //parameters[0][0] = 1.;
        //parameters[1][0] = 2.;

        std::vector <double> xx(nrData);
        xx[0] = 0;
        xx[1] = 1;
        xx[2] = 2;
        xx[3] = 5;
        //parametrizzazione per spezzata

        parametersMin[0][0] = -100;
        parametersMax[0][0]= 1500;
        parametersMin[0][1]= -40;
        parametersMax[0][1]= 55;
        parametersMin[0][2]= -0.1;
        parametersMax[0][2]= 1500;
        parametersMin[0][3]= -40;
        parametersMax[0][3]= 55;
        parametersMin[0][4]= -0.05;
        parametersMax[0][4]= 0.001;

 /*       parametersMin[0][0] = -100;
        parametersMax[0][0] = 100;
*/

        //parametersMin[1][0] = -1000;
        //parametersMax[1][0] = 1000;

/*
        //parametrizzazione per Frei
        parametersMin[0]= -50;
        parametersMax[0]= 50;
        parametersMin[1]= -0.05;
        parametersMax[1]= 0.001;
        parametersMin[2]= -1;
        parametersMax[2]= 30;
        parametersMin[3]= -20;
        parametersMax[3]= 5000;
        parametersMin[4]= -20;
        parametersMax[4]= 5000;
*/
        //van Genuchten
/*        parametersMin[0][0] = 0.1;
        parametersMax[0][0] = 0.6;
        parametersMin[0][1] = 0.09;
        parametersMax[0][1] = 0.59;
        parametersMin[0][2] = 5;
        parametersMax[0][2] = 50;
        parametersMin[0][3] = 0.002;
        parametersMax[0][3] = 0.5;
        parametersMin[0][4] = 1.0001;
        parametersMax[0][4] = 1.7;
*/
        /*for (int i=0;i<nrPredictors;i++)
        {
            for (int j=0;j<nrParameters[i];j++)
            {
                parameters[i][j]= 0.5*(parametersMax[i][j]+parametersMin[i][j]);
                parametersDelta[i][j]= 0.00001;
            }
        }*/

        clock_t startTime = clock();
        //nrSteps = interpolation::bestFittingMarquardt_nDimension(&tempVsHeightPiecewise,maxIterationsNr,nrMinima,parametersMin,parametersMax,parameters,nrParameters,parametersDelta,maxIterationsNr,myEpsilon,height,value,nrData,1,false,weights);
        //const double RATIO_DELTA = 1000;

        std::vector<std::function<double(double, std::vector<double>&)>> myFunc;

        myFunc.push_back(lapseRatePiecewise_three);
        //myFunc.push_back(modifiedVanGenuchtenRestricted_nDim);
        //myFunc.push_back(functionLinear);
        //double result = functionSum(myFunc,xx,parameters);
        //printf("risultato %f \n",result);
        //nrSteps = interpolation::bestFittingMarquardt_nDimension(&functionSum, myFunc, 10000, nrMinima, parametersMin, parametersMax, parameters, parametersDelta,
                                        //100, EPSILON, 0.01, predictors, value, false, weights);
        nrSteps = interpolation::bestFittingMarquardt_nDimension(&functionSum, myFunc, 10000, nrMinima, parametersMin, parametersMax, parameters, parametersDelta,
                                        100, EPSILON, 0.01, predictors, value, weights);

        clock_t endTime = clock();
        double deltaTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
        printf("Tempo impiegato: %f secondi\n", deltaTime);
        printf("nrsteps %d\n",nrSteps);

        for (int i=0;i<nrPredictors;i++)
        {
            for (int j=0;j<nrParameters[i];j++)
            {
                std::cout << parameters[i][j] << std::endl;
            }
            printf("\n");
        }


        printf("\n");
        //double result = functionSum(myFunc,xx,parameters);
        //printf("risultatoR2Weights %f \n",interpolation::computeWeighted_R2(value,xx,weights));
        //printf("risultatoR2 %f \n",interpolation::computeR2(value,xx));
        //double valueFunc;
        //startTime = clock();
        //for (int i=0;i<100;i++)
        //{
            //std::vector <double> xx(1);
            //xx[0] = -10. + i*1;
            //valueFunc = tempVsHeightPiecewise(&xx,parameters,1,5);
            //valueFunc = lapseRatePiecewise(xx,parameters);
            //printf("%.1f\t%.1f\n",xx[0],valueFunc);
        //}
        //endTime = clock();
        //deltaTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
        //printf("Tempo impiegato: %f secondi\n", deltaTime);


        //double minimum,maximum;
        //find_min_max(&tempVsHeightPiecewise,parameters,1,5,-10,2000,&minimum,&maximum);
        //printf("%.2f\t%.2f\n",minimum,maximum);
        //double x1,y1,x2,y2;
        //x1 = parameters[0];
        //y1 = parameters[1];
        //x2 = x1 + parameters[2];
        //y2 = parameters[3];
        //printf("%.2f\t%.2f\t%.2f\t%.2f\n",x1,y1,x2,y2);
        printf("\n");


    }

    nrPredictors = 2;
    nrData = 5;
    std::vector<float> value2;
    value2.resize(nrData);
    std::vector<float> weights2;
    weights2.resize(nrData);
    std::vector<float> slope;
    //slope.resize(nrPredictors);
    float q;
    std::vector <std::vector <float>> predictors2(nrData);
    for (int i=0;i<nrData;i++)
    {
        predictors2[i].resize(nrPredictors);
    }
    for (int i=0;i<nrData;i++)
    {
        weights2[i] = 1.0;
    }


    // prova doppia lineare
    /*
    predictors2[0][0] = 0;
    predictors2[1][0] = 1;
    predictors2[2][0] = 2;
    predictors2[3][0] = 3;

    predictors2[0][1] = 0;
    predictors2[1][1] = -0.5;
    predictors2[2][1] = 1;
    predictors2[3][1] = 1.2;


    value2[0] = 0;
    value2[1] = 1.;
    value2[2] = 2.;
    value2[3] = 3.+6;
    */
    predictors2[0][0] = 2;
    predictors2[1][0] = 3;
    predictors2[2][0] = 5;
    predictors2[3][0] = 7;
    predictors2[4][0] = 8;

    predictors2[0][1] = 1;
    predictors2[1][1] = 5;
    predictors2[2][1] = 3;
    predictors2[3][1] = 6;
    predictors2[4][1] = 7;


    value2[0] = 3;
    value2[1] = 2.;
    value2[2] = 4.;
    value2[3] = 5;
    value2[4] = 8;

    for (int i=0;i<nrData;i++)
    {
        weights2[i] = 1.0;
    }
    //weights2[3] = 0.01;

    statistics::weightedMultiRegressionLinear(predictors2,value2,weights2,nrData,&q,slope,nrPredictors);
    printf("linear regression %f\t%f\n",slope[0],q);
    float qSE,R2,stdError;
    std::vector<float> slopeSE;
    //slopeSE.resize(nrPredictors);
    // prova doppia lineare
    /*
    predictors2[0][0] = 0;
    predictors2[1][0] = 1;
    predictors2[2][0] = 2;
    predictors2[3][0] = 3;
    predictors2[0][1] = 0;
    predictors2[1][1] = -0.5;
    predictors2[2][1] = 1;
    predictors2[3][1] = 1.2;

    value2[0] = 0;
    value2[1] = 1.;
    value2[2] = 2.;
    value2[3] = 3.+6;
    for (int i=0;i<nrData;i++)
    {
        weights2[i] = 1.0;
    }
    weights2[3] = 0.001;
    */
    predictors2[0][0] = 2;
    predictors2[1][0] = 3;
    predictors2[2][0] = 5;
    predictors2[3][0] = 7;
    predictors2[4][0] = 8;

    predictors2[0][1] = 1;
    predictors2[1][1] = 5;
    predictors2[2][1] = 3;
    predictors2[3][1] = 6;
    predictors2[4][1] = 7;


    value2[0] = 3;
    value2[1] = 2.;
    value2[2] = 4.;
    value2[3] = 5;
    value2[4] = 8;

    std::vector<float> slope2;
    float q2;
    statistics::weightedMultiRegressionLinearWithStats(predictors2,value2,weights2,&q2,slope2,true,true,&R2,&stdError,&qSE,slopeSE);
    printf("linear regression with stats %f\t%f\n",slope[1],q);
    printf("R2, m, q %f\t%f\t%f\n",R2,slopeSE[1],qSE);




    return 0;
}




