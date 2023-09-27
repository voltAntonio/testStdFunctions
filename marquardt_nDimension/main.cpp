#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <functional>

#include "furtherMathFunctions.h"
#include "commonConstants.h"
//#include "inputOutput.h"
//#include "functions.h"

using namespace std;

//#include <stdio.h>
//#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1000

//#include <stdio.h>
//#include <math.h>



// Function to find the minimum and maximum using the bisection method
void find_min_max(double (*func)(double*,double*,int,int),double* parameters,int xDim, int nrParameters,double a, double b, double *min, double *max) {
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
    const char *filename = "../dataset/erg5_dailytmin_20230201.csv";
    //const char *filename = "../dataset/erg5_dailytmin_20230724.csv";
    //const char *filename = "../dataset/erg5_dailytmax_20230201.csv";
    //const char *filename = "../dataset/erg5_dailytmax_20230724.csv";
    CSVData *csv_data = readCSV(filename);
    double** height;
    double* value;
    double* weights;
    int nrParameters =5; // to be parameterized
    int maxIterationsNr = 10000; // to be parameterized
    int nrMinima = 5; // to be parameterized
    double myEpsilon = EPSILON;
    std::vector<double> parametersMin(nrParameters);
    std::vector<double> parametersMax(nrParameters);
    std::vector<double> parameters(nrParameters);
    std::vector<double> parametersDelta(nrParameters);
    int nrSteps=0;
    int nrData;
    if (csv_data != NULL)
    {
        printf("Number of lines: %d\n", csv_data->line_count);

        // Access and print values of fields 13 and 14 for each line
        for (int i = 0; i < csv_data->line_count; i++)
        {
            //printf("Line %d - Field 13: %.2lf, Field 14: %.2lf\n", i+1, csv_data->fields[i].field13, csv_data->fields[i].field14);
        }
        nrData = csv_data->line_count;
        //printf("Lines nr %d\n",csv_data->line_count);
        std::vector<double> value;
        value.resize(nrData);
        std::vector<double> weights;
        weights.resize(nrData);
        std::vector <std::vector <double>> height(nrData) ;
        for (int i=0;i<nrData;i++)
        {
            height[i].resize(1);
        }
        for (int i=0;i<nrData;i++)
        {
            height[i][0] = csv_data->fields[i].field13;
            value[i] = csv_data->fields[i].field14;
        }
            // Free allocated memory
        freeCSVData(csv_data);


        //parametrizzazione per spezzata
        parametersMin[0]= -0;
        parametersMax[0]= 1500;
        parametersMin[1]= -40;
        parametersMax[1]= 55;
        parametersMin[2]= -0.1;
        parametersMax[2]= 1500;
        parametersMin[3]= -40;
        parametersMax[3]= 55;
        parametersMin[4]= -0.05;
        parametersMax[4]= 0.001;

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

        for (int i=0;i<nrParameters;i++)
        {
            //parametersMin[i]= -1000;
            //parametersMax[i]= 1000;
            parameters[i]= 0.5*(parametersMax[i]+parametersMin[i]);
            parametersDelta[i] = 0.00001;
        }

        clock_t startTime = clock();
        //nrSteps = interpolation::bestFittingMarquardt_nDimension(&tempVsHeightPiecewise,maxIterationsNr,nrMinima,parametersMin,parametersMax,parameters,nrParameters,parametersDelta,maxIterationsNr,myEpsilon,height,value,nrData,1,false,weights);
        //const double RATIO_DELTA = 1000;

        std::vector<std::function<double(std::vector<double>&, std::vector<double>&)>> myFunc;
        //myFunc[0] = lapseRatePiecewise;
        myFunc.push_back(lapseRatePiecewise);
        nrSteps = interpolation::bestFittingMarquardt_nDimension(&functionSum, myFunc, 10000, 5, parametersMin, parametersMax, parameters, parametersDelta,
                                        100, EPSILON, 0.01, height, value, value.size(), 1, false, weights);
        clock_t endTime = clock();
        double deltaTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
        printf("Tempo impiegato: %f secondi\n", deltaTime);
        printf("nrsteps %d\n",nrSteps);

        for (int i=0;i<nrParameters;i++)
        {
            printf("%f\t",parameters[i]);
        }
        printf("\n");
        double valueFunc;
        startTime = clock();
        for (int i=0;i<1000000;i++)
        {
            std::vector <double> xx(1);
            xx[0] = -10. + i*1;
            //valueFunc = tempVsHeightPiecewise(&xx,parameters,1,5);
            valueFunc = lapseRatePiecewise(xx,parameters);
            //printf("%.1f\t%.1f\n",xx[0],valueFunc);
        }
        endTime = clock();
        deltaTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
        printf("Tempo impiegato: %f secondi\n", deltaTime);


        //double minimum,maximum;
        //find_min_max(&tempVsHeightPiecewise,parameters,1,5,-10,2000,&minimum,&maximum);
        //printf("%.2f\t%.2f\n",minimum,maximum);
        double x1,y1,x2,y2;
        x1 = parameters[0];
        y1 = parameters[1];
        x2 = x1 + parameters[2];
        y2 = parameters[3];
        //printf("%.2f\t%.2f\t%.2f\t%.2f\n",x1,y1,x2,y2);
        printf("\n");
        //free(height);
        //free(value);
        //free(weights);
        //free(parameters);
        //free(parametersDelta);
        //free(parametersMax);
        //free(parametersMin);
    }


    return 0;
}
