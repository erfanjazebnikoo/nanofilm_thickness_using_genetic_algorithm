/* 
 * File:   GeneticAlgorithem.hpp
 * Author: Erfan Jazeb Nikoo
 *
 */

#ifndef GENETICALGORITHEM_HPP
#define	GENETICALGORITHEM_HPP

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <string>
#include <limits>
#include <wtypes.h>
#include <algorithm>
#include <vector>

#define POPSIZE         50
#define MAXGENS         1000
#define NVARS           12
#define PXOVER          0.8
#define PMUTATION       0.5
#define LANDASIZE       601

#define RANGE           0.1

#define PI		3.14159265358979323846

double nInf;
double Efb;
double A1;
double A2;
double B1;
double B2;
double C1;
double C2;
double Ed;
double Egamma0;
double s;
double d;

double landa[LANDASIZE];
double Tmeas[LANDASIZE];
double Tcalc[LANDASIZE];

static double varMin[NVARS] = {1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.5, 80.0};
static double varMax[NVARS] = {3.0, 4.5, 1.0, 0.01, 200.0, 150.0, 50.0, 20.0, 4.5, 1.0, 2, 300.0};

struct genotype
{
    double gene[NVARS];
    double fitness;
    double upper[NVARS];
    double lower[NVARS];
    double rfitness;
    double cfitness;
};

struct genotype population[POPSIZE + 1];
struct genotype newPopulation[POPSIZE * 2];
struct genotype tempPopulation[POPSIZE * 4];

std::vector<genotype> bestPopulation;

struct genotype newpopulation[POPSIZE + 1];

int main();
void evaluate();
void initialize(int &seed);
void report(int generation);
void Xover(int one, int two, int &seed);
void keepTheBest();

////////////////////////////////////////////// implement for NANO question
void excelReader();
void TcalcCalculator(int i);
double fitnessCalculator();

void crossover(int &seed);
void mutate(int &seed);
void newEvaluate();
void sortPopulation();
void selector();
bool sortByFitness(const genotype &a, const genotype &b);

int newPopulationSize;
int bestPopulationSize;

#endif	/* GENETICALGORITHEM_HPP */

