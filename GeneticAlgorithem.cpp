/* 
 * File:   GeneticAlgorithem.hpp
 * Author: Erfan Jazeb Nikoo
 *
 */

#include "GeneticAlgorithem.hpp"
#include "Utilities.hpp"
#include "ExcelReader.hpp"

using namespace std;

int main()
{
    timestamp();
    const clock_t beginTime = clock();
    cout << "\n";

    if (NVARS < 2)
    {
        cout << "\n";
        cout << "  The crossover modification will not be available,\n";
        cout << "  since it requires 2 <= NVARS.\n";
    }

    int seed = 13700207;

    initialize(seed);

    excelReader();

    evaluate();

    keepTheBest();

    for (int generation = 0; generation < MAXGENS; generation++)
    {
        crossover(seed);
        mutate(seed);
        newEvaluate();
        sortPopulation();
        selector();
        report(generation);
    }

    for (int i = 0; i < NVARS; i++)
    {
        cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
    }

    ofstream o("Tcalc.dat");

    for (int i = 0; i < LANDASIZE; i++)
    {
        o << landa[i] << " " << Tcalc[i] << endl;
    }

    cout << endl << "  Best fitness = " << population[0].fitness << endl;
    cout <<"Total Time: "<< float( clock() - beginTime) / CLOCKS_PER_SEC << endl;
    timestamp();

    return 0;
}

//****************************************************************************

void evaluate()
{

    for (int member = 0; member < POPSIZE; member++)
    {
        nInf = population[member].gene[0];
        Efb = population[member].gene[1];
        A1 = population[member].gene[2];
        A2 = population[member].gene[3];
        C1 = population[member].gene[4];
        C2 = population[member].gene[5];
        B1 = population[member].gene[6];
        B2 = population[member].gene[7];
        Ed = population[member].gene[8];
        Egamma0 = population[member].gene[9];
        s = population[member].gene[10];
        d = population[member].gene[11];

        population[member].fitness = fitnessCalculator();
    }
}

//****************************************************************************

void initialize(int &seed)
{
    double lbound;
    double ubound;
    int B1Pos = 6;
    int B2Pos = 7;
    int C1Pos = 4;
    int C2Pos = 5;
    double randNum;

    ///////////////////////////// new
    newPopulationSize = 0;
    bestPopulationSize = 0;
    /////////////////////////////

    // 
    //  Initialize variables within the bounds 
    //
    for (int i = 0; i < NVARS; i++)
    {
        lbound = varMin[i];
        ubound = varMax[i];
        for (int j = 0; j < POPSIZE; j++)
        {
            population[j].fitness = 0;
            population[j].rfitness = 0;
            population[j].cfitness = 0;
            population[j].lower[i] = lbound;
            population[j].upper[i] = ubound;
            randNum = r8_uniform_ab(lbound, ubound, seed);
            if (i == B1Pos)
            {
                while (4 * population[j].gene[C1Pos] <= pow(randNum, 2))
                {
                    randNum = r8_uniform_ab(lbound, ubound, seed);
                }
            }
            else if (i == B2Pos)
            {
                while (4 * population[j].gene[C2Pos] <= pow(randNum, 2))
                {
                    randNum = r8_uniform_ab(lbound, ubound, seed);
                }
            }
            population[j].gene[i] = randNum;
        }
    }
}

//****************************************************************************

void report(int generation)
{
    double avg;
    double best_val;
    int i;
    double square_sum;
    double stddev;
    double sum;
    double sum_square;

    //    if (generation % 10 == 0)
    cout << generation << ": " << population[0].fitness << endl;

    //    if (generation == 0)
    //    {
    //        cout << "\n";
    //        cout << "  Generation       Best            Average       Standard \n";
    //        cout << "  number           value           fitness       deviation \n";
    //        cout << "\n";
    //    }
    //
    //    sum = 0.0;
    //    sum_square = 0.0;
    //
    //    for (i = 0; i < POPSIZE; i++)
    //    {
    //        sum = sum + population[i].fitness;
    //        sum_square = sum_square + population[i].fitness * population[i].fitness;
    //    }
    //
    //    avg = sum / (double) POPSIZE;
    //    square_sum = avg * avg * POPSIZE;
    //    stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));
    //    best_val = population[POPSIZE].fitness;
    //
    //    cout << "  " << setw(8) << generation
    //            << "  " << setw(14) << best_val
    //            << "  " << setw(14) << avg
    //            << "  " << setw(14) << stddev << "\n";
}

//****************************************************************************

void Xover(int one, int two, int &seed)
{
    //    int point = NVARS; //i4_uniform_ab(0, NVARS - 1, seed);

    for (int i = 0; i < NVARS; i++)
    {
        swap(population[one].gene[i], population[two].gene[i]);
    }
}

//****************************************************************************

void keepTheBest()
{
    int cur_best = 0;

    for (int i = 0; i < POPSIZE; i++)
    {
        if (population[POPSIZE].fitness < population[i].fitness)
        {
            cur_best = i;
            population[POPSIZE].fitness = population[i].fitness;
        }
    }
    for (int j = 0; j < NVARS; j++)
    {
        population[POPSIZE].gene[j] = population[cur_best].gene[j];
    }
}

//****************************************************************************

void excelReader()
{
    ExcelReader excelReader("Landa.csv", 2);

    excelReader.readTable(0, landa);
    excelReader.readTable(1, Tmeas);

    ofstream o("Tmeas.dat");

    for (int i = 0; i < LANDASIZE; i++)
    {
        o << landa[i] << " " << Tmeas[i] << endl;
    }
}

//****************************************************************************

void TcalcCalculator(int i)
{
    const static double n0 = 1;
    const static double eInf = 5.3;
    const static double ns = 1.5;

    Tcalc[i] = 0;

    double E = (double) 1242.0 / landa[i];

    double Egamma = Egamma0 * pow((Ed / sqrt(eInf)), -s) * pow(E, s);

    double e1D = -1.0 * (Ed / (pow(E, 2) + pow(Egamma, 2)));

    double e2D = (pow(Ed, 2) * Egamma) / (pow(E, 3)+(E * pow(Egamma, 2)));

    double nD = sqrt((sqrt(pow(e1D, 2) + pow(e2D, 2)) + e1D) / 2.0);

    double kD = sqrt((sqrt(pow(e1D, 2) + pow(e2D, 2)) - e1D) / 2.0);

    double Q1 = 0.5 * sqrt((4.0 * C1) - pow(B1, 2));

    double Q2 = 0.5 * sqrt((4.0 * C2) - pow(B2, 2));

    double Bo1 = (A1 / Q1)*((Efb * B1) - pow(Efb, 2) + C1 - (pow(B1, 2) / 2.0));

    double Bo2 = (A2 / Q2)*((Efb * B2) - pow(Efb, 2) + C2 - (pow(B2, 2) / 2.0));

    double Co1 = (A1 / Q1)*(((pow(Efb, 2) + C1)*(B1 / 2.0))-(2.0 * Efb * C1));

    double Co2 = (A2 / Q2)*(((pow(Efb, 2) + C2)*(B2 / 2.0))-(2.0 * Efb * C2));

    double kFB = ((A1 / (pow(E, 2)-(B1 * E) + C1))+(A2 / (pow(E, 2)-(B2 * E) + C2))) * pow((E - Efb), 2);

    double nFB = nInf + ((((Bo1 * E) + Co1) / (pow(E, 2)-(B1 * E) + C1))+(((Bo2 * E) + Co2) / (pow(E, 2)-(B2 * E) + C2)));

    double e1 = (pow(nFB, 2) + pow(nD, 2))-(pow(kFB, 2) + pow(kD, 2));

    double e2 = 2.0 * ((nFB * kFB)+(nD * kD));

    double n = sqrt((sqrt(pow(e1, 2) + pow(e2, 2)) + e1) / 2.0);

    double k = sqrt((sqrt(pow(e1, 2) + pow(e2, 2)) - e1) / 2.0);

    double g1 = (pow(n0, 2) - pow(n, 2) - pow(k, 2)) / (pow((n0 + n), 2) + pow(k, 2));

    double g2 = (pow(n, 2) - pow(ns, 2) + pow(k, 2)) / (pow((n + ns), 2) + pow(k, 2));

    double h1 = (2.0 * n0 * k) / (pow((n0 + n), 2) + pow(k, 2));

    double h2 = (-2.0 * ns * k) / (pow((n + ns), 2) + pow(k, 2));

    double C = 2.0 * ((g1 * g2)-(h1 * h2));

    double D = 2.0 * ((g1 * h2)+(g2 * h1));

    double beta = (2.0 * PI * k * d) / landa[i];

    double gamma = (2.0 * PI * n * d) / landa[i];

    double A = (pow((1 + g1), 2) + pow(h1, 2))*(pow((1.0 + g2), 2) + pow(h2, 2));

    double B = exp(2.0 * beta)+((pow(g1, 2) + pow(h1, 2))*(pow(g2, 2) + pow(h2, 2)) * exp(-2.0 * beta))+(C * cos(2.0 * gamma))+(D * sin(2.0 * gamma));

    Tcalc[i] = ((ns * A) / (n0 * B))*((4 * ns * n0) / (pow((ns + n0), 2)));
}

//****************************************************************************

double fitnessCalculator()
{
    double sum = 0;
    for (int i = 0; i < LANDASIZE; i++)
    {
        TcalcCalculator(i);

        sum += pow((Tcalc[i] - Tmeas[i]), 2);
    }

    if (isnan(sum))
    {
        return 100000.0;
    }

    return sum;
}

//****************************************************************************

void crossover(int &seed)
{
    const double a = 0.0;
    const double b = 1.0;
    double x;
    double delta = 0.02;
    double alpha;
    double y;

    int divider = POPSIZE / 2;

    for (int i = 0; i < POPSIZE; i++)
    {
        x = r8_uniform_ab(a, b, seed);

        if (x < PXOVER)
        {
            newPopulation[newPopulationSize].fitness = 0;
            //            newPopulation[newPopulationSize].rfitness = 0;
            //            newPopulation[newPopulationSize].cfitness = 0;
            alpha = r8_uniform_ab(-1.0 * delta, 1.0 + delta, seed);
            for (int j = 0; j < NVARS; j++)
            {
                newPopulation[newPopulationSize].lower[j] = varMin[j];
                newPopulation[newPopulationSize].upper[j] = varMax[j];
                if (i < POPSIZE / 2)
                {
                    y = (alpha * population[i].gene[j])+((1 - alpha) * population[i + divider].gene[j]);
                    y = min(max(y, varMin[j]), varMax[j]);
                }
                else
                {
                    y = (alpha * population[i + divider].gene[j])+((1 - alpha) * population[i].gene[j]);
                    y = min(max(y, varMin[j]), varMax[j]);
                }

                if (y != 0)
                {
                    newPopulation[newPopulationSize].gene[j] = y;
                }
                else
                {
                    newPopulation[newPopulationSize].gene[j] = population[i].gene[j];
                }
            }
            newPopulationSize++;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

void mutate(int& seed)
{
    const double a = 0.0;
    const double b = 1.0;
    double lbound;
    double ubound;
    double x;
    int B1Pos = 6;
    int B2Pos = 7;
    int C1Pos = 4;
    int C2Pos = 5;
    double randNum;

    for (int i = 0; i < POPSIZE; i++)
    {
        x = r8_uniform_ab(a, b, seed);
        if (x < PMUTATION)
        {
            newPopulation[newPopulationSize].fitness = 0;
            //            newPopulation[newPopulationSize].rfitness = 0;
            //            newPopulation[newPopulationSize].cfitness = 0;
            for (int j = 0; j < NVARS; j++)
            {
                newPopulation[newPopulationSize].lower[j] = varMin[j];
                newPopulation[newPopulationSize].upper[j] = varMax[j];
                x = r8_uniform_ab(a, b, seed);
                if (x < PMUTATION)
                {
                    lbound = population[i].lower[j];
                    ubound = population[i].upper[j];
                    randNum = r8_uniform_ab(lbound, ubound, seed);
                    if (j == B1Pos)
                    {
                        while (4 * population[i].gene[C1Pos] <= pow(randNum, 2))
                        {
                            randNum = r8_uniform_ab(lbound, ubound, seed);
                        }
                    }
                    else if (j == B2Pos)
                    {
                        while (4 * population[i].gene[C2Pos] <= pow(randNum, 2))
                        {
                            randNum = r8_uniform_ab(lbound, ubound, seed);
                        }
                    }

                    newPopulation[newPopulationSize].gene[j] = randNum;
                }
                else
                {
                    newPopulation[newPopulationSize].gene[j] = population[i].gene[j];
                }
            }
            newPopulationSize++;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

void newEvaluate()
{
    for (int i = 0; i < POPSIZE; i++)
    {
        nInf = population[i].gene[0];
        Efb = population[i].gene[1];
        A1 = population[i].gene[2];
        A2 = population[i].gene[3];
        C1 = population[i].gene[4];
        C2 = population[i].gene[5];
        B1 = population[i].gene[6];
        B2 = population[i].gene[7];
        Ed = population[i].gene[8];
        Egamma0 = population[i].gene[9];
        s = population[i].gene[10];
        d = population[i].gene[11];

        population[i].fitness = fitnessCalculator();
    }

    for (int i = 0; i < newPopulationSize; i++)
    {
        nInf = newPopulation[i].gene[0];
        Efb = newPopulation[i].gene[1];
        A1 = newPopulation[i].gene[2];
        A2 = newPopulation[i].gene[3];
        C1 = newPopulation[i].gene[4];
        C2 = newPopulation[i].gene[5];
        B1 = newPopulation[i].gene[6];
        B2 = newPopulation[i].gene[7];
        Ed = newPopulation[i].gene[8];
        Egamma0 = newPopulation[i].gene[9];
        s = newPopulation[i].gene[10];
        d = newPopulation[i].gene[11];

        newPopulation[i].fitness = fitnessCalculator();
    }
}

//////////////////////////////////////////////////////////////////////////////

void sortPopulation()
{
    for (int i = 0; i < POPSIZE; i++)
    {
        for (int j = 0; j < NVARS; j++)
        {
            tempPopulation[i].gene[j] = population[i].gene[j];
        }
        tempPopulation[i].fitness = population[i].fitness;
        bestPopulationSize++;
    }

    bestPopulationSize = POPSIZE;
    int counter = POPSIZE;

    for (int i = 0; i < newPopulationSize; i++)
    {
        for (int j = 0; j < NVARS; j++)
        {
            tempPopulation[counter + i].gene[j] = newPopulation[i].gene[j];
        }
        tempPopulation[counter + i].fitness = newPopulation[i].fitness;
        bestPopulationSize++;
    }

    bestPopulation.clear();

    for (int i = 0; i < bestPopulationSize; i++)
    {
        bestPopulation.push_back(genotype());
        bestPopulation[i].fitness = tempPopulation[i].fitness;
        for (int j = 0; j < NVARS; j++)
        {
            bestPopulation[i].lower[j] = varMin[j];
            bestPopulation[i].upper[j] = varMax[j];
            bestPopulation[i].gene[j] = tempPopulation[i].gene[j];
        }
    }

    //    cout << "+++++++++++++++++++++++++++++++++++++++++++ START BEST" << endl;
    //    for (int i = 0; i < bestPopulationSize; i++)
    //    {
    //        cout << bestPopulation[i].fitness << endl;
    //    }

    //    cout << "------------------------------------------- START SORT" << endl;

    sort(bestPopulation.begin(), bestPopulation.end(), sortByFitness);

    //    for (int i = 0; i < bestPopulationSize; i++)
    //    {
    //        cout << bestPopulation[i].fitness << endl;
    //    }
    //
    //    cout << "******************************************* END SORT" << endl;
}

//////////////////////////////////////////////////////////////////////////////

void selector()
{
    for (int i = 0; i < POPSIZE; i++)
    {
        population[i].fitness = bestPopulation[i].fitness;
        for (int j = 0; j < NVARS; j++)
        {
            population[i].lower[j] = bestPopulation[i].lower[j];
            population[i].upper[j] = bestPopulation[i].upper[j];
            population[i].gene[j] = bestPopulation[i].gene[j];
        }
    }

    newPopulationSize = 0;
    bestPopulationSize = 0;
}

//////////////////////////////////////////////////////////////////////////////

bool sortByFitness(const genotype &a, const genotype &b)
{
    return a.fitness < b.fitness;
}