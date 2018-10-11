/* 
 * File:   GeneticAlgorithem.hpp
 * Author: Erfan Jazeb Nikoo
 *
 */

#include "Utilities.hpp"

using namespace std;

//****************************************************************************80

int i4_uniform_ab(int a, int b, int &seed)

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
    int c;
    const int i4_huge = 2147483647;
    int k;
    float r;
    int value;

    if (seed == 0)
    {
        cerr << "\n";
        cerr << "I4_UNIFORM_AB - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit(1);
    }
    //
    //  Guarantee A <= B.
    //
    if (b < a)
    {
        c = a;
        a = b;
        b = c;
    }

    k = seed / 127773;

    seed = 16807 * (seed - k * 127773) - k * 2836;

    if (seed < 0)
    {
        seed = seed + i4_huge;
    }

    r = (float) (seed) * 4.656612875E-10;
    //
    //  Scale R to lie between A-0.5 and B+0.5.
    //
    r = (1.0 - r) * ((float) a - 0.5)
            + r * ((float) b + 0.5);
    //
    //  Use rounding to convert R to an integer between A and B.
    //
    value = round(r);
    //
    //  Guarantee A <= VALUE <= B.
    //
    if (value < a)
    {
        value = a;
    }
    if (b < value)
    {
        value = b;
    }

    return value;
}
//****************************************************************************80

double r8_uniform_ab(double a, double b, int &seed)

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
    int i4_huge = 2147483647;
    int k;
    double value;

    if (seed == 0)
    {
        cerr << "\n";
        cerr << "R8_UNIFORM_AB - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit(1);
    }

    k = seed / 127773;

    seed = 16807 * (seed - k * 127773) - k * 2836;

    if (seed < 0)
    {
        seed = seed + i4_huge;
    }

    value = (double) (seed) * 4.656612875E-10;

    value = a + (b - a) * value;

    return value;
}
//****************************************************************************80

void timestamp()

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
{
#define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    size_t len;
    time_t now;

    now = time(NULL);
    tm = localtime(&now);

    len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

    cout << time_buffer << "\n";

    return;
#undef TIME_SIZE
}
//****************************************************************************80

