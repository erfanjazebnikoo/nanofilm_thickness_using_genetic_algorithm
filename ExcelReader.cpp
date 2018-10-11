/* 
 * File:   ExcelReader.cpp
 * Author: Erfan Jazeb Nikoo
 */

#include "ExcelReader.hpp"

using namespace std;

ExcelReader::ExcelReader(char* filename, int row)
: filename(filename), row(row)
{
}

int ExcelReader::readTable(int selectRow, double *output)
{
    int p = 0;
    int counter = 0;
    ifstream file(this->filename);
    string rowValue, value;
    while (file.good())
    {
        getline(file, rowValue);
        istringstream tmp(rowValue);
        while (tmp.good())
        {
            getline(tmp, value, ',');
            if (p == 0 || p == 1)
            {
                p++;
                continue;
            }
            if (p % this->row == selectRow)
            {
                output[counter] = atof(value.c_str());
                counter++;
            }
            p++;
        }
    }
    return counter;
}

ExcelReader::~ExcelReader()
{
}

