/* 
 * File:   ExcelReader.hpp
 * Author: Erfan Jazeb Nikoo
 */

#ifndef EXCELREADER_HPP
#define	EXCELREADER_HPP

#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<sstream>
#include <math.h>

class ExcelReader
{
public:

    ExcelReader(char* filename, int row);
    int readTable(int selectRow, double *output);
    virtual ~ExcelReader();
private:
    int row;
    char* filename;
};

#endif	/* EXCELREADER_HPP */

