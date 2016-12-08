
/********************************************
*  Type Definitions for Rigid Body Class    *
*  Version : 0.1                            *
*  Author: Emant                            *
*  Year : 2014                              *
*********************************************/

#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED
#include <array>

// use double precision
#define _DOUBLE

#ifdef _DOUBLE
    #define FLOAT double
#else
    #define FLOAT float
#endif // _DOUBLE

#define PI 3.14159265

typedef std::array<FLOAT,3> vect;
typedef std::array<FLOAT,6> vect6;
typedef std::array<FLOAT,4> quat;
typedef std::array<FLOAT,13> vstate;
typedef std::array<FLOAT,16> matrix;


#endif // DEFS_H_INCLUDED
