#define _USE_MATH_DEFINES
#include <string>
#include <sstream>
#include "stdio.h"
#include <fstream>
#include <math.h>

double myRound(double x, double mode = 1){ //round to nearest mode
	x = x/mode;
	double upper = ceil(x);
	double lower = floor(x);
	if (x+0.5 >= upper){
		x = upper;
	} else{
		x = lower;
	}
	x = x*mode;
	return(x);
}

std::string doubleToString(double inDouble){
	return static_cast<std::ostringstream*>( &(std::ostringstream() << inDouble) )->str();
}

double stringToDouble( const std::string& s ){
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

double deltaphi(double a, double b) {
	double dphi = a-b;
	if (a >= M_PI) {
		dphi -= 2*M_PI;
	} else if (dphi <= -M_PI) {
		dphi += 2*M_PI;
	}
	return dphi;
}
