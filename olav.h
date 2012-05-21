// Created: olav, Mon Nov 17 16:56:02 CET 2008
// some utility functions to make things a little bit nicer to use

#include <iostream>
#include <cstring>
#include <cmath>
#include <sstream>
#include <ctime>
#include <limits>

void status(unsigned long x, unsigned long R, unsigned long N);

std::string sec2string(double seconds);
std::string sec2string(long seconds);

std::string ETAstring(long completed, long total, double elapsedtime);
std::string ETAstring(long completed, long total, time_t starttime);
