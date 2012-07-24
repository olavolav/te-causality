// Copyright 2012, Olav Stetter
// 
// This file is part of TE-Causality.
// 
// TE-Causality is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// TE-Causality is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with TE-Causality.  If not, see <http://www.gnu.org/licenses/>.

// Some utility functions to make things a little bit nicer to use.

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
