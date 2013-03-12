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

// This is a very simple (and probably not very accurate) time
// measurement object to see which tasks of a given program are taking
// how much time to compute.

#include <cstdlib>
#include <iostream>
#include <cstring>
#include <sstream>
#include <time.h>
#include <vector>

#define __MINIPROFILER_H

struct MiniProfilerTask {
  std::string name;
  clock_t elapsed_clock_ticks;
  bool currently_active;
  clock_t clock_tick_of_activation;
};


class MiniProfiler {
  
private:
  std::vector<MiniProfilerTask> tasks;
  int search(const std::string& t_name) const;
  
public:
  MiniProfiler();
  ~MiniProfiler();
  int number_of_registered_tasks() const;
  void register_task(const std::string& t_name);
  void resuming_task(const std::string& t_name);
  void stopping_task(const std::string& t_name);
  
  std::string summary();
};
