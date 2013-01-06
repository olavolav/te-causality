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

#include "miniprofiler.h"

MiniProfiler::MiniProfiler() {
  // std::cout <<"DEBUG: MiniProfiler#init."<<std::endl;
}

MiniProfiler::~MiniProfiler() {
  // std::cout <<"DEBUG: MiniProfiler#delete."<<std::endl;
}

int MiniProfiler::number_of_registered_tasks() {
  return tasks.size();
}

void MiniProfiler::register_task(std::string t_name) {
  struct MiniProfilerTask* new_task = new struct MiniProfilerTask;
  new_task->name = t_name;
  new_task->elapsed_clock_ticks = 0;
  new_task->currently_active = false;

  tasks.push_back( *new_task );
}

int MiniProfiler::search(std::string t_name) {
  for(int i = 0; i < tasks.size(); i++) {
    if( tasks[i].name == t_name ) {
      // std::cout <<"DEBUG: found the task."<<std::endl;
      return i;
    }
  }
  return -1;
}

void MiniProfiler::resuming_task(std::string t_name) {
  int t_index = search(t_name);
  if( t_index == -1 ) {
    std::cout <<"Error in MiniProfiler#resuming_task: Invalid task key."<<std::endl;
    exit(1);
  }
  if( tasks[t_index].currently_active ) {
    std::cout <<"Warning in MiniProfiler#resuming_task: task '"<<t_name<<"' is already running!"<<std::endl;
  } else {
    tasks[t_index].clock_tick_of_activation = clock();
    tasks[t_index].currently_active = true;
  }
}

void MiniProfiler::stopping_task(std::string t_name) {
  int t_index = search(t_name);
  if( t_index == -1 ) {
    std::cout <<"Error in MiniProfiler#stopping_task: Invalid task key."<<std::endl;
    exit(1);
  }
  if( !tasks[t_index].currently_active ) {
    std::cout <<"Warning in MiniProfiler#stopping_task: task '"<<t_name<<"' is not running!"<<std::endl;
  } else {
    tasks[t_index].elapsed_clock_ticks += clock() - tasks[t_index].clock_tick_of_activation;
    tasks[t_index].currently_active = false;
  }
}

std::string MiniProfiler::summary() {
  clock_t now = clock();
  std::stringstream summ;
  for(int i = 0; i < tasks.size(); i++) {
    summ  <<"- "<<tasks[i].name<<": ";
    if(tasks[i].currently_active) {
      summ <<float(tasks[i].elapsed_clock_ticks + (now - tasks[i].clock_tick_of_activation))/CLOCKS_PER_SEC;
    } else {
      summ  <<float(tasks[i].elapsed_clock_ticks)/CLOCKS_PER_SEC;
    }
    summ  <<"s\n";
  }
  return summ.str();
}
