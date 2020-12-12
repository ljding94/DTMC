#ifndef _PT_H
#define _PT_H
#include "triangulation.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

int PT_swap(triangulation &singlemesh_a, triangulation &singlemesh_b);

void PT_Thermal(std::vector<triangulation> multimesh, int sweeps_p_PT,
                int MC_sweeps, int step_p_sweep, double ds);

void PT_O_MC_measure(std::vector<triangulation> multimesh, int sweeps_p_PT,
                     int MC_sweeps, int step_p_sweep, double ds,
                     std::string folder);
#endif