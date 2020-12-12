#include "PT.h"
#include "triangulation.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

int PT_swap(triangulation &singlemesh_a, triangulation &singlemesh_b) {
    double beta_a, beta_b;
    double lam_a, lam_b;
    double Le_a, Le_b;
    double E_a, E_b;
    beta_a = singlemesh_a.beta;
    beta_b = singlemesh_b.beta;
    lam_a = singlemesh_a.lam;
    lam_b = singlemesh_b.lam;
    Le_a = singlemesh_a.Les[0];
    Le_b = singlemesh_b.Les[0];
    E_a = singlemesh_a.E;
    E_b = singlemesh_b.E;

    // only swap lam_a lam_b!!

    // W_init = std::exp(-beta_a * E_a - beta_b * E_b);
    // W_swap = std::exp(-beta_a * E_a - beta_b * E_b -(lam_b - lam_a) * (Le_a
    // * beta_a - Le_b * beta_b));
    // W_init = Exp[-E_a]*Exp[-E_b]
    // W_swap = Exp[-E_a-(lam_b-lam_a)Le_a]*Exp[-E_b-(lam_a-lam_b)Le_b]
    if (singlemesh_a.rand_uni(singlemesh_a.gen) <
        std::exp((lam_b - lam_a) * (Le_b * beta_b - Le_a * beta_a))) {
        // accept swap! do it!
        singlemesh_a.lam = lam_b;
        singlemesh_b.lam = lam_a;
        singlemesh_a.E = E_a + (lam_b - lam_a) * Le_a;
        singlemesh_b.E = E_b + (lam_a - lam_b) * Le_b;
        triangulation exchange = singlemesh_a;
        singlemesh_a = singlemesh_b;
        singlemesh_b = exchange;
        return 1;
    }
    return 0;
}

void PT_Thermal(std::vector<triangulation> multimesh, int sweeps_p_PT,
                int MC_sweeps, int step_p_sweep, double ds) {
    // thermalisation utilizing PT
    int N2swap;
    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++) {
        // std::cout << sweep_n * 1.0 / MC_sweeps * 100 << "\%\n";
        for (int j = 0; j < multimesh.size(); j++) {
            for (int i = 0; i < step_p_sweep; i++) {
                multimesh[j].vertex_metropolis(ds);
                multimesh[j].bond_metropolis();
                if (i % int(std::sqrt(multimesh[j].N)) == 0) {
                    multimesh[j].edge_metropolis();
                }
            }
        }
        if (sweep_n % sweeps_p_PT == 0) {
            N2swap = int((multimesh.size() - 1) *
                         multimesh[0].rand_uni(multimesh[0].gen));
            PT_swap(multimesh[N2swap], multimesh[N2swap + 1]);
        }
    }
}

void PT_O_MC_measure(std::vector<triangulation> multimesh, int sweeps_p_PT,
                     int MC_sweeps, int step_p_sweep, double ds,
                     std::string folder) {
    // need to put data into multimesh.size() files

    int mesh_num = multimesh.size();
    std::vector<std::vector<double>> PT_E_all(mesh_num);
    std::vector<std::vector<double>> PT_Le_all(mesh_num);
    int N2swap;
    std::string filename;
    std::clock_t c_start = std::clock();
    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++) {
        for (int j = 0; j < multimesh.size(); j++) {
            // update every mesh during each sweep
            for (int i = 0; i < step_p_sweep; i++) {
                multimesh[j].vertex_metropolis(ds);
                multimesh[j].bond_metropolis();
                if (i % int(std::sqrt(multimesh[j].N)) == 0) {
                    multimesh[j].edge_metropolis();
                }
            }
            PT_E_all[j].push_back(multimesh[j].E);
            PT_Le_all[j].push_back(multimesh[j].Les[0]);
        }
        if (sweep_n % sweeps_p_PT == 0) {
            N2swap = int((multimesh.size() - 1) *
                         multimesh[0].rand_uni(multimesh[0].gen));
            PT_swap(multimesh[N2swap], multimesh[N2swap + 1]);
        }
    }
    std::clock_t c_end = std::clock();

    for (int j = 0; j < multimesh.size(); j++) {
        std::string kappa_s = std::to_string(multimesh[j].kar);
        int pos = kappa_s.find(".");
        kappa_s = kappa_s.substr(0, pos + 2);
        std::string lam_s = std::to_string(multimesh[j].lam);
        pos = lam_s.find(".");
        lam_s = lam_s.substr(0, pos + 2);
        filename =
            folder + "/PT_O_MC_kappa" + kappa_s + "_lam" + lam_s + ".txt";
        std::ofstream f(filename);
        if (f.is_open()) {
            f << "N: " << multimesh[j].N << "\n";
            f << "Lr: " << multimesh[j].Lr << "\n";
            f << "l0: " << multimesh[j].l0 << "\n";
            f << "step_p_sweep: " << step_p_sweep << "\n";
            f << "MC_sweeps: " << MC_sweeps << "\n";
            f << "CPUtime: " << double(c_end - c_start) * 1.0 / CLOCKS_PER_SEC
              << "\n";
            f << "E_total,Le\n";
            for (int i = 0; i < PT_E_all[j].size(); i++) {
                f << PT_E_all[j][i] << "," << PT_Le_all[j][i] << "\n";
            }
        }
        f.close();
    }
}