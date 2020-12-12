// Copyright[2018] [Lijie Ding]
//#include "PT.h"
#include "triangulation.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char const *argv[]) {
    std::clock_t c_start = std::clock();
    double beta = 1;
    int N;      // number of beads
    int Ne = 1; // number of edges
    int Lr;     // number of bead per line
    int L;      // number of bead on vertical direction, used for cylinder
                // initialization
    // also help adjusting the fixed distance of two beads
    double d0 = 1.2;
    // double l0 = std::sqrt(3); // constrained by self-avoidance limit<sqrt(3),
    double l0 = 1.56;
    double l1 = 1.95; // only constrained by triangle limit l1<l0*(4-l0^2)
    double kar;
    double C0;
    double karg;
    double lam;
    double B;
    double Bc;
    double tau_0;
    double ds = 0.2;

    std::string folder;

    N = std::atoi(argv[1]);
    Ne = std::atoi(argv[2]);
    L = std::atoi(argv[3]);

    kar = std::atof(argv[4]);
    C0 = std::atof(argv[5]);
    karg = std::atof(argv[6]);
    lam = std::atof(argv[7]);
    B = std::atof(argv[8]);
    Bc = std::atof(argv[9]);
    tau_0 = std::atof(argv[10]);
    bool fix_bead_on = 1;
    if (L == 0) {
        fix_bead_on = 0;
        L = int(std::sqrt(N));
        if (Ne == 2) {
            L = 10;
            if (N % L != 0) {
                std::cout << "N % L != 0 \n";
            }
            fix_bead_on = 0;
            Lr = N / L;
        }
    }

    std::string finfo;
    finfo = "N" + std::string(argv[1]) + "_Ne" + std::string(argv[2]) + "_L" +
            std::string(argv[3]) + "_kar" + std::string(argv[4]) + "_C0" +
            std::string(argv[5]) + "_karg" + std::string(argv[6]) + "_lam" +
            std::string(argv[7]) + "_B" + std::string(argv[8]) + "_Bc" +
            std::string(argv[9]) + "_tau0" + std::string(argv[10]);
    if (argc == 12) {
        // use "triangulation kar lam local" for local running
        // used for local running!
        std::cout << "running on local machine\n";
        folder = "../data/scratch_local";
        // singlemesh.State_write(folder + "/State_" + finfo + "_init.txt");
        // singlemesh.State_load(folder + "/Init_State_kar8.0_lam8.0.txt");
        // singlemesh.shape_set_helix(0.1);
        // singlemesh.O_reset();
        triangulation singlemesh(beta, N, Ne, L, d0, l0, l1, kar, C0, karg, lam,
                                 B, Bc, tau_0);

        if (fix_bead_on) {
            int fix_bead0 = int(N / (2 * L)) * L - 1;
            if (std::abs((singlemesh.mesh[fix_bead0 + L - 1].pos[0] -
                          singlemesh.mesh[fix_bead0].pos[0]) -
                         (L * d0 - d0)) > 0.1) {
                std::cout << "not exact fixed length\n";
            }
            singlemesh.fixed_beads = {fix_bead0, fix_bead0 + L - 1};
        }
        singlemesh.State_write(folder + "/State_" + finfo + "_init.txt");
        singlemesh.Thermal(0, int(N / (ds * ds)) + 1, 2, ds);
        // singlemesh.O_MC_measure(5, 1, int(N / (ds * ds)) + 1, ds,
        singlemesh.O_MC_measure(10, 11, int(N / (ds * ds)) + 1, ds,
                                folder + "/O_MC_" + finfo + ".txt",
                                folder + "/Cncnc_MC_" + finfo + ".txt",
                                folder + "/Crr_MC_" + finfo + ".txt",
                                folder + "/Iij_MC_" + finfo + ".txt");
        singlemesh.State_write(folder + "/State_" + finfo + ".txt");

        return 0;
    } else if (argc == 11) {
        // ccv running
        folder = "/users/lding3/scratch";

        triangulation singlemesh(beta, N, Ne, L, d0, l0, l1, kar, C0, karg, lam,
                                 B, Bc, tau_0); // kar~6000, kar_g~100, lam
        if (fix_bead_on) {
            int fix_bead0 = int(N / (2 * L)) * L - 1;
            if ((singlemesh.mesh[fix_bead0 + L - 1].pos[0] -
                 singlemesh.mesh[fix_bead0].pos[0]) != L * d0 - d0) {
                std::cout << "not exact fixed length\n";
            }
            singlemesh.fixed_beads = {fix_bead0, fix_bead0 + L - 1};
        }
        singlemesh.Thermal(20000, N / (ds * ds), 10, ds);
        singlemesh.Thermal(20000, N / (ds * ds), 1, ds);
        singlemesh.O_MC_measure(40000, 100, N / (ds * ds) + 1, ds,
                                folder + "/O_MC_" + finfo + ".txt",
                                folder + "/Cncnc_MC_" + finfo + ".txt",
                                folder + "/Crr_MC_" + finfo + ".txt",
                                folder + "/Iij_MC_" + finfo + ".txt");
        singlemesh.State_write(folder + "/State_" + finfo + ".txt");

        return 0;
    }
}