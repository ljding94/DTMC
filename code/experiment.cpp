#include "triangulation.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

void triangulation::State_write(std::string filename) {
    std::ofstream f(filename);
    if (f.is_open()) {
        // find maxim nei.size()
        int max_nei_size = 0;
        for (int i = 0; i < mesh.size(); i++) {
            if (mesh[i].nei.size() > max_nei_size) {
                max_nei_size = mesh[i].nei.size();
            }
        }
        f << "beta=" << beta << "\n";
        f << "N=" << N << "\n";
        f << "L=" << L << "\n";
        f << "l0=" << l0 << "\n";
        f << "l1=" << l1 << "\n";
        f << "kar=" << kar << "\n";
        f << "kar_g=" << karg << "\n";
        f << "B=" << B << "\n";
        f << "Bc=" << Bc << "\n";
        f << "tau_0=" << tau0 << "\n";
        f << "lam=" << lam << "\n";
        f << "E=" << E << "\n";
        f << "Les=";
        for (int e = 0; e < Ne; e++) {
            f << Les[e] << ",";
        }
        f << "\n";
        f << "I2H2=" << I2H2 << "\n";
        f << "Ikg=" << Ikg << "\n";
        f << "Ik2=" << Ik2 << "\n";
        f << "Itau=" << Itau << "\n";
        f << "max_nei_size=" << max_nei_size << "\n";
        f << "ds,dA,2H,dskg,dsk2,dstau,dE,x,y,z,edge_num,edge_neibs,neibs";
        for (int i = 0; i < mesh.size(); i++) {
            f << "\n"
              << mesh[i].ds << "," << mesh[i].dAn2H[0] << ","
              << mesh[i].dAn2H[1] << "," << mesh[i].dskg << "," << mesh[i].dsk2
              << "," << mesh[i].dstau << "," << mesh[i].dE;
            f << "," << mesh[i].pos[0] << "," << mesh[i].pos[1] << ","
              << mesh[i].pos[2] << "," << mesh[i].edge_num;
            for (int j = 0; j < mesh[i].edge_nei.size(); j++) {
                f << "," << mesh[i].edge_nei[j];
            }
            for (int j = 0; j < 2 - mesh[i].edge_nei.size(); j++) {
                f << ",-1";
            }
            for (int j = 0; j < mesh[i].nei.size(); j++) {
                f << "," << mesh[i].nei[j];
            }
            for (int j = 0; j < max_nei_size - mesh[i].nei.size(); j++) {
                f << ",-1";
            }
        }
    }
    f.close();
} // end State_write

void triangulation::State_load(std::string state_file) {
    std::ifstream f(state_file);
    std::string buff;
    int max_nei_size;
    while (f.good()) {
        // skip 8 lines
        for (int i = 0; i < 13; i++) {
            std::getline(f, buff);
            // std::cout << buff << "\n";
        }
        // load geometric obsevables
        std::getline(f, buff, '=');
        std::getline(f, buff);
        Les[0] = std::stof(buff);
        // kinda don't wanna use it any more
        std::getline(f, buff, '=');
        std::getline(f, buff);
        I2H2 = std::stof(buff);
        std::getline(f, buff, '=');
        std::getline(f, buff);
        Ikg = std::stof(buff);
        std::getline(f, buff, '=');
        std::getline(f, buff);
        Ik2 = std::stof(buff);
        std::getline(f, buff, '=');
        std::getline(f, buff);
        Itau = std::stof(buff);
        std::getline(f, buff, '=');
        std::getline(f, buff);
        max_nei_size = std::stof(buff);
        // skipp another line
        std::getline(f, buff);
        // std::cout << buff << "\n";

        int mcount = 0;

        for (mcount = 0; mcount < mesh.size(); mcount++) {
            // load observables
            std::getline(f, buff, ',');
            mesh[mcount].ds = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dAn2H[0] = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dAn2H[1] = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dskg = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dsk2 = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dstau = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dE = std::stof(buff);
            // load pos
            for (int k = 0; k < 3; k++) {
                std::getline(f, buff, ',');
                mesh[mcount].pos[k] = std::stof(buff);
            }
            // load edge_nei
            mesh[mcount].edge_nei.clear();
            for (int k = 0; k < 2; k++) {
                std::getline(f, buff, ',');
                if (std::stof(buff) != -1) {
                    mesh[mcount].edge_nei.push_back(std::stof(buff));
                }
            }
            // load neibs
            // should read to the end of line, but how?
            mesh[mcount].nei.clear();
            for (int k = 0; k < max_nei_size - 1; k++) {
                std::getline(f, buff, ',');
                if (std::stof(buff) != -1) {
                    mesh[mcount].nei.push_back(std::stof(buff));
                }
            }
            std::getline(f, buff);
            if (std::stof(buff) != -1) {
                mesh[mcount].nei.push_back(std::stof(buff));
            }
        }
    }

    // also need to update v2flip, edge list
    edge_lists.clear();
    v2flip_list.clear();
    E = 0;
    for (int ind = 0; ind < mesh.size(); ind++) {
        update_v2flip_local(ind);
        mesh[ind].dE = dE_m(ind);
        E += mesh[ind].dE;
    }
    // and the other list
    v2add_list.clear();
    renew_v2add_list();
    ve2remove_list.clear();
    renew_ve2remove_list();
}

void triangulation::State_write_seq(std::string filename, int MC_sweeps,
                                    int step_p_sweep, double ds_) {
    std::string savename;
    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++) {
        Thermal(1, step_p_sweep, 1, ds_);
        savename = filename;
        savename.insert(savename.size() - 4, "_" + std::to_string(sweep_n));
        State_write(savename);
    }
}

void triangulation::shape_set(double theta) {
    // set a spherical cap with angle θ
    // R = r/cos(θ)
    double R2;
    R2 = std::pow(mesh[edge_lists[0][0]].pos[0], 2) +
         std::pow(mesh[edge_lists[0][0]].pos[1], 2);
    R2 = R2 / std::pow(std::sin(theta), 2);
    for (int i = 0; i < mesh.size(); i++) {
        mesh[i].pos[2] = std::sqrt(R2 - std::pow(mesh[i].pos[0], 2) -
                                   std::pow(mesh[i].pos[1], 2)) -
                         std::sqrt(R2) * std::cos(theta);
    }
    E = 0;
    for (int e = 0; e < Ne; e++) {
        Les[e] = 0;
    }
    for (int i = 0; i < mesh.size(); i++) {
        mesh[i].ds = ds_m(i);
        mesh[i].dAn2H = dAn2H_m(i);
        mesh[i].dskg = dskg_m(i);
        mesh[i].dsk2 = dsk2_m(i);
        mesh[i].dstau = dstau_m(i);
        mesh[i].dE = dE_m(i);
        E += mesh[i].dE;
        Les[mesh[i].edge_num] += mesh[i].ds;
    }
}
void triangulation::shape_set_helix(double q) {
    for (int i = 0; i < N; i++) {
        mesh[i].pos[2] = mesh[i].pos[1] * std::sin(q * mesh[i].pos[0]);
        mesh[i].pos[1] = mesh[i].pos[1] * std::cos(q * mesh[i].pos[0]);
    }
}

void triangulation::Thermal(int MC_sweeps, int step_p_sweep, int beta_steps,
                            double delta_s) {
    beta = 0.0;
    for (int n_beta = 0; n_beta < beta_steps; n_beta++) {
        beta += 1.0 / beta_steps;
        for (int sweep_n = 0; sweep_n < MC_sweeps / beta_steps; sweep_n++) {
            // std::cout << sweep_n << "/" << MC_sweeps << "\n";
            for (int i = 0; i < step_p_sweep; i++) {
                vertex_metropolis(delta_s);
                bond_metropolis();
                if (i % int(std::sqrt(N)) == 0) {
                    edge_metropolis();
                }
            }
            // std::cout << "thermo, beta=" << beta << "," << sweep_n << "/"<<
            // MC_sweeps << "\n";
        }
    }
}
void triangulation::O_MC_measure(int MC_sweeps, int step_p_Cm, int step_p_sweep,
                                 double delta_s, std::string filename,
                                 std::string filename_C,
                                 std::string filename_Cr,
                                 std::string filename_I) {
    std::vector<double> E_all;
    std::vector<std::vector<double>> Les_all;
    Les_all.resize(Ne);
    std::vector<double> IdA_all;
    std::vector<double> I2H_all;
    std::vector<double> I2H2_all;
    std::vector<double> Ikg_all;
    std::vector<double> Ik2_all;
    std::vector<double> Itau_all;

    std::vector<std::vector<double>> Iij_all;

    std::vector<std::vector<double>> Cncnc0_all;
    std::vector<std::vector<double>> Cncnc1_all;
    std::vector<std::vector<double>> Cncnc2_all;

    std::vector<std::vector<double>> Crr0_all;
    std::vector<std::vector<double>> Crr1_all;
    std::vector<std::vector<double>> Crr2_all;

    double vertex_accept = 0;
    double bond_accept = 0;
    double edge_accept = 0;
    std::clock_t c_start = std::clock();
    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++) {
        for (int i = 0; i < step_p_sweep; i++) {
            vertex_accept += vertex_metropolis(delta_s);
            bond_accept += bond_metropolis();
            // edge_accept += edge_metropolis();

            if (i % int(std::sqrt(N)) == 0) {
                edge_accept += edge_metropolis();
            }
        }
        // std::cout << sweep_n << "/" << MC_sweeps << "\n";
        // std::cout << "E,Le,I2H2,Ikg,Ik2,Itau=" << E << "," << Le << "," <<
        // I2H2
        //<< "," << Ikg << "," << Ik2 << "," << Itau << "\n";
        E_all.push_back(E);
        for (int e = 0; e < Ne; e++) {
            Les_all[e].push_back(Les[e]);
        }
        IdA_all.push_back(IdA);
        I2H_all.push_back(I2H);
        I2H2_all.push_back(I2H2);
        Ikg_all.push_back(Ikg);
        Ik2_all.push_back(Ik2);
        Itau_all.push_back(Itau);

        if (sweep_n % step_p_Cm == 0) {
            Iij_all.push_back(Iij_m());
            // Cncnc0_all.push_back(Cncnc_m(0));
            // Crr0_all.push_back(Crr_m(0));
            if (Ne >= 2) {
                // Cncnc1_all.push_back(Cncnc_m(1));
                // Crr1_all.push_back(Crr_m(1));
            }
            if (Ne >= 3) {
                // Cncnc2_all.push_back(Cncnc_m(2));
                // Crr2_all.push_back(Crr_m(2));
            }
        }
        /*
        for (int j = 0; j < Cncnc_all.back().size(); j++) {
            std::cout << "," << Cncnc_all.back()[j];
        }
        std::cout << "\n";
        */
        // std::cout << "MC_m" << sweep_n << "/" << MC_sweeps << "\n";
    }
    vertex_accept /= MC_sweeps * step_p_sweep;
    bond_accept /= MC_sweeps * step_p_sweep;
    edge_accept /= MC_sweeps * (step_p_sweep / int(std::sqrt(N)));
    // edge_accept /= MC_sweeps * step_p_sweep;
    std::clock_t c_end = std::clock();
    std::ofstream f(filename);
    std::ofstream f_I(filename_I);
    std::ofstream f_C0(filename_C.substr(0, filename_C.length() - 4) +
                       "_e0.txt");
    std::ofstream f_Cr0(filename_Cr.substr(0, filename_Cr.length() - 4) +
                        "_e0.txt");

    if (f.is_open()) {
        f << "beta=" << beta << "\n";
        f << "N=" << N << "\n";
        f << "Ne=" << Ne << "\n";
        f << "L=" << L << "\n";
        f << "l0=" << l0 << "\n";
        f << "l1=" << l1 << "\n";
        f << "step_p_sweep=" << step_p_sweep << "\n";
        f << "MC_sweeps=" << MC_sweeps << "\n";
        f << "CPUtime=" << double(c_end - c_start) * 1.0 / CLOCKS_PER_SEC
          << "\n";
        f << "vertex_accept," << vertex_accept << "\n";
        f << "bond_accept," << bond_accept << "\n";
        f << "edge_accept," << edge_accept << "\n";
        f << "E,";
        for (int e = 0; e < Ne; e++) {
            f << "Les[" << e << "],";
        }
        f << "IdA,I2H,I2H2,Ikg,Ik2,Itau\n";
        for (int i = 0; i < E_all.size(); i++) {
            f << E_all[i] << ",";
            for (int e = 0; e < Ne; e++) {
                f << Les_all[e][i] << ",";
            }
            f << IdA_all[i] << "," << I2H_all[i] << "," << I2H2_all[i] << ","
              << Ikg_all[i] << "," << Ik2_all[i] << "," << Itau_all[i] << "\n";
        }
    }
    f.close();
    if (f_I.is_open()) {
        f_I << "Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz\n";
        for (int i = 0; i < Iij_all.size(); i++) {
            f_I << Iij_all[i][0];
            for (int j = 1; j < Iij_all[i].size(); j++) {
                f_I << "," << Iij_all[i][j];
            }
            f_I << "\n";
        }
    }
    f_I.close();
    if (f_C0.is_open()) {
        f_C0 << "Len, |<nc>|,Cncnc0\n";
        for (int i = 0; i < Cncnc0_all.size(); i++) {
            f_C0 << Cncnc0_all[i].size() - 1;
            for (int j = 0; j < Cncnc0_all[i].size(); j++) {
                f_C0 << "," << Cncnc0_all[i][j];
            }
            f_C0 << "\n";
        }
    }
    f_C0.close();
    if (f_Cr0.is_open()) {
        f_Cr0 << "Len,<r>,Crr0\n";
        for (int i = 0; i < Crr0_all.size(); i++) {
            f_Cr0 << Crr0_all[i].size() - 1;
            for (int j = 0; j < Crr0_all[i].size(); j++) {
                f_Cr0 << "," << Crr0_all[i][j];
            }
            f_Cr0 << "\n";
        }
    }
    f_Cr0.close();

    if (Ne >= 2) {
        std::ofstream f_C1(filename_C.substr(0, filename_C.length() - 4) +
                           "_e1.txt");
        if (f_C1.is_open()) {
            f_C1 << "Len, |<nc>| ,Cncnc1\n";
            for (int i = 0; i < Cncnc1_all.size(); i++) {
                f_C1 << Cncnc1_all[i].size() - 1;
                for (int j = 0; j < Cncnc1_all[i].size(); j++) {
                    f_C1 << "," << Cncnc1_all[i][j];
                }
                f_C1 << "\n";
            }
        }
        f_C1.close();

        std::ofstream f_Cr1(filename_Cr.substr(0, filename_Cr.length() - 4) +
                            "_e1.txt");
        if (f_Cr1.is_open()) {
            f_Cr1 << "Len,<r>,Crr1\n";
            for (int i = 0; i < Crr1_all.size(); i++) {
                f_Cr1 << Crr1_all[i].size() - 1;
                for (int j = 0; j < Crr1_all[i].size(); j++) {
                    f_Cr1 << "," << Crr1_all[i][j];
                }
                f_Cr1 << "\n";
            }
        }
        f_Cr1.close();
    }
    if (Ne >= 3) {
        std::ofstream f_C2(filename_C.substr(0, filename_C.length() - 4) +
                           "_e2.txt");
        if (f_C2.is_open()) {
            f_C2 << "Len, |<nc>|,Cncnc2\n";
            for (int i = 0; i < Cncnc2_all.size(); i++) {
                f_C2 << Cncnc2_all[i].size() - 1;
                for (int j = 0; j < Cncnc2_all[i].size(); j++) {
                    f_C2 << "," << Cncnc2_all[i][j];
                }
                f_C2 << "\n";
            }
        }
        f_C2.close();

        std::ofstream f_Cr2(filename_Cr.substr(0, filename_Cr.length() - 4) +
                            "_e2.txt");
        if (f_Cr2.is_open()) {
            f_Cr2 << "Len,<r>,Crr2\n";
            for (int i = 0; i < Crr2_all.size(); i++) {
                f_Cr2 << Crr2_all[i].size() - 1;
                for (int j = 0; j < Crr2_all[i].size(); j++) {
                    f_Cr2 << "," << Crr2_all[i][j];
                }
                f_Cr2 << "\n";
            }
        }
        f_Cr2.close();
    }
}
