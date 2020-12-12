#include "triangulation.h"
#include <cmath>
#include <iostream>
#include <random>
#define PI 3.14159265358979323846

// initialization
triangulation::triangulation(double beta_, int N_, int Ne_, int L_, double d0_,
                             double l0_, double l1_, double kar_, double C0_,
                             double karg_, double lam_, double B_, double Bc_,
                             double tau0_) {
    // system related
    beta = beta_;
    N = N_;
    Ne = Ne_;
    L = L_;
    l0 = l0_; // sigma0 is always 1;
    l1 = l1_;

    // energy related
    kar = kar_;
    C0 = C0_;
    karg = karg_;
    lam = lam_;
    B = B_;
    Bc = Bc_;
    tau0 = tau0_;

    edge_lists.resize(Ne);
    for (int n = 0; n < Ne; n++) {
        edge_lists[n].clear();
    }
    bond_list.clear();

    // new way, assume (N+2)%L != 0
    // if L=sqrt(N), N can't be 700
    if (Ne == 1) {
        init_rhombus_shape(d0_);
    } else if (Ne == 2) {
        init_cylinder_shape(d0_);
    }
    /*
    if (Ne > 1) {
        // two edges
        int i0, i1, i2;
        // i0 = int(N / 2);
        i0 = 3 * L + 1;
        // i0 = int(N / 3);
        i1 = i0 + 1;
        i2 = i0 + L;
        mesh[i0].edge_num = 1;
        push_eneis_back(i0, {1, L});
        edge_lists[1].push_back(i0);
        mesh[i1].edge_num = 1;
        push_eneis_back(i1, {L - 1, -1});
        edge_lists[1].push_back(i1);
        mesh[i2].edge_num = 1;
        push_eneis_back(i2, {-L, -L + 1});
        edge_lists[1].push_back(i2);

        // delete these three edge bond from the bulk bond_list
        delete_bond_list(i0, i1);
        delete_bond_list(i0, i2);
        delete_bond_list(i1, i2);
    }
    if (Ne > 2) {
        // 3 edges
        // std::cout << "not implemented but planing to do!\n";
        // some how it workd for N=200
        int i0, i1, i2;
        // i0 = int(2 * N / 3);
        i0 = 6 * L + 1;
        i1 = i0 + 1;
        i2 = i0 + L;
        mesh[i0].edge_num = 2;
        push_eneis_back(i0, {1, L});
        edge_lists[2].push_back(i0);
        mesh[i1].edge_num = 2;
        push_eneis_back(i1, {L - 1, -1});
        edge_lists[2].push_back(i1);
        mesh[i2].edge_num = 2;
        push_eneis_back(i2, {-L, -L + 1});
        edge_lists[2].push_back(i2);
        // might work? this way?

        // again, delete these three edge bond from the bulk bond_list
        delete_bond_list(i0, i1);
        delete_bond_list(i0, i2);
        delete_bond_list(i1, i2);
    }
    if (Ne > 3) {
        std::cout << "not implemented yet\n";
    }*/

    /*
        // the old way (which works for any N)
        for (int i = 0; i < N; i++) {
            // assign position
            x_n = i % L;
            y_n = i / L;
            mesh[i].pos[0] = d0_ * (x_n + 0.5 * y_n);
            mesh[i].pos[1] = d0_ * 0.5 * std::sqrt(3) * y_n;
            mesh[i].pos[2] = 0;

            // put bonds
            if (y_n == ((N - 1) / L - 1) && x_n > ((N - 1) % L)) {
                edge_list.push_back(i);
                if (x_n == (N - 1) % L + 1) {
                    push_neis_back(i, {L - 1, -1, -L, -L + 1, 1});
                    push_eneis_back(i, {L - 1, 1});
                } else if (x_n == L - 1) {
                    push_neis_back(i, {-1, -L});
                    push_eneis_back(i, {-1, -L});
                } else {
                    push_neis_back(i, {-1, -L, -L + 1, 1});
                    push_eneis_back(i, {-1, 1});
                }
            } else if (y_n == (N - 1) / L) {
                edge_list.push_back(i);
                if (x_n == (N - 1) % L) {
                    if (N % L != 0) {
                        // the corner on the middle of top line exist
                        if (x_n == 0) {
                            // only one bead on the top line
                            push_neis_back(i, {-L, -L + 1});
                            push_eneis_back(i, {-L, -L + 1});
                        } else {
                            push_neis_back(i, {-1, -L, -L + 1});
                            push_eneis_back(i, {-1, -L + 1});
                        }
                    } else {
                        // no such corner
                        push_neis_back(i, {-1, -L});
                        push_eneis_back(i, {-1, -L});
                    }
                } else if (x_n == 0) {
                    push_neis_back(i, {-L, -L + 1, 1});
                    push_eneis_back(i, {-L, 1});
                } else {
                    push_neis_back(i, {-1, -L, -L + 1, +1});
                    push_eneis_back(i, {-1, 1});
                }
            } else if (x_n == 0) {
                edge_list.push_back(i);
                if (y_n == 0) {
                    push_neis_back(i, {1, L});
                    push_eneis_back(i, {1, L});
                } else {
                    push_neis_back(i, {-L, -L + 1, +1, L});
                    push_eneis_back(i, {-L, L});
                }
            } else if (x_n == L - 1) {
                edge_list.push_back(i);
                if (y_n == 0) {
                    push_neis_back(i, {L, L - 1, -1});
                    push_eneis_back(i, {L, -1});
                } else {
                    push_neis_back(i, {L, L - 1, -1, -L});
                    push_eneis_back(i, {L, -L});
                }
            } else if (y_n == 0) {
                edge_list.push_back(i);
                push_neis_back(i, {1, L, L - 1, -1});
                push_eneis_back(i, {1, -1});
            } else {
                push_neis_back(i, {1, L, L - 1, -1, -L, -L + 1});
            }
        }
    */
    /*
        // sort edge direction and edge list.
        // initla
        // edge_nei[0]---edge_list[0]---edge_nei[1]
        int edge_i = edge_list[0];             // edge initial
        int edge_c = edge_i;                   // current sorted bead
        int edge_n = mesh[edge_c].edge_nei[1]; // next edge bead
        edge_list.clear();
        edge_list.push_back(edge_c);
        while (edge_n != edge_i) {
            if (mesh[edge_n].edge_nei[0] != edge_c) {
                mesh[edge_n].edge_nei[1] = mesh[edge_n].edge_nei[0];
                mesh[edge_n].edge_nei[0] = edge_c;
                std::cout << "it happens\n";
            } // sorted now
            edge_c = edge_n;
            edge_n = mesh[edge_c].edge_nei[1];
            edge_list.push_back(edge_c);
        }
        // so, it's actually not necessary since it's already sorted when
        generating the edge_nei manually
    */

    // add nei2flip list
    // instead just maintain bond_list that contains all bond in bulk, twice
    // somehow
    /*
        int ind_o, i_nei_a, i_nei_b, ind_a, ind_b;
        for (int ind_i = 0; ind_i < mesh.size(); ind_i++) {
            if (mesh[ind_i].nei.size() < 4) {
                continue; // ind_i has to have more than 3 neighbours
            }
            for (int i_nei_o = 0; i_nei_o < mesh[ind_i].nei.size();
       i_nei_o++) {
                // enumerate ind_i's neighbours
                ind_o = mesh[ind_i].nei[i_nei_o]; // the current neighbour
                if (list_a_nei_b(mesh[ind_i].edge_nei, ind_o) == -1 &&
                    mesh[ind_o].nei.size() >= 4) {
                    i_nei_a = (i_nei_o + 1) % mesh[ind_i].nei.size();
                    i_nei_b = (i_nei_o - 1 + mesh[ind_i].nei.size()) %
                              mesh[ind_i].nei.size();
                    ind_a = mesh[ind_i].nei[i_nei_a];
                    ind_b = mesh[ind_i].nei[i_nei_b];
                    // clearly, a-b can't be connected
                    if (list_a_nei_b(mesh[ind_a].nei, ind_b) == -1) {
                        mesh[ind_i].nei2flip.push_back(i_nei_o);
                    }
                }
            }
            if (mesh[ind_i].nei2flip.size() != 0) {
                v2flip_list.push_back(ind_i);
            }
        }
    */
    // checking on bond_list
    /*
    for (int i = 0; i < bond_list.size(); i++) {
        int ind_i = bond_list[i].first;
        int ind_o = bond_list[i].second;
        if (list_a_nei_b(mesh[ind_i].edge_nei, ind_o) != -1) {
            std::cout << "an edge bond?!\n";
        }
        if (list_a_nei_b(mesh[ind_i].nei, ind_o) == -1) {
            std::cout << "i-o not connected at the beginning?!" << L <<
    "\n"; std::cout << "ind_i/ind_o" << ind_i << "/" << ind_o << "\n"; for
    (int j = 0; j < mesh[ind_i].nei.size(); j++) { std::cout <<
    mesh[ind_i].nei[j] << ",";
            }
            std::cout << "\n";
        }
    }
    */
    // set inital observable value
    E = 0;
    // Le = 0;
    Les.resize(Ne);
    for (int e = 0; e < Ne; e++) {
        Les[e] = 0;
    }
    IdA = 0;
    I2H = 0;
    I2H2 = 0;
    Ikg = 0;
    Ik2 = 0;
    Itau = 0;
    // Rg2 = 0;
    for (int i = 0; i < mesh.size(); i++) {
        mesh[i].ds = ds_m(i);
        mesh[i].dAn2H = dAn2H_m(i);
        // mesh[i].dA2H2 = dA2H2_m(i);
        mesh[i].dskg = dskg_m(i);
        mesh[i].dsk2 = dsk2_m(i);
        mesh[i].nc = nc_m(i);
        mesh[i].dstau = dstau_m(i);
        mesh[i].dE = dE_m(i);
        // mesh[i].dRg2 = dRg2_m(i);
        IdA += mesh[i].dAn2H[0];
        I2H += mesh[i].dAn2H[1] * mesh[i].dAn2H[0];
        I2H2 += mesh[i].dAn2H[1] * mesh[i].dAn2H[1] * mesh[i].dAn2H[0];
        if (mesh[i].edge_num != -1) {
            Les[mesh[i].edge_num] += mesh[i].ds;
        }
        Ikg += mesh[i].dskg;
        Ik2 += mesh[i].dsk2;
        Itau += mesh[i].dstau;
        // Rg2 += mesh[i].dRg2;
        // if (mesh[i].edge_nei.size() != 0) {
        // Itau2 += mesh[i].dstau * mesh[i].dstau / mesh[i].ds;
        //}
        E += mesh[i].dE;
    }
    // std::cout << "Le vs Le_m = " << d0_ * edge_list.size() << " vs " <<
    // Le; std::cout << "after intialization\n"; std::cout << "E = " << E <<
    // "\n";

    // set random number generators
    std::random_device rd;
    std::mt19937 gen_set(rd());
    std::uniform_int_distribution<> rand_pos_set(0,
                                                 mesh.size() - 1); // random pos
    std::uniform_real_distribution<> rand_uni_set(0, 1);
    gen = gen_set;
    rand_pos = rand_pos_set;
    rand_uni = rand_uni_set;

} // end of triangulation

void triangulation::init_rhombus_shape(double d0_) {
    int x_n, y_n; // position of the vertex in the two vector coordinate
    mesh.resize(N);
    if ((N + 2) % L == 1 || (N + 2) % L >= (L - 1)) {
        std::cout << "invalid N!\n";
    }
    for (int i = 0; i < N; i++) {
        // assign position
        if (i < N - N % L - 2) {
            x_n = (i + 1) % L;
            y_n = (i + 1) / L;
        } else {
            // ignore upper-right corner due to #-of-nei limit
            x_n = (i + 2) % L;
            y_n = (i + 2) / L;
        }

        mesh[i].pos[0] = d0_ * (x_n + 0.5 * y_n);
        mesh[i].pos[1] = d0_ * 0.5 * std::sqrt(3) * y_n;
        mesh[i].pos[2] = 0;

        // put bonds
        if (x_n == L - 1) {
            // right line
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (y_n == 0) {
                // lower-right corner
                push_neis_back(i, {L, L - 1, -1});
                push_eneis_back(i, {L, -1});
                push_bneis_list(i, {L - 1});
            } else if (y_n == (N + 1) / L - 2) {
                // upper-right corner
                push_neis_back(i, {L - 1, -1, -L});
                push_eneis_back(i, {L - 1, -L});
                push_bneis_list(i, {-1});
            } else {
                push_neis_back(i, {L, L - 1, -1, -L});
                push_eneis_back(i, {L, -L});
                push_bneis_list(i, {L - 1, -1});
            }
        } else if (y_n == (N + 1) / L - 1 && x_n != 0) {
            if (x_n > (N + 1) % L) {
                edge_lists[0].push_back(i);
                mesh[i].edge_num = 0;
                if (x_n == L - 2) {
                    if ((N + 2) % L < (L - 2)) {
                        push_neis_back(i, {-1, -L, -L + 1});
                        push_eneis_back(i, {-1, -L + 1});
                        push_bneis_list(i, {-L});
                    } else {
                        push_neis_back(i, {L - 2, -1, -L, -L + 1});
                        push_eneis_back(i, {L - 2, -L + 1});
                        push_bneis_list(i, {-1, -L});
                    }

                } else if (x_n == (N + 1) % L + 1) {
                    push_neis_back(i, {L - 2, -1, -L, -L + 1, 1});
                    push_eneis_back(i, {L - 2, 1});
                    push_bneis_list(i, {-1, -L, -L + 1});
                } else {
                    push_neis_back(i, {-1, -L, -L + 1, 1});
                    push_eneis_back(i, {-1, 1});
                    push_bneis_list(i, {-L, -L + 1});
                }
            } else {
                // in the bulk
                mesh[i].edge_num = -1;
                push_neis_back(i, {1, L - 1, L - 2, -1, -L, -L + 1});
                push_bneis_list(i, {1, L - 1, L - 2, -1, -L, -L + 1});
            }
        } else if (y_n == (N + 1) / L) {
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (x_n == (N + 1) % L) {
                push_neis_back(i, {-1, -L + 1, -L + 2});
                push_eneis_back(i, {-1, -L + 2});
                push_bneis_list(i, {-L + 1});
            } else if (x_n == 0) {
                push_neis_back(i, {-L + 1, -L + 2, 1});
                push_eneis_back(i, {-L + 1, 1});
                push_bneis_list(i, {-L + 2});
            } else {
                push_neis_back(i, {-1, -L + 1, -L + 2, 1});
                push_eneis_back(i, {-1, 1});
                push_bneis_list(i, {-L + 1, -L + 2});
            }
        } else if (x_n == 0) {
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (y_n == 1) {
                push_neis_back(i, {-L + 1, 1, L});
                push_eneis_back(i, {-L + 1, L});
                push_bneis_list(i, {1});
            } else if (y_n == (N + 1) / L - 1) {
                push_neis_back(i, {-L, -L + 1, 1, L - 1});
                push_eneis_back(i, {-L, L - 1});
                push_bneis_list(i, {-L + 1, 1});
            } else {
                push_neis_back(i, {-L, -L + 1, 1, L});
                push_eneis_back(i, {-L, L});
                push_bneis_list(i, {-L + 1, 1});
            }
        } else if (y_n == 0) {
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (x_n == 1) {
                push_neis_back(i, {1, L, L - 1});
                push_eneis_back(i, {1, L - 1});
                push_bneis_list(i, {L});
            } else {
                push_neis_back(i, {1, L, L - 1, -1});
                push_eneis_back(i, {1, -1});
                push_bneis_list(i, {L, L - 1});
            }
        } else {
            mesh[i].edge_num = -1; // in the bulk
            push_neis_back(i, {1, L, L - 1, -1, -L, -L + 1});
            push_bneis_list(i, {1, L, L - 1, -1, -L, -L + 1});
        }
    }
}
void triangulation::init_cylinder_shape(double d0_) {
    int x_n, y_n; // position of the vertex in the two vector coordinate
    // cylinder initial shape
    Lr = N / L;
    double R = d0_ / (2 * std::sin(PI / Lr));
    mesh.resize(N);
    for (int i = 0; i < N; i++) {
        // assign position
        x_n = i % Lr;
        y_n = i / Lr;
        mesh[i].pos[0] = R * std::cos(2 * PI * (x_n + 0.5 * y_n) / Lr);
        mesh[i].pos[1] = R * std::sin(2 * PI * (x_n + 0.5 * y_n) / Lr);
        mesh[i].pos[2] = d0_ * 0.5 * std::sqrt(3) * y_n;

        // put bonds
        if (y_n == 0) {
            // lower edge
            mesh[i].edge_num = 0;
            edge_lists[0].push_back(i);
            if (x_n == 0) {
                push_neis_back(i, {1, Lr, 2 * Lr - 1, Lr - 1});
                push_eneis_back(i, {1, Lr - 1});
                push_bneis_list(i, {Lr, 2 * Lr - 1});
            } else if (x_n == Lr - 1) {
                push_neis_back(i, {-Lr + 1, Lr, Lr - 1, -1});
                push_eneis_back(i, {-Lr + 1, -1});
                push_bneis_list(i, {Lr, Lr - 1});
            } else {
                push_neis_back(i, {1, Lr, Lr - 1, -1});
                push_eneis_back(i, {1, -1});
                push_bneis_list(i, {Lr, Lr - 1});
            }
        } else if (y_n == L - 1) {
            // upper edge
            mesh[i].edge_num = 1;
            edge_lists[1].push_back(i);
            if (x_n == 0) {
                push_neis_back(i, {Lr - 1, -Lr, -Lr + 1, 1});
                push_eneis_back(i, {Lr - 1, 1});
                push_bneis_list(i, {-Lr, -Lr + 1});
            } else if (x_n == Lr - 1) {
                push_neis_back(i, {-1, -Lr, -2 * Lr + 1, -Lr + 1});
                push_eneis_back(i, {-1, -Lr + 1});
                push_bneis_list(i, {-Lr, -2 * Lr + 1});
            } else {
                push_neis_back(i, {-1, -Lr, -Lr + 1, 1});
                push_eneis_back(i, {-1, 1});
                push_bneis_list(i, {-Lr, -Lr + 1});
            }
        } else {
            // internal
            mesh[i].edge_num = -1;
            if (x_n == 0) {
                push_neis_back(i, {1, Lr, 2 * Lr - 1, Lr - 1, -Lr, -Lr + 1});
                push_bneis_list(i, {1, Lr, 2 * Lr - 1, Lr - 1, -Lr, -Lr + 1});
            } else if (x_n == Lr - 1) {
                push_neis_back(i, {-Lr + 1, Lr, Lr - 1, -1, -Lr, -2 * Lr + 1});
                push_bneis_list(i, {-Lr + 1, Lr, Lr - 1, -1, -Lr, -2 * Lr + 1});
            } else {
                push_neis_back(i, {1, Lr, Lr - 1, -1, -Lr, -Lr + 1});
                push_bneis_list(i, {1, Lr, Lr - 1, -1, -Lr, -Lr + 1});
            }
        }
    }
}

void triangulation::push_neis_back(int i, std::vector<int> nei_dist) {
    for (int j = 0; j < nei_dist.size(); j++) {
        mesh[i].nei.push_back(i + nei_dist[j]);
    }
}
void triangulation::push_eneis_back(int i, std::vector<int> enei_dist) {
    for (int j = 0; j < enei_dist.size(); j++) {
        mesh[i].edge_nei.push_back(i + enei_dist[j]);
    }
}
void triangulation::push_bneis_list(int i, std::vector<int> bnei_dist) {
    std::pair<int, int> bond;
    for (int j = 0; j < bnei_dist.size(); j++) {
        bond.first = i;
        bond.second = i + bnei_dist[j];
        bond_list.push_back(bond);
    }
}
void triangulation::delete_bond_list(int ind_i, int ind_j) {
    std::pair<int, int> bond0, bond1;
    bond0.first = ind_i;
    bond0.second = ind_j;
    bond1.first = ind_j;
    bond1.second = ind_i;
    for (int i = 0; i < bond_list.size(); i++) {
        if (bond_list[i] == bond0 || bond_list[i] == bond1) {
            bond_list.erase(bond_list.begin() + i);
            i--;
        }
    }
}