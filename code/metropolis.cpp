#include "triangulation.h"
#include <ctime>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>

int print_vec(std::vector<int> vec) {
    for (int j = 0; j < vec.size(); j++) {
        std::cout << vec[j] << ",";
    }
    std::cout << "\n";
    return 0;
}

int triangulation::vertex_metropolis(double delta_s) {
    int index, nei_ind;
    std::vector<double> dsr_old; // ds related
    std::vector<double> Lers_old(Ne, 0);
    std::vector<double> Lers_new(Ne, 0);
    std::vector<std::vector<double>> dAn2Hr_old; // dAn2H related
    double IdAr_old, IdAr_new;
    double I2Hr_old, I2Hr_new;
    double I2H2r_old, I2H2r_new;
    std::vector<double> dskgr_old; // dskg related
    double Ikgr_old, Ikgr_new;
    std::vector<double> dsk2r_old; // dsk2 related
    double Ik2r_old, Ik2r_new;
    std::vector<double> dstaur_old; // dstau related
    double Itaur_old, Itaur_new;
    std::vector<double> dEr_old; // dE related
    double Er_old, Er_new;       // summation of related dr
    double distance_buff;
    std::set<int, std::greater<int>>
        afc_beads; // beads that energy are affacted by the update
    std::set<int, std::greater<int>> afc_ebeads; // accordingly, but on the edge
    std::set<int, std::greater<int>>::iterator itr;
    do {
        index = rand_pos(gen);
    } while (list_a_nei_b(fixed_beads, index) != -1);
    // old local energy summation
    // find affacted beads
    afc_beads.insert(index);
    for (int j = 0; j < mesh[index].nei.size(); j++) {
        // std::cout<<"Er_old="<<Er_old<<"\n";
        nei_ind = mesh[index].nei[j];
        if (mesh[nei_ind].edge_nei.size() != 0) {
            afc_ebeads.insert(nei_ind);
        }
        afc_beads.insert(nei_ind);
    }
    for (itr = afc_ebeads.begin(); itr != afc_ebeads.end(); itr++) {
        for (int j = 0; j < mesh[*itr].edge_nei.size(); j++) {
            afc_beads.insert(mesh[*itr].edge_nei[j]);
        }
    }
    // store observable of affected beads
    dAn2Hr_old.clear();
    IdAr_old = 0;
    I2Hr_old = 0;
    I2H2r_old = 0;
    dsr_old.clear();
    for (int e = 0; e < Ne; e++) {
        Lers_old[e] = 0;
    }
    dskgr_old.clear();
    Ikgr_old = 0;
    dsk2r_old.clear();
    Ik2r_old = 0;
    dstaur_old.clear();
    Itaur_old = 0;
    dEr_old.clear();
    Er_old = 0;
    for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
        // kappa term
        dAn2Hr_old.push_back(mesh[*itr].dAn2H);
        IdAr_old += mesh[*itr].dAn2H[0];
        I2Hr_old += mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1];
        I2H2r_old +=
            mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1] * mesh[*itr].dAn2H[1];
        // lambda term
        dsr_old.push_back(mesh[*itr].ds);
        Lers_old[mesh[*itr].edge_num] += mesh[*itr].ds;
        // karg term
        dskgr_old.push_back(mesh[*itr].dskg);
        Ikgr_old += mesh[*itr].dskg;
        // B term
        dsk2r_old.push_back(mesh[*itr].dsk2);
        Ik2r_old += mesh[*itr].dsk2;
        // Bc term
        dstaur_old.push_back(mesh[*itr].dstau);
        Itaur_old += mesh[*itr].dstau;
        // energy
        dEr_old.push_back(mesh[*itr].dE);
        Er_old += mesh[*itr].dE;
    }

    // Monte Carlo update proposal
    std::vector<double> delta_pos{0, 0, 0};
    for (int k = 0; k < 3; k++) {
        delta_pos[k] = 2 * delta_s * rand_uni(gen) - delta_s;
        mesh[index].pos[k] += delta_pos[k];
    }

    // hard bead potential between all beads
    for (int nei_ind = 0; nei_ind < mesh.size(); nei_ind++) {
        if ((std::abs(mesh[index].pos[0] - mesh[nei_ind].pos[0])) < 1 &&
            (std::abs(mesh[index].pos[1] - mesh[nei_ind].pos[1])) < 1 &&
            (std::abs(mesh[index].pos[2] - mesh[nei_ind].pos[2])) < 1) {
            distance_buff = distance(index, nei_ind);
            if ((distance_buff < 1) && (nei_ind != index)) {
                for (int k = 0; k < 3; k++) {
                    mesh[index].pos[k] -= delta_pos[k];
                } // return previous position
                return 0;
            }
        }
    }

    // tether potential between linked beads
    for (int j = 0; j < mesh[index].nei.size(); j++) {
        nei_ind = mesh[index].nei[j];
        distance_buff = distance(index, nei_ind);
        if (mesh[index].edge_nei.size() != 0 &&
            (nei_ind == mesh[index].edge_nei[0] ||
             nei_ind == mesh[index].edge_nei[1])) {
            if (distance_buff >= l1) {
                for (int k = 0; k < 3; k++) {
                    mesh[index].pos[k] -= delta_pos[k];
                } // return previous position
                return 0;
            }
        } else if (distance_buff >= l0) {
            for (int k = 0; k < 3; k++) {
                mesh[index].pos[k] -= delta_pos[k];
            } // return previous position
            return 0;
        }
    }

    // after-update observables
    IdAr_new = 0;
    I2Hr_new = 0;
    I2H2r_new = 0;
    for (int e = 0; e < Ne; e++) {
        Lers_new[e] = 0;
    }
    Ikgr_new = 0;
    Ik2r_new = 0;
    Itaur_new = 0;
    Er_new = 0;
    for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
        // kappa term
        mesh[*itr].dAn2H = dAn2H_m(*itr);
        IdAr_new += mesh[*itr].dAn2H[0];
        I2Hr_new += mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1];
        I2H2r_new +=
            mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1] * mesh[*itr].dAn2H[1];
        // lambda term (needed by other edge terms)
        mesh[*itr].ds = ds_m(*itr);
        Lers_new[mesh[*itr].edge_num] += mesh[*itr].ds;
        // karg term
        mesh[*itr].dskg = dskg_m(*itr);
        Ikgr_new += mesh[*itr].dskg;
        // B term
        mesh[*itr].dsk2 = dsk2_m(*itr);
        Ik2r_new += mesh[*itr].dsk2;
        // Bc term
        mesh[*itr].dstau = dstau_m(*itr);
        Itaur_new += mesh[*itr].dstau;
        // energy
        mesh[*itr].dE = dE_m(*itr);
        Er_new += mesh[*itr].dE;
    }
    // [Metropolis]
    // std::cout << "Er_new-Er_old=" << Er_new - Er_old << "\n";
    if (rand_uni(gen) <= std::exp(-beta * (Er_new - Er_old))) {
        // [accepted]
        IdA += IdAr_new - IdAr_old;
        I2H += I2Hr_new - I2Hr_old;
        I2H2 += I2H2r_new - I2H2r_old;
        for (int e = 0; e < Ne; e++) {
            Les[e] += Lers_new[e] - Lers_old[e];
        }
        Ikg += Ikgr_new - Ikgr_old;
        Ik2 += Ik2r_new - Ik2r_old;
        Itau += Itaur_new - Itaur_old;
        E += Er_new - Er_old;
        return 1;
    } else {
        // [rejected]
        for (int k = 0; k < 3; k++) {
            // return previous position
            mesh[index].pos[k] -= delta_pos[k];
        }
        int count = 0;
        for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
            mesh[*itr].ds = dsr_old[count];
            mesh[*itr].dAn2H = dAn2Hr_old[count];
            mesh[*itr].dskg = dskgr_old[count];
            mesh[*itr].dsk2 = dsk2r_old[count];
            mesh[*itr].dstau = dstaur_old[count];
            mesh[*itr].dE = dEr_old[count];
            count++;
        }

        return 0;
    }
}

int triangulation::bond_metropolis()
/*
           - (o) -
         -    |       -
       -      |       - (b)
    (a)       |     -
        -     |   -
           - (i)

    to:
              - (o) -
         -            -
       -      - - - - - (b)
    (a) - -         -
        -         -
           - (i)

*/
{
    // int v2flip_ind_i, v2flip_ind_i_nei_o;
    int bondlist_ind_i;
    int ind_i, ind_a, ind_b, ind_o;
    int i_nei_o, i_nei_a, i_nei_b; // mesh[ind_i].nei[i_nei_o] = ind_o;
    int o_nei_i, o_nei_a, o_nei_b;
    int a_nei_i, a_nei_o, a_nei_b;
    int b_nei_i, b_nei_o, b_nei_a;

    // std::vector<int> i_nei, o_nei, a_nei, b_nei;
    // vertex i_v, o_v, a_v, b_v;

    std::set<int, std::greater<int>>
        afc_beads; // beads that energy are affacted by the update
    std::set<int, std::greater<int>> afc_ebeads; // accordingly, but on the edge
    std::set<int, std::greater<int>>::iterator itr;

    // observables
    std::vector<double> dsr_old; // ds related
    std::vector<double> Lers_old(Ne, 0);
    std::vector<double> Lers_new(Ne, 0);
    std::vector<std::vector<double>> dAn2Hr_old; // dAn2H related
    double IdAr_old, IdAr_new;
    double I2Hr_old, I2Hr_new;
    double I2H2r_old, I2H2r_new;
    std::vector<double> dskgr_old; // dskg related
    double Ikgr_old, Ikgr_new;
    std::vector<double> dsk2r_old; // dsk2 related
    double Ik2r_old, Ik2r_new;
    std::vector<double> dstaur_old; // dstau related
    double Itaur_old, Itaur_new;
    std::vector<double> dEr_old; // dE related
    double Er_old, Er_new;       // summation of related dr

    // randomize ind_i and i_nei_io
    // renew_v2flip_list();

    // throw away v2flip method
    // if (v2flip_list.size() == 0) {
    // return 0;
    //}
    if (bond_list.size() == 0) {
        return 0;
    }
    bondlist_ind_i = int(bond_list.size() * rand_uni(gen));
    // v2flip_ind_i = int(v2flip_list.size() * rand_uni(gen));
    // ind_i = v2flip_list[v2flip_ind_i];
    ind_i = bond_list[bondlist_ind_i].first;
    ind_o = bond_list[bondlist_ind_i].second;
    /*
    if (list_a_nei_b(mesh[ind_i].nei, ind_o) == -1) {
        std::cout << "i-o not connected?!\n";
    }
    if (list_a_nei_b(mesh[ind_i].edge_nei, ind_o) != -1) {
        std::cout << "an edge bond?!\n";
    }
    */
    // v2flip_ind_i_nei_o = int(mesh[ind_i].nei2flip.size() * rand_uni(gen));
    // this actually gonna break the symmetry, since bead with less bond will be
    // more likely to have their bond flipped
    // i_nei_o = mesh[ind_i].nei2flip[v2flip_ind_i_nei_o];
    i_nei_o = list_a_nei_b(mesh[ind_i].nei, ind_o);
    // ind_o = mesh[ind_i].nei[i_nei_o];
    i_nei_a = (i_nei_o + 1) % mesh[ind_i].nei.size();
    i_nei_b = (i_nei_o - 1 + mesh[ind_i].nei.size()) % mesh[ind_i].nei.size();
    ind_a = mesh[ind_i].nei[i_nei_a];
    ind_b = mesh[ind_i].nei[i_nei_b];
    if (list_a_nei_b(mesh[ind_a].nei, ind_b) != -1) {
        // a and b can't be connected
        return 0;
    }

    // check the a-b distance
    if (distance(ind_a, ind_b) > l0) {
        return 0;
    }
    // check # of nei
    // can't be greater than 9, otherwise will have more than 9 after flip
    if (mesh[ind_a].nei.size() >= 9 || mesh[ind_b].nei.size() >= 9) {
        return 0;
    }
    // can't be less than 3
    // both has to have at least 4 neighbours, after flip, will have at least
    // 3;
    if (mesh[ind_i].nei.size() <= 3 || mesh[ind_o].nei.size() <= 3) {
        return 0;
    }
    // this is for debuging

    // i_nei = mesh[ind_i].nei;
    // o_nei = mesh[ind_o].nei;
    // a_nei = mesh[ind_a].nei;
    // b_nei = mesh[ind_b].nei;

    // i_v = mesh[ind_i];
    // o_v = mesh[ind_o];
    // a_v = mesh[ind_a];
    // b_v = mesh[ind_b];

    // find the rest of bonds
    o_nei_i = list_a_nei_b(mesh[ind_o].nei, ind_i);
    o_nei_a = list_a_nei_b(mesh[ind_o].nei, ind_a);
    o_nei_b = list_a_nei_b(mesh[ind_o].nei, ind_b);
    a_nei_i = list_a_nei_b(mesh[ind_a].nei, ind_i);
    a_nei_o = list_a_nei_b(mesh[ind_a].nei, ind_o);
    b_nei_i = list_a_nei_b(mesh[ind_b].nei, ind_i);
    b_nei_o = list_a_nei_b(mesh[ind_b].nei, ind_o);

    /*
    if (((a_nei_i - a_nei_o + mesh[ind_a].nei.size()) %
             mesh[ind_a].nei.size() !=
         1) &&
        ((a_nei_o - a_nei_i + mesh[ind_a].nei.size()) %
             mesh[ind_a].nei.size() !=
         1)) {
        std::cout << "i-o not neib in a anei=\n";
    }
    if (((b_nei_i - b_nei_o + mesh[ind_b].nei.size()) %
             mesh[ind_b].nei.size() !=
         1) &&
        ((b_nei_o - b_nei_i + mesh[ind_b].nei.size()) %
             mesh[ind_b].nei.size() !=
         1)) {
        std::cout << "i-o not neib in b bnei=\n";
    }
    */

    // energy before bond update
    // find related beads
    afc_beads.insert(ind_i);
    afc_beads.insert(ind_o);
    afc_beads.insert(ind_a);
    afc_beads.insert(ind_b);
    if (mesh[ind_i].edge_nei.size() != 0) {
        for (int j = 0; j < mesh[ind_i].edge_nei.size(); j++) {
            afc_beads.insert(mesh[ind_i].edge_nei[j]);
        }
    }
    if (mesh[ind_o].edge_nei.size() != 0) {
        for (int j = 0; j < mesh[ind_o].edge_nei.size(); j++) {
            afc_beads.insert(mesh[ind_o].edge_nei[j]);
        }
    }
    if (mesh[ind_a].edge_nei.size() != 0) {
        for (int j = 0; j < mesh[ind_a].edge_nei.size(); j++) {
            afc_beads.insert(mesh[ind_a].edge_nei[j]);
        }
    }
    if (mesh[ind_b].edge_nei.size() != 0) {
        for (int j = 0; j < mesh[ind_b].edge_nei.size(); j++) {
            afc_beads.insert(mesh[ind_b].edge_nei[j]);
        }
    }
    // store observable of affected beads
    dAn2Hr_old.clear();
    IdAr_old = 0;
    I2Hr_old = 0;
    I2H2r_old = 0;
    dsr_old.clear();
    for (int e = 0; e < Ne; e++) {
        Lers_old[e] = 0;
    }
    dskgr_old.clear();
    Ikgr_old = 0;
    dsk2r_old.clear();
    Ik2r_old = 0;
    dstaur_old.clear();
    Itaur_old = 0;
    dEr_old.clear();
    Er_old = 0;
    for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
        // kappa term
        dAn2Hr_old.push_back(mesh[*itr].dAn2H);
        IdAr_old += mesh[*itr].dAn2H[0];
        I2Hr_old += mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1];
        I2H2r_old +=
            mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1] * mesh[*itr].dAn2H[1];
        // lambda term
        dsr_old.push_back(mesh[*itr].ds);
        Lers_old[mesh[*itr].edge_num] += mesh[*itr].ds;
        // karg term
        dskgr_old.push_back(mesh[*itr].dskg);
        Ikgr_old += mesh[*itr].dskg;
        // B term
        dsk2r_old.push_back(mesh[*itr].dsk2);
        Ik2r_old += mesh[*itr].dsk2;
        // Bc term
        dstaur_old.push_back(mesh[*itr].dstau);
        Itaur_old += mesh[*itr].dstau;
        // energy
        dEr_old.push_back(mesh[*itr].dE);
        Er_old += mesh[*itr].dE;
    }

    // [update]
    // remove i-o bond
    mesh[ind_i].nei.erase(mesh[ind_i].nei.begin() + i_nei_o);
    mesh[ind_o].nei.erase(mesh[ind_o].nei.begin() + o_nei_i);
    // one of i_nei_a, i_nei_b changes
    // one of o_nei_a, o_nei_b changes

    // add a-b bond
    if ((a_nei_i - a_nei_o + mesh[ind_a].nei.size()) % mesh[ind_a].nei.size() ==
        1) {
        mesh[ind_a].nei.insert(mesh[ind_a].nei.begin() + a_nei_i, ind_b);
        a_nei_b = a_nei_i;
        a_nei_i++;
    } else if ((a_nei_o - a_nei_i + mesh[ind_a].nei.size()) %
                   mesh[ind_a].nei.size() ==
               1) {
        mesh[ind_a].nei.insert(mesh[ind_a].nei.begin() + a_nei_o, ind_b);
        a_nei_b = a_nei_o;
        a_nei_o++;
    } else {
        std::cout << "i-o not next in a_nei!"
                  << "a_nei_i,a_nei_o=" << a_nei_i << "," << a_nei_o
                  << " a.nei.size()=" << mesh[ind_a].nei.size() << "\n";
        std::cout << "miao?";
    } // error code 2
    if ((b_nei_i - b_nei_o + mesh[ind_b].nei.size()) % mesh[ind_b].nei.size() ==
        1) {
        mesh[ind_b].nei.insert(mesh[ind_b].nei.begin() + b_nei_i, ind_a);
        b_nei_a = b_nei_i;
        b_nei_i++;
    } else if ((b_nei_o - b_nei_i + mesh[ind_b].nei.size()) %
                   mesh[ind_b].nei.size() ==
               1) {
        mesh[ind_b].nei.insert(mesh[ind_b].nei.begin() + b_nei_o, ind_a);
        b_nei_a = b_nei_o;
        b_nei_o++;
    } else {
        std::cout << "i-o not next in b_nei!"
                  << "b_nei_i,b_nei_o=" << b_nei_i << "," << b_nei_o
                  << " b.nei.size()=" << mesh[ind_b].nei.size() << "\n";
        std::cout << "huh?";
    } // error code 2

    // after-update observables
    IdAr_new = 0;
    I2Hr_new = 0;
    I2H2r_new = 0;
    for (int e = 0; e < Ne; e++) {
        Lers_new[e] = 0;
    }
    Ikgr_new = 0;
    Ik2r_new = 0;
    Itaur_new = 0;
    Er_new = 0;
    for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
        // kappa term
        mesh[*itr].dAn2H = dAn2H_m(*itr);
        IdAr_new += mesh[*itr].dAn2H[0];
        I2Hr_new += mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1];
        I2H2r_new +=
            mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1] * mesh[*itr].dAn2H[1];
        // lambda term (needed by other edge terms)
        mesh[*itr].ds = ds_m(*itr);
        Lers_new[mesh[*itr].edge_num] += mesh[*itr].ds;
        // karg term
        mesh[*itr].dskg = dskg_m(*itr);
        Ikgr_new += mesh[*itr].dskg;
        // B term
        mesh[*itr].dsk2 = dsk2_m(*itr);
        Ik2r_new += mesh[*itr].dsk2;
        // Bc term
        mesh[*itr].dstau = dstau_m(*itr);
        Itaur_new += mesh[*itr].dstau;
        // energy
        mesh[*itr].dE = dE_m(*itr);
        Er_new += mesh[*itr].dE;
    }

    // [Metropolis]
    if (rand_uni(gen) <= std::exp(-beta * (Er_new - Er_old))) {
        // [accepted]
        E += Er_new - Er_old;
        for (int e = 0; e < Ne; e++) {
            Les[e] += Lers_new[e] - Lers_old[e];
            if (Les[e] < 0) {
                std::cout << "Les[e]<0\n";
            }
        }
        IdA += IdAr_new - IdAr_old;
        I2H += I2Hr_new - I2Hr_old;
        I2H2 += I2H2r_new - I2H2r_old;
        Ikg += Ikgr_new - Ikgr_old;
        Ik2 += Ik2r_new - Ik2r_old;
        Itau += Itaur_new - Itaur_old;

        // update bond_list
        delete_bond_list(ind_i, ind_o);
        std::pair<int, int> bond0, bond1;
        bond0.first = ind_a;
        bond0.second = ind_b;
        bond1.first = ind_b;
        bond1.second = ind_a;
        bond_list.push_back(bond0);
        bond_list.push_back(bond1);

        // update bond 2 flip
        /*
        update_v2flip_local(ind_i);
        update_v2flip_group(mesh[ind_i].nei);
        update_v2flip_local(ind_o);
        update_v2flip_group(mesh[ind_o].nei);
        update_v2flip_local(ind_a);
        update_v2flip_group(mesh[ind_a].nei);
        update_v2flip_local(ind_b);
        update_v2flip_group(mesh[ind_b].nei);
        */
        // check duplication
        /*
        check_duplication(ind_i);
        check_duplication(ind_o);
        check_duplication(ind_a);
        check_duplication(ind_b);

        if (std::abs(Er_new - Er_old) > 1000) {
            std::cout << "energy change (bond): " << (Er_new - Er_old) << "\n";
        }
        */
        return 1;
    } else {
        // [rejected]
        // put everything back (so sad T.T)
        mesh[ind_a].nei.erase(mesh[ind_a].nei.begin() + a_nei_b);
        mesh[ind_b].nei.erase(mesh[ind_b].nei.begin() + b_nei_a);
        mesh[ind_i].nei.insert(mesh[ind_i].nei.begin() + i_nei_o, ind_o);
        mesh[ind_o].nei.insert(mesh[ind_o].nei.begin() + o_nei_i, ind_i);

        // if (mesh[ind_i].nei != i_nei || mesh[ind_o].nei != o_nei ||
        //   mesh[ind_a].nei != a_nei || mesh[ind_b].nei != b_nei) {
        //  std::cout << "put back wrong(bond flip)!\n";
        //}
        int count = 0;
        for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
            mesh[*itr].ds = dsr_old[count];
            mesh[*itr].dAn2H = dAn2Hr_old[count];
            mesh[*itr].dskg = dskgr_old[count];
            mesh[*itr].dsk2 = dsk2r_old[count];
            mesh[*itr].dstau = dstaur_old[count];
            mesh[*itr].dE = dEr_old[count];
            count++;
        }
        return 0;
    }
}

int triangulation::edge_metropolis() {

    // observables
    std::vector<double> dsr_old; // ds related
    std::vector<double> Lers_old(Ne, 0);
    std::vector<double> Lers_new(Ne, 0);
    std::vector<std::vector<double>> dAn2Hr_old; // dAn2H related
    double IdAr_old, IdAr_new;
    double I2Hr_old, I2Hr_new;
    double I2H2r_old, I2H2r_new;
    std::vector<double> dskgr_old; // dskg related
    double Ikgr_old, Ikgr_new;
    std::vector<double> dsk2r_old; // dsk2 related
    double Ik2r_old, Ik2r_new;
    std::vector<double> dstaur_old; // dstau related
    double Itaur_old, Itaur_new;
    std::vector<double> dEr_old; // dE related
    double Er_old, Er_new;       // summation of related dr

    std::set<int, std::greater<int>>
        afc_beads; // beads that energy are affacted by the update
    std::set<int, std::greater<int>> afc_ebeads; // accordingly, but on the edge
    std::set<int, std::greater<int>>::iterator itr;
    double N_shrink, N_extend;
    std::vector<int> fedge_list;

    fedge_list.clear();
    for (int i = 0; i < edge_lists.size(); i++) {
        fedge_list.insert(fedge_list.end(), edge_lists[i].begin(),
                          edge_lists[i].end());
    } // flatten edge_lists and put in to one int vector
    // [shrink] (add a bond:2 bulk bond-1 edge bond)
    if (rand_uni(gen) < 0.5) {
        int ve_ind_i, e_ind_i, ind_i;
        /*
        renew_ve2remove_list();
        N_shrink = ve2remove_list.size();
        if (N_shrink == 0) {
            // std::cout << "no edge to remove\n";
            return 0;
        }
        ve_ind_i = int(N_shrink * rand_uni(gen));
        ind_i = ve2remove_list[ve_ind_i];
        */
        ve_ind_i = int(fedge_list.size() * rand_uni(gen));
        ind_i = fedge_list[ve_ind_i];
        e_ind_i = list_a_nei_b(edge_lists[mesh[ind_i].edge_num], ind_i);

        // remove ind_i
        /*
                        (j)
        (k)             =
            =         =
            = (i)

        to:
                = =   (j)
        (k)  =        -
            -       -
            - (i)
        */
        // notice the bond limit changes from l1 to l0 for i during this update,
        // vive versa for j,k
        int ind_j, ind_k;
        int j_e_nei_i, k_e_nei_i;
        int j_nei_i, k_nei_i;
        int j_nei_k, k_nei_j;
        int j_nei_i_next, j_nei_i_previous;
        int k_nei_i_next, k_nei_i_previous;

        std::vector<int> i_nei, j_nei, k_nei;
        // std::vector<int> i_e_nei, j_e_nei, k_e_nei;

        ind_j = mesh[ind_i].edge_nei[0];
        ind_k = mesh[ind_i].edge_nei[1];
        // check j k connection
        if (list_a_nei_b(mesh[ind_j].nei, ind_k) != -1) {
            return 0;
        }
        // check distance
        if (distance(ind_j, ind_k) >= l1 || distance(ind_i, ind_j) >= l0 ||
            distance(ind_i, ind_k) >= l0) {
            return 0;
        }
        // check # nei, can't be greater than 9
        if (mesh[ind_j].nei.size() >= 9 || mesh[ind_k].nei.size() >= 9) {
            return 0;
        }

        // for debug

        k_nei = mesh[ind_k].nei;
        // k_e_nei = mesh[ind_k].edge_nei;
        j_nei = mesh[ind_j].nei;
        // j_e_nei = mesh[ind_j].edge_nei;
        i_nei = mesh[ind_i].nei;
        // i_e_nei = mesh[ind_i].edge_nei;

        // find affected beads
        afc_beads.insert(ind_i);
        afc_beads.insert(ind_j);
        afc_beads.insert(ind_k);
        for (int j = 0; j < mesh[ind_j].edge_nei.size(); j++) {
            afc_beads.insert(mesh[ind_j].edge_nei[j]);
        }
        for (int j = 0; j < mesh[ind_k].edge_nei.size(); j++) {
            afc_beads.insert(mesh[ind_k].edge_nei[j]);
        }
        // store observable of affected beads
        dAn2Hr_old.clear();
        IdAr_old = 0;
        I2Hr_old = 0;
        I2H2r_old = 0;
        dsr_old.clear();
        for (int e = 0; e < Ne; e++) {
            Lers_old[e] = 0;
        }
        dskgr_old.clear();
        Ikgr_old = 0;
        dsk2r_old.clear();
        Ik2r_old = 0;
        dstaur_old.clear();
        Itaur_old = 0;
        dEr_old.clear();
        Er_old = 0;
        for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
            // kappa term
            dAn2Hr_old.push_back(mesh[*itr].dAn2H);
            IdAr_old += mesh[*itr].dAn2H[0];
            I2Hr_old += mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1];
            I2H2r_old +=
                mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1] * mesh[*itr].dAn2H[1];
            // lambda term
            dsr_old.push_back(mesh[*itr].ds);
            Lers_old[mesh[*itr].edge_num] += mesh[*itr].ds;
            // karg term
            dskgr_old.push_back(mesh[*itr].dskg);
            Ikgr_old += mesh[*itr].dskg;
            // B term
            dsk2r_old.push_back(mesh[*itr].dsk2);
            Ik2r_old += mesh[*itr].dsk2;
            // Bc term
            dstaur_old.push_back(mesh[*itr].dstau);
            Itaur_old += mesh[*itr].dstau;
            // energy
            dEr_old.push_back(mesh[*itr].dE);
            Er_old += mesh[*itr].dE;
        }

        // neighbors are also sorted during the energy calculation

        // [update]
        mesh[ind_i].edge_nei.clear();

        j_e_nei_i = list_a_nei_b(mesh[ind_j].edge_nei, ind_i);
        if (j_e_nei_i != 1) {
            std::cout << "j->i order?\n";
        }
        j_nei_i = list_a_nei_b(mesh[ind_j].nei, ind_i);
        if (list_a_nei_b(mesh[ind_k].nei, ind_j) != -1) {
            std::cout << "j-k connected!"
                      << "\n";
        }
        mesh[ind_j].edge_nei[j_e_nei_i] = ind_k;

        j_nei_i_next = (j_nei_i + 1) % mesh[ind_j].nei.size();
        j_nei_i_previous =
            (j_nei_i - 1 + mesh[ind_j].nei.size()) % mesh[ind_j].nei.size();
        if (mesh[ind_j].nei[j_nei_i_next] ==
            mesh[ind_j].edge_nei[1 - j_e_nei_i]) {
            // j_nei_i_next is the other edge, can push i forward
            mesh[ind_j].nei.insert(mesh[ind_j].nei.begin() + j_nei_i_next,
                                   ind_k);
            j_nei_k = j_nei_i_next;
        } else {
            mesh[ind_j].nei.insert(mesh[ind_j].nei.begin() + j_nei_i, ind_k);
            j_nei_k = j_nei_i;
        } // okie!

        k_e_nei_i = list_a_nei_b(mesh[ind_k].edge_nei, ind_i);
        if (k_e_nei_i != 0) {
            std::cout << "i->k order?\n";
        }
        k_nei_i = list_a_nei_b(mesh[ind_k].nei, ind_i);
        if (list_a_nei_b(mesh[ind_k].nei, ind_j) != -1) {
            std::cout << "k-j connected!"
                      << "\n";
        }
        mesh[ind_k].edge_nei[k_e_nei_i] = ind_j;
        k_nei_i_next = (k_nei_i + 1) % mesh[ind_k].nei.size();
        k_nei_i_previous =
            (k_nei_i - 1 + mesh[ind_k].nei.size()) % mesh[ind_k].nei.size();
        if (mesh[ind_k].nei[k_nei_i_next] ==
            mesh[ind_k].edge_nei[1 - k_e_nei_i]) {
            // k_nei_i_next in side, can push i forward
            mesh[ind_k].nei.insert(mesh[ind_k].nei.begin() + k_nei_i_next,
                                   ind_j);
            k_nei_j = k_nei_i_next;
        } else {
            // j_nei_i_next is edge, can push i backward
            mesh[ind_k].nei.insert(mesh[ind_k].nei.begin() + k_nei_i, ind_j);
            k_nei_j = k_nei_i;
        }

        if (mesh[ind_i].edge_nei.size() != 0) {
            std::cout << "ind_i edge not clear"
                      << "\n";
        }
        // renew_v2add_list();
        // N_extend = v2add_list.size();

        // after-update observables
        IdAr_new = 0;
        I2Hr_new = 0;
        I2H2r_new = 0;
        for (int e = 0; e < Ne; e++) {
            Lers_new[e] = 0;
        }
        Ikgr_new = 0;
        Ik2r_new = 0;
        Itaur_new = 0;
        Er_new = 0;
        for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
            // kappa term
            mesh[*itr].dAn2H = dAn2H_m(*itr);
            IdAr_new += mesh[*itr].dAn2H[0];
            I2Hr_new += mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1];
            I2H2r_new +=
                mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1] * mesh[*itr].dAn2H[1];
            // lambda term (needed by other edge terms)
            mesh[*itr].ds = ds_m(*itr);
            Lers_new[mesh[*itr].edge_num] += mesh[*itr].ds;
            // karg term
            mesh[*itr].dskg = dskg_m(*itr);
            Ikgr_new += mesh[*itr].dskg;
            // B term
            mesh[*itr].dsk2 = dsk2_m(*itr);
            Ik2r_new += mesh[*itr].dsk2;
            // Bc term
            mesh[*itr].dstau = dstau_m(*itr);
            Itaur_new += mesh[*itr].dstau;
            // energy
            mesh[*itr].dE = dE_m(*itr);
            Er_new += mesh[*itr].dE;
        }

        // [Metropolis]
        if (rand_uni(gen) <= 1.0 * fedge_list.size() / (fedge_list.size() - 1) *
                                 std::exp(-beta * (Er_new - Er_old))) {
            // [accepted]
            // std::cout << "shrink: (Er_new - Er_old)=" << (Er_new - Er_old)<<
            // "\n";
            // std::cout << "N_edge/N_shrink/N_extend=" << edge_lists[0].size()
            // << "/" << N_shrink << "/" << N_extend << "\n";
            edge_lists[mesh[ind_i].edge_num].erase(
                edge_lists[mesh[ind_i].edge_num].begin() + e_ind_i);
            mesh[ind_i].edge_num = -1;

            E += Er_new - Er_old;
            for (int e = 0; e < Ne; e++) {
                Les[e] += Lers_new[e] - Lers_old[e];
            }
            IdA += IdAr_new - IdAr_old;
            I2H += I2Hr_new - I2Hr_old;
            I2H2 += I2H2r_new - I2H2r_old;
            Ikg += Ikgr_new - Ikgr_old;
            Ik2 += Ik2r_new - Ik2r_old;
            Itau += Itaur_new - Itaur_old;

            // update bond_list, add i-j i-k as bulk bond
            std::pair<int, int> bond1, bond2;
            bond1.first = ind_i;
            bond1.second = ind_j;
            bond2.first = ind_j;
            bond2.second = ind_i;
            bond_list.push_back(bond1);
            bond_list.push_back(bond2);
            bond1.first = ind_i;
            bond1.second = ind_k;
            bond2.first = ind_k;
            bond2.second = ind_i;
            bond_list.push_back(bond1);
            bond_list.push_back(bond2);

            /*
            update_v2flip_local(ind_i);
            update_v2flip_group(mesh[ind_i].nei);
            update_v2flip_local(ind_j);
            update_v2flip_group(mesh[ind_j].nei);
            update_v2flip_local(ind_k);
            update_v2flip_group(mesh[ind_k].nei);
            */
            // check duplication
            // check_duplication(ind_i);
            // check_duplication(ind_j);
            // check_duplication(ind_k);

            // std::cout << "edge_list.size()=" << edge_list.size() << "\n";
            // if (std::abs(Er_new - Er_old) > 100) {
            //   std::cout << "energy change (edge remove): "
            //            << (Er_new - Er_old) << "\n";
            //}
            return 1;
        } else {
            // [rejected]
            // put everything back (so sad)
            // notice nei of j and k are sorted! not anymore

            mesh[ind_k].nei.erase(mesh[ind_k].nei.begin() + k_nei_j);
            mesh[ind_k].edge_nei[k_e_nei_i] = ind_i;

            mesh[ind_j].nei.erase(mesh[ind_j].nei.begin() + j_nei_k);
            mesh[ind_j].edge_nei[j_e_nei_i] = ind_i;

            mesh[ind_i].edge_nei.push_back(ind_j);
            mesh[ind_i].edge_nei.push_back(ind_k);
            if (mesh[ind_i].nei != i_nei || mesh[ind_j].nei != j_nei ||
                mesh[ind_k].nei != k_nei) {
                std::cout << "put back wrong()!\n";
            }
            int count = 0;
            for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
                mesh[*itr].ds = dsr_old[count];
                mesh[*itr].dAn2H = dAn2Hr_old[count];
                mesh[*itr].dskg = dskgr_old[count];
                mesh[*itr].dsk2 = dsk2r_old[count];
                // mesh[*itr].nc = ncr_old[count];
                mesh[*itr].dstau = dstaur_old[count];
                mesh[*itr].dE = dEr_old[count];
                count++;
            }
            return 0;
        }
    } // end remove ind_i
    // [extend] (remove a bond: 1 edge bond - 2 bulk bond)
    else {
        // add i (from the bulk) to edge

        std::vector<int> i_nei, j_nei, k_nei, l_nei;
        // std::vector<int> i_e_nei, j_e_nei, k_e_nei;
        int ind_i, e_ind_j, ind_j, ind_k;
        int k_nei_j, k_nei_i, j_nei_k;
        int k_e_nei_j, j_e_nei_k;
        int v2a_ind_i;
        int i_nei_j, i_nei_k;
        // renew_v2add_list();
        // N_extend = v2add_list.size();
        // if (N_extend == 0) {
        // return 0;
        //}
        // v2a_ind_i = int(v2add_list.size() * rand_uni(gen));
        e_ind_j = int(fedge_list.size() * rand_uni(gen));
        // ind_i = v2add_list[v2a_ind_i];
        ind_j = fedge_list[e_ind_j];
        ind_k = mesh[ind_j].edge_nei[1];
        if (mesh[ind_k].edge_nei[0] != ind_j) {
            std::cout << "j-k order wrong here!\n";
        }
        // check # nei, can't be less than 3
        if (mesh[ind_j].nei.size() <= 3 || mesh[ind_k].nei.size() <= 3) {
            return 0;
        }
        // let's find i
        j_nei = mesh[ind_j].nei;
        k_nei = mesh[ind_k].nei;
        ind_i = -1;
        int j_next, j_previous;
        int ind_l, l_nei_j;
        for (int l = 0; l < mesh[ind_j].nei.size(); l++) {
            ind_l = mesh[ind_j].nei[l];
            l_nei = mesh[ind_l].nei;
            if ((list_a_nei_b(mesh[ind_k].nei, ind_l) != -1) &&
                (mesh[ind_l].edge_nei.size() == 0)) {
                // if (ind_i != -1 && ind_i != mesh[ind_j].nei[l]) {
                // State_write("State_error.txt");
                // std::cout << "must be kidding ?\n";
                //}

                l_nei_j = list_a_nei_b(mesh[ind_l].nei, ind_j);
                j_next = (l_nei_j + 1) % mesh[ind_l].nei.size();
                j_previous = (l_nei_j - 1 + mesh[ind_l].nei.size()) %
                             mesh[ind_l].nei.size();
                if (mesh[ind_l].nei[j_next] == ind_k ||
                    mesh[ind_l].nei[j_previous] == ind_k) {
                    ind_i = ind_l;
                    break; // k is the next of j in i_nei
                }
            }
        }

        if (ind_i == -1) {
            // std::cout << "couldn't find ind_i!!?\n";
            return 0;
        }
        i_nei = mesh[ind_i].nei;

        /*
        int j_next;
        for (int j = 0; j < mesh[ind_i].nei.size(); j++) {
            j_next = (j + 1) % mesh[ind_i].nei.size();
            if (mesh[mesh[ind_i].nei[j]].edge_nei.size() != 0 &&
                mesh[mesh[ind_i].nei[j_next]].edge_nei.size() != 0 &&
                (mesh[mesh[ind_i].nei[j]].edge_nei[0] ==
                     mesh[ind_i].nei[j_next] ||
                 mesh[mesh[ind_i].nei[j]].edge_nei[1] ==
                     mesh[ind_i].nei[j_next])) {
                i_nei_j = j;
                i_nei_k = j_next;
                break;
            }
        }


        ind_j = mesh[ind_i].nei[i_nei_j];
        ind_k = mesh[ind_i].nei[i_nei_k];
        */
        // order of j k is unclear, unclear now

        // add i to the edge

        afc_beads.insert(ind_i);
        afc_beads.insert(ind_j);
        afc_beads.insert(ind_k);
        for (int j = 0; j < mesh[ind_j].edge_nei.size(); j++) {
            afc_beads.insert(mesh[ind_j].edge_nei[j]);
        }
        for (int j = 0; j < mesh[ind_k].edge_nei.size(); j++) {
            afc_beads.insert(mesh[ind_k].edge_nei[j]);
        }
        // store observable of affected beads
        dAn2Hr_old.clear();
        IdAr_old = 0;
        I2Hr_old = 0;
        I2H2r_old = 0;
        dsr_old.clear();
        for (int e = 0; e < Ne; e++) {
            Lers_old[e] = 0;
        }
        dskgr_old.clear();
        Ikgr_old = 0;
        dsk2r_old.clear();
        Ik2r_old = 0;
        dstaur_old.clear();
        Itaur_old = 0;
        dEr_old.clear();
        Er_old = 0;
        for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
            // kappa term
            dAn2Hr_old.push_back(mesh[*itr].dAn2H);
            IdAr_old += mesh[*itr].dAn2H[0];
            I2Hr_old += mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1];
            I2H2r_old +=
                mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1] * mesh[*itr].dAn2H[1];
            // lambda term
            dsr_old.push_back(mesh[*itr].ds);
            Lers_old[mesh[*itr].edge_num] += mesh[*itr].ds;
            // karg term
            dskgr_old.push_back(mesh[*itr].dskg);
            Ikgr_old += mesh[*itr].dskg;
            // B term
            dsk2r_old.push_back(mesh[*itr].dsk2);
            Ik2r_old += mesh[*itr].dsk2;
            // Bc term
            dstaur_old.push_back(mesh[*itr].dstau);
            Itaur_old += mesh[*itr].dstau;
            // energy
            dEr_old.push_back(mesh[*itr].dE);
            Er_old += mesh[*itr].dE;
        }
        // sort_nei(ind_i);

        // add energy abd sorting k

        // just for debugging

        k_nei = mesh[ind_k].nei;
        // k_e_nei = mesh[ind_k].edge_nei;
        j_nei = mesh[ind_j].nei;
        // j_e_nei = mesh[ind_j].edge_nei;
        i_nei = mesh[ind_i].nei;
        // i_e_nei = mesh[ind_i].edge_nei;

        // [update]
        // will remove k_j link and use j-i k_i as new edge
        // need to follow the order of J-K link
        if (mesh[ind_j].edge_nei[0] == ind_k) {
            // j<-k
            std::cout << "oh?\n";
            mesh[ind_i].edge_nei.push_back(ind_k);
            mesh[ind_i].edge_nei.push_back(ind_j);
        } else if (mesh[ind_j].edge_nei[1] == ind_k) {
            // j->k
            mesh[ind_i].edge_nei.push_back(ind_j);
            mesh[ind_i].edge_nei.push_back(ind_k);
        } else {
            std::cout << "j-k not connected!\n";
        }

        k_nei_j = list_a_nei_b(mesh[ind_k].nei, ind_j);
        k_e_nei_j = list_a_nei_b(mesh[ind_k].edge_nei, ind_j);
        if (k_e_nei_j != 0) {
            std::cout << "k_e_nei_j!=0\n";
        }

        if (k_nei_j == -1 || k_e_nei_j == -1) {
            std::cout << "k_j not connected here (1)"
                      << "\n";
        }
        // ind_k break up with j
        mesh[ind_k].edge_nei[k_e_nei_j] = ind_i;
        mesh[ind_k].nei.erase(mesh[ind_k].nei.begin() + k_nei_j);

        j_e_nei_k = list_a_nei_b(mesh[ind_j].edge_nei, ind_k);
        j_nei_k = list_a_nei_b(mesh[ind_j].nei, ind_k);
        if (j_nei_k == -1 || j_e_nei_k == -1) {
            std::cout << "j_k not connected here (1)"
                      << "\n";
        }

        // ind_j break up with k
        // recall j_e_nei_k = 1s
        mesh[ind_j].edge_nei[j_e_nei_k] = ind_i;
        mesh[ind_j].nei.erase(mesh[ind_j].nei.begin() + j_nei_k);

        mesh[ind_i].edge_num = mesh[ind_j].edge_num;
        // renew_ve2remove_list();
        // N_shrink = ve2remove_list.size();

        // after-update observables
        IdAr_new = 0;
        I2Hr_new = 0;
        I2H2r_new = 0;
        for (int e = 0; e < Ne; e++) {
            Lers_new[e] = 0;
        }
        Ikgr_new = 0;
        Ik2r_new = 0;
        Itaur_new = 0;
        Er_new = 0;
        for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
            // kappa term
            mesh[*itr].dAn2H = dAn2H_m(*itr);
            IdAr_new += mesh[*itr].dAn2H[0];
            I2Hr_new += mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1];
            I2H2r_new +=
                mesh[*itr].dAn2H[0] * mesh[*itr].dAn2H[1] * mesh[*itr].dAn2H[1];
            // lambda term (needed by other edge terms)
            mesh[*itr].ds = ds_m(*itr);
            Lers_new[mesh[*itr].edge_num] += mesh[*itr].ds;
            // karg term
            mesh[*itr].dskg = dskg_m(*itr);
            Ikgr_new += mesh[*itr].dskg;
            // B term
            mesh[*itr].dsk2 = dsk2_m(*itr);
            Ik2r_new += mesh[*itr].dsk2;
            // Bc term
            mesh[*itr].dstau = dstau_m(*itr);
            Itaur_new += mesh[*itr].dstau;
            // energy
            mesh[*itr].dE = dE_m(*itr);
            Er_new += mesh[*itr].dE;
        }

        // [ Metropolis]

        // if (Er_new - Er_old < 5) {
        //}
        if (rand_uni(gen) <= 1.0 * fedge_list.size() / (fedge_list.size() + 1) *
                                 std::exp(-beta * (Er_new - Er_old))) {
            // [accepted]
            // std::cout << "extend: (Er_new - Er_old)=" << (Er_new - Er_old)<<
            // "\n";
            // "\n"; std::cout << "N_extend / N_shrink=" << N_extend / N_shrink
            // <<
            // "\n";
            // mesh[ind_i].edge_num = mesh[ind_j].edge_num;
            if (mesh[ind_i].edge_num == -1) {
                std::cout << "somewhat edge_num wrong\n";
            }
            edge_lists[mesh[ind_i].edge_num].push_back(ind_i);
            // I2H_new = 0;
            // i j k are sorted again!
            E += Er_new - Er_old;
            for (int e = 0; e < Ne; e++) {
                Les[e] += Lers_new[e] - Lers_old[e];
            }
            IdA += IdAr_new - IdAr_old;
            I2H += I2Hr_new - I2Hr_old;
            I2H2 += I2H2r_new - I2H2r_old;
            Ikg += Ikgr_new - Ikgr_old;
            Ik2 += Ik2r_new - Ik2r_old;
            Itau += Itaur_new - Itaur_old;

            // these are edge bond now
            delete_bond_list(ind_i, ind_j);
            delete_bond_list(ind_i, ind_k);
            // update nei2flip
            /*
            update_v2flip_local(ind_i);
            update_v2flip_group(mesh[ind_i].nei);
            update_v2flip_local(ind_j);
            update_v2flip_group(mesh[ind_j].nei);
            update_v2flip_local(ind_k);
            update_v2flip_group(mesh[ind_k].nei);
            */
            // check duplication
            /*
            check_duplication(ind_i);
            check_duplication(ind_j);
            check_duplication(ind_k);
            */
            return 1;
        } else {
            // [rejected]
            // rejected put everything back (so sad)
            mesh[ind_i].edge_num = -1;
            mesh[ind_j].nei.insert(mesh[ind_j].nei.begin() + j_nei_k, ind_k);
            mesh[ind_j].edge_nei[j_e_nei_k] = ind_k;

            mesh[ind_k].nei.insert(mesh[ind_k].nei.begin() + k_nei_j, ind_j);
            mesh[ind_k].edge_nei[k_e_nei_j] = ind_j;

            mesh[ind_i].edge_nei.clear();
            if (mesh[ind_i].nei != i_nei || mesh[ind_j].nei != j_nei ||
                mesh[ind_k].nei != k_nei) {
                std::cout << "put back wrong!\n";
            }
            int count = 0;
            for (itr = afc_beads.begin(); itr != afc_beads.end(); itr++) {
                mesh[*itr].ds = dsr_old[count];
                mesh[*itr].dAn2H = dAn2Hr_old[count];
                mesh[*itr].dskg = dskgr_old[count];
                mesh[*itr].dsk2 = dsk2r_old[count];
                // mesh[*itr].nc = ncr_old[count];
                mesh[*itr].dstau = dstaur_old[count];
                mesh[*itr].dE = dEr_old[count];
                count++;
            }

            return 0;
        }
    } // end add ind_i
}
