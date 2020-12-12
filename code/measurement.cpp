#include "triangulation.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#define PI 3.14159265358979323846

double innerproduct(std::vector<double> a, std::vector<double> b) {
    double result = 0;
    for (int k = 0; k < a.size(); k++) {
        result += a[k] * b[k];
    }
    return result;
}
std::vector<double> crossproduct(std::vector<double> a, std::vector<double> b) {
    std::vector<double> cp{0, 0, 0};
    cp[0] = a[1] * b[2] - a[2] * b[1];
    cp[1] = a[2] * b[0] - a[0] * b[2];
    cp[2] = a[0] * b[1] - a[1] * b[0];
    return cp;
}
double triangulation::distance(int ind_1, int ind_2) {
    double result = 0;
    for (int k = 0; k < 3; k++) {
        result += std::pow(mesh[ind_2].pos[k] - mesh[ind_1].pos[k], 2);
    }
    result = std::sqrt(result);
    return result;
}

int triangulation::sort_nei(int index) {
    std::vector<int> index_nei = mesh[index].nei;
    std::vector<int> index_enei = mesh[index].edge_nei;
    if (mesh[index].edge_nei.size() == 0) {
        return 0;
    }
    if (mesh[index].nei.size() == 3) {
        return 0;
    }
    while (mesh[index].nei[0] != mesh[index].edge_nei[0]) {
        mesh[index].nei.push_back(mesh[index].nei[0]);
        mesh[index].nei.erase(mesh[index].nei.begin());
    }
    if (mesh[index].nei[1] == mesh[index].edge_nei[1]) {
        mesh[index].nei.push_back(mesh[index].nei[0]);
        mesh[index].nei.erase(mesh[index].nei.begin());
    } else if (mesh[index].nei.back() != mesh[index].edge_nei[1]) {
        std::cout << "wrong, can't sort neighbors!\n";
        return 0;
    }
    return 1;
}

double triangulation::ds_m(int index) {
    if (mesh[index].edge_nei.size() == 0) {
        return 0;
    }
    double Le_l = 0;
    int nei_ind;
    double distance;
    for (int i = 0; i < mesh[index].edge_nei.size(); i++) {
        nei_ind = mesh[index].edge_nei[i];
        distance = 0;
        for (int j = 0; j < 3; j++) {
            // Rij[j] = mesh[index].pos[j] - mesh[nei_ind].pos[j];
            distance += std::pow(mesh[index].pos[j] - mesh[nei_ind].pos[j], 2);
        }
        distance = std::sqrt(distance);
        Le_l += 0.5 * distance;
    }
    if (Le_l > 3) {
        std::cout << "Le_l=" << Le_l << "\n";
    }
    return Le_l;
}
double triangulation::dA2H2_m(int index) {
    if (mesh[index].edge_nei.size() != 0) {
        return 0;
    }
    double l_ij, l_ij0, l_ij2, l_jj0, l_jj2;
    double dot_ij0j, dot_ij2j, cot0, cot2, theta0, theta2;
    double sigma_ij, sigma_i;
    int ind_j0, ind_j, ind_j2;

    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ij0{0, 0, 0};
    std::vector<double> r_ij2{0, 0, 0};
    std::vector<double> r_jj0{0, 0, 0};
    std::vector<double> r_jj2{0, 0, 0};
    /* vertex configuration:
               - (j2)
             -
        (i)- - - (j)
             -
               - (j0)
    */
    sigma_i = 0;
    std::vector<double> dH_i{0, 0, 0}; // 2H_i, depends on r_ij
    for (int j = 0; j < mesh[index].nei.size(); j++) {
        // get nei index
        ind_j = mesh[index].nei[j];
        if (j == 0) {
            ind_j0 = mesh[index].nei[mesh[index].nei.size() - 1];
            ind_j2 = mesh[index].nei[j + 1];
        } else if (j == mesh[index].nei.size() - 1) {
            ind_j0 = mesh[index].nei[j - 1];
            ind_j2 = mesh[index].nei[0];
        } else {
            ind_j0 = mesh[index].nei[j - 1];
            ind_j2 = mesh[index].nei[j + 1];
        }
        // end set nei_ind
        // triangulation site to site vectors
        for (int k = 0; k < 3; k++) {
            r_ij[k] = mesh[ind_j].pos[k] - mesh[index].pos[k];
            r_ij0[k] = mesh[ind_j0].pos[k] - mesh[index].pos[k];
            r_ij2[k] = mesh[ind_j2].pos[k] - mesh[index].pos[k];
            r_jj0[k] = mesh[ind_j0].pos[k] - mesh[ind_j].pos[k];
            r_jj2[k] = mesh[ind_j2].pos[k] - mesh[ind_j].pos[k];
        }
        // site-site distance
        l_ij = std::sqrt(innerproduct(r_ij, r_ij));
        l_ij0 = std::sqrt(innerproduct(r_ij0, r_ij0));
        l_ij2 = std::sqrt(innerproduct(r_ij2, r_ij2));
        l_jj0 = std::sqrt(innerproduct(r_jj0, r_jj0));
        l_jj2 = std::sqrt(innerproduct(r_jj2, r_jj2));
        // inner product for cot calculation
        dot_ij0j = innerproduct(r_ij0, r_jj0);
        theta0 = std::acos(dot_ij0j / (l_ij0 * l_jj0));
        dot_ij2j = innerproduct(r_ij2, r_jj2);
        theta2 = std::acos(dot_ij2j / (l_ij2 * l_jj2));
        // get cot
        /*
        cot0 = dot_ij0j / std::sqrt(std::pow(l_ij0 * l_jj0, 2) -
                                    std::pow(dot_ij0j, 2));
        cot2 = dot_ij2j / std::sqrt(std::pow(l_ij2 * l_jj2, 2) -
                                    std::pow(dot_ij2j, 2));
        */
        cot0 = std::cos(theta0) / std::sin(theta0);
        cot2 = std::cos(theta2) / std::sin(theta2);
        sigma_ij = 0.5 * l_ij * (cot0 + cot2);
        sigma_i += 0.25 * sigma_ij * l_ij;
        for (int k = 0; k < 3; k++) {
            dH_i[k] += sigma_ij / l_ij * r_ij[k];
        }
    } // end for nei
    for (int k = 0; k < 3; k++) {
        dH_i[k] = dH_i[k] / sigma_i - C0;
    }
    return sigma_i * innerproduct(dH_i, dH_i);
}

std::vector<double> triangulation::dAn2H_m(int index) {
    std::vector<double> dAn2H{0.0, 0.0};
    if (mesh[index].edge_nei.size() != 0) {
        return dAn2H;
    }
    double l_ij, l_ij0, l_ij2, l_jj0, l_jj2;
    double dot_ij0j, dot_ij2j, cot0, cot2, theta0, theta2;
    double sigma_ij, sigma_i;
    int ind_j0, ind_j, ind_j2;

    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ij0{0, 0, 0};
    std::vector<double> r_ij2{0, 0, 0};
    std::vector<double> r_jj0{0, 0, 0};
    std::vector<double> r_jj2{0, 0, 0};
    /* vertex configuration:
               - (j2)
             -
        (i)- - - (j)
             -
               - (j0)
    */
    sigma_i = 0;
    std::vector<double> dH_i{0, 0, 0}; // 2H_i, depends on r_ij
    for (int j = 0; j < mesh[index].nei.size(); j++) {
        // get nei index
        ind_j = mesh[index].nei[j];
        if (j == 0) {
            ind_j0 = mesh[index].nei[mesh[index].nei.size() - 1];
            ind_j2 = mesh[index].nei[j + 1];
        } else if (j == mesh[index].nei.size() - 1) {
            ind_j0 = mesh[index].nei[j - 1];
            ind_j2 = mesh[index].nei[0];
        } else {
            ind_j0 = mesh[index].nei[j - 1];
            ind_j2 = mesh[index].nei[j + 1];
        }
        // end set nei_ind
        // triangulation site to site vectors
        for (int k = 0; k < 3; k++) {
            r_ij[k] = mesh[ind_j].pos[k] - mesh[index].pos[k];
            r_ij0[k] = mesh[ind_j0].pos[k] - mesh[index].pos[k];
            r_ij2[k] = mesh[ind_j2].pos[k] - mesh[index].pos[k];
            r_jj0[k] = mesh[ind_j0].pos[k] - mesh[ind_j].pos[k];
            r_jj2[k] = mesh[ind_j2].pos[k] - mesh[ind_j].pos[k];
        }
        // site-site distance
        l_ij = std::sqrt(innerproduct(r_ij, r_ij));
        l_ij0 = std::sqrt(innerproduct(r_ij0, r_ij0));
        l_ij2 = std::sqrt(innerproduct(r_ij2, r_ij2));
        l_jj0 = std::sqrt(innerproduct(r_jj0, r_jj0));
        l_jj2 = std::sqrt(innerproduct(r_jj2, r_jj2));
        // inner product for cot calculation
        dot_ij0j = innerproduct(r_ij0, r_jj0);
        theta0 = std::acos(dot_ij0j / (l_ij0 * l_jj0));
        dot_ij2j = innerproduct(r_ij2, r_jj2);
        theta2 = std::acos(dot_ij2j / (l_ij2 * l_jj2));
        // get cot
        /*
        cot0 = dot_ij0j / std::sqrt(std::pow(l_ij0 * l_jj0, 2) -
                                    std::pow(dot_ij0j, 2));
        cot2 = dot_ij2j / std::sqrt(std::pow(l_ij2 * l_jj2, 2) -
                                    std::pow(dot_ij2j, 2));
        */
        cot0 = std::cos(theta0) / std::sin(theta0);
        cot2 = std::cos(theta2) / std::sin(theta2);
        sigma_ij = 0.5 * l_ij * (cot0 + cot2);
        sigma_i += 0.25 * sigma_ij * l_ij;
        for (int k = 0; k < 3; k++) {
            dH_i[k] += sigma_ij / l_ij * r_ij[k];
        }
    } // end for nei
    for (int k = 0; k < 3; k++) {
        dH_i[k] = dH_i[k] / sigma_i;
    }
    dAn2H[0] = sigma_i;
    dAn2H[1] = std::sqrt(innerproduct(dH_i, dH_i));
    return dAn2H;
}

double triangulation::dskg_m(int index) {
    if (mesh[index].edge_nei.size() == 0) {
        return 0;
    }
    double kg;
    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    double l_ij, l_ik;
    int ind_j, ind_k;
    double cos_jk;
    /* vertex configuration:
               - (k)
             -
        (i)- - - (j)
    */
    std::vector<int> nei_original;
    // sort neighbors
    nei_original = mesh[index].nei;
    sort_nei(index);
    // get thetai
    // kg = PI - sum(theta_i)
    kg = PI;
    for (int j = 0; j < mesh[index].nei.size() - 1; j++) {
        cos_jk = 0;
        l_ij = 0;
        l_ik = 0;
        ind_j = mesh[index].nei[j];
        ind_k = mesh[index].nei[j + 1];
        for (int k = 0; k < 3; k++) {
            r_ij[k] = mesh[ind_j].pos[k] - mesh[index].pos[k];
            l_ij += r_ij[k] * r_ij[k];
            r_ik[k] = mesh[ind_k].pos[k] - mesh[index].pos[k];
            l_ik += r_ik[k] * r_ik[k];
            cos_jk += r_ij[k] * r_ik[k];
        }
        cos_jk = cos_jk / (std::sqrt(l_ij * l_ik));
        kg -= std::acos(cos_jk);
    }
    mesh[index].nei = nei_original; // put nei sort back;
    return kg;
}

double triangulation::dsk2_m(int index) {
    if (mesh[index].edge_nei.size() == 0) {
        return 0;
    }
    /*
    (j)===(i)===(k)
    */
    double cos_jk, theta_jk;
    double l2_ji, l2_ik;
    int ind_j, ind_k;
    std::vector<double> r_ji{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    cos_jk = 0;
    l2_ji = 0;
    l2_ik = 0;
    ind_j = mesh[index].edge_nei[0];
    ind_k = mesh[index].edge_nei[1];
    for (int k = 0; k < 3; k++) {
        r_ji[k] = mesh[index].pos[k] - mesh[ind_j].pos[k];
        l2_ji += r_ji[k] * r_ji[k];
        r_ik[k] = mesh[ind_k].pos[k] - mesh[index].pos[k];
        l2_ik += r_ik[k] * r_ik[k];
        cos_jk += r_ji[k] * r_ik[k];
    }
    cos_jk = cos_jk / (std::sqrt(l2_ji * l2_ik));
    cos_jk = std::min(cos_jk, 1.0);  // make sure cos is not over bound
    cos_jk = std::max(cos_jk, -1.0); // make sure cos is not over bound
    theta_jk = std::acos(cos_jk);
    return theta_jk * theta_jk / ds_m(index);
}

std::vector<double> triangulation::nc_m(int index) {
    std::vector<double> nc{0, 0, 0};
    if (mesh[index].edge_nei.size() == 0) {
        return nc; // not valid for beads in the bulk
    }
    double nc_abs;
    int ind_j, ind_k;
    double l_ij, l_ik;
    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    std::vector<double> n_jk{0, 0, 0};
    double n_jk_abs;
    double theta_jk;
    std::vector<int> nei_original = mesh[index].nei;
    sort_nei(index);
    for (int j = 0; j < mesh[index].nei.size() - 1; j++) {
        ind_j = mesh[index].nei[j];
        ind_k = mesh[index].nei[j + 1];
        for (int k = 0; k < 3; k++) {
            r_ij[k] = mesh[ind_j].pos[k] - mesh[index].pos[k];
            r_ik[k] = mesh[ind_k].pos[k] - mesh[index].pos[k];
        }
        l_ij = std::sqrt(innerproduct(r_ij, r_ij));
        l_ik = std::sqrt(innerproduct(r_ik, r_ik));
        n_jk = crossproduct(r_ij, r_ik);
        n_jk_abs = std::sqrt(innerproduct(n_jk, n_jk));
        theta_jk = std::asin(n_jk_abs / (l_ij * l_ik));
        for (int k = 0; k < 3; k++) {
            n_jk[k] = n_jk[k] / n_jk_abs;
            nc[k] += theta_jk * n_jk[k];
        }
    }
    nc_abs = std::sqrt(innerproduct(nc, nc));
    for (int k = 0; k < 3; k++) {
        nc[k] = nc[k] / nc_abs;
    }
    mesh[index].nei = nei_original; // put nei sort back;
    return nc;
}

double triangulation::dstau_m(int index) {
    if (mesh[index].edge_nei.size() == 0) {
        return 0;
    }
    std::vector<double> t{0, 0, 0};    // unit tangent vector
    std::vector<double> nc_0{0, 0, 0}; // unit surface normal on edge_nei[0]
    std::vector<double> nc_1{0, 0, 0}; // unit surface normal on edge_nei[1]

    std::vector<double> nc_i{0, 0, 0}; // nc(i)
    std::vector<double> dnc{0, 0, 0};  // 1/2 (nc_1-nc_0)

    int ind_0 = mesh[index].edge_nei[0];
    int ind_1 = mesh[index].edge_nei[1];

    // calculate tangent vector
    for (int k = 0; k < 3; k++) {
        t[k] = mesh[ind_1].pos[k] - mesh[ind_0].pos[k];
    }
    double t_abs = std::sqrt(innerproduct(t, t));
    for (int k = 0; k < 3; k++) {
        t[k] = t[k] / t_abs;
    }
    // calculate surface normal
    nc_0 = nc_m(ind_0);
    nc_1 = nc_m(ind_1);
    nc_i = nc_m(index);
    for (int k = 0; k < 3; k++) {
        // nc_i[k] = (nc_1[k] + nc_0[k]) / 2;
        dnc[k] = (nc_1[k] - nc_0[k]) / 2;
    }
    return innerproduct(t, crossproduct(nc_i, dnc));
}
double triangulation::dE_m(int index) {
    double dE_i = 0;
    double dstau_i;
    double ds_i;
    if (mesh[index].edge_nei.size() == 0) {
        dE_i += 0.5 * kar * mesh[index].dAn2H[0] * (mesh[index].dAn2H[1] - C0) *
                (mesh[index].dAn2H[1] - C0);
    } else {
        dE_i += -karg * mesh[index].dskg;
        dE_i += lam * mesh[index].ds;
        dE_i += 0.5 * B * mesh[index].dsk2;
        dstau_i = mesh[index].dstau;
        ds_i = mesh[index].ds;
        dE_i +=
            0.5 * Bc * ds_i * (dstau_i / ds_i - tau0) * (dstau_i / ds_i - tau0);
    }
    return dE_i;
}
/*
double triangulation::dRg2_m(int index) {
    double dRg2_i = 0;
    for (int j = 0; j < N; j++) {
        dRg2_i += distance(j, index);
    }
    dRg2_i /= 2 * N * N;
}
 */

std::vector<double> triangulation::Cncnc_m(int ednum) {
    // sort edge
    std::vector<double> Cncnc;
    std::vector<int> sort_edge_list;
    std::vector<double> ave_nc(3, 0);
    int next_i;

    sort_edge_list.push_back(edge_lists[ednum][0]);
    for (int i = 0; i < edge_lists[ednum].size() - 1; i++) {
        next_i = mesh[sort_edge_list.back()].edge_nei[0];
        sort_edge_list.push_back(next_i);
    }
    if (mesh[sort_edge_list.back()].edge_nei[0] != sort_edge_list[0]) {
        std::cout << "not in a circle?\n";
    }
    // measure nc of all edge beads
    for (int i = 0; i < sort_edge_list.size(); i++) {
        mesh[sort_edge_list[i]].nc = nc_m(sort_edge_list[i]);
        for (int k = 0; k < 3; k++) {
            ave_nc[k] += mesh[sort_edge_list[i]].nc[k];
        }
    }
    for (int k = 0; k < 3; k++) {
        ave_nc[k] /= sort_edge_list.size();
    }

    if (sort_edge_list.size() != edge_lists[ednum].size()) {
        std::cout << "not same size after sorted edge\n";
    }
    // calculate correlation
    for (int s = 0; s < (sort_edge_list.size() + 1); s++) {
        Cncnc.push_back(0.0);
        for (int i = 0; i < sort_edge_list.size(); i++) {
            int j = (i + s) % sort_edge_list.size();
            Cncnc[s] += innerproduct(mesh[sort_edge_list[i]].nc,
                                     mesh[sort_edge_list[j]].nc);
        }
        Cncnc[s] /= sort_edge_list.size();
    }
    Cncnc.insert(Cncnc.begin(), std::sqrt(innerproduct(ave_nc, ave_nc)));
    return Cncnc;
}
std::vector<double> triangulation::Iij_m() {
    // I_ij  = Jij - N*(\delta_ij r_c*r_c - xc_i*x_c*j)
    std::vector<double> Jij(9, 0);
    std::vector<double> Iij(9, 0);
    double xc = 0, yc = 0, zc = 0;
    for (int n = 0; n < mesh.size(); n++) {
        // Jxx
        Jij[0] +=
            mesh[n].pos[1] * mesh[n].pos[1] + mesh[n].pos[2] * mesh[n].pos[2];
        // Jxy
        Jij[1] += -mesh[n].pos[0] * mesh[n].pos[1];
        // Jxz
        Jij[2] += -mesh[n].pos[0] * mesh[n].pos[2];

        // Jyx
        Jij[3] += -mesh[n].pos[1] * mesh[n].pos[0];
        // Jyy
        Jij[4] +=
            mesh[n].pos[0] * mesh[n].pos[0] + mesh[n].pos[2] * mesh[n].pos[2];
        // Jyz
        Jij[5] += -mesh[n].pos[1] * mesh[n].pos[2];

        // Jzx
        Jij[6] += -mesh[n].pos[2] * mesh[n].pos[0];
        // Jzy
        Jij[7] += -mesh[n].pos[2] * mesh[n].pos[1];
        // Jzz
        Jij[8] +=
            mesh[n].pos[0] * mesh[n].pos[0] + mesh[n].pos[1] * mesh[n].pos[1];

        // xc,yc,zc
        xc += mesh[n].pos[0];
        yc += mesh[n].pos[1];
        zc += mesh[n].pos[2];
    }
    xc /= N;
    yc /= N;
    zc /= N;
    Iij[0] = Jij[0] - N * (yc * yc + zc * zc);
    Iij[1] = Jij[1] - N * (-xc * yc);
    Iij[2] = Jij[2] - N * (-xc * zc);
    Iij[3] = Jij[3] - N * (-yc * xc);
    Iij[4] = Jij[4] - N * (xc * xc + zc * zc);
    Iij[5] = Jij[5] - N * (-yc * zc);
    Iij[6] = Jij[6] - N * (-zc * xc);
    Iij[7] = Jij[7] - N * (-zc * yc);
    Iij[8] = Jij[8] - N * (xc * xc + yc * yc);
    return Iij;
}
std::vector<double> triangulation::Crr_m(int ednum) {
    // sort edge
    std::vector<double> Crr;
    std::vector<int> sort_edge_list;
    std::vector<double> posc(3, 0);
    std::vector<double> r1(3, 0), r2(3, 0);
    double r_ave = 0;
    double l_r1, l_r2;
    int next_i;

    sort_edge_list.push_back(edge_lists[ednum][0]);
    for (int i = 0; i < edge_lists[ednum].size() - 1; i++) {
        next_i = mesh[sort_edge_list.back()].edge_nei[0];
        sort_edge_list.push_back(next_i);
    }
    if (mesh[sort_edge_list.back()].edge_nei[0] != sort_edge_list[0]) {
        std::cout << "not in a circle?\n";
    }
    // measure posc: sx,yc,zc
    for (int i = 0; i < sort_edge_list.size(); i++) {
        for (int k = 0; k < 3; k++) {
            posc[k] += mesh[sort_edge_list[i]].pos[k];
        }
    }
    for (int k = 0; k < 3; k++) {
        posc[k] /= sort_edge_list.size();
    }
    // measure correlation
    Crr.clear();
    r_ave = 0;
    for (int s = 0; s < (sort_edge_list.size() + 1); s++) {
        Crr.push_back(0.0);
        for (int i = 0; i < sort_edge_list.size(); i++) {
            int j = (i + s) % sort_edge_list.size();
            // measure correlation betwee i and j, then
            for (int k = 0; k < 3; k++) {
                r1[k] = mesh[sort_edge_list[i]].pos[k] - posc[k];
                r2[k] = mesh[sort_edge_list[j]].pos[k] - posc[k];
            }
            Crr[s] += std::sqrt(innerproduct(r1, r1) * innerproduct(r2, r2));
            if (s == 0) {
                r_ave += std::sqrt(innerproduct(r1, r1));
            }
        }

        Crr[s] /= sort_edge_list.size();
    }
    r_ave /= sort_edge_list.size();
    Crr.insert(Crr.begin(), r_ave);
    return Crr;
}

void triangulation::O_reset() {
    E = 0;
    for (int e = 0; e < Ne; e++) {
        Les[e] = 0;
    }
    IdA = 0;
    I2H = 0;
    I2H2 = 0;
    Ikg = 0;
    Ik2 = 0;
    Itau = 0;
    for (int i = 0; i < mesh.size(); i++) {
        mesh[i].ds = ds_m(i);
        mesh[i].dAn2H = dAn2H_m(i);
        mesh[i].dskg = dskg_m(i);
        mesh[i].dsk2 = dsk2_m(i);
        mesh[i].dstau = dstau_m(i);
        mesh[i].dE = dE_m(i);
        Les[mesh[i].edge_num] += mesh[i].ds;
        IdA += mesh[i].dAn2H[0];
        I2H += mesh[i].dAn2H[1] * mesh[i].dAn2H[0];
        I2H2 += mesh[i].dAn2H[1] * mesh[i].dAn2H[1] * mesh[i].dAn2H[0];
        Ikg += mesh[i].dskg;
        Ik2 += mesh[i].dsk2;
        Itau += mesh[i].dstau;
        E += mesh[i].dE;
    }
    std::cout << "after reset\n";
    std::cout << "E = " << E << "\n";
}

int triangulation::update_v2flip_local(int ind_i) {
    int ind_o, i_nei_a, i_nei_b, ind_a, ind_b;
    int v2flip_ind_i = list_a_nei_b(v2flip_list, ind_i);
    if (v2flip_ind_i != -1) {
        v2flip_list.erase(v2flip_list.begin() + v2flip_ind_i);
    }
    mesh[ind_i].nei2flip.clear();
    if (mesh[ind_i].nei.size() < 4) {
        return 0;
    }
    for (int i_nei_o = 0; i_nei_o < mesh[ind_i].nei.size(); i_nei_o++) {
        ind_o = mesh[ind_i].nei[i_nei_o];
        if (list_a_nei_b(mesh[ind_i].edge_nei, ind_o) == -1 &&
            mesh[ind_o].nei.size() > 3) {
            i_nei_a = (i_nei_o + 1) % mesh[ind_i].nei.size();
            i_nei_b =
                (i_nei_o - 1 + mesh[ind_i].nei.size()) % mesh[ind_i].nei.size();
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
        return 1;
    }
    return 0;
}
int triangulation::update_v2flip_group(std::vector<int> ind_i_nei) {
    for (int i = 0; i < ind_i_nei.size(); i++) {
        update_v2flip_local(ind_i_nei[i]);
    }
    return 0;
}

void triangulation::renew_v2add_list() {
    v2add_list.clear();
    int j, j_next;
    std::vector<int> i_nei, j_nei, j_next_nei;
    std::vector<int> j_e_nei, j_next_e_nei;
    for (int i = 0; i < mesh.size(); i++) {
        if (mesh[i].edge_nei.size() == 0 && mesh[i].nei.size() >= 3) {
            i_nei = mesh[i].nei;
            for (j = 0; j < mesh[i].nei.size(); j++) {
                j_next = (j + 1) % mesh[i].nei.size();
                // if j j_next are on the edge and connected, and have more than
                // 3 nei
                if (mesh[mesh[i].nei[j]].edge_nei.size() != 0 &&
                    mesh[mesh[i].nei[j_next]].edge_nei.size() != 0 &&
                    (mesh[mesh[i].nei[j]].edge_nei[0] == mesh[i].nei[j_next] ||
                     mesh[mesh[i].nei[j]].edge_nei[1] == mesh[i].nei[j_next]) &&
                    mesh[mesh[i].nei[j]].nei.size() >= 3 &&
                    mesh[mesh[i].nei[j_next]].nei.size() >= 3) {
                    j_nei = mesh[mesh[i].nei[j]].nei;
                    j_next_nei = mesh[mesh[i].nei[j_next]].nei;
                    j_e_nei = mesh[mesh[i].nei[j]].edge_nei;
                    j_next_e_nei = mesh[mesh[i].nei[j_next]].edge_nei;
                    if (mesh[mesh[i].nei[j]].edge_nei[0] !=
                            mesh[i].nei[j_next] &&
                        mesh[mesh[i].nei[j]].edge_nei[1] !=
                            mesh[i].nei[j_next]) {
                        std::cout << "j-j_next not connected!\n";
                    }
                    v2add_list.push_back(i);
                    break;
                }
            }
        }
    }
}

void triangulation::renew_ve2remove_list() {
    ve2remove_list.clear();
    int i_e_nei0, i_e_nei1;
    for (int i = 0; i < edge_lists.size(); i++) {
        for (int j = 0; j < edge_lists[i].size(); j++) {
            i_e_nei0 = mesh[edge_lists[i][j]].edge_nei[0];
            i_e_nei1 = mesh[edge_lists[i][j]].edge_nei[1];
            if (mesh[edge_lists[i][j]].nei.size() > 2 &&
                list_a_nei_b(mesh[i_e_nei0].nei, i_e_nei1) == -1 &&
                list_a_nei_b(mesh[i_e_nei1].nei, i_e_nei0) == -1) {
                ve2remove_list.push_back(edge_lists[i][j]);
            }
        }
    }
}

int triangulation::list_a_nei_b(std::vector<int> a, int b) {
    for (int j = 0; j < a.size(); j++) {
        if (a[j] == b) {
            return j;
        }
    }
    return -1;
}

int triangulation::check_nei_connect() {
    int j_next;
    for (int i = 0; i < mesh.size(); i++) {
        for (int j = 0; j < mesh[i].nei.size(); j++) {
            j_next = 0;
        }
    }
    return 0;
}
int triangulation::check_duplication(int ind_i) {
    for (int j = 0; j < mesh[ind_i].nei.size() - 1; j++) {
        for (int k = j + 1; k < mesh[ind_i].nei.size(); k++) {
            if (mesh[ind_i].nei[j] == mesh[ind_i].nei[k]) {
                std::cout << "duplicated nei!\n";
                return 1;
            }
        }
    }
    return 0;
}