#ifndef _TRIANGULATION_H
#define _TRIANGULATION_H

#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>
class triangulation {
  public:
    // eternal parameters
    double beta; // system temperature
    int N;       // number of beads
    int Ne;      // number of edges
    int Lr;      // initial number of bead per line
    int L;       // vertical height of cylinder initialization
    // also used to set fixed distance betwen two beads
    double l0;   // tether maximum length
    double l1;   //  edge tether maximum length
    double kar;  // mean curvature bending stiffness
    double C0;   // spontaneous mean curvature
    double karg; // gaussian curvature bending stiffness
    double lam;  // line tension coefficient
    double B;    // line bending stiffness
    double Bc;   // geodesic torsion stiffness
    double tau0; // nature torsion

    // system configuration
    struct vertex {
        double ds;
        // double dA2H2;
        // seperate into |2H| and dA
        std::vector<double> dAn2H;
        double dskg;
        double dsk2;
        double dstau;
        double dE;
        // double dRg2; // radius of gyration square, has nothing to do with the
        // energy
        std::vector<double> pos{0, 0, 0}; // position (x,y,z)
        std::vector<double> nc{0, 0, 0};  // normal for edge beads
        std::vector<int> nei;             // index of neighbors
        // have to make sure nei[k+1],nei[k] are connected!!!
        // except for the boundary?
        std::vector<int> edge_nei;
        // neighbors form edge with this one (if exist)
        int edge_num;              // which edge
        std::vector<int> nei2flip; // index in nei that could be flipped with
    };

    std::vector<vertex> mesh; // the mesh
    std::vector<std::vector<int>>
        edge_lists; // list of the edge beads sorted, linked
    std::vector<std::pair<int, int>> bond_list;
    // bond_list[i] is a pair of two bead connected in [bulk!}
    // notice only bond in the bulk are included, and there will be on
    // repetition due to symmetry,
    std::vector<int> fixed_beads; // beads can't be moved
    /*
    struct v2flip {
        int mesh_index;                // index of the vertex in mesh
        std::vector<int> flip_nei_ind; // index of nei that can flip bond with
        // flip_nei_ind +-1 in mesh_index.nei are a, b
    };
    std::vector<v2flip> v2flip_list; // list of vertex having bond to flip
    */
    std::vector<int> v2flip_list;
    std::vector<int> v2add_list;     // list of vertex can be added to the edge
    std::vector<int> ve2remove_list; // list of edge vertex can be removed

    // observable
    double E;                // system energy
    std::vector<double> Les; // length of the boundary (edge)
    double IdA;              // integral of dA
    double I2H;              // integral of dA(2H)
    double I2H2;             // integral of dA(2h)^2
    double Ikg;              // integral of ds kg
    double Ik2;              // integral of ds k^2
    double Itau;             // integral of ds tau

    // randomnumber generators
    std::mt19937 gen;
    std::uniform_int_distribution<> rand_pos;  // random position
    std::uniform_real_distribution<> rand_uni; // uniform distribution

    // initialization
    triangulation(double beta_, int N_, int Ne_, int L_, double d0_, double l0_,
                  double l1_, double kar_, double C0_, double karg_,
                  double lam_, double B_, double Bc_, double tau0_);

    // put the beads and bonds in to position accordingly
    void init_rhombus_shape(double d0_);
    void init_cylinder_shape(double d0_);

    void reset_config();
    void push_neis_back(int i, std::vector<int> nei_dist);
    void push_eneis_back(int i, std::vector<int> enei_dist);
    void push_bneis_list(int i, std::vector<int> bnei_dist);

    // bond in bulk to bond_list
    // d0_ initial vertex-vertex distance
    // other parameters are the same as in the class

    // local energy-related measurement
    // _m stand for measurement
    double ds_m(int index);                 // length of the local edge index
    double dA2H2_m(int index);              // (2H)^2 dA for mean curvature
    std::vector<double> dAn2H_m(int index); // measure and set dA and |2H|
    double dskg_m(int index);               // kg*ds, for gaussian
    double dsk2_m(int index);               // k^2*ds for edge bending
    std::vector<double> nc_m(int index);    // surface normal
    double dstau_m(int index);              // geodesic torsion * ds
    double dE_m(int index);                 // energy of the local vertex
    // double dRg2_m(int index); // R_g^2, just for the initial measurement

    // correlation function measurement
    std::vector<double> Cncnc_m(int ednum);
    // return correlation of edge(enum) surface normal
    // moment of inertia matrix Iij measurement
    std::vector<double> Iij_m();
    // edge displacement correlation measurement
    std::vector<double> Crr_m(int ednum);

    void O_reset();
    double distance(int ind_1, int ind_2);
    int energy_check(); // check if the energy set are correct, for debug use

    void delete_bond_list(int ind_i, int ind_j);
    // delete bond i-j, including both <i,j> and <j,i> in the bulk bond_list
    void renew_v2add_list();
    void renew_ve2remove_list();
    int update_v2flip_local(int ind_i);
    int update_v2flip_group(std::vector<int> ind_i_nei);

    // Monte Carlo updates
    int vertex_metropolis(double delta_s);
    // execute vertex move update in [-ds,ds]^3
    // return 1: accepted, 0: rejected

    int bond_metropolis();
    // execute bond switch update
    // return 1: accepted, 0: rejected

    int edge_metropolis();
    // execute bond remove/add update
    // return 1: accepted, 0: rejected

    // experiment
    void State_write(std::string filename);
    // write system state to file
    void shape_set(double theta);
    // set cap-like shape for testing
    void State_write_seq(std::string filename, int MC_sweeps, int step_p_sweep,
                         double ds_);

    void shape_set_helix(double q);
    // set helical shape for testing

    void State_load(std::string state_file);
    // load system state from file
    void Thermal(int MC_sweeps, int step_p_sweep, int beta_steps, double ds_);
    // thermalisation of the system, starting from beta=0 (hot start)
    void O_MC_measure(int MC_sweeps, int step_p_Cm, int step_p_sweep, double ds,
                      std::string filename, std::string filename_C,
                      std::string filename_Cr, std::string filename_I);
    // measure the obserables
    void E_vs_lam(std::vector<double> lam_measure, int MC_steps_init,
                  int MC_steps_measure, std::string filename);
    // energy versus lambda curve testing

    // little tool
    int sort_nei(int index);
    int list_a_nei_b(std::vector<int> a, int b);
    int check_nei_connect();
    int check_duplication(int ind_i);
};
#endif
