//
//  tbmd.hpp
//  GQMD-Hubbard-liquid
//
//  Created by Gia-Wei Chern on 1/11/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#ifndef tbmd_hpp
#define tbmd_hpp

#include <iostream>
#include <fstream>

#include "util.hpp"
#include "potential.hpp"
#include "gutzwiller.hpp"

const vec3 p0 = {0, 0, 0};

class SOrbitalSystem {
public:
    int numAtoms;
    
    Vec<vec3> position;
    Vec<vec3> velocity;
    Vec<vec3> force;
    Vec<vec3> force_elec;
    Vec<Vec<vec3>> force_ij;  //force_ij[i][j] means force on j from i
    Vec<Vec<vec3>> force_elec_ij;
    
    Vec<vec3> pos_init;
    Vec<vec3> v_init;
    
    Vec<double> mu_n;
    Vec<double> n0;
    
    vec3 cm_velocity;
    
    double mass;
    double onsite_U;
    double kT;
    double gamma;   // Langevin damping
    
    double filling_fraction;
    int num_filling;
    double mu;              // chemical potential
    
    double kT_elec;
    
    arma::vec eigE;
    arma::cx_mat eig_mode;
    

    GA_S_Orbital *potential;
    
    GA_S_Orbital::GA_Data atom(int i) {return potential->atom[i];};
    
    Vec<arma::cx_vec> f;
    
    RNG rng;
    std::random_device seed;
    
    int dimH;
    
    arma::cx_mat Hamiltonian;
    arma::cx_mat rho;
    
    vec3 system_box;
    double volume;

    int solu_tag;
    
    int max_GA_iter;
    double accuracy_GA_iter;
    double cp_mix_factor;
    double GA_kT_start, GA_kT_stop;
    int GA_n_annealing;
    double GA_r_annealing;
    double GA_r_lambda;
    
    // constructor:
    SOrbitalSystem(int nAtoms, double atomMass, double U_, double filling, double temperature, double gm, int solu_tag_pass);
    
    // initialization:
    void init_potential(void);
    
    // ============================================================================
    // initialization of various lattices:
    void init_bcc(double rs, int M, vec3& bdsLo, vec3& bdsHi);
    void init_sc(double rs, int M, vec3& bdsLo, vec3& bdsHi);
    void init_linear_chain(double rnn, int M, vec3& bdsLo, vec3& bdsHi);
    void init_random(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi, double rmin = 0);
    void init_random_dimers(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi);

    double lat_a;
    int lat_num_k_1D;
    
    // ============================================================================
    // set parameters:
    int boundary_type;   // 0 : Open BC,  1 : Periodic BC
    int langevin_dym;
    
    void set_BD_type(int t) {boundary_type = t;};   // default = 1
    void set_LD_dynm(int t) {langevin_dym = t;};    // default = 1
    
    // ============================================================================
    Vec<AtomPair> pairs;
    
    void find_all_pairs(void);
    void find_all_pairs(double cutoff);

    void wrap_PBC_position();

    vector<vector<double>> atom_dist, atom_dist2;
    vector<double> atom_nearest_dist;   //min distance of each atom to other atoms
    double ave_atom_nearest_dist;  //average of atom_nearest_dist
    double nearest_dist;   //min distance of all atoms

    void compute_atom_distance();
    
    // ============================================================================
    void built_Hamiltonian(void);
    void built_bare_Hamiltonian(void);
    void compute_uncorreted_densities(void);
    
    // ============================================================================
    // for single-particle density matrix:
    double compute_chemical_potential(void);
    double compute_avg_density(double x);
    arma::cx_mat compute_density_matrix(arma::cx_mat & H_elec);
    void compute_density_matrix(void);
    arma::vec fd_factor;
    
    // ============================================================================
    double copy_phi_matrix(void);
    
    // ============================================================================
    // GA iteration routines
    void set_GA_param(int max_iter, double accu, double T_start, double T_end, int n_annealing, double r_annealing, double factor, double mix);
    void iteration_GA(int & iter, double & err1, double & err2);

    void GA_solver(int & iter, double & err1, double & err2);

    void save_GA(string const filename);
    void read_GA(string const filename);
    
    // ============================================================================
    // Solving bare Hamiltonian
    void solve_bare_TB(void);
    
    // ============================================================================
    // MD-related routines
    void move_atoms(double dt);
    void compute_forces(void);
    void integrate_forces(double dt);
    void integrate_Langevin(double dt);
    void step_NVT(double dt, int & iter, double & err1, double & err2);
    
    // ============================================================================
    // compute energies
    double ee, ep, ek;
    
    double e_kin(void);
    double e_elec(void);
    double e_pair(void);
    
    // ============================================================================
    // measurements:
    double renormalization;
    double db_occupancy;
    double force_avg;
    
    double compute_avg_renormalization(void);
    double compute_avg_db_occupancy(void);
    double compute_avg_force();
    
    double avg_n, sgm_n, avg_rho_d;
    void compute_standard_dev_density(void);
    
    void compute_cm_velocity(void);
    void subtract_cm_velocity(void);
    
    double virial;
    double pressure;
    double virial_elec;
    arma::mat P_tensor;
    void compute_pressure_tensor(void);
    double compute_virial(void);
    double compute_pressure(void);
    
    Vec<double> tnn_eff;
    Vec<double> tnn_bare;
    double tnn_min, tnn_max;
    void compute_tnn_eff(void);
    
    double dist_pair, D2_pair;
    void compute_avg_pair_distance(void);
    
    
    inline double nf_ga(int i) {
        return real(conj(atom(i).f(1)) * atom(i).f(1) + conj(atom(i).f(2)) * atom(i).f(2));
    };
    inline double db_ga(int i) {
        return norm(atom(i).f(2));
    };

    double fe, s_qp, s_sb;
    double compute_elec_fe(void);
    
    // ============================================================================
    // timing routines:
    double time_diagonalization;
    double time_gw_iteration;
    void reset_time(void);
    
    // ============================================================================
    void save_configuration(string const filename, string const filename_1, double db_condition);
    void save_configuration(string const filename);
    void read_configuration(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi, string const filename);
    
    // ============================================================================
    void plot_binding_curve(string const filename);
    
    // ============================================================================
    void print_spectrum(const string filename);
    
    void print_ga_data(const string filename);
    
    void compute_forces_elec(void);

    double compute_metallic_fraction(void);

};


#endif /* tbmd_hpp */
