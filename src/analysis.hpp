//
//  analysis.hpp
//  GQMD-Hubbard-liquid
//
//  Created by Gia-Wei Chern on 1/11/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#ifndef analysis_hpp
#define analysis_hpp

#include <stdio.h>
#include "util.hpp"
#include "potential.hpp"
#include "gutzwiller.hpp"
#include "tbmd.hpp"

using namespace std;

class cluster_site{
    public:

    int parent;
    int rank;
};

class MD_param {
public:
    
    double U;
    double r_s;
    double kT;
    
    double atom_mass;
    double Langevin_damping;
    
    double filling_fraction;
    int N_atoms;
    
    int solu_tag;

    int GA_max_iter;
    double GA_accuracy;
    double GA_r_lambda;
    double GA_cp_mix;
    
    
    double GA_kT_start;
    double GA_kT_stop;
    double ratio_T;
    
    int GA_n_annealing;
    double GA_r_annealing;
    
    double dlt_t;
    int N_steps;

    int wait_steps;  //steps waiting for equilibrium
    
    int quench_tag;
    double quench_temp;

    int autocorr_numSteps;
    int transp_atomic_tag;
    int transp_elec_tag;
};


class MD_Data {
public:
    
    double M1_ee, M2_ee;
    double M1_ep, M2_ep;
    double M1_ek, M2_ek;
    double M1_rq, M2_rq;
    double M1_dc, M2_dc;
    double M1_force, M2_force;
    double M1_e_tot, M2_e_tot;
    
    double M1_pr, M2_pr;        // pressure
    double M1_vl, M2_vl;
    
    double M1_vol, M2_vol;
    double M1_D, M2_D;
    double M1_K_atomic, M2_K_atomic;

    double M1_ave_atom_nearest_dist, M2_ave_atom_nearest_dist;

    vector<double> M1_beta_L11;
    vector<double> M1_beta_L22;
    vector<double> M1_beta_L12;
    vector<double> M1_K_elec;
    vector<double> M1_Seeback;
    
    double M1_ne, M2_ne;
    double M1_rho, M2_rho;
    
    double M1_fc, M2_fc;
    
    MD_Data(void) {
        reset();
    };
    
    ~MD_Data(void) {
        
    };
    
    void reset(void) {
        M1_ee = M2_ee = M1_ep = M2_ep = M1_ek = M2_ek = 0;
        M1_rq = M2_rq = M1_dc = M2_dc = M1_force = M2_force = 0;
        M1_pr = M2_pr = M1_vl = M2_vl = 0;
        M1_vol = M2_vol = 0;
        M1_D = M2_D = 0;
        M1_ne = M2_ne = M1_rho = M2_rho = 0;
        M1_fc = M2_fc = 0;
        M1_K_atomic = M2_K_atomic = 0;

        M1_ave_atom_nearest_dist = M2_ave_atom_nearest_dist = 0;

        M1_e_tot = M2_e_tot = 0;

        ratio_dimer = ratio_dimer2 = 0;
        ratio_nonmott_atomic = ratio_nonmott_atomic2 = 0;
        p_A = p_A2 = 0;
        p_B = p_B2 = 0;
        p_AB = 0;
        N_atomA = N_atomB = 0;

    };



    //====================== ratio of different atoms ====================
    double ratio_dimer, ratio_dimer2;
    double ratio_nonmott_atomic, ratio_nonmott_atomic2; 
    double p_A, p_A2;   //ratio of nonmott
    double p_B, p_B2;   //ratio of mott
    double p_AB;
    int N_atomA, N_atomB;
    double ratio_dimer_by_dist;

    vector<double> atom_nearest_rho;  //max rho_ij of each atom to other atoms
    
    void compute_atom_ratio(SOrbitalSystem & system, double db_condition);

    //==========================================================
    Vec<double> g_r;
    int num_r;
    double dr;

    Vec<double> gr_AA, gr_AB, gr_BB;
    Vec<double> gr_AA_1, gr_AB_1, gr_BB_1;
    
    void init_g_r(int size, double dlt_r) {
        num_r = size;
        dr = dlt_r;

        g_r = Vec<double>(size, 0);
        gr_AA = Vec<double>(size, 0);
        gr_AB = Vec<double>(size, 0);
        gr_BB = Vec<double>(size, 0);

        gr_AA_1 = Vec<double>(size, 0);
        gr_AB_1 = Vec<double>(size, 0);
        gr_BB_1 = Vec<double>(size, 0);
    };

    void compute_g_r(SOrbitalSystem & system, double db_condition, double r_min, double r_max);
    
    void print_g_r(string const filename) {
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m<g_r.size(); m++)
            fs << (m + 0.5) * dr << '\t' << g_r[m] << '\t' << gr_AA[m] << '\t' << gr_AB[m] << '\t' << gr_BB[m] << '\t' << gr_AA_1[m] << '\t' << gr_AB_1[m] << '\t' << gr_BB_1[m] << endl;
        fs.close();
    };

    Vec<double> hs_R;
    int num_bin_R;
    double delta_R;
    double R_min;

    void init_hist_R(int size, double delta_R0, double R_min0) {
        num_bin_R = size;
        delta_R = delta_R0;
        R_min = R_min0;
        
        hs_R = Vec<double>(num_bin_R, 0);
    };

    void compute_hist_renormalization(SOrbitalSystem & system);
    
    void print_hist_R(string const filename) {
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m<hs_R.size(); m++) fs << R_min + (m + 0.5) * delta_R << '\t' << hs_R[m] << endl;
        fs.close();
        
    };

    Vec<double> hs_d;
    int num_bin_d;
    double delta_d;
    double d_min;
    
    void init_hist_d(int size, double delta_d0, double d_min0) {
        num_bin_d = size;
        delta_d = delta_d0;
        d_min = d_min0;

        hs_d = Vec<double>(num_bin_d, 0);
    };
    
    void compute_hist_double_occupancy(SOrbitalSystem & system);

    void print_hist_d(string const filename) {
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m< hs_d.size(); m++) fs << d_min + (m + 0.5) * delta_d << '\t' << hs_d[m] << endl;
        fs.close();
    };

    Vec<double> hs_force;
    int num_bin_force;
    double delta_force;
    double force_min;
    
    void init_hist_force(int size, double delta_force0, double force_min0) {
        num_bin_force = size;
        delta_force = delta_force0;
        force_min = force_min0;

        hs_force = Vec<double>(num_bin_force, 0);
    };
    
    void compute_hist_force(SOrbitalSystem & system);

    void print_hist_force(string const filename) {
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m< hs_force.size(); m++) fs << force_min + (m + 0.5) * delta_force << '\t' << hs_force[m] << endl;
        fs.close();
    };

    vector<vector<double>> hs2D_force_d;
    int num_bin_d_2D;
    int num_bin_force_2D;
    double delta_d_2D;
    double delta_force_2D;

    void init_hist2D_force_d(int num_bin_d_2D0, int num_bin_force_2D0, double delta_d_2D0, double delta_force_2D0){
        num_bin_d_2D = num_bin_d_2D0;
        num_bin_force_2D = num_bin_force_2D0;
        delta_d_2D = delta_d_2D0;
        delta_force_2D = delta_force_2D0;

        hs2D_force_d = vector<vector<double>>(num_bin_d_2D, vector<double>(num_bin_force_2D, 0));
    };
    
    void compute_hist2D_force_d(SOrbitalSystem & system);

    void print_hist2D_force_d(string const filename) {
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m < hs2D_force_d.size(); m++){
        for (unsigned n=0; n < hs2D_force_d[m].size(); n++){
            fs << d_min + (m + 0.5) * delta_d_2D << '\t' << force_min + (n + 0.5) * delta_force_2D << '\t' << hs2D_force_d[m][n] << endl;
        }
        }

        fs.close();
    };

    Vec<double> hs_rho;
    int num_bin_rho;
    double delta_rho;
    double rho_min;

    void init_hist_rho(int size, double delta_rho0, double rho_min0) {
        num_bin_rho = size;
        delta_rho = delta_rho0;
        rho_min = rho_min0;

        hs_rho = Vec<double>(num_bin_rho, 0);
    };

    void compute_hist_rho(SOrbitalSystem & system);
    
    void print_hist_rho(string const filename) {
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m<hs_rho.size(); m++) fs << rho_min + (m + 0.5) * delta_rho << '\t' << hs_rho[m] << endl;
        fs.close();
    };

    vector<vector<double>> hs2D_rho_distance;
    int num_bin_distance_2D;
    int num_bin_rho_2D;
    double delta_distance_2D;
    double delta_rho_2D;
    double distance_min;

    void init_hist2D_rho_distance(int num_bin_distance_2D0, int num_bin_rho_2D0, double delta_distance_2D0, double delta_rho_2D0, double distance_min0){
        num_bin_distance_2D = num_bin_distance_2D0;
        num_bin_rho_2D = num_bin_rho_2D0;
        delta_distance_2D = delta_distance_2D0;
        delta_rho_2D = delta_rho_2D0;
        distance_min = distance_min0;

        hs2D_rho_distance = vector<vector<double>>(num_bin_distance_2D, vector<double>(num_bin_rho_2D, 0));
    };

    void compute_hist2D_rho_distance(SOrbitalSystem & system);

    void print_hist2D_rho_distance(string const filename) {
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        
        for (int n=0; n < num_bin_rho_2D; n++){
        for (int m=0; m < num_bin_distance_2D; m++){
            fs << distance_min + (m + 0.5) * delta_distance_2D << '\t' << rho_min + (n + 0.5) * delta_rho_2D << '\t' << hs2D_rho_distance[m][n] << endl;
        }
        }

        fs.close();
    };


    vector<vector<double>> hs2D_fij_distance;
    
    void init_hist2D_fij_distance(){
        hs2D_fij_distance = vector<vector<double>>(num_bin_distance_2D, vector<double>(num_bin_force, 0));
    };

    void compute_hist2D_fij_distance(SOrbitalSystem & system);

    void print_hist2D_fij_distance(string const filename){
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        
        for (int n=0; n < num_bin_force; n++){
        for (int m=0; m < num_bin_distance_2D; m++){
            fs << distance_min + (m + 0.5) * delta_distance_2D << '\t' << force_min + (n + 0.5) * delta_force << '\t' << hs2D_fij_distance[m][n] << endl;
        }
        }

        fs.close();
    };

    Vec<double> hs_clusterSize;
    int num_bin_clusterSize;
    int clusterSize_min;
    double num_clusters, num_clusters2;   //equals num_roots, but specifically for accu_data in average_MD_Data(), so it is a double, not an int
    double num_singular_roots, num_singular_roots2;

    void init_hist_clusterSize(int num_bin_clusterSize0, int clusterSize_min0) {
        num_bin_clusterSize = num_bin_clusterSize0;
        clusterSize_min = clusterSize_min0;

        hs_clusterSize = Vec<double>(num_bin_clusterSize, 0);

        num_clusters = num_clusters2 = 0;
        
        num_singular_roots = num_singular_roots2 = 0;
    };

    void compute_hist_clusterSize(SOrbitalSystem & system, double cluster_rho_condition);

    void print_hist_clusterSize(string const filename){
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m<hs_clusterSize.size(); m++) fs << clusterSize_min + m << '\t' << hs_clusterSize[m] << endl;
        fs.close();
    };
    void cluster_analysis(SOrbitalSystem & system, double cluster_rho_condition);

    //============ union-find algorithm to find metallic clusters ==========
    int num_sites;
    Vec<cluster_site> site;
    
    int num_roots;   //equals number of clusters
    vector<int> site_of_root;
    vector<int> num_in_root;
    vector<vector<int>> clusters;

    void init_sites(int num_sites0);
    int uf_find(int i);
    void uf_union(int i, int j);
    void union_find(arma::mat rho, double cluster_rho_condition);

    //===============================================================
    int num_bin_numAtomAB;
    Vec<double> hs_numAtomA, hs_numAtomB;

    void init_hist_numAtomAB(int N_atoms) {
        num_bin_numAtomAB = N_atoms + 1;
        hs_numAtomA = Vec<double>(num_bin_numAtomAB, 0);
        hs_numAtomB = Vec<double>(num_bin_numAtomAB, 0);
    };

    //this is dependent on compute_atom_ratio() of N_atomA and N_atomB
    void print_hist_numAtomAB(string const filename){
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m<hs_numAtomA.size(); m++)
            fs << m << '\t' << hs_numAtomA[m] << '\t' << hs_numAtomB[m] << endl;
        fs.close();
    };

    vector<double> hs_nearest_distance, hs_nearest_distance_mott;
    vector<double> hs_nearest_rho, hs_nearest_rho_mott;
    
    //used the bins from init_hist2D_rho_distance()
    void init_hist_nearest_neighbor(){
        hs_nearest_distance = vector<double>(num_bin_distance_2D, 0);
        hs_nearest_distance_mott = vector<double>(num_bin_distance_2D, 0);
        hs_nearest_rho = vector<double>(num_bin_rho_2D, 0);
        hs_nearest_rho_mott = vector<double>(num_bin_rho_2D, 0);
    }
    
    void compute_hist_nearest_neighbor(SOrbitalSystem & system, double db_condition);

    void print_hist_nearest_neighbor(string const filename_1, string const filename_2){
        std::ofstream fs;
        
        fs.open(filename_1.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m<hs_nearest_distance.size(); m++)
            fs << distance_min + (m + 0.5) * delta_distance_2D << '\t' << hs_nearest_distance[m] << '\t' << hs_nearest_distance_mott[m] << endl;
        fs.close();

        fs.open(filename_2.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m<hs_nearest_rho.size(); m++)
            fs << rho_min + (m + 0.5) * delta_rho_2D << '\t' << hs_nearest_rho[m] << '\t' << hs_nearest_rho_mott[m] << endl;
        fs.close();
    };

    Vec<double> hs_E;
    int num_bin_E;
    double delta_E;
    double E_min;
    
    void init_hist_E(int size, double delta_E_p, double E_min_p) {
        num_bin_E = size;
        delta_E = delta_E_p;
        E_min = E_min_p;

        hs_E = Vec<double>(num_bin_E, 0);
    };
    
    void compute_hist_E(SOrbitalSystem & system);

    void print_hist_E(string const filename) {
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for (unsigned m=0; m< hs_E.size(); m++) 
            fs << E_min + (m + 0.5) * delta_E << '\t' << hs_E[m] << endl;
        fs.close();
    };


    void print_data(string const filename, double param) {
        
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        //fs.precision(12);
        fs << setprecision(6);
        
        fs << param << '\t';

        fs << M1_ee << '\t' << sqrt(M2_ee - pow(M1_ee, 2)) << '\t'; //NO.2
        fs << M1_ep << '\t' << sqrt(M2_ep - pow(M1_ep, 2)) << '\t';
        fs << M1_ek << '\t' << sqrt(M2_ek - pow(M1_ek, 2)) << '\t';
        fs << M1_rq << '\t' << sqrt(M2_rq - pow(M1_rq, 2)) << '\t';
        fs << M1_dc << '\t' << sqrt(M2_dc - pow(M1_dc, 2)) << '\t';

        fs << M1_force << '\t' << sqrt(M2_force - pow(M1_force, 2)) << '\t'; //NO. 12
        fs << M1_pr << '\t' << sqrt(M2_pr - pow(M1_pr, 2)) << '\t';
        fs << M1_vol << '\t' << sqrt(M2_vol - pow(M1_vol, 2)) << '\t';
        fs << M1_fc << '\t' << sqrt(M2_fc - pow(M1_fc, 2)) << '\t';
        fs << M1_ne << '\t' << sqrt(M2_ne - pow(M1_ne, 2)) << '\t';

        fs << endl;

        fs.close();
    };

    void print_elec_transp(string const filename, vector<double> param){
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs << setprecision(6);

        for (int i=0; i < param.size(); i++)
            fs << param[i] << '\t' << M1_beta_L11[i] << '\t' << M1_beta_L22[i] << '\t' << M1_beta_L12[i] << '\t' << M1_K_elec[i] << '\t' << M1_Seeback[i] << endl;

        fs.close();
    }
    
    void basic_measurements(SOrbitalSystem & system);
};

void average_MD_Data(MD_Data & accu_data, MD_Data & mddata, int nc);

#endif /* analysis_hpp */
