//
//  analysis.cpp
//  GQMD-Hubbard-liquid
//
//  Created by Gia-Wei Chern on 1/11/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "util.hpp"
#include "potential.hpp"
#include "gutzwiller.hpp"
#include "tbmd.hpp"
#include "analysis.hpp"

void average_MD_Data(MD_Data & accu_data, MD_Data & mddata, int nc) {
    
    if (nc == 0) {
        cout << "resetting accumulated MD_data." << endl;
        accu_data.reset();

        accu_data.init_g_r(mddata.num_r, mddata.dr);

        accu_data.init_hist_R(mddata.num_bin_R, mddata.delta_R, mddata.R_min);

        accu_data.init_hist_d(mddata.num_bin_d, mddata.delta_d, mddata.d_min);

        accu_data.init_hist_force(mddata.num_bin_force, mddata.delta_force, mddata.force_min);
        
        //accu_data.init_hist2D_force_d(mddata.num_bin_d_2D, mddata.num_bin_force_2D, mddata.delta_d_2D, mddata.delta_force_2D);

        accu_data.init_hist_rho(mddata.num_bin_rho, mddata.delta_rho, mddata.rho_min);

        accu_data.init_hist2D_rho_distance(mddata.num_bin_distance_2D, mddata.num_bin_rho_2D, mddata.delta_distance_2D, mddata.delta_rho_2D, mddata.distance_min);

        accu_data.init_hist2D_fij_distance();

        accu_data.init_hist_clusterSize(mddata.num_bin_clusterSize, mddata.clusterSize_min);

        accu_data.init_hist_nearest_neighbor();

        accu_data.init_hist_numAtomAB(mddata.num_bin_numAtomAB);

        accu_data.init_hist_E(mddata.num_bin_E, mddata.delta_E, mddata.E_min);
    }

    //nc means there are number of counts before this average
    accu_data.M1_ee = (accu_data.M1_ee * nc + mddata.M1_ee) / (nc + 1.);
    accu_data.M2_ee = (accu_data.M2_ee * nc + mddata.M2_ee) / (nc + 1.);
    
    accu_data.M1_ep = (accu_data.M1_ep * nc + mddata.M1_ep) / (nc + 1.);
    accu_data.M2_ep = (accu_data.M2_ep * nc + mddata.M2_ep) / (nc + 1.);
    
    accu_data.M1_ek = (accu_data.M1_ek * nc + mddata.M1_ek) / (nc + 1.);
    accu_data.M2_ek = (accu_data.M2_ek * nc + mddata.M2_ek) / (nc + 1.);

    accu_data.M1_e_tot = (accu_data.M1_e_tot * nc + mddata.M1_e_tot) / (nc + 1.);
    accu_data.M2_e_tot = (accu_data.M2_e_tot * nc + mddata.M2_e_tot) / (nc + 1.);
    
    accu_data.M1_rq = (accu_data.M1_rq * nc + mddata.M1_rq) / (nc + 1.);
    accu_data.M2_rq = (accu_data.M2_rq * nc + mddata.M2_rq) / (nc + 1.);
    
    accu_data.M1_dc = (accu_data.M1_dc * nc + mddata.M1_dc) / (nc + 1.);
    accu_data.M2_dc = (accu_data.M2_dc * nc + mddata.M2_dc) / (nc + 1.);

    accu_data.M1_force = (accu_data.M1_force * nc + mddata.M1_force) / (nc + 1.);
    accu_data.M2_force = (accu_data.M2_force * nc + mddata.M2_force) / (nc + 1.);
    
    accu_data.M1_pr = (accu_data.M1_pr * nc + mddata.M1_pr) / (nc + 1.);
    accu_data.M2_pr = (accu_data.M2_pr * nc + mddata.M2_pr) / (nc + 1.);
    
    accu_data.M1_vl = (accu_data.M1_vl * nc + mddata.M1_vl) / (nc + 1.);
    accu_data.M2_vl = (accu_data.M2_vl * nc + mddata.M2_vl) / (nc + 1.);
    
    accu_data.M1_vol = (accu_data.M1_vol * nc + mddata.M1_vol) / (nc + 1.);
    accu_data.M2_vol = (accu_data.M2_vol * nc + mddata.M2_vol) / (nc + 1.);

    accu_data.M1_fc = (accu_data.M1_fc * nc + mddata.M1_fc) / (nc + 1.);
    accu_data.M2_fc = (accu_data.M2_fc * nc + mddata.M2_fc) / (nc + 1.);
        
    accu_data.M1_ne = (accu_data.M1_ne * nc + mddata.M1_ne) / (nc + 1.);
    accu_data.M2_ne = (accu_data.M2_ne * nc + mddata.M2_ne) / (nc + 1.);
    
    accu_data.M1_rho = (accu_data.M1_rho * nc + mddata.M1_rho) / (nc + 1.);
    accu_data.M2_rho = (accu_data.M2_rho * nc + mddata.M2_rho) / (nc + 1.);

    accu_data.M1_ave_atom_nearest_dist = (accu_data.M1_ave_atom_nearest_dist * nc + mddata.M1_ave_atom_nearest_dist) / (nc + 1.);
    accu_data.M2_ave_atom_nearest_dist = (accu_data.M2_ave_atom_nearest_dist * nc + mddata.M2_ave_atom_nearest_dist) / (nc + 1.);

    for (unsigned m=0; m<accu_data.g_r.size(); m++){
        accu_data.g_r[m] = (accu_data.g_r[m] * nc + mddata.g_r[m]) / (nc + 1.);

        // Special weighted average for partial pair correlation function
        if (pow((nc * accu_data.p_A + mddata.p_A), 2) / (nc + 1.) > 1.E-8)
            accu_data.gr_AA[m] = (accu_data.gr_AA[m] * nc * pow(accu_data.p_A, 2) + mddata.gr_AA[m] * mddata.p_A2) / (nc * pow(accu_data.p_A, 2) + mddata.p_A2); //(pow((nc * accu_data.p_A + mddata.p_A), 2) / (nc + 1.));
        
        if ((nc * accu_data.p_A + mddata.p_A) * (nc * accu_data.p_B + mddata.p_B) / (nc + 1.) > 1.E-8)
            accu_data.gr_AB[m] = (accu_data.gr_AB[m] * nc * accu_data.p_A * accu_data.p_B + mddata.gr_AB[m] * mddata.p_AB) / (nc * accu_data.p_A * accu_data.p_B + mddata.p_AB); // ((nc * accu_data.p_A + mddata.p_A) * (nc * accu_data.p_B + mddata.p_B) / (nc + 1.));
        
        if (pow((nc * accu_data.p_B + mddata.p_B), 2) / (nc + 1.) > 1.E-8)
            accu_data.gr_BB[m] = (accu_data.gr_BB[m] * nc * pow(accu_data.p_B, 2) + mddata.gr_BB[m] * mddata.p_B2) / (nc * pow(accu_data.p_B, 2) + mddata.p_B2); // (pow((nc * accu_data.p_B + mddata.p_B), 2) / (nc + 1.));

        //different way of normalization, test them
        if (nc * accu_data.p_A2 + mddata.p_A2 > 1.E-8)
            accu_data.gr_AA_1[m] = (accu_data.gr_AA_1[m] * nc * accu_data.p_A2 + mddata.gr_AA[m] * mddata.p_A2) / (nc * accu_data.p_A2 + mddata.p_A2);
        
        if (nc * accu_data.p_AB + mddata.p_AB > 1.E-8)
            accu_data.gr_AB_1[m] = (accu_data.gr_AB_1[m] * nc * accu_data.p_AB + mddata.gr_AB[m] * mddata.p_AB) / (nc * accu_data.p_AB + mddata.p_AB);
        
        if (nc * accu_data.p_B2 + mddata.p_B2 > 1.E-8)
            accu_data.gr_BB_1[m] = (accu_data.gr_BB_1[m] * nc * accu_data.p_B2 + mddata.gr_BB[m] * mddata.p_B2) / (nc * accu_data.p_B2 + mddata.p_B2);
    }
    //accu_data.p_A and accu_data.p_B is updated in the buttom, they still have usage
    accu_data.p_A2 = (nc * accu_data.p_A2 + mddata.p_A2) / (nc + 1.);
    accu_data.p_AB = (nc * accu_data.p_AB + mddata.p_AB) / (nc + 1.);
    accu_data.p_B2 = (nc * accu_data.p_B2 + mddata.p_B2) / (nc + 1.);
    
    for (unsigned m=0; m<accu_data.hs_R.size(); m++)
        accu_data.hs_R[m] = (accu_data.hs_R[m] * nc + mddata.hs_R[m]) / (nc + 1.);

    for (unsigned m=0; m<accu_data.hs_d.size(); m++)
        accu_data.hs_d[m] = (accu_data.hs_d[m] * nc + mddata.hs_d[m]) / (nc + 1.);

    for (unsigned m=0; m<accu_data.hs_force.size(); m++) {
        accu_data.hs_force[m] = (accu_data.hs_force[m] * nc + mddata.hs_force[m]) / (nc + 1.);
    }

    for (unsigned m=0; m<accu_data.hs_E.size(); m++) {
        accu_data.hs_E[m] = (accu_data.hs_E[m] * nc + mddata.hs_E[m]) / (nc + 1.);
    }

    /* for (unsigned m=0; m<accu_data.hs2D_force_d.size(); m++) {
    for (unsigned n=0; n<accu_data.hs2D_force_d[m].size(); n++) {
        accu_data.hs2D_force_d[m][n] = (accu_data.hs2D_force_d[m][n] * nc + mddata.hs2D_force_d[m][n]) / (nc + 1.);
    }
    } */

    for (unsigned m=0; m<accu_data.hs_rho.size(); m++) {
        accu_data.hs_rho[m] = (accu_data.hs_rho[m] * nc + mddata.hs_rho[m]) / (nc + 1.);
    }

    /* for (unsigned m=0; m<accu_data.hs2D_rho_distance.size(); m++){
    for (unsigned n=0; n<accu_data.hs2D_rho_distance[m].size(); n++){
        accu_data.hs2D_rho_distance[m][n] = (accu_data.hs2D_rho_distance[m][n] * nc + mddata.hs2D_rho_distance[m][n]) / (nc + 1.);
    }
    } */

    for (unsigned m=0; m<accu_data.hs2D_fij_distance.size(); m++){
    for (unsigned n=0; n<accu_data.hs2D_fij_distance[m].size(); n++){
        accu_data.hs2D_fij_distance[m][n] = (accu_data.hs2D_fij_distance[m][n] * nc + mddata.hs2D_fij_distance[m][n]) / (nc + 1.);
    }
    }

    for (unsigned m=0; m<accu_data.hs_clusterSize.size(); m++){
        accu_data.hs_clusterSize[m] = (accu_data.hs_clusterSize[m] * nc + mddata.hs_clusterSize[m]) / (nc + 1.);   //for unnormalized hs_clusterSize
        
        //if (mddata.num_roots > 0)   //for normalized hs_clusterSize
        //accu_data.hs_clusterSize[m] = (accu_data.hs_clusterSize[m] * nc * accu_data.num_clusters + mddata.hs_clusterSize[m] * mddata.num_roots) / (nc * accu_data.num_clusters + mddata.num_roots);
    }
    accu_data.num_clusters = (accu_data.num_clusters * nc + mddata.num_roots) / (nc + 1.);
    accu_data.num_clusters2 = (accu_data.num_clusters2 * nc + mddata.num_roots * mddata.num_roots) / (nc + 1.);

    accu_data.num_singular_roots = (accu_data.num_singular_roots * nc + mddata.num_singular_roots) / (nc + 1.);
    accu_data.num_singular_roots2 = (accu_data.num_singular_roots2 * nc + mddata.num_singular_roots2) / (nc + 1.);

    /* //===============================================
    for (unsigned m=0; m<accu_data.hs_numAtomA.size(); m++){
        accu_data.hs_numAtomA[m] = accu_data.hs_numAtomA[m] * nc / (nc + 1.);
        accu_data.hs_numAtomB[m] = accu_data.hs_numAtomB[m] * nc / (nc + 1.);
    }
    accu_data.hs_numAtomA[mddata.N_atomA] += 1. / (nc + 1.);
    accu_data.hs_numAtomB[mddata.N_atomB] += 1. / (nc + 1.); */

    //================================================
    accu_data.ratio_dimer = (accu_data.ratio_dimer * nc + mddata.ratio_dimer) / (nc + 1.);
    accu_data.ratio_dimer2 = (accu_data.ratio_dimer2 * nc + mddata.ratio_dimer2) / (nc + 1.);
    
    accu_data.ratio_nonmott_atomic = (accu_data.ratio_nonmott_atomic * nc + mddata.ratio_nonmott_atomic) / (nc + 1.);
    accu_data.ratio_nonmott_atomic2 = (accu_data.ratio_nonmott_atomic2 * nc + mddata.ratio_nonmott_atomic2) / (nc + 1.);

    /* for (unsigned m=0; m<accu_data.hs_nearest_distance.size(); m++){
        if (nc * accu_data.p_A + mddata.p_A > 1.E-8)
            accu_data.hs_nearest_distance[m] = (accu_data.hs_nearest_distance[m] * nc * accu_data.p_A + mddata.hs_nearest_distance[m] * mddata.p_A) / (nc * accu_data.p_A + mddata.p_A);
        if (nc * accu_data.p_B + mddata.p_B > 1.E-8)
            accu_data.hs_nearest_distance_mott[m] = (accu_data.hs_nearest_distance_mott[m] * nc * accu_data.p_B + mddata.hs_nearest_distance_mott[m] * mddata.p_B) / (nc * accu_data.p_B + mddata.p_B);
    }
    for (unsigned m=0; m<accu_data.hs_nearest_rho.size(); m++){
        if (nc * accu_data.p_A + mddata.p_A > 1.E-8)
            accu_data.hs_nearest_rho[m] = (accu_data.hs_nearest_rho[m] * nc * accu_data.p_A + mddata.hs_nearest_rho[m] * mddata.p_A) / (nc * accu_data.p_A + mddata.p_A);
        if (nc * accu_data.p_B + mddata.p_B > 1.E-8)
            accu_data.hs_nearest_rho_mott[m] = (accu_data.hs_nearest_rho_mott[m] * nc * accu_data.p_B + mddata.hs_nearest_rho_mott[m] * mddata.p_B) / (nc * accu_data.p_B + mddata.p_B);
    } */

    accu_data.p_A = (accu_data.p_A * nc + mddata.p_A) / (nc + 1.);
    accu_data.p_B = (accu_data.p_B * nc + mddata.p_B) / (nc + 1.);
}

void MD_Data::basic_measurements(SOrbitalSystem & system) {
    
    double ee = system.e_elec();
    double ep = system.e_pair();
    double ek = system.e_kin();
    double rq = system.compute_avg_renormalization();
    double dc = system.compute_avg_db_occupancy();
    double force = system.compute_avg_force();
    
    double pr = system.compute_pressure();
    double vl = system.virial;
    
    double vol = system.volume;
    
    system.compute_standard_dev_density();
    double ne = system.avg_n;
    double rho = system.avg_rho_d;
    
    double frac = system.compute_metallic_fraction();

    M1_ave_atom_nearest_dist = system.ave_atom_nearest_dist;
    M2_ave_atom_nearest_dist = pow(system.ave_atom_nearest_dist, 2);


    M1_ee = ee;
    M2_ee = pow(ee, 2);
    
    M1_ep = ep;
    M2_ep = pow(ep, 2);
    
    M1_ek = ek;
    M2_ek = pow(ek, 2);

    M1_e_tot = ee + ep + ek;
    M2_e_tot = pow(M1_e_tot, 2);
    
    M1_rq = rq;
    M2_rq = pow(rq, 2);
    
    M1_dc = dc;
    M2_dc = pow(dc, 2);

    M1_force = force;
    M2_force = pow(force, 2);
    
    M1_pr = pr;
    M2_pr = pow(pr, 2);
    
    M1_vl = vl;
    M2_vl = pow(vl, 2);
    
    M1_vol = vol;
    M2_vol = pow(vol, 2);
    
    M1_fc = frac;
    M2_fc = pow(frac, 2);
    
    M1_ne = ne;
    M2_ne = pow(ne, 2);
    
    M1_rho = rho;
    M2_rho = pow(rho, 2);
    
}

void MD_Data::compute_atom_ratio(SOrbitalSystem & system, double db_condition){
    ratio_dimer = ratio_dimer2 = 0;
    ratio_nonmott_atomic = ratio_nonmott_atomic2 = 0;
    ratio_dimer_by_dist = 0;

    atom_nearest_rho = vector<double>(system.numAtoms, 0);

    for (int i=0; i<system.numAtoms; i++){
    for (int j=0; j<system.numAtoms; j++){
        if (j!=i && abs(system.rho(i, j)) > atom_nearest_rho[i])
            atom_nearest_rho[i] = abs(system.rho(i, j));
    }
    }

    for (int i=0; i<system.numAtoms; i++){
    if (system.db_ga(i) > db_condition){
        int count = 0;
        for (int j=0; j<system.numAtoms; j++)
            if (j!=i && abs(system.rho(i, j)) > atom_nearest_rho[i] - 0.24) count++;

        if (count == 1) ratio_dimer++;
        else ratio_nonmott_atomic++;
    }
    }

    for (int i=0; i<system.numAtoms; i++){
    if (system.db_ga(i) > db_condition){
        if (system.atom_nearest_dist[i] < 1.25) ratio_dimer_by_dist++;
    }
    }
    ratio_dimer_by_dist /= system.numAtoms;

    N_atomA = ratio_dimer + ratio_nonmott_atomic;
    N_atomB = system.numAtoms - N_atomA;

    ratio_dimer /= system.numAtoms;
    ratio_nonmott_atomic /= system.numAtoms;

    p_A = ratio_dimer + ratio_nonmott_atomic;
    p_B = 1. - p_A;

    ratio_dimer2 = ratio_dimer * ratio_dimer;
    ratio_nonmott_atomic2 = ratio_nonmott_atomic * ratio_nonmott_atomic;
    p_A2 = p_A * p_A;
    p_B2 = p_B * p_B;
    p_AB = p_A * p_B;
}

void MD_Data::compute_g_r(SOrbitalSystem & system, double db_condition, double r_min, double r_max){
    
    for (int m=0; m<num_r; m++){
        g_r[m] = 0;
        gr_AA[m] = gr_AB[m] = gr_BB[m] = 0;
    }
    
    for (int i=0; i<system.numAtoms; i++){
    for (int j = i+1; j<system.numAtoms; j++){
        double rij = system.atom_dist[i][j];
        if (rij < r_max && rij > r_min) {
            int m = (int) ((rij - r_min)/dr);
            
            if(m<num_r && m>=0){   
                g_r[m]++;
                
                if (system.db_ga(i) > db_condition && system.db_ga(j) > db_condition) gr_AA[m]++;
                else if (system.db_ga(i) <= db_condition && system.db_ga(j) <= db_condition) gr_BB[m]++;
                else gr_AB[m]++;
            }
        }
    }
    }
    
    for (int m=0; m<num_r; m++) {
        double rm = (m + 0.5) * dr;
        g_r[m] *= system.volume / (2. * 3.1415926 * dr * pow(system.numAtoms * rm, 2));

        // Be careful with zero cases
        // If one species vanishes, then the pair correlation function involve that species should also vanish,
        // therefore no need for normalization.
        // This kind of case can be due phase transition or finite size effect.
        if (N_atomA != 0)
            gr_AA[m] *= system.volume / (2. * 3.1415926 * dr * pow(N_atomA * rm, 2));
        if (N_atomA * N_atomB != 0)
            gr_AB[m] *= system.volume / (4. * 3.1415926 * dr * N_atomA * N_atomB * pow(rm, 2));
        if (N_atomB != 0)
            gr_BB[m] *= system.volume / (2. * 3.1415926 * dr * pow(N_atomB * rm, 2));
    }
}

void MD_Data::compute_hist_renormalization(SOrbitalSystem & system) {
    for (int m=0; m<num_bin_R; m++) hs_R[m] = 0;
    int N_tot = 0;
    
    for (int i=0; i<system.numAtoms; i++) {
        N_tot++;
        double R_m = norm(system.atom(i).Rq);   //in C++ and armadillo, norm() actually gives norm^2 
        int m = (int) ((R_m - R_min) / delta_R);
        
        if (m>=0 && m<num_bin_R) {
            hs_R[m]++;
        }
    }
    
    for(int m=0; m<num_bin_R; m++) {
        hs_R[m] /= (N_tot * delta_R);
    }
}

void MD_Data::compute_hist_double_occupancy(SOrbitalSystem & system) {
    for (int m=0; m<num_bin_d; m++) hs_d[m] = 0;
    int N_tot = 0;
    
    for (int i=0; i<system.numAtoms; i++) {
        N_tot++;
        double d_m = norm(system.atom(i).f(2));
        int m = (int) ((d_m - d_min) / delta_d);
        
        if (m>=0 && m<num_bin_d) {
            hs_d[m]++;
        }
    }
    
    for (int m=0; m<num_bin_d; m++) {
        hs_d[m] /= ((N_tot + 1.e-5) * delta_d);
    }
}

void MD_Data::compute_hist_force(SOrbitalSystem & system) {
    for (int m=0; m<num_bin_force; m++) hs_force[m] = 0;
    int N_tot = 0;
    
    for (int i=0; i<system.numAtoms; i++) {
        N_tot++;
        double force_m = system.force_elec[i].norm();
        int m = (int) ((force_m - force_min) / delta_force);
        
        if (m>=0 && m<num_bin_force) {
            hs_force[m]++;
        }
    }
    
    for (int m=0; m<num_bin_force; m++) {
        hs_force[m] /= (N_tot * delta_force);
    }
}

void MD_Data::compute_hist_E(SOrbitalSystem & system) {
    for (int m=0; m<num_bin_E; m++) hs_E[m] = 0;
    int N_tot = 0;

    for (int i=0; i<system.dimH; i++){
        N_tot++;
        double E_m = system.eigE(i) - system.mu;
        int m = (int) ((E_m - E_min) / delta_E);
        
        if (m>=0 && m<num_bin_E) {
            hs_E[m]++;
        }
    }

    for (int m=0; m<num_bin_E; m++) {
        hs_E[m] /= ((N_tot + 1.e-8) * delta_E);
    }
}

void MD_Data::compute_hist2D_force_d(SOrbitalSystem & system){
    hs2D_force_d = vector<vector<double>>(num_bin_d_2D, vector<double>(num_bin_force_2D, 0));
    int N_tot = 0;

    for (int i=0; i<system.numAtoms; i++){
        N_tot++;
        double d_m = system.db_ga(i);
        double force_n = system.force_elec[i].norm();

        int m = (int) ((d_m - d_min) / delta_d_2D);
        int n = (int) ((force_n - force_min) / delta_force_2D);
        
        if (m>=0 && m<num_bin_d_2D){
        if (n>=0 && n<num_bin_force_2D){
            hs2D_force_d[m][n]++;
        }
        }
    }

    for (int m=0; m<num_bin_d_2D; m++){
    for (int n=0; n<num_bin_force_2D; n++){
        hs2D_force_d[m][n] /= ((N_tot + 1.e-6) * delta_d_2D * delta_force_2D);
    }
    }

}

void MD_Data::compute_hist_rho(SOrbitalSystem & system) {
    
    for (int m=0; m<num_bin_rho; m++) hs_rho[m] = 0;
    int N_tot = 0;
    
    for (int i=0; i<system.numAtoms; i++){
    for (int j=0; j<system.numAtoms; j++){
    if (i != j){
        N_tot++;

        double rho_m = abs(system.rho(i, j));
        int m = (int) ((rho_m - rho_min) / delta_rho);

        if(m>=0 && m<num_bin_rho) {    
            hs_rho[m]++;
        }
    }
    }
    }
    
    for (int m=0; m<num_bin_rho; m++) {
        hs_rho[m] /= (N_tot * delta_rho);
    }
}

void MD_Data::compute_hist2D_rho_distance(SOrbitalSystem & system){
    hs2D_rho_distance = vector<vector<double>>(num_bin_distance_2D, vector<double>(num_bin_rho_2D, 0));
    int N_tot = 0;

    for (int i=0; i<system.numAtoms; i++){
    for (int j=i+1; j<system.numAtoms; j++){
        N_tot++;
        double distance_m = system.atom_dist[i][j];
        double rho_n = abs(system.rho(i, j));

        int m = (int) ((distance_m - distance_min) / delta_distance_2D);
        int n = (int) ((rho_n - rho_min) / delta_rho_2D);
        
        if (m>=0 && m < num_bin_distance_2D){
        if (n>=0 && n < num_bin_rho_2D){
            hs2D_rho_distance[m][n]++;
        }
        }
    }
    }

    for (int m=0; m < num_bin_distance_2D; m++){
    for (int n=0; n < num_bin_rho_2D; n++){
        hs2D_rho_distance[m][n] /= ((N_tot + 1.e-6) * delta_distance_2D * delta_rho_2D);
    }
    }
}

void MD_Data::compute_hist2D_fij_distance(SOrbitalSystem & system){
    hs2D_fij_distance = vector<vector<double>>(num_bin_distance_2D, vector<double>(num_bin_force, 0));
    int N_tot = 0;

    for (int i=0; i<system.numAtoms; i++){
    for (int j=i+1; j<system.numAtoms; j++){
        N_tot++;
        double distance_m = system.atom_dist[i][j];
        double fij_n = system.force_elec_ij[i][j].norm();

        int m = (int) ((distance_m - distance_min) / delta_distance_2D);
        int n = (int) ((fij_n - force_min) / delta_force);
        
        if (m>=0 && m < num_bin_distance_2D){
        if (n>=0 && n < num_bin_force){
            hs2D_fij_distance[m][n]++;
        }
        }
    }
    }

    for (int m=0; m < num_bin_distance_2D; m++){
    for (int n=0; n < num_bin_force; n++){
        hs2D_fij_distance[m][n] /= ((N_tot + 1.e-6) * delta_distance_2D * delta_force);
    }
    }
}

void MD_Data::compute_hist_clusterSize(SOrbitalSystem & system, double cluster_rho_condition){

    init_sites(system.numAtoms);
    union_find(abs(system.rho), cluster_rho_condition);

    for (int i=0; i<num_bin_clusterSize; i++) hs_clusterSize[i] = 0;

    if (num_roots > 0){
        for (int i=0; i<num_roots; i++)
            hs_clusterSize[ num_in_root[i] - clusterSize_min ]++;

        //for (int m=0; m<num_bin_clusterSize; m++) hs_clusterSize[m] /= (double) num_roots;  //normalizing hs_clusterSize, but do not forget to change the accumulation data in average_MD_Data()
    }
}

void MD_Data::init_sites(int num_sites0){
    num_sites = num_sites0;
    site = Vec<cluster_site>(num_sites);

    for (int i=0; i<num_sites; i++){
        site[i].parent = i;
        site[i].rank = 0;
    }
}

int MD_Data::uf_find(int i){
    if (site[i].parent != i){
        site[i].parent = uf_find(site[i].parent);
    }
    return site[i].parent;
}

void MD_Data::uf_union(int i, int j){
    int i_root = uf_find(i);
    int j_root = uf_find(j);

    if (i_root != j_root) {
        if (site[i_root].rank < site[j_root].rank) {
            site[i_root].parent = j_root;
        } 
        else if (site[i_root].rank > site[j_root].rank) {
            site[j_root].parent = i_root;
        } 
        else {
            site[j_root].parent = i_root;
            site[i_root].rank++;
        }
    }
}

void MD_Data::union_find(arma::mat rho, double cluster_rho_condition){
    for (int i=0; i<num_sites; i++){
    for (int j = i+1; j<num_sites; j++){
    if (rho(i, j) > cluster_rho_condition){
        uf_union(i, j);
    }
    }
    }

    num_roots = 0;
    site_of_root = vector<int>(num_sites, -1);
    
    for (int i=0; i<num_sites; i++){
        if (uf_find(i) == i){
            num_roots++;
            site_of_root[num_roots - 1] = i;
        }
    }
    site_of_root.resize(num_roots);

    num_in_root = vector<int>(num_roots, 0);
    clusters = vector<vector<int>>(num_roots, vector<int>(0, 0));
    for (int r=0; r<num_roots; r++){
    for (int i=0; i<num_sites; i++){
        if (uf_find(i) == site_of_root[r]){
            num_in_root[r]++;
            clusters[r].push_back(i);
        }
    }
    }

    //===== exclude roots that only contain themself =====
    num_singular_roots = 0;
    num_singular_roots2 = 0;
    for (int r = num_roots - 1; r >= 0; r--){  //erase reverse so that erase one elements does not change the positions of unscanned elements
    if (num_in_root[r] == 1){
        num_singular_roots++;
        num_roots--;
        site_of_root.erase(site_of_root.begin() + r);
        num_in_root.erase(num_in_root.begin() + r);
        clusters.erase(clusters.begin() + r);
    }
    }
    num_singular_roots /= (double) num_sites;
    num_singular_roots2 = num_singular_roots * num_singular_roots;

    /* cout << "Number of clusters: " << num_roots << ", roots are:";
    for (unsigned i=0; i<site_of_root.size(); i++){
        cout << " " << site_of_root[i];
    }
    cout << endl;

    for (unsigned i=0; i<clusters.size(); i++){
        cout << "cluster " << i << ", " << clusters[i].size() << " sites: ";
        for (unsigned j=0; j<clusters[i].size(); j++){
            cout << " " << clusters[i][j];
        }
        cout << endl;
    } */

}

//dependent on compute_atom_ratio() of atom_nearest_rho[i] and p_A, p_B
void MD_Data::compute_hist_nearest_neighbor(SOrbitalSystem & system, double db_condition){
    hs_nearest_distance = vector<double>(num_bin_distance_2D, 0);
    hs_nearest_distance_mott = vector<double>(num_bin_distance_2D, 0);
    hs_nearest_rho = vector<double>(num_bin_rho_2D, 0);
    hs_nearest_rho_mott = vector<double>(num_bin_rho_2D, 0);

    for (int i=0; i<system.numAtoms; i++){
        int m = (int) ((system.atom_nearest_dist[i] - distance_min) / delta_distance_2D);
        if (m>=0 && m<num_bin_distance_2D){
            if (system.db_ga(i) > db_condition) hs_nearest_distance[m]++;
            else hs_nearest_distance_mott[m]++;
        }
    }
    for (int m=0; m<num_bin_distance_2D; m++){
        if (p_A > 1.E-8) hs_nearest_distance[m] /= (system.numAtoms * p_A * delta_distance_2D);
        if (p_B > 1.E-8) hs_nearest_distance_mott[m] /= (system.numAtoms * p_B * delta_distance_2D);
    }

    //=========================================
    for (int i=0; i<system.numAtoms; i++){
        int m = (int) ((atom_nearest_rho[i] - rho_min) / delta_rho_2D);

        if (m>=0 && m<num_bin_rho_2D){
            if (system.db_ga(i) > db_condition) hs_nearest_rho[m]++;
            else hs_nearest_rho_mott[m]++;
        }
    }
    for (int m=0; m<num_bin_rho_2D; m++){
        if (p_A > 1.E-8) hs_nearest_rho[m] /= (system.numAtoms * p_A * delta_rho_2D);
        if (p_B > 1.E-8) hs_nearest_rho_mott[m] /= (system.numAtoms * p_B * delta_rho_2D);
    }
}




