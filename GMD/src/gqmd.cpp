//
//  main.cpp
//  GA-QMD-Hubbard-liquid
//
//  Created by Gia-Wei Chern on 1/31/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#include <iostream>
#include <sstream>
#include "util.hpp"
#include "potential.hpp"
#include "gutzwiller.hpp"
#include "tbmd.hpp"
#include "analysis.hpp"

using namespace std;

void GQMD_simulation(MD_param param, MD_Data &accu_data, int restart_tag){
    
    double l_box = pow((4. * M_PI / 3.) * ((double) param.N_atoms), 1./3.) * param.r_s;
    int N_atoms = param.N_atoms;
    double dt = param.dlt_t;
    
    SOrbitalSystem atoms(N_atoms, param.atom_mass, param.U, param.filling_fraction, param.kT, param.Langevin_damping, param.solu_tag);
    atoms.set_BD_type(1);
    
    atoms.set_GA_param(param.GA_max_iter, param.GA_accuracy, param.GA_kT_start, param.GA_kT_stop, param.GA_n_annealing, param.GA_r_annealing, param.GA_r_lambda, param.GA_cp_mix);

    atoms.plot_binding_curve("eb.dat");

    cout << "l_box = " << l_box << endl;
    
    double lt = l_box;
    double lz = l_box;
    vec3 bd_lo, bd_hi;

    if (restart_tag == 1)
        atoms.read_configuration(lt, lz, bd_lo, bd_hi, "c.dat");  //we did not compute the force here
    if (restart_tag == 0)
        atoms.init_random(lt, lz, bd_lo, bd_hi, 0.5);
    if (restart_tag == 2)
        atoms.init_random_dimers(lt, lz, bd_lo, bd_hi);
    
    int iter;
    double err_n, err_f;
    
    atoms.init_potential();
    
    //======================================================
    MD_Data mddata;

    double r_min = 0.000001;
    double r_max = 0.25 * (lt + lz);
    double dr = 0.01;
    int num_r = (int) ((r_max - r_min)/dr);
    mddata.init_g_r(num_r, dr);

    double R_min = 0;
    double R_max = 1.002;
    int num_bin_R = 501;
    double delta_R = (R_max - R_min) / (double) num_bin_R;
    mddata.init_hist_R(num_bin_R, delta_R, R_min);

    double d_min = 0;
    double d_max = 0.30;
    int num_bin_d = 300;
    double delta_d = (d_max - d_min) / (double) num_bin_d;
    mddata.init_hist_d(num_bin_d, delta_d, d_min);

    double force_min = 0;  //electronic force
    double force_max = 15.;    //30
    int num_bin_force = 100;   //300
    double delta_force = (force_max - force_min) / (double) num_bin_force;
    mddata.init_hist_force(num_bin_force, delta_force, force_min);

    //===== d_min and force_min are dependent on init_hist_d() and init_hist_force()
    int num_bin_d_2D = 100;
    int num_bin_force_2D = 100;
    double delta_d_2D = (d_max - d_min) / (double) num_bin_d_2D;
    double delta_force_2D = (force_max - force_min) / (double) num_bin_force_2D;
    //mddata.init_hist2D_force_d(num_bin_d_2D, num_bin_force_2D, delta_d_2D, delta_force_2D);

    double rho_min = 0.;
    double rho_max = 0.6;
    int num_bin_rho = 300;
    double delta_rho = (rho_max - rho_min) / (double) num_bin_rho;
    mddata.init_hist_rho(num_bin_rho, delta_rho, rho_min);

    //===== initialization of rho_min are dependent on init_hist_rho()
    double distance_min = 0.;
    double distance_max = floor(l_box / 2 * sqrt(3));
    int num_bin_distance_2D = 100;  //300
    int num_bin_rho_2D = 100;       //300
    double delta_distance_2D = (distance_max - distance_min) / (double) num_bin_distance_2D;
    double delta_rho_2D = (rho_max - rho_min) / (double) num_bin_rho_2D;
    mddata.init_hist2D_rho_distance(num_bin_distance_2D, num_bin_rho_2D, delta_distance_2D, delta_rho_2D, distance_min);

    //===== dependent on init_hist2D_rho_distance() and init_hist_force()
    mddata.init_hist2D_fij_distance();

    int clusterSize_min = 2;
    int clusterSize_max = N_atoms;
    int num_bin_clusterSize = clusterSize_max - clusterSize_min + 1;
    mddata.init_hist_clusterSize(num_bin_clusterSize, clusterSize_min);

    mddata.init_hist_numAtomAB(N_atoms);

    //===== use the bins from init_hist2D_rho_distance()
    mddata.init_hist_nearest_neighbor();

    //========== electronic energy spectrum =========
    double E_min = -60;
    double E_max = 30;
    int num_bin_E = 900;
    double delta_E = (E_max - E_min) / num_bin_E;
    mddata.init_hist_E(num_bin_E, delta_E, E_min);

    //===============================================================
    ofstream fs("r.dat");
    
    int n_save = 20;
    int nx = 0;

    double db_condition = 0.001;
    double cluster_rho_condition = 0.1;

    int n_config = 0;

    //===========================================================
    for (int r=0; r < param.N_steps; r++) {

        atoms.step_NVT(dt, iter, err_n, err_f);
        
        double ee = atoms.e_elec();
        double ep = atoms.e_pair();
        double ek = atoms.e_kin();
        
        double rq = atoms.compute_avg_renormalization();
        double db = atoms.compute_avg_db_occupancy();
        
        mddata.compute_atom_ratio(atoms, db_condition);
        
        fs << setprecision(6);
        fs << r * dt << '\t' << ee << '\t' << ep << '\t' << ek << '\t';
        fs << rq << '\t' << db << '\t';  //No. 5
        fs << iter << '\t' << err_n << '\t' << err_f;
        fs << endl;

        //============= temperature quench ================
        if (param.quench_tag == 1){
            if (r == param.wait_steps)
                atoms.kT = param.quench_temp;

            if (r >= param.wait_steps){
                atoms.kT_elec = 2./3. * ek;
                atoms.GA_kT_start = 2./3. * ek * param.ratio_T;
                atoms.GA_kT_stop = 2./3. * ek;
            }
        }


        if (r % n_save == 0)
            cout << "r = " << r << endl;
            
        if (r % n_save == 0 && r >= param.wait_steps){

            mddata.basic_measurements(atoms);

            mddata.compute_g_r(atoms, db_condition, r_min, r_max);

            mddata.compute_hist_renormalization(atoms);

            mddata.compute_hist_double_occupancy(atoms);

            mddata.compute_hist_force(atoms);

            //mddata.compute_hist2D_force_d(atoms);

            mddata.compute_hist_E(atoms);

            mddata.compute_hist_rho(atoms);
            
            //mddata.compute_hist2D_rho_distance(atoms);

            mddata.compute_hist2D_fij_distance(atoms);

            mddata.compute_hist_clusterSize(atoms, cluster_rho_condition);

            //dependent on compute_atom_ratio() of atom_nearest_rho[i] and p_A, p_B
            //mddata.compute_hist_nearest_neighbor(atoms, db_condition);

            average_MD_Data(accu_data, mddata, nx);

            //================================================================================
            accu_data.print_g_r("gr.dat");

            accu_data.print_hist_R("hist_R.dat");

            accu_data.print_hist_d("hist_d.dat");

            accu_data.print_hist_force("hist_force.dat");

            //accu_data.print_hist2D_force_d("hist2D_force_d.dat");

            accu_data.print_hist_E("hist_E.dat");

            accu_data.print_hist_rho("hist_rho.dat");

            //accu_data.print_hist2D_rho_distance("hist2D_rho_distance.dat");

            accu_data.print_hist2D_fij_distance("hist2D_fij_distance.dat");

            accu_data.print_hist_clusterSize("hist_clusterSize.dat");

            //accu_data.print_hist_nearest_neighbor("hist_nearest_distance.dat", "hist_nearest_rho.dat");

            //dependent on compute_atom_ratio()
            //accu_data.print_hist_numAtomAB("hist_numAtomAB.dat");

            accu_data.print_data("a.dat", atoms.potential->onsiteU);

            nx++;

            atoms.save_configuration("c.dat");

            atoms.print_ga_data("q.dat");

            atoms.print_spectrum("em.dat");
        }

    }
    
    fs.close();
}

int main(int argc, const char * argv[]) {

    MD_param sim1;
    
    sim1.U    = argc > 1 ? atof(argv[1]) : 0;
    sim1.kT   = argc > 2 ? atof(argv[2]) : 0.2; 
    double ratio_T      = argc > 3 ? atof(argv[3]) : 1.15;
    sim1.GA_r_lambda    = argc > 4 ? atof(argv[4]) : 0.175;
    sim1.GA_cp_mix      = argc > 5 ? atof(argv[5]) : 0.87;
    
    sim1.N_atoms    = argc > 6 ? atoi(argv[6]) : 100;
    sim1.r_s        = argc > 7 ? atof(argv[7]) : 1.0;
    sim1.N_steps    = argc > 8 ? atoi(argv[8]) : 100;
    sim1.wait_steps = argc > 9 ? atoi(argv[9]) : 40;     //steps waiting for equilibrium
    
    int restart_tag  = argc > 10 ? atoi(argv[10]) : 0;         //turn on: 1
    sim1.quench_tag  = argc > 11 ? atoi(argv[11]) : 0;         //turn on: 1
    sim1.quench_temp = argc > 12 ? atof(argv[12]) : 0.02;      //useless if quench_switch is off

    cout << "U = " << sim1.U << endl;
    cout << "T_ion = " << sim1.kT << endl;
    cout << "r_lambda = " << sim1.GA_r_lambda << endl;
    cout << "N_atoms = " << sim1.N_atoms << endl;
    
    sim1.filling_fraction = 0.5;
    sim1.atom_mass = 5.0;
    sim1.Langevin_damping = 0.05;
    sim1.dlt_t = 0.02;
    sim1.solu_tag = 1;   //1 for GA solution, 0 for bare-t solution, 2 for classical bare-phi potential
    
    sim1.GA_kT_start = sim1.kT * ratio_T; // 1500.;
    sim1.GA_kT_stop  = sim1.kT;           // 1400.;
    sim1.ratio_T = ratio_T;
    sim1.GA_max_iter = 2000;       //1000;
    sim1.GA_accuracy = 1.e-5;      //1.e-5; //this is for err1, err2 = GA_accuracy * 1.e-N is stricter
    sim1.GA_n_annealing = 50;
    sim1.GA_r_annealing = 0.98;    //0.98
    
    MD_Data mdt;
    GQMD_simulation(sim1, mdt, restart_tag);

    return 0;
}


