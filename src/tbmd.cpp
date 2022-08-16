//
//  tbmd.cpp
//  GQMD-Hubbard-liquid
//
//  Created by Gia-Wei Chern on 1/11/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cassert>

#include "util.hpp"
#include "potential.hpp"
#include "gutzwiller.hpp"
#include "tbmd.hpp"

SOrbitalSystem::SOrbitalSystem(int nAtoms, double atomMass, double _U, double filling, double temperature, double gm, int solu_tag_pass) {
    numAtoms = nAtoms;
    position = Vec<vec3>(nAtoms);
    velocity = Vec<vec3>(nAtoms);
    force = Vec<vec3>(nAtoms);
    force_elec = Vec<vec3>(nAtoms);
    force_ij = Vec<Vec<vec3>>(nAtoms, Vec<vec3>(nAtoms, p0));
    force_elec_ij = force_ij;
    
    mu_n = Vec<double>(nAtoms);
    n0 = Vec<double>(nAtoms);
    
    f = Vec<arma::cx_vec>(nAtoms);
    for (int i=0; i<nAtoms; i++) 
        f[i] = arma::cx_vec(potential->numSOrbitalGAParam);
    
    mass = atomMass;
    kT = temperature;
    gamma = gm;
    onsite_U = _U;

    solu_tag = solu_tag_pass;
    
    ee = ep = ek = 0;
    renormalization = 0;
    db_occupancy = 0;
    
    kT_elec = kT;
    
    boundary_type = 1;     // 0 : Open BC,  1 : Periodic BC
    langevin_dym = 1;
    
    filling_fraction = filling;
    num_filling = (int) (filling_fraction * nAtoms);
    
    dimH = nAtoms;

    eigE = arma::vec(dimH);
    eig_mode = arma::cx_mat(dimH, dimH);
    fd_factor = arma::vec(dimH);
    
    //    Hamiltonian = arma::sp_cx_mat(nAtoms, nAtoms);
    //    rho = arma::sp_cx_mat(nAtoms, nAtoms);
    
    
    potential = new(GA_S_Orbital)(nAtoms, onsite_U);
        
    //potential->kT = temperature;
    
    rng = RNG(seed());
    
    P_tensor.resize(3,3);
    
}

void SOrbitalSystem::init_bcc(double rnn, int M, vec3& bdsLo, vec3& bdsHi) {
    position.resize(numAtoms);
    
    double a0 = (2./sqrt(3.)) * rnn;
    
    bdsLo = {0, 0, 0};
    bdsHi = {M*a0, M*a0, M*a0};
    
    int ir = 0;
    for (int z=0; z<M; z++){
    for (int y=0; y<M; y++){
    for (int x=0; x<M; x++){
                
        position[ir] = vec3 {x*a0, y*a0, z*a0};
        position[ir+1] = position[ir] + vec3 {0.5*a0, 0.5*a0, 0.5*a0};
        
        ir+=2;
    }
    }
    }
    
    
    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box.x * system_box.y * system_box.z;
    cout << "volumn = " << volume << endl;
    
    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }
    
    double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    for(int i=0; i<numAtoms; i++) {
        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
    }
    
    pos_init = Vec<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::init_sc(double rnn, int M, vec3& bdsLo, vec3& bdsHi) {
    position.resize(numAtoms);
    
    lat_a = rnn;
    lat_num_k_1D = M;
    
    bdsLo = {0, 0, 0};
    bdsHi = {M * lat_a, M * lat_a, M * lat_a};
    
    int ir = 0;
    for (int z=0; z<M; z++){
    for (int y=0; y<M; y++){
    for (int x=0; x<M; x++){
                
        position[ir] = vec3 {x * lat_a, y * lat_a, z * lat_a};
        
        ir++;
    }
    }
    }
    
    
    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box(0) * system_box(1) * system_box(2);
    cout << "volumn = " << volume << endl;
    
    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }
    
    double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;
    for(int i=0; i<numAtoms; i++) {
        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
    }
    
    pos_init = Vec<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::init_linear_chain(double rnn, int L, vec3& bdsLo, vec3& bdsHi) {
    position.resize(numAtoms);
    
    double a0 = rnn;
    
    bdsLo = {0, 0, 0};
    bdsHi = {L*a0, a0, a0};
    
    for(int x=0; x<L; x++) position[x] = vec3 {x * a0, 0, 0};
    
    std::normal_distribution<double> rd;
    
    double dm = 0.0001 * a0;
    for(int i=0; i<L; i++) position[i] += vec3 {dm * rd(rng), dm * rd(rng), dm * rd(rng) };
    
    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box(0) * system_box(1) * system_box(2);
    cout << "volumn = " << volume << endl;
    
    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }
    
    double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;
    for(int i=0; i<numAtoms; i++) {
        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
    }
    
    pos_init = Vec<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::init_random(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi, double rmin) {
    std::uniform_real_distribution<double> rd;
    
    bdsLo = {0, 0, 0};
    bdsHi = {Lt, Lt, Lz};
    position.resize(numAtoms);

    system_box = bdsHi - bdsLo;
    //cout << "system box = " << system_box << endl;
    volume = system_box(0) * system_box(1) * system_box(2);
    //cout << "volumn = " << volume << endl;
    
    for (int i=0; i<numAtoms; i++){
        double rij = Lt;
        do {
            position[i] = vec3 {rd(rng) * Lt, rd(rng) * Lt, rd(rng) * Lz};

            rij = Lt;
            for (int j=0; j<i && rij >= rmin; j++) {
                if (boundary_type == 0) rij = (position[j] - position[i]).norm();
                if (boundary_type == 1) rij = (wrapDelta(position[j] - position[i], system_box)).norm();
            }                
        }
        while (rij < rmin);
    }

    for (int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }
    
    /* double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    for(int i=0; i<numAtoms; i++) {
        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
    } */
    
    pos_init = Vec<vec3>(numAtoms);
    for (int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::init_random_dimers(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi) {
    std::uniform_real_distribution<double> rd;
    
    bdsLo = {0, 0, 0};
    bdsHi = {Lt, Lt, Lz};
    position.resize(numAtoms);
    
    double bond_length = 0.95;
    double bond_x = rd(rng);
    double bond_y = rd(rng);
    double bond_z = rd(rng);
    double norm_xyz = sqrt(norm(bond_x) + norm(bond_y) + norm(bond_z));

    bond_x = bond_x / norm_xyz * bond_length;
    bond_y = bond_y / norm_xyz * bond_length;
    bond_z = bond_z / norm_xyz * bond_length;

    for (int i=0; i<numAtoms; i+=2){
        position[i] = vec3 {rd(rng) * Lt, rd(rng) * Lt, rd(rng) * Lz};
        position[i+1] = position[i] + vec3{bond_x, bond_y, bond_z};
    }
    
    system_box = bdsHi - bdsLo;
    volume = system_box(0) * system_box(1) * system_box(2);
    
    for (int i=0; i<numAtoms; i++){
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }
    
    pos_init = Vec<vec3>(numAtoms);
    for (int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::find_all_pairs(void) {
    if (boundary_type == 0)
        pairs = allPairs(position, potential->rcut());
    else
        pairs = allPairsPeriodic(position, potential->rcut(), system_box);
}

void SOrbitalSystem::find_all_pairs(double cutoff) {
    if (boundary_type == 0)
        pairs = allPairs(position, cutoff);
    else
        pairs = allPairsPeriodic(position, cutoff, system_box);
}

void SOrbitalSystem::wrap_PBC_position(){
    for (int i=0; i<numAtoms; i++) 
        position[i] = wrapPosition(position[i], p0, system_box);
}

void SOrbitalSystem::compute_atom_distance(){
    atom_dist = atom_dist2 = vector<vector<double>>(numAtoms, vector<double>(numAtoms, 0));

    for (int i=0; i<numAtoms; i++){
    for (int j=0; j<numAtoms; j++){
        if (boundary_type == 0) atom_dist2[i][j] = (position[j] - position[i]).norm2();
        if (boundary_type == 1) atom_dist2[i][j] = (wrapDelta(position[j] - position[i], system_box)).norm2();

        atom_dist[i][j] = sqrt(atom_dist2[i][j]);
    }
    }

    atom_nearest_dist = vector<double>(numAtoms, 0.1);
    ave_atom_nearest_dist = 0;

    for (int i=0; i<numAtoms; i++){
        if (i!=0)
            atom_nearest_dist[i] = atom_dist[i][0];
        else
            atom_nearest_dist[i] = atom_dist[i][1];

        for (int j=0; j<numAtoms; j++){
            if (j!=i && atom_dist[i][j] < atom_nearest_dist[i])
                atom_nearest_dist[i] = atom_dist[i][j];
        }

        ave_atom_nearest_dist += atom_nearest_dist[i] / numAtoms;
    }

    nearest_dist = *min_element(atom_nearest_dist.begin(), atom_nearest_dist.end());
}

void SOrbitalSystem::built_Hamiltonian(void) {
    Hamiltonian = potential->build_Hamiltonian(numAtoms, pairs);
    //Hamiltonian = potential->buildHamiltonian(numAtoms, pairs, system_box);
}

void SOrbitalSystem::built_bare_Hamiltonian(void) {
    Hamiltonian = potential->barePotential.build_Hamiltonian(numAtoms, pairs);
}

double SOrbitalSystem::compute_avg_density(double x) {
    double sum = 0;
    for (unsigned i=0; i<eigE.size(); i++) {
        sum += fermi_density(eigE(i), kT_elec, x);
    }
    return sum / ((double) numAtoms);
}

// this requires eigE in ascending order!!
double SOrbitalSystem::compute_chemical_potential(void) {
    
    double x1 = eigE(0);
    double x2 = eigE(eigE.size() - 1);
    
    int max_bisection = 100;
    double eps_bisection = 1.e-12;
    
    int iter = 0;
    while (iter < max_bisection || fabs(x2 - x1) > eps_bisection) {
        
        double xm = 0.5*(x1 + x2);
        double density = compute_avg_density(xm);
        
        if (density <= filling_fraction) x1 = xm;
        else x2 = xm;
        
        iter++;
    }
    
    return 0.5*(x1 + x2);
}

arma::cx_mat SOrbitalSystem::compute_density_matrix(arma::cx_mat & H_elec) {
    
    arma::cx_mat rho_elec(dimH, dimH);
    
    arma::cx_mat eigvec;
    auto Hd = H_elec; //fkpm::sparseToDense(H_elec);
    arma::eig_sym(eigE, eigvec, Hd);

    eig_mode = eigvec;
    
    mu = compute_chemical_potential();
    
    for (int i=0; i<dimH; i++) fd_factor(i) = fermi_density(eigE(i), kT_elec, mu);

    rho_elec.zeros();
    for(int a=0; a<dimH; a++)
    for(int b=a; b<dimH; b++) {

        cx_double sum = 0;
        for(int m=0; m<eigE.size(); m++) {
            sum += fd_factor(m) * conj(eigvec(a, m)) * eigvec(b, m);
        }
        rho_elec(a, b) = sum;
        if(a != b) rho_elec(b, a) = conj(sum);
    }
    
    return rho_elec;
}

void SOrbitalSystem::compute_density_matrix(void) {
    rho = compute_density_matrix(Hamiltonian);
}

void SOrbitalSystem::compute_forces(void) {

    potential->force(rho, pairs, force, force_ij, virial);
    potential->force_elec(rho, pairs, force_elec, force_elec_ij, virial_elec);
    
    //potential->force(potential->numSpins() * rho, pairs, force, virial, system_box);  // old version CHECK !!!!!!
    //the numSpins() is absorbed into the calculation in potential->force(...) !!!
    
}

void SOrbitalSystem::compute_forces_elec(void) {
    potential->force_elec(rho, pairs, force_elec, force_elec_ij, virial_elec);
}


void SOrbitalSystem::set_GA_param(int max_iter, double accu, double T_start, double T_stop, int n_annealing, double r_annealing, double factor, double mix) {
    max_GA_iter = max_iter;
    accuracy_GA_iter = accu;
    GA_kT_start = T_start;
    GA_kT_stop  = T_stop;
    GA_n_annealing = n_annealing;
    GA_r_annealing = r_annealing;
    GA_r_lambda = factor;
    
    cp_mix_factor = mix;
}

void SOrbitalSystem::compute_uncorreted_densities(void) {
    
    arma::cx_mat H0 = potential->barePotential.build_Hamiltonian(numAtoms, pairs);
    
    arma::cx_mat rho0 = compute_density_matrix(H0);
    
    for (int a=0; a<numAtoms; a++) {
        std::complex<double> tmp = rho0(a, a);
        n0[a] = real(tmp);
    }
    
    double sum1 = 0;
    
    for (int i=0; i<dimH; i++) {
        sum1 += eigE(i) * fermi_density(eigE(i), kT_elec, mu);
    }
}

void SOrbitalSystem::init_potential(void) {
    find_all_pairs();
    compute_atom_distance();

    potential->compute_neighbor_list(pairs);
    potential->reset_GA_parameters();
    built_Hamiltonian();
    
    compute_density_matrix();

    potential->init_var_param(rho);

}

void SOrbitalSystem::solve_bare_TB(void) {
    potential->reset_GA_parameters();
    built_Hamiltonian();
    compute_density_matrix();
}

void SOrbitalSystem::compute_tnn_eff(void) {
    int numOrbs = potential->numOrbitalsPerSite();
    arma::mat h(numOrbs, numOrbs);
    
    tnn_eff.clear();
    tnn_bare.clear();
    
    tnn_min = 1.e10;
    tnn_max = -1.e10;
    
    for (auto const& p : pairs) {
        
        if (p.delta.norm2() < pow(potential->rcut(), 2)) {
            int m = p.index1;
            int n = p.index2;
            assert(m < n);
            
            h = potential->hoppingMatrix(p.delta.norm());
            double tnn = real(h(0,0) * atom(m).Rq * atom(n).Rq);
            double tnn_b = h(0,0);
            
            tnn_eff.push_back(tnn);
            tnn_bare.push_back(tnn_b);
            
            if(tnn < tnn_min) tnn_min = tnn;
            if(tnn > tnn_max) tnn_max = tnn;
        }
    }
}

double SOrbitalSystem::copy_phi_matrix(void) {
    double sum = 0;
    
    for(int i=0; i<numAtoms; i++) {
        
        double sm = 0;
        for(int k=0; k<potential->numSOrbitalGAParam; k++) {
            sm += norm(f[i](k)) - norm(atom(i).f(k));
        }
        
        sum += fabs(sm);
        
        f[i] = atom(i).f;
    }
    return sum / ((double) numAtoms);
}

void SOrbitalSystem::iteration_GA(int & iter, double & err1, double & err2) {
    
    err1 = 100;
    err2 = 100;

    iter = 0;
    
    //Crucial: whether to re-start the GA iteration from uncorrelated state:
    
    potential->reset_GA_parameters();
    built_Hamiltonian();
    compute_density_matrix();
    potential->init_var_param(rho);
    
    kT_elec = GA_kT_start;
    
    while ((iter < max_GA_iter && (err1 > accuracy_GA_iter || err2 > accuracy_GA_iter * 1.e-3)) || (kT_elec > GA_kT_stop)){
        
        potential->compute_renormalizations(rho);
        built_Hamiltonian();
        compute_density_matrix();

        potential->compute_Deltas(rho);
        
        potential->compute_var_param(rho, cp_mix_factor, kT_elec);

        err1 = potential->adjust_lambda(rho, GA_r_lambda);
        err2 = copy_phi_matrix();

        iter++;

        if (iter % GA_n_annealing == 0 && kT_elec > GA_kT_stop) {
            kT_elec *= GA_r_annealing;
        }

        //potential->compute_renormalizations(rho);
        
        //====================== for testing ===========================
        /* if(iter % 50 == 0) {
            ofstream fs("t.dat");
            fs.precision(15);
            for(int i=0; i<numAtoms; i++) {
                //fs << i << '\t' << real(atom(i).Rq) << '\t';
                //fs << real(atom(i).f(0)) << '\t' << real(atom(i).f(1)) << '\t' << real(atom(i).f(2)) << '\t';
                //fs << imag(atom(i).f(0)) << '\t' << imag(atom(i).f(1)) << '\t' << imag(atom(i).f(2)) << '\t';
                fs << real(rho(i, i)) << '\t' << real(conj(atom(i).f(1)) * atom(i).f(1) + conj(atom(i).f(2)) * atom(i).f(2)) << '\t' << atom(i).lambda << '\t' << norm(atom(i).f(2)) << endl;
            }
            fs.close();
        } */
       
    }

    kT_elec = GA_kT_stop;
}

void SOrbitalSystem::save_GA(string const filename){
    std::ofstream fs;
    
    fs.open(filename.c_str(), ios::out);
    fs.precision(6);

    for (int i=0; i<numAtoms; i++){
        fs << real(atom(i).Rq) << '\t' << atom(i).lambda << endl;
    }

    fs.close();
}


void SOrbitalSystem::read_GA(string const filename){
    std::ifstream fs;
    fs.open(filename.c_str(), ios::in);

    for (int i=0; i<numAtoms; i++) {
        fs >> potential->atom[i].Rq >> potential->atom[i].lambda;
    }
}


void SOrbitalSystem::move_atoms(double dt) {
    for(int i=0; i<numAtoms; i++) {
        vec3 dlt = 0.5 * (force[i]/mass) * dt;
        dlt += velocity[i];
        velocity[i] = dlt;  // velocity at t + 0.5*dt
        dlt *= dt;
        position[i] += dlt;
    }
}

void SOrbitalSystem::integrate_forces(double dt) {
    for(int i=0; i<numAtoms; i++) {
        velocity[i] += 0.5 * (force[i]/mass) * dt;
    }
}

void SOrbitalSystem::integrate_Langevin(double dt) {
    std::normal_distribution<double> rd;    // default mean = 0, var = 1
    
    double alpha2 = exp(-gamma * dt);
    double sigma2 = sqrt((1 - pow(alpha2, 2)) * kT / mass);
    
    for (int i=0; i<numAtoms; i++) {
        velocity[i] = alpha2 * velocity[i] + sigma2 * vec3(rd(rng), rd(rng), rd(rng));
    }
}

void SOrbitalSystem::GA_solver(int &iter, double &err1, double &err2) {
    
    find_all_pairs();
    potential->compute_neighbor_list(pairs);
    
    iteration_GA(iter, err1, err2);
    
}

void SOrbitalSystem::step_NVT(double dt, int & iter, double & err1, double & err2) {
    
    // velocity verlet integration with Langevin damping
    move_atoms(dt);

    if (boundary_type == 1) wrap_PBC_position();
    
    find_all_pairs();
    compute_atom_distance();
    potential->compute_neighbor_list(pairs);

    if (solu_tag == 0) solve_bare_TB();
    if (solu_tag == 1) iteration_GA(iter, err1, err2);
    if (solu_tag == 2) rho.zeros();
    
    compute_forces();
    
    integrate_forces(dt);
    
    if (langevin_dym == 1)
        integrate_Langevin(dt);
}

// MUST run only after "e_elec" is called
double EPS_fe = 1.e-10;
double SOrbitalSystem::compute_elec_fe(void) {
    
    double sum1 = 0;
    for(int m=0; m<dimH; m++) {
        sum1 += -(fd_factor(m) * log(fabs(fd_factor(m) + EPS_fe)) + (1. - fd_factor(m)) * log(fabs(1. - fd_factor(m)) + EPS_fe));
    }
    
    s_qp = 2. * sum1 / ((double) numAtoms);
    
    double sum2 = 0;
    for(int i=0; i<numAtoms; i++) {
        double n0 = real(conj(atom(i).f(1)) * atom(i).f(1) + conj(atom(i).f(2)) * atom(i).f(2));

        double P[3] = {pow(1. - n0, 2), fabs(n0 * (1. - n0)), pow(n0, 2)};
        
        double sm = 0;
        for(int k=0; k<3; k++) {
            int deg = (k == 1) ? 2 : 1;
            
            double tmp = real(conj(atom(i).f(k)) * atom(i).f(k));
            
            sm += deg * tmp * log(tmp / (P[k] + EPS_fe));
        }
    
        sum2 += -sm;
    }
    
    s_sb = sum2 / ((double) numAtoms);
    
    fe = ee - kT_elec * (s_qp + s_sb);
    
    return fe;
}

double SOrbitalSystem::e_elec(void) {
    double sum1 = 0;
    double sum2 = 0;
    
    auto H_elec = Hamiltonian;
    auto rho_elec = rho;
    
    auto rr = H_elec * rho_elec;
    
    sum1 = trace(rr).real();
    
    for(int m=0; m<numAtoms; m++) sum2 += norm(atom(m).f(2));
    
    double sum3 = 0;
    for(int i=0; i<numAtoms; i++) sum3 += potential->atom[i].lambda * real(rho_elec(i, i));
    
    ee = (potential->numSpins() * (sum1 - sum3) + onsite_U * sum2) / ((double) numAtoms);
    
    return ee;
}

double SOrbitalSystem::e_kin(void) {
    double sum = 0;
    for (int i=0; i<numAtoms; i++) sum += velocity[i].norm2();
    ek = (0.5 * mass * sum) / ((double) numAtoms);
    return ek;
}

double SOrbitalSystem::e_pair(void) {
    ep = potential->pair_energy(pairs) / ((double) numAtoms);
    return ep;
}

double SOrbitalSystem::compute_avg_renormalization(void) {
    double sum = 0;
    for(int i=0; i<numAtoms; i++) {
        sum += norm(atom(i).Rq);
    }
    renormalization = sum / ((double) numAtoms);
    
    return renormalization;
}

double SOrbitalSystem::compute_avg_db_occupancy(void) {
    double sum = 0;
    for(int i=0; i<numAtoms; i++) {
        sum += norm(atom(i).f(2));  //in arma norm(X,P), if P is not specified the default values is 2, which gives norm^2
    }
    db_occupancy = sum / ((double) numAtoms);
    return db_occupancy;
}

double SOrbitalSystem::compute_avg_force(void) {
    double sum = 0;
    for(int i=0; i<numAtoms; i++) {
        sum += force[i].norm();
    }
    force_avg = sum / ((double) numAtoms);
    return force_avg;
}

void SOrbitalSystem::compute_standard_dev_density(void){
    double sum1 = 0;
    double sum1b = 0;
    double sum2 = 0;
    for(int i=0; i<numAtoms; i++) {
        double tmp = real(conj(atom(i).f(1)) * atom(i).f(1) + conj(atom(i).f(2)) * atom(i).f(2));
        sum1 += tmp;
        sum2 += pow(tmp, 2);
        
        sum1b += real(rho(i, i));
        
    }
    avg_n = sum1 / ((double) numAtoms);
    sgm_n = sqrt((sum2/(double) numAtoms) - pow(avg_n, 2));
    
    avg_rho_d = sum1b / ((double) numAtoms);
}

void SOrbitalSystem::compute_cm_velocity(void) {
    cm_velocity = {0, 0, 0};
    for(int i=0; i<numAtoms; i++) cm_velocity += velocity[i];
    cm_velocity /= ((double) numAtoms);
}

void SOrbitalSystem::subtract_cm_velocity(void) {
    for(int i=0; i<numAtoms; i++) velocity[i] -= cm_velocity;
}

void SOrbitalSystem::compute_pressure_tensor(void) {
    
    P_tensor.zeros();
    
    for(int i=0; i<numAtoms; i++) {
        for(int x=0; x<3; x++)
            for(int y=0; y<3; y++) P_tensor(x, y) += mass * velocity[i](x) * velocity[i](y);
    }
    
    compute_forces();
    
    for (auto & p : pairs) {
        if (p.delta.norm2() < pow(potential->rcut(),2)) {
            int i = p.index1;
            int j = p.index2;
            
            for(int x=0; x<3; x++)
                for(int y=0; y<3; y++) P_tensor(x, y) += 0.5 * p.delta(x) * (force[j](y) - force[i](y));
            
        }
    }
    
    P_tensor /= volume;
}

double SOrbitalSystem::compute_virial(void) {
    
    double EPS_V = 1.e-3;
    Vec<AtomPair> pairs_x = pairs;
    
    for (unsigned i=0; i<pairs.size(); i++)
        pairs_x[i].delta = (1. + EPS_V) * pairs[i].delta;
    
    Hamiltonian = potential->build_Hamiltonian(numAtoms, pairs_x);
    compute_density_matrix();
    
    double U1 = e_elec() * numAtoms + potential->pair_energy(pairs_x);
    
    for (unsigned i=0; i<pairs.size(); i++)
        pairs_x[i].delta = (1. - EPS_V) * pairs[i].delta;
    
    Hamiltonian = potential->build_Hamiltonian(numAtoms, pairs_x);
    compute_density_matrix();
    
    double U2 = e_elec() * numAtoms + potential->pair_energy(pairs_x);
    
    virial = -(U1 - U2) / (pow(1.+EPS_V, 3) - pow(1.-EPS_V, 3));
    
    return virial;
}

double SOrbitalSystem::compute_pressure(void){
    
    //        compute_virial();
    //        pressure = (virial + (2./3.) * ek * numAtoms) / volume;     // 3-dimensions
    
    pressure = (virial/3. + (2./3.) * ek * numAtoms) / volume;     // 3-dimensions
    
    return pressure;
}

void SOrbitalSystem::compute_avg_pair_distance(void){
    
    dist_pair = 0;
    D2_pair = 0;
    int n_tot = 0;
    
    for (auto const& p : pairs) {
        
        double dl2 = p.delta.norm2();
        
        if (dl2 < pow(potential->rcut(), 2)) {
            
            dist_pair += sqrt(dl2);
            D2_pair += dl2;
            n_tot++;
        }
    }
    
    dist_pair /= ((double) n_tot);
    D2_pair /= ((double) n_tot);
}

void SOrbitalSystem::reset_time(void) {
    time_diagonalization = 0;
    time_gw_iteration = 0;
}


void SOrbitalSystem::save_configuration(string const filename, string const filename_1, double db_condition) {
    
    std::ofstream fs;
    std::ofstream fs_1;
    
    fs.open(filename.c_str(), ios::out);
    fs.precision(12);

    fs_1.open(filename_1.c_str(), ios::out);
    fs_1.precision(12);

    for (int i=0; i<numAtoms; i++){
        if (db_ga(i) >= db_condition){
            fs << position[i].x << '\t' << position[i].y << '\t' << position[i].z << '\t';
            fs << velocity[i].x << '\t' << velocity[i].y << '\t' << velocity[i].z << endl;
        }
        else {
            fs_1 << position[i].x << '\t' << position[i].y << '\t' << position[i].z << '\t';
            fs_1 << velocity[i].x << '\t' << velocity[i].y << '\t' << velocity[i].z << endl;
        }
    }

    fs.close();
    fs_1.close();
}

void SOrbitalSystem::save_configuration(string const filename) {
    std::ofstream fs;
    
    fs.open(filename.c_str(), ios::out);
    fs.precision(12);

    for (int i=0; i<numAtoms; i++){
        fs << position[i].x << '\t' << position[i].y << '\t' << position[i].z << '\t';
        fs << velocity[i].x << '\t' << velocity[i].y << '\t' << velocity[i].z << endl;
    }

    fs.close();
}

void SOrbitalSystem::read_configuration(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi, const string filename) {

    position.resize(numAtoms);

    std::ifstream fs;
    fs.open(filename.c_str(), ios::in);
    for (int i=0; i<numAtoms; i++) {
        fs >> position[i].x >> position[i].y >> position[i].z >> velocity[i].x >> velocity[i].y >> velocity[i].z;
    }
    
    bdsLo = {0, 0, 0};
    bdsHi = {Lt, Lt, Lz};
    
    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box(0) * system_box(1) * system_box(2);
    
    for (int i=0; i<numAtoms; i++) {
        force[i] = vec3 {0, 0, 0};
    }
    
    
    pos_init = Vec<vec3>(numAtoms);
    for (int i=0; i<numAtoms; i++) pos_init[i] = position[i];

}

void SOrbitalSystem::plot_binding_curve(const string filename) {
    
    std::ofstream fs;
    
    fs.open(filename.c_str(), ios::out);
    fs.precision(12);
    
    double dr = 0.01;
    for(double r = 1.e-5; r<=potential->barePotential.rcut(); r+= dr) {
        double e_t = potential->barePotential.hopping(r)(0,0);
        double e_p = potential->barePotential.phi(r);
        fs << r << '\t' << e_t << '\t' << e_p << '\t' << 2 * e_t + e_p << endl;
    }
    
    fs.close();
}

void SOrbitalSystem::print_spectrum(const string filename) {
    
    std::ofstream fs;
    
    fs.open(filename.c_str(), ios::out);
    fs.precision(4);

    for (int m=0; m<dimH; m++){
        fs << m << '\t' 
        << setw(9) << eigE(m) << '\t' 
        << setw(9) << fd_factor(m) << '\t' 
        << setw(9) << 1 / kT_elec * exp((eigE(m) - mu) / kT_elec) / norm(exp((eigE(m) - mu) / kT_elec) + 1) 
        << endl;
    }
    fs.close();
}

void SOrbitalSystem::print_ga_data(const string filename) {
    
    std::ofstream fs;
    
    fs.open(filename.c_str(), ios::out);
    fs.precision(6); //12
    
    for (int i=0; i<numAtoms; i++){
        fs << rho(i, i).real() << '\t' << nf_ga(i) << '\t' << db_ga(i) << '\t' << norm(atom(i).Rq) << '\t';
        
        /* fs << rho(i, i).real() << '\t' << conj(atom(i).f(2)) * atom(i).f(1) + conj(atom(i).f(1)) * atom(i).f(0) << '\t' << sqrt(fabs(rho(i, i).real() * (1. - rho(i, i).real()))) << '\t' 
        << norm( ( conj(atom(i).f(2)) * atom(i).f(1) + conj(atom(i).f(1)) * atom(i).f(0) ) / sqrt(fabs(rho(i, i).real() * (1. - rho(i, i).real())) + 1.e-7) ) << '\t' << norm(atom(i).Rq) << '\t';
        fs << endl; */
    }
    
    fs.close();
}

double EPS_mm = 1.e-8;
double SOrbitalSystem::compute_metallic_fraction(void) {
    
    int sum = 0;
    for (int i=0; i<numAtoms; i++) {
        if (norm(atom(i).Rq) > EPS_mm) sum++;
    }
    
    return ((double) sum) / ((double) numAtoms);
}