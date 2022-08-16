//
//  gutzwiller.cpp
//  GQMD-Hubbard-liquid
//
//  Created by Gia-Wei Chern on 1/11/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include "util.hpp"
#include "potential.hpp"
#include "gutzwiller.hpp"


GA_S_Orbital::SD GA_S_Orbital::sd;

GA_S_Orbital::GA_S_Orbital(int natoms, double _U) {
    numAtoms = natoms;
    atom = Vec<GA_Data>(natoms);
    for(int i=0; i<natoms; i++) atom[i].Rq = 1;
    
    onsiteU = _U;
    
    std::random_device seed;
    rng = RNG(seed());
}

// ===== tight-binding and pair-potential model for s-orbital

arma::mat GA_S_Orbital::hoppingMatrix(double r) {
    return barePotential.hopping(r);
}

double GA_S_Orbital::phi(double r) {
    return barePotential.phi(r);
    //return (r > rcut()) ? 0 : ; 
}

double GA_S_Orbital::dphi_dr(double r) {
    return barePotential.dphi_dr(r);
    //return (r > rcut()) ? 0 : ;
}

void GA_S_Orbital::fill_TB_hoppings(AtomPair pair,
                                    arma::mat& h,
                                    arma::mat& dh_dx,
                                    arma::mat& dh_dy,
                                    arma::mat& dh_dz) {
    
    barePotential.fill_TB_hoppings(pair, h, dh_dx, dh_dy, dh_dz);
    
    if (pair.index1 != pair.index2){

        // !!! Important !!! assume "R" are always real-valued:  (IF NOT, need to modify the definition of fill_TB_hoppings)
        
        double renormFactor = real( atom[pair.index1].Rq * atom[pair.index2].Rq );
        
        h *= renormFactor;
        dh_dx *= renormFactor;
        dh_dy *= renormFactor;
        dh_dz *= renormFactor;
    }
    else {    // on-site term:
        h(0,0) = atom[pair.index1].lambda;
    }
}

// ===========================================================================
// Gutzwiller optimization:

void GA_S_Orbital::init_var_param(arma::cx_mat const& rho) {
    for (int i=0; i<numAtoms; i++){
        double n0 = real(rho(i, i));
        atom[i].f(0) = 1. - n0;
        atom[i].f(1) = sqrt(n0 * (1. - n0));
        atom[i].f(2) = n0;
        
        //atom[i].f(0) = 0.5;
        //atom[i].f(1) = 1./sqrt(2.);
        //atom[i].f(2) = 0.5;

        //atom[i].f(0) = sqrt(1. - 2 * n0);
        //atom[i].f(1) = sqrt(n0);
        //atom[i].f(2) = 0;
        
    }
}

void GA_S_Orbital::compute_neighbor_list(const Vec<AtomPair> &pairs) {
    
    int numOrbs = numOrbitalsPerSite();
    arma::mat h(numOrbs, numOrbs);
    
    for(int i=0; i<numAtoms; i++) atom[i].neighbors.clear();
    
    for (auto const& p : pairs) {
        
        if (p.delta.norm2() < rcut()*rcut()) {
            int m = p.index1;
            int n = p.index2;
            assert(m < n);
            
            h = hoppingMatrix(p.delta.norm());
            Neighbor nb;
            nb.hopping = h;
            nb.idx = n;
            atom[m].neighbors.push_back(nb);
            
            nb.idx = m;
            atom[n].neighbors.push_back(nb);
        }
    }
}

void GA_S_Orbital::reset_GA_parameters(void) {
    for(int i=0; i<numAtoms; i++) {
        atom[i].Rq = 1;
        //atom[i].lambda = onsiteU/2.;
        atom[i].lambda = 0;
    }
}

double EPS_Rq = 1.e-7;

void GA_S_Orbital::compute_renormalizations(arma::cx_mat const& rho) {
    for (int i=0; i<numAtoms; i++) {

        double n0 = rho(i, i).real();
        //double n0 = real(conj(atom[i].f(1)) * atom[i].f(1) + conj(atom[i].f(2)) * atom[i].f(2));
        
        atom[i].Rq = (conj(atom[i].f(2)) * atom[i].f(1) + conj(atom[i].f(1)) * atom[i].f(0)) / sqrt(fabs(n0 * (1. - n0)) + EPS_Rq);
        
        if (norm(atom[i].Rq) > 1) atom[i].Rq = 1.;  //prevent |Rq| > 1
    }

}

void GA_S_Orbital::compute_Deltas(arma::cx_mat const& rho) {
    
    for(int i=0; i<numAtoms; i++) {
        atom[i].Delta = 0;
        
        auto num_neighbors = atom[i].neighbors.size();
        for(int r=0; r<num_neighbors; r++) {
            int j = atom[i].neighbors[r].idx;
            
            atom[i].Delta += atom[i].neighbors[r].hopping(0,0) * conj(atom[j].Rq) * rho(j, i);
        }
    }
    
}

double EPS_nf = 1.e-10;
arma::mat GA_S_Orbital::compute_entropy_term(arma::cx_mat const& rho, int i) {
    arma::mat M(3, 3);
    M.zeros();
    //double n0 = rho(i, i).real();
    double n0 = real(conj(atom[i].f(1)) * atom[i].f(1) + conj(atom[i].f(2)) * atom[i].f(2));

    double P[3] = {pow(1. - n0, 2), fabs(n0 * (1. - n0)), pow(n0, 2)};
    for(int k=0; k<3; k++) {
        M(k, k) = 1. + log(real(conj(atom[i].f(k)) * atom[i].f(k)) / (P[k] + EPS_nf));
    }
    return M;
}

void GA_S_Orbital::compute_var_param(arma::cx_mat const& rho, double mix, double T_e) {
    
    //compute_renormalizations(Dm, Phi);
    
    compute_Deltas(rho);
    
    for(int i=0; i<numAtoms; i++) {
        arma::cx_mat Hpp(3, 3);
        double n0 = rho(i, i).real();
        double nf = real(conj(atom[i].f(1)) * atom[i].f(1) + conj(atom[i].f(2)) * atom[i].f(2));
                
        double Lmd = -4. * fabs(real(atom[i].Rq * atom[i].Delta)) * (nf - 0.5) / (fabs(nf * (1. - nf)) + EPS_nf);
        
        /*
        double Lmd_S = -2. * T_e * real(-conj(atom[i].f(0)) * atom[i].f(0) / ((1. - nf) + EPS_nf)
                                        + conj(atom[i].f(1)) * atom[i].f(1) * (1. - 2. * nf) / (fabs(nf * (1. - nf)) + EPS_nf)
                                        + conj(atom[i].f(2)) * atom[i].f(2) / (nf + EPS_nf));
        
        arma::mat Spp = compute_entropy_term(rho, i);
        
        Hpp = (atom[i].Delta / sqrt(fabs(n0 * (1. - n0))) + EPS_Rq) * sd.Mpp + onsiteU * sd.Upp + (Lmd - 2. * atom[i].lambda + Lmd_S) * sd.Npp + T_e * Spp;
        */
        
        
        Hpp = (atom[i].Delta / sqrt(fabs(n0 * (1. - n0))) + EPS_Rq) * sd.Mpp + onsiteU * sd.Upp + (Lmd - 2. * atom[i].lambda) * sd.Npp;

        /*
        cout << Spp << endl;
        cout << n0 << '\t' << nx << '\t' << atom[i].Delta << '\t' << atom[i].Rq << '\t' << atom[i].lambda << '\t' << Lmd << endl;
        cout << Hpp << endl;
        */
         
         
        arma::vec E(3);
        arma::cx_mat v(3, 3);
        
        arma::eig_sym(E, v, Hpp);
        
        int idx_min = 0;
        
        double e_min = E(0);
        for(int k=1; k<3; k++) {
            if(E(k) < e_min) {
                e_min = E(k);
                idx_min = k;
            }
        }
        
        for(int k=0; k<numSOrbitalGAParam; k++) atom[i].f_tmp(k) = atom[i].f(k);
        
        arma::cx_vec _f = arma::cx_vec(numSOrbitalGAParam);
        
        _f(0) = v(0, idx_min);
        _f(1) = v(1, idx_min) / sqrt(2.);
        _f(2) = v(2, idx_min);
        
        double _Q[numSOrbitalGAParam];
        for(int k=0; k<numSOrbitalGAParam; k++) {
            
            _Q[k] = mix * real(conj(atom[i].f_tmp(k)) * atom[i].f_tmp(k)) + (1. - mix) * real(conj(_f(k)) * _f(k));
            
            atom[i].f(k) = sqrt(_Q[k]);
        }


        /*
        atom[i].f(0) = v(0, idx_min);
        atom[i].f(1) = v(1, idx_min) / sqrt(2.);
        atom[i].f(2) = v(2, idx_min);
        */
    }
}

double GA_S_Orbital::adjust_lambda(arma::cx_mat const& rho, double r) {
    //std::ofstream fs("GA_iter");

    double sum = 0;
    for(int i=0; i<numAtoms; i++) {

        double n0 = rho(i, i).real();
        double npp = real(conj(atom[i].f(1)) * atom[i].f(1) + conj(atom[i].f(2)) * atom[i].f(2));

        //        int sgn = (n0 > npp) ? +1 : -1;
        //        atom[i].lambda += sgn * 0.1 * pow(fabs(n0 - npp), 0.3); // 0.01 * sgn * pow(fabs(n0 - npp), 0.8);
        
        atom[i].lambda += r * (n0 - npp);

        //fs << n0 << '\t' << npp << '\t' << endl;
        
        sum += fabs(n0 - npp);
    }
    
    return sum / ((double) numAtoms);
}

void GA_S_Orbital::print_GA_data(std::ofstream &fs) {
    for(int i=0; i<numAtoms; i++) {
        fs << "atom-" << i << "\t Rq = " << atom[i].Rq << "\t Delta = " << atom[i].Delta << '\t';
        fs << "cp = ";
        for(int a=0; a<numSOrbitalGAParam; a++) fs << atom[i].f(a) << '\t';
        fs << std::endl;
    }
}

void GA_S_Orbital::print_neighbor_list(std::ofstream &fs) {
    fs << "Neighbor list:\n";
    for(int i=0; i<numAtoms; i++) {
        fs << "atom-" << i << '\t';
        auto numNeighbors = atom[i].neighbors.size();
        for(int j=0; j<numNeighbors; j++) {
            fs << atom[i].neighbors[j].idx << '[' << atom[i].neighbors[j].hopping(0,0)/rydberg << ']' << '\t';
        }
        fs << std::endl;
    }
}

