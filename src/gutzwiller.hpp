//
//  gutzwiller.hpp
//  GQMD-Hubbard-liquid
//
//  Created by Gia-Wei Chern on 1/11/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#ifndef gutzwiller_hpp
#define gutzwiller_hpp

#include <iostream>
#include <fstream>

#include "util.hpp"
#include "potential.hpp"


class GA_S_Orbital : public Potential {
public:
    
    static constexpr int numOrbitals = 1;
    static constexpr int numSOrbitalGAParam = 3;
    
    struct SD {
        arma::mat Mpp = arma::mat(numSOrbitalGAParam, numSOrbitalGAParam);
        arma::mat Npp = arma::mat(numSOrbitalGAParam, numSOrbitalGAParam);
        arma::mat Upp = arma::mat(numSOrbitalGAParam, numSOrbitalGAParam);
        SD() {
            
            Mpp(0, 1) = Mpp(1, 0) = Mpp(1, 2) = Mpp(2, 1) = sqrt(2.);
            Npp(1, 1) = 0.5;
            Npp(2, 2) = 1;
            Upp(2, 2) = 1;
        }
    };
    static SD sd;
    
    int numAtoms;
    double onsiteU;
    
    //double kT;
    
    //ModelSOribtal barePotential;
    ToyModelSOribtal barePotential;
    
    RNG rng;
    
    struct Neighbor {
        int idx;
        arma::mat hopping = arma::mat(numOrbitals, numOrbitals);
    };
    
    struct GA_Data {
    public:
        arma::cx_vec f = arma::cx_vec(numSOrbitalGAParam);
        //arma::cx_vec cp = arma::cx_vec(numSOrbitalGAParam);
        arma::cx_mat eigvec = arma::cx_mat(numSOrbitalGAParam, numSOrbitalGAParam);
        arma::vec eigval = arma::vec(numSOrbitalGAParam);
        cx_double Rq;   // hopping renormalization R
        cx_double Delta;
        
        arma::cx_vec f_tmp = arma::cx_vec(numSOrbitalGAParam);
        
        double lambda;     // lagrangin multiplier
        
        Vec<Neighbor> neighbors;
    };
    
    Vec<GA_Data> atom;
    
    // constructor:
    GA_S_Orbital(int natoms, double _U);
    
    // ====================================================
    // basic routines
    int numSpins() {return 2;};
    double rcut() {return 5.001;}   //5.001 //5.29
    int numOrbitalsPerSite() {return numOrbitals;}
    
    // ====================================================
    // initialization of variational parameters cp[]
    void init_var_param(arma::cx_mat const& rho);
    
    // ====================================================
    // setup neighbor list for Dm
    void compute_neighbor_list(Vec<AtomPair> const& pairs);
    
    // ====================================================
    // renormalization factor Rq:
    void compute_renormalizations(arma::cx_mat const& rho);
    void compute_Deltas(arma::cx_mat const& rho);

    void reset_GA_parameters(void);
    double adjust_lambda(arma::cx_mat const& rho, double r);

    // ====================================================
    // entropy contribution:
    
    arma::mat compute_entropy_term(arma::cx_mat const& rho, int i);
    
    // ====================================================
    // main GA root-finding routine:
    void compute_var_param(arma::cx_mat const& rho1, double r, double T_e);
    
    // ====================================================
    // output GA data
    void print_GA_data(std::ofstream &fs);
    void print_neighbor_list(std::ofstream &fs);
    
    // ====================================================
    // compute renormalized TB matrix:
    arma::mat hoppingMatrix(double r);
    double phi(double r);
    double dphi_dr(double r);
    void fill_TB_hoppings(AtomPair pair,
                          arma::mat& h,
                          arma::mat& dh_dx,
                          arma::mat& dh_dy,
                          arma::mat& dh_dz);
};


#endif /* gutzwiller_hpp */
