//
//  potential.hpp
//  GA_solver-TBMD
//
//  Created by Gia-Wei Chern on 9/29/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#ifndef potential_hpp
#define potential_hpp

#include <stdio.h>
#include <vector>
#include "util.hpp"


using namespace std;

inline int positive_mod(int i, int n) {
    return (i%n + n) % n;
}

// Polynomial that splices together functions f and g between points x0 and x1, of the form
//     a0 + a1 (x - x0) + a2 (x - x0)^2 + a3 (x - x0)^3
// Coefficients a0..a3 are selected so that the spliced function and its derivatives are continuous.
struct SplicePolynomial {
    double a0, a1, a2, a3, x0;
    SplicePolynomial(double f0, double df0, double x0, double x1, double g1, double dg1);
    double eval(double x);
    double deriv(double x);
};

double wrapDelta1d(double dx, double boxLength);
vec3 wrapDelta(vec3 delta, vec3 boxLengths);

double wrapPosition1d(double x, double bdsLo, double bdsHi);
vec3 wrapPosition(vec3 p, vec3 bdsLo, vec3 bdsHi);

struct AtomPair {
    int index1;
    int index2;
    vec3 delta;  // r_2 - r_1
};
std::ostream& operator<< (std::ostream& os, AtomPair const& p);

// Replace with fast neighbor finding algorithm
Vec<AtomPair> allPairs(Vec<vec3> const& pos, double cutoff);
Vec<AtomPair> allPairsPeriodic(Vec<vec3> const& pos, double cutoff, vec3 boxLengths);


class Potential {
protected:
    Vec<int> colors_to_groups(Vec<int> const& colors);
    
public:
    
    // Build Hamiltonian matrix
    arma::cx_mat build_Hamiltonian(int numAtoms, Vec<AtomPair> const& pairs);
    
    // Classical pair potential
    double pair_energy(Vec<AtomPair> const& pairs);
    
    // Electronic and classical force
    void force(arma::cx_mat const& dE_dH, Vec<AtomPair> const& pairs, Vec<vec3>& forces, Vec<Vec<vec3>>& force_ij, double& pressure);
    void force_elec(arma::cx_mat const& dE_dH, Vec<AtomPair> const& pairs, Vec<vec3>& forces, Vec<Vec<vec3>>& force_elec_ij, double& pressure);

    //void force(fkpm::SpMatBsr<cx_double> const& dE_dH, Vec<AtomPair> const& pairs, Vec<vec3>& forces, double& pressure, vec3 const& box);
    
    virtual int numSpins() = 0;
    virtual double rcut() = 0;
    virtual int numOrbitalsPerSite() = 0;
    
    virtual double phi(double r) = 0;
    virtual double dphi_dr(double r) = 0;
    virtual void fill_TB_hoppings(AtomPair pair,
                                  arma::mat& h,
                                  arma::mat& dh_dx,
                                  arma::mat& dh_dy,
                                  arma::mat& dh_dz) = 0;
};

class ToyModelSOribtal : public Potential {
public:

    // ss hopping
    static constexpr double h0 = 24;   // should be positive !!!
    static constexpr double alpha = 1.9;

    // pair potential v0 * exp(-beta * (pow(r, 4) + r)), equilibrium at ~1 Angstrom, bonding anti-bonding gap is 7.193 eV
    static constexpr double v0 = 100;
    static constexpr double beta = 2.2;

    
    
    static constexpr double rmax = 5.001;   //5.001 //12.0
    
    ToyModelSOribtal() {};
    int numSpins() {return 2;};
    double rcut() {return rmax;}
    int numOrbitalsPerSite() {return 1;}
    
    arma::mat hopping(double r);
    vec3 dt_dr(vec3 delta);
    double phi(double r);
    double dphi_dr(double r);
    void fill_TB_hoppings(AtomPair pair,
                          arma::mat& h,
                          arma::mat& dh_dx,
                          arma::mat& dh_dy,
                          arma::mat& dh_dz);
    
};

#endif /* potential_hpp */
