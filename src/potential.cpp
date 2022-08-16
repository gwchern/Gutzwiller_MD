//
//  potential.cpp
//  GA_solver-TBMD
//
//  Created by Gia-Wei Chern on 9/29/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include <cassert>
#include "potential.hpp"

SplicePolynomial::SplicePolynomial(double f0, double df0, double x0, double x1, double g1, double dg1) {
    this->x0 = x0;
    a0 = f0;
    a1 = df0;
    double dx = x1 - x0;
    double c1 = a0 + a1*dx - g1;
    double c2 = dx*dx;
    double c3 = dx*dx*dx;
    double d1 = a1 - dg1;
    double d2 = 2*dx;
    double d3 = 3*dx*dx;
    double det = c2*d3 - c3*d2;
    a2 = (c3*d1 - c1*d3) / det;
    a3 = (c1*d2 - c2*d1) / det;
}

double SplicePolynomial::eval(double x) {
    double dx = x - x0;
    return a0 + dx*(a1 + dx*(a2 + dx*a3));
}

double SplicePolynomial::deriv(double x) {
    double dx = x - x0;
    return a1 + dx*(2*a2 + dx*3*a3);
}


std::ostream& operator<< (std::ostream& os, AtomPair const& p) {
    os << "AtomPair { " << p.index1 << ", " << p.index2 << ", " << p.delta << " }";
    return os;
}

Vec<AtomPair> allPairs(Vec<vec3> const& pos, double cutoff) {
    Vec<AtomPair> pairs;
    for (int i = 0; i < pos.size(); i++) {
        for (int j = i+1; j < pos.size(); j++) {
            vec3 delta = pos[j]-pos[i];
            if (delta.norm2() < cutoff*cutoff)
                pairs.push_back(AtomPair {i, j, delta});
        }
    }
    return pairs;
}

double wrapDelta1d(double dx, double boxLength) {
    while (dx > boxLength/2)
        dx -= boxLength;
    while (dx < -boxLength/2)
        dx += boxLength;
    return dx;
}

vec3 wrapDelta(vec3 delta, vec3 boxLengths) {
    delta.x = wrapDelta1d(delta.x, boxLengths.x);
    delta.y = wrapDelta1d(delta.y, boxLengths.y);
    delta.z = wrapDelta1d(delta.z, boxLengths.z);
    return delta;
}

double wrapPosition1d(double x, double bdsLo, double bdsHi) {
    double length = bdsHi - bdsLo;
    while (x > bdsHi) {
        x -= length;
    }
    while (x < bdsLo) {
        x += length;
    }
    return x;
}

vec3 wrapPosition(vec3 p, vec3 bdsLo, vec3 bdsHi) {
    p.x = wrapPosition1d(p.x, bdsLo.x, bdsHi.x);
    p.y = wrapPosition1d(p.y, bdsLo.y, bdsHi.y);
    p.z = wrapPosition1d(p.z, bdsLo.z, bdsHi.z);
    return p;
}

Vec<AtomPair> allPairsPeriodic(Vec<vec3> const& pos, double cutoff, vec3 boxLengths) {
#ifdef WITH_TBB
    tbb::concurrent_vector<AtomPair> pairs;
    auto fn = [&](size_t i) {
        for (int j = i+1; j < pos.size(); j++) {
            vec3 delta = wrapDelta(pos[j]-pos[i], boxLengths);
            if (delta.norm2() < cutoff*cutoff)
                pairs.push_back(AtomPair {(int)i, j, delta});
        }
    };
    tbb::parallel_for(size_t(0), size_t(pos.size()), fn);
    return Vec<AtomPair>(pairs.begin(), pairs.end());
#else
    Vec<AtomPair> pairs;
    for (int i = 0; i < pos.size(); i++) {
        for (int j = i+1; j < pos.size(); j++) {
            vec3 delta = wrapDelta(pos[j]-pos[i], boxLengths);
            if (delta.norm2() < cutoff*cutoff)
                pairs.push_back(AtomPair {i, j, delta});
        }
    }
    return pairs;
#endif
}

arma::cx_mat Potential::build_Hamiltonian(int numAtoms, Vec<AtomPair> const& pairs) {
    
    int numOrbs = numOrbitalsPerSite();
    arma::mat h(numOrbs, numOrbs);
    arma::mat tmp(numOrbs, numOrbs);
    
    int dim = numAtoms*numOrbs;
    
    arma::cx_mat H(dim, dim);
    H.zeros();
    
    // same site orbital overlap
    for (int i = 0; i < numAtoms; i++) {
        fill_TB_hoppings({i, i, vec3{0,0,0}}, h, tmp, tmp, tmp);
        for (int o1 = 0; o1 < numOrbs; o1++) {
            for (int o2 = 0; o2 < numOrbs; o2++) {
                cx_flt v = h(o1, o2);
                // ensure diagonal is non-zero
                if (o1 == o2 && std::abs(v) == 0.0) {
                    v = std::numeric_limits<double>::epsilon();
                }
                H(i*numOrbs+o1, i*numOrbs+o2) = v;
            }
        }
    }
    
    // neighbor overlap
    for (auto const& p : pairs) {
        if (p.delta.norm2() < rcut()*rcut()) {
            int i = p.index1;
            int j = p.index2;
            assert(i < j);
            fill_TB_hoppings(p, h, tmp, tmp, tmp);
            for (int o1 = 0; o1 < numOrbs; o1++) {
                for (int o2 = 0; o2 < numOrbs; o2++) {
                    
                    cx_flt v = h(o1, o2);
                    
                    H(i*numOrbs+o1, j*numOrbs+o2) = v;
                    H(j*numOrbs+o2, i*numOrbs+o1) = v;
                }
            }
        }
    }
    
    return H;
}

double Potential::pair_energy(Vec<AtomPair> const& pairs) {
    double E_pair = 0.0;
    for (auto const& p : pairs) {
        if (p.delta.norm2() < rcut()*rcut()) {
            E_pair += phi(p.delta.norm());
        }
    }
    return E_pair;
}

void Potential::force(arma::cx_mat const& dE_dH, Vec<AtomPair> const& pairs, Vec<vec3>& forces, Vec<Vec<vec3>>& force_ij, double &virial) {
    std::fill(forces.begin(), forces.end(), vec3{0, 0, 0});

    std::fill(force_ij.begin(), force_ij.end(), forces);
    
    virial = 0;
    
    int nOrbs = numOrbitalsPerSite();
    arma::mat tmp(nOrbs, nOrbs);
    arma::mat dh_dx(nOrbs, nOrbs);
    arma::mat dh_dy(nOrbs, nOrbs);
    arma::mat dh_dz(nOrbs, nOrbs);
    
    for (auto const& p : pairs) {
        if (p.delta.norm2() < rcut()*rcut()) {
            int i = p.index1;
            int j = p.index2;
            vec3 f_ij {0, 0, 0};
            
            // electronic part
            fill_TB_hoppings(p, tmp, dh_dx, dh_dy, dh_dz);
            for (int o1 = 0; o1 < nOrbs; o1++) {
                for (int o2 = 0; o2 < nOrbs; o2++) {
                    cx_double dE_dH_ij = dE_dH(i*nOrbs + o1, j*nOrbs + o2);
                    cx_double dE_dH_ji = (i == j) ? 0.0 : dE_dH(j*nOrbs + o2, i*nOrbs + o1);
                    vec3 dh_dr { dh_dx(o1, o2), dh_dy(o1, o2), dh_dz(o1, o2) };
                    f_ij += - dh_dr * numSpins() * real(dE_dH_ij + dE_dH_ji);
                }
            }
            
            // classical part
            double r = p.delta.norm();
            f_ij += p.delta * (-dphi_dr(r) / r);   //f_ij is force on atom j due to i
            
            forces[j] += f_ij;  
            forces[i] -= f_ij;

            force_ij[i][j] = f_ij;
            force_ij[j][i] = -f_ij;

            virial += f_ij.dot(p.delta);
        }
    }

}

void Potential::force_elec(arma::cx_mat const& dE_dH, Vec<AtomPair> const& pairs, Vec<vec3>& forces, Vec<Vec<vec3>>& force_elec_ij, double &virial) {
    std::fill(forces.begin(), forces.end(), vec3{0, 0, 0});

    std::fill(force_elec_ij.begin(), force_elec_ij.end(), forces);
    
    virial = 0;
    
    int nOrbs = numOrbitalsPerSite();
    arma::mat tmp(nOrbs, nOrbs);
    arma::mat dh_dx(nOrbs, nOrbs);
    arma::mat dh_dy(nOrbs, nOrbs);
    arma::mat dh_dz(nOrbs, nOrbs);
    
    for (auto const& p : pairs) {
        if (p.delta.norm2() < rcut()*rcut()) {
            int i = p.index1;
            int j = p.index2;
            vec3 f_ij {0, 0, 0};
            
            // electronic part
            fill_TB_hoppings(p, tmp, dh_dx, dh_dy, dh_dz);
            for (int o1 = 0; o1 < nOrbs; o1++) {
                for (int o2 = 0; o2 < nOrbs; o2++) {
                    cx_double dE_dH_ij = dE_dH(i*nOrbs + o1, j*nOrbs + o2);
                    cx_double dE_dH_ji = (i == j) ? 0.0 : dE_dH(j*nOrbs + o2, i*nOrbs + o1);
                    vec3 dh_dr { dh_dx(o1, o2), dh_dy(o1, o2), dh_dz(o1, o2) };
                    f_ij += - dh_dr * numSpins() * real(dE_dH_ij + dE_dH_ji);
                }
            }
                        
            forces[j] += f_ij;
            forces[i] -= f_ij;

            force_elec_ij[i][j] = f_ij;
            force_elec_ij[j][i] = -f_ij;

            virial += f_ij.dot(p.delta);
        }
    }
}

//===================================================================================================
// codes for toy-model s-orbital
//===================================================================================================

double ToyModelSOribtal::phi(double r) {
    return (r > rmax) ? 0 : v0 * exp(-beta * (pow(r, 4) + r)); 
}

double ToyModelSOribtal::dphi_dr(double r) {
    return (r > rmax) ? 0 : -v0 * beta * (4 * pow(r, 3) + 1) * exp(-beta * (pow(r, 4) + r));
}

arma::mat ToyModelSOribtal::hopping(double r) {
    arma::mat h(1,1);
    double t = (r > rmax) ? 0 : -h0 * exp(-alpha * r);
    
    h(0, 0) = t;
    return h;
}

vec3 ToyModelSOribtal::dt_dr(vec3 delta) {
    double r = delta.norm();

    return (r > rmax) ? vec3(0, 0, 0) : h0 * alpha * exp(-alpha * r) * delta / r;
}

void ToyModelSOribtal::fill_TB_hoppings(AtomPair pair,
                                        arma::mat& h,
                                        arma::mat& dh_dx,
                                        arma::mat& dh_dy,
                                        arma::mat& dh_dz) {
    h.fill(0);
    dh_dx.fill(0);
    dh_dy.fill(0);
    dh_dz.fill(0);
    if (pair.index1 == pair.index2) {
        h(0, 0) = 0;
        return;
    }
    
    double r = pair.delta.norm();
    
    // ss hopping
    double t     = (r > rmax) ? 0 : -h0 * exp(-alpha * r);
    vec3   dt_dr = (r > rmax) ? vec3(0, 0, 0) : h0 * alpha * exp(-alpha * r) * pair.delta / r;
    
    h(0, 0) = t;
    dh_dx(0, 0) = dt_dr.x;
    dh_dy(0, 0) = dt_dr.y;
    dh_dz(0, 0) = dt_dr.z;
}





