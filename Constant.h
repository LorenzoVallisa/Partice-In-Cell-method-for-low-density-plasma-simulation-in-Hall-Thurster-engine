#ifndef CONSTANT_H
#define CONSTANT_H

#include <cmath>
#include <cstddef>

namespace constant {
const double initial_particles=0;
const double mpq = 1.0;
const double AMU = 1.661e-27;
const double Kb = 1.381e-23;			
const double vdrift=7e3;
const double spwt=3.6846e+05;
const std::size_t nodes_number = 160;
const std::size_t nodes=nodes_number;
const std::size_t x_nodes = 16;
const std::size_t y_nodes = 10;
const double mp_charge=mpq;
const double m_ion=32*AMU;
const double mp_mass=m_ion*spwt;
const double T=10.;
const std::size_t maxit=200;
const std::size_t n_inflow=135;
const std::size_t n_plots=4;
const double n0 = 1e12;
const double Te = 1.;
const double Ti = 0.1;
const double Qe = 1.602e-19;
const double Eps0= 8.854e-12;
const double lambda = std::sqrt(Eps0*Te/(n0*Qe));
const double Lx=15*lambda;
const double Ly=9*lambda;
const double surface_cell=lambda*lambda;
const double x_domain = lambda*(x_nodes-1);
const double y_domain = lambda*(y_nodes-1);
const double vth = std::sqrt(2*Qe*Ti/m_ion);
};

#endif
