#ifndef SOLVE_H
#define SOLVE_H

#include <utility>
#include <tuple>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <cstddef>

#include"Constant.h"
#include"Matrix.h"



class Solve{

	static constexpr bool debug = false;

private:
	double tol = 1e-2;
	std::size_t get_nodes_index(std::size_t i,std::size_t j)const;

	Matrix A; 	//Inizializza una matrice di zeri di dimensioni nodes_number x nodes_number, quindi una 160x160



	std::array<double,constant::nodes_number> convert_row_col(const std::array<double,constant::nodes_number> &)const;

	std::array<double,constant::nodes_number> convert_col_row(const std::array<double,constant::nodes_number> & input)const;

	double compute_residual(const std::array<double,constant::nodes_number> & bb,
													const std::array<double,constant::nodes_number> & xx)const;


public:

	//Costruttore: (crea la matrice A, default per il resto)

	Solve();

	//Metodo LeapFrog per la risoluzione delle equazioni di Newton-Lorentz, caso elettrostatico
	//Vantaggi: metodo esplicito, ordine 2, condizione di stabilità rispetto alla frequenza di plasma
	//nota dalla letteratura. Inoltre è un metodo simplettico e il sistema da risolvere è hamiltoniano, quindi
	//preserva l'energia del sistema dinamico --> utile se si volessero calolare quantità legate all'energia come output
	//Svantaggio principale: v e x devono essere integrate "sfalsate" in tempo di mezzo passo temporale,
	//è necessario quindi calcolare v(-dt/2) prima di iniziare il time-loop
	std::tuple<double,double,double,double>
	leapfrog(const double& vx0, const double& vy0, const double& Ex, const double& Ey,
			   const double& x0, const double& y0, const double& dt, const double& mult_fact);
	
	//Eulero avanti, stessa firma di leapfrog per facilitare lo scambio in system::motion_equation
	std::tuple<double,double,double,double>
	eulero_avanti(const double& vx0, const double& vy0, const double& Ex, const double& Ey,
			   const double& x0, const double& y0, const double& dt, const double& mult_fact);
	

	std::array<std::pair<double,double>,constant::nodes_number> compute_electric_field( 
																const std::array<double,constant::nodes_number> & phi);



	std::array<double,constant::nodes_number> poisson(const std::array<double,constant::nodes_number> & rho_ions,
																										const std::array<double,constant::nodes_number> & phi_start);

};

#endif
