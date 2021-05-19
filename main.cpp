#include "Solve.h"
#include "System.h"
#include <fstream>
#include <memory>
#include <algorithm>
#include <iostream>

using std::ofstream;

template <typename phi_out>
ofstream & print_phi(ofstream & os, const phi_out & phi){
	for(const auto & r: phi)
		os<<r<<"\n";
	return os;
}

int main(){

	//Inizializzazione
	Solve sol;
	System sys; 

	ofstream os_phi("Phi.csv");
	ofstream os_vel("NodesVelocity.csv");
	//ofstream os_rho_out("Rho_out.csv");
	std::array<double,constant::nodes_number> phi_start;
	phi_start.fill(0);
	
	//time loop
	const double dt=0.1*constant::lambda/constant::vdrift;
	for(size_t nt=0; nt<constant::maxit; ++nt){
	
		std::cout<<"Iteration number: "<<nt+1<<"/"<<constant::maxit<<std::endl;
		//sys.print_particles();
		
		//0. Assegna ad ogni nodo inizializzato la sua cella 
		
		sys.assign_cells();
		
		//1. Deposita carica delle particelle sui nodi della mesh
		sys.update_charge();
		
		//2. Solve the Poisson problem for the scalar potential phi
		auto rho_ions = sys.get_rho_ions(); //l'output è un vettore con la carica in ogni nodo
		
		auto phi = sol.poisson(rho_ions,phi_start);
		if(nt%(constant::maxit/constant::n_plots-1)==0) {
			print_phi(os_phi,phi); 
			//print_phi(os_rho_out,rho_ions); //unused: we plot only (vx,vy) and phi
		}

		//3. Compute E on each node, knowing phi
		auto E_nodes=sol.compute_electric_field(phi);
		sys.set_nodes_field(E_nodes);
		
		//4. Interpola E sulle particelle
		sys.update_particles_field();

		//5. Solve the Newton-Lorentz equation of motion
	
	  //Initial condition for the LeapFrog method. 
	  //To be commented if another integrator is used
	 	if(nt==0){	sys.initial_velocity(dt); }
	
		//Calculate and assign update v_i and x_i to the particles
		sys.motion_equation(sol,dt); 
	
		//6. Interpola v sui nodi e stampa il vettore
	 	if(nt%(constant::maxit/constant::n_plots)==0){
	 		//std::cout<<"nt="<<nt<<std::endl;
	 		sys.interp_vel(); // ora ogni nodo ha una media pesata delle velocità delle
	 											// particelle che si trovano nelle celle circostanti
	 		sys.print_vel(os_vel); //chiama un metodo di mesh che cicla su mesh_nodes e per ogni nodo stampa vx,vy su file
	 		sys.reset_node_velocities();
	 	}
	
		//7. Remove lost particles, add new ones
		sys.check_position();
		sys.emplace_front(constant::n_inflow);
	
		//8. Azzera la carica sui nodi e aggiorna phi_start
		sys.reset_ion_charge();
		//std::swap(phi,phi_start);
		phi_start=phi;
		

	}//end time loop
	
	return 0;
}
