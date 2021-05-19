#ifndef SYSTEM_H
#define SYSTEM_H

#include <forward_list>
#include <list>
#include <memory>
#include <tuple>
#include <cstddef>
#include <iostream>

#include "MacroParticle.h"
#include "Mesh.h"
#include "Solve.h"


//Necessario per poter usare std::forward_list::remove_if
struct Node_terminator{
	bool operator()(const MacroParticle & mp){
		return (!(mp.get_y() < constant::y_domain && mp.get_x() < constant::x_domain && mp.get_y()>0 && mp.get_x()>0));

	}
};



class System{


	static constexpr bool debug = false;

private:
	std::forward_list<MacroParticle> mplist;
	Mesh mesh;

public:
	System():mplist(constant::initial_particles), mesh() {}
	
	void print_particles() const; //for debug

	//Seguono i metodi usati in main.cpp
	//Stesso ordine e stessa suddivisione in sezioni del main

	//0. Assegna ad ogni nodo inizializzato la sua cella 
	
	//(NOTA: ogni cella viene riconosciuta dal suo nodo in BASSO A SINISTRA,
	//quindi la particella cercherà quella coppia di inidici)

	void assign_cells();

	
	//1. Deposita carica delle particelle sui nodi della mesh

	//Ogni particella chiama un metodo della cella in cui si trova che effettua
	//la media pesata della carica e la ripartisce tra i suoi 4 vertici
	void update_charge();


	//2. Solve the Poisson problem for the scalar potential phi

	//Alias per un metodo di mesh che cicla sui nodi e restituisce un vettore
	//contenente il valore di rho -densità di carica- per ogni nodo
	std::array<double,constant::nodes_number> get_rho_ions(){

		/////////////////////////////////////////
		if(debug)
			std::cout<<" Calling get_rho_ions from system"<<std::endl;
		/////////////////////////////////////////

		return mesh.get_rho_ions();
	}

	//3. Compute E on each node, knowing phi

	//Alias per un metodo di mesh che cicla sui nodi e per ogni nodo
	//chiama nodo.set_electric_field(E)
	void set_nodes_field(std::array<std::pair<double,double>,constant::nodes_number> E);


	//4. Interpola E sulle particelle

	//Ogni particella mp chiama il metodo della cella in cui si trova "update_field"
	//che prende la posizione della particella e restituisce il valore del vettore
	//campo elettrico in quel punto, effettuando una media pesata sui valori del
	//campo nei vertici della cella.
	//I pesi della media sono gli stessi del metodo "cella.update_charge" del puto 1
	void update_particles_field();

	//5. Solve the Newton-Lorentz equation of motion v' = q * E / m, x'=v

	//Usando il metodo Leapfrog, è necessario conoscere la velocità
	//all'istante temporale t-dt/2. Per questo alla prima iterazione si
	//calcola v(-dt/2) usando il campo E(t) appena calcolato.
	void initial_velocity(const double& dt);


	//Chiama il metodo Leapfrog della classe Solve che risolve localmente l'equazione del moto
	void motion_equation(Solve & solref, const double& dt);


	//6. Interpola v sui nodi e stampa il vettore --> 1 metodo da scrivere in Cell.h

	//Come system.update_charge() del punto 1. Necessita di cell.update_velocity(),
	//analogo al cell.update_charge() già implementato ma che interpola vx,vy invece che la carica
	void interp_vel();

	//Alias per un metodo di Mesh che
	//scorre sui nodi e stampa in un file il valore di vx,vy per ogni nodo
	void print_vel(std::ofstream& os);

	void reset_node_velocities();


	//7. Remove lost particles, add new ones

	//Cicla sulle particelle ed elimina da System quelle che sono uscite dal dominio
	void check_position();

	//Aggiunge n nuove particelle a System. Le particelle vengono direttamente
	//construite usando il costruttore MacroParticle(double t), evitando quindi copie
	void emplace_front(std::size_t n);


	//8. Azzera la carica sui nodi
	void reset_ion_charge();
};

#endif
