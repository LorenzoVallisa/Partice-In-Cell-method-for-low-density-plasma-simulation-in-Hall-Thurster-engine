#include "System.h"

void System::print_particles() const{
	for(auto & mp : mplist){
		std::cout<<"vx=";
		std::cout<<mp.get_vx();
		std::cout<<" vy=";
		std::cout<<mp.get_vy();
		std::cout<<std::endl;
		
		std::cout<<"x=";
		std::cout<<mp.indexes().first;
		std::cout<<" y=";
		std::cout<<mp.indexes().second;
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}


void System::assign_cells(){

  //Viene chiamato dopo aver inizializzato le particelle
  //per associarle alla cella che le contiene

  /////////////////////////////////////////
  if(debug)
    std::cout<<"Assigning Cells"<<std::endl;
    /////////////////////////////////////////


  for(auto & mp : mplist){
    const auto & ij = mp.indexes(); 	//restituisce una pair

    if(debug)
      std::cout<<" Found indexes "<<ij.first<<" - "<<ij.second<<std::endl;

    auto c = mesh.find_my_cell(ij);	//restituisce uno shared_ptr
    mp.set_cell(c);

    /////////////////////////////////////////
    if(debug){
      std::cout<<"To particle whose position is "<<mp.get_x()<<" - "<<mp.get_y()<<std::endl;
      std::cout<<" the cell whose x and y are     "<<std::endl;
      c->print_my_nodes();
      std::cout<<" have been assigned"<<std::endl;
    }
  /////////////////////////////////////////

}
}






void System::print_vel(std::ofstream& os){
  mesh.print_nodes_velocity(os);
}



void System::check_position(){

  /////////////////////////////////////////
  if(debug){
    std::cout<<"Particles AT THE MOMENT in the system "<<std::endl;
    for(const auto & it:mplist){
      std::cout<<it.get_x()<<"  -  "<<it.get_y()<<std::endl;
    }
  std::cout<<"Preparing to delete those OUT OF BOUNDARIES "<<std::endl;
  }
  /////////////////////////////////////////


  mplist.remove_if(Node_terminator());

  /////////////////////////////////////////
  if(debug){
    std::cout<<"Particles LEFT in the system "<<std::endl;
    for(const auto & it:mplist)
      std::cout<<it.get_x()<<"  -  "<<it.get_y()<<std::endl;
    }
  /////////////////////////////////////////



}




void System::reset_node_velocities(){

  /////////////////////////////////////////
  if(debug)
  std::cout<<"Resetting all nodes velocities from system"<<std::endl;
  /////////////////////////////////////////

  mesh.reset_velocities_nodes();
}




void System::update_charge(){

  /////////////////////////////////////////
  if(debug)
    std::cout<<"Updating charge in system"<<std::endl;
  /////////////////////////////////////////

  for(auto & mp : mplist)
    mp.get_mycell()->update_charge(mp.get_x(),mp.get_y());


}






	void System::update_particles_field(){

    /////////////////////////////////////////
    if(debug)
      std::cout<<"Updating Particle field from system"<<std::endl;
    /////////////////////////////////////////

  for(auto & mp : mplist){
    auto campo=mp.get_mycell() -> update_field(mp.get_x(),mp.get_y());
    mp.set_E( campo.first, campo.second );
  }

  /////////////////////////////////////////
  if(debug){
    std::cout<<"Electric Field of each particle "<<std::endl;
    for(auto & mp : mplist)
      std::cout<<"    ++   "<<mp.get_Ex()<<"   --    "<<mp.get_Ey()<<"   ++   "<<std::endl;
  }
  /////////////////////////////////////////

}



void System::set_nodes_field(std::array<std::pair<double,double>,constant::nodes_number> E){

  /////////////////////////////////////////
  if(debug)
    std::cout<<"Setting nodes field from system"<<std::endl;
  /////////////////////////////////////////

  mesh.update_electric_field(E);
}


void System::initial_velocity(const double& dt){

    /////////////////////////////////////////
    if(debug)
      std::cout<<"Setting initial velocities"<<std::endl;
      /////////////////////////////////////////

  for(auto & mp:mplist){
  	double vx=mp.get_vx();
  	double vy=mp.get_vy();
  	vx+=constant::mp_charge / constant::mp_mass * mp.get_Ex() * (-dt/2);
  	vy+=constant::mp_charge / constant::mp_mass * mp.get_Ey() * (-dt/2);
		mp.set_v(vx,vy);
  }
}



void System::motion_equation(Solve & solref, const double& dt){

  /////////////////////////////////////////
  if(debug)
    std::cout<<"Starting motion equation solver"<<std::endl;
  /////////////////////////////////////////

  for(auto & mp : mplist){
    auto tupla=solref.leapfrog(mp.get_vx(),mp.get_vy(),mp.get_Ex(),mp.get_Ey(),
                      mp.get_x(),mp.get_y(),dt,constant::Qe / constant::mp_mass);
    mp.set_pos(std::get<0>(tupla),std::get<1>(tupla));
    mp.set_v(std::get<2>(tupla),std::get<3>(tupla));
  }
}



void System::interp_vel(){

  for(auto & mp : mplist)
    mp.get_mycell()->update_velocity(mp.get_x(),mp.get_y(),mp.get_vx(),mp.get_vy());
}



void System::emplace_front(size_t n){

  //Metodo che aggiunge particelle al bordo sinistro chiamando
  //il costruttore di MacroParticle che prende un double in ingresso
  for(size_t i=0; i<n; ++i)
    mplist.emplace_front(constant::vdrift);
}




void System::reset_ion_charge(){ mesh.reset_mesh_charge();}
