#include "Mesh.h"



Mesh::indexes Mesh::get_lower_left_indexes_pair(const std::size_t & xx)const{

  // Metodo ausiliario usato in costruttore: data il numero della cella, trova gli indici
  // in basso a sinistra della cella in questione


  std::size_t j = xx/(constant::x_nodes-1);
  std::size_t i = xx%(constant::x_nodes-1);
  indexes ret(i,j);
  return ret;

}







Mesh::Mesh(){

  // Mesh costruisce tutti i nodi della griglia, facendo in modo da assegnare i puntatori alle celle corrispettive,
  // dopo avendone individuato l'indice di riferimento (quello in basso a sinistra).
  // Puntatori a Node perchè viene usato anche nell'array di Nodes all'interno di Mesh e perchè Cell diverse possono avere gli stessi Node.
  // Costruisco tutti i shared_ptr a nodi e li inserisco in una mappa: ad ogni shared_ptr a nodo
  // faccio corrispondere una coppia di indici, di modo da facilitare il find

  /////////////////////////////////////////
  if(debug)
    std::cout<<"Initializing Mesh"<<std::endl;
  /////////////////////////////////////////

  std::map<std::pair<int,int>,std::shared_ptr<Node>> inter_map;
  for(std::size_t i = 0;i<constant::x_nodes;++i){
    for(std::size_t j = 0;j<constant::y_nodes;++j){
      Node nd(i,j);
      std::shared_ptr<Node> nd_ptr=std::make_shared<Node>(nd);
      inter_map.insert(std::make_pair(std::make_pair(i,j),nd_ptr));
    }
  }

  //Ora che ho tutti i shared_ptr a nodi in una mappa, costruisco un ciclo sul numero di celle:
  //per ogni cella: 

  // 1 trovo la sua coppia di indici di riferimento
  // 2 cerco il puntatore al nodo di riferimeneto
  // 3 costruisco un puntatore allo stesso nodo
  // 4 faccio lo stesso con gli altri 4 nodi relativi ad una stessa cella (usando metodi implementati su Node)
  // 5 assemblo array
  // 6 chiamo add element

  std::size_t cell_number =(constant::x_nodes-1)*(constant::y_nodes-1);

  for(std::size_t ii = 0;ii<cell_number;++ii){
    indexes cell_pair = get_lower_left_indexes_pair(ii);

    auto p1 = inter_map.find(cell_pair);
    std::shared_ptr<Node> first_ptr(p1->second);
    auto p2 = inter_map.find(p1->second->next_first());
    std::shared_ptr<Node> second_ptr(p2->second);
    auto p3 = inter_map.find(p1->second->next_second());
    std::shared_ptr<Node> third_ptr(p3->second);
    auto p4 = inter_map.find(p1->second->next_third());
    std::shared_ptr<Node> fourth_ptr(p4->second);

    cell_nodes cn{first_ptr,second_ptr,third_ptr,fourth_ptr};

    /////////////////////////////////////////
    if (debug){
    std::cout<<" I insert INTO the UNORDERED_MAP a shared_ptr to a Cell with shared_pointer to the Nodes whose indexes are:  "<<std::endl;
    std::cout<<first_ptr->get_i()<<" - "<<first_ptr->get_j()<<std::endl;
    std::cout<<second_ptr->get_i()<<" - "<<second_ptr->get_j()<<std::endl;
    std::cout<<third_ptr->get_i()<<" - "<<third_ptr->get_j()<<std::endl;
    std::cout<<fourth_ptr->get_i()<<" - "<<fourth_ptr->get_j()<<std::endl;
    }
    /////////////////////////////////////////


    this->add_element(cn,get_lower_left_indexes_pair(ii));
}


std::size_t count = 0;
for (const auto & it : inter_map){
  // da inserire in ordine
  nodes[count]=it.second;
  count++;
}

/////////////////////////////////////////
if(debug){
std::cout<<" I have initialized an array with all following Nodes "<<std::endl;
/////////////////////////////////////////

for(const auto & nd : nodes)
  std::cout<<nd->get_i()<<" - "<<nd->get_j()<<std::endl;
}


}






void Mesh::reset_mesh_charge(){

  // Ho bisogno di questo metodo dal momento in cui
  // ogni volta che finisco un iterazione temporale le particelle si spostano
  // e la carica sui nodi va ricalcolata, ma NON è cumulativa quindi resetto i valori a zero

  for(const auto & it:nodes)
    it->reset_charge();

}





void Mesh::update_electric_field(const std::array<std::pair<double,double>,constant::nodes_number> & v_e){

  // In questo metodo si assume ovviamente che l'ordine dei valori, che viene inizialmente
  // estratto appunto da nodes, non vari e che quindi posso reinserire ciclando con un semplice
  // ciclo for

  /////////////////////////////////////////
  if(debug)
    std::cout<<" Updating Elecrtic Field from Mesh"<<std::endl;
  /////////////////////////////////////////

  for(std::size_t i=0;i<v_e.size();++i){
    nodes[i]->set_electric_field(v_e[i].first,v_e[i].second);
  }

}






void Mesh::print_my_cell()const{
  for(const auto & i:mesh_map){


    std::cout<<" Index referring to the following displayed cell :  "<<i.first.first<<" -  "<<i.first.second<<std::endl;
    (i.second)->print_my_nodes();
  }
  std::cout<<"Total "<<mesh_map.size()<<std::endl;
}










bool Mesh::add_element(const cell_nodes & nodi,const indexes & p){

  // nodes è array di shared_ptr a nodi: mesh lo prende come input insieme alla coppia di indici a cui fa riferimento,
  // come già dichiarato nell'analisi è quella relativa al nodo in alto a sinistra della cella

  if (nodi.size()==4){
  Cell new_cell(nodi);
  std::shared_ptr<Cell> ptr_cell = std::make_shared<Cell>(new_cell);

  // shared_ptr a Cell perchè abbiamo usato puntatori a Cell anche all'interno di Macroparticle

  auto ret = mesh_map.insert(std::make_pair(p,ptr_cell));
  return ret.second;
  }
  else
  return false;
}






std::shared_ptr<Cell> Mesh::find_my_cell(const indexes & p) const{

  return mesh_map.find(p)->second;
}





Mesh::rho_out Mesh::get_rho_ions() const{

  /////////////////////////////////////////
  if(debug)
    std::cout<<"Returning the density charge (ions) on whole nodes!"<<std::endl;
  /////////////////////////////////////////

  rho_out out;
  for(std::size_t i = 0; i<out.size(); ++i){
    out[i] = nodes[i]->get_density_charge();

    /////////////////////////////////////////
    if(debug)
      std::cout<<out[i]<<std::endl;
    ////////////////////////////////////////

    }
  return out;
}





void Mesh::reset_velocities_nodes(){
  for(auto & it:nodes)
    it->reset_velocities();

}






std::ofstream & Mesh::print_nodes_velocity(std::ofstream & os){
	for(const auto & n: nodes)
		os<<n->get_vx()<<","<<n->get_vy()<<"\n";
	return os;
}
