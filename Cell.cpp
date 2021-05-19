#include "Cell.h"

void Cell::print_my_nodes()const{
  std::cout<<" This cell has the following nodes:  "<<std::endl;
  for(const auto & i:cell_nodes){
    i->print_me();
  }

}


void Cell::update_charge(const double x, const double y){

  // In questo metodo devo calcolare il peso che ogni nodo ha in base alla posizione della particella
  // contenuta nella cella. Come carica da spartire uso l'approssimazione di Carica delle Macroparticelle.

  /////////////////////////////////////////
    if(debug){
    std::cout<< " Updating nodes contained in the Cell whose indexes are (LLI) "<<std::endl;
    cell_nodes[0]->print_me();
    }
  /////////////////////////////////////////

  for(auto & it : cell_nodes){

    double weight = compute_surface(it->get_x(), x,it->get_y(), y);

    /////////////////////////////////////////
    if(debug){
    std::cout<<" Adding density charge of node "<<std::endl;
    it->print_me();
    std::cout<<" up by "<<constant::mpq*weight<<std::endl;
    std::cout<<" With weight "<<weight<<std::endl;
    }
    /////////////////////////////////////////


    it->set_density_charge(constant::mpq*weight*constant::spwt/(constant::lambda*constant::lambda));

  }
}



double Cell::compute_surface(double x1, double x2,double y1, double y2)const{

  // Metodo ausiliario usato per costruire la weighting function

  return (1-std::fabs(x1-x2)/constant::lambda)*(1-std::fabs(y1-y2)/constant::lambda);
}


void Cell::update_velocity(double x,double y,double vx,double vy){

  for(auto & it : cell_nodes){
    double weight = compute_surface(it->get_x(), x,it->get_y(), y);

    it->set_vx(vx*weight);
    it->set_vy(vy*weight);
    //I nodi in base al loro peso si ripartiscono il valore della velocità sulla particella

}

}




std::pair<double,double> Cell::update_field(const double x, const double y)const{

  // In questo metodo devo aggiornare il campo eletrico dai nodi alle particelle.
  // Questo metodo fa parte di un a chiamata di un metodo del sistema, che restituisce
  // campi elettrici (x ed y) di una singola particella che si trova nella cella chiamata

  std::pair<double,double> out;
  std::vector<double> ex_vect;
  std::vector<double> ey_vect;

  /////////////////////////////////////////
  if(debug){
  std::cout<< " Updating particle electric filed contained in the Cell whose indexes are (LLI) "<<std::endl;
  cell_nodes[0]->print_me();
  }
  /////////////////////////////////////////


  for(const auto & it : cell_nodes){
    double weight = compute_surface(it->get_x(), x,it->get_y(), y);
    ex_vect.push_back(weight*(it->get_Ex()));
    ey_vect.push_back(weight*(it->get_Ey()));
  }
  out.first = std::accumulate(ex_vect.begin(),ex_vect.end(),0);
  out.second = std::accumulate(ey_vect.begin(),ey_vect.end(),0);

  /////////////////////////////////////////
  if(debug){
  std::cout<<" Total contribution of x Electric Field "<<out.first<<std::endl;
  std::cout<<" Total contribution of y Electric Field "<<out.second<<std::endl;
  }
  /////////////////////////////////////////

  return out;
}
