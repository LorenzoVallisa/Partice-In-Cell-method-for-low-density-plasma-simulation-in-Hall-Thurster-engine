#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <utility>
#include <cstddef>

#include "Constant.h"



class Node{

  static constexpr bool debug = false;
private:
  const std::size_t i;
  const std::size_t j;
  double rho=0; // Densità di carica nel nodo della mesh
  double Ex=0;  // Valore campo elettrico in quel nodo per x
  double Ey=0;  // Valore campo elettrico in quel nodo per x
  double vx=0;  // Velocità media delle particelle nelle 4 celle circostanti lungo l'asse x
  double vy=0;  // Velocità media delle particelle nelle 4 celle circostanti lungo l'asse y



public:
  Node(std::size_t ii, std::size_t jj):i(ii),j(jj){};
  // Niente Default perchè la posizione del nodo DEVE essere data da subito e non può cambiare

  Node(std::pair<std::size_t,std::size_t> ind_pair):i(ind_pair.first),j(ind_pair.second) {};

  void print_me()const;

  void reset_charge();
  // Alla fine di ogni ciclo temporale ho bisogno di azzerare la carica sul nodo

  double get_x()const{
    return constant::lambda*static_cast<double>(i);
  }

  double get_y()const{
    return constant::lambda*static_cast<double>(j);
  }
  // Calcola la posizione di x e y usando lambda(lunghezza singola mesh quadratica)

  std::size_t get_i()const {
    return i;
  }

  std::size_t get_j()const{
    return j;
  }


  double get_density_charge()const{
    return rho;
  }

  double get_Ex()const{
    return Ex;
  }
  double get_Ey()const{
    return Ey;
  }


  double get_vx()const{
    return vx;
  }
  double get_vy()const{
    return vy;
  }

  void set_vx( double vxx){
    vx+= vxx;
  }

  void set_vy( double vyy){
    vy+= vyy;
  }

  void reset_velocities();

  void set_density_charge(double new_rho){
    rho+=new_rho;
  }

  // Metodo usato per aggioranre densità di carica del nodo:
  // se piu volte viene chiamato in un ciclo (piu particelle lo hanno in comune)
  // allora sommo i valori

  void set_electric_field(double Exx, double Eyy);

  std::pair<std::size_t,std::size_t> next_first()const{
    return std::make_pair((i+1)%constant::x_nodes,j);
  }
  std::pair<std::size_t,std::size_t> next_second()const{
    return std::make_pair((i)%constant::x_nodes,j+1);
  }
  std::pair<std::size_t,std::size_t> next_third()const{
    return std::make_pair((i+1)%constant::x_nodes,j+1);
  }

  // Metodi usati per trovare i 4 nodi limitrofi alla cella


  bool operator<(const Node & l)const;

};

#endif
