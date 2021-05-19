#ifndef MESH_H
#define MESH_H

#include <unordered_map>
#include <array>
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <map>
#include <fstream>
#include <cstddef>

#include "Node.h"
#include "Constant.h"
#include "Cell.h"


struct my_hash{
  size_t operator()(const std::pair<std::size_t,std::size_t> & coppia)const{
    return (coppia.first*constant::x_nodes+coppia.second);
  }
};


class Mesh{

  static constexpr bool debug = false;

  using mesh_nodes =std::array<std::shared_ptr<Node>,constant::nodes_number>;
  typedef std::array<std::shared_ptr<Node>,4> cell_nodes;
  typedef std::array<double,constant::nodes_number> rho_out;
  typedef std::pair<std::size_t,std::size_t> indexes;

private:

  std::unordered_map<indexes,std::shared_ptr<Cell>,my_hash> mesh_map;
  mesh_nodes nodes;
  indexes get_upper_left_indexes_pair(const std::size_t & xx)const;
  indexes get_lower_left_indexes_pair(const std::size_t & xx)const;


public:
  Mesh();

  bool add_element(const cell_nodes & nodi,const indexes & p);
  // Prende dal main un array di 4 nodi, controlla che siano effettivamente 4,
  // costruisce la cella e verifica che la mesh abbia inserito la cella

  std::shared_ptr<Cell> find_my_cell(const indexes & p) const;
  // Ogni volta che particella si sposta e quando viene inizializzata troverà il suo indice in basso a sinistra,
  // questo serve per poi potergli assegnare la cella che contiene i nodi che la circondano

  void print_my_cell()const;

  void reset_velocities_nodes();


  void reset_mesh_charge();
  // Alla fine di ogni iterazione temporale resetto i valori della densità di carica

  rho_out get_rho_ions() const;
  // All'interno di solve ad un certo punto nella risoluzione del Poisson devo ricavare un vettore del valore della densità di carica in ogni nodo.

  void update_electric_field(const std::array<std::pair<double,double>,constant::nodes_number> &);
  // Metodo che dato in ingresso un vettore di coppie (Ex ed Ey) aggiorna il campo elettrico sui nodi

  std::ofstream & print_nodes_velocity(std::ofstream & os);
  //Metodo che stampa su file le velocità medie memorizzate nei nodi

};

#endif
