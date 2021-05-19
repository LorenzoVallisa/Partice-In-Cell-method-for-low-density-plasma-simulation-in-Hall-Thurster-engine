#ifndef CELL_H
#define CELL_H

#include <memory>
#include <array>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>

#include "Constant.h"
#include "Node.h"

class Cell{

  static constexpr bool debug = false;

private:
  std::array<std::shared_ptr<Node>,4> cell_nodes;
  // Contiene i shared_ptr ai suoi 4 nodi

  double compute_surface(double x1, double x2,double y1, double y2)const;
  // Metodo privato che calcola calcola area relativa a posizione particella e nodo

public:
  Cell(std::array<std::shared_ptr<Node>,4> nodes):cell_nodes(nodes){}
  // Non c'è defualt constructor

  void print_my_nodes()const;
  // Metodo per controllare quali nodi ha la cella

  void update_charge(const double x, const double y);
  // Data posizione particella fa update carica (weighting function)

  std::pair<double,double> update_field(const double x, const double y)const;

  void update_velocity(double x,double y,double vx,double vy);
  // Per fini di stampa di velocità sui nodi, vado a calcolarmi la velocità sui 4 nodi
  // della cella dat la posizione della particella


};

#endif
