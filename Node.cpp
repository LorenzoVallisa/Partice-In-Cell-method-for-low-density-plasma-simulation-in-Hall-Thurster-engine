#include "Node.h"


  void Node::print_me()const{
    std::cout<< " i----->      "<<i<<" j----->    "<<j<<std::endl;
    std::cout<< " x----->      "<<get_x()<<" y----->    "<<get_y()<<std::endl;
  }

  void Node::reset_charge(){
    /////////////////////////////////////////
    if(debug)
      std::cout<<" Resetting density charge of Node whose indexes are  "<<i<<" - "<<j<<std::endl;
    /////////////////////////////////////////

    rho = 0;
  }// Alla fine di ogni ciclo temporale ho bisogno di azzerare la carica sul nodo


  void Node::reset_velocities(){

    /////////////////////////////////////////
    if(debug)
      std::cout<<" Resetting velocities at Node whose indexes are  "<<i<<" - "<<j<<std::endl;
    /////////////////////////////////////////

    vx=0;
    vy=0;
  }



  void Node::set_electric_field(double Exx, double Eyy){
    Ex=Exx;
    Ey=Eyy;

    /////////////////////////////////////////
    if(debug)
      std::cout<<" Setting Ex = "<<Exx<<"  and Ey = "<<Eyy<<" on Node whose indexes are  "<<i<<"  and  "<<j<<std::endl;
    /////////////////////////////////////////

  }



  bool Node::operator<(const Node & l)const{
    if (i<l.i){
      return true;
    }
    else if (i==l.i){
      if(j<l.j){
        return true;
      }
      else
      return false;
    }
    else
    return false;
    }


