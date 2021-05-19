#ifndef MPARTICLE_H
#define MPARTICLE_H

#include <utility>
#include <memory>
#include <random>
#include "Constant.h"
#include "Cell.h"


class MacroParticle{

	static constexpr bool debug = false;

private:

	//Posizione e velocità della macroparticella:
	double x{0.};
	double y{0.};
	double vx{0.};
	double vy{0.};
	//Campo elettrico calcolato nella posizione della particella:
	double Ex{0.};
	double Ey{0.};

	//Shared pointer alla cella corrente nella mesh
	std::shared_ptr<Cell> my_cell;

	//helper routines
	double rand_pos(const double , const double);

	double rand_vel(const double mu, const double sigma);

public:
	//Costruttore iniziale
	//importante che sia di default perché System usa quello di list di tipo list(std::size_t),
	//che a sua volta usa quello di default dei propri membri
	MacroParticle();
	//Costruttore usato nel loop per le particelle entranti da sinistra
	MacroParticle(double vdriftx);

	//Getters:
	double get_x() const {return x;}
	double get_y() const {return y;}
	double get_vx() const {return vx;}
	double get_vy() const {return vy;}
	double get_Ex() const {return Ex;}
	double get_Ey() const {return Ey;}
	std::pair<double,double> get_E() const {return std::make_pair(Ex,Ey);}
	std::shared_ptr<Cell> get_mycell() const {return my_cell;}

	//Setters:
	void set_pos(const double xx, const double yy){x=xx; y=yy;}
	void set_pos(const std::pair<double,double> p){x=p.first; y=p.second;}
	void set_v(const double vxx, const double vyy){vx=vxx; vy=vyy;}
	void set_E(const double Exx, const double Eyy){Ex=Exx; Ey=Eyy;}
	void set_cell(const std::shared_ptr<Cell> cella){my_cell=cella;}

	//compute the cell indexes of the particle
	std::pair<size_t,size_t> indexes() const;

};


#endif
