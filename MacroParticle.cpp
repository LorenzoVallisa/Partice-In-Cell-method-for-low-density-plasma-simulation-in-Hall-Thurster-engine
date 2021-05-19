#include "MacroParticle.h"

//routine per generare un numero casuale
//tra 0 ed L, distribuzione uniforme
double MacroParticle::rand_pos(const double max,const double min=0.){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> distrib(min,max);

	////////////////////////////////
	if(debug)
		std::cout<<" Position set "<<distrib(gen)<<std::endl;
	////////////////////////////////

	return distrib(gen) ;
}

//routine per generare un numero casuale
//da una distribuzione normale (non standard)
double MacroParticle::rand_vel(const double mu, const double sigma){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<double> distrib(mu,sigma);
	return distrib(gen);
}

//Costruttore di default:
//genera una macroparticella con posizione casuale,
//distribuzione uniforme nella mesh, e velocità casuale,
//distribuzione di Maxwell Boltzmann (situazione di equilibrio)
MacroParticle::MacroParticle(){

	////////////////////////////////
	if(debug)
		std::cout<<"Initializing Macroparticle with no velocity "<<std::endl;
////////////////////////////////

	x=rand_pos(constant::Lx);
	y=rand_pos(constant::Ly);
	vx=rand_vel(0.,constant::vth/2);
	vy=rand_vel(0.,constant::vth/2);
}

//Costruttore:
//genera una particella con y casuale-uniforme sul bordo sinistro,
//velocità longitudinale perturbata da una componente di drift,
//velocità assiale maxwelliana di default, (condizione di inflow)
MacroParticle::MacroParticle(double vdriftx):x(0.){

	////////////////////////////
	if(debug)
			std::cout<<"Initializing Macroparticle"<<std::endl;
	////////////////////////////


		y=rand_pos(constant::Ly);
		vx=rand_vel(0.,constant::vth)+vdriftx;
		vy=rand_vel(0.,constant::vth/2);
}

//compute the cell indexes of the particle
std::pair<std::size_t,std::size_t> MacroParticle::indexes() const {
	std::size_t i=static_cast<std::size_t> (x/constant::lambda);
	std::size_t j=static_cast<std::size_t> (y/constant::lambda);
	return std::make_pair(i,j);
}
