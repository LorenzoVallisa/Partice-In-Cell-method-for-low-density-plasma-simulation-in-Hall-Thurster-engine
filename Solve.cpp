#include "Solve.h"
#include <fstream>
#include <cmath>
#include <numeric>

//Costruttore:

Solve::Solve():A(constant::nodes){

	/////////////////////////////////////////
	if(debug)
		std::cout<<"Initializing A matrix"<<std::endl;
	/////////////////////////////////////////

	using constant::x_nodes;
	using constant::y_nodes;
	using constant::lambda;
	using constant::nodes;

	size_t k{0}; //helper index

	//Nodi interni:
	for(size_t j=1; j<y_nodes-1; ++j){
		for(size_t i=1; i<x_nodes-1; ++i){
			k=j*x_nodes+i; //--> max k= (ny-1)*nx +nx-1=nn-1
			A(k,k) = -4/(lambda*lambda);
      A(k,k-1)=1/(lambda*lambda);
      if(k+1<nodes)
      	A(k,k+1)=1/(lambda*lambda);
      A(k,k-x_nodes)=1/(lambda*lambda);
      if(k+x_nodes<nodes)
      	A(k,k+x_nodes)=1/(lambda*lambda);
		}
	}

	/////////////////////////////////////////
	if(debug)
		std::cout<<"Initializing Boundary Conditions of A"<<std::endl;
	/////////////////////////////////////////

	//Neumann bc sotto (y=0)
	for(size_t i=0; i<x_nodes; ++i){
		A(i,i)=-1/lambda;
		if(i==0 || i==x_nodes-1)
			A(i,i+x_nodes)=0.;
		else
			A(i,i+x_nodes)=1/lambda;
	}

	//Neumann BC sopra (y=Ly)
	for(size_t i=0; i<x_nodes; ++i){
		k=(y_nodes-1)*x_nodes+i; //max k=(ny-1)*nx+nx>nn
		if(i==0 || i==x_nodes-1)
			A(k,k-x_nodes)=0.;
		else
			A(k,k-x_nodes)=1/lambda;
		A(k,k)=-1/lambda;
	}

	//Neumann BC a destra (x=Lx)
	for(size_t i=0; i<y_nodes; ++i){
		k=(i+1)*x_nodes-1;		//max k=ny*nx-1=nn-1 ok
		for(size_t j=0; j<y_nodes; ++j)
			A(k,j)=0.;
		A(k,k-1)=1/lambda;
		A(k,k)=-1/lambda;
	}

	/////////////////////////////////////////
	if(debug)
		std::cout<<"End Neumann "<<std::endl;
	/////////////////////////////////////////

	//Dirichlet bc a sinistra (x=Ly)
	for(size_t i=0; i<y_nodes; ++i){
		k=i*x_nodes;	//max k=(ny-1)*nx<nn ok
		for(size_t j=0; j<y_nodes; ++j)
			A(k,j)=0.;
		A(k,k)=1.;
	}

		/////////////////////////////////////////
		if(debug){
			std::cout<<"Writing Matrix"<<std::endl;
		std::ofstream os("PoissonMatrix.csv");
		A.print(os);
		}
  	/////////////////////////////////////////
	
}

std::size_t Solve::get_nodes_index(std::size_t i,std::size_t j)const{

	//Metodo ausiliario che data una coppia di indici (di un nodo) trova il corrispettivo indice
	//del vettore di double del campo elettrostatico

	return i*(constant::y_nodes)+j;

}


double Solve::compute_residual(const std::array<double,constant::nodes_number> & bb,const std::array<double,constant::nodes_number> & xx)const{

	std::vector<double> x(xx.begin(),xx.end());
	std::vector<double> x1 = A*x;
	std::vector<double> res_sq(160);
	std::transform(bb.begin(),bb.end(),x1.begin(),res_sq.begin(),std::minus<double>());
	double norm_squared = std::inner_product(res_sq.begin(),res_sq.end(),res_sq.begin(),0);

	/////////////////////////////////////////
	if(debug)
		std::cout<<" Norm Squared "<<norm_squared<<std::endl;
	/////////////////////////////////////////

	double res = std::sqrt(norm_squared);
	return res;
}



std::array<double,constant::nodes_number> Solve::convert_row_col(const std::array<double,constant::nodes_number> & input)const{

	std::array<double,constant::nodes_number> inter;
	for(std::size_t i = 0;i<constant::nodes_number;++i)
		inter[i] = input[(i%constant::x_nodes)*constant::y_nodes + i/constant::x_nodes];
	return inter;
}



std::array<double,constant::nodes_number> Solve::convert_col_row(const std::array<double,constant::nodes_number> & input)const{

	std::array<double,constant::nodes_number> inter;
	for(std::size_t i = 0;i<constant::nodes_number;++i)
		inter[i] = input[(i%constant::y_nodes)*constant::x_nodes + i/constant::y_nodes];
	return inter;
}


std::tuple<double,double,double,double>
Solve::eulero_avanti(const double& vx0, const double& vy0, const double& Ex, const double& Ey,
			   const double& x0, const double& y0, const double& dt, const double& mult_fact){
	double vx{vx0},vy{vy0},x{x0},y{y0};
	vx+=mult_fact * Ex * dt;
  vy+=mult_fact * Ey * dt;
  x+=vx*dt;
  y+=vy*dt;
	return std::make_tuple(x,y,vx,vy);
}


std::tuple<double,double,double,double>
Solve::leapfrog(const double& vx0, const double& vy0, const double& Ex, const double& Ey,
       const double& x0, const double& y0, const double& dt, const double& mult_fact=1.){
  double vx{vx0},vy{vy0},x{x0},y{y0};
  //Eulero Avanti per le velocità, half time step
  vx+=mult_fact * Ex * dt/2;
  vy+=mult_fact * Ey * dt/2;
  //Eulero Avanti per le posizioni, usando la nuova velocità
  x+=vx*dt;
  y+=vy*dt;
  //Eulero Avanti per le velocità, half time step
  vx+=mult_fact * Ex * dt/2;
  vy+=mult_fact * Ey* dt/2;

  return std::make_tuple(x,y,vx,vy);
}



std::array<std::pair<double,double>,constant::nodes_number> Solve::compute_electric_field( const std::array<double,constant::nodes_number> & phi){

	// Method that takes as input a vector of double containing values of electrostatic potential for every node.
	// Goal is to implement algorithm for 2d gradient through spatial finite differences (centred)
	// Notare che phi è un vettore lungo nx*ny

	std::array<std::pair<double,double>,constant::nodes_number> Ef;

			/////////////////////////////////////////
	if(debug)
			std::cout<<" Initializing vector of Electric Field along x "<<std::endl;
			/////////////////////////////////////////



	for(std::size_t j=0;j<constant::y_nodes;++j){
		for(std::size_t i=0;i<constant::x_nodes-2;++i){
			// Nodi interni
			Ef[get_nodes_index(i+1,j)].first=1/(2*constant::lambda)*(phi[get_nodes_index(i,j)]-phi[get_nodes_index(i+2,j)]);
		}
		// x = 0
		Ef[get_nodes_index(0,j)].first = 1/(constant::lambda)*(phi[get_nodes_index(0,j)]-phi[get_nodes_index(1,j)]);
		// x = Lx
		Ef[get_nodes_index(constant::x_nodes-1,j)].first = 1/(constant::lambda)*(phi[get_nodes_index(constant::x_nodes-2,j)]-phi[get_nodes_index(constant::x_nodes-1,j)]);
	}




			/////////////////////////////////////////
	if(debug){
			std::cout<<" Initializing vector of Electric Field along y "<<std::endl;
	}
			/////////////////////////////////////////


			for(std::size_t i=0;i<constant::x_nodes;++i){
				for(std::size_t j=0;j<constant::y_nodes-2;++j){
					//Internal nodes
					Ef[get_nodes_index(i,j+1)].second =1/(2*constant::lambda)*(phi[get_nodes_index(i,j)]-phi[get_nodes_index(i,j+2)]);
				}
				// y=0
				Ef[get_nodes_index(i,0)].second = 1/(constant::lambda)*(phi[get_nodes_index(i,0)]-phi[get_nodes_index(i,1)]);
				//y = Ly
				Ef[get_nodes_index(i,constant::y_nodes-1)].second = 1/(constant::lambda)*(phi[get_nodes_index(i,constant::y_nodes-2)]-phi[get_nodes_index(i,constant::y_nodes-1)]);
			}


			/////////////////////////////////////////
			if(debug){
				std::cout<<"X component of Electric Field    "<<std::endl;
				for(const auto & it : Ef)
					std::cout<<" - "<<it.first<<" - "<<std::endl;
				std::cout<<"Y component of Electric Field    "<<std::endl;
				for(const auto & it : Ef)
					std::cout<<" - "<<it.second<<" - "<<std::endl;
				}
			/////////////////////////////////////////


	return Ef;

	// Giustificazione uso due cicli for: usando functional (minus) e algoritmo per inizializzare subranges di vettori avrei creato tante coppie inutili

}


std::array<double,constant::nodes_number> Solve::poisson(const std::array<double,constant::nodes_number> & rho_ions,
	const std::array<double,constant::nodes_number> & phi_start){

/////////////////////////////////////////
if(debug)
		std::cout<<"Calling Poisson Solver"<<std::endl;
/////////////////////////////////////////

	//1.Converto i vettori

	std::array<double,constant::nodes_number> rho_ions_ready(convert_row_col(rho_ions));

	std::array<double,constant::nodes_number> phi_start_ready(convert_row_col(phi_start));

	//2.Costruisco b iniziale: per quello mi serve poter far differenza tra array, oltre che appunto poter valutare un array in una funzione.



	std::array<double,constant::nodes_number> b_inter,b,x_electron,x(phi_start_ready);



	auto val_exp=[](double alfa){ return constant::n0*std::exp(alfa/constant::Te);};

	auto val_b=[](double beta){ return -constant::Qe*beta/constant::Eps0;};

	// Valuto phi_start_ready nell'espressione della forzante di Poisson (formula per densità elettronica
	// secondo distribuzione Boltzmann) e lo salvo in x_0


	//INIZIO CICLO ITERAZIONI PER CALCOLARMI PHI



	for (std::size_t nit=0;nit<2000;++nit){

		/////////////////////////////////////////
		if(debug)
		std::cout<<"Entering cycle number nit "<<nit<<std::endl;
		/////////////////////////////////////////

		std::transform(x.begin(),x.end(),x_electron.begin(),val_exp);


		//Termino rhs includendo termine ioni: ora tutto è salvato in b

		std::transform(rho_ions_ready.begin(),rho_ions_ready.end(),x_electron.begin(),b_inter.begin(),std::minus<double>());


		// Ultima normalizzazione di b per lambda

		std::transform(b_inter.begin(),b_inter.end(),b.begin(),val_b);



				for(std::size_t i = 0;i<constant::x_nodes;++i){
					//y = 0
					b[i]=0;
					//y=L
					b[constant::nodes_number+1-constant::x_nodes + i]=0;

					/////////////////////////////////////////
					if(debug)
					 	std::cout<<"Setting "<<i<< "  and  "<<constant::nodes_number+1-constant::x_nodes + i<<" indexes equal to zero on RHS of Poisson Solver(condition on y)"<<std::endl;
					/////////////////////////////////////////

				}

				for(std::size_t j = 0;j<constant::y_nodes;++j){
					//x = L
					b[(j+1)*constant::x_nodes -1 ] = 0;
					//x = 0
					b[j*constant::x_nodes] = 0;

					/////////////////////////////////////////
					if(debug)
						std::cout<<"Setting "<<(j+1)*constant::x_nodes -1<< "  and  "<<j*constant::x_nodes<<" indexes equal to zero on RHS of Poisson Solver(condition on x)"<<std::endl;
					/////////////////////////////////////////

				}

				/////////////////////////////////////////
				if(debug){
					std::cout<<"Size of b is    "<<b.size()<<std::endl;
					for(const auto & it : b)
						std::cout<<" - "<<it<<" - "<<std::endl;
					}
				/////////////////////////////////////////


	//Calcolo x (phi) con Gauss_Seidel (metodo risoluzione sistemi lineari)


		/////////////////////////////////////////
		if(debug)
			std::cout<<"Entering the cycle to compute phi"<<std::endl;
		/////////////////////////////////////////

			for (std::size_t jj = 0;jj<constant::nodes_number;++jj){


				if(jj==0){


					x[jj] = (b[jj] - A.scalar_product(jj,jj+1,constant::nodes_number-1,&x[jj+1]))/A(jj,jj);

				}

				else if (jj==constant::nodes_number-1){

				x[jj] = (b[jj] - A.scalar_product(jj,0,jj-1,&x[0]))/A(jj,jj);


				}

			else

			{

				x[jj] = (b[jj] - A.scalar_product(jj,0,jj-1,&x[0])- A.scalar_product(jj,jj+1,constant::nodes_number-1,&x[jj+1]))/A(jj,jj);


			}

		}


		/////////////////////////////////////////
		if(debug){
			std::cout<<"Size of x after cicle number "<<nit<<" before entering residual check is    "<<x.size()<<std::endl;
			for(const auto & it : x)
				std::cout<<" - "<<it<<" - "<<std::endl;

			//Further debug to check every cycle
			// bool val;
			// std::cout<<" Continue 1 - Abort 0"<<std::endl;
			// std::cin>>val;
			// if(!val)
			// 	break;

			}

		/////////////////////////////////////////





		if(nit%10==0){

			/////////////////////////////////////////
			if(debug)
			std::cout<<"Residual = "<<compute_residual(b,x)<<"  at iteration   "<<nit<<std::endl;
			/////////////////////////////////////////


			if(tol>compute_residual(b,x)){

				/////////////////////////////////////////
				if(debug){
					std::cout<<"Returning final x"<<std::endl;
					for(const auto & it : x)
						std::cout<<" - "<<it<<" - "<<std::endl;
					}
				/////////////////////////////////////////

				std::array<double,constant::nodes_number> final_x(convert_col_row(x));
				return final_x;
			}


				/////////////////////////////////////////
				if(debug)
					std::cout<<"Residual checked"<<std::endl;
				/////////////////////////////////////////


		}




	}

	/////////////////////////////////////////
	if(debug){
		std::cout<<"Returning final x"<<std::endl;
		for(const auto & it : x)
			std::cout<<" - "<<it<<" - "<<std::endl;
		}
	/////////////////////////////////////////

		std::array<double,constant::nodes_number> final_x(convert_col_row(x));
		return final_x;
}
