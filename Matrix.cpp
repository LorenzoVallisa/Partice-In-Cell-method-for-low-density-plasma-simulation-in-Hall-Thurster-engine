#include "Matrix.h"


std::vector<double> operator * (const Matrix & A  ,std::vector<double> const & B){
  using size_type = std::vector<double>::size_type;

  std::vector<double> C (A.rows);

  for (size_type i = 0; i < A.rows; ++i)
    for (size_type k = 0; k < A.cols; ++k)
        C[i] += A(i, k) * B[k];

  return C;
}


std::ofstream & Matrix::print(std::ofstream & os){
  for(std::size_t i=0; i<rows; ++i){
		for(std::size_t j=0; j<cols; ++j){
			os<<const_index(i,j);
			if(j<cols-1) os<<",";
		}
		if(i<rows-1) os<<"\n";
  }
  return os;
}




double Matrix::scalar_product(size_t i, size_t j1, size_t j2, 
									const std::array<double,constant::nodes_number>::iterator ptr)const{
  std::vector<double> v;
  for(size_t j=j1; j<=j2; ++j)
    v.push_back(const_index(i,j));
  return std::inner_product(v.begin(),v.end(),ptr,0.);
}
