#ifndef PIC_MATRIX_H
#define PIC_MATRIX_H

#include "Constant.h"

#include <vector>
#include <fstream>
#include <numeric>
#include <iostream>
#include <array>
#include <cstddef>

template <typename T>
class Matrix {

friend std::vector<T> operator * (const Matrix & A  ,std::vector<T> const & B);

private :

  std::vector<T> elements;
  const std::size_t rows;
  const std::size_t cols;

  inline	std::size_t	sub2ind (const std::size_t i, const std::size_t j) const
  { return (i + j * rows); }

  double & index (std::size_t i, std::size_t j)
  { return elements[sub2ind (i, j)]; }

  const double & const_index (std::size_t i, std::size_t j) const
  { return elements[sub2ind (i, j)]; }

public :

  Matrix (std::size_t sz): rows (sz), cols (sz)
  { elements.resize (rows * cols, 0.0); };

  Matrix (std::size_t r, std::size_t c): rows (r), cols (c)
  { elements.resize (rows * cols, 0.0); };

  Matrix (Matrix const &) = default;

	std::size_t get_rows () const { return rows; }

	std::size_t get_cols () const { return cols; }

  T & operator() (std::size_t i, std::size_t j)
  { return index (i, j); };

  const T & operator()  (std::size_t i, std::size_t j) const
  { return const_index (i, j); };

  const T * get_elements () const { return &(elements[0]); };

  T * get_elements () { return &(elements[0]); };

  std::ofstream & print(std::ofstream & os);

  double scalar_product(size_t i, size_t j1, size_t j2,
  											const std::array<T,constant::nodes_number>::iterator ptr)const;


};


#endif
