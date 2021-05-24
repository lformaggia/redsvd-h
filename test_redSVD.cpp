#include "Eigen/Core"
#include "Eigen/Sparse"
#include "RedSVD/RedSVD.hpp"
#include <algorithm>
#include <iostream>
int
main ()
{
  int constexpr nrow = 10000;
  int constexpr ncol = 1000;
  // Eigen::MatrixXd A(nrow,ncol);
  // A matrix
  // Testing a sparse matrix
  std::cout<<" Constructing the matrix"<<std::endl;
  Eigen::SparseMatrix<double> A (nrow, ncol);
  A.reserve (10 * nrow);
  for (int i = 0; i < nrow; ++i)
    for (int j = std::max (0, i - 4);
         j < std::min (i + 4, static_cast<int> (ncol)); ++j)
      A.coeffRef (i, j) = 1. / (1 + j + i);
  A.makeCompressed ();
  std::cout<<"Done"<<std::endl;
  // try to change the rank
  unsigned int   rank = 20;
  RedSVD::RedSVD redsvd (A, rank);
  // never use auto with Eigen matrices!
  Eigen::MatrixXd U = redsvd.matrixU ();
  Eigen::MatrixXd V = redsvd.matrixV ().transpose ();
  Eigen::MatrixXd S = redsvd.singularValues ().asDiagonal ();
  auto            actual_rank = S.rows ();
  std::cout << "U=" << U.rows () << "X" << U.cols () << ", S= " << S.rows ()
            << "X" << S.cols ();
  std::cout << " Vt=" << V.rows () << "X" << V.cols () << std::endl;
  auto As = U * S * V;
  auto diff = (A - As).norm ()/A.norm();
  std::cout << "Error  " << diff << std::endl;
  std::cout << "Actual rank " << actual_rank << std::endl;
  std::cout << "First " << actual_rank
            << " approximated singular values:" << std::endl;
  for (unsigned int i = 0; i < actual_rank; ++i)
    {
      std::cout << S (i, i) << " ";
      if (i + 1 % 8 == 0)
        std::cout << std::endl;
    }
  if (actual_rank % 8 != 0)
    std::cout << std::endl;
}
