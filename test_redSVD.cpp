#include "Eigen/Core"
#include "RedSVD/RedSVD.hpp"
#include<iostream>
int main()
{
  unsigned int constexpr nrow=1000;
  unsigned int constexpr ncol=100;
  Eigen::MatrixXd A(nrow,ncol);
  // A matrix
  for(unsigned int i=0; i<nrow; ++i)
    for(unsigned int j=0; j<ncol; ++j)
      A(i,j)=1./(1+j+i);
  // try to change the rank
  unsigned int rank=10;
  RedSVD::RedSVD redsvd(A,rank);
  // never use auto with Eigen matrices!
  Eigen::MatrixXd  U=redsvd.matrixU();
  Eigen::MatrixXd  V=redsvd.matrixV().transpose();
  Eigen::MatrixXd  S=redsvd.singularValues().asDiagonal();
  auto actual_rank=S.rows();
  std::cout<<U.rows()<<" "<<U.cols()<<" "<<S.rows()<<" "<<S.cols();
  std::cout<<V.rows()<<" "<<V.cols()<<std::endl;
  auto As = U*S*V;
  auto diff = (A-As).norm();
  std::cout<<"Error  "<<diff<<std::endl;
  std::cout<<"Actual rank "<<actual_rank<<std::endl;
  std::cout<<"First "<<actual_rank<<" approximated singular values:"<<std::endl;
  for (unsigned int i=0;i<actual_rank;++i)
    {
      std::cout<<S(i,i)<<" ";
      if(i+1 % 8 == 0)std::cout<<std::endl;
    }
  if(actual_rank %8 !=0)std::cout<<std::endl;
}
