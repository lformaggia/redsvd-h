/* 
 * A header-only version of RedSVD
 * 
 * Copyright (c) 2014 Nicolas Tessore
 * 
 * based on RedSVD
 * 
 * Copyright (c) 2010 Daisuke Okanohara
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 * 1. Redistributions of source code must retain the above Copyright
 *    notice, this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above Copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the authors nor the names of its contributors
 *    may be uses to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 *    Minor modifications by Luca Formaggia.
 */

#ifndef REDSVD_MODULE_H
#define REDSVD_MODULE_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cstdlib>
#include <cmath>
#include <random>

namespace RedSVD
{
  //! It generates a (pseudo) x and y extracted from a (pseudo) random Gaussian distribution
  /*!
   * It uses the Box-Muller transform, see 
   * [this Wikipedia page](https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform),
   * to generate a couple of variables from a normally distribution
   * with zero mean and unitary variance.  
   *
   * This version uses rand() to generate the uniform distribution from
   * which the Gaussian one is derived.  It may be changed to implement
   * the new random distribution.
   *
   * \note Compared with the original version, it is no more used to generate
   * random Gaussian matrices, since we use the random number generator
   * of the standard library.
   */
  template<typename Scalar>
  inline void sample_gaussian(Scalar& x, Scalar& y)
  {
    using std::sqrt;
    using std::log;
    using std::cos;
    using std::sin;
		
    constexpr Scalar PI(3.1415926535897932384626433832795028841971693993751);
		
    Scalar v1 = (Scalar)(std::rand() + Scalar(1)) / ((Scalar)RAND_MAX+Scalar(2));
    Scalar v2 = (Scalar)(std::rand() + Scalar(1)) / ((Scalar)RAND_MAX+Scalar(2));
    Scalar len = sqrt(Scalar(-2) * log(v1));
    x = len * cos(Scalar(2) * PI * v2);
    y = len * sin(Scalar(2) * PI * v2);
  }

  //! Generates a random Gaussian matrix
  /*! A matrix whose elements are drawn from Normal distribution with
   * zero mean and unitary variance.
   *
   * \note Modified from the original code by Luca Formaggia: using standard library random number distribution.
  */
  template<typename MatrixType>
  inline void sample_gaussian(MatrixType& mat)
  {
    // generate random seed
    std::random_device rd{};
    std::mt19937 gen{rd()};
    // define the distribution
    typedef typename MatrixType::Scalar Scalar;
    std::normal_distribution<Scalar> d;

    typedef typename MatrixType::Index Index;
          
    for(Index i = 0; i < mat.rows(); ++i)
      for(Index j = 0; j < mat.cols(); ++j)
        mat(i, j)=d(gen);
  }
  //! Performs Gram Schmidt on the matrix columns
  /*!
   * \tparam MatrixType An Eigen Dense Matrix
   * \param mat The matrix on which to operate. The orthonormalisation is performed in-place
   * \param EPS tolerance to discard small columns. Defaulted to 1e-6
   */
  template<typename MatrixType>
  inline void gram_schmidt(MatrixType& mat,
                           typename MatrixType::Scalar const EPS=1.E-6)
  {
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
		
    for(Index i = 0; i < mat.cols(); ++i)
      {
        // c_i = c_i - sum_{j<i} (c_i * c_j)c_j
        for(Index j = 0; j < i; ++j)
          {
            Scalar r = mat.col(i).dot(mat.col(j));
            mat.col(i) -= r * mat.col(j);
          }
			
        Scalar norm = mat.col(i).norm();
        // If the norm is too small it means that rank A = i-1
        // so we can put=0 all columns k with k>= i.
        if(norm < EPS)
          {
            for(Index k = i; k < mat.cols(); ++k)
              mat.col(k).setZero();
            return;
          }
        // Orhonormalization
        mat.col(i) /= norm;
      }
  }

  //! Performs reduced Singular Value Decomposition
  /*!
   * Given a nxm matrix \f$A\f$ it computes matrices \f$U\f$ (nxk, orthogonal), \f$\Sigma\f$ (kxk, diagonal) and
   * \f$V\f$ (mxk, orthogonal), so that \f$A^*=U\Sigma V^T\f$ is a rank-k matrix approximation of \f$A\f$.
   * Matrix \f$\Sigma\f$ is returned as vector of size k containing the diagonal elements, which are an approximation
   * of the first singular values of \f$A\f$.
   *
   * This algorithm uses is a probabilistic algorithm. More details on this type of algorithms may be found in
   * <em>Finding structures with randomness:..., Halko, N. and Martinsson, P. G. and Tropp, J. A. (2009)</em>
   */
  template<typename _MatrixType>
  class RedSVD
  {
  public:
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ScalarVector;
    //! Default constructor
    RedSVD()=default;

    //! A may be a dense or sparse Eigen matrix.
    /*!
      Here rank is taken as min between n. rows and n. columns
      No great memory savings.
    */
    RedSVD(const MatrixType& A)
    {
      int r = (A.rows() < A.cols()) ? A.rows() : A.cols();
      compute(A, r);
    }
    //! Specifying the rank
    /*!
      \param A a sparse or dense Eigen matrix
      \param rank the desired rank for the reduced decomposition
    */
    RedSVD(const MatrixType& A, const Index rank)
    {
      compute(A, rank);
    }
    /*!
      \brief Computes the actual reduced SVD
      It is automatically called by constructors that take a matrix in input.
    */
    void compute(const MatrixType& A, const Index rank)
    {
      if(A.cols() == 0 || A.rows() == 0)
        return;
			
      Index r = (rank < A.cols()) ? rank : A.cols();
			
      r = (r < A.rows()) ? r : A.rows();
			
      // Gaussian Random Matrix for A^T
      DenseMatrix O(A.rows(), r);
      sample_gaussian(O);
			
      // Compute Sample Matrix of A^T
      DenseMatrix Y = A.transpose() * O;
			
      // Orthonormalize Y
      gram_schmidt(Y);
			
      // Range(B) = Range(A^T)
      DenseMatrix B = A * Y;
			
      // Gaussian Random Matrix
      DenseMatrix P(B.cols(), r);
      sample_gaussian(P);
			
      // Compute Sample Matrix of B
      DenseMatrix Z = B * P;
			
      // Orthonormalize Z
      gram_schmidt(Z);
			
      // Range(C) = Range(B)
      DenseMatrix C = Z.transpose() * B; 
      // Thin SVD. 			
      Eigen::JacobiSVD<DenseMatrix> svdOfC(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
			
      // C = USV^T
      // A^* = (Z * U) * S * (Y * V)^T (low rank approx of A).
      m_matrixU = Z * svdOfC.matrixU();
      m_vectorS = svdOfC.singularValues();
      m_matrixV = Y * svdOfC.matrixV();
    }

    //! Returns matrix U
    DenseMatrix matrixU() const
    {
      return m_matrixU;
    }
    //! Returns vector with first k singular values    
    ScalarVector singularValues() const
    {
      return m_vectorS;
    }
    //! Returns matrix V
    DenseMatrix matrixV() const
    {
      return m_matrixV;
    }
		
  private:
    DenseMatrix m_matrixU;
    ScalarVector m_vectorS;
    DenseMatrix m_matrixV;
  };

  //! Like RedSVD but for symmetric matrices
  /*!
   * In this case, \f$U=V\f$ is the matrix with the first (approximate) eigenvectors of \f$A\f$ and \f$\Sigma\f$ contains
   * the (approximate) first k eigenvalues.
   */
  template<typename _MatrixType>
  class RedSymEigen
  {
  public:
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ScalarVector;
		
    RedSymEigen()=default;
		
    RedSymEigen(const MatrixType& A)
    {
      int r = (A.rows() < A.cols()) ? A.rows() : A.cols();
      compute(A, r);
    }
		
    RedSymEigen(const MatrixType& A, const Index rank)
    {
      compute(A, rank);
    }  
		
    void compute(const MatrixType& A, const Index rank)
    {
      if(A.cols() == 0 || A.rows() == 0)
        return;
			
      Index r = (rank < A.cols()) ? rank : A.cols();
			
      r = (r < A.rows()) ? r : A.rows();
			
      // Gaussian Random Matrix
      DenseMatrix O(A.rows(), r);
      sample_gaussian(O);
			
      // Compute Sample Matrix of A
      DenseMatrix Y = A.transpose() * O;
			
      // Orthonormalize Y
      gram_schmidt(Y);
			
      DenseMatrix B = Y.transpose() * A * Y;
      Eigen::SelfAdjointEigenSolver<DenseMatrix> eigenOfB(B);
			
      m_eigenvalues = eigenOfB.eigenvalues();
      m_eigenvectors = Y * eigenOfB.eigenvectors();
    }
    //! Returns a vector with the first k approximate eigenvalues
    /*!
      They are the diagonal elements of \f$\Sigma\f$.
     */
    ScalarVector eigenvalues() const
    {
      return m_eigenvalues;
    }
    //! Returns a nxk orthogonal matrix whose columns are the first k approximate eigenvectors
    DenseMatrix eigenvectors() const
    {
      return m_eigenvectors;
    }
		
  private:
    ScalarVector m_eigenvalues;
    DenseMatrix m_eigenvectors;
  };

  //! Performs the principal component analysis
  /*!
   * It works by contructing a RedSVD object to perform an low-rank SVD decomposition and then
   * it returns the nxk matrix \f$V\f$ that represent the components of the PCA and the nxk matrix
   * \f$U\Sigma\f$, which are the scores.
   */
  template<typename _MatrixType>
  class RedPCA
  {
  public:
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ScalarVector;
		
    RedPCA() {}
		
    RedPCA(const MatrixType& A)
    {
      int r = (A.rows() < A.cols()) ? A.rows() : A.cols();
      compute(A, r);
    }
		
    RedPCA(const MatrixType& A, const Index rank)
    {
      compute(A, rank);
    }  
		
    void compute(const DenseMatrix& A, const Index rank)
    {
      RedSVD<MatrixType> redsvd(A, rank);
			
      ScalarVector S = redsvd.singularValues();
			
      m_components = redsvd.matrixV();
      m_scores = redsvd.matrixU() * S.asDiagonal();
    }
    //! Returns the components
    DenseMatrix components() const
    {
      return m_components;
    }
    //! Returns the scores		
    DenseMatrix scores() const
    {
      return m_scores;
    }
		
  private:
    DenseMatrix m_components;
    DenseMatrix m_scores;
  };
}

#endif
