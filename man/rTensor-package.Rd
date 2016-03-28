\docType{package}
\name{rTensor-package}
\alias{rTensor}
\alias{rTensor-package}
\title{Tools for tensor analysis and decomposition}
\description{
This package is centered around the \code{\link{Tensor-class}}, which defines a S4 class for tensors of arbitrary number of modes. A vignette and/or a possible paper will be included in a future release of this package.
}
\details{
This page will summarize the full functionality of this package. Note that since all the methods associated with S4 class \code{\link{Tensor-class}} are documented there, we will not duplicate it here.
 
The remaining functions can be split into two groups: the first is a set of tensor decompositions, and the second is a set of helper functions that are useful in tensor manipulation.

rTensor implements the following tensor decompositions: \describe{
	\item{\code{\link{cp}}}{Canonical Polyadic (CP) decomposition}
	\item{\code{\link{tucker}}}{General Tucker decomposition}
	\item{\code{\link{mpca}}}{Multilinear Principal Component Analysis; note that for 3-Tensors this is also known as Generalized Low Rank Approximation of Matrices(GLRAM)}
	\item{\code{\link{hosvd}}}{(Truncated-)Higher-order singular value decomposition}
	\item{\code{\link{t_svd}}}{Tensor singular value decomposition; 3-Tensors only; also note that there is an asociated reconstruction function \code{\link{t_svd_reconstruct}}}
	\item{\code{\link{pvd}}}{Population value decomposition of images; 3-Tensors only}
}

rTensor also provides a set functions for tensors multiplication: \describe{
	\item{\code{\link{ttm}}}{Tensor times matrix, aka m-mode product}
	\item{\code{\link{ttl}}}{Tensor times list (of matrices)}
	\item{\code{\link{t_mult}}}{Tensor product based on block circulant unfolding; only implemented for a pair of 3-Tensors}
}

...as well as for matrices: \describe{
	\item{\code{\link{hamadard_list}}}{Computes the Hamadard (element-wise) product of a list of matrices}
	\item{\code{\link{kronecker_list}}}{Computes the Kronecker product of a list of matrices}
	\item{\code{\link{khatri_rao}}}{Computes the Khatri-Rao product of two matrices}
	\item{\code{\link{khatri_rao_list}}}{Computes the Khatri-Rao product of a list of matrices}
	\item{\code{\link{fold}}}{General folding of a matrix into a tensor}
	\item{\code{\link{k_fold}}}{Inverse operation for \code{\link{k_unfold}}}
	\item{\code{\link{unmatvec}}}{Inverse operation for \code{\link{matvec}}}
}

For more information on any of the functions, please consult the individual man pages.
}
\author{
  James Li \email{jamesyili@gmail.com}, Jacob Bien, and Martin T. Wells
}
