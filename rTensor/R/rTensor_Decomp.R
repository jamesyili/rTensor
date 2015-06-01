###Tensor Decompositions

#'(Truncated-)Higher-order SVD
#'
#'Higher-order SVD of a K-Tensor. Write the K-Tensor as a (m-mode) product of a core Tensor (possibly smaller modes) and K orthogonal factor matrices. Truncations can be specified via \code{ranks} (making them smaller than the original modes of the K-Tensor will result in a truncation). For the mathematical details on HOSVD, consult Lathauwer et. al. (2000).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure. A progress bar is included to help monitor operations on large tensors.
#'@name hosvd
#'@rdname hosvd
#'@aliases hosvd
#'@param tnsr Tensor with K modes
#'@param ranks a vector of desired modes in the output core tensor, default is \code{tnsr@@modes}
#'@return a list containing the following:\describe{
#'\item{\code{Z}}{core tensor with modes speficied by \code{ranks}}
#'\item{\code{U}}{a list of orthogonal matrices, one for each mode}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)} - if there was no truncation, then this is O(mach_eps) }
#'}
#'@seealso \code{\link{tucker}}
#'@references L. Lathauwer, B.Moor, J. Vanderwalle "A multilinear singular value decomposition". Journal of Matrix Analysis and Applications 2000.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes}.
#'@examples
#'tnsr <- rand_tensor(c(6,7,8))
#'hosvdD <-hosvd(tnsr)
#'hosvdD$fnorm_resid
#'hosvdD2 <-hosvd(tnsr,ranks=c(3,3,4))
#'hosvdD2$fnorm_resid
hosvd <- function(tnsr,ranks=NULL){
	#stopifnot(is(tnsr,"Tensor"))
	num_modes <- tnsr@num_modes
	#no truncation if ranks not provided
	if(is.null(ranks)){
		ranks <- tnsr@modes
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=num_modes,style=3)
	#loops through and performs SVD on mode-m matricization of tnsr
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		temp_mat <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
		setTxtProgressBar(pb,m)
	}
	close(pb)
	#computes the core tensor
	Z <- ttl(tnsr,lapply(U_list,t),ms=1:num_modes)
	est <- ttl(Z,U_list,ms=1:num_modes)
	resid <- fnorm(est-tnsr)
	#put together the return list, and returns
	list(Z=Z,U=U_list,est=est,fnorm_resid=resid)	
}

#'Canonical Polyadic Decomposition
#'
#'Canonical Polyadic (CP) decomposition of a tensor, aka CANDECOMP/PARAFRAC. Approximate a K-Tensor using a sum of \code{num_components} rank-1 K-Tensors. A rank-1 K-Tensor can be written as an outer product of K vectors. There are a total of \code{num_compoents *tnsr@@num_modes} vectors in the output, stored in \code{tnsr@@num_modes} matrices, each with \code{num_components} columns. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on CP decomposition, consult Kolda and Bader (2009).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure. A progress bar is included to help monitor operations on large tensors.
#'@name cp
#'@rdname cp
#'@aliases cp
#'@param tnsr Tensor with K modes
#'@param num_components the number of rank-1 K-Tensors to use in approximation
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following \describe{
#'\item{\code{lambdas}}{a vector of normalizing constants, one for each component}
#'\item{\code{U}}{a list of matrices - one for each mode - each matrix with \code{num_components} columns}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{tucker}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@examples
#'tnsr <- rand_tensor(c(6,7,8))
#'cpD <- cp(tnsr,num_components=5) 
#'cpD$conv 
#'cpD$norm_percent 
#'plot(cpD$all_resids) 
cp <- function(tnsr, num_components=NULL,max_iter=25, tol=1e-5){
	if(is.null(num_components)) stop("num_components must be specified")
	stopifnot(is(tnsr,"Tensor"))
	#initialization via truncated hosvd
	num_modes <- tnsr@num_modes
	modes <- tnsr@modes
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	tnsr_norm <- fnorm(tnsr)
	for(m in 1:num_modes){
		unfolded_mat[[m]] <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
	}
	est <- tnsr
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resid <- rep(0, max_iter)
	CHECK_CONV <- function(est){
		curr_resid <- fnorm(est - tnsr)
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{ return(FALSE)}
	}	
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	norm_vec <- function(vec){
	norm(as.matrix(vec))
	}
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)
		for(m in 1:num_modes){
			V <- hamadard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
			V_inv <- solve(V)			
			tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],reverse=TRUE)%*%V_inv
			lambdas <- apply(tmp,2,norm_vec)
			U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
			Z <- .superdiagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
			est <- ttl(Z,U_list,ms=1:num_modes)
		}
		#checks convergence
		if(CHECK_CONV(est)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
			 }
	}
	if(!converged){setTxtProgressBar(pb,max_iter)}
	close(pb)
	#end of main loop
	#put together return list, and returns
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent<- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(lambdas=lambdas, U=U_list, conv=converged, est=est, norm_percent=norm_percent, fnorm_resid = tail(fnorm_resid,1),all_resids=fnorm_resid))
}

#'Tucker Decomposition
#'
#'The Tucker decomposition of a tensor. Approximates a K-Tensor using a n-mode product of a core tensor (with modes specified by \code{ranks}) with orthogonal factor matrices. If there is no truncation in one of the modes, then this is the same as the MPCA, \code{\link{mpca}}. If there is no truncation in all the modes (i.e. \code{ranks = tnsr@@modes}), then this is the same as the HOSVD, \code{\link{hosvd}}. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on the Tucker decomposition, consult Kolda and Bader (2009).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure also known as Higher-Order Orthogonal Iteration (HOOI). Intialized using a (Truncated-)HOSVD. A progress bar is included to help monitor operations on large tensors.
#'@name tucker
#'@rdname tucker
#'@aliases tucker
#'@param tnsr Tensor with K modes
#'@param ranks a vector of the modes of the output core Tensor
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{Z}}{the core tensor, with modes specified by \code{ranks}}
#'\item{\code{U}}{a list of orthgonal factor matrices - one for each mode, with the number of columns of the matrices given by \code{ranks}}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{hosvd}}, \code{\link{mpca}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes}.
#'@examples
#'tnsr <- rand_tensor(c(6,7,8))
#'tuckerD <- tucker(tnsr,ranks=c(3,3,4))
#'tuckerD$conv 
#'tuckerD$norm_percent
#'plot(tuckerD$all_resids)
tucker <- function(tnsr,ranks=NULL,max_iter=25,tol=1e-5){
	stopifnot(is(tnsr,"Tensor"))
	if(is.null(ranks)) stop("ranks must be specified")
	#initialization via truncated hosvd
	num_modes <- tnsr@num_modes
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		temp_mat <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
	}
	tnsr_norm <- fnorm(tnsr)
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resid <- rep(0, max_iter)
	CHECK_CONV <- function(Z,U_list){
		est <- ttl(Z,U_list,ms=1:num_modes)
		curr_resid <- fnorm(tnsr - est)
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)	
	modes <- tnsr@modes
	modes_seq <- 1:num_modes
		for(m in modes_seq){
			#core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-m],t),ms=modes_seq[-m])
			#truncated SVD of X
			#U_list[[m]] <- (svd(rs_unfold(X,m=m)@data,nu=ranks[m],nv=prod(modes[-m]))$u)[,1:ranks[m]]
			U_list[[m]] <- svd(rs_unfold(X,m=m)@data,nu=ranks[m])$u
		}
		#compute core tensor Z
		Z <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)

		#checks convergence
		if(CHECK_CONV(Z, U_list)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)	
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	close(pb)
	#end of main loop
	#put together return list, and returns
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent<-(1-(tail(fnorm_resid,1)/tnsr_norm))*100
	est <- ttl(Z,U_list,ms=1:num_modes)
	invisible(list(Z=Z, U=U_list, conv=converged, est=est, norm_percent = norm_percent, fnorm_resid=tail(fnorm_resid,1), all_resids=fnorm_resid))
}

#'Multilinear Principal Components Analysis
#'
#'This is basically the Tucker decomposition of a K-Tensor, \code{\link{tucker}}, with one of the modes uncompressed. If K = 3, then this is also known as the Generalized Low Rank Approximation of Matrices (GLRAM). This implementation assumes that the last mode is the measurement mode and hence uncompressed. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on the MPCA of tensors, consult Lu et al. (2008).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure. A progress bar is included to help monitor operations on large tensors.
#'@name mpca
#'@rdname mpca
#'@aliases mpca
#'@param tnsr Tensor with K modes
#'@param ranks a vector of the compressed modes of the output core Tensor, this has length K-1
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{Z_ext}}{the extended core tensor, with the first K-1 modes given by \code{ranks}}
#'\item{\code{U}}{a list of K-1 orthgonal factor matrices - one for each compressed mode, with the number of columns of the matrices given by \code{ranks}}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{tucker}}, \code{\link{hosvd}}
#'@references H. Lu, K. Plataniotis, A. Venetsanopoulos, "Mpca: Multilinear principal component analysis of tensor objects". IEEE Trans. Neural networks, 2008.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes-1}.
#'@examples
#'tnsr <-rand_tensor(c(100,10,10))
#'mpcaD <- mpca(tnsr,ranks=c(30,5))
#'mpcaD$conv
#'mpcaD$norm_percent
#'plot(mpcaD$all_resids)
mpca <- function(tnsr, ranks = NULL, max_iter = 25, tol=1e-5){
	if(is.null(ranks)) stop("ranks must be specified")
	stopifnot(is(tnsr,"Tensor"))
	#initialization via hosvd of M-1 modes
	num_modes <- tnsr@num_modes
	stopifnot(length(ranks)==(num_modes-1))
	ranks <- c(ranks,1)
	modes <- tnsr@modes
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	for(m in 1:(num_modes-1)){
		unfolded_mat <- rs_unfold(tnsr,m=m)@data
		mode_m_cov <- unfolded_mat%*%t(unfolded_mat)
		U_list[[m]] <- svd(mode_m_cov, nu=ranks[m])$u
	}
	Z_ext <- ttl(tnsr,lapply(U_list[-num_modes],t),ms=1:(num_modes-1))
	tnsr_norm <- fnorm(tnsr)
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resid <- rep(0, max_iter)
	CHECK_CONV <- function(Z_ext,U_list){
		est <- ttl(Z_ext,U_list[-num_modes],ms=1:(num_modes-1))
		curr_resid <- fnorm(tnsr - est)
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)
	modes <-tnsr@modes
	modes_seq <- 1:(num_modes-1)
		for(m in modes_seq){
			#extended core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-c(m,num_modes)],t),ms=modes_seq[-m])
			#truncated SVD of X
			U_list[[m]] <- svd(rs_unfold(X,m=m)@data,nu=ranks[m])$u
		}
		#compute core tensor Z_ext
		Z_ext <- ttm(X,mat=t(U_list[[num_modes-1]]),m=num_modes-1)
		#checks convergence
		if(CHECK_CONV(Z_ext, U_list)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	close(pb)
	#end of main loop
	#put together return list, and returns
	est <- ttl(Z_ext,U_list[-num_modes],ms=1:(num_modes-1))
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent<-(1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(Z_ext=Z_ext, U=U_list, conv=converged, est=est, norm_percent = norm_percent, fnorm_resid=tail(fnorm_resid,1), all_resids=fnorm_resid))
}

#'Population Value Decomposition
#'
#'The default Population Value Decomposition (PVD) of a series of 2D images. Constructs population-level matrices P, V, and D to account for variances within as well as across the images. Structurally similar to Tucker (\code{\link{tucker}}) and GLRAM (\code{\link{mpca}}), but retains crucial differences. Requires \code{2*n3 + 2} parameters to specified the final ranks of P, V, and D, where n3 is the third mode (how many images are in the set). Consult Crainiceanu et al. (2013) for the construction and rationale behind the PVD model.
#'@export
#'@details The PVD is not an iterative method, but instead relies on \code{n3 + 2}separate PCA decompositions. The third mode is for how many images are in the set.
#'@name pvd
#'@rdname pvd
#'@aliases pvd
#'@param tnsr 3-Tensor with the third mode being the measurement mode
#'@param uranks ranks of the U matrices
#'@param wranks ranks of the W matrices
#'@param a rank of \code{P = U\%*\%t(U)}
#'@param b rank of \code{D = W\%*\%t(W)}
#'@return a list containing the following:\describe{
#'\item{\code{P}}{population-level matrix \code{P = U\%*\%t(U)}, where U is constructed by stacking the truncated left eigenvectors of slicewise PCA along the third mode}
#'\item{\code{V}}{a list of image-level core matrices}
#'\item{\code{D}}{population-leve matrix \code{D = W\%*\%t(W)}, where W is constructed by stacking the truncated right eigenvectors of slicewise PCA along the third mode}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'}
#'@references C. Crainiceanu, B. Caffo, S. Luo, V. Zipunnikov, N. Punjabi, "Population value decomposition: a framework for the analysis of image populations". Journal of the American Statistical Association, 2013.
#'@examples
#'tnsr <- rand_tensor(c(10,5,100))
#'pvdD<-pvd(tnsr,uranks=rep(8,100),wranks=rep(4,100),a=8,b=4)
pvd <- function(tnsr,uranks=NULL,wranks=NULL,a=NULL,b=NULL){
	if(tnsr@num_modes!=3) stop("PVD only for 3D")
	if(is.null(uranks)||is.null(wranks)) stop("U and V ranks must be specified")
	if(is.null(a)||is.null(b)) stop("a and b must be specified")
	modes <- tnsr@modes
	n <- modes[3]
	if(length(uranks)!=n||length(wranks)!=n) stop("ranks must be of length n3")
	pb <- txtProgressBar(min=0,max=(n+3),style=3)
	x <- tnsr@data
	Us <- vector('list',n)
	Vs <- vector('list',n)
	S <- vector('list',n)
	for(i in 1:n){
		svdz <- svd(x[,,i],nu=uranks[i],nv=wranks[i])
		Us[[i]] <- svdz$u
		Vs[[i]] <- svdz$v
		S[[i]] <- svdz$d[1:min(uranks[i],wranks[i])]
		setTxtProgressBar(pb,i)
	}
	U <- matrix(unlist(Us),nrow=modes[1],ncol=sum(uranks)*n)
	#eigenU <- eigen(U%*%t(U))
	P <- eigen(U%*%t(U))$vectors[,1:a] #E-vecs of UU^T
	setTxtProgressBar(pb,n+1)
	V <- matrix(unlist(Vs),nrow=modes[2],ncol=sum(wranks)*n)
	#eigenV <- eigen(V%*%t(V))
	Dt <- eigen(V%*%t(V))$vectors[,1:b] #E-vecs of VV^T
	D <- t(Dt)
	setTxtProgressBar(pb,n+2)
	V2 <- vector('list',n)
	est <- array(0,dim=modes)
	for(i in 1:n){
		V2[[i]] <- (t(P)%*%Us[[i]])%*%diag(S[[i]],nrow=uranks[i],ncol=wranks[i])%*%(t(Vs[[i]])%*%Dt)
		est[,,i] <- P%*%V2[[i]]%*%D
	}
	est <- as.tensor(est)
	fnorm_resid <- fnorm(est-tnsr)	
	setTxtProgressBar(pb,n+3)
	norm_percent<-(1-(fnorm_resid/fnorm(tnsr)))*100
	invisible(list(P=P,D=D,V=V2,est=est,norm_percent=norm_percent,fnorm_resid=fnorm_resid))
}

#'Tensor Singular Value Decomposition
#'
#'TSVD for a 3-Tensor. Constructs 3-Tensors \code{U, S, V} such that \code{tnsr = t_mult(t_mult(U,S),t(V))}. \code{U} and \code{V} are orthgonal 3-Tensors with orthogonality defined in Kilmer et al. (2013), and \code{S} is a 3-Tensor consists of facewise diagonal matrices. For more details on the TSVD, consult Kilmer et al. (2013).
#'@export
#'@name t_svd
#'@rdname t_svd
#'@aliases t_svd
#'@param tnsr 3-Tensor to decompose via TSVD
#'@return a list containing the following:\describe{
#'\item{\code{U}}{the left orthgonal 3-Tensor}
#'\item{\code{V}}{the right orthgonal 3-Tensor}
#'\item{\code{S}}{the middle 3-Tensor consisting of face-wise diagonal matrices}
#'}
#'@seealso \code{\link{t_mult}}, \code{\link{t_svd_reconstruct}}
#'@references M. Kilmer, K. Braman, N. Hao, and R. Hoover, "Third-order tensors as operators on matrices: a theoretical and computational framework with applications in imaging". SIAM Journal on Matrix Analysis and Applications 2013.
#'@note Computation involves complex values, but if the inputs are real, then the outputs are also real. Some loss of precision occurs in the truncation of the imaginary components during the FFT and inverse FFT.
#'@examples
#'tnsr <- rand_tensor()
#'tsvdD <- t_svd(tnsr)
t_svd<-function(tnsr){
	if(tnsr@num_modes!=3) stop("T-SVD only implemented for 3d so far")
	modes <- tnsr@modes
	n1 <- modes[1]
	n2 <- modes[2]
	n3 <- modes[3]
	#progress bar
	pb <- txtProgressBar(min=0,max=n3,style=3)
	#define ifft
	ifft <- function(x){suppressWarnings(as.numeric(fft(x,inverse=TRUE))/length(x))}
	#fft for each of the n1n2 vectors (of length n3) along mode 3
	fftz <- aperm(apply(tnsr@data,MARGIN=1:2,fft),c(2,3,1))
	#svd for each face (svdz is a list of the results)
	U_arr <- array(0,dim=c(n1,n1,n3))
	V_arr <- array(0,dim=c(n2,n2,n3))
	m <- min(n1,n2)		
	S_arr <- array(0,dim=c(n1,n2,n3))
	#Think of a way to avoid a loop in the beginning
	#Problem is that svd returns a list but ideally we want 3 arrays
	#Even with unlist this doesn't seem possible
	for (j in 1:n3){
		setTxtProgressBar(pb,j)
		decomp <- svd(fftz[,,j],nu=n1,nv=n2)
		U_arr[,,j] <- decomp$u
		V_arr[,,j] <- decomp$v
		S_arr[,,j] <- diag(decomp$d,nrow=n1,ncol=n2) #length is min(n1,n2)
	}	
	close(pb)
	#for each svd result, we want to apply ifft
	U <- as.tensor(aperm(apply(U_arr,MARGIN=1:2,ifft),c(2,3,1)))
	V <- as.tensor(aperm(apply(V_arr,MARGIN=1:2,ifft),c(2,3,1)))
	S <- as.tensor(aperm(apply(S_arr,MARGIN=1:2,ifft),c(2,3,1)))
	invisible(list(U=U,V=V,S=S))
}

#'Reconstruct Tensor From TSVD
#'
#'Reconstruct the original 3-Tensor after it has been decomposed into \code{U, S, V} via \code{\link{t_svd}}.
#'@export
#'@name t_svd_reconstruct
#'@rdname t_svd_reconstruct
#'@aliases t_svd_reconstruct
#'@param L list that is an output from \code{\link{t_svd}}
#'@return a 3-Tensor 
#'@seealso \code{\link{t_svd}}
#'@examples
#'tnsr <- rand_tensor(c(10,10,10))
#'tsvdD <- t_svd(tnsr)
#'1 - fnorm(t_svd_reconstruct(tsvdD)-tnsr)/fnorm(tnsr)
t_svd_reconstruct <- function(L){
	t_mult(t_mult(L$U,L$S),t(L$V))
}

###t-compress (Not Supported)
.t_compress <- function(tnsr,k){
	modes <- tnsr@modes
	n1 <- modes[1]
	n2 <- modes[2]
	n3 <- modes[3]
	#progress bar
	pb <- txtProgressBar(min=0,max=n3,style=3)
	#define ifft
	ifft <- function(x){suppressWarnings(as.numeric(fft(x,inverse=TRUE))/length(x))}
	#fft for each of the n1n2 vectors (of length n3) along mode 3
	fftz <- aperm(apply(tnsr@data,MARGIN=1:2,fft),c(2,3,1))
	#svd for each face (svdz is a list of the results)
	U_arr <- array(0,dim=c(n1,n1,n3))
	V_arr <- array(0,dim=c(n2,n2,n3))
	m <- min(n1,n2)		
	S_arr <- array(0,dim=c(n1,n2,n3))
	#Think of a way to avoid a loop in the beginning
	#Problem is that svd returns a list but ideally we want 3 arrays
	#Even with unlist this doesn't seem possible
	for (j in 1:n3){
		setTxtProgressBar(pb,j)
		decomp <- svd(fftz[,,j],nu=n1,nv=n2)
		U_arr[,,j] <- decomp$u
		V_arr[,,j] <- decomp$v
		S_arr[,,j] <- diag(decomp$d,nrow=n1,ncol=n2) #length is min(n1,n2)
	}	
	close(pb)
	#for each svd result, we want to apply ifft
	U <- as.tensor(aperm(apply(U_arr,MARGIN=1:2,ifft),c(2,3,1)))
	V <- as.tensor(aperm(apply(V_arr,MARGIN=1:2,ifft),c(2,3,1)))
	S <- as.tensor(aperm(apply(S_arr,MARGIN=1:2,ifft),c(2,3,1)))
	
	est <- as.tensor(array(0,dim=modes))
	for (i in 1:k){
		est <- est + t_mult(t_mult(U[,i,,drop=FALSE],S[i,i,,drop=FALSE]),t(V[,i,,drop=FALSE]))
	}
	resid <- fnorm(est-tnsr)
	invisible(list(est=est, fnorm_resid = resid, norm_percent = (1-resid/fnorm(tnsr))*100))
}

###t-compress2 (Not Supported)
.t_compress2 <- function(tnsr,k1,k2){
	A = modeSum(tnsr,m=3,drop=TRUE)
	svdz <- svd(A@data,nu=k1,nv=k2)
	Util <- svdz$u
	Vtil <- svdz$v
	modes <- tnsr@modes
	n3 <- modes[3]
	core <- array(0,dim=c(k1,k2,n3))
	for(i in 1:n3){
	core[,,i]<-t(Util)%*%tnsr[,,i]@data%*%Vtil
	}
	est <- array(0,dim=modes)
	for(i in 1:k1){
		for (j in 1:k2){
			est = est + Util[,i] %o% Vtil[,j] %o% core[i,j,]
		}	
	}
	resid <- fnorm(tnsr - est)
	invisible(list(core = as.tensor(core), est=est, fnorm_resid = resid, norm_percent = (1-resid/fnorm(tnsr))*100))
}

#' sparse nonnegative Tucker decomposition
#'
#' Decomposes nonnegative tensor \code{tnsr} into core nonnegative tensor \code{Z} and sparse nonnegative factor matrices \code{U[n]}.
#'@export
#'@param tnsr nonnegative tensor with \code{K} modes
#'@param ranks an integer vector of length \code{K} specifying the modes sizes for the output core tensor \code{Z}
#'@param tol relative Frobenius norm error tolerance
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param max_time max running time
#'@param lambda \code{K+1} vector of sparsity regularizer coefficients for the factor matrices and the core tensor
#'@param L_min lower bound for Lipschitz constant for the gradients of residual error \eqn{l(Z,U) = fnorm(tnsr - ttl(Z, U))} by \code{Z} and each \code{U}
#'@param rw controls the extrapolation weight
#'@param bound upper bound for the elements of \code{Z} and \code{U[[n]]} (the ones that have zero regularization coefficient \code{lambda})
#'@param U0 initial factor matrices, defaults to nonnegative Gaussian random matrices
#'@param Z0 initial core tensor \code{Z}, defaults to nonnegative Gaussian random tensor
#'@param verbose more output algorithm progress
#'@param unfold_tnsr precalculate \code{tnsr} to matrix unfolding by every mode (speeds up calculation, but may require lots of memory)
#'@return a list:\describe{
#'\item{\code{U}}{nonnegative factor matrices}
#'\item{\code{Z}}{nonnegative core tensor}
#'\item{\code{est}}{estimate \eqn{Z \times_1 U_1 \ldots \times_K U_K}}
#'\item{\code{conv}}{method convergence indicator}
#'\item{\code{resid}}{the Frobenius norm of the residual error \code{l(Z,U)} plus regularization penalty (if any)}
#'\item{\code{n_iter}}{number of iterations}
#'\item{\code{n_redo}}{number of times Z and U were recalculated to avoid the increase in objective function}
#'\item{\code{diag}}{convergence info for each iteration\describe{
#'\item{\code{all_resids}}{residues}
#'\item{\code{all_rel_resid_deltas}}{residue delta relative to the current residue}
#'\item{\code{all_rel_resids}}{residue relative to the \code{sqrt(||tnsr||)}}
#'}}}
#'
#'@details The function uses the alternating proximal gradient method to solve the following optimization problem:
#' \deqn{\min 0.5 \|tnsr - Z \times_1 U_1 \ldots \times_K U_K \|_{F^2} + 
#' \sum_{n=1}^{K} \lambda_n \|U_n\|_1 + \lambda_{K+1} \|Z\|_1, \;\text{where}\; Z \geq 0, \, U_i \geq 0.}
#' The method stops if either the relative improvement of the error is below the tolerance \code{tol} for 3 consequitive iterations or
#' both the relative error improvement and relative error (wrt the \code{tnsr} norm) are below the tolerance.
#' Otherwise it stops if the maximal number of iterations or the time limit were reached.
#'
#'@note The implementation is based on ntds() MATLAB code by Yangyang Xu and Wotao Yin.
#'@references Y. Xu, "Alternating proximal gradient method for sparse nonnegative Tucker decomposition", Math. Prog. Comp., 7, 39-70, 2013.
#'@seealso \code{\link{tucker}}
#'@seealso \url{http://www.caam.rice.edu/~optimization/bcu/}
tucker.nonneg <- function( tnsr, ranks,
                           tol=1e-4, hosvd=FALSE,
                           max_iter = 500, max_time=0,
                           lambda = rep.int( 0, length(ranks)+1 ), L_min = 1, rw=0.9999,
                           bound = Inf,
                           U0=NULL, Z0=NULL, verbose=FALSE,
                           unfold_tnsr=length(dim(tnsr))*prod(dim(tnsr)) < 4000^2 )
{
  #progress bar
  start_time <- proc.time()
  pb <- txtProgressBar(min=0,max=max_iter,style=3)

  make_nonneg.tnsr <- function( tnsr )
  {
    tnsr@data[ tnsr@data < 0 ] <- 0
    return ( tnsr )
  }
  make_nonneg.mtx <- function( mtx )
  {
    mtx[ mtx < 0 ] <- 0
    return ( mtx )
  }
  # update core tensor
  # return new residual error
  makeZStep <- function( curZ )
  {
    gradZ <- ttl(curZ, Usq, seq_len(K)) - TtU
    # update core vector
    Z <<- curZ - gradZ/L[[K+1]]
    if ( lambda[[K+1]] > 0 ) {
      Z <<- Z - lambda[[K+1]]/L[[K+1]]
    }
    Z <<- make_nonneg.tnsr( Z )
    # do projection
    if ( doproj[[K+1]] ) {
      mask <- abs(Z@data) > bound
      Z@data[mask] <<- sign(Z@data[mask]) * bound
    }
    return ( invisible() )
  }

  # update n-th factor matrix (U[[n]])
  # return new residual error
  makeUnStep <- function( curU, n )
  {
    if ( !is.null(Tmtx) ) {
      B <- unfold( ttl( Z, U[-n], seq_len(K)[-n] ), n, seq_len(K)[-n] )
      Bsq <- tcrossprod(B@data)
      TB <- tcrossprod(Tmtx[[n]], B@data)
    } else {
      B <- unfold( ttl( Z, Usq[-n], seq_len(K)[-n] ), n, seq_len(K)[-n] )
      TB <- unfold( ttl( tnsr, U[-n], seq_len(K)[-n], transpose=TRUE ), n, seq_len(K)[-n] )
      Zn <- unfold( Z, n, seq_len(K)[-n] )

      Bsq <- tcrossprod( B@data, Zn@data )
      TB <- tcrossprod( TB@data, Zn@data )
    }
    # compute the gradient
    gradU <- curU %*% Bsq - TB
    # update Lipschitz constant
    L0[[n]] <<- L[[n]]
    L[[n]] <<- max( L_min, norm(Bsq, '2') )
    # update n-th factor matrix
    newU <- make_nonneg.mtx( curU - (gradU+lambda[[n]])/L[[n]] )
    if ( doproj[[n]] ) newU[ newU > bound ] <- bound

    # update U[[n]]
    U[[n]] <<- newU
    Usq[[n]] <<- crossprod( U[[n]] )
    nrmUsq[[n]] <<- norm( Usq[[n]], '2' )

    # --- diagnostics, reporting, stopping checks ---
    newResid <- 0.5*(sum(Usq[[n]]*Bsq)-2*sum(U[[n]]*TB)+Tnrm^2)
    if (sparse.reg) {
      newResid <- newResid + lambda %*% c( sapply( U, sum ), sum(abs(Z@data)) )
    }
    return ( newResid )
  }

  Kway <- dim(tnsr) # dimension of tnsr
  K <- length(Kway) # tnsr is an K-way tensor

  if ( is.null(U0) ) {
    if ( verbose ) message( 'Generating random initial factor matrices estimates...' )
    U0 <- lapply( seq_len(K), function(n) make_nonneg.mtx( matrix( rnorm( Kway[[n]]*ranks[[n]] ), ncol = ranks[[n]] ) ) )
  }
  if ( is.null(Z0) ) {
    if ( verbose ) message( 'Generating random initial core tensor estimate...' )
    Z0 <- make_nonneg.tnsr( rand_tensor( modes = ranks, drop = FALSE) )
  }

  # pre-process the starting point
  if (hosvd) {
    if ( verbose ) message( 'Applying High Order SVD to improve initial U and Z...' )
    # "solve" Z = tnsr x_1 U_1' ... x_K U_K'
    U0 <- lapply( seq_len(K), function(n) {
      U0n_tilde <- unfold( ttl(tnsr, U0[-n], seq_len(K)[-n], transpose=TRUE ),
                           row_idx = n, col_idx = seq_len(K)[-n] )@data
      U0n_vecs <- svd( U0n_tilde, nu = ranks[[n]], nv = 0 )$u
      U0n <- matrix( unlist( lapply( U0n_vecs, function( Uvec ) {
        # make the largest absolute element positive
        i <- which.max( abs(Uvec) )
        if ( Uvec[[i]] < 0 ) Uvec <- -Uvec
        # project to > 0
        Uvec <- pmax( .Machine$double.eps, Uvec )
      } ) ), ncol=ranks[[n]] )
      return ( U0n/sum(U0n) )
    } )
    Z0 <- ttl( tnsr, U0, seq_len(K), transpose=TRUE )
  }
  # check the existence of sparseness regularizer
  sparse.reg <- any(lambda>0)
  # add bound constraint for well-definedness
  doproj <- lambda == 0 & is.finite(bound)

  Tnrm <- fnorm(tnsr)

  # rescale the initial point according to the number of elements
  Knum <- Kway * ranks
  totalNum <- prod(ranks) + sum(Knum)
  U0 <- lapply( seq_along(U0), function(n) U0[[n]]/norm(U0[[n]],"F")*Tnrm^(Knum[[n]]/totalNum) )
  Usq0 <- lapply( U0, crossprod )
  nrmUsq <- sapply( Usq0, norm, '2' )
  Z0 <- Z0/fnorm(Z0)*Tnrm^(prod(ranks)/totalNum)

  resid0 <- 0.5*fnorm( tnsr-ttl(Z0,U0,seq_len(K)) )^2
  if (sparse.reg) resid0 <- resid0 + lambda %*% c( sapply( U0, sum ), sum(abs(Z0@data)) )
  resid <- resid0

  # Precompute matrix unfoldings of input tensor to save computing time if it is not too large
  if (unfold_tnsr) {
    if ( verbose ) message( 'Precomputing input tensor unfoldings...' )
    Tmtx <- lapply( seq_len(K), function(n) unfold( tnsr, row_idx = n, col_idx = seq_len(K)[-n] )@data )
  } else {
    if ( verbose ) message( 'No precomputing of tensor unfoldings' )
    Tmtx <- NULL
  }

  # Iterations of block-coordinate update
  # iteratively updated variables:
  # GradU: gradients with respect to each component matrix of U
  # GradZ: gradient with respect to Z
  # U,Z: new updates
  # U0,Z0: old updates
  # Um,Zm: extrapolations of U
  # L, L0: current and previous Lipschitz bounds
  # resid, resid0: current and previous residual error
  U <- U0
  Um <- U0
  Usq <- Usq0
  Z <- Z0
  Zm <- Z0

  t0 <- rep.int( 1, K+1 )
  t <- t0
  wU <- rep.int( 0, K+1 )
  L0 <- rep.int( 1, K+1 )
  L <- L0

  all_resids <- numeric(0)
  all_rel_resid_deltas <- numeric(0)
  all_rel_resids <- numeric(0)
  n_stall <- 0
  n_redo <- 0
  conv <- FALSE

  # do the iterations
  if ( verbose ) message( 'Starting iterations...' )
  for (n_iter in seq_len(max_iter)) {
    setTxtProgressBar(pb, n_iter)

    residn0 <- resid
    TtU0 <- list( ttm( tnsr, U0[[1]], 1, transpose=TRUE ) )
    for (n in 2:K) {
      TtU0[[n]] <- ttm( TtU0[[n-1]], U0[[n]], n, transpose=TRUE )
    }

    for (n in seq.int(from=K,to=1) ) {
      # -- update the core tensor Z --
      L0[[K+1]] <- L[[K+1]]
      L[[K+1]] <- pmax( L_min, prod( nrmUsq ) )

      # try to make a step using extrapolated decompositon (Zm,Um)
      TtU <- if ( n < K ) ttl( TtU0[[n]], U[(n+1):K], (n+1):K, transpose=TRUE ) else TtU0[[K]]
      makeZStep( Zm )
      residn <- makeUnStep( Um[[n]], n )
      if ( residn>residn0 ) {
        # extrapolated Zm,Um decomposition lead to residual norm increase,
        # revert extrapolation and make a step using Z0,U0 to ensure
        # objective function is decreased
        n_redo <- n_redo + 1
        # re-update to make objective nonincreasing
        Usq[[n]] <- Usq0[[n]] # Z update needs it
        makeZStep( Z0 )
        residn <- makeUnStep( U0[[n]], n )
        if ( residn>residn0 ) warning( n_iter, ': residue increase at redo step' )
      }
      # --- correction and extrapolation ---
      t[[n]] <- (1+sqrt(1+4*t0[[n]]^2))/2
      # choose smaller weight of U[[n]] for convergence
      wU[[n]] <- min( (t0[[n]]-1)/t[[n]], rw*sqrt(L0[[n]]/L[[n]]) )
      Um[[n]] <- U[[n]] + wU[[n]]*( U[[n]]-U0[[n]] )
      t[[K+1]] <- (1+sqrt(1+4*t0[[K+1]]^2))/2
      # choose smaller weight of Z for convergence
      wU[[K+1]] <- min( (t0[[K+1]]-1)/t[[K+1]], rw*sqrt(L0[[K+1]]/L[[K+1]]) )
      Zm <- Z + wU[K+1]*(Z-Z0)

      # store the current update
      Z0 <- Z
      U0[[n]] <- U[[n]]
      Usq0[[n]] <- Usq[[n]]
      t0[c(n,K+1)] <- t[c(n,K+1)]
      residn0 <- residn
    }

    # --- diagnostics, reporting, stopping checks ---
    resid0 <- resid
    resid <- max( 0, residn )

    rel_resid_delta <- abs(resid-resid0)/(resid0+1)
    rel_resid <- sqrt(2*resid)/Tnrm

    # reporting
    all_resids <- append( all_resids, resid )
    all_rel_resid_deltas <- append( all_rel_resid_deltas, rel_resid_delta )
    all_rel_resids <- append( all_rel_resids, rel_resid )

    # check stopping criterion
    crit <- rel_resid_delta < tol
    n_stall <- ifelse( crit, n_stall+1, 0 )
    if ( n_stall >= 3 || rel_resid < tol ) {
      if ( verbose ) {
        if ( rel_resid == 0 ) message( 'Residue is zero. Exact decomposition was found' )
        if ( n_stall >= 3 ) message( 'Residue relative delta below ', tol, ' ', n_stall, ' times in a row' )
        if ( rel_resid < tol ) message( 'Residue is ', rel_resid, ' times below input tensor norm' )
        message( 'tucker.nonneg() converged in ', n_iter, ' iteration(s), ', n_redo, ' redo steps' )
      }
      conv = TRUE
      break
    }
    if ( max_time > 0 && ( proc.time() - start_time )[[3]] > max_time ) {
      warning( "Maximal time exceeded, might be not an optimal solution")
      break
    }
  }
  setTxtProgressBar(pb,max_iter)
  close(pb)
  if ( !conv && n_iter == max_iter ) {
    warning( "Maximal number of iterations reached, might be not an optimal solution")
  }

  return ( invisible( list(U=U, Z=Z, est=ttl(Z,U,seq_len(K)),
                           n_iter = n_iter,
                           n_redo = n_redo,
                           conv = conv,
                           resid = resid,
                           diag = list( all_resids = all_resids,
                                        all_rel_resid_deltas = all_rel_resid_deltas,
                                        all_rel_resids = all_rel_resids ) ) ) )
}
