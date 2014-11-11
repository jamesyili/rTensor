pkgname <- "rTensor"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rTensor')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Ops-methods")
### * Ops-methods

flush(stderr()); flush(stdout())

### Name: Ops-methods
### Title: Conformable elementwise operators for Tensor
### Aliases: Ops,Tensor,Tensor-method Ops,Tensor,array-method
###   Ops,Tensor,numeric-method Ops,array,Tensor-method
###   Ops,numeric,Tensor-method Ops-methods

### ** Examples

tnsr <- rand_tensor(c(3,4,5))
tnsr2 <- rand_tensor(c(3,4,5))
tnsrsum <- tnsr + tnsr2
tnsrdiff <- tnsr - tnsr2
tnsrelemprod <- tnsr * tnsr2
tnsrelemquot <- tnsr / tnsr2
for (i in 1:3L){
	for (j in 1:4L){
		for (k in 1:5L){
			stopifnot(tnsrsum@data[i,j,k]==tnsr@data[i,j,k]+tnsr2@data[i,j,k])
			stopifnot(tnsrdiff@data[i,j,k]==(tnsr@data[i,j,k]-tnsr2@data[i,j,k]))
			stopifnot(tnsrelemprod@data[i,j,k]==tnsr@data[i,j,k]*tnsr2@data[i,j,k])
			stopifnot(tnsrelemquot@data[i,j,k]==tnsr@data[i,j,k]/tnsr2@data[i,j,k])
}
}
}



cleanEx()
nameEx("Tensor-class")
### * Tensor-class

flush(stderr()); flush(stdout())

### Name: Tensor-class
### Title: S4 Class for a Tensor
### Aliases: Tensor Tensor-class

### ** Examples

tnsr <- rand_tensor()
class(tnsr)
tnsr
print(tnsr)
dim(tnsr)
tnsr@num_modes
tnsr@data



cleanEx()
nameEx("as.tensor")
### * as.tensor

flush(stderr()); flush(stdout())

### Name: as.tensor
### Title: Tensor Conversion
### Aliases: as.tensor

### ** Examples

#From vector
vec <- runif(100); vecT <- as.tensor(vec); vecT
#From matrix
mat <- matrix(runif(1000),nrow=100,ncol=10)
matT <- as.tensor(mat); matT
#From array
indices <- c(10,20,30,40)
arr <- array(runif(prod(indices)), dim = indices)
arrT <- as.tensor(arr); arrT



cleanEx()
nameEx("cp")
### * cp

flush(stderr()); flush(stdout())

### Name: cp
### Title: Canonical Polyadic Decomposition
### Aliases: cp

### ** Examples

tnsr <- rand_tensor(c(6,7,8))
cpD <- cp(tnsr,num_components=5)
cpD$conv
cpD$norm_percent
plot(cpD$all_resids)



cleanEx()
nameEx("dim-methods")
### * dim-methods

flush(stderr()); flush(stdout())

### Name: dim-methods
### Title: Mode Getter for Tensor
### Aliases: dim,Tensor-method dim-methods

### ** Examples

tnsr <- rand_tensor()
dim(tnsr)



cleanEx()
nameEx("extract-methods")
### * extract-methods

flush(stderr()); flush(stdout())

### Name: [-methods
### Title: Extract or Replace Subtensors
### Aliases: [,Tensor-method [-methods [<-,Tensor-method
###   extract,Tensor-method

### ** Examples

tnsr <- rand_tensor()
tnsr[1,2,3]
tnsr[3,1,]
tnsr[,,5]
tnsr[,,5,drop=FALSE]

tnsr[1,2,3] <- 3; tnsr[1,2,3]
tnsr[3,1,] <- rep(0,5); tnsr[3,1,]
tnsr[,2,] <- matrix(0,nrow=3,ncol=5); tnsr[,2,]



cleanEx()
nameEx("fnorm-methods")
### * fnorm-methods

flush(stderr()); flush(stdout())

### Name: fnorm-methods
### Title: Tensor Frobenius Norm
### Aliases: fnorm fnorm,Tensor-method fnorm-methods

### ** Examples

tnsr <- rand_tensor()
fnorm(tnsr)



cleanEx()
nameEx("fold")
### * fold

flush(stderr()); flush(stdout())

### Name: fold
### Title: General Folding of Matrix
### Aliases: fold

### ** Examples

tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
matT3<-unfold(tnsr,row_idx=2,col_idx=c(3,1))
identical(fold(matT3,row_idx=2,col_idx=c(3,1),modes=c(3,4,5)),tnsr)



cleanEx()
nameEx("hamadard_list")
### * hamadard_list

flush(stderr()); flush(stdout())

### Name: hamadard_list
### Title: List Hamadard Product
### Aliases: hamadard_list

### ** Examples

lizt <- list('mat1' = matrix(runif(40),ncol=4),
'mat2' = matrix(runif(40),ncol=4),
'mat3' = matrix(runif(40),ncol=4))
dim(hamadard_list(lizt))



cleanEx()
nameEx("head-methods")
### * head-methods

flush(stderr()); flush(stdout())

### Name: head-methods
### Title: Head for Tensor
### Aliases: head,Tensor-method head-methods

### ** Examples

tnsr <- rand_tensor()
head(tnsr)



cleanEx()
nameEx("hosvd")
### * hosvd

flush(stderr()); flush(stdout())

### Name: hosvd
### Title: (Truncated-)Higher-order SVD
### Aliases: hosvd

### ** Examples

tnsr <- rand_tensor(c(6,7,8))
hosvdD <-hosvd(tnsr)
hosvdD$fnorm_resid
hosvdD2 <-hosvd(tnsr,ranks=c(3,3,4))
hosvdD2$fnorm_resid



cleanEx()
nameEx("innerProd-methods")
### * innerProd-methods

flush(stderr()); flush(stdout())

### Name: innerProd-methods
### Title: Tensors Inner Product
### Aliases: innerProd innerProd,Tensor,Tensor-method innerProd-methods

### ** Examples

tnsr1 <- rand_tensor()
tnsr2 <- rand_tensor()
innerProd(tnsr1,tnsr2)



cleanEx()
nameEx("k_fold")
### * k_fold

flush(stderr()); flush(stdout())

### Name: k_fold
### Title: k-mode Folding of Matrix
### Aliases: k_fold

### ** Examples

tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
matT2<-k_unfold(tnsr,m=2)
identical(k_fold(matT2,m=2,modes=c(3,4,5)),tnsr)



cleanEx()
nameEx("k_unfold-methods")
### * k_unfold-methods

flush(stderr()); flush(stdout())

### Name: k_unfold-methods
### Title: Tensor k-mode Unfolding
### Aliases: k_unfold k_unfold,Tensor-method k_unfold-methods

### ** Examples

tnsr <- rand_tensor()
matT2<-rs_unfold(tnsr,m=2)



cleanEx()
nameEx("khatri_rao")
### * khatri_rao

flush(stderr()); flush(stdout())

### Name: khatri_rao
### Title: Khatri-Rao Product
### Aliases: khatri_rao

### ** Examples

dim(khatri_rao(matrix(runif(12),ncol=4),matrix(runif(12),ncol=4)))



cleanEx()
nameEx("khatri_rao_list")
### * khatri_rao_list

flush(stderr()); flush(stdout())

### Name: khatri_rao_list
### Title: List Khatri-Rao Product
### Aliases: khatri_rao_list

### ** Examples

smalllizt <- list('mat1' = matrix(runif(12),ncol=4),
'mat2' = matrix(runif(12),ncol=4),
'mat3' = matrix(runif(12),ncol=4))
dim(khatri_rao_list(smalllizt))



cleanEx()
nameEx("kronecker_list")
### * kronecker_list

flush(stderr()); flush(stdout())

### Name: kronecker_list
### Title: List Kronecker Product
### Aliases: kronecker_list

### ** Examples

smalllizt <- list('mat1' = matrix(runif(12),ncol=4),
'mat2' = matrix(runif(12),ncol=4),
'mat3' = matrix(runif(12),ncol=4))
dim(kronecker_list(smalllizt))



cleanEx()
nameEx("matvec-methods")
### * matvec-methods

flush(stderr()); flush(stdout())

### Name: matvec-methods
### Title: Tensor Matvec Unfolding
### Aliases: matvec matvec,Tensor-method matvec-methods

### ** Examples

tnsr <- rand_tensor(c(2,3,4))
matT1<- matvec(tnsr)



cleanEx()
nameEx("modeMean-methods")
### * modeMean-methods

flush(stderr()); flush(stdout())

### Name: modeMean-methods
### Title: Tensor Mean Across Single Mode
### Aliases: modeMean modeMean,Tensor-method modeMean-methods

### ** Examples

tnsr <- rand_tensor()
modeMean(tnsr,1,drop=TRUE)



cleanEx()
nameEx("modeSum-methods")
### * modeSum-methods

flush(stderr()); flush(stdout())

### Name: modeSum-methods
### Title: Tensor Sum Across Single Mode
### Aliases: modeSum modeSum,Tensor-method modeSum-methods

### ** Examples

tnsr <- rand_tensor()
modeSum(tnsr,3,drop=TRUE)



cleanEx()
nameEx("mpca")
### * mpca

flush(stderr()); flush(stdout())

### Name: mpca
### Title: Multilinear Principal Components Analysis
### Aliases: mpca

### ** Examples

tnsr <-rand_tensor(c(100,10,10))
mpcaD <- mpca(tnsr,ranks=c(30,5))
mpcaD$conv
mpcaD$norm_percent
plot(mpcaD$all_resids)



cleanEx()
nameEx("print-methods")
### * print-methods

flush(stderr()); flush(stdout())

### Name: print-methods
### Title: Print for Tensor
### Aliases: print,Tensor-method print-methods

### ** Examples

tnsr <- rand_tensor()
print(tnsr)



cleanEx()
nameEx("pvd")
### * pvd

flush(stderr()); flush(stdout())

### Name: pvd
### Title: Population Value Decomposition
### Aliases: pvd

### ** Examples

tnsr <- rand_tensor(c(10,5,100))
pvdD<-pvd(tnsr,uranks=rep(8,100),wranks=rep(4,100),a=8,b=4)



cleanEx()
nameEx("rand_tensor")
### * rand_tensor

flush(stderr()); flush(stdout())

### Name: rand_tensor
### Title: Tensor with Random Entries
### Aliases: rand_tensor

### ** Examples

rand_tensor()
rand_tensor(c(4,4,4))
rand_tensor(c(10,2,1),TRUE)



cleanEx()
nameEx("show-methods")
### * show-methods

flush(stderr()); flush(stdout())

### Name: show-methods
### Title: Show for Tensor
### Aliases: show,Tensor-method show-methods

### ** Examples

tnsr <- rand_tensor()
tnsr



cleanEx()
nameEx("t-methods")
### * t-methods

flush(stderr()); flush(stdout())

### Name: t-methods
### Title: Tensor Transpose
### Aliases: t,Tensor-method t-methods

### ** Examples

tnsr <- rand_tensor()
identical(t(tnsr)@data[,,1],t(tnsr@data[,,1]))
identical(t(tnsr)@data[,,2],t(tnsr@data[,,5]))
identical(t(t(tnsr)),tnsr)



cleanEx()
nameEx("t_mult")
### * t_mult

flush(stderr()); flush(stdout())

### Name: t_mult
### Title: Tensor Multiplication (T-MULT)
### Aliases: t_mult

### ** Examples

tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
tnsr2 <- new("Tensor",3L,c(4L,3L,5L),data=runif(60))
t_mult(tnsr, tnsr2)



cleanEx()
nameEx("t_svd")
### * t_svd

flush(stderr()); flush(stdout())

### Name: t_svd
### Title: Tensor Singular Value Decomposition
### Aliases: t_svd

### ** Examples

tnsr <- rand_tensor()
tsvdD <- t_svd(tnsr)



cleanEx()
nameEx("t_svd_reconstruct")
### * t_svd_reconstruct

flush(stderr()); flush(stdout())

### Name: t_svd_reconstruct
### Title: Reconstruct Tensor From TSVD
### Aliases: t_svd_reconstruct

### ** Examples

tnsr <- rand_tensor(c(10,10,10))
tsvdD <- t_svd(tnsr)
1 - fnorm(t_svd_reconstruct(tsvdD)-tnsr)/fnorm(tnsr)



cleanEx()
nameEx("tail-methods")
### * tail-methods

flush(stderr()); flush(stdout())

### Name: tail-methods
### Title: Tail for Tensor
### Aliases: tail,Tensor-method tail-methods

### ** Examples

tnsr <- rand_tensor()
tail(tnsr)



cleanEx()
nameEx("tperm-methods")
### * tperm-methods

flush(stderr()); flush(stdout())

### Name: tperm-methods
### Title: Mode Permutation for Tensor
### Aliases: tperm tperm,Tensor-method tperm-methods

### ** Examples

tnsr <- rand_tensor(c(3,4,5))
dim(tperm(tnsr,perm=c(2,1,3)))
dim(tperm(tnsr,perm=c(1,3,2)))



cleanEx()
nameEx("ttl")
### * ttl

flush(stderr()); flush(stdout())

### Name: ttl
### Title: Tensor Times List
### Aliases: ttl

### ** Examples

tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
lizt <- list('mat1' = matrix(runif(30),ncol=3),
'mat2' = matrix(runif(40),ncol=4),
'mat3' = matrix(runif(50),ncol=5))
ttl(tnsr,lizt,ms=c(1,2,3))



cleanEx()
nameEx("ttm")
### * ttm

flush(stderr()); flush(stdout())

### Name: ttm
### Title: Tensor Times Matrix (m-Mode Product)
### Aliases: ttm

### ** Examples

tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
mat <- matrix(runif(50),ncol=5)
ttm(tnsr,mat,m=3)



cleanEx()
nameEx("tucker")
### * tucker

flush(stderr()); flush(stdout())

### Name: tucker
### Title: Tucker Decomposition
### Aliases: tucker

### ** Examples

tnsr <- rand_tensor(c(6,7,8))
tuckerD <- tucker(tnsr,ranks=c(3,3,4))
tuckerD$conv
tuckerD$norm_percent
plot(tuckerD$all_resids)



cleanEx()
nameEx("unfold-methods")
### * unfold-methods

flush(stderr()); flush(stdout())

### Name: unfold-methods
### Title: Tensor Unfolding
### Aliases: unfold unfold,Tensor-method unfold-methods

### ** Examples

tnsr <- rand_tensor()
matT3<-unfold(tnsr,row_idx=2,col_idx=c(3,1))



cleanEx()
nameEx("unmatvec")
### * unmatvec

flush(stderr()); flush(stdout())

### Name: unmatvec
### Title: Unmatvec Folding of Matrix
### Aliases: unmatvec

### ** Examples

tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
matT1<-matvec(tnsr)
identical(unmatvec(matT1,modes=c(3,4,5)),tnsr)



cleanEx()
nameEx("vec-methods")
### * vec-methods

flush(stderr()); flush(stdout())

### Name: vec-methods
### Title: Tensor Vec
### Aliases: vec vec,Tensor-method vec-methods

### ** Examples

tnsr <- rand_tensor(c(4,5,6,7))
vec(tnsr)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
