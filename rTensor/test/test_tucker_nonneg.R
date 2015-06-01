# An example of nonnegative and semi-nonnegative Tucker decomposition

require( rTensor )
## Generate nonnegative synthetic 3-order tensor

T.dim <- c( 40, 50, 60 )
# randomly generate core tensor
Z.dim <- c( 4, 5, 6 )
Z <- rand_tensor( Z.dim )
Z@data <- pmax( Z@data, 0 )

# randomly generate factor matrices
U <- lapply( seq_along(T.dim), function( mode_ix ) {
  matrix( pmax( 0, rnorm( T.dim[[mode_ix]]*Z.dim[[mode_ix]] ) ),
          ncol = Z.dim[[mode_ix]] )
} )

# generate tensor
T.orig <- ttl( Z, U, seq_along(Z.dim) )
Z <- Z/max(T.orig@data)
T.orig@data <- T.orig@data / max(T.orig@data)

sn.ratio <- 0.6
# -- add noise --
T.noise <- rand_tensor( T.dim )
T <- T.orig + 10^(-sn.ratio/0.2)*fnorm(T.orig)/fnorm(T.noise)*T.noise;

# Solve the problem
T.tucker_nonneg <- tucker.nonneg( T, Z.dim, tol = 1E-4, hosvd = TRUE,
                                  max_iter = 1000, verbose = TRUE, lambda = rep.int(0.1,4) )

# Reporting
rel_err <- fnorm(T.tucker_nonneg$est - T.orig)/fnorm(T.orig)
message('Relative error of decomposition = ', rel_err)

# test semi-nonnegative Tucker decomposition

Z.SNTD <- rand_tensor( Z.dim )

# generate tensor
T.SNTD_orig <- ttl( Z.SNTD, U, seq_along(Z.dim) )
T.SNTD_orig@data <- T.SNTD_orig@data / max(T.SNTD_orig@data)

# -- add noise --
T.SNTD <- T.SNTD_orig + 10^(-sn.ratio/0.2)*fnorm(T.SNTD_orig)/fnorm(T.noise)*T.noise;

# Solve the problem
T.SNTD_tucker_nonneg <- tucker.nonneg( T.SNTD, Z.dim, tol = 1E-4, hosvd = TRUE, core_nonneg = FALSE,
                                  max_iter = 1000, verbose = TRUE, lambda = rep.int(0.5,4) )

# Reporting
rel_err.SNTD <- fnorm(T.SNTD_tucker_nonneg$est - T.SNTD_orig)/fnorm(T.SNTD_orig)
message('Relative error of decomposition = ', rel_err.SNTD)
