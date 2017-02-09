library(mccbn)
library(doMC)

############################### INPUT OPTIONS ################################
L = 100                             # number of repetitions
p = seq(4, 12, 2)                   # number of events
lambda_s = 1                        # sampling rate
N = sapply(50*p, min, 1000)         # number of observations / genotypes
eps = 0.05
thrds = 4
hcbn_path = "/Users/susanap/Documents/software/ct-cbn-0.1.04b/"
datadir = "/Users/susanap/Documents/hivX/CBN/hcbn_sampling/testdata/"
mccbn_path = "/Users/susanap/Documents/software/MC-CBN"
###############################################################################

source(file.path(mccbn_path, "hcbn_functions.r"))

write_poset <- function(p, poset, filename, datadir) {
  
  outfile = file.path(datadir, paste(filename, ".poset", sep=""))
  #  write header 
  write(p, outfile)
  
  relations = which(poset == 1)
  i = relations %% p
  j = (relations %/% p) + 1
  write.table(cbind(i, j), outfile,row.names=FALSE, col.names=FALSE,
              append=TRUE)

  #  write footer 
  write(0, outfile, append=TRUE)
}

# Set seed for reproducibility
set.seed(47)
registerDoMC(thrds)
OMPthreads = paste("OMP_NUM_THREADS=", thrds, " ", sep="")

obs_llhood_mc           = matrix(0, nrow=length(p), ncol=L)
relative_abs_error_mc   = matrix(0, nrow=length(p), ncol=L)
runtime_mc              = matrix(0, nrow=length(p), ncol=L)
obs_llhood_hcbn         = matrix(0, nrow=length(p), ncol=L)
relative_abs_error_hcbn = matrix(0, nrow=length(p), ncol=L)
runtime_hcbn            = matrix(0, nrow=length(p), ncol=L)

for (i in 1:length(p)) {
  res = foreach(j = 1:L, .combine=rbind) %dopar% {
    poset = random_poset(p[i]) 
    lambdas = runif(p[i], 1/3*lambda_s, 3*lambda_s)
    simulated_obs = sample_genotypes(N[i], poset, sampling_param=lambda_s, 
                                     lambdas=lambdas, eps=eps)
    t0 <- Sys.time()
    ret = MCMC_hcbn(poset, simulated_obs$obs_events)
    run_t_mc = as.numeric(difftime(Sys.time(), t0))
    llhood_mc = obs_log_likelihood(simulated_obs$obs_events, poset, ret$lambdas,
                                lambda_s, ret$eps, L=10000)
    error_mc = mean(abs(ret$lambdas - lambdas))/mean(lambdas)
    
    # save observations and poset for h-cbn
    filename = paste("simulated_obs_n", N[i], "_p", p[i], "_j", j, sep="")
    write(c(N[i], p[i]+1), file.path(datadir, paste(filename, ".pat", sep="")))
    write.table(cbind(rep(1, N[i]), simulated_obs$obs_events), 
                file.path(datadir, paste(filename, ".pat", sep="")),
                row.names=FALSE, col.names=FALSE, append=TRUE)
    
    write_poset(p[i], poset, filename, datadir)
    
    if (!dir.exists(file.path(datadir, filename))) {
      dir.create(file.path(datadir, filename))
    }
    
    t0 <- Sys.time()
    system(paste(OMPthreads, hcbn_path, "h-cbn -f", datadir, filename, " -w > ", datadir,
                 filename, ".out.txt", sep=""))
    run_t_hcbn = as.numeric(difftime(Sys.time(), t0))
    lambdas_hcbn = read.csv(file.path(datadir, 
                                      paste(filename, ".lambda", sep="")))
    lambdas_hcbn = as.vector(t(lambdas_hcbn))
    
    res_hcbn = read.table(file.path(datadir, paste(filename, ".out.txt", sep="")),
                          skip=1)
    llhood_hcbn = as.numeric(res_hcbn[3])
    error_hcbn = mean(abs(lambdas_hcbn - lambdas))/mean(lambdas)
    
    return(c(llhood_mc, error_mc, run_t_mc, llhood_hcbn, error_hcbn, run_t_hcbn))
  }
  obs_llhood_mc[i, ] = res[, 1]
  relative_abs_error_mc[i, ] = res[, 2]
  runtime_mc[i, ] = res[, 3]
  obs_llhood_hcbn[i, ] = res[, 4]
  relative_abs_error_hcbn[i, ] = res[, 5]
  runtime_hcbn[i, ] = res[, 6]
}

save(obs_llhood_mc, file=file.path(datadir,"obs_llhood_mc.Rdata"))
save(relative_abs_error_mc, 
     file=file.path(datadir, "relative_abs_error_mc.Rdata"))
save(runtime_mc, file=file.path(datadir, "runtime_mc.Rdata"))

save(obs_llhood_hcbn, file=file.path(datadir,"obs_llhood_hcbn.Rdata"))
save(relative_abs_error_hcbn, 
     file=file.path(datadir, "relative_abs_error_hcbn.Rdata"))
save(runtime_hcbn, file=file.path(datadir, "runtime_hcbn.Rdata"))
