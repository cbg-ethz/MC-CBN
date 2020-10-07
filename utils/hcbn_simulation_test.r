library(mccbn)
library(doMC)

############################### INPUT OPTIONS ################################
K = 100                      # number of repetitions
p = seq(4, 12, 2)            # number of events
lambda_s = 1                 # sampling rate
N = sapply(50*p, min, 1000)  # number of observations / genotypes
eps = 0.05
thrds = 4
run_hcbn = FALSE             # indicate whether or not H-CBN should be executed
save_output = FALSE          # indicate whether or not save output

# The following options are NOT MANDATORY to specify/modify
# If you want to perform comparisons to h-cbn, specify where to find hcbn 
# executables (if so, option 'run_hcbn' should be set to TRUE)
hcbn_path = "/Users/susanap/Documents/software/ct-cbn-0.1.04b/"

# Specify the directory where output files are to be saved. It is needed, if 
# you want to perform comparisons to h-cbn (save .pat and .poset files) and if 
# you want to save output files (if so, option 'save_output' should be set to
# TRUE). If path doesn't exit, it will set to the working directory
datadir = "/Users/susanap/Documents/hivX/CBN/hcbn_sampling/testdata/"

# Specify where to find 'hcbn_functions.r' (directory where MC-CBN library is
# saved). If path doesn't exit, it will set to the working directory
mccbn_path = "/Users/susanap/Documents/software/MC-CBN"
###############################################################################

### Check that required paths are given and exist
if (!dir.exists(mccbn_path)) {
  cat("Specified directory doesn't exist. Setting 'mccbn_path' to working" ,
      " directory, \'", getwd(), "\'\n", sep="")
  mccbn_path = getwd()
} 
if (!dir.exists(hcbn_path) & run_hcbn) {
  stop("Path to H-CBN executables should be specified")
}
if (!dir.exists(datadir) & (run_hcbn | save_output)) {
  cat("Specified directory doesn't exist. Setting 'datadir' to working" ,
      " directory, \'", getwd(), "\'\n", sep="")
  datadir = getwd()
}

### Load functions 
if (file.exists(file.path(mccbn_path, "hcbn_functions.r"))) {
  cat("Loading functions from \'" , mccbn_path, "\'\n", sep="")
  source(file.path(mccbn_path, "hcbn_functions.r"))
} else {
  stop("Script 'hcbn_functions.r' cannot be found in \'", mccbn_path, "\'")
}

###############################################################################
### MAIN PROGRAM
###############################################################################
# Set seed for reproducibility
set.seed(47)
registerDoMC(thrds)
OMPthreads = paste("OMP_NUM_THREADS=", thrds, " ", sep="")

obs_llhood_mc           = matrix(0, nrow=length(p), ncol=K)
relative_abs_error_mc   = matrix(0, nrow=length(p), ncol=K)
runtime_mc              = matrix(0, nrow=length(p), ncol=K)
if (run_hcbn) {
  obs_llhood_hcbn         = matrix(0, nrow=length(p), ncol=K)
  relative_abs_error_hcbn = matrix(0, nrow=length(p), ncol=K)
  runtime_hcbn            = matrix(0, nrow=length(p), ncol=K)
}

for (i in 1:length(p)) {
  res = foreach(k = 1:K, .combine=rbind) %dopar% {
    poset = random_poset(p[i]) 
    lambdas = runif(p[i], 1/3*lambda_s, 3*lambda_s)
    simulated_obs = sample_genotypes(N[i], poset=poset, sampling_param=lambda_s, 
                                     lambdas=lambdas, eps=eps)
    t0 = Sys.time()
    ret = MCEM_hcbn(poset, simulated_obs$obs_events)
    run_t_mc = as.numeric(difftime(Sys.time(), t0, units='mins'))
    llhood_mc = obs_log_likelihood(simulated_obs$obs_events, poset, ret$lambdas, 
                                   lambda_s, ret$eps, L=10000)
    error_mc = mean(abs(ret$lambdas - lambdas))/mean(lambdas)
    
    if (run_hcbn) {
      # save observations and poset for h-cbn
      filename = paste("simulated_obs_n", N[i], "_p", p[i], "_k", k, sep="")
      write(c(N[i], p[i] + 1), 
            file.path(datadir, paste(filename, ".pat", sep="")))
      write.table(cbind(rep(1, N[i]), simulated_obs$obs_events), 
                  file.path(datadir, paste(filename, ".pat", sep="")),
                  row.names=FALSE, col.names=FALSE, append=TRUE)
      
      write_poset(p[i], poset, filename, datadir)
      
      if (!dir.exists(file.path(datadir, filename))) {
        dir.create(file.path(datadir, filename))
      }
      
      t0 = Sys.time()
      system(paste(OMPthreads, hcbn_path, "h-cbn -f", datadir, filename,
                   " -w > ", datadir, filename, ".out.txt", sep=""))
      run_t_hcbn = as.numeric(difftime(Sys.time(), t0, units='mins'))
      lambdas_hcbn = read.csv(file.path(datadir, 
                                        paste(filename, ".lambda", sep="")))
      lambdas_hcbn = as.vector(t(lambdas_hcbn))
      
      res_hcbn = read.table(file.path(datadir, paste(filename, ".out.txt", 
                                                     sep="")), skip=1)
      llhood_hcbn = as.numeric(res_hcbn[3])
      error_hcbn = mean(abs(lambdas_hcbn - lambdas))/mean(lambdas)
      
      return(c(llhood_mc, error_mc, run_t_mc, llhood_hcbn, error_hcbn, 
               run_t_hcbn))
    } else {
      return(c(llhood_mc, error_mc, run_t_mc))
    }
  }
  obs_llhood_mc[i, ]         = res[, 1]
  relative_abs_error_mc[i, ] = res[, 2]
  runtime_mc[i, ]            = res[, 3]
  if (run_hcbn) {
    obs_llhood_hcbn[i, ]         = res[, 4]
    relative_abs_error_hcbn[i, ] = res[, 5]
    runtime_hcbn[i, ]            = res[, 6]
  }
}

if (save_output) {
  save(obs_llhood_mc, file=file.path(datadir,"obs_llhood_mc.Rdata"))
  save(relative_abs_error_mc, 
       file=file.path(datadir, "relative_abs_error_mc.Rdata"))
  save(runtime_mc, file=file.path(datadir, "runtime_mc.Rdata"))
  
  if (run_hcbn) {
    save(obs_llhood_hcbn, file=file.path(datadir,"obs_llhood_hcbn.Rdata"))
    save(relative_abs_error_hcbn, 
         file=file.path(datadir, "relative_abs_error_hcbn.Rdata"))
    save(runtime_hcbn, file=file.path(datadir, "runtime_hcbn.Rdata"))
  }
}
