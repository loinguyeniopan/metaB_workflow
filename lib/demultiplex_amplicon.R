# demultiplex amplicon libraries

install.packages("remotes")
remotes::install_github("trtcrd/fastqUtils")

library(fastqUtils)

demultRun(params = "run_params.csv", t2s = "t2s.csv",
          accountMixedOrientation = F, cores = 4)
