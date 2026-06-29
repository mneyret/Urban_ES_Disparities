# config.R
path_data_raw   <- "/Users/margotneyret/Desktop/Research/Data/"
path_intermdediate_results    <- "/Users/margotneyret/Desktop/Research/Data /Users/margotneyret/Desktop/Research/Projets_intermediaire/EquitES_UrbanES/Intermediate\ data"
path_figures    <- "/Users/marie/local-research/cool-project/figures"

if (nchar(path_data_raw) == 0)
  stop("path_data_raw is not set or does not exist — copy config.R.template to config.R and fill in your local paths")

if (nchar(path_intermdediate_results) == 0)
  stop("path_intermdediate_results is not set or does not exist — copy config.R.template to config.R and fill in your local paths")

if (nchar(path_figures) == 0)
  stop("path_figures is not set or does not exist — copy config.R.template to config.R and fill in your local paths")