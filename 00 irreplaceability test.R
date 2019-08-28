library(raster)
library(prioritizr)
library(marxan)
library(foreach)
library(doParallel)
library(uuid)
library(here)
library(tidyverse)

pkg_list <- c("raster", "prioritizr", "marxan", "uuid",  "here", "tidyverse")
select <- dplyr::select
walk(list.files("R", full.names = TRUE), source)
prioritizr_timed <- add_timer(prioritizr::solve)
# parallelization
n_cores <- 12
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# load nplcc data ----

# species list
species <- here("data", "nplcc_species.csv") %>% 
  read_csv() %>% 
  mutate(id = as.integer(id))

# cost and features
cost_r <- raster(here("data/cost.tif"))
feat_st <- stack(list.files(here("data/features/"), full.names = TRUE))

# setup runs ----

# define run matrix
marxan_runs <- expand.grid(
  marxan_iterations = 1e8,
  spf = 25
)
runs <- expand.grid(target = seq(0.1, 0.9, by = 0.1),
                    n_features = 72,
                    n_pu = 50625,
                    blm = 0.0) %>%
  # add marxan specific parameters
  mutate(marxan = list(marxan_runs),
         run_id = 300 + row_number()) %>%
  select(run_id, everything())

# fixed run parameters
ilp_gap <- 0.001
marxan_reps <- 100
random_subset <- FALSE
sysname <- tolower(Sys.info()[["sysname"]])
marxan_path <- switch(sysname, 
                      windows = here("marxan", "Marxan_x64.exe"), 
                      darwin = here("marxan", "MarOpt_v243_Mac64"), 
                      linux = here("marxan", "MarOpt_v243_Linux64")
)
stopifnot(file.exists(marxan_path))

# iterate over runs ----

# clean up old files
gurobi_dir <- here("output_blm", "gurobi")
#unlink(gurobi_dir, recursive = TRUE)
dir.create(gurobi_dir)
rsymphony_dir <- here("output_blm", "rsymphony")
#unlink(rsymphony_dir, recursive = TRUE)
dir.create(rsymphony_dir)
marxan_dir <- here("output_blm", "marxan")
#unlink(marxan_dir, recursive = TRUE)
dir.create(marxan_dir)
runs_dir <- here("output_blm", "runs")
#unlink(runs_dir, recursive = TRUE)
dir.create(runs_dir)

set.seed(1)

e <- extent(560000, 560000 + 22500, 5300000 - 22500, 5300000)
cost_crop <- crop(cost_r, e)
feat_crop <- crop(feat_st, e)

rij <- rij_matrix(cost_crop, feat_crop)

bnd_mat <- boundary_matrix(cost_crop)
smm_mat <- summary(bnd_mat)
# df <- as.data.frame(as.matrix(bnd_mat))

bnd_df <- data.frame(id1 = smm_mat$i,
                     id2 = smm_mat$j,
                     amount = round(smm_mat$x,0))

run <- 3

r <- runs[run, ]
str_glue_data(r, "Run ", run, 
              ": Target {target}; Features {n_features}; PUs {n_pu}; BLM {blm}") %>% 
  message()

# ilp 
p <- problem(cost_crop, 
             feat_crop) %>% 
  add_min_set_objective() %>%
  add_relative_targets(r$target) %>%
  add_proportion_decisions() 
# add_binary_decisions() 
# gurobi
s_gur <- p %>% 
  add_gurobi_solver(gap = ilp_gap) %>%
  prioritizr_timed(force = TRUE)

spplot(s_gur$result, "cost", main = "Solution",  at = c(0, 0.99),
       col.regions = c("#440154", "#FDE725"))

system.time(rc <- replacement_cost(p, s_gur$result, force = TRUE, threads = n_cores))

# rw <- rarity_weighted_richness(p, s_gur$result)
# plot(rw)
