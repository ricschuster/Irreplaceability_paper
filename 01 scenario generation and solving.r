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
gurobi_dir <- here("output", "gurobi")
# unlink(gurobi_dir, recursive = TRUE)
dir.create(gurobi_dir)
marxan_dir <- here("output", "marxan")
# unlink(marxan_dir, recursive = TRUE)
dir.create(marxan_dir)
runs_dir <- here("output", "runs")
# unlink(runs_dir, recursive = TRUE)
dir.create(runs_dir)

set.seed(1)

# e <- extent(560000, 560000 + 22500, 5300000 - 22500, 5300000)
e <- extent(560000, 560000 + 2000, 5300000 - 2000, 5300000)
cost_crop <- crop(cost_r, e)
feat_crop <- crop(feat_st, e)

rij <- rij_matrix(cost_crop, feat_crop)

bnd_mat <- boundary_matrix(cost_crop)
smm_mat <- summary(bnd_mat)
# df <- as.data.frame(as.matrix(bnd_mat))

bnd_df <- data.frame(id1 = smm_mat$i,
                     id2 = smm_mat$j,
                     amount = round(smm_mat$x,0))


runs <- foreach(run = seq_len(nrow(runs)), .combine = bind_rows) %do% {
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
  
  rc <- replacement_cost(p, s_gur$result, force = TRUE, threads = n_cores)
  
  # spplot(rc, "cost", main = "Replacement cost",  at = c(seq(0, 0.9, 0.1), 1.01, 1.1),
  #        col.regions = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
  #                        "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725",
  #                        "#FF0000"))
  
  rw <- rarity_weighted_richness(p, s_gur$result)
  # plot(rw)
  # spplot(rw, "cost", main = "Rarity weigthed richness",  at = c(seq(0, 0.9, 0.1), 1.01, 1.1),
  #        col.regions = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
  #                        "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725",
  #                        "#FF0000"))
  
  # solution summary
  cost_gurobi <- attr(s_gur$result, "objective")
  r$gurobi <- list(tibble(n_solutions = 1,
                          cost = cost_gurobi, 
                          time = s_gur$time[["elapsed"]]))
  # save solution
  s_gur <- "gurobi_target-{target}_solution.tif" %>% 
    str_glue_data(r, .) %>% 
    file.path(gurobi_dir, .) %>% 
    writeRaster(s_gur$result, overwrite = TRUE, .)
  
  
  rc <- "gurobi_target-{target}_replacement_cost.tif" %>% 
    str_glue_data(r, .) %>% 
    file.path(gurobi_dir, .) %>% 
    writeRaster(rc, overwrite = TRUE, .)
  
  rw <- "gurobi_target-{target}_rarity_weighted_richness.tif" %>% 
    str_glue_data(r, .) %>% 
    file.path(gurobi_dir, .) %>% 
    writeRaster(rw, overwrite = TRUE, .)
  
  rm(s_gur, rc, rw)
  
                           
  #marxan
  pu <- data.frame(id = 1:ncell(cost_crop),
                   cost = cost_crop[],
                   status = 0L)
  
  spec <- data.frame(id = 1:nlayers(feat_crop),
                     target = r$target,
                     spf = r$marxan[[1]]$spf, 
                     name = names(feat_crop),
                     stringsAsFactors = FALSE)
  puvspecies <- data.frame(species = rep(1:nlayers(feat_crop),1, each = ncell(cost_crop)),
                           pu = rep(1:ncell(cost_crop), nlayers(feat_crop)),
                           amount = as.numeric(unlist(as.data.frame(feat_crop))),
                           stringsAsFactors = FALSE)
  puvspecies <- puvspecies[order(puvspecies$pu, puvspecies$species),]
  puvspecies <- puvspecies[puvspecies$amount > 0, ]
  
  m_data <- MarxanData(pu = pu,
                       species = spec,
                       puvspecies = puvspecies, 
                       boundary = bnd_df, skipchecks = TRUE)
  
  # options
  m_opts <- MarxanOpts(BLM = 0.0, NCORES = 1L, VERBOSITY = 3L)
  m_opts@NUMREPS <- as.integer(marxan_reps)
  m_opts@NUMITNS <- as.integer(r$marxan[[1]]$marxan_iterations)
  m_opts@NUMTEMP <- as.integer(ceiling(m_opts@NUMITNS * 0.2))
  m_unsolved <- MarxanUnsolved(opts = m_opts, data = m_data)
  # solve
  td <- file.path(tempdir(), paste0("marxan-run_", UUIDgenerate()))
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  write.MarxanUnsolved(m_unsolved, td)
  file.copy(marxan_path, td)
  system(paste0("chmod +x ", file.path(td, basename(marxan_path))))
  setwd(td)
  m_time <- system.time(system2(marxan_path, args = c("input.dat", "-s")))
  m_results <- safely(read.MarxanResults)(td)$result
  setwd(here())
  unlink(td, recursive = TRUE)
  # save
  if (!is.null(m_results)) {
    str_glue_data(cbind(r, r_marxan),
                  "marxan_target-{target}_features-{n_features}_pu-{n_pu}_blm-{blm}_",
                  "spf-{spf}_iters-{marxan_iterations}.csv") %>%
      file.path(marxan_dir, .) %>%
      write_csv(m_results@summary, .)
    r_marxan$n_solutions <- sum(m_results@summary$Shortfall == 0)
    if (r_marxan$n_solutions > 0) {
      best <- filter(m_results@summary, Shortfall == 0) %>%
        arrange(Score) %>%
        slice(1)
      r_marxan$cost <- best$Score
      # raster solution
      s_mar <- m_results@selections[best$Run_Number,] %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        setNames(c("id", "solution_1")) %>%
        mutate(id = str_replace(id, "^P", "") %>% as.numeric()) %>%
        solution_to_raster(pus)
      str_glue_data(cbind(r, r_marxan),
                    "marxan_target-{target}_features-{n_features}_pu-{n_pu}_blm-{blm}_",
                    "spf-{spf}_iters-{marxan_iterations}.tif") %>%
        file.path(marxan_dir, .) %>%
        writeRaster(s_mar, overwrite = TRUE, .)
      rm(s_mar)
    } else {
      r_marxan$cost <- NA_real_
    }
    r_marxan$time <- m_time["elapsed"]
  } else {
    r_marxan$n_solutions <- NA_integer_
    r_marxan$cost <- NA_real_
    r_marxan$time <- NA_real_
  }
  select(r_marxan, marxan_iterations, spf, n_solutions, cost, time)
}
  r$marxan <- list(r_marxan)
  # save this iteration in case of crashing
  str_glue_data(r, "run-", run,
                "_target-{target}_features-{n_features}_pu-{n_pu}_blm-{blm}.rds") %>%
    file.path(runs_dir, .) %>%
    saveRDS(r, .)
  r
}

# unnest
runs_g <- runs %>% 
  mutate(solver = "gurobi") %>% 
  select(run_id, solver, target, n_features, n_pu, species, gurobi) %>% 
  unnest()
runs_s <- runs %>% 
  mutate(solver = "rsymphony") %>% 
  select(run_id, solver, target, n_features, n_pu, species, rsymphony) %>% 
  unnest()
# runs_m <- runs %>% 
#   mutate(solver = "marxan") %>% 
#   select(run_id, solver, target, n_features, n_pu, species, marxan) %>% 
#   unnest()
runs_long <- bind_rows(runs_g, runs_s)
# runs_long <- bind_rows(runs_g, runs_s, runs_m)
write_csv(runs_long, here("output_blm", "ilp-comparison-runs_no_marx.csv"))

# clean up
stopCluster(cl)



