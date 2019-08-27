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

# cost and occupancy
nplcc_file <- here("data", "nplcc_cost_occupancy.zip")
if (!file.exists(nplcc_file)) {
  "https://s3.amazonaws.com/marxan-vs-ilp/nplcc_cost_occupancy.zip" %>% 
    download.file(destfile = nplcc_file)
}
cost_occ <- read_csv(nplcc_file, 
                     col_types = cols(.default = col_double(),
                                      pu = col_integer()))

# split out cost and occupancy
cost <- select(cost_occ, id = pu, cost) %>% 
  arrange(id)
occ <- select(cost_occ, -cost) %>% 
  gather("species_code", "amount", -pu) %>% 
  inner_join(species %>% select(species_code, species = id), 
             by = "species_code") %>% 
  select(pu, species, amount, name = species_code) %>% 
  arrange(pu, species)
rm(cost_occ)

# planning unit raster
pus <- here("data", "nplcc_planning-units.tif") %>% 
  raster() %>% 
  # create an empty template
  raster()

# setup runs ----

# define run matrix
marxan_runs <- expand.grid(
  marxan_iterations = 1e8,
  spf = 25
)
runs <- expand.grid(target = seq(0.1, 0.9, by = 0.1),
                    n_features = 72,
                    n_pu = 50625,
                    blm = c(0.1, 1, 10, 100, 1000)) %>%
  # add marxan specific parameters
  mutate(marxan = list(marxan_runs),
         run_id = 300 + row_number()) %>%
  select(run_id, everything())

# fixed run parameters
ilp_gap <- 0.001
marxan_reps <- 10
random_subset <- TRUE
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

# convert vector of selected units to raster using template
solution_to_raster <- function(x, y) {
  x <- filter(x, solution_1 == 1) %>% 
    pull(id)
  y[x] <- 1
  return(y)
}

set.seed(1)


tt <- here("data", "nplcc_planning-units.tif") %>% 
  raster()

tt[] <- 1:ncell(tt)

e <- extent(560000, 560000 + 22500, 5300000 - 22500, 5300000)
tmp.r <- crop(tt, e)

cost_ss <- cost[cost$id %in% tmp.r[], ] %>% 
  arrange(id)

bnd_mat <- boundary_matrix(tmp.r)
smm_mat <- summary(bnd_mat)
# df <- as.data.frame(as.matrix(bnd_mat))

bnd_df <- data.frame(id1 = tmp.r[][smm_mat$i],
                     id2 = tmp.r[][smm_mat$j],
                     amount = round(smm_mat$x,0))

runs <- foreach(run = seq_len(nrow(runs)), .combine = bind_rows) %do% {
  r <- runs[run, ]
  str_glue_data(r, "Run ", run, 
                ": Target {target}; Features {n_features}; PUs {n_pu}; BLM {blm}") %>% 
    message()
  
  
  
  # sample species and planning units
  features <- species %>% 
    select(id, name = species_code) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  
  r$species <- paste(features$name, collapse = ",")
  pu_ss <- pus
  pu_ss[cost_ss$id] <- 0
  rij <- filter(occ, species %in% features$id, pu %in% cost_ss$id) %>% 
    arrange(pu, species)
  targets <- group_by(rij, species) %>% 
    summarize(tot_amount = sum(amount, na.rm = TRUE)) %>% 
    ungroup() %>% 
    transmute(id = species, amount = r$target * tot_amount)
  features <- inner_join(features, targets, by = "id")
  
  # ilp 
  p <- problem(cost_ss, 
               features = features %>% select(id, name), 
               rij = rij, 
               cost_column = "cost") %>% 
    add_min_set_objective() %>%
    add_relative_targets(r$target) %>%
    add_binary_decisions() %>%
    add_boundary_penalties(penalty = r$blm, edge_factor = 0.5, data = bnd_mat)
  # gurobi
  s_gur <- p %>% 
    add_gurobi_solver(gap = ilp_gap) %>%
    prioritizr_timed(force = TRUE)
  # solution summary
  cost_gurobi <- attr(s_gur$result, "objective")
  r$gurobi <- list(tibble(n_solutions = 1,
                          cost = cost_gurobi, 
                          time = s_gur$time[["elapsed"]]))
  # save solution
  s_gur <- "gurobi_target-{target}_features-{n_features}_pu-{n_pu}_blm-{blm}.tif" %>% 
    str_glue_data(r, .) %>% 
    file.path(gurobi_dir, .) %>% 
    writeRaster(solution_to_raster(s_gur$result, pus), overwrite = TRUE, .)
  rm(s_gur)
  
  # symphony
  s_sym <- p %>% 
    add_rsymphony_solver(gap = ilp_gap * cost_gurobi) %>% 
    prioritizr_timed(force = TRUE)
  # solution summary
  r$rsymphony <- list(tibble(n_solutions = 1,
                             cost = attr(s_sym$result, "objective"), 
                             time = s_sym$time[["elapsed"]]))
  # save solution
  s_sym <- "rsymphony_target-{target}_features-{n_features}_pu-{n_pu}_blm-{blm}.tif" %>% 
    str_glue_data(r, .) %>% 
    file.path(rsymphony_dir, .) %>% 
    writeRaster(solution_to_raster(s_sym$result, pus), overwrite = TRUE, .)
  rm(s_sym)
  
  #marxan
  m_data <- MarxanData(pu = mutate(cost_ss, status = 0L),
                       species = features %>%
                         mutate(spf = 1) %>%
                         select(id, target = amount, spf, name),
                       puvspecies = select(rij, species, pu, amount),
                       boundary = bnd_df, skipchecks = TRUE)
  #loop over marxan iterations
  r_marxan <- foreach(i_marxan = seq_len(nrow(r$marxan[[1]])),
                      .combine = bind_rows, .packages = pkg_list) %dopar% {
                        r_marxan <- r$marxan[[1]][i_marxan, , drop = FALSE]
                        message(paste0("  Marxan: ",
                                       "SPF = ", r_marxan$spf, "; ",
                                       r_marxan$marxan_iterations, " iterations"))
                        # data
                        m_data@species$spf <- r_marxan$spf
                        # options
                        m_opts <- MarxanOpts(BLM = r$blm, NCORES = 1L, VERBOSITY = 3L)
                        m_opts@NUMREPS <- as.integer(marxan_reps)
                        m_opts@NUMITNS <- as.integer(r_marxan$marxan_iterations)
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



