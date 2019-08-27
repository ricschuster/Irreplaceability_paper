nms <- names(cost_occ)

for(ii in 1:ncol(cost_occ)){
  tmp_r <- here("data", "nplcc_planning-units.tif") %>% 
    raster()
  vv <- tmp_r[]
  vv[!is.na(vv)] <- as.vector(unlist(cost_occ[,ii]))
  tmp_r[] <- vv
  writeRaster(tmp_r, file = paste0(here("data/raster//"), nms[ii], ".tif"), overwrite = TRUE)
}

# load packages
library(prioritizrdata)
library(prioritizr)

# load planning unit data
data(tas_pu)

# load conservation feature data
data(tas_features)

spplot(tas_pu, zcol = "cost", names.attr = "Cost",
       main = "Planning unit costs")

# build problem
p1 <- problem(tas_pu, tas_features, cost_column = "cost", run_checks = FALSE) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.3) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions() %>%
  add_boundary_penalties(penalty = 0.0005, edge_factor = 1)

# print the problem
s1 <- solve(p1)

spplot(s1, "solution_1", col.regions = c("grey90", "darkgreen"),
       main = "problem solution")

rc <- replacement_cost(p1, s1[,"solution_1"])

rc$rc[rc$rc > 100] <- 1.09


spplot(rc, "rc", main = "Irreplaceability",  at = c(seq(0, 0.9, 0.1), 1.01, 1.1),
       col.regions = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
                       "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725",
                       "#FF0000"))

rw <- rarity_weighted_richness(p1, s1[,"solution_1"])

spplot(rw, "rwr", main = "Irreplaceability",  at = c(seq(0, 0.9, 0.1), 1.01, 1.1),
       col.regions = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
                       "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725",
                       "#FF0000"))
