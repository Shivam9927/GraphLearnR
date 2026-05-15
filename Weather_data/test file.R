library(dplyr)
library(GraphLearnR)
library(tidyr)
data <- read.csv("4066002.csv")
station_counts <- data %>%
  group_by(STATION) %>%
  summarise(n_obs = sum(!is.na(TMAX))) %>%
  arrange(desc(n_obs))

head(station_counts, 30)
top20_stations <- station_counts %>%
  slice(1:20) %>%
  pull(STATION)

top20_stations

noaa20 <- data %>%
  filter(STATION %in% top20_stations) %>%
  select(STATION, NAME, LATITUDE, LONGITUDE, ELEVATION, DATE, TMAX, TMIN) %>%
  arrange(STATION, DATE)

str(noaa20)
head(noaa20)

TMAX_wide <- noaa20 %>%
  select(STATION, DATE, TMAX) %>%
  pivot_wider(names_from = DATE, values_from = TMAX)

TMIN_wide <- noaa20 %>%
  select(STATION, DATE, TMIN) %>%
  pivot_wider(names_from = DATE, values_from = TMIN)

sum(is.na(TMAX_wide))
sum(is.na(TMIN_wide))

# remove station column temporarily
TMAX_df <- TMAX_wide[, -1]
TMIN_df <- TMIN_wide[, -1]

# combine for joint imputation
combined_df <- cbind(TMAX_df, TMIN_df)


elevation <- noaa20 %>%
  distinct(STATION, ELEVATION) %>%
  arrange(STATION) %>%
  pull(ELEVATION)

N <- nrow(X_noisy)
n <- ncol(X_noisy)

P1 <- matrix(rep(elevation, n), nrow = N)
P2 <- TMIN_mat
P3 <- matrix(1, nrow = N, ncol = n)

P <- list(P1, P2, P3)


noaa20$DATE <- as.Date(noaa20$DATE)

dates <- sort(unique(noaa20$DATE))
selected_dates <- dates[1:90]

noaa90 <- noaa20 %>%
  filter(DATE %in% selected_dates)

length(unique(noaa90$DATE))

########################################## Tmax signal##################
TMAX_wide <- noaa90 %>%
  select(STATION, DATE, TMAX) %>%
  pivot_wider(names_from = DATE, values_from = TMAX) %>%
  arrange(STATION)
################################Pivot TMIN (regressor)#########
TMIN_wide <- noaa90 %>%
  select(STATION, DATE, TMIN) %>%
  pivot_wider(names_from = DATE, values_from = TMIN) %>%
  arrange(STATION)

###############################
X_noisy <- as.matrix(TMAX_wide[, -1])
TMIN_mat <- as.matrix(TMIN_wide[, -1])

elevation <- noaa90 %>%
  distinct(STATION, ELEVATION) %>%
  arrange(STATION) %>%
  pull(ELEVATION)

N <- nrow(X_noisy)
n <- ncol(X_noisy)

P1 <- matrix(rep(elevation, n), nrow = N)   # elevation
P2 <- TMIN_mat                             # TMIN
P3 <- matrix(1, nrow = N, ncol = n)        # intercept

P <- list(P1, P2, P3)

dim(X_noisy)
dim(P1)
dim(P2)
dim(P3)
dim(X_noisy)

fit <- gsp_regression(
  X_noisy = X_noisy,
  P = P,
  w1 = 0.001,
  w2 = 0.99,
  max_iter = 20,
  threshold = 0.005
)


png("NOAA_graph.png", width = 800, height = 800)
plot_laplacian(fit$L, main = "NOAA Learned Graph")
dev.off()

fit_a <- gsp_regression(X_noisy, P, w1 = 0.01, w2 = 2, max_iter = 20, threshold = 0.005)
fit_b <- gsp_regression(X_noisy, P, w1 = 0.02, w2 = 3, max_iter = 20, threshold = 0.005)
fit_c <- gsp_regression(X_noisy, P, w1 = 0.05, w2 = 5, max_iter = 20, threshold = 0.005)

sum(fit_a$L[lower.tri(fit_a$L)] != 0)
sum(fit_b$L[lower.tri(fit_b$L)] != 0)
sum(fit_c$L[lower.tri(fit_c$L)] != 0)


fit_final <- gsp_regression(
  X_noisy = X_noisy,
  P = P,
  w1 = 0.02,
  w2 = 3,
  max_iter = 20,
  threshold = 0.005
)

png("NOAA_graph_circular.png", width = 800, height = 800)

par(mar = c(2,2,2,2))
plot_laplacian(fit_final$L, main = "NOAA Learned Graph (Circular)")

dev.off()
png("NOAA_graph_map.png", width = 900, height = 800)

coords <- noaa90 %>%
  distinct(STATION, LATITUDE, LONGITUDE) %>%
  arrange(STATION)

plot(coords$LONGITUDE, coords$LATITUDE,
     xlab = "Longitude",
     ylab = "Latitude",
     main = "NOAA Learned Graph (Map)",
     pch = 19)

L <- fit_final$L

for (i in 1:(nrow(L)-1)) {
  for (j in (i+1):nrow(L)) {
    
    if (abs(L[i, j]) > 0.01) {
      
      segments(
        coords$LONGITUDE[i], coords$LATITUDE[i],
        coords$LONGITUDE[j], coords$LATITUDE[j],
        col = "blue"
      )
    }
  }
}

dev.off()

station_info <- noaa90 %>%
  distinct(STATION, NAME, LATITUDE, LONGITUDE, ELEVATION) %>%
  arrange(STATION)

adj_binary <- (fit_final$L != 0) * 1
diag(adj_binary) <- 0

degree_counts <- rowSums(adj_binary)

connectivity_df <- data.frame(
  Node = 1:nrow(fit_final$L),
  STATION = station_info$STATION,
  NAME = station_info$NAME,
  LATITUDE = station_info$LATITUDE,
  LONGITUDE = station_info$LONGITUDE,
  ELEVATION = station_info$ELEVATION,
  Degree = degree_counts
)

connectivity_df <- connectivity_df[order(-connectivity_df$Degree), ]
print(connectivity_df)

png("NOAA_node_connectivity.png", width = 1000, height = 700)

barplot(
  connectivity_df$Degree,
  names.arg = connectivity_df$Node,
  las = 2,
  main = "NOAA Node Connectivity (Degree)",
  xlab = "Node",
  ylab = "Number of Connections"
)

dev.off()

library(igraph)

# Build adjacency from learned Laplacian
adj_mat <- abs(fit_final$L)
diag(adj_mat) <- 0
adj_mat[adj_mat > 0] <- abs(fit_final$L[adj_mat > 0])

# Create graph object
g_noaa <- graph_from_adjacency_matrix(
  adj_mat,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

# Degree for node connectivity
deg <- degree(g_noaa)

# Save figure
png("NOAA_node_connectivity_network.png", width = 900, height = 800)

plot(
  g_noaa,
  layout = layout_in_circle(g_noaa),
  vertex.size = 6 + deg * 4,
  vertex.label = NA,
  vertex.color = "steelblue",
  vertex.frame.color = NA,
  edge.width = 1 + E(g_noaa)$weight * 8,
  edge.color = "gray30",
  main = "Estimated Graph from Real Data"
)

dev.off()


library(igraph)

# 1. Build weighted adjacency from learned Laplacian
adj_mat <- abs(fit_final$L)
diag(adj_mat) <- 0

# remove tiny numerical noise
adj_mat[adj_mat < 1e-10] <- 0

# 2. Create graph
g_noaa <- graph_from_adjacency_matrix(
  adj_mat,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

# 3. Degree
deg <- degree(g_noaa)

# 4. Edge widths
edge_w <- E(g_noaa)$weight
if (is.null(edge_w) || length(edge_w) == 0) {
  edge_w <- rep(1, gsize(g_noaa))
} else {
  edge_w <- 1 + 6 * edge_w / max(edge_w)
}

# 5. Save figure
png("NOAA_node_connectivity_network.png", width = 900, height = 800)

plot(
  g_noaa,
  layout = layout_in_circle(g_noaa),
  vertex.size = 8 + 4 * deg,
  vertex.label = 1:vcount(g_noaa),
  vertex.label.cex = 0.8,
  vertex.color = "lightblue",
  vertex.frame.color = "navy",
  edge.width = edge_w,
  edge.color = "gray30",
  main = "Estimated Graph from Real Data"
)

dev.off()

station_info <- noaa90 %>%
  distinct(STATION, NAME, LATITUDE, LONGITUDE, ELEVATION) %>%
  arrange(STATION)

station_info$Node <- 1:nrow(station_info)

station_info

coords_mat <- as.matrix(station_info[, c("LONGITUDE", "LATITUDE")])

dist_mat_m <- geosphere::distm(coords_mat, fun = geosphere::distGeo)

dim(dist_mat_m)

close_thr <- 40000
far_thr   <- 193121

L_partial <- matrix(999, nrow(dist_mat_m), ncol(dist_mat_m))
L_partial[dist_mat_m > 0 & dist_mat_m < close_thr] <- -1
L_partial[dist_mat_m > far_thr] <- 0
diag(L_partial) <- 0

ref_edges <- data.frame()

for (i in 1:(nrow(L_partial) - 1)) {
  for (j in (i + 1):nrow(L_partial)) {
    if (L_partial[i, j] == -1) {
      ref_edges <- rbind(
        ref_edges,
        data.frame(
          from = i,
          to = j,
          from_station = station_info$STATION[i],
          to_station = station_info$STATION[j],
          from_name = station_info$NAME[i],
          to_name = station_info$NAME[j],
          from_lat = station_info$LATITUDE[i],
          from_lon = station_info$LONGITUDE[i],
          to_lat = station_info$LATITUDE[j],
          to_lon = station_info$LONGITUDE[j],
          dist_km = round(dist_mat_m[i, j] / 1000, 2)
        )
      )
    }
  }
}

nrow(ref_edges)
head(ref_edges)

library(htmlwidgets)
library(geosphere)
library(dplyr)
library(leaflet)
m_ref <- leaflet(station_info) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(
    lng = ~LONGITUDE,
    lat = ~LATITUDE,
    radius = 6,
    color = "red",
    stroke = TRUE,
    fillOpacity = 0.9,
    popup = ~paste0(
      "<b>Node:</b> ", Node, "<br>",
      "<b>Station:</b> ", STATION, "<br>",
      "<b>Name:</b> ", NAME, "<br>",
      "<b>Elevation:</b> ", ELEVATION
    )
  )

if (nrow(ref_edges) > 0) {
  for (k in 1:nrow(ref_edges)) {
    m_ref <- m_ref %>%
      addPolylines(
        lng = c(ref_edges$from_lon[k], ref_edges$to_lon[k]),
        lat = c(ref_edges$from_lat[k], ref_edges$to_lat[k]),
        color = "blue",
        weight = 2,
        opacity = 0.7,
        popup = paste0(
          "<b>Reference edge:</b> ", ref_edges$from[k], " - ", ref_edges$to[k], "<br>",
          "<b>From:</b> ", ref_edges$from_name[k], "<br>",
          "<b>To:</b> ", ref_edges$to_name[k], "<br>",
          "<b>Distance (km):</b> ", ref_edges$dist_km[k]
        )
      )
  }
}

m_ref
saveWidget(m_ref, "NOAA_reference_graph_map.html", selfcontained = TRUE)

L_learned <- fit_final$L

learned_edges <- data.frame()

for (i in 1:(nrow(L_learned) - 1)) {
  for (j in (i + 1):nrow(L_learned)) {
    if (abs(L_learned[i, j]) > 0.01) {
      learned_edges <- rbind(
        learned_edges,
        data.frame(
          from = i,
          to = j,
          weight = abs(L_learned[i, j]),
          from_station = station_info$STATION[i],
          to_station = station_info$STATION[j],
          from_name = station_info$NAME[i],
          to_name = station_info$NAME[j],
          from_lat = station_info$LATITUDE[i],
          from_lon = station_info$LONGITUDE[i],
          to_lat = station_info$LATITUDE[j],
          to_lon = station_info$LONGITUDE[j]
        )
      )
    }
  }
}

nrow(learned_edges)
head(learned_edges)


m_learned <- leaflet(station_info) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(
    lng = ~LONGITUDE,
    lat = ~LATITUDE,
    radius = 6,
    color = "red",
    stroke = TRUE,
    fillOpacity = 0.9,
    popup = ~paste0(
      "<b>Node:</b> ", Node, "<br>",
      "<b>Station:</b> ", STATION, "<br>",
      "<b>Name:</b> ", NAME, "<br>",
      "<b>Elevation:</b> ", ELEVATION
    )
  )

if (nrow(learned_edges) > 0) {
  max_w <- max(learned_edges$weight)
  
  for (k in 1:nrow(learned_edges)) {
    m_learned <- m_learned %>%
      addPolylines(
        lng = c(learned_edges$from_lon[k], learned_edges$to_lon[k]),
        lat = c(learned_edges$from_lat[k], learned_edges$to_lat[k]),
        color = "darkgreen",
        weight = 2 + 8 * learned_edges$weight[k] / max_w,
        opacity = 0.75,
        popup = paste0(
          "<b>Learned edge:</b> ", learned_edges$from[k], " - ", learned_edges$to[k], "<br>",
          "<b>Weight:</b> ", round(learned_edges$weight[k], 4), "<br>",
          "<b>From:</b> ", learned_edges$from_name[k], "<br>",
          "<b>To:</b> ", learned_edges$to_name[k]
        )
      )
  }
}

m_learned

isolated_nodes <- connectivity_df %>%
  filter(Degree == 0)

isolated_nodes



library(leaflet)
library(dplyr)
library(htmlwidgets)
library(GraphLearnR)

station_info <- noaa90 %>%
  distinct(STATION, NAME, LATITUDE, LONGITUDE, ELEVATION) %>%
  arrange(STATION)

station_info$Node <- 1:nrow(station_info)

#degree / connectivity from learned graph
adj_binary <- (fit_final$L != 0) * 1
diag(adj_binary) <- 0
deg <- rowSums(adj_binary)

station_info$Degree <- deg
station_info$NodeColor <- ifelse(station_info$Degree == 0, "black", "red")
station_info$NodeRadius <- 5 + 2 * station_info$Degree

#learned edge list
L_learned <- fit_final$L
learned_edges <- data.frame()

for (i in 1:(nrow(L_learned) - 1)) {
  for (j in (i + 1):nrow(L_learned)) {
    if (abs(L_learned[i, j]) > 0.01) {
      learned_edges <- rbind(
        learned_edges,
        data.frame(
          from = i,
          to = j,
          weight = abs(L_learned[i, j]),
          from_station = station_info$STATION[i],
          to_station = station_info$STATION[j],
          from_name = station_info$NAME[i],
          to_name = station_info$NAME[j],
          from_lat = station_info$LATITUDE[i],
          from_lon = station_info$LONGITUDE[i],
          to_lat = station_info$LATITUDE[j],
          to_lon = station_info$LONGITUDE[j]
        )
      )
    }
  }
}
#satellite map
m_learned_sat <- leaflet(station_info) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addCircleMarkers(
    lng = ~LONGITUDE,
    lat = ~LATITUDE,
    radius = ~NodeRadius,
    color = ~NodeColor,
    stroke = TRUE,
    weight = 2,
    fillOpacity = 0.9,
    popup = ~paste0(
      "<b>Node:</b> ", Node, "<br>",
      "<b>Station:</b> ", STATION, "<br>",
      "<b>Name:</b> ", NAME, "<br>",
      "<b>Elevation:</b> ", ELEVATION, "<br>",
      "<b>Degree:</b> ", Degree
    )
  )

if (nrow(learned_edges) > 0) {
  max_w <- max(learned_edges$weight)
  
  for (k in 1:nrow(learned_edges)) {
    m_learned_sat <- m_learned_sat %>%
      addPolylines(
        lng = c(learned_edges$from_lon[k], learned_edges$to_lon[k]),
        lat = c(learned_edges$from_lat[k], learned_edges$to_lat[k]),
        color = "cyan",
        weight = 2 + 8 * learned_edges$weight[k] / max_w,
        opacity = 0.8,
        popup = paste0(
          "<b>Learned edge:</b> ", learned_edges$from[k], " - ", learned_edges$to[k], "<br>",
          "<b>Weight:</b> ", round(learned_edges$weight[k], 4), "<br>",
          "<b>From:</b> ", learned_edges$from_name[k], "<br>",
          "<b>To:</b> ", learned_edges$to_name[k]
        )
      )
  }
}

m_learned_sat
saveWidget(m_learned_sat, "NOAA_learned_graph_satellite.html", selfcontained = TRUE)
saveWidget(m_learned_sat, "NOAA_learned_graph_satellite.html", selfcontained = TRUE)

plot(X_noisy[15, ], type = "l", col = "blue", lwd = 2,
     main = "Node 15 vs Node 16",
     ylab = "TMAX")

lines(X_noisy[16, ], col = "red", lwd = 2)

legend("topright", legend = c("Node 15", "Node 16"),
       col = c("blue", "red"), lwd = 2)


residual <- fit_final$hatS

plot(residual[15, ], type = "l", col = "blue", lwd = 2,
     main = "Residual Comparison: Node 15 vs 16",
     ylab = "Residual")

lines(residual[16, ], col = "red", lwd = 2)

legend("topright", legend = c("Node 15", "Node 16"),
       col = c("blue", "red"), lwd = 2)


reg_effect <- P[[1]] * beta[1] +
  P[[2]] * beta[2] +
  P[[3]] * beta[3]

residual_true <- fit_final$hatS - reg_effect

plot(residual_true[15, ], type = "l", col = "blue", lwd = 2,
     main = "True Residual: Node 15 vs 16",
     ylab = "Residual")

lines(residual_true[16, ], col = "red", lwd = 2)

legend("topright", legend = c("Node 15", "Node 16"),
       col = c("blue", "red"), lwd = 2)


beta_hat <- fit_final$beta

reg_effect <- P[[1]] * beta_hat[1] +
  P[[2]] * beta_hat[2] +
  P[[3]] * beta_hat[3]

residual_true <- fit_final$hatS - reg_effect

plot(residual_true[15, ], type = "l", col = "blue", lwd = 2,
     main = "True Residual: Node 15 vs 16",
     ylab = "Residual")

lines(residual_true[16, ], col = "red", lwd = 2)

legend("topright",
       legend = c("Node 15", "Node 16"),
       col = c("blue", "red"),
       lwd = 2)


library(igraph)
library(graphics)
adj_mat <- abs(fit_final$L)
diag(adj_mat) <- 0
adj_mat[adj_mat < 1e-10] <- 0

g_noaa <- graph_from_adjacency_matrix(
  adj_mat,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

deg <- degree(g_noaa)
edge_w <- E(g_noaa)$weight
edge_w <- 1 + 6 * edge_w / max(edge_w)

png("NOAA_laplacian_network.png", width = 900, height = 800)

plot(
  g_noaa,
  layout = layout_in_circle(g_noaa),
  vertex.size = 8 + 4 * deg,
  vertex.label = 1:vcount(g_noaa),
  vertex.label.cex = 0.8,
  vertex.color = "lightblue",
  vertex.frame.color = "navy",
  edge.width = edge_w,
  edge.color = "gray30",
  main = "NOAA Learned Graph from Laplacian"
)

dev.off()

png("NOAA_laplacian_heatmap.png", width = 900, height = 800)

image(
  1:nrow(L_learned), 1:ncol(L_learned),
  t(L_learned[nrow(L_learned):1, ]),
  col = heat.colors(100),
  xlab = "Node",
  ylab = "Node",
  main = "Learned Laplacian Matrix"
)

axis(1, at = 1:nrow(L_learned), labels = 1:nrow(L_learned), las = 2)
axis(2, at = 1:ncol(L_learned), labels = rev(1:ncol(L_learned)), las = 2)

dev.off()

library(dplyr)
library(geosphere)

station_info <- noaa90 %>%
  distinct(STATION, NAME, LATITUDE, LONGITUDE, ELEVATION) %>%
  arrange(STATION)

coords_mat <- as.matrix(station_info[, c("LONGITUDE", "LATITUDE")])
dist_mat_m <- geosphere::distm(coords_mat, fun = geosphere::distGeo)

# thresholds
close_thr <- 40000
far_thr   <- 193121

L_partial <- matrix(999, nrow(dist_mat_m), ncol(dist_mat_m))
L_partial[dist_mat_m > 0 & dist_mat_m < close_thr] <- -1
L_partial[dist_mat_m > far_thr] <- 0
diag(L_partial) <- 0


table(L_partial)

library(leaflet)
library(htmlwidgets)

ref_edges <- data.frame()

for (i in 1:(nrow(L_partial) - 1)) {
  for (j in (i + 1):nrow(L_partial)) {
    if (L_partial[i, j] == -1) {
      ref_edges <- rbind(
        ref_edges,
        data.frame(
          from = i,
          to = j,
          from_name = station_info$NAME[i],
          to_name = station_info$NAME[j],
          from_lat = station_info$LATITUDE[i],
          from_lon = station_info$LONGITUDE[i],
          to_lat = station_info$LATITUDE[j],
          to_lon = station_info$LONGITUDE[j],
          dist_km = round(dist_mat_m[i, j] / 1000, 2)
        )
      )
    }
  }
}

m_ref <- leaflet(station_info) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addCircleMarkers(
    lng = ~LONGITUDE,
    lat = ~LATITUDE,
    radius = 6,
    color = "yellow",
    stroke = TRUE,
    fillOpacity = 0.9,
    popup = ~paste0(
      "<b>Station:</b> ", STATION, "<br>",
      "<b>Name:</b> ", NAME
    )
  )

if (nrow(ref_edges) > 0) {
  for (k in 1:nrow(ref_edges)) {
    m_ref <- m_ref %>%
      addPolylines(
        lng = c(ref_edges$from_lon[k], ref_edges$to_lon[k]),
        lat = c(ref_edges$from_lat[k], ref_edges$to_lat[k]),
        color = "yellow",
        weight = 2,
        opacity = 0.8,
        popup = paste0(
          "<b>Reference edge</b><br>",
          ref_edges$from_name[k], "<br>",
          ref_edges$to_name[k], "<br>",
          "<b>Distance (km):</b> ", ref_edges$dist_km[k]
        )
      )
  }
}

m_ref
saveWidget(m_ref, "NOAA_partial_reference_map.html", selfcontained = TRUE)

metricsprf_partial <- function(L_partial, L_learned) {
  true_vals <- L_partial[lower.tri(L_partial)]
  learned_vals <- L_learned[lower.tri(L_learned)]
  
  keep <- true_vals != 999
  
  true_bin <- ifelse(true_vals[keep] == -1, 1, 0)
  learned_bin <- ifelse(abs(learned_vals[keep]) > 0, 1, 0)
  
  tp <- sum(true_bin == 1 & learned_bin == 1)
  fp <- sum(true_bin == 0 & learned_bin == 1)
  fn <- sum(true_bin == 1 & learned_bin == 0)
  tn <- sum(true_bin == 0 & learned_bin == 0)
  
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  
  data.frame(
    Precision = precision,
    Recall = recall,
    F1 = f1,
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn,
    Known_Entries = sum(keep)
  )
}

partial_metrics <- metricsprf_partial(L_partial, fit_final$L)
print(partial_metrics)
