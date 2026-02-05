library(RColorBrewer)
library(colorspace)
library(ggplot2)

show_col(viridis(6))
show_col(inferno(6))
show_col(turbo(6))

# Extract the default ggplot2 color scheme
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# Save the default color scheme for 12 color palette
palette <- ggplotColours(n=12)
palette <- RColorBrewer::brewer.pal(n=8, name="Dark2")
palette <- RColorBrewer::brewer.pal(n=12, name="Paired")
palette <- RColorBrewer::brewer.pal(n=11, name="RdYlBu")
palette <- RColorBrewer::brewer.pal(11, "RdYlBu")[c(1,11)]

n <- 25
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

# Visualize the colors
pie(rep(1,length(palette)), col=palette)

# Print the color codes
my_palette <- c("#FF1F5B", "#F28522", "#009ADE", "#AF58BA", "#3C005A",
                "#00B000", "#FFC61E", "#808080", "#A6761D", "#F6D2E0",
                "#E75480", "#C8E7F5", "#2E5984", "#FFFF99", "#B15928",
                "#A8A9AD", "#000000", "#BF812D", "#35978F", "#C51B7D",
                "#7FBC41", "#762A83", "#D6604D", "#4393C3", "#FFFFBF",
                "#9E0142", "#E41A1C", "#4DAF4A", "#FF7F00", "#FFFF33",
                "#A65628", "#F781BF", "#66C2A5", "#FC8D62", "#1F77B4",
                "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
                "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8",
                "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94",
                "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5", "#333333",
                "#B8E80C", "#DE77AE", "#F1B6DA", "#FDE0EF", "#F7F7F7",
                "#E6F5D0", "#B8E186", "#FFC3FF", "#FFFF72", "#E08214",
                "#542788", "#878787", "#1A1A1A", "#377EB8", "#984EA3",
                "#999999")

scanpy_default_102 <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
  "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
  "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
  "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#e6550d", "#31a354",
  "#756bb1", "#636363", "#9e9ac8", "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4",
  "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec", "#f2f2f2", "#a6cee3", "#1f78b4",
  "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
  "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072",
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5",
  "#ffed6f", "#a65628", "#ffff99", "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  "#66a61e", "#e6ab02", "#a6761d", "#666666", "#e41a1c", "#377eb8", "#4daf4a",
  "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5",
  "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d",
  "#666666")

glasbey_100 <- c(
  "#0000ff", "#ff0000", "#00ff00", "#000033", "#ff00b6", "#005300", "#ffd300", "#b6ff00", "#4f00ff", "#ff6eff",
  "#6200ff", "#00ff90", "#ffb6c1", "#0077ff", "#ff4f00", "#aaff00", "#00ffd0", "#ff0090", "#004dff", "#ffa600",
  "#b600ff", "#00b6ff", "#ffd0b6", "#8cff00", "#ff0033", "#00ff4f", "#0033ff", "#d000ff", "#b6b6ff", "#a6ff00",
  "#00aaff", "#ff00d0", "#007777", "#ffb66e", "#3300ff", "#ff0077", "#4fff00", "#00d6ff", "#ffb600", "#7200ff",
  "#00ffb6", "#ff0066", "#003366", "#ffb600", "#00ff33", "#d6ff00", "#ff0033", "#00b66e", "#ff6e00", "#3300cc",
  "#00ffd6", "#ff3399", "#6600ff", "#00ff66", "#ff00ff", "#ffcc00", "#6699ff", "#cc0066", "#00cc00", "#cc6600",
  "#9966ff", "#66cc00", "#0066cc", "#cc0099", "#009933", "#996600", "#660066", "#0099cc", "#666600", "#cc00cc",
  "#00cc99", "#ff6600", "#00cccc", "#ff9999", "#ccff00", "#006666", "#9900cc", "#cc3333", "#999900", "#660000",
  "#00ffcc", "#6600cc", "#ffcccc", "#9999ff", "#cc0000", "#339900", "#ffcc99", "#006600", "#cc99cc", "#336600",
  "#993333", "#00cc66", "#ff99cc", "#333366", "#ccffcc", "#999966", "#66cccc", "#ff0066", "#660066", "#9999cc")


iwanthue_100 <- c(
  "#b12a90", "#3acb30", "#3f51cb", "#d03239", "#27c8a1", "#ad45d6", "#cb9a36", "#8a4fa3", "#658d33", "#ca4aa2",
  "#4369da", "#78c13b", "#2d3c92", "#cf3b7f", "#38c76a", "#bc3669", "#50c9c2", "#5c3fa3", "#9b9735", "#5272b4",
  "#ce6854", "#3aaf4a", "#3e50a1", "#cb517e", "#9dc55b", "#5b60b3", "#3f947a", "#6a46c3", "#c85c36", "#349da5",
  "#6a9634", "#865ac2", "#bf3850", "#4a8a2e", "#c94f96", "#3685cb", "#d0472e", "#51b06f", "#8938a4", "#7cbc42",
  "#b84f35", "#2c86c4", "#9e4f93", "#3a9e56", "#a25acd", "#7c9e3c", "#c6365e", "#49c192", "#ba3b97", "#349353",
  "#6751c6", "#c13b3f", "#2fa77b", "#b045c6", "#43a03e", "#bb497c", "#3e8fcc", "#d36031", "#48c07b", "#8d41cb",
  "#779f30", "#c73543", "#50a697", "#b944ab", "#3cb753", "#8b4fb6", "#cf3e64", "#52b16d", "#4f46a4", "#b96833",
  "#46c3b3", "#6f3d98", "#92ac3b", "#a348ce", "#53a15e", "#c13a9d", "#458a3e", "#d54b4d", "#3f88b7", "#bd3680",
  "#4cb17c", "#8455c9", "#a34b39", "#55be5b", "#6247c2", "#b14957", "#2da66e", "#c13688", "#3ea747", "#b63f95",
  "#55a06f", "#7c52c5", "#d6462f", "#3db383", "#7a47ab", "#b43d35", "#34a4b7", "#9761c6", "#61a633", "#cb3653")



my_palette <- c(my_palette, 
                 colorspace::adjust_transparency(col = my_palette, alpha = 0.6), 
                 colorspace::adjust_transparency(col = my_palette, alpha = 0.3))
my_palette <-c("grey", viridis(n = 10, option = "C", direction = -1))
hcl_palettes(plot = TRUE)


# Visualize the colors
pie(rep(1,length(scanpy_default_102)), col=scanpy_default_102)
