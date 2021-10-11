# Prep for NephoMonograph
source(here::here("src", "init.R"))

# drviz
hach_mds <- read_tsv(here("assets", "hachelijk.mds.tsv"), col_types = cols())
hach_umap <- read_tsv(here("assets", "hachelijk.umap.tsv"), col_types = cols())                                                                                    
top_grid <- plot_grid(                                                                                                                                             
     plotBasic(hach_mds),
     plotBasic(d$hachelijk$medoidCoords[[1]]),
     labels = c("NMDS", "t-SNE"), label_fontfamily = "Bookman Old Style")
plot_grid(
     top_grid, plotBasic(hach_umap),
     labels = c("", "UMAP"),
     label_fontfamily = "Bookman Old Style", ncol = 1, hjust = -4)
ggsave("../NephoMonograph/data/img/drviz.png")


# Popular plots
popular_plots <- map(
   setNames(to_show, to_show),
   ~plotWithCws(.x, str_glue("{.x}.{popular_medoid}"), colorguide = "none"))
plotPopular <- function(...) plot_grid(plotlist = popular_plots[c(...)], labels = c(...),
                                       label_fontfamily = "Bookman Old Style")
plotPopular("heet", "stof")
ggsave("../NephoMonograph/data/img/popular_heet_stof.png")

plotPopular("dof", "huldigen")
ggsave("../NephoMonograph/data/img/popular_dof_huldigen.png")

plotPopular("haten", "hoop")
ggsave("../NephoMonograph/data/img/popular_haten_hoop.png")

# Semantic interpretation
saveSemInt <- function(l, n){
  plotWithCws(l, n)
  ggsave(str_glue("../NephoMonograph/data/img/{l}-{n}.png"))
}
saveSemInt("heilzaam", 1)
saveSemInt("schaal", 1)
saveSemInt("heffen", 1)
saveSemInt("hachelijk", 1)
saveSemInt("stof", 3)
saveSemInt("herstructureren", 7)
saveSemInt("herhalen", 7)
saveSemInt("herinneren", 3)
saveSemInt("horde", 4)
saveSemInt("grijs", 4)
saveSemInt("herroepen", 6)
saveSemInt("haken", 1)
plotWithCws('heet', 3) +
  scale_shape_manual(values = c(16, 17, 18, 0, 2, 4, 8), guide = "none")
ggsave("../NephoMonograph/data/img/heet-3.png")
saveSemInt("geldig", 1)
saveSemInt("blik", 3)
saveSemInt("huldigen", 6)
network_cws(d$huldigen$medoidCoords[[6]]$coords %>% filter(cluster == 1))
ggsave("../NephoMonograph/data/img/huldigen-net.png")
