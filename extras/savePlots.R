source("src/init.R")
iwalk(d, function(ldata, lname){
  iwalk(unname(ldata$medoidCoords), function(mdata, mnum){
    fname <- str_glue("C:/Users/u0118974/Pictures/Clouds/{lname}.{mnum}.png")
    mplot <- plotBasic(mdata) +
      theme(
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA)
      )
    ggsave(fname, plot = mplot, dpi = 600, bg = "transparent")
  })
})
