library(tidyverse)
library(gggenes)
library(ggnewscale)
library(ggside)

db_amanda <- readxl::read_xlsx(path = "interacoes.xlsx",
                               col_names = T)

db_amanda %>% 
  ggplot(aes(xmin = start, xmax = end, y = lncrna)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_subgene_arrow(aes(xmin = start, xmax = end,
                         y = lncrna, fill = mirna,
                         xsubmin = int_start, xsubmax = int_end),
                     color="black", alpha=.7) +
  theme_genes() +
  labs(y = "")


# Create sub db's ---------------------------------------------------------



db_lncrna <-
  db_amanda %>% 
  select(lncrna:mirna) %>% 
  rename(molecule = lncrna)

(
db_mirna <-
    db_amanda %>% 
    mutate(mir_pos_start = int_start - int_mir_start + 1,
           mir_pos_end = mir_pos_start + mir_end) %>% 
    select(mirna:mir_pos_end) %>%
    mutate(mir_pos_int_start = mir_pos_start + int_mir_start -1,
           mir_pos_int_end = mir_pos_start + int_mir_end -1) %>% 
    rename(molecule = mirna,
           start = mir_pos_start,
           end = mir_pos_end,
           int_start = mir_pos_int_start,
           int_end = mir_pos_int_end,
           "Binding Energy" = int) %>% 
    select(molecule, `Binding Energy`:int_end) %>% 
    mutate(
      #mirna = molecule,
      #molecule = "miRNA",
      molecule = fct_reorder(molecule, -`Binding Energy`))
)



db_lncrna %>% 
  ggplot(aes(xmin = start, xmax = end, y = molecule)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_subgene_arrow(aes(xmin = start, xmax = end,
                         y = molecule,
                         xsubmin = int_start, xsubmax = int_end),
                     fill = "dodgerblue",
                     color="black", alpha=.4) +
  geom_xsidehistogram(data = coverage,
                      aes(x = pos),
                      binwidth = 5,
                      inherit.aes = F,
                      fill = "dodgerblue") +
  new_scale_fill() +
  geom_gene_arrow(data = db_mirna,
                  aes(xmin = start, xmax = end,
                      y = molecule),
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(.5, "mm")) +
  geom_subgene_arrow(data = db_mirna,
                     aes(xmin = start, xmax = end,
                         y = molecule,
                         fill = `Binding Energy`,
                         xsubmin = int_start, xsubmax = int_end),
                     arrowhead_height = unit(3, "mm"),
                     arrowhead_width = unit(.5, "mm")) + 
  #scale_y_discrete(expand = c(0,2)) +
  scale_fill_gradient(low = "red",high = "white",
                      limits = c(min(db_amanda$int),0)) +
  #theme_genes() +
  labs(y = "",x = "") +
  theme(legend.position = "bottom",
        line = element_blank(),
        panel.grid.major.y = element_line(linewidth = .05,colour = "grey50"),
        panel.background = element_rect(fill = "white",
                                        colour = "black"))
  
ggsave(filename = "coord_interactions.png",dpi = 600,bg = "white",
       width = 12,
       height = 10)

# db_lncrna %>% 
#   ggplot(aes(xmin = start, xmax = end, y = molecule)) +
#   geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
#   geom_subgene_arrow(aes(xmin = start, xmax = end,
#                          y = molecule, fill = mirna,
#                          xsubmin = int_start, xsubmax = int_end),
#                      color="black", alpha=.7) +
#   theme_genes() +
#   labs(y = "")

# db_mirna %>% 
#   ggplot(aes(xmin = start, xmax = end, y = mirna)) +
#   geom_gene_arrow(data = db_mirna,
#                   aes(xmin = start, xmax = end,
#                       y = molecule)) +
#   geom_subgene_arrow(data = db_mirna,
#                      aes(xmin = start, xmax = end,
#                          y = molecule,
#                          fill = `Binding Energy`,
#                          xsubmin = int_start, xsubmax = int_end)) + 
#   scale_y_discrete(expand = c(0,6)) +
#   theme_genes() +
#   labs(y = "")


limits <- select(db_amanda,int_start,int_end)

seq_pos <- list()
for(i in 1:nrow(limits)){
  
  seq_pos[[i]] <- seq(limits$int_start[i],limits$int_end[i],1)
}

coverage <-
  tibble(seq_pos) %>% 
  unnest() %>% 
  rename(pos = seq_pos)

# Only mirs ---------------------------------------------------------------


db_amanda %>% 
  ggplot(aes(xmin = mir_start, xmax = mir_end, y = mirna)) +
  geom_gene_arrow() +
  geom_subgene_arrow(aes(xmin = mir_start, xmax = mir_end,
                         y = mirna, fill = int,
                         xsubmin = int_mir_start, xsubmax = int_mir_end),
                     color="black", alpha=.7) +
  scale_fill_gradient(low = "red",high = "white",
                      limits = c(min(db_amanda$int),0)) +
  scale_y_discrete(expand = c(0,6)) +
  theme_genes() +
  labs(y = "")
