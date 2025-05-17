library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

# Rds format
path_data <- "../data"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
df_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)
path_images <-"../images"

# colnames(df)[colnames(df) == 'oldName'] <- 'newName'
colnames(df_core)[colnames(df_core) == 'entity_name'] <- 'Gene'
df_core$id <- as.numeric(df_core$id)


# Plots ----
# # Count gene per confidence level ----
# p1 <- df_core %>%
#   select(id, Gene, confidence_level) %>%
#   distinct() %>%
#   count(confidence_level) %>%
#   ggplot(aes(x = confidence_level, y = n) ) +
#   geom_bar(stat = "identity", fill = "orange", color = "black") + 
#   labs(y = "Number of\ngenes per\nconfidence level")
# p1
# 
# # No. panels where the same gene is included ----
# p2 <- df_core %>%
#   select(id, Gene) %>%
#   distinct() %>%
#   count(Gene) %>%
#   arrange(n) %>%
#   mutate(gene_index = row_number()) %>%
#   ggplot(aes(x = gene_index, y = n)) +
#   geom_bar(stat = "identity", fill = "blue") +
#   labs(title = "No. panels where the same gene is included",
#        x = "Gene index number",
#        y = "Number of\npanels with gene")
# p2
# 
# # No. genes per panel index ----
# p3 <- df_core %>%
#   select(id, Gene) %>%
#   distinct() %>%
#   select(id) %>%
#   count(id) %>%
#   arrange(n) %>%
#   mutate(id_index = row_number()) %>%
#   ggplot(aes(x = id_index, y = n)) +
#   geom_bar(stat = "identity", fill = "purple") +
#   labs(title = "Number of genes per panel ID", 
#        x = "Panel index number",
#        y = "Number of\npanels with gene")
# 
# library(patchwork)
# patch <- p1 / p2 / p3
# patch 

# ggsave(patch, file = paste0(path_images, "/plot_patch1.pdf") )

# Count gene per confidence level ----
p1 <- df_core %>%
  select(id, Gene, confidence_level) %>%
  distinct() %>%
  count(confidence_level) %>%
  ggplot(aes(x = confidence_level, y = n, fill = confidence_level) ) +
  guides(fill = "none") +
  geom_bar(stat = "identity", 
           # fill = "orange", 
           color = "black") +
  scale_fill_brewer(palette = "YlGnBu") +
  labs(y = "Number of\ngenes per\nconfidence level",
       x = "Confidence level")
p1

p2_data <- df_core %>%
  select(id, Gene) %>%
  distinct() %>%
  count(Gene) %>%
  arrange(n) %>%
  mutate(gene_index = row_number()) 


p2b <- p2_data %>%
  ggplot(aes(x = gene_index, y = n, fill =n)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low="orange1", high="red3") +
  # scale_fill_gradient2(low = "white", high = "red3", mid = "blue", midpoint = 20) +
  guides(fill = "none") +
  labs(
    # title = "No. panels where the same gene is included",
       x = "Gene index number",
       y = "Number of panels\nwith a given gene") + 
ggrepel::geom_label_repel(
  data = filter(p2_data, Gene == "RAG1"),  # Filter df_core to include only ID 467
  aes(label = stringr::str_wrap(Gene, width = 30)), 
  fill = "white",
  box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50',
  size = 3, color = "black", nudge_x = 0, nudge_y = 20
)

p2b

# example named
p3_data <- df_core %>%
  select(id, Gene, name) %>%
  distinct() %>%
  select(id, name) %>%
  count(id, name) %>%
  arrange(n) %>%
  mutate(id_index = row_number())

p3b <- p3_data %>%
  ggplot(aes(x = id_index, y = n, fill = n)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low="red3", high="purple4") +
  labs(#title = "Number of genes per panel ID", 
       x = "Panel index number",
       y = "Number of genes\nper panel", 
       ) +
  guides(fill = "none") + 
  ggrepel::geom_label_repel(
    data = filter(p3_data, id == 467 | id == 398),  # Filter df_core to include only ID 467
    aes(label = stringr::str_wrap(name, width = 30)), 
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50',
    fill = "white",
    size = 3, color = "black", nudge_x = -200, nudge_y = 2000
  )

p3b 

patch2 <- p1 / p2b / p3b + plot_layout(guides = 'collect', axis = "collect")  + plot_annotation(tag_levels = 'A')
patch2

ggsave(patch2, file = paste0(path_images, "/plot_patch2_annotated_example.pdf"), width = 5, height = 6)
ggsave(patch2, file = paste0(path_images, "/plot_patch2_annotated_example.png"), width = 5, height = 6)
# example info ----
print("An example is panel 398 with 572 PID genes which are well established as consensus in the community.")
print("An example is panel 1220 with 1675 genes which are associated with unexplained death in infancy and sudden unexplained death in childhood.")

