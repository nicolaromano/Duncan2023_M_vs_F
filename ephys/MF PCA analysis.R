library(ggplot2)
library(ggrepel)
library(car)
library(xtable)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggcorrplot)
library(gridExtra)

# Read the data in
mf <- read.csv("MF_ephys.csv")

mf %>%
  separate(Group, sep = "-", into = c("Sex", "Treatment")) %>%
  mutate(SexShort = substr(Sex, 1, 1), .after = Sex) %>%
  # Convert NAs to 0
  replace(is.na(.), 0) -> mf

mf %>%  
  # Remove unused variables
  select(-c(Cell, Treatment, Sex, SexShort, Treatment, Days.in.Culture)) -> mf.2

# z-score the values
mf.2.norm <- apply(mf.2, 2, function(x) {
  (x - mean(x)) / diff(range(x))
})

# Check correlation between variables
var.cor <- cor(mf.2)
colnames(var.cor) <- rownames(var.cor) <- sub("\\.", " ", colnames(var.cor))

ggcorrplot(var.cor,
  hc.order = TRUE, type = "lower",
  lab = TRUE, lab_size = 5, colors = c("#6D9EC1", "white", "#f1732f")
)

# Run PCA
mf.pca <- princomp(mf.2.norm)

# Remove variables with high correlation
# Burst.Factor, Active.Time
var.cor <- cor(mf.2 %>% select(-Burst.Factor, -Active.Time, -Event.Frequency, -Event.Duration))
colnames(var.cor) <- rownames(var.cor) <- sub("\\.", " ", colnames(var.cor))

ggcorrplot(var.cor,
  hc.order = TRUE, type = "lower",
  lab = TRUE, lab_size = 5, colors = c("#6D9EC1", "white", "#f1732f")
)

# Redo PCA
mf.3 <- mf.2 %>% select(-Burst.Factor, -Active.Time, -Event.Frequency, -Event.Duration)

mf.3.norm <- apply(mf.3, 2, function(x) {
  (x - mean(x)) / diff(range(x))
})

mf.pca.no_corr <- princomp(mf.3.norm, scale = TRUE)

# Screeplots
p1 <- mf.pca$sdev %>%
  as.data.frame() %>%
  rename_at(vars(1), ~"Variance explained") %>%
  mutate(PC = 1:nrow(.)) %>%
  ggplot(aes(x = PC, y = mf.pca$sdev)) +
  geom_col(fill = "#39aee0") +
  scale_x_continuous(breaks = seq(1, length(mf.pca$sdev), 1)) +
  geom_text(aes(label = round(mf.pca$sdev, 2)), vjust = -0.5) +
  xlab("Principal dimension") +
  ylab("% variance explained") +
  ggtitle("All variables") +
  theme_bw() +
  theme(
    title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

p2 <- mf.pca.no_corr$sdev %>%
  as.data.frame() %>%
  rename_at(vars(1), ~"Variance explained") %>%
  mutate(PC = 1:nrow(.)) %>%
  ggplot(aes(x = PC, y = mf.pca.no_corr$sdev)) +
  geom_col(fill = "#39aee0") +
  scale_x_continuous(breaks = seq(1, length(mf.pca.no_corr$sdev), 1)) +
  geom_text(aes(label = round(mf.pca.no_corr$sdev, 2)), vjust = -0.5) +
  xlab("Principal dimension") +
  ylab("% variance explained") +
  ggtitle("Removed highly correlated") +
  theme_bw() +
  theme(
    title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

grid.arrange(p1, p2, ncol = 2)

pretty_pca <- function(pca, plot_direction = FALSE) {
  #' Pretty PCA plot
  #' Params
  #' ------
  #' pca: PCA object
  #' plot_direction: logical, whether to plot the average Control->CRH direction

  my.pca <- data.frame(
    PC1 = pca$scores[, 1],
    PC2 = pca$scores[, 2],
    Cell = mf$Cell,
    Sex = mf$Sex,
    Treatment = mf$Treatment
  )

  my.pca$PC1_CRH <- NULL
  my.pca$PC1_CRH[my.pca$Treatment == "Base"] <- my.pca$PC1[my.pca$Treatment == "CRH"]
  my.pca$PC2_CRH <- NULL
  my.pca$PC2_CRH[my.pca$Treatment == "Base"] <- my.pca$PC2[my.pca$Treatment == "CRH"]

  cells <- data.frame(aggregate(row(mf)[, 1],
    by = list(
      Cell = mf$Cell,
      Sex = mf$Sex
    ),
    FUN = identity
  ))

  # Aggregate returns a funny data frame in this case,
  # with two columns appearing as nested into a single one...
  # Converting to matrix and back to dataframe seems
  # to solve the issue.
  cells <- as.matrix(cells)
  cells <- as.data.frame(cells, stringsAsFactors = F)
  colnames(cells)[3:4] <- c("Base", "CRH")
  cells$Base <- as.integer(cells$Base)
  cells$CRH <- as.integer(cells$CRH)

  cells$PC1.start <- mf.pca$scores[cells$Base, "Comp.1"]
  cells$PC1.end <- mf.pca$scores[cells$CRH, "Comp.1"]
  cells$PC2.start <- mf.pca$scores[cells$Base, "Comp.2"]
  cells$PC2.end <- mf.pca$scores[cells$CRH, "Comp.2"]

  centers <- aggregate(my.pca[, 1:2],
    by = list(Sex = my.pca$Sex, Treatment = my.pca$Treatment),
    FUN = mean
  )
  centers$PC1_CRH <- NA
  centers$PC2_CRH <- NA
  centers$PC1_CRH[centers$Treatment == "Base"] <- centers$PC1[centers$Treatment == "CRH"]
  centers$PC2_CRH[centers$Treatment == "Base"] <- centers$PC2[centers$Treatment == "CRH"]

  # PCA plot
  g <- ggplot(my.pca, aes(PC1, PC2)) +
    geom_point(aes(pch = Treatment), size = 5, col = "#4D4D4D") +
    scale_shape_manual(values = 21:20) +
    geom_segment(aes(x = PC1, xend = PC1_CRH, y = PC2, yend = PC2_CRH),
      col = "lightgray", alpha = 0.5
    ) +
    xlab(expression(PC[1])) +
    ylab(expression(PC[2]))

  if (plot_direction) {
    g <- g +
      geom_segment(
        data = centers,
        aes(
          x = PC1, y = PC2, xend = PC1_CRH, yend = PC2_CRH,
          col = Sex
        ),
        arrow = arrow(type = "closed", length = unit(0.1, "inches"))
      ) +
      scale_color_discrete(guide = "none")
  }

  g <- g +
    facet_wrap(~Sex) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(size = .5),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      strip.text = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    )

  print(g)

  return(g)
}

pretty_loadings_plot <- function(pca) {
  data.frame(pca$loadings[, 1:2]) %>%
    rename(PC1 = Comp.1, PC2 = Comp.2) %>%
    rownames_to_column("Property") %>%
    mutate(Property = str_replace(Property, "\\.", " ")) %>%
    mutate(Magnitude = sqrt(PC1^2 + PC2^2)) -> loadings

  ggplot(loadings) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2, col = Magnitude),
      arrow = arrow(
        angle = 30,
        type = "closed",
        length = unit(0.1, "inches")
      )
    ) +
    scale_color_gradientn(colors = c("#3093e4", "#f1be66", "#ff4800")) +
    geom_text_repel(aes(x = PC1, y = PC2, label = Property),
      size = 6,
      min.segment.length = 2,
      hjust = 0,
      nudge_y = sign(loadings$PC2) * 0.02
    ) +
    xlab(expression(PC[1])) +
    ylab(expression(PC[2])) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    )
}

p1 <- pretty_pca(mf.pca, plot_direction = TRUE)
p2 <- pretty_pca(mf.pca.no_corr, plot_direction = TRUE)

pdf("pca_allvars.pdf", width = 18, height = 10)
print(p1)
dev.off()

pdf("pca_no_high_corr.pdf", width = 18, height = 12)
print(p2)
dev.off()

p1 <- pretty_loadings_plot(mf.pca)
p2 <- pretty_loadings_plot(mf.pca.no_corr)

pdf("loadings_allvars.pdf", width = 12, height = 10)
print(p1)
dev.off()

pdf("loadings_no_high_corr.pdf", width = 12, height = 10)
print(p2)
dev.off()