### Preprocessing and Cleaning ADNI data

# Load library and functions 

setwd("Desktop/Neurotract")
options(digits = 3)
library(knitr)
library(ADNIMERGE)
library(ggplot2)
library(Hmisc)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
library(kableExtra)
library(consort)
library(gtsummary)

source("https://adni.bitbucket.io/myfunctions.R")
theme_set(theme_bw())
dxpal <- c( "#9CCDBA", "#F1A84A", "#825C9F")
scale_colour_discrete <- function(...) scale_colour_manual(..., values = dxpal)
scale_fill_discrete <- function(...) scale_fill_manual(..., values = dxpal)

# Define the datasets

dfall <- adnimerge %>%
  select(RID, ORIGPROT, COLPROT, VISCODE, DX, DX.bl, AGE, Years.bl, PTGENDER, 
         PTEDUCAT, APOE4, ADAS11, ADAS11.bl, ADAS13, ADAS13.bl, ADASQ4, 
         ADASQ4.bl, CDRSB, CDRSB.bl, MOCA, MOCA.bl, MMSE, MMSE.bl, ABETA, 
         ABETA.bl, TAU, TAU.bl, PTAU, PTAU.bl)

phs <- desikanlab
gap43 <- blennow_lab_csf_gap_43
nfl.bl <- blennowplasmanfl 
nfl.lng <- blennowplasmanfllong

# Create a unique Sample ID (SID) for every sample of each individual

dfall$SID <- paste0(dfall$RID, "_", dfall$VISCODE)
gap43$SID <- paste0(gap43$RID, "_", gap43$VISCODE)
nfl.bl$SID <- paste0(nfl.bl$RID, "_", nfl.bl$VISCODE)
nfl.lng$SID <- paste0(nfl.lng$RID, "_", nfl.lng$VISCODE)


# Merge all data in the single data frame  

dfall$RID <- as.numeric(dfall$RID)
phs <- phs[, -c(1,3)]
gap43 <- gap43[, c(7,8)]
nfl.bl <- nfl.bl[, c(7,9)]
nfl.lng <- nfl.lng[, c(9,11)]
nfl <- rbind(nfl.bl, nfl.lng)

mdf <- dfall %>%
  full_join(phs, by = "RID") %>%
  full_join(gap43, by = "SID") %>%
  full_join(nfl, by = "SID")

## Assessment eligibility of all RIDs for the downstream analysis

consort.id.raw <- sort(unique(mdf$RID))
n_subjects <- length(consort.id.raw)

# Filter subjects with only one visit

ids_one_visit <- mdf %>%
  group_by(RID) %>%
  summarise(n_visit = n()) %>%
  filter(n_visit == 1) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects with single visit only =", length(ids_one_visit))

is_single_visit <- consort.id.raw %in% ids_one_visit
consort.id.singlevisit <- rep(NA, n_subjects) 
consort.id.singlevisit[is_single_visit] <- "One visit only" 
consort.id.interim1 <- consort.id.raw 
consort.id.interim1[is_single_visit] <- NA

paste("[Count] Subjects removed =", sum(!is.na(consort.id.singlevisit)))

# Filter subjects without any diagnosis at baseline

ids_DXbl_NA <- mdf %>%
  filter(VISCODE == "bl") %>%
  filter(is.na(DX)) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects with DX==NA at baseline =", length(ids_DXbl_NA))

is_DXbl_NA <- consort.id.interim1 %in% ids_DXbl_NA
consort.id.DXbl_NA <- rep(NA, n_subjects) 
consort.id.DXbl_NA[is_DXbl_NA] <- "No baseline diagnosis"
consort.id.interim2 <- consort.id.interim1
consort.id.interim2[is_DXbl_NA] <- NA

# Filter subjects without any diagnosis after baseline

ids_DX_allNA <- mdf %>%
  filter(VISCODE != "bl") %>%
  group_by(RID) %>%
  summarise(
    all_DX_NA = all(is.na(DX))
  ) %>% 
  filter(all_DX_NA == TRUE) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects without any DX after baseline", length(ids_DX_allNA))

is_DX_allNA <- consort.id.interim2 %in% ids_DX_allNA
consort.id.DX_allNA <- rep(NA, n_subjects) 
consort.id.DX_allNA[is_DX_allNA] <- "No diagnosis after baseline" 
consort.id.interim3 <- consort.id.interim2
consort.id.interim3[is_DX_allNA] <- NA

# Filter subjects with NA or unknown value in baseline covariates

vars_baseline <- c("AGE", "PTGENDER", "PTEDUCAT", "APOE4")

ids_base_NA <- lapply(vars_baseline, function(x) {
  ids_NA <- mdf %>%
    filter(VISCODE == "bl") %>%
    filter(is.na(get(x)) | get(x) == "Unknown") %>%
    select(RID) %>%
    unlist(use.names = FALSE)
  print(paste("[Count] Subjects with NA or unknown in baseline covariate", 
              x, " =", length(ids_NA)))
  return (ids_NA)
})

is_no_basecov <- consort.id.interim3 %in% Reduce(union, ids_base_NA)
consort.id.no_basecov <- rep(NA, n_subjects) 
consort.id.no_basecov[is_no_basecov] <- "NA in baseline covariates"
consort.id.interim4 <- consort.id.interim3
consort.id.interim4[is_no_basecov] <- NA


# Remove missing data

ids_exclude <- Reduce(union, 
                      c(ids_base_NA, list(ids_one_visit, ids_DXbl_NA, ids_DX_allNA)))

paste("[Count] Subjects to be excluded due to single visit or without diagnosis at baseline or no DX after baseline =", 
      length(ids_exclude))

mdf.clean <- mdf %>%
  filter(!(RID %in% ids_exclude))

mdf.clean[, "PTGENDER"] <- as.factor(mdf.clean[, "PTGENDER"])
mdf.clean[, "APOE4"] <- as.factor(mdf.clean[, "APOE4"])

# Consort diagram to describe exclusion of subjects                                  

df.consort <- data.frame(
  consort.id.raw,
  consort.id.singlevisit, 
  consort.id.interim1,
  consort.id.DXbl_NA,
  consort.id.interim2,
  consort.id.DX_allNA,
  consort.id.interim3,
  consort.id.no_basecov,
  consort.id.interim4)

df.consort <- df.consort %>% 
  unite(
    "consort.id.singlevisit", 
    "consort.id.DXbl_NA", 
    "consort.id.DX_allNA",
    "consort.id.no_basecov", 
    sep = "", remove = TRUE, 
    na.rm = TRUE) %>%
  rename(consort.id.excl = "consort.id.singlevisit")

df.consort$consort.id.excl[df.consort$consort.id.excl == ""] <- NA

g.consort <- consort_plot(data = df.consort,
                          orders = c(consort.id.raw = "ADNI participants",
                                     consort.id.excl = "Excluded",
                                     consort.id.interim4 = "Subjects for analysis"),
                          side_box = c("consort.id.excl"),
                          cex = 0.7)

ggsave("consort_diagram.pdf", plot = g.consort, device = "pdf",
       path = "Results/",
       width = 6, height = 2)

save(mdf.clean, 
     file = "Data/adni_cleaned.RData")

# Filter subjects without any PHS

ids_PHS_allNA <- mdf.clean %>%
  group_by(RID) %>%
  summarise(
    all_PHS_NA = all(is.na(PHS))
  ) %>% 
  filter(all_PHS_NA == TRUE) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects without any PHS", length(ids_PHS_allNA))

is_PHS_allNA <- consort.id.interim4 %in% ids_PHS_allNA
consort.id.PHS_allNA <- rep(NA, n_subjects) 
consort.id.PHS_allNA[is_PHS_allNA] <- "NA in PHS" 
consort.id.interim5 <- consort.id.interim4
consort.id.interim5[is_PHS_allNA] <- NA

# Filter subjects without any GAP.43

ids_GAP.43_allNA <- mdf.clean %>%
  group_by(RID) %>%
  summarise(
    all_GAP.43_NA = all(is.na(GAP.43))
  ) %>% 
  filter(all_GAP.43_NA == TRUE) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects without any GAP.43", length(ids_GAP.43_allNA))

is_GAP.43_allNA <- consort.id.interim5 %in% ids_GAP.43_allNA
consort.id.GAP.43_allNA <- rep(NA, n_subjects) 
consort.id.GAP.43_allNA[is_GAP.43_allNA] <- "NA in GAP-43" 
consort.id.interim6 <- consort.id.interim5
consort.id.interim6[is_GAP.43_allNA] <- NA

# Filter subjects without any PLASMA_NFL

ids_PLASMA_NFL_allNA <- mdf.clean %>%
  group_by(RID) %>%
  summarise(
    all_PLASMA_NFL_NA = all(is.na(PLASMA_NFL))
  ) %>% 
  filter(all_PLASMA_NFL_NA == TRUE) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects without any PLASMA_NFL", length(ids_PLASMA_NFL_allNA))

is_PLASMA_NFL_allNA <- consort.id.interim6 %in% ids_PLASMA_NFL_allNA
consort.id.PLASMA_NFL_allNA <- rep(NA, n_subjects) 
consort.id.PLASMA_NFL_allNA[is_PLASMA_NFL_allNA] <- "NA in PLASMA_NFL" 
consort.id.interim7 <- consort.id.interim6
consort.id.interim7[is_PLASMA_NFL_allNA] <- NA

# Remove missing data

ids_exclude <- Reduce(union, 
                      c(list(ids_PHS_allNA, ids_GAP.43_allNA, ids_PLASMA_NFL_allNA)))

paste("[Count] Subjects to be excluded due to without PHS, GAP-43 and NFL =", 
      length(ids_exclude))

df <- mdf.clean %>%
  filter(!(RID %in% ids_exclude))

df$ABETA[df$ABETA == ">1700"] <- "1700"
df$ABETA <- as.numeric(df$ABETA)
df$TAU <- as.numeric(df$TAU)
df$PTAU <- as.numeric(df$PTAU)
df$ABETA.bl[df$ABETA.bl == ">1700"] <- "1700"
df$ABETA.bl <- as.numeric(df$ABETA.bl)
df$TAU.bl <- as.numeric(df$TAU.bl)
df$PTAU.bl <- as.numeric(df$PTAU.bl)
df$DX <- gsub("Dementia", "AD", df$DX)

df <- merge(df, df.bl.cp[, c("RID", "GAP.43", "PLASMA_NFL")], by = "RID", all.x = TRUE)
colnames(df)[colnames(df) == "GAP.43.x"] <- "GAP.43"
colnames(df)[colnames(df) == "PLASMA_NFL.x"] <- "PLASMA_NFL"
colnames(df)[colnames(df) == "GAP.43.y"] <- "GAP.43.bl"
colnames(df)[colnames(df) == "PLASMA_NFL.y"] <- "PLASMA_NFL.bl"

# Remove technical replicate problems of GAP-43 and NFL data
rownames(df)[duplicated(df$SID)]
df <- df[!(rownames(df) %in% 
             rownames(df)[duplicated(df$SID)]), ]


df.bl <- df[df$VISCODE == "bl",]
df.bl <- na.omit(df.bl[, 33])
df.bl <- df.bl[complete.cases(df.bl$GAP.43), , drop = FALSE]
df.bl <- df.bl[complete.cases(df.bl$PLASMA_NFL), , drop = FALSE]
colSums(is.na(df.bl))
dff <- df[df$RID %in% df.bl$RID, ]

# Chcek distribution normality

columns <- c("GAP.43", "PLASMA_NFL", "MMSE.bl", "CDRSB.bl", "ADAS11.bl", 
             "ADAS13.bl", "ADASQ4.bl", "MOCA.bl", "PHS", "ABETA.bl", 
             "TAU.bl", "PTAU.bl", "PTEDUCAT", "AGE")
##Kolmogorov-Smirnov normality test
##ks_results <- list()
##for (col in columns) {
  test_result <- ks.test(ddf[[col]], "pnorm", mean = mean(ddf[[col]], na.rm = TRUE), sd = sd(ddf[[col]], na.rm = TRUE))
  ks_results[[col]] <- test_result
}
##ks_results

# Shapiro-Wilk normality test

shapiro_results <- list()
for (col in columns) {
  test_result <- shapiro.test(ddf[[col]])
  shapiro_results[[col]] <- test_result
}
shapiro_results

p_values <- data.frame(Column = columns, 
                       p_value = sapply(shapiro_results, 
                                        function(result) result$p.value))

plots <- list()
for (col in columns) {
  median_val <- median(ddf[[col]], na.rm = TRUE)
  mean_val <- mean(ddf[[col]], na.rm = TRUE)
  
  p <- ggplot(ddf, aes_string(x = col)) +
    geom_density(color = "black", fill = "grey") +
    geom_vline(xintercept = median_val, color = "red", linetype = "dashed", size = 0.75, show.legend = TRUE) +
    geom_vline(xintercept = mean_val, color = "blue", linetype = "dashed", size = 0.75, show.legend = TRUE) +
    labs(title = paste("Shapiro-Wilk p-value = ", p_values[col,2]),
         x = col,
         y = "Density") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(size = 10)
    )
  
  plots[[col]] <- p
}

do.call("grid.arrange", c(plots, ncol = 2))

mdf.clean.bl <- mdf.clean[mdf.clean$VISCODE == 0,]
columns <- c( "MMSE.bl", "CDRSB.bl", "ADAS11.bl", 
             "ADAS13.bl", "ADASQ4.bl", "MOCA.bl", "PHS", 
             "PTEDUCAT", "AGE")

shapiro_results <- list()
for (col in columns) {
  test_result <- shapiro.test( mdf.clean.bl[[col]])
  shapiro_results[[col]] <- test_result
}
shapiro_results

p_values <- data.frame(Column = columns, 
                       p_value = sapply(shapiro_results, 
                                        function(result) result$p.value))

plots <- list()
for (col in columns) {
  median_val <- median( mdf.clean.bl[[col]], na.rm = TRUE)
  mean_val <- mean( mdf.clean.bl[[col]], na.rm = TRUE)
  
  p <- ggplot( mdf.clean.bl, aes_string(x = col)) +
    geom_density(color = "black", fill = "grey") +
    geom_vline(xintercept = median_val, color = "red", linetype = "dashed", size = 0.75, show.legend = TRUE) +
    geom_vline(xintercept = mean_val, color = "blue", linetype = "dashed", size = 0.75, show.legend = TRUE) +
    labs(title = paste("Shapiro-Wilk p-value = ", p_values[col,2]),
         x = col,
         y = "Density") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(size = 10)
    )
  
  plots[[col]] <- p
}

do.call("grid.arrange", c(plots, ncol = 2))

# Basic demography

ddf %>% 
  select(c(AGE,PTGENDER,PTEDUCAT,APOE4,PHS,ABETA.bl,TAU.bl,PTAU.bl,
           GAP.43,PLASMA_NFL,MMSE.bl,CDRSB.bl,ADAS11.bl,ADAS13.bl,
           ADASQ4.bl,MOCA.bl,DX)) %>% 
  tbl_summary(by = DX, statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n} / {N} ({p}%)"
  )) %>% 
  add_p()

mdf.clean %>% 
  select(c(AGE,PTGENDER,PTEDUCAT,APOE4,PHS,
           MMSE.bl,CDRSB.bl,ADAS11.bl,ADAS13.bl,
           ADASQ4.bl,MOCA.bl,DX)) %>% 
  tbl_summary(by = DX, statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n} / {N} ({p}%)"
  )) %>% 
  add_p()

# Saved cleaned data
write.table(ddf, 
            file = "~/drive/MyDrive/df_bl.csv",
            quote = F, sep = ",", row.names = F, col.names = T)

write.table(mdf.clean, 
            file = "~/drive/MyDrive/long_sub.csv",
            quote = F, sep = ",", row.names = F, col.names = T)

save(dff, 
     file = "~/drive/MyDrive/df_bl.RData")

save(mdf.clean,
     file = "~/drive/MyDrive/long_sub.RData")
