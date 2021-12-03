#packages
################
library(tidyverse)
library(patchwork)
library(broom)
library(ggbeeswarm)
library(modelr)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_grey()+theme(panel.background = element_rect(fill = "grey96")))
################

#functions
#################

#summarises mean, sd, se, ci for bar charts
mean_ci_summary <- function(data, groups, variable, log10 = FALSE){
  
  data <- data|>ungroup()
  
  for (i in 1:length(groups)) {
    
    data <- data|>
      group_by(.data[[ groups[i] ]], .add = TRUE)
    
  }
  
  if(log10 == TRUE){
    
    data|>
      mutate(varna = case_when(is.na(.data[[variable]]) | is.nan(.data[[variable]]) | is.infinite(.data[[variable]]) ~ as.integer(1),
                               TRUE ~ NA_integer_))|>
      summarise(mean = mean(log10(.data[[variable]]), na.rm = TRUE),
                sd = sd(log10(.data[[variable]]), na.rm = TRUE),
                n = n() - sum(varna, na.rm = TRUE))|>
      mutate(se = sd / sqrt(n),
             ci = se * qt(1 - 0.05 / 2, n-1))
    
  }else {
    
    data|>
      mutate(varna = case_when(is.na(.data[[variable]]) | is.nan(.data[[variable]]) | is.infinite(.data[[variable]]) ~ as.integer(1),
                               TRUE ~ NA_integer_))|>
      summarise(mean = mean(.data[[variable]], na.rm = TRUE),
                sd = sd(.data[[variable]], na.rm = TRUE),
                n = n()  - sum(varna, na.rm = TRUE))|>
      mutate(se = sd / sqrt(n),
             ci = se * qt(1 - 0.05 / 2, n-1))
    
  }
}

pvalue_fun <- function(p.value){
  if(is.na(p.value)){
    NA_character_
  }else if(p.value > 0.05){
    "> 0.05"
  }else if(p.value < 0.05 & p.value > 0.01){
    "< 0.05"
  }else if(p.value < 0.01 & p.value > 0.001){
    "< 0.01"
  }else if(p.value < 0.001){
    "< 0.001"
  }else {
    "something went wrong"
  }
}

pvalue_star <- function(p.value){
  if(is.na(p.value)){
    NA_character_
  }else if(p.value > 0.05){
    "."
  }else if(p.value < 0.05 & p.value > 0.01){
    "*"
  }else if(p.value < 0.01 & p.value > 0.001){
    "**"
  }else if(p.value < 0.001){
    "***"
  }else {
    "something went wrong"
  }
}

stat_paste_fun <- function(stat, log = FALSE){
  if(log != TRUE){
  paste0(str_remove_all(stat[["method"]][[1]], "\n\t"), ": ", 
         str_remove_all(names(stat[["statistic"]]), "Kruskal-Wallis "), " = ",
         round(stat[["statistic"]][[1]],1), ", ",
         if(length(stat[["parameter"]]) == 0)
         {""
         }else if(is.na(stat[["parameter"]])){
           ""
         }else if(length(stat[["parameter"]]) > 0){
           paste0("df = ", round(stat[["parameter"]], 0), ", ")
         } ,
         "p-value ",  pvalue_fun(stat[["p.value"]]))
  }else{
    paste0(str_remove_all(stat[["method"]][[1]], "\n\t"), " (log y): ", 
           str_remove_all(names(stat[["statistic"]]), "Kruskal-Wallis "), " = ",
           round(stat[["statistic"]][[1]],1), ", ",
           if(length(stat[["parameter"]]) == 0)
           {""
           }else if(is.na(stat[["parameter"]])){
             ""
           }else if(length(stat[["parameter"]]) > 0){
             paste0("df = ", round(stat[["parameter"]], 0), ", ")
           } ,
           "p-value ",  pvalue_fun(stat[["p.value"]]))
  }
}

paste_anova_fun <- function(anova, log = FALSE){
  
  if(log != TRUE){
  paste0("ANOVA: ", "F value = ", round((anova$`F value`)[1], 2), ", ", "p-value ", pvalue_fun(anova$`Pr(>F)`[1]))
  }else {
    paste0("ANOVA (log y): ", "F value = ", round((anova$`F value`)[1], 2), ", ", "p-value ", pvalue_fun(anova$`Pr(>F)`[1])) 
  }
}

lmR2adj <- function(lm){
  
  R2 <- 1-((sum(residuals(lm)^2))/
             (sum((lm$model[[1]] - mean(lm$model[[1]]))^2)))
  
  1- ((1 - R2) * (length(lm$model[[1]] ) - 1)/
        (length(lm$model[[1]] ) - length(lm$coefficients) - 1))
}

lmR2 <- function(lm){
  
  1-((sum(residuals(lm)^2))/
       (sum((lm$model[[1]] - mean(lm$model[[1]]))^2)))
}

pairwise_diff_comps <- function(pairwise_pv_matrix){
  
  pairwise_diff <- as_tibble(pairwise_pv_matrix, rownames = "vars")|>
    pivot_longer(2:last_col(), names_to = "var.y", values_to = "p.value")|>
    na.omit()|>
    filter(p.value < 0.05)|>
    select(1,2)|>
    unite("comp", c(1,2), sep = " ~ ")
  
  comp_diffs <- ""
  
  for (i in 1:nrow(pairwise_diff)) {
    
    if(i == nrow(pairwise_diff)){
      comp_diffs <- paste0(c(comp_diffs, pairwise_diff[[i,1]]), collapse = ", and ")
    }else if(i != 1){
      comp_diffs <- paste0(c(comp_diffs, pairwise_diff[[i,1]]), collapse = ", ")
    }else{
      comp_diffs <- pairwise_diff[[i,1]]
    }
    
  }
  
  return(comp_diffs)
}

plot_assumptions <- function(model, obs, bins = 30){
  #create tibble for use in ggplot of model assumptions
  
  
  df <- tibble(
    observations = model,
    fitted = fitted(model),
    residuals = resid(model),
    standardised_residuals = sqrt(resid(model)^2),
    RMSresiduals = sqrt(mean(resid(model)^2)))
  
  
  a <- ggplot(df, aes(fitted, residuals))+
    geom_point()+
    geom_abline(aes(slope = 0, intercept = 0))
  #labs(title = "Heteroscedasticity", subtitle = "Check trends not captured by model")
  
  b <- ggplot(df, aes(residuals))+
    geom_histogram(bins = bins)
  
  c <- ggplot(df, aes(sample = residuals))+
    geom_qq()+
    geom_qq_line()#+
  #labs(title = "Quantile-quantile plot", subtitle = "check normality of residuals")
  
  d <- ggplot(df, aes(observations, fitted))+
    geom_point()+
    geom_abline()#+
  #labs(subtitle = "Check accuracy of model across whole dataset")
  
  plot <- a+b+c+d+plot_layout(ncol = 2, nrow = 2)+plot_annotation(tag_levels = "a")
  return(plot)
} 
##############################

#import and clean
#####################

data <- read_csv("Spiking_data.csv")

virus_char <- read_csv("virus_characteristics.csv")|>
  pivot_longer(2:12, names_to = "virus", values_to = "value")|>
  pivot_wider(names_from = "X1", values_from = "value")|>
  mutate(virus = str_replace_all(virus, c("SARS" = "N1")))

spike <- data|>
  select(-Sample_code)|>
  rename("method" = Method, "spike" = Spike,
         "water_type" = Water_type,
         "volume" = Volume, "replicate" = Replicate)|>
  pivot_longer(cols = contains("rec"), names_to = "virus", values_to = "rec")|>
  mutate(virus = str_remove_all(virus, "_Rec%| Rec%"))|>
  mutate(rec = case_when((spike == "U" & virus %in% c("N1", "Phi6", "NoVGII") )~ NA_real_,
                         spike == "S" & water_type == "DW" & virus %in% c("AdV", "PMMoV", "CrAss") ~ NA_real_,
                         TRUE ~ rec))|>
  filter(!method %in% c("Pellet", "Pellet_BE"))|>
  left_join(virus_char, by = "virus")|>
  mutate(rec = case_when(round(rec,2) <= 0 ~ NA_real_,
                         TRUE ~ round(rec,2)),
         shape = factor(shape),
         water_type = factor(water_type),
         enveloped = factor(enveloped),
         size_nm = as.double(size_nm),
         method = factor(method))|>
  na.omit()


pellet <- data|>
  select(-Sample_code)|>
  rename("method" = Method, "spike" = Spike,
         "water_type" = Water_type,
         "volume" = Volume, "replicate" = Replicate)|>
  pivot_longer(cols = contains("rec"), names_to = "virus", values_to = "rec")|>
  mutate(virus = str_remove_all(virus, "_Rec%| Rec%"))|>
  mutate(rec = case_when((spike == "U" & virus %in% c("N1", "Phi6", "NoVGII") )~ NA_real_,
                         spike == "S" & water_type == "DW" & virus %in% c("AdV", "PMMoV", "CrAss") ~ NA_real_,
                         TRUE ~ rec))|>
  filter(method %in% c("Pellet", "Pellet_BE", "PEG", "BE-PEG"))|>
  left_join(virus_char, by = "virus")|>
  mutate(rec = case_when(round(rec,2) <= 0 ~ NA_real_,
                         TRUE ~ round(rec,2)),
         shape = factor(shape),
         water_type = factor(water_type),
         enveloped = factor(enveloped),
         size_nm = as.double(size_nm),
         method = factor(method))|>
  na.omit()|>
  mutate(BE = case_when(str_detect(method, "BE") ~ "Beef extract",
                        TRUE ~ "No beef extract"))
##########################

#all plot
########################
ww_all <- spike|>
  filter(water_type == "WW")|>
  mutate(volume = fct_reorder(.f = factor(paste0(volume, "ml")), .x = volume, .fun = median),
         shape1 = case_when(method == "Amicon" & volume == "20ml" ~ as.integer(32),
                            virus %in% c("AdV", "PMMoV", "CrAss") ~ as.integer(32),
                            TRUE ~ as.integer(16)),
         shape2 = case_when(method == "Amicon" & volume == "20ml" ~ as.integer(16),
                            TRUE ~ as.integer(32)),
         shape3 = case_when(virus %in% c("AdV", "PMMoV", "CrAss") & volume != "20ml" ~ as.integer(16),
                            TRUE ~ as.integer(32)))|>
  ggplot(aes(y = rec, x =method))+
  geom_beeswarm(aes(shape = shape1), priority = "ascending", cex = 4, show.legend = F)+
  geom_beeswarm(aes(shape = shape2), priority = "ascending", cex = 12, show.legend = F)+
  geom_beeswarm(aes(shape = shape3), priority = "ascending", cex = 6, show.legend = F)+
  geom_beeswarm(aes(shape = shape1), priority = "ascending", cex = 4, show.legend = F, size = 0.1, colour = "white")+
  geom_beeswarm(aes(shape = shape2), priority = "ascending", cex = 12, show.legend = F, size = 0.1, colour = "white")+
  geom_beeswarm(aes(shape = shape3), priority = "ascending", cex = 6, show.legend = F, size = 0.1, colour = "white")+
  scale_shape_identity()+
  scale_y_continuous(sec.axis = sec_axis( trans=~.*1, name="Virus", labels = NULL, breaks = NULL))+
  facet_grid(spike + virus~ water_type + volume, scale = "free", space = "free_x")+
  labs(x = "Method", y = "Recovery (%)")

dw_all <- spike|>
  filter(water_type == "DW")|>
  mutate(volume = fct_reorder(.f = factor(paste0(volume, "ml")), .x = volume, .fun = median),
         shape1 = case_when(method == "Amicon" & volume == "20ml" ~ as.integer(32),
                            virus %in% c("AdV", "PMMoV", "CrAss") ~ as.integer(32),
                            TRUE ~ as.integer(16)),
         shape2 = case_when(method == "Amicon" & volume == "20ml" ~ as.integer(16),
                            TRUE ~ as.integer(32)),
         shape3 = case_when(virus %in% c("AdV", "PMMoV", "CrAss") & volume != "20ml" ~ as.integer(16),
                            TRUE ~ as.integer(32)))|>
  ggplot(aes(y = rec, x =method))+
  geom_beeswarm(aes(shape = shape1), priority = "ascending", cex = 4, show.legend = F)+
  geom_beeswarm(aes(shape = shape2), priority = "ascending", cex = 12, show.legend = F)+
  geom_beeswarm(aes(shape = shape3), priority = "ascending", cex = 6, show.legend = F)+
  geom_beeswarm(aes(shape = shape1), priority = "ascending", cex = 4, show.legend = F, size = 0.1, colour = "white")+
  geom_beeswarm(aes(shape = shape2), priority = "ascending", cex = 12, show.legend = F, size = 0.1, colour = "white")+
  geom_beeswarm(aes(shape = shape3), priority = "ascending", cex = 6, show.legend = F, size = 0.1, colour = "white")+
  scale_shape_identity()+
  scale_y_continuous(sec.axis = sec_axis( trans=~.*1, name="Virus", labels = NULL, breaks = NULL))+
  facet_grid(spike + virus~ water_type + volume, scale = "free", space = "free_x")+
  labs(x = "Method", y = "Recovery (%)")

ww_all / dw_all
#########################

#overview glm
##########################
glm <- glm(rec ~ volume + method + water_type + shape + enveloped + genome, data = spike, family = Gamma(link = "log"))
summary(glm)

glmR2(glm, adj = TRUE)
glm$model[[1]]
plot_assumptions(glm, glm[["y"]])
##########################

#water type
########################
wt_comp <- spike|>
  filter(spike == "S" & samples_spiked == "yes")

wt_comp_t.t <- tibble(label = stat_paste_fun(t.test(log(rec) ~ water_type, data = wt_comp), log = TRUE),
       row = 1)

wt_comp|>
  filter(rec < 120)|>
  mutate(row = row_number())|>
  left_join(wt_comp_t.t, by = "row")|>
  ggplot(aes(water_type, rec))+
  geom_boxplot(width = 2, outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.5)+
  geom_text(aes(x = 1.5, y = max(rec) + (0.07 * (max(rec) - min(rec))),
                label = label))+
  labs(x = "Water type",
       y = "Recovery (%)")

#######################

#volume
########################
spike_ww <- spike|>
  filter(water_type == "WW")

glm_ww <- glm(rec ~ volume, data = filter(spike_ww, !method %in% c("Amicon", "IP")), family = Gamma(link = "log"))
lm_ww <- lm(log(rec) ~ as.factor(volume), data = filter(spike_ww, !method %in% c("Amicon", "IP")))
summary(lm_ww)
plot_assumptions(glm_ww, glm_ww[["y"]])
plot_assumptions(lm_ww, lm_ww[["model"]][["log(rec)"]])

vol_comp_aov <- tibble(label = paste_anova_fun(anova = anova(lm_ww), log = TRUE),
                       row = 1)

spike_ww|>
  filter(!method %in% c("Amicon", "IP"))|>
  mutate(volume = fct_reorder(factor(volume), volume),
         row = row_number())|>
  left_join(vol_comp_aov, by = "row")|>
  ggplot(aes(volume, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.5)+
  geom_text(aes(x = 2, y = max(rec) + (0.07 * (max(rec) - min(rec))),
                label = label), size = 3)+
  labs(x = "Volume (ml)",
       y = "Recovery (%)")
#########################

#method
########################
spike_ww_15 <- spike_ww|>
  filter(volume == 15)

lm_meth <- lm(log(rec) ~ method, data = spike_ww_15)
meth_comp_aov <- tibble(label = paste_anova_fun(anova = anova(lm_meth), log = TRUE),
                       row = 1)

pairwise.t.test(log(clean_ww_15$rec), clean_ww_15$method)[["p.value"]]|>
  as_tibble(rownames = "rownames")|>
  pivot_longer(2:last_col(), names_to = "voly", values_to = "p.v")|>
  na.omit()|>
  rowwise()|>
  mutate(p.v = pvalue_star(p.v))


#pairwise plot fun
pairwise_plot_layer <- function(pairwise.test, ymin, ymax, xmin, xmax, x, y, scale, data){

pairwise <- pairwise.test[["p.value"]]|>
  as_tibble(rownames = "rownames")|>
  rowwise()|>
  mutate(across(2:last_col(), pvalue_star),
         across(2:last_col(), ~replace_na(.x, "")))

title = tibble(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, 
               x = x, y = y, scale = scale,
               ncol = (ncol(pairwise)-1)/2,
               labels = paste0("Pairwise ",
                              pairwise.test[["method"]],
                              " (", pairwise.test[["p.adjust.method"]],
                              " adj.)"))|>
  mutate(labely = y + ((ymax-ymin)* 1/(16/scale)),
         labelx = x + ((xmax-xmin)* ncol/(16/scale)))
  

ylabel <- tibble(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, 
       x = x, y = y, scale = scale,
       nrow = 1:nrow(pairwise),
       labels = pull(pairwise, rownames))|>
  mutate(labely = y - ((ymax-ymin)* nrow/(16/scale)),
         labelx = x)

xlabel <- tibble(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, 
                 x = x, y = y, scale = scale,
       ncol = 1:(ncol(pairwise)-1),
       labels = colnames(pairwise)[-1])|>
  mutate(labely = y,
         labelx = x + ((xmax-xmin)* ncol/(16/scale)))


xypvalue <- tibble(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, 
                   x = x, y = y, scale = scale,
       nrow = rep(1:nrow(pairwise), each = ncol(pairwise)-1),
       ncol = rep(1:(ncol(pairwise)-1), nrow(pairwise)),
       labels = NA_character_)|>
  mutate(labely = y- ((ymax-ymin)* nrow/(16/scale)),
         labelx = x + ((xmax-xmin)* ncol/(16/scale)))


for (i in 1:nrow(xypvalue)) {

 nrow <-  (xypvalue[i,])|>
    pull(nrow)
 ncol <-  (xypvalue[i,])|>
   pull(ncol)
 
 xypvalue[[i, 10]] <- pairwise[[nrow, (ncol+1)]]
}

pairwise_label <- ylabel|>bind_rows(xlabel, xypvalue, title)

ex.pairwise_label <- pairwise_label|>
  bind_rows(tibble(
    ymin = rep(NA_real_, (nrow(data)-nrow(pairwise_label)))
    ))

return(ex.pairwise_label)
}
#
pairwise.test <- pairwise.t.test(log(clean_ww_15$rec), clean_ww_15$method)

ex.pairwise_label <- pairwise_plot_layer(pairwise.test,
                                         ymin = 0, ymax = 40, xmin = 0, xmax = 6,
                                         x = 0.8, y = 35.5, scale = 0.6,
                                         data = spike_ww_15)

spike_ww_15|>
  bind_cols(ex.pairwise_label)|>
  mutate(method = fct_reorder(method, rec),
         row = row_number())|>
  left_join(meth_comp_aov, by = "row")|>
  ggplot(aes(method, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.5)+
  geom_tile(aes(x = median(labelx, na.rm = TRUE),
                y = median(labely, na.rm = TRUE)), fill = "white", colour = "black", width = 1.4, height = 10)+
  geom_text(aes(x = 3, y = max(rec) + (0.07 * (max(rec) - min(rec))),
                label = label), size = 3)+
  geom_text(aes(x = labelx, y = labely,
                label = labels), size = 3)+
  labs(x = "Concentration method",
       y = "Recovery (%)")


########################

#pellet
########################
facets <- pellet|>ungroup()|>count(BE)|>pull(BE)
t.test.pellet <- tibble(BE = facets, label = NA_character_, row = 2)
for (i in 1:length(facets)) {
  
  t.test.pellet[i,2] <- stat_paste_fun(t.test(log10(rec) ~ method, data = (pellet|>filter(BE == facets[i])) ), log = TRUE)
  
}



pellet|>
  group_by(BE)|>
  mutate(row = row_number(),
         max = max(rec),
         min = min(rec))|>
  left_join(t.test.pellet, by = c("BE", "row"))|>
  ggplot(aes(method, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.3)+
  geom_text(aes(x = 1.5, y = max + (0.07 * (max - min)),
                label = label), size = 3)+
  facet_wrap(~BE, scales = "free")+
  labs(x = "Concentration method",
       y = "Recovery (%)")
#########################

#virus characteristics 
###########################
glm_vc <- glm(rec ~ shape + enveloped + genome, data = spike_ww, family = Gamma(link = "log"))
summary(glm_vc)
plot_assumptions(glm_vc, glm_vc[["y"]])

shape.anova <- tibble(label = paste_anova_fun(anova(lm(log(rec) ~ shape, spike_ww)), log = TRUE),
                      row = 1)

shape <- spike_ww|>
  mutate(row = row_number(),
         shape = fct_reorder(shape, rec))|>
  left_join(shape.anova, by = "row")|>
  ggplot(aes(shape, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.2)+
  geom_text(aes(x = 2.5, y = max(rec) + (0.07 * (max(rec) - min(rec))),
                label = label), size = 3)+
  labs(x = "Virus shape",
       y = "Recovery (%)")

enveloped.t.test <- tibble(label = stat_paste_fun(t.test(log(rec) ~ enveloped, spike_ww), log = TRUE),
                           row = 1)

enveloped <- spike_ww|>
  mutate(row = row_number())|>
  left_join(enveloped.t.test, by = "row")|>
  ggplot(aes(enveloped, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.2)+
  geom_text(aes(x = 1.5, y = max(rec) + (0.07 * (max(rec) - min(rec))),
                label = label), size = 3)+
  labs(x = "Virus enveloped",
       y = "Recovery (%)")

genome.t.test <- tibble(label = stat_paste_fun(t.test(log(rec) ~ genome, spike_ww), log = TRUE),
                        row = 1)

genome <- spike_ww|>
  mutate(row = row_number())|>
  left_join(genome.t.test, by = "row")|>
  ggplot(aes(genome, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.2)+
  geom_text(aes(x = 1.5, y = max(rec) + (0.07 * (max(rec) - min(rec))),
                label = label), size = 3)+
  labs(x = "Viral genome",
       y = "Recovery (%)")

(enveloped + genome) / shape
###########################



summary(glm(rec ~ method, data = spike_ww_15, family = Gamma(link = "log")))

pairwise.t.test(log(spike_ww_15$rec), spike_ww_15$method)


glm_vc <- glm(rec ~ 0 + shape + enveloped + genome, data = spike_ww, family = Gamma(link = "log"))
summary(glm_vc)
plot_assumptions(glm_vc, glm_vc[["y"]])

spike|>
  filter(rec < 0.02)
