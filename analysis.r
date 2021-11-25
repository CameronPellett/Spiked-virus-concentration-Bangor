remove(list = ls())

#packages
################
library(tidyverse)
library(patchwork)
################

#functions
###############
#summarises mean, sd, se, ci for bar charts
mean_ci_summary <- function(data, groups, variable, log10 = FALSE){
  
  data <- data|>ungroup()
  
  for (i in 1:length(groups)) {
    
    data <- data|>
      group_by(.data[[ groups[i] ]], .add = TRUE)
    
  }
  
  if(log10 == TRUE){
    
    data|>
      summarise(mean = mean(log10(.data[[variable]]), na.rm = TRUE),
                sd = sd(log10(.data[[variable]]), na.rm = TRUE),
                n = n())|>
      mutate(se = sd / sqrt(n),
             ci = se * qt(1 - 0.05 / 2, n-1))
    
  }else {
    
    data|>
      summarise(mean = mean(.data[[variable]], na.rm = TRUE),
                sd = sd(.data[[variable]], na.rm = TRUE),
                n = n())|>
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

stat_paste_fun <- function(stat){
  paste0(str_remove_all(stat[["method"]][[1]], "\n\t"), ": ", 
         names(stat[["statistic"]]), " = ",
         round(stat[["statistic"]][[1]],2), ", ",
         if(length(stat[["parameter"]]) == 0)
         {""
         }else if(is.na(stat[["parameter"]])){
           ""
         }else if(length(stat[["parameter"]]) > 0){
           paste0("df = ", round(stat[["parameter"]], 2), ", ")
         } ,
         "p-value ",  pvalue_fun(stat[["p.value"]]))
}
###############

#cleanup
################
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
  filter(!virus %in% c("AdV", "PMMoV", "CrAss"))|>
  filter(!method %in% c("Amicon", "Pellet", "Pellet_BE"))|>
  left_join(virus_char, by = "virus")|>
  mutate(spike = str_replace_all(spike, c("U" = "unspiked",
                                          "S" = "spiked")))|>
  pivot_wider(names_from = "spike", values_from = "rec", names_prefix = "rec_")|>
  rowwise()|>
  mutate(rec = case_when(is.na(rec_spiked) | is.na(rec_unspiked) ~ rec_spiked,
                                !is.na(rec_spiked) & !is.na(rec_unspiked) ~ rec_spiked - rec_unspiked))|>
  group_by(method, water_type, volume, virus)|>
  summarise(rec = mean(rec, na.rm = TRUE),
            host = host[[1]],
            enveloped = enveloped[[1]],
            genome = genome[[1]],
            size_nm = as.double(size_nm[[1]]),
            shape = shape[[1]])|>
  mutate(rec = case_when(rec <= 0 ~ NA_real_,
                         TRUE ~ rec))


spike_ami_pel <- data|>
  select(-Sample_code)|>
  rename("method" = Method, "spike" = Spike,
         "water_type" = Water_type,
         "volume" = Volume, "replicate" = Replicate)|>
  pivot_longer(cols = contains("rec"), names_to = "virus", values_to = "rec")|>
  mutate(virus = str_remove_all(virus, "_Rec%"))|>
  filter(!virus %in% c("AdV", "PMMoV", "CrAss"))|>
  left_join(virus_char, by = "virus")|>
  mutate(spike = str_replace_all(spike, c("U" = "unspiked",
                                          "S" = "spiked")))|>
  pivot_wider(names_from = "spike", values_from = "rec", names_prefix = "rec_")|>
  rowwise()|>
  mutate(rec = case_when(is.na(rec_spiked) | is.na(rec_unspiked) ~ rec_spiked,
                         !is.na(rec_spiked) & !is.na(rec_unspiked) ~ rec_spiked - rec_unspiked))|>
  group_by(method, water_type, volume, virus)|>
  summarise(rec = mean(rec, na.rm = TRUE),
            host = host[[1]],
            enveloped = enveloped[[1]],
            genome = genome[[1]],
            size_nm = as.double(size_nm[[1]]),
            shape = shape[[1]])|>
  mutate(rec = case_when(rec <= 0 ~ NA_real_,
                         TRUE ~ rec))

#################

#overview stats
####################
anova(lm(log10(rec) ~ volume + method + enveloped + size_nm + shape + water_type, spike))
#summary(lm(log10(rec) ~  volume  + enveloped  + method + size_nm + shape + water_type, spike))

summary(lm(log10(rec) ~  volume + method + size_nm + enveloped + water_type + shape, spike))

####################

#DW WW
#############################################
spike|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~virus + water_type, scales = "free_y")

spike|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~water_type, scales = "free_y")

#t.test p.value for all virus comp
all.t.p.v <- tibble(p.value = paste0("p-value ",
                                pvalue_fun(t.test(rec ~ water_type, data = spike)[["p.value"]])),
                    water_type = "WW"
)

#t.test p.value comp for individual virus comp
viruses <- spike|>ungroup()|>count(virus)|>pull(virus)
t.p.v <- tibble(virus = viruses, p.value = NA_character_, water_type = "WW")
for (i in 1:length(viruses)) {
  
t.p.v[i,2] <- paste0("p-value ",
                     pvalue_fun(t.test(rec ~ water_type, data = (spike|>filter(virus == viruses[i])) )[["p.value"]])
                     )
  
}

#t.test p.value comp for individual shapes 
shapes <- spike|>ungroup()|>count(shape)|>pull(shape)
t.p.v.s <- tibble(shape = shapes, p.value = NA_character_, water_type = "WW")
for (i in 1:length(shapes)) {
  
  t.p.v.s[i,2] <- paste0("p-value ",
                       pvalue_fun(t.test(rec ~ water_type, data = (spike|>filter(shape == shapes[i])) )[["p.value"]])
  )
  
}

#t.test p.value comp for individual envelopes
envelopes <- spike|>ungroup()|>count(enveloped)|>pull(enveloped)
t.p.v.e <- tibble(enveloped = envelopes, p.value = NA_character_, water_type = "WW")
for (i in 1:length(envelopes)) {
  
  t.p.v.e[i,2] <- paste0("p-value ",
                         pvalue_fun(t.test(rec ~ water_type, data = (spike|>filter(enveloped == envelopes[i])) )[["p.value"]])
  )
  
}

#t.test p.value comp for individual methods
methods <- spike|>ungroup()|>count(method)|>pull(method)
t.p.v.m <- tibble(method = methods, p.value = NA_character_, water_type = "WW")
for (i in 1:length(methods)) {
  
  t.p.v.m[i,2] <- paste0("p-value ",
                         pvalue_fun(t.test(rec ~ water_type, data = (spike|>filter(method == methods[i])) )[["p.value"]])
  )
  
}


water_type <- spike|>
    na.omit()|>
    mean_ci_summary(groups = "water_type", variable = "rec", log10 = TRUE)|>
    bind_cols(virus = "All viruses")|>
    left_join(all.t.p.v, by = "water_type")|>
  ggplot(aes(water_type, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
    geom_text(aes(label = p.value, x = 1.5, y = 1.05), size = 3)+
    geom_segment(aes(x = 1, xend = 2, y = 1, yend = 1))+
    geom_segment(aes(x = 1, xend = 1, y = 1, yend = .95))+
    geom_segment(aes(x = 2, xend = 2, y = 1, yend = .95))+
  labs(y = "Mean log10 recovery",
       x = "Water type")+
    theme(axis.title.x = element_blank())

water_type_shape <- spike|>
  na.omit()|>
  mean_ci_summary(groups = c("water_type", "shape"), variable = "rec", log10 = TRUE)|>
  left_join(t.p.v.s, by = c("shape", "water_type"))|>
  ggplot(aes(water_type, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  geom_text(aes(label = p.value, x = 1.5, y = 1.2), size = 3)+
  geom_segment(aes(x = 1, xend = 2, y = 1.15, yend = 1.15))+
  geom_segment(aes(x = 1, xend = 1, y = 1.15, yend = 1.1))+
  geom_segment(aes(x = 2, xend = 2, y = 1.15, yend = 1.1))+
  facet_grid(~shape, scales = "free_y")+
  labs(y = "Mean log10 recovery",
       x = "Water type")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())

water_type_envelope <- spike|>
  na.omit()|>
  mean_ci_summary(groups = c("water_type", "enveloped"), variable = "rec", log10 = TRUE)|>
  left_join(t.p.v.e, by = c("enveloped", "water_type"))|>
  mutate(enveloped = str_replace_all(enveloped, c("no" = "Not enveloped",
                                                  "yes" = "Enveloped")))|>
  ggplot(aes(water_type, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  geom_text(aes(label = p.value, x = 1.5, y = 1.2), size = 3)+
  geom_segment(aes(x = 1, xend = 2, y = 1.15, yend = 1.15))+
  geom_segment(aes(x = 1, xend = 1, y = 1.15, yend = 1.1))+
  geom_segment(aes(x = 2, xend = 2, y = 1.15, yend = 1.1))+
  facet_grid(~enveloped, scales = "free_y")+
  labs(y = "Mean log10 recovery",
       x = "Water type")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())

water_type_method <- spike|>
  na.omit()|>
  mean_ci_summary(groups = c("water_type", "method"), variable = "rec", log10 = TRUE)|>
  left_join(t.p.v.m, by = c("method", "water_type"))|>
  ggplot(aes(water_type, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  geom_text(aes(label = p.value, x = 1.5, y = 1.5), size = 3)+
  geom_segment(aes(x = 1, xend = 2, y = 1.45, yend = 1.45))+
  geom_segment(aes(x = 1, xend = 1, y = 1.45, yend = 1.4))+
  geom_segment(aes(x = 2, xend = 2, y = 1.45, yend = 1.4))+
  facet_grid(~method, scales = "free_y")+
  labs(y = "Mean log10 recovery",
       x = "Water type")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())

water_type_virus_2 <- spike|>
    filter(virus %in% c("Flu-A", "Flu-B"))|>
    na.omit()|>
    mean_ci_summary(groups = c("water_type", "virus"), variable = "rec", log10 = TRUE)|>
    left_join(t.p.v, by = c("virus", "water_type"))|>
    ggplot(aes(water_type, mean))+
    geom_col()+
    geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
    geom_text(aes(label = p.value, x = 1.5, y = 2.0), size = 3)+
    geom_segment(aes(x = 1, xend = 2, y = 1.9, yend = 1.9))+
    geom_segment(aes(x = 1, xend = 1, y = 1.9, yend = 1.85))+
    geom_segment(aes(x = 2, xend = 2, y = 1.9, yend = 1.85))+
    facet_grid(~virus, scales = "free_y")+
    labs(y = "Mean log10 recovery",
         x = "Water type")+
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())

water_type_virus_7 <- spike|>
    filter(!virus %in% c("Flu-A", "Flu-B"))|>
    na.omit()|>
    mean_ci_summary(groups = c("water_type", "virus"), variable = "rec", log10 = TRUE)|>
    left_join(t.p.v, by = c("virus", "water_type"))|>
    ggplot(aes(water_type, mean))+
    geom_col()+
    geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
    geom_text(aes(label = p.value, x = 1.5, y = 2.0), size = 3)+
    geom_segment(aes(x = 1, xend = 2, y = 1.9, yend = 1.9))+
    geom_segment(aes(x = 1, xend = 1, y = 1.9, yend = 1.85))+
    geom_segment(aes(x = 2, xend = 2, y = 1.9, yend = 1.85))+
    facet_grid(~virus, scales = "free_y")+
    labs(y = "Mean log10 recovery",
         x = "Water type")

water_type /
  (water_type_shape + water_type_envelope)
#############################################

#volume
######################
spike|>
  mutate(volume = factor(as.character(volume), levels = c("15", "20", "37.5", "50", "150")))|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~volume + method, scales = "free")

volume <- spike|>
  mutate(volume = factor(as.character(volume), levels = c("15", "20", "37.5", "50", "150")))|>
  mean_ci_summary(groups = c("volume"), variable = "rec", log10 = TRUE)|> 
  ggplot(aes(volume, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  labs(y = "Mean log10 recovery",
       x = "Volume (ml)")

volume_method <- spike|>
  mutate(volume = factor(as.character(volume), levels = c("15", "20", "37.5", "50", "150")))|>
  mean_ci_summary(groups = c("volume", "method"), variable = "rec", log10 = TRUE)|> #, "water_type" # , "virus"
  ggplot(aes(volume, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  #facet_wrap(~virus, scales = "free")+
  #facet_wrap(~water_type, scales = "free")+
  facet_wrap(~method, scales = "free")+
  labs(y = "Mean log10 recovery",
       x = "Volume (ml)")

volume_virus <- spike|>
  mutate(volume = factor(as.character(volume), levels = c("15", "20", "37.5", "50", "150")))|>
  mean_ci_summary(groups = c("volume", "virus"), variable = "rec", log10 = TRUE)|> # , "virus"
  ggplot(aes(volume, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  facet_wrap(~virus, scales = "free", nrow = 2)+
  labs(y = "Mean log10 recovery",
       x = "Volume (ml)")

(volume + volume_method)
###############################

#method
######################
spike_ami_pel|>
  filter(volume == 15)|>
  mutate(method = fct_reorder(.f = method, .x = rec, .fun = median))|>
  ggplot(aes(method, log10(rec)))+
  geom_boxplot()
######################

#size
#############################

spike|>
  mutate(size = fct_reorder(as.character(size_nm), size_nm))|>
  ggplot(aes(size, log10(rec)))+
  geom_boxplot()

##############################

