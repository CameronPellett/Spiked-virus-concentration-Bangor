remove(list = ls())

#packages
################
library(tidyverse)
library(patchwork)


theme_set(theme_grey()+theme(panel.background = element_rect(fill = "grey96")))
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

paste_anova_fun <- function(anova){
  paste0("ANOVA: ", "F value = ", round((anova$`F value`)[1], 2), ", ", "p-value ", pvalue_fun(anova$`Pr(>F)`[1]))
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
  filter(!virus %in% c("AdV", "PMMoV", "CrAss"), spike == "S")|>
  filter(!method %in% c("Amicon", "Pellet", "Pellet_BE"))|>
  left_join(virus_char, by = "virus")|>
  mutate(rec = case_when(rec <= 0 ~ NA_real_,
                         TRUE ~ rec),
         size_nm = as.double(size_nm))


spike_ami_pel <- data|>
  select(-Sample_code)|>
  rename("method" = Method, "spike" = Spike,
         "water_type" = Water_type,
         "volume" = Volume, "replicate" = Replicate)|>
  pivot_longer(cols = contains("rec"), names_to = "virus", values_to = "rec")|>
  mutate(virus = str_remove_all(virus, "_Rec%| Rec%"))|>
  filter(!virus %in% c("AdV", "PMMoV", "CrAss"), spike == "S")|>
  left_join(virus_char, by = "virus")|>
  mutate(rec = case_when(rec <= 0 ~ NA_real_,
                         TRUE ~ rec),
         size_nm = as.double(size_nm))


unspiked <- data|>
  select(-Sample_code)|>
  rename("method" = Method, "spike" = Spike,
         "water_type" = Water_type,
         "volume" = Volume, "replicate" = Replicate)|>
  pivot_longer(cols = contains("rec"), names_to = "virus", values_to = "rec")|>
  mutate(virus = str_remove_all(virus, "_Rec%| Rec%"))|>
  filter(virus %in% c("NoVGII", "AdV", "PMMoV", "CrAss", "N1"), spike == "U", water_type == "WW")|>
  filter(!method %in% c("Amicon", "Pellet", "Pellet_BE"))|>
  left_join(virus_char, by = "virus")|>
  mutate(rec = case_when(rec <= 0 ~ NA_real_,
                         TRUE ~ rec),
         size_nm = as.double(size_nm))
  

unspiked_ami_pel <- data|>
  select(-Sample_code)|>
  rename("method" = Method, "spike" = Spike,
         "water_type" = Water_type,
         "volume" = Volume, "replicate" = Replicate)|>
  pivot_longer(cols = contains("rec"), names_to = "virus", values_to = "rec")|>
  mutate(virus = str_remove_all(virus, "_Rec%| Rec%"))|>
  filter(virus %in% c("NoVGII", "AdV", "PMMoV", "CrAss", "N1"), spike == "U", water_type == "WW")|>
  #filter(!method %in% c("Amicon", "Pellet", "Pellet_BE"))|>
  left_join(virus_char, by = "virus")|>
  mutate(rec = case_when(rec <= 0 ~ NA_real_,
                         TRUE ~ rec),
         size_nm = as.double(size_nm))
#################

#spiked
#overview stats
####################
anova(lm(log10(rec) ~ volume + method + size_nm + enveloped + shape + water_type, spike))
#summary(lm(log10(rec) ~  volume  + enveloped  + method + size_nm + shape + water_type, spike))

summary(lm(log10(rec) ~  volume + method + size_nm + water_type + shape, spike))

####################

#volume
######################
qq_volume <- spike|>
  mutate(volume = factor(as.character(volume), levels = c("15", "20", "37.5", "50", "150")))|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~volume, scales = "free")


anova_volume_lbl <- tibble(label = paste_anova_fun(anova(lm(log10(rec) ~ volume, spike))),
       volume = 37.5)

volume <- spike|>
  mean_ci_summary(groups = c("volume"), variable = "rec", log10 = TRUE)|>
  left_join(anova_volume_lbl, by = "volume")|>
  mutate(volume = fct_reorder(as.character(volume), volume))|>
  ggplot(aes(volume, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  geom_text(aes(y = 1.1, label = label), size = 3)+
  labs(y = "log10 recovery",
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
  labs(y = "log10 recovery",
       x = "Volume (ml)")

volume_virus <- spike|>
  mutate(volume = factor(as.character(volume), levels = c("15", "20", "37.5", "50", "150")))|>
  mean_ci_summary(groups = c("volume", "virus"), variable = "rec", log10 = TRUE)|> # , "virus"
  ggplot(aes(volume, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  facet_wrap(~virus, scales = "free", nrow = 2)+
  labs(y = "log10 recovery",
       x = "Volume (ml)")

(volume + volume_method)
###############################

#method
######################
spike_method_data <- spike_ami_pel|>
  filter(volume == 15)

qq_method <- spike_method_data|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~method, scales = "free_y")

stat <- anova(lm(log10(rec) ~ method, spike_method_data))
pairwise.t.test(log10(spike_method_data$rec), spike_method_data$method, paired = TRUE)

anova_lbl <- tibble(label = paste_anova_fun(stat), method = "AS")

method <- spike_method_data|>
  mean_ci_summary(c("method"), "rec", TRUE)|>
  left_join(anova_lbl, by = "method")|>
  mutate(method = fct_reorder(.f = method, .x = mean, .fun = median))|>
  ggplot(aes(method, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  geom_text(aes(label = label, y = 1.75), size = 3)+
  labs(y = "log10 recovery",
       x = "Method")

######################

#shape
##########################
qq_shape <- spike|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~shape, scales = "free_y")

shape.t.p.v <- tibble(p.value = paste0(
                                     #"p-value ",
                                     #pvalue_fun(t.test(log10(rec) ~ shape, data = spike)[["p.value"]])),
                                     stat_paste_fun(t.test(log10(rec) ~ shape, data = spike))),
                    shape = "spherical"
)

stat_paste_fun(t.test(log10(rec) ~ shape, data = spike))

shape <- spike|>
  mean_ci_summary(c("shape"), "rec", TRUE)|>
  left_join(shape.t.p.v, by = "shape")|>
  ggplot(aes(shape, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  geom_text(aes(label = p.value, x = 1.5, y = max(mean + ci) + max(ci * 1.6)), size = 3)+
  geom_segment(aes(x = 1, xend = 2, y = max(mean + ci) + max(ci) * 0.7, yend = max(mean + ci) + (max(ci) * 0.7)))+
  geom_segment(aes(x = 1, xend = 1, y = max(mean + ci) + max(ci) * 0.7, yend = max(mean + ci) + (max(ci) * 0.4)))+
  geom_segment(aes(x = 2, xend = 2, y = max(mean + ci) + max(ci) * 0.7, yend = max(mean + ci) + (max(ci) * 0.4)))+
  labs(y = "log10 recovery",
       x = "Shape")

spike_method_data|>
  mean_ci_summary(c("shape", "method"), "rec", TRUE)|>
  left_join(shape.t.p.v, by = "shape")|>
  ggplot(aes(shape, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  facet_wrap(~method, scales = "free_y")

#########################

#water type
#############################################
qq_water_type <- spike|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~water_type, scales = "free_y")



#t.test p.value for all virus comp
all.t.p.v <- tibble(p.value = paste0(
                                #"p-value ",
                                #pvalue_fun(t.test(log10(rec) ~ water_type, data = spike)[["p.value"]])),
                                stat_paste_fun(t.test(log10(rec) ~ water_type, data = spike))),
                    water_type = "WW"
)

#t.test p.value comp for individual virus comp
viruses <- spike|>ungroup()|>count(virus)|>pull(virus)
t.p.v <- tibble(virus = viruses, p.value = NA_character_, water_type = "WW")
for (i in 1:length(viruses)) {
  
t.p.v[i,2] <- paste0("p-value ",
                     pvalue_fun(t.test(log10(rec) ~ water_type, data = (spike|>filter(virus == viruses[i])) )[["p.value"]])
                     )
  
}

#t.test p.value comp for individual shapes 
shapes <- spike|>ungroup()|>count(shape)|>pull(shape)
t.p.v.s <- tibble(shape = shapes, p.value = NA_character_, water_type = "WW")
for (i in 1:length(shapes)) {
  
  t.p.v.s[i,2] <- paste0("p-value ",
                       pvalue_fun(t.test(log10(rec) ~ water_type, data = (spike|>filter(shape == shapes[i])) )[["p.value"]])
  )
  
}

#t.test p.value comp for individual envelopes
envelopes <- spike|>ungroup()|>count(enveloped)|>pull(enveloped)
t.p.v.e <- tibble(enveloped = envelopes, p.value = NA_character_, water_type = "WW")
for (i in 1:length(envelopes)) {
  
  t.p.v.e[i,2] <- paste0("p-value ",
                         pvalue_fun(t.test(log10(rec) ~ water_type, data = (spike|>filter(enveloped == envelopes[i])) )[["p.value"]])
  )
  
}

#t.test p.value comp for individual methods
methods <- spike|>ungroup()|>count(method)|>pull(method)
t.p.v.m <- tibble(method = methods, p.value = NA_character_, water_type = "WW")
for (i in 1:length(methods)) {
  
  t.p.v.m[i,2] <- paste0("p-value ",
                         pvalue_fun(t.test(log10(rec) ~ water_type, data = (spike|>filter(method == methods[i])) )[["p.value"]])
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
  labs(y = "log10 recovery",
       x = "Water type")#+
    #theme(axis.title.x = element_blank())

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
  labs(y = "log10 recovery",
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
  labs(y = "log10 recovery",
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
  labs(y = "log10 recovery",
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
    labs(y = "log10 recovery",
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
    labs(y = "log10 recovery",
         x = "Water type")

water_type /
  (water_type_shape + water_type_envelope)
#############################################

(volume + method + theme(axis.title.y = element_blank())) / 
  (shape + water_type + theme(axis.title.y = element_blank()))

#size
#############################

anova(lm(log10(rec) ~ as.character(size_nm), data = spike))

qq_size <- spike|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~size_nm, scales = "free_y")


genome_size <- spike|>
  mutate(size = as.character(size_nm))|>
  mean_ci_summary(groups = c("size", "size_nm"), variable = "rec", log10 = TRUE)|>
  ungroup()|>
  mutate(size = fct_reorder(.f = size, .x = size_nm, .fun = min))|>
  ggplot(aes(size, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean+ci), width = 0.2)


spike|>
  ggplot(aes(size_nm, log10(rec)))+
  geom_jitter()

summary(lm(log10(rec) ~ size_nm, spike))

##############################

#envelope
#############################
qq_envelope <- spike|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~enveloped, scales = "free_y")


enveloped <- spike|>
  mean_ci_summary(groups = c("enveloped"), variable = "rec", log10 = TRUE)|>
  ggplot(aes(enveloped, mean))+
  geom_col()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean+ci), width = 0.2)

##############################

genome_size + enveloped

#pellet
############################
pellet <- spike_ami_pel|>
 filter(method %in% c("PEG", "BE-PEG", "Pellet", "Pellet_BE"))|>
  mutate(BE = str_detect(method, "BE"),
         BE = case_when(BE == TRUE ~ "Beef extract",
                        BE == FALSE ~ "No beef extract"))
  
pellet|>
  ggplot(aes(sample =log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~method, scales = "free")

facets <- pellet|>ungroup()|>count(BE)|>pull(BE)
t.test.pellet <- tibble(BE = facets, stat = NA_character_, row = 2)
for (i in 1:length(facets)) {
  
  t.test.pellet[i,2] <- stat_paste_fun(t.test(log10(rec) ~ method, data = (pellet|>filter(BE == facets[i])) ))
  
}

pellet_fig <- pellet|>
  mean_ci_summary(groups = c("method", "BE"), variable = "rec", log10 = TRUE)|>
  mutate(method = str_replace_all(method, c("Pellet_BE" = "Pellet")))|>
  group_by(BE)|>
  mutate(row = row_number(),
         max = max(mean + ci) * (row - (max(row)-1)),
         max = case_when(max == 0 ~ NA_real_,
                         TRUE ~ max))|>
  left_join(t.test.pellet, by = c("BE", "row"))|>
  ggplot(aes(method, mean))+
  geom_boxplot()+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2)+
  facet_wrap(~BE, scales = "free")+
  geom_text(aes(label = stat, x = 1.5, y = max + ci * 1.2), size = 3)+
  geom_segment(aes(x = 1, xend = 2, y = max + (ci * 0.7), yend = max + (ci * 0.7)))+
  geom_segment(aes(x = 1, xend = 1, y = max + ci * 0.7, yend = max + (ci * 0.4)))+
  geom_segment(aes(x = 2, xend = 2, y = max + ci * 0.7, yend = max + (ci * 0.4)))+
  labs(y = "log10 recovery",
       x = "Method")

###########################

pellet_fig

#unspiked
#overview stats
####################
anova(lm(log10(rec) ~ volume  + method, unspiked))
#summary(lm(log10(rec) ~  volume  + enveloped  + method + size_nm + shape + water_type, unspiked))

summary(lm(log10(rec) ~  volume + method, unspiked))
####################

#volume
######################
qq_Uvolume <- unspiked|>
  mutate(volume = factor(as.character(volume), levels = c("15", "20", "37.5", "50", "150")))|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~volume, scales = "free")



kruskal_volume_lbl <- tibble(label = stat_paste_fun(kruskal.test(log10(rec) ~ volume, unspiked)),
                           row = 1)

Uvolume <- unspiked|>
  mutate(row = row_number())|>
  left_join(kruskal_volume_lbl, by = c("row"))|>
  mutate(volume = fct_reorder(as.character(volume), volume))|>
  ggplot(aes(volume, log10(rec)))+
  geom_boxplot()+
  geom_text(aes(x = 2, y = 1.5, label = label), size = 3)+
  labs(y = "log10 recovery",
       x = "Volume (ml)")

volume_method <- unspiked|>
  mutate(volume = fct_reorder(as.character(volume), volume))|>
  ggplot(aes(volume, log10(rec)))+
  geom_boxplot()+
  facet_wrap(~method, scales = "free_y")+
  labs(y = "log10 recovery",
       x = "Volume (ml)")

Uvolume_virus <- unspiked|>
  mutate(volume = fct_reorder(as.character(volume), volume))|>
  ggplot(aes(volume, log10(rec)))+
  geom_boxplot()+
  facet_wrap(~virus, scales = "free_y")+
  labs(y = "log10 recovery",
       x = "Volume (ml)")

(Uvolume / Uvolume_virus)
###############################

#method
######################
unspiked_method_data <- unspiked_ami_pel|>
  filter(volume == 15)

qq_Umethod <- unspiked_method_data|>
  ggplot(aes(sample = log10(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~method, scales = "free_y")

stat <- kruskal.test(log10(rec) ~ method, unspiked_method_data)
pairwise.wilcox.test(log10(unspiked_method_data$rec), unspiked_method_data$method, paired = TRUE)

kruskal_lbl <- tibble(label = stat_paste_fun(stat), row = 1)

Umethod <- unspiked_method_data|>
  mutate(row = row_number())|>
  left_join(kruskal_lbl, by = "row")|>
  ungroup()|>
  mutate(method = fct_reorder(.f = method, .x = log10(rec), .fun = function(.x){median(.x, na.rm = TRUE)}))|>
  ggplot(aes(method, log10(rec)))+
  geom_boxplot()+
  geom_text(aes(label = label, y = 1.5, x = 2.5), size = 3)+
  labs(y = "log10 recovery",
       x = "Method")

Umethod_virus <- unspiked_method_data|>
  ungroup()|>
  mutate(method = fct_reorder(.f = method, .x = log10(rec), .fun = function(.x){median(.x, na.rm = TRUE)}))|>
  ggplot(aes(method, log10(rec)))+
  geom_boxplot()+
  facet_wrap(~virus, scales = "free_y")+
  labs(y = "log10 recovery",
       x = "Method")

Umethod / Umethod_virus
######################

Uvolume / Umethod 

Uvolume + (Uvolume_virus + theme(axis.text.x = element_text(angle = 90), axis.title = element_blank())) +
  Umethod + (Umethod_virus + theme(axis.text.x = element_text(angle = 90), axis.title = element_blank())) +
  plot_layout(ncol = 2, nrow = 2, widths = c(1.2, 0.8, 1.2, 0.8))


#report
###########################

#figures
#1 spiked
(volume + method + theme(axis.title.y = element_blank())) / 
  (shape + water_type + theme(axis.title.y = element_blank()))
#-lower volume = better (low inhibitors) #-DW better than WW (low inhibitors)
#-amicon method best followed by BE-PEG
#-spherical viruses have lower recovery


#2 unspiked
Uvolume + (Uvolume_virus + theme(axis.text.x = element_text(angle = 90), axis.title = element_blank())) +
  Umethod + (Umethod_virus + theme(axis.text.x = element_text(angle = 90), axis.title = element_blank())) +
  plot_layout(ncol = 2, nrow = 2, widths = c(1.2, 0.8, 1.2, 0.8))
#


#3 pellet
pellet_fig
#-removal of the pellet will improve recovery

#############################


#supplimentary
########################
genome_size + enveloped

(qq_volume + qq_method)/
  (qq_shape + qq_water_type)

qq_Uvolume / 
  qq_Umethod

qq_size + qq_envelope
######################