Comparative assessment of filtration- and precipitation-based methods
for the concentration of SARS-CoV-2 and other viruses from wastewater
================
Cameron Pellett
17/03/2022

This document and github repo provides the r code for results (3.1)
produced in *Comparative assessment of filtration- and
precipitation-based methods for the concentration of SARS-CoV-2 and
other viruses from wastewater*. The full report with code is included
below. The source code for this document can be found in README.rmd. For
a cleaner report with a html output, Spiked-virus-data-report.rmd can be
used. qPCR recovery data can be found in Spiking_data.csv, virus
characteristics can be found in virus_characteristics.csv, and some
additional untidy notes and exploration can be found in analysis.r.
Results in section 3.2. can be found in another repo:
<https://github.com/CameronPellett/Chem-Con-Bangor-Bath>.

<br> <br>

``` r
#packages
library(tidyverse)
library(patchwork)
library(broom)
library(ggbeeswarm)

theme_set(theme_grey()+theme(panel.background = element_rect(fill = "grey96")))
```

``` r
#functions


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

#converts p-values to easily pasted groups
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

#converts p-values to star based indicators, similar to those seen in summary(lm())
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

#extracts results from a statistical test and reforms into one neat line
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

#as above but for an anova
paste_anova_fun <- function(anova, log = FALSE){
  
  if(log != TRUE){
  paste0("ANOVA: ", "F value = ", round((anova$`F value`)[1], 2), ", ", "p-value ", pvalue_fun(anova$`Pr(>F)`[1]))
  }else {
    paste0("ANOVA (log y): ", "F value = ", round((anova$`F value`)[1], 2), ", ", "p-value ", pvalue_fun(anova$`Pr(>F)`[1])) 
  }
}


#next three pull R2 and adjR2
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

glmR2 <- function(glm, adj = FALSE){
  
  residuals <- glm[["y"]] - glm[["fitted.values"]]
  
  R2 <- 1-((sum(residuals^2))/
             (sum((glm$model[[1]] - mean(glm$model[[1]]))^2)))
  
  if(adj != TRUE){
    
  return(R2)
    
  }else {
    
    1- ((1 - R2) * (length(glm$model[[1]] ) - 1)/
                 (length(glm$model[[1]] ) - length(glm$coefficients) - 1))
  }
}


#pastes significant pairwaise comparisons (not used in final edit)
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

#plots a model residuals to assess assumptions
plot_assumptions <- function(model, obs, bins = 30){
        #create tibble for use in ggplot of model assumptions
        df <- tibble(
                observations = obs,
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
```

``` r
#data import and clean

data <- read_csv("Spiking_data.csv")

virus_char <- read_csv("virus_characteristics.csv")|>
  pivot_longer(2:12, names_to = "virus", values_to = "value")|>
  pivot_wider(names_from = "...1", values_from = "value")|>
  mutate(virus = str_replace_all(virus, c("SARS" = "N1")))

clean <- data|>
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
         method = factor(method),
         genome = factor(genome),
         samples_spiked = str_replace_all(samples_spiked, c("yes" = "S",
                                                            "no" = "U")),
         method = str_replace_all(method, "Amicon", "AM"),
         virus = str_replace_all(virus, c("N1" = "SARS-CoV-2",
                                          "CrAss" = "crAss")))|>
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
         method = factor(method),
         virus = str_replace_all(virus, c("N1" = "SARS",
                                          "CrAss" = "crAss")))|>
  na.omit()|>
  mutate(BE = case_when(str_detect(method, "BE") ~ "Beef extract",
                        TRUE ~ "No beef extract"))
```

``` r
notes <- (tibble(estimate.x = "p-value: < 0.001 [***]; < 0.01 [**]; < 0.05 [*]; > 0.05 [.]"))
```

``` r
#supplimentary figure full dataset plot

ww_all <- clean|>
  filter(water_type == "WW")|>
  mutate(volume = fct_reorder(.f = factor(paste0(volume, "ml WW")), .x = volume, .fun = median),
         shape1 = case_when(method == "AM" & volume == "20ml" ~ as.integer(32),
                            virus %in% c("AdV", "PMMoV", "crAss") ~ as.integer(32),
                            TRUE ~ as.integer(16)),
         shape2 = case_when(method == "AM" & volume == "20ml" ~ as.integer(16),
                            TRUE ~ as.integer(32)),
         shape3 = case_when(virus %in% c("AdV", "PMMoV", "crAss") & volume != "20ml" ~ as.integer(16),
                            TRUE ~ as.integer(32)),
         virus = str_replace_all(virus, "SARS-CoV-2", "SARS"))|>
  ggplot(aes(y = rec, x =method))+
  geom_beeswarm(aes(shape = shape1), priority = "ascending", cex = 4, show.legend = F)+
  geom_beeswarm(aes(shape = shape2), priority = "ascending", cex = 12, show.legend = F)+
  geom_beeswarm(aes(shape = shape3), priority = "ascending", cex = 6, show.legend = F)+
  geom_beeswarm(aes(shape = shape1), priority = "ascending", cex = 4, show.legend = F, size = 0.1, colour = "white")+
  geom_beeswarm(aes(shape = shape2), priority = "ascending", cex = 12, show.legend = F, size = 0.1, colour = "white")+
  geom_beeswarm(aes(shape = shape3), priority = "ascending", cex = 6, show.legend = F, size = 0.1, colour = "white")+
  scale_shape_identity()+
  scale_y_continuous(sec.axis = sec_axis( trans=~.*1, name="Virus", labels = NULL, breaks = NULL))+
  facet_grid(samples_spiked + virus~ volume, scale = "free", space = "free_x")+
  labs(x = "Concentration method", y = "Recovery (%)")

dw_all <- clean|>
  filter(water_type == "DW")|>
  mutate(volume = fct_reorder(.f = factor(paste0(volume, "ml DW")), .x = volume, .fun = median),
         shape1 = case_when(method == "AM" & volume == "20ml" ~ as.integer(32),
                            virus %in% c("AdV", "PMMoV", "crAss") ~ as.integer(32),
                            TRUE ~ as.integer(16)),
         shape2 = case_when(method == "AM" & volume == "20ml" ~ as.integer(16),
                            TRUE ~ as.integer(32)),
         shape3 = case_when(virus %in% c("AdV", "PMMoV", "crAss") & volume != "20ml" ~ as.integer(16),
                            TRUE ~ as.integer(32)),
         virus = str_replace_all(virus, "SARS-CoV-2", "SARS"))|>
  ggplot(aes(y = rec, x =method))+
  geom_beeswarm(aes(shape = shape1), priority = "ascending", cex = 4, show.legend = F)+
  geom_beeswarm(aes(shape = shape2), priority = "ascending", cex = 12, show.legend = F)+
  geom_beeswarm(aes(shape = shape3), priority = "ascending", cex = 6, show.legend = F)+
  geom_beeswarm(aes(shape = shape1), priority = "ascending", cex = 4, show.legend = F, size = 0.1, colour = "white")+
  geom_beeswarm(aes(shape = shape2), priority = "ascending", cex = 12, show.legend = F, size = 0.1, colour = "white")+
  geom_beeswarm(aes(shape = shape3), priority = "ascending", cex = 6, show.legend = F, size = 0.1, colour = "white")+
  scale_shape_identity()+
  scale_y_continuous(sec.axis = sec_axis( trans=~.*1, name="Virus", labels = NULL, breaks = NULL))+
  facet_grid(samples_spiked + virus~ volume, scale = "free", space = "free_x")+
  labs(x = "Concentration method", y = "Recovery (%)")

WWDW_all <- (ww_all + theme(axis.title.x = element_blank())) / dw_all + plot_layout(heights = c(1.1, 0.8))
```

## 2. Methods

### 2.1 Statistical methodology

Data analysis and statistical tests were carried out in R (R Core Team,
2020); the full script and data are provided in a dedicated repository
(<https://github.com/CameronPellett/Spiked-virus-concentration-Bangor>).

Sample starting volume (15ml, 30ml, 37.5ml, and 150ml), concentration
method (PEG, AS, IP, BE-PEG, and AM), and water type (wastewater WW, and
deionised water DW) were selected as factors and co-variates of viral
recovery. The selected features were assessed in combination as
predictors of viral recovery using a generalised linear model (glm) with
the response variable modelled as a Gamma distribution with a
logarithmic link function (Eq. 1; see supplementary materials Figure S2
for residuals plots).

$$  
Recovery\_{i} = exp\[\\beta\_{0} + \\beta\_{1} \\space Volume\_{i} + \\beta\_{2} \\space Method\_{i} + \\beta\_{3} \\space Water \\space Type\_{i} + \\epsilon\] \\quad \\quad (Eq.1) 
$$

Following the glm, comparisons of individual features were visualised,
and statistical tests were carried out. For statistical tests the
recovery percentile was log transformed to meet assumptions of a
Gaussian distribution (see supplementary materials Figure S3-S6 for
quantile-quantile plots). Equality of variances were tested with F
tests. Statistical comparisons of features with two levels and non-equal
variance were made with Welch two sample t-tests. Comparisons with three
or more levels with equal variance were made with a one-way ANOVA,
followed by pairwise two sample t-tests with pooled standard deviations
(SD), adjusting p-values with the holm method. Comparisons with three or
more levels and non-equal variance were made with a Welch ANOVA (one way
comparison of means), followed by pairwise two sample t-tests without
pooled SD, adjusting p-values with the holm method. Paired tests were
not selected due to missing data created by removal of undetermined
results and sample removal during qPCR quality control. Water type was
compared with only spiked viruses due to the lack of naturally present
virus in deionised water. Starting volume was compared using only
wastewater to simulate typical interactions seen in wastewater projects
and not using the AM or IP methods due to differing volumes. Method
comparison was made using 15ml wastewater sample volumes which were used
in all concentration methods, including AM and IP.

<br>

## 3. Results

The generalised linear model identified sample starting volume,
concentration method, and water type as significant predictors of viral
recovery (Table 1). The model coefficients reveal that sample starting
volume is a significant negative correlate of viral recovery; the AM
method has the highest mean recovery; and wastewater has reduced
recovery over deionised water.

<br>

Table 1: generalised linear model with Gamma residuals (link = log)
predicting recovery of 11 viruses. Pellet methods were excluded from the
model due to non-standard starting volumes. The model results include
the variable coefficient on a log scale indicating its effect (positive
numbers indicating increased recovery), followed by the standard error
in brackets, and a significance code ( `` `r notes` `` p-value: \< 0.001
\[\*\*\*\]; \< 0.01 \[\*\*\]; \< 0.05 \[\*\]; \> 0.05 \[.\] ) rounded to
three decimal places. The Amicon (AM) method and deionised water (DW)
are included in the baseline intercept. Compared methods include amonium
sulphate (AS), polyethene glycol (PEG), beef extract PEG (BE-PEG), and
innova prep (IP).

``` r
#overview glm

glm <- glm(rec ~ volume + method + water_type, data = clean, family = Gamma(link = "log"))

lm_table <- tidy(glm)|>
  na.omit()|>
  mutate(term = str_replace_all(term, c("[(]Intercept[)]" = "Intercept",
                                        "volume" = "Volume",
                                        "methodAS" = "Method: AS",
                                        "methodBE-PEG" = "Method: BE-PEG",
                                        "methodIP" = "Method: IP",
                                        "methodPEG" = "Method: PEG",
                                        "size_nm" = "Virus genome size (nm)",
                                        "water_typeWW" = "Water type: WW",
                                        "shapespherical:" = "Virus shape-genome: spherical-",
                                        "shapehelical:" = "Virus shape-genome: helical-",
                                        "shapeicosahedral:" = "Virus shape-genome: icosahedral-",
                                        "shapeHead-tail:" = "Virus shape-genome: Head-tail-",
                                        "envelopedyes" = "Enveloped: TRUE",
                                        "genomeRNA" = "RNA",
                                        "genomeDNA" = "DNA")),
         across(2:4, ~round(.x, 3)))|>
  rowwise()|>
  mutate(p.value = pvalue_star(p.value))|>
  unite(col = "estimate", c("estimate", "std.error"), sep = " (")|>
  unite(col = "estimate", c(2,4), sep = ") ")|>
  select(-statistic)

aicr2 <- tibble(
  term = c("AIC", "R-squared", "Adj. R-squared"), 
  estimate = c(AIC(glm), glmR2(glm), glmR2(glm, adj = TRUE)),
  #estimate.y = c(AIC(glm2), glmR2(glm2), glmR2(glm2, adj = TRUE))
)|>
  mutate(across(contains("estimate"), ~as.character(round(.x, 3))))



lm_table_fin <- lm_table|>
  #left_join(lm_table, by = "term")|>
  bind_rows(tibble(term = ""))|>
  bind_rows(aicr2)|>
  mutate(across(everything(), ~replace_na(.x, "")))

knitr::kable(lm_table_fin, col.names = c("Variables", "Model estimates"))
```

| Variables      | Model estimates       |
|:---------------|:----------------------|
| Intercept      | 3.759 (0.078) \*\*\*  |
| Volume         | -0.007 (0.001) \*\*\* |
| Method: AS     | -1.245 (0.092) \*\*\* |
| Method: BE-PEG | -0.368 (0.092) \*\*\* |
| Method: IP     | -1.55 (0.101) \*\*\*  |
| Method: PEG    | -1.569 (0.094) \*\*\* |
| Water type: WW | -0.848 (0.057) \*\*\* |
|                |                       |
| AIC            | 4748.24               |
| R-squared      | 0.597                 |
| Adj. R-squared | 0.594                 |

<br>

### 3.1. Greater recovery of spiked virus in deionised water compared to watewater

``` r
#water type comp

wt_comp <- clean|>
  filter(spike == "S" & samples_spiked == "S")

wt_comp_t.t <- tibble(label = stat_paste_fun(t.test(log(rec) ~ water_type, data = wt_comp), log = TRUE),
       row = 1)

f <- var.test(log(rec) ~ water_type, data = wt_comp)

Figure1 <- (wt_comp|>
  mutate(row = row_number())|>
  left_join(wt_comp_t.t, by = "row")|>
  ggplot(aes(water_type, rec))+
  geom_boxplot(width = 2, outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.5)+
    geom_beeswarm(shape = ".", colour = "white")+
  geom_text(aes(x = 1.5, y = max(rec) + (0.07 * (max(rec) - min(rec))),
                label = label), size = 3)+
  labs(x = "Water type",
       y = "Recovery (%)")) +
  
  (wt_comp|>
  ggplot(aes(water_type, rec))+
  geom_boxplot(width = 2, outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.3, cex = 3, size = 1)+
    geom_beeswarm(cex = 3, shape = ".", colour = "white")+
  facet_wrap(~virus, ncol = 4)+
  labs(x = "Water type",
       y = "Recovery (%)") +
  theme(axis.title.y = element_blank())) +
  
  plot_layout(widths = c(1.3, 0.7), nrow = 2, ncol = 1) +
    plot_annotation(tag_levels = "a")
```

Deionised water (DW) had greater median recovery of spiked virus
(`` `r median(wt_comp|>filter(water_type == "DW")|>pull(rec))` `` 7.64%
) compared to wastewater
(`` `r median(wt_comp|>filter(water_type == "WW")|>pull(rec))` `` 5.23%
), which was found to be significant when comparing mean log transformed
viral recovery (Figure 1a;
`` `r stat_paste_fun(t.test(log(rec) ~ water_type, data = wt_comp), log = TRUE)` ``
Welch Two Sample t-test (log y): t = 5.5, df = 578, p-value \< 0.001 ).
Greater recovery from deionised water was seen in all viruses except
Phi6 and RoV (Figure 1b). These results suggest wastewater contains
other chemicals or materials which reduce efficiency of concentration,
extraction, and qPCR.

<br>

``` r
Figure1
```

![](README_files/figure-gfm/fig1-1.png)<!-- -->

Figure 1: Greater recovery of spiked virus in deionised water (DW)
compared to wastewater (WW). Panel ???a??? combines all spiked virus
recovery results and panel ???b??? separates recovery by individual virus.
Comparisons were made with a Welch Two Sample t-test after log
transformation of recovery.

<br>

### 3.2. Viral recovery improved with a reduced starting volume of wastewater

``` r
#volume assess

clean_ww <- clean|>
  filter(water_type == "WW")

vol_ww <- clean_ww|>
  filter(!method %in% c("AM", "IP"))

glm_ww <- glm(rec ~ volume, data = vol_ww, family = Gamma(link = "log"))
lm_ww <- lm(log(rec) ~ as.factor(volume), data = vol_ww)

vol_comp_aov <- tibble(label = paste_anova_fun(anova = anova(lm_ww), log = TRUE),
                       row = 1)

f <- var.test(log(rec) ~ as.factor(volume), data = filter(vol_ww, volume != 150))
f <- var.test(log(rec) ~ as.factor(volume), data = filter(vol_ww, volume != 15))
f <- var.test(log(rec) ~ as.factor(volume), data = filter(vol_ww, volume != 37.5))

pairwise.vol.ttest <- pairwise.t.test(log(vol_ww$rec), as.factor(vol_ww$volume))[["p.value"]]|>
  as_tibble(rownames = "volx")|>
  pivot_longer(2:3, names_to = "voly", values_to = "p.v")

fig2a <- vol_ww|>
  mutate(volume = fct_reorder(factor(volume), volume),
         row = row_number())|>
  left_join(vol_comp_aov, by = "row")|>
  ggplot(aes(volume, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.5)+
  geom_beeswarm(shape = ".", colour = "white")+
  geom_text(aes(x = 2, y = max(rec) + (0.07 * (max(rec) - min(rec))),
                label = label), size = 3)+
  labs(x = "Volume (ml)",
       y = "Recovery (%)")
  
fig2b <- vol_ww|>
  mutate(volume = fct_reorder(factor(volume), volume))|>
  ggplot(aes(volume, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.3, cex = 3, size = 1)+
  geom_beeswarm(shape = ".", cex = 3, colour = "white")+
  facet_wrap(~virus, scales = "free")+
  labs(x = "Volume (ml)",
       y = "Recovery (%)")+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90)) 

fig2c <- pairwise.vol.ttest|>
  rowwise()|>
  mutate(p.v = pvalue_star(p.v))|>
  ggplot(aes(volx, voly))+
  geom_tile(fill = "grey97", colour = "white")+
  geom_text(aes(label = p.v))+
  labs(x = "Volume (ml)")+
  theme(axis.title.y = element_blank())

fig2 <- fig2a/ ((fig2b + (fig2c/plot_spacer()/plot_spacer())) +
  plot_layout(widths =  c(1.35, 0.25))) + plot_layout( heights = c(1.3, 0.7)) + plot_annotation(tag_levels = "a")
```

The lowest starting volume (15ml) had the greatest median viral recovery
( `` `r vol_ww|>filter(volume == 15)|>pull(rec)|>median()` `` 4.85% )
followed by 37.5ml (
`` `r vol_ww|>filter(volume == 37.5)|>pull(rec)|>median()` `` 3.84% ),
then the largest volume ( 150ml;
`` `r vol_ww|>filter(volume == 150)|>pull(rec)|>median()` `` 1.75%;
Figure 2a ). The mean log transformed viral recovery was significantly
different between groups (
`` `r paste_anova_fun(anova = anova(lm_ww), log = TRUE)` `` ANOVA (log
y): F value = 36.25, p-value \< 0.001 ), as well as all pairwise
comparisons (Figure 2c; Pairwise t-tests with pooled SD; p-value \<
0.05; holm adjustment method). The negative trend between volume and
recovery was seen in the recovery of all individual viruses; Although,
MNV and RoV did have higher median recovery for 37.5ml samples compared
to 15ml (Figure 2b). These results support conclusions in 3.1, that
wastewater contains inhibitors of viral recovery, which are increased
with a greater volume.

<br>

``` r
fig2
```

![](README_files/figure-gfm/fig2-1.png)<!-- -->

Figure 2: Viral recovery improved with a reduced starting volume of
wastewater. Panel ???a??? combines all recovery results, panel ???b??? separates
recovery by virus with a variable y scale, and panel ???c??? shows p-values
(holm adjustment method) of pairwise comparisons with two sample t-tests
with pooled standard deviations ( `` `r notes` `` p-value: \< 0.001
\[\*\*\*\]; \< 0.01 \[\*\*\]; \< 0.05 \[\*\]; \> 0.05 \[.\] ).
Comparisons were made with an ANOVA after log transformation of
recovery, followed by pairwise t tests; Pairwise comparisons found
significant differences between all volumes with the holm adjustment
method (p-value \< 0.05).

<br>

### 3.3. BE-PEG and Amicon concentration methods have the greatest viral recovery

``` r
#method comp

clean_ww_15 <- clean_ww|>
  filter(volume == 15)

lm_meth <- lm(log(rec) ~ method, data = clean_ww_15)

f <- var.test(log(rec) ~ method, data = filter(clean_ww_15, method %in% c("BE-PEG", "AM")))
f2 <- var.test(log(rec) ~ method, data = filter(clean_ww_15, method %in% c("AS", "AM")))
f3 <- var.test(log(rec) ~ method, data = filter(clean_ww_15, method %in% c("IP", "AM")))
f4 <- var.test(log(rec) ~ method, data = filter(clean_ww_15, method %in% c("PEG", "AM")))


meth_comp_aov <- tibble(label = stat_paste_fun(oneway.test(log(rec) ~ method, data = clean_ww_15), log = TRUE)[1],
                       row = 1)

meth_pairwise.t.test <- pairwise.t.test(log(clean_ww_15$rec), as.factor(clean_ww_15$method), pool.sd = FALSE)[["p.value"]]|>
  as_tibble(rownames = "varx")|>
  pivot_longer(2:5, names_to = "vary", values_to = "p.v")

fig3a <- clean_ww_15|>
  mutate(method = fct_reorder(method, rec),
         row = row_number())|>
  left_join(meth_comp_aov, by = "row")|>
  ggplot(aes(method, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.5)+
  geom_beeswarm(shape = ".", colour = "white")+
  geom_text(aes(x = 3, y = max(rec) + (0.07 * (max(rec) - min(rec))),
                label = label), size = 3)+
  labs(x = "Concentration method",
       y = "Recovery (%)") 
  
fig3b <- clean_ww_15|>
  mutate(method = fct_reorder(str_replace_all(method, c("BE-PEG" = "BE-P")), rec))|>
  ggplot(aes(method, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.3, cex = 3, size = 1)+
  geom_beeswarm(shape = ".", cex = 3, colour = "white")+
  facet_wrap(~virus, scale = "free")+
  labs(x = "Concentration method",
       y = "Recovery (%)") 
  
fig3c <- meth_pairwise.t.test|>
  rowwise()|>
  mutate(p.v = pvalue_star(p.v),
         varx = str_replace_all(varx, c("BE-PEG" = "BE-P")),
         vary = str_replace_all(vary, c("BE-PEG" = "BE-P")))|>
  ggplot(aes(varx, vary))+
  geom_tile(fill = "grey97", colour = "white")+
  geom_text(aes(label = p.v))+
  xlab("Method")+
  theme(axis.title.y = element_blank())

fig3 <- fig3a / (fig3b + (fig3c/ plot_spacer()) + plot_layout(widths = c(1.35, 0.25))) + plot_layout(heights = c(1.3, 0.7), nrow = 2, ncol = 1) + plot_annotation(tag_levels = "a")
```

Different concentration methods did have an effect of viral recovery
(Figure 3); AM and BE-PEG methods had the highest median recovery (
`` `r clean_ww_15|>filter(method == "AM")|>pull(rec)|>median()` ``
12.18% and
`` `r clean_ww_15|>filter(method == "BE-PEG")|>pull(rec)|>median()` ``
10.885%, respectively ), followed by IP (
`` `r clean_ww_15|>filter(method == "IP")|>pull(rec)|>median()` `` 5.1%
), then AS (
`` `r clean_ww_15|>filter(method == "AS")|>pull(rec)|>median()` `` 5.04%
) and PEG (
`` `r clean_ww_15|>filter(method == "PEG")|>pull(rec)|>median()` ``
2.29% ). The variance of the AM method was significantly different to
all other methods (p-value \< 0.05), so a Welch ANOVA was selected
(one-way analysis of means not assuming equal variances) which found
significant differences between the method???s mean log transformed
recovery (
`` `r stat_paste_fun(oneway.test(log(rec) ~ method, data = clean_ww_15), log = TRUE)[1]` ``
One-way analysis of means (not assuming equal variances) (log y): F =
13.7, df = 4, p-value \< 0.001 ). Post-hoc pairwise comparisons with
t-tests found that the AM method???s viral recovery was not significantly
different to the IP method due to AM???s greater variance, whereas BE-PEG
was significantly different to all other methods except AM (Figure 3c).
AM???s performance compared to other methods varied highly between
viruses, performing well for Flu, MNV and PMMoV, whilst poorly for AdV
and crAss (Figure 3b). BE-PEG, on the other hand, was consistently
better than other methods (excluding AM) for individual virus recovery.
These results suggest the BE-PEG method would have the greatest and most
consistent recovery of virus.

<br>

``` r
fig3
```

![](README_files/figure-gfm/fig3-1.png)<!-- -->

Figure 3: Concentration methods influence virus recovery at a starting
volume of 15ml. Panel ???a??? combines all recovery results with starting
volumes of 15ml, panel ???b??? separates recovery by virus with a variable y
scale, and panel ???c??? shows p-values of pairwise comparisons of method
recovery using t-tests without pooled standard deviations (
`` `r notes` `` p-value: \< 0.001 \[\*\*\*\]; \< 0.01 \[\*\*\]; \< 0.05
\[\*\]; \> 0.05 \[.\]; p-value holm adjustment method ). Comparisons
were made with an ANOVA after log transformation of recovery, followed
by pairwise t tests (c); BE-PEG and AM methods had the highest median
recovery, but due to AM???s greater variance, pairwise comparisons with IP
(third highest median recovery) were only significantly different for
BE-PEG (p-value \< 0.05; c).

<br>

### 3.4. Removal of the sample solid fraction via centrifugation is likely to increase viral recovery

``` r
#pellet assess

facets <- pellet|>ungroup()|>count(BE)|>pull(BE)
t.test.pellet <- tibble(BE = facets, label = NA_character_, row = 2)
for (i in 1:length(facets)) {
  
  t.test.pellet[i,2] <- stat_paste_fun(t.test(log10(rec) ~ method, data = (pellet|>filter(BE == facets[i])) ), log = TRUE)
  
}



fig4 <- pellet|>
  group_by(BE)|>
  mutate(row = row_number(),
         max = max(rec),
         min = min(rec),
         method = str_replace_all(method, "Pellet_BE", "BE-Pellet"))|>
  left_join(t.test.pellet, by = c("BE", "row"))|>
  ggplot(aes(method, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.5)+
  geom_beeswarm(shape = ".", colour = "white")+
  geom_text(aes(x = 1.5, y = max + (0.07 * (max - min)),
                label = label), size = 3)+
  facet_wrap(~BE, scales = "free")+
  labs(x = "Concentration method",
       y = "Recovery (%)") +
  theme(axis.title.x = element_blank())+ 
  
  pellet|>
  mutate(method = str_remove_all(method, "BE-|_BE"))|>
  ggplot(aes(method, rec))+
  geom_boxplot(outlier.alpha = 0)+
  geom_beeswarm(alpha = 0.3, cex = 3, size = 1)+
  geom_beeswarm(shape = ".", cex = 3, colour = "white")+
  facet_grid(BE~virus, scales = "free")+
  labs(x = "Concentration method",
       y = "Recovery (%)")+
  theme(axis.text.x = element_text(angle = 90)) +
  
  plot_layout(nrow = 2, heights = c(1.3,0.7)) + plot_annotation(tag_level = "a")
```

The solid fraction of wastewater samples removed as the initial pellet
after spiking had significantly lower mean log transformed viral
recovery compared to the concentrate while using BE-PEG and PEG
concentration methods ( Figure 4a; (BE-PEG)
`` `r t.test.pellet[[1,2]]` `` Welch Two Sample t-test (log y): t =
16.9, df = 147, p-value \< 0.001; (PEG) `` `r t.test.pellet[[2,2]]` ``
Welch Two Sample t-test (log y): t = 7, df = 55, p-value \< 0.001 ). The
pellet had consistently lower recovery of all individual viruses (Figure
4b). These results suggest a greater proportion of virus is suspended in
the liquid of a sample compared to the solid fraction. It is, thus,
likely that removal of the pellet will improve overall virus recovery.

<br>

``` r
fig4
```

![](README_files/figure-gfm/fig4-1.png)<!-- -->

Figure 4: The sample solid fraction (Pellet) has significantly lower
viral recovery than the concentrated sample; thus, removal via
centrifugation will likely increase the median viral recovery of a
concentrated sample. Panel ???a??? combines all recovery results with
starting volumes of 50ml and panel ???b??? separates recovery by virus with
a variable y scale. Comparisons were made with a Welch Two Sample t-test
after log transformation of recovery.

<br>

## 5. Supplementary materials

``` r
WWDW_all
```

![](README_files/figure-gfm/plot%20S1%20dwww_all-1.png)<!-- -->

Figure S1: recovery percentage for all 11 spiked (S) and unspiked (U)
viruses grouped by wastewater (WW), deionised water (DW), sample volume
(ml), and concentration method. Points have been offset on the
horizontal axis to avoid over plotting.

<br>

``` r
plot_assumptions(glm, glm[["y"]])
```

![](README_files/figure-gfm/S2%20glm%20residual%20plots-1.png)<!-- -->

Figure S2: Residual plot of the generalised linear model fitted as a
Gamma distribution with a logarithmic link function. Residual plots of
(a) homoscedasticity, (b) histogram normality, (c) quantile-quantile
normality, and (d) true observations against fitted values.

<br>

``` r
wt_comp|>
  ggplot(aes(sample = log(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~water_type)+
  labs(y = "log recovery",
       x = "theoretical")
```

![](README_files/figure-gfm/S3%20-1.png)<!-- -->

Figure S3: quantile-quantile plot of log transformed virus recovery
grouped by water type (De-ionised water DW and wastewater WW).

<br>

``` r
vol_ww|>
  mutate(volume = fct_reorder(paste0(volume, "ml"), volume))|>
  ggplot(aes(sample = log(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~volume)+
  labs(y = "log recovery",
       x = "theoretical")
```

![](README_files/figure-gfm/S4-1.png)<!-- -->

Figure S4: quantile-quantile plot of log transformed virus recovery
grouped starting volume of wastewater.

<br>

``` r
clean_ww_15|>
  ggplot(aes(sample = log(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~method)+
  labs(y = "log recovery",
       x = "theoretical")
```

![](README_files/figure-gfm/s5-1.png)<!-- -->

Figure S5: quantile-quantile plot of log transformed virus recovery
grouped by concentration method.

<br>

``` r
pellet|>
  ggplot(aes(sample = log(rec)))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~method)+
  labs(y = "log recovery",
       x = "theoretical")
```

![](README_files/figure-gfm/s6-1.png)<!-- -->

Figure S6: quantile-quantile plot of log transformed virus recovery
grouped by pellet and concentrated sample.
