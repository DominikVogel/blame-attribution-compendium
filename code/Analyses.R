# Encoding: UTF-8
set.seed(20170824)


# 0.1 Packages ################################################################

# Load packages
library(broman) # for myround()
library(stargazer) # for tables
library(yarrr) # for pirateplot
library(data.table) # for setattr()
library(psych) # for Cronbach's Alpha
library(data.table) # for row names in correlation table
library(compute.es) # for effect sizes
library(weights) # for rounding p values
library(here)
library(tidyverse) # for everything


# 0.2 Functions ###############################################################
# format p values in APA format
pformat <- function(p) {
  paste0("p ",
         ifelse(p <=1 & p >= 0.1,
                paste0("= ",
                       weights::rd(p, 2)),
                ifelse(p < .1 & p >= 0.001,
                       paste0("= ",
                              weights::rd(p, 3)),
                       paste0("< .001"))))
}

# ANOVA output
anova_output <- function(anova) {
  anova <- summary(anova)
  output <- paste0("F(",
                   as.integer(anova[[1]]$Df[1]),
                   ",",
                   as.integer(anova[[1]]$Df[2]),
                   ") = ", 
                   myround(anova[[1]]$`F value`[1], 2),
                   ", ",
                   pformat(anova[[1]]$`Pr(>F)`[1]))
  return(output)
}

# Chi2 output
chi_output <- function(chi_result) {
  paste0("Chi2(",
         as.integer(chi_result$parameter),
         ") = ",
         myround(chi_result$statistic, 2),
         ", ",
         pformat(chi_result$p.value))
}

# Return Hedge's g and its CI
hedge_function <- function(treatment) {
  tc <- df %>% 
    dplyr::filter(group == "Control") %>% 
    dplyr::select(blame_mean)
  mc <- mean(tc$blame_mean, na.rm = TRUE)
  sdc <- sd(tc$blame_mean, na.rm = TRUE)
  nc <- nrow(tc)
  t1 <- df %>% 
    dplyr::filter(group == treatment) %>% 
    dplyr::select(blame_mean)
  m1 <- mean(t1$blame_mean, na.rm = TRUE)
  sd1 <- sd(t1$blame_mean, na.rm = TRUE)
  n1 <- nrow(t1)
  effect <- compute.es::mes(m1, mc, sd1, sdc, n1, nc, verbose = FALSE)
  g <- effect$g
  lg <- effect$l.g
  ug <- effect$u.g
  return(c(effect$g, effect$l.g, effect$u.g))
}

# 1. Read in data #############################################################
df <- read_csv(here("data", "blame_data.csv"))


# 2. Generate Variables #######################################################

# 2.1 Treatment as factor -----------------------------------------------------
df <- df %>% 
  mutate(group = factor(group,
                        labels = c("Control", 
                                   "Public opinion in favor of blaming",
                                   "Public opinion not in favor of blaming",
                                   "High intensity of pre-existing blame",
                                   "No pre-existing blame")))

# 2.2 Age ---------------------------------------------------------------------
df <- df %>% mutate(age = 2017 - yearbirth)



# 2.3 Dummy and factor Gender -------------------------------------------------
df <- df %>% mutate(gender_num = ifelse(gender == "Female", 1, NA),
                    gender_num = ifelse(gender == "Male", 0, gender_num),
                    gender_f = factor(gender_num, 
                                      labels = c("Male", "Female")))




# 2.4 Quality summary ---------------------------------------------------------
M1 <-  df %>% 
  select(quality_time_management, 
         quality_service_on_time,
         quality_amount_of_staff, 
         quality_overall)
M2 <- rowMeans(M1, na.rm = TRUE) 
df <- cbind(df, as.numeric(M2))
df <- df %>% 
  rename(quality_summary = `as.numeric(M2)`) # rename variable

# Cronbach's Alpha
alpha_quality_summary <- psych::alpha(M1)
alpha_quality_summary <- alpha_quality_summary$total$raw_alpha
print(myround(alpha_quality_summary, 2))

rm(M1, M2) # delete matrix



# 2.5 Registration Office -----------------------------------------------------
df$visit_service_center[df$visit_service_center == "No"] <- 0
df$visit_service_center[df$visit_service_center == "Yes"] <- 1
df$visit_service_center[df$visit_service_center == "N/A"] <- NA
df$visit_service_center <- as.numeric(df$visit_service_center)






# 2.6 Blame -------------------------------------------------------------------

# Mean index
M1 <- df %>%
  select(blame_responsible,
         blame_fault,
         blame_time_management,
         blame_substantial_service,
         blame_employees,
         blame_overall,
         blame_willing)
M2 <- rowMeans(M1, na.rm = TRUE)
df <- bind_cols(df, tibble::enframe(M2))
df <- df %>% 
  rename(blame_mean = value) # rename variable


# Cronbach's Alpha 
alpha_acbs <- psych::alpha(M1)
alpha_acbs <- alpha_acbs$total$raw_alpha
print(myround(alpha_acbs, 2))

rm(M1, M2) # delete matrix



# 3. Restrict sample  #########################################################
# (no missing values on dv and iv)

# Get participants per group who received treatment
df %>% group_by(group) %>% summarise(mean = n())


# number of missings per variable
sapply(df, function(x) sum(is.na(x)))

df <- df %>%
  filter(!is.na(blame_mean) & #!is.na(income) &
           !is.na(age) & !is.na(gender_num) &
           !is.na(visit_service_center) & 
           !is.na(quality_summary) & 
           !is.na(political_orientation))

df %>% group_by(group) %>% summarise(mean = n())

# 4. Exclude obs failing attention check ######################################
df_full <- df # keep full sample
df_full <- df_full %>% 
  mutate(attention_failed = ifelse(attention_check == "6 - 8 weeks", 0, 1))

df <- filter(df, attention_check == "6 - 8 weeks")

df %>% group_by(group) %>% summarise(mean = n())




# 5. Exploratory Factor Analysis ##############################################

# 5.1 Quality Summary ---------------------------------------------------------
# Parallel Analysis to determine number of factors
para_qs <- fa.parallel(df %>%
                         select(quality_time_management, 
                                quality_service_on_time,
                                quality_amount_of_staff, 
                                quality_overall),
                       fa = "fa", 
                       n.iter = 100)
para_qs

fa.qs <- fa(df %>%
              select(quality_time_management, 
                     quality_service_on_time,
                     quality_amount_of_staff, 
                     quality_overall),
            nfactors = 1, 
            fm = "minres", 
            rotate = "oblimin",
            scores = "regression")
print(fa.qs, digits = 2, sort=TRUE)


fac_qs <- fa.qs$scores # safe factor scores in a new c
fac_qs <- data.frame(fac_qs)
names(fac_qs) <- "fac_qs"

df <- bind_cols(df, fac_qs) 
# add factor scores to data frame and bind by rownames
df <- as_tibble(df)
rm(fac_qs)


# 5.2 Blame -------------------------------------------------------------------
# Parallel Analysis to determine number of factors
para_blame <- fa.parallel(df %>%
                            select(blame_responsible,
                                   blame_fault,
                                   blame_time_management,
                                   blame_substantial_service,
                                   blame_employees,
                                   blame_overall,
                                   blame_willing),
                          fa = "fa", 
                          n.iter = 100)
para_blame

fa_blame <- fa(df %>%
                select(blame_responsible,
                       blame_fault,
                       blame_time_management,
                       blame_substantial_service,
                       blame_employees,
                       blame_overall,
                       blame_willing),
              nfactors =  1, 
              fm = "minres", 
              rotate = "oblimin",
              scores = "regression")

print(fa_blame, digits = 2, sort=TRUE)


fac_blame <- fa_blame$scores # safe factor scores in a new c
fac_blame <- data.frame(fac_blame)
names(fac_blame) <- "fac_blame"

df <- bind_cols(df, fac_blame) 
# add factor scores to data frame and bind by rownames
df <- as_tibble(df)
rm(fac_blame)



# 6. Descriptive Statistics ###################################################

# 6.1 Table 2: Descriptive Statistics =========================================
descr <- df %>%
  select(blame_mean, age, gender, #income, 
         visit_service_center, 
         quality_summary, political_orientation) %>%
  as.data.frame() %>%
  stargazer(type = "html", header = FALSE,
            title = "Table 2: Descriptive Statistics",
            digits = 2, digits.extra = 2,
            covariate.labels = c("Citizens' blame",
                                 "Age",
                                 "Gender (1 = female)",
                                 "Visited municipal service center",
                                 "Perception of quality of service",
                                 "Political orientation (0 = extreme left)"))
descr <- gsub("<td>N</td>", "<td>n</td>", descr)
cat(paste(descr, sep = "\n", collapse = "\n"), "\n", 
    file = here("output", "Table2_descriptives.html"))











# 7. Center age ###############################################################
df$age_center <- df$age - mean(df$age, na.rm = T)



# 8. Regression analysis (Table 3) ############################################

regModel1 <- lm(blame_mean ~ group, 
                data = df)
summary(regModel1)

regModel2 <- lm(blame_mean ~ group + 
                  age_center + gender_f + visit_service_center + 
                  quality_summary + political_orientation, 
                data = df)
summary(regModel2)

star <- stargazer(regModel1, regModel2, 
                  type = "html",
                  report = c("vc*s"),
                  align = TRUE, 
                  title="Effects of treatments on blame attribution", 
                  dep.var.caption = "", 
                  dep.var.labels.include = FALSE,
                  covariate.labels = c("Public opinion in favor of blaming",
                                       "Public opinion not in favor of blaming",
                                       "High intensity of pre-existing blame",
                                       "No pre-existing blame",
                                       "Age (centered)",
                                       "Gender (1 = female)",
                                       "Visited municipal service center",
                                       "Perception of quality of service",
                                       "Political orientation (0 = extreme left)"),
                  omit.stat=c("ser","f"), 
                  star.cutoffs = c(0.05, 0.01, 0.001), 
                  notes = c("<sup>*</sup>p<0.5; <sup>**</sup>p<0.01; <sup>***</sup>p<0.001; standard errors in parentheses"),
                  notes.append = FALSE
                  )
star
star <- c(star[1:36], '<tr><td colspan=\"3\" style=\"border-bottom: 1px solid black\"></td></tr><tr><td colspan=\"3\" style=\"text-align:left\"><em>Note:</em> Ordinary least squares; <sup>*</sup>p&lt;0.05; <sup>**</sup>p&lt;0.01; <sup>***</sup>p&lt;0.001; Standard errors in parentheses.</td></tr>', star[38])
cat(star, sep = "\n", file = here("output", "Table3_Regression.html"))


# 9 ANOVA =====================================================================
anova_blame <- aov(blame_mean ~ group, data = df)
summary(anova_blame)

# 10. Graphs ##################################################################

# 10.1 Pirateplot (Figure 2) ==================================================
df$group_graph <- factor(df$group, 
                            labels = c("Control", 
                                       "Public opinion\n in favor\n of blaming",
                                       "Public opinion\n not in favor\n of blaming",
                                       "High intensity\n of pre-existing\n blame",
                                       "No\n pre-existing\n blame"))
df$group_graph <- relevel(df$group_graph, ref = 1)

pirateplot <- pirateplot(formula = blame_mean ~ group_graph,
                         data = df,
                         theme = 1,
                         inf.method = "ci",
                         pal = gray(.1),
                         main = "Blame attribution in the five experimental groups",
                         ylab = "Blame attribution",
                         xlab = "")

dev.print(pdf, here("output", "Figure2_Pirateplot.pdf"))




# 11. Randomization Check (Table 1) ###########################################

treat_sum <- df %>%
  group_by(group) %>%
  mutate(n = n()) %>%
  select(age, 
         gender_num, 
         visit_service_center, 
         quality_summary, 
         political_orientation, n) %>%
  summarise_all(mean) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%
  filter(rowname != "group") %>% 
  as_tibble() %>%
  mutate(V1 = as.numeric(levels(V1))[V1],
         V2 = as.numeric(levels(V2))[V2],
         V3 = as.numeric(levels(V3))[V3],
         V4 = as.numeric(levels(V4))[V4],
         V5 = as.numeric(levels(V5))[V5]) %>%
  select(rowname, V2, V3, V4, V5, V1) %>%
  rename("Treatment 1" = V2, "Treatment 2" = V3,
         "Treatment 3" = V4, "Treatment 4" = V5,
         "Control" = V1) %>%
  mutate(rowname = ifelse(rowname == "age", "Age", rowname),
         rowname = ifelse(rowname == "gender_num", "Gender", rowname),
         rowname = ifelse(rowname == "visit_service_center", 
                          "Visited municipal service center", rowname),
         rowname = ifelse(rowname == "quality_summary", 
                          "Perception quality of PS", rowname),
         rowname = ifelse(rowname == "political_orientation", 
                          "Political orientation", rowname))
treat_sum


# Checks
chi.gender <- chisq.test(table(df$gender_num, df$group))
chi.reg <- chisq.test(table(df$visit_service_center, df$group))
anova.age <- aov(age ~ group, data = df)
anova.pol <- aov(political_orientation ~ group, data = df)
anova.qs <- aov(quality_summary ~ group, data = df)



# Write results to treat_sum
treat_sum$test <-  NA
treat_sum <- treat_sum %>% 
  mutate(test = ifelse(rowname == "Gender", 
                       chi_output(chi.gender), test),
         test = ifelse(rowname == "Visited municipal service center", 
                       chi_output(chi.reg), test),
         test = ifelse(rowname == "Age", 
                       anova_output(anova.age), test),
         test = ifelse(rowname == "Perception quality of PS", 
                       anova_output(anova.qs), test),
         test = ifelse(rowname == "Political orientation", 
                       anova_output(anova.pol), test)) %>%
  rename(`Randomization Check` = test)
  


write.csv(format(as.data.frame(treat_sum), 
                  digits=2, 
                  na.encode = FALSE),
           here("output", "Table1_Randomization_check.csv"),
           na = "", 
           row.names = FALSE) 



# 12. Correlations (Table 4) ##################################################

cor_df <- df %>% 
  select(blame_mean, 
         age, 
         gender_num, 
         visit_service_center, 
         quality_summary, 
         political_orientation) %>% 
  cor()


# Hide upper triangle
upper <- cor_df # new object
upper[upper.tri(cor_df)] <- NA # remove upper triangle

# remove last columns
limit <- ncol(upper) - 1 # determine last column
upper <- upper[,1:limit] # remove last columns

# Transform to data frame
upper <- data.frame(upper) 

upper <- (setattr(upper, "row.names",  
                  c("(1) Citizens' blame",
                    "(2) Age",
                    "(3) Gender (1 = female)",
                    "(4) Visited municipal service center",
                    "(5) Perception of quality of service",
                    "(6) Political orientation (0 = extreme left)")))

stargazer(upper, 
          type = "html", 
          out = here("output", "Table4_Correlations.html"), 
          summary = FALSE, 
          digits = 2, 
          digits.extra = 2, 
          covariate.labels = c("", "(1)", "(2)", "(3)", "(4)", "(5)")
          )



# 13. Effect Sizes ############################################################
g_t1 <- hedge_function("Public opinion in favor of blaming")
g_t1

g_t2 <- hedge_function("Public opinion not in favor of blaming")
g_t2

g_t3 <- hedge_function("High intensity of pre-existing blame")
g_t3

g_t4 <- hedge_function("No pre-existing blame")
g_t4




# 14. Attention robustness check ##############################################

# 14.1 Descriptive Statistics -------------------------------------------------
descr <- df_full %>%
  select(blame_mean, 
         age, 
         gender_num,
         visit_service_center, 
         quality_summary, 
         political_orientation) %>%
  as.data.frame() %>%
  stargazer(type = "html", 
            header = FALSE,
            title = "Table 2: Descriptive Statistics",
            digits = 2, 
            digits.extra = 2,
            covariate.labels = c("Citizens' blame",
                                 "Age",
                                 "Gender (1 = female)",
                                 "Visited municipal service center",
                                 "Perception of quality of service",
                                 "Political orientation (0 = extreme left)"))
descr <- gsub("<td>N</td>", "<td>n</td>", descr)
cat(paste(descr, sep = "\n", collapse = "\n"), "\n", 
    file = here("output", "Table2_Descriptives_full.html"))






# 14.2 Center age =============================================================
df_full$age_center <- df_full$age - mean(df_full$age, na.rm = T)

# 14.3 Regression Analysis ====================================================

# groups only
regModel3 <- lm(blame_mean ~ group, 
                data = df_full)
summary(regModel3)

# groups + attention dummy
regModel4 <- lm(blame_mean ~ group + attention_failed, 
                data = df_full)
summary(regModel4)


# groups # control
regModel5 <- lm(blame_mean ~ group + 
                  age_center + gender_num + visit_service_center + 
                  quality_summary + political_orientation, 
                data = df_full)
summary(regModel5)


regModel6 <- lm(blame_mean ~ group + attention_failed +
                  age_center + gender_num + visit_service_center + 
                  quality_summary + political_orientation, 
                data = df_full)
summary(regModel6)

star_full <- stargazer(regModel1, regModel3, regModel4, 
                       regModel2, regModel5, regModel6,
                       type = "html",
                       report = c("vc*s"),
                       align = TRUE,
                       title="Effects of treatments on blame attribution",
                       dep.var.caption = "", 
                       dep.var.labels.include = FALSE,
                       covariate.labels = c("Public opinion in favor of blaming",
                                            "Public opinion not in favor of blaming",
                                            "High intensity of pre-existing blame",
                                            "No pre-existing blame",
                                            "Attention dummy",
                                            "Age (centered)",
                                            "Gender (1 = female)",
                                            "Visited municipal service center",
                                            "Perception of quality of service",
                                            "Political orientation (0 = extreme left)"),
                       omit.stat=c("ser","f"), 
                       star.cutoffs = c(0.05, 0.01, 0.001), 
                       notes = c("<sup>*</sup>p<0.5; <sup>**</sup>p<0.01; <sup>***</sup>p<0.001; standard errors in parentheses"), 
                       notes.append = FALSE
)
star_full
star_full <- c(star_full[1:39], '<tr><td colspan=\"7\" style=\"border-bottom: 1px solid black\"></td></tr><tr><td colspan=\"7\" style=\"text-align:left\"><em>Note:</em> Ordinary least squares; <sup>*</sup>p&lt;0.05; <sup>**</sup>p&lt;0.01; <sup>***</sup>p&lt;0.001; Standard errors in parentheses.</td></tr>', star[41])
cat(star_full, sep = "\n", file = here("output", "Table3_Regression_full.html"))



# 14.4 Randomization Check ----------------------------------------------------

# 14.1.1 per treatment group 
treat_sum_full <- df_full %>%
  group_by(group) %>%
  mutate(n = n()) %>%
  select(age, 
         gender_num, 
         visit_service_center, 
         quality_summary, 
         political_orientation, 
         n) %>%
  summarise_all(mean) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%
  filter(rowname != "group") %>% 
  as_tibble() %>%
  mutate(V1 = as.numeric(levels(V1))[V1],
         V2 = as.numeric(levels(V2))[V2],
         V3 = as.numeric(levels(V3))[V3],
         V4 = as.numeric(levels(V4))[V4],
         V5 = as.numeric(levels(V5))[V5]) %>%
  select(rowname, V2, V3, V4, V5, V1) %>%
  rename("Treatment 1" = V2, "Treatment 2" = V3,
         "Treatment 3" = V4, "Treatment 4" = V5,
         "Control" = V1) %>%
  mutate(rowname = ifelse(rowname == "age", "Age", rowname),
         rowname = ifelse(rowname == "gender_num", "Gender", rowname),
         rowname = ifelse(rowname == "visit_service_center", 
                          "Visited municipal service center", rowname),
         rowname = ifelse(rowname == "quality_summary", 
                          "Perception quality of PS", rowname),
         rowname = ifelse(rowname == "political_orientation", 
                          "Political orientation", rowname))
treat_sum_full

# Checks
chi.gender_full <- chisq.test(table(df_full$gender_num, df_full$group))
chi.reg_full <- chisq.test(table(df_full$visit_service_center, df_full$group))
anova.age_full <- aov(age ~ group, data = df_full)
anova.pol_full <- aov(political_orientation ~ group, data = df_full)
anova.qs_full <- aov(quality_summary ~ group, data = df_full)




# Write results to treat_sum
treat_sum_full$test <-  NA
treat_sum_full <- treat_sum_full %>% 
  mutate(test = ifelse(rowname == "Gender", 
                       chi_output(chi.gender_full), test),
         test = ifelse(rowname == "Visited municipal service center", 
                       chi_output(chi.reg_full), test),
         test = ifelse(rowname == "Age", 
                       anova_output(anova.age_full), test),
         test = ifelse(rowname == "Perception quality of PS", 
                       anova_output(anova.qs_full), test),
         test = ifelse(rowname == "Political orientation", 
                       anova_output(anova.pol_full), test)) %>%
  rename(`Randomization Check` = test)



write.csv(format(as.data.frame(treat_sum_full), 
                 digits=2, 
                 na.encode = FALSE), 
          here("output", "randomization_check_full.csv"),
          na = "",
          row.names = FALSE) 
