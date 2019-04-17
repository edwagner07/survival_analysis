
# Part 1: Building the data set


```{r, cache = T, results = "hide", echo = F}
library(dplyr)
library(cowplot)
library(dat.table)

## Log-odds and response probability functions
loglink <- function(z) 1 / (1 + exp(-z))
invlog  <- function(z) log(z / (1 - z))
purchase_prob <- function(b_0, b_1, x_1, b_2, x_2, b_3 = 1, x_3) {
  loglink(b_0 + 
            (b_1 * x_1) + 
            (b_2 * x_2) + 
            (b_3 * x_3))
}

## Functions for generating random age, income, and time function
agefunc <- function() {
  rnorm(n = 1, mean = 45, sd = 20)}

incfunc <- function(age) {
  30000 + 1000 * ((age - 20)/5 + rnorm(1, mean = 0, sd = 10))
}

gamma_log  <- invlog(dgamma(1:20, 2, 0.35))
gamma_func <- function(x) invlog(dgamma(x + 1, 2, 0.35))


## Set beta values
b_0 <- -6
b_1 <- 0.02
b_2 <- 2.5e-05

```


```{r, cache = T, results = "hide", echo = F}
set.seed(123)
subcount = 1000000
subject_list = as.list(seq_len(subcount))


## Reset counting variables; build per-subject data frames
for (i in 1:subcount) {
  j      <- 1
  k      <- sample(1:20, 1)
  subage <- agefunc()
  subinc <- incfunc(subage)

  subject_weeks <- data.frame(
    income   = rep(subinc, 20),
    age      = rep(subage, 20),
    id       = rep(i, 20),
    week     = 0:19,
    purchase = NA)

  subject_weeks$z <- purchase_prob(
    b_0,
    b_1,
    subage,
    b_2,
    subinc,
    x_3 = gamma_log)

### Loop to calculate whether purchase occurs for each week
  while (j < k + 1 && sum(subject_weeks$purchase, na.rm = T) == 0) {
    subject_weeks$purchase[j]  <- rbinom(1, 1, subject_weeks$z[j])
    j  <- j + 1
  }

subject_list[[i]] <- subject_weeks

}

## bind all individual frames and remove NA values
week_frame <- data.table::rbindlist(subject_list)
week_frame <- week_frame[complete.cases(week_frame), ]


```


## Calculating weekly hazard rates
```{r, cache = T, results = "hide", echo = F}
sub_counts = data.frame(week_frame %>%
  group_by(week) %>%
  summarise("total_subs" = n())
)

weekdata <- data.frame(week_frame$week[week_frame$purchase == 1])
names(weekdata) <- "weeknum"
week_summary <- data.frame(weekdata %>%
  group_by(weeknum) %>%
  summarise(count = n() / nrow(weekdata)))
week_summary$gamma <- dgamma(week_summary$weeknum + 1, 2, 0.35)

hazards <- week_frame %>%
  group_by(week) %>%
  summarise(hazard = sum(purchase) / n())
hazards$d_hazard <- hazards$hazard / sum(hazards$hazard)
hazards$gamma <- dgamma(hazards$week + 1, 2, 0.35)
hazards = hazards[1:18, ]


```


# Part 2: Demonstrating the hazard function
## Plotting response counts and hazard rates
```{r, echo = F}
## Histogram of all response weeks versus original distribution function
plot1 <- ggplot(data = data.frame(week_summary), aes(x = weeknum, y = count)) +
  geom_bar(stat = "identity", col = "black", fill = "grey50") +
  geom_line(aes(x = weeknum, y = gamma), size = 1) +
  ylab("Response probability density") +
  xlab("week") +
  ylim(0, 0.17) +
  coord_fixed(ratio = 120) +
  ggtitle("Simple response counts")

plot2 <- ggplot(data = data.frame(hazards), aes(x = week, y = d_hazard)) +
  geom_bar(stat = "identity", col = "black", fill = "grey50") +
  geom_line(aes(x = week, y = gamma), size = 1) +
  ylab("") + 
  ylim(0, 0.17) +
  coord_fixed(ratio = 120) +
  ggtitle("Hazard Rates")

fig2 <- cowplot::plot_grid(plot1, plot2, labels = c("A", "B"))
title2 <- ggdraw() + 
  draw_label("Predicting the response distribution", fontface = "bold")

cowplot::plot_grid(title2, fig2, ncol=1, rel_heights=c(0.1, 1))

```


# Part 4: Demonstrating the stratified sample
## Defining the functions to build SRS and stratified sample models
```{r, cache = T, results = "hide", echo = F, message = F, warning = F}

##Function for 95% confidence interval half-width
ci95 = function(x) 1.96 * (sd(x) / sqrt(length(x)))

## Functions and purpose
### strat_pred_func calculates RMSE for predictions built on stratified samples
### srs_pred_func calculates RMSE for predictions built on SRS
### strat_c_func creates table of coefficients for stratified models
### srs_c_func creates table of coefficients for SRS models

## NOTE: for this report, aggregate tables were created in advance and pre-
## loaded to save processing time.

## Stratified Prediction Function ----------------------------------------------
strat_pred_func = function(k = 1, C, nboot = 1) {
  
  rmse = function(predicted, actual) {
    sqrt(mean(actual - predicted)^2)
  }

  all_cval_list = as.list(seq_len(length(C)))
  all_strat_rmse_list = as.list(seq_len(length(C)))
  for (Cval in 1:length(C)) {
    kval = 1
    print(paste("round", Cval, "of", length(C), ": initializing"))
    weekly_offset <- 
      data.frame(
        week_frame_1 %>%
          group_by(week) %>%
          summarise(
            N_0 = n() - sum(purchase),
            N_1 = sum(purchase),
            n_0 = round(k[kval] * (n() / C[Cval] - sum(purchase))),
            n_1 = round(sum(purchase) * k[kval])
          )
    )
    weekly_offset$offset = NA
    weekly_offset$offset = with(weekly_offset, log((n_1 * N_0) / (n_0 * N_1)))
    weekly_offset$offset[is.nan(weekly_offset$offset)] = log(C[Cval])
    weekly_offset$offset[weekly_offset$offset == -Inf] = log(C[Cval])
    
    week_frame = week_frame_1[, 1:7]
    week_frame = merge(week_frame, weekly_offset[, c(1, 6)], by = "week")
    week_0 = week_frame[week_frame$purchase == 0, ]
    week_1 = week_frame[week_frame$purchase == 1, ]
    sample_counts <-
      data.frame(week_frame %>% #total non-responses needed per week
                   group_by(week) %>%
                   summarize(n_1 = round(sum(purchase) * k[kval]),
                             n_0 = round(k[kval] * (
                               n() / C[Cval] - sum(purchase)
                             ))))
    
    
    coef_list = as.list(seq_len(nboot))
    prediction_list = as.list(seq_len(nboot))
    rmse_list = as.list(seq_len(nboot))
    week_sample = sample_n(week_frame, 1000)
    
    #Full loop for parametric bootstrap
    print(paste("round", Cval, "of", length(C), ": resampling"))    
    for (j in 1:nboot) {
      
      if (j %% 5 == 0) {
        print(paste(
          "round",
          Cval,
          "of",
          length(C),
          "resampling",
          round(j / nboot * 100),
          "% complete"
        ))
      }
      control_list <- as.list(seq_len(20))
      case_list    <- as.list(seq_len(20))
      
      ## Creating collection of non-responses
      for (i in 1:20) {
        control_list[[i]] <- sample_n(week_0[week_0$week == (i - 1), ],
                                      size = weekly_offset$n_0[i])
      }
      
      if (k[kval] != 1) {
        ## Creating collection of responses
        for (i in 1:20) {
          case_list[[i]] <- sample_n(week_1[week_1$week == (i - 1), ],
                                     size = weekly_offset$n_1[i])
        }
        case_control <- rbind(data.table::rbindlist(case_list),
                              data.table::rbindlist(control_list))
      } else {
        case_control <- rbind(week_frame[week_frame$purchase == 1, ],
                              data.table::rbindlist(control_list))
      }
      
      stratreg <-  glm(
        purchase ~ age +
          income +
          gamma_func(week) +
          offset(offset),
        data = case_control,
        family = binomial
      )
      
      stratpred = 
        invlog(predict(stratreg, week_sample, type = "response")) - 
        week_sample$offset
      
      strat_pred_diff = data.frame(invlog(week_sample$z2) - stratpred)
      names(strat_pred_diff) = "difference"
      strat_pred_diff$condition = 
        as.factor(rep("Stratified", nrow(strat_pred_diff)))
      
      difference_rmse = data.frame(rmse(invlog(week_sample$z2), stratpred))
      names(difference_rmse) = "rmse"
      difference_rmse$cval = C[Cval]
      difference_rmse$condition = "Stratified"
      
      rmse_list[[j]] = difference_rmse
      strat_pred_diff$cval = rep(C[Cval], nrow(strat_pred_diff))
      prediction_list[[j]] <- data.frame(strat_pred_diff)
      
    }
    
    all_cval_list[[Cval]] = data.table::rbindlist(prediction_list)
    all_strat_rmse_list[[Cval]] = data.table::rbindlist(rmse_list)

  }
  print("Finishing up...")
  all_cval_table = data.table::rbindlist(all_cval_list)
  all_cval_table$cval = factor(all_cval_table$cval)
  all_cval_table$cval = 
    factor(all_cval_table$cval, levels = rev(levels(all_cval_table$cval)))
  
  all_strat_rmse_table = data.table::rbindlist(all_strat_rmse_list)
  all_strat_rmse_table$cval = factor(all_strat_rmse_table$cval)
  all_strat_rmse_table$cval = 
    factor(all_strat_rmse_table$cval, 
           levels = rev(levels(all_strat_rmse_table$cval)))
  
  saveRDS(all_strat_rmse_table, "all_strat_rmse_table_100.rds")
  saveRDS(all_cval_table, "strat_prediction_table.rds")
}



## SRS Prediction Function ----------------------------------------------------

srs_pred_func = function(k = 1, C, nboot = 1) {
  
  rmse = function(predicted, actual) {
    sqrt(mean(actual - predicted)^2)
  }
  
  srs_all_cval = as.list(seq_len(length(C)))
  all_srs_rmse_list = as.list(seq_len(length(C)))
  for (Cval in 1:length(C)) {
    kval = 1
    print(paste("round", Cval, "of", length(C), ": initializing"))
    week_frame = week_frame_1[, 1:7]
    n_0 = (nrow(week_frame) - sum(week_frame$purchase)) / C[Cval] - 
      sum(week_frame$purchase)
    N_0 = nrow(week_frame) - sum(week_frame$purchase)
    week_frame$offset = rep(log(N_0 / n_0), nrow(week_frame))
    
    srs_prediction_list = as.list(seq_len(nboot))
    srs_rmse_list         = as.list(seq_len(nboot))
    week_sample = sample_n(week_frame, 1000)
    #Full loop for parametric bootstrap
    print(paste("round", Cval, "of", length(C), ": resampling"))    
    for (j in 1:nboot) {
      
      if (j %% 5 == 0) {
        print(paste(
          "round",
          Cval,
          "of",
          length(C),
          "resampling",
          round(j / nboot * 100),
          "% complete"
        ))
      }
      
      ## Creating case-control set
      
      case_control <- 
        rbind(week_frame[week_frame$purchase == 1, ],
              sample_n(week_frame[week_frame$purchase == 0, ], n_0))
      
      srsreg <-  glm(
        purchase ~ age +
          income +
          gamma_func(week) +
          offset(offset),
        data = case_control,
        family = binomial
      )
      
      srspred = 
        invlog(predict(srsreg, week_sample, type = "response")) - 
        week_sample$offset
      
      srs_pred_diff = data.frame(invlog(week_sample$z2) - srspred)
      names(srs_pred_diff) = "difference"
      srs_pred_diff$condition = 
        as.factor(rep("SRS", nrow(srs_pred_diff)))
      
      difference_rmse = data.frame(rmse(invlog(week_sample$z2), srspred))
      names(difference_rmse) = "rmse"
      difference_rmse$cval = C[Cval]
      difference_rmse$condition = "SRS"
      
      srs_rmse_list[[j]] = difference_rmse
      
      srs_pred_diff$cval = rep(C[Cval], nrow(srs_pred_diff))
      srs_prediction_list[[j]] <- data.frame(srs_pred_diff)
      
      srs_prediction_list[[j]] <- data.frame(srs_pred_diff)
      
    }
    srs_all_cval[[Cval]] = data.table::rbindlist(srs_prediction_list)
    all_srs_rmse_list[[Cval]] = data.table::rbindlist(srs_rmse_list)

  }
  srs_difference_table = data.table::rbindlist(srs_all_cval)
  srs_difference_table$cval = factor(srs_difference_table$cval)
  srs_difference_table$cval = 
    factor(srs_difference_table$cval, 
           levels = rev(levels(srs_difference_table$cval)))
  
  all_srs_rmse_table = data.table::rbindlist(all_srs_rmse_list)
  all_srs_rmse_table$cval = factor(all_srs_rmse_table$cval)
  all_srs_rmse_table$cval = 
    factor(all_srs_rmse_table$cval, 
           levels = rev(levels(all_srs_rmse_table$cval)))
  
  saveRDS(all_srs_rmse_table, "all_srs_rmse_table_100.rds")  
  
  
  saveRDS(srs_difference_table, "srs_prediction_table.rds")
}

```


## Building the plots
```{r, cache = T, results = "hide", echo = F, message = F, warning = F}
## Processing -----------------------------------------------------------------

strat_rmse_100 = readRDS("all_strat_rmse_table_100.rds")
srs_rmse_100 = readRDS("all_srs_rmse_table_100.rds")
all_rmse = rbind(strat_rmse_100, srs_rmse_100)
rmse_summary = all_rmse%>%
  group_by(cval, condition) %>%
  summarise(mean_rmse = mean(rmse), ci95 = ci95(rmse))

fullcoef    <- readRDS("fullcoef.rds")
strat_table <- readRDS("strat_table.rds")
srs_table   <- readRDS("srs_table.rds")
names(strat_table)[4] = "timevar"
names(srs_table)[4] = "timevar"

## Plotting -------------------------------------------------------------------
pd = position_dodge(0.2)
ggplot(data = rmse_summary, aes(x = cval, group = condition)) +
  geom_line(aes(x = cval, y = mean_rmse, colour = condition, group = condition),
            position = pd, size = 1) +
  geom_errorbar(aes(ymin = mean_rmse - ci95, ymax = mean_rmse + ci95,
                    colour = condition),
                position = pd, width = 0.4, size = 1) +
  geom_point(aes(y = mean_rmse, colour = condition), size = 2, position = pd) +
  scale_colour_manual(values = c("grey50", "black")) +
  ggtitle("Figure 4: Error of stratified and SRS model predictions") +
  xlab("Compression factor") +
  ylab("Root Mean Squared Error") +
  theme(axis.text.x = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10))

```
# Code snippets for examples

logreg_object <- glm(response ~ age + income + gamma_func(week), 
                     data = case_control, family = "binomial")

logreg_object <- glm(response ~ age + income + factor(week), 
                     data = case_control, family = "binomial")
