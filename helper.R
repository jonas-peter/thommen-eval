library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(emmeans)

# https://www.datanovia.com/en/lessons/ancova-in-r/

################################################################################
#Functions
################################################################################
plot2setts <- function(x1,y1,x2,y2, xlable, ylable, col1, col2, lg1, lg2, titel)
{
  xlim <- range(c(x1,x2))
  ylim <- range(c(0,y1,y2))

  par(mar=c(5,4,2,7))
  # Plot x1 - y1 with regression line

  plot(x1, y1, col=col1, xlim=xlim, ylim=ylim,
       main = titel,
       pch = 4,
       xlab = xlable,
       ylab = ylable)


  model1 <- lm(y1 ~ x1)
  abline(model1, col=col1)
  print(paste('Lin Reg. of Model1: ', lg1))
  print(summary(model1))

  points(x2, y2, col=col2,
         pch=4)
  model2 <- lm(y2 ~ x2)
  abline(model2, col=col2)
  print(paste('Lin Reg. of Model2: ', lg2))
  print(summary(model2))

  legend("left", inset=c(1,0),
         legend=c(paste(lg1, '\nR^2:', round(summary(model1)$adj.r.square, 3)),
                  # paste('p-value:', round(summary(model1)$coefficients[2,4],3)),
                  '',
                  paste(lg2, '\nR^2:', round(summary(model2)$adj.r.square, 3))),
         # paste('p-value:', round(summary(model2)$coefficients[2,4],3)),

         col=c(col1,
               col1,
               col2,
               col2),
         lty=c(0,-1,0,-1,0),
         pch=c(4,-1,4,-1,4), cex=0.8,
         box.lty=0,
         xpd=TRUE, bty="n")
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  return(model1)
}

ancova_analysis <- function(daten, gruppe, dependent, independent)
{
  daten$dependent_fit<-(daten[[dependent]])
  daten$independent_fit<-(daten[[independent]])
  daten$gruppe_fit <- (daten[[gruppe]])

  # Check for Linearity assumption:
  plot1 <- ggscatter(
    daten, x = independent, y = dependent,
    color = gruppe, add = "reg.line"
  ) +
    stat_regline_equation(
      aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = gruppe_fit)
    )
  ## Check for Homogeneity of regression slopes:
  table1 <- daten %>% anova_test(dependent_fit ~ gruppe_fit*independent_fit)

  ##  Check for Normality of residuals
  # Fit the model, the covariate goes first
  model <- lm(dependent_fit ~ independent_fit + gruppe_fit, data = daten)
  # Inspect the model diagnostic metrics
  model.metrics <- augment(model) %>% select(-.hat, -.sigma, -.fitted) # Remove details
  table2 <- head(model.metrics, 3)
  # Assess normality of residuals using shapiro wilk test
  table3 <- shapiro_test(model.metrics$.resid)

  # ## Check for Homogenity of varainces:
  table4 <- model.metrics %>% levene_test(.resid ~ gruppe_fit)
  # table4 <- 1

  ## Check for Outliers
  table5 <- model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame()

  ## Computation of ANCOVA
  res.aov <- daten %>% anova_test(dependent_fit ~ independent_fit + gruppe_fit)
  table6 <- get_anova_table(res.aov)
  # table6 <- Anova(res.aov, type="III")

  # Pairwise comparisons
  library(emmeans)
  pwc <- daten %>%
    emmeans_test(
      dependent_fit ~ gruppe_fit,
      covariate = independent_fit,
      p.adjust.method = "bonferroni"
    )
  pwc1 <- pwc

  # Display the adjusted means of each group
  # Also called as the estimated marginal means (emmeans)
  table7 <- get_emmeans(pwc)

  # Visualization: line plots with p-values
  pwc <- pwc %>% add_xy_position(x = "gruppe_fit", fun = "mean_se")
  plot2 <-
    ggline(get_emmeans(pwc), x = "gruppe_fit", y = "emmean") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
    stat_pvalue_manual(pwc, hide.ns = TRUE, tip.length = FALSE) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE),
      caption = get_pwc_label(pwc)
    )
  return(list(plot1, table1, table2, table3, table4, table5, table6, pwc1, table7, plot2))
}
