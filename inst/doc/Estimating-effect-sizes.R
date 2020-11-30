## ---- message = FALSE---------------------------------------------------------
library(scdhlm)

## ----Lambert------------------------------------------------------------------
data(Lambert)

Lambert_ES <- effect_size_ABk(outcome = outcome, treatment = treatment, 
                              id = case, phase = phase, time = time, 
                              data = Lambert)

print(Lambert_ES, 5)

## ----Lambert-summary----------------------------------------------------------
summary(Lambert_ES)

## ----Lambert_sensitivity------------------------------------------------------
Lambert_ES2 <- effect_size_ABk(outcome = outcome, treatment = treatment, 
                id = case, phase = phase, time = time, 
                data = Lambert, phi = 0.10)
print(Lambert_ES2, digits = 5)

Lambert_ES3 <- effect_size_ABk(outcome = outcome, treatment = treatment, 
                id = case, phase = phase, time = time, 
                data = Lambert, phi = 0.35)
print(Lambert_ES3, digits = 5)

## ----Anglesea-----------------------------------------------------------------
data(Anglesea)
Anglesea_ES <- effect_size_ABk(outcome, condition, case, phase, session, data = Anglesea)
Anglesea_ES

## ----Saddler------------------------------------------------------------------
data(Saddler)

quality_ES <- effect_size_MB(outcome, treatment, case, time, 
                             data = subset(Saddler, measure=="writing quality"))
complexity_ES <- effect_size_MB(outcome, treatment, case, time , 
                                data = subset(Saddler, measure=="T-unit length"))
construction_ES <- effect_size_MB(outcome, treatment, case, time, 
                                  data = subset(Saddler, measure=="number of constructions"))

data.frame(
  quality = unlist(quality_ES),
  complexity = unlist(complexity_ES), 
  construction = unlist(construction_ES) 
)[c("delta_hat","V_delta_hat","nu","phi","rho"),]


## ----quality-RML--------------------------------------------------------------
quality_RML <- lme(fixed = outcome ~ treatment, 
                   random = ~ 1 | case, 
                   correlation = corAR1(0, ~ time | case), 
                   data = subset(Saddler, measure=="writing quality"))
summary(quality_RML)

## ----quality-g-mlm, eval=TRUE-------------------------------------------------
quality_ES_RML <- g_mlm(quality_RML, p_const = c(0,1), 
                         r_const = c(1,0,1), returnModel = TRUE)
summary(quality_ES_RML)

## ----Laski-RML, eval = TRUE---------------------------------------------------
data(Laski)

# Hedges, Pustejovsky, & Shadish (2013)
Laski_ES_HPS <- effect_size_MB(outcome, treatment, case, time, data= Laski)

# Pustejovsky, Hedges, & Shadish (2014)
Laski_RML <- lme(fixed = outcome ~ treatment,
                 random = ~ 1 | case, 
                 correlation = corAR1(0, ~ time | case), 
                 data = Laski)
summary(Laski_RML)

## ----Laski-ES-RML, eval = TRUE------------------------------------------------
Laski_ES_RML <- g_mlm(Laski_RML, p_const = c(0,1), r_const = c(1,0,1))

# compare the estimates
data.frame(
  HPS = with(Laski_ES_HPS, c(SMD = delta_hat, SE = sqrt(V_delta_hat), 
                             phi = phi, rho = rho, nu = nu)),
  RML = with(Laski_ES_RML, c(g_AB, SE_g_AB, theta$cor_params, 
                             theta$Tau$case / (theta$Tau$case + theta$sigma_sq), nu))
)

## ----Laski-RML2---------------------------------------------------------------
Laski_RML2 <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case, 
                 correlation = corAR1(0, ~ time | case), 
                 data = Laski)
summary(Laski_RML2)
anova(Laski_RML, Laski_RML2)

## ----Laski-ES-RML2, eval=TRUE-------------------------------------------------
Laski_ES_RML2 <- g_mlm(Laski_RML2, p_const = c(0,1), r_const = c(1,0,0,0,1))
Laski_ES_RML2

## ----Schutte-graph, warning = FALSE, message = FALSE, fig.width = 7, fig.height = 7----
library(ggplot2)

data(Schutte)

graph_SCD(case = case, phase = treatment,
          session = week, outcome = fatigue,
          design = "MB", data = Schutte) + 
  facet_wrap(~ case, ncol = 3) + 
  theme(legend.position = "bottom")


## ----Schutte-center-----------------------------------------------------------
# time-point constants
A <- 2
B <- 9
# center at follow-up time

# create time-by-trt interaction
Schutte <- 
  preprocess_SCD(case = case, phase = treatment,
                 session = week, outcome = fatigue,
                 design = "MB", center = B, data = Schutte)


## ----Schutte-Model3-----------------------------------------------------------
hlm1 <- lme(fixed = fatigue ~ week + treatment + week_trt, 
            random = ~ 1 | case, 
            correlation = corAR1(0, ~ week | case),
            data = Schutte,
            method = "REML")
summary(hlm1)

## ----Schutte-ES-Model3, eval = TRUE-------------------------------------------
Schutte_g1 <- g_mlm(hlm1, p_const = c(0,0,1, B - A), r_const = c(1,0,1))
Schutte_g1

## ----Schutte-Model4-----------------------------------------------------------
hlm2 <- update(hlm1, random = ~ week | case, 
               control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))
summary(hlm2)
anova(hlm1, hlm2)

## ----Schutte-ES-Model4, eval = TRUE-------------------------------------------
Schutte_g2 <- g_mlm(hlm2, p_const = c(0,0,1, B - A), r_const = c(1,0,0,0,1))
Schutte_g2

## -----------------------------------------------------------------------------
unlist(extract_varcomp(hlm2))

## ----Schutte-Model5, warning = FALSE------------------------------------------
hlm3 <- update(hlm2, random = ~ week + week_trt | case, 
               control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))
summary(hlm3)
anova(hlm2, hlm3)

## ----Schutte-ES-Model5, eval = TRUE-------------------------------------------
Schutte_g3 <- g_mlm(hlm3, p_const = c(0,0,1,B - A), r_const = c(1,0,0,0,0,0,0,1))
Schutte_g3

# compare effect size estimates
data.frame(
  MB3 = round(unlist(Schutte_g1[c("g_AB","SE_g_AB","nu")]), 3), 
  MB4 = round(unlist(Schutte_g2[c("g_AB","SE_g_AB","nu")]), 3),
  MB5 = round(unlist(Schutte_g3[c("g_AB","SE_g_AB","nu")]), 3)
)

