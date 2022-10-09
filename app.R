#
# Strength of Evidence
# Reproducible R Script for Binder
# Paper: Why Most Results of Socio-Technical Security User Studies Are False
#
# Note: The paper was directly compiled from the R sourcecode and raw data using knitR.
#       This environment uses the same functions, but constitutes a separate code-base w/o knitR integration.
#

#
# Copyright (c) 2022 Anonymized Author of this work.
#
# This software is made as available  as is for the purpose of anonymous peer-review.
# The software can be be executed and evaluated by peer-reviewers.
# Once the paper and software is published, the overall work will be made available under
# an appropriate open source license.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#


#
# Setup of Environment
#

# Loading Libraries

library(plyr)
library(dplyr)
library(pwr)
library(ggplot2)
library(ggExtra)
library(ggridges)
library(viridis)
library(scico)
library(weights)
library(car)
library(shiny)


# Reading Sample Sizes of Coded Papers

dfSampleSizes <- read.csv("./Classification_Assignments-Paper_anonymized.csv",
                          colClasses=c("character", "character", "character", "character",
                                       "character", "character", "character", "character", 
                                       "character", "character", "character", "character",
                                       "character", "character", "character", 
                                       "numeric", "numeric",
                                       "numeric", "numeric",
                                       NULL, NULL, NULL, NULL))
dfSampleSizes$Sample <- as.numeric(dfSampleSizes$Sample)


# Reading Effect Size Table

dfStatcheck_SLR_ES_N <- read.csv("./SLR_ES_Master_Table_anonymized.csv", as.is = T)
dfStatcheck_SLR_ES_N <- dfStatcheck_SLR_ES_N[which(is.na(dfStatcheck_SLR_ES_N$Excluded) | dfStatcheck_SLR_ES_N$Excluded == FALSE),]


#
# ES/Power Simulation Utility Functions
#

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


#
# Computes the PPV including considering biases
# that is computes the posteriori probability, given prior, alpha and power.
# 
# Inputs:
#   prior: Prior probability
#   aplha: Significance level
#   power=1-beta: Power of the study
#
computePPVBias <- function(prior, alpha = 0.05, power = 0.8, bias = 0.2) {
  #beta <- 1-power
  #alpha_power_ratio <- alpha/power
  if (missing(prior) || is.null(prior)) {
    stop('The prior must be specified, received NULL.')
  }
  if (missing(power) || is.null(power)) {
    stop('The power must be specified, received NULL.')
  }
  
  # Formula by Ioannidis
  ppv = tryCatch({
    (power*prior + bias*(1-power)*prior) / (prior + alpha - (1-power)*prior + bias - bias*alpha + bias*(1-power)*prior)
  }, warning = function(w) {
    stop(w)
  }, error = function(e) {
    stop(e)
  }, finally = {
    
  })
  
  names(ppv) <- "ppvb"
  ppv
}

#
# Computes the PPV without considering biases
# that is computes the posteriori probability, given prior, alpha and power.
# 
# Inputs:
#   prior: Prior probability
#   aplha: Significance level
#   power=1-beta: Power of the study
#
computePPV <- function(prior, alpha = 0.05, power = 0.8) {
  #beta <- 1-power
  #alpha_power_ratio <- alpha/power
  if (missing(prior) || is.null(prior)) {
    stop('The prior must be specified, received NULL.')
  }
  if (missing(power) || is.null(power)) {
    stop('The power must be specified, received NULL.')
  }
  
  ### ppv is P(T|R) = P(R|T)*P(T) / (P(R|T)P(T)   +   P(R|not T))P(not T)
  ###               = power*prior / (power*prior + alpha*(1-prior))
  ppv = tryCatch({
    power*prior / (power*prior + alpha*(1-prior))
  }, warning = function(w) {
    stop(w)
  }, error = function(e) {
    stop(e)
  }, finally = {
    
  })
  
  names(ppv) <- "ppv"
  ppv
}


computeSEfromCI <- function(UL, LL) {
  return((UL-LL)/3.92)
}

convertCItoCI.MCC <- function(es, UL, LL, NMC = 1, sig.level = 0.05) {
  se <- computeSEfromCI(UL, LL)
  conf.level <- 1-(sig.level/NMC)
  z <- qnorm(1-(sig.level/(2*NMC)))
  UL.MCC = es + z*se
  LL.MCC = es - z*se
  return(list(es = es, UL = UL, LL = LL, NMC, orig.sig.level = sig.level,
              UL.MCC = UL.MCC, LL.MCC = LL.MCC, sig.level = sig.level/NMC,
              z = z,
              conf.level.MCC = conf.level))
}

convertCItoCI.CL <- function(es, UL, LL, conf.level = 0.95) {
  se <- computeSEfromCI(UL, LL)
  sig.level <- 1-conf.level
  z <- qnorm(1-(sig.level/2))
  UL.MCC = es + z*se
  LL.MCC = es - z*se
  return(list(es = es, UL = UL, LL = LL, orig.sig.level = sig.level,
              UL.MCC = UL.MCC, LL.MCC = LL.MCC, sig.level = sig.level,
              conf.level.MCC = conf.level))
}

# Based on the article by Altman and Bland (2011) - How to obtain the P value from a confidence interval (BMJ 2011;343:d2304)
# Tested against his example
computePfromCI <- function(est, LL, UL, confidence.level=.95) {
  z.level <- qnorm((1-confidence.level)/2,lower.tail=FALSE)
  SE <-  (UL - LL)/(2*z.level)
  z <- abs(est/SE)
  p.value <- 2*pnorm(q=z, lower.tail=FALSE)
  p.value.Altman <- exp(-0.717*z - 0.416*z^2)
  return(list("conf.level" = confidence.level, 
              "z.level" = z.level, 
              "SE" = SE,
              "z" = z,
              "p.value" = p.value,
              "p.value.A" = p.value.Altman))
}




#
# Analysis of Actual Power wrt. to Effect Sizes
#

#
# ES Extraction utility functions. 
# Meant to extract standard effect sizes from different scenarions
# Typically these effect sizes are then converted to log odds and their confidence interval as standard measure.
# The output dfESCombined contains all effect sizes coded.
# This is the list of all ES with their ***observed*** power (incl. under MCC)
#
  
  esc_from_r <- function(r, n, study) {
    require(compute.es)
    es.res <- compute.es::res(r = r, n = n, dig = 10)
    return(as.data.frame(list(study = study, 
                              es = es.res$lOR, 
                              weight = NA, 
                              sample.size = es.res$N.total,
                              se = NA, 
                              var = (pi^2*es.res$var.d)/(3), 
                              ci.lo = es.res$l.lor, 
                              ci.hi = es.res$u.lor,
                              measure = "logit")))
  }

esc_from_msd <- function(m1i, m2i, sd1i, sd2i, n1i, n2i, study) {
  require(compute.es)
  
  es.mes <- compute.es::mes(m.1 = m1i, m.2 = m2i, sd.1 = sd1i, sd.2 = sd2i, n.1 = n1i, n.2 = n2i, dig = 10)
  return(as.data.frame(list(study = study, 
                            es = es.mes$lOR, 
                            weight = NA, 
                            sample.size = es.mes$N.total,
                            se = NA, 
                            var = (pi^2*es.mes$var.d)/(3),
                            ci.lo = es.mes$l.lor, 
                            ci.hi = es.mes$u.lor,
                            measure = "logit")))
  # Source for variance compuration for log odds ratio: compute.es documentation
  # v_(log(o))= (pi^2v_(d))/ (3)
}

esc_from_Z <- function(zi, n1i, n2i, pvalue = NULL, study, alternative = "two") {
  require(compute.es)
  
  if (is.null(pvalue)) {
    pvalue = 2*pnorm(abs(zi), lower.tail = F)
  }
  
  es.pes <- compute.es::pes(p = pvalue, n.1 = n1i, n.2 = n2i, dig = 10, tail = alternative)
  return(as.data.frame(list(study = study, 
                            es = es.pes$lOR, 
                            weight = NA, 
                            sample.size = es.pes$N.total,
                            se = NA, 
                            var = (pi^2*es.pes$var.d)/(3),
                            ci.lo = es.pes$l.lor, 
                            ci.hi = es.pes$u.lor,
                            measure = "logit")))
  # Source for variance compuration for log odds ratio: compute.es documentation
  # v_(log(o))= (pi^2v_(d))/ (3)
}


# Transform known effect sizes

dfES <- quiet(ldply(1:nrow(dfStatcheck_SLR_ES_N), function(x) {
  require(esc)
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "Chi2" & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) & !is.na(dfStatcheck_SLR_ES_N[[x, "ni"]]) & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "df1"]]) & (dfStatcheck_SLR_ES_N[[x, "df1"]] == 1) &
     (dfStatcheck_SLR_ES_N[[x, "ni"]]-dfStatcheck_SLR_ES_N[[x, "Value"]]) > 0) {
    return(as.data.frame(esc::esc_chisq(chisq = dfStatcheck_SLR_ES_N[[x, "Value"]],
                                        totaln = dfStatcheck_SLR_ES_N[[x, "ni"]],
                                        study = dfStatcheck_SLR_ES_N[[x, "ATag"]],
                                        es.type = "logit"))) 
  }
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "t" & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "n1i"]]) & !is.na(dfStatcheck_SLR_ES_N[[x, "n2i"]])) {
    return(as.data.frame(esc::esc_t(t = dfStatcheck_SLR_ES_N[[x, "Value"]],
                                    grp1n = dfStatcheck_SLR_ES_N[[x, "n1i"]],
                                    grp2n = dfStatcheck_SLR_ES_N[[x, "n2i"]],
                                    study = dfStatcheck_SLR_ES_N[[x, "ATag"]],
                                    es.type = "logit")))
  }
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "t" & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "ni"]])) {
    return(as.data.frame(esc::esc_t(t = dfStatcheck_SLR_ES_N[[x, "Value"]],
                                    totaln = dfStatcheck_SLR_ES_N[[x, "ni"]],
                                    study = dfStatcheck_SLR_ES_N[[x, "ATag"]],
                                    es.type = "logit")))
  }
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "t" & 
     is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "n1i"]]) & !is.na(dfStatcheck_SLR_ES_N[[x, "n2i"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "m1i"]]) & !is.na(dfStatcheck_SLR_ES_N[[x, "m2i"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "sd1i"]]) & !is.na(dfStatcheck_SLR_ES_N[[x, "sd2i"]])) {
    return(as.data.frame(esc_from_msd(m1i = dfStatcheck_SLR_ES_N[[x, "m1i"]],
                                      m2i = dfStatcheck_SLR_ES_N[[x, "m2i"]],
                                      sd1i = dfStatcheck_SLR_ES_N[[x, "sd1i"]],
                                      sd2i = dfStatcheck_SLR_ES_N[[x, "sd2i"]],
                                      n1i = dfStatcheck_SLR_ES_N[[x, "n1i"]],
                                      n2i = dfStatcheck_SLR_ES_N[[x, "n2i"]],
                                      study = dfStatcheck_SLR_ES_N[[x, "ATag"]])))
  }
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "t" & 
     is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) &
     is.na(dfStatcheck_SLR_ES_N[[x, "n1i"]]) & is.na(dfStatcheck_SLR_ES_N[[x, "n2i"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "ni"]]) & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "m1i"]]) & !is.na(dfStatcheck_SLR_ES_N[[x, "m2i"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "sd1i"]]) & !is.na(dfStatcheck_SLR_ES_N[[x, "sd2i"]])) {
    return(as.data.frame(esc_from_msd(m1i = dfStatcheck_SLR_ES_N[[x, "m1i"]],
                                      m2i = dfStatcheck_SLR_ES_N[[x, "m2i"]],
                                      sd1i = dfStatcheck_SLR_ES_N[[x, "sd1i"]],
                                      sd2i = dfStatcheck_SLR_ES_N[[x, "sd2i"]],
                                      n1i = floor(dfStatcheck_SLR_ES_N[[x, "ni"]]/2),
                                      n2i = floor(dfStatcheck_SLR_ES_N[[x, "ni"]]/2),
                                      study = dfStatcheck_SLR_ES_N[[x, "ATag"]])))
  }
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "F" & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "n1i"]]) & !is.na(dfStatcheck_SLR_ES_N[[x, "n2i"]])) {
    return(as.data.frame(esc::esc_f(f = dfStatcheck_SLR_ES_N[[x, "Value"]],
                                    grp1n = dfStatcheck_SLR_ES_N[[x, "n1i"]],
                                    grp2n = dfStatcheck_SLR_ES_N[[x, "n2i"]],
                                    study = dfStatcheck_SLR_ES_N[[x, "ATag"]],
                                    es.type = "logit")))
  }
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "F" & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "ni"]])) {
    return(as.data.frame(esc::esc_f(f = dfStatcheck_SLR_ES_N[[x, "Value"]],
                                    totaln = dfStatcheck_SLR_ES_N[[x, "ni"]],
                                    study = dfStatcheck_SLR_ES_N[[x, "ATag"]],
                                    es.type = "logit")))
  }
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "r" &
     !is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "ni"]])) {
    return(as.data.frame(esc_from_r(r = dfStatcheck_SLR_ES_N[[x, "Value"]],
                                    n = dfStatcheck_SLR_ES_N[[x, "ni"]],
                                    study = dfStatcheck_SLR_ES_N[[x, "ATag"]])))
  }
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "Z" & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) &
     !is.na(dfStatcheck_SLR_ES_N[[x, "n1i"]]) & !is.na(dfStatcheck_SLR_ES_N[[x, "n1i"]])) {
    return(as.data.frame(esc_from_Z(zi = dfStatcheck_SLR_ES_N[[x, "Value"]],
                                    n1i = dfStatcheck_SLR_ES_N[[x, "n1i"]],
                                    n2i = dfStatcheck_SLR_ES_N[[x, "n2i"]],
                                    study = dfStatcheck_SLR_ES_N[[x, "ATag"]])))
  }
  if(dfStatcheck_SLR_ES_N[[x, "Statistic"]] == "Z" & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "Value"]]) &
     is.na(dfStatcheck_SLR_ES_N[[x, "n1i"]]) & is.na(dfStatcheck_SLR_ES_N[[x, "n1i"]]) & 
     !is.na(dfStatcheck_SLR_ES_N[[x, "ni"]])) {
    return(as.data.frame(esc_from_Z(zi = dfStatcheck_SLR_ES_N[[x, "Value"]],
                                    n1i = floor(dfStatcheck_SLR_ES_N[[x, "ni"]]/2),
                                    n2i = floor(dfStatcheck_SLR_ES_N[[x, "ni"]]/2),
                                    study = dfStatcheck_SLR_ES_N[[x, "ATag"]])))
  }
  
  return(as.data.frame(list(study = dfStatcheck_SLR_ES_N[[x, "ATag"]], 
                            es = NA, 
                            weight = NA, 
                            sample.size = NA,
                            se = NA, 
                            var = NA, 
                            ci.lo = NA, 
                            ci.hi = NA,
                            measure = NA)))
}))


#
# The final table contains the extracted effect sizes in standard format for Metafor
#
dfESCombined <- dplyr::bind_cols(dfStatcheck_SLR_ES_N, dfES)


#
# Compute likelihood ratios
# 

# This analysis is based on dfESCombined

compute.LR.RBP.t <- function(t, m1i, m2i, sd1i, sd2i, n1i, n2i, FPR = 0.05) {
  df <- n1i + n2i - 2
  y0 <- dt(x = t, df = df, ncp = 0)
  
  sdpooled <- sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/(df))
  ncp <- (m2i-m1i)/sdpooled
  y1 <- dt(x = t, df = df, ncp = ncp)
  p1 = y1
  p0 = 2*y0
  
  RBP = p0*(1-FPR)/(p0*(1-FPR) + p1*FPR)
  
  return(list(df = df, ncp = ncp, y0 = y0, y1 = y1, p0 = 2*y0, p1 = y1, LR = p1/p0, RBP = RBP))
}

compute.LR.RBP.d <- function(t, d, n1i, n2i, FPR = 0.05) {
  df <- n1i + n2i - 2
  y0 <- dt(x = t, df = df, ncp = 0)
  
  ncp <- d*sqrt((n1i*n2i)/(n1i+n2i))
  y1 <- dt(x = t, df = df, ncp = ncp)
  p1 = y1
  p0 = 2*y0
  
  RBP = p0*(1-FPR)/(p0*(1-FPR) + p1*FPR)
  
  return(list(df = df, ncp = ncp, y0 = y0, y1 = y1, p0 = 2*y0, p1 = y1, LR = p1/p0, RBP = RBP))
}


#
# Analysis of Power wrt. actual observed effect sizes and MCC test families
#

# Lookup table for effect size classification as small, medium and large according to Cohen.
createThresholds <- function() {
  statistic <- c("t", "r", "Chi2", "F")
  es.type <- c("d", "r", "w", "f")
  small <- c(.20, .10, .10, .10)
  medium <- c(.50, .30, .30, .25)
  large <- c(.80, .50, .50, .40)
  df <- data.frame(statistic, es.type, small, medium, large)
  rownames(df) <- c("t", "r", "Chi2", "F")
  return (df)
}
dfThresholds <- createThresholds()


indexFunTTest2N <- function(x, threshold, sig_level) {
  require(pwr)
  pwr_t <- pwr.t2n.test(n1 = dfESCombined[x,]$n1i, n2 = dfESCombined[x,]$n2i, d = dfThresholds["t", threshold], sig.level = sig_level, alternative="two.sided")
  pwr_t_MCC <- pwr.t2n.test(n1 = dfESCombined[x,]$n1i, n2 = dfESCombined[x,]$n2i, d = dfThresholds["t", threshold], sig.level = (sig_level/dfESCombined[x,]$NMC), alternative="two.sided")
  pwr_real <- list("pwr.N" = pwr_t$n1+pwr_t$n2, "pwr.n1i" = pwr_t$n1, "pwr.n2i" = pwr_t$n2, 
                   "pwr.k" = 2, "pwr.es" = pwr_t$d, "pwr.es.type" = "d", "pwr.es.threshold" = threshold,
                   "df" = NA, "pwr.power" = pwr_t$power, "pwr.NMC" = dfESCombined[x,]$NMC, 
                   "pwr.power.MCC" = ifelse(!is.null(pwr_t_MCC), pwr_t_MCC$power, NA), 
                   "pwr.method" = pwr_t$method, "pwr.alternative" = pwr_t$alternative, "sig.level" = pwr_t$sig.level)
  return(pwr_real)
}

indexFunTTest <- function(x, threshold, sig_level) {
  require(pwr)
  pwr_t <- pwr.t.test(n = floor(dfESCombined[x,]$ni/2), d = dfThresholds["t", threshold], sig.level = sig_level, type = "two.sample", alternative="two.sided")
  pwr_t_MCC <- pwr.t.test(n = floor(dfESCombined[x,]$ni/2), d = dfThresholds["t", threshold], sig.level = (sig_level/dfESCombined[x,]$NMC), type = "two.sample", alternative="two.sided")
  pwr_real <- list("pwr.N" = pwr_t$n*2, "pwr.n1i" = pwr_t$n, "pwr.n2i" = pwr_t$n, 
                   "pwr.k" = 2, "pwr.es" = pwr_t$d, "pwr.es.type" = "d", "pwr.es.threshold" = threshold,
                   "df" = NA, "pwr.power" = pwr_t$power, "pwr.NMC" = dfESCombined[x,]$NMC, 
                   "pwr.power.MCC" = ifelse(!is.null(pwr_t_MCC), pwr_t_MCC$power, NA),
                   "pwr.method" = pwr_t$method, "pwr.alternative" = pwr_t$alternative, "sig.level" = pwr_t$sig.level)
  return(pwr_real)
}

indexFunChisq <- function(x, threshold, sig_level) {
  require(pwr)
  pwr_chisq <- pwr.chisq.test(N = dfESCombined[x,]$ni, w = dfThresholds["Chi2", threshold], df = dfESCombined[x,]$df1, sig.level = sig_level)
  pwr_chisq_MCC <- pwr.chisq.test(N = dfESCombined[x,]$ni, w = dfThresholds["Chi2", threshold], df = dfESCombined[x,]$df1, sig.level = (sig_level/dfESCombined[x,]$NMC))
  pwr_real <- list("pwr.N" = pwr_chisq$N, "pwr.n1i" = NA, "pwr.n2i" = NA, 
                   "pwr.k" = 2, "pwr.es" = pwr_chisq$w, "pwr.es.type" = "w", "pwr.es.threshold" = threshold,
                   "df" = pwr_chisq$df, "pwr.power" = pwr_chisq$power,  "pwr.NMC" = dfESCombined[x,]$NMC, 
                   "pwr.power.MCC"  = ifelse(!is.null(pwr_chisq_MCC), pwr_chisq_MCC$power, NA),
                   "pwr.method" = pwr_chisq$method, "pwr.alternative" = NA, "sig.level" = pwr_chisq$sig.level)
  return(pwr_real)
}

indexFunR <- function(x, threshold, sig_level) {
  require(pwr)
  pwr_r <- pwr.r.test(n = dfESCombined[x,]$ni, r = dfThresholds["r", threshold], sig.level = sig_level, alternative = "two.sided")
  pwr_r_MCC <- pwr.r.test(n = dfESCombined[x,]$ni, r = dfThresholds["r", threshold], sig.level = (sig_level/dfESCombined[x,]$NMC), alternative = "two.sided")
  pwr_real <- list("pwr.N" = pwr_r$n, "pwr.n1i" = NA, "pwr.n2i" = NA, 
                   "pwr.k" = 2, "pwr.es" = pwr_r$r, "pwr.es.type" = "r", "pwr.es.threshold" = threshold,
                   "df" = NA, "pwr.power" = pwr_r$power, "pwr.NMC" = dfESCombined[x,]$NMC, 
                   "pwr.power.MCC"  = ifelse(!is.null(pwr_r_MCC), pwr_r_MCC$power, NA),
                   "pwr.method" = pwr_r$method, "pwr.alternative" = pwr_r$alternative, "sig.level" = pwr_r$sig.level)
  return(pwr_real)
}

indexFunF <- function(x, threshold, sig_level) {
  require(pwr)
  pwr_F <- pwr.anova.test(n = floor(dfESCombined[x,]$ni/2), k =2, f = dfThresholds["F", threshold], sig.level = sig_level)
  pwr_F_MCC <- pwr.anova.test(n = floor(dfESCombined[x,]$ni/2), k =2, f = dfThresholds["F", threshold], sig.level = (sig_level/dfESCombined[x,]$NMC))
  pwr_real <- list("pwr.N" = pwr_F$n, "pwr.n1i" = NA, "pwr.n2i" = NA, 
                   "pwr.k" = pwr_F$k, "pwr.es" = pwr_F$f, "pwr.es.type" = "f", "pwr.es.threshold" = threshold,
                   "df" = NA, "pwr.power" = pwr_F$power, "pwr.NMC" = dfESCombined[x,]$NMC, 
                   "pwr.power.MCC"= ifelse(!is.null(pwr_F_MCC),  pwr_F_MCC$power, NA),
                   "pwr.method" = pwr_F$method, "pwr.alternative" = NA, "sig.level" = pwr_F$sig.level)
  return(pwr_real)
}

indexFunThresholdESParam <- function(x, threshold = "medium", prior = 0.5, bias = 0.2, sig_level = 0.05) {
  require(pwr)
  ID <- as.character(dfESCombined[x,]$ID)
  names(ID) <- "ID"
  
  type <- "simulated"
  names(type) <- "type"
  names(prior) <- "prior"
  names(bias) <- "bias"
  names(sig_level) <- "sig_level"
  
  
  pwr_real <- NULL
  if (dfESCombined[[x,"Statistic"]] == "t" &
      !is.na(dfESCombined[[x, "n1i"]]) & !is.na(dfESCombined[[x,"n2i"]])) {
    pwr_real <- indexFunTTest2N(x, threshold, sig_level)
  } else if (dfESCombined[[x,"Statistic"]] == "t" &
             !is.na(dfESCombined[[x,"ni"]])) {
    pwr_real <- indexFunTTest(x, threshold, sig_level)
  } else if (dfESCombined[[x,"Statistic"]] == "Chi2" &
             !is.na(dfESCombined[[x,"ni"]]) & !is.na(dfESCombined[[x,"df1"]])) {
    pwr_real <- indexFunChisq(x, threshold, sig_level)
  } else if (dfESCombined[[x,"Statistic"]] == "r" &
             !is.na(dfESCombined[[x,"ni"]])) {
    pwr_real <- indexFunR(x, threshold, sig_level)
  } else if (dfESCombined[[x,"Statistic"]] == "F" &
             !is.na(dfESCombined[[x,"ni"]])) {
    pwr_real <- indexFunF(x, threshold, sig_level)
  } else {
    
  }
  
  if (!is.null(pwr_real)) {
    if (bias == 0) {
      ppv_real <- computePPV(prior = prior, alpha = sig_level, power = pwr_real$pwr.power)
      ppv_real_MCC <- computePPV(prior = prior, alpha = (sig_level/dfESCombined[x,]$NMC), power = pwr_real$pwr.power)
    } else {
      ppv_real <- computePPVBias(prior = prior, alpha = sig_level, power = pwr_real$pwr.power, bias = bias)
      ppv_real_MCC <- computePPVBias(prior = prior, alpha = (sig_level/dfESCombined[x,]$NMC), power = pwr_real$pwr.power, bias = bias)
    }
    names(ppv_real) <- "ppv"
    
    return(list("ID" = ID, "ATag" = dfESCombined[[x,"ATag"]],
                "type" = type, 
                "N" = pwr_real$pwr.N,
                "statistic" = dfESCombined[[x,"Statistic"]], "es.type" = dfThresholds[[dfESCombined[[x,"Statistic"]], "es.type"]], 
                "es" = dfThresholds[[dfESCombined[[x,"Statistic"]], threshold]], 
                "threshold" = threshold,
                "sig_level" = sig_level, "prior" = prior, "bias" = bias, 
                "power" = pwr_real$pwr.power, 
                "power.MCC" = pwr_real$pwr.power.MCC, 
                "ppv" = ppv_real,
                "ppv.MCC" = ppv_real_MCC))
  } else {
    statistic <- dfESCombined[[x,"Statistic"]]
    return(list("ID" = ID, "ATag" = dfESCombined[[x,"ATag"]],
                "type" = type, 
                "N" = NA,
                "statistic" = statistic, "es.type" = NA, 
                "es" = NA, 
                "threshold" = threshold,
                "sig_level" = sig_level, "prior" = prior, "bias" = bias, 
                "power" = NA, "power.MCC" = NA, 
                "ppv" = NA, "ppv.MCC" = NA))
  }
}


dfESCombined <- dfESCombined %>% dplyr::group_by(ATag, Study, TestFamily) %>% dplyr::mutate(NMC = n()) %>% ungroup()

dfESMCC <- ldply(1:nrow(dfESCombined), function(x) {
  if (!is.na(dfESCombined[[x, "es"]]) & !is.na(dfESCombined[[x, "ci.lo"]]) & !is.na(dfESCombined[[x, "ci.hi"]]) & !is.na(dfESCombined[[x, "NMC"]]) ) { 
    adjCI <- convertCItoCI.MCC(es = dfESCombined[[x, "es"]], UL = dfESCombined[[x, "ci.hi"]], LL = dfESCombined[[x, "ci.lo"]], NMC = dfESCombined[[x, "NMC"]])
    return(as.data.frame(list(ci.hi.MCC = adjCI$UL.MCC, ci.lo.MCC = adjCI$LL.MCC, sig.level.MCC = adjCI$sig.level)))
  } else {
    return(as.data.frame(list(ci.hi.MCC = NA, ci.lo.MCC = NA, sig.level.MCC = NA)))
  }
})

dfESCombined <- dplyr::bind_cols(dfESCombined, dfESMCC)

dfESp <- ldply(1:nrow(dfESCombined), function(x) {
  if (!is.na(dfESCombined[[x, "es"]]) & !is.na(dfESCombined[[x, "ci.lo"]]) & !is.na(dfESCombined[[x, "ci.hi"]]) & !is.na(dfESCombined[[x, "ci.lo.MCC"]]) & !is.na(dfESCombined[[x, "ci.hi.MCC"]])) {
    est.p <- computePfromCI(est = dfESCombined[[x, "es"]], LL = dfESCombined[[x, "ci.lo"]], UL = dfESCombined[[x, "ci.hi"]])$p.value
    est.p.MCC <- computePfromCI(est = dfESCombined[[x, "es"]], LL = dfESCombined[[x, "ci.lo.MCC"]], UL = dfESCombined[[x, "ci.hi.MCC"]])$p.value
    if (!is.na(dfESCombined[[x, "Computed"]])) {
      adj.p <- p.adjust(p = dfESCombined[[x, "Computed"]], method = "bonferroni", n=dfESCombined[[x, "NMC"]])
    } else {
      adj.p <- NA
    }
    return(as.data.frame(list(adj.p = adj.p, est.p = est.p, est.p.MCC = est.p.MCC)))
  } else {
    return(as.data.frame(list(adj.p = NA, est.p = NA, est.p.MCC = NA)))
  }
})

dfESCombined <- dplyr::bind_cols(dfESCombined, dfESp)



#
# Simulation of Strength of Evidence (PPV etc.) based on observed ES
#

#
# Subsequent analysis operates on dfESCombined
# Computes a df with PPV for assumed ES thresholds.
#


assumed_thresholds <- c("small", "medium", "large")
assumed_priors <- c(0.5, 0.2, 0.1)
assumed_sig_level <- 0.05
assumed_biases <- c(0, 0.2, 0.3, 0.8)

if (exists("analysis")) rm("analysis")
for(y in assumed_thresholds) {
  for(z in assumed_priors) {
    for(u in assumed_biases) {
      newanalysis <- sapply(1:nrow(dfESCombined), function(x) {
        indexFunThresholdESParam(x, threshold = y, prior = z, bias = u, sig_level = assumed_sig_level)
      })
      if (!exists("analysis")) {
        analysis <- as.data.frame(t(newanalysis), stringsAsFactors = FALSE)
      } else {
        analysis <- rbind(analysis, as.data.frame(t(newanalysis), stringsAsFactors = FALSE))
      }
    }
  }
}

analysis2 <- as.data.frame(lapply(analysis, unlist))
analysis2$ID <- as.character(analysis2$ID)


dfESCombinedWPower <- dplyr::left_join(dfESCombined, analysis2, by = c("ID", "ATag"))
rm(analysis)
rm(analysis2)

### Likelihood-Ratio (vs. simulated power)
dfESCombinedWPower$LR <- dfESCombinedWPower$power/dfESCombinedWPower$est.p
dfESCombinedWPower$LR.MCC2 <- dfESCombinedWPower$power.MCC/dfESCombinedWPower$adj.p

FPR = 0.05
dfESCombinedWPower$RBP <- dfESCombinedWPower$est.p*(1-FPR)/(dfESCombinedWPower$est.p*(1-FPR) + dfESCombinedWPower$power*FPR)
dfESCombinedWPower$RBP.MCC <- dfESCombinedWPower$est.p.MCC*(1-FPR)/(dfESCombinedWPower$est.p.MCC*(1-FPR) + dfESCombinedWPower$power.MCC*FPR)
dfESCombinedWPower$RBP.MCC2 <- dfESCombinedWPower$est.p.MCC*(1-FPR)/(dfESCombinedWPower$est.p.MCC*(1-FPR) + dfESCombinedWPower$power*FPR)




server <- function(input, output, session) {
  visual_N <- 1250
  
  plotScatter <- function(df, xvar = "N", yvar = "1-ppv.MCC", ylabel = {if(input$graph=="ppv" || input$graph=="ppv.MCC"){"Positive Predictive Value (PPV)"}else{"False Positive Risk (FPR)"}}) {
    plot_scatter <- ggplot(df, aes_string(x=xvar, y=yvar, color="threshold", fill="threshold", "shape=as.factor(bias)", group="as.factor(threshold)")) +
      geom_jitter() +
      geom_smooth(formula = y ~ x, method = "loess") +
      geom_hline(yintercept=c(0.5), linetype="dotted") +
      scale_y_continuous(labels = scales::percent, breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits = c(0, 1)) +
      scale_color_viridis(discrete = TRUE, option = "D", name = "Effect Size Threshold") + 
      scale_fill_viridis(discrete = TRUE, option = "D", name = "Effect Size Threshold") + 
      xlab("Sample Size") + ylab(ylabel) +
      theme_bw() +
      theme(legend.position='bottom', 
            plot.background = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'))
    ggMarginal(plot_scatter, type = "density", groupColour = T, groupFill = T, alpha = 0.2, margins = "y")
  }
  
  plotHeatmap <- function(df) {
    df %>% dplyr::filter(N < visual_N) %>%
      dplyr::group_by(ATag, threshold, prior, bias) %>%
      dplyr::summarise(max.ppv.MCC = max(ppv.MCC)) %>%
      ggplot(., aes(x=as.factor(prior), y=reorder(ATag, max.ppv.MCC), fill=max.ppv.MCC, group = threshold, color = max.ppv.MCC)) +
      geom_tile() +
      scale_color_scico(palette = "vik", name = "Maximal PPV", midpoint = 0.5) +
      scale_fill_scico(palette = "vik", name = "Maximal PPV", midpoint = 0.5) +
      xlab("Prior") + ylab("") +
      theme_bw() +
      theme(legend.position='bottom',
            plot.background = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black')) +
      facet_grid(
        factor(bias) ~   
          threshold)
  }
  
  output$ppvPlot <- renderPlot({
    dfESCombinedWPower %>% dplyr::filter(N < visual_N) %>%
          dplyr::filter(threshold %in% input$threshold) %>% 
          dplyr::filter(bias %in% input$bias) %>% 
          dplyr::filter(prior %in% input$prior) %>% plotScatter(., yvar=input$graph)
  })
  
  output$heatPlot <- renderPlot({
    dfESCombinedWPower %>% dplyr::filter(N < visual_N) %>%
      dplyr::filter(threshold %in% input$thresholdHeat) %>% 
      dplyr::filter(bias %in% input$biasHeat) %>% 
      plotHeatmap()
  })
  
  analysisPlan <- reactive({
    for (x in seq(from=input$power[1], to=input$power[2], by=0.01)) {
      for (z in c(0.01, seq(from=0.05, to=0.5, by=0.05))) {
        for (u in seq(from=input$biasPlanning[1], to=input$biasPlanning[2], by = 0.1)) {
          if (!exists("analysisP")) {
            analysisP <- c(z, 
                              x,
                              u,
                              computePPVBias(prior = z, alpha = assumed_sig_level, power = x, bias = u))
          } else {
            analysisP <- rbind(analysisP, c(z, 
                                              x,
                                              u,
                                              computePPVBias(prior = z, alpha = assumed_sig_level, power = x, bias = u)))
          }
        }
      }
    }
    colnames(analysisP) <- c("prior", "power", "bias", "ppv")
    return(as.data.frame(analysisP))
  })
  
  output$planPlot <- renderPlot({
      ggplot(analysisPlan(), aes(x = prior, y=ppv, alpha=power, color=as.factor(bias), fill=as.factor(bias), group = as.factor(bias))) +
      geom_jitter() +
      geom_smooth(formula = y ~ x, method = "loess") +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      geom_hline(yintercept=c(0.5), linetype="dotted") +
      geom_vline(xintercept=as.numeric((input$priorPlanningTrue/input$priorPlanningFalse)/(input$priorPlanningTrue/input$priorPlanningFalse+1)), linetype="dashed") +
      scale_y_continuous(labels = scales::percent, breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits = c(0, 1)) +
      scale_x_continuous(labels = scales::percent, breaks = c(0, .1, .2, .3, .4, .5), limits = c(0, .6)) +
      scale_alpha_continuous(labels = scales::percent, name="Power") +
      scale_color_viridis(discrete = TRUE, option = "D", name = "Bias") + 
      scale_fill_viridis(discrete = TRUE, option = "D", name = "Bias") + 
      annotate(geom="text", x=.15, y=.05, label="No Information Threshold", color="red") + 
      xlab("Prior") + ylab("Positive Predictive Value") +
      theme_bw() +
      theme(legend.position='bottom', 
            plot.background = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'))
  })
  
  output$plannedR <- renderText({
    paste0("R = ", input$priorPlanningTrue, ":", input$priorPlanningFalse)
  })
  output$plannedPrior <- renderText({
    paste0("prior = ", weights::rd((input$priorPlanningTrue/input$priorPlanningFalse)/(input$priorPlanningTrue/input$priorPlanningFalse+1), digits = 3))
  })
}

ui <- fluidPage(
  # useShinyjs(),
  
  titlePanel("Strength of Evidence Analysis"),
  
  tabsetPanel(
    tabPanel("False Positive Risk", fluid = TRUE,
             h2("False Positive Risk"),
             p("Distribution of False Positive Risk (FPR) across studies dependent on assumed situation in reality."),
  sidebarLayout(
    sidebarPanel(
      
      selectInput(
        "graph",
        "Select the type of graph displayed:",
        choices = 
          list(`graph type` = list('FPR'='1-ppv.MCC',
                              'PPV'='ppv.MCC')),
        selected = 'FPR'
      ),
 
      h3("Parameters:"),
      
      selectInput(
        "threshold",
        "Select the effect-size theshold:",
        choices =
          list(`threshold` = list('small', 'medium', 'large')),
        multiple = TRUE,
        selected = c('small', 'medium', 'large')
      ),
      
      selectInput(
        "prior",
        "Select the anticipated prior probabiliy:",
        choices = list(`prior` = list('Confirmatory (.50)'=.50, 'Intermediate (.20)'=.20, 'Exploratory (.10)'=.10)),
        selected = c(.20)
      ),
      
      selectInput(
        "bias",
        "Select the anticipated bias:",
        choices =
          list(`bias` = list('Well-run RCT (.2)'='0.2', 'Weak RCT (.3)'='0.3', 'Biased study (.8)'='0.8')),
          selected = c(0.3)
      ),
    ),
    
    mainPanel(
      plotOutput("ppvPlot"),
      br(),
      conditionalPanel(
        condition = "input.prior == 0.10 && input.bias == 0.8",
        p(strong("Note: "), "When an exploratory study investigating many relations (prior = .10) 
          is executed in a biased fashion (u = .80), 
          then we have can have little confidence in the outcome,
          irrespective of the effect size thresholds present in reality
          and statistical power exerted by the study.")
      ),
      conditionalPanel(
        condition = "input.prior == 0.10 && input.bias == 0.3",
        p(strong("Note: "), "When an exploratory study investigating many relations (prior = .10)
          is realized as a random-controlled trial (RCT) with some biases (u = .30),
          we obsesrve a differentiation depending on effect size thresholds in reality 
          and power achieved by the study.")
      ),
      conditionalPanel(
        condition = "input.prior == 0.10 && input.bias == 0.2",
        p(strong("Note: "), "An exploratory study investigating many relations (prior = .10)
          in the form of a well-done random-controlled trial (RCT, u = .20)
          still bears a considerable false positive risk (> 50%) 
          even against large effect sizes in reality and with sufficient power achieved.")
      ),
      conditionalPanel(
        condition = "input.prior == 0.20 && input.bias == 0.8",
        p(strong("Note: "), "In an intermediate case of a study with a prior of .20,
          we find that an experiment setup with considerable biases (u = .80) diminishes
          any gains from large effect sizes or great statistical power achieved.")
      ),
      conditionalPanel(
        condition = "input.prior == 0.20 && input.bias == 0.3",
        p(strong("Note: "), "An intermediate-case study (prior = .20)
          run as a random-controlled trial (RCT) with some biases (u = .30)
          can achieve a 60% false positive risk against large effect sizes, 
          given sufficient statistical power.")
      ),
      conditionalPanel(
        condition = "input.prior == 0.20 && input.bias == 0.2",
        p(strong("Note: "), "If an intermediate-case study (prior = .20)
          is executed as well-done random-controlled trial (RCT) with a bias of .20,
          we observe that the false positive risk reaches the 50% mark for 
          medium and large effect size thresholds and sufficient statistical power.")
      ),
      conditionalPanel(
        condition = "input.prior == 0.50 && input.bias == 0.8",
        p(strong("Note: "), "In a confirmatory study with a prior of .50
          with considerable biases (u = .80), we observe that the 
          false positive risk after the study is greater than before the study,
          irrespective of the effect size threshold in reality and the statistical power achieved.
          The biases of the study muddle the waters.")
      ),
      conditionalPanel(
        condition = "input.prior == 0.50 && input.bias == 0.3",
        p(strong("Note: "), "Considering a confirmatory study (prior = .50)
          done as a random-controlled trial with some biases (u = .30),
          we find that false positive risks less than 40% can be achieved
          against medium and large effect sizes with suffient statistical power, 
          against small effect sizes with large samples.")
      ),
      conditionalPanel(
        condition = "input.prior == 0.50 && input.bias == 0.2",
        p(strong("Note: "), "In a confirmatory study (prior = .50) 
          executed as a well-done random-controlled trial (RCT) with a bias of u = .20,
          we have the strongest gain of information reaching a false positive risk 
          of 30% against medium and large effect size thresholds, given sufficient statistical power.")
      )
    )
  )), # tab panel False Positive Risk
  tabPanel("Upper-Bound PPV", fluid = TRUE,
           h2("Upper-Bound PPV"),
           p("Maximal upper-bound on the Positive Predictive Value (PPV) the studies could have achieved."),
           
           sidebarLayout(
             sidebarPanel(

               h3("Parameters:"),

               selectInput(
                 "thresholdHeat",
                 "Select the effect-size theshold:",
                 choices =
                   list(`threshold` = list('small', 'medium', 'large')),
                 multiple = TRUE,
                 selected = c('small', 'medium', 'large')
               ),

               selectInput(
                 "biasHeat",
                 "Select the anticipated bias:",
                 choices =
                   list(`bias` = list('Well-run RCT (.2)'='0.2', 'Weak RCT (.3)'='0.3', 'Biased study (.8)'='0.8')),
                 selected = c(0.3)
               ),

             ),

             mainPanel(
               plotOutput("heatPlot"),
               br(),
               conditionalPanel(
                 condition = "input.thresholdHeat == 'small' && input.biasHeat == 0.8",
                 p(strong("Note: "), "Biased studies (u = .80) cannot reach a 50% Positive Predictive Value (PPV)
                   in face of small effect size thresholds present in reality.")
               ),
               conditionalPanel(
                 condition = "input.thresholdHeat == 'small' && input.biasHeat == 0.3",
                 p(strong("Note: "), "For weak random-controlled trials (RCT) with some biases (u = .30) 
                 against small effect sizes present in reality,
                   we observe that they can achieve positive predictive values (PPV) greater than 50%
                   from a prior of .35.")
               ),
               conditionalPanel(
                 condition = "input.thresholdHeat == 'small' && input.biasHeat == 0.2",
                 p(strong("Note: "), "Well-done random controlled trials (RCT, u =.20) 
                   against small effect sizes present in reality
                   can achieve positive predictive values (PPV) greater than 50% 
                   from priors of .25.")
               ),
               conditionalPanel(
                 condition = "input.thresholdHeat == 'medium' && input.biasHeat == 0.8",
                 p(strong("Note: "), "Against medium effect size thresholds,
                   biased studies (u = .80) would not reach a positive predictive value (PPV) of 50%.")
               ),
               conditionalPanel(
                 condition = "input.thresholdHeat == 'medium' && input.biasHeat == 0.3",
                 p(strong("Note: "), "If we assume studies to be weak random-controlled trials (RCTs) with a bias of u = .30
                   executed against medium effect sizes, we observe that a number of studies could yield
                   a positive predictive value (PPV) greater than 50% from a prior of .30.")
               ),
               conditionalPanel(
                 condition = "input.thresholdHeat == 'medium' && input.biasHeat == 0.2",
                 p(strong("Note: "), "If studies in the field were well-done 
                   random-controlled trials (RCTs) with a bias of u = .20
                   pitted against medium effect sizes present in reality,
                   we would observe that studies could obtain a PPC greater than 50%
                   from priors of .20."
                 ),
                 conditionalPanel(
                   condition = "input.thresholdHeat == 'large' && input.biasHeat == 0.8",
                   p(strong("Note: "), "Even against large effect sizes, biased studies
                     could not obtain a PPV greater than 50%.")
                 ),
                 conditionalPanel(
                   condition = "input.thresholdHeat == 'large' && input.biasHeat == 0.3",
                   p(strong("Note: "), "Weak random-controlled trials (RCTs) with a bias of u = .30
                     would fare well against large effect sizes in reality,
                     reaching a PPV greater than 50% from a prior of .30.")
                 ),
                 conditionalPanel(
                   condition = "input.thresholdHeat == 'large' && input.biasHeat == 0.2",
                   p(strong("Note: "), "Well-done random controlled trials (RCTs) with a bias of u = .20
                     pitted against large effect size thresholds in reality yield positive predictive values 
                     greater than 50% from priors of .20")
                 )
               )
             )
           )
  ), # tab panel Max PPV Heatmap
  tabPanel("Planning", fluid = TRUE,
           h2("Planning"),
           p("Planning for a Positive Predictive Value."),
           sidebarLayout(
             sidebarPanel(
               
               sliderInput("power", "A priori power:",
                           min = 0, max = 1,
                           value = c(.80,.95),
                           step = 0.01),
               
               # sliderInput("es", "Effect size:",
               #             min = 0, max = 2,
               #             value = c(.20,.50),
               #             step = 0.10),
               
               sliderInput("biasPlanning", "Bias:",
                           min = 0, max = 1,
                           value = c(.20,.30),
                           step = 0.10),
               
               # selectInput(
               #   "priorPlanning",
               #   "Select the reference prior:",
               #   choices = list(`prior` = list('Confirmatory (.50)'=.50, 'Intermediate (.20)'=.20, 'Exploratory (.10)'=.10)),
               #   selected = c(.20)
               # ),
               
               h5("Anticipated Prior"), 
               
               sliderInput("priorPlanningTrue", "Likely number of true relations:",
                           min = 1, max = 100,
                           value = 1,
                           step = 1
                           ),
               
               sliderInput("priorPlanningFalse", "Likely number of false relations:",
                           min = 1, max = 100,
                           value = 1,
                           step = 1
                           ),
               
               textOutput("plannedR"), textOutput("plannedPrior"),
             ),
             
             mainPanel(
               plotOutput("planPlot")
             )
           )
  ) # tab panel Planning
  ) # tabset panel
) # fluid page

shinyApp(ui = ui, server = server)



