# Author: Khang Tong, as part of the project "Understanding the Efficacy of Phishing Training in Practice"

library(lme4)
library(dplyr)

# Formatting functions
pvalFormat = function(p.values, method = 'none', replace = FALSE, math = TRUE){
  ## Formats p-values for reports, can report adjusted pvalues
  ##    Inputs:
  ##       - p.value: numeric p-value
  ##       - method: pvalue adjustment, passed to p.adjust.methods
  ##       - replace: if TRUE, replaces p-values with their adjusted value
  ##    Outputs:
  ##       - out: formatted p-value
  
  p.values <- suppressWarnings(as.numeric(p.values))
  out      <- rep(NA, length(p.values))
  sig      <- p.adjust(p.values, method)
  if(replace) p.values <- sig
  
  for(i in 1:length(p.values)){
    if(is.na(p.values[i])){out[i] <- NA}else{
      if(p.values[i] >= .001){
        out[i] <- paste('$', formatC(p.values[i], format = 'f', digits = 3), '$', sep = '')
      }
      
      if(p.values[i] < .001){
        out[i] <- '< $0.001$'
      }
      
      if(sig[i] > 0.01 & sig[i] <= 0.05){
        out[i] <- paste(out[i], '*', sep = '')
      }
      
      if(sig[i] > 0.001 & sig[i] <= 0.01) {
        out[i] <- paste(out[i], '**', sep = '')
      }
      
      if(sig[i] <= 0.001){
        out[i] <- paste(out[i], '***', sep = '')
      }}
  }
  
  out[is.na(out)] <- '-'
  return(out)
}

printMod = function(mod, mod.null = NULL, row.names = NULL, caption = NULL, digits = 3, tex = F, gee.se = 'sandwich',
         overall.test = F, ci = F, ci.profile = T, level = .95, d.f. = F, cat.cov.test = F, lmer.norm.p = F,
         csv = NULL, ...){
  ## A function to print regression tables in LaTeX and to csv (truncated to only include glmerMod needed for this file). 
  ##    Inputs:
  ##       - mod: regression model
  ##       - row.names: a character vector, giving row names for the final printed table
  ##       - caption: a character string, caption passed to latex(...)     
  ##       - digits: number of digits to print
  ##       - tex: logical, should the latex(...) function be run, or should a matrix be returned
  ##       - gee.se: a character string, 'sandwich' or 'naive'
  ##       - overall.test: logical, likelihood ratio / conditional F test for overall model fit
  ##       - ci: logical, should confidence intervals be printed for each parameter
  ##       - ci.profile: logical, should profile likelihood be used for computing confidence intervals
  ##       - level: level of confidence interval, ignored unless ci = T
  ##       - d.f.: logical, should the degress of freedom be printed in the regression table - currently unemplemented
  ##       - cat.cov.test: logical, should categorical covariates be tested with something like ANOVA / LRT
  ##       - lmer.norm.p: logical, should a normal approximation be used for p-values in lmer() models
  ##       - ...: additional parameters passed to latex(...)
  
  ## Cleaning and defining unique model classes
  if(length(class(mod)) == 1){
    mC <- class(mod)
  }

  if(mC == 'glmerMod'){
    out <- summary(mod)$coefficients
    out.t <- out[, 3]
    out[, c(1:3)] <- formatC(out[, c(1:3), drop = F], format = 'f', digits = digits)
    out[, 1:3] <- paste('$', out[, 1:3, drop = F], '$', sep = '')
    out[, 4] <- pvalFormat(out[, 4])
    colnames(out) <- c('Estimate', 'Std. Error', 'z-value', 'p-value')
    
    if(summary(mod)$family == 'binomial'){
      if(summary(mod)$link == 'logit'){
        ccs <- summary(mod)$coefficients
        ints <- cbind(ccs[, 1] - abs(qnorm((1 - level) / 2)) * ccs[, 2], ccs[, 1] + abs(qnorm((1 - level) / 2)) * ccs[, 2])
        CI <- exp(cbind(OR = ccs[, 1], ints))
        CI <- cbind(paste('$', formatC(CI[, 1], format = 'f', digits = 3), '$', sep = ''),
                    paste('($', formatC(CI[, 2], format = 'f', digits = 3), '$, $',
                          formatC(CI[, 3], format = 'f', digits = 3), '$)', sep = ''))
        out <- cbind(out, CI); out[1, 5:6] <- c('-', '-')
        colnames(out)[5:6] <- c('Odds Ratio', paste(round(level * 100), '\\% CI', sep = ''))
      } else {
        if(ci) {
          ccs <- summary(mod)$coefficients
          ints <- cbind(ccs[, 1] - abs(qnorm((1 - level) / 2)) * ccs[, 2], ccs[, 1] + abs(qnorm((1 - level) / 2)) * ccs[, 2])
          ints <- paste('($', formatC(ints[, 1], digits = digits, format = 'f'), '$, $', 
                        formatC(ints[, 2], digits = digits, format = 'f'), '$)', sep = '')
          out <- cbind(out[, 1:2], ints, out[, 3:ncol(out)])
          colnames(out)[3] <- paste(round(level * 100), '\\% CI', sep = '')
        }
      }
    } else {
      if(ci){
        ccs <- summary(mod)$coefficients
        ints <- cbind(ccs[, 1] - abs(qnorm((1 - level) / 2)) * ccs[, 2], ccs[, 1] + abs(qnorm((1 - level) / 2)) * ccs[, 2])
        ints <- paste('($', formatC(ints[, 1], digits = digits, format = 'f'), '$, $', 
                      formatC(ints[, 2], digits = digits, format = 'f'), '$)', sep = '')
        out <- cbind(out[, 1:2], ints, out[, 3:ncol(out), drop = F])
        colnames(out)[3] <- paste(round(level * 100), '\\% CI', sep = '')
      }
    }
    
    if(cat.cov.test){
      tt1 <- rep('', nrow(out))
      tt2 <- rep('', nrow(out))
      tAOV <- anova(mod)
      int.cor <- ifelse(rownames(out)[1] == '(Intercept)', 2, 1)
      for(i in which(tAOV$npar > 1)){
        tAOV. <- anova(update(mod), update(mod, as.formula(paste('. ~ . -', rownames(tAOV)[i]))))
        tt1[sum(tAOV$npar[1:i]) + int.cor - tAOV$npar[i]] <- paste('$\\chi^{2}_{', tAOV.$'Chi Df'[2], 
                                                                   '} = ', formatC(tAOV.$Chisq[2], format = 'f', digits = digits), '$', sep = '')
        tt2[sum(tAOV$npar[1:i]) + int.cor - tAOV$npar[i]] <- pvalFormat(tAOV.[2, 'Pr(>Chisq)'])
      }
      out <- cbind(out, tt1, tt2)
      colnames(out)[(ncol(out) - 1):ncol(out)] <- c('$\\chi^2$ value', 'p-value')
    }
    
    if(!is.null(mod.null)){
      lrt_stat <- anova(mod, mod.null)
      out <- cbind(out, c('', paste('$\\chi^2_{', abs(diff(lrt_stat$Df)), '} = ', formatC(lrt_stat$Chisq[2], format = 'f', digits = digits), 
                                    '$', sep = ''), rep('', nrow(out) - 2)))
      out <- cbind(out, c('', pvalFormat(lrt_stat$'Pr(>Chisq)'[2]),
                          rep('', nrow(out) - 2)))
      colnames(out)[c(ncol(out) - 1, ncol(out))] <- c('$\\Delta$Dev.', 'p-value')
    }
    
    if(overall.test){
      LRT <- anova(update(mod), update(mod, as.formula(paste('. ~ . - (', paste(rownames(anova(mod)), collapse = ' + '), ')'))))
      out <- cbind(out, c('', paste('$\\chi^2_{', abs(diff(LRT$Df)), '} = ', 
                                    formatC(LRT$Chisq[2], format = 'f', digits = 2), '$', sep = ''),
                          rep('', nrow(out) - 2)))
      out <- cbind(out, c('', pvalFormat(LRT$'Pr(>Chisq)'[2]),
                          rep('', nrow(out) - 2)))
      colnames(out)[(ncol(out) - 1):ncol(out)] <- c('Overall $\\chi^2$', 'Pr($>$$\\chi^2$)')
    }  
  }
  if(!is.null(row.names)) rownames(out) <- row.names
  if(tex){
    latex(out,
          file = '',
          title = '',
          where = '!htp',
          col.just = rep('c', ncol(out)),
          caption = caption,
          insert.bottom = "Significance codes: ***$0.001$, **$0.01$, *$0.05$.",
          ...)
  } else{return(out)}
  if(!is.null(csv)) write.csv(gsub('$', '', out, fixed = T), file = csv, quote = F)
}

# 4. Annual Security Awareness Training
fit.glmer.q1 = glmer(fail ~ days.since.annual.training + cumul.prev.fail + 
                            month + phish.label + track + (1 | user), 
                     data = df, family = binomial)

# 5. Embedded Phishing Training
fit.glmer.q2 = glmer(fail ~ treatment + cumul.prev.fail + month + phish.label + 
                            track + (1 | user), 
                     data=df, family = binomial)

# 6. Embedded Training Engagement

## OVERALL

### Cumulative number of acknowledgements
fit.glmer.q3.cumul.ack = glmer(fail ~ cumul.ack + cumul.prev.fail + month + 
                                      phish.label + track + (1 | user), 
                               data = df, family = binomial)

### Completed at least 1 previous training
fit.glmer.q3.prev.ack = glmer(fail ~ prev.ack + cumul.prev.fail + month + 
                                     phish.label + track + (1 | user), 
                              data = df, family = binomial)

### Cumulative time training (30s)
fit.glmer.q3.train.time = glmer(fail ~ cumul.training.time30 + cumul.prev.fail + 
                                       month + phish.label + track + (1 | user), 
                                data = df, family = binomial)


## Models with interaction term between training type and engagement variable
fit.glmer.interaction.cumul.ack = glmer(fail ~ cumul.ack*training.type2 + cumul.prev.fail + 
                                               month + phish.label + track + (1 | user), 
                                        data = df, family = binomial)
tbl.fit.glmer.interaction.cumul.ack = printMod(fit.glmer.interaction.cumul.ack, cat.cov.test=T)

fit.glmer.interaction.prev.ack = glmer(fail ~ prev.ack*training.type2 + cumul.prev.fail + 
                                              month + phish.label + track + (1 | user), 
                                       data = df, family = binomial)
tbl.fit.glmer.interaction.prev.ack = printMod(fit.glmer.interaction.prev.ack, cat.cov.test=T)

fit.glmer.interaction.train.time = glmer(fail ~ cumul.training.time30*training.type2 + cumul.prev.fail + 
                                                month + phish.label + track + (1 | user), 
                                         data = df, family = binomial)
tbl.fit.glmer.interaction.train.time = printMod(fit.glmer.interaction.train.time, cat.cov.test=T)

### Generate summary table for interaction models
est.interactive.cumul.ack = summary(fit.glmer.interaction.cumul.ack)$coefficients[2,1] + summary(fit.glmer.interaction.cumul.ack)$coefficients[21,1]
or.interactive.cumul.ack = paste0("$", exp(est.interactive.cumul.ack) %>% round(3) %>% format(nsmall=2), "$")
se.interactive.cumul.ack = sqrt(vcov(fit.glmer.interaction.cumul.ack)[2,2] + vcov(fit.glmer.interaction.cumul.ack)[21,21] + 2*vcov(fit.glmer.interaction.cumul.ack)[2,21])
ci.interactive.cumul.ack = paste0("($", exp(est.interactive.cumul.ack - 1.965*se.interactive.cumul.ack) %>% round(3) %>% format(nsmall=2), "$, $", exp(est.interactive.cumul.ack + 1.965*se.interactive.cumul.ack) %>% round(3) %>% format(nsmall=2), "$)")

est.interactive.prev.ack = summary(fit.glmer.interaction.prev.ack)$coefficients[2,1] + summary(fit.glmer.interaction.prev.ack)$coefficients[21,1]
or.interactive.prev.ack = paste0("$", exp(est.interactive.prev.ack) %>% round(3) %>% format(nsmall=2), "$")
se.interactive.prev.ack = sqrt(vcov(fit.glmer.interaction.prev.ack)[2,2] + vcov(fit.glmer.interaction.prev.ack)[21,21] + 2*vcov(fit.glmer.interaction.prev.ack)[2,21])
ci.interactive.prev.ack = paste0("($", exp(est.interactive.prev.ack - 1.965*se.interactive.prev.ack) %>% round(3) %>% format(nsmall=2), "$, $", exp(est.interactive.prev.ack + 1.965*se.interactive.prev.ack) %>% round(3) %>% format(nsmall=2), "$)")

est.interactive.train.time = summary(fit.glmer.interaction.train.time)$coefficients[2,1] + summary(fit.glmer.interaction.train.time)$coefficients[21,1]
or.interactive.train.time = paste0("$", exp(est.interactive.train.time) %>% round(3) %>% format(nsmall=2), "$")
se.interactive.train.time = sqrt(vcov(fit.glmer.interaction.train.time)[2,2] + vcov(fit.glmer.interaction.train.time)[21,21] + 2*vcov(fit.glmer.interaction.train.time)[2,21])
ci.interactive.train.time = paste0("($", exp(est.interactive.train.time - 1.965*se.interactive.train.time) %>% round(3) %>% format(nsmall=2), "$, $", exp(est.interactive.train.time + 1.965*se.interactive.train.time) %>% round(3) %>% format(nsmall=2), "$)")


tbl.summary = data.frame(
  or.overall = c(tbl.fit.glmer.q3.prev.ack["prev.ackYes", "Odds Ratio"],
                 tbl.fit.glmer.q3.cumul.ack["cum.ack", "Odds Ratio"],
                 tbl.fit.glmer.q3.train.time["cum.training.time30", "Odds Ratio"]),
  ci.overall = c(tbl.fit.glmer.q3.prev.ack["prev.ackYes", "95\\% CI"],
                 tbl.fit.glmer.q3.cumul.ack["cum.ack", "95\\% CI"],
                 tbl.fit.glmer.q3.train.time["cum.training.time30", "95\\% CI"]),
  p.overall = c(tbl.fit.glmer.q3.prev.ack["prev.ackYes", "p-value"],
                tbl.fit.glmer.q3.cumul.ack["cum.ack", "p-value"],
                tbl.fit.glmer.q3.train.time["cum.training.time30", "p-value"]),
  or.static = c(tbl.fit.glmer.interaction.prev.ack["prev.ackYes", "Odds Ratio"],
                tbl.fit.glmer.interaction.cumul.ack["cumul.ack", "Odds Ratio"],
                tbl.fit.glmer.interaction.train.time["cumul.train.time30", "Odds Ratio"]),
  ci.static = c(tbl.fit.glmer.interaction.prev.ack["prev.ackYes", "95\\% CI"],
                tbl.fit.glmer.interaction.cumul.ack["cumul.ack", "95\\% CI"],
                tbl.fit.glmer.interaction.train.time["cumul.train.time30", "95\\% CI"]),
  or.interactive = c(or.interactive.prev.ack,
                     or.interactive.cumul.ack,
                     or.interactive.cumul.train),
  ci.interactive = c(ci.interactive.prev.ack,
                     ci.interactive.cumul.ack,
                     ci.interactive.cumul.train),
  p.interactive = c(tbl.fit.glmer.interaction.prev.ack["prev.ackYes:training.type2interactive", "p-value"],
                    tbl.fit.glmer.interaction.cumul.ack["cumul.ack:training.type2interactive", "p-value"],
                    tbl.fit.glmer.interaction.train.time["cumul.train.time30:training.type2interactive", "p-value"]))
rownames(tbl.summary) = c("prev.ackYes", "cumul.ack", "cumul.train.time30")
colnames(tbl.summary) = c("OR","95\\% CI", "P-value", "OR", "95\\% CI", "OR", "95\\% CI", "P-value for difference in OR")


df.n = data.frame(training.type = c("Overall", "Static", "Interactive"), 
                  n.recs = 
                    c(dim(df)[1],
                      dim(df %>% filter(training.type %in% c("generic.static", "contextual.static")))[1], 
                      dim(df %>% filter(training.type %in% c("generic.interactive", "contextual.interactive")))[1]), 
                  n.id = 
                    c(dim(unique(df%>%select(user)))[1], 
                      dim(unique(df%>% filter(training.type %in% c("generic.static", "contextual.static")) %>%select(user)))[1], 
                      dim(unique(df%>% filter(training.type %in% c("generic.interactive", "contextual.interactive")) %>%select(user)))[1]))

clabel2 = c(
  paste0("Overall (n = ", df.n$n.id[1], ", ", df.n$n.recs[1], ")"),
  paste0("Static (n = ", df.n$n.id[2], ", ", df.n$n.recs[2], ")"),
  paste0("Interactive (n = ", df.n$n.id[3], ", ", df.n$n.recs[3], ")"))

# Power analysis
p1 = seq(0.1, 0.5, by=0.05) #Baseline failure rate
delta = 0.05 # Minimum detectable effect
alpha = 0.05 
power = 0.9

n=c()
for (i in 1:length(p1)){
  out = power.prop.test(p1=p1[i], p2=p1[i]-delta, sig.level = alpha, power=power, alternative="two.sided")
  n = c(n, out$n)
}

table.power = data.frame(p = p1, n = n %>% ceiling())
colnames(table.power) = c("Baseline failure rate", "Sample size per group")

n=c()
for (i in 1:length(p1)){
  out = power.prop.test(p1=p1[i], p2=p1[i]-delta, sig.level = alpha/3, power=power, alternative="two.sided")
  n = c(n, out$n)
}

table.power2 = data.frame(p = p1, n = n %>% ceiling())
colnames(table.power2) = c("Baseline failure rate", "Sample size per group")
