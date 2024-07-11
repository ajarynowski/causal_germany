library(tidyverse)
library(lubridate)
library(bnlearn)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(networkD3)
library(Hmisc)
library(MASS)
library(visNetwork)
library(caret)
library(pheatmap)
library(lubridate)
library(rtweet)

#l=ts_plot(twi2, "days", bins = 1)
#twi2$date <- as.Date(as.POSIXct(as.numeric(twi2$timestamp)))
# as.POSIXct(twi2$timestamp[424232])

#zrost= left_join(germanysel,dsel, by=c("date"="date"))
#zrost$day[c(3:16)]=floor((zrost$icu_patients_per_million[c(3:16)]+3*c(16:3))*40)
set.seed(1020)

plot.network <- function(structure, ht = "400px"){
  nodes.uniq <- unique(c(structure$arc_str[,1], structure$arc_str[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      shadow = TRUE)
  edges <- data.frame(from = structure$arc_str[,1],
                      to = structure$arc_str[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black",
                      width = abs(structure$arc_str[,3])/5)
  return(visNetwork(nodes, edges, height = ht, width = "100%"))
}


set.seed(1020)






#poprzednie dane juz nie potrzebne bo na nowo ladowanie
data_sel=read.csv2("data_casual_dsds.csv")
#data_sel$Twitter=zrost$day
matrix_corr=cor(data_sel[,-c(1,2)])
pheatmap(matrix_corr, display_numbers = T)

ccfvalues = ccf(data_sel$R_0, data_sel$COVIMOD, na.action = na.pass, lag.max=4)
ccfvalues = ccf(data_sel$INCIDENCE, data_sel$COVIMOD, na.action = na.pass, lag.max=4)
ccfvalues = ccf(data_sel$DEATHS, data_sel$COVIMOD, na.action = na.pass, lag.max=4)
ccfvalues = ccf(data_sel$DEATHS, data_sel$NETCHECK, na.action = na.pass, lag.max=4)

ccfvalues = ccf(data.ext$DEATHS_, data.ext$NETCHECK, na.action = na.pass, lag.max=4)
ccfvalues = ccf(data.ext$DEATHS_, data.ext$COVIMOD, na.action = na.pass, lag.max=4)
ccfvalues = ccf(data.ext$DEATHS, data.ext$COVIMOD, na.action = na.pass, lag.max=4)
ccfvalues = ccf(data.ext$DEATHS, data.ext$GOOGLE, na.action = na.pass, lag.max=4)
ccfvalues = ccf(data.ext$DEATHS_, data.ext$COVIMOD, na.action = na.pass, lag.max=4)


process <- preProcess(as.data.frame(data_sel), method=c("range"))

norm_scale <- predict(process, as.data.frame(data_sel))


shift <- function(d, k) rbind( tail(d,k), head(d,-k), deparse.level = 0 )

sh=shift(norm_scale$R_0,2)
data.del <- data.frame(R_0_=shift(norm_scale$R_0,2)[2,], DEATHS_=shift(norm_scale$DEATHS,2)[2,] )

data.del <- data.frame(R_0_=shift(norm_scale$R_0,2)[2,],INCIDENCE_=shift(norm_scale$INCIDENCE,2)[2,], DEATHS_=shift(norm_scale$DEATHS,2)[2,] )

data.ext <- cbind(norm_scale[c(2:32),], data.del)

ccfvalues = ccf(data.ext$R_0_, data.ext$COVIMOD, na.action = na.pass, lag.max=4)
ccfvalues = ccf(data.ext$R_0, data.ext$COVIMOD, na.action = na.pass, lag.max=4)


data <- data.ext[,-2]
data <- data.ext

vars <- c("COVIMOD", "NETCHECK")
output <- c("R_0", "DEATHS")
base.model <- "[COVIMOD][NETCHECK][R_0|COVIMOD:NETCHECK][DEATHS|COVIMOD:NETCHECK]"

calc.model <- function(data, 
                       vars,
                       output,
                       base.model)
{
  
  data_sel <- data[,c(vars, output)]
  structure_base <- empty.graph(c(vars, output))
  modelstring(structure_base) <- base.model
  
  blacklist = data.frame(from = c(c("R_0", "DEATHS"),rep("R_0", length(vars)),rep("DEATHS", length(vars))),
                         to =  c(c("DEATHS", "R_0"),vars,vars))
  
  all_models <- list(base = list(struct = structure_base),
                     hc = list(struct = hc(data_sel, score = "bic-g", blacklist = blacklist))
                     #iamb_cor = list(struct = iamb(x = data_sel, test = "cor")),
                     #iamb_zf = list(struct = iamb(x = data_sel, test = "zf")),
                     #pc_cor = list(struct = pc.stable(x = data_sel, test = "cor")),
                     #pc_zf = list(struct = pc.stable(x = data_sel, test = "zf"))
  )
  
  for( name in names(all_models)){
    all_models[[name]][["fit"]] <- bn.fit(all_models[[name]][["struct"]], data = data_sel)
    all_models[[name]][["bic_score"]] <- bnlearn::score(all_models[[name]][["struct"]], data = data_sel, type = "bic-g")
    all_models[[name]][["arc_str"]] <- arc.strength(all_models[[name]][["struct"]], data = data_sel, criterion = "bic-g")
    all_models[[name]][["R_0_R2"]] <- 1 - (var(all_models[[name]][["fit"]]$R_0$residuals)/var(data_sel$R_0))
    all_models[[name]][["DEATHS_R2"]] <- 1 - (var(all_models[[name]][["fit"]]$DEATHS$residuals)/var(data_sel$DEATHS))
  }
  
  lm_models <- list(
    R_0 = lm(R_0 ~. ,data[c(vars,"R_0")]),
    DEATHS = lm(DEATHS ~. ,data[c(vars,"DEATHS")])
  )
  
  summary <- data.frame(model = c("lm", names(all_models)),
                        R_0_R2 = c(summary(lm_models$R_0)$r.squared, sapply(names(all_models), FUN = function(name){all_models[[name]][["R_0_R2"]]})),
                        DEATHS_R2 = c(summary(lm_models$DEATHS)$r.squared, sapply(names(all_models), FUN = function(name){all_models[[name]][["DEATHS_R2"]]}))
  )
  
  return(list(bn=all_models,
              lm=lm_models,
              summary=summary))
}



vars <- c("COVIMOD", "NETCHECK")
output <- c("R_0", "DEATHS")
base.model <- "[COVIMOD][NETCHECK][R_0|COVIMOD:NETCHECK][DEATHS|COVIMOD:NETCHECK]"

calc.model <- function(data, 
                       vars,
                       output,
                       base.model)
{
  
  data_sel <- data[,c(vars, output)]
  structure_base <- empty.graph(c(vars, output))
  modelstring(structure_base) <- base.model
  
  blacklist = data.frame(from = c(c("R_0", "DEATHS"),rep("R_0", length(vars)),rep("DEATHS", length(vars))),
                         to =  c(c("DEATHS", "R_0"),vars,vars))
  
  all_models <- list(base = list(struct = structure_base),
                     hc = list(struct = hc(data_sel, score = "bic-g", blacklist = blacklist))
                     #iamb_cor = list(struct = iamb(x = data_sel, test = "cor")),
                     #iamb_zf = list(struct = iamb(x = data_sel, test = "zf")),
                     #pc_cor = list(struct = pc.stable(x = data_sel, test = "cor")),
                     #pc_zf = list(struct = pc.stable(x = data_sel, test = "zf"))
  )
  
  for( name in names(all_models)){
    all_models[[name]][["fit"]] <- bn.fit(all_models[[name]][["struct"]], data = data_sel)
    all_models[[name]][["bic_score"]] <- bnlearn::score(all_models[[name]][["struct"]], data = data_sel, type = "bic-g")
    all_models[[name]][["arc_str"]] <- arc.strength(all_models[[name]][["struct"]], data = data_sel, criterion = "bic-g")
    all_models[[name]][["R_0_R2"]] <- 1 - (var(all_models[[name]][["fit"]]$R_0$residuals)/var(data_sel$R_0))
    all_models[[name]][["DEATHS_R2"]] <- 1 - (var(all_models[[name]][["fit"]]$DEATHS$residuals)/var(data_sel$DEATHS))
  }
  
  lm_models <- list(
    R_0 = lm(R_0 ~. ,data[c(vars,"R_0")]),
    DEATHS = lm(DEATHS ~. ,data[c(vars,"DEATHS")])
  )
  
  summary <- data.frame(model = c("lm", names(all_models)),
                        R_0_R2 = c(summary(lm_models$R_0)$r.squared, sapply(names(all_models), FUN = function(name){all_models[[name]][["R_0_R2"]]})),
                        DEATHS_R2 = c(summary(lm_models$DEATHS)$r.squared, sapply(names(all_models), FUN = function(name){all_models[[name]][["DEATHS_R2"]]}))
  )
  
  return(list(bn=all_models,
              lm=lm_models,
              summary=summary))
}

models_H0 <- calc.model(data.ext,
                        c("COVIMOD", "NETCHECK"),
                        c("R_0", "DEATHS"),
                        "[COVIMOD][NETCHECK][R_0|COVIMOD:NETCHECK][DEATHS|COVIMOD:NETCHECK]")

models_H1 <- calc.model(data.ext,
                        c("COVIMOD", "NETCHECK", "FACEBOOK", "GOOGLE"),
                        c("R_0", "DEATHS"),
                        "[COVIMOD][NETCHECK][FACEBOOK][GOOGLE][R_0|COVIMOD:NETCHECK:FACEBOOK:GOOGLE][DEATHS|COVIMOD:NETCHECK:FACEBOOK:GOOGLE]")


#dobry model do analizy!
models_H1 <- calc.model(data.ext,
                        c("COVIMOD", "NETCHECK", "FACEBOOK", "GOOGLE","INCIDENCE","Stri","FEAR_in","FEAR_sev","INCIDENCE_","Twitter"),
                        c("R_0", "DEATHS"),
                        "[COVIMOD][Twitter][FEAR_in][FEAR_sev][Stri][INCIDENCE][NETCHECK][FACEBOOK][INCIDENCE_][GOOGLE][R_0|COVIMOD:NETCHECK:FACEBOOK:GOOGLE:Twitter:INCIDENCE_:INCIDENCE:Stri:FEAR_in:FEAR_sev][DEATHS|COVIMOD:NETCHECK:INCIDENCE_:FACEBOOK:GOOGLE:INCIDENCE:Twitter:Stri:FEAR_in:FEAR_sev]")



models_H0$summary
models_H0$bn$hc
models_H0$bn$hc$arc_str
plot.network(models_H0$bn$hc$struct)
plot.network(models_H0$bn$hc)

models_H1$summary
models_H1$bn$hc
models_H1$bn$hc$arc_str
#plot.network(models_H1$bn$hc$struct)

plot.network(models_H1$bn$hc)
