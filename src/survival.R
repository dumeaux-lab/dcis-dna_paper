require(rlang)
require(maftools)
require(tidyr)
require(forcats)
require(dplyr)
require(broom)
require(patchwork)


# make a function to perform survival group analysys with different cohorts
survival.analysis <- function(maf_data, time="time.10", event="event.10", 
                              geneSetSize = 1, genecombo = 'intersection', minMutations = 1,
                              genes=NULL, 
                              clinsubset = NULL, top = 50, pvalue = 0.05, output = "") {
  
  
  if (!is.null(clinsubset)) {
  maf <- subsetMaf(maf_data, clinQuery = clinsubset)} else {maf <- maf_data}
  
  # Apply survGroup function
  sg_result <- mysurvGroup(maf = maf, top = top, geneSetSize = geneSetSize , genecombo = genecombo, genes=genes, time = time, Status = event, verbose = FALSE)
  genes <- sg_result$restb$Gene_combination[sg_result$restb$P_value < pvalue & sg_result$restb$Mutant >= min(c(5, nrow(maf@clinical.data)*0.10))]
  
  if (length(genes)>0){
     pdf(paste0("../../figures/exp07/survplot_top", top, "_", genecombo, "_set", geneSetSize, "_sig_", clinsubset, ".pdf"))
     for (i in genes){
       km <- ggsurvplot(
         sg_result$res[[i]]$fit.km,# survfit object with calculated statistics.
         data=sg_result$res[[i]]$data,
         pval = TRUE,             # show p-value of log-rank test.
         conf.int = TRUE,         # show confidence intervals for
         # point estimaes of survival curves.
         conf.int.style = "step",  # customize style of confidence intervals
         xlab = "Time in years",   # customize X axis label.
         break.time.by = 3,     # break X axis in time intervals by 200.
         ggtheme = theme_light(), # customize plot and risk table with a theme.
         risk.table = "abs_pct",  # absolute number and percentage at risk.
         risk.table.y.text.col = T,# colour risk table text annotations.
         risk.table.y.text = FALSE,# show bars instead of names in text annotations
         # in legend of risk table.
         ncensor.plot = TRUE,      # plot the number of censored subjects at time t
         surv.median.line = "hv",  # add the median survival pointer.
         title = i)
       print(km)
       }
     dev.off()
     
     pp <- list()
     pdf(paste0("../../figures/exp07/hr_", top, "_", genecombo, "_set", geneSetSize, "_sig_", clinsubset, ".pdf"), height = 40, width = 20)
     #for (j in genes){
       # Create the forest plot
      pp <- ggplot(sg_result$restb[sg_result$restb$Gene_combination %in% genes,], aes(y = Gene_combination)) +
         geom_point(aes(x = hr), size = 3) +
         geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
         geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
         scale_x_log10() +
         theme_bw() +
         theme(
           panel.grid.minor = element_blank(),
           panel.grid.major.y = element_blank()
         ) +
         labs(
           x = "",
           y = "",
           title = ""
         ) +
         # Add HR and CI values as text
         geom_text(aes(x = max(conf.high) * 1.2, 
                       label = sprintf("%.2f (%.2f-%.2f)", hr, conf.low, conf.high)),
                   hjust = 0) +
         # Add p-values
         geom_text(aes(x = max(conf.high) * 2, 
                       label = sprintf("p = %.3f", P_value)),
                   hjust = 0)
      
      print(pp)
    dev.off()
     
     } else {pp <- NULL}
  
  return(list(surv_result=sg_result, maf=maf, genes=genes, pp=pp))
}

my_run_surv = function(cd, tsbs){
  groupNames = c("Mutant", "WT")
  col = c('maroon', 'royalblue')
  cd$Group = ifelse(test = cd$Tumor_Sample_Barcode %in% tsbs, yes = groupNames[1], no = groupNames[2])
  cd <- cd %>% 
    mutate(Group = Group %>%
             fct_relevel("WT"))
  
  surv.km = survival::survfit(formula = survival::Surv(time = Time, event = Status) ~ Group, data = cd, conf.type = "log-log")

  res = summary(surv.km)

  surv.diff = survival::survdiff(formula = survival::Surv(time = Time, event = Status) ~ Group, data = cd)
  surv.diff.pval = signif(1 - pchisq(surv.diff$chisq, length(surv.diff$n) - 1), digits = 3)

  surv.cox = coxphf::coxphf(formula = survival::Surv(time = Time, event = Status) ~ Group, data = cd)
  surv.cox.tb <- data.frame(
    term = names(surv.cox$coefficients),
    estimate = exp(surv.cox$coefficients),
    conf.low = surv.cox$ci.lower,
    conf.high = surv.cox$ci.upper,
    p.value = surv.cox$prob
  ) |>
    dplyr::rename(Group = term, hr = estimate, P_value = p.value) |>
    mutate(WT = nrow(cd[cd$Group == "WT",]),
           Mutant = nrow(cd[cd$Group == "Mutant",]))
  
  # First, ensure your data frame is properly formatted
  surv.cox.tb <- surv.cox.tb |>
    mutate(
      # Create label for the plot
      # label = Gene_combination,
      # Calculate standard error
      se = (log(conf.high) - log(conf.low))/(2*1.96)
    )


  surv.dat <- data.table::data.table(surv.cox.tb)
  surv.dat$Group = gsub(pattern = 'Group', replacement = '', x = surv.dat$Group)
  
  return(list(fit.km=surv.km, fit.cox= surv.cox, surv.diff= surv.diff, data=cd, surv.tb=surv.dat))
}


mysurvGroup = function(maf=maf, top = top, genes = genes, geneSetSize = geneSetSize, genecombo = genecombo, minMutations = minMutations
                       minSamples = 5, clinicalData = NULL, time=time, Status=event, verbose = TRUE, pvalue=0.05){
  
  if(is.null(genes)){
    genes = getGeneSummary(x = maf)[1:top, Hugo_Symbol]
  }
  

  genesCombn = combn(x = genes, m = geneSetSize)
  
  if(verbose){
    cat("------------------\n")
    cat(paste0("genes: ", length(genes), "\n"))
    cat(paste0("geneset size: ", geneSetSize, "\n"))
    cat(paste0(ncol(genesCombn), " combinations\n"))
  }
  
  if(is.null(clinicalData)){
    if(verbose){
      message("Looking for clinical data in annoatation slot of MAF..")
    }
    
    clinicalData = data.table::copy(getClinicalData(x = maf))
    clinicalData = data.table::setDT(clinicalData)
  }else{
    clinicalData = data.table::setDT(clinicalData)
  }
  
  if(!"Tumor_Sample_Barcode" %in% colnames(clinicalData)){
    print(colnames(clinicalData))
    stop("Column Tumo_Sample_Barcode not found in clinical data. Check column names and rename it to Tumo_Sample_Barcode if necessary.")
  }
  
  if(length(colnames(clinicalData)[colnames(clinicalData) %in% time]) == 0){
    print(colnames(clinicalData))
    stop(paste0(time, " not found in clinicalData. Use argument time to povide column name containing time to event."))
  }else{
    colnames(clinicalData)[colnames(clinicalData) %in% time] = 'Time'
  }
  
  if(length(colnames(clinicalData)[colnames(clinicalData) %in% Status]) == 0){
    print(colnames(clinicalData))
    stop(paste0(Status, " not found in clinicalData. Use argument Status to povide column name containing events (Dead or Alive)."))
  }else{
    colnames(clinicalData)[colnames(clinicalData) %in% Status] = 'Status'
  }
  
  clinicalData$Time = suppressWarnings(as.numeric(as.character(clinicalData$Time)) )
  clinicalData$Status = suppressWarnings(as.integer(clinicalData$Status))
  clinicalData$Time = ifelse(test = is.infinite(clinicalData$Time), yes = NA, no = clinicalData$Time)
  if(nrow(clinicalData[!is.na(Time)][!is.na(Status)]) < nrow(clinicalData)){
    message(paste0("Removed ", nrow(clinicalData) - nrow(clinicalData[!is.na(Time)][!is.na(Status)]),
                   " samples with NA's"))
    clinicalData = clinicalData[!is.na(Time)][!is.na(Status)]
  }
  
  om = maftools:::createOncoMatrix(m = maf, g = genes)
  all.tsbs = levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])
  
  mutMat = t(om$numericMatrix)
  missing.tsbs = all.tsbs[!all.tsbs %in% rownames(mutMat)]
  
  mutMat[mutMat > 0 ] = 1
  
  if(genecombo == 'union'){
    res = lapply(seq_along(1:ncol(genesCombn)), function(i){
      x = genesCombn[,i]
      mm = mutMat[,x, drop = FALSE]
      genesTSB = names(which(rowSums(mm) >= minMutations))
      if(length(genesTSB) >= minSamples){
        if(verbose){
          cat("Union Geneset: ", paste0(x, collapse = ","), "[N=", length(genesTSB),"]\n")
        }
        surv.dat = my_run_surv(cd = clinicalData, tsbs = genesTSB)
      }else{
        surv.dat = NULL
      }
      surv.dat
    })
  } else if (genecombo == 'intersection'){
    res = lapply(seq_along(1:ncol(genesCombn)), function(i){
    x = genesCombn[,i]
    mm = mutMat[,x, drop = FALSE]
    genesTSB = names(which(rowSums(mm) == geneSetSize))
    if(length(genesTSB) >= minSamples){
      if(verbose){
        cat("Intersection Geneset: ", paste0(x, collapse = ","), "[N=", length(genesTSB),"]\n")
      }
      surv.dat = my_run_surv(cd = clinicalData, tsbs = genesTSB)
    }else{
      surv.dat = NULL
    }
    surv.dat
  })
  } else {
    stop("Invalid value for genecombo: must be 'union' or 'intersection'")
  }
  
  names(res) = apply(genesCombn, 2, paste, collapse = '_')
  resfit <- lapply(res, "[[", "surv.tb")
  restb = data.table::rbindlist(l = resfit, idcol = "Gene_combination")
  restb = restb[order(P_value, decreasing = FALSE)]
  
  
  return(list(res=res, restb=restb))
}
