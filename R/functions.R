######################### ------------------------------------------------------
# Subset and droplevels #
#########################
  subsetD <- function(...) droplevels(subset(...))
  
######################################################## -----------------------
# Read raw data of enzyme assay from microplate reader #
########################################################  

  # Thie function removes unnecessary elemetns, adds plate numbers and
  # measurement dates taken from file names.
  
  read_enzyme_rawdata <- function(file_names, with_plate){
    col_names <- c("WellRow", "WellCol", "Content", "Group", "value")
    
    # Read data
    rawdd <- read.xlsx(file_names, 
                       sheetName  = "End point", 
                       header     = FALSE,
                       colIndex   = 1:5,
                       stringsAsFactors = FALSE)
    
    # identify rows to be used; find the row including "Content", below which 96 raws
    # are data to be used
    row_val <- which(rawdd == "Content", arr.ind = TRUE)[1, 1]
    rawdd_data <- rawdd[(row_val + 1):(row_val + 96), 1:5]
    colnames(rawdd_data) <- col_names
    
    # organise data frame
    rawdd_data <- within(rawdd_data, {
      WellCol  <- factor(WellCol)
      WellRow  <- factor(WellRow)
      date     <- ymd(substring(file_names, 16, 23))
      value    <- as.numeric(value)
      Group    <- factor(Group)
      Content <- factor(Content)
    })
    rawdd_data$plate <- factor(
      ifelse(with_plate, substring(file_names, 30, 30), "1")
      )
    
    return(rawdd_data)
  }
 
###################### ---------------------------------------------------------   
# Calibration curves #
######################
  
  # Plot result
  plot_lm_results <- function(x){
    m  <- lm(value ~ MUB, data = x)
    r2 <- round(summary(m)$r.squared, 4)
    plot(value ~ MUB, 
         data  = x, 
         main  = paste0("Sample", unique(x$sample), "\n r2 =", r2))
    abline(m)
  }
  
  # return summary results
  return_lm_results <- function(x) {
    m             <- lm(value ~ MUB, data = x)
    coefm         <- coef(m)
    r2            <- summary(m)$r.squared
    msd           <- data.frame(t(coef(m)), r2)
    colnames(msd) <- c("intercept", "slope", "r2")
    return(msd)
  }
  
  # re-calculate r2 for by dropping some observations when r2 > 0.98
  recalculate_r2 <- function(x, remove_element, ...) {
    # x              ; origianl data frame
    # remove_element ; number of element to be removed
    # combs          ; all combinations of elements after removing element(s) as defined above
    combs       <- combn(1:nrow(x), nrow(x) - remove_element)
    newx_list   <- alply(combs, 2, function(y) x[y, ])
    new_results <- ldply(newx_list, 
                         function(y) return_lm_results(y))
    max_r2      <- which(new_results$r2 == max(new_results$r2))
    new_res     <- new_results[max_r2, -1]
    points(value ~ MUB, data = newx_list[[max_r2]], pch = 19, ...)
    abline(lm(value ~ MUB, data = newx_list[[max_r2]]), lty = 2, ...)
    mtext(paste0("r2 = ", round(max(new_results$r2), 4)), 3, col = ...)
    return(new_res)
  } 
  
  #  summarise parameters for calibration curves using above functions to obtain
  #  a curve with r2 > 0.98
  correct_standard_curve <- function(x){
    lm_results <- return_lm_results(x)
    plot_lm_results(x)
    if(lm_results$r2 > 0.98) {
      return(lm_results)
    } else {
      new_res <- recalculate_r2(x, remove_element = 1, col = "red")
      if(new_res$r2 > 0.98){
        return(new_res)
      } else {
        new_res2 <- recalculate_r2(x, remove_element = 2, col = "blue", line = -1)
        return(new_res2)
      }
    }
  }
  
  
  
  
  