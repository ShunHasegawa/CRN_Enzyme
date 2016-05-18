# File names with and without plate
file_names <- dir("Data/raw_data/", full.names = TRUE)

file_wt_plates    <- file_names[grepl("plate", file_names)]
file_wtout_plates <- setdiff(file_names, file_wt_plates)

# read data
dd1 <- ldply(file_wt_plates, read_enzyme_rawdata, with_plate = TRUE)
dd2 <- ldply(file_wtout_plates, read_enzyme_rawdata, with_plate = FALSE)
enzyme_raw_dd <- rbind(dd1, dd2)

# merge with sample ids
sample_id_dd <- read.csv("Data/sample_number.csv")
sample_id_dd <- within(sample_id_dd, {
  date   <- dmy(date)
  sample <- factor(sample)
})

enzyme_raw_dd <- merge(enzyme_raw_dd, sample_id_dd, 
                       by    = c("date", "Group", "plate"), 
                       all.x = TRUE)
enzyme_raw_dd <- enzyme_raw_dd[complete.cases(enzyme_raw_dd), ]

# correction
  # 25/4/2016: column12 was wrong so remove
  enzyme_raw_dd <- subsetD(enzyme_raw_dd, 
                           !(date == as.Date("2016-4-25") & WellCol == 12))

# Split to measurement and standard data frames
measure_dd <- subsetD(enzyme_raw_dd, Content == "Sample X1")
std_dd     <- subsetD(enzyme_raw_dd, grepl("Standard", as.character(Content)))

# Inspect standard--------------------------------------------------------------

# Differnet ranges of standard concentrations (unit = umol) were used
 # standard1 -> before 6/4/2016; standard2 -> after 6/4/2016
  sds <- data.frame(Content   = rep(paste0("Standard S", 1:7), 2), 
                    st_ranges = rep(paste0("standard", c(1, 2)), each = 7), 
                    MUB       = c(0, 0.000125, 0.00025, 0.0005, 0.00125, 0.0025, 0.005, 
                                  0, 0.00025, 0.0005, 0.00125, 0.0025, 0.005, 0.01))

  std_dd$st_ranges <- ifelse(std_dd$date < as.Date("2016-4-6"), "standard1", "standard2")
  
  std_dd <- merge(std_dd, sds, by = c("Content", "st_ranges"), all.x = TRUE)
  
# generate and plot summary results
  pdf("Output/Figs/Enzyme_standard_curves.pdf", width  = 4, height = 4, onefile = TRUE)
  lm_summary_dd <- ddply(std_dd, .(sample), function(x) correct_standard_curve(x))
  dev.off()

# inspect sample measurements---------------------------------------------------
  boxplot(value ~ sample, data = measure_dd)
  boxplot(value ~ sample, data = corrected_measure_dd)
  original_measure_dd <- measure_dd
  corrected_measure_dd <- ddply(measure_dd, .(sample), function(x){
    combns <- combn(4, 3)
    newx_list <- alply(combns, 2, function(y) x[y, ])
    vars <- ldply(newx_list, function(y) var(y$value))
    min_vars <- which(vars$V1 == min(vars$V1))
    new_value <- mean(newx_list[[min_vars]]$value)
    # return(data.frame(new_value))
    return(newx_list[[min_vars]])
  })
  
  sm_dd <- merge(lm_summary_dd, corrected_measure_dd, by = "sample")
  sm_dd$MUB <- with(sm_dd, (value - intercept)/slope)
  plot(MUB ~ sample, data = sm_dd)
  
  bset(measure_dd, sample == 1)
