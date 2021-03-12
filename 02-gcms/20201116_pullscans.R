scans2pull <- c(
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200214/JWL61_blanked.mzxml" = 1183,
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200609/JWL97_blanked.mzxml" = 1148,
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200612/JWL115_blanked.mzxml" = 898,
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200212/JWL29_blanked.mzxml" = 878,
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200212/JWL16_blanked.mzxml" = 1182,
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200212/JWL17_blanked.mzxml" = 841,
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200608/JWL69_blanked.mzxml" = 1160
)

scans2pull %>%
  names() %>%
  pblapply(
    cl = 8L,
    .,
    function(file_data){
      scannum <- scans2pull[[file_data]]
      file_data %>%
        read_tidymass() %>%
        filter(scan == scannum)
    }
  ) %>%
  do.call(rbind, .) %>%
  write_tidymass("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/ocfas/ocfas.mzxml")

