# Clean data

dat <- dat[,c("{user_parms}")]
names(dat) <- c("case","session","phase","outcome")

dat <- preprocess_SCD(design = "{user_design}", 
                      case = case, 
                      phase = phase, 
                      session = session, 
                      outcome = outcome, 
                      data = dat)

