# Clean data
dat <- preprocess_SCD(design = "{user_design}", 
                      case = {user_caseID}, 
                      phase = {user_phaseID}, 
                      session = {user_session}, 
                      outcome = {user_outcome}, 
                      round_session = {user_round}, 
                      treatment_name = "{user_treatment}", 
                      data = dat)
