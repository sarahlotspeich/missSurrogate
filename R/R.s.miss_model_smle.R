# Sieve maximum likelihood estimator -- Only type = "model"
R.s.miss_model_smle = function(sone, szero, yone, yzero, conv.res, max.it = 1E4, tol = 1E-3, full.output = FALSE, orig.smle) {
  if (orig.smle) {
    R.s.miss_model_smle_original(sone = sone, 
                                 szero = szero, 
                                 yone = yone, 
                                 yzero = yzero, 
                                 conv.res = conv.res, 
                                 max.it = max.it, 
                                 tol = tol, 
                                 full.output = full.output)
  } else {
    R.s.miss_model_smle_chatgpt(sone = sone, 
                                szero = szero, 
                                yone = yone, 
                                yzero = yzero, 
                                conv.res = conv.res, 
                                max.it = max.it, 
                                tol = tol, 
                                full.output = full.output)
  }
}