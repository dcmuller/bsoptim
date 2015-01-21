*! bsoptim
*! v0.3
*! Author:      David C Muller 
*! Email:       <davidmull@gmail.com> 
*! Repository:  <https://github.com/dcmuller/bsoptim>
** 
** Command to calculate the bootstrap optimism corrected c-statistic
** after fitting a logistic regression model, or a survival model

program define bsoptim, rclass
version 11.0
syntax, estname(name) classvar(varname) [reps(integer 200) cluster(varname) predopts(string) complement]
tempvar pred  // variable to hold predictions from the models
tempvar insample // variable to flag the estimation sample
tempvar idx // variable to flag the number of clusters (if cluster specified)
tempname c o corrected model bsmodel clstid


// initialise scalar to hold the sum of optimism
scalar `o' = 0

// store estimates frin the model
if "`estname'" != "" {
  qui est restore `estname'
}
est store `model'
gen `insample' = e(sample)
local sampsize = e(N)
local call `e(cmdline)'
local yvar "`classvar'"
local wgt ""
if "`e(wtype)'" != "" {
  local wgt "[`e(wtype)' `e(wexp)']"
}

if "`cluster'" != "" {
  local clstopt "cluster(`cluster') idcluster(`clstid')"
  bys `cluster': gen `idx'=_n==1
  qui count if `idx' & `insample'
  local nclust = r(N)
}


di "{txt} Bootstrap optimism corrected c-statistic"
di "{txt} Numer of bootstrap replications = `reps'"

// apparent c from full model
qui predict double `pred' if `insample', `predopts'
if "`complement'" != "" {
  qui replace `pred' = 1 - `pred'
}
qui somersd `yvar' `pred' if `insample' `wgt', transf(c)
mat `c' = e(b)
local c_app = `c'[1,1]
drop `pred'
*di `c_app'

// loop over bootstrap samples
local i=1
while `i' <= `reps' {
  preserve
  if "`cluster'" !="" {
    local ndraw `nclust'
  }
  else {
    local ndraw `sampsize'
  }
  bsample `ndraw' if `insample', `clstopt'
  if "`clstopt'" != "" {
    qui replace `cluster'=`clstid'
  }

  qui `call' 
  est store `bsmodel'
  *assert e(N)==`sampsize'
  predict double `pred' if e(sample), `predopts'
  if "`complement'" != "" {
    qui replace `pred' = 1 - `pred'
  }
  qui somersd `yvar' `pred' if e(sample) `wgt', transf(c)
  mat `c' = e(b)
  local c_boot = `c'[1,1]
  restore
  qui est restore `bsmodel'
  qui predict double `pred' if `insample', `predopts'
  if "`complement'" != "" {
    qui replace `pred' = 1 - `pred'
  }
  qui somersd `yvar' `pred' `wgt', transf(c)
  mat `c' = e(b)
  local c_orig = `c'[1,1]
  scalar `o' = `o' + (`c_boot' - `c_orig')
  *di scalar(`o')
  drop `pred'
  qui est restore `model'
  local `++i'
}
scalar `o' = `o'/`reps'
scalar `corrected' = `c_app' - `o'
di "{txt}{hline 45}"
di "{ralign 15: corrected}{ralign 15: optimism}{ralign 15: apparent}"
di "{txt}{hline 45}"
di  "{ralign 10:}{res}" %04.3f `corrected' ///
    "{ralign 10:}{res}" %04.3f `o' ///  
    "{ralign 10:}{res}" %04.3f `c_app'  
di "{txt}{hline 45}"

return scalar c_apparent  = `c_app'
return scalar optimism    = `o'
return scalar c_corrected = `corrected'
return local  N_reps      = `reps'
end


