*! bsoptim
*! v0.1
*! Author:      David C Muller 
*! Email:       <davidmull@gmail.com> 
*! Repository:  <https://github.com/dcmuller/bsoptim>
** 
** Command to calculate the bootstrap optimism corrected c-statistic
** after fitting a logistic regression model

program define bsoptim, rclass
version 11.0
syntax [, reps(integer 200) estname(name)]
tempvar pred  // variable to hold predictions from the models
tempvar insample // variable to flag the estimation sample
tempname o corrected model

// initialise scalar to hold the sum of optimism
scalar `o' = 0

// store estimates frin the model
if "`estname'" != "" {
  qui est restore `estname'
}
est store `model'
gen `insample' = e(sample)
local sampsize = e(N)
local call = e(cmdline)
local yvar = e(depvar)

di "{txt} Bootstrap optimism corrected c-statistic"
di "{txt} Numer of bootstrap replications = `reps'"

// apparent c from full model
qui lroc, nograph
local c_app = r(area)
*di `c_app'

// loop over bootstrap samples
local i=1
while `i' <= `reps' {
  preserve
  bsample `sampsize' if `insample'
  qui `call' 
  assert e(N)==`sampsize'
  qui lroc, nograph
  local c_boot = r(area)
  restore
  predict `pred' if `insample', pr
  qui roctab `yvar' `pred'
  local c_orig = r(area)
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


