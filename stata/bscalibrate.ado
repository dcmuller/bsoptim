*! bscalibrate
*! v0.1
*! Author:      David C Muller 
*! Email:       <davidmull@gmail.com> 
*! Repository:  <https://github.com/dcmuller/bsoptim>
** 
** Command to calculate the bootstrap optimism corrected c-statistic
** after fitting a logistic regression model, or a survival model

program define bscalibrate, rclass
version 11.0
*syntax, estname(name) classvar(varname) ngroups(integer) [reps(integer 200) predopts(string) at(real) complement]
syntax, estname(name) classvar(varname) ngroups(integer) at(real) [reps(integer 200) predopts(string) complement]
tempvar pred  // variable to hold predictions from the models
tempvar grp // variable to hold quantiles of predicted quantity
tempvar obs // variable to hold 1 - KM estimates
tempvar insample // variable to flag the estimation sample
tempfile km // file to hold the KM estimates
tempname  result predicted observed apparent corrected ///
          boot_in_err boot_out_err o model bsmodel

// initialise vector to hold the sum of optimism
matrix `o' = J(`ngroups', 1, 0)

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


di "{txt} Bootstrap optimism corrected calibration"
di "{txt} Numer of bootstrap replications = `reps'"

// apparent errors from full model
qui predict double `pred' if `insample', `predopts'
if "`complement'" != "" {
  qui replace `pred' = 1 - `pred'
}
_pctile `pred' `wgt', nquantile(`ngroups')
qui gen `grp' = .
local endi = `ngroups' - 1
forval i=1/`endi' {
  local p`i' = r(r`i')
}
forval i=1/`endi' {
  qui recode `grp' .=`i' if `pred' <  `p`i''
}
qui recode `grp' .=`ngroups' if `pred' < .
qui sts list, by(`grp') ///
              at(0 `at') ///
              saving(`km', replace)
preserve
collapse (mean) `pred' `wgt', by(`grp')
qui drop if missing(`grp')
qui merge 1:m `grp' using `km'
assert _merge==3
drop _merge
qui keep if time == `at'
qui gen double `obs'=.
if "`complement'" != "" {
  qui replace `obs' = 1 - survivor
}
else {
  qui replace `obs' = survivor
}
gen double app_err = `pred' - `obs'
mkmat `pred', mat(`predicted')
mkmat `obs', mat(`observed')
mkmat app_err, mat(`apparent')
restore
drop `pred' `grp'

// loop over bootstrap samples
local b=1
while `b' <= `reps' {
  preserve
  bsample `sampsize' if `insample'
  qui `call' 
  est store `bsmodel'
  assert e(N)==`sampsize'
  predict double `pred' if e(sample), `predopts'
  if "`complement'" != "" {
    qui replace `pred' = 1 - `pred'
  }
  qui gen `grp' = .
  forval i=1/`endi' {
    qui recode `grp' .=`i' if `pred' <  `p`i''
  }
  qui recode `grp' .=`ngroups' if `pred' < .
  qui sts list, by(`grp') ///
                at(0 `at') ///
                saving(`km', replace)
  collapse (mean) `pred' `wgt', by(`grp')
  qui drop if missing(`grp')
  qui merge 1:m `grp' using `km'
  assert _merge==3
  drop _merge
  qui keep if time == `at'
  qui gen double `obs'=.
  if "`complement'" != "" {
    qui replace `obs' = 1 - survivor
  }
  else {
    qui replace `obs' = survivor
  }
  gen double app_err = `pred' - `obs'
  mkmat app_err, mat(`boot_in_err')
  restore
  preserve
  qui est restore `bsmodel'
  qui predict double `pred' if `insample', `predopts'
  if "`complement'" != "" {
    qui replace `pred' = 1 - `pred'
  }
  qui gen `grp' = .
  forval i=1/`endi' {
    qui recode `grp' .=`i' if `pred' <  `p`i''
  }
  qui recode `grp' .=`ngroups' if `pred' < .
  qui sts list, by(`grp') ///
                at(0 `at') ///
                saving(`km', replace)
  collapse (mean) `pred' `wgt', by(`grp')
  qui drop if missing(`grp')
  qui merge 1:m `grp' using `km'
  assert _merge==3
  drop _merge
  qui keep if time == `at'
  qui gen double `obs'=.
  if "`complement'" != "" {
    qui replace `obs' = 1 - survivor
  }
  else {
    qui replace `obs' = survivor
  }
  qui gen double app_err = `pred' - `obs'
  mkmat app_err, mat(`boot_out_err')
  restore
  mat `o' = `o' + (`boot_in_err' - `boot_out_err' )
  qui est restore `model'
  local `++b'
}
mat `o' = `o'/`reps'
mat `corrected' = `apparent' + `o'
mat `result' = `predicted', `observed', `apparent', `o', `corrected'
mat colnames `result' = predicted observed apparent optimism corrected
local rn ""
forval i=1/`ngroups' {
  local rn "`rn' `i'"
}
mat rownames `result' = `rn'
matlist `result', rowtitle("Pr. group")

return matrix result=`result'
end


