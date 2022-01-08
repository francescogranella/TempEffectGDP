set scheme Tab10
cd "C:\Users\Granella\Dropbox (CMCC)\PhD\Research\TempEffectGDP"
import delimited "Data\WB_UDel.csv", numericcols(4 5 6 7) clear

drop v1
encode countrycode, gen(id)
xtset id year

rename (udel_pop_preci udel_pop_temp) (prec temp)

drop if temp==.

gen year2 = year*year
replace prec = prec/1000

bys id: gegen tempmean = mean(temp)
bys id: gegen precmean = mean(prec)

reghdfe temp, abs(id id#c.year id#c.year2) resid(tempres)
reghdfe prec, abs(id id#c.year id#c.year2) resid(precres)

lab var temp "Temp"
lab var prec "Prec"


drop if growth==. | temp==. | prec==.
bys id: gen N = _N
keep if N>=30

*	BHM
eststo clear
eststo bhm: reghdfe growth c.temp##c.temp c.prec##c.prec, abs(year id id#c.year id#c.year2) vce(clus id)
estimates save "estimates/bhm.ster", replace
* Auxiliary regression for adjustment factor
qui eststo m0: reghdfe temp c.temp#c.temp c.prec##c.prec, abs(year id id#c.year id#c.year2) vce(clus id) resid(m0)

*	Filter
foreach i in 3 5 10 15 {

	* Filter temp
	bys id: butterworth tempres, freq(`i') stub(x)
	egen temp_`i' = rowtotal(x*)
	replace temp_`i' = tempres - temp_`i'
	replace temp_`i' = tempmean + temp_`i'
	drop x*
	* Filter prec
	bys id: butterworth precres, freq(`i') stub(x)
	egen prec_`i' = rowtotal(x*)  // What?
	replace prec_`i' = precres - prec_`i'
	replace prec_`i' = precmean + prec_`i'
	drop x*

	* Auxiliary regressions for adjustment factor
	// qui reghdfe temp_`i' c.temp_`i'#c.temp_`i' c.prec##c.prec, abs(year id id#c.year id#c.year2) vce(clus id) resid(m`i')
	* Regress
	qui eststo filter_`i': ///
		reghdfe growth ///
			c.temp_`i'##c.temp_`i' ///
			c.prec_`i'##c.prec_`i' ///
			if N>=2*`i', ///
		abs(year id id#c.year id#c.year2) vce(clus id)
	estimates save "estimates/filter_`i'.ster", replace
	// gen x = m0/m`i'
	// sum x, d
	// estadd scalar R = r(p50)
	// drop x
}
*	Sanity check
// preserve
// 	gcollapse (mean) temp temp_3 temp_5 temp_10 temp_15, by(year)
// 	line temp temp_3 temp_5 temp_10 temp_15 year
// restore

estimates use "estimates/bhm.ster"
estimates store bhm
foreach i in 3 5 10 15 {
	estimates use "estimates/filter_`i'.ster"
	estimates store filter_`i'
	loc estlist `estlist' filter_`i'
	loc namelist `namelist' ///
	temp_`i' temp ///
	prec_`i' prec ///
	c.temp_`i'#c.temp_`i' "c.temp#c.temp" ///
	c.prec_`i'#c.prec_`i' "c.prec#c.prec"
}

estfe bhm `estlist', labels(year "Year FE" id "Country FE" id#c.year "Country linear trend" id#c.year2 "Country quadratic trend")
return list
esttab bhm filter_3 filter_5 filter_10 filter_15 using "bhm.txt", label replace nogap compress keep(*temp* *prec*) ///
	indicate(`r(indicate_fe)') varwidth(25) ///
	rename(`namelist') mti("BHM" "3" "5" "10" "15") star(* 0.1 ** 0.05 *** 0.01) se addnote("Standard errors clustered by country") /* stats(R) */

estfe bhm `estlist', labels(year "Year FE" id "Country FE" id#c.year "Country linear trend" id#c.year2 "Country quadratic trend")
return list
esttab bhm filter_3 filter_5 filter_10 filter_15 using "bhm.rtf", label replace nogap compress keep(*temp* *prec*) ///
	indicate(`r(indicate_fe)') varwidth(25) ///
	rename(`namelist') mti("BHM" "3" "5" "10" "15") star(* 0.1 ** 0.05 *** 0.01) se addnote("Standard errors clustered by country") /* stats(R) */

estfe bhm `estlist', labels(year "Year FE" id "Country FE" id#c.year "Country linear trend" id#c.year2 "Country quadratic trend")
return list
esttab bhm filter_3 filter_5 filter_10 filter_15 , label nogap compress keep(*temp* *prec*) ///
	indicate(`r(indicate_fe)') varwidth(25) ///
	rename(`namelist') mti("BHM" "3" "5" "10" "15") star(* 0.1 ** 0.05 *** 0.01) se addnote("Standard errors clustered by country") /* stats(R) */

save "data/bhm", replace

* --------------- *
*	Plot
* --------------- *
estimates restore bhm
margins, dydx(temp) at(temp=(-6(1)30)) vsquish post 
// marginsplot, yline(0) recast(line) ciopts(recast(rarea) lw(none) color(%20))
est store bhm_

loc mlist (bhm_)
loc mlist2 bhm_

foreach i in 3 5 10 15 {
	estimates restore filter_`i'
	// qui margins, at(temp = (-6(1)30)) post
	qui margins, dydx(temp_`i') at(temp_`i'=(-6(1)30)) vsquish post 
	est store filter_`i'_
	loc mlist `mlist' (filter_`i'_)
	loc mlist2 `mlist2'  || filter_`i'_
}

coefplot `mlist', ///
	recast(line) ciopts(recast(rarea) lw(none) color(%20)) at vert xtitle("Temperature") /// leg(lab(2 "BHM") lab(4 "Filter 15")) 
	yline(0)
graph export "panel_quadratic.png", replace

coefplot `mlist2', ///
	recast(line) ciopts(recast(rarea) lw(none) color(%20)) at vert xtitle("Temperature") /// leg(lab(2 "BHM") lab(4 "Filter 15")) 
	yline(0)
graph export "panel_quadratic_by.png", replace

* --------------- *
*	Lags
* --------------- *

global plot

 1 3 5

foreach p of loc periods {
	foreach l of loc lags {
		di "Filter `p',  Lags `l'"
		qui eststo lags_`p'_`l': reghdfe growth ///
					c.temp_`p'##c.temp_`p' ///
					L(0/`l').(c.temp_`p'##c.temp_`p') ///
					c.prec_`p'##c.prec_`p' ///
					L(0/`l').(c.prec_`p'##c.prec_`p') ///
					if N>=2*`p', ///
					abs(year id id#c.year id#c.year2) vce(clus id)
		estimates save "estimates/lags_`p'_`l'.ster", replace			
	}
}

clear
set obs 30
gen t = _n

local periods 0 3 5 10 15
local lags 0 1 3 5
set graphics off

foreach p of loc periods {
	foreach l of loc lags {
		
		estimates use "estimates/lags_`p'_`l'.ster"
		estimates store lags_`p'_`l'
		// estimates replay lags_`p'_`l'
		
		gen Y_`p'_`l'_est = .
		gen Y_`p'_`l'_se = .
		
		forvalues i = 1/30 {

			if `l'==0 {
				qui lincom temp_`p' +  2*`i'*c.temp_`p'#c.temp_`p'
				}
			if `l'==1 {
				qui lincom temp_`p' + L1.temp_`p' +  2*`i'*(c.temp_`p'#c.temp_`p' + cL1.temp_`p'#cL1.temp_`p')
				}
			if `l'==3 {
				qui lincom temp_`p' + L1.temp_`p' + L2.temp_`p' + L3.temp_`p' + 2*`i'*(c.temp_`p'#c.temp_`p' + cL1.temp_`p'#cL1.temp_`p' + cL2.temp_`p'#cL2.temp_`p' + cL3.temp_`p'#cL3.temp_`p')
				}
			if `l'==5 {
				qui lincom temp_`p' + L1.temp_`p' + L2.temp_`p' + L3.temp_`p' + L4.temp_`p' + L5.temp_`p' + ///
					2*`i'*(c.temp_`p'#c.temp_`p' + cL.temp_`p'#cL.temp_`p' + cL2.temp_`p'#cL2.temp_`p' + cL3.temp_`p'#cL3.temp_`p' + cL4.temp_`p'#cL4.temp_`p' + cL5.temp_`p'#cL5.temp_`p')
				}
			qui replace Y_`p'_`l'_est= r(estimate) in `i'
			qui replace Y_`p'_`l'_se = r(se) in `i'
		}
		gen cilo_`p'_`l' = Y_`p'_`l'_est - 1.96*Y_`p'_`l'_se
		gen cihi_`p'_`l' = Y_`p'_`l'_est + 1.96*Y_`p'_`l'_se
		tw (rarea cilo_`p'_`l' cihi_`p'_`l' t, color(%20) blw(none) ti("Filter: `p' Lags: `l'") legend(off) ylabel(,angle(horizontal))) ///
			(line Y_`p'_`l'_est t , lc(black) yline(0, lp(dash))), name(_`p'_`l', replace) // saving($plots/g`j'`p'.gph, replace)) 
		loc glist `glist' _`p'_`l'
	}
}
di "`glist'"
set graphics on
graph combine `glist', rows(5) xcommon ycommon imargin(b=1 t=1)
graph export "img/BHM_lags.png", replace


* --------------- *
*	Compare with Bernie's
* --------------- *
import delimited using "data/Panel_data.csv", clear
drop v1
keep if econdata=="wb"
rename years year
// gen year2 = year^2
// encode countrycode, gen(id)
// reghdfe growth_wdi ///
// 	c.temp##c.temp ///
// 	c.prec##c.prec ///
// 	if filter==15, ///
// 	abs(year id id#c.year id#c.year2) vce(clus id)
keep countrycode year filter growth temp preci
rename (growth temp preci) (growth_b temp_b prec_b)
tempfile bernie
save `bernie'

use "data/bhm", clear
keep if growth!=.
keep countrycode year growth prec temp temp_3 prec_3 temp_5 prec_5 temp_10 prec_10 temp_15 prec_15
rename (temp prec) (temp_0 prec_0)
reshape long temp_ prec_, i(countrycode year) j(f)
rename (temp_ prec_ f) (temp prec filter)
gen year2 = year^2
encode countrycode, gen(id)
bys id filter: gen N = _N
foreach i in 0 3 5 10 15 {
	loc i2 = `i'*2
	qui eststo m_`i': reghdfe growth ///
		c.temp##c.temp ///
		c.prec##c.prec ///
		if filter==`i'  & N>=`i2', ///
		abs(year id id#c.year id#c.year2) vce(clus id)
}

merge 1:1 countrycode year filter using `bernie'
reghdfe growth ///
	c.temp_b##c.temp_b ///
	c.prec_b##c.prec_b ///
	if filter==15, ///
	abs(year id id#c.year id#c.year2) vce(clus id)

replace prec_b = prec_b/1000

reghdfe temp temp_b if filter==15 & N>=35, abs(id)
reghdfe prec prec_b, abs(id)



		reghdfe growth ///
			c.temp_15##c.temp_15 ///
			c.prec_15##c.prec_15 ///
			if N>=2*15, ///
		abs(year id id#c.year id#c.year2) vce(clus id)
