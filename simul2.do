set scheme Tab10
cd "C:\Users\Granella\Desktop\repli"

clear
set obs 350
gen time = _n
gen T = sin(time/5) + cos(time/10) + sin(time/50) + rnormal()
tsset time
foreach i in 10 100 1000 10000 {
	hprescott T, stub(hp_`i') smooth(`i')
	gen inv_r_`i' = T / hp_`i'_T_sm_1
}
line T hp*sm* time, name(bw, replace)

graph export "temp.png", replace

set seed 1234
capture program drop myreg
program define myreg, rclass
	clear
	set obs 300
	gen time = _n
	tsset time
	gen T = 20 + 0.01*time + sin(time/5) + cos(time/10) - sin(time/50) + rnormal()
	// gen T = 20 + 0.01*time + sin(time/(rnormal(1,1)*5)) + cos(time/(rnormal(1,1)*10)) - sin(time/(rnormal(1,1)*50)) + rnormal()
	* DGP: g_t = g + \gamma T, g = 0, gamma=1
	gen g1 = T - 0.5*L1.T - 0.25*L2.T - 0.25*L3.T + rnormal()
	* DGP: g_t = g + \beta \Delta T, g = 0, beta=1
	gen dT = T-L1.T
	gen g2 = dT + rnormal()
	* Detrend and center
	reg T time
	predict res, res
	drop T
	rename res T
	* 
	tsset time
	* Baseline
	gen hp_0_T_sm_1 = T
	gen inv_r_0 = T / hp_0_T_sm_1
	foreach i in 10 25 50 {
		// Low-pass filter
		butterworth T, freq(`i') stub(bw_`i')
		replace bw_`i'_T_1 = T - bw_`i'_T_1
		gen inv_r_`i' = T / bw_`i'_T_1
		// winsor2 inv_r_`i', replace
	}

	foreach i in 10 25 50 {
		// reg T g1
		// dfbeta, stub(dfbeta_`i')
		// reg inv_r_`i'[w=dfbeta_i] 
		// gen w_`i' = abs(dfbeta_`i'1)
		sum inv_r_`i', d // [iw=w_`i']
		return scalar inv_r_`i' =  r(p50)
		reg g1 bw_`i'_T_1, nohead
		return scalar b_1_`i' = r(table)[1,1]
		// return scalar ll_1_`i' = r(table)[5,1]
		// return scalar ul_1_`i' = r(table)[6,1]
		reg g2 bw_`i'_T_1, nohead
		return scalar b_2_`i' = r(table)[1,1]
		// return scalar ll_2_`i' = r(table)[5,1]
		// return scalar ul_2_`i' = r(table)[6,1]
	}
end

simulate b_1_50=r(b_1_50) b_2_50=r(b_2_50) inv_r_50=r(inv_r_50) /// 
		b_1_25=r(b_1_25) b_2_25=r(b_2_25) inv_r_25=r(inv_r_25) /// 
		b_1_10=r(b_1_10) b_2_10=r(b_2_10) inv_r_10=r(inv_r_10) /// 
		, reps(100) seed(1234): myreg

// simulate ul_2_10000=r(ul_2_10000) ll_2_10000=r(ll_2_10000) b_2_10000=r(b_2_10000) ///
// 		 ul_1_10000=r(ul_1_10000) ll_1_10000=r(ll_1_10000) b_1_10000=r(b_1_10000) ///
// 		 ul_2_1000=r(ul_2_1000) ll_2_1000=r(ll_2_1000) b_2_1000=r(b_2_1000) ///
// 		 ul_1_1000=r(ul_1_1000) ll_1_1000=r(ll_1_1000) b_1_1000=r(b_1_1000) ///
// 		 ul_2_100=r(ul_2_100) ll_2_100=r(ll_2_100) b_2_100=r(b_2_100) ///
// 		 ul_1_100=r(ul_1_100) ll_1_100=r(ll_1_100) b_1_100=r(b_1_100) ///
// 		 ul_2_10=r(ul_2_10) ll_2_10=r(ll_2_10) b_2_10=r(b_2_10) ///
// 		 ul_1_10=r(ul_1_10) ll_1_10=r(ll_1_10) b_1_10=r(b_1_10) /// 
// 		 ul_2_0=r(ul_2_0) ll_2_0=r(ll_2_0) b_2_0=r(b_2_0) ///
// 		 ul_1_0=r(ul_1_0) ll_1_0=r(ll_1_0) b_1_0=r(b_1_0) /// 
// 		 inv_r_10000=r(inv_r_10000) ///
// 		 inv_r_1000=r(inv_r_1000) ///
// 		 inv_r_100=r(inv_r_100) ///
// 		 inv_r_10=r(inv_r_10) ///
// 		 inv_r_0=r(inv_r_0) ///
// 	, reps(50) seed(1234): myreg

// drop ul* ll*

// foreach i in 10 25 50 {
// 	tw (kdensity b_1_`i') (kdensity b_2_`i'), xline(1) name(b_`i', replace) ti(`i') leg(lab(1 "{&gamma} model") lab(2 "{&beta} model")) xti({&gamma})
// }
// grc1leg b_0 b_10 b_100 b_1000 b_10000, imargin(zero) name(b, replace) xcommon ycommon
// graph export "repli.png", replace
set graphics off
foreach i in 10 25 50 {
	cap gen adj_b_1_`i' = b_1_`i' / inv_r_`i'
	tw (kdensity b_1_`i') (kdensity adj_b_1_`i'), xline(1) name(k_`i', replace) ti(`i') leg(lab(1 "{&gamma} model unadjusted") lab(2 "{&gamma} model adjusted"))  xti({&gamma})
}
set graphics on
grc1leg k_10 k_25 k_50, imargin(zero) name(adj, replace) xcommon ycommon
graph export "repli_adj.png", replace
