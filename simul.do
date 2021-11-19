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
foreach i in 3 6 9 12 {
	butterworth T, freq(`i') stub(bw`i')
	replace bw`i'_T_1 = T - bw`i'_T_1
}
line T hp_*_T_sm_1 time, name(line, replace)
graph export "temp.png", replace

capture program drop myreg
program define myreg, rclass
	clear
	set obs 350
	gen time = _n
	gen T = sin(time/5) + cos(time/10) + sin(time/50) + rnormal()
	sum T
	replace T = T-r(mean)
	tsset time
	gen hp_0_T_sm_1 = T
	gen inv_r_0 = T / hp_0_T_sm_1
	foreach i in 10 100 1000 10000 {
		hprescott T, stub(hp_`i') smooth(`i')
		drop hp_`i'_T_1	
		winsor2 hp_`i'_T_sm_1, replace
		gen inv_r_`i' = T / hp_`i'_T_sm_1
	}
	// DGP: g_t = g + \gamma T, g = 0, gamma=1
	gen g1 = T + rnormal()
	// DGP: g_t = g + \beta \Delta T, g = 0, beta=1
	gen dT = T-L1.T
	gen g2 = dT + rnormal()

	foreach i in 0 10 100 1000 10000 {
		sum inv_r_`i', d
		return scalar inv_r_`i' = r(p50)
		reg g1 hp_`i'_T_sm_1
		return scalar b_1_`i' = r(table)[1,1]
		return scalar ll_1_`i' = r(table)[5,1]
		return scalar ul_1_`i' = r(table)[6,1]
		reg g2 hp_`i'_T_sm_1
		return scalar b_2_`i' = r(table)[1,1]
		return scalar ll_2_`i' = r(table)[5,1]
		return scalar ul_2_`i' = r(table)[6,1]
	}
end

simulate ul_2_10000=r(ul_2_10000) ll_2_10000=r(ll_2_10000) b_2_10000=r(b_2_10000) ///
		 ul_1_10000=r(ul_1_10000) ll_1_10000=r(ll_1_10000) b_1_10000=r(b_1_10000) ///
		 ul_2_1000=r(ul_2_1000) ll_2_1000=r(ll_2_1000) b_2_1000=r(b_2_1000) ///
		 ul_1_1000=r(ul_1_1000) ll_1_1000=r(ll_1_1000) b_1_1000=r(b_1_1000) ///
		 ul_2_100=r(ul_2_100) ll_2_100=r(ll_2_100) b_2_100=r(b_2_100) ///
		 ul_1_100=r(ul_1_100) ll_1_100=r(ll_1_100) b_1_100=r(b_1_100) ///
		 ul_2_10=r(ul_2_10) ll_2_10=r(ll_2_10) b_2_10=r(b_2_10) ///
		 ul_1_10=r(ul_1_10) ll_1_10=r(ll_1_10) b_1_10=r(b_1_10) /// 
		 ul_2_0=r(ul_2_0) ll_2_0=r(ll_2_0) b_2_0=r(b_2_0) ///
		 ul_1_0=r(ul_1_0) ll_1_0=r(ll_1_0) b_1_0=r(b_1_0) /// 
		 inv_r_10000=r(inv_r_10000) ///
		 inv_r_1000=r(inv_r_1000) ///
		 inv_r_100=r(inv_r_100) ///
		 inv_r_10=r(inv_r_10) ///
		 inv_r_0=r(inv_r_0) ///
	, reps(100) seed(1234): myreg

drop ul* ll*

foreach i in 0 10 100 1000 10000 {
	tw (kdensity b_1_`i') (kdensity b_2_`i'), xline(1) name(b_`i', replace) ti(`i') leg(lab(1 "{&gamma} model") lab(2 "{&beta} model")) xti({&gamma})
}
grc1leg b_0 b_10 b_100 b_1000 b_10000, imargin(zero) name(b, replace) xcommon ycommon
graph export "repli.png", replace

foreach i in 0 10 100 1000 10000 {
	cap gen adj_b_1_`i' = b_1_`i' / inv_r_`i'
	tw (kdensity b_1_`i') (kdensity adj_b_1_`i'), xline(1) name(k_`i', replace) ti(`i') leg(lab(1 "{&gamma} model unadjusted") lab(2 "{&gamma} model adjusted"))  xti({&gamma})
}
grc1leg k_0 k_10 k_100 k_1000 k_10000, imargin(zero) name(adj, replace) xcommon ycommon
graph export "repli_adj.png", replace
