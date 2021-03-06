	# Bias-independent calculations
	if (TYPE == ntype) 
		devsign = 1
	else 
		devsign = -1
	end
	
	# Constants
	epssi    = EPSRSUB * EPS0
	epsox    = EPSROX * EPS0
	Cox      = EPSROX * EPS0 / TOXE
	epsratio = EPSRSUB / EPSROX
	
	# Physical oxide thickness
	if (!$param_given(TOXP)) 
		BSIMBULKTOXP = (TOXE * EPSROX / 3.9) - DTOX
	else 
		BSIMBULKTOXP = TOXP
	end
	L_mult = L * LMLT
	W_mult = W * WMLT
	Lnew = L_mult + XL
	if (Lnew <= 0.0) 
		println("Fatal: Ldrawn * LMLT + XL = %e for %M is non-positive", Lnew)
		
	end
	W_by_NF = W_mult / NF
	Wnew    = W_by_NF + XW
	if (Wnew <= 0.0) 
		println("Fatal: W / NF * WMLT + XW = %e for %M is non-positive", Wnew)
		
	end
	
	# Leff and Weff for I-V
	L_LLN      = pow(Lnew, -LLN)
	W_LWN      = pow(Wnew, -LWN)
	LW_LLN_LWN = L_LLN * W_LWN
	dLIV       = LINT + LL * L_LLN + LW * W_LWN + LWL * LW_LLN_LWN
	L_WLN      = pow(Lnew, -WLN)
	W_WWN      = pow(Wnew, -WWN)
	LW_WLN_WWN = L_WLN * W_WWN
	dWIV       = WINT + WL * L_WLN + WW * W_WWN + WWL * LW_WLN_WWN
	Leff       = Lnew - 2.0 * dLIV
	if (Leff <= 0.0) 
		println("Fatal: Effective channel length = %e for  %M is non-positive", Leff)
		
	elseif (Leff <= 1.0e-9) 
		println("Warning: Effective channel length = %e for %M is <= 1.0e-9. Recommended Leff >= 1e-8", Leff)
	end
	Weff = Wnew - 2.0 * dWIV
	if (Weff <= 0.0) 
		 println("Fatal: Effective channel Width = %e for %M is non-positive", Weff)
		
	elseif (Weff <= 1.0e-9) 
		println("Warning: Effective channel width = %e for %M is <= 1.0e-9. Recommended Weff >= 1e-8", Weff)
	end
	
	# Leff and Weff for C-V
	dLCV = DLC + LLC * L_LLN + LWC * W_LWN + LWLC * LW_LLN_LWN
	dWCV = DWC + WLC * L_WLN + WWC * W_WWN + WWLC * LW_WLN_WWN
	Lact = Lnew - 2.0 * dLCV
	if (Lact <= 0.0) 
		 println("Fatal: Effective channel length for C-V = %e for %M is non-positive", Lact)
		
	elseif (Lact <= 1.0e-9) 
		println("Warning: Effective channel length for C-V = %e for %M is <= 1.0e-9. Recommended Lact >= 1e-8", Lact)
	end
	Wact = Wnew - 2.0 * dWCV
	if (Wact <= 0.0) 
		println("Fatal: Effective channel width for C-V = %e for %M is non-positive", Wact)
		
	elseif (Wact <= 1.0e-9) 
		println("Warning: Effective channel width for C-V = %e for %M is <= 1.0e-9. Recommended Wact >= 1e-8", Wact)
	end
	
	# Weffcj for diode, GIDL etc.
	dWJ    = DWJ + WLC / pow(Lnew, WLN) + WWC / pow(Wnew, WWN) + WWLC / pow(Lnew, WLN) / pow(Wnew, WWN)
	Weffcj = Wnew - 2.0 * dWJ
	if (Weffcj <= 0.0) 
		println("Fatal: Effective channel width for S/D junctions = %e for %M is non-positive", Weffcj)
		
	 end
	Inv_L     = 1.0e-6 / Leff
	Inv_W     = 1.0e-6 / Weff
	Inv_Lact  = 1.0e-6 / Lact
	Inv_Wact  = 1.0e-6 / Wact
	Inv_Llong = 1.0e-6 / LLONG
	Inv_Wwide = 1.0e-6 / WWIDE
	Inv_WL    = Inv_L * Inv_W
	
	# Effective length and width for binning
	L_LLN1 = L_LLN
	L_WLN1 = L_WLN
	if (DLBIN != 0.0) 
		if (DLBIN <= -Lnew) 
			println("Fatal: DLBIN for %M = %e is <= -Ldrawn * LMLT", DLBIN)
			
		else 
			L_LLN1 = pow(Lnew + DLBIN, -LLN)
			L_WLN1 = pow(Lnew + DLBIN, -WLN)
		end
	end
	W_LWN1 = W_LWN
	W_WWN1 = W_WWN
	if (DWBIN != 0.0) 
		if (DWBIN <= -Wnew) 
			println("Fatal: DWBIN for %M = %e is <= -Wdrawn * WMLT", DWBIN)
			
		else 
			W_LWN1 = pow(Wnew + DWBIN, -LWN)
			W_WWN1 = pow(Wnew + DWBIN, -WWN)
		end
	end
	LW_LLN_LWN1 = L_LLN1 * W_LWN1
	dLB         = LINT + LL * L_LLN1 + LW * W_LWN1 + LWL * LW_LLN_LWN1
	LW_WLN_WWN1 = L_WLN1 * W_WWN1
	dWB         = WINT + WL * L_WLN1 + WW * W_WWN1 + WWL * LW_WLN_WWN1
	Leff1 = Lnew - 2.0 * dLB + DLBIN
	if (Leff1 <= 0.0) 
		println("Fatal: Effective channel length for binning = %e for %M is non-positive", Leff1)
		
	end
	Weff1 = Wnew - 2.0 * dWB + DWBIN
	if (Weff1 <= 0.0) 
		println("Fatal: Effective channel width for binning = %e for %M is non-positive", Weff1)
		
	end
	if (BINUNIT == 1) 
		BIN_L = 1.0e-6 / Leff1
		BIN_W = 1.0e-6 / Weff1
	else 
		BIN_L = 1.0 / Leff1
		BIN_W = 1.0 / Weff1
	end
	BIN_WL         = BIN_L * BIN_W
	VFB_i          = VFB + BIN_L * LVFB + BIN_W * WVFB + BIN_WL * PVFB
	VFBCV_i        = VFBCV + BIN_L * LVFBCV + BIN_W * WVFBCV + BIN_WL * PVFBCV
	NSD_i          = NSD + BIN_L * LNSD + BIN_W * WNSD + BIN_WL * PNSD
	NDEP_i         = NDEP + BIN_L * LNDEP + BIN_W * WNDEP + BIN_WL * PNDEP
	NDEPCV_i       = NDEPCV + BIN_L * LNDEPCV + BIN_W * WNDEPCV + BIN_WL * PNDEPCV
	NGATE_i        = NGATE + BIN_L * LNGATE + BIN_W * WNGATE + BIN_WL * PNGATE
	CIT_i          = CIT + BIN_L * LCIT + BIN_W * WCIT + BIN_WL * PCIT
	NFACTOR_i      = NFACTOR + BIN_L * LNFACTOR + BIN_W * WNFACTOR + BIN_WL * PNFACTOR
	CDSCD_i        = CDSCD + BIN_L * LCDSCD + BIN_W * WCDSCD + BIN_WL * PCDSCD
	CDSCB_i        = CDSCB + BIN_L * LCDSCB + BIN_W * WCDSCB + BIN_WL * PCDSCB
	DVTP0_i        = DVTP0 + BIN_L * LDVTP0 + BIN_W * WDVTP0 + BIN_WL * PDVTP0
	DVTP1_i        = DVTP1 + BIN_L * LDVTP1 + BIN_W * WDVTP1 + BIN_WL * PDVTP1
	DVTP2_i        = DVTP2 + BIN_L * LDVTP2 + BIN_W * WDVTP2 + BIN_WL * PDVTP2
	DVTP3_i        = DVTP3 + BIN_L * LDVTP3 + BIN_W * WDVTP3 + BIN_WL * PDVTP3
	DVTP4_i        = DVTP4 + BIN_L * LDVTP4 + BIN_W * WDVTP4 + BIN_WL * PDVTP4
	DVTP5_i        = DVTP5 + BIN_L * LDVTP5 + BIN_W * WDVTP5 + BIN_WL * PDVTP5
	K2_i           = K2 + BIN_L * LK2 + BIN_W * WK2 + BIN_WL * PK2
	K1_i           = K1 + BIN_L * LK1 + BIN_W * WK1 + BIN_WL * PK1
	XJ_i           = XJ + BIN_L * LXJ + BIN_W * WXJ + BIN_WL * PXJ
	PHIN_i         = PHIN + BIN_L * LPHIN + BIN_W * WPHIN + BIN_WL * PPHIN
	ETA0_i         = ETA0 + BIN_L * LETA0 + BIN_W * WETA0 + BIN_WL * PETA0
	ETAB_i         = ETAB + BIN_L * LETAB + BIN_W * WETAB + BIN_WL * PETAB
	DELTA_i        = DELTA + BIN_L * LDELTA + BIN_W * WDELTA + BIN_WL * PDELTA
	U0_i           = U0 + BIN_L * LU0 + BIN_W * WU0 + BIN_WL * PU0
	UA_i           = UA + BIN_L * LUA + BIN_W * WUA + BIN_WL * PUA
	UD_i           = UD + BIN_L * LUD + BIN_W * WUD + BIN_WL * PUD
	EU_i           = EU + BIN_L * LEU + BIN_W * WEU + BIN_WL * PEU
	UCS_i          = UCS + BIN_L * LUCS + BIN_W * WUCS + BIN_WL * PUCS
	UC_i           = UC + BIN_L * LUC + BIN_W * WUC + BIN_WL * PUC
	PCLM_i         = PCLM + BIN_L * LPCLM + BIN_W * WPCLM + BIN_WL * PPCLM
	PCLMCV_i       = PCLMCV + BIN_L * LPCLMCV + BIN_W * WPCLMCV + BIN_WL * PPCLMCV
	RSW_i          = RSW + BIN_L * LRSW + BIN_W * WRSW + BIN_WL * PRSW
	RDW_i          = RDW + BIN_L * LRDW + BIN_W * WRDW + BIN_WL * PRDW
	PRWG_i         = PRWG + BIN_L * LPRWG + BIN_W * WPRWG + BIN_WL * PPRWG
	PRWB_i         = PRWB + BIN_L * LPRWB + BIN_W * WPRWB + BIN_WL * PPRWB
	WR_i           = WR + BIN_L * LWR + BIN_W * WWR + BIN_WL * PWR
	RSWMIN_i       = RSWMIN + BIN_L * LRSWMIN + BIN_W * WRSWMIN + BIN_WL * PRSWMIN
	RDWMIN_i       = RDWMIN + BIN_L * LRDWMIN + BIN_W * WRDWMIN + BIN_WL * PRDWMIN
	RDSW_i         = RDSW + BIN_L * LRDSW + BIN_W * WRDSW + BIN_WL * PRDSW
	RDSWMIN_i      = RDSWMIN + BIN_L * LRDSWMIN + BIN_W * WRDSWMIN + BIN_WL * PRDSWMIN
	PTWG_i         = PTWG + BIN_L * LPTWG + BIN_W * WPTWG + BIN_WL * PPTWG
	PDIBLC_i       = PDIBLC + BIN_L * LPDIBLC + BIN_W * WPDIBLC + BIN_WL * PPDIBLC
	PDIBLCB_i      = PDIBLCB + BIN_L * LPDIBLCB + BIN_W * WPDIBLCB + BIN_WL * PPDIBLCB
	PSCBE1_i       = PSCBE1 + BIN_L * LPSCBE1 + BIN_W * WPSCBE1 + BIN_WL * PPSCBE1
	PSCBE2_i       = PSCBE2 + BIN_L * LPSCBE2 + BIN_W * WPSCBE2 + BIN_WL * PPSCBE2
	PDITS_i        = PDITS + BIN_L * LPDITS + BIN_W * WPDITS + BIN_WL * PPDITS
	PDITSD_i       = PDITSD + BIN_L * LPDITSD + BIN_W * WPDITSD + BIN_WL * PPDITSD
	FPROUT_i       = FPROUT + BIN_L * LFPROUT + BIN_W * WFPROUT + BIN_WL * PFPROUT
	PVAG_i         = PVAG + BIN_L * LPVAG + BIN_W * WPVAG + BIN_WL * PPVAG
	VSAT_i         = VSAT + BIN_L * LVSAT + BIN_W * WVSAT + BIN_WL * PVSAT
	PSAT_i         = PSAT + BIN_L * LPSAT + BIN_W * WPSAT + BIN_WL * PPSAT
	VSATCV_i       = VSATCV + BIN_L * LVSATCV + BIN_W * WVSATCV + BIN_WL * PVSATCV
	CF_i           = CF + BIN_L * LCF + BIN_W * WCF + BIN_WL * PCF
	CGSL_i         = CGSL + BIN_L * LCGSL + BIN_W * WCGSL + BIN_WL * PCGSL
	CGDL_i         = CGDL + BIN_L * LCGDL + BIN_W * WCGDL + BIN_WL * PCGDL
	CKAPPAS_i      = CKAPPAS + BIN_L * LCKAPPAS + BIN_W * WCKAPPAS + BIN_WL * PCKAPPAS
	CKAPPAD_i      = CKAPPAD + BIN_L * LCKAPPAD + BIN_W * WCKAPPAD + BIN_WL * PCKAPPAD
	ALPHA0_i       = ALPHA0 + BIN_L * LALPHA0 + BIN_W * WALPHA0 + BIN_WL * PALPHA0
	BETA0_i        = BETA0 + BIN_L * LBETA0 + BIN_W * WBETA0 + BIN_WL * PBETA0
	KVTH0WE_i      = KVTH0WE + BIN_L * LKVTH0WE  + BIN_W * WKVTH0WE + BIN_WL * PKVTH0WE
	K2WE_i         = K2WE + BIN_L * LK2WE + BIN_W * WK2WE + BIN_WL * PK2WE
	KU0WE_i        = KU0WE + BIN_L * LKU0WE + BIN_W * WKU0WE + BIN_WL * PKU0WE
	AGIDL_i        = AGIDL + BIN_L * LAGIDL + BIN_W * WAGIDL + BIN_WL * PAGIDL
	BGIDL_i        = BGIDL + BIN_L * LBGIDL + BIN_W * WBGIDL + BIN_WL * PBGIDL
	CGIDL_i        = CGIDL + BIN_L * LCGIDL + BIN_W * WCGIDL + BIN_WL * PCGIDL
	EGIDL_i        = EGIDL + BIN_L * LEGIDL + BIN_W * WEGIDL + BIN_WL * PEGIDL
	AGISL_i        = AGISL + BIN_L * LAGISL + BIN_W * WAGISL + BIN_WL * PAGISL
	BGISL_i        = BGISL + BIN_L * LBGISL + BIN_W * WBGISL + BIN_WL * PBGISL
	CGISL_i        = CGISL + BIN_L * LCGISL + BIN_W * WCGISL + BIN_WL * PCGISL
	EGISL_i        = EGISL + BIN_L * LEGISL + BIN_W * WEGISL + BIN_WL * PEGISL
	UTE_i          = UTE + BIN_L * LUTE + BIN_W * WUTE + BIN_WL * PUTE
	UA1_i          = UA1 + BIN_L * LUA1 + BIN_W * WUA1 + BIN_WL * PUA1
	UC1_i          = UC1 + BIN_L * LUC1 + BIN_W * WUC1 + BIN_WL * PUC1
	UD1_i          = UD1 + BIN_L * LUD1 + BIN_W * WUD1 + BIN_WL * PUD1
	EU1_i          = EU1 + BIN_L * LEU1 + BIN_W * WEU1 + BIN_WL * PEU1
	UCSTE_i        = UCSTE + BIN_L * LUCSTE + BIN_W * WUCSTE + BIN_WL * PUCSTE
	PRT_i          = PRT + BIN_L * LPRT + BIN_W * WPRT + BIN_WL * PPRT
	AT_i           = AT + BIN_L * LAT + BIN_W * WAT + BIN_WL * PAT
	PTWGT_i        = PTWGT + BIN_L * LPTWGT + BIN_W * WPTWGT + BIN_WL * PPTWGT
	IIT_i          = IIT + BIN_L * LIIT + BIN_W * WIIT + BIN_WL * PIIT
	TGIDL_i        = TGIDL + BIN_L * LTGIDL + BIN_W * WTGIDL + BIN_WL * PTGIDL
	IGT_i          = IGT + BIN_L * LIGT + BIN_W * WIGT + BIN_WL * PIGT
	AIGBINV_i      = AIGBINV + BIN_L * LAIGBINV + BIN_W * WAIGBINV + BIN_WL * PAIGBINV
	BIGBINV_i      = BIGBINV + BIN_L * LBIGBINV + BIN_W * WBIGBINV + BIN_WL * PBIGBINV
	CIGBINV_i      = CIGBINV + BIN_L * LCIGBINV + BIN_W * WCIGBINV + BIN_WL * PCIGBINV
	EIGBINV_i      = EIGBINV + BIN_L * LEIGBINV + BIN_W * WEIGBINV + BIN_WL * PEIGBINV
	NIGBINV_i      = NIGBINV + BIN_L * LNIGBINV + BIN_W * WNIGBINV + BIN_WL * PNIGBINV
	AIGBACC_i      = AIGBACC + BIN_L * LAIGBACC + BIN_W * WAIGBACC + BIN_WL * PAIGBACC
	BIGBACC_i      = BIGBACC + BIN_L * LBIGBACC + BIN_W * WBIGBACC + BIN_WL * PBIGBACC
	CIGBACC_i      = CIGBACC + BIN_L * LCIGBACC + BIN_W * WCIGBACC + BIN_WL * PCIGBACC
	NIGBACC_i      = NIGBACC + BIN_L * LNIGBACC + BIN_W * WNIGBACC + BIN_WL * PNIGBACC
	AIGC_i         = AIGC + BIN_L * LAIGC + BIN_W * WAIGC + BIN_WL * PAIGC
	BIGC_i         = BIGC + BIN_L * LBIGC + BIN_W * WBIGC + BIN_WL * PBIGC
	CIGC_i         = CIGC + BIN_L * LCIGC + BIN_W * WCIGC + BIN_WL * PCIGC
	AIGS_i         = AIGS + BIN_L * LAIGS + BIN_W * WAIGS + BIN_WL * PAIGS
	BIGS_i         = BIGS + BIN_L * LBIGS + BIN_W * WBIGS + BIN_WL * PBIGS
	CIGS_i         = CIGS + BIN_L * LCIGS + BIN_W * WCIGS + BIN_WL * PCIGS
	AIGD_i         = AIGD + BIN_L * LAIGD + BIN_W * WAIGD + BIN_WL * PAIGD
	BIGD_i         = BIGD + BIN_L * LBIGD + BIN_W * WBIGD + BIN_WL * PBIGD
	CIGD_i         = CIGD + BIN_L * LCIGD + BIN_W * WCIGD + BIN_WL * PCIGD
	POXEDGE_i      = POXEDGE + BIN_L * LPOXEDGE + BIN_W * WPOXEDGE + BIN_WL * PPOXEDGE
	DLCIG_i        = DLCIG + BIN_L * LDLCIG + BIN_W * WDLCIG + BIN_WL * PDLCIG
	DLCIGD_i       = DLCIGD + BIN_L * LDLCIGD + BIN_W * WDLCIGD + BIN_WL * PDLCIGD
	NTOX_i         = NTOX + BIN_L * LNTOX + BIN_W * WNTOX + BIN_WL * PNTOX
	KT1_i          = KT1 + BIN_L * LKT1 + BIN_W * WKT1 + BIN_WL * PKT1
	KT2_i          = KT2 + BIN_L * LKT2 + BIN_W * WKT2 + BIN_WL * PKT2
	PSATB_i        = PSATB + BIN_L * LPSATB + BIN_W * WPSATB + BIN_WL * PPSATB
	A1_i           = A1 + BIN_L * LA1 + BIN_W * WA1 + BIN_WL * PA1
	A11_i          = A11 + BIN_L * LA11 + BIN_W * WA11 + BIN_WL * PA11
	A2_i           = A2 + BIN_L * LA2 + BIN_W * WA2 + BIN_WL * PA2
	A21_i          = A21 + BIN_L * LA21 + BIN_W * WA21 + BIN_WL * PA21
	K0_i           = K0 + BIN_L * LK0 + BIN_W * WK0 + BIN_WL * PK0
	M0_i           = M0 + BIN_L * LM0 + BIN_W * WM0 + BIN_WL * PM0
	K01_i          = K01 + BIN_L * LK01 + BIN_W * WK01 + BIN_WL * PK01
	M01_i          = M01 + BIN_L * LM01 + BIN_W * WM01 + BIN_WL * PM01
	NFACTOREDGE_i  = NFACTOREDGE + BIN_L * LNFACTOREDGE + BIN_W * WNFACTOREDGE + BIN_WL * PNFACTOREDGE
	NDEPEDGE_i     = NDEPEDGE + BIN_L * LNDEPEDGE + BIN_W * WNDEPEDGE + BIN_WL * PNDEPEDGE
	CITEDGE_i      = CITEDGE + BIN_L * LCITEDGE + BIN_W * WCITEDGE + BIN_WL * PCITEDGE
	CDSCDEDGE_i    = CDSCDEDGE + BIN_L * LCDSCDEDGE + BIN_W * WCDSCDEDGE + BIN_WL * PCDSCDEDGE
	CDSCBEDGE_i    = CDSCBEDGE + BIN_L * LCDSCBEDGE + BIN_W * WCDSCBEDGE + BIN_WL * PCDSCBEDGE
	ETA0EDGE_i     = ETA0EDGE + BIN_L * LETA0EDGE + BIN_W * WETA0EDGE + BIN_WL * PETA0EDGE
	ETABEDGE_i     = ETABEDGE + BIN_L * LETABEDGE + BIN_W * WETABEDGE + BIN_WL * PETABEDGE
	KT1EDGE_i      = KT1EDGE + BIN_L * LKT1EDGE + BIN_W * WKT1EDGE + BIN_WL * PKT1EDGE
	KT1LEDGE_i     = KT1LEDGE + BIN_L * LKT1LEDGE + BIN_W * WKT1LEDGE + BIN_WL * PKT1LEDGE
	KT2EDGE_i      = KT2EDGE + BIN_L * LKT2EDGE + BIN_W * WKT2EDGE + BIN_WL * PKT2EDGE
	KT1EXPEDGE_i   = KT1EXPEDGE + BIN_L * LKT1EXPEDGE + BIN_W * WKT1EXPEDGE + BIN_WL * PKT1EXPEDGE
	TNFACTOREDGE_i = TNFACTOREDGE + BIN_L * LTNFACTOREDGE + BIN_W * WTNFACTOREDGE + BIN_WL * PTNFACTOREDGE
	TETA0EDGE_i    = TETA0EDGE + BIN_L * LTETA0EDGE + BIN_W * WTETA0EDGE + BIN_WL * PTETA0EDGE
	K2EDGE_i       = K2EDGE + BIN_L * LK2EDGE + BIN_W * WK2EDGE + BIN_WL * PK2EDGE
	KVTH0EDGE_i    = KVTH0EDGE + BIN_L * LKVTH0EDGE + BIN_W * WKVTH0EDGE + BIN_WL * PKVTH0EDGE
	STK2EDGE_i     = STK2EDGE + BIN_L * LSTK2EDGE + BIN_W * WSTK2EDGE + BIN_WL * PSTK2EDGE
	STETA0EDGE_i   = STETA0EDGE + BIN_L * LSTETA0EDGE + BIN_W * WSTETA0EDGE + BIN_WL * PSTETA0EDGE
	C0_i           = C0 + BIN_L * LC0 + BIN_W * WC0 + BIN_WL * PC0
	C01_i          = C01 + BIN_L * LC01 + BIN_W * WC01 + BIN_WL * PC01
	C0SI_i         = C0SI + BIN_L * LC0SI + BIN_W * WC0SI + BIN_WL * PC0SI
	C0SI1_i        = C0SI1 + BIN_L * LC0SI1 + BIN_W * WC0SI1 + BIN_WL * PC0SI1
	C0SISAT_i      = C0SISAT + BIN_L * LC0SISAT + BIN_W * WC0SISAT + BIN_WL * PC0SISAT
	C0SISAT1_i     = C0SISAT1 + BIN_L * LC0SISAT1 + BIN_W * WC0SISAT1 + BIN_WL * PC0SISAT1
	
	if (ASYMMOD != 0) 
		CDSCDR_i  = CDSCDR + BIN_L * LCDSCDR + BIN_W * WCDSCDR + BIN_WL * PCDSCDR
		ETA0R_i   = ETA0R + BIN_L * LETA0R + BIN_W * WETA0R + BIN_WL * PETA0R
		U0R_i     = U0R + BIN_L * LU0R + BIN_W * WU0R + BIN_WL * PU0R
		UAR_i     = UAR + BIN_L * LUAR + BIN_W * WUAR + BIN_WL * PUAR
		UDR_i     = UDR + BIN_L * LUDR + BIN_W * WUDR + BIN_WL * PUDR
		UCSR_i    = UCSR + BIN_L * LUCSR + BIN_W * WUCSR + BIN_WL * PUCSR
		UCR_i     = UCR + BIN_L * LUCR + BIN_W * WUCR + BIN_WL * PUCR
		PCLMR_i   = PCLMR + BIN_L * LPCLMR + BIN_W * WPCLMR + BIN_WL * PPCLMR
		PDIBLCR_i = PDIBLCR + BIN_L * LPDIBLCR + BIN_W * WPDIBLCR + BIN_WL * PPDIBLCR
		VSATR_i   = VSATR + BIN_L * LVSATR + BIN_W * WVSATR + BIN_WL * PVSATR
		PSATR_i   = PSATR + BIN_L * LPSATR + BIN_W * WPSATR + BIN_WL * PPSATR
		PTWGR_i   = PTWGR + BIN_L * LPTWGR + BIN_W * WPTWGR + BIN_WL * PPTWGR
	end
	
	# Geometrical scaling
	T0        = NDEPL1 * max(pow(Inv_L, NDEPLEXP1) - pow(Inv_Llong, NDEPLEXP1), 0.0) + NDEPL2 * max(pow(Inv_L, NDEPLEXP2) - pow(Inv_Llong, NDEPLEXP2), 0.0)
	T1        = NDEPW * max(pow(Inv_W, NDEPWEXP) - pow(Inv_Wwide, NDEPWEXP), 0.0) + NDEPWL * pow(Inv_W * Inv_L, NDEPWLEXP)
	NDEP_i    = NDEP_i * (1.0 + T0 + T1)
	T0        = NFACTORL * max( pow(Inv_L, NFACTORLEXP) - pow(Inv_Llong, NFACTORLEXP), 0.0)
	T1        = NFACTORW * max( pow(Inv_W, NFACTORWEXP) - pow(Inv_Wwide, NFACTORWEXP), 0.0) + NFACTORWL * pow(Inv_WL, NFACTORWLEXP)
	NFACTOR_i = NFACTOR_i * (1.0 + T0 + T1)
	T0        = (1.0 + CDSCDL * max(pow(Inv_L, CDSCDLEXP) - pow(Inv_Llong, CDSCDLEXP), 0.0))
	CDSCD_i   = CDSCD_i * T0
	if (ASYMMOD != 0) 
		CDSCDR_i = CDSCDR_i * T0
	end
	CDSCB_i = CDSCB_i * (1.0 + CDSCBL * max(pow(Inv_L, CDSCBLEXP) - pow(Inv_Llong, CDSCBLEXP), 0.0))
	U0_i    = MULU0 * U0_i
	if (MOBSCALE != 1) 
		if (U0LEXP > 0.0) 
			U0_i = U0_i * (1.0 - U0L * max(pow(Inv_L, U0LEXP) - pow(Inv_Llong, U0LEXP), 0.0))
			if (ASYMMOD != 0) 
				U0R_i = U0R_i * (1.0 - U0L * max(pow(Inv_L, U0LEXP) - pow(Inv_Llong, U0LEXP), 0.0))
			end
		else 
			U0_i = U0_i * (1.0 - U0L)
			if (ASYMMOD != 0) 
				U0R_i = U0R_i * (1.0 - U0L)
			end
		end
	else 
		U0_i = U0_i * (1.0 - (UP1 * lexp(-Leff / LP1)) - (UP2 * lexp(-Leff / LP2)))
		if (ASYMMOD != 0) 
			U0R_i = U0R_i * (1.0 - (UP1 * lexp(-Leff / LP1)) - (UP2 * lexp(-Leff / LP2)))
		end
	end
	T0   = UAL * max(pow(Inv_L, UALEXP) - pow(Inv_Llong, UALEXP), 0.0)
	T1   = UAW * max(pow(Inv_W, UAWEXP) - pow(Inv_Wwide, UAWEXP), 0.0) + UAWL * pow(Inv_WL, UAWLEXP)
	UA_i = UA_i * (1.0 + T0 + T1)
	if (ASYMMOD != 0) 
		UAR_i = UAR_i * (1.0 + T0 + T1)
	end
	T0   = EUL * max(pow(Inv_L, EULEXP) - pow(Inv_Llong, EULEXP), 0.0)
	T1   = EUW * max(pow(Inv_W, EUWEXP) - pow(Inv_Wwide, EUWEXP), 0.0) + EUWL * pow(Inv_WL, EUWLEXP)
	EU_i = EU_i * (1.0 + T0 + T1)
	T0   = 1.0 + UDL * max(pow(Inv_L, UDLEXP) - pow(Inv_Llong, UDLEXP), 0.0)
	UD_i = UD_i * T0
	if (ASYMMOD != 0) 
		UDR_i = UDR_i * T0
	end
	T0   = UCL * max(pow(Inv_L, UCLEXP) - pow(Inv_Llong, UCLEXP), 0.0)
	T1   = UCW * max(pow(Inv_W, UCWEXP) - pow(Inv_Wwide, UCWEXP), 0.0) + UCWL * pow(Inv_WL, UCWLEXP)
	UC_i = UC_i * (1.0 + T0 + T1)
	if (ASYMMOD != 0) 
		UCR_i = UCR_i * (1.0 + T0 + T1)
	end
	T0     = max(pow(Inv_L, DSUB) - pow(Inv_Llong, DSUB), 0.0)
	ETA0_i = ETA0_i * T0
	if (ASYMMOD != 0) 
		ETA0R_i = ETA0R_i * T0
	end
	ETAB_i   = ETAB_i * max(pow(Inv_L, ETABEXP) - pow(Inv_Llong, ETABEXP), 0.0)
	T0       = 1.0 + PDIBLCL * max(pow(Inv_L, PDIBLCLEXP) - pow(Inv_Llong, PDIBLCLEXP), 0.0)
	PDIBLC_i = PDIBLC_i * T0
	if (ASYMMOD != 0) 
		PDIBLCR_i = PDIBLCR_i * T0
	end
	T0       = DELTA_i * (1.0 + DELTAL * max(pow(Inv_L, DELTALEXP) - pow(Inv_Llong, DELTALEXP), 0.0))
	DELTA_i  = min(T0, 0.5)
	FPROUT_i = FPROUT_i * (1.0 + FPROUTL * max(pow(Inv_L, FPROUTLEXP) - pow(Inv_Llong, FPROUTLEXP), 0.0))
	T0       = (1.0 + PCLML * max(pow(Inv_L, PCLMLEXP) - pow(Inv_Llong, PCLMLEXP), 0.0))
	PCLM_i   = PCLM_i * T0
	PCLM_i   = max(PCLM_i, 0.0)
	if (ASYMMOD != 0) 
		PCLMR_i = PCLMR_i * T0
		PCLMR_i = max(PCLMR_i, 0.0)
	end
	T0     = VSATL * max(pow(Inv_L, VSATLEXP) - pow(Inv_Llong, VSATLEXP), 0.0)
	T1     = VSATW * max(pow(Inv_W, VSATWEXP) - pow(Inv_Wwide, VSATWEXP), 0.0) + VSATWL * pow(Inv_WL, VSATWLEXP)
	VSAT_i = VSAT_i * (1.0 + T0 + T1)
	if (ASYMMOD != 0) 
		VSATR_i = VSATR_i * (1.0 + T0 + T1)
	end
	PSAT_i = max(PSAT_i * (1.0 + PSATL * max(pow(Inv_L, PSATLEXP) - pow(Inv_Llong, PSATLEXP), 0.0)), 0.25)
	if (ASYMMOD != 0) 
		PSATR_i = max(PSATR_i * (1.0 + PSATL * max(pow(Inv_L, PSATLEXP) - pow(Inv_Llong, PSATLEXP), 0.0)), 0.25)
	end
	T0     = (1.0 + PTWGL * max(pow(Inv_L, PTWGLEXP) - pow(Inv_Llong, PTWGLEXP), 0.0))
	PTWG_i = PTWG_i * T0
	if (ASYMMOD != 0) 
		PTWGR_i = PTWGR_i * T0
	end
	ALPHA0_i = ALPHA0_i * (1.0 + ALPHA0L * max(pow(Inv_L, ALPHA0LEXP) - pow(Inv_Llong, ALPHA0LEXP), 0.0))
	AGIDL_i  = AGIDL_i * (1.0 + AGIDLL * Inv_L + AGIDLW * Inv_W)
	AGISL_i  = AGISL_i * (1.0 + AGISLL * Inv_L + AGISLW * Inv_W)
	AIGC_i   = AIGC_i * (1.0 + AIGCL * Inv_L + AIGCW * Inv_W)
	AIGS_i   = AIGS_i * (1.0 + AIGSL * Inv_L + AIGSW * Inv_W)
	AIGD_i   = AIGD_i * (1.0 + AIGDL * Inv_L + AIGDW * Inv_W)
	PIGCD_i  = PIGCD * (1.0 + PIGCDL * Inv_L)
	T0       = NDEPCVL1 * max(pow(Inv_Lact, NDEPCVLEXP1) - pow(Inv_Llong, NDEPCVLEXP1), 0.0) + NDEPCVL2 * max( pow(Inv_Lact, NDEPCVLEXP2) - pow(Inv_Llong, NDEPCVLEXP2), 0.0)
	T1       = NDEPCVW * max(pow(Inv_Wact, NDEPCVWEXP) - pow(Inv_Wwide, NDEPCVWEXP), 0.0) + NDEPCVWL * pow(Inv_Wact * Inv_Lact, NDEPCVWLEXP)
	NDEPCV_i = NDEPCV_i * (1.0 + T0 + T1)
	T0       = VFBCVL * max(pow(Inv_Lact, VFBCVLEXP) - pow(Inv_Llong, VFBCVLEXP), 0.0)
	T1       = VFBCVW * max(pow(Inv_Wact, VFBCVWEXP) - pow(Inv_Wwide, VFBCVWEXP), 0.0) + VFBCVWL * pow(Inv_WL, VFBCVWLEXP)
	VFBCV_i  = VFBCV_i * (1.0 + T0 + T1)
	T0       = VSATCVL * max(pow(Inv_Lact, VSATCVLEXP) - pow(Inv_Llong, VSATCVLEXP), 0.0)
	T1       = VSATCVW * max(pow(Inv_W, VSATCVWEXP) - pow(Inv_Wwide, VSATCVWEXP), 0.0) + VSATCVWL * pow(Inv_WL, VSATCVWLEXP)
	VSATCV_i = VSATCV_i * (1.0 + T0 + T1)
	PCLMCV_i = PCLMCV_i * (1.0 + PCLMCVL * max(pow(Inv_Lact, PCLMCVLEXP) - pow(Inv_Llong, PCLMCVLEXP), 0.0))
	PCLMCV_i = max(PCLMCV_i, 0.0)
	T0       = K1L * max(pow(Inv_L, K1LEXP) - pow(Inv_Llong, K1LEXP), 0.0)
	T1       = K1W * max(pow(Inv_W, K1WEXP) - pow(Inv_Wwide, K1WEXP), 0.0) + K1WL * pow(Inv_WL, K1WLEXP)
	K1_i     = K1_i * (1.0 + T0 + T1)
	T0       = K2L * max(pow(Inv_L, K2LEXP) - pow(Inv_Llong, K2LEXP), 0.0)
	T1       = K2W * max(pow(Inv_W, K2WEXP) - pow(Inv_Wwide, K2WEXP), 0.0) + K2WL * pow(Inv_WL, K2WLEXP)
	K2_i     = K2_i * (1.0 + T0 + T1)
	PRWB_i   = PRWB_i * (1.0 + PRWBL * max( pow(Inv_L, PRWBLEXP) - pow(Inv_Llong, PRWBLEXP), 0))
	
	# Global scaling parameters for temperature
	UTE_i   = UTE_i * (1.0 + Inv_L * UTEL)
	UA1_i   = UA1_i * (1.0 + Inv_L * UA1L)
	UD1_i   = UD1_i * (1.0 + Inv_L * UD1L)
	AT_i    = AT_i * (1.0 + Inv_L * ATL)
	PTWGT_i = PTWGT_i * (1.0 + Inv_L * PTWGTL)
	if ($port_connected(t) == 0) 
		if (SHMOD == 0 || RTH0 == 0.0) 
			Temp(t) <+ 0.0
		else 
			println("5 terminal Module, while 't' node is not connected, SH is activated.")
		end
	end
	if (RDSMOD == 1) 
		RSW_i = RSW_i * (1.0 + RSWL * max(pow(Inv_L, RSWLEXP) - pow(Inv_Llong, RSWLEXP), 0.0))
		RDW_i = RDW_i * (1.0 + RDWL * max(pow(Inv_L, RDWLEXP) - pow(Inv_Llong, RDWLEXP), 0.0))
	else 
		RDSW_i = RDSW_i * (1.0 + RDSWL * max(pow(Inv_L, RDSWLEXP) - pow(Inv_Llong, RDSWLEXP), 0.0))
	end
	
	# Parameter checking
	if (UCS_i < 1.0) 
		UCS_i = 1.0
	elseif (UCS_i > 2.0) 
		UCS_i = 2.0
	end
	if (ASYMMOD != 0) 
		if (UCSR_i < 1.0) 
			UCSR_i = 1.0
		elseif (UCSR_i > 2.0) 
			UCSR_i = 2.0
		end
	end
	if (CGIDL_i < 0.0) 
		println("Fatal: CGIDL_i = %e is negative.", CGIDL_i)
		
	end
	if (CGISL_i < 0.0) 
		println("Fatal: CGISL_i = %e is negative.", CGISL_i)
		
	end
	if (CKAPPAD_i <= 0.0) 
		println("Fatal: CKAPPAD_i = %e is non-positive.", CKAPPAD_i)
		
	end
	if (CKAPPAS_i <= 0.0) 
		println("Fatal: CKAPPAS_i = %e is non-positive.", CKAPPAS_i)
		
	end
	if (PDITS_i < 0.0) 
		println("Fatal: PDITS_i = %e is negative.", PDITS_i)
		
	end
	if (CIT_i < 0.0) 
		println("Fatal: CIT_i = %e is negative.", CIT_i)
		
	end
	if (NFACTOR_i < 0.0) 
		println("Fatal: NFACTOR_i = %e is negative.", NFACTOR_i)
		
	end
	if (K1_i < 0.0) 
		println("Fatal: K1_i = %e is negative.", K1_i)
		
	end
	
	if (NSD_i <= 0.0) 
		println("Fatal: NSD_i = %e is non-positive.", NSD_i)
		
	end
	if (NDEP_i <= 0.0) 
		println("Fatal: NDEP_i = %e is non-positive.", NDEP_i)
		
	end
	if (NDEPCV_i <= 0.0) 
		println("Fatal: NDEPCV_i = %e is non-positive.", NDEPCV_i)
		
	end
	if (IGBMOD != 0) 
		if (NIGBINV_i <= 0.0) 
			println("Fatal: NIGBINV_i = %e is non-positive.", NIGBINV_i)
			
		end
		if (NIGBACC_i <= 0.0) 
			println("Fatal: NIGBACC_i = %e is non-positive.", NIGBACC_i)
			
		end
	end
	if (IGCMOD != 0) 
		if (POXEDGE_i <= 0.0) 
			println("Fatal: POXEDGE_i = %e is non-positive.", POXEDGE_i)
			
		end
	end
	if (CDSCD_i < 0.0) 
		println("Fatal: CDSCD_i = %e is negative.", CDSCD_i)
		
	end
	if (ASYMMOD != 0) 
		if (CDSCDR_i < 0.0) 
			println("Fatal: CDSCDR_i = %e is negative.", CDSCDR_i)
			
		end
	end
	if (DLCIG_i < 0.0) 
		println("Warning: DLCIG = %e is negative, setting it to 0.", DLCIG_i)
		DLCIG_i = 0.0
	end
	if (DLCIGD_i < 0.0) 
		println("Warning: DLCIGD = %e is negative, setting it to 0.", DLCIGD_i)
		DLCIGD_i = 0.0
	end
	if (M0_i < 0.0) 
		println("Warning: M0_i = %e is negative, setting it to 0.", M0_i)
		M0_i = 0.0
	end
	if (U0_i <= 0.0) 
		println("Warning: U0_i = %e is non-positive, setting it to the default value.", U0_i)
		U0_i = 0.067
	end
	if (UA_i < 0.0) 
		println("Warning: UA_i = %e is negative, setting it to 0.", UA_i)
		UA_i = 0.0
	end
	if (EU_i < 0.0) 
		println("Warning: EU_i = %e is negative, setting it to 0.", EU_i)
		EU_i = 0.0
	end
	if (UD_i < 0.0) 
		println("Warning: UD_i = %e is negative, setting it to 0.", UD_i)
		UD_i = 0.0
	end
	if (UCS_i < 0.0) 
		println("Warning: UCS_i = %e is negative, setting it to 0.", UCS_i)
		UCS_i = 0.0
	end
	
	# Initialize variables used in geometry macros
	nuEndD = 0.0 nuEndS = 0.0 nuIntD = 0.0 nuIntS = 0.0 Rend = 0.0 Rint = 0.0
	
	# Process drain series resistance
	DMCGeff = DMCG - DMCGT
	DMCIeff = DMCI
	DMDGeff = DMDG - DMCGT
	
	# Processing S/D resistances and conductances
	if($param_given(NRS)) 
		RSourceGeo = RSH * NRS
	elseif (RGEOMOD > 0 && RSH > 0.0) 
		BSIMBULKRdseffGeo(NF, GEOMOD, RGEOMOD, MINZ, Weff, RSH, DMCGeff, DMCIeff, DMDGeff, 1, RSourceGeo)
	else 
		RSourceGeo = 0.0
	end
	
	if ($param_given(NRD)) 
		RDrainGeo = RSH * NRD
	elseif (RGEOMOD > 0 && RSH > 0.0) 
		BSIMBULKRdseffGeo(NF, GEOMOD, RGEOMOD, MINZ, Weff, RSH, DMCGeff, DMCIeff, DMDGeff, 0, RDrainGeo)
	else 
		RDrainGeo = 0.0
	end
	# Clamping of S/D resistances 
	
	if (RDSMOD == 0) 
		if (RSourceGeo < minr) 
			RSourceGeo = 0
		end
		if (RDrainGeo < minr) 
			RDrainGeo = 0
		end
	else 
		if (RSourceGeo <= minr) 
			RSourceGeo = minr
		end
		if (RDrainGeo <= minr) 
			RDrainGeo = minr
		end
	end
	
	
	if (RDSMOD == 1) 
		if (RSWMIN_i <= 0.0) 
			RSWMIN_i = 0.0
		end
		if (RDWMIN_i <= 0.0) 
			RDWMIN_i = 0.0
		end
		if (RSW_i <= 0.0) 
			RSW_i = 0.0
		end
		if (RDW_i <= 0.0) 
			RDW_i = 0.0
		end
	else 
		if (RDSWMIN_i <= 0.0) 
			RDSWMIN_i = 0.0
		end
		if (RDSW_i <= 0.0) 
			RDSW_i = 0.0
		end
	end
	
	# Body resistance network
	Grbsb = 0.0 Grbdb = 0.0 Grbpb = 0.0 Grbps = 0.0 Grbpd = 0.0
	
	if (RBODYMOD != 0) 
		Lnl  = lln(Leff * 1.0e6)
		Lnw  = lln(Weff * 1.0e6)
		Lnnf = lln(NF)
		Bodymode = 5
		Rbpb = RBPB
		Rbpd = RBPD
		Rbps = RBPS
		Rbdb = RBDB
		Rbsb = RBSB
		if (!$param_given(RBPS0) || !$param_given(RBPD0)) 
			Bodymode = 1
		end
		elseif (!$param_given(RBSBX0) && !$param_given(RBSBY0) || !$param_given(RBDBX0) && !$param_given(RBDBY0)) 
			Bodymode = 3
		end
		if (RBODYMOD == 2) 
			if (Bodymode == 5) 
				Rbsbx = RBSBX0 * lexp(RBSDBXL * Lnl + RBSDBXW * Lnw + RBSDBXNF * Lnnf)
				Rbsby = RBSBY0 * lexp(RBSDBYL * Lnl + RBSDBYW * Lnw + RBSDBYNF * Lnnf)
				Rbsb  = Rbsbx * Rbsby / (Rbsbx + Rbsby)
				Rbdbx = RBDBX0 * lexp(RBSDBXL * Lnl + RBSDBXW * Lnw + RBSDBXNF * Lnnf)
				Rbdby = RBDBY0 * lexp(RBSDBYL * Lnl + RBSDBYW * Lnw + RBSDBYNF * Lnnf)
				Rbdb  = Rbdbx * Rbdby / (Rbdbx + Rbdby)
			end
			if (Bodymode == 3 || Bodymode == 5) 
				Rbps = RBPS0 * lexp(RBPSL * Lnl + RBPSW * Lnw + RBPSNF * Lnnf)
				Rbpd = RBPD0 * lexp(RBPDL * Lnl + RBPDW * Lnw + RBPDNF * Lnnf)
			end
			Rbpbx = RBPBX0 * lexp(RBPBXL * Lnl + RBPBXW * Lnw + RBPBXNF * Lnnf)
			Rbpby = RBPBY0 * lexp(RBPBYL * Lnl + RBPBYW * Lnw + RBPBYNF * Lnnf)
			Rbpb  = Rbpbx * Rbpby / (Rbpbx + Rbpby)
		end
		if (RBODYMOD == 1 || (RBODYMOD == 2 && Bodymode == 5)) 
			if (Rbdb < 1.0e-3) 
				Grbdb = 1.0e3  # in mho
			else 
				Grbdb = GBMIN + 1.0 / Rbdb
			end
			if (Rbpb < 1.0e-3) 
				Grbpb = 1.0e3
			else 
				Grbpb = GBMIN + 1.0 / Rbpb
			end
			if (Rbps < 1.0e-3) 
				Grbps = 1.0e3
			else 
				Grbps = GBMIN + 1.0 / Rbps
			end
			if (Rbsb < 1.0e-3) 
				Grbsb = 1.0e3
			else 
				Grbsb = GBMIN + 1.0 / Rbsb
			end
			if (Rbpd < 1.0e-3) 
				Grbpd = 1.0e3
			else 
				Grbpd = GBMIN + 1.0 / Rbpd
			end
		elseif (RBODYMOD == 2 && Bodymode == 3) 
			Grbdb = GBMIN
			Grbsb = GBMIN
			if (Rbpb < 1.0e-3) 
				Grbpb = 1.0e3
			else 
				Grbpb = GBMIN + 1.0 / Rbpb
			end
			if (Rbps < 1.0e-3) 
				Grbps = 1.0e3
			else 
				Grbps = GBMIN + 1.0 / Rbps
			end
			if (Rbpd < 1.0e-3) 
				Grbpd = 1.0e3
			else 
				Grbpd = GBMIN + 1.0 / Rbpd
			end
		elseif (RBODYMOD == 2 && Bodymode == 1) 
			Grbdb = GBMIN
			Grbsb = GBMIN
			Grbps = 1.0e3
			Grbpd = 1.0e3
			if (Rbpb < 1.0e-3) 
				Grbpb = 1.0e3
			else 
				Grbpb = GBMIN + 1.0 / Rbpb
			end
		end
	end
	
	# Gate process resistance
	Grgeltd = RSHG * (XGW + Weffcj / 3.0 / NGCON) / (NGCON * NF * (Lnew - XGL))
	if (Grgeltd > 0.0) 
		Grgeltd = 1.0 / Grgeltd
	else 
		Grgeltd = 1.0e3
		if (RGATEMOD != 0) 
			STROBE("Warning: (instance %M) The gate conductance reset to 1.0e3 mho.")
		end
	end
	T0           = TOXE * TOXE
	T1           = TOXE * POXEDGE_i
	T2           = T1 * T1
	ToxRatio     = lexp(NTOX_i * lln(TOXREF / TOXE)) / T0
	ToxRatioEdge = lexp(NTOX_i * lln(TOXREF / T1)) / T2
	Aechvb       = (TYPE == ntype) ? 4.97232e-7 : 3.42537e-7
	Bechvb       = (TYPE == ntype) ? 7.45669e11 : 1.16645e12
	AechvbEdge   = Aechvb * Weff * ToxRatioEdge
	BechvbEdge   = -Bechvb * TOXE * POXEDGE_i
	Aechvb       = Aechvb * (Weff * Leff * ToxRatio)
	Bechvb       = -Bechvb * TOXE
	Weff_SH      = WTH0 + Weff
	
	# Parameters for self-heating effects
	if((SHMOD != 0) && (RTH0 > 0.0) && (Weff_SH > 0.0)) 
		gth = Weff_SH * NF / RTH0
		cth = CTH0 * Weff_SH * NF
	else 
		# Set gth to some value to prevent a singular G matrix
		gth = 1.0
		cth = 0.0
	end
	
	# Temperature-dependent calculations
	if (TNOM <= -P_CELSIUS0) 
		T0 = REFTEMP - P_CELSIUS0
		println("Warning: TNOM = %e C <= %e C. Setting TNOM to %e C.", TNOM, -P_CELSIUS0, T0)
		Tnom = REFTEMP
	else 
		Tnom = TNOM + P_CELSIUS0
	end
	DevTemp = $temperature + DTEMP
	
	# Calculate temperature dependent values for self-heating effects
	if ((SHMOD != 0) && (RTH0 > 0.0) && (Weff_SH > 0.0)) 
		delTemp1 = Temp(t)
	else 
		delTemp1 = 0.0
	end
	DevTemp    = delTemp1 + DevTemp
	T_DELTA_SH = Temp(t)
	T_TOTAL_K  = DevTemp
	T_TOTAL_C  = DevTemp - P_CELSIUS0
	Vt         = KboQ * DevTemp
	inv_Vt     = 1.0 / Vt
	TRatio     = DevTemp / Tnom
	delTemp    = DevTemp - Tnom
	Vtm        = KboQ * DevTemp
	Vtm0       = KboQ * Tnom
	Eg         = BG0SUB - TBGASUB * DevTemp * DevTemp / (DevTemp + TBGBSUB)
	Eg0        = BG0SUB - TBGASUB * Tnom * Tnom / (Tnom + TBGBSUB)
	T1         = (DevTemp / Tnom) * sqrt(DevTemp / Tnom)
	ni         = NI0SUB * T1 * lexp(Eg / (2.0 * Vtm0) - Eg / (2.0 * Vtm))
	if ((SHMOD != 0) && (RTH0 > 0.0) && (Weff_SH > 0.0)) 
		T0   = lln(NDEP_i / ni)
		phib = sqrt(T0 * T0 + 1.0e-6)
	else 
		phib = lln(NDEP_i / ni)
	end
	if ((SHMOD != 0) && (RTH0 > 0.0) && (Weff_SH > 0.0)) 
		T0  = lln(NDEPEDGE_i * NSD_i / (ni * ni))
		Vbi_edge = sqrt(T0 * T0 + 1.0e-6)
	else 
		Vbi_edge = lln(NDEPEDGE_i * NSD_i / (ni * ni))
	end
	if (NGATE_i > 0.0) 
		Vfbsdr = -devsign * Vt * lln(NGATE_i / NSD_i) + VFBSDOFF
	else 
		Vfbsdr = 0.0
	end
	
	# Short channel effects
	Phist     = max(0.4 + Vt * phib + PHIN_i, 0.4)
	sqrtPhist = sqrt(Phist)
	T1DEP     = sqrt(2.0 * epssi / (q * NDEP_i))
	litl      = sqrt((epssi / epsox) * TOXE * XJ_i)
	NFACTOR_t = NFACTOR_i * hypsmooth((1.0 + TNFACTOR * (TRatio - 1.0)), 1e-3)
	ETA0_t    = ETA0_i * (1.0 + TETA0 * (TRatio - 1.0))
	if (ASYMMOD != 0) 
		ETA0R_t = ETA0R_i * (1.0 + TETA0 * (TRatio - 1.0))
	end
	
	# Mobility degradation
	eta_mu = (TYPE != ntype) ? (Oneby3 * ETAMOB) : (0.5 * ETAMOB)
	U0_t   = U0_i * pow(TRatio, UTE_i)
	UA_t   = UA_i * hypsmooth(1.0 + UA1_i * delTemp - 1.0e-6, 1.0e-3)
	UC_t   = UC_i * hypsmooth(1.0 + UC1_i * delTemp - 1.0e-6, 1.0e-3)
	UD_t   = UD_i * pow(TRatio, UD1_i)
	UCS_t  = UCS_i * pow(TRatio, UCSTE_i)
	EU_t   = EU_i * hypsmooth((1.0 + EU1_i * (TRatio - 1.0)), 1e-3)
	if (ASYMMOD != 0) 
		U0R_t  = U0R_i * pow(TRatio, UTE_i)
		UAR_t  = UAR_i * hypsmooth(1.0 + UA1_i * delTemp - 1.0e-6, 1.0e-3)
		UCR_t  = UCR_i * hypsmooth(1.0 + UC1_i * delTemp - 1.0e-6, 1.0e-3)
		UDR_t  = UDR_i * pow(TRatio, UD1_i)
		UCSR_t = UCSR_i * pow(TRatio, UCSTE_i)
	end
	rdstemp = pow(TRatio, PRT_i)
	VSAT_t  = VSAT_i * pow(TRatio, -AT_i)
	if (VSAT_t < 100.0) 
		println("Warning: VSAT(%f) = %e is less than 100, setting it to 100.", DevTemp, VSAT_t)
		VSAT_t = 100.0
	end
	if (HVMOD == 1) 
		rdstemphv = pow(TRatio, PRTHV)
		VDRIFT_t  = VDRIFT * pow(TRatio, -ATHV)
	end
	if (ASYMMOD != 0) 
		VSATR_t = VSATR_i * pow(TRatio, -AT_i)
		if(VSATR_t < 100.0) 
			println("Warning: VSATR(%f) = %e is less than 100, setting it to 100.", DevTemp, VSATR_t)
			VSATR_t = 100.0
		end
	end
	VSATCV_t = VSATCV_i * pow(TRatio, -AT_i)
	if (VSATCV_t < 100.0) 
		println("Warning: VSATCV(%f) = %e is less than 100, setting it to 100.", DevTemp, VSATCV_t)
		VSATCV_t = 100.0
	end
	DELTA_t = 1.0 / ( hypsmooth((1.0 / DELTA_i) * (1.0 + TDELTA * delTemp) - 2.0 , 1.0e-3) + 2.0)
	PTWG_t  = PTWG_i * hypsmooth(1.0 - PTWGT_i * delTemp - 1.0e-6, 1.0e-3)
	if (ASYMMOD != 0) 
		PTWGR_t = PTWGR_i * hypsmooth(1.0 - PTWGT_i * delTemp - 1.0e-6, 1.0e-3)
	end
	A1_t    = A1_i * hypsmooth(1.0 + A11_i * delTemp - 1.0e-6, 1.0e-3)
	A2_t    = A2_i * hypsmooth(1.0 + A21_i * delTemp - 1.0e-6, 1.0e-3)
	BETA0_t = BETA0_i * pow(TRatio, IIT_i)
	BGIDL_t = BGIDL_i * hypsmooth(1.0 + TGIDL_i * delTemp - 1.0e-6, 1.0e-3)
	BGISL_t = BGISL_i * hypsmooth(1.0 + TGIDL_i * delTemp - 1.0e-6, 1.0e-3)
	igtemp  = lexp(IGT_i * lln(TRatio))
	K0_t    = K0_i * hypsmooth(1.0 + K01_i * delTemp - 1.0e-6, 1.0e-3)
	M0_t    = M0_i * hypsmooth(1.0 + M01_i * delTemp - 1.0e-6, 1.0e-3)
	C0_t    = C0_i * hypsmooth(1.0 + C01_i * delTemp - 1.0e-6, 1.0e-3)
	C0SI_t  = C0SI_i * hypsmooth(1.0 + C0SI1_i * delTemp - 1.0e-6, 1.0e-3)
	C0SISAT_t = C0SISAT_i * hypsmooth(1.0 + C0SISAT1_i * delTemp - 1.0e-6, 1.0e-3)
	
	# Diode model temperature effects
	CJS_t     = CJS * hypsmooth(1.0 + TCJ * delTemp - 1.0e-6, 1.0e-3)
	CJD_t     = CJD * hypsmooth(1.0 + TCJ * delTemp - 1.0e-6, 1.0e-3)
	CJSWS_t   = CJSWS * hypsmooth(1.0 + TCJSW * delTemp - 1.0e-6, 1.0e-3)
	CJSWD_t   = CJSWD * hypsmooth(1.0 + TCJSW * delTemp - 1.0e-6, 1.0e-3)
	CJSWGS_t  = CJSWGS * hypsmooth(1.0 + TCJSWG * delTemp - 1.0e-6, 1.0e-3)
	CJSWGD_t  = CJSWGD * hypsmooth(1.0 + TCJSWG * delTemp - 1.0e-6, 1.0e-3)
	PBS_t     = hypsmooth(PBS - TPB * delTemp - 0.01, 1.0e-3) + 0.01
	PBD_t     = hypsmooth(PBD - TPB * delTemp - 0.01, 1.0e-3) + 0.01
	PBSWS_t   = hypsmooth(PBSWS - TPBSW * delTemp - 0.01, 1.0e-3) + 0.01
	PBSWD_t   = hypsmooth(PBSWD - TPBSW * delTemp - 0.01, 1.0e-3) + 0.01
	PBSWGS_t  = hypsmooth(PBSWGS - TPBSWG * delTemp - 0.01, 1.0e-3) + 0.01
	PBSWGD_t  = hypsmooth(PBSWGD - TPBSWG * delTemp - 0.01, 1.0e-3) + 0.01
	T0        = Eg0 / Vtm0 - Eg / Vtm
	T1        = lln(TRatio)
	T3        = lexp((T0 + XTIS * T1) / NJS)
	JSS_t     = JSS * T3
	JSWS_t    = JSWS * T3
	JSWGS_t   = JSWGS * T3
	T3        = lexp((T0 + XTID * T1) / NJD)
	JSD_t     = JSD * T3
	JSWD_t    = JSWD * T3
	JSWGD_t   = JSWGD * T3
	JTSS_t    = JTSS * lexp(Eg0 * XTSS * (TRatio - 1.0) / Vtm)
	JTSSWS_t  = JTSSWS * lexp(Eg0 * XTSSWS * (TRatio - 1.0) / Vtm)
	JTSSWGS_t = JTSSWGS * (sqrt(JTWEFF / Weffcj) + 1.0) * lexp(Eg0 * XTSSWGS * (TRatio - 1) / Vtm)
	JTSD_t    = JTSD * lexp(Eg0 * XTSD * (TRatio - 1.0) / Vtm)
	JTSSWD_t  = JTSSWD * lexp(Eg0 * XTSSWD * (TRatio - 1.0) / Vtm)
	JTSSWGD_t = JTSSWGD * (sqrt(JTWEFF / Weffcj) + 1.0) * lexp(Eg0 * XTSSWGD * (TRatio - 1) / Vtm)
	
	# All NJT*'s smoothed to 0.01 to prevent divide by zero/negative values
	NJTS_t     = hypsmooth(NJTS * (1.0 + TNJTS * (TRatio - 1.0)) - 0.01, 1.0e-3) + 0.01
	NJTSSW_t   = hypsmooth(NJTSSW * (1.0 + TNJTSSW * (TRatio - 1.0)) - 0.01, 1.0e-3) + 0.01
	NJTSSWG_t  = hypsmooth(NJTSSWG * (1.0 + TNJTSSWG * (TRatio - 1.0)) - 0.01, 1.0e-3) + 0.01
	NJTSD_t    = hypsmooth(NJTSD * (1.0 + TNJTSD * (TRatio - 1.0)) - 0.01, 1.0e-3) + 0.01
	NJTSSWD_t  = hypsmooth(NJTSSWD * (1.0 + TNJTSSWD * (TRatio - 1.0)) - 0.01, 1.0e-3) + 0.01
	NJTSSWGD_t = hypsmooth(NJTSSWGD * (1.0 + TNJTSSWGD * (TRatio - 1.0)) - 0.01, 1.0e-3) + 0.01
	
	# Effective S/D junction area and perimeters
	BSIMBULKPAeffGeo(NF, GEOMOD, MINZ, Weffcj, DMCGeff, DMCIeff, DMDGeff, temp_PSeff, temp_PDeff, temp_ASeff, temp_ADeff)
	if ($param_given(AS)) 
		ASeff = AS * WMLT * LMLT
	else 
		ASeff = temp_ASeff
	end
	if (ASeff < 0.0) 
		println("Warning: (instance %M) ASeff = %e is negative. Set to 0.0.", ASeff)
		ASeff = 0.0
	end
	if ($param_given(AD)) 
		ADeff = AD * WMLT * LMLT
	else 
		ADeff = temp_ADeff
	end
	if (ADeff < 0.0) 
		println("Warning: (instance %M) ADeff = %e is negative. Set to 0.0.", ADeff)
		ADeff = 0.0
	end
	if ($param_given(PS)) 
		if (PERMOD == 0) 
			# PS does not include gate-edge perimeters
			PSeff = PS * WMLT
		else 
			# PS includes gate-edge perimeters
			PSeff = max(PS * WMLT - Weffcj * NF, 0.0)
		end
	else 
		PSeff = temp_PSeff
		if (PSeff < 0.0) 
			println("Warning: (instance %M) PSeff = %e is negative. Set to 0.0.", PSeff)
			PSeff = 0.0
		end
	end
	if ($param_given(PD)) 
		if (PERMOD == 0) 
			# PD does not include gate-edge perimeters
			PDeff = PD * WMLT
		else 
			# PD includes gate-edge perimeters
			PDeff = max(PD * WMLT - Weffcj * NF, 0.0)
		end
	else 
		PDeff = temp_PDeff
		if (PDeff < 0.0) 
			println("Warning: (instance %M) PDeff = %e is negative. Set to 0.0.", PDeff)
			PDeff = 0.0
		end
	end
	Isbs = ASeff * JSS_t + PSeff * JSWS_t + Weffcj * NF * JSWGS_t
	if (Isbs > 0.0) 
		Nvtms    = Vtm * NJS
		XExpBVS  = lexp(-BVS / Nvtms) * XJBVS
		T2       = max(IJTHSFWD / Isbs, 10.0)
		Tb       = 1.0 + T2 - XExpBVS
		VjsmFwd  = Nvtms * lln(0.5 * (Tb + sqrt(Tb * Tb + 4.0 * XExpBVS)))
		T0       = lexp(VjsmFwd / Nvtms)
		IVjsmFwd = Isbs * (T0 - XExpBVS / T0 + XExpBVS - 1.0)
		SslpFwd  = Isbs * (T0 + XExpBVS / T0) / Nvtms
		T2       = hypsmooth(IJTHSREV / Isbs - 10.0, 1.0e-3) + 10.0
		VjsmRev  = -BVS - Nvtms * lln((T2 - 1.0) / XJBVS)
		T1       = XJBVS * lexp(-(BVS + VjsmRev) / Nvtms)
		IVjsmRev = Isbs * (1.0 + T1)
		SslpRev  = -Isbs * T1 / Nvtms
	else 
		Nvtms    = 0.0
		XExpBVS  = 0.0
		VjsmFwd  = 0.0
		IVjsmFwd = 0.0
		SslpFwd  = 0.0
		VjsmRev  = 0.0
		IVjsmRev = 0.0
		SslpRev  = 0.0
	end
	
	# Drain-side junction currents
	Isbd = ADeff * JSD_t + PDeff * JSWD_t + Weffcj * NF * JSWGD_t
	if (Isbd > 0.0) 
		Nvtmd    = Vtm * NJD
		XExpBVD  = lexp(-BVD / Nvtmd) * XJBVD
		T2       = max(IJTHDFWD / Isbd, 10.0)
		Tb       = 1.0 + T2 - XExpBVD
		VjdmFwd  = Nvtmd * lln(0.5 * (Tb + sqrt(Tb * Tb + 4.0 * XExpBVD)))
		T0       = lexp(VjdmFwd / Nvtmd)
		IVjdmFwd = Isbd * (T0 - XExpBVD / T0 + XExpBVD - 1.0)
		DslpFwd  = Isbd * (T0 + XExpBVD / T0) / Nvtmd
		T2       = hypsmooth(IJTHDREV / Isbd - 10.0, 1.0e-3) + 10.0
		VjdmRev  = -BVD - Nvtmd * lln((T2 - 1.0) / XJBVD)
		T1       = XJBVD * lexp(-(BVD + VjdmRev) / Nvtmd)
		IVjdmRev = Isbd * (1.0 + T1)
		DslpRev  = -Isbd * T1 / Nvtmd
	else 
		Nvtmd    = 0.0
		XExpBVD  = 0.0
		VjdmFwd  = 0.0
		IVjdmFwd = 0.0
		DslpFwd  = 0.0
		VjdmRev  = 0.0
		IVjdmRev = 0.0
		DslpRev  = 0.0
	end
	
	# STI stress equations
	if((SA > 0.0) && (SB > 0.0) && ((NF == 1.0) || ((NF > 1.0) && (SD > 0.0)))) 
		T0              = pow(Lnew, LLODKU0)
		W_tmp_stress    = Wnew + WLOD
		T1              = pow(W_tmp_stress, WLODKU0)
		tmp1_stress     = LKU0 / T0 + WKU0 / T1 + PKU0 / (T0 * T1)
		kstress_u0      = 1.0 + tmp1_stress
		T0              = pow(Lnew, LLODVTH)
		T1              = pow(W_tmp_stress, WLODVTH)
		tmp1_stress_vth = LKVTH0 / T0 + WKVTH0 / T1 + PKVTH0 / (T0 * T1)
		kstress_vth0    = 1.0 + tmp1_stress_vth
		T0              = TRatio - 1.0
		ku0_temp        = kstress_u0 * (1.0 + TKU0 * T0) + 1.0e-9
		for (i = 0 i < NF i = i + 1) : forloop
			T0     = 1.0 / NF / (SA + 0.5 * L_mult + i * (SD + L_mult))
			T1     = 1.0 / NF / (SB + 0.5 * L_mult + i * (SD + L_mult))
			Inv_sa = Inv_sa + T0
			Inv_sb = Inv_sb + T1
		end
		Inv_saref   = 1.0 / (SAREF + 0.5 * L_mult)
		Inv_sbref   = 1.0 / (SBREF + 0.5 * L_mult)
		Inv_odref   = Inv_saref + Inv_sbref
		rho_ref     = (KU0 / ku0_temp) * Inv_odref
		Inv_od      = Inv_sa + Inv_sb
		rho         = (KU0 / ku0_temp) * Inv_od
		mu0_mult    = (1.0 + rho) / (1.0 + rho_ref)
		vsat_mult   = (1.0 + rho * KVSAT) / (1.0 + rho_ref * KVSAT)
		vth0_stress = (KVTH0 / kstress_vth0) * (Inv_od - Inv_odref)
		k2_stress   = (STK2 / pow(kstress_vth0, LODK2)) * (Inv_od - Inv_odref)
		eta_stress  = (STETA0 / pow(kstress_vth0, LODETA0)) * (Inv_od - Inv_odref)
		U0_t        = U0_t * mu0_mult
		VSAT_t      = VSAT_t * vsat_mult
		K2_i        = K2_i + k2_stress
		ETA0_t      = ETA0_t + eta_stress
		if (EDGEFET == 1) 
			vth0_stress_EDGE = (KVTH0EDGE_i / kstress_vth0) * (Inv_od - Inv_odref)
			k2_stress_EDGE   = (STK2EDGE_i / pow(kstress_vth0, LODK2)) * (Inv_od - Inv_odref)
			eta_stress_EDGE  = (STETA0EDGE_i / pow(kstress_vth0, LODETA0)) * (Inv_od - Inv_odref)
		end
		K2EDGE_i   = K2EDGE_i + k2_stress_EDGE
		ETA0EDGE_i = ETA0EDGE_i + eta_stress_EDGE
	else 
		vth0_stress = 0.0
		vth0_stress_EDGE = 0.0
	end
	
	# Well proximity effect
	if (WPEMOD == 1) 
		Wdrn      = W / NF
		local_sca = SCA
		local_scb = SCB
		local_scc = SCC
		if (!$param_given(SCA) && !$param_given(SCB) && !$param_given(SCC)) 
			if($param_given(SC) && SC > 0.0) 
				T1        = SC + Wdrn
				T2        = 1.0 / SCREF
				local_sca = SCREF * SCREF / (SC * T1)
				local_scb = ((0.1 * SC + 0.01 * SCREF) * lexp(-10.0 * SC * T2)  - (0.1 * T1 + 0.01 * SCREF) *
							lexp(-10.0 * T1 * T2)) / Wdrn
				local_scc = ((0.05 * SC + 0.0025 * SCREF) * lexp(-20.0 * SC * T2)  - (0.05 * T1 + 0.0025 * SCREF) *
							lexp(-20.0 * T1 * T2)) / Wdrn
			else 
				STROBE("Warning: (Instance %M) No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive.")
			end
		end
	end
	vth0_well = KVTH0WE_i * (local_sca + WEB * local_scb + WEC * local_scc)
	k2_well   = K2WE_i * (local_sca + WEB * local_scb + WEC * local_scc)
	mu_well   = 1.0 + KU0WE_i * (local_sca + WEB * local_scb + WEC * local_scc)
	U0_t      = U0_t * mu_well
	K2_i      = K2_i + k2_well
	
	# Load terminal voltages
	Vg            = devsign * V(gi, bi)
	Vd            = devsign * V(di, bi)
	Vs            = devsign * V(si, bi)
	Vds           = Vd - Vs
	Vds_noswap    = Vds
	Vsb_noswap    = Vs
	Vdb_noswap    = Vd
	Vbs_jct       = devsign * V(sbulk, si)
	Vbd_jct       = devsign * V(dbulk, di)
	Vgd_noswap    = Vg - Vd
	Vgs_noswap    = Vg - Vs
	Vgd_ov_noswap = devsign * V(gm, di)
	Vgs_ov_noswap = devsign * V(gm, si)
	
	# Terminal voltage conditioning
	# Source-drain interchange
	sigvds = 1.0
	if (Vds < 0.0) 
		sigvds = -1.0
		Vd = devsign * V(si, bi)
		Vs = devsign * V(di, bi)
	end
	Vds  = Vd - Vs
	T0   = AVDSX * Vds
	if (T0 > EXPL_THRESHOLD) 
	   T1 = T0
	else 
	   T1 = ln(1.0 + exp(T0))
	end
	Vdsx = ((2.0 / AVDSX) * T1) - Vds - ((2.0 / AVDSX) * ln(2.0))
	Vbsx = -(Vs + 0.5 * (Vds - Vdsx))
	
	# Asymmetry model
	T0 = tanh(0.6 * Vds_noswap / Vtm)
	wf = 0.5 + 0.5 * T0
	wr = 1.0 - wf
	if (ASYMMOD != 0) 
		CDSCD_a  = CDSCDR_i * wr + CDSCD_i * wf
		ETA0_a   = ETA0R_t * wr + ETA0_t * wf
		PDIBLC_a = PDIBLCR_i * wr + PDIBLC_i * wf
		PCLM_a   = PCLMR_i * wr + PCLM_i * wf
		PSAT_a   = PSATR_i * wr + PSAT_i * wf
		VSAT_a   = VSATR_t * wr + VSAT_t * wf
		PTWG_a   = PTWGR_t * wr + PTWG_t * wf
		U0_a     = U0R_t * wr + U0_t * wf
		UA_a     = UAR_t * wr + UA_t * wf
		UC_a     = UCR_t * wr + UC_t * wf
		UD_a     = UDR_t * wr + UD_t * wf
		UCS_a    = UCSR_t * wr + UCS_t * wf
	else 
		CDSCD_a  = CDSCD_i
		ETA0_a   = ETA0_t
		PDIBLC_a = PDIBLC_i
		PCLM_a   = PCLM_i
		PSAT_a   = PSAT_i
		VSAT_a   = VSAT_t
		PTWG_a   = PTWG_t
		U0_a     = U0_t
		UA_a     = UA_t
		UC_a     = UC_t
		UD_a     = UD_t
		UCS_a    = UCS_t
	end
	
	# SCE, DIBL, SS degradation effects, Ref: BSIM4
	Smooth(Phist - Vbsx, 0.05, 0.1, PhistVbs)
	sqrtPhistVbs = sqrt(PhistVbs)
	Xdep         = T1DEP * sqrtPhistVbs
	Cdep         = epssi / Xdep
	cdsc         = CIT_i + NFACTOR_t + CDSCD_a * Vdsx - CDSCB_i * Vbsx
	T1           = 1.0 + cdsc/Cox
	Smooth(T1, 1.0, 0.05, n)
	nVt     = n * Vt
	inv_nVt = 1.0 / nVt
	
	# Vth shift for DIBL
	dVth_dibl = -(ETA0_a + ETAB_i * Vbsx) * Vdsx
	Smooth2(dVth_dibl, 0.0, 5.0e-5, dVth_dibl)
	
	# Vth shift with temperature
	dvth_temp = (KT1_i + KT1L / Leff + KT2_i * Vbsx) * (pow(TRatio, KT1EXP) - 1.0)
	
	
	# Vth correction for pocket implants
	if (DVTP0_i > 0.0) 
		T0 = -DVTP1_i * Vdsx
		if (T0 < -EXPL_THRESHOLD) 
			T2 = MIN_EXPL
		else 
			T2 = lexp(T0)
		end
		T3        = Leff + DVTP0_i * (1.0 + T2)
		dVth_ldop = -nVt * lln(Leff / T3)
	else 
		dVth_ldop = 0.0
	end
	T4        = DVTP5_i + DVTP2_i / pow(Leff, DVTP3_i)
	dVth_ldop = dVth_ldop - T4 * tanh(DVTP4_i * Vdsx)
	
	# Normalization of terminal and flatband voltage by nVt
	VFB_i = VFB_i + DELVTO
	vg    = Vg * inv_nVt
	vs    = Vs * inv_nVt
	vfb   = VFB_i * inv_nVt
	
	# Compute dVth_VNUD with "first-order" and "second-order" body-bias effect
	dVth_VNUD = K1_i * (sqrtPhistVbs - sqrtPhist) - K2_i * Vbsx
	Vth_shift = dVth_dibl + dVth_ldop + dVth_VNUD - dvth_temp + vth0_stress + vth0_well
	vgfb      = vg - vfb - Vth_shift * inv_nVt
	
	# Threshold voltage for operating point information
	gam     = sqrt(2.0 * q * epssi * NDEP_i * inv_Vt) / Cox
	q_vth   = 0.5
	T0      = hypsmooth((2.0 * phib + Vs * inv_Vt), 1.0e-3)
	nq      = 1.0 + gam / (2.0 * sqrt(T0))
	psip_th = hypsmooth((Vs * inv_Vt + 2.0 * phib + lln(q_vth) + 2.0 * q_vth + lln(2.0 * nq / gam * (2.0 * q_vth * nq / gam + 2.0 * sqrt(T0)))), 1.0e-3)
	VTH     = devsign * (VFB_i + (psip_th - Vs * inv_Vt) * Vt + Vt * gam * sqrt(psip_th) + Vth_shift)
	
	# Normalized body factor
	gam     = sqrt(2.0 * q * epssi * NDEP_i * inv_nVt) / Cox
	inv_gam = 1.0 / gam
	
	# psip: pinch-off voltage
	phib_n = phib / n
	PO_psip(vgfb, gam, 0.0, phib_n, psip)
	
	# Normalized inversion charge at source end of channel
	BSIM_q(psip, phib_n, vs, gam, qs)
	
	# Average charge-surface potential slope, Ref: Charge-based MOS Transistor Modeling by C. Enz & E. Vittoz
	Smooth(psip, 1.0, 2.0, psipclamp)
	sqrtpsip = sqrt(psipclamp)
	
	# Source side surface potential
	psiavg = psip - 2.0 * qs
	Smooth(psiavg, 1.0, 2.0, T0)
	nq = 1.0 + gam / (sqrtpsip + sqrt(T0))
	
	# Drain saturation voltage
	EeffFactor = 1.0e-8 / (epsratio * TOXE)
	T0 = nVt * (vgfb - psip - 2.0 * qs * (nq - 1.0))
	Smooth(T0, 0.0, 0.1, qbs)
	
	# Source side qi and qb for Vdsat- normalized to Cox
	qis = 2.0 * nq * nVt * qs
	Eeffs = EeffFactor * (qbs + eta_mu * qis)
	
	# Ref: BSIM4 mobility model
	T2 = pow(0.5 * (1.0 + (qis / qbs)), UCS_a)
	T3 = (UA_a + UC_a * Vbsx) * pow(Eeffs, EU_t) + UD_a / T2
	T4 = 1.0 + T3
	Smooth(T4, 1.0, 0.0015, Dmobs)
	WeffWRFactor = 1.0 / (pow(Weff * 1.0e6, WR_i) * NF)
	
	if (RDSMOD == 1) 
		Rdss = 0.0
	else 
		T0   = 1.0 + PRWG_i * qis
		T1   = PRWB_i * (sqrtPhistVbs - sqrtPhist)
		T2   = 1.0 / T0 + T1
		T3   = T2 + sqrt(T2 * T2 + 0.01)
		Rdss = (RDSWMIN_i + RDSW_i * T3) * WeffWRFactor * NF * rdstemp
		if (RDSMOD == 2) 
			Rdss = (RSourceGeo + (RDSWMIN_i + RDSW_i * T3) * WeffWRFactor * NF + RDrainGeo) * rdstemp
		end
	end
	T0  = pow(Dmobs, 1.0 / PSAT_a)
	T11 = PSATB_i * Vbsx
	T12 = sqrt(0.1+T11*T11)
	T1  = 0.5*(1-T11+sqrt((1-T11)*(1-T11)+T12))
	T2  = 10.0 * PSATX * qs * T1 / (10.0 * PSATX + qs * T1)
	if (PTWG_a < 0.0) 
		LambdaC = 2.0 * ((U0_a / T0) * nVt / (VSAT_a * Leff)) * (1.0 / (1.0 - PTWG_a * T2))
	else 
		LambdaC = 2.0 * ((U0_a / T0) * nVt / (VSAT_a * Leff)) * (1.0 + PTWG_a * T2)
	end
	
	# qdsat for external Rds
	if (Rdss == 0) 
		# Accurate qdsat derived from consistent I-V
		T0 = 0.5 * LambdaC * (qs * qs + qs) / (1.0 + 0.5 * LambdaC * (1.0 + qs))
		T1 = 2.0 * LambdaC * (qs - T0)
		T2 = sqrt(1.0 + T1 * T1)
		ln_T1_T2 = asinh(T1)
		if (T1 != 0.0) 
			T3 = T2 + (1.0 / T1) * ln_T1_T2
		else 
			T3 = T2 + (1.0 / T2)
		end
		T4 = T0 * T3 - LambdaC * ((qs * qs + qs) - (T0 * T0 + T0))
		if (T1 != 0.0) 
			T5 = -2.0 * LambdaC * (T1 * T2 - ln_T1_T2) / (T1 * T1)
		else 
			T5 = -2.0 * LambdaC * (T1/T2) * (T1/T2) *(T1/T2)
		end
		T6 = T0 * T5 + T3 + LambdaC * (2.0 * T0 + 1.0)
		T0 = T0 - (T4 / T6)
		T1 = 2.0 * LambdaC * (qs - T0)
		T2 = sqrt(1.0 + T1 * T1)
		ln_T1_T2 = asinh(T1)
		if (T1 != 0.0) 
			T3 = T2 + (1.0 / T1) * ln_T1_T2
		else 
			T3 = T2 + (1.0 / T2)
		end
		T4 = T0 * T3 - LambdaC * ((qs * qs + qs) - (T0 * T0 + T0))
		if (T1 != 0.0) 
			T5 = -2.0 * LambdaC * (T1 * T2 - ln_T1_T2) / (T1 * T1)
		else 
			T5 = (T1 / T2) * (T1 / T2) * (T1 / T2)
		end
		T6    = T0 * T5 + T3 + LambdaC * (2.0 * T0 + 1.0)
		qdsat = T0 - (T4/T6)
	# qdsat for internal Rds, Ref: BSIM4
	else 
		# Accurate qdsat derived from consistent I-V
		T11 = Weff * 2.0 * nq * Cox * nVt * VSAT_a
		T12 = T11 * LambdaC * Rdss / (2.0 * nVt)
		T0  = 0.5 * LambdaC * (qs * qs + qs) / (1.0 + 0.5 * LambdaC * (1.0 + qs))
		T1  = 2.0 * LambdaC * (qs - T0)
		T2  = sqrt(1.0 + T1 * T1)
		ln_T1_T2 = asinh(T1)
		if (T1 != 0.0) 
			T3 = T2 + (1.0 / T1) * ln_T1_T2
		else 
			T3 = T2 + (1.0 / T2)
		end
		T4 = T0 * T3 + T12 * T0 * (qs + T0 + 1.0) - LambdaC * ((qs * qs + qs) - (T0 * T0 + T0))
		if (T1 != 0.0) 
			T5 = -2.0 * LambdaC * (T1 * T2 - ln_T1_T2) / (T1 * T1)
		else 
			T5 = -2.0 * LambdaC * (T1 / T2) * (T1 / T2) * (T1 / T2)
		end
		T6 = T0 * T5 + T3 + T12 * (qs + 2.0 * T0 + 1.0) + LambdaC * (2.0 * T0 + 1.0)
		T0 = T0 - T4 / T6
		T1 = 2.0 * LambdaC * (qs - T0)
		T2 = sqrt(1.0 + T1 * T1)
		ln_T1_T2 = asinh(T1)
		if (T1 != 0) 
			T3 = T2 + (1.0 / T1) * ln_T1_T2
		else 
			T3 = T2 + (1.0 / T2)
		end
		T4 = T0 * T3 + T12 * T0 * (qs + T0 + 1.0) - LambdaC * ((qs * qs + qs) - (T0 * T0 + T0))
		if (T1 != 0.0) 
			T5 = -2.0 * LambdaC * (T1 * T2 - ln_T1_T2) / (T1 * T1)
		else 
			T5 = -2.0 * LambdaC * (T1 / T2) * (T1 / T2) * (T1 / T2)
		end
		T6    = T0 * T5 + T3 + T12 * (qs + 2.0 * T0 + 1.0) + LambdaC * (2.0 * T0 + 1.0)
		qdsat = T0 - T4 / T6
	end
	vdsat = psip - 2.0 * phib_n - (2.0 * qdsat + lln((qdsat * 2.0 * nq * inv_gam) * ((qdsat * 2.0 * nq * inv_gam) + (gam / (nq - 1.0)))))
	Vdsat = vdsat * nVt
	
	# Normalized charge qdeff at drain end of channel
	# Vdssat clamped to avoid negative values during transient simulation
	Smooth(Vdsat - Vs, 0.0, 1.0e-3, Vdssat)
	T7      = pow(Vds / Vdssat , 1.0 / DELTA_t)
	T8      = pow(1.0 + T7, -DELTA_t)
	Vdseff  = Vds * T8
	vdeff   = (Vdseff + Vs) * inv_nVt
	BSIM_q(psip, phib_n, vdeff, gam, qdeff)
	
	# Reevaluation of nq to include qdeff
	psiavg = psip - qs - qdeff -1.0
	Smooth(psiavg, 1.0, 2.0, T0)
	T2 = sqrt(T0)
	nq = 1.0 + gam / (sqrtpsip + T2)
	
	# Inversion and bulk charge
	DQSD2 = (qs - qdeff) * (qs - qdeff)
	T0    = 1.0 / (1.0 + qs + qdeff)
	T1    = DQSD2 * T0
	Qb    = vgfb - psip - (nq - 1.0) * (qs + qdeff + Oneby3 * T1)
	T2    = Oneby3 * nq
	T3    = T1 * T0
	Qs    = T2 * (2.0 * qs + qdeff + 0.5 * (1.0 + 0.8 * qs + 1.2 * qdeff) * T3)
	Qd    = T2 * (qs + 2.0 * qdeff + 0.5 * (1.0 + 1.2 * qs + 0.8 * qdeff) * T3)
	
	# Mobility degradation, Ref: BSIM4
	# Average charges (qba and qia) - normalized to Cox
	Smooth(nVt * Qb, 0.0, 0.1, qba)
	qia   = nVt * (Qs + Qd)
	
	Eeffm = EeffFactor * (qba + eta_mu * qia)
	T2    = pow(0.5 * (1.0 + (qia / qba)), UCS_a)
	T3    = (UA_a + UC_a * Vbsx) * pow(Eeffm, EU_t) + UD_a / T2
	T4    = 1.0 + T3
	Smooth(T4, 1.0, 0.0015, Dmob)
	
	# Output conductance
	Esat  = 2.0 * VSAT_a / (U0_a / Dmob)
	EsatL = Esat * Leff
	if (PVAG_i > 0.0) 
		PVAGfactor = 1.0 + PVAG_i * qia / EsatL
	else 
		PVAGfactor = 1.0 / (1.0 - PVAG_i * qia / EsatL)
	end
	
	# Output conductance due to DIBL, Ref: BSIM4
	DIBLfactor = PDIBLC_a
	diffVds    = Vds - Vdseff
	Vgst2Vtm   = qia + 2.0 * nVt
	if (DIBLfactor > 0.0) 
		T3     = Vgst2Vtm / (Vdssat + Vgst2Vtm)
		T4     = hypsmooth((1.0 + PDIBLCB_i * Vbsx), 1.0e-3)
		T5     = 1.0 / T4
		VaDIBL = Vgst2Vtm / DIBLfactor * T3 * PVAGfactor * T5
		Moc    = 1.0 + diffVds / VaDIBL
	else 
		Moc = 1.0
	end
	
	# Degradation factor due to pocket implants, Ref: BSIM4
	if (FPROUT_i <= 0.0) 
		Fp = 1.0
	else 
		T9 = FPROUT_i * sqrt(Leff) / Vgst2Vtm
		Fp = 1.0 / (1.0 + T9)
	end
	
	# Channel length modulation, Ref: BSIM4
	Vasat = Vdssat + EsatL
	if(PCLM_a != 0.0) 
		if (PCLMG < 0.0) 
			T1 = PCLM_a / (1.0 - PCLMG * qia / EsatL) / Fp
		else 
			T1 = PCLM_a * (1.0 + PCLMG * qia / EsatL) / Fp
		end
		MdL = 1.0 + T1 * lln(1.0 + diffVds / T1 / Vasat)
	else 
		MdL = 1.0
	end
	Moc = Moc * MdL
	
	# Calculate Va_DITS, Ref: BSIM4
	T1 = lexp(PDITSD_i * Vds)
	if (PDITS_i > 0.0) 
		T2      = 1.0 + PDITSL * Leff
		VaDITS  = (1.0 + T2 * T1) / PDITS_i
		VaDITS  = VaDITS * Fp
	else 
		VaDITS  = MAX_EXPL
	end
	T4  = diffVds / VaDITS
	T0  = 1.0 + T4
	Moc = Moc * T0
	
	# Calculate Va_SCBE, Ref: BSIM4
	if (PSCBE2_i > 0.0) 
		if (diffVds > PSCBE1_i * litl / EXPL_THRESHOLD) 
			T0     = PSCBE1_i * litl / diffVds
			VaSCBE = Leff * lexp(T0) / PSCBE2_i
		else 
			VaSCBE = MAX_EXPL * Leff/PSCBE2_i
		end
	else 
		VaSCBE = MAX_EXPL
	end
	Mscbe = 1.0 + (diffVds / VaSCBE)
	Moc   = Moc * Mscbe
	
	# Velocity saturation
	T0 = pow(Dmob, 1.0 / PSAT_a)
	T11 = PSATB_i * Vbsx
	T12 = sqrt(0.1+T11*T11)
	T1  = 0.5*(1-T11+sqrt((1-T11)*(1-T11)+T12))
	T2  = 10.0 * PSATX * qia * T1 / (10.0 * PSATX + qia * T1)
	if (PTWG_a < 0.0) 
		LambdaC = 2.0 * ((U0_a / T0) * nVt / (VSAT_a * Leff)) * (1.0 / (1.0 - PTWG_a * T2))
	else 
		LambdaC = 2.0 * ((U0_a / T0) * nVt / (VSAT_a * Leff)) * (1.0 + PTWG_a * T2)
	end
	T1 = 2.0 * LambdaC * (qs - qdeff)
	T2 = sqrt(1.0 + T1 * T1)
	if (T1 != 0.0) 
		Dvsat = 0.5 * (T2 + (1.0 / T1) * asinh(T1))
	else 
		Dvsat = 0.5 * (T2 + (1.0 / T2))
	end
	Dptwg = Dvsat
	
	# S/D series resistances, Ref: BSIM4
	Rsource = 0.0
	Rdrain  = 0.0
	if (RDSMOD == 1) 
		Rdsi = 0.0
		Dr   = 1.0
		# Rs (Source side resistance for all fingers)
		T2      = Vgs_noswap - Vfbsdr
		T3      = sqrt(T2 * T2 + 0.01)
		Vgs_eff = 0.5 * (T2 + T3)
		T5      = 1.0 + PRWG_i * Vgs_eff
		T6      = (1.0 / T5) + PRWB_i * Vsb_noswap
		T4      = 0.5 * (T6 + sqrt(T6 * T6 + 0.01))
		Rsource = rdstemp * (RSourceGeo + (RSWMIN_i + RSW_i * T4) * WeffWRFactor)
		# Rd (Drain side resistance for all fingers)
		T2      = Vgd_noswap - Vfbsdr
		T3      = sqrt(T2 * T2 + 0.01)
		Vgd_eff = 0.5 * (T2 + T3)
		T5      = 1.0 + PRWG_i * Vgd_eff
		T6      = (1.0 / T5) + PRWB_i * Vdb_noswap
		T4      = 0.5 * (T6 + sqrt(T6 * T6 + 0.01))
		Rdrain  = rdstemp * (RDrainGeo + (RDWMIN_i + RDW_i * T4) * WeffWRFactor)
	else 
		# Ref: (1) BSIM4 (2) "Operation and Modeling of the MOS Transistor" by Yannis Tsividis
		T0      = 1.0 + PRWG_i * qia
		T1      = PRWB_i * (sqrtPhistVbs - sqrtPhist)
		T2      = 1.0 / T0 + T1
		T3      = 0.5 * (T2 + sqrt(T2 * T2 + 0.01))
		Rdsi    = rdstemp * (RDSWMIN_i + RDSW_i * T3) * WeffWRFactor * NF
		Rdrain  = RDrainGeo
		Rsource = RSourceGeo
		Dr      = 1.0 + U0_a /(Dvsat * Dmob) * Cox * Weff / Leff * qia * Rdsi
		if (RDSMOD == 2) 
			Rdsi    = rdstemp * (RSourceGeo + (RDSWMIN_i + RDSW_i * T3) * WeffWRFactor * NF + RDrainGeo)
			Rdrain  = 0.0
			Rsource = 0.0
			Dr      = 1.0 + U0_a /(Dvsat * Dmob) * Cox * Weff / Leff * qia * Rdsi
		end
	end
	
	# Non-saturation effect
	T0   = A1_t + A2_t / (qia + 2.0 * n * Vtm)
	DQSD = qs - qdeff
	T1   = T0 * DQSD * DQSD
	T2   = T1 + 1.0 - 0.001
	T3   = -1.0 + 0.5 * (T2 + sqrt(T2 * T2 + 0.004))
	Nsat = 0.5 * (1.0 + sqrt(1.0 + T3))
	
	# MNUD model to enhance Id-Vd fitting flexibility
	T0   = (qs + qdeff)
	T1   = (qs - qdeff)
	T2   = T1 / (T0 + M0_t)
	T3   = K0_t * T2 * T2
	Mnud = 1.0 + T3
	
	# MNUD1 to enhance the fitting flexiblity for the gm/Id - similar approach used in BSIM-CMG
	T9    = C0_t / (max(0, C0SI_t + C0SISAT_t * T1 * T1) * T0 + 2.0 * n * Vtm)
	Mnud1 = lexp(-T9)
	Dtot  = Dmob * Dvsat * Dr
	
	# Effective mobility including mobility degradation
	ueff = U0_a / Dtot
	
	# I-V
	ids  = 2.0 * NF * nq * ueff * Weff / Leff * Cox * nVt * nVt * ((qs - qdeff) * (1.0 + qs + qdeff)) * Moc / Nsat * Mnud * Mnud1
	ids  = ids * IDS0MULT
	
	# High-voltage model s: Ref. - Harshit Agarwal et.al., IEEE TED vol. 66, issue 10, pp. 4258, 2019
	
	if (RDSMOD == 1 && HVMOD == 1) 
		T4  = 1 + PDRWB * Vbsx
		T0  = ids 
		T11 = NF * Weff * q  * VDRIFT_t 
		if (RDLCW != 0) 
			idrift_sat_d = T11 * NDRIFTD  
			delta_hv = pow(ids,4-MDRIFT) / (pow(ids,4-MDRIFT)+ HVFACTOR * pow(idrift_sat_d,4-MDRIFT))
			T5  = T0/idrift_sat_d 
			if (T5 >= 0.99) 
				T5  = 0.5 * ((T5 + 0.99) - sqrt( (T5 - 0.99) * (T5 - 0.99) + 1.0e-6) + 0.001 )
			end
			T0D = delta_hv * pow(T5,MDRIFT) 
			T1D = 1.0 - T0D 
			T2D = pow(T1D,1.0/MDRIFT) 
			rdrift_d = rdstemphv * RDLCW * WeffWRFactor/T2D * T4
			IDRIFTSATD = idrift_sat_d
			if (rdrift_d < 0) 
				rdrift_d = 0
			end
		end
	
		if (RSLCW != 0) 
			idrift_sat_s = T11 * NDRIFTS  
			delta_hv = pow(ids,4-MDRIFT) / (pow(ids,4-MDRIFT)+ HVFACTOR * pow(idrift_sat_s,4-MDRIFT))
			T5  = T0/idrift_sat_s 
			if (T5 >= 0.99) 
				T5  = 0.5 * ((T5 + 0.99) - sqrt( (T5 - 0.99) * (T5 - 0.99) + 1.0e-6) + 0.001 )
			end
			T0S = delta_hv * pow(T5,MDRIFT) 
			T1S = 1.0 - T0S 
			T2S = pow(T1S,1.0/MDRIFT) 
			rdrift_s = rdstemphv * RSLCW * WeffWRFactor/T2S * T4
			if (rdrift_s < 0) 
				rdrift_s = 0
			end
		end
	
		Rdrain  = Rdrain + rdrift_d  
		Rsource = Rsource + rdrift_s  
	end
	
	QIOV  = 0
	QBOV  = 0
	QIOVS = 0
	QBOVS = 0
	
	# CV calculations for HVMOD
	if (RDSMOD == 1 && HVCAP == 1 && HVMOD == 1) 
		vgfbdrift = -devsign * V(gm,di) - VFBOV 
		vgfbdrift = vgfbdrift/Vt
		gamhv     = sqrt(2.0 * q * epssi * NDR * inv_Vt) / Cox
		phibHV    = lln(NDR / ni)
		PO_psip(vgfbdrift,gamhv,0,phibHV,psip_k)
		BSIM_q(psip_k, phibHV, devsign *V(di,bi)/Vt, gamhv, q_k)
	
		# calculate nq for the drift region
		Smooth(psip_k, 1.0, 2.0, psipclamp_hv)
		sqrtpsip_k = sqrt(psipclamp_hv)
		psiavg_hv = psip_k - 2.0 * q_k
		Smooth(psiavg_hv, 1.0, 2.0, T0)
		nq_hv = 1.0 + gamhv / (sqrtpsip_k + sqrt(T0))
		psi_k = psip_k - 2 * q_k
	
		# contribution due to accumulation of the overlap region
		QBOV = NF * Wact * LOVER * EPS0 * EPSROX / BSIMBULKTOXP * Vt * (vgfbdrift - psi_k - 2 * nq_hv * q_k)
	
		# contribution due to inversion of the overlap region
		if (SLHV > 0) 
			T1 = 1 + q_k / SLHV1 
			T2 = SLHV * 1.9e-9 / T1
			T0 = 3.9 * EPS0 / (BSIMBULKTOXP * 3.9 / EPSROX + T2 / epsratio)
		else 
			T0 = EPS0 * EPSROX / BSIMBULKTOXP
		end
		QIOV = NF * Wact * LOVERACC * 2 * nq_hv * Vt * T0 * q_k 
	
		# For symmetric device, adding contribution of the source side drift region
	
		if (HVCAPS == 1) 
			vgfbdrift = -devsign * V(gm,si) - VFBOV 
			vgfbdrift = vgfbdrift/Vt
			PO_psip(vgfbdrift,gamhv,0,phibHV,psip_k)
			BSIM_q(psip_k, phibHV, devsign * V(si,bi)/Vt, gamhv, q_k)           
	
			Smooth(psip_k, 1.0, 2.0, psipclamp_hv)
			sqrtpsip_k = sqrt(psipclamp_hv)
			psiavg_hv = psip_k - 2.0 * q_k
			Smooth(psiavg_hv, 1.0, 2.0, T0)
			nq_hv = 1.0 + gamhv / (sqrtpsip_k + sqrt(T0))
			psi_k = psip_k - 2 * q_k
	
			QBOVS = NF * Wact * LOVER * EPS0 * EPSROX / BSIMBULKTOXP * Vt * (vgfbdrift - psi_k - 2 * nq_hv * q_k)
	
			if (SLHV > 0) 
				T1 = 1 + q_k / SLHV1 
				T2 = SLHV * 1.9e-9 / T1
				T0 = 3.9 * EPS0 / (BSIMBULKTOXP * 3.9 / EPSROX + T2 / epsratio)
			else 
				T0 = EPS0 * EPSROX / BSIMBULKTOXP
			end
	
		   QIOVS = NF * Wact * LOVERACC * 2 * nq_hv * Vt * T0 * q_k 
		end
	end 
	
	
	Gcrg = 0.0
	if (RGATEMOD > 1) 
		idsovvds = ueff * Weff / Leff * Cox * qia
		T9       = XRCRG2 * Vt
		T0       = T9 * ueff * Weff / Leff * Cox
		Gcrg     = XRCRG1 * NF * (T0 + idsovvds)
		if (RGATEMOD == 2) 
			T11  = Grgeltd + Gcrg
			Gcrg = Grgeltd * Gcrg / T11
		end
	end
	
	# Impact ionization currents, Ref: BSIM4
	if ((ALPHA0_i <= 0.0) || (BETA0_t <= 0.0)) 
		Iii = 0.0
	elseif (diffVds > BETA0_t / EXPL_THRESHOLD) 
		T1  = -BETA0_t / diffVds
		Iii = ALPHA0_i * diffVds * ids * lexp(T1) / Mscbe
	else 
		Iii = ALPHA0_i * diffVds * ids * MIN_EXPL / Mscbe
	end
	
	# Secondary impact ionization in the drift region
	if (HVMOD == 1 && IIMOD == 1) 
		Ntot = DRII1 * ids/(NF * Weff * q  * VDRIFT_t )   
		Nextra = Ntot/NDRIFTD - 1
		Smooth(Nextra, 0, DELTAII, Nextra)
		Nextra = NDRIFTD * Nextra
	
		Smooth(devsign * V(d,bi) - Vdseff - DRII2, 0, 0.05, T2)
		T3 = 2.0 * q /(EPSRSUB * EPS0) * Nextra
		T3 = T3 * T2                             
	
		if (T3 > BETADR / EXPL_THRESHOLD) 
			T1  = -BETADR/T3
			IsubDR = ALPHADR * T2 * ids * lexp(T1)
		else 
			IsubDR = ALPHADR * T2 * ids * MIN_EXPL
		end
		Iii = Iii + IsubDR
	end
	
	ISUB = Iii * devsign
	
	# Gate currents, Ref: BSIM4
	igbinv = 0.0 igbacc = 0.0 igb = 0.0 igcs = 0.0
	igcd   = 0.0 igs    = 0.0 igd = 0.0
	
	if ((IGCMOD != 0) || (IGBMOD != 0)) 
		Voxm    = nVt * (vgfb - psip + qs + qdeff)
		T1      = sqrt(Voxm * Voxm + 1.0e-4)
		Voxmacc = 0.5 * (-Voxm + T1)
		Voxminv = 0.5 * (Voxm + T1)
	# Igbinv
	if (IGBMOD != 0) 
		T1     = Voxm / NIGBACC_i / Vt
		Vaux_Igbacc = NIGBACC_i * Vt * lln(1.0 + lexp(-T1))
		T2     = AIGBACC_i - BIGBACC_i * Voxmacc
		T3     = 1.0 + CIGBACC_i * Voxmacc
		T4     = -7.45669e11 * TOXE * T2 * T3
		T5     = lexp(T4)
		T6     = 4.97232e-7
		igbacc = NF * Weff * Leff * T6 * ToxRatio * Vg * Vaux_Igbacc * T5
		igbacc = igbacc * igtemp
		T1     = (Voxm - EIGBINV_i) / NIGBINV_i / Vt
		Vaux_Igbinv = NIGBINV_i * Vt * lln(1.0 + lexp(T1))
		T2     = AIGBINV_i - BIGBINV_i * Voxminv
		T3     = 1.0 + CIGBINV_i * Voxminv
		T4     = -9.82222e11 * TOXE * T2 * T3
		T5     = lexp (T4)
		T6     = 3.75956e-7
		igbinv = NF * Weff * Leff * T6 * ToxRatio * Vg * Vaux_Igbinv * T5
		igbinv = igbinv * igtemp
		igb    = igbacc + igbinv
	end
	
	if (IGCMOD != 0) 
		# Igcinv
		T1   = AIGC_i - BIGC_i * Voxminv
		T2   = 1.0 + CIGC_i * Voxminv
		T3   = Bechvb * T1 * T2
		T4   = nq * nVt * (qs + qdeff) * lexp(T3)
		igc0 = NF * Aechvb * T4 * (Vg + 0.5 * Vdsx - 0.5 * (Vs + Vd)) * igtemp
		# Gate-current partitioning
		Vdseffx = sqrt(Vdseff * Vdseff + 0.01) - 0.1
		T1      = PIGCD_i * Vdseffx
		T1_exp  = lexp(-T1)
		T3      = T1 + T1_exp -1.0 + 1.0e-4
		T4      = 1.0 - (T1 + 1.0) * T1_exp + 1.0e-4
		T5      = T1 * T1 + 2.0e-4
		if (sigvds > 0) 
			igcd = igc0 * T4 / T5
			igcs = igc0 * T3 / T5
		else 
			igcs = igc0 * T4 / T5
			igcd = igc0 * T3 / T5
		end
		# Igs
		T2      = Vgs_noswap - Vfbsdr
		Vgs_eff = sqrt(T2 * T2 + 1.0e-4)
		if (IGCLAMP == 1) 
			T1 = hypsmooth((AIGS_i - BIGS_i * Vgs_eff), 1.0e-6)
			if (CIGS_i < 0.01) 
				CIGS_i = 0.01
			end
		else 
			T1 = AIGS_i - BIGS_i * Vgs_eff
		end
		T2       = 1.0 + CIGS_i * Vgs_eff
		T3       = BechvbEdge * T1 * T2
		T4       = lexp(T3)
		igs_mult = igtemp * NF * AechvbEdge * DLCIG_i
		igs      = igs_mult * Vgs_noswap * Vgs_eff * T4
		# Igd
		T2      = Vgd_noswap - Vfbsdr
		Vgd_eff = sqrt(T2 * T2 + 1.0e-4)
		if (IGCLAMP == 1) 
			T1 = hypsmooth((AIGD_i - BIGD_i * Vgd_eff), 1.0e-6)
			if (CIGD_i < 0.01) 
				CIGD_i = 0.01
			end
		else 
			T1 = AIGD_i - BIGD_i * Vgd_eff
		end
		T2       = 1.0 + CIGD_i * Vgd_eff
		T3       = BechvbEdge * T1 * T2
		T4       = lexp(T3)
		igd_mult = igtemp * NF * AechvbEdge * DLCIGD_i
		igd      = igd_mult * Vgd_noswap * Vgd_eff * T4
	end
	end
	IGS  = devsign * igs
	IGD  = devsign * igd
	IGB  = devsign * igb
	IGCS = devsign * igcs
	IGCD = devsign * igcd
	
	# GIDL and GISL currents, Ref: BSIM4
	igisl = 0.0
	igidl = 0.0
	if (GIDLMOD != 0) 
		T0 = epsratio * TOXE
		# GIDL
		if ((AGIDL_i <= 0.0) || (BGIDL_t <= 0.0) || (CGIDL_i < 0.0)) 
			T6 = 0.0
		else 
			T1 = (-Vgd_noswap - EGIDL_i + Vfbsdr) / T0
			T1 = hypsmooth(T1, 1.0e-2)
			T2 = BGIDL_t / (T1 + 1.0e-3)
			if (CGIDL_i != 0.0) 
				T3 = Vdb_noswap * Vdb_noswap * Vdb_noswap
				T4 = CGIDL_i + abs(T3) + 1.0e-4
				T5 = hypsmooth(T3 / T4, 1.0e-6) - 1.0e-6
			else 
				T5 = 1.0
			end
			T6 = AGIDL_i * Weff * T1 * lexp(-T2) * T5
		end
		igidl = T6
		# GISL
		if ((AGISL_i <= 0.0) || (BGISL_t <= 0.0) || (CGISL_i < 0.0)) 
			T6 = 0.0
		else 
			T1 = (-Vgs_noswap - EGISL_i + Vfbsdr) / T0
			T1 = hypsmooth(T1, 1.0e-2)
			T2 = BGISL_t / (T1 + 1.0e-3)
			if (CGISL_i != 0.0) 
				T3 = Vsb_noswap * Vsb_noswap * Vsb_noswap
				T4 = CGISL_i + abs(T3) + 1.0e-4
				T5 = hypsmooth(T3 / T4, 1.0e-6) - 1.0e-6
			else 
				T5 = 1.0
			end
			T6 = AGISL_i * Weff * T1 * lexp(-T2) * T5
		end
		igisl = T6
	end
	IGIDL = devsign * NF * igidl
	IGISL = devsign * NF * igisl
	
	# Junction currents and capacitances
	# Source-side junction currents
	if (Isbs > 0.0) 
		if (Vbs_jct < VjsmRev) 
			T0  = Vbs_jct / Nvtms
			T1  = lexp(T0) - 1.0
			T2  = IVjsmRev + SslpRev * (Vbs_jct - VjsmRev)
			Ibs = T1 * T2
		elseif (Vbs_jct <= VjsmFwd) 
			T0  = Vbs_jct / Nvtms
			T1  = (BVS + Vbs_jct) / Nvtms
			T2  = lexp(-T1)
			Ibs = Isbs * (lexp(T0) + XExpBVS - 1.0 - XJBVS * T2)
		else 
			Ibs = IVjsmFwd + SslpFwd * (Vbs_jct - VjsmFwd)
		end
	else 
		Ibs = 0.0
	end
	
	# Source-side junction tunneling currents
	if (JTSS_t > 0.0) 
		if ((VTSS - Vbs_jct) < (VTSS * 1.0e-3)) 
			T0  = -Vbs_jct / Vtm0 / NJTS_t
			T1  = lexp(T0 * 1.0e3) - 1.0
			Ibs = Ibs - ASeff * JTSS_t * T1
		else 
			T0  = -Vbs_jct / Vtm0 / NJTS_t
			T1  = lexp(T0 * VTSS / (VTSS - Vbs_jct)) - 1.0
			Ibs = Ibs - ASeff * JTSS_t * T1
		end
	end
	if (JTSSWS_t > 0.0) 
		if ((VTSSWS - Vbs_jct) < (VTSSWS * 1.0e-3)) 
			T0  = -Vbs_jct / Vtm0 / NJTSSW_t
			T1  = lexp(T0 * 1.0e3) - 1.0
			Ibs = Ibs - PSeff * JTSSWS_t * T1
		else 
			T0  = -Vbs_jct / Vtm0 / NJTSSW_t
			T1  = lexp(T0 * VTSSWS / (VTSSWS - Vbs_jct)) - 1.0
			Ibs = Ibs - PSeff * JTSSWS_t * T1
		end
	end
	if (JTSSWGS_t > 0.0) 
		if((VTSSWGS - Vbs_jct) < (VTSSWGS * 1.0e-3)) 
			T0  = -Vbs_jct / Vtm0 / NJTSSWG_t
			T1  = lexp(T0 * 1.0e3) - 1.0
			Ibs = Ibs - Weffcj * NF * JTSSWGS_t * T1
		else 
			T0  = -Vbs_jct / Vtm0 / NJTSSWG_t
			T1  = lexp(T0 * VTSSWGS / (VTSSWGS - Vbs_jct)) - 1.0
			Ibs = Ibs - Weffcj * NF * JTSSWGS_t * T1
		end
	end
	
	# Drain-side junction currents
	if (Isbd > 0.0) 
		if (Vbd_jct < VjdmRev) 
			T0  = Vbd_jct / Nvtmd
			T1  = lexp(T0) - 1.0
			T2  = IVjdmRev + DslpRev * (Vbd_jct - VjdmRev)
			Ibd = T1 * T2
		elseif (Vbd_jct <= VjdmFwd) 
			T0  = Vbd_jct / Nvtmd
			T1  = (BVD + Vbd_jct) / Nvtmd
			T2  = lexp(-T1)
			Ibd = Isbd * (lexp(T0) + XExpBVD - 1.0 - XJBVD * T2)
		else 
			Ibd = IVjdmFwd + DslpFwd * (Vbd_jct - VjdmFwd)
		end
	else 
		Ibd = 0.0
	end
	
	# Drain-side junction tunneling currents
	if (JTSD_t > 0.0) 
		if ((VTSD - Vbd_jct) < (VTSD * 1.0e-3)) 
			T0  = -Vbd_jct / Vtm0 / NJTSD_t
			T1  = lexp(T0 * 1.0e3) - 1.0
			Ibd = Ibd - ADeff * JTSD_t * T1
		else 
			T0  = -Vbd_jct / Vtm0 / NJTSD_t
			T1  = lexp(T0 * VTSD/ (VTSD - Vbd_jct)) - 1.0
			Ibd = Ibd - ADeff * JTSD_t * T1
		end
	end
	if (JTSSWD_t > 0.0) 
		if ((VTSSWD - Vbd_jct) < (VTSSWD * 1.0e-3)) 
			T0  = -Vbd_jct / Vtm0 / NJTSSWD_t
			T1  = lexp(T0 * 1.0e3) - 1.0
			Ibd = Ibd - PDeff * JTSSWD_t * T1
		else 
			T0  = -Vbd_jct / Vtm0 / NJTSSWD_t
			T1  = lexp(T0 * VTSSWD / (VTSSWD - Vbd_jct)) - 1.0
			Ibd = Ibd - PDeff * JTSSWD_t * T1
		end
	end
	if (JTSSWGD_t > 0.0) 
		if ((VTSSWGD - Vbd_jct) < (VTSSWGD * 1.0e-3)) 
			T0  = -Vbd_jct / Vtm0 / NJTSSWGD_t
			T1  = lexp(T0 * 1.0e3) - 1.0
			Ibd = Ibd - Weffcj * NF * JTSSWGD_t * T1
		else 
			T0  = -Vbd_jct / Vtm0 / NJTSSWGD_t
			T1  = lexp(T0 * VTSSWGD / (VTSSWGD - Vbd_jct)) - 1.0
			Ibd = Ibd - Weffcj * NF * JTSSWGD_t * T1
		end
	end
	
	# Junction capacitances (no swapping)
	# Source-to-bulk junction
	Czbs       = CJS_t * ASeff
	Czbssw     = CJSWS_t * PSeff
	Czbsswg    = CJSWGS_t * Weffcj * NF
	czbs_p1    = pow(0.1, -MJS)
	czbs_p2    = 1.0 / (1.0 - MJS) * (1.0 - 0.05 * MJS * (1.0 + MJS) * czbs_p1)
	czbssw_p1  = pow(0.1, -MJSWS)
	czbssw_p2  = 1.0 / (1.0 - MJSWS) * (1.0 - 0.05 * MJSWS * (1.0 + MJSWS) * czbssw_p1)
	czbsswg_p1 = pow(0.1, -MJSWGS)
	czbsswg_p2 = 1.0 / (1.0 - MJSWGS) * (1.0 - 0.05 * MJSWGS * (1.0 + MJSWGS) * czbsswg_p1)
	JunCap(Czbs, Vbs_jct, PBS_t, MJS, czbs_p1, czbs_p2, Qbsj1)
	JunCap(Czbssw, Vbs_jct, PBSWS_t, MJSWS, czbssw_p1, czbssw_p2, Qbsj2)
	JunCap(Czbsswg, Vbs_jct, PBSWGS_t, MJSWGS, czbsswg_p1, czbsswg_p2, Qbsj3)
	Qbsj = Qbsj1 + Qbsj2 + Qbsj3
	
	# Drain-to-bulk junction
	Czbd       = CJD_t * ADeff
	Czbdsw     = CJSWD_t * PDeff
	Czbdswg    = CJSWGD_t * Weffcj * NF
	czbd_p1    = pow(0.1, -MJD)
	czbd_p2    = 1.0 / (1.0 - MJD) * (1.0 - 0.05 * MJD * (1.0 + MJD) * czbd_p1)
	czbdsw_p1  = pow(0.1, -MJSWD)
	czbdsw_p2  = 1.0 / (1.0 - MJSWD) * (1.0 - 0.05 * MJSWD * (1.0 + MJSWD) * czbdsw_p1)
	czbdswg_p1 = pow(0.1, -MJSWGD)
	czbdswg_p2 = 1.0 / (1.0 - MJSWGD) * (1.0 - 0.05 * MJSWGD * (1.0 + MJSWGD) * czbdswg_p1)
	JunCap(Czbd, Vbd_jct, PBD_t, MJD, czbd_p1, czbd_p2, Qbdj1)
	JunCap(Czbdsw, Vbd_jct, PBSWD_t, MJSWD, czbdsw_p1, czbdsw_p2, Qbdj2)
	JunCap(Czbdswg, Vbd_jct, PBSWGD_t, MJSWGD, czbdswg_p1, czbdswg_p2, Qbdj3)
	Qbdj = Qbdj1 + Qbdj2 + Qbdj3
	
	# Sub-surface leakage drain current
	if (SSLMOD != 0) 
		T1 = pow(NDEP_i / 1.0e23, SSLEXP1)
		T2 = pow(300.0 / DevTemp, SSLEXP2)
		T3 = (devsign*SSL5 * V(bi,si)) / Vt         
		SSL0_NT  = SSL0 * lexp(-T1 * T2)
		SSL1_NT  = SSL1 * T2 * T1
		PHIB_SSL = SSL3 * tanh(lexp(devsign * SSL4 * (V(gi, bi) - VTH - V(si,bi))))
		Issl     = sigvds * NF * Weff * SSL0_NT * lexp(T3) * lexp(-SSL1_NT * Leff) * lexp(PHIB_SSL / Vt) * (lexp(SSL2 * Vdsx / Vt) - 1.0)
		I(di, si) <+ devsign * Issl
	end
	
	# Harshit's new flicker noise model. Ref: H. Agarwal et. al., IEEE J-EDS, vol. 3, no. 4, April 2015.
	Nt      = 4.0 * Vt * q
	Esatnoi = 2.0 * VSAT_a / ueff
	if (EM <= 0.0) 
	   DelClm = 0.0
	else 
		T0     = (diffVds / litl + EM) / Esatnoi
		DelClm = litl * lln(T0)
		if (DelClm < 0.0) 
		   DelClm = 0.0
		end
	end
	Nstar = Vt / q * (Cox + Cdep + CIT_i)
	Nl    = 2.0 * nq * Cox * Vt * qdeff * Mnud1 * Mnud / q
	T0a   = q * q * q * Vt * abs(ids) * ueff
	T0b   = q * Vt * ids * ids
	T0c   = NOIA + NOIB * Nl + NOIC * Nl * Nl
	T0d   = (Nl + Nstar) * (Nl + Nstar)
	T0e   = NOIA * q * Vt
	if (FNOIMOD == 1) 
		LH1 = LH
		if (Leff > LH1) 
			T0 = (Leff - LH1)
		else 
			LH1 = Leff
			T0 = LH1
		end            
		if(LINTNOI >= T0 / 2.0) 
			println("Warning: LINTNOI = %e is too large - Leff for noise is negative. Re-setting LINTNOI = 0.", LINTNOI)
			LINTNOI_i = 0.0
		else 
			LINTNOI_i = LINTNOI
		end
		LeffnoiH = Leff
		vgfbh  = (Vg - VFB_i) / Vt
		gam_h  = sqrt(2.0 * q * epssi * HNDEP / Vt) / Cox
		phib_h = ln(HNDEP / ni)
	
		# Pinch-Off potential for halo region
		PO_psip(vgfbh, gam_h, 0.0, phib_h, psiph)
	
		# Normalized inversion charge at source end of halo MOSFET
		BSIM_q(psiph, phib_h, vs, gam_h, qsh)
		nq_h = 1.0 + gam_h / (2.0 * sqrt(psiph))
	
		# Setting mobility of halo region equal to the mobility of the channel. In general, U0H < ueff
		U0_i_h  = ueff
		beta_h  = U0_i_h * Cox * Weff
		beta_ch = ueff * Cox * Weff
	
		# Normalized drain current for halo transistor. Eq. (14) of the paper
		i1 = ids * LH1 / (2.0 * nq_h * beta_h * Vt * Vt)
	
		# Normalized drain current for channel transistor. Eq. (15) of the paper
		i2 = ids * (LeffnoiH - LH1) / (2.0 * nq * beta_ch * nVt * nVt)
		T0 = (1.0 + 4.0 * (qsh * qsh + qsh - i1))
		if (T0 < 1.0) 
			qdh = 0.0
		else 
			# Drain charge of halo transistor. Eq. (16) of the paper
			qdh = -0.5 + 0.5 * sqrt(T0)
		end
	
		# Source charge of channel transistor. Eq. (17) of the paper
		qsch   = -0.5 + 0.5 * sqrt(1.0 + 4.0 * (qdeff * qdeff + qdeff + i2))
		gds_h  = 2.0 * nq_h * beta_h * Vt * qdh
		gds_ch = 2.0 * nq * beta_ch * Vt * qdeff
		gm_ch  = 2.0 * beta_ch * Vt * (qsch - qdeff)
		R_ch   = gds_h * (LeffnoiH - LH1)
		R_h    = gm_ch * LH1 + gds_ch * LH1
		t_tot  = 1.0 / (R_ch + R_h) / (R_ch + R_h)
		CF_ch  = R_ch * R_ch * t_tot
		CF_h   = R_h * R_h * t_tot
	
		# Local noise source
		if (Leff != LH1) 
			Np2       = 2.0 * nq * Cox * Vt * qsch / q
			Leffnoi   = LeffnoiH - 2.0 * LINTNOI_i-LH1
			Leffnoisq = Leffnoi * Leffnoi
			# Channel transistor LNS
			T1     = 1.0e10 * Cox * Leffnoisq
			T2     = NOIA * lln((Np2 + Nstar) / (Nl + Nstar))
			T3     = NOIB * (Np2 - Nl)
			T4     = 0.5 * NOIC * (Np2 * Np2 - Nl * Nl)
			T5     = 1.0e10 * Leffnoisq * Weff * NF
			Ssi_ch = T0a / T1 * (T2 + T3 + T4) + T0b / T5 * DelClm * T0c / T0d
			T6     = Weff * NF * Leffnoi * 1.0e10 * Nstar * Nstar
			Swi_ch = T0e / T6 * ids * ids
			T7 = Swi_ch + Ssi_ch
			if (T7 > 0.0)  
				FNPowerAt1Hz_ch = (Ssi_ch * Swi_ch) / T7
			else 
				FNPowerAt1Hz_ch = 0.0
			end
		else 
			FNPowerAt1Hz_ch = 0.0
		end
	
		# Halo transistor LNS
		T8    = NOIA2 * q * Vt
		T9    = Weff * NF * LH1 * 1.0e10 * Nstar * Nstar
		Swi_h = T8 / T9 * ids * ids
		T10   = Swi_h
		if (T10 > 0.0) 
			FNPowerAt1Hz_h = Swi_h
		else 
			FNPowerAt1Hz_h = 0.0
		end
		# Overall noise
		FNPowerAt1Hz = FNPowerAt1Hz_ch * CF_ch + FNPowerAt1Hz_h * CF_h
		I(di, si) <+ flicker_noise(sigvds*FNPowerAt1Hz, EF, "1overf")
	else 
		# Parameter checking
		if (LINTNOI >= Leff/2.0) 
			println("Warning: LINTNOI = %e is too large - Leff for noise is negative. Re-setting LINTNOI = 0.", LINTNOI)
			LINTNOI_i = 0.0
		else 
			LINTNOI_i = LINTNOI
		end
		if (NOIA > 0.0 || NOIB > 0.0 || NOIC > 0.0) 
			Leffnoi   = Leff - 2.0 * LINTNOI_i
			Leffnoisq = Leffnoi * Leffnoi
			T0        = 1.0e10 * Cox * Leffnoisq
			N0        = 2.0 * nq * Cox * Vt * qs * Mnud1 * Mnud / q
			T1        = NOIA * lln((N0 + Nstar) / (Nl + Nstar))
			T2        = NOIB * (N0 - Nl)
			T3        = 0.5 * NOIC * (N0 * N0 - Nl * Nl)
			T4        = 1.0e10 * Leffnoisq * Weff * NF
			Ssi       = T0a / T0 * (T1 + T2 + T3) + T0b / T4 * DelClm * T0c / T0d
			T5        = Weff * NF * Leffnoi * 1.0e10 * Nstar * Nstar
			Swi       = T0e / T5 * ids * ids
			T6        = Swi + Ssi
			if (T6 > 0.0) 
				FNPowerAt1Hz = (Ssi * Swi) / T6/(1+NOIA1*pow((qs-qdeff),NOIAX))
			else 
				FNPowerAt1Hz = 0.0
			end
		else 
			FNPowerAt1Hz = 0.0
		end
		I(di, si) <+ flicker_noise(sigvds*FNPowerAt1Hz, EF, "1overf")
	end
	
	T0         = qia / Esatnoi / Leff
	T1         = T0 * T0
	T3         = RNOIA * (1.0 + TNOIA * Leff * T1)
	T4         = RNOIB * (1.0 + TNOIB * Leff * T1)
	T5         = RNOIK * (1.0 + TNOIK * Leff * T1)
	ctnoi      = RNOIC * (1.0 + TNOIC * Leff * T1)
	betanoisq  = 3.0 * T3 * T3
	betanoisq  = (betanoisq - 1.0) * exp(-Leff / LP) + 1.0
	betaLowId  = T5 * T5
	thetanoisq = T4 * T4
	cm_igid    = 0.0
	case (TNOIMOD)
		0: 
			QSi   = -NF * Weff * Leff * Cox * Vt * Qs
			QDi   = -NF * Weff * Leff * Cox * Vt * Qd
			T0    = ueff * abs(QSi + QDi)
			T1    = T0 * Rdsi + Leff * Leff
			Gtnoi = (T0 / T1) * NTNOI
			sidn  = Nt * Gtnoi
			I(di, si) <+ white_noise(sidn, "id")
			V(N1)     <+ 0.0
		end
		1: 
			Vtn   = 2.0 * nq * nVt
			T0    = ueff * Dptwg * Moc * Cox * Vtn
			T1    = 0.5 * (qs + qdeff)
			T3    = T1 + 0.5
			T4    = T3 * T3
			T5    = T4 * T3
			T6    = (qs - qdeff)
			T7    = T6 * T6
			T8    = T7 * T6
			T9    = (6.0 * T1 + 0.5) * T7
			Lvsat = Leff * Dptwg
			T10   = Lvsat / Leff
			T12   = 1.0 + (betaLowId * (Vdseff / Vdssat) / (TNOIK2 + qia))
			T12   = ((T12 - 1.0) * exp(-Leff / LP)) + 1.0
			Smooth(T12, 0.0, 1.0e-1, T12)
			mid   = T0 * NF * Weff / Lvsat * (T1 * T12 + T7 * betanoisq / (12.0 * T3))
			mig   = Lvsat * T10 * T10 * (T1 / T4 - T9 / (60.0 * T4 * T4) + T7 * T7 / (144.0 * T4 * T5)) * 15.0 / 4.0 * thetanoisq / (NF * Weff * 12.0 * T0)
			migid = T10 * (T6 / (12.0 * T3) - T8 / (144.0 * T5)) * ctnoi / 0.395
			sqid  = sqrt(Nt * mid)
			if (mig == 0.0) 
				sqig    = 0.0
				cm_igid = 0.0
			else 
				sqig = sqrt(Nt / mig)
				if (sqid == 0.0) 
					cm_igid = 0.0
				else 
					cm_igid = migid * sqig / sqid
				end
			end
			I(N2)     <+ white_noise(cm_igid, "corl")
			I(NI)     <+ white_noise(sqig * sqig * (1.0 - cm_igid), "corl")
			I(NI)     <+ -sqig * V(N2)
			I(NC)     <+ ddt(mig * Cox * Weff * NF * Leff * V(NC))
			I(di, si) <+ white_noise(sqid * sqid * (1.0 - cm_igid), "id")
			I(di, si) <+ sqid * V(N2)
			I(gi, si) <+ ddt(0.5 * ((1.0 + sigvds) * mig * Cox * Weff * NF * Leff * V(NC)))
			I(gi, di) <+ ddt(0.5 * ((1.0 - sigvds) * mig * Cox * Weff * NF * Leff * V(NC)))
		end
	endcase
	I(N2) <+ V(N2)
	I(NR) <+ V(NR)
	
	# Gate current shot noise
	if (IGCMOD != 0) 
		I(gi, si) <+ white_noise(2.0 * q * abs(igcs + igs), "igs")
		I(gi, di) <+ white_noise(2.0 * q * abs(igcd + igd), "igd")
	end
	if (IGBMOD != 0) 
		I(gi, bi) <+ white_noise(2.0 * q * abs(igb), "igb")
	end
	
	# C-V model
	vgfbCV   = vgfb
	gamg2    = (2.0 * q * epssi * NGATE_i) / (Cox * Cox * Vt)
	invgamg2 = 0.0
	if (CVMOD == 1) 
		VFBCV_i = VFBCV_i + DELVTO
		vg      = Vg * inv_Vt
		vs      = Vs * inv_Vt
		vfb     = VFBCV_i * inv_Vt
		vgfbCV  = vg - vfb
		phibCV    = lln(NDEPCV_i / ni)
		# Normalized body factor
		gamCV      = sqrt(2.0 * q * epssi * NDEPCV_i * inv_Vt) / Cox
		inv_gam  = 1.0 / gamCV
		gamg2    = (2.0 * q * epssi * NGATE_i) / (Cox * Cox * Vt)
		invgamg2 = (NGATE_i > 0.0) ? (1.0 / gamg2) : 0.0
		DPD      = (NGATE_i > 0.0) ? (NDEPCV_i / NGATE_i) : 0.0
	
		# psip: pinch-off voltage
		PO_psip(vgfbCV, gamCV, DPD, phibCV, psip)
	
		# Normalized inversion charge at source end of channel
		BSIM_q(psip, phibCV, vs, gamCV, qs)
		Smooth(psip, 1.0, 2.0, psipclamp)
		sqrtpsip = sqrt(psipclamp)
	
		# Source side surface potential
		psiavg = psip - 2.0 * qs
		Smooth(psiavg, 1.0, 2.0, T0)
		nq = 1.0 + gamCV / (sqrtpsip + sqrt(T0))
	
		# Drain saturation voltage
		T0 = Vt * (vgfbCV - psip - 2.0 * qs * (nq - 1.0))
		Smooth(T0, 0.0, 0.1, qbs)
	
		# Source side qi and qb for Vdsat (normalized to Cox)
		qis = 2.0 * nq * Vt * qs
		Eeffs = EeffFactor * (qbs + eta_mu * qis)
	
		# Ref: BSIM4 mobility model
		T3 = (UA_a + UC_a * Vbsx) * pow(Eeffs, EU_t)
		T4 = 1.0 + T3
		Smooth(T4, 1.0, 0.0015, Dmobs)
		LambdaC_by2 = (U0_a / Dmobs) * Vt / (VSATCV_t * Lact)
		qdsat       = LambdaC_by2 * (qs * qs + qs) / (1.0 + LambdaC_by2 * (1.0 + qs))
		vdsatcv     = psip - 2.0 * phibCV - (2.0 * qdsat + lln((qdsat * 2.0 * nq * inv_gam) * ((qdsat * 2.0 * nq * inv_gam) + (gam / (nq - 1.0)))))
		VdsatCV     = vdsatcv * Vt
		
		# Normalized charge qdeff at drain end of channel
		Smooth(VdsatCV - Vs, 0.0, 1e-3, VdssatCV)
		VdssatCV     = VdssatCV / ABULK
		T7     = pow(Vds / VdssatCV , 1.0 / DELTA_t)
		T8     = pow(1.0 + T7, -DELTA_t)
		Vdseff = Vds * T8
		vdeff  = (Vdseff + Vs) * inv_Vt
		BSIM_q(psip, phibCV, vdeff, gamCV, qdeff)
	
		# Reevaluation of nq to include qdeff needed for gummel symmetry
		psiavg = psip - qs - qdeff - 1.0
		Smooth(psiavg, 1.0, 2.0, T0)
		T2 = sqrt(T0)
		T3 = 1.0 + DPD + gamCV / (sqrtpsip + T2)
		T4 = 0.5 + DPD * T2 * inv_gam
		T5 = sqrt(T4 * T4 + T3 * (qs + qdeff) * invgamg2)
		nq = T3 / (T4 + T5)
	
		# C-V expressions including velocity saturation and CLM
		# Velocity saturation for C-V
		T0  = Vt * (vgfbCV - psip - 2.0 * qs * (nq - 1.0))
		Smooth(T0, 0.0, 0.1, qbs)
		T1  = Vt * (vgfbCV - psip - 2.0 * qdeff * (nq - 1.0))
		Smooth(T1, 0.0, 0.1, qbd)
		qb  = 0.5 * (qbs + qbd)
		qia = nq * Vt * (qs + qdeff)
		Eeffm = EeffFactor * (qb + eta_mu * qia)
	   PO_psip((vgfbCV + DELVFBACC * inv_Vt), gamCV, DPD, phibCV, psip)
		T3    = (UA_a + UC_a * Vbsx) * pow(Eeffm, EU_t)
		T4    = 1.0 + T3
		Smooth(T4, 1.0, 0.0015, Dmob)
		LambdaC = 2.0 * (U0_a / Dmob) * Vt / (VSATCV_t * Lact)
		dps     = qs - qdeff
		T1      = 2.0 * (LambdaC * dps) * (LambdaC * dps)
		zsat    = sqrt(1.0 + T1)
		Dvsat   = 0.5 * (1.0 + zsat)
		# CLM for C-V
		Esat    = 2.0 * VSATCV_t / (U0_a / Dmob)
		EsatL   = Esat * Lact
		Vasat   = VdssatCV + EsatL
		diffVds = Vds - Vdseff
	end
	if (PCLMCV_i != 0.0) 
		MdL = 1.0 + PCLMCV_i * lln(1.0 + diffVds / PCLMCV_i / Vasat)
	else 
		MdL = 1.0
	end
	MdL_2       = MdL * MdL
	inv_MdL     = 1.0 / MdL
	inv_MdL_2   = 1.0 / MdL_2
	MdL_less_1  = MdL - 1.0
	vgpqm = vgfbCV - psip
	DQSD  = (qs - qdeff)
	DQSD2 = (qs - qdeff) * (qs - qdeff)
	sis   = vgpqm + 2.0 * qs
	sid   = vgpqm + 2.0 * qdeff
	Smooth(sis, 0.0, 0.5, T1)
	Smooth(sid, 0.0, 0.5, T2)
	Temps = sqrt(0.25 + T1 * invgamg2)
	Tempd = sqrt(0.25 + T2 * invgamg2)
	T1 = sis / (1.0 + 2.0 * Temps)
	T2 = sid / (1.0 + 2.0 * Tempd)
	T3 = Temps + Tempd
	T4 = Oneby3 * (DQSD2 / (T3 * T3 * T3))
	T5 = (ABULK*Dvsat * inv_MdL) / (1.0 + qs + qdeff)
	T6 = 0.8 * (T3 * T3 + Temps * Tempd) * T5
	T7 = T6 + (2.0 * invgamg2)
	T8 = Oneby3 * DQSD2 * T5
	dqgeff = sid * (2.0 * Tempd - 1.0) / (2.0 * Tempd + 1.0)
	qbeff  = vgpqm - 2.0 * (nq - 1.0) * qdeff + dqgeff
	Qb  = inv_MdL * (T1 + T2 + (T4 * T7 - nq * (qs + qdeff + T8))) + MdL_less_1 * qbeff
	T9  = qs + qdeff
	T10 = DQSD2 * T5 * T5
	Qi  = nq * inv_MdL * (T9 + Oneby3 * DQSD2 * T5) + 2.0 * nq * MdL_less_1 * qdeff
	Qd1 = nq * inv_MdL_2 * (0.5 * T9 - (DQSD / 6.0) * (1.0 - DQSD * T5 - 0.2 * T10))
	Qd2 = nq * (MdL - inv_MdL) * qdeff
	Qd  = Qd1 + Qd2
	Qs  = Qi - Qd
	
	# Quantum mechanical effects
	Smooth(Vt * Qb, 0.0, 0.1, qbaCV)
	qiaCV      = Vt * (Qs + Qd)
	T0         = (qiaCV + ETAQM * qbaCV) / QM0
	T1         = 1.0 + pow(T0, 0.7 * BDOS)
	XDCinv     = ADOS * 1.9e-9 / T1
	Coxeffinv  = 3.9 * EPS0 / (BSIMBULKTOXP * 3.9 / EPSROX + XDCinv / epsratio)
	QBi        = -NF * Wact * Lact * (EPS0 * EPSROX / BSIMBULKTOXP) * Vt * Qb
	WLCOXVtinv = NF * Wact * Lact * Coxeffinv * Vt
	QSi        = -WLCOXVtinv * Qs
	QDi        = -WLCOXVtinv * Qd
	QGi        = -(QBi + QSi + QDi)
	
	# Outer fringing capacitances
	if (!$param_given(CF)) 
		CF_i = 2.0 * EPSROX * EPS0 / M_PI * lln(CFRCOEFF * (1.0 + 0.4e-6 / TOXE))
	end
	Cgsof = CGSO + CF_i
	Cgdof = CGDO + CF_i
	
	# Overlap capacitances
	if (COVMOD == 0) 
		Qovs = -Wact * NF * Cgsof * Vgs_ov_noswap
		Qovd = -Wact * NF * Cgdof * Vgd_ov_noswap
	else 
		T0    = sqrt((Vgs_ov_noswap - Vfbsdr + DELTA_1) * (Vgs_ov_noswap - Vfbsdr + DELTA_1) + 4.0 * DELTA_1)
		Vgsov = 0.5 * (Vgs_ov_noswap - Vfbsdr + DELTA_1 - T0)
		T1    = sqrt(1.0 - 4.0 * Vgsov / CKAPPAS_i)
		Qovs  = -Wact * NF * (Cgsof * Vgs_ov_noswap + CGSL_i * (Vgs_ov_noswap - Vfbsdr - Vgsov - 0.5 * CKAPPAS_i * (-1.0 + T1)))
		T0    = sqrt((Vgd_ov_noswap - Vfbsdr + DELTA_1) * (Vgd_ov_noswap - Vfbsdr + DELTA_1) + 4.0 * DELTA_1)
		Vgdov = 0.5 * (Vgd_ov_noswap - Vfbsdr + DELTA_1 - T0)
		T2    = sqrt(1.0 - 4.0 * Vgdov / CKAPPAD_i)
		Qovd  = -Wact * NF * (Cgdof * Vgd_ov_noswap + CGDL_i * (Vgd_ov_noswap - Vfbsdr - Vgdov - 0.5 * CKAPPAD_i * (-1.0 + T2)))
	end
	Qovb = -devsign * NF * Lact * CGBO * V(gm, bi)
	Qovg = -(Qovs + Qovd + Qovb)
		
	# Edge FET model
	if (EDGEFET == 1) 
		phib_edge     = lln(NDEPEDGE_i / ni)
		Phist         = max(0.4 + Vt * phib_edge + PHIN_i, 0.4)
		sqrtPhist     = sqrt(Phist)
		T1DEP         = sqrt(2.0 * epssi / (q * NDEPEDGE_i))
		NFACTOREDGE_t = NFACTOREDGE_i * hypsmooth((1.0 + TNFACTOREDGE_i * (TRatio - 1.0)), 1e-3)
		ETA0EDGE_t    = ETA0EDGE_i * (1.0 + TETA0EDGE_i * (TRatio - 1.0))
		Smooth(Phist - Vbsx, 0.05, 0.1, PhistVbs)
		sqrtPhistVbs  = sqrt(PhistVbs)
		Xdep          = T1DEP * sqrtPhistVbs
		Cdep          = epssi / Xdep
		cdsc          = CITEDGE_i + NFACTOREDGE_t + CDSCDEDGE_i * Vdsx - CDSCBEDGE_i * Vbsx
		T1            = 1.0 + cdsc/Cox
		Smooth(T1, 1.0, 0.05, n)
		nVt       = n * Vt
		inv_nVt   = 1.0 / nVt
		vg        = Vg * inv_nVt
		vs        = Vs * inv_nVt
		vfb       = VFB_i * inv_nVt
		dvth_dibl = -(ETA0EDGE_t + ETABEDGE_i * Vbsx) * Vdsx
		dvth_temp = (KT1EDGE_i + KT1LEDGE_i / Leff + KT2EDGE_i * Vbsx) * (pow(TRatio, KT1EXPEDGE_i) - 1.0)
		litl_edge = litl * (1.0 + DVT2EDGE * Vbsx)
		T0        = DVT1EDGE * Leff / litl_edge
		if (T0 < 40.0) 
			theta_sce_edge = 0.5 * DVT0EDGE / (cosh(T0) - 1.0)
		else 
			theta_sce_edge = DVT0EDGE * lexp(-T0)
		end
		dvth_sce  = theta_sce_edge * (Vbi_edge - Phist)
		Vth_shift = dvth_dibl - dvth_temp + dvth_sce + DVTEDGE + vth0_stress_EDGE - K2EDGE_i * Vbsx
		vgfb      = vg - vfb - Vth_shift * inv_nVt
	
		# Normalized body factor
		DGAMMAEDGE_i = DGAMMAEDGE * (1.0 + DGAMMAEDGEL * pow(Leff, -DGAMMAEDGELEXP))
		gam_edge          = sqrt(2.0 * q * epssi * NDEPEDGE_i * inv_nVt) / Cox
		gam_edge          = gam_edge * (1.0 + DGAMMAEDGE_i)
		inv_gam           = 1.0 / gam_edge
		
		# psip: pinch-off voltage
		phib_n_edge  = phib_edge / n
		PO_psip(vgfb, gam_edge, 0.0, phib_n_edge, psip)
		
		
		BSIM_q(psip, phib_n_edge, vs, gam_edge, qs_edge)
	
		# Approximate pinch-off voltage
		vdsatedge = 2.0 * nVt * qs_edge + 2.0 * nVt
		Vdsatedge = vdsatedge
		Vdsatedge = Vdsatedge + Vs
	
		# Vdssat clamped to avoid negative values during transient simulation
		Smooth(Vdsatedge - Vs, 0.0, 1.0e-3, Vdssate)
		T7     = pow(Vds / Vdssate, 1.0 / DELTA_t)
		T8     = pow(1.0 + T7, -DELTA_t)
		Vdseff = Vds * T8
		vdeff  = (Vdseff + Vs) * inv_nVt
		BSIM_q(psip, phib_n_edge, vdeff, gam_edge, qdeff_edge)
	
		# Nq calculation for Edge FET
		Smooth(psip, 1.0, 2.0, psipclamp)
		sqrtpsip = sqrt(psipclamp)
		psiavg   = psip - qs_edge - qdeff_edge -1.0
		Smooth(psiavg, 1.0, 2.0, T0)
		T2       = sqrt(T0)
		nq_edge  = 1.0 + gam_edge / (sqrtpsip + T2)
		ids_edge = 2.0 * NF * nq_edge * ueff * WEDGE / Leff * Cox * nVt * nVt * ((qs_edge - qdeff_edge) * (1.0 + qs_edge + qdeff_edge)) * Moc
		ids      = ids_edge + ids