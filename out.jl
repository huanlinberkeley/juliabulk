DELTA_1 = 0.02
Oneby3 = 0.33333333333333333
REFTEMP = 300.15       # 27 degrees C

# Clamped exponential function
function lexp(x)::Float64
	EXPL_THRESHOLD::Float64 = 80.0
	MAX_EXPL::Float64 = 5.540622384e34
	MIN_EXPL::Float64 = 1.804851387e-35
	if (x > EXPL_THRESHOLD)
		return MAX_EXPL * (1.0 + x - EXPL_THRESHOLD)
	elseif (x < -EXPL_THRESHOLD)
		return MIN_EXPL
	else
		return exp(x)
	end
end

# Clamped log function
function lln(x)::Float64
	N_MINLOG::Float64 = 1.0e-38
	return log(max(x, N_MINLOG))
end

# Hyperbolic smoothing function
function hypsmooth(x, c)::Float64
    return 0.5 * (x + sqrt(x * x + 4.0 * c * c))
end

function bsimbulk()
	ntype::Int8 = 1
	ptype::Int8 = -1
	q::Float64 = 1.60219e-19
	EPS0::Float64 = 8.85418e-12
	KboQ::Float64 = 8.617087e-5      # Joule/degree

	# Pure instance parameters
	L::Float64 = 1.0e-5
	W::Float64 = 1.0e-5
	NF::UInt16 = 1
	NRS::Float64 = 1.0
	NRD::Float64 = 1.0
	VFBSDOFF::Float64 = 0.0
	MINZ::UInt16 = 0
	RGATEMOD::UInt16 = 0
	RBODYMOD::UInt16 = 0
	GEOMOD::UInt16 = 0
	RGEOMOD::UInt16 = 0
	RBPB::Float64 = 50.0
	RBPD::Float64 = 50.0
	RBPS::Float64 = 50.0
	RBDB::Float64 = 50.0
	RBSB::Float64 = 50.0
	SA::Float64 = 0.0
	SB::Float64 = 0.0
	SD::Float64 = 0.0
	SCA::Float64 = 0.0
	SCB::Float64 = 0.0
	SCC::Float64 = 0.0
	SC::Float64 = 0.0
	AS::Float64 = 0.0
	AD::Float64 = 0.0
	PS::Float64 = 0.0
	PD::Float64 = 0.0

	# Both model and instance parameters
	XGW::Float64 = 0.0
	NGCON::UInt16 = 1
	DTEMP::Float64 = 0.0
	MULU0::Float64 = 1.0
	DELVTO::Float64 = 0.0
	IDS0MULT::Float64 = 1.0
	EDGEFET::UInt16 = 0
	SSLMOD::UInt16 = 0

	# Pure model parameters
	TYPE::UInt16 = 1
	CVMOD::UInt16 = 0
	COVMOD::UInt16 = 0
	RDSMOD::UInt16 = 0
	WPEMOD::UInt16 = 0
	ASYMMOD::UInt16 = 0
	GIDLMOD::UInt16 = 0
	IGCMOD::UInt16 = 0
	IGBMOD::UInt16 = 0
	TNOIMOD::UInt16 = 0
	SHMOD::UInt16 = 0
	MOBSCALE::UInt16 = 0

	# Device parameters
	LLONG::Float64 = 1.0e-5
	LMLT::Float64 = 1.0
	WMLT::Float64 = 1.0
	XL::Float64 = 0.0
	WWIDE::Float64 = 1.0e-5
	XW::Float64 = 0.0
	LINT::Float64 = 0.0
	LL::Float64 = 0.0
	LW::Float64 = 0.0
	LWL::Float64 = 0.0
	LLN::Float64 = 1.0
	LWN::Float64 = 1.0
	WINT::Float64 = 0.0
	WL::Float64 = 0.0
	WW::Float64 = 0.0
	WWL::Float64 = 0.0
	WLN::Float64 = 1.0
	WWN::Float64 = 1.0
	DLC::Float64 = 0.0
	LLC::Float64 = 0.0
	LWC::Float64 = 0.0
	LWLC::Float64 = 0.0
	DWC::Float64 = 0.0
	WLC::Float64 = 0.0
	WWC::Float64 = 0.0
	WWLC::Float64 = 0.0
	TOXE::Float64 = 3.0e-9
	TOXP::Float64 = TOXE
	DTOX::Float64 = 0.0
	NDEP::Float64 = 1e24
	NDEPL1::Float64 = 0.0
	NDEPLEXP1::Float64 = 1.0
	NDEPL2::Float64 = 0.0
	NDEPLEXP2::Float64 = 2.0
	NDEPW::Float64 = 0.0
	NDEPWEXP::Float64 = 1.0
	NDEPWL::Float64 = 0.0
	NDEPWLEXP::Float64 = 1.0
	LNDEP::Float64 = 0.0
	WNDEP::Float64 = 0.0
	PNDEP::Float64 = 0.0
	NDEPCV::Float64 = 1e24
	NDEPCVL1::Float64 = 0.0
	NDEPCVLEXP1::Float64 = 1.0
	NDEPCVL2::Float64 = 0.0
	NDEPCVLEXP2::Float64 = 2.0
	NDEPCVW::Float64 = 0.0
	NDEPCVWEXP::Float64 = 1.0
	NDEPCVWL::Float64 = 0.0
	NDEPCVWLEXP::Float64 = 1.0
	LNDEPCV::Float64 = 0.0
	WNDEPCV::Float64 = 0.0
	PNDEPCV::Float64 = 0.0
	NGATE::Float64 = 5e25
	LNGATE::Float64 = 0.0
	WNGATE::Float64 = 0.0
	PNGATE::Float64 = 0.0
	EASUB::Float64 = 4.05
	NI0SUB::Float64 = 1.1e16
	BG0SUB::Float64 = 1.17
	EPSRSUB::Float64 = 11.9
	EPSROX::Float64 = 3.9
	XJ::Float64 = 1.5e-7
	LXJ::Float64 = 0.0
	WXJ::Float64 = 0.0
	PXJ::Float64 = 0.0
	VFB::Float64 = -0.5
	LVFB::Float64 = 0.0
	WVFB::Float64 = 0.0
	PVFB::Float64 = 0.0
	VFBCV::Float64 = -0.5
	LVFBCV::Float64 = 0.0
	WVFBCV::Float64 = 0.0
	PVFBCV::Float64 = 0.0
	VFBCVL::Float64 = 0.0
	VFBCVLEXP::Float64 = 1.0
	VFBCVW::Float64 = 0.0
	VFBCVWEXP::Float64 = 1.0
	VFBCVWL::Float64 = 0.0
	VFBCVWLEXP::Float64 = 1.0
	DELVFBACC::Float64 = 0.0

	# Diode parameters
	PERMOD::UInt16 = 1
	DWJ::Float64 = DWC

	# Short channel effects
	NSD::Float64 = 1e26
	LNSD::Float64 = 0.0
	WNSD::Float64 = 0.0
	PNSD::Float64 = 0.0
	DVTP0::Float64 = 0.0
	LDVTP0::Float64 = 0
	WDVTP0::Float64 = 0
	PDVTP0::Float64 = 0
	DVTP1::Float64 = 0.0
	LDVTP1::Float64 = 0
	WDVTP1::Float64 = 0
	PDVTP1::Float64 = 0
	DVTP2::Float64 = 0.0
	LDVTP2::Float64 = 0
	WDVTP2::Float64 = 0
	PDVTP2::Float64 = 0
	DVTP3::Float64 = 0.0
	LDVTP3::Float64 = 0
	WDVTP3::Float64 = 0
	PDVTP3::Float64 = 0
	DVTP4::Float64 = 0.0
	LDVTP4::Float64 = 0
	WDVTP4::Float64 = 0
	PDVTP4::Float64 = 0
	DVTP5::Float64 = 0.0
	LDVTP5::Float64 = 0
	WDVTP5::Float64 = 0
	PDVTP5::Float64 = 0
	PHIN::Float64 = 0.045
	LPHIN::Float64 = 0.0
	WPHIN::Float64 = 0.0
	PPHIN::Float64 = 0.0
	ETA0::Float64 = 0.08
	LETA0::Float64 = 0.0
	WETA0::Float64 = 0.0
	PETA0::Float64 = 0.0
	ETA0R::Float64 = ETA0
	LETA0R::Float64 = LETA0
	WETA0R::Float64 = WETA0
	PETA0R::Float64 = PETA0
	DSUB::Float64 = 1.0
	ETAB::Float64 = -0.07
	ETABEXP::Float64 = 1.0
	LETAB::Float64 = 0.0
	WETAB::Float64 = 0.0
	PETAB::Float64 = 0.0
	K1::Float64 = 0.0
	K1L::Float64 = 0.0
	K1LEXP::Float64 = 1.0
	K1W::Float64 = 0.0
	K1WEXP::Float64 = 1.0
	K1WL::Float64 = 0.0
	K1WLEXP::Float64 = 1.0
	LK1::Float64 = 0.0
	WK1::Float64 = 0.0
	PK1::Float64 = 0.0
	K2::Float64 = 0.0
	K2L::Float64 = 0.0
	K2LEXP::Float64 = 1.0
	K2W::Float64 = 0.0
	K2WEXP::Float64 = 1.0
	K2WL::Float64 = 0.0
	K2WLEXP::Float64 = 1.0
	LK2::Float64 = 0.0
	WK2::Float64 = 0.0
	PK2::Float64 = 0.0

	# Quantum mechanical effects
	ADOS::Float64 = 0.0
	BDOS::Float64 = 1.0
	QM0::Float64 = 1.0e-3
	ETAQM::Float64 = 0.54

	# Sub-threshold swing factor
	CIT::Float64 = 0.0
	LCIT::Float64 = 0.0
	WCIT::Float64 = 0.0
	PCIT::Float64 = 0.0
	NFACTOR::Float64 = 0.0
	NFACTORL::Float64 = 0.0
	NFACTORLEXP::Float64 = 1.0
	NFACTORW::Float64 = 0.0
	NFACTORWEXP::Float64 = 1.0
	NFACTORWL::Float64 = 0.0
	NFACTORWLEXP::Float64 = 1.0
	LNFACTOR::Float64 = 0.0
	WNFACTOR::Float64 = 0.0
	PNFACTOR::Float64 = 0.0
	CDSCD::Float64 = 1e-9
	CDSCDL::Float64 = 0.0
	CDSCDLEXP::Float64 = 1.0
	LCDSCD::Float64 = 0.0
	WCDSCD::Float64 = 0.0
	PCDSCD::Float64 = 0.0
	CDSCDR::Float64 = CDSCD
	CDSCDLR::Float64 = CDSCDL
	LCDSCDR::Float64 = LCDSCD
	WCDSCDR::Float64 = WCDSCD
	PCDSCDR::Float64 = PCDSCD
	CDSCB::Float64 = 0.0
	CDSCBL::Float64 = 0.0
	CDSCBLEXP::Float64 = 1.0
	LCDSCB::Float64 = 0.0
	WCDSCB::Float64 = 0.0
	PCDSCB::Float64 = 0.0

	# Drain saturation voltage
	VSAT::Float64 = 1e5
	LVSAT::Float64 = 0.0
	WVSAT::Float64 = 0.0
	PVSAT::Float64 = 0.0
	VSATL::Float64 = 0.0
	VSATLEXP::Float64 = 1.0
	VSATW::Float64 = 0.0
	VSATWEXP::Float64 = 1.0
	VSATWL::Float64 = 0.0
	VSATWLEXP::Float64 = 1.0
	VSATR::Float64 = VSAT
	LVSATR::Float64 = LVSAT
	WVSATR::Float64 = WVSAT
	PVSATR::Float64 = PVSAT
	DELTA::Float64 = 0.125
	LDELTA::Float64 = 0.0
	WDELTA::Float64 = 0.0
	PDELTA::Float64 = 0.0
	DELTAL::Float64 = 0.0
	DELTALEXP::Float64 = 1.0
	VSATCV::Float64 = 1e5
	LVSATCV::Float64 = 0.0
	WVSATCV::Float64 = 0.0
	PVSATCV::Float64 = 0.0
	VSATCVL::Float64 = 0.0
	VSATCVLEXP::Float64 = 1.0
	VSATCVW::Float64 = 0.0
	VSATCVWEXP::Float64 = 1.0
	VSATCVWL::Float64 = 0.0
	VSATCVWLEXP::Float64 = 1.0

	# Mobility degradation
	UP1::Float64 = 0.0
	LP1::Float64 = 1.0e-8
	UP2::Float64 = 0.0
	LP2::Float64 = 1.0e-8
	U0::Float64 = 67.0e-3
	U0L::Float64 = 0.0
	U0LEXP::Float64 = 1.0
	LU0::Float64 = 0.0
	WU0::Float64 = 0.0
	PU0::Float64 = 0.0
	U0R::Float64 = U0
	LU0R::Float64 = LU0
	WU0R::Float64 = WU0
	PU0R::Float64 = PU0
	ETAMOB::Float64 = 1.0
	UA::Float64 = 0.001
	UAL::Float64 = 0.0
	UALEXP::Float64 = 1.0
	UAW::Float64 = 0.0
	UAWEXP::Float64 = 1.0
	UAWL::Float64 = 0.0
	UAWLEXP::Float64 = 1.0
	LUA::Float64 = 0.0
	WUA::Float64 = 0.0
	PUA::Float64 = 0.0
	UAR::Float64 = UA
	LUAR::Float64 = LUA
	WUAR::Float64 = WUA
	PUAR::Float64 = PUA
	EU::Float64 = 1.5
	LEU::Float64 = 0.0
	WEU::Float64 = 0.0
	PEU::Float64 = 0.0
	EUL::Float64 = 0.0
	EULEXP::Float64 = 1.0
	EUW::Float64 = 0.0
	EUWEXP::Float64 = 1.0
	EUWL::Float64 = 0.0
	EUWLEXP::Float64 = 1.0
	UD::Float64 = 0.001
	UDL::Float64 = 0.0
	UDLEXP::Float64 = 1.0
	LUD::Float64 = 0.0
	WUD::Float64 = 0.0
	PUD::Float64 = 0.0
	UDR::Float64 = UD
	LUDR::Float64 = LUD
	WUDR::Float64 = WUD
	PUDR::Float64 = PUD
	UCS::Float64 = 2.0
	LUCS::Float64 = 0.0
	WUCS::Float64 = 0.0
	PUCS::Float64 = 0.0
	UCSR::Float64 = UCS
	LUCSR::Float64 = LUCS
	WUCSR::Float64 = WUCS
	PUCSR::Float64 = PUCS
	UC::Float64 = 0.0
	UCL::Float64 = 0.0
	UCLEXP::Float64 = 1.0
	UCW::Float64 = 0.0
	UCWEXP::Float64 = 1.0
	UCWL::Float64 = 0.0
	UCWLEXP::Float64 = 1.0
	LUC::Float64 = 0.0
	WUC::Float64 = 0.0
	PUC::Float64 = 0.0
	UCR::Float64 = UC
	LUCR::Float64 = LUC
	WUCR::Float64 = WUC
	PUCR::Float64 = PUC

	# Channel length modulation
	PCLM::Float64 = 0.0
	PCLML::Float64 = 0.0
	PCLMLEXP::Float64 = 1.0
	LPCLM::Float64 = 0.0
	WPCLM::Float64 = 0.0
	PPCLM::Float64 = 0.0
	PCLMR::Float64 = PCLM
	LPCLMR::Float64 = LPCLM
	WPCLMR::Float64 = WPCLM
	PPCLMR::Float64 = PPCLM
	PCLMG::Float64 = 0.0
	PCLMCV::Float64 = PCLM
	PCLMCVL::Float64 = PCLML
	PCLMCVLEXP::Float64 = PCLMLEXP
	LPCLMCV::Float64 = LPCLM
	WPCLMCV::Float64 = WPCLM
	PPCLMCV::Float64 = PPCLM
	PSCBE1::Float64 = 4.24e8
	LPSCBE1::Float64 = 0.0
	WPSCBE1::Float64 = 0.0
	PPSCBE1::Float64 = 0.0
	PSCBE2::Float64 = 1.0e-8
	LPSCBE2::Float64 = 0.0
	WPSCBE2::Float64 = 0.0
	PPSCBE2::Float64 = 0.0
	PDITS::Float64 = 0.0
	LPDITS::Float64 = 0.0
	WPDITS::Float64 = 0.0
	PPDITS::Float64 = 0.0
	PDITSL::Float64 = 0.0
	PDITSD::Float64 = 0.0
	LPDITSD::Float64 = 0.0
	WPDITSD::Float64 = 0.0
	PPDITSD::Float64 = 0.0

	# S/D series resistances
	RSH::Float64 = 0.0
	PRWG::Float64 = 1.0
	LPRWG::Float64 = 0.0
	WPRWG::Float64 = 0.0
	PPRWG::Float64 = 0.0
	PRWB::Float64 = 0.0
	LPRWB::Float64 = 0.0
	WPRWB::Float64 = 0.0
	PPRWB::Float64 = 0.0
	PRWBL::Float64 = 0.0
	PRWBLEXP::Float64 = 1.0
	WR::Float64 = 1.0
	LWR::Float64 = 0.0
	WWR::Float64 = 0.0
	PWR::Float64 = 0.0
	RSWMIN::Float64 = 0.0
	LRSWMIN::Float64 = 0.0
	WRSWMIN::Float64 = 0.0
	PRSWMIN::Float64 = 0.0
	RSW::Float64 = 10.0
	LRSW::Float64 = 0.0
	WRSW::Float64 = 0.0
	PRSW::Float64 = 0.0
	RSWL::Float64 = 0.0
	RSWLEXP::Float64 = 1.0
	RDWMIN::Float64 = RSWMIN
	LRDWMIN::Float64 = LRSWMIN
	WRDWMIN::Float64 = WRSWMIN
	PRDWMIN::Float64 = PRSWMIN
	RDW::Float64 = RSW
	LRDW::Float64 = LRSW
	WRDW::Float64 = WRSW
	PRDW::Float64 = PRSW
	RDWL::Float64 = RSWL
	RDWLEXP::Float64 = RSWLEXP
	RDSWMIN::Float64 = 0.0
	LRDSWMIN::Float64 = 0.0
	WRDSWMIN::Float64 = 0.0
	PRDSWMIN::Float64 = 0.0
	RDSW::Float64 = 20.0
	RDSWL::Float64 = 0.0
	RDSWLEXP::Float64 = 1.0
	LRDSW::Float64 = 0.0
	WRDSW::Float64 = 0.0
	PRDSW::Float64 = 0.0

	# Velocity saturation
	PSAT::Float64 = 1.0
	LPSAT::Float64 = 0.0
	WPSAT::Float64 = 0.0
	PPSAT::Float64 = 0.0
	PSATL::Float64 = 0.0
	PSATLEXP::Float64 = 1.0
	PSATB::Float64 = 0.0
	PSATR::Float64 = PSAT
	LPSATR::Float64 = LPSAT
	WPSATR::Float64 = WPSAT
	PPSATR::Float64 = PPSAT
	LPSATB::Float64 = 0.0
	WPSATB::Float64 = 0.0
	PPSATB::Float64 = 0.0
	PSATX::Float64 = 1.0
	PTWG::Float64 = 0.0
	LPTWG::Float64 = 0.0
	WPTWG::Float64 = 0.0
	PPTWG::Float64 = 0.0
	PTWGL::Float64 = 0.0
	PTWGLEXP::Float64 = 1.0
	PTWGR::Float64 = PTWG
	LPTWGR::Float64 = LPTWG
	WPTWGR::Float64 = WPTWG
	PPTWGR::Float64 = PPTWG
	PTWGLR::Float64 = PTWGL
	PTWGLEXPR::Float64 = PTWGLEXP

	# Velocity non-saturation effect
	A1::Float64 = 0.0
	LA1::Float64 = 0.0
	WA1::Float64 = 0.0
	PA1::Float64 = 0.0
	A11::Float64 = 0.0
	LA11::Float64 = 0.0
	WA11::Float64 = 0.0
	PA11::Float64 = 0.0
	A2::Float64 = 0.0
	LA2::Float64 = 0.0
	WA2::Float64 = 0.0
	PA2::Float64 = 0.0
	A21::Float64 = 0.0
	LA21::Float64 = 0.0
	WA21::Float64 = 0.0
	PA21::Float64 = 0.0

	# Output conductance
	PDIBLC::Float64 = 0.0
	PDIBLCL::Float64 = 0.0
	PDIBLCLEXP::Float64 = 1.0
	LPDIBLC::Float64 = 0.0
	WPDIBLC::Float64 = 0.0
	PPDIBLC::Float64 = 0.0
	PDIBLCR::Float64 = PDIBLC
	PDIBLCLR::Float64 = PDIBLCL
	PDIBLCLEXPR::Float64 = PDIBLCLEXP
	LPDIBLCR::Float64 = LPDIBLC
	WPDIBLCR::Float64 = WPDIBLC
	PPDIBLCR::Float64 = PPDIBLC
	PDIBLCB::Float64 = 0.0
	LPDIBLCB::Float64 = 0.0
	WPDIBLCB::Float64 = 0.0
	PPDIBLCB::Float64 = 0.0
	PVAG::Float64 = 1.0
	LPVAG::Float64 = 0.0
	WPVAG::Float64 = 0.0
	PPVAG::Float64 = 0.0
	FPROUT::Float64 = 0.0
	FPROUTL::Float64 = 0.0
	FPROUTLEXP::Float64 = 1.0
	LFPROUT::Float64 = 0.0
	WFPROUT::Float64 = 0.0
	PFPROUT::Float64 = 0.0

	# Impact ionization current
	ALPHA0::Float64 = 0.0
	ALPHA0L::Float64 = 0.0
	ALPHA0LEXP::Float64 = 1.0
	LALPHA0::Float64 = 0.0
	WALPHA0::Float64 = 0.0
	PALPHA0::Float64 = 0.0
	BETA0::Float64 = 0.0
	LBETA0::Float64 = 0.0
	WBETA0::Float64 = 0.0
	PBETA0::Float64 = 0.0

	# Gate dielectric tunnelling current model parameters
	AIGBACC::Float64 = 1.36e-2
	BIGBACC::Float64 = 1.71e-3
	CIGBACC::Float64 = 0.075
	NIGBACC::Float64 = 1.0
	AIGBINV::Float64 = 1.11e-2
	BIGBINV::Float64 = 9.49e-4
	CIGBINV::Float64 = 0.006
	EIGBINV::Float64 = 1.1
	NIGBINV::Float64 = 3.0
	AIGC::Float64 = ((TYPE == 1) ? 1.36e-2 : 9.8e-3)
	BIGC::Float64 = ((TYPE == 1) ? 1.71e-3 : 7.59e-4)
	CIGC::Float64 = ((TYPE == 1) ? 0.075 : 0.03)
	AIGS::Float64 = ((TYPE == 1) ? 1.36e-2 : 9.8e-3)
	BIGS::Float64 = ((TYPE == 1) ? 1.71e-3 : 7.59e-4)
	CIGS::Float64 = ((TYPE == 1) ? 0.075 : 0.03)
	AIGD::Float64 = ((TYPE == 1) ? 1.36e-2 : 9.8e-3)
	BIGD::Float64 = ((TYPE == 1) ? 1.71e-3 : 7.59e-4)
	CIGD::Float64 = ((TYPE == 1) ? 0.075 : 0.03)
	DLCIG::Float64 = LINT
	DLCIGD::Float64 = DLCIG
	POXEDGE::Float64 = 1.0
	NTOX::Float64 = 1.0
	TOXREF::Float64 = 3.0e-9
	PIGCD::Float64 = 1.0
	AIGCL::Float64 = 0.0
	AIGCW::Float64 = 0.0
	AIGSL::Float64 = 0.0
	AIGSW::Float64 = 0.0
	AIGDL::Float64 = 0.0
	AIGDW::Float64 = 0.0
	PIGCDL::Float64 = 0.0
	LAIGBINV::Float64 = 0.0
	WAIGBINV::Float64 = 0.0
	PAIGBINV::Float64 = 0.0
	LBIGBINV::Float64 = 0.0
	WBIGBINV::Float64 = 0.0
	PBIGBINV::Float64 = 0.0
	LCIGBINV::Float64 = 0.0
	WCIGBINV::Float64 = 0.0
	PCIGBINV::Float64 = 0.0
	LEIGBINV::Float64 = 0.0
	WEIGBINV::Float64 = 0.0
	PEIGBINV::Float64 = 0.0
	LNIGBINV::Float64 = 0.0
	WNIGBINV::Float64 = 0.0
	PNIGBINV::Float64 = 0.0
	LAIGBACC::Float64 = 0.0
	WAIGBACC::Float64 = 0.0
	PAIGBACC::Float64 = 0.0
	LBIGBACC::Float64 = 0.0
	WBIGBACC::Float64 = 0.0
	PBIGBACC::Float64 = 0.0
	LCIGBACC::Float64 = 0.0
	WCIGBACC::Float64 = 0.0
	PCIGBACC::Float64 = 0.0
	LNIGBACC::Float64 = 0.0
	WNIGBACC::Float64 = 0.0
	PNIGBACC::Float64 = 0.0
	LAIGC::Float64 = 0.0
	WAIGC::Float64 = 0.0
	PAIGC::Float64 = 0.0
	LBIGC::Float64 = 0.0
	WBIGC::Float64 = 0.0
	PBIGC::Float64 = 0.0
	LCIGC::Float64 = 0.0
	WCIGC::Float64 = 0.0
	PCIGC::Float64 = 0.0
	LAIGS::Float64 = 0.0
	WAIGS::Float64 = 0.0
	PAIGS::Float64 = 0.0
	LBIGS::Float64 = 0.0
	WBIGS::Float64 = 0.0
	PBIGS::Float64 = 0.0
	LCIGS::Float64 = 0.0
	WCIGS::Float64 = 0.0
	PCIGS::Float64 = 0.0
	LAIGD::Float64 = 0.0
	WAIGD::Float64 = 0.0
	PAIGD::Float64 = 0.0
	LBIGD::Float64 = 0.0
	WBIGD::Float64 = 0.0
	PBIGD::Float64 = 0.0
	LCIGD::Float64 = 0.0
	WCIGD::Float64 = 0.0
	PCIGD::Float64 = 0.0
	LPOXEDGE::Float64 = 0.0
	WPOXEDGE::Float64 = 0.0
	PPOXEDGE::Float64 = 0.0
	LDLCIG::Float64 = 0.0
	WDLCIG::Float64 = 0.0
	PDLCIG::Float64 = 0.0
	LDLCIGD::Float64 = 0.0
	WDLCIGD::Float64 = 0.0
	PDLCIGD::Float64 = 0.0
	LNTOX::Float64 = 0.0
	WNTOX::Float64 = 0.0
	PNTOX::Float64 = 0.0

	# GIDL and GISL currents
	AGIDL::Float64 = 0.0
	AGIDLL::Float64 = 0.0
	AGIDLW::Float64 = 0.0
	LAGIDL::Float64 = 0.0
	WAGIDL::Float64 = 0.0
	PAGIDL::Float64 = 0.0
	BGIDL::Float64 = 2.3e9
	LBGIDL::Float64 = 0.0
	WBGIDL::Float64 = 0.0
	PBGIDL::Float64 = 0.0
	CGIDL::Float64 = 0.5
	LCGIDL::Float64 = 0.0
	WCGIDL::Float64 = 0.0
	PCGIDL::Float64 = 0.0
	EGIDL::Float64 = 0.8
	LEGIDL::Float64 = 0.0
	WEGIDL::Float64 = 0.0
	PEGIDL::Float64 = 0.0
	AGISL::Float64 = AGIDL
	AGISLL::Float64 = AGIDLL
	AGISLW::Float64 = AGIDLW
	LAGISL::Float64 = LAGIDL
	WAGISL::Float64 = WAGIDL
	PAGISL::Float64 = PAGIDL
	BGISL::Float64 = BGIDL
	LBGISL::Float64 = LBGIDL
	WBGISL::Float64 = WBGIDL
	PBGISL::Float64 = PBGIDL
	CGISL::Float64 = CGIDL
	LCGISL::Float64 = LCGIDL
	WCGISL::Float64 = WCGIDL
	PCGISL::Float64 = PCGIDL
	EGISL::Float64 = EGIDL
	LEGISL::Float64 = LEGIDL
	WEGISL::Float64 = WEGIDL
	PEGISL::Float64 = PEGIDL

	# Overlap capacitance and fringing capacitance
	CF::Float64 = 0.0
	LCF::Float64 = 0.0
	WCF::Float64 = 0.0
	PCF::Float64 = 0.0
	CFRCOEFF::Float64 = 1.0
	CGSO::Float64 = 0.0
	CGDO::Float64 = 0.0
	CGBO::Float64 = 0.0
	CGSL::Float64 = 0.0
	LCGSL::Float64 = 0.0
	WCGSL::Float64 = 0.0
	PCGSL::Float64 = 0.0
	CGDL::Float64 = 0.0
	LCGDL::Float64 = 0.0
	WCGDL::Float64 = 0.0
	PCGDL::Float64 = 0.0
	CKAPPAS::Float64 = 0.6
	LCKAPPAS::Float64 = 0.0
	WCKAPPAS::Float64 = 0.0
	PCKAPPAS::Float64 = 0.0
	CKAPPAD::Float64 = 0.6
	LCKAPPAD::Float64 = 0.0
	WCKAPPAD::Float64 = 0.0
	PCKAPPAD::Float64 = 0.0

	# Layout-dependent parasitics model parameters (resistance only)
	DMCG::Float64 = 0.0
	DMCI::Float64 = DMCG
	DMDG::Float64 = 0.0
	DMCGT::Float64 = 0.0
	XGL::Float64 = 0.0
	RSHG::Float64 = 0.1

	# Junction capacitance
	CJS::Float64 = 5.0e-4
	CJD::Float64 = CJS
	CJSWS::Float64 = 5.0e-10
	CJSWD::Float64 = CJSWS
	CJSWGS::Float64 = 0.0
	CJSWGD::Float64 = CJSWGS
	PBS::Float64 = 1.0
	PBD::Float64 = PBS
	PBSWS::Float64 = 1.0
	PBSWD::Float64 = PBSWS
	PBSWGS::Float64 = PBSWS
	PBSWGD::Float64 = PBSWGS
	MJS::Float64 = 0.5
	MJD::Float64 = MJS
	MJSWS::Float64 = 0.33
	MJSWD::Float64 = MJSWS
	MJSWGS::Float64 = MJSWS
	MJSWGD::Float64 = MJSWGS

	# Junction current
	JSS::Float64 = 1.0e-4
	JSD::Float64 = JSS
	JSWS::Float64 = 0.0
	JSWD::Float64 = JSWS
	JSWGS::Float64 = 0.0
	JSWGD::Float64 = JSWGS
	NJS::Float64 = 1.0
	NJD::Float64 = NJS
	IJTHSFWD::Float64 = 0.1
	IJTHDFWD::Float64 = IJTHSFWD
	IJTHSREV::Float64 = 0.1
	IJTHDREV::Float64 = IJTHSREV
	BVS::Float64 = 10.0
	BVD::Float64 = BVS
	XJBVS::Float64 = 1.0
	XJBVD::Float64 = XJBVS

	# Tunneling component of junction current
	JTSS::Float64 = 0.0
	JTSD::Float64 = JTSS
	JTSSWS::Float64 = 0.0
	JTSSWD::Float64 = JTSSWS
	JTSSWGS::Float64 = 0.0
	JTSSWGD::Float64 = JTSSWGS
	JTWEFF::Float64 = 0.0
	NJTS::Float64 = 20.0
	NJTSD::Float64 = NJTS
	NJTSSW::Float64 = 20.0
	NJTSSWD::Float64 = NJTSSW
	NJTSSWG::Float64 = 20.0
	NJTSSWGD::Float64 = NJTSSWG
	VTSS::Float64 = 10.0
	VTSD::Float64 = VTSS
	VTSSWS::Float64 = 10.0
	VTSSWD::Float64 = VTSSWS
	VTSSWGS::Float64 = 10.0
	VTSSWGD::Float64 = VTSSWGS

	# High-speed/RF model parameters
	XRCRG1::Float64 = 12.0
	XRCRG2::Float64 = 1.0
	GBMIN::Float64 = 1.0e-12
	RBPS0::Float64 = 50.0
	RBPSL::Float64 = 0.0
	RBPSW::Float64 = 0.0
	RBPSNF::Float64 = 0.0
	RBPD0::Float64 = 50.0
	RBPDL::Float64 = 0.0
	RBPDW::Float64 = 0.0
	RBPDNF::Float64 = 0.0
	RBPBX0::Float64 = 100.0
	RBPBXL::Float64 = 0.0
	RBPBXW::Float64 = 0.0
	RBPBXNF::Float64 = 0.0
	RBPBY0::Float64 = 100.0
	RBPBYL::Float64 = 0.0
	RBPBYW::Float64 = 0.0
	RBPBYNF::Float64 = 0.0
	RBSBX0::Float64 = 100.0
	RBSBY0::Float64 = 100.0
	RBDBX0::Float64 = 100.0
	RBDBY0::Float64 = 100.0
	RBSDBXL::Float64 = 0.0
	RBSDBXW::Float64 = 0.0
	RBSDBXNF::Float64 = 0.0
	RBSDBYL::Float64 = 0.0
	RBSDBYW::Float64 = 0.0
	RBSDBYNF::Float64 = 0.0

	# Flicker noise
	EF::Float64 = 1.0
	EM::Float64 = 4.1e7
	NOIA::Float64 = 6.250e40
	NOIB::Float64 = 3.125e25
	NOIC::Float64 = 8.750e8
	LINTNOI::Float64 = 0.0
	NOIA1::Float64 = 0.0
	NOIAX::Float64 = 1.0

	# Thermal noise
	NTNOI::Float64 = 1.0
	RNOIA::Float64 = 0.577
	RNOIB::Float64 = 0.5164
	RNOIC::Float64 = 0.395
	TNOIA::Float64 = 1.5
	TNOIB::Float64 = 3.5
	TNOIC::Float64 = 0.0

	# Binning parameters
	BINUNIT::UInt16 = 1
	DLBIN::Float64 = 0.0
	DWBIN::Float64 = 0.0

	# Temperature dependence parameters
	TNOM::Float64 = 27.0
	TBGASUB::Float64 = 4.73e-4
	TBGBSUB::Float64 = 636.0
	TNFACTOR::Float64 = 0.0
	UTE::Float64 = -1.5
	LUTE::Float64 = 0.0
	WUTE::Float64 = 0.0
	PUTE::Float64 = 0.0
	UTEL::Float64 = 0.0
	UA1::Float64 = 1.0e-3
	LUA1::Float64 = 0.0
	WUA1::Float64 = 0.0
	PUA1::Float64 = 0.0
	UA1L::Float64 = 0.0
	UC1::Float64 = 0.056e-9
	LUC1::Float64 = 0.0
	WUC1::Float64 = 0.0
	PUC1::Float64 = 0.0
	UD1::Float64 = 0.0
	LUD1::Float64 = 0.0
	WUD1::Float64 = 0.0
	PUD1::Float64 = 0.0
	UD1L::Float64 = 0.0
	EU1::Float64 = 0.0
	LEU1::Float64 = 0.0
	WEU1::Float64 = 0.0
	PEU1::Float64 = 0.0
	UCSTE::Float64 = -4.775e-3
	LUCSTE::Float64 = 0.0
	WUCSTE::Float64 = 0.0
	PUCSTE::Float64 = 0.0
	TETA0::Float64 = 0.0
	PRT::Float64 = 0.0
	LPRT::Float64 = 0.0
	WPRT::Float64 = 0.0
	PPRT::Float64 = 0.0
	AT::Float64 = -1.56e-3
	LAT::Float64 = 0.0
	WAT::Float64 = 0.0
	PAT::Float64 = 0.0
	ATL::Float64 = 0.0
	TDELTA::Float64 = 0.0
	PTWGT::Float64 = 0.0
	LPTWGT::Float64 = 0.0
	WPTWGT::Float64 = 0.0
	PPTWGT::Float64 = 0.0
	PTWGTL::Float64 = 0.0
	KT1::Float64 = -0.11
	KT1EXP::Float64 = 1.0
	KT1L::Float64 = 0.0
	LKT1::Float64 = 0.0
	WKT1::Float64 = 0.0
	PKT1::Float64 = 0.0
	KT2::Float64 = 0.022
	LKT2::Float64 = 0.0
	WKT2::Float64 = 0.0
	PKT2::Float64 = 0.0
	IIT::Float64 = 0.0
	LIIT::Float64 = 0.0
	WIIT::Float64 = 0.0
	PIIT::Float64 = 0.0
	IGT::Float64 = 2.5
	LIGT::Float64 = 0.0
	WIGT::Float64 = 0.0
	PIGT::Float64 = 0.0
	TGIDL::Float64 = 0.0
	LTGIDL::Float64 = 0.0
	WTGIDL::Float64 = 0.0
	PTGIDL::Float64 = 0.0
	TCJ::Float64 = 0.0
	TCJSW::Float64 = 0.0
	TCJSWG::Float64 = 0.0
	TPB::Float64 = 0.0
	TPBSW::Float64 = 0.0
	TPBSWG::Float64 = 0.0
	XTIS::Float64 = 3.0
	XTID::Float64 = XTIS
	XTSS::Float64 = 0.02
	XTSD::Float64 = XTSS
	XTSSWS::Float64 = 0.02
	XTSSWD::Float64 = XTSSWS
	XTSSWGS::Float64 = 0.02
	XTSSWGD::Float64 = XTSSWGS
	TNJTS::Float64 = 0.0
	TNJTSD::Float64 = TNJTS
	TNJTSSW::Float64 = 0.0
	TNJTSSWD::Float64 = TNJTSSW
	TNJTSSWG::Float64 = 0.0
	TNJTSSWGD::Float64 = TNJTSSWG

	# Self heating parameters
	RTH0::Float64 = 0.0
	CTH0::Float64 = 1.0e-5
	WTH0::Float64 = 0.0

	# Stress related parameters
	SAREF::Float64 = 1.0e-6
	SBREF::Float64 = 1.0e-6
	WLOD::Float64 = 0.0
	KU0::Float64 = 0.0
	KVSAT::Float64 = 0.0
	TKU0::Float64 = 0.0
	LKU0::Float64 = 0.0
	WKU0::Float64 = 0.0
	PKU0::Float64 = 0.0
	LLODKU0::Float64 = 0.0
	WLODKU0::Float64 = 0.0
	KVTH0::Float64 = 0.0
	LKVTH0::Float64 = 0.0
	WKVTH0::Float64 = 0.0
	PKVTH0::Float64 = 0.0
	LLODVTH::Float64 = 0.0
	WLODVTH::Float64 = 0.0
	STK2::Float64 = 0.0
	LODK2::Float64 = 0.0
	STETA0::Float64 = 0.0
	LODETA0::Float64 = 0.0

	# Well proximity parameters
	WEB::Float64 = 0.0
	WEC::Float64 = 0.0
	KVTH0WE::Float64 = 0.0
	LKVTH0WE::Float64 = 0.0
	WKVTH0WE::Float64 = 0.0
	PKVTH0WE::Float64 = 0.0
	K2WE::Float64 = 0.0
	LK2WE::Float64 = 0.0
	WK2WE::Float64 = 0.0
	PK2WE::Float64 = 0.0
	KU0WE::Float64 = 0.0
	LKU0WE::Float64 = 0.0
	WKU0WE::Float64 = 0.0
	PKU0WE::Float64 = 0.0
	SCREF::Float64 = 1.0e-6

	# Sub-surface leakage drain current
	SSL0::Float64 = 4.0e2
	SSL1::Float64 = 3.36e8
	SSL2::Float64 = 0.185
	SSL3::Float64 = 0.3
	SSL4::Float64 = 1.4
	SSL5::Float64 = 0
	SSLEXP1::Float64 = 0.490
	SSLEXP2::Float64 = 1.42

	# Vdsx smoothing
	AVDSX::Float64 = 20.0

	# STI edge FET device parameters
	WEDGE::Float64 = 10.0e-9
	DGAMMAEDGE::Float64 = 0.0
	DGAMMAEDGEL::Float64 = 0.0
	DGAMMAEDGELEXP::Float64 = 1.0
	DVTEDGE::Float64 = 0.0
	NDEPEDGE::Float64 = 1e24
	LNDEPEDGE::Float64 = 0.0
	WNDEPEDGE::Float64 = 0.0
	PNDEPEDGE::Float64 = 0.0
	NFACTOREDGE::Float64 = 0.0
	LNFACTOREDGE::Float64 = 0.0
	WNFACTOREDGE::Float64 = 0.0
	PNFACTOREDGE::Float64 = 0.0
	CITEDGE::Float64 = 0.0
	LCITEDGE::Float64 = 0.0
	WCITEDGE::Float64 = 0.0
	PCITEDGE::Float64 = 0.0
	CDSCDEDGE::Float64 = 1e-9
	LCDSCDEDGE::Float64 = 0.0
	WCDSCDEDGE::Float64 = 0.0
	PCDSCDEDGE::Float64 = 0.0
	CDSCBEDGE::Float64 = 0.0
	LCDSCBEDGE::Float64 = 0.0
	WCDSCBEDGE::Float64 = 0.0
	PCDSCBEDGE::Float64 = 0.0
	ETA0EDGE::Float64 = 0.08
	LETA0EDGE::Float64 = 0.0
	WETA0EDGE::Float64 = 0.0
	PETA0EDGE::Float64 = 0.0
	ETABEDGE::Float64 = -0.07
	LETABEDGE::Float64 = 0.0
	WETABEDGE::Float64 = 0.0
	PETABEDGE::Float64 = 0.0
	KT1EDGE::Float64 = -0.11
	LKT1EDGE::Float64 = 0.0
	WKT1EDGE::Float64 = 0.0
	PKT1EDGE::Float64 = 0.0
	KT1LEDGE::Float64 = 0.0
	LKT1LEDGE::Float64 = 0.0
	WKT1LEDGE::Float64 = 0.0
	PKT1LEDGE::Float64 = 0.0
	KT2EDGE::Float64 = 0.022
	LKT2EDGE::Float64 = 0.0
	WKT2EDGE::Float64 = 0.0
	PKT2EDGE::Float64 = 0.0
	KT1EXPEDGE::Float64 = 1.0
	LKT1EXPEDGE::Float64 = 0.0
	WKT1EXPEDGE::Float64 = 0.0
	PKT1EXPEDGE::Float64 = 0.0
	TNFACTOREDGE::Float64 = 0.0
	LTNFACTOREDGE::Float64 = 0.0
	WTNFACTOREDGE::Float64 = 0.0
	PTNFACTOREDGE::Float64 = 0.0
	TETA0EDGE::Float64 = 0.0
	LTETA0EDGE::Float64 = 0.0
	WTETA0EDGE::Float64 = 0.0
	PTETA0EDGE::Float64 = 0.0
	DVT0EDGE::Float64 = 2.2
	DVT1EDGE::Float64 = 0.53
	DVT2EDGE::Float64 = 0.0
	K2EDGE::Float64 = 0.0
	LK2EDGE::Float64 = 0.0
	WK2EDGE::Float64 = 0.0
	PK2EDGE::Float64 = 0.0
	KVTH0EDGE::Float64 = 0.0
	LKVTH0EDGE::Float64 = 0.0
	WKVTH0EDGE::Float64 = 0.0
	PKVTH0EDGE::Float64 = 0.0
	STK2EDGE::Float64 = 0.0
	LSTK2EDGE::Float64 = 0.0
	WSTK2EDGE::Float64 = 0.0
	PSTK2EDGE::Float64 = 0.0
	STETA0EDGE::Float64 = 0.0
	LSTETA0EDGE::Float64 = 0.0
	WSTETA0EDGE::Float64 = 0.0
	PSTETA0EDGE::Float64 = 0.0
	IGCLAMP::UInt16 = 1
	LP::Float64 = 1.0e-5
	RNOIK::Float64 = 0.0
	TNOIK::Float64 = 0.0
	TNOIK2::Float64 = 0.1
	K0::Float64 = 0.0
	LK0::Float64 = 0.0
	WK0::Float64 = 0.0
	PK0::Float64 = 0.0
	K01::Float64 = 0.0
	LK01::Float64 = 0.0
	WK01::Float64 = 0.0
	PK01::Float64 = 0.0
	M0::Float64 = 1.0
	LM0::Float64 = 0.0
	WM0::Float64 = 0.0
	PM0::Float64 = 0.0
	M01::Float64 = 0.0
	LM01::Float64 = 0.0
	WM01::Float64 = 0.0
	PM01::Float64 = 0.0

	# Flicker noise model parameter for EDGE FET transistor
	NEDGE::Float64 = 1
	NOIA1_EDGE::Float64 = 0.0
	NOIAX_EDGE::Float64 = 1.0


	# Flicker noise model parameter for Halo transistor
	FNOIMOD::UInt16 = 0
	LH::Float64 = 1.0e-8
	NOIA2::Float64 = NOIA
	HNDEP::Float64 = NDEP

	# Flexibility of tuning Cgg in strong inversion
	ABULK::Float64 = 1.0

	# To enhance the fitting flexiblity for the gm/Id
	C0::Float64 = 0.0
	LC0::Float64 = 0.0
	WC0::Float64 = 0.0
	PC0::Float64 = 0.0
	C01::Float64 = 0.0
	LC01::Float64 = 0.0
	WC01::Float64 = 0.0
	PC01::Float64 = 0.0
	C0SI::Float64 = 1.0
	LC0SI::Float64 = 0.0
	WC0SI::Float64 = 0.0
	PC0SI::Float64 = 0.0
	C0SI1::Float64 = 0.0
	LC0SI1::Float64 = 0.0
	WC0SI1::Float64 = 0.0
	PC0SI1::Float64 = 0.0
	C0SISAT::Float64 = 0.0
	LC0SISAT::Float64 = 0.0
	WC0SISAT::Float64 = 0.0
	PC0SISAT::Float64 = 0.0
	C0SISAT1::Float64 = 0.0
	LC0SISAT1::Float64 = 0.0
	WC0SISAT1::Float64 = 0.0
	PC0SISAT1::Float64 = 0.0

	# Minimum resistance value
	minr::Float64 = 1.0e-3
	# High Voltage Model Parameters

	# --- Mod selectors -----
	HVMOD::UInt16 = 0
	HVCAP::UInt16 = 0
	HVCAPS::UInt16 = 0
	IIMOD::UInt16 = 0

	# --- Other parameters -----
	NDRIFTD::Float64 = 5.0e16
	VDRIFT::Float64 = 1.0e5
	MDRIFT::Float64 = 1.0
	NDRIFTS::Float64 = NDRIFTD
	RDLCW::Float64 = 100.0
	RSLCW::Float64 = 0
	PDRWB::Float64 = 0
	VFBDRIFT::Float64 = -1
	VFBOV::Float64 = -1
	LOVER::Float64 = 500e-9
	LOVERACC::Float64 = LOVER
	NDR::Float64 = NDEP
	SLHV::Float64 = 0
	SLHV1::Float64 = 1.0
	ALPHADR::Float64 = ALPHA0
	BETADR::Float64 = BETA0
	PRTHV::Float64 = 0.0
	ATHV::Float64 = 0
	HVFACTOR::Float64 = 1e-3
	DRII1::Float64 = 1.0
	DRII2::Float64 = 5
	DELTAII::Float64 = 0.5


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
    #if not_given(TOXP) 
    #    BSIMBULKTOXP = (TOXE * EPSROX / 3.9) - DTOX
    #else 
        BSIMBULKTOXP = TOXP
    #end

    L_mult = L * LMLT
    W_mult = W * WMLT
    Lnew = L_mult + XL
    W_by_NF = W_mult / NF
    Wnew    = W_by_NF + XW



	AIGC
end

yy=bsimbulk()
print(yy)