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

# Power function
function pow(x, y)::Float64
	return x^y
end

# Smoothing function for (max of x, x0 with deltax)
function Smooth(x, x0, deltax)::Float64
    return 0.5 * (x + x0 + sqrt((x - x0) * (x - x0) + 0.25 * deltax * deltax))
end

# Smoothing function for (min of x, x0 with deltax)
function Smooth2(x, x0, deltax)::Float64
    return 0.5 * (x + x0 - sqrt((x - x0) * (x - x0) + 0.25 * deltax * deltax)) + 0.25 * deltax
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

# Junction capacitance macro between S/D and bulk
function JunCap(Czbx, Vbx_jct, PBX_t, MJX, czbx_p1, czbx_p2)::Float64
    if (Czbx > 0.0)
        T1 = Vbx_jct / PBX_t
        if (T1 < 0.9)
            arg = 1.0 - T1
            if (MJX == 0.5)
                sarg = 1.0 / sqrt(arg)
            else
                sarg = lexp(-MJX * lln(arg))
            end
            return PBX_t * Czbx * (1.0 - arg * sarg) / (1.0 - MJX)
        else
            T2 = czbx_p1 * (T1 - 1.0) * (5.0 * MJX * (T1 - 1.0) + (1.0 + MJX))
            return PBX_t * Czbx * (T2 + czbx_p2)
        end
    else
        return 0.0
	end
end

# Normalized pinch-off voltage including PD
function PO_psip(vg_vfb, gamma, DPD, phif)::Float64
	M_SQRT2::Float64 = 1.41421356237309504880
    T1       = 1.0 + DPD
    vgfbPD   = vg_vfb / T1
    gammaPD  = gamma / T1
    T1       = 0.5 * vgfbPD - 3.0 * (1.0 + gammaPD / M_SQRT2)
    T2       = T1 + sqrt(T1 * T1 + 6.0 * vgfbPD)
    if (vgfbPD < 0.0)
        T3   = (vgfbPD - T2) / gammaPD
        return -lln(1.0 - T2 + T3 * T3)
    else
        T3   = lexp(-T2)
        T1   = 0.5 * gammaPD
        T2   = sqrt(vgfbPD - 1.0 + T3 + T1 * T1) - T1
        return T2 * T2 + 1.0 - T3
    end
end

# Normalized charge-voltage relationship
function BSIM_q(psip, phib, vch, gam)::Float64
    T8 = 0.5 * (psip + 1.0 + sqrt((psip - 1.0) * (psip - 1.0) + 0.25 * 2.0 * 2.0))
    sqrtpsip = sqrt(T8)
    T9 = 1.0 + gam / (2.0 * sqrtpsip)
    T0 = (1.0 + (gam / (2.0 * sqrtpsip))) / gam
    T1 = psip - 2.0 * phib - vch
    T2 = T1 - lln(4.0 * T0 * sqrtpsip)
    T8 = 0.5 * (T2 - 0.201491 - sqrt(T2 * (T2 + 0.402982) + 2.446562))
    sqrtpsisa = sqrtpsip
    if (T8 <= -68.0)
        T4 = -100.0
        T5 = 20.0
        if (T8 < T4 - 0.5 * T5)
            T3 = lexp(T4)
        else
            if (T8 > T4 + 0.5 * T5)
                T3 = lexp(T8)
            else
                T2 = (T8 - T4) / T5
                T6 = T2 * T2
                T3 = lexp(T4 + T5 * ((5.0 / 64.0) + 0.5 * T2 + T6 * ((15.0 / 16.0) - T6 * (1.25 - T6))))
            end
        end
        return T3 * (1.0 + T1 - T8 - lln(2.0 * T0 * (T3 * 2.0 * T0 + 2.0 * sqrtpsisa)))
    else
        T3 = lexp(T8)
        sqrtpsisainv = 1.0 / sqrtpsisa
        T4 = 2.0 * T3 + lln(T3 * 2.0 * T0 * (T3 *  2.0 * T0 + 2.0 * sqrtpsisa)) - T1
        T5 = 2.0 + (1.0 / T3) + (T0 + sqrtpsisainv) / (T0 * T3 + sqrtpsisa)
        T3 = T3 - T4 / T5
        T4 = 2.0 * T3 + lln(T3 * 2.0 * T0 * (T3 * 2.0 * T0 + 2.0 * sqrtpsisa)) - T1
        T5 = 2.0 + (1.0 / T3) + (T0 + sqrtpsisainv) / (T0 * T3 + sqrtpsisa)
        T6 = ((T0 + sqrtpsisainv) / (T0 * T3 + sqrtpsisa)) * ((T0 + sqrtpsisainv) / (T0 * T3 + sqrtpsisa))
        T7 = -((1.0 / T3) * (1.0 / T3)) - (1.0 / (sqrtpsisa * sqrtpsisa * sqrtpsisa * (T0 * T3 + sqrtpsisa))) - T6
        return T3 - (T4 / T5) * (1.0 + T4 * T7 / (2.0 * T5 * T5))
    end
end

# Define GEOMOD and RGEOMOD in the modelcard
function BSIMBULKNumFingerDiff(nf, minSD)
    if (nf % 2) != 0
        nuEndD = 1
        nuIntD = 2 * max((nf - 1) / 2, 0)
        nuEndS = 1
        nuIntS = nuIntD
    else
        if (minSD == 1)
            nuEndD = 2
            nuIntD = 2 * max((nf / 2 - 1), 0.0)
            nuEndS = 0
            nuIntS = nf
        else
            nuEndD = 0
            nuIntD = nf
            nuEndS = 2
            nuIntS = 2 * max((nf / 2 - 1), 0)
        end
    end
	return nuEndD, nuIntD, nuEndS, nuIntS
end

function BSIMBULKRdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEnd, rgeo, SRCFLAG)::Float64
    if (SRCFLAG == 1)
		if (rgeo == 1) || (rgeo == 2) || (rgeo == 5)
            if (nuEnd == 0)
                Rend = 0.0
            else
                Rend = Rsh * DMCG / (Weffcj * nuEnd)
            end
		elseif (rgeo == 3) || (rgeo == 4) || (rgeo == 6)
            if ((DMCG + DMCI) == 0.0)
                println("(DMCG + DMCI) can not be equal to zero")
            end
            if (nuEnd == 0)
                Rend = 0.0
            else
                Rend = Rsh * Weffcj / (3.0 * nuEnd * (DMCG + DMCI))
            end
		else
            println("Warning: (instance %M) Specified RGEO = %d not matched (BSIMBULKRdsEndIso), Rend is set to zero.", rgeo)
            Rend = 0.0
        end
    else
		if (rgeo == 1) || (rgeo == 3) || (rgeo == 7)
            if (nuEnd == 0)
                Rend = 0.0
            else
                Rend = Rsh * DMCG / (Weffcj * nuEnd)
            end
		elseif (rgeo == 2) || (rgeo == 4) || (rgeo == 8)
            if ((DMCG + DMCI) == 0.0)
                println("(DMCG + DMCI) can not be equal to zero")
            end
            if (nuEnd == 0)
                Rend = 0.0
            else
                Rend = Rsh * Weffcj / (3.0 * nuEnd * (DMCG + DMCI))
            end
		else
			println("Warning: (instance %M) Specified RGEO=%d not matched (BSIMBULKRdsEndIso type 2), Rend is set to zero.", rgeo)
            Rend = 0.0
        end
    end
	return Rend
end

function BSIMBULKRdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEnd, rgeo, SRCFLAG)::Float64
	if (SRCFLAG == 1)
		if (rgeo == 1) || (rgeo == 2) || (rgeo == 5)
			if (nuEnd == 0)
				Rend = 0.0
			else
				Rend = Rsh * DMCG / (Weffcj * nuEnd)
			end
		elseif (rgeo == 3) || (rgeo == 4) || (rgeo == 6)
			if (DMCG == 0.0)
				println("DMCG can not be equal to zero")
			end
			if (nuEnd == 0)
				Rend = 0.0
			else
				Rend = Rsh * Weffcj / (6.0 * nuEnd * DMCG)
			end
		else
			println("Warning: (instance %M) Specified RGEO = %d not matched (BSIMBULKRdsEndSha), Rend is set to zero.", rgeo)
			Rend = 0.0
		end
	else
		if (rgeo == 1) || (rgeo == 3) || (rgeo == 7)
			if (nuEnd == 0)
				Rend = 0.0
			else
				Rend = Rsh * DMCG / (Weffcj * nuEnd)
			end
		elseif (rgeo == 2) || (rgeo == 4) || (rgeo == 8)
			if (DMCG == 0.0)
				println("DMCG can not be equal to zero")
			end
			if (nuEnd == 0)
				Rend = 0.0
			else
				Rend = Rsh * Weffcj / (6.0 * nuEnd * DMCG)
			end
		else
			println("Warning: (instance %M) Specified RGEO=%d not matched (BSIMBULKRdsEndSha type 2), Rend is set to zero.", rgeo)
			Rend = 0.0
		end
	end
	return Rend
end

function BSIMBULKRdseffGeo(nf, geo, rgeo, minSD, Weffcj, Rsh, DMCG, DMCI, DMDG, SRCFLAG)::Float64
	if (geo < 9)
		nuEndD, nuIntD, nuEndS, nuIntS = BSIMBULKNumFingerDiff(nf, minSD)
		if (SRCFLAG == 1)
			if (nuIntS == 0)
				Rint = 0.0
			else
				Rint = Rsh * DMCG / ( Weffcj * nuIntS)
			end
		else
			if (nuIntD == 0)
				Rint = 0.0
			else
				Rint = Rsh * DMCG / ( Weffcj * nuIntD)
			end
		end
	end
	if (geo == 0)
		if (SRCFLAG == 1)
			Rend = BSIMBULKRdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndS, rgeo, 1)
		else
			Rend = BSIMBULKRdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndD, rgeo, 0)
		end
	elseif (geo == 1)
		if (SRCFLAG == 1)
			Rend = BSIMBULKRdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndS, rgeo, 1)
		else
			Rend = BSIMBULKRdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndD, rgeo, 0)
		end
	elseif (geo == 2)
		if (SRCFLAG == 1)
			Rend = BSIMBULKRdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndS, rgeo, 1)
		else
			Rend = BSIMBULKRdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndD, rgeo, 0)
		end
	elseif (geo == 3)
		if (SRCFLAG == 1)
			Rend = BSIMBULKRdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndS, rgeo, 1)
		else
			Rend = BSIMBULKRdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndD, rgeo, 0)
		end
	elseif (geo == 4)
		if (SRCFLAG == 1)
			Rend = BSIMBULKRdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndS, rgeo, 1)
		else
			Rend = Rsh * DMDG / Weffcj
		end
	elseif (geo == 5)
		if (SRCFLAG == 1)
			Rend = BSIMBULKRdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndS, rgeo, 1)
		else
			if (nuEndD == 0)
				Rend = 0.0
			else
				Rend = Rsh * DMDG / (Weffcj * nuEndD)
			end
		end
	elseif (geo == 6)
		if (SRCFLAG == 1)
			Rend = Rsh * DMDG / Weffcj
		else
			Rend = BSIMBULKRdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndD, rgeo, 0)
		end
	elseif (geo == 7)
		if (SRCFLAG == 1)
			if (nuEndS == 0)
				Rend = 0.0
			else
				Rend = Rsh * DMDG / (Weffcj * nuEndS)
			end
		else
			Rend = BSIMBULKRdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndD, rgeo, 0)
		end
	elseif (geo == 8)
		Rend = Rsh * DMDG / Weffcj
	elseif (geo == 9)		# all wide contacts assumed for geo = 9 and 10
		if (SRCFLAG == 1)
			Rend = 0.5 * Rsh * DMCG / Weffcj
			if (nf == 2)
				Rint = 0.0
			else
				Rint = Rsh * DMCG / (Weffcj * (nf - 2.0))
			end
		else
			Rend = 0.0
			Rint = Rsh * DMCG / (Weffcj * nf)
		end
	elseif (geo == 10)
		if (SRCFLAG == 1)
			Rend = 0.0
			Rint = Rsh * DMCG / (Weffcj * nf)
		else
			Rend = 0.5 * Rsh * DMCG / Weffcj
			if (nf == 2)
				Rint = 0.0
			else
				Rint = Rsh * DMCG / (Weffcj * (nf - 2.0))
			end
		end
	else
		println("Warning: (instance %M) Specified GEO=%d not matched (BSIMBULKRdseffGeo), Rint is set to zero.", geo)
		Rint = 0.0
	end
	if (Rint <= 0.0)
		Rtot = Rend
	elseif (Rend <= 0.0)
		Rtot = Rint
	else
		Rtot = Rint * Rend / (Rint + Rend)
	end
	if (Rtot == 0.0)
		println("Warning: (instance %M) Zero resistance returned from RdseffGeo")
	end
	return Rtot
end

# Effective PS, PD, AS, AD calculation, Ref: BSIM4
function BSIMBULKPAeffGeo(nf, geo, minSD, Weffcj, DMCG, DMCI, DMDG)
    if (geo < 9)
		nuEndD, nuIntD, nuEndS, nuIntS = BSIMBULKNumFingerDiff(nf, minSD)
		T0 = DMCG + DMCI
		T1 = DMCG + DMCG
		T2 = DMDG + DMDG
		PSiso = T0 + T0 + Weffcj
		PDiso = T0 + T0 + Weffcj
		PSsha = T1
		PDsha = T1
		PSmer = T2
		PDmer = T2
		ASiso = T0 * Weffcj
		ADiso = T0 * Weffcj
		ASsha = DMCG * Weffcj
		ADsha = DMCG * Weffcj
		ASmer = DMDG * Weffcj
		ADmer = DMDG * Weffcj
	end
	if (geo == 0)
		Ps = nuEndS * PSiso + nuIntS * PSsha
		Pd = nuEndD * PDiso + nuIntD * PDsha
		As = nuEndS * ASiso + nuIntS * ASsha
		Ad = nuEndD * ADiso + nuIntD * ADsha
	elseif (geo == 1)
		Ps = nuEndS * PSiso + nuIntS * PSsha
		Pd = (nuEndD + nuIntD) * PDsha
		As = nuEndS * ASiso + nuIntS * ASsha
		Ad = (nuEndD + nuIntD) * ADsha
	elseif (geo == 2)
		Ps = (nuEndS + nuIntS) * PSsha
		Pd = nuEndD * PDiso + nuIntD * PDsha
		As = (nuEndS + nuIntS) * ASsha
		Ad = nuEndD * ADiso + nuIntD * ADsha
	elseif (geo == 3)
		Ps = (nuEndS + nuIntS) * PSsha
		Pd = (nuEndD + nuIntD) * PDsha
		As = (nuEndS + nuIntS) * ASsha
		Ad = (nuEndD + nuIntD) * ADsha
	elseif (geo == 4)
		Ps = nuEndS * PSiso + nuIntS * PSsha
		Pd = nuEndD * PDmer + nuIntD * PDsha
		As = nuEndS * ASiso + nuIntS * ASsha
		Ad = nuEndD * ADmer + nuIntD * ADsha
	elseif (geo == 5)
		Ps = (nuEndS + nuIntS) * PSsha
		Pd = nuEndD * PDmer + nuIntD * PDsha
		As = (nuEndS + nuIntS) * ASsha
		Ad = nuEndD * ADmer + nuIntD * ADsha
	elseif (geo == 6)
		Ps = nuEndS * PSmer + nuIntS * PSsha
		Pd = nuEndD * PDiso + nuIntD * PDsha
		As = nuEndS * ASmer + nuIntS * ASsha
		Ad = nuEndD * ADiso + nuIntD * ADsha
	elseif (geo == 7)
		Ps = nuEndS * PSmer + nuIntS * PSsha
		Pd = (nuEndD + nuIntD) * PDsha
		As = nuEndS * ASmer + nuIntS * ASsha
		Ad = (nuEndD + nuIntD) * ADsha
	elseif (geo == 8)
		Ps = nuEndS * PSmer + nuIntS * PSsha
		Pd = nuEndD * PDmer + nuIntD * PDsha
		As = nuEndS * ASmer + nuIntS * ASsha
		Ad = nuEndD * ADmer + nuIntD * ADsha
	elseif (geo == 9)
		Ps = PSiso + (nf - 1.0) * PSsha
		Pd = nf * PDsha
		As = ASiso + (nf - 1.0) * ASsha
		Ad = nf * ADsha
	elseif (geo == 10)
		Ps = nf * PSsha
		Pd = PDiso + (nf - 1.0) * PDsha
		As = nf * ASsha
		Ad = ADiso + (nf - 1.0) * ADsha
	else
		println("Warning: (instance %M) Specified GEO=%d not matched (BSIMBULKPAeffGeo), PS,PD,AS,AD set to zero.", geo)
		Ps = 0
		Pd = 0
		As = 0
		Ad = 0
	end
	return Ps, Pd, As, Ad
end

function bsimbulk(param)
	ntype::Int8 = 1
	ptype::Int8 = -1
	q::Float64 = 1.60219e-19
	EPS0::Float64 = 8.85418e-12
	KboQ::Float64 = 8.617087e-5      # Joule/degree
	P_CELSIUS0::Float64 = 273.15
	Oneby3::Float64 = 0.33333333333333333
	DELTA_1::Float64 = 0.02
	REFTEMP::Float64 = 300.15       # 27 degrees C
	EXPL_THRESHOLD::Float64 = 80.0
	MAX_EXPL::Float64 = 5.540622384e34
	M_PI::Float64 = 3.14159265358979323846

	# Bias and temperature
	Vd::Float64 = get(param, "Vd", 1.0)
	Vg::Float64 = get(param, "Vg", 1.0)
	Vs::Float64 = get(param, "Vs", 0.0)
	Vb::Float64 = get(param, "Vb", 0.0)
	Temp::Float64 = get(param, "Temp", 25.0)

	# Pure instance parameters
	L::Float64 = get(param, "L", 1.0e-5)
	W::Float64 = get(param, "W", 1.0e-5)
	NF::UInt8 = get(param, "NF", 1)
	NRS::Float64 = get(param, "NRS", 1.0)
	NRD::Float64 = get(param, "NRD", 1.0)
	VFBSDOFF::Float64 = get(param, "VFBSDOFF", 0.0)
	MINZ::UInt8 = get(param, "MINZ", 0)
	RGATEMOD::UInt8 = get(param, "RGATEMOD", 0)
	RBODYMOD::UInt8 = get(param, "RBODYMOD", 0)
	GEOMOD::UInt8 = get(param, "GEOMOD", 0)
	RGEOMOD::UInt8 = get(param, "RGEOMOD", 0)
	RBPB::Float64 = get(param, "RBPB", 50.0)
	RBPD::Float64 = get(param, "RBPD", 50.0)
	RBPS::Float64 = get(param, "RBPS", 50.0)
	RBDB::Float64 = get(param, "RBDB", 50.0)
	RBSB::Float64 = get(param, "RBSB", 50.0)
	SA::Float64 = get(param, "SA", 0.0)
	SB::Float64 = get(param, "SB", 0.0)
	SD::Float64 = get(param, "SD", 0.0)
	SCA::Float64 = get(param, "SCA", 0.0)
	SCB::Float64 = get(param, "SCB", 0.0)
	SCC::Float64 = get(param, "SCC", 0.0)
	SC::Float64 = get(param, "SC", 0.0)
	AS::Float64 = get(param, "AS", 0.0)
	AD::Float64 = get(param, "AD", 0.0)
	PS::Float64 = get(param, "PS", 0.0)
	PD::Float64 = get(param, "PD", 0.0)

	# Both model and instance parameters
	XGW::Float64 = get(param, "XGW", 0.0)
	NGCON::UInt8 = get(param, "NGCON", 1)
	DTEMP::Float64 = get(param, "DTEMP", 0.0)
	MULU0::Float64 = get(param, "MULU0", 1.0)
	DELVTO::Float64 = get(param, "DELVTO", 0.0)
	IDS0MULT::Float64 = get(param, "IDS0MULT", 1.0)
	EDGEFET::UInt8 = get(param, "EDGEFET", 0)
	SSLMOD::UInt8 = get(param, "SSLMOD", 0)

	# Pure model parameters
	TYPE::UInt8 = get(param, "TYPE", ntype)
	CVMOD::UInt8 = get(param, "CVMOD", 0)
	COVMOD::UInt8 = get(param, "COVMOD", 0)
	RDSMOD::UInt8 = get(param, "RDSMOD", 0)
	WPEMOD::UInt8 = get(param, "WPEMOD", 0)
	ASYMMOD::UInt8 = get(param, "ASYMMOD", 0)
	GIDLMOD::UInt8 = get(param, "GIDLMOD", 0)
	IGCMOD::UInt8 = get(param, "IGCMOD", 0)
	IGBMOD::UInt8 = get(param, "IGBMOD", 0)
	TNOIMOD::UInt8 = get(param, "TNOIMOD", 0)
	SHMOD::UInt8 = get(param, "SHMOD", 0)
	MOBSCALE::UInt8 = get(param, "MOBSCALE", 0)

	# Device parameters
	LLONG::Float64 = get(param, "LLONG", 1.0e-5)
	LMLT::Float64 = get(param, "LMLT", 1.0)
	WMLT::Float64 = get(param, "WMLT", 1.0)
	XL::Float64 = get(param, "XL", 0.0)
	WWIDE::Float64 = get(param, "WWIDE", 1.0e-5)
	XW::Float64 = get(param, "XW", 0.0)
	LINT::Float64 = get(param, "LINT", 0.0)
	LL::Float64 = get(param, "LL", 0.0)
	LW::Float64 = get(param, "LW", 0.0)
	LWL::Float64 = get(param, "LWL", 0.0)
	LLN::Float64 = get(param, "LLN", 1.0)
	LWN::Float64 = get(param, "LWN", 1.0)
	WINT::Float64 = get(param, "WINT", 0.0)
	WL::Float64 = get(param, "WL", 0.0)
	WW::Float64 = get(param, "WW", 0.0)
	WWL::Float64 = get(param, "WWL", 0.0)
	WLN::Float64 = get(param, "WLN", 1.0)
	WWN::Float64 = get(param, "WWN", 1.0)
	DLC::Float64 = get(param, "DLC", 0.0)
	LLC::Float64 = get(param, "LLC", 0.0)
	LWC::Float64 = get(param, "LWC", 0.0)
	LWLC::Float64 = get(param, "LWLC", 0.0)
	DWC::Float64 = get(param, "DWC", 0.0)
	WLC::Float64 = get(param, "WLC", 0.0)
	WWC::Float64 = get(param, "WWC", 0.0)
	WWLC::Float64 = get(param, "WWLC", 0.0)
	TOXE::Float64 = get(param, "TOXE", 3.0e-9)
	TOXP::Float64 = get(param, "TOXP", TOXE)
	DTOX::Float64 = get(param, "DTOX", 0.0)
	NDEP::Float64 = get(param, "NDEP", 1e24)
	NDEPL1::Float64 = get(param, "NDEPL1", 0.0)
	NDEPLEXP1::Float64 = get(param, "NDEPLEXP1", 1.0)
	NDEPL2::Float64 = get(param, "NDEPL2", 0.0)
	NDEPLEXP2::Float64 = get(param, "NDEPLEXP2", 2.0)
	NDEPW::Float64 = get(param, "NDEPW", 0.0)
	NDEPWEXP::Float64 = get(param, "NDEPWEXP", 1.0)
	NDEPWL::Float64 = get(param, "NDEPWL", 0.0)
	NDEPWLEXP::Float64 = get(param, "NDEPWLEXP", 1.0)
	LNDEP::Float64 = get(param, "LNDEP", 0.0)
	WNDEP::Float64 = get(param, "WNDEP", 0.0)
	PNDEP::Float64 = get(param, "PNDEP", 0.0)
	NDEPCV::Float64 = get(param, "NDEPCV", 1e24)
	NDEPCVL1::Float64 = get(param, "NDEPCVL1", 0.0)
	NDEPCVLEXP1::Float64 = get(param, "NDEPCVLEXP1", 1.0)
	NDEPCVL2::Float64 = get(param, "NDEPCVL2", 0.0)
	NDEPCVLEXP2::Float64 = get(param, "NDEPCVLEXP2", 2.0)
	NDEPCVW::Float64 = get(param, "NDEPCVW", 0.0)
	NDEPCVWEXP::Float64 = get(param, "NDEPCVWEXP", 1.0)
	NDEPCVWL::Float64 = get(param, "NDEPCVWL", 0.0)
	NDEPCVWLEXP::Float64 = get(param, "NDEPCVWLEXP", 1.0)
	LNDEPCV::Float64 = get(param, "LNDEPCV", 0.0)
	WNDEPCV::Float64 = get(param, "WNDEPCV", 0.0)
	PNDEPCV::Float64 = get(param, "PNDEPCV", 0.0)
	NGATE::Float64 = get(param, "NGATE", 5e25)
	LNGATE::Float64 = get(param, "LNGATE", 0.0)
	WNGATE::Float64 = get(param, "WNGATE", 0.0)
	PNGATE::Float64 = get(param, "PNGATE", 0.0)
	EASUB::Float64 = get(param, "EASUB", 4.05)
	NI0SUB::Float64 = get(param, "NI0SUB", 1.1e16)
	BG0SUB::Float64 = get(param, "BG0SUB", 1.17)
	EPSRSUB::Float64 = get(param, "EPSRSUB", 11.9)
	EPSROX::Float64 = get(param, "EPSROX", 3.9)
	XJ::Float64 = get(param, "XJ", 1.5e-7)
	LXJ::Float64 = get(param, "LXJ", 0.0)
	WXJ::Float64 = get(param, "WXJ", 0.0)
	PXJ::Float64 = get(param, "PXJ", 0.0)
	VFB::Float64 = get(param, "VFB", -0.5)
	LVFB::Float64 = get(param, "LVFB", 0.0)
	WVFB::Float64 = get(param, "WVFB", 0.0)
	PVFB::Float64 = get(param, "PVFB", 0.0)
	VFBCV::Float64 = get(param, "VFBCV", -0.5)
	LVFBCV::Float64 = get(param, "LVFBCV", 0.0)
	WVFBCV::Float64 = get(param, "WVFBCV", 0.0)
	PVFBCV::Float64 = get(param, "PVFBCV", 0.0)
	VFBCVL::Float64 = get(param, "VFBCVL", 0.0)
	VFBCVLEXP::Float64 = get(param, "VFBCVLEXP", 1.0)
	VFBCVW::Float64 = get(param, "VFBCVW", 0.0)
	VFBCVWEXP::Float64 = get(param, "VFBCVWEXP", 1.0)
	VFBCVWL::Float64 = get(param, "VFBCVWL", 0.0)
	VFBCVWLEXP::Float64 = get(param, "VFBCVWLEXP", 1.0)
	DELVFBACC::Float64 = get(param, "DELVFBACC", 0.0)

	# Diode parameters
	PERMOD::UInt8 = get(param, "PERMOD", 1)
	DWJ::Float64 = get(param, "DWJ", DWC)

	# Short channel effects
	NSD::Float64 = get(param, "NSD", 1e26)
	LNSD::Float64 = get(param, "LNSD", 0.0)
	WNSD::Float64 = get(param, "WNSD", 0.0)
	PNSD::Float64 = get(param, "PNSD", 0.0)
	DVTP0::Float64 = get(param, "DVTP0", 0.0)
	LDVTP0::Float64 = get(param, "LDVTP0", 0)
	WDVTP0::Float64 = get(param, "WDVTP0", 0)
	PDVTP0::Float64 = get(param, "PDVTP0", 0)
	DVTP1::Float64 = get(param, "DVTP1", 0.0)
	LDVTP1::Float64 = get(param, "LDVTP1", 0)
	WDVTP1::Float64 = get(param, "WDVTP1", 0)
	PDVTP1::Float64 = get(param, "PDVTP1", 0)
	DVTP2::Float64 = get(param, "DVTP2", 0.0)
	LDVTP2::Float64 = get(param, "LDVTP2", 0)
	WDVTP2::Float64 = get(param, "WDVTP2", 0)
	PDVTP2::Float64 = get(param, "PDVTP2", 0)
	DVTP3::Float64 = get(param, "DVTP3", 0.0)
	LDVTP3::Float64 = get(param, "LDVTP3", 0)
	WDVTP3::Float64 = get(param, "WDVTP3", 0)
	PDVTP3::Float64 = get(param, "PDVTP3", 0)
	DVTP4::Float64 = get(param, "DVTP4", 0.0)
	LDVTP4::Float64 = get(param, "LDVTP4", 0)
	WDVTP4::Float64 = get(param, "WDVTP4", 0)
	PDVTP4::Float64 = get(param, "PDVTP4", 0)
	DVTP5::Float64 = get(param, "DVTP5", 0.0)
	LDVTP5::Float64 = get(param, "LDVTP5", 0)
	WDVTP5::Float64 = get(param, "WDVTP5", 0)
	PDVTP5::Float64 = get(param, "PDVTP5", 0)
	PHIN::Float64 = get(param, "PHIN", 0.045)
	LPHIN::Float64 = get(param, "LPHIN", 0.0)
	WPHIN::Float64 = get(param, "WPHIN", 0.0)
	PPHIN::Float64 = get(param, "PPHIN", 0.0)
	ETA0::Float64 = get(param, "ETA0", 0.08)
	LETA0::Float64 = get(param, "LETA0", 0.0)
	WETA0::Float64 = get(param, "WETA0", 0.0)
	PETA0::Float64 = get(param, "PETA0", 0.0)
	ETA0R::Float64 = get(param, "ETA0R", ETA0)
	LETA0R::Float64 = get(param, "LETA0R", LETA0)
	WETA0R::Float64 = get(param, "WETA0R", WETA0)
	PETA0R::Float64 = get(param, "PETA0R", PETA0)
	DSUB::Float64 = get(param, "DSUB", 1.0)
	ETAB::Float64 = get(param, "ETAB", -0.07)
	ETABEXP::Float64 = get(param, "ETABEXP", 1.0)
	LETAB::Float64 = get(param, "LETAB", 0.0)
	WETAB::Float64 = get(param, "WETAB", 0.0)
	PETAB::Float64 = get(param, "PETAB", 0.0)
	K1::Float64 = get(param, "K1", 0.0)
	K1L::Float64 = get(param, "K1L", 0.0)
	K1LEXP::Float64 = get(param, "K1LEXP", 1.0)
	K1W::Float64 = get(param, "K1W", 0.0)
	K1WEXP::Float64 = get(param, "K1WEXP", 1.0)
	K1WL::Float64 = get(param, "K1WL", 0.0)
	K1WLEXP::Float64 = get(param, "K1WLEXP", 1.0)
	LK1::Float64 = get(param, "LK1", 0.0)
	WK1::Float64 = get(param, "WK1", 0.0)
	PK1::Float64 = get(param, "PK1", 0.0)
	K2::Float64 = get(param, "K2", 0.0)
	K2L::Float64 = get(param, "K2L", 0.0)
	K2LEXP::Float64 = get(param, "K2LEXP", 1.0)
	K2W::Float64 = get(param, "K2W", 0.0)
	K2WEXP::Float64 = get(param, "K2WEXP", 1.0)
	K2WL::Float64 = get(param, "K2WL", 0.0)
	K2WLEXP::Float64 = get(param, "K2WLEXP", 1.0)
	LK2::Float64 = get(param, "LK2", 0.0)
	WK2::Float64 = get(param, "WK2", 0.0)
	PK2::Float64 = get(param, "PK2", 0.0)

	# Quantum mechanical effects
	ADOS::Float64 = get(param, "ADOS", 0.0)
	BDOS::Float64 = get(param, "BDOS", 1.0)
	QM0::Float64 = get(param, "QM0", 1.0e-3)
	ETAQM::Float64 = get(param, "ETAQM", 0.54)

	# Sub-threshold swing factor
	CIT::Float64 = get(param, "CIT", 0.0)
	LCIT::Float64 = get(param, "LCIT", 0.0)
	WCIT::Float64 = get(param, "WCIT", 0.0)
	PCIT::Float64 = get(param, "PCIT", 0.0)
	NFACTOR::Float64 = get(param, "NFACTOR", 0.0)
	NFACTORL::Float64 = get(param, "NFACTORL", 0.0)
	NFACTORLEXP::Float64 = get(param, "NFACTORLEXP", 1.0)
	NFACTORW::Float64 = get(param, "NFACTORW", 0.0)
	NFACTORWEXP::Float64 = get(param, "NFACTORWEXP", 1.0)
	NFACTORWL::Float64 = get(param, "NFACTORWL", 0.0)
	NFACTORWLEXP::Float64 = get(param, "NFACTORWLEXP", 1.0)
	LNFACTOR::Float64 = get(param, "LNFACTOR", 0.0)
	WNFACTOR::Float64 = get(param, "WNFACTOR", 0.0)
	PNFACTOR::Float64 = get(param, "PNFACTOR", 0.0)
	CDSCD::Float64 = get(param, "CDSCD", 1e-9)
	CDSCDL::Float64 = get(param, "CDSCDL", 0.0)
	CDSCDLEXP::Float64 = get(param, "CDSCDLEXP", 1.0)
	LCDSCD::Float64 = get(param, "LCDSCD", 0.0)
	WCDSCD::Float64 = get(param, "WCDSCD", 0.0)
	PCDSCD::Float64 = get(param, "PCDSCD", 0.0)
	CDSCDR::Float64 = get(param, "CDSCDR", CDSCD)
	CDSCDLR::Float64 = get(param, "CDSCDLR", CDSCDL)
	LCDSCDR::Float64 = get(param, "LCDSCDR", LCDSCD)
	WCDSCDR::Float64 = get(param, "WCDSCDR", WCDSCD)
	PCDSCDR::Float64 = get(param, "PCDSCDR", PCDSCD)
	CDSCB::Float64 = get(param, "CDSCB", 0.0)
	CDSCBL::Float64 = get(param, "CDSCBL", 0.0)
	CDSCBLEXP::Float64 = get(param, "CDSCBLEXP", 1.0)
	LCDSCB::Float64 = get(param, "LCDSCB", 0.0)
	WCDSCB::Float64 = get(param, "WCDSCB", 0.0)
	PCDSCB::Float64 = get(param, "PCDSCB", 0.0)

	# Drain saturation voltage
	VSAT::Float64 = get(param, "VSAT", 1e5)
	LVSAT::Float64 = get(param, "LVSAT", 0.0)
	WVSAT::Float64 = get(param, "WVSAT", 0.0)
	PVSAT::Float64 = get(param, "PVSAT", 0.0)
	VSATL::Float64 = get(param, "VSATL", 0.0)
	VSATLEXP::Float64 = get(param, "VSATLEXP", 1.0)
	VSATW::Float64 = get(param, "VSATW", 0.0)
	VSATWEXP::Float64 = get(param, "VSATWEXP", 1.0)
	VSATWL::Float64 = get(param, "VSATWL", 0.0)
	VSATWLEXP::Float64 = get(param, "VSATWLEXP", 1.0)
	VSATR::Float64 = get(param, "VSATR", VSAT)
	LVSATR::Float64 = get(param, "LVSATR", LVSAT)
	WVSATR::Float64 = get(param, "WVSATR", WVSAT)
	PVSATR::Float64 = get(param, "PVSATR", PVSAT)
	DELTA::Float64 = get(param, "DELTA", 0.125)
	LDELTA::Float64 = get(param, "LDELTA", 0.0)
	WDELTA::Float64 = get(param, "WDELTA", 0.0)
	PDELTA::Float64 = get(param, "PDELTA", 0.0)
	DELTAL::Float64 = get(param, "DELTAL", 0.0)
	DELTALEXP::Float64 = get(param, "DELTALEXP", 1.0)
	VSATCV::Float64 = get(param, "VSATCV", 1e5)
	LVSATCV::Float64 = get(param, "LVSATCV", 0.0)
	WVSATCV::Float64 = get(param, "WVSATCV", 0.0)
	PVSATCV::Float64 = get(param, "PVSATCV", 0.0)
	VSATCVL::Float64 = get(param, "VSATCVL", 0.0)
	VSATCVLEXP::Float64 = get(param, "VSATCVLEXP", 1.0)
	VSATCVW::Float64 = get(param, "VSATCVW", 0.0)
	VSATCVWEXP::Float64 = get(param, "VSATCVWEXP", 1.0)
	VSATCVWL::Float64 = get(param, "VSATCVWL", 0.0)
	VSATCVWLEXP::Float64 = get(param, "VSATCVWLEXP", 1.0)

	# Mobility degradation
	UP1::Float64 = get(param, "UP1", 0.0)
	LP1::Float64 = get(param, "LP1", 1.0e-8)
	UP2::Float64 = get(param, "UP2", 0.0)
	LP2::Float64 = get(param, "LP2", 1.0e-8)
	U0::Float64 = get(param, "U0", 67.0e-3)
	U0L::Float64 = get(param, "U0L", 0.0)
	U0LEXP::Float64 = get(param, "U0LEXP", 1.0)
	LU0::Float64 = get(param, "LU0", 0.0)
	WU0::Float64 = get(param, "WU0", 0.0)
	PU0::Float64 = get(param, "PU0", 0.0)
	U0R::Float64 = get(param, "U0R", U0)
	LU0R::Float64 = get(param, "LU0R", LU0)
	WU0R::Float64 = get(param, "WU0R", WU0)
	PU0R::Float64 = get(param, "PU0R", PU0)
	ETAMOB::Float64 = get(param, "ETAMOB", 1.0)
	UA::Float64 = get(param, "UA", 0.001)
	UAL::Float64 = get(param, "UAL", 0.0)
	UALEXP::Float64 = get(param, "UALEXP", 1.0)
	UAW::Float64 = get(param, "UAW", 0.0)
	UAWEXP::Float64 = get(param, "UAWEXP", 1.0)
	UAWL::Float64 = get(param, "UAWL", 0.0)
	UAWLEXP::Float64 = get(param, "UAWLEXP", 1.0)
	LUA::Float64 = get(param, "LUA", 0.0)
	WUA::Float64 = get(param, "WUA", 0.0)
	PUA::Float64 = get(param, "PUA", 0.0)
	UAR::Float64 = get(param, "UAR", UA)
	LUAR::Float64 = get(param, "LUAR", LUA)
	WUAR::Float64 = get(param, "WUAR", WUA)
	PUAR::Float64 = get(param, "PUAR", PUA)
	EU::Float64 = get(param, "EU", 1.5)
	LEU::Float64 = get(param, "LEU", 0.0)
	WEU::Float64 = get(param, "WEU", 0.0)
	PEU::Float64 = get(param, "PEU", 0.0)
	EUL::Float64 = get(param, "EUL", 0.0)
	EULEXP::Float64 = get(param, "EULEXP", 1.0)
	EUW::Float64 = get(param, "EUW", 0.0)
	EUWEXP::Float64 = get(param, "EUWEXP", 1.0)
	EUWL::Float64 = get(param, "EUWL", 0.0)
	EUWLEXP::Float64 = get(param, "EUWLEXP", 1.0)
	UD::Float64 = get(param, "UD", 0.001)
	UDL::Float64 = get(param, "UDL", 0.0)
	UDLEXP::Float64 = get(param, "UDLEXP", 1.0)
	LUD::Float64 = get(param, "LUD", 0.0)
	WUD::Float64 = get(param, "WUD", 0.0)
	PUD::Float64 = get(param, "PUD", 0.0)
	UDR::Float64 = get(param, "UDR", UD)
	LUDR::Float64 = get(param, "LUDR", LUD)
	WUDR::Float64 = get(param, "WUDR", WUD)
	PUDR::Float64 = get(param, "PUDR", PUD)
	UCS::Float64 = get(param, "UCS", 2.0)
	LUCS::Float64 = get(param, "LUCS", 0.0)
	WUCS::Float64 = get(param, "WUCS", 0.0)
	PUCS::Float64 = get(param, "PUCS", 0.0)
	UCSR::Float64 = get(param, "UCSR", UCS)
	LUCSR::Float64 = get(param, "LUCSR", LUCS)
	WUCSR::Float64 = get(param, "WUCSR", WUCS)
	PUCSR::Float64 = get(param, "PUCSR", PUCS)
	UC::Float64 = get(param, "UC", 0.0)
	UCL::Float64 = get(param, "UCL", 0.0)
	UCLEXP::Float64 = get(param, "UCLEXP", 1.0)
	UCW::Float64 = get(param, "UCW", 0.0)
	UCWEXP::Float64 = get(param, "UCWEXP", 1.0)
	UCWL::Float64 = get(param, "UCWL", 0.0)
	UCWLEXP::Float64 = get(param, "UCWLEXP", 1.0)
	LUC::Float64 = get(param, "LUC", 0.0)
	WUC::Float64 = get(param, "WUC", 0.0)
	PUC::Float64 = get(param, "PUC", 0.0)
	UCR::Float64 = get(param, "UCR", UC)
	LUCR::Float64 = get(param, "LUCR", LUC)
	WUCR::Float64 = get(param, "WUCR", WUC)
	PUCR::Float64 = get(param, "PUCR", PUC)

	# Channel length modulation
	PCLM::Float64 = get(param, "PCLM", 0.0)
	PCLML::Float64 = get(param, "PCLML", 0.0)
	PCLMLEXP::Float64 = get(param, "PCLMLEXP", 1.0)
	LPCLM::Float64 = get(param, "LPCLM", 0.0)
	WPCLM::Float64 = get(param, "WPCLM", 0.0)
	PPCLM::Float64 = get(param, "PPCLM", 0.0)
	PCLMR::Float64 = get(param, "PCLMR", PCLM)
	LPCLMR::Float64 = get(param, "LPCLMR", LPCLM)
	WPCLMR::Float64 = get(param, "WPCLMR", WPCLM)
	PPCLMR::Float64 = get(param, "PPCLMR", PPCLM)
	PCLMG::Float64 = get(param, "PCLMG", 0.0)
	PCLMCV::Float64 = get(param, "PCLMCV", PCLM)
	PCLMCVL::Float64 = get(param, "PCLMCVL", PCLML)
	PCLMCVLEXP::Float64 = get(param, "PCLMCVLEXP", PCLMLEXP)
	LPCLMCV::Float64 = get(param, "LPCLMCV", LPCLM)
	WPCLMCV::Float64 = get(param, "WPCLMCV", WPCLM)
	PPCLMCV::Float64 = get(param, "PPCLMCV", PPCLM)
	PSCBE1::Float64 = get(param, "PSCBE1", 4.24e8)
	LPSCBE1::Float64 = get(param, "LPSCBE1", 0.0)
	WPSCBE1::Float64 = get(param, "WPSCBE1", 0.0)
	PPSCBE1::Float64 = get(param, "PPSCBE1", 0.0)
	PSCBE2::Float64 = get(param, "PSCBE2", 1.0e-8)
	LPSCBE2::Float64 = get(param, "LPSCBE2", 0.0)
	WPSCBE2::Float64 = get(param, "WPSCBE2", 0.0)
	PPSCBE2::Float64 = get(param, "PPSCBE2", 0.0)
	PDITS::Float64 = get(param, "PDITS", 0.0)
	LPDITS::Float64 = get(param, "LPDITS", 0.0)
	WPDITS::Float64 = get(param, "WPDITS", 0.0)
	PPDITS::Float64 = get(param, "PPDITS", 0.0)
	PDITSL::Float64 = get(param, "PDITSL", 0.0)
	PDITSD::Float64 = get(param, "PDITSD", 0.0)
	LPDITSD::Float64 = get(param, "LPDITSD", 0.0)
	WPDITSD::Float64 = get(param, "WPDITSD", 0.0)
	PPDITSD::Float64 = get(param, "PPDITSD", 0.0)

	# S/D series resistances
	RSH::Float64 = get(param, "RSH", 0.0)
	PRWG::Float64 = get(param, "PRWG", 1.0)
	LPRWG::Float64 = get(param, "LPRWG", 0.0)
	WPRWG::Float64 = get(param, "WPRWG", 0.0)
	PPRWG::Float64 = get(param, "PPRWG", 0.0)
	PRWB::Float64 = get(param, "PRWB", 0.0)
	LPRWB::Float64 = get(param, "LPRWB", 0.0)
	WPRWB::Float64 = get(param, "WPRWB", 0.0)
	PPRWB::Float64 = get(param, "PPRWB", 0.0)
	PRWBL::Float64 = get(param, "PRWBL", 0.0)
	PRWBLEXP::Float64 = get(param, "PRWBLEXP", 1.0)
	WR::Float64 = get(param, "WR", 1.0)
	LWR::Float64 = get(param, "LWR", 0.0)
	WWR::Float64 = get(param, "WWR", 0.0)
	PWR::Float64 = get(param, "PWR", 0.0)
	RSWMIN::Float64 = get(param, "RSWMIN", 0.0)
	LRSWMIN::Float64 = get(param, "LRSWMIN", 0.0)
	WRSWMIN::Float64 = get(param, "WRSWMIN", 0.0)
	PRSWMIN::Float64 = get(param, "PRSWMIN", 0.0)
	RSW::Float64 = get(param, "RSW", 10.0)
	LRSW::Float64 = get(param, "LRSW", 0.0)
	WRSW::Float64 = get(param, "WRSW", 0.0)
	PRSW::Float64 = get(param, "PRSW", 0.0)
	RSWL::Float64 = get(param, "RSWL", 0.0)
	RSWLEXP::Float64 = get(param, "RSWLEXP", 1.0)
	RDWMIN::Float64 = get(param, "RDWMIN", RSWMIN)
	LRDWMIN::Float64 = get(param, "LRDWMIN", LRSWMIN)
	WRDWMIN::Float64 = get(param, "WRDWMIN", WRSWMIN)
	PRDWMIN::Float64 = get(param, "PRDWMIN", PRSWMIN)
	RDW::Float64 = get(param, "RDW", RSW)
	LRDW::Float64 = get(param, "LRDW", LRSW)
	WRDW::Float64 = get(param, "WRDW", WRSW)
	PRDW::Float64 = get(param, "PRDW", PRSW)
	RDWL::Float64 = get(param, "RDWL", RSWL)
	RDWLEXP::Float64 = get(param, "RDWLEXP", RSWLEXP)
	RDSWMIN::Float64 = get(param, "RDSWMIN", 0.0)
	LRDSWMIN::Float64 = get(param, "LRDSWMIN", 0.0)
	WRDSWMIN::Float64 = get(param, "WRDSWMIN", 0.0)
	PRDSWMIN::Float64 = get(param, "PRDSWMIN", 0.0)
	RDSW::Float64 = get(param, "RDSW", 20.0)
	RDSWL::Float64 = get(param, "RDSWL", 0.0)
	RDSWLEXP::Float64 = get(param, "RDSWLEXP", 1.0)
	LRDSW::Float64 = get(param, "LRDSW", 0.0)
	WRDSW::Float64 = get(param, "WRDSW", 0.0)
	PRDSW::Float64 = get(param, "PRDSW", 0.0)

	# Velocity saturation
	PSAT::Float64 = get(param, "PSAT", 1.0)
	LPSAT::Float64 = get(param, "LPSAT", 0.0)
	WPSAT::Float64 = get(param, "WPSAT", 0.0)
	PPSAT::Float64 = get(param, "PPSAT", 0.0)
	PSATL::Float64 = get(param, "PSATL", 0.0)
	PSATLEXP::Float64 = get(param, "PSATLEXP", 1.0)
	PSATB::Float64 = get(param, "PSATB", 0.0)
	PSATR::Float64 = get(param, "PSATR", PSAT)
	LPSATR::Float64 = get(param, "LPSATR", LPSAT)
	WPSATR::Float64 = get(param, "WPSATR", WPSAT)
	PPSATR::Float64 = get(param, "PPSATR", PPSAT)
	LPSATB::Float64 = get(param, "LPSATB", 0.0)
	WPSATB::Float64 = get(param, "WPSATB", 0.0)
	PPSATB::Float64 = get(param, "PPSATB", 0.0)
	PSATX::Float64 = get(param, "PSATX", 1.0)
	PTWG::Float64 = get(param, "PTWG", 0.0)
	LPTWG::Float64 = get(param, "LPTWG", 0.0)
	WPTWG::Float64 = get(param, "WPTWG", 0.0)
	PPTWG::Float64 = get(param, "PPTWG", 0.0)
	PTWGL::Float64 = get(param, "PTWGL", 0.0)
	PTWGLEXP::Float64 = get(param, "PTWGLEXP", 1.0)
	PTWGR::Float64 = get(param, "PTWGR", PTWG)
	LPTWGR::Float64 = get(param, "LPTWGR", LPTWG)
	WPTWGR::Float64 = get(param, "WPTWGR", WPTWG)
	PPTWGR::Float64 = get(param, "PPTWGR", PPTWG)
	PTWGLR::Float64 = get(param, "PTWGLR", PTWGL)
	PTWGLEXPR::Float64 = get(param, "PTWGLEXPR", PTWGLEXP)

	# Velocity non-saturation effect
	A1::Float64 = get(param, "A1", 0.0)
	LA1::Float64 = get(param, "LA1", 0.0)
	WA1::Float64 = get(param, "WA1", 0.0)
	PA1::Float64 = get(param, "PA1", 0.0)
	A11::Float64 = get(param, "A11", 0.0)
	LA11::Float64 = get(param, "LA11", 0.0)
	WA11::Float64 = get(param, "WA11", 0.0)
	PA11::Float64 = get(param, "PA11", 0.0)
	A2::Float64 = get(param, "A2", 0.0)
	LA2::Float64 = get(param, "LA2", 0.0)
	WA2::Float64 = get(param, "WA2", 0.0)
	PA2::Float64 = get(param, "PA2", 0.0)
	A21::Float64 = get(param, "A21", 0.0)
	LA21::Float64 = get(param, "LA21", 0.0)
	WA21::Float64 = get(param, "WA21", 0.0)
	PA21::Float64 = get(param, "PA21", 0.0)

	# Output conductance
	PDIBLC::Float64 = get(param, "PDIBLC", 0.0)
	PDIBLCL::Float64 = get(param, "PDIBLCL", 0.0)
	PDIBLCLEXP::Float64 = get(param, "PDIBLCLEXP", 1.0)
	LPDIBLC::Float64 = get(param, "LPDIBLC", 0.0)
	WPDIBLC::Float64 = get(param, "WPDIBLC", 0.0)
	PPDIBLC::Float64 = get(param, "PPDIBLC", 0.0)
	PDIBLCR::Float64 = get(param, "PDIBLCR", PDIBLC)
	PDIBLCLR::Float64 = get(param, "PDIBLCLR", PDIBLCL)
	PDIBLCLEXPR::Float64 = get(param, "PDIBLCLEXPR", PDIBLCLEXP)
	LPDIBLCR::Float64 = get(param, "LPDIBLCR", LPDIBLC)
	WPDIBLCR::Float64 = get(param, "WPDIBLCR", WPDIBLC)
	PPDIBLCR::Float64 = get(param, "PPDIBLCR", PPDIBLC)
	PDIBLCB::Float64 = get(param, "PDIBLCB", 0.0)
	LPDIBLCB::Float64 = get(param, "LPDIBLCB", 0.0)
	WPDIBLCB::Float64 = get(param, "WPDIBLCB", 0.0)
	PPDIBLCB::Float64 = get(param, "PPDIBLCB", 0.0)
	PVAG::Float64 = get(param, "PVAG", 1.0)
	LPVAG::Float64 = get(param, "LPVAG", 0.0)
	WPVAG::Float64 = get(param, "WPVAG", 0.0)
	PPVAG::Float64 = get(param, "PPVAG", 0.0)
	FPROUT::Float64 = get(param, "FPROUT", 0.0)
	FPROUTL::Float64 = get(param, "FPROUTL", 0.0)
	FPROUTLEXP::Float64 = get(param, "FPROUTLEXP", 1.0)
	LFPROUT::Float64 = get(param, "LFPROUT", 0.0)
	WFPROUT::Float64 = get(param, "WFPROUT", 0.0)
	PFPROUT::Float64 = get(param, "PFPROUT", 0.0)

	# Impact ionization current
	ALPHA0::Float64 = get(param, "ALPHA0", 0.0)
	ALPHA0L::Float64 = get(param, "ALPHA0L", 0.0)
	ALPHA0LEXP::Float64 = get(param, "ALPHA0LEXP", 1.0)
	LALPHA0::Float64 = get(param, "LALPHA0", 0.0)
	WALPHA0::Float64 = get(param, "WALPHA0", 0.0)
	PALPHA0::Float64 = get(param, "PALPHA0", 0.0)
	BETA0::Float64 = get(param, "BETA0", 0.0)
	LBETA0::Float64 = get(param, "LBETA0", 0.0)
	WBETA0::Float64 = get(param, "WBETA0", 0.0)
	PBETA0::Float64 = get(param, "PBETA0", 0.0)

	# Gate dielectric tunnelling current model parameters
	AIGBACC::Float64 = get(param, "AIGBACC", 1.36e-2)
	BIGBACC::Float64 = get(param, "BIGBACC", 1.71e-3)
	CIGBACC::Float64 = get(param, "CIGBACC", 0.075)
	NIGBACC::Float64 = get(param, "NIGBACC", 1.0)
	AIGBINV::Float64 = get(param, "AIGBINV", 1.11e-2)
	BIGBINV::Float64 = get(param, "BIGBINV", 9.49e-4)
	CIGBINV::Float64 = get(param, "CIGBINV", 0.006)
	EIGBINV::Float64 = get(param, "EIGBINV", 1.1)
	NIGBINV::Float64 = get(param, "NIGBINV", 3.0)
	AIGC::Float64 = get(param, "AIGC", ((TYPE == ntype) ? 1.36e-2 : 9.8e-3))
	BIGC::Float64 = get(param, "BIGC", ((TYPE == ntype) ? 1.71e-3 : 7.59e-4))
	CIGC::Float64 = get(param, "CIGC", ((TYPE == ntype) ? 0.075 : 0.03))
	AIGS::Float64 = get(param, "AIGS", ((TYPE == ntype) ? 1.36e-2 : 9.8e-3))
	BIGS::Float64 = get(param, "BIGS", ((TYPE == ntype) ? 1.71e-3 : 7.59e-4))
	CIGS::Float64 = get(param, "CIGS", ((TYPE == ntype) ? 0.075 : 0.03))
	AIGD::Float64 = get(param, "AIGD", ((TYPE == ntype) ? 1.36e-2 : 9.8e-3))
	BIGD::Float64 = get(param, "BIGD", ((TYPE == ntype) ? 1.71e-3 : 7.59e-4))
	CIGD::Float64 = get(param, "CIGD", ((TYPE == ntype) ? 0.075 : 0.03))
	DLCIG::Float64 = get(param, "DLCIG", LINT)
	DLCIGD::Float64 = get(param, "DLCIGD", DLCIG)
	POXEDGE::Float64 = get(param, "POXEDGE", 1.0)
	NTOX::Float64 = get(param, "NTOX", 1.0)
	TOXREF::Float64 = get(param, "TOXREF", 3.0e-9)
	PIGCD::Float64 = get(param, "PIGCD", 1.0)
	AIGCL::Float64 = get(param, "AIGCL", 0.0)
	AIGCW::Float64 = get(param, "AIGCW", 0.0)
	AIGSL::Float64 = get(param, "AIGSL", 0.0)
	AIGSW::Float64 = get(param, "AIGSW", 0.0)
	AIGDL::Float64 = get(param, "AIGDL", 0.0)
	AIGDW::Float64 = get(param, "AIGDW", 0.0)
	PIGCDL::Float64 = get(param, "PIGCDL", 0.0)
	LAIGBINV::Float64 = get(param, "LAIGBINV", 0.0)
	WAIGBINV::Float64 = get(param, "WAIGBINV", 0.0)
	PAIGBINV::Float64 = get(param, "PAIGBINV", 0.0)
	LBIGBINV::Float64 = get(param, "LBIGBINV", 0.0)
	WBIGBINV::Float64 = get(param, "WBIGBINV", 0.0)
	PBIGBINV::Float64 = get(param, "PBIGBINV", 0.0)
	LCIGBINV::Float64 = get(param, "LCIGBINV", 0.0)
	WCIGBINV::Float64 = get(param, "WCIGBINV", 0.0)
	PCIGBINV::Float64 = get(param, "PCIGBINV", 0.0)
	LEIGBINV::Float64 = get(param, "LEIGBINV", 0.0)
	WEIGBINV::Float64 = get(param, "WEIGBINV", 0.0)
	PEIGBINV::Float64 = get(param, "PEIGBINV", 0.0)
	LNIGBINV::Float64 = get(param, "LNIGBINV", 0.0)
	WNIGBINV::Float64 = get(param, "WNIGBINV", 0.0)
	PNIGBINV::Float64 = get(param, "PNIGBINV", 0.0)
	LAIGBACC::Float64 = get(param, "LAIGBACC", 0.0)
	WAIGBACC::Float64 = get(param, "WAIGBACC", 0.0)
	PAIGBACC::Float64 = get(param, "PAIGBACC", 0.0)
	LBIGBACC::Float64 = get(param, "LBIGBACC", 0.0)
	WBIGBACC::Float64 = get(param, "WBIGBACC", 0.0)
	PBIGBACC::Float64 = get(param, "PBIGBACC", 0.0)
	LCIGBACC::Float64 = get(param, "LCIGBACC", 0.0)
	WCIGBACC::Float64 = get(param, "WCIGBACC", 0.0)
	PCIGBACC::Float64 = get(param, "PCIGBACC", 0.0)
	LNIGBACC::Float64 = get(param, "LNIGBACC", 0.0)
	WNIGBACC::Float64 = get(param, "WNIGBACC", 0.0)
	PNIGBACC::Float64 = get(param, "PNIGBACC", 0.0)
	LAIGC::Float64 = get(param, "LAIGC", 0.0)
	WAIGC::Float64 = get(param, "WAIGC", 0.0)
	PAIGC::Float64 = get(param, "PAIGC", 0.0)
	LBIGC::Float64 = get(param, "LBIGC", 0.0)
	WBIGC::Float64 = get(param, "WBIGC", 0.0)
	PBIGC::Float64 = get(param, "PBIGC", 0.0)
	LCIGC::Float64 = get(param, "LCIGC", 0.0)
	WCIGC::Float64 = get(param, "WCIGC", 0.0)
	PCIGC::Float64 = get(param, "PCIGC", 0.0)
	LAIGS::Float64 = get(param, "LAIGS", 0.0)
	WAIGS::Float64 = get(param, "WAIGS", 0.0)
	PAIGS::Float64 = get(param, "PAIGS", 0.0)
	LBIGS::Float64 = get(param, "LBIGS", 0.0)
	WBIGS::Float64 = get(param, "WBIGS", 0.0)
	PBIGS::Float64 = get(param, "PBIGS", 0.0)
	LCIGS::Float64 = get(param, "LCIGS", 0.0)
	WCIGS::Float64 = get(param, "WCIGS", 0.0)
	PCIGS::Float64 = get(param, "PCIGS", 0.0)
	LAIGD::Float64 = get(param, "LAIGD", 0.0)
	WAIGD::Float64 = get(param, "WAIGD", 0.0)
	PAIGD::Float64 = get(param, "PAIGD", 0.0)
	LBIGD::Float64 = get(param, "LBIGD", 0.0)
	WBIGD::Float64 = get(param, "WBIGD", 0.0)
	PBIGD::Float64 = get(param, "PBIGD", 0.0)
	LCIGD::Float64 = get(param, "LCIGD", 0.0)
	WCIGD::Float64 = get(param, "WCIGD", 0.0)
	PCIGD::Float64 = get(param, "PCIGD", 0.0)
	LPOXEDGE::Float64 = get(param, "LPOXEDGE", 0.0)
	WPOXEDGE::Float64 = get(param, "WPOXEDGE", 0.0)
	PPOXEDGE::Float64 = get(param, "PPOXEDGE", 0.0)
	LDLCIG::Float64 = get(param, "LDLCIG", 0.0)
	WDLCIG::Float64 = get(param, "WDLCIG", 0.0)
	PDLCIG::Float64 = get(param, "PDLCIG", 0.0)
	LDLCIGD::Float64 = get(param, "LDLCIGD", 0.0)
	WDLCIGD::Float64 = get(param, "WDLCIGD", 0.0)
	PDLCIGD::Float64 = get(param, "PDLCIGD", 0.0)
	LNTOX::Float64 = get(param, "LNTOX", 0.0)
	WNTOX::Float64 = get(param, "WNTOX", 0.0)
	PNTOX::Float64 = get(param, "PNTOX", 0.0)

	# GIDL and GISL currents
	AGIDL::Float64 = get(param, "AGIDL", 0.0)
	AGIDLL::Float64 = get(param, "AGIDLL", 0.0)
	AGIDLW::Float64 = get(param, "AGIDLW", 0.0)
	LAGIDL::Float64 = get(param, "LAGIDL", 0.0)
	WAGIDL::Float64 = get(param, "WAGIDL", 0.0)
	PAGIDL::Float64 = get(param, "PAGIDL", 0.0)
	BGIDL::Float64 = get(param, "BGIDL", 2.3e9)
	LBGIDL::Float64 = get(param, "LBGIDL", 0.0)
	WBGIDL::Float64 = get(param, "WBGIDL", 0.0)
	PBGIDL::Float64 = get(param, "PBGIDL", 0.0)
	CGIDL::Float64 = get(param, "CGIDL", 0.5)
	LCGIDL::Float64 = get(param, "LCGIDL", 0.0)
	WCGIDL::Float64 = get(param, "WCGIDL", 0.0)
	PCGIDL::Float64 = get(param, "PCGIDL", 0.0)
	EGIDL::Float64 = get(param, "EGIDL", 0.8)
	LEGIDL::Float64 = get(param, "LEGIDL", 0.0)
	WEGIDL::Float64 = get(param, "WEGIDL", 0.0)
	PEGIDL::Float64 = get(param, "PEGIDL", 0.0)
	AGISL::Float64 = get(param, "AGISL", AGIDL)
	AGISLL::Float64 = get(param, "AGISLL", AGIDLL)
	AGISLW::Float64 = get(param, "AGISLW", AGIDLW)
	LAGISL::Float64 = get(param, "LAGISL", LAGIDL)
	WAGISL::Float64 = get(param, "WAGISL", WAGIDL)
	PAGISL::Float64 = get(param, "PAGISL", PAGIDL)
	BGISL::Float64 = get(param, "BGISL", BGIDL)
	LBGISL::Float64 = get(param, "LBGISL", LBGIDL)
	WBGISL::Float64 = get(param, "WBGISL", WBGIDL)
	PBGISL::Float64 = get(param, "PBGISL", PBGIDL)
	CGISL::Float64 = get(param, "CGISL", CGIDL)
	LCGISL::Float64 = get(param, "LCGISL", LCGIDL)
	WCGISL::Float64 = get(param, "WCGISL", WCGIDL)
	PCGISL::Float64 = get(param, "PCGISL", PCGIDL)
	EGISL::Float64 = get(param, "EGISL", EGIDL)
	LEGISL::Float64 = get(param, "LEGISL", LEGIDL)
	WEGISL::Float64 = get(param, "WEGISL", WEGIDL)
	PEGISL::Float64 = get(param, "PEGISL", PEGIDL)

	# Overlap capacitance and fringing capacitance
	CF::Float64 = get(param, "CF", 0.0)
	LCF::Float64 = get(param, "LCF", 0.0)
	WCF::Float64 = get(param, "WCF", 0.0)
	PCF::Float64 = get(param, "PCF", 0.0)
	CFRCOEFF::Float64 = get(param, "CFRCOEFF", 1.0)
	CGSO::Float64 = get(param, "CGSO", 0.0)
	CGDO::Float64 = get(param, "CGDO", 0.0)
	CGBO::Float64 = get(param, "CGBO", 0.0)
	CGSL::Float64 = get(param, "CGSL", 0.0)
	LCGSL::Float64 = get(param, "LCGSL", 0.0)
	WCGSL::Float64 = get(param, "WCGSL", 0.0)
	PCGSL::Float64 = get(param, "PCGSL", 0.0)
	CGDL::Float64 = get(param, "CGDL", 0.0)
	LCGDL::Float64 = get(param, "LCGDL", 0.0)
	WCGDL::Float64 = get(param, "WCGDL", 0.0)
	PCGDL::Float64 = get(param, "PCGDL", 0.0)
	CKAPPAS::Float64 = get(param, "CKAPPAS", 0.6)
	LCKAPPAS::Float64 = get(param, "LCKAPPAS", 0.0)
	WCKAPPAS::Float64 = get(param, "WCKAPPAS", 0.0)
	PCKAPPAS::Float64 = get(param, "PCKAPPAS", 0.0)
	CKAPPAD::Float64 = get(param, "CKAPPAD", 0.6)
	LCKAPPAD::Float64 = get(param, "LCKAPPAD", 0.0)
	WCKAPPAD::Float64 = get(param, "WCKAPPAD", 0.0)
	PCKAPPAD::Float64 = get(param, "PCKAPPAD", 0.0)

	# Layout-dependent parasitics model parameters (resistance only)
	DMCG::Float64 = get(param, "DMCG", 0.0)
	DMCI::Float64 = get(param, "DMCI", DMCG)
	DMDG::Float64 = get(param, "DMDG", 0.0)
	DMCGT::Float64 = get(param, "DMCGT", 0.0)
	XGL::Float64 = get(param, "XGL", 0.0)
	RSHG::Float64 = get(param, "RSHG", 0.1)

	# Junction capacitance
	CJS::Float64 = get(param, "CJS", 5.0e-4)
	CJD::Float64 = get(param, "CJD", CJS)
	CJSWS::Float64 = get(param, "CJSWS", 5.0e-10)
	CJSWD::Float64 = get(param, "CJSWD", CJSWS)
	CJSWGS::Float64 = get(param, "CJSWGS", 0.0)
	CJSWGD::Float64 = get(param, "CJSWGD", CJSWGS)
	PBS::Float64 = get(param, "PBS", 1.0)
	PBD::Float64 = get(param, "PBD", PBS)
	PBSWS::Float64 = get(param, "PBSWS", 1.0)
	PBSWD::Float64 = get(param, "PBSWD", PBSWS)
	PBSWGS::Float64 = get(param, "PBSWGS", PBSWS)
	PBSWGD::Float64 = get(param, "PBSWGD", PBSWGS)
	MJS::Float64 = get(param, "MJS", 0.5)
	MJD::Float64 = get(param, "MJD", MJS)
	MJSWS::Float64 = get(param, "MJSWS", 0.33)
	MJSWD::Float64 = get(param, "MJSWD", MJSWS)
	MJSWGS::Float64 = get(param, "MJSWGS", MJSWS)
	MJSWGD::Float64 = get(param, "MJSWGD", MJSWGS)

	# Junction current
	JSS::Float64 = get(param, "JSS", 1.0e-4)
	JSD::Float64 = get(param, "JSD", JSS)
	JSWS::Float64 = get(param, "JSWS", 0.0)
	JSWD::Float64 = get(param, "JSWD", JSWS)
	JSWGS::Float64 = get(param, "JSWGS", 0.0)
	JSWGD::Float64 = get(param, "JSWGD", JSWGS)
	NJS::Float64 = get(param, "NJS", 1.0)
	NJD::Float64 = get(param, "NJD", NJS)
	IJTHSFWD::Float64 = get(param, "IJTHSFWD", 0.1)
	IJTHDFWD::Float64 = get(param, "IJTHDFWD", IJTHSFWD)
	IJTHSREV::Float64 = get(param, "IJTHSREV", 0.1)
	IJTHDREV::Float64 = get(param, "IJTHDREV", IJTHSREV)
	BVS::Float64 = get(param, "BVS", 10.0)
	BVD::Float64 = get(param, "BVD", BVS)
	XJBVS::Float64 = get(param, "XJBVS", 1.0)
	XJBVD::Float64 = get(param, "XJBVD", XJBVS)

	# Tunneling component of junction current
	JTSS::Float64 = get(param, "JTSS", 0.0)
	JTSD::Float64 = get(param, "JTSD", JTSS)
	JTSSWS::Float64 = get(param, "JTSSWS", 0.0)
	JTSSWD::Float64 = get(param, "JTSSWD", JTSSWS)
	JTSSWGS::Float64 = get(param, "JTSSWGS", 0.0)
	JTSSWGD::Float64 = get(param, "JTSSWGD", JTSSWGS)
	JTWEFF::Float64 = get(param, "JTWEFF", 0.0)
	NJTS::Float64 = get(param, "NJTS", 20.0)
	NJTSD::Float64 = get(param, "NJTSD", NJTS)
	NJTSSW::Float64 = get(param, "NJTSSW", 20.0)
	NJTSSWD::Float64 = get(param, "NJTSSWD", NJTSSW)
	NJTSSWG::Float64 = get(param, "NJTSSWG", 20.0)
	NJTSSWGD::Float64 = get(param, "NJTSSWGD", NJTSSWG)
	VTSS::Float64 = get(param, "VTSS", 10.0)
	VTSD::Float64 = get(param, "VTSD", VTSS)
	VTSSWS::Float64 = get(param, "VTSSWS", 10.0)
	VTSSWD::Float64 = get(param, "VTSSWD", VTSSWS)
	VTSSWGS::Float64 = get(param, "VTSSWGS", 10.0)
	VTSSWGD::Float64 = get(param, "VTSSWGD", VTSSWGS)

	# High-speed/RF model parameters
	XRCRG1::Float64 = get(param, "XRCRG1", 12.0)
	XRCRG2::Float64 = get(param, "XRCRG2", 1.0)
	GBMIN::Float64 = get(param, "GBMIN", 1.0e-12)
	RBPS0::Float64 = get(param, "RBPS0", 50.0)
	RBPSL::Float64 = get(param, "RBPSL", 0.0)
	RBPSW::Float64 = get(param, "RBPSW", 0.0)
	RBPSNF::Float64 = get(param, "RBPSNF", 0.0)
	RBPD0::Float64 = get(param, "RBPD0", 50.0)
	RBPDL::Float64 = get(param, "RBPDL", 0.0)
	RBPDW::Float64 = get(param, "RBPDW", 0.0)
	RBPDNF::Float64 = get(param, "RBPDNF", 0.0)
	RBPBX0::Float64 = get(param, "RBPBX0", 100.0)
	RBPBXL::Float64 = get(param, "RBPBXL", 0.0)
	RBPBXW::Float64 = get(param, "RBPBXW", 0.0)
	RBPBXNF::Float64 = get(param, "RBPBXNF", 0.0)
	RBPBY0::Float64 = get(param, "RBPBY0", 100.0)
	RBPBYL::Float64 = get(param, "RBPBYL", 0.0)
	RBPBYW::Float64 = get(param, "RBPBYW", 0.0)
	RBPBYNF::Float64 = get(param, "RBPBYNF", 0.0)
	RBSBX0::Float64 = get(param, "RBSBX0", 100.0)
	RBSBY0::Float64 = get(param, "RBSBY0", 100.0)
	RBDBX0::Float64 = get(param, "RBDBX0", 100.0)
	RBDBY0::Float64 = get(param, "RBDBY0", 100.0)
	RBSDBXL::Float64 = get(param, "RBSDBXL", 0.0)
	RBSDBXW::Float64 = get(param, "RBSDBXW", 0.0)
	RBSDBXNF::Float64 = get(param, "RBSDBXNF", 0.0)
	RBSDBYL::Float64 = get(param, "RBSDBYL", 0.0)
	RBSDBYW::Float64 = get(param, "RBSDBYW", 0.0)
	RBSDBYNF::Float64 = get(param, "RBSDBYNF", 0.0)

	# Flicker noise
	EF::Float64 = get(param, "EF", 1.0)
	EM::Float64 = get(param, "EM", 4.1e7)
	NOIA::Float64 = get(param, "NOIA", 6.250e40)
	NOIB::Float64 = get(param, "NOIB", 3.125e25)
	NOIC::Float64 = get(param, "NOIC", 8.750e8)
	LINTNOI::Float64 = get(param, "LINTNOI", 0.0)
	NOIA1::Float64 = get(param, "NOIA1", 0.0)
	NOIAX::Float64 = get(param, "NOIAX", 1.0)

	# Thermal noise
	NTNOI::Float64 = get(param, "NTNOI", 1.0)
	RNOIA::Float64 = get(param, "RNOIA", 0.577)
	RNOIB::Float64 = get(param, "RNOIB", 0.5164)
	RNOIC::Float64 = get(param, "RNOIC", 0.395)
	TNOIA::Float64 = get(param, "TNOIA", 1.5)
	TNOIB::Float64 = get(param, "TNOIB", 3.5)
	TNOIC::Float64 = get(param, "TNOIC", 0.0)

	# Binning parameters
	BINUNIT::UInt8 = get(param, "BINUNIT", 1)
	DLBIN::Float64 = get(param, "DLBIN", 0.0)
	DWBIN::Float64 = get(param, "DWBIN", 0.0)

	# Temperature dependence parameters
	TNOM::Float64 = get(param, "TNOM", 27.0)
	TBGASUB::Float64 = get(param, "TBGASUB", 4.73e-4)
	TBGBSUB::Float64 = get(param, "TBGBSUB", 636.0)
	TNFACTOR::Float64 = get(param, "TNFACTOR", 0.0)
	UTE::Float64 = get(param, "UTE", -1.5)
	LUTE::Float64 = get(param, "LUTE", 0.0)
	WUTE::Float64 = get(param, "WUTE", 0.0)
	PUTE::Float64 = get(param, "PUTE", 0.0)
	UTEL::Float64 = get(param, "UTEL", 0.0)
	UA1::Float64 = get(param, "UA1", 1.0e-3)
	LUA1::Float64 = get(param, "LUA1", 0.0)
	WUA1::Float64 = get(param, "WUA1", 0.0)
	PUA1::Float64 = get(param, "PUA1", 0.0)
	UA1L::Float64 = get(param, "UA1L", 0.0)
	UC1::Float64 = get(param, "UC1", 0.056e-9)
	LUC1::Float64 = get(param, "LUC1", 0.0)
	WUC1::Float64 = get(param, "WUC1", 0.0)
	PUC1::Float64 = get(param, "PUC1", 0.0)
	UD1::Float64 = get(param, "UD1", 0.0)
	LUD1::Float64 = get(param, "LUD1", 0.0)
	WUD1::Float64 = get(param, "WUD1", 0.0)
	PUD1::Float64 = get(param, "PUD1", 0.0)
	UD1L::Float64 = get(param, "UD1L", 0.0)
	EU1::Float64 = get(param, "EU1", 0.0)
	LEU1::Float64 = get(param, "LEU1", 0.0)
	WEU1::Float64 = get(param, "WEU1", 0.0)
	PEU1::Float64 = get(param, "PEU1", 0.0)
	UCSTE::Float64 = get(param, "UCSTE", -4.775e-3)
	LUCSTE::Float64 = get(param, "LUCSTE", 0.0)
	WUCSTE::Float64 = get(param, "WUCSTE", 0.0)
	PUCSTE::Float64 = get(param, "PUCSTE", 0.0)
	TETA0::Float64 = get(param, "TETA0", 0.0)
	PRT::Float64 = get(param, "PRT", 0.0)
	LPRT::Float64 = get(param, "LPRT", 0.0)
	WPRT::Float64 = get(param, "WPRT", 0.0)
	PPRT::Float64 = get(param, "PPRT", 0.0)
	AT::Float64 = get(param, "AT", -1.56e-3)
	LAT::Float64 = get(param, "LAT", 0.0)
	WAT::Float64 = get(param, "WAT", 0.0)
	PAT::Float64 = get(param, "PAT", 0.0)
	ATL::Float64 = get(param, "ATL", 0.0)
	TDELTA::Float64 = get(param, "TDELTA", 0.0)
	PTWGT::Float64 = get(param, "PTWGT", 0.0)
	LPTWGT::Float64 = get(param, "LPTWGT", 0.0)
	WPTWGT::Float64 = get(param, "WPTWGT", 0.0)
	PPTWGT::Float64 = get(param, "PPTWGT", 0.0)
	PTWGTL::Float64 = get(param, "PTWGTL", 0.0)
	KT1::Float64 = get(param, "KT1", -0.11)
	KT1EXP::Float64 = get(param, "KT1EXP", 1.0)
	KT1L::Float64 = get(param, "KT1L", 0.0)
	LKT1::Float64 = get(param, "LKT1", 0.0)
	WKT1::Float64 = get(param, "WKT1", 0.0)
	PKT1::Float64 = get(param, "PKT1", 0.0)
	KT2::Float64 = get(param, "KT2", 0.022)
	LKT2::Float64 = get(param, "LKT2", 0.0)
	WKT2::Float64 = get(param, "WKT2", 0.0)
	PKT2::Float64 = get(param, "PKT2", 0.0)
	IIT::Float64 = get(param, "IIT", 0.0)
	LIIT::Float64 = get(param, "LIIT", 0.0)
	WIIT::Float64 = get(param, "WIIT", 0.0)
	PIIT::Float64 = get(param, "PIIT", 0.0)
	IGT::Float64 = get(param, "IGT", 2.5)
	LIGT::Float64 = get(param, "LIGT", 0.0)
	WIGT::Float64 = get(param, "WIGT", 0.0)
	PIGT::Float64 = get(param, "PIGT", 0.0)
	TGIDL::Float64 = get(param, "TGIDL", 0.0)
	LTGIDL::Float64 = get(param, "LTGIDL", 0.0)
	WTGIDL::Float64 = get(param, "WTGIDL", 0.0)
	PTGIDL::Float64 = get(param, "PTGIDL", 0.0)
	TCJ::Float64 = get(param, "TCJ", 0.0)
	TCJSW::Float64 = get(param, "TCJSW", 0.0)
	TCJSWG::Float64 = get(param, "TCJSWG", 0.0)
	TPB::Float64 = get(param, "TPB", 0.0)
	TPBSW::Float64 = get(param, "TPBSW", 0.0)
	TPBSWG::Float64 = get(param, "TPBSWG", 0.0)
	XTIS::Float64 = get(param, "XTIS", 3.0)
	XTID::Float64 = get(param, "XTID", XTIS)
	XTSS::Float64 = get(param, "XTSS", 0.02)
	XTSD::Float64 = get(param, "XTSD", XTSS)
	XTSSWS::Float64 = get(param, "XTSSWS", 0.02)
	XTSSWD::Float64 = get(param, "XTSSWD", XTSSWS)
	XTSSWGS::Float64 = get(param, "XTSSWGS", 0.02)
	XTSSWGD::Float64 = get(param, "XTSSWGD", XTSSWGS)
	TNJTS::Float64 = get(param, "TNJTS", 0.0)
	TNJTSD::Float64 = get(param, "TNJTSD", TNJTS)
	TNJTSSW::Float64 = get(param, "TNJTSSW", 0.0)
	TNJTSSWD::Float64 = get(param, "TNJTSSWD", TNJTSSW)
	TNJTSSWG::Float64 = get(param, "TNJTSSWG", 0.0)
	TNJTSSWGD::Float64 = get(param, "TNJTSSWGD", TNJTSSWG)

	# Self heating parameters
	RTH0::Float64 = get(param, "RTH0", 0.0)
	CTH0::Float64 = get(param, "CTH0", 1.0e-5)
	WTH0::Float64 = get(param, "WTH0", 0.0)

	# Stress related parameters
	SAREF::Float64 = get(param, "SAREF", 1.0e-6)
	SBREF::Float64 = get(param, "SBREF", 1.0e-6)
	WLOD::Float64 = get(param, "WLOD", 0.0)
	KU0::Float64 = get(param, "KU0", 0.0)
	KVSAT::Float64 = get(param, "KVSAT", 0.0)
	TKU0::Float64 = get(param, "TKU0", 0.0)
	LKU0::Float64 = get(param, "LKU0", 0.0)
	WKU0::Float64 = get(param, "WKU0", 0.0)
	PKU0::Float64 = get(param, "PKU0", 0.0)
	LLODKU0::Float64 = get(param, "LLODKU0", 0.0)
	WLODKU0::Float64 = get(param, "WLODKU0", 0.0)
	KVTH0::Float64 = get(param, "KVTH0", 0.0)
	LKVTH0::Float64 = get(param, "LKVTH0", 0.0)
	WKVTH0::Float64 = get(param, "WKVTH0", 0.0)
	PKVTH0::Float64 = get(param, "PKVTH0", 0.0)
	LLODVTH::Float64 = get(param, "LLODVTH", 0.0)
	WLODVTH::Float64 = get(param, "WLODVTH", 0.0)
	STK2::Float64 = get(param, "STK2", 0.0)
	LODK2::Float64 = get(param, "LODK2", 0.0)
	STETA0::Float64 = get(param, "STETA0", 0.0)
	LODETA0::Float64 = get(param, "LODETA0", 0.0)

	# Well proximity parameters
	WEB::Float64 = get(param, "WEB", 0.0)
	WEC::Float64 = get(param, "WEC", 0.0)
	KVTH0WE::Float64 = get(param, "KVTH0WE", 0.0)
	LKVTH0WE::Float64 = get(param, "LKVTH0WE", 0.0)
	WKVTH0WE::Float64 = get(param, "WKVTH0WE", 0.0)
	PKVTH0WE::Float64 = get(param, "PKVTH0WE", 0.0)
	K2WE::Float64 = get(param, "K2WE", 0.0)
	LK2WE::Float64 = get(param, "LK2WE", 0.0)
	WK2WE::Float64 = get(param, "WK2WE", 0.0)
	PK2WE::Float64 = get(param, "PK2WE", 0.0)
	KU0WE::Float64 = get(param, "KU0WE", 0.0)
	LKU0WE::Float64 = get(param, "LKU0WE", 0.0)
	WKU0WE::Float64 = get(param, "WKU0WE", 0.0)
	PKU0WE::Float64 = get(param, "PKU0WE", 0.0)
	SCREF::Float64 = get(param, "SCREF", 1.0e-6)

	# Sub-surface leakage drain current
	SSL0::Float64 = get(param, "SSL0", 4.0e2)
	SSL1::Float64 = get(param, "SSL1", 3.36e8)
	SSL2::Float64 = get(param, "SSL2", 0.185)
	SSL3::Float64 = get(param, "SSL3", 0.3)
	SSL4::Float64 = get(param, "SSL4", 1.4)
	SSL5::Float64 = get(param, "SSL5", 0)
	SSLEXP1::Float64 = get(param, "SSLEXP1", 0.490)
	SSLEXP2::Float64 = get(param, "SSLEXP2", 1.42)

	# Vdsx smoothing
	AVDSX::Float64 = get(param, "AVDSX", 20.0)

	# STI edge FET device parameters
	WEDGE::Float64 = get(param, "WEDGE", 10.0e-9)
	DGAMMAEDGE::Float64 = get(param, "DGAMMAEDGE", 0.0)
	DGAMMAEDGEL::Float64 = get(param, "DGAMMAEDGEL", 0.0)
	DGAMMAEDGELEXP::Float64 = get(param, "DGAMMAEDGELEXP", 1.0)
	DVTEDGE::Float64 = get(param, "DVTEDGE", 0.0)
	NDEPEDGE::Float64 = get(param, "NDEPEDGE", 1e24)
	LNDEPEDGE::Float64 = get(param, "LNDEPEDGE", 0.0)
	WNDEPEDGE::Float64 = get(param, "WNDEPEDGE", 0.0)
	PNDEPEDGE::Float64 = get(param, "PNDEPEDGE", 0.0)
	NFACTOREDGE::Float64 = get(param, "NFACTOREDGE", 0.0)
	LNFACTOREDGE::Float64 = get(param, "LNFACTOREDGE", 0.0)
	WNFACTOREDGE::Float64 = get(param, "WNFACTOREDGE", 0.0)
	PNFACTOREDGE::Float64 = get(param, "PNFACTOREDGE", 0.0)
	CITEDGE::Float64 = get(param, "CITEDGE", 0.0)
	LCITEDGE::Float64 = get(param, "LCITEDGE", 0.0)
	WCITEDGE::Float64 = get(param, "WCITEDGE", 0.0)
	PCITEDGE::Float64 = get(param, "PCITEDGE", 0.0)
	CDSCDEDGE::Float64 = get(param, "CDSCDEDGE", 1e-9)
	LCDSCDEDGE::Float64 = get(param, "LCDSCDEDGE", 0.0)
	WCDSCDEDGE::Float64 = get(param, "WCDSCDEDGE", 0.0)
	PCDSCDEDGE::Float64 = get(param, "PCDSCDEDGE", 0.0)
	CDSCBEDGE::Float64 = get(param, "CDSCBEDGE", 0.0)
	LCDSCBEDGE::Float64 = get(param, "LCDSCBEDGE", 0.0)
	WCDSCBEDGE::Float64 = get(param, "WCDSCBEDGE", 0.0)
	PCDSCBEDGE::Float64 = get(param, "PCDSCBEDGE", 0.0)
	ETA0EDGE::Float64 = get(param, "ETA0EDGE", 0.08)
	LETA0EDGE::Float64 = get(param, "LETA0EDGE", 0.0)
	WETA0EDGE::Float64 = get(param, "WETA0EDGE", 0.0)
	PETA0EDGE::Float64 = get(param, "PETA0EDGE", 0.0)
	ETABEDGE::Float64 = get(param, "ETABEDGE", -0.07)
	LETABEDGE::Float64 = get(param, "LETABEDGE", 0.0)
	WETABEDGE::Float64 = get(param, "WETABEDGE", 0.0)
	PETABEDGE::Float64 = get(param, "PETABEDGE", 0.0)
	KT1EDGE::Float64 = get(param, "KT1EDGE", -0.11)
	LKT1EDGE::Float64 = get(param, "LKT1EDGE", 0.0)
	WKT1EDGE::Float64 = get(param, "WKT1EDGE", 0.0)
	PKT1EDGE::Float64 = get(param, "PKT1EDGE", 0.0)
	KT1LEDGE::Float64 = get(param, "KT1LEDGE", 0.0)
	LKT1LEDGE::Float64 = get(param, "LKT1LEDGE", 0.0)
	WKT1LEDGE::Float64 = get(param, "WKT1LEDGE", 0.0)
	PKT1LEDGE::Float64 = get(param, "PKT1LEDGE", 0.0)
	KT2EDGE::Float64 = get(param, "KT2EDGE", 0.022)
	LKT2EDGE::Float64 = get(param, "LKT2EDGE", 0.0)
	WKT2EDGE::Float64 = get(param, "WKT2EDGE", 0.0)
	PKT2EDGE::Float64 = get(param, "PKT2EDGE", 0.0)
	KT1EXPEDGE::Float64 = get(param, "KT1EXPEDGE", 1.0)
	LKT1EXPEDGE::Float64 = get(param, "LKT1EXPEDGE", 0.0)
	WKT1EXPEDGE::Float64 = get(param, "WKT1EXPEDGE", 0.0)
	PKT1EXPEDGE::Float64 = get(param, "PKT1EXPEDGE", 0.0)
	TNFACTOREDGE::Float64 = get(param, "TNFACTOREDGE", 0.0)
	LTNFACTOREDGE::Float64 = get(param, "LTNFACTOREDGE", 0.0)
	WTNFACTOREDGE::Float64 = get(param, "WTNFACTOREDGE", 0.0)
	PTNFACTOREDGE::Float64 = get(param, "PTNFACTOREDGE", 0.0)
	TETA0EDGE::Float64 = get(param, "TETA0EDGE", 0.0)
	LTETA0EDGE::Float64 = get(param, "LTETA0EDGE", 0.0)
	WTETA0EDGE::Float64 = get(param, "WTETA0EDGE", 0.0)
	PTETA0EDGE::Float64 = get(param, "PTETA0EDGE", 0.0)
	DVT0EDGE::Float64 = get(param, "DVT0EDGE", 2.2)
	DVT1EDGE::Float64 = get(param, "DVT1EDGE", 0.53)
	DVT2EDGE::Float64 = get(param, "DVT2EDGE", 0.0)
	K2EDGE::Float64 = get(param, "K2EDGE", 0.0)
	LK2EDGE::Float64 = get(param, "LK2EDGE", 0.0)
	WK2EDGE::Float64 = get(param, "WK2EDGE", 0.0)
	PK2EDGE::Float64 = get(param, "PK2EDGE", 0.0)
	KVTH0EDGE::Float64 = get(param, "KVTH0EDGE", 0.0)
	LKVTH0EDGE::Float64 = get(param, "LKVTH0EDGE", 0.0)
	WKVTH0EDGE::Float64 = get(param, "WKVTH0EDGE", 0.0)
	PKVTH0EDGE::Float64 = get(param, "PKVTH0EDGE", 0.0)
	STK2EDGE::Float64 = get(param, "STK2EDGE", 0.0)
	LSTK2EDGE::Float64 = get(param, "LSTK2EDGE", 0.0)
	WSTK2EDGE::Float64 = get(param, "WSTK2EDGE", 0.0)
	PSTK2EDGE::Float64 = get(param, "PSTK2EDGE", 0.0)
	STETA0EDGE::Float64 = get(param, "STETA0EDGE", 0.0)
	LSTETA0EDGE::Float64 = get(param, "LSTETA0EDGE", 0.0)
	WSTETA0EDGE::Float64 = get(param, "WSTETA0EDGE", 0.0)
	PSTETA0EDGE::Float64 = get(param, "PSTETA0EDGE", 0.0)
	IGCLAMP::UInt8 = get(param, "IGCLAMP", 1)
	LP::Float64 = get(param, "LP", 1.0e-5)
	RNOIK::Float64 = get(param, "RNOIK", 0.0)
	TNOIK::Float64 = get(param, "TNOIK", 0.0)
	TNOIK2::Float64 = get(param, "TNOIK2", 0.1)
	K0::Float64 = get(param, "K0", 0.0)
	LK0::Float64 = get(param, "LK0", 0.0)
	WK0::Float64 = get(param, "WK0", 0.0)
	PK0::Float64 = get(param, "PK0", 0.0)
	K01::Float64 = get(param, "K01", 0.0)
	LK01::Float64 = get(param, "LK01", 0.0)
	WK01::Float64 = get(param, "WK01", 0.0)
	PK01::Float64 = get(param, "PK01", 0.0)
	M0::Float64 = get(param, "M0", 1.0)
	LM0::Float64 = get(param, "LM0", 0.0)
	WM0::Float64 = get(param, "WM0", 0.0)
	PM0::Float64 = get(param, "PM0", 0.0)
	M01::Float64 = get(param, "M01", 0.0)
	LM01::Float64 = get(param, "LM01", 0.0)
	WM01::Float64 = get(param, "WM01", 0.0)
	PM01::Float64 = get(param, "PM01", 0.0)

	# Flicker noise model parameter for EDGE FET transistor
	NEDGE::Float64 = get(param, "NEDGE", 1)
	NOIA1_EDGE::Float64 = get(param, "NOIA1_EDGE", 0.0)
	NOIAX_EDGE::Float64 = get(param, "NOIAX_EDGE", 1.0)


	# Flicker noise model parameter for Halo transistor
	FNOIMOD::UInt8 = get(param, "FNOIMOD", 0)
	LH::Float64 = get(param, "LH", 1.0e-8)
	NOIA2::Float64 = get(param, "NOIA2", NOIA)
	HNDEP::Float64 = get(param, "HNDEP", NDEP)

	# Flexibility of tuning Cgg in strong inversion
	ABULK::Float64 = get(param, "ABULK", 1.0)

	# To enhance the fitting flexiblity for the gm/Id
	C0::Float64 = get(param, "C0", 0.0)
	LC0::Float64 = get(param, "LC0", 0.0)
	WC0::Float64 = get(param, "WC0", 0.0)
	PC0::Float64 = get(param, "PC0", 0.0)
	C01::Float64 = get(param, "C01", 0.0)
	LC01::Float64 = get(param, "LC01", 0.0)
	WC01::Float64 = get(param, "WC01", 0.0)
	PC01::Float64 = get(param, "PC01", 0.0)
	C0SI::Float64 = get(param, "C0SI", 1.0)
	LC0SI::Float64 = get(param, "LC0SI", 0.0)
	WC0SI::Float64 = get(param, "WC0SI", 0.0)
	PC0SI::Float64 = get(param, "PC0SI", 0.0)
	C0SI1::Float64 = get(param, "C0SI1", 0.0)
	LC0SI1::Float64 = get(param, "LC0SI1", 0.0)
	WC0SI1::Float64 = get(param, "WC0SI1", 0.0)
	PC0SI1::Float64 = get(param, "PC0SI1", 0.0)
	C0SISAT::Float64 = get(param, "C0SISAT", 0.0)
	LC0SISAT::Float64 = get(param, "LC0SISAT", 0.0)
	WC0SISAT::Float64 = get(param, "WC0SISAT", 0.0)
	PC0SISAT::Float64 = get(param, "PC0SISAT", 0.0)
	C0SISAT1::Float64 = get(param, "C0SISAT1", 0.0)
	LC0SISAT1::Float64 = get(param, "LC0SISAT1", 0.0)
	WC0SISAT1::Float64 = get(param, "WC0SISAT1", 0.0)
	PC0SISAT1::Float64 = get(param, "PC0SISAT1", 0.0)

	# Minimum resistance value
	minr::Float64 = get(param, "minr", 1.0e-3)
	# High Voltage Model Parameters

	# --- Mod selectors -----
	HVMOD::UInt8 = get(param, "HVMOD", 0)
	HVCAP::UInt8 = get(param, "HVCAP", 0)
	HVCAPS::UInt8 = get(param, "HVCAPS", 0)
	IIMOD::UInt8 = get(param, "IIMOD", 0)

	# --- Other parameters -----
	NDRIFTD::Float64 = get(param, "NDRIFTD", 5.0e16)
	VDRIFT::Float64 = get(param, "VDRIFT", 1.0e5)
	MDRIFT::Float64 = get(param, "MDRIFT", 1.0)
	NDRIFTS::Float64 = get(param, "NDRIFTS", NDRIFTD)
	RDLCW::Float64 = get(param, "RDLCW", 100.0)
	RSLCW::Float64 = get(param, "RSLCW", 0)
	PDRWB::Float64 = get(param, "PDRWB", 0)
	VFBDRIFT::Float64 = get(param, "VFBDRIFT", -1)
	VFBOV::Float64 = get(param, "VFBOV", -1)
	LOVER::Float64 = get(param, "LOVER", 500e-9)
	LOVERACC::Float64 = get(param, "LOVERACC", LOVER)
	NDR::Float64 = get(param, "NDR", NDEP)
	SLHV::Float64 = get(param, "SLHV", 0)
	SLHV1::Float64 = get(param, "SLHV1", 1.0)
	ALPHADR::Float64 = get(param, "ALPHADR", ALPHA0)
	BETADR::Float64 = get(param, "BETADR", BETA0)
	PRTHV::Float64 = get(param, "PRTHV", 0.0)
	ATHV::Float64 = get(param, "ATHV", 0)
	HVFACTOR::Float64 = get(param, "HVFACTOR", 1e-3)
	DRII1::Float64 = get(param, "DRII1", 1.0)
	DRII2::Float64 = get(param, "DRII2", 5)
	DELTAII::Float64 = get(param, "DELTAII", 0.5)

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
	if !("TOXP" in keys(param))
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

	# Process drain series resistance
	DMCGeff = DMCG - DMCGT
	DMCIeff = DMCI
	DMDGeff = DMDG - DMCGT

	# Processing S/D resistances and conductances
	if "NRS" in keys(param)
		RSourceGeo = RSH * NRS
	elseif (RGEOMOD > 0 && RSH > 0.0)
		RSourceGeo = BSIMBULKRdseffGeo(NF, GEOMOD, RGEOMOD, MINZ, Weff, RSH, DMCGeff, DMCIeff, DMDGeff, 1)
	else
		RSourceGeo = 0.0
	end

	if "NRD" in keys(param)
		RDrainGeo = RSH * NRD
	elseif (RGEOMOD > 0 && RSH > 0.0)
		RDrainGeo = BSIMBULKRdseffGeo(NF, GEOMOD, RGEOMOD, MINZ, Weff, RSH, DMCGeff, DMCIeff, DMDGeff, 0)
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
		if !("RBPS0" in keys(param)) || !("RBPD0" in keys(param))
			Bodymode = 1
		elseif !("RBSBX0" in keys(param)) && !("RBSBY0" in keys(param)) || !("RBDBX0" in keys(param)) && !("RBDBY0" in keys(param))
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
			println("Warning: (instance %M) The gate conductance reset to 1.0e3 mho.")
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
	DevTemp = Temp + P_CELSIUS0 + DTEMP

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
	temp_PSeff, temp_PDeff, temp_ASeff, temp_ADeff = BSIMBULKPAeffGeo(NF, GEOMOD, MINZ, Weffcj, DMCGeff, DMCIeff, DMDGeff)
	if "AS" in keys(param)
		ASeff = AS * WMLT * LMLT
	else
		ASeff = temp_ASeff
	end
	if (ASeff < 0.0)
		println("Warning: (instance %M) ASeff = %e is negative. Set to 0.0.", ASeff)
		ASeff = 0.0
	end
	if "AD" in keys(param)
		ADeff = AD * WMLT * LMLT
	else
		ADeff = temp_ADeff
	end
	if (ADeff < 0.0)
		println("Warning: (instance %M) ADeff = %e is negative. Set to 0.0.", ADeff)
		ADeff = 0.0
	end
	if "PS" in keys(param)
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
	if "PD" in keys(param)
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
		for i = 0 : NF-1
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
		if !("SCA" in keys(param)) && !("SCB" in keys(param)) && !("SCC" in keys(param))
			if ("SC" in keys(param)) && (SC > 0.0)
				T1        = SC + Wdrn
				T2        = 1.0 / SCREF
				local_sca = SCREF * SCREF / (SC * T1)
				local_scb = ((0.1 * SC + 0.01 * SCREF) * lexp(-10.0 * SC * T2)  - (0.1 * T1 + 0.01 * SCREF) * lexp(-10.0 * T1 * T2)) / Wdrn
				local_scc = ((0.05 * SC + 0.0025 * SCREF) * lexp(-20.0 * SC * T2)  - (0.05 * T1 + 0.0025 * SCREF) * lexp(-20.0 * T1 * T2)) / Wdrn
			else
				println("Warning: (Instance %M) No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive.")
			end
		end
	else
		local_sca = 0
		local_scb = 0
		local_scc = 0
	end

	vth0_well = KVTH0WE_i * (local_sca + WEB * local_scb + WEC * local_scc)
	k2_well   = K2WE_i * (local_sca + WEB * local_scb + WEC * local_scc)
	mu_well   = 1.0 + KU0WE_i * (local_sca + WEB * local_scb + WEC * local_scc)
	U0_t      = U0_t * mu_well
	K2_i      = K2_i + k2_well

	# Load terminal voltages
	Vds           = Vd - Vs
	Vds_noswap    = Vds
	Vsb_noswap    = Vs
	Vdb_noswap    = Vd
	Vbs_jct       = devsign * (Vb - Vs)
	Vbd_jct       = devsign * (Vb - Vd)
	Vgd_noswap    = Vg - Vd
	Vgs_noswap    = Vg - Vs
	Vgd_ov_noswap = devsign * (Vg - Vd)
	Vgs_ov_noswap = devsign * (Vg - Vs)

	# Terminal voltage conditioning
	# Source-drain interchange
	sigvds = 1.0
	if (Vds < 0.0)
		sigvds = -1.0
		Vd = devsign * Vs
		Vs = devsign * Vd
	end
	Vds  = Vd - Vs
	T0   = AVDSX * Vds
	if (T0 > EXPL_THRESHOLD)
	   T1 = T0
	else
	   T1 = log(1.0 + exp(T0))
	end
	Vdsx = ((2.0 / AVDSX) * T1) - Vds - ((2.0 / AVDSX) * log(2.0))
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
	PhistVbs = Smooth(Phist - Vbsx, 0.05, 0.1)
	sqrtPhistVbs = sqrt(PhistVbs)
	Xdep         = T1DEP * sqrtPhistVbs
	Cdep         = epssi / Xdep
	cdsc         = CIT_i + NFACTOR_t + CDSCD_a * Vdsx - CDSCB_i * Vbsx
	T1           = 1.0 + cdsc/Cox
	n = Smooth(T1, 1.0, 0.05)
	nVt     = n * Vt
	inv_nVt = 1.0 / nVt

	# Vth shift for DIBL
	dVth_dibl = -(ETA0_a + ETAB_i * Vbsx) * Vdsx
	dVth_dibl = Smooth2(dVth_dibl, 0.0, 5.0e-5)

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
	psip = PO_psip(vgfb, gam, 0.0, phib_n)

	# Normalized inversion charge at source end of channel
	qs = BSIM_q(psip, phib_n, vs, gam)

	# Average charge-surface potential slope, Ref: Charge-based MOS Transistor Modeling by C. Enz & E. Vittoz
	psipclamp = Smooth(psip, 1.0, 2.0)
	sqrtpsip = sqrt(psipclamp)

	# Source side surface potential
	psiavg = psip - 2.0 * qs
	T0 = Smooth(psiavg, 1.0, 2.0)
	nq = 1.0 + gam / (sqrtpsip + sqrt(T0))

	# Drain saturation voltage
	EeffFactor = 1.0e-8 / (epsratio * TOXE)
	T0 = nVt * (vgfb - psip - 2.0 * qs * (nq - 1.0))
	qbs = Smooth(T0, 0.0, 0.1)

	# Source side qi and qb for Vdsat- normalized to Cox
	qis = 2.0 * nq * nVt * qs
	Eeffs = EeffFactor * (qbs + eta_mu * qis)

	# Ref: BSIM4 mobility model
	T2 = pow(0.5 * (1.0 + (qis / qbs)), UCS_a)
	T3 = (UA_a + UC_a * Vbsx) * pow(Eeffs, EU_t) + UD_a / T2
	T4 = 1.0 + T3
	Dmobs = Smooth(T4, 1.0, 0.0015)
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
	Vdssat = Smooth(Vdsat - Vs, 0.0, 1.0e-3)
	T7      = pow(Vds / Vdssat , 1.0 / DELTA_t)
	T8      = pow(1.0 + T7, -DELTA_t)
	Vdseff  = Vds * T8
	vdeff   = (Vdseff + Vs) * inv_nVt
	qdeff = BSIM_q(psip, phib_n, vdeff, gam)

	# Reevaluation of nq to include qdeff
	psiavg = psip - qs - qdeff -1.0
	T0 = Smooth(psiavg, 1.0, 2.0)
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
	qba = Smooth(nVt * Qb, 0.0, 0.1)
	qia   = nVt * (Qs + Qd)

	Eeffm = EeffFactor * (qba + eta_mu * qia)
	T2    = pow(0.5 * (1.0 + (qia / qba)), UCS_a)
	T3    = (UA_a + UC_a * Vbsx) * pow(Eeffm, EU_t) + UD_a / T2
	T4    = 1.0 + T3
	Dmob = Smooth(T4, 1.0, 0.0015)

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

	# CV calculations for HVMOD
	if (RDSMOD == 1 && HVCAP == 1 && HVMOD == 1)
		vgfbdrift = -devsign * (Vg - Vd) - VFBOV
		vgfbdrift = vgfbdrift/Vt
		gamhv     = sqrt(2.0 * q * epssi * NDR * inv_Vt) / Cox
		phibHV    = lln(NDR / ni)
		psip_k = PO_psip(vgfbdrift, gamhv, 0, phibHV)
		q_k = BSIM_q(psip_k, phibHV, devsign * (Vd - Vb) / Vt, gamhv)

		# calculate nq for the drift region
		psipclamp_hv = Smooth(psip_k, 1.0, 2.0)
		sqrtpsip_k = sqrt(psipclamp_hv)
		psiavg_hv = psip_k - 2.0 * q_k
		T0 = Smooth(psiavg_hv, 1.0, 2.0)
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
			vgfbdrift = -devsign * (Vg - Vs) - VFBOV
			vgfbdrift = vgfbdrift/Vt
			psip_k = PO_psip(vgfbdrift, gamhv, 0, phibHV)
			q_k = BSIM_q(psip_k, phibHV, devsign * (Vs - Vb) / Vt, gamhv)
			psipclamp_hv = Smooth(psip_k, 1.0, 2.0)
			sqrtpsip_k = sqrt(psipclamp_hv)
			psiavg_hv = psip_k - 2.0 * q_k
			T0 = Smooth(psiavg_hv, 1.0, 2.0)
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
		Nextra = Smooth(Nextra, 0, DELTAII)
		Nextra = NDRIFTD * Nextra

		T2 = Smooth(devsign * (Vd - Vb) - Vdseff - DRII2, 0, 0.05)
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

	# Gate currents, Ref: BSIM4
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
			T5     = lexp(T4)
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

	# GIDL and GISL currents, Ref: BSIM4
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
	Qbsj1 = JunCap(Czbs, Vbs_jct, PBS_t, MJS, czbs_p1, czbs_p2)
	Qbsj2 = JunCap(Czbssw, Vbs_jct, PBSWS_t, MJSWS, czbssw_p1, czbssw_p2)
	Qbsj3 = JunCap(Czbsswg, Vbs_jct, PBSWGS_t, MJSWGS, czbsswg_p1, czbsswg_p2)
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
	Qbdj1 = JunCap(Czbd, Vbd_jct, PBD_t, MJD, czbd_p1, czbd_p2)
	Qbdj2 = JunCap(Czbdsw, Vbd_jct, PBSWD_t, MJSWD, czbdsw_p1, czbdsw_p2)
	Qbdj3 = JunCap(Czbdswg, Vbd_jct, PBSWGD_t, MJSWGD, czbdswg_p1, czbdswg_p2)
	Qbdj = Qbdj1 + Qbdj2 + Qbdj3

	# Sub-surface leakage drain current
	if (SSLMOD != 0)
		T1 = pow(NDEP_i / 1.0e23, SSLEXP1)
		T2 = pow(300.0 / DevTemp, SSLEXP2)
		T3 = (devsign * SSL5 * (Vb - Vs)) / Vt
		SSL0_NT  = SSL0 * lexp(-T1 * T2)
		SSL1_NT  = SSL1 * T2 * T1
		PHIB_SSL = SSL3 * tanh(lexp(devsign * SSL4 * ((Vg - Vb) - VTH - (Vs - Vb))))
		Issl     = sigvds * NF * Weff * SSL0_NT * lexp(T3) * lexp(-SSL1_NT * Leff) * lexp(PHIB_SSL / Vt) * (lexp(SSL2 * Vdsx / Vt) - 1.0)
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
		phib_h = log(HNDEP / ni)

		# Pinch-Off potential for halo region
		psiph = PO_psip(vgfbh, gam_h, 0.0, phib_h)

		# Normalized inversion charge at source end of halo MOSFET
		qsh = BSIM_q(psiph, phib_h, vs, gam_h)
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
		psip = PO_psip(vgfbCV, gamCV, DPD, phibCV)

		# Normalized inversion charge at source end of channel
		qs = BSIM_q(psip, phibCV, vs, gamCV)
		psipclamp = Smooth(psip, 1.0, 2.0)
		sqrtpsip = sqrt(psipclamp)

		# Source side surface potential
		psiavg = psip - 2.0 * qs
		T0 = Smooth(psiavg, 1.0, 2.0)
		nq = 1.0 + gamCV / (sqrtpsip + sqrt(T0))

		# Drain saturation voltage
		T0 = Vt * (vgfbCV - psip - 2.0 * qs * (nq - 1.0))
		qbs = Smooth(T0, 0.0, 0.1)

		# Source side qi and qb for Vdsat (normalized to Cox)
		qis = 2.0 * nq * Vt * qs
		Eeffs = EeffFactor * (qbs + eta_mu * qis)

		# Ref: BSIM4 mobility model
		T3 = (UA_a + UC_a * Vbsx) * pow(Eeffs, EU_t)
		T4 = 1.0 + T3
		Dmobs = Smooth(T4, 1.0, 0.0015)
		LambdaC_by2 = (U0_a / Dmobs) * Vt / (VSATCV_t * Lact)
		qdsat       = LambdaC_by2 * (qs * qs + qs) / (1.0 + LambdaC_by2 * (1.0 + qs))
		vdsatcv     = psip - 2.0 * phibCV - (2.0 * qdsat + lln((qdsat * 2.0 * nq * inv_gam) * ((qdsat * 2.0 * nq * inv_gam) + (gam / (nq - 1.0)))))
		VdsatCV     = vdsatcv * Vt

		# Normalized charge qdeff at drain end of channel
		VdssatCV = Smooth(VdsatCV - Vs, 0.0, 1e-3)
		VdssatCV     = VdssatCV / ABULK
		T7     = pow(Vds / VdssatCV , 1.0 / DELTA_t)
		T8     = pow(1.0 + T7, -DELTA_t)
		Vdseff = Vds * T8
		vdeff  = (Vdseff + Vs) * inv_Vt
		qdeff = BSIM_q(psip, phibCV, vdeff, gamCV)

		# Reevaluation of nq to include qdeff needed for gummel symmetry
		psiavg = psip - qs - qdeff - 1.0
		T0 = Smooth(psiavg, 1.0, 2.0)
		T2 = sqrt(T0)
		T3 = 1.0 + DPD + gamCV / (sqrtpsip + T2)
		T4 = 0.5 + DPD * T2 * inv_gam
		T5 = sqrt(T4 * T4 + T3 * (qs + qdeff) * invgamg2)
		nq = T3 / (T4 + T5)

		# C-V expressions including velocity saturation and CLM
		# Velocity saturation for C-V
		T0  = Vt * (vgfbCV - psip - 2.0 * qs * (nq - 1.0))
		qbs = Smooth(T0, 0.0, 0.1)
		T1  = Vt * (vgfbCV - psip - 2.0 * qdeff * (nq - 1.0))
		qbd = Smooth(T1, 0.0, 0.1)
		qb  = 0.5 * (qbs + qbd)
		qia = nq * Vt * (qs + qdeff)
		Eeffm = EeffFactor * (qb + eta_mu * qia)
	    psip = PO_psip((vgfbCV + DELVFBACC * inv_Vt), gamCV, DPD, phibCV)
		T3    = (UA_a + UC_a * Vbsx) * pow(Eeffm, EU_t)
		T4    = 1.0 + T3
		Dmob = Smooth(T4, 1.0, 0.0015)
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
	T1 = Smooth(sis, 0.0, 0.5)
	T2 = Smooth(sid, 0.0, 0.5)
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
	qbaCV = Smooth(Vt * Qb, 0.0, 0.1)
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
	if !("CF" in keys(param))
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
	Qovb = -devsign * NF * Lact * CGBO * (Vg - Vb)
	Qovg = -(Qovs + Qovd + Qovb)

	# Edge FET model
	if (EDGEFET == 1)
		phib_edge     = lln(NDEPEDGE_i / ni)
		Phist         = max(0.4 + Vt * phib_edge + PHIN_i, 0.4)
		sqrtPhist     = sqrt(Phist)
		T1DEP         = sqrt(2.0 * epssi / (q * NDEPEDGE_i))
		NFACTOREDGE_t = NFACTOREDGE_i * hypsmooth((1.0 + TNFACTOREDGE_i * (TRatio - 1.0)), 1e-3)
		ETA0EDGE_t    = ETA0EDGE_i * (1.0 + TETA0EDGE_i * (TRatio - 1.0))
		PhistVbs = Smooth(Phist - Vbsx, 0.05, 0.1)
		sqrtPhistVbs  = sqrt(PhistVbs)
		Xdep          = T1DEP * sqrtPhistVbs
		Cdep          = epssi / Xdep
		cdsc          = CITEDGE_i + NFACTOREDGE_t + CDSCDEDGE_i * Vdsx - CDSCBEDGE_i * Vbsx
		T1            = 1.0 + cdsc/Cox
		n = Smooth(T1, 1.0, 0.05)
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
		psip = PO_psip(vgfb, gam_edge, 0.0, phib_n_edge)
		qs_edge = BSIM_q(psip, phib_n_edge, vs, gam_edge)

		# Approximate pinch-off voltage
		vdsatedge = 2.0 * nVt * qs_edge + 2.0 * nVt
		Vdsatedge = vdsatedge
		Vdsatedge = Vdsatedge + Vs

		# Vdssat clamped to avoid negative values during transient simulation
		Vdssate = Smooth(Vdsatedge - Vs, 0.0, 1.0e-3)
		T7     = pow(Vds / Vdssate, 1.0 / DELTA_t)
		T8     = pow(1.0 + T7, -DELTA_t)
		Vdseff = Vds * T8
		vdeff  = (Vdseff + Vs) * inv_nVt
		qdeff_edge = BSIM_q(psip, phib_n_edge, vdeff, gam_edge)

		# Nq calculation for Edge FET
		psipclamp = Smooth(psip, 1.0, 2.0)
		sqrtpsip = sqrt(psipclamp)
		psiavg   = psip - qs_edge - qdeff_edge -1.0
		T0 = Smooth(psiavg, 1.0, 2.0)
		T2       = sqrt(T0)
		nq_edge  = 1.0 + gam_edge / (sqrtpsip + T2)
		ids_edge = 2.0 * NF * nq_edge * ueff * WEDGE / Leff * Cox * nVt * nVt * ((qs_edge - qdeff_edge) * (1.0 + qs_edge + qdeff_edge)) * Moc
		ids      = ids_edge + ids
	end
	return ids
end

param = Dict(
	"L" => 1e-6,
	"TOXP" => 1e-9,
	"Temp" => 125.0
	)

kk=bsimbulk(param)
println(kk)