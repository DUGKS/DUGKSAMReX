#ifndef ARK_DVM_DmVnEq_H
#define ARK_DVM_DmVnEq_H

void DUGKSAmr::Update_DVDF_EqSo_hDt(int const lev)
{	
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_EqSo_hDt()", Update_DVDF_EqSo_hDt);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase(ld.hDVDFxF[Ark::EQF], ld.hDVDFxF[Ark::SOF], 
											   			  ld.msv_xF);

	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase(ld.hDVDFyF[Ark::EQF], ld.hDVDFyF[Ark::SOF], 
											   			  ld.msv_yF);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase(ld.hDVDFzF[Ark::EQF], ld.hDVDFzF[Ark::SOF], 
											   			  ld.msv_zF);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum(ld.fDVDFxF[Ark::EQF], ld.fDVDFxF[Ark::SOF], 
												     			ld.msv_xF); //, {Alpha()}

	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum(ld.fDVDFyF[Ark::EQF], ld.fDVDFyF[Ark::SOF], 
												     			ld.msv_yF); //, {Alpha()}
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum(ld.fDVDFzF[Ark::EQF], ld.fDVDFzF[Ark::SOF], 
												     			ld.msv_zF); //, {Alpha()}
#endif
#endif

#ifdef ARK_USE_Energy
	DUGKSCore<Ark::DmVnEnergy>::Update_Eq_Energy(ld.gDVDFxF[Ark::EQF], ld.gDVDFxF[Ark::SOF], 
												     		ld.msv_xF);

	DUGKSCore<Ark::DmVnEnergy>::Update_Eq_Energy(ld.gDVDFyF[Ark::EQF], ld.gDVDFyF[Ark::SOF], 
												     		ld.msv_yF);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Eq_Energy(ld.gDVDFzF[Ark::EQF], ld.gDVDFzF[Ark::SOF], 
												     		ld.msv_zF);
#endif
#endif	
}

void DUGKSAmr::Update_DVDF_Eq_hDt(int const lev)
{	
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Eq_hDt()", Update_DVDF_Eq_hDt);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase<false> (ld.hDVDFxF[Ark::EQF], ld.hDVDFxF[Ark::SOF], 
											   			  ld.msv_xF);

	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase<false> (ld.hDVDFyF[Ark::EQF], ld.hDVDFyF[Ark::SOF], 
											   			  ld.msv_yF);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase<false> (ld.hDVDFzF[Ark::EQF], ld.hDVDFzF[Ark::SOF], 
											   			  ld.msv_zF);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum<false> (ld.fDVDFxF[Ark::EQF], ld.fDVDFxF[Ark::SOF], 
												     			ld.msv_xF); //, {Alpha()}

	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum<false> (ld.fDVDFyF[Ark::EQF], ld.fDVDFyF[Ark::SOF], 
												     			ld.msv_yF); //, {Alpha()}
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum<false> (ld.fDVDFzF[Ark::EQF], ld.fDVDFzF[Ark::SOF], 
												     			ld.msv_zF); //, {Alpha()}
#endif
#endif

#ifdef ARK_USE_Energy
	DUGKSCore<Ark::DmVnEnergy>::Update_Eq_Energy<false> (ld.gDVDFxF[Ark::EQF], ld.gDVDFxF[Ark::SOF], 
												     		ld.msv_xF);

	DUGKSCore<Ark::DmVnEnergy>::Update_Eq_Energy<false> (ld.gDVDFyF[Ark::EQF], ld.gDVDFyF[Ark::SOF], 
												     		ld.msv_yF);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Eq_Energy<false> (ld.gDVDFzF[Ark::EQF], ld.gDVDFzF[Ark::SOF], 
												     		ld.msv_zF);
#endif
#endif	
}

void DUGKSAmr::Update_DVDF_EqSo (int const lev)
{	
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_EqSo()", Update_DVDF_EqSo);

	auto &ld = Data(lev);
	
#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase (ld.hDVDF[Ark::EQ], ld.hDVDF[Ark::SO],ld.msv);
#endif

#ifdef ARK_USE_MOMENTUM	
	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum (ld.fDVDF[Ark::EQ], ld.fDVDF[Ark::SO], 
																ld.msv); //, {Alpha()}
#endif

#ifdef ARK_USE_ENERGY	
	DUGKSCore<Ark::DmVnEnergy>::Update_Eq_Energy (ld.gDVDF[Ark::EQ], ld.gDVDF[Ark::SO], ld.msv);
#endif
}

void DUGKSAmr::Update_DVDF_Eq (int const lev)
{	
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Eq()", Update_DVDF_Eq);

	auto &ld = Data(lev);
	
#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase<false> (ld.hDVDF[Ark::EQ], ld.hDVDF[Ark::SO],ld.msv);
#endif

#ifdef ARK_USE_MOMENTUM	
	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum<false> (ld.fDVDF[Ark::EQ], ld.fDVDF[Ark::SO], 
																ld.msv); //, {Alpha()}
#endif

#ifdef ARK_USE_ENERGY	
	DUGKSCore<Ark::DmVnEnergy>::Update_Eq_Energy<false> (ld.gDVDF[Ark::EQ], ld.gDVDF[Ark::SO], ld.msv);
#endif
}

#endif
