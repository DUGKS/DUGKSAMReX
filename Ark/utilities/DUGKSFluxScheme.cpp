#include <DUGKSAmr.H>

void DUGKSAmr::Update_DVDF_Grad_BarPlus_Null (int const lev)
{
	return;
}

void DUGKSAmr::Update_DVDF_Grad_BarPlus_Default (int const lev)
{
BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Grad_BarPlus_Default()", Update_DVDF_Grad_BarPlus_Default);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Grad_BarPlus (ld.hDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Grad_BarPlus (ld.fDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Grad_BarPlus (ld.gDVDF, Geom(lev).CellSize());
#endif	
}

void DUGKSAmr::Update_DVDF_Grad_BarPlus_Third (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Grad_BarPlus_Third()", Update_DVDF_Grad_BarPlus_Third);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Grad_BarPlus_3rd(ld.hDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Grad_BarPlus_3rd(ld.fDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Grad_BarPlus_3rd(ld.gDVDF, Geom(lev).CellSize());
#endif	
}

void DUGKSAmr::Update_DVDF_Grad_BarPlus_Fourth (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Grad_BarPlus_Fourth()", Update_DVDF_Grad_BarPlus_Fourth);	

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Grad_BarPlus_Fourth_PHASE()", Grad_BarPlus_Fourth_PHASE);
	DUGKSCore<Ark::DmVnPhase>::Update_Grad_BarPlus_4th(ld.hDVDF, Geom(lev).CellSize());
BL_PROFILE_VAR_STOP(Grad_BarPlus_Fourth_PHASE);

#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Grad_BarPlus_4th(ld.fDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Grad_BarPlus_4th(ld.gDVDF, Geom(lev).CellSize());
#endif	
}

void DUGKSAmr::Update_DVDF_Grad_BarPlus_Fifth (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Grad_BarPlus_Fifth()", Update_DVDF_Grad_BarPlus_Fifth);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Grad_BarPlus_5th(ld.hDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Grad_BarPlus_5th(ld.fDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Grad_BarPlus_5th(ld.gDVDF, Geom(lev).CellSize());
#endif	
}

void DUGKSAmr::Update_DVDF_Grad_BarPlus_Fourth_Third (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Grad_BarPlus_Fourth_Third()", Grad_BarPlus_Fourth_Third);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Grad_BarPlus_4th(ld.hDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Grad_BarPlus_3rd(ld.fDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Grad_BarPlus_3rd(ld.gDVDF, Geom(lev).CellSize());
#endif	
}

void DUGKSAmr::Update_DVDF_Grad_BarPlus_Default_Second (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Grad_BarPlus_Default_Second()", Grad_BarPlus_Default_Second);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Grad_BarPlus(ld.hDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Grad_BarPlus_2nd(ld.fDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Grad_BarPlus_2nd(ld.gDVDF, Geom(lev).CellSize());
#endif	
}

void DUGKSAmr::Update_DVDF_Grad_BarPlus_Default_VanLeer (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Grad_BarPlus_Default_VanLeer()", Grad_BarPlus_VanLeer);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Grad_BarPlus(ld.hDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Grad_BarPlus_VanLeer(ld.fDVDF, Geom(lev).CellSize());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Grad_BarPlus_VanLeer(ld.gDVDF, Geom(lev).CellSize());
#endif	
}

void DUGKSAmr::Update_DVDF_Barh_UW_Default (int const lev)
{
	Update_DVDF_Barh_UW_Third (lev);
}

void DUGKSAmr::Update_DVDF_Barh_UW_Third (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Barh_UW()", Update_DVDF_Barh_UW);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_3rd (ld.hDVDFxF[Ark::FLUX], ld.hDVDF,
						 				  			Geom(lev).CellSize(), hDt(lev),
						 				  			amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_3rd (ld.hDVDFyF[Ark::FLUX], ld.hDVDF,
						 				  			Geom(lev).CellSize(), hDt(lev),
						 				  			amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_3rd (ld.hDVDFzF[Ark::FLUX], ld.hDVDF,
						 				  			Geom(lev).CellSize(), hDt(lev),
						 				  			amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_3rd (ld.fDVDFxF[Ark::FLUX], ld.fDVDF,
						 				  	 	  	  Geom(lev).CellSize(), hDt(lev),
						 				  	 	  	  amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_3rd (ld.fDVDFyF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_3rd (ld.fDVDFzF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_3rd (ld.gDVDFxF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_3rd (ld.gDVDFyF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_3rd (ld.gDVDFzF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::z);
#endif

#endif
}

void DUGKSAmr::Update_DVDF_Barh_UW_Fourth (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Barh_UW_Fourth()", Update_DVDF_Barh_UW_Fourth);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_4th (ld.hDVDFxF[Ark::FLUX], ld.hDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_4th (ld.hDVDFyF[Ark::FLUX], ld.hDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  		 	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_4th (ld.hDVDFzF[Ark::FLUX], ld.hDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  		 	amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_4th (ld.fDVDFxF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_4th (ld.fDVDFyF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_4th (ld.fDVDFzF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_4th (ld.gDVDFxF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_4th (ld.gDVDFyF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_4th (ld.gDVDFzF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::z);
#endif

#endif	
}
void DUGKSAmr::Update_DVDF_Barh_UW_Fourth_Cross (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Barh_UW_Fourth_Cross()", Update_DVDF_Barh_UW_Fourth_Cross);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_4th_cross (ld.hDVDFxF[Ark::FLUX], ld.hDVDF,
						 				  	 			Geom(lev).CellSize(), hDt(lev),
						 				  	 		  	amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_4th_cross (ld.hDVDFyF[Ark::FLUX], ld.hDVDF,
						 				  	 		  	Geom(lev).CellSize(), hDt(lev),
						 				  		 	  	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_4th_cross (ld.hDVDFzF[Ark::FLUX], ld.hDVDF,
						 				  	 		  	Geom(lev).CellSize(), hDt(lev),
						 				  		 	  	amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_4th_cross (ld.fDVDFxF[Ark::FLUX], ld.fDVDF,
						 				  	 				Geom(lev).CellSize(), hDt(lev),
						 				  	 				amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_4th_cross (ld.fDVDFyF[Ark::FLUX], ld.fDVDF,
						 				  	 				Geom(lev).CellSize(), hDt(lev),
						 				  	 				amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_4th_cross (ld.fDVDFzF[Ark::FLUX], ld.fDVDF,
						 				  	 				Geom(lev).CellSize(), hDt(lev),
						 				  	 				amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_4th_cross (ld.gDVDFxF[Ark::FLUX], ld.gDVDF,
						 				  	 			  Geom(lev).CellSize(), hDt(lev),
						 				  	 			  amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_4th_cross (ld.gDVDFyF[Ark::FLUX], ld.gDVDF,
						 				  	 			  Geom(lev).CellSize(), hDt(lev),
						 				  	 			  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_4th_cross (ld.gDVDFzF[Ark::FLUX], ld.gDVDF,
						 				  	 			  Geom(lev).CellSize(), hDt(lev),
						 				  	 			  amrex::Direction::z);
#endif

#endif	
}
void DUGKSAmr::Update_DVDF_Barh_UW_Fifth (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Barh_UW_Fifth()", Update_DVDF_Barh_UW_Fifth);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_5th (ld.hDVDFxF[Ark::FLUX], ld.hDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_5th (ld.hDVDFyF[Ark::FLUX], ld.hDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  		 	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_5th (ld.hDVDFzF[Ark::FLUX], ld.hDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  		 	amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_5th (ld.fDVDFxF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_5th (ld.fDVDFyF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)		
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_5th (ld.fDVDFzF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_5th (ld.gDVDFxF[Ark::FLUX], ld.gDVDF,
						 				  			Geom(lev).CellSize(), hDt(lev),
						 				  			amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_5th (ld.gDVDFyF[Ark::FLUX], ld.gDVDF,
						 				  			Geom(lev).CellSize(), hDt(lev),
						 				  			amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_5th (ld.gDVDFzF[Ark::FLUX], ld.gDVDF,
						 				  			Geom(lev).CellSize(), hDt(lev),
						 				  			amrex::Direction::z);
#endif

#endif	
}
void DUGKSAmr::Update_DVDF_Barh_UW_Fourth_Cross_Third (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Barh_UW_Fourth_Cross_Third()", UW_Fourth_Cross_Third);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_4th_cross (ld.hDVDFxF[Ark::FLUX], ld.hDVDF,
						 				  	 		  	Geom(lev).CellSize(), hDt(lev),
						 				  	 		  	amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_4th_cross (ld.hDVDFyF[Ark::FLUX], ld.hDVDF,
						 				  	 		  	Geom(lev).CellSize(), hDt(lev),
						 				  		 	  	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_4th_cross (ld.hDVDFzF[Ark::FLUX], ld.hDVDF,
						 				  	 		  	Geom(lev).CellSize(), hDt(lev),
						 				  		 	  	amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_3rd (ld.fDVDFxF[Ark::FLUX], ld.fDVDF,
						 				  	 	  	  Geom(lev).CellSize(), hDt(lev),
						 				  	 	  	  amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_3rd (ld.fDVDFyF[Ark::FLUX], ld.fDVDF,
						 				  	 	  	  Geom(lev).CellSize(), hDt(lev),
						 				  	 	  	  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_3rd (ld.fDVDFzF[Ark::FLUX], ld.fDVDF,
						 				  	 	  	  Geom(lev).CellSize(), hDt(lev),
						 				  	 	  	  amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_3rd (ld.gDVDFxF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_3rd (ld.gDVDFyF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_3rd (ld.gDVDFzF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::z);
#endif

#endif		
}

void DUGKSAmr::Update_DVDF_Barh_UW_Default_Second (int const lev)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Barh_UW_Fourth_Cross_Third()", UW_Fourth_Cross_Third);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_3rd (ld.hDVDFxF[Ark::FLUX], ld.hDVDF,
						 				  	   		Geom(lev).CellSize(), hDt(lev),
						 				  	   		amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_3rd (ld.hDVDFyF[Ark::FLUX], ld.hDVDF,
						 				  	   		Geom(lev).CellSize(), hDt(lev),
						 				  	   		amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW_3rd (ld.hDVDFzF[Ark::FLUX], ld.hDVDF,
						 				  	   		Geom(lev).CellSize(), hDt(lev),
						 				  	   		amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_2nd (ld.fDVDFxF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_2nd (ld.fDVDFyF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW_2nd (ld.fDVDFzF[Ark::FLUX], ld.fDVDF,
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_2nd (ld.gDVDFxF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_2nd (ld.gDVDFyF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW_2nd (ld.gDVDFzF[Ark::FLUX], ld.gDVDF,
						 				  	 		Geom(lev).CellSize(), hDt(lev),
						 				  	 		amrex::Direction::z);
#endif

#endif
}
void DUGKSAmr::Update_DVDF_Barh_CD (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Barh_CD()", Update_DVDF_Barh_CD);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_CD (ld.hDVDFxF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  		Geom(lev).CellSize(), hDt(lev),
						 				  		amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_CD (ld.hDVDFyF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  		Geom(lev).CellSize(), hDt(lev),
						 				  		amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_CD (ld.hDVDFzF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  		Geom(lev).CellSize(), hDt(lev),
						 				  		amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_CD (ld.fDVDFxF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 	  Geom(lev).CellSize(), hDt(lev),
						 				  	 	  amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_CD (ld.fDVDFyF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 	  Geom(lev).CellSize(), hDt(lev),
						 				  	 	  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_CD (ld.fDVDFzF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 	  Geom(lev).CellSize(), hDt(lev),
						 				  	 	  amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_CD (ld.gDVDFxF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 	Geom(lev).CellSize(), hDt(lev),
						 				  	 	amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_CD (ld.gDVDFyF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 	Geom(lev).CellSize(), hDt(lev),
						 				  	 	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_CD (ld.gDVDFzF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 	Geom(lev).CellSize(), hDt(lev),
						 				  	 	amrex::Direction::z);
#endif

#endif
}

void DUGKSAmr::Update_DVDF_Barh_WENO_Third (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Barh_WENO_Third()", Update_DVDF_Barh_WENO_Third);

	auto &ld = Data(lev);
	
#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_WENO_3rd (ld.hDVDFxF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  	   		Geom(lev).CellSize(), hDt(lev),
						 				  	   		amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_WENO_3rd (ld.hDVDFyF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  	   		Geom(lev).CellSize(), hDt(lev),
						 				  	   		amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_WENO_3rd (ld.hDVDFzF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  	   		Geom(lev).CellSize(), hDt(lev),
						 				  	   		amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_WENO_3rd (ld.fDVDFxF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 			Geom(lev).CellSize(), hDt(lev),
						 				  	 			amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_WENO_3rd (ld.fDVDFyF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 			Geom(lev).CellSize(), hDt(lev),
						 				  	 			amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_WENO_3rd (ld.fDVDFzF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 			Geom(lev).CellSize(), hDt(lev),
						 				  	 			amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_WENO_3rd (ld.gDVDFxF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_WENO_3rd (ld.gDVDFyF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  			  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_WENO_3rd (ld.gDVDFzF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::z);
#endif

#endif
}

void DUGKSAmr::Update_DVDF_Barh_WENO_Fifth (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Barh_WENO_Fifth()", Update_DVDF_Barh_WENO_Fifth);

	auto &ld = Data(lev);
	
#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_WENO_5th (ld.hDVDFxF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  	   		Geom(lev).CellSize(), hDt(lev),
						 				  	   		amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_WENO_5th (ld.hDVDFyF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  	   		Geom(lev).CellSize(), hDt(lev),
						 				  	   		amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_WENO_5th (ld.hDVDFzF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  	   		Geom(lev).CellSize(), hDt(lev),
						 				  	   		amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_WENO_5th (ld.fDVDFxF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 			Geom(lev).CellSize(), hDt(lev),
						 				  	 			amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_WENO_5th (ld.fDVDFyF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 			Geom(lev).CellSize(), hDt(lev),
						 				  	 			amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_WENO_5th (ld.fDVDFzF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 			Geom(lev).CellSize(), hDt(lev),
						 				  	 			amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_WENO_5th (ld.gDVDFxF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_WENO_5th (ld.gDVDFyF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  			  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_WENO_5th (ld.gDVDFzF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 		  Geom(lev).CellSize(), hDt(lev),
						 				  	 		  amrex::Direction::z);
	#endif

#endif
}

void DUGKSAmr::Update_FD_DVDF_Barh_WENO_Third (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_FD_DVDF_Barh_WENO_Third()", Update_FD_DVDF_Barh_WENO_Third);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_FD_Barh_WENO_3rd (ld.hDVDFxF[Ark::FLUX], 
														ld.hDVDF[Ark::WAVE],
						 				  		 	  	amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_FD_Barh_WENO_3rd (ld.hDVDFyF[Ark::FLUX], 
														ld.hDVDF[Ark::WAVE],
						 				  		 	  	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_FD_Barh_WENO_3rd (ld.hDVDFzF[Ark::FLUX], 
														ld.hDVDF[Ark::WAVE],
						 				  		 	  	amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_3rd (ld.fDVDFxF[Ark::FLUX], 
															ld.fDVDF[Ark::WAVE],
						 				  	 		     	amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_3rd (ld.fDVDFyF[Ark::FLUX], 
															ld.fDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_3rd (ld.fDVDFzF[Ark::FLUX], 
															ld.fDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_3rd (ld.gDVDFxF[Ark::FLUX], 
															ld.gDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_3rd (ld.gDVDFyF[Ark::FLUX], 
															ld.gDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_3rd (ld.gDVDFzF[Ark::FLUX], 
															ld.gDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::z);
#endif

#endif
}

void DUGKSAmr::Update_FD_DVDF_Barh_WENO_Fifth (int const lev)
{

BL_PROFILE_VAR("DUGKSAmr::Update_FD_DVDF_Barh_WENO_Fifth()", Update_FD_DVDF_Barh_WENO_Fifth);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_FD_Barh_WENO_5th (ld.hDVDFxF[Ark::FLUX],
														ld.hDVDF[Ark::WAVE],
						 				  		 	  	amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_FD_Barh_WENO_5th (ld.hDVDFyF[Ark::FLUX], 
														ld.hDVDF[Ark::WAVE],
						 				  		 	    amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_FD_Barh_WENO_5th (ld.hDVDFzF[Ark::FLUX], 
														ld.hDVDF[Ark::WAVE],
						 				  		 	  	amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_5th (ld.fDVDFxF[Ark::FLUX], 
															ld.fDVDF[Ark::WAVE],
						 				  	 		     	amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_5th (ld.fDVDFyF[Ark::FLUX], 
															ld.fDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_5th (ld.fDVDFzF[Ark::FLUX], 
															ld.fDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_5th (ld.gDVDFxF[Ark::FLUX], 
															ld.gDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_5th (ld.gDVDFyF[Ark::FLUX], 
															ld.gDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_FD_Barh_WENO_5th (ld.gDVDFzF[Ark::FLUX], 
															ld.gDVDF[Ark::WAVE],
						 				  	 			 	amrex::Direction::z);
#endif

#endif
}
