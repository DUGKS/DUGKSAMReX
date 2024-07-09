#ifndef ARK_DUGKS_CORE_H
#define ARK_DUGKS_CORE_H

#include <DUGKSNameSpace.H>

template<typename DVM>
class DUGKSCore
{

public:
	template<bool ENABLESOURCE>
	static void Construct_Prime (Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hdTime);

	template<bool ENABLESOURCE>
	static void Construct_NonEq (Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hdTime);

	template<bool ENABLESOURCE>
	static void Construct_Wave (Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hdTime);

    static void Construct_WavePlus (MultiFab &xWaveFab, const MultiFab &xBarPlusFab,
									const Real hTime, const Real dTime);

	static void Transform (MultiFab &TargetFab, const double (&M)[DVM::Q][DVM::Q]);

	// - construct \bar{f}^{+}
	static void Collision_BarPlus (Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hTime, Real const dTime);

	// - construct \bar{f}^{+} from f^{NEQ}
	static void Collision_BarPlus (Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hTime);

	static void Update_Grad_BarPlus (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real *dx);

	static void Update_Grad_BarPlus_2nd (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real *dx);

	static void Update_Grad_BarPlus_3rd (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real *dx);

	static void Update_Grad_BarPlus_4th (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real *dx);

	static void Update_Grad_BarPlus_5th (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real *dx);

	static void Update_Grad_BarPlus_VanLeer (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real *dx);

	static void Update_Barh_UW_1st (MultiFab &xBarhFabnF,
									Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
									const Real* dx, const Real hTime,
									amrex::Direction FaceID);

	static void Update_Barh_UW_2nd (MultiFab &xBarhFabnF,
									Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
									const Real* dx, const Real hTime,
									amrex::Direction FaceID);

	static void Update_Barh_UW_3rd (MultiFab &xBarhFab,
									Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
						 			const Real* dx, const Real hTime,
						 			amrex::Direction FaceID);

	static void Update_Barh_UW_4th (MultiFab &xBarhFab,
									Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
						 			const Real* dx, const Real hTime,
						 			amrex::Direction FaceID);

	static void Update_Barh_UW_4th_cross (MultiFab &xBarhFab,
										  Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
						 				  const Real* dx, const Real hTime,
						 				  amrex::Direction FaceID);

	static void Update_Barh_UW_5th (MultiFab &xBarhFab,
									Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
						 			const Real* dx, const Real hTime,
						 			amrex::Direction FaceID);

	static void Update_Barh_UW_Limiter (MultiFab &xBarhFab,
										Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
						 				amrex::Direction FaceID);

	static void Update_Barh_CD (MultiFab &xBarhFab, const MultiFab &xBarPlusFab,
						 		const Real* dx, const Real hTime,
						 		amrex::Direction FaceID);

	static void Update_Barh_WENO_3rd (MultiFab &xBarhFab, const MultiFab &xBarPlusFab,
						 			  const Real* dx, const Real hTime,
						 			  amrex::Direction FaceID);

	static void Update_Barh_WENO_5th (MultiFab &xBarhFab, const MultiFab &xBarPlusFab,
						 			  const Real* dx, const Real hTime,
						 			  amrex::Direction FaceID);

	static void Update_FD_Barh_WENO_3rd (MultiFab &xBarhFab, const MultiFab &xBarPlusFab,
								  	 	 amrex::Direction FaceID);

	static void Update_FD_Barh_WENO_5th (MultiFab &xBarhFab, const MultiFab &xBarPlusFab,
								  	 	 amrex::Direction FaceID);

	template<bool ENABLESOURCE>
	static void Collision_xFluxX (Array<MultiFab, Ark::FACEDVDFNUM> &xDVDFnF, const Real hTime);

	static void Update_Wave (MultiFab &xWaveFab, const MultiFab &xBarPlusFab,
							 const Vector<MultiFab const *> &xFluxFabnF,
							 const Real *dx, const Real hTime, const Real dTime);

	static void Update_Wave (MultiFab &xWaveFab, const Vector<MultiFab const *> &xFluxFabnF,
							 const Real *dx, const Real dTime);

	static void Update_Wave_SS (MultiFab &xWaveFab, const MultiFab &xSoFab,
							 	const Real hTime);

	//../discrete_velocity_model/DUGKSDmVn***.cpp
	template<bool ENABLESOURCE=true>
	static void Update_Eq_Phase (MultiFab &hEq, MultiFab &hSo, const Array<MultiFab, Ark::MSVNUM> &msv);

	template<bool ENABLESOURCE=true>
	static void Update_Eq_Momentum (MultiFab &fEq, MultiFab &fSo, const Array<MultiFab, Ark::MSVNUM> &msv);

	//used for strang-splitting
	template<bool EnableSource=true>
	static void Update_Eq_Momentum (MultiFab &fEq, MultiFab &fSo, const Array<MultiFab, Ark::MSVNUM> &msv,
									const amrex::Vector<double>& magic);

	template<bool EnableSource=true>
	static void Update_Eq_Energy (MultiFab &gEq, MultiFab &gSo, const Array<MultiFab, Ark::MSVNUM> &msv);

	static void Integral_MacroBaseVars_Phase (Array<MultiFab, Ark::MSVNUM> &msv,
										      const MultiFab &xDVDF);

	static void Integral_MacroBaseVars_Momentum (Array<MultiFab, Ark::MSVNUM> &msv,
										   	     const MultiFab &xDVDF);
	static void Integral_MacroBaseVars_Momentum_Face (Array<MultiFab, Ark::MSVNUM> &msv,
										   	     const MultiFab &xDVDF);
	static void Integral_MacroBaseVars_Energy (Array<MultiFab, Ark::MSVNUM> &msv,
										   	   const MultiFab &xDVDF);
};
template<typename DVM>
class DUGKSBoundary
{
public:
	static void FillPhysicalBoundaryFlux_Eq  (const Vector<MultiFab*> &xFluxFabPtr,
											  const Vector<const MultiFab*> &xEqFabPtr,
											  int const lo_min, int const hi_max, int const dir);

	static void FillPhysicalBoundaryFlux_NEE (const Vector<MultiFab*> &xFluxFabPtr,
											  const Vector<const MultiFab*> &xEqFabPtr,
											  int const lo_min, int const hi_max, int const dir);


	static void FillPhysicalBoundaryFlux_NEECF (const Vector<MultiFab*> &xFluxFabPtr,
												const Vector<const MultiFab*> &xEqFabPtr,
												const Array<MultiFab, Ark::DVDFNUM> &xDVDF,
												int const lo_min, int const hi_max,
												Real const dTime, int const dir);

	static void FillPhysicalBoundaryFlux_BB (const Vector<MultiFab*> &xFluxFabPtr,
									  		 const Vector<const MultiFab*> &xVarFabPtr,
									  		 const Vector<const MultiFab*> &VelFabPtr,
									  		 int const lo_min, int const hi_max,
									  		 int const dir);
private:
	//This function is only valid for the equilibrium function based upon Taylor expansion,
	//which means the equilibrium distribution function should contain the typical
	//first order term of LaTex:\frac{(\xi\cdot{u})}{RT}
	static void BounceBack (Array4<Real> const &xFlux,
                 			Array4<const Real> const &/*conservative variable : phi,rho,temp */,
                 			Array4<const Real> const &/*velocity*/,
                 			int const i, int const j, int const k, int const n);
};

template<typename DVM>
class DVMCore
{
public:
	static void Update_Wave (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real dTime,
							 const Real minL);
};

template<typename DVM>
class DVMBoundary
{
public:

	static void FillPhysicalBoundaryWave (Array<MultiFab, Ark::DVDFNUM> &xDVDF,
										  const MultiFab& PhiFab,
										  double const dTime, double const minL,
										  int const lo_min, int const hi_max,
									  	  int const dir);
};

#include <DUGKSCoreImp.H>
#endif