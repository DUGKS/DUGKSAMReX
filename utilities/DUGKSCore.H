#ifndef ARK_DUGKS_CORE_H
#define ARK_DUGKS_CORE_H

#include <DUGKSNameSpace.H>

template<typename DVM>
class DUGKSCore
{

public:
	static constexpr double aT = 4.0/3.0; 

	static constexpr double bT = 1.0/3.0;

	//DUGKSCoreImp.H
	static void Update_BarPlus (Array<MultiFab, Ark::DVDFNUM> &xDVDF, 
								Real const dTime, Real const hTime);

	static void Update_Grad_BarPlus (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real *dx);

	static void Update_xFluxX (Array<MultiFab, Ark::FACEDVDFNUM> &xDVDFnF, const Real hTime);

	static void Update_Barh_UW (MultiFab &xBarhFab, 
								Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
						 		const Real* dx, const Real hTime,
						 		amrex::Direction FaceID);
	
	static void Update_Barh_CD (MultiFab &xBarhFab, const MultiFab &xBarPlusFab,
						 		const Real* dx, const Real hTime,
						 		amrex::Direction FaceID);
	
	static void Update_Wave (MultiFab &xWaveFab, const MultiFab &xBarPlusFab,
							 const Vector<MultiFab const *> &xFluxFabnF,
							 const Real *dx, const Real dTime);

	//../discrete_velocity_model/DUGKSDmVn***.cpp
	template<long unsigned NUM>
	static void Update_Eq_Phase (MultiFab &hEq, MultiFab &hSo, 
								 const Array<MultiFab, NUM> &msv);

	template<long unsigned NUM>
	static void Update_Eq_Momentum (MultiFab &fEq, MultiFab &fSo,
								 	const Array<MultiFab, NUM> &msv);

	template<long unsigned NUM>
	static void Update_Eq_Energy (MultiFab &gEq, MultiFab &gSo, 
								  const Array<MultiFab, NUM> &msv);

	static void Update_MacroBaseVars_Phase(Array<MultiFab, Ark::MSVNUM> &msv, 
										   const MultiFab &xDVDF);

	static void Update_MacroBaseVars_Momentum(Array<MultiFab, Ark::MSVNUM> &msv, 
										   	  const MultiFab &xDVDF);

	static void Update_MacroBaseVars_Energy(Array<MultiFab, Ark::MSVNUM> &msv, 
										   	  const MultiFab &xDVDF);
};
template<typename DVM>
class DUGKSBoundary
{
public:
	static void FillPhysicalBoundaryFlux_NEE (const Vector<MultiFab*> &xFluxFabPtr,
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
#include <DUGKSCoreImp.H>
#endif