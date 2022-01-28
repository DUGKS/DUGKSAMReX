#include <DUGKSAmr.H>
#include <DmVnAC.H>

#error "Aborted due to inconsistency between PHASE and MOMENTUM kinetic equation"

#if(AMREX_SPACEDIM != 2)
#error "AMREX_SPACEDIM is expected to be 2!!!"
#endif

template <typename DVM>
template <long unsigned NUM>
void DUGKSCore<DVM>::Update_Eq_Phase (MultiFab &xEq, MultiFab &xSo,
									  const Array<MultiFab, NUM> &msv)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xEq, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const amrex::Box& tbx = mfi.tilebox();

		Array4<Real const> const &phi_array     = msv[Ark::PHI].const_array(mfi);
		Array4<Real const> const &vel_array     = msv[Ark::VEL].const_array(mfi);
		Array4<Real const> const &gradphi_array = msv[Ark::GradPHI].const_array(mfi);

		Array4<Real> const &hEq = xEq.array(mfi);
		Array4<Real> const &hSo = xSo.array(mfi);

		ParallelFor(tbx, DVM::Q,
		[=] AMREX_GPU_DEVICE (int i,int j,int k,int n) noexcept
		{
			const Real phi  = phi_array(i,j,k);
			const Real u    = vel_array(i,j,k, Ark::X);
			const Real v    = vel_array(i,j,k, Ark::Y);
			const Real phix = gradphi_array(i,j,k, Ark::X);
			const Real phiy = gradphi_array(i,j,k, Ark::Y);

			hEq(i,j,k,n) = DVM::omega[n]*phi*(1.0 + (DVM::xi_u[n]*u + DVM::xi_v[n]*v)/DVM::RT);
			hSo(i,j,k,n) = DVM::omega[n]*(DVM::xi_u[n]*phix + DVM::xi_v[n]*phiy);
		});
	}
}

template <typename DVM>
template <long unsigned NUM>
void DUGKSCore<DVM>::Update_Eq_Momentum (MultiFab &xEq, MultiFab &xSo, 
										 const Array<MultiFab, NUM> &msv)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xEq, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const amrex::Box& tbx = mfi.tilebox();
	//
		Array4<Real const> const &Den     = msv[Ark::RHO].const_array(mfi);
		Array4<Real const> const &Vel     = msv[Ark::VEL].const_array(mfi);
		Array4<Real const> const &Pres    = msv[Ark::P].const_array(mfi);
		Array4<Real const> const &GradDen = msv[Ark::GradRHO].const_array(mfi);
		Array4<Real const> const &Force   = msv[Ark::FORCE].const_array(mfi);

		Array4<Real> const &fEq = xEq.array(mfi);
		Array4<Real> const &fSo = xSo.array(mfi);

		ParallelFor(tbx, DVM::Q,
		[=] AMREX_GPU_DEVICE (int i,int j,int k,int n) noexcept
		{
       		const double rho  = Den(i,j,k);
       		const double p    = Pres(i,j,k);
			const double u    = Vel(i,j,k,Ark::X);
     		const double v    = Vel(i,j,k,Ark::Y);
       		const double rhox = GradDen(i,j,k,Ark::X);
       		const double rhoy = GradDen(i,j,k,Ark::Y);
       		const double fx   = Force(i,j,k,Ark::X);
       		const double fy   = Force(i,j,k,Ark::Y);

			const double uu = Ark::Square(u,v);

			const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v)/DVM::RT;
			const double Gamma = u1 + 0.5*u1*u1 - uu*DVM::Lambda;

			fEq(i,j,k,n) = DVM::omega[n]*(p/DVM::RT + rho*Gamma);

			#ifdef _ARK_M_LSGC_FLIP
			fEq(i,j,k,n) -= static_cast<int>(n == 0)*p/DVM::RT;
			#endif

			// fSo(i,j,k,n) = DVM::omega[n]*
			// 		   (
			// 		   (1+Gamma)*((DVM::xi_u[n] - u)*fx   + (DVM::xi_v[n] - v)*fy)/DVM::RT
			// 		   +
			// 		   Gamma*((DVM::xi_u[n] - u)*rhox + (DVM::xi_v[n] - v)*rhoy)
			// 		   );

			// Me
			fSo(i,j,k,n) = DVM::omega[n]*
			(
			  ((DVM::xi_u[n] - u)*fx + (DVM::xi_v[n] - v)*fy)
			+ (DVM::xi_u[n]*(fx + rhox*DVM::RT) + DVM::xi_v[n]*(fy + rhoy*DVM::RT))*u1
			)/DVM::RT;

			//Liang-2018
			// fSo(i,j,k,n) = DVM::omega[n]*
			// (
			//   (DVM::xi_u[n]*fx + DVM::xi_v[n]*fy)
			// + (DVM::xi_u[n]*rhox + DVM::xi_v[n]*rhoy)*u1*DVM::RT
			// )/DVM::RT;
		});
	}	
}

void DUGKSAmr::Update_MacroBaseVars (Array<MultiFab, Ark::MSVNUM> &msv)
{
#ifdef ARK_USE_MOMENTUM
	const MultiFab &RhoFab   = msv[Ark::RHO];
	
	MultiFab &VelFab   = msv[Ark::VEL];
	MultiFab &PresFab  = msv[Ark::P];

	PresFab.mult(Ark::DmVnMomentum::RT, 0);

	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &rho = RhoFab.array(mfi);

		Array4<Real> const &vel   = VelFab.array(mfi);

		amrex::ParallelFor(bx, AMREX_SPACEDIM,
		[=] AMREX_GPU_DEVICE (int i,int j,int k, int n) noexcept
		{
			vel(i,j,k,n) /= rho(i,j,k);
		});
	}
#endif
}
void DUGKSAmr::Update_MacroSourceVars (Array<MultiFab, Ark::MSVNUM> &msv, const Real interval)
{
#ifdef ARK_USE_MOMENTUM

	MultiFab const &RhoFab     = msv[Ark::RHO];
	MultiFab const &GradRhoFab = msv[Ark::GradRHO];
	MultiFab const &ForceFab   = msv[Ark::FORCE];

	MultiFab &VelFab  = msv[Ark::VEL];
	MultiFab &PresFab = msv[Ark::P];

	for(amrex::MFIter mfi(PresFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &rho 	   = RhoFab.array(mfi);
		Array4<Real const> const &grad_rho = GradRhoFab.const_array(mfi);
		Array4<Real const> const &force    = ForceFab.array(mfi);

		Array4<Real> const &vel   = VelFab.array(mfi);
		Array4<Real> const &pres  = PresFab.array(mfi);

		amrex::ParallelFor(bx, AMREX_SPACEDIM,[=] AMREX_GPU_DEVICE (int i,int j,int k, int n) noexcept
		{
			vel(i,j,k,n) += interval*force(i,j,k,n)/rho(i,j,k);
		});

		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
		{
			pres(i,j,k) +=  interval*Ark::DmVnMomentum::RT*
						   	(
								vel(i,j,k,Ark::X)*grad_rho(i,j,k,Ark::X)
							+	vel(i,j,k,Ark::Y)*grad_rho(i,j,k,Ark::Y)
							);
		});		
	}

#endif
}

void DUGKSAmr::Update_MacroVars_hDt(int const lev, double const hTime)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Integral_MacroBaseVars_Phase(ld.msv_xF, ld.hDVDFxF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnPhase>::Integral_MacroBaseVars_Phase(ld.msv_yF, ld.hDVDFyF[Ark::FLUX]);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum(ld.msv_xF, ld.fDVDFxF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum(ld.msv_yF, ld.fDVDFyF[Ark::FLUX]);
#endif

	// * * * * * * * * * * * *Interface Elongation * * * * * * * * * * * *
	// if(nstep == end_step/2)
	// {
	// 	ElongationVelocityField (ld.msv_xF[Ark::VEL], Geom(lev).CellSize(), 
	// 						  	 Geom(lev).ProbDomain(), Ark::X);
	// 	ElongationVelocityField (ld.msv_yF[Ark::VEL], Geom(lev).CellSize(), 
	// 						  	 Geom(lev).ProbDomain(), Ark::Y);
	// }
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	//
	// * * * * * * * * * * * *Interface Smooth Deformation * * * * * * * *
	// DeformationVelocityField (ld.msv_xF[Ark::VEL], Geom(lev).CellSize(), 
	// 					   		 cur_time+hDt(lev), Ark::X);
	// DeformationVelocityField (ld.msv_yF[Ark::VEL], Geom(lev).CellSize(), 
	// 					   		 cur_time+hDt(lev), Ark::Y);
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	//
	// * * * * * * * * * * * *Interface Smooth Deformation * * * * * * * *
	// ReversedVortexVelocityField (ld.msv_xF[Ark::VEL], Geom(lev).CellSize(), 
	// 					   		 cur_time+hDt(lev), Ark::X);
	// ReversedVortexVelocityField (ld.msv_yF[Ark::VEL], Geom(lev).CellSize(), 
	// 					   		 cur_time+hDt(lev), Ark::Y);
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	AllenCahnSource_hDt(ld.msv_xF, ld.msv, IntVect::TheDimensionVector(Ark::X));
	AllenCahnSource_hDt(ld.msv_yF, ld.msv, IntVect::TheDimensionVector(Ark::Y));

	Update_MacroBaseVars(ld.msv_xF);
	Update_MacroBaseVars(ld.msv_yF);

// #ifdef ARK_USE_MOMENTUM
// 	Update_Tau_Momentum(ld.fDVDFxF[Ark::TAUF], ld.msv_xF[Ark::PHI], ld.msv_xF[Ark::RHO]);
// 	Update_Tau_Momentum(ld.fDVDFyF[Ark::TAUF], ld.msv_yF[Ark::PHI], ld.msv_yF[Ark::RHO]);
// #endif

	FillPhysicalBoundaryMacroVars(lev);
}

void DUGKSAmr::Update_MacroVars(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Integral_MacroBaseVars_Phase(ld.msv, ld.hDVDF[Ark::WAVE]);
	// ReinitializePhase(ld.msv[Ark::PHI]);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum(ld.msv, ld.fDVDF[Ark::WAVE]);
#endif

	// * * * * * * * * * * * *Interface Elongation * * * * * * * * * * * *
	// if(nstep == end_step/2)
	// {
	// 	ElongationVelocityField (ld.msv[Ark::VEL], Geom(lev).CellSize(), 
	// 						  	 Geom(lev).ProbDomain(), -1);
	// }

	// * * * * * * * * * * * *Interface Smooth Deformation * * * * * * * *
	// DeformationVelocityField (ld.msv[Ark::VEL], Geom(lev).CellSize(), cur_time+Dt(), -1);

	// * * * * * * * * * * * *Reversed Single Vortex * * * * * * * *
	// ReversedVortexVelocityField(ld.msv[Ark::VEL], Geom(lev).CellSize(), cur_time+Dt(), -1);
	
	FillPatchPhase(lev, cur_time, ld.msv[Ark::PHI]);

	AllenCahnSource(lev);

	Update_MacroBaseVars(ld.msv);

// #ifdef ARK_USE_MOMENTUM
// 	Update_Tau_Momentum(ld.fDVDF[Ark::TAU], ld.msv);
// #endif

	AllenCahnTemporalPhiU(ld);
}