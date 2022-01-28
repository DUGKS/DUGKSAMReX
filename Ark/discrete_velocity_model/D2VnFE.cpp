#include <DUGKSAmr.H>

#if(AMREX_SPACEDIM != 2)
#error "AMREX_SPACEDIM is expected to be 2!!!"
#endif

#if ARK_USE_PHASE
#error "wrong velocity model!!!"
#endif

template <typename DVM>
template <long unsigned NUM>
void DUGKSCore<DVM>::Update_Eq_Phase (MultiFab &xEq, MultiFab &xSo, const Array<MultiFab, NUM> &msv)
{
	
}


//Rho_0 == 0
template <typename DVM>
template <long unsigned NUM>
void DUGKSCore<DVM>::Update_Eq_Momentum (MultiFab &xEq, MultiFab &xSo, const Array<MultiFab, NUM> &msv)
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
		Array4<Real const> const &GradDen = msv[Ark::GradRHO].const_array(mfi);
		Array4<Real const> const &Force   = msv[Ark::FORCE].const_array(mfi);

		Array4<Real> const &fEq = xEq.array(mfi);
		Array4<Real> const &fSo = xSo.array(mfi);

		ParallelFor(tbx, DVM::Q, [=] AMREX_GPU_DEVICE (int i,int j,int k,int n) noexcept
		{
       		const double rho  = Den(i,j,k);
       		// const double p    = Pres(i,j,k);
			const double u    = Vel(i,j,k,Ark::X);
     		const double v    = Vel(i,j,k,Ark::Y);
     		const double rhox = GradDen(i,j,k,Ark::X);
       		const double rhoy = GradDen(i,j,k,Ark::Y);
       		const double fx   = Force(i,j,k,Ark::X);
       		const double fy   = Force(i,j,k,Ark::Y);

			const double uu = Ark::Square(u,v);

			const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v)/DVM::RT;
			const double Gamma = u1 + 0.5*u1*u1 - uu*DVM::Lambda;

			fEq(i,j,k,n) = DVM::omega[n]*rho*Gamma + static_cast<int>(n == 0)*rho;

			fSo(i,j,k,n) = DVM::omega[n]*
			(
				(
				((DVM::xi_u[n] - u)*fx + (DVM::xi_v[n] - v)*fy)
				+ (DVM::xi_u[n]*(fx + rhox*DVM::RT) + DVM::xi_v[n]*(fy + rhoy*DVM::RT))*u1
				)/DVM::RT
			+   ((DVM::xi_u[n]*DVM::xi_u[n] + DVM::xi_v[n]*DVM::xi_v[n])*DVM::Lambda -2)
				*(u*rhox + v*rhoy)
			);
		});
	}	
}

template <typename DVM>
void DUGKSCore<DVM>::Integral_MacroBaseVars_Phase(Array<MultiFab, Ark::MSVNUM> &msv, 
										   		const MultiFab &xDVDF)
{

}

template <typename DVM>
void DUGKSCore<DVM>::Integral_MacroBaseVars_Momentum(Array<MultiFab, Ark::MSVNUM> &msv, 
										   		   const MultiFab &xDVDF)
{
	msv[Ark::RHO].setVal(0.0);
	msv[Ark::VEL].setVal(0.0);
	msv[Ark::T].setVal(0.0);

	#ifdef _OPENMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif
	for(amrex::MFIter mfi(msv[Ark::RHO], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real> const &rho = msv[Ark::RHO].array(mfi);
		Array4<Real> const &vel = msv[Ark::VEL].array(mfi);

		Array4<Real const> const &fdvdf = xDVDF.const_array(mfi);

		amrex::For(bx, DVM::Q,
		[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
		{
			rho(i,j,k) 	  += fdvdf(i,j,k,n);
			vel(i,j,k,Ark::X) += fdvdf(i,j,k,n)*DVM::xi_u[n];
			vel(i,j,k,Ark::Y) += fdvdf(i,j,k,n)*DVM::xi_v[n];
		});
	}
}
void DUGKSAmr::Update_MacroBaseVars (Array<MultiFab, Ark::MSVNUM> &msv)
{
#ifdef ARK_USE_MOMENTUM
	const MultiFab &RhoFab   = msv[Ark::RHO];

	MultiFab &VelFab   = msv[Ark::VEL];

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
	MultiFab const &ForceFab = msv[Ark::FORCE];
	MultiFab const &RhoFab   = msv[Ark::RHO];

	MultiFab &VelFab  = msv[Ark::VEL];

	for(amrex::MFIter mfi(ForceFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &rho 	   = RhoFab.array(mfi);
		Array4<Real const> const &force    = ForceFab.array(mfi);

		Array4<Real> const &vel   = VelFab.array(mfi);

		amrex::ParallelFor(bx, AMREX_SPACEDIM,[=] AMREX_GPU_DEVICE (int i,int j,int k, int n) noexcept
		{
			vel(i,j,k,n) += interval*force(i,j,k,n)/rho(i,j,k);
		});	
	}
#endif
}

void DUGKSAmr::Update_MacroVars(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum(ld.msv, ld.fDVDF[Ark::WAVE]);
#endif

	FillPatchDensity(lev, cur_time, ld.msv[Ark::RHO]);

	FreeEnergySource(lev);

	Update_MacroBaseVars(ld.msv);
}
void DUGKSAmr::Update_MacroVars_hDt(int const lev, const double hTime)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum (ld.msv_xF, ld.fDVDFxF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum (ld.msv_yF, ld.fDVDFyF[Ark::FLUX]);
#endif

	Update_MacroBaseVars(ld.msv_xF);
	Update_MacroBaseVars(ld.msv_yF);

	FillPhysicalBoundaryMacroVars(lev);
}

#include "DmVnEq.H"