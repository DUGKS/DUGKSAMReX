#include <DUGKSAmr.H>

#if(AMREX_SPACEDIM != 2)
#error "AMREX_SPACEDIM is expected to be 2!!!"
#endif

#if ARK_USE_PHASE
#error "wrong velocity model!!!"
#endif

template <typename DVM>
template <long unsigned NUM>
void DUGKSCore<DVM>::Update_Eq_Phase (MultiFab &xEq, MultiFab &xSo,
									  const Array<MultiFab, NUM> &msv)
{
	
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
		Array4<Real const> const &Force   = msv[Ark::FORCE].const_array(mfi);

		Array4<Real> const &fEq = xEq.array(mfi);
		Array4<Real> const &fSo = xSo.array(mfi);

		ParallelFor(tbx, DVM::Q,
		[=] AMREX_GPU_DEVICE (int i,int j,int k,int n) noexcept
		{
       		const double rho  = Den(i,j,k);
       		// const double p    = Pres(i,j,k);
			const double u    = Vel(i,j,k,Ark::X);
     		const double v    = Vel(i,j,k,Ark::Y);
       		const double fx   = Force(i,j,k,Ark::X);
       		const double fy   = Force(i,j,k,Ark::Y);

			const double uu = Ark::Square(u,v);

			const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v)/DVM::RT;
			const double Gamma = u1 + 0.5*u1*u1 - uu*DVM::Lambda;

			fEq(i,j,k,n) = DVM::omega[n]*(rho + Gamma);
			fSo(i,j,k,n) = ( (DVM::xi_u[n]-u)*fx + (DVM::xi_v[n]-v)*fy )*fEq(i,j,k,n)/DVM::RT;
		});
	}	
}

template <typename DVM>
void DUGKSCore<DVM>::Update_MacroBaseVars_Phase(Array<MultiFab, Ark::MSVNUM> &msv, 
										   		const MultiFab &xDVDF)
{

}

template <typename DVM>
void DUGKSCore<DVM>::Update_MacroBaseVars_Momentum(Array<MultiFab, Ark::MSVNUM> &msv, 
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

void DUGKSAmr::Update_MacroSourceVars (Array<MultiFab, Ark::MSVNUM> &msv, const Real hTime)
{
#ifdef ARK_USE_MOMENTUM
	//density
	MultiFab &RhoFab  = msv[Ark::RHO];
	MultiFab &VelFab  = msv[Ark::VEL];
	MultiFab &PresFab = msv[Ark::P];
	MultiFab &ForceFab   = msv[Ark::FORCE];
	
	// MultiFab::Divide(VelFab, RhoFab, 0, 0, 1, 0);
	// MultiFab::Divide(VelFab, RhoFab, 0, 1, 1, 0);
	MultiFab::Saxpy(VelFab, 0.5*hTime, ForceFab, 0, 0, AMREX_SPACEDIM, 0);

	MultiFab::Copy(PresFab, RhoFab, 0, 0, 1, zero_ghost());
	PresFab.plus(-1.0, 0);
	PresFab.mult(Ark::DmVnMomentum::RT, 0);
#endif
}
void DUGKSAmr::Update_MacroVars(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv, ld.fDVDF[Ark::WAVE]);
#endif

	Update_MacroSourceVars(ld.msv, Dt());
}
void DUGKSAmr::Update_MacroVars_hDt(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv_xF, ld.fDVDFxF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv_yF, ld.fDVDFyF[Ark::FLUX]);
#endif

	// AllenCahnSource_hDt(ld.msv_xF, ld.msv, IntVect::TheDimensionVector(Ark::X));
	// AllenCahnSource_hDt(ld.msv_yF, ld.msv, IntVect::TheDimensionVector(Ark::Y));

	Update_MacroSourceVars(ld.msv_xF, hDt(lev));
	Update_MacroSourceVars(ld.msv_yF, hDt(lev));

	FillPhysicalBoundaryMacroVars(lev);
}