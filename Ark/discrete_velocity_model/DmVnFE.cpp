#include <DUGKSAmr.H>

#if ARK_USE_PHASE
#error "wrong velocity model!!!"
#endif

const std::string DUGKSAmr::DmVnModel = "DmVnFE";

template <typename DVM>
template <bool ENABLESO>
void DUGKSCore<DVM>::Update_Eq_Phase (MultiFab &xEq, MultiFab &xSo, const Array<MultiFab, Ark::MSVNUM> &msv)
{
	amrex::Abort("Free Energy model disables PHASE distribution.");
}


//Rho_0 == 0
template <typename DVM>
template <bool ENABLESO>
void DUGKSCore<DVM>::Update_Eq_Momentum (MultiFab &xEq, MultiFab &xSo,const Array<MultiFab, Ark::MSVNUM> &msv)
{
	xSo.setVal(0);

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xEq, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const amrex::Box& tbx = mfi.tilebox();
	//
		Array4<Real const> const &Den     = msv[Ark::RHO].const_array(mfi);
		Array4<Real const> const &Vel     = msv[Ark::VEL].const_array(mfi);
		Array4<Real const> const &ChemiP  = msv[Ark::CHEMIP].const_array(mfi);
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
#if(AMREX_SPACEDIM == 3)
     		const double w    = Vel(i,j,k,Ark::Z);
       		const double rhoz = GradDen(i,j,k,Ark::Z);
       		const double fz   = Force(i,j,k,Ark::Z);
			const double uGradrho = u*rhox + v*rhoy + w*rhoz;
			const double uu   = Ark::Square(u,v,w);
			const double u1   = (DVM::xi_u[n]*u + DVM::xi_v[n]*v + DVM::xi_w[n]*w)/DVM::RT;
#else
			const double uu = Ark::Square(u,v);
			const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v)/DVM::RT;
			const double uGradrho = u*rhox + v*rhoy;
#endif
			const double Gamma = u1 + 0.5*u1*u1 - uu*DVM::Lambda;

// // * * * * * * * * * * * * * * * primitive model * * * * * * * * * * * * * * * * * * * * 
			fEq(i,j,k,n) = DVM::omega[n]*rho*(1.0 + Gamma);	
		
		if constexpr (ENABLESO) {
			fSo(i,j,k,n) = DVM::omega[n]/DVM::RT*
			(
				(DVM::xi_u[n] - u)*fx 
			  + (DVM::xi_v[n] - v)*fy
#if(AMREX_SPACEDIM == 3) 
			  + (DVM::xi_w[n] - w)*fz
#endif
			  + (
				  DVM::xi_u[n]*fx 
			    + DVM::xi_v[n]*fy
#if(AMREX_SPACEDIM == 3) 
			    + DVM::xi_w[n]*fz
#endif
			  	)*u1
			);
		}

// * * * * * * * * * * * * * * * well-balanced model * * * * * * * * * * * * * * * * * * * * 
// 		fEq(i,j,k,n) = DVM::omega[n]*rho*Gamma + static_cast<int>(n == 0)*rho;

// 		if constexpr (ENABLESO) {
// 			fSo(i,j,k,n) = DVM::omega[n]/DVM::RT*
// 			(
// 				(DVM::xi_u[n] - u)*fx 
// 			  + (DVM::xi_v[n] - v)*fy
// #if(AMREX_SPACEDIM == 3) 
// 			  + (DVM::xi_w[n] - w)*fz
// #endif
// 			  + (
// 				  DVM::xi_u[n]*(fx + rhox*DVM::RT) 
// 			    + DVM::xi_v[n]*(fy + rhoy*DVM::RT)
// #if(AMREX_SPACEDIM == 3) 
// 			    + DVM::xi_w[n]*(fz + rhoz*DVM::RT)
// #endif
// 			  	)*u1
// 			)
// 			+ DVM::omega[n]*uGradrho*0.5*(
// #if(AMREX_SPACEDIM == 3)
//                 Ark::Square(DVM::xi_u[n],DVM::xi_v[n],DVM::xi_w[n])/DVM::RT-(AMREX_SPACEDIM)-2
// #elif(AMREX_SPACEDIM == 2)
//                 Ark::Square(DVM::xi_u[n],DVM::xi_v[n])/DVM::RT-(AMREX_SPACEDIM)-2
// #endif 
// 			);
// 		}

// * * * * * * * * * * * * * * * Lee2006PRE Deprecated* * * * * * * * * * * * * * * * * * * *
// 		fEq(i,j,k,n) = DVM::omega[n]*rho*(1.0 + Gamma);	

// 		if constexpr (ENABLESO) {
// 			const int ni  = static_cast<int>(DVM::xi_u[n]);
// 			const int nj  = static_cast<int>(DVM::xi_v[n]);
// 			const int nii = ni + ni; 
// 			const int njj = nj + nj;
// #if(AMREX_SPACEDIM == 3) 			
// 			const int nk  = static_cast<int>(DVM::xi_w[n]);
// 			const int nkk = nk + nk;
// #else
// 			const int nk  = 0;
// 			const int nkk = 0;
// #endif

// 			const double biasedDirectionalDensity 
// 			= (-Den(i+nii,j+njj,k+nkk) + 4*Den(i+ni,j+nj,k+nk) - 3*rho)*0.5;

// 			const double biasedDirectionalChemiP  
// 			= (-ChemiP(i+nii,j+njj,k+nkk) + 4*ChemiP(i+ni,j+nj,k+nk) - 3*ChemiP(i,j,k))*0.5;

// 			const double centralDirectionalDensity = (Den(i+ni,j+nj,k+nk) - Den(i-ni,j-nj,k-nk))*0.5;
// 			const double centralDirectionalChemiP  = (ChemiP(i+ni,j+nj,k+nk) - ChemiP(i-ni,j-nj,k-nk))*0.5;

// 			const double directionalDensity = (biasedDirectionalDensity + centralDirectionalDensity)*0.5;
// 			const double directionalChemiP  = (biasedDirectionalChemiP + centralDirectionalChemiP)*0.5;

// 			// const double directionalDensity = centralDirectionalDensity;
// 			// const double directionalChemiP  = centralDirectionalChemiP;

// 			const double grad_rho_x = (Den(i+1,j,k) - Den(i-1,j,k))*0.5;
// 			const double grad_rho_y = (Den(i,j+1,k) - Den(i,j-1,k))*0.5;

// 			const double grad_chemip_x = (ChemiP(i+1,j,k) - ChemiP(i-1,j,k))*0.5;
// 			const double grad_chemip_y = (ChemiP(i,j+1,k) - ChemiP(i,j-1,k))*0.5;

// 			const double fxNow = (DVM::RT*grad_rho_x - rho*grad_chemip_x);
// 			const double fyNow = (DVM::RT*grad_rho_y - rho*grad_chemip_y);  								 

// 			fSo(i,j,k,n) = fEq(i,j,k,n)*
// 			(
// 				// DVM::RT*directionalDensity - rho*directionalChemiP -
// 				(
// 				  (DVM::xi_u[n])*fx 
// 			    + (DVM::xi_v[n])*fy
// #if(AMREX_SPACEDIM == 3) 
// 			    + (DVM::xi_w[n])*fz
// #endif
// 			  	) -
				
// 			   (
// 				  (u)*fx 
// 			    + (v)*fy
// #if(AMREX_SPACEDIM == 3) 
// 			    + (w)*fz
// #endif
// 			  	)
// 			)/(DVM::RT*rho);
// 		}
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
#if(AMREX_SPACEDIM == 3) 
			vel(i,j,k,Ark::Z) += fdvdf(i,j,k,n)*DVM::xi_w[n];
#endif
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

	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum(ld.msv, ld.fDVDF[Ark::WAVE]);

	FillPatchDensity(lev, cur_time, ld.msv[Ark::RHO]);

	FreeEnergySource(lev);
	// PseudoPotential(lev);

	Update_MacroBaseVars(ld.msv);
}
void DUGKSAmr::Update_MacroVars_hDt(int const lev, const double hTime)
{
	auto &ld = Data(lev);

	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum (ld.msv_xF, ld.fDVDFxF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum (ld.msv_yF, ld.fDVDFyF[Ark::FLUX]);
#if(AMREX_SPACEDIM == 3) 
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum (ld.msv_zF, ld.fDVDFzF[Ark::FLUX]);
#endif

	Update_MacroBaseVars(ld.msv_xF);
	Update_MacroBaseVars(ld.msv_yF);
#if(AMREX_SPACEDIM == 3) 
	Update_MacroBaseVars(ld.msv_zF);
#endif

	FillPhysicalBoundaryMacroVars(lev);
}

void DUGKSAmr::Update_Tau_Momentum (MultiFab &tauFab, Array<MultiFab, Ark::MSVNUM> &msv)
{
	tauFab.setVal(TauMomentum());
}

#include <DUGKSDmVnEq.H>