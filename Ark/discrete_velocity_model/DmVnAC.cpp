#include <DUGKSAmr.H>

const std::string DUGKSAmr::DmVnModel = "DmVnAC";

template <typename DVM>
template <long unsigned NUM>
void DUGKSCore<DVM>::Update_Eq_Phase (MultiFab &xEq, MultiFab &xSo, const Array<MultiFab, NUM> &msv)
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

		ParallelFor(tbx, DVM::Q, [=] AMREX_GPU_DEVICE (int i,int j,int k,int n) noexcept
		{
			const double phi  = phi_array(i,j,k);
			const double u    = vel_array(i,j,k, Ark::X);
			const double v    = vel_array(i,j,k, Ark::Y);
			const double phix = gradphi_array(i,j,k, Ark::X);
			const double phiy = gradphi_array(i,j,k, Ark::Y);
#if(AMREX_SPACEDIM == 3)
			const double w    = vel_array(i,j,k, Ark::Z);
			const double phiz = gradphi_array(i,j,k, Ark::Z);
			const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v + DVM::xi_w[n]*w)/DVM::RT;
			const double so = DVM::xi_u[n]*phix + DVM::xi_v[n]*phiy + DVM::xi_w[n]*phiz;
#elif(AMREX_SPACEDIM == 2)
			const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v)/DVM::RT;
			const double so = DVM::xi_u[n]*phix + DVM::xi_v[n]*phiy;
#endif
			hEq(i,j,k,n) = DVM::omega[n]*phi*(1.0 + u1);
			hSo(i,j,k,n) = DVM::omega[n]*so;
		});
	}
}

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
		Array4<Real const> const &Pres    = msv[Ark::P].const_array(mfi);
		Array4<Real const> const &GradDen = msv[Ark::GradRHO].const_array(mfi);
		Array4<Real const> const &Force   = msv[Ark::FORCE].const_array(mfi);

		Array4<Real> const &fEq = xEq.array(mfi);
		Array4<Real> const &fSo = xSo.array(mfi);

		ParallelFor(tbx, DVM::Q, [=] AMREX_GPU_DEVICE (int i,int j,int k,int n) noexcept
		{
       		const double rho  = Den(i,j,k);
       		const double p    = Pres(i,j,k);
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
			const double uu   = Ark::Square(u,v,w);
			const double u1   = (DVM::xi_u[n]*u + DVM::xi_v[n]*v + DVM::xi_w[n]*w)/DVM::RT;
#else
			const double uu = Ark::Square(u,v);
			const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v)/DVM::RT;
#endif
			const double Gamma = u1 + 0.5*u1*u1 - uu*DVM::Lambda;

			fEq(i,j,k,n) = DVM::omega[n]*(p/DVM::RT + rho*Gamma);

			#ifdef _ARK_M_LSGC_FLIP
			fEq(i,j,k,n) -= static_cast<int>(n == 0)*p;
			#endif

			fSo(i,j,k,n) = DVM::omega[n]/DVM::RT*
			(
			  (
			  	(DVM::xi_u[n] - u)*fx 
			  + (DVM::xi_v[n] - v)*fy
#if(AMREX_SPACEDIM == 3) 
			  + (DVM::xi_w[n] - w)*fz
#endif
			  )
			+ (
				DVM::xi_u[n]*(fx + rhox*DVM::RT) 
			  + DVM::xi_v[n]*(fy + rhoy*DVM::RT)
#if(AMREX_SPACEDIM == 3) 
			  + DVM::xi_w[n]*(fz + rhoz*DVM::RT)
#endif
			  )*u1
			);
		});
	}	
}

template <typename DVM>
void DUGKSCore<DVM>::Integral_MacroBaseVars_Phase (Array<MultiFab, Ark::MSVNUM> &msv, 
										   		   const MultiFab &xDVDF)
{
	msv[Ark::PHI].setVal(0.0);

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(msv[Ark::PHI], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real> const &phi = msv[Ark::PHI].array(mfi);

		Array4<Real const> const &hdvdf = xDVDF.const_array(mfi);

		amrex::For(bx, DVM::Q, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
		{
			phi(i,j,k) += hdvdf(i,j,k,n);
		});
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Integral_MacroBaseVars_Momentum (Array<MultiFab, Ark::MSVNUM> &msv, 
										   		   	  const MultiFab &xDVDF)
{
	msv[Ark::RHO].setVal(0.0);
	msv[Ark::VEL].setVal(0.0);
	msv[Ark::P].setVal(0.0);
	msv[Ark::T].setVal(0.0);

	#ifdef _OPENMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif
	for(amrex::MFIter mfi(msv[Ark::RHO], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real> const &vel  = msv[Ark::VEL].array(mfi);
		Array4<Real> const &pres = msv[Ark::P].array(mfi);

		Array4<Real const> const &fdvdf = xDVDF.const_array(mfi);

		amrex::For(bx, DVM::Q, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
		{
			pres(i,j,k) 	  += fdvdf(i,j,k,n);
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
	amrex::Abort("DimVnACSS model is supposed to be used");
}

void DUGKSAmr::Update_MacroSourceVars (Array<MultiFab, Ark::MSVNUM> &msv, const Real hTime)
{
	// for(amrex::MFIter mfi(msv[Ark::VEL], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	// {
	// 	const amrex::Box& bx = mfi.validbox();

	// 	FArrayBox& VelFab = msv[Ark::VEL][mfi];

	// 	VelFab.setVal(0,bx,Ark::X);
	// 	VelFab.setVal(U0(),bx,Ark::Y);
	// }
		
#ifdef ARK_USE_MOMENTUM
	const MultiFab &PhiFab     = msv[Ark::PHI];
	const MultiFab &GradRhoFab = msv[Ark::GradRHO];

	//density
	MultiFab &RhoFab  = msv[Ark::RHO];
	MultiFab &VelFab  = msv[Ark::VEL];
	MultiFab &PresFab = msv[Ark::P];
	MultiFab &ForceFab   = msv[Ark::FORCE];

	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &phi 	   = PhiFab.const_array(mfi);
		Array4<Real const> const &grad_rho = GradRhoFab.const_array(mfi);

		Array4<Real> const &rho   = RhoFab.array(mfi);
		Array4<Real> const &vel   = VelFab.array(mfi);
		Array4<Real> const &pres  = PresFab.array(mfi);
		Array4<Real> const &force    = ForceFab.array(mfi);

		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
		{
			rho(i,j,k) = APhi()*phi(i,j,k) + BPhi();
		});

		// Do NOT move this part elsewhere
		// Both CELL values and FACE values should perform this operation
		amrex::ParallelFor(bx, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i,int j,int k, int n) noexcept
		{
			// force(i,j,k,n) += (rho(i,j,k) - RhoL())*Gravity(n);
			force(i,j,k,n) += (rho(i,j,k) - RhoV())*Gravity(n);
			vel(i,j,k,n) += hTime*force(i,j,k,n);
			vel(i,j,k,n) /= rho(i,j,k);
		});

		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
		{
			pres(i,j,k) +=  hTime*
						   	(
								vel(i,j,k,Ark::X)*grad_rho(i,j,k,Ark::X)
							+	vel(i,j,k,Ark::Y)*grad_rho(i,j,k,Ark::Y)
#if(AMREX_SPACEDIM == 3)
							+	vel(i,j,k,Ark::Z)*grad_rho(i,j,k,Ark::Z)
#endif
							);

			#ifdef _ARK_M_LSGC_FLIP
			msq(i,j,k,Ark::P) -= 2.0/9.0*msq(i,j,k,Ark::RHO)*msq(i,j,k,Ark::SQUV);
			msq(i,j,k,Ark::P) *= 9.0/5.0;
			#endif
		});
	}
	PresFab.mult(Ark::DmVnMomentum::RT, 0);
#endif
}

void DUGKSAmr::Update_MacroVars(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Integral_MacroBaseVars_Phase(ld.msv, ld.hDVDF[Ark::WAVE]);
	if(lev != finest_level) 
	{
		BoundPhase(ld.msv[Ark::PHI]);
	}
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum(ld.msv, ld.fDVDF[Ark::WAVE]);
#endif

	// * * * * * * * * * * * *Interface Elongation * * * * * * * * * * * *
	// if(nstep == end_step/2)
	{
// #if(AMREX_SPACEDIM == 2)
// 		ElongationVelocityField (ld.msv[Ark::VEL], Geom(lev).CellSize(), Geom(lev).ProbDomain(), -1);
// #elif(AMREX_SPACEDIM == 3)
// 		ElongationVelocityField3D (ld.msv[Ark::VEL], Geom(lev).CellSize(), Geom(lev).ProbDomain(), -1);
// #endif
	}

	// * * * * * * * * * * * *Interface Smooth Deformation * * * * * * * *
	// DeformationVelocityField (ld.msv[Ark::VEL], Geom(lev).CellSize(), cur_time+Dt(), -1);

	// * * * * * * * * * * * *Reversed Single Vortex * * * * * * * *
	// ReversedVortexVelocityField(ld.msv[Ark::VEL], Geom(lev).CellSize(), cur_time+Dt(), -1);
	
	FillPatchPhase(lev, cur_time, ld.msv[Ark::PHI]);

	AllenCahnSource(lev);

	Update_MacroSourceVars(ld.msv, 0.5*Dt());

	AllenCahnTemporalPhiU(ld);
}
void DUGKSAmr::Update_MacroVars_hDt(int const lev, double const hTime)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Integral_MacroBaseVars_Phase(ld.msv_xF, ld.hDVDFxF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnPhase>::Integral_MacroBaseVars_Phase(ld.msv_yF, ld.hDVDFyF[Ark::FLUX]);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Integral_MacroBaseVars_Phase(ld.msv_zF, ld.hDVDFzF[Ark::FLUX]);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum(ld.msv_xF, ld.fDVDFxF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum(ld.msv_yF, ld.fDVDFyF[Ark::FLUX]);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Integral_MacroBaseVars_Momentum(ld.msv_zF, ld.fDVDFzF[Ark::FLUX]);
#endif
#endif

	// * * * * * * * * * * * *Interface Elongation * * * * * * * * * * * *
	// if(nstep == end_step/2) 
	{
// #if(AMREX_SPACEDIM == 2)
// 	ElongationVelocityField (ld.msv_xF[Ark::VEL], Geom(lev).CellSize(), Geom(lev).ProbDomain(), Ark::X);
// 	ElongationVelocityField (ld.msv_yF[Ark::VEL], Geom(lev).CellSize(), Geom(lev).ProbDomain(), Ark::Y);
// #elif(AMREX_SPACEDIM == 3)
// 	ElongationVelocityField3D (ld.msv_xF[Ark::VEL], Geom(lev).CellSize(), Geom(lev).ProbDomain(), Ark::X);
// 	ElongationVelocityField3D (ld.msv_yF[Ark::VEL], Geom(lev).CellSize(), Geom(lev).ProbDomain(), Ark::Y);
// 	ElongationVelocityField3D (ld.msv_zF[Ark::VEL], Geom(lev).CellSize(), Geom(lev).ProbDomain(), Ark::Z);
// #endif
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	//
	// * * * * * * * * * * * *Interface Smooth Deformation * * * * * * * *
	// DeformationVelocityField (ld.msv_xF[Ark::VEL], Geom(lev).CellSize(),cur_time+hDt(lev), Ark::X);
	// DeformationVelocityField (ld.msv_yF[Ark::VEL], Geom(lev).CellSize(),cur_time+hDt(lev), Ark::Y);
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	//
	// * * * * * * * * * * * *Interface Smooth Deformation * * * * * * * *
	// ReversedVortexVelocityField(ld.msv_xF[Ark::VEL], Geom(lev).CellSize(),cur_time+hDt(lev), Ark::X);
	// ReversedVortexVelocityField(ld.msv_yF[Ark::VEL], Geom(lev).CellSize(),cur_time+hDt(lev), Ark::Y);
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	AllenCahnSource_hDt(ld.msv_xF, ld.msv, IntVect::TheDimensionVector(Ark::X));
	AllenCahnSource_hDt(ld.msv_yF, ld.msv, IntVect::TheDimensionVector(Ark::Y));
#if(AMREX_SPACEDIM == 3)
	AllenCahnSource_hDt(ld.msv_zF, ld.msv, IntVect::TheDimensionVector(Ark::Z));
#endif

	Update_MacroSourceVars(ld.msv_xF, hTime);
	Update_MacroSourceVars(ld.msv_yF, hTime);
#if(AMREX_SPACEDIM == 3)
	Update_MacroSourceVars(ld.msv_zF, hTime);
#endif

// #ifdef ARK_USE_MOMENTUM
// 	Update_Tau_Momentum(ld.fDVDFxF[Ark::TAUF], ld.msv_xF);
// 	Update_Tau_Momentum(ld.fDVDFyF[Ark::TAUF], ld.msv_yF);
// #if(AMREX_SPACEDIM == 3)
// 	Update_Tau_Momentum(ld.fDVDFzF[Ark::TAUF], ld.msv_zF);
// #endif
// #endif

	FillPhysicalBoundaryMacroVars(lev);
}
void DUGKSAmr::Update_Tau_Phase(MultiFab &tauFab, Array<MultiFab, Ark::MSVNUM> &msv)
{
	tauFab.setVal(TauPhase());
}
void DUGKSAmr::Update_Tau_Momentum (MultiFab &tauFab, Array<MultiFab, Ark::MSVNUM> &msv)
{
	const MultiFab &PhiFab = msv[Ark::PHI];
	const MultiFab &RhoFab = msv[Ark::RHO];
	
	// const Real AMu = (MuL() - MuV())/(PhiL() - PhiV());
	// const Real BMu = MuV() - AMu*PhiV();

	// for(amrex::MFIter mfi(tauFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	// {
	// 	const amrex::Box& bx = mfi.tilebox();

	// 	Array4<Real const> const &phi = PhiFab.const_array(mfi);
	// 	// Array4<Real const> const &rho = RhoFab.const_array(mfi);

	// 	Array4<Real> const &tau = tauFab.array(mfi);

	// 	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
	// 	{
	// 		//linear
	// 		// tau(i,j,k) = AMu*phi(i,j,k) + BMu;

	// 		// -reciprocal linear
	// 		tau(i,j,k) = MuL()*MuV()*(PhiL() - PhiV())
	// 		/((phi(i,j,k)-PhiV())*MuV() - (phi(i,j,k)-PhiL())*MuL());
			
	// 	});
	// }

	// tauFab.mult(2*Ark::DmVnMomentum::Lambda,0);
	// MultiFab::Divide(tauFab, RhoFab, 0, 0, 1, 0);



	const Real ANu = (NuL() - NuV())/(PhiL() - PhiV());
	const Real BNu = NuV() - ANu*PhiV();

	for(amrex::MFIter mfi(tauFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &phi = PhiFab.const_array(mfi);
		// Array4<Real const> const &rho = RhoFab.const_array(mfi);

		Array4<Real> const &tau = tauFab.array(mfi);

		amrex::ParallelFor(bx,
		[=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
		{
			// linear
			tau(i,j,k) = ANu*phi(i,j,k) + BNu;

			//reciprocal linear
			// tau(i,j,k) = NuL()*NuV()*(PhiL() - PhiV())
			// /((phi(i,j,k)-PhiV())*NuV() - (phi(i,j,k)-PhiL())*NuL());
			
		});
	}

	tauFab.mult(2*Ark::DmVnMomentum::Lambda, 0);

	// tauFab.setVal(TauMomentum());
}

#include <DmVnEq.H>