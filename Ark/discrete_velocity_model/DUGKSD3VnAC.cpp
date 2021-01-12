#include <DUGKSAmr.H>

#if(AMREX_SPACEDIM != 3)
#error "AMREX_SPACEDIM is expected to be 3!!!"
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
			const Real w    = vel_array(i,j,k, Ark::Z);
			const Real phix = gradphi_array(i,j,k, Ark::X);
			const Real phiy = gradphi_array(i,j,k, Ark::Y);
			const Real phiz = gradphi_array(i,j,k, Ark::Z);

			const Real u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v + DVM::xi_w[n]*w)/DVM::RT;

			hEq(i,j,k,n) = DVM::omega[n]*phi*(1.0 + u1);
			hSo(i,j,k,n) = DVM::omega[n]*(DVM::xi_u[n]*phix + DVM::xi_v[n]*phiy
															+ DVM::xi_w[n]*phiz);
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
       		const double &rho  = Den(i,j,k);
       		const double &p    = Pres(i,j,k);
			const double &u    = Vel(i,j,k,Ark::X);
     		const double &v    = Vel(i,j,k,Ark::Y);
     		const double &w    = Vel(i,j,k,Ark::Z);
       		const double &rhox = GradDen(i,j,k,Ark::X);
       		const double &rhoy = GradDen(i,j,k,Ark::Y);
       		const double &rhoz = GradDen(i,j,k,Ark::Z);
       		const double &fx   = Force(i,j,k,Ark::X);
       		const double &fy   = Force(i,j,k,Ark::Y);
       		const double &fz   = Force(i,j,k,Ark::Z);

			const double uu = Ark::Square(u,v,w);

			const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v + DVM::xi_w[n]*w)/DVM::RT;
			const double Gamma = u1 + 0.5*u1*u1 - uu*DVM::Lambda;

			fEq(i,j,k,n) = DVM::omega[n]*(p/DVM::RT + rho*Gamma);

			#ifdef _ARK_M_LSGC_FLIP
			fEq(i,j,k,n) -= static_cast<int>(n == 0)*p;
			#endif

			// fSo(i,j,k,n) = DVM::omega[n]*
			// 		   (
			// 		   (1+Gamma)*((DVM::xi_u[n] - u)*fx + (DVM::xi_v[n] - v)*fy 
			// 		   									+ (DVM::xi_w[n] - w)*fz)/DVM::RT
			// 		   +
			// 		   Gamma*((DVM::xi_u[n] - u)*rhox + (DVM::xi_v[n] - v)*rhoy
			// 		   								  + (DVM::xi_w[n] - w)*rhoz)
			// 		   );
			fSo(i,j,k,n) = DVM::omega[n]*
			(
			  ((DVM::xi_u[n] - u)*fx + (DVM::xi_v[n] - v)*fy + (DVM::xi_w[n] - w)*fz)
			+ (
				DVM::xi_u[n]*(fx + rhox*DVM::RT) 
			  + DVM::xi_v[n]*(fy + rhoy*DVM::RT)
			  + DVM::xi_w[n]*(fz + rhoz*DVM::RT)
			  )*u1
			)/DVM::RT;
		});
	}	
}

template <typename DVM>
void DUGKSCore<DVM>::Update_MacroBaseVars_Phase(Array<MultiFab, Ark::MSVNUM> &msv, 
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

		amrex::For(bx, DVM::Q,
		[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
		{
			phi(i,j,k) += hdvdf(i,j,k,n);
		});
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Update_MacroBaseVars_Momentum(Array<MultiFab, Ark::MSVNUM> &msv, 
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

		amrex::For(bx, DVM::Q,
		[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
		{
			pres(i,j,k) 	  += fdvdf(i,j,k,n);
			vel(i,j,k,Ark::X) += fdvdf(i,j,k,n)*DVM::xi_u[n];
			vel(i,j,k,Ark::Y) += fdvdf(i,j,k,n)*DVM::xi_v[n];
			vel(i,j,k,Ark::Z) += fdvdf(i,j,k,n)*DVM::xi_w[n];
		});
	}
}
void DUGKSAmr::Update_MacroSourceVars (Array<MultiFab, Ark::MSVNUM> &msv, const Real hTime)
{
#ifdef ARK_USE_MOMENTUM
	const MultiFab &PhiFab     = msv[Ark::PHI];
	const MultiFab &GradRhoFab = msv[Ark::GradRHO];

	//density
	MultiFab &RhoFab  = msv[Ark::RHO];
	MultiFab &VelFab  = msv[Ark::VEL];
	MultiFab &PresFab = msv[Ark::P];
	MultiFab &ForceFab   = msv[Ark::FORCE];

	MultiFab::Copy(RhoFab, PhiFab, 0, 0, 1, zero_ghost());
	RhoFab.mult(APhi(), 0);
	RhoFab.plus(BPhi(), 0);

	// MultiFab::Saxpy(VelFab, 0.5*hTime, ForceFab, 0, 0, AMREX_SPACEDIM, 0);
	// MultiFab::Divide(VelFab, RhoFab, 0, 0, 1, 0);
	// MultiFab::Divide(VelFab, RhoFab, 0, 1, 1, 0);
	// MultiFab::Divide(VelFab, RhoFab, 0, 2, 1, 0);


	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &rho 	   = RhoFab.array(mfi);
		Array4<Real const> const &grad_rho = GradRhoFab.const_array(mfi);

		Array4<Real> const &vel   = VelFab.array(mfi);
		Array4<Real> const &pres  = PresFab.array(mfi);
		Array4<Real> const &force = ForceFab.array(mfi);

		amrex::ParallelFor(bx, AMREX_SPACEDIM,
		[=] AMREX_GPU_DEVICE (int i,int j,int k, int n) noexcept
		{
			force(i,j,k,n) += rho(i,j,k)*Gravity(n);
			vel(i,j,k,n) += 0.5*hTime*(force(i,j,k,n));
			vel(i,j,k,n) /= rho(i,j,k);
		});

		amrex::ParallelFor(bx,
		[=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
		{
			pres(i,j,k) +=  0.5*hTime*
						   	(
								vel(i,j,k,Ark::X)*grad_rho(i,j,k,Ark::X)
							+	vel(i,j,k,Ark::Y)*grad_rho(i,j,k,Ark::Y)
							+	vel(i,j,k,Ark::Z)*grad_rho(i,j,k,Ark::Z)
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
void DUGKSAmr::Update_Tau_Momentum (MultiFab &tauFab, MultiFab const &PhiFab, 
													  MultiFab const &RhoFab)
{
	// const Real AMu = (MuL() - MuV())/(PhiL() - PhiV());
	// const Real BMu = MuV() - AMu*PhiV();

	// MultiFab::Copy(tauFab, PhiFab, 0, 0, 1, zero_ghost());
	// tauFab.mult(AMu, 0);
	// tauFab.plus(BMu, 0);

	// tauFab.mult(2*Ark::DmVnMomentum::Lambda,0);
	// MultiFab::Divide(tauFab, RhoFab, 0, 0, 1, 0);

	//***************************
	//     Inverse Linear       *
	//***************************
	// for(amrex::MFIter mfi(tauFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	// {
	// 	const amrex::Box& bx = mfi.tilebox();

	// 	Array4<Real const> const &phi = PhiFab.const_array(mfi);
	// 	Array4<Real const> const &rho = RhoFab.const_array(mfi);

	// 	Array4<Real> const &tau = tauFab.array(mfi);

	// 	amrex::ParallelFor(bx,
	// 	[=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
	// 	{
	// 		tau(i,j,k) = MuL()*MuV()/
	// 		((phi(i,j,k)-PhiV())*MuV() - (phi(i,j,k)-PhiL())*MuL());
	// 	});
	// }
}
void DUGKSAmr::Update_DVDF_Eq(int const lev)
{	
	auto &ld = Data(lev);
#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase(ld.hDVDF[Ark::EQ], ld.hDVDF[Ark::SO], ld.msv);
#endif

#ifdef ARK_USE_MOMENTUM	
	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum(ld.fDVDF[Ark::EQ], ld.fDVDF[Ark::SO], 
													 ld.msv);
#endif
}
void DUGKSAmr::Update_DVDF_Eq_hDt(int const lev)
{	
	auto &ld = Data(lev);

	#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase(ld.hDVDFxF[Ark::EQF], ld.hDVDFxF[Ark::SOF], 
											   ld.msv_xF);

	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase(ld.hDVDFyF[Ark::EQF], ld.hDVDFyF[Ark::SOF], 
											   ld.msv_yF);

	DUGKSCore<Ark::DmVnPhase>::Update_Eq_Phase(ld.hDVDFzF[Ark::EQF], ld.hDVDFzF[Ark::SOF], 
											   ld.msv_zF);
	#endif

	#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum(ld.fDVDFxF[Ark::EQF], ld.fDVDFxF[Ark::SOF], 
												     ld.msv_xF);

	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum(ld.fDVDFyF[Ark::EQF], ld.fDVDFyF[Ark::SOF], 
												     ld.msv_yF);

	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum(ld.fDVDFzF[Ark::EQF], ld.fDVDFzF[Ark::SOF], 
												     ld.msv_zF);
	#endif
}
void DUGKSAmr::Update_MacroVars(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_MacroBaseVars_Phase(ld.msv, ld.hDVDF[Ark::WAVE]);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv, ld.fDVDF[Ark::WAVE]);
#endif
	
	FillPatchPhase(lev, cur_time, ld.msv[Ark::PHI]);

	AllenCahnSource(lev);

	Update_MacroSourceVars(ld.msv, Dt());

#ifdef ARK_USE_MOMENTUM
	Update_Tau_Momentum(ld.fDVDF[Ark::TAU], ld.msv[Ark::PHI], ld.msv[Ark::RHO]);
#endif

	AllenCahnTemporalPhiU(ld);
}
void DUGKSAmr::Update_MacroVars_hDt(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_MacroBaseVars_Phase(ld.msv_xF, ld.hDVDFxF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnPhase>::Update_MacroBaseVars_Phase(ld.msv_yF, ld.hDVDFyF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnPhase>::Update_MacroBaseVars_Phase(ld.msv_zF, ld.hDVDFzF[Ark::FLUX]);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv_xF, ld.fDVDFxF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv_yF, ld.fDVDFyF[Ark::FLUX]);
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv_zF, ld.fDVDFzF[Ark::FLUX]);
#endif

	AllenCahnSource_hDt(ld.msv_xF, ld.msv, IntVect::TheDimensionVector(Ark::X));
	AllenCahnSource_hDt(ld.msv_yF, ld.msv, IntVect::TheDimensionVector(Ark::Y));
	AllenCahnSource_hDt(ld.msv_zF, ld.msv, IntVect::TheDimensionVector(Ark::Z));

	Update_MacroSourceVars(ld.msv_xF, hDt(lev));
	Update_MacroSourceVars(ld.msv_yF, hDt(lev));
	Update_MacroSourceVars(ld.msv_zF, hDt(lev));
#ifdef ARK_USE_MOMENTUM
	Update_Tau_Momentum(ld.fDVDFxF[Ark::TAUF], ld.msv_xF[Ark::PHI], ld.msv_xF[Ark::RHO]);
	Update_Tau_Momentum(ld.fDVDFyF[Ark::TAUF], ld.msv_yF[Ark::PHI], ld.msv_yF[Ark::RHO]);
	Update_Tau_Momentum(ld.fDVDFzF[Ark::TAUF], ld.msv_zF[Ark::PHI], ld.msv_zF[Ark::RHO]);
#endif
	FillPhysicalBoundaryMacroVars(lev);
}