#include <DUGKSAmr.H>

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
			const Real u    = vel_array(i,j,k, Ark::dir_x);
			const Real v    = vel_array(i,j,k, Ark::dir_y);
			const Real w    = vel_array(i,j,k, Ark::dir_z);
			const Real phix = gradphi_array(i,j,k, Ark::dir_x);
			const Real phiy = gradphi_array(i,j,k, Ark::dir_y);
			const Real phiz = gradphi_array(i,j,k, Ark::dir_z);

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
			const double &u    = Vel(i,j,k,Ark::dir_x);
     		const double &v    = Vel(i,j,k,Ark::dir_y);
     		const double &w    = Vel(i,j,k,Ark::dir_z);
       		const double &rhox = GradDen(i,j,k,Ark::dir_x);
       		const double &rhoy = GradDen(i,j,k,Ark::dir_y);
       		const double &rhoz = GradDen(i,j,k,Ark::dir_z);
       		const double &fx   = Force(i,j,k,Ark::dir_x);
       		const double &fy   = Force(i,j,k,Ark::dir_y);
       		const double &fz   = Force(i,j,k,Ark::dir_z);

			const double uu = Ark::Square(u,v,w);

			const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v + DVM::xi_w[n]*w)/DVM::RT;
			const double Gamma = u1 + 0.5*u1*u1 - uu*DVM::Lambda;

			fEq(i,j,k,n) = DVM::omega[n]*(p/DVM::RT + rho*Gamma);

			#ifdef _ARK_M_LSGC_FLIP
			fEq(i,j,k,n) -= static_cast<int>(n == 0)*p;
			#endif

			fSo(i,j,k,n) = DVM::omega[n]*
					   (
					   (1+Gamma)*((DVM::xi_u[n] - u)*fx + (DVM::xi_v[n] - v)*fy 
					   									+ (DVM::xi_w[n] - w)*fz)/DVM::RT
					   +
					   Gamma*((DVM::xi_u[n] - u)*rhox + (DVM::xi_v[n] - v)*rhoy
					   								  + (DVM::xi_w[n] - w)*rhoz)
					   );
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
			pres(i,j,k) 		  += fdvdf(i,j,k,n);
			vel(i,j,k,Ark::dir_x) += fdvdf(i,j,k,n)*DVM::xi_u[n];
			vel(i,j,k,Ark::dir_y) += fdvdf(i,j,k,n)*DVM::xi_v[n];
			vel(i,j,k,Ark::dir_z) += fdvdf(i,j,k,n)*DVM::xi_w[n];
		});
	}
}
void DUGKSAmr::Update_MacroSourceVars(Array<MultiFab, Ark::MSVNUM> &msv, 
									  const Real hTime, const Real RT)
{
#if ARK_USE_MOMENTUM
	const MultiFab &PhiFab     = msv[Ark::PHI];
	const MultiFab &ForceFab   = msv[Ark::FORCE];
	const MultiFab &GradRhoFab = msv[Ark::GradRHO];

	//density
	MultiFab &RhoFab  = msv[Ark::RHO];
	MultiFab &PresFab = msv[Ark::P];
	MultiFab::Saxpy(RhoFab, APhi(), PhiFab,0,0,1,0);
	RhoFab.plus(BPhi(), 0);

	MultiFab &VelFab = msv[Ark::VEL];

	MultiFab::Saxpy(VelFab, 0.5*hTime, ForceFab, 0, 0, AMREX_SPACEDIM, 0);
	MultiFab::Divide(VelFab, RhoFab, 0, 0, 1, 0);
	MultiFab::Divide(VelFab, RhoFab, 0, 1, 1, 0);
	MultiFab::Divide(VelFab, RhoFab, 0, 2, 1, 0);


	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &grad_rho = GradRhoFab.const_array(mfi);

		Array4<Real const> const &vel  = VelFab.array(mfi);
		Array4<Real> const &pres = PresFab.array(mfi);

		amrex::ParallelFor(bx,
		[=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
		{
			pres(i,j,k) +=  0.5*hTime*
						   	(
								vel(i,j,k,Ark::dir_x)*grad_rho(i,j,k,Ark::dir_x)
							+	vel(i,j,k,Ark::dir_y)*grad_rho(i,j,k,Ark::dir_y)
							+	vel(i,j,k,Ark::dir_z)*grad_rho(i,j,k,Ark::dir_z)
							);
			#ifdef _ARK_M_LSGC_FLIP
			msq(i,j,k,Ark::P) -= 2.0/9.0*msq(i,j,k,Ark::RHO)*msq(i,j,k,Ark::SQUV);
			msq(i,j,k,Ark::P) *= 9.0/5.0;
			#endif
		});
	}

	PresFab.mult(RT, 0);
#endif
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

	Update_MacroSourceVars(ld.msv, hDt(lev), Ark::DmVnMomentum::RT);

	AllenCahnTemporalPhiU(ld, Ark::DmVnPhase::RT);
}
void DUGKSAmr::Update_MacroVars_hDt(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_MacroBaseVars_Phase(ld.msv_xF, ld.hDVDFxF[Ark::BHF]);
	DUGKSCore<Ark::DmVnPhase>::Update_MacroBaseVars_Phase(ld.msv_yF, ld.hDVDFyF[Ark::BHF]);
	DUGKSCore<Ark::DmVnPhase>::Update_MacroBaseVars_Phase(ld.msv_zF, ld.hDVDFzF[Ark::BHF]);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv_xF, ld.fDVDFxF[Ark::BHF]);
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv_yF, ld.fDVDFyF[Ark::BHF]);
	DUGKSCore<Ark::DmVnMomentum>::Update_MacroBaseVars_Momentum(ld.msv_zF, ld.fDVDFzF[Ark::BHF]);
#endif

	AllenCahnSource_hDt(ld.msv_xF, ld.msv, IntVect::TheDimensionVector(Ark::dir_x));
	AllenCahnSource_hDt(ld.msv_yF, ld.msv, IntVect::TheDimensionVector(Ark::dir_y));
	AllenCahnSource_hDt(ld.msv_zF, ld.msv, IntVect::TheDimensionVector(Ark::dir_z));

	Update_MacroSourceVars(ld.msv_xF, hDt(lev), Ark::DmVnMomentum::RT);
	Update_MacroSourceVars(ld.msv_yF, hDt(lev), Ark::DmVnMomentum::RT);
	Update_MacroSourceVars(ld.msv_zF, hDt(lev), Ark::DmVnMomentum::RT);

	FillPhysicalBoundaryMacroVars(lev);
}