#include <DUGKSAmr.H>

inline
double DUGKSAmr::SourcePhi (double const phi)
{	
	return -4*(phi - PhiL())*(phi - PhiV())/(WidthI()*(PhiL() - PhiV()));
}
inline
double DUGKSAmr::ChemicalPotential (const double phi,const double DivGradPhi)
{
	return
	(
		2*Beta()*(phi - PhiL())*(phi - PhiV())*(2*phi - PhiL() - PhiV()) - Kappa()*DivGradPhi
	);
}

void DUGKSAmr::AllenCahnSource (int const lev)
{
	auto &ld = Data(lev);

	const Real *dl = Geom(lev).CellSize();

	MultiFab const &PhiFab = ld.msv[Ark::PHI];
	MultiFab &GradPhiFab   = ld.msv[Ark::GradPHI];

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(PhiFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		//source
		Array4<Real const> const &phi = PhiFab.array(mfi);

		//target
		Array4<Real> const &grad_phi  = GradPhiFab.array(mfi);

		GradFunc::gradient_operator_isotropy(bx, dl, phi, grad_phi);


	}
	if(lev != finest_level) 
	{
		BoundGradient(GradPhiFab);
	}

#ifdef ARK_USE_MOMENTUM

	MultiFab &GradRhoFab = ld.msv[Ark::GradRHO];
	MultiFab &ForceFab   = ld.msv[Ark::FORCE];
	MultiFab &ChemiPFab  = ld.chemiP;

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(PhiFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &phi      = PhiFab.array(mfi);
		Array4<Real const> const &grad_phi = GradPhiFab.array(mfi);

		Array4<Real> const &grad_rho  = GradRhoFab.array(mfi);
		Array4<Real> const &chemi_p   = ChemiPFab.array(mfi);
		Array4<Real> const &force     = ForceFab.array(mfi);

		// Kim Analytical mu\nabla\phi
		GradFunc::laplacian_operator_isotropy(bx, dl, phi, chemi_p);

		// Kim discrete
		// GradFunc::laplacian_normal_operator_isotropy(bx, dl, phi, grad_phi, chemi_p);

		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
			const Real order = phi(i,j,k);

			grad_rho(i,j,k,Ark::X) = APhi()*grad_phi(i,j,k,Ark::X);
			grad_rho(i,j,k,Ark::Y) = APhi()*grad_phi(i,j,k,Ark::Y);
#if(AMREX_SPACEDIM == 3)
			grad_rho(i,j,k,Ark::Z) = APhi()*grad_phi(i,j,k,Ark::Z);
#endif
			// Kim Analytical
			const Real mu_phi = ChemicalPotential(order,chemi_p(i,j,k));

			// Kim discrete
			// const Real mu_phi = Kappa()*chemi_p(i,j,k);

			chemi_p(i,j,k) = mu_phi;
						
			//F = Fs
			force(i,j,k,Ark::X) = mu_phi*grad_phi(i,j,k,Ark::X);
			force(i,j,k,Ark::Y) = mu_phi*grad_phi(i,j,k,Ark::Y);
#if(AMREX_SPACEDIM == 3)
			force(i,j,k,Ark::Z) = mu_phi*grad_phi(i,j,k,Ark::Z);
#endif
		});
	}
#endif
}

// F_s = \phi\nabla\mu
// void DUGKSAmr::AllenCahnSource (int const lev)
// {
// 	auto &ld = Data(lev);

// 	const Real *dl = Geom(lev).CellSize();

// 	MultiFab const &PhiFab = ld.msv[Ark::PHI];
// 	MultiFab &GradPhiFab   = ld.msv[Ark::GradPHI];
// 	MultiFab &ChemiPFab    = ld.chemiP;

// 	#ifdef _OPENMP
// 	#pragma omp parallel
// 	#endif
// 	for(amrex::MFIter mfi(PhiFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
// 	{
// 		const amrex::Box& bx = mfi.tilebox();

// 		Array4<Real const> const &phi  = PhiFab.array(mfi);

// 		Array4<Real> const &grad_phi   = GradPhiFab.array(mfi);
// 		Array4<Real> const &chemi_p    = ChemiPFab.array(mfi);

// 		GradFunc::gradient_operator_isotropy (bx, dl, phi, grad_phi);
// 		GradFunc::laplacian_operator_isotropy (bx, dl, phi, chemi_p);
// 		// GradFunc::laplacian_normal_operator_isotropy(bx, dl, phi, grad_phi, chemi_p);

// 		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
// 		{
// 			const Real order = phi(i,j,k);
// 			chemi_p(i,j,k)   = ChemicalPotential(order,chemi_p(i,j,k));
// 		});
// 	}

// 	if(lev != finest_level) 
// 	{
// 		BoundGradient(GradPhiFab);
// 	}
// #ifdef ARK_USE_MOMENTUM
// 	FillPatchChemiP(lev, cur_time, ChemiPFab);

// 	MultiFab &GradRhoFab = ld.msv[Ark::GradRHO];
// 	MultiFab &ForceFab = ld.msv[Ark::FORCE];

// 	#ifdef _OPENMP
// 	#pragma omp parallel
// 	#endif
// 	for(amrex::MFIter mfi(PhiFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
// 	{
// 		const amrex::Box& bx = mfi.tilebox();

// 		Array4<Real const> const &phi      = PhiFab.array(mfi);
// 		Array4<Real const> const &grad_phi = GradPhiFab.array(mfi);

// 		Array4<Real> const &grad_rho  = GradRhoFab.array(mfi);
// 		Array4<Real> const &chemi_p   = ChemiPFab.array(mfi);
// 		Array4<Real> const &force     = ForceFab.array(mfi);

// 		GradFunc::gradient_operator_isotropy(bx, dl, chemi_p, force);

// 		amrex::ParallelFor(bx, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
// 		{
// 			const Real order = phi(i,j,k);

// 			grad_rho(i,j,k,Ark::X) = APhi()*grad_phi(i,j,k,Ark::X);
// 			grad_rho(i,j,k,Ark::Y) = APhi()*grad_phi(i,j,k,Ark::Y);
// #if(AMREX_SPACEDIM == 3)
// 			grad_rho(i,j,k,Ark::Z) = APhi()*grad_phi(i,j,k,Ark::Z);
// #endif			
// 			force(i,j,k,n) *= -phi(i,j,k);
// 		});
// 	}
// #endif
// }

void DUGKSAmr::AllenCahnSource_hDt (Array<MultiFab, Ark::MSVNUM> &msv_nF,
							 		Array<MultiFab, Ark::MSVNUM> const &msv,
							 		const IntVect &faceIndex)
{
	int const di  = faceIndex[0];
	int const di2 = di  + di;
	int const di3 = di2 + di;

	int const dj  = faceIndex[1];
	int const dj2 = dj  + dj;
	int const dj3 = dj2 + dj;

#if(AMREX_SPACEDIM == 3)
	int const dk  = faceIndex[2];
	int const dk2 = dk  + dk;
	int const dk3 = dk2 + dk;
#else
	int const dk = 0;
	int const dk2 = 0;
	int const dk3 = 0;
#endif

// #ifdef ARK_USE_MOMENTUM
// 	MultiFab &PhiFabnF = msv_nF[Ark::PHI];
// 	MultiFab &RhoFabnF = msv_nF[Ark::RHO];

// 	#ifdef _OPENMP
// 	#pragma omp parallel
// 	#endif
// 	for(amrex::MFIter mfi(PhiFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
// 	{
// 		const amrex::Box& bx = mfi.tilebox();

// 		//source
// 		Array4<Real const> const &phi = PhiFabnF.const_array(mfi);

// 		//destination
// 		Array4<Real> const &rho = RhoFabnF.array(mfi);

// 		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
// 		{
// 			rho(i,j,k) = APhi()*phi(i,j,k) + BPhi();
// 		});
// 	}

// #endif

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(msv_nF[Ark::GradPHI], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &grad_phi = msv[Ark::GradPHI].const_array(mfi);
		Array4<Real const> const &grad_rho = msv[Ark::GradRHO].const_array(mfi);
		Array4<Real const> const &force    = msv[Ark::FORCE].const_array(mfi);

		Array4<Real> const &grad_phi_nF = msv_nF[Ark::GradPHI].array(mfi);
		Array4<Real> const &grad_rho_nF = msv_nF[Ark::GradRHO].array(mfi);
		Array4<Real> const &force_nF    = msv_nF[Ark::FORCE].array(mfi);
		

		amrex::ParallelFor(bx, AMREX_SPACEDIM,
		[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
		{
			grad_phi_nF(i,j,k,n) = 
			GradFunc::double_interpolation(grad_phi(i-di2,j-dj2,k-dk2,n),
										    grad_phi(i-di ,j-dj ,k-dk ,n),
										    grad_phi(i    ,j    ,k    ,n),
										    grad_phi(i+di ,j+dj ,k+dk ,n),
											SchemeIndex());
			

			// grad_phi_nF(i,j,k,n) = 
			// GradFunc::triple_interpolation (grad_phi(i-di3,j-dj3,k-dk3,n),
			// 								grad_phi(i-di2,j-dj2,k-dk2,n),
			// 								grad_phi(i-di ,j-dj ,k-dk ,n),
			// 								grad_phi(i    ,j    ,k    ,n),
			// 								grad_phi(i+di ,j+dj ,k+dk ,n),
			// 								grad_phi(i+di2,j+dj2,k+dk2,n),
			// 								SchemeIndex());
			

#ifdef ARK_USE_MOMENTUM
			grad_rho_nF(i,j,k,n) = 0.5*(grad_rho(i,j,k,n) + grad_rho(i-di, j-dj, k-dk, n));
			force_nF(i,j,k,n)    = 0.5*(force(i,j,k,n)    + force(i-di, j-dj, k-dk, n));
			// grad_rho_nF(i,j,k,n) = 
			// GradFunc::central_interpolation(grad_rho(i-di2,j-dj2,k-dk2,n),
			// 							    grad_rho(i-di ,j-dj ,k-dk ,n),
			// 							    grad_rho(i    ,j    ,k    ,n),
			// 							    grad_rho(i+di ,j+dj ,k+dk ,n));

			// force_nF(i,j,k,n) =
			// GradFunc::central_interpolation(force(i-di2,j-dj2,k-dk2,n),
			// 							    force(i-di ,j-dj ,k-dk ,n),
			// 							    force(i    ,j    ,k    ,n),
			// 							    force(i+di ,j+dj ,k+dk ,n),
			// 							    SchemeIndex());
#endif
		});
	}
}

void DUGKSAmr::AllenCahnTemporalPhiU (DUGKSAmr::LevelData &ld)
{
	const Real RT = Ark::DmVnPhase::RT;

	MultiFab const& PhaseFab    = ld.msv[Ark::PHI];
	MultiFab const& VelocityFab = ld.msv[Ark::VEL];
//
	MultiFab &PhasePrevFab    = ld.phase_prev;
	MultiFab &VelocityPrevFab = ld.velocity_prev;
	MultiFab &GradPhiFab      = ld.msv[Ark::GradPHI];
//
	#ifdef _OPENMP
	#pragma omp parallel
  	#endif
  	for(amrex::MFIter mfi(PhaseFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
  	{
    	const amrex::Box& bx = mfi.tilebox();

    	Array4<Real const> const &phi = PhaseFab.const_array(mfi);
    	Array4<Real const> const &vel = VelocityFab.const_array(mfi);

    	Array4<Real> const &phi_prev  = PhasePrevFab.array(mfi);
    	Array4<Real> const &vel_prev  = VelocityPrevFab.array(mfi);
    	Array4<Real> const &grad_phi = GradPhiFab.array(mfi);
    	//
    	amrex::For(bx,
    	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    	{
    		Real const &phi_value = phi(i,j,k);
    		Real &phi_prev_value = phi_prev(i,j,k);

    		Real const &vel_x = vel(i,j,k,Ark::X);
    		Real const &vel_y = vel(i,j,k,Ark::Y);

    		Real &vel_prev_x = vel_prev(i,j,k,Ark::X);
    		Real &vel_prev_y = vel_prev(i,j,k,Ark::Y);

    		Real const dPhiU = (phi_value*vel_x - phi_prev_value*vel_prev_x)/(RT*Dt());
    		Real const dPhiV = (phi_value*vel_y - phi_prev_value*vel_prev_y)/(RT*Dt());

#if(AMREX_SPACEDIM == 3)
    		Real const &vel_z      = vel(i,j,k,Ark::Z);
    		Real       &vel_prev_z = vel_prev(i,j,k,Ark::Z);
    		Real const dPhiW = (phi_value*vel_z - phi_prev_value*vel_prev_z)/(RT*Dt());

    		vel_prev_z = vel_z;
#endif
    		vel_prev_x = vel_x;
    		vel_prev_y = vel_y;

    		Real &grad_phi_x = grad_phi(i,j,k,Ark::X);
    		Real &grad_phi_y = grad_phi(i,j,k,Ark::Y); 

#if(AMREX_SPACEDIM == 3)
    		Real &grad_phi_z = grad_phi(i,j,k,Ark::Z);
    		
    		//c++ 17
    		//Real const mod_phi = std::hypot(grad_phi_x, grad_phi_y, grad_phi_z);
    		Real const mod_phi = Ark::Hypot(grad_phi_x, grad_phi_y, grad_phi_z);
#else
    		Real const mod_phi = std::hypot(grad_phi_x, grad_phi_y);
#endif		
    		if(fabs(mod_phi) > Ark::tiny)
    		{
    			Real const mod_phi_eq = SourcePhi(phi_value);
    			grad_phi_x *= mod_phi_eq/mod_phi;
    			grad_phi_y *= mod_phi_eq/mod_phi;
#if(AMREX_SPACEDIM == 3)
    			grad_phi_z *= mod_phi_eq/mod_phi;
#endif
    		}
    			grad_phi_x += dPhiU;
    			grad_phi_y += dPhiV;

#if(AMREX_SPACEDIM == 3)
    			grad_phi_z += dPhiW;
#endif
    		phi_prev_value = phi_value;
    	});
  	}
}
