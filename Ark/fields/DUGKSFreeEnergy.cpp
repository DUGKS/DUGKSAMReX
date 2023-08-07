#include <DUGKSAmr.H>

// void DUGKSAmr::FreeEnergySource (int const lev)
// {
// 	auto &ld = Data(lev);

// 	const Real *dl = Geom(lev).CellSize();

// 	MultiFab const &RhoFab = ld.msv[Ark::RHO];
// 	MultiFab &GradRhoFab   = ld.msv[Ark::GradRHO];

// 	#ifdef _OPENMP
// 	#pragma omp parallel
// 	#endif
// 	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
// 	{
// 		const amrex::Box& bx = mfi.tilebox();

// 		//source
// 		Array4<Real const> const &rho = RhoFab.array(mfi);

// 		//target
// 		Array4<Real> const &grad_rho  = GradRhoFab.array(mfi);

// 		GradFunc::gradient_operator_isotropy(bx, dl, rho, grad_rho);
// 	}

// 	MultiFab &ForceFab   = ld.msv[Ark::FORCE];
// 	MultiFab &ChemiPFab  = ld.msv[Ark::CHEMIP];

// 	#ifdef _OPENMP
// 	#pragma omp parallel
// 	#endif
// 	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
// 	{
// 		const amrex::Box& bx = mfi.tilebox();

// 		Array4<Real const> const &rho       = RhoFab.array(mfi);
// 		Array4<Real const> const &grad_rho  = GradRhoFab.array(mfi);
// 		Array4<Real> const &chemi_p   = ChemiPFab.array(mfi);
// 		Array4<Real> const &force     = ForceFab.array(mfi);

// 		GradFunc::laplacian_operator_isotropy(bx, dl, rho, chemi_p);

// 		// GradFunc::laplacian_normal_chemip_isotropy(bx, dl, rho, grad_rho, chemi_p);

// 		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
// 		{
// 			const Real order = rho(i,j,k);

// 			const Real mu_rho = 
// 			2*Beta()*(order - RhoL())*(order - RhoV())*(2*order - RhoL() - RhoV())
// 			- Kappa()*chemi_p(i,j,k);

// 			// const Real mu_rho = Kappa()*chemi_p(i,j,k);

// 			chemi_p(i,j,k) = mu_rho;
						
// 			//F = Fs + G
// 			force(i,j,k,Ark::X) = mu_rho*grad_rho(i,j,k,Ark::X) + rho(i,j,k)*Gravity(Ark::X);
// 			force(i,j,k,Ark::Y) = mu_rho*grad_rho(i,j,k,Ark::Y) + rho(i,j,k)*Gravity(Ark::Y);
// #if(AMREX_SPACEDIM == 3)
// 			force(i,j,k,Ark::Z) = mu_rho*grad_rho(i,j,k,Ark::Z) + rho(i,j,k)*Gravity(Ark::Z);
// #endif
// 		});
// 	}
// }

// void DUGKSAmr::FreeEnergySource (int const lev)
// {
// 	auto &ld = Data(lev);

// 	const Real *dl = Geom(lev).CellSize();

// 	MultiFab const &RhoFab = ld.msv[Ark::RHO];
// 	MultiFab &GradRhoFab   = ld.msv[Ark::GradRHO];

// 	#ifdef _OPENMP
// 	#pragma omp parallel
// 	#endif
// 	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
// 	{
// 		const amrex::Box& bx = mfi.tilebox();

// 		//source
// 		Array4<Real const> const &rho = RhoFab.array(mfi);

// 		//target
// 		Array4<Real> const &grad_rho  = GradRhoFab.array(mfi);

// 		GradFunc::gradient_operator_isotropy(bx, dl, rho, grad_rho);
// 	}

// 	MultiFab &ForceFab   = ld.msv[Ark::FORCE];
// 	MultiFab &ChemiPFab  = ld.msv[Ark::CHEMIP];

// 	#ifdef _OPENMP
// 	#pragma omp parallel
// 	#endif
// 	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
// 	{
// 		const amrex::Box& bx = mfi.tilebox();

// 		Array4<Real const> const &rho = RhoFab.array(mfi);
// 		Array4<Real> const &chemi_p   = ChemiPFab.array(mfi);

// 		//continuous
// 		GradFunc::laplacian_operator_isotropy(bx, dl, rho, chemi_p);

// 		//discrete
// 		// GradFunc::laplacian_normal_chemip_isotropy(bx, dl, rho, grad_rho, chemi_p);

// 		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
// 		{
// 			const Real order = rho(i,j,k);

// 			//continuous double well
// 			// const Real mu_rho = 2*Beta()*EoS::ChemicalPotentialDW(order, RhoL(), RhoV())
// 			// 					- Kappa()*chemi_p(i,j,k);

// 			//discrete double well
// 			// const Real mu_rho = Kappa()*chemi_p(i,j,k);

// 			const Real mu_rho = EoS::ChemicalPotentialVdW(order, T0()) - Kappa()*chemi_p(i,j,k);

// 			chemi_p(i,j,k) = mu_rho;
// 		});
// 	}

// 	if(lev != finest_level) { AverageDownLevelPotential(lev); }
// 	FillPatchChemicalP(lev, cur_time, ChemiPFab);

// 	#ifdef _OPENMP
// 	#pragma omp parallel
// 	#endif
// 	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
// 	{
// 		const amrex::Box& bx = mfi.tilebox();

// 		Array4<Real const> const &rho  = RhoFab.array(mfi);

// 		Array4<Real> const &grad_rho = GradRhoFab.array(mfi);
// 		Array4<Real> const &chemi_p  = ChemiPFab.array(mfi);
// 		Array4<Real> const &force    = ForceFab.array(mfi);

// 		GradFunc::gradient_operator_isotropy(bx, dl, chemi_p, force);

// 		amrex::ParallelFor(bx, 
// 		AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
// 		{
// 			// well-balanced force for chemical potential model
// 			// force(i,j,k,n) = rho(i,j,k)*(Gravity(n) - force(i,j,k,n));

// 			// original force for chemical potential model
// 			force(i,j,k,n) =  -rho(i,j,k)*force(i,j,k,n) + Ark::DmVnMomentum::RT*grad_rho(i,j,k,n);
// 		});
// 	}
// }

void DUGKSAmr::FreeEnergySource (int const lev)
{
	auto &ld = Data(lev);

	const Real *dl = Geom(lev).CellSize();

	MultiFab const &RhoFab = ld.msv[Ark::RHO];
	MultiFab &GradRhoFab   = ld.msv[Ark::GradRHO];

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		//source
		Array4<Real const> const &rho = RhoFab.array(mfi);

		//target
		Array4<Real> const &grad_rho  = GradRhoFab.array(mfi);

		GradFunc::gradient_operator_isotropy(bx, dl, rho, grad_rho);
		// GradFunc::gradient_operator_isotropy_blend(bx, dl, rho, grad_rho);
		// GradFunc::gradient_operator_central(bx, dl, rho, grad_rho);
	}

	MultiFab &ForceFab   = ld.msv[Ark::FORCE];
	MultiFab &ChemiPFab  = ld.msv[Ark::CHEMIP];

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &rho = RhoFab.array(mfi);
		Array4<Real> const &chemi_p   = ChemiPFab.array(mfi);

		//continuous laplacian density
		GradFunc::laplacian_operator_isotropy(bx, dl, rho, chemi_p);

		//discrete
		// GradFunc::laplacian_normal_chemip_isotropy(bx, dl, rho, grad_rho, chemi_p);

		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
			const Real order = rho(i,j,k);

			// - continuous double well
			// const Real mu_rho = 2*Beta()*EoS::ChemicalPotentialDW(order, RhoL(), RhoV()) - Kappa()*chemi_p(i,j,k);				

			//discrete double well
			// const Real mu_rho = Kappa()*chemi_p(i,j,k);

			const Real mu_rho = EoS::ChemicalPotentialVdW(order, T0()) - Kappa()*chemi_p(i,j,k);

			chemi_p(i,j,k) = mu_rho;
		});
	}

	if(lev != finest_level) { AverageDownLevelPotential(lev); }
	FillPatchChemicalP(lev, cur_time, ChemiPFab);

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &rho  = RhoFab.array(mfi);

		Array4<Real> const &grad_rho = GradRhoFab.array(mfi);
		Array4<Real> const &chemi_p  = ChemiPFab.array(mfi);
		Array4<Real> const &force    = ForceFab.array(mfi);

		GradFunc::gradient_operator_isotropy(bx, dl, chemi_p, force);
		// GradFunc::gradient_operator_isotropy_blend(bx, dl, chemi_p, force);
		// GradFunc::gradient_operator_central(bx, dl, chemi_p, force);

		amrex::ParallelFor(bx, 
		AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
		{
			// well-balanced force for chemical potential model
			force(i,j,k,n) = rho(i,j,k)*(Gravity(n) - force(i,j,k,n));

			// original force for chemical potential model
			// force(i,j,k,n) =  Ark::DmVnMomentum::RT*grad_rho(i,j,k,n) - rho(i,j,k)*force(i,j,k,n);
		});
	}
}