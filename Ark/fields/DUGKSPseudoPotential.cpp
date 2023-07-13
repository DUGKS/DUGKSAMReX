#include <DUGKSAmr.H>

void DUGKSAmr::PseudoPotential (int const lev)
{
    auto &ld = Data(lev);

    const Real *dl = Geom(lev).CellSize();

    MultiFab const &RhoFab = ld.msv[Ark::RHO];
    MultiFab const &xEqFab = ld.fDVDF[Ark::EQ];
    MultiFab &PresFab      = ld.msv[Ark::P];
    MultiFab &GradRhoFab   = ld.msv[Ark::GradRHO];
    MultiFab &PseudoPFab   = ld.msv[Ark::CHEMIP];

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        // const amrex::Box& bx = mfi.tilebox();
        const amrex::Box& bx = mfi.growntilebox(1);

        //source
        Array4<Real const> const &rho = RhoFab.array(mfi);

        //target
        Array4<Real> const &pres      = PresFab.array(mfi);
        Array4<Real> const &pseudoP   = PseudoPFab.array(mfi);
        Array4<Real> const &grad_rho  = GradRhoFab.array(mfi);

        // used as the refinement criterion
        GradFunc::gradient_operator_isotropy(bx, dl, rho, grad_rho);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real order  = rho(i,j,k);

            pres(i,j,k) = EoS::PressureVdW(order, T0());

            //Shan-ChenPRE1994
            pseudoP(i,j,k) = std::sqrt(2*std::fabs(pres(i,j,k) - order*Ark::DmVnMomentum::RT));

            // ZhangPRE2003
            // pseudoP(i,j,k) = std::fabs(pres(i,j,k) - order*Ark::DmVnMomentum::RT);

            //well-balance PP
            // pseudoP(i,j,k) = std::sqrt(2*std::fabs(pres(i,j,k)));
        });
    }

    // if(lev != finest_level) { AverageDownLevelPotential(lev); }
    // FillPatchChemicalP(lev, cur_time, PseudoPFab);

    MultiFab &ForceFab   = ld.msv[Ark::FORCE];

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        //source
        Array4<Real const> const &rho = RhoFab.array(mfi);
        
        Array4<Real> const &pseudoP = PseudoPFab.array(mfi);
        Array4<Real> const &force   = ForceFab.array(mfi);
        Array4<Real> const &pres    = PresFab.array(mfi);

        GradFunc::gradient_operator_isotropy(bx, dl, pseudoP, force);

        amrex::ParallelFor(bx, 
        AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const Real order  = rho(i,j,k);

            // well-balanced Shan-Chen
            // force(i,j,k,n) *= -pseudoP(i,j,k)*Ark::sign(pres(i,j,k));

            //Shan-Chen
            force(i,j,k,n) *= -pseudoP(i,j,k)*Ark::sign(pres(i,j,k)-rho(i,j,k)*Ark::DmVnMomentum::RT);

            // Zhang
            // force(i,j,k,n) *= -Ark::sign(pres(i,j,k) - order*Ark::DmVnMomentum::RT);
        });
    }
}