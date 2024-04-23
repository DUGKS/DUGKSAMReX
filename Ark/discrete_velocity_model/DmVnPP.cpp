#include <DUGKSAmr.H>

#if ARK_USE_PHASE
#error "wrong velocity model!!!"
#endif

const std::string DUGKSAmr::DmVnModel = "D2Mn";

template <typename DVM>
template <long unsigned NUM>
void DUGKSCore<DVM>::Update_Eq_Phase (MultiFab &xEq, MultiFab &xSo, const Array<MultiFab, NUM> &msv)
{
    amrex::Abort("D2MnPP model involves no PHASE distribution.");
}

//Rho_0 == 0
template <typename DVM>
template <bool ENABLESOURCE>
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

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
        {
            const double rho  = Den(i,j,k);
            // const double p    = Pres(i,j,k);
            const double u    = Vel(i,j,k,Ark::X);
            const double v    = Vel(i,j,k,Ark::Y);

            const double fx   = Force(i,j,k,Ark::X);
            const double fy   = Force(i,j,k,Ark::Y);

            const double uu = u*u  + v*v;
            const double uF = u*fx + v*fy;

            const double MEq[DVM::Q] = {1,3*uu-2,magic[0]-3*uu,u,-u,v,-v,u*u-v*v,u*v};
            const double MSo[DVM::Q] = {0,6.0*uF,-6.0*uF,fx,-fx,fy,-fy,2.0*(fx*u-fy*v),fx*v+fy*u};

            for(int nn = 0;nn < DVM::Q;++nn) {
                fEq(i,j,k,nn) = rho*MEq[nn];
            }
            if constexpr (ENABLESOURCE) {
                for(int nn = 0;nn < DVM::Q;++nn) {
                    fSo(i,j,k,nn) = MSo[nn];
                }
            }
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
            rho(i,j,k)    += fdvdf(i,j,k,n);
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

        Array4<Real const> const &rho      = RhoFab.array(mfi);
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

    // FreeEnergySource(lev);
    PseudoPotential(lev);

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