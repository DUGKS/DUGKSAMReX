#include <DUGKSAmr.H>

const std::string DUGKSAmr::DmVnModel = "DmVnCHSS";

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

        Array4<Real const> const &phi_array = msv[Ark::PHI].const_array(mfi);
        Array4<Real const> const &vel_array = msv[Ark::VEL].const_array(mfi);
        Array4<Real const> const &gradphi_array = msv[Ark::GradPHI].const_array(mfi);
        Array4<Real const> const &chemi_potential = msv[Ark::CHEMIP].const_array(mfi);

        Array4<Real> const &hEq = xEq.array(mfi);
        Array4<Real> const &hSo = xSo.array(mfi);

        ParallelFor(tbx, DVM::Q, [=] AMREX_GPU_DEVICE (int i,int j,int k,int n) noexcept
        {
            const double phi  = phi_array(i,j,k);
            const double chemiP = chemi_potential(i,j,k);
            const double u    = vel_array(i,j,k, Ark::X);
            const double v    = vel_array(i,j,k, Ark::Y);
            const double phix = gradphi_array(i,j,k, Ark::X);
            const double phiy = gradphi_array(i,j,k, Ark::Y);
#if(AMREX_SPACEDIM == 3)
            const double w    = vel_array(i,j,k, Ark::Z);
            const double phiz = gradphi_array(i,j,k, Ark::Z);
            const double uu   = Ark::Square(u,v,w);
            const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v + DVM::xi_w[n]*w)/DVM::RT;
            const double so = DVM::xi_u[n]*phix + DVM::xi_v[n]*phiy + DVM::xi_w[n]*phiz;
#elif(AMREX_SPACEDIM == 2)
            const double uu = Ark::Square(u,v);
            const double u1 = (DVM::xi_u[n]*u + DVM::xi_v[n]*v)/DVM::RT;
            const double so = DVM::xi_u[n]*phix + DVM::xi_v[n]*phiy;
#endif

            hEq(i,j,k,n) = DVM::omega[n]*(chemiP + phi*u1) + static_cast<int>(n == 0)*(phi-chemiP);
            // hSo(i,j,k,n) = DVM::omega[n]*so;
            hSo(i,j,k,n) = 0;
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

void DUGKSAmr::Update_MacroBaseVars (Array<MultiFab, Ark::MSVNUM> &msv)
{
#ifdef ARK_USE_MOMENTUM
    const MultiFab &RhoFab   = msv[Ark::RHO];
    
    MultiFab &VelFab   = msv[Ark::VEL];
    MultiFab &PresFab  = msv[Ark::P];

    PresFab.mult(Ark::DmVnMomentum::RT, 0);

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

    MultiFab const &RhoFab     = msv[Ark::RHO];
    MultiFab const &GradRhoFab = msv[Ark::GradRHO];
    MultiFab const &ForceFab   = msv[Ark::FORCE];

    MultiFab &VelFab  = msv[Ark::VEL];
    MultiFab &PresFab = msv[Ark::P];

    for(amrex::MFIter mfi(PresFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        Array4<Real const> const &rho      = RhoFab.array(mfi);
        Array4<Real const> const &grad_rho = GradRhoFab.const_array(mfi);
        Array4<Real const> const &force    = ForceFab.array(mfi);

        Array4<Real> const &vel   = VelFab.array(mfi);
        Array4<Real> const &pres  = PresFab.array(mfi);

        amrex::ParallelFor(bx, AMREX_SPACEDIM,[=] AMREX_GPU_DEVICE (int i,int j,int k, int n) noexcept
        {
            vel(i,j,k,n) += interval*force(i,j,k,n)/rho(i,j,k);
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept
        {
            pres(i,j,k) +=  interval*Ark::DmVnMomentum::RT*
                            (
                                vel(i,j,k,Ark::X)*grad_rho(i,j,k,Ark::X)
                            +   vel(i,j,k,Ark::Y)*grad_rho(i,j,k,Ark::Y)
                            );
        });     
    }

#endif
}