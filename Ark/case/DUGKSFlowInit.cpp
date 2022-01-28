#include <DUGKSAmr.H>

void DUGKSAmr::BasicCase(int const lev)
{
    AMREX_ALWAYS_ASSERT(CaseIndex() == -666 && CaseName() == "BasicCase");

#ifdef AMREX_DEBUG

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();

    for(amrex::MFIter mfi(ld.msv[Ark::RHO], true);mfi.isValid();++mfi)
    {
        const Box& tbx = mfi.tilebox();
        // const Box& tbx = mfi.growntilebox(1);
        // const Box& tbx = mfi.fabbox();
        //
        Array4<Real> const &phi  = ld.msv[Ark::PHI].array(mfi);
        Array4<Real> const &rho  = ld.msv[Ark::RHO].array(mfi);
        Array4<Real> const &vel  = ld.msv[Ark::VEL].array(mfi);
        Array4<Real> const &pres = ld.msv[Ark::P].array(mfi);
        Array4<Real> const &temp = ld.msv[Ark::T].array(mfi);
        Array4<Real> const &grad_phi = ld.msv[Ark::GradPHI].array(mfi);
        //
        Array4<Real> const &hWave = ld.hDVDF[Ark::WAVE].array(mfi);
        Array4<Real> const &fWave = ld.fDVDF[Ark::WAVE].array(mfi);

        Array4<Real> const &hBP = ld.hDVDF[Ark::BP].array(mfi);
        Array4<Real> const &fBP = ld.fDVDF[Ark::BP].array(mfi);



        ParallelFor(tbx, Ark::DVM_PHASE, [=](int i,int j,int k,int n) noexcept
        {
            const double x = dx[Ark::X]*(i+0.5) + rb.lo(Ark::X);
            const double y = dx[Ark::Y]*(j+0.5) + rb.lo(Ark::Y);
#if(AMREX_SPACEDIM == 3)
            const double z = dx[Ark::Z]*(k+0.5) + rb.lo(Ark::Z);
#else
            const double z = 0;
#endif
            phi(i,j,k) = 2*x*x + y*y;
            rho(i,j,k) = 2*x + x*y + y*z;

            vel(i,j,k,0) = 0.1*x + x*y;
            vel(i,j,k,1) = 1.1*x + x*y;

            pres(i,j,k) = 2*x + 3*y;
            temp(i,j,k) = 2*x + 4*y;

            hWave(i,j,k,n) = (2.0  *x + (0.1 + n*0.01)*y);
            hBP(i,j,k,n)   = (2.0 * x + (0.1 + n*0.01)* y + 1);

            fWave(i,j,k,n) = (2.0 * x *x + (0.1 + n*0.01)*x* y + 2);
            fBP(i,j,k,n)   = (2.0  *x*x + (0.1 + n*0.01)*x* y + 3);
        });
    }

    for(amrex::MFIter mfi(ld.hDVDFxF[Ark::FLUX],true);mfi.isValid();++mfi)
    {
        const Box& tbx = mfi.tilebox();
        //
        Array4<Real> const &hFluxX = ld.hDVDFxF[Ark::FLUX].array(mfi);
        Array4<Real> const &fFluxX = ld.fDVDFxF[Ark::FLUX].array(mfi);


        ParallelFor(tbx, Ark::DVM_PHASE, [=](int i,int j,int k,int n) noexcept
        {
            const double x = dx[Ark::X]*(i);
            const double y = dx[Ark::Y]*(j+0.5);

            hFluxX(i,j,k,n) = (2.0 * x + (2 + n)*y);
            fFluxX(i,j,k,n) = ((2.0 +n)*x + n*0.1*y);
        });
    }
#endif
}
//
void DUGKSAmr::UniformFlow(int const lev)
{
	AMREX_ALWAYS_ASSERT((CaseIndex() == 0 && CaseName() == "UniformFlow"));

    LevelData &ld = Data(lev);

	ld.msv[Ark::PHI].setVal(Rho0());
	ld.msv[Ark::RHO].setVal(Rho0());
	ld.msv[Ark::VEL].setVal(0);
	ld.msv[Ark::P].setVal(P0());
	ld.msv[Ark::T].setVal(T0());
}
void DUGKSAmr::UnsteadyTaylorGreenVortex(int const lev)
{
    AMREX_ALWAYS_ASSERT((CaseIndex() == 1 && CaseName() == "UnsteadyTaylorGreenVortex"));

    auto &ld = Data(lev);

    ld.msv[Ark::RHO].setVal(Rho0());

    FillLevelUnsteadyTaylorGreenVortex(lev, ld.msv[Ark::VEL], ld.msv[Ark::P]);
    
    amrex::MultiFab::Saxpy(ld.msv[Ark::RHO], 2*Ark::DmVnMomentum::Lambda, ld.msv[Ark::P], 0, 0, 1, 0);
}
void DUGKSAmr::ForceTaylorGreenVortex(int const lev)
{
    AMREX_ALWAYS_ASSERT((CaseIndex() == 5 && CaseName() == "ForceTaylorGreenVortex"));

    auto &ld = Data(lev);

    ld.msv[Ark::RHO].setVal(Rho0());

    FillLevelForceTaylorGreenVortex(lev, ld.msv[Ark::VEL], ld.msv[Ark::P], ld.msv[Ark::FORCE],-1);

    FillLevelForceTaylorGreenVortex(lev, ld.msv_xF[Ark::VEL], ld.msv_xF[Ark::P], 
                                    ld.msv_xF[Ark::FORCE], Ark::X);
                                    
    FillLevelForceTaylorGreenVortex(lev, ld.msv_yF[Ark::VEL], ld.msv_yF[Ark::P],
                                    ld.msv_yF[Ark::FORCE], Ark::Y);

    FillLevelForceTaylorGreenVortex(lev, ld.velocity_prev, ld.phase_prev, ld.msv[Ark::FORCE], -1);

    amrex::MultiFab::Saxpy(ld.msv[Ark::RHO], 2*Ark::DmVnMomentum::Lambda, ld.msv[Ark::P], 0, 0, 1, 0);
}
//**********************two phase flows*****************************
void DUGKSAmr::ReversedSingleVortex(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 219 && CaseName() == "ReversedSingleVortex"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();

    FillLevelPhaseSphere(ld.msv[Ark::PHI], dx, LV::drop);
    FillLevelPhaseSphere(ld.density_old, dx, LV::drop);

    MultiFab::Copy(ld.msv[Ark::RHO], ld.msv[Ark::PHI], 0, 0, 1, zero_ghost());
    ld.msv[Ark::RHO].mult(APhi(), 0);
    ld.msv[Ark::RHO].plus(BPhi(), 0);

    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    {
    MultiFab &VelFab    = ld.msv[Ark::VEL];
    for(amrex::MFIter mfi(VelFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &vel = VelFab.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x = dx[0]*(i+0.5)/ChLength();
            const double y = dx[1]*(j+0.5)/ChLength();

            vel(i,j,k,Ark::X) =  U0()*sin(PI*x)*sin(PI*x)*sin(2.0*PI*y);
            vel(i,j,k,Ark::Y) = -U0()*sin(PI*y)*sin(PI*y)*sin(2.0*PI*x);           
        });
    }}
}
void DUGKSAmr::QuiescentDrop(int const lev)
{
    AMREX_ALWAYS_ASSERT((CaseIndex() == 220 && CaseName() == "QuiescentDrop"));

    LevelData &ld = Data(lev);
    const Real *dl = Geom(lev).CellSize();

	FillLevelPhaseSphere(ld.msv[Ark::PHI], dl, LV::drop);
    FillLevelPhaseSphere(ld.density_old, dl, LV::drop);

	MultiFab::Copy(ld.msv[Ark::RHO], ld.msv[Ark::PHI], 0, 0, 1, zero_ghost());
	ld.msv[Ark::RHO].mult(APhi(), 0);
    ld.msv[Ark::RHO].plus(BPhi(), 0);
    // ld.msv[Ark::VEL].setVal(0);
	ld.msv[Ark::P].setVal(P0());
	ld.msv[Ark::T].setVal(T0());

    ld.msv[Ark::VEL].setVal(U0());
    ld.msv_xF[Ark::VEL].setVal(U0());
    ld.msv_yF[Ark::VEL].setVal(U0());
#if(AMREX_SPACEDIM == 3)
    ld.msv_zF[Ark::VEL].setVal(U0());
#endif
}
void DUGKSAmr::TransDrop(int const lev)
{
    AMREX_ALWAYS_ASSERT((CaseIndex() == 221 && CaseName() == "TransDrop"));

    LevelData &ld = Data(lev);
    const Real *dl = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();

    FillLevelPhaseSphere(ld.msv[Ark::PHI], dl, LV::drop);
    FillLevelPhaseSphere(ld.density_old, dl, LV::drop);

    MultiFab::Copy(ld.msv[Ark::RHO], ld.msv[Ark::PHI], 0, 0, 1, zero_ghost());
    ld.msv[Ark::RHO].mult(APhi(), 0);
    ld.msv[Ark::RHO].plus(BPhi(), 0);

    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    ld.msv[Ark::VEL].setVal(U0());
    ld.msv_xF[Ark::VEL].setVal(U0());
    ld.msv_yF[Ark::VEL].setVal(U0());
#if(AMREX_SPACEDIM == 3)
    ld.msv_zF[Ark::VEL].setVal(U0());
#endif
}
void DUGKSAmr::ZalesakDisk(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 222 && CaseName() == "ZalesakDisk"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();

    FillLevelPhaseZalesakDisk(ld.msv[Ark::PHI], dx);
    FillLevelPhaseZalesakDisk(ld.density_old, dx);

    MultiFab::Copy(ld.msv[Ark::RHO], ld.msv[Ark::PHI], 0, 0, 1, zero_ghost());
    ld.msv[Ark::RHO].mult(APhi(), 0);
    ld.msv[Ark::RHO].plus(BPhi(), 0);

    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    ld.msv[Ark::VEL].setVal(0);

    {
    MultiFab &VelFab    = ld.msv[Ark::VEL];
    for(amrex::MFIter mfi(VelFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &vel = VelFab.array(mfi);
        //
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x = dx[0]*(i+0.5);
            const double y = dx[1]*(j+0.5);

            vel(i,j,k,Ark::X) = -U0()*pi*(y-CenterPF(Ark::Y))/ChLength();
            vel(i,j,k,Ark::Y) =  U0()*pi*(x-CenterPF(Ark::X))/ChLength();           
        });
    }}
    {
    MultiFab &VelFab_xF = ld.msv_xF[Ark::VEL];
    for(amrex::MFIter mfi(VelFab_xF);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &vel = VelFab_xF.array(mfi);
        //
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x = dx[0]*(i);
            const double y = dx[1]*(j+0.5);

            vel(i,j,k,Ark::X) = -U0()*pi*(y-CenterPF(Ark::Y))/ChLength();
            vel(i,j,k,Ark::Y) =  U0()*pi*(x-CenterPF(Ark::X))/ChLength();
            
        });
    }}
    {
        MultiFab &VelFab_yF = ld.msv_yF[Ark::VEL];
        for(amrex::MFIter mfi(VelFab_yF);mfi.isValid();++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
        //
            Array4<Real> const &vel = VelFab_yF.array(mfi);
        //
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
            {
                const double x = dx[0]*(i+0.5);
                const double y = dx[1]*(j);

                vel(i,j,k,Ark::X) = -U0()*pi*(y-CenterPF(Ark::Y))/ChLength();
                vel(i,j,k,Ark::Y) =  U0()*pi*(x-CenterPF(Ark::X))/ChLength();            
            });
    }}
#if(AMREX_SPACEDIM == 3)
    {
        MultiFab &VelFab_zF = ld.msv_zF[Ark::VEL];
        for(amrex::MFIter mfi(VelFab_zF);mfi.isValid();++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
        //
            Array4<Real> const &vel = VelFab_zF.array(mfi);
        //
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
            {
                const double x = dx[0]*(i+0.5);
                const double y = dx[1]*(j+0.5);

                vel(i,j,k,Ark::X) = -U0()*pi*(y-CenterPF(Ark::Y))/ChLength();
                vel(i,j,k,Ark::Y) =  U0()*pi*(x-CenterPF(Ark::X))/ChLength();         
            });
    }}
#endif
}
void DUGKSAmr::ElongatedInterface(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() ==223 && CaseName() == "ElongatedInterface"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();

    FillLevelPhaseSphere(ld.msv[Ark::PHI], dx, LV::drop);
    FillLevelPhaseSphere(ld.density_old, dx, LV::drop);

    MultiFab::Copy(ld.msv[Ark::RHO], ld.msv[Ark::PHI], 0, 0, 1, zero_ghost());
    ld.msv[Ark::RHO].mult(APhi(), 0);
    ld.msv[Ark::RHO].plus(BPhi(), 0);

    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

#if(AMREX_SPACEDIM == 2)
    ElongationVelocityField (ld.msv[Ark::VEL],    dx, rb, -1);
    ElongationVelocityField (ld.msv_xF[Ark::VEL], dx, rb, Ark::X);
    ElongationVelocityField (ld.msv_yF[Ark::VEL], dx, rb, Ark::Y);
#elif(AMREX_SPACEDIM == 3)
    ElongationVelocityField3D (ld.msv[Ark::VEL],    dx, rb, -1);
    ElongationVelocityField3D (ld.msv_xF[Ark::VEL], dx, rb, Ark::X);
    ElongationVelocityField3D (ld.msv_yF[Ark::VEL], dx, rb, Ark::Y);
    ElongationVelocityField3D (ld.msv_zF[Ark::VEL], dx, rb, Ark::Z);
#endif

}
void DUGKSAmr::SmoothDeformation(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 224 && CaseName() == "SmoothDeformation"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();

    FillLevelPhaseSphere(ld.msv[Ark::PHI], dx, LV::drop);
    FillLevelPhaseSphere(ld.density_old, dx, LV::drop);

    MultiFab::Copy(ld.msv[Ark::RHO], ld.msv[Ark::PHI], 0, 0, 1, zero_ghost());
    ld.msv[Ark::RHO].mult(APhi(), 0);
    ld.msv[Ark::RHO].plus(BPhi(), 0);

    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    {
    MultiFab &VelFab    = ld.msv[Ark::VEL];
    for(amrex::MFIter mfi(VelFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &vel = VelFab.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x = dx[0]*(i+0.5)/ChLength();
            const double y = dx[1]*(j+0.5)/ChLength();

            vel(i,j,k,Ark::X) = -U0()*sin(4.0*PI*x)*sin(4.0*PI*y);
            vel(i,j,k,Ark::Y) = -U0()*cos(4.0*PI*x)*cos(4.0*PI*y);            
        });
    }}
}
void DUGKSAmr::LayeredPoiseuille (int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 225 && CaseName() == "LayeredPoiseuille")); 

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();

    ld.msv[Ark::VEL].setVal(0);
    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    MultiFab &PhiFab = ld.msv[Ark::PHI];
    MultiFab &RhoFab = ld.msv[Ark::RHO];

    for(amrex::MFIter mfi(PhiFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &phi = PhiFab.array(mfi);
        Array4<Real> const &rho = RhoFab.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double y = dx[Ark::Y]*(j+0.5) + rb.lo(Ark::Y);

            phi(i,j,k) = 0.5*(PhiL() + PhiV()) + 0.5*(PhiL() - PhiV())
                         *tanh((CenterPF(Ark::Y) - 2.*y)/WidthI());
            rho(i,j,k) = 0.5*(RhoL() + RhoV()) + 0.5*(RhoL() - RhoV())
                         *tanh((CenterPF(Ark::Y) - 2.*y)/WidthI());          
        });
    }
}
void DUGKSAmr::SpinodalDecomposition(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 226 && CaseName() == __func__)); 

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();

    ld.msv[Ark::VEL].setVal(0);
    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    MultiFab &PhiFab = ld.msv[Ark::PHI];
    MultiFab &RhoFab = ld.msv[Ark::RHO];

    for(amrex::MFIter mfi(PhiFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &phi = PhiFab.array(mfi);
        Array4<Real> const &rho = RhoFab.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double y = dx[Ark::Y]*(j+0.5) + rb.lo(Ark::Y);

            phi(i,j,k) = 0.6 + (amrex::Random())*0.01;
            rho(i,j,k) = APhi()*phi(i,j,k) + BPhi();          
        });
    }
}
void DUGKSAmr::RayleighTaylorInstability(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 227 && CaseName() == "RayleighTaylorInstability")); 

    _Peclet = std::sqrt(std::fabs(Gravity(Ark::Y))*ChLength())*ChLength()/MPhi();

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();

    ld.msv[Ark::VEL].setVal(0);
    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    MultiFab &PhiFab = ld.msv[Ark::PHI];
    MultiFab &RhoFab = ld.msv[Ark::RHO];

    for(amrex::MFIter mfi(PhiFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &phi = PhiFab.array(mfi);
        Array4<Real> const &rho = RhoFab.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x = dx[Ark::X]*(i+0.5) + rb.lo(Ark::X);
            const double y = dx[Ark::Y]*(j+0.5) + rb.lo(Ark::Y);
#if(AMREX_SPACEDIM == 3)
            const double z = dx[Ark::Z]*(k+0.5) + rb.lo(Ark::Z);
            const double midY = CenterPF(Ark::Y) 
            + 0.05*ChLength()*(cos(2.*PI*x/ChLength()) + cos(2.*PI*z/ChLength()));
#else
            const double midY = CenterPF(Ark::Y) + 0.05*ChLength()*cos(2.*PI*x/ChLength());
#endif
            phi(i,j,k) = 0.5*(PhiL() + PhiV()) + 0.5*(PhiL() - PhiV())
                         *tanh((y - midY)*2/WidthI());
            rho(i,j,k) = APhi()*phi(i,j,k) + BPhi();          
        });
    }
}
void DUGKSAmr::FlatInterface(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 228 && CaseName() == "FlatInterface"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();
    const Box &domain = Geom(lev).Domain();

    ld.msv[Ark::VEL].setVal(0);
    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    MultiFab &PhiFab = ld.msv[Ark::PHI];
    MultiFab &RhoFab = ld.msv[Ark::RHO];
    MultiFab &RhoOldFab = ld.density_old;

    for(amrex::MFIter mfi(PhiFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &phi = PhiFab.array(mfi);
        Array4<Real> const &rho = RhoFab.array(mfi);
        Array4<Real> const &rho_old = RhoOldFab.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double y = dx[Ark::Y]*(j+0.5) + rb.lo(Ark::Y);

            const double yA = CenterPF(Ark::Y) - 0.25*rb.length(Ark::Y);
            const double yB = CenterPF(Ark::Y) + 0.25*rb.length(Ark::Y);

            phi(i,j,k) = PhiV() + 0.5*(PhiL() - PhiV())
                       *(tanh(2.*(y-yA)/WidthI()) - tanh(2.*(y-yB)/WidthI()));
            rho(i,j,k) = RhoV() + 0.5*(RhoL() - RhoV())
                       *(tanh(2.*(y-yA)/WidthI()) - tanh(2.*(y-yB)/WidthI()));
            rho_old(i,j,k) = rho(i,j,k);          
        });
    }
}

void DUGKSAmr::DropletImpact(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 229 && CaseName() == "DropletImpact"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();
    const Box &domain = Geom(lev).Domain();

    ld.msv[Ark::VEL].setVal(0);
    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    MultiFab &PhiFab = ld.msv[Ark::PHI];
    MultiFab &RhoFab = ld.msv[Ark::RHO];
    MultiFab &VelFab = ld.msv[Ark::VEL];

    for(amrex::MFIter mfi(PhiFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &phi = PhiFab.array(mfi);
        Array4<Real> const &rho = RhoFab.array(mfi);
        Array4<Real> const &vel = VelFab.array(mfi);

        const double Hw = 50.0;
        const double wI = WidthI();

        const Real centerX = CenterPF(Ark::X);
        const Real centerY = Hw + RadiusPF() + 2*wI;
#if(AMREX_SPACEDIM == 3)
        const Real centerZ = CenterPF(Ark::Z);
#endif
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x  = dx[Ark::X]*(i+0.5) + rb.lo(Ark::X);
            const double y  = dx[Ark::Y]*(j+0.5) + rb.lo(Ark::Y);
#if(AMREX_SPACEDIM == 3)
            const double z  = dx[Ark::Z]*(k+0.5) + rb.lo(Ark::Z);
            const double distance = std::sqrt(Ark::Square(x-centerX, y-centerY, z-centerZ)) - RadiusPF();
#elif(AMREX_SPACEDIM == 2)
            const double distance = std::hypot(x-centerX, y-centerY) - RadiusPF();
#endif
            phi(i,j,k) = (PhiV() + PhiL())/2.0 + (PhiV() - PhiL())*tanh(distance*2.0/WidthI())/2.0;

            vel(i,j,k,Ark::Y) = -phi(i,j,k)*U0();

            if(y < (Hw + wI)) {
                phi(i,j,k) = (PhiV() + PhiL())/2.0 + (PhiV() - PhiL())/2.0*tanh(2.*(y-Hw)/WidthI());
                vel(i,j,k,Ark::Y) = 0;
            }
            
            rho(i,j,k) = APhi()*phi(i,j,k) + BPhi();
        });
    }
}

void DUGKSAmr::TaylorBubble(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 230 && CaseName() == "TaylorBubble"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();
    const Box &domain = Geom(lev).Domain();

    ld.msv[Ark::VEL].setVal(0);
    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    MultiFab &PhiFab = ld.msv[Ark::PHI];
    MultiFab &RhoFab = ld.msv[Ark::RHO];
    // MultiFab &VelFab = ld.msv[Ark::VEL];

    for(amrex::MFIter mfi(PhiFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        Array4<Real> const &phi = PhiFab.array(mfi);
        Array4<Real> const &rho = RhoFab.array(mfi);
        // Array4<Real> const &vel = VelFab.array(mfi);

        const Real centerX = CenterPF(Ark::X);
        const Real centerY = CenterPF(Ark::Y);
#if(AMREX_SPACEDIM == 3)
        const Real centerZ = CenterPF(Ark::Z);
#endif
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x  = dx[Ark::X]*(i+0.5) + rb.lo(Ark::X);
            const double y  = dx[Ark::Y]*(j+0.5) + rb.lo(Ark::Y);
#if(AMREX_SPACEDIM == 3)
            const double z  = dx[Ark::Z]*(k+0.5) + rb.lo(Ark::Z);
            const double distance = Ark::Hypot(x-centerX, y-centerY, z-centerZ) - RadiusPF();
#elif(AMREX_SPACEDIM == 2)
            const double distance = std::hypot(x-centerX, y-centerY) - RadiusPF();
#endif
            phi(i,j,k) = PhaseProfile(distance, LV::bubble);

            const double tailY = RadiusPF();

            if(y < centerY + rb.lo(Ark::Y))
            {
#if(AMREX_SPACEDIM == 3)
                const double distanceCylinder = std::hypot(x-centerX,z-centerZ) - RadiusPF();
#elif(AMREX_SPACEDIM == 2)
                const double distanceCylinder = std::hypot(x-centerX,0) - RadiusPF();
#endif
                const double distanceTail = tailY - y;

                double distanceMin = 0;

                if(distanceCylinder <= 0 && distanceTail <= 0) 
                {
                    distanceMin = 
                    distanceCylinder > distanceTail ? distanceCylinder : distanceTail;
                }
                else if(distanceCylinder > 0 && distanceTail <= 0)
                {
                    distanceMin = distanceCylinder;
                }
                else if(distanceCylinder <= 0 && distanceTail > 0)
                {
                    distanceMin = distanceTail;
                }
                else
                {
                    distanceMin = std::hypot(distanceCylinder, distanceTail);
                }

                phi(i,j,k) = PhaseProfile(distanceMin, LV::bubble);
            }
        });
    }
}


void DUGKSAmr::RisingBubble(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 231 && CaseName() == "RisingBubble"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();
    const Box &domain = Geom(lev).Domain();

    ld.msv[Ark::VEL].setVal(0);
    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    MultiFab &PhiFab = ld.msv[Ark::PHI];
    MultiFab &RhoFab = ld.msv[Ark::RHO];
    // MultiFab &VelFab = ld.msv[Ark::VEL];

    for(amrex::MFIter mfi(PhiFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        Array4<Real> const &phi = PhiFab.array(mfi);
        Array4<Real> const &rho = RhoFab.array(mfi);
        // Array4<Real> const &vel = VelFab.array(mfi);

        const Real centerX = CenterPF(Ark::X);
        const Real centerY = CenterPF(Ark::Y);
#if(AMREX_SPACEDIM == 3)
        const Real centerZ = CenterPF(Ark::Z);
#endif
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x  = dx[Ark::X]*(i+0.5) + rb.lo(Ark::X);
            const double y  = dx[Ark::Y]*(j+0.5) + rb.lo(Ark::Y);
#if(AMREX_SPACEDIM == 3)
            const double z  = dx[Ark::Z]*(k+0.5) + rb.lo(Ark::Z);
            const double distance = Ark::Hypot(x-centerX, y-centerY, z-centerZ) - RadiusPF();
#elif(AMREX_SPACEDIM == 2)
            const double distance = std::hypot(x-centerX, y-centerY) - RadiusPF();
#endif
            phi(i,j,k) = PhaseProfile(distance, LV::bubble);
        });
    }
}

void DUGKSAmr::FallingDrop(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 232 && CaseName() == "FallingDrop"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    // const RealBox &rb = Geom(lev).ProbDomain();
    const Real* prob_lo = Geom(lev).ProbDomain().lo();
    const Box &domain = Geom(lev).Domain();

    ld.msv[Ark::VEL].setVal(0);
    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    MultiFab &PhiFab = ld.msv[Ark::PHI];
    MultiFab &RhoFab = ld.msv[Ark::RHO];
    // MultiFab &VelFab = ld.msv[Ark::VEL];

    for(amrex::MFIter mfi(PhiFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        Array4<Real> const &phi = PhiFab.array(mfi);
        Array4<Real> const &rho = RhoFab.array(mfi);
        // Array4<Real> const &vel = VelFab.array(mfi);

        const Real centerX = CenterPF(Ark::X);
        const Real centerY = CenterPF(Ark::Y);
#if(AMREX_SPACEDIM == 3)
        const Real centerZ = CenterPF(Ark::Z);
#endif
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x  = dx[Ark::X]*(i+0.5) + prob_lo[Ark::X];
            const double y  = dx[Ark::Y]*(j+0.5) + prob_lo[Ark::Y];
#if(AMREX_SPACEDIM == 3)
            const double z  = dx[Ark::Z]*(k+0.5) + prob_lo[Ark::Z];
            const double distance = Ark::Hypot(x-centerX, y-centerY, z-centerZ) - RadiusPF();
#elif(AMREX_SPACEDIM == 2)
            const double distance = std::hypot(x-centerX, y-centerY) - RadiusPF();
#endif
            phi(i,j,k) = PhaseProfile(distance, LV::drop);
        });
    }
}