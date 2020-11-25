#include <DUGKSAmr.H>

void DUGKSAmr::BasicCase(int const lev)
{
    AMREX_ALWAYS_ASSERT(CaseIndex() == -666 && CaseName() == "BasicCase");
#ifdef AMREX_DEBUG
    auto &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();

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



        ParallelFor(tbx, Ark::DVM_MOMENTUM, [=](int i,int j,int k,int n) noexcept
        {
            const double x = dx[Ark::X]*(i+0.5);
            const double y = dx[Ark::Y]*(j+0.5);
            phi(i,j,k) = 2*x*x + y*y;
            rho(i,j,k) = 2*x + x*y;

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


        ParallelFor(tbx, Ark::DVM_MOMENTUM, [=](int i,int j,int k,int n) noexcept
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
void DUGKSAmr::QuiescentDrop(int const lev)
{
    AMREX_ALWAYS_ASSERT((CaseIndex() == 10 && CaseName() == "QuiescentDrop"));

    LevelData &ld = Data(lev);
    const Real *dl = Geom(lev).CellSize();

	FillLevelPhaseSphere(ld.msv[Ark::PHI], dl, LV::drop);

	MultiFab::Copy(ld.msv[Ark::RHO], ld.msv[Ark::PHI], 0, 0, 1, zero_ghost());
	ld.msv[Ark::RHO].mult(APhi(), 0);
    ld.msv[Ark::RHO].plus(BPhi(), 0);
    ld.msv[Ark::VEL].setVal(0);
	ld.msv[Ark::P].setVal(P0());
	ld.msv[Ark::T].setVal(T0());
}
void DUGKSAmr::TransDrop(int const lev)
{
    AMREX_ALWAYS_ASSERT((CaseIndex() == 11 && CaseName() == "TransDrop"));

    LevelData &ld = Data(lev);
    const Real *dl = Geom(lev).CellSize();

    FillLevelPhaseSphere(ld.msv[Ark::PHI], dl, LV::drop);

    MultiFab::Copy(ld.msv[Ark::RHO], ld.msv[Ark::PHI], 0, 0, 1, zero_ghost());
    ld.msv[Ark::RHO].mult(APhi(), 0);
    ld.msv[Ark::RHO].plus(BPhi(), 0);

    ld.msv[Ark::VEL].setVal(U0());
    ld.msv[Ark::P].setVal(P0());
    ld.msv[Ark::T].setVal(T0());

    ld.msv_xF[Ark::VEL].setVal(U0());
    ld.msv_yF[Ark::VEL].setVal(U0());
#if(AMREX_SPACEDIM == 3)
    ld.msv_zF[Ark::VEL].setVal(U0());
#endif
}
void DUGKSAmr::ZalesakDisk(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 12 && CaseName() == "ZalesakDisk"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();

    FillLevelPhaseZalesakDisk(ld.msv[Ark::PHI], dx);

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
        //
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x = dx[0]*(i+0.5);
            const double y = dx[1]*(j+0.5);

            vel(i,j,k,Ark::X) = -U0()*pi/ChLength()*(y-CenterPF(Ark::Y));
            vel(i,j,k,Ark::Y) =  U0()*pi/ChLength()*(x-CenterPF(Ark::X));
            
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
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x = dx[0]*(i);
            const double y = dx[1]*(j+0.5);

            vel(i,j,k,Ark::X) = -U0()*pi/ChLength()*(y-CenterPF(Ark::Y));
            vel(i,j,k,Ark::Y) =  U0()*pi/ChLength()*(x-CenterPF(Ark::X));
            
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
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i,int j,int k)
            {
                const double x = dx[0]*(i+0.5);
                const double y = dx[1]*(j);

                vel(i,j,k,Ark::X) = -U0()*pi/ChLength()*(y-CenterPF(Ark::Y));
                vel(i,j,k,Ark::Y) =  U0()*pi/ChLength()*(x-CenterPF(Ark::X));            
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
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i,int j,int k)
            {
                const double x = dx[0]*(i+0.5);
                const double y = dx[1]*(j+0.5);

                vel(i,j,k,Ark::X) = -U0()*pi/ChLength()*(y-CenterPF(Ark::Y));
                vel(i,j,k,Ark::Y) =  U0()*pi/ChLength()*(x-CenterPF(Ark::X));         
            });
    }}
#endif
}
void DUGKSAmr::InterfaceElongation(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 13 && CaseName() == "InterfaceElongation"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();

    FillLevelPhaseSphere(ld.msv[Ark::PHI], dx, LV::drop);

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
            const double x = dx[0]*(i+0.5);
            const double y = dx[1]*(j+0.5);

            vel(i,j,k,Ark::X) = U0()*pi*sin(pi*x/ChLength())*cos(pi*y/ChLength());
            vel(i,j,k,Ark::Y) = -U0()*pi*cos(pi*x/ChLength())*sin(pi*y/ChLength());
            
        });
    }}

    {
    const double real_lo_x = rb.lo(Ark::X);
    const double real_hi_x = rb.hi(Ark::X);

    MultiFab &VelFab_xF = ld.msv_xF[Ark::VEL];
    for(amrex::MFIter mfi(VelFab_xF);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &vel = VelFab_xF.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x = dx[0]*(i);
            const double y = dx[1]*(j+0.5);

            vel(i,j,k,Ark::X) =  U0()*pi*sin(pi*x/ChLength())*cos(pi*y/ChLength());
            vel(i,j,k,Ark::Y) = -U0()*pi*cos(pi*x/ChLength())*sin(pi*y/ChLength());

            if(x == real_lo_x or x == real_hi_x) {
                vel(i,j,k,Ark::X) = 0;
                vel(i,j,k,Ark::Y) = 0;
            }
            
        });
    }}

    {
    const double real_lo_y = rb.lo(Ark::Y);
    const double real_hi_y = rb.hi(Ark::Y);

    MultiFab &VelFab_yF = ld.msv_yF[Ark::VEL];
    for(amrex::MFIter mfi(VelFab_yF);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &vel = VelFab_yF.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            const double x = dx[0]*(i+0.5);
            const double y = dx[1]*(j);

            vel(i,j,k,Ark::X) =  U0()*pi*sin(pi*x/ChLength())*cos(pi*y/ChLength());
            vel(i,j,k,Ark::Y) = -U0()*pi*cos(pi*x/ChLength())*sin(pi*y/ChLength());

            if(y == real_lo_y or y == real_hi_y) {
                vel(i,j,k,Ark::X) = 0;
                vel(i,j,k,Ark::Y) = 0;
            }
            
        });    
    }}
}
void DUGKSAmr::SmoothDeformation(int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 14 && CaseName() == "SmoothDeformation"));

    LevelData &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();

    FillLevelPhaseSphere(ld.msv[Ark::PHI], dx, LV::drop);

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
            const double x = dx[0]*(i+0.5);
            const double y = dx[1]*(j+0.5);

            vel(i,j,k,Ark::X) = 
            -U0()*sin(4.0*PI*x/ChLength())*sin(4.0*PI*y/ChLength());
            vel(i,j,k,Ark::Y) = 
            -U0()*cos(4.0*PI*x/ChLength())*cos(4.0*PI*y/ChLength());            
        });
    }}
}
void DUGKSAmr::LayeredPoiseuille (int const lev)
{
    using namespace Ark;

    AMREX_ALWAYS_ASSERT((CaseIndex() == 15 && CaseName() == "LayeredPoiseuille")); 

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