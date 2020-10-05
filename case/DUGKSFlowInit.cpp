#include <DUGKSAmr.H>

void DUGKSAmr::BasicCase(LevelData &ld, const Real *dx)
{
    AMREX_ALWAYS_ASSERT(CaseIndex() == -666 && CaseName() == "BasicCase");
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
            const double x = dx[Ark::dir_x]*(i+0.5);
            const double y = dx[Ark::dir_y]*(j+0.5);
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
            const double x = dx[Ark::dir_x]*(i);
            const double y = dx[Ark::dir_y]*(j+0.5);

            hFluxX(i,j,k,n) = (2.0 * x + (2 + n)*y);
            fFluxX(i,j,k,n) = ((2.0 +n)*x + n*0.1*y);

            });
    }
}
void DUGKSAmr::UniformFlow(LevelData &ld)
{
	AMREX_ALWAYS_ASSERT((CaseIndex() == 0 && CaseName() == "UniformFlow"));

	ld.msv[Ark::PHI].setVal(Rho0());
	ld.msv[Ark::RHO].setVal(Rho0());
	ld.msv[Ark::VEL].setVal(0);
	ld.msv[Ark::P].setVal(P0());
	ld.msv[Ark::T].setVal(T0());
}
void DUGKSAmr::QuiescentDrop(LevelData &ld, const Real *dl)
{
	AMREX_ALWAYS_ASSERT((CaseIndex() == 10 && CaseName() == "QuiescentDrop"));

	FillLevelPhaseAnalytical(ld.msv[Ark::PHI], dl, LV::drop);
	MultiFab::Copy(ld.msv[Ark::RHO], ld.msv[Ark::PHI], 0, 0, 1, zero_ghost());
	ld.msv[Ark::RHO].mult(APhi(),0);
	ld.msv[Ark::P].setVal(P0());
	ld.msv[Ark::T].setVal(T0());
    
    ld.msv[Ark::PHI].FillBoundary(Geom(0).periodicity());

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for(amrex::MFIter mfi(ld.msv[Ark::PHI], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        Array4<Real const> const &phi  = ld.msv[Ark::PHI].array(mfi);
        Array4<Real> const &grad_phi   = ld.msv[Ark::GradPHI].array(mfi);

        GradFunc::gradient_operator_isotropy(bx, dl, phi, grad_phi);
    }
}