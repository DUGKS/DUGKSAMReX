#include <DUGKSAmr.H>

//
// Print maximum values (useful for tracking evolution)
//
void DUGKSAmr::PrintMaxValues(Real time_in)
{
    return; // xxxxx TODO

#if 0
    ComputeDivU(time_in);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        amrex::Print() << "Level " << lev << std::endl;
        PrintMaxVel(lev);
        PrintMaxGp(lev);
    }
    amrex::Print() << std::endl;
#endif
}

//
// Print the maximum values of the velocity components and velocity divergence
//
void DUGKSAmr::PrintMaxVel(int lev)
{
#if 0
    // xxxxx TODO
    amrex::Print() << "max(abs(u/v/w/divu))  = "
                   << Norm(vel , lev, 0, 0) << "  "
		   << Norm(vel , lev, 1, 0) << "  "
                   << Norm(vel , lev, 2, 0) << "  "
                   << Norm(divu, lev, 0, 0) << "  " << std::endl;
    if (ntrac > 0)
    {
       for (int i = 0; i < ntrac; i++)
          amrex::Print() << "max tracer" << i << " = " << Norm(tracer,lev,i,0) << std::endl;;
    }
#endif
}

//
// Print the maximum values of the pressure gradient components and pressure
//
void DUGKSAmr::PrintMaxGp(int lev)
{
#if 0
    // xxxxx TODO
    amrex::Print() << "max(abs(gpx/gpy/gpz/p))  = "
                   << Norm(gp, lev, 0, 0) << "  "
		   << Norm(gp, lev, 1, 0) << "  "
                   << Norm(gp, lev, 2, 0) << "  "
		   << Norm(p , lev, 0, 0) << "  " << std::endl;
#endif
}

void DUGKSAmr::CheckForNans(int lev)
{
#if 0
    // xxxxx
    bool ro_has_nans = density[lev]->contains_nan(0);
    bool ug_has_nans = vel[lev]->contains_nan(0);
    bool vg_has_nans = vel[lev]->contains_nan(1);
    bool wg_has_nans = vel[lev]->contains_nan(2);
    bool pg_has_nans = p[lev]->contains_nan(0);

    if (ro_has_nans)
	amrex::Print() << "WARNING: ro contains NaNs!!!";

    if (ug_has_nans)
	amrex::Print() << "WARNING: u contains NaNs!!!";

    if(vg_has_nans)
	amrex::Print() << "WARNING: v contains NaNs!!!";

    if(wg_has_nans)
	amrex::Print() << "WARNING: w contains NaNs!!!";

    if(pg_has_nans)
	amrex::Print() << "WARNING: p contains NaNs!!!";
#endif
}
void DUGKSAmr::CheckDUGKSAmr()
{
	amrex::Print()<<"Checking variables in DUGKSAmr ..."<<std::endl;

//
    if(Ark::DVM_MOMENTUM < Ark::DVM_PHASE or Ark::DVM_MOMENTUM < Ark::DVM_ENERGY)
    {
        std::ostringstream ss;   
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  << "DVM_MOMENTUM is expected to be the maximum," 
            << " otherwise the boundary conditions for DVDF would be ill-posed"
            << '\n';
        amrex::Abort(ss.str());
    }

    if(steady_state && end_step > 0) 
    {
        std::ostringstream ss;   
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  << "end_step is expected to be negative if steady_state is used!" 
            << '\n';
        amrex::Abort(ss.str());
    }
    if(!steady_state && end_step <= 0)
    {
        std::ostringstream ss;   
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  << "either steady_state or end_step should be used!" 
            << '\n';
        amrex::Abort(ss.str());
    }
    if( Update_DVDF_Barh != &DUGKSAmr::Update_DVDF_Barh_UW 
        and 
        Update_DVDF_Barh != &DUGKSAmr::Update_DVDF_Barh_CD)
    {
        std::ostringstream ss;   
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  << "Invalid function for Update_DVDF_Barh!!!" 
            << '\n';
        amrex::Abort(ss.str());
    }
	amrex::Print()<<"Checking Done."<<std::endl;
}
void DUGKSAmr::Verify(int const lev)
{
	auto &ld = *m_leveldata[lev];

	const Real *dx = Geom(lev).CellSize();

	for(amrex::MFIter mfi(ld.msv[Ark::PHI]);mfi.isValid();++mfi)
    {
        // const Box& bx = mfi.validbox();
        const Box& bx = mfi.growntilebox(1);
        // const Box& bx = mfi.fabbox();
        //
        Array4<Real> const &phi  = ld.msv[Ark::PHI].array(mfi);
        Array4<Real> const &rho  = ld.msv[Ark::RHO].array(mfi);
        Array4<Real> const &vel  = ld.msv[Ark::VEL].array(mfi);
        Array4<Real> const &pres = ld.msv[Ark::P].array(mfi);
        Array4<Real> const &temp = ld.msv[Ark::T].array(mfi);
        Array4<Real> const &grad_phi = ld.msv[Ark::GradPHI].array(mfi);

        // Array4<Real> const &phi_old  = ld.phase_old.array(mfi);
        // Array4<Real> const &rho_old  = ld.density_old.array(mfi);
        // Array4<Real> const &vel_old  = ld.velocity_old.array(mfi);

        Array4<Real> const &hWave = ld.hDVDF[Ark::WAVE].array(mfi);
        Array4<Real> const &hBP = ld.hDVDF[Ark::BP].array(mfi);

        Array4<Real> const &fBP = ld.fDVDF[Ark::BP].array(mfi);
        Array4<Real> const &fWave = ld.fDVDF[Ark::WAVE].array(mfi);

        int const Num = 1;//Ark::DVM_MOMENTUM;
        // 

        std::shared_ptr<int> nn = std::make_shared<int>(0);

        ParallelFor(bx,Num,[=](int i,int j,int k,int n) noexcept
        {
            const double x = dx[Ark::dir_x]*(i+0.5);
            const double y = dx[Ark::dir_y]*(j+0.5);

            Real value_real = 2*x + x*y;
            Real value_average = rho(i,j,k);

            if(fabs(value_average - value_real) > 1E-12)
            {
            	// if(not (i == -1 || i == 16 || j == -1 || j == 16 ) )
                // if((j == 31 or j == 0) or fabs(value_average - value_real) > 1E-12)
                {
                amrex::AllPrint().SetPrecision(17)<<"i : "<<i<<"    "<<"j : "<<j<<"    "
                                      <<"phi : "<<value_average<<"    "<<"Real phi : "
                                      <<value_real<<" nn : "<<++(*nn)<<std::endl;
                }
            }
        });
    }

    // for(amrex::MFIter mfi(ld.msv_yF[Ark::RHO], true);mfi.isValid();++mfi)
    // {
    //     const Box& tbx = mfi.tilebox();
    //     //
    //     Array4<Real> const &rho  = ld.msv_yF[Ark::RHO].array(mfi);
    //     Array4<Real> const &pres = ld.msv_yF[Ark::P].array(mfi);

    //     std::shared_ptr<int> nn = std::make_shared<int>(0);
    //     amrex::ParallelFor(tbx, 1,
    //     [=](int i,int j,int k,int n) noexcept
    //     {
    //         const Real x = dx[Ark::dir_x]*(i+0.5);
    //         const Real y = dx[Ark::dir_y]*(j);

    //         Real value_real = 2*x + 3*y;
    //         // Real value_real = (2.0 * x + 4.0*y);
    //         Real value_average = pres(i,j,k);
    //         // if(fabs(value_average - value_real) < 1E-14)
    //         if(value_average != 1)
    //         {
    //         	amrex::Print().SetPrecision(17)<<"i : "<<i<<"    "<<"j : "<<j<<"    "
    //                                   <<"phi : "<<value_average<<"    "<<"Real phi : "
    //                                   <<value_real<<" nn : "<<++(*nn)<<std::endl;
    //         }
    //     });
    // }
    //
    if (m_verbose > 0 and amrex::ParallelDescriptor::IOProcessor()) {
        printGridSummary(amrex::OutStream(), 0, finest_level);
    }
}
void DUGKSAmr::SetUnit(int const lev)
{
    auto &ld = *m_leveldata[lev];

    for(auto &mf : ld.msv) {
        mf.setVal(1);
    }

    for(auto &mf : ld.msv_xF) {
        mf.setVal(1);
    }

    for(auto &mf : ld.msv_yF) {
        mf.setVal(1);
    }

    for(auto &hdvdf : ld.hDVDF) {
        hdvdf.setVal(1);
    }

    for(auto &fdvdf : ld.fDVDF) {
        fdvdf.setVal(1);
    }

    ld.hDVDFxF[Ark::FLUX].setVal(1);
    ld.hDVDFyF[Ark::FLUX].setVal(1);

    ld.fDVDFxF[Ark::FLUX].setVal(1);
    ld.fDVDFyF[Ark::FLUX].setVal(1);

    // ld.velocity_old.setVal(3,0,1);
    // ld.velocity_old.setVal(5,1,1);
    // ld.phase_old.setVal(9,0,1);
}
void DUGKSAmr::Test()
{
    int const lev = 0;
    auto &ld = Data(lev);
    SetUnit(lev);
    // ld.msv[Ark::RHO].FillBoundary(0, 1, IntVect(2), Geom(lev).periodicity());
    // FillPatchNew(lev,cur_time,ld);
    // FillCoarsePatchDensity(lev,cur_time,ld.density_old);
    // ld.density_old.FillBoundary(0,1,IntVect(1), Geom(lev).periodicity());
    // amrex::Print()<<ld.msv[Ark::PHI].norm2();
    L2NormPhase(0);
    // Verify(lev);
    WritePlotFile();
}
// Real DUGKSAmr::L2NormDiff(const MultiFab& xFabA, const MultiFab& xFabB, int const comp)
// {
//     std::unique_ptr<MultiFab> diff_xFab;
// #ifdef AMREX_USE_EB
//     (*diff_xFab).reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, 
//                         MFInfo(), EBFactory(lev)));
// #else
//     (*diff_xFab).reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, 
//                         MFInfo(), Factory(lev)));
// #endif
//     MultiFab::LinComb(*diff_vel[lev], 
//                       1.0, xFabA, 0, 
//                       -1.0, xFabB, 0, 
//                       0, AMREX_SPACEDIM, 0);

//     amrex::Print()<<(*diff_vel[lev]).norm2(0) << std::endl;
//     amrex::Print()<<(*diff_vel[lev]).norm2(1) << std::endl;

//     return (*diff_vel[lev]).norm2();
// }

