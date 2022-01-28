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
void DUGKSAmr::CheckForPhysicalBoundary (int const lev)
{
    auto &ld = Data(lev);
    const Real *dx = Geom(lev).CellSize();
    {

    const Real real_lo_x = Geom(lev).ProbDomain().lo(Ark::X);
    const Real real_hi_x = Geom(lev).ProbDomain().hi(Ark::X);
    const auto &msvFaceX = ld.msv_xF;
    for(amrex::MFIter mfi(msvFaceX[Ark::VEL], true);mfi.isValid();++mfi)
    {
        const Box &bx = mfi.tilebox();

        const Array4<const Real> &vel = msvFaceX[Ark::VEL].const_array(mfi);

        amrex::ParallelFor(bx, [=] (int i,int j,int k)
        {
            const double x = real_lo_x + i*dx[Ark::X];

            if(x == real_lo_x or x == real_hi_x)
            {
                // if(vel(i,j,k,Ark::X) != 0)
                // {
                //     amrex::AllPrint()   << "x : "<<x<<"    "
                //                         << "i : "<<i<<"    "
                //                         << "j : "<<j<<"    "
                //                         << "vel_x : "<<vel(i,j,k,Ark::X)<<std::endl;
                // }

                // if(vel(i,j,k,Ark::Y) != 0)
                // {
                //     amrex::AllPrint()   << "x : "<<x<<"    "
                //                         << "i : "<<i<<"    "
                //                         << "j : "<<j<<"    "
                //                         << "vel_y : "<<vel(i,j,k,Ark::Y)<<std::endl;
                // }
            }
        });
    }
    }
    //
    {
        const auto &msvFaceY = ld.msv_yF;
        const Real real_lo_y = Geom(lev).ProbDomain().lo(Ark::Y);
        const Real real_hi_y = Geom(lev).ProbDomain().hi(Ark::Y);
        if(real_hi_y != 100 or real_lo_y != -100) 
            amrex::Abort("if(real_hi_y != 100 or real_lo_y != -100) amrex::Abort();");
        for(amrex::MFIter mfi(msvFaceY[Ark::VEL], true);mfi.isValid();++mfi)
        {
            const Box &bx = mfi.tilebox();

            const Array4<const Real> &vel = msvFaceY[Ark::VEL].const_array(mfi);

            amrex::ParallelFor(bx, [=] (int i,int j,int k)
            {
                const double y = real_lo_y + j*dx[Ark::Y];

                if(y == real_lo_y or y == real_hi_y)
                {
                    if(vel(i,j,k,Ark::X) != 0)
                    {
                        amrex::AllPrint()   << "y : "<<y<<"    "
                                            << "i : "<<i<<"    "
                                            << "j : "<<j<<"    "
                                            << "vel_x : "<<vel(i,j,k,Ark::X)<<std::endl;
                    }

                    if(vel(i,j,k,Ark::Y) != 0)
                    {
                        amrex::AllPrint()   << "y : "<<y<<"    "
                                            << "i : "<<i<<"    "
                                            << "j : "<<j<<"    "
                                            << "vel_y : "<<vel(i,j,k,Ark::Y)<<std::endl;
                    }
                }
            });
        }
    }
    //
    // {
    //     amrex::Print()<<"--------------"<<nstep<<"--------------\n";
    //     const auto &VelFab = ld.msv[Ark::VEL];

    //     const int lo_min = Geom(lev).Domain().smallEnd(Ark::Y);
    //     const int hi_max = Geom(lev).Domain().bigEnd(Ark::Y);
    //     amrex::Print()<<lo_min<<"  "<<hi_max<<std::endl;
    //     for(amrex::MFIter mfi(VelFab, true);mfi.isValid();++mfi)
    //     {
    //         const Box &bx = mfi.tilebox();

    //         auto const lo = lbound(bx);
    //         auto const hi = ubound(bx);

    //         const Array4<const Real> &vel = VelFab.const_array(mfi);

    //         amrex::ParallelFor(bx, [=] (int i,int j,int k)
    //         {
    //             if(j == lo_min) 
    //             {
    //                 amrex::Print().SetPrecision(18)
    //                 <<"i "<<i<<"  "
    //                 <<"j "<<j<<"  " 
    //                 <<"lo_min "<<lo_min<<"  "
    //                 <<"vel_y : "<< vel(i,j,k,Ark::Y) <<std::endl;
    //             }
    //             if(j == hi_max)
    //             {
    //                 amrex::Print().SetPrecision(18)
    //                 <<"i "<<i<<"  "
    //                 <<"j "<<j<<"  " 
    //                 <<"hi_max "<<hi_max<<"  "
    //                 <<"vel_y : "<< vel(i,j,k,Ark::Y) <<std::endl;
    //             }
    //         });
    //     }
    // }
    {
        const auto &msvFaceY = ld.msv_yF;
        const MultiFab& fFluxFabyF = ld.fDVDFyF[Ark::FLUX];
        const MultiFab& PhiFab  = ld.msv[Ark::PHI];
        const MultiFab& PresFab = ld.msv[Ark::P];

        const Real real_lo_y = Geom(lev).ProbDomain().lo(Ark::Y);
        const Real real_hi_y = Geom(lev).ProbDomain().hi(Ark::Y);
        if(real_hi_y != 100 or real_lo_y != -100) amrex::Abort(__func__);

        for(amrex::MFIter mfi(msvFaceY[Ark::P], true);mfi.isValid();++mfi)
        {
            const Box &bx = mfi.tilebox();

            const Array4<const Real> &phi  = PhiFab.const_array(mfi);
            const Array4<const Real> &phi_yF  = msvFaceY[Ark::PHI].const_array(mfi);
            const Array4<const Real> &pres = PresFab.const_array(mfi);
            const Array4<const Real> &pres_yF = msvFaceY[Ark::P].const_array(mfi);
            const Array4<const Real> &fFlux_yF = fFluxFabyF.const_array(mfi);

            amrex::ParallelFor(bx, [=] (int i,int j,int k)
            {
                const double y = real_lo_y + j*dx[Ark::Y];
                // amrex::Print()  << get_dirichlet_pressure()[0]<<"    "
                //                 << get_dirichlet_pressure()[1]<<"    "
                //                 << get_dirichlet_pressure()[2]<<"    "
                //                 << get_dirichlet_pressure()[3]<<"    \n";
                if(y == real_lo_y)
                {
                    if(pres(i,j,k) != pres_yF(i,j,k) or pres(i,j-1,k) != pres_yF(i,j,k))
                    // if(pres_yF(i,j,k) != 1)
                    {
                        amrex::Print().SetPrecision(18)
                        <<"i "<<i<<"  "
                        <<"j "<<j<<"  "
                        <<"cell : "<<pres(i,j,k)<<"  "
                        <<"face : "<<pres_yF(i,j,k)<<std::endl;
                    }
                    if(phi_yF(i,j,k) != phi(i,j,k) or phi_yF(i,j,k) != phi(i,j-1,k) )
                    {
                        amrex::Print().SetPrecision(18)
                        <<"i "<<i<<"  "
                        <<"j "<<j<<"  "
                        <<"cellphi : "<<phi(i,j,k)<<"  "
                        <<"facephi : "<<phi_yF(i,j,k)<<std::endl;
                    }
                    // for(int n = 0;n < 9;++n){
                    // if(
                    //     (fabs(fFlux_yF(i,j,k,n) - Ark::DmVnMomentum::omega[n]*1*3.0)) 
                    //     > 1E-16
                    //   )
                    // {
                    //     amrex::Print().SetPrecision(18)
                    //     <<"i "<<i<<"  "
                    //     <<"j "<<j<<"  "
                    //     <<"n "<<n<<"  "
                    //     <<"flux : "<<fFlux_yF(i,j,k,n)<<"  "
                    //     <<"true : "<<Ark::DmVnMomentum::omega[n]*3.0<<std::endl;
                    // }}                 
                }

                if(y == real_hi_y)
                {
                    if(pres(i,j-1,k) != pres_yF(i,j,k) or pres(i,j,k) != pres_yF(i,j,k) )
                    // if(pres_yF(i,j,k) != 1)
                    {
                        amrex::Print().SetPrecision(18)
                        <<"i "<<i<<"  "
                        <<"j "<<j<<"  "
                        <<"cell : "<<pres(i,j,k)<<"  "
                        <<"face : "<<pres_yF(i,j,k)<<std::endl;
                    }
                    if(phi_yF(i,j,k) != phi(i,j,k) or phi_yF(i,j,k) != phi(i,j-1,k) )
                    {
                        amrex::Print().SetPrecision(18)
                        <<"i "<<i<<"  "
                        <<"j "<<j<<"  "
                        <<"cellphi : "<<phi(i,j,k)<<"  "
                        <<"facephi : "<<phi_yF(i,j,k)<<std::endl;
                    }
                    // for(int n = 0;n < 9;++n){
                    // if(
                    //     (fabs(fFlux_yF(i,j,k,n) - Ark::DmVnMomentum::omega[n]*1*3.0)) 
                    //     > 1E-16
                    //   )
                    // {
                    //     amrex::Print().SetPrecision(18)
                    //     <<"i "<<i<<"  "
                    //     <<"j "<<j<<"  "
                    //     <<"n "<<n<<"  "
                    //     <<"flux : "<<fFlux_yF(i,j,k,n)<<"  "
                    //     <<"true : "<<Ark::DmVnMomentum::omega[n]*3.0<<std::endl;
                    // }}
                }
            });
        }
    }
}
void DUGKSAmr::CheckDUGKSAmr()
{
	amrex::Print()<<"Checking variables in DUGKSAmr ..."<<std::endl;

#ifdef ARK_USE_PHASE
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM == Ark::DmVnPhase::D,
        "Fatal Error : mismatch between DVM_PHASE and DmVnPhase");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ark::approximatelyEqual(MPhi(), TauPhase()*Ark::DmVnPhase::RT),
    "Fatal Error : mismatch between mobility_phase and tau_phase");
#endif

#ifdef ARK_USE_MOMENTUM
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM == Ark::DmVnMomentum::D,
        "Fatal Error : mismatch between DVM_MOMENTUM and DmVnMomentum");
#endif

#ifdef ARK_USE_ENERGY
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM == Ark::DmVnEnergy::D,
        "Fatal Error : mismatch between DVM_ENERGY and DmVnEnergy");
#endif

    if( DmVnModel != DmVnName())
    {
        std::ostringstream ss;   
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  << "wrong discrete velocity model!!!" 
            << '\n';
        amrex::Abort(ss.str());
    }

    if (Geom(0).Domain().length(Ark::X) != ChNum() 
        and 
        Geom(0).Domain().length(Ark::Y) != ChNum()
#if(AMREX_SPACEDIM == 3)
        and 
        Geom(0).Domain().length(Ark::Z) != ChNum()
#endif
        )
    {
        std::ostringstream ss;   
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  << "incosistence between ch_num:" << ChNum() 
            << " and ncell_x:" << Geom(0).Domain().length(Ark::X) 
            << '\n';
        amrex::Abort(ss.str());
    }    
    if (static_cast<int>(ChLength()/MinL()) != ChNum())
    {
        std::ostringstream ss;   
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  << "incosistence among ch_length:" << ChLength() 
            << " ch_num:" << ChNum()  
            <<" and min_length:" << MinL()
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

    {
        if (!Ark::tinyEqual(CFL()*MinL()/std::pow(2,max_level), Dt(), 1E-10))
        {
           std::ostringstream ss;   
           ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
           ss  << "incosistence between time step and CFL number!!!" 
               << '\n';
           amrex::Abort(ss.str()); 
        }
    }

	amrex::Print()<<"Checking Done."<<std::endl;
}
//
#ifdef AMREX_DEBUG

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
    // SetUnit(lev);
    // ld.msv[Ark::RHO].FillBoundary(0, 1, IntVect(2), Geom(lev).periodicity());
    FillPatchNew(lev,cur_time,ld);
    // FillCoarsePatchDensity(lev,cur_time,ld.density_old);
    // ld.density_old.FillBoundary(0,1,IntVect(1), Geom(lev).periodicity());
    // amrex::Print()<<ld.msv[Ark::PHI].norm2();
    // L2NormPhase(0);
    Verify(lev);
    WritePlotFile();
}

void DUGKSAmr::Verify(int const lev)
{
    LevelData &ld = Data(lev);
    const Box &domain_box = Geom(lev).Domain();
    const Real *dx = Geom(lev).CellSize();
    const RealBox &rb = Geom(lev).ProbDomain();

    const auto& domain_hi = domain_box.hiVect();
    const auto& domain_lo = domain_box.loVect();

    const int ilo = domain_lo[0];
    const int ihi = domain_hi[0];

    for(amrex::MFIter mfi(ld.msv[Ark::PHI]);mfi.isValid();++mfi)
    {
        // const Box& bx = mfi.validbox();
        const Box& bx = mfi.growntilebox(3);
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
            const double x = dx[Ark::X]*(i+0.5) + rb.lo(Ark::X);
            const double y = dx[Ark::Y]*(j+0.5) + rb.lo(Ark::Y);
#if(AMREX_SPACEDIM == 3)
            const double z = dx[Ark::Z]*(k+0.5) + rb.lo(Ark::Z);
#else
            const double z = 0;
#endif

            Real value_real = 2*x + x*y + y*z;
            Real value_average = rho(i,j,k);

            if(fabs(value_average - value_real) > 1E-12)
            {
                // if(not (i == -1 || i == 16 || j == -1 || j == 16 ) )
                // if((j == 31 or j == 0) or fabs(value_average - value_real) > 1E-12)
                {
                amrex::AllPrint().SetPrecision(17)<<"i : "<<i<<"    "<<"j : "<<j<<"    "
                                    <<"k : "<<k<<"    "
                                    <<"phi : "<<value_average<<"    "
                                    <<"Real phi : "<<value_real
                                    <<" nn : "<<++(*nn)<<std::endl;
                }
            }
            
            // if(i < ilo) {
            // if(rho(i,j,k) != rho(ihi+(i-ilo+1),j,k))
            // {
            // amrex::AllPrint().SetPrecision(17)<<"i : "<<i<<"    "<<"j : "<<j<<"    "
            //                     <<"k : "<<k<<"    "
            //                     <<"phi : "<<rho(i,j,k)<<"    "
            //                     <<"Real phi : "<<rho(ihi-(i-ilo+1),j,k)
            //                     <<" nn : "<<++(*nn)<<std::endl;
            // }
            // }
            // //
            // if(i > ihi) {
            // if(rho(i,j,k) != rho(ilo+(i-ihi-1),j,k))
            // {
            // amrex::AllPrint().SetPrecision(17)<<"i : "<<i<<"    "<<"j : "<<j<<"    "
            //                     <<"k : "<<k<<"    "
            //                     <<"phi : "<<rho(i,j,k)<<"    "
            //                     <<"Real phi : "<<rho(ilo+(i-ihi-1),j,k)
            //                     <<" nn : "<<++(*nn)<<std::endl;
            // }
            // }
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
    //         const Real x = dx[Ark::X]*(i+0.5);
    //         const Real y = dx[Ark::Y]*(j);

    //         Real value_real = 2*x + 3*y;
    //         // Real value_real = (2.0 * x + 4.0*y);
    //         Real value_average = pres(i,j,k);
    //         // if(fabs(value_average - value_real) < 1E-14)
    //         if(value_average != 1)
    //         {
    //          amrex::Print().SetPrecision(17)<<"i : "<<i<<"    "<<"j : "<<j<<"    "
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
#endif
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

