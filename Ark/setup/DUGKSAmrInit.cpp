#include <DUGKSAmr.H>

// Resize all arrays when instance of DUGKSAmr class is constructed.
// This is only done at the very start of the simulation.
void DUGKSAmr::ResizeArrays ()
{
    // Time holders for fillpatch stuff
    m_t_new.resize(max_level + 1, 0.);
    m_t_old.resize(max_level + 1, 0.);

    trace_time_step.resize(max_level + 1, 0.);

    m_leveldata.resize(max_level+1);
    m_factory.resize(max_level+1);

    m_dirichlet_velocity_x.fill(0);
    m_dirichlet_velocity_y.fill(0);
    m_dirichlet_velocity_z.fill(0);
    m_dirichlet_pressure.fill(0);
    m_dirichlet_temperature.fill(0);
}
void DUGKSAmr::ReadParametersDUGKSAmr()
{
    ReadIOParameters();
    // ReadRheologyParameters();

    { // Prefix amr
 	ParmParse pp("amr");

	pp.query("regrid_int", m_regrid_int);
#ifdef AMREX_USE_EB
        pp.query("refine_cutcells", m_refine_cutcells);
#endif
    } // end prefix amr

    {
    amrex::ParmParse ppTime("time");
    ppTime.get("constant_trace", use_constant_trace);
    ppTime.get("CFL", Courant_Friedrichs_Lewy);
    ppTime.get("advance_time_step", advance_time_step);

    double hTime = 0.5*advance_time_step;
    ppTime.query("trace_time_step", hTime);

    if(use_constant_trace) {
        for(auto &d : trace_time_step) { d = hTime; }
    }
    else {
        SetMutableTraceTimeStep(hTime);
    }

    }

    { // Prefix DUGKSAmr
	ParmParse pp("DUGKSAmr");

    pp.get("verbose", m_verbose);
    pp.get("n_ghost", n_ghost);

    pp.queryarr("phi_error", phi_error);
    // if (phi_error.size() > 0) {
    //     Real last = phi_error.back();
    //     phi_error.resize(max_level+1, last);
    // }

    pp.queryarr("grad_phi_error", grad_phi_error);
    // if (grad_phi_error.size() > 0) {
    //     Real last = grad_phi_error.back();
    //     grad_phi_error.resize(max_level+1, last);
    // }
    pp.queryarr("vel_error", vel_error);

    tag_region_lo.resize(3);
    tag_region_hi.resize(3);

    pp.queryarr("tag_region",tag_region);
    pp.queryarr("tag_region_lo", tag_region_lo);
    pp.queryarr("tag_region_hi", tag_region_hi);
    }// end prefix DUGKSamr
}

void DUGKSAmr::InitData()
{
    BL_PROFILE("DUGKSAmr::InitData()");

    if (name_restart_file.empty())
    {
        // This tells the AmrMesh class not to iterate when creating the initial
        // grid hierarchy
        // SetIterateToFalse();

        // This tells the Cluster routine to use the new chopping routine which
        // rejects cuts if they don't improve the efficiency
        SetUseNewChop();

        // This is an AmrCore member function which recursively makes new levels
        // with MakeNewLevelFromScratch.
        InitFromScratch(cur_time);

        // AverageDownAll();

        // Set m_nstep to 0 before entering time loop
        nstep = 0;

        // this is rather stupid.
        if (write_check_interval > 0) { WriteCheckPointFile(); }

        // Plot initial distribution
        if (write_plot_interval > 0)  { WritePlotFile();       }
    }
    else
    {
        // Read starting configuration from chk file.
        ReadCheckpointFile();

        if (write_plot_on_restart)
        {
            WritePlotFile();
        }
    }

#ifdef AMREX_USE_EB
    ParmParse pp("DUGKSAmr");
    bool write_eb_surface = false;
    pp.query("write_eb_surface", write_eb_surface);
    if (write_eb_surface) WriteMyEBSurface();
#endif

    if (m_verbose > 0 and amrex::ParallelDescriptor::IOProcessor()) {
        printGridSummary(amrex::OutStream(), 0, finest_level);
    }
}

void DUGKSAmr::ReadIOParameters()
{
    {// Prefix control
    ParmParse pp("control");
        //control parameters
    pp.get("end_step", end_step);
    pp.get("steady_state", steady_state);
    if(steady_state) { pp.get("residual_tolerance", residual_tolerance); }
    pp.get("adapt_tolerance", adapt_tolerance);

    pp.query("check_file", name_check_file);
    pp.query("check_interval", write_check_interval);

    pp.query("restart_file", name_restart_file);
    pp.query("plot_on_restart", write_plot_on_restart);

    pp.query("plot_file", name_plot_file);
    pp.query("plot_interval" , write_plot_interval);

    pp.query("conv_interval" , write_conv_interval);
    pp.query("residual_interval" , write_residual_interval);
    }

    {// Prefix plot
    ParmParse pp("plot");
    // The plt_ccse_regtest resets the defaults,
    int plot_ccse_regtest = 0;
    pp.get("ccse_regtest", plot_ccse_regtest);

    // Which variables to write to plotfile
    if (plot_ccse_regtest) {
        pp.get("phi"       , m_plt_phi);
        pp.get("rho"       , m_plt_rho);
        pp.get("velx"      , m_plt_velx);
        pp.get("vely"      , m_plt_vely);
        pp.get("velz"      , m_plt_velz);
        pp.get("pres"      , m_plt_pres);
        pp.get("temp"      , m_plt_temp);
        pp.get("gphix"     , m_plt_gphix);
        pp.get("gphiy"     , m_plt_gphiy);
        pp.get("gphiz"     , m_plt_gphiz);
        pp.get("forcex"    , m_plt_forcex);
        pp.get("forcey"    , m_plt_forcey);
        pp.get("forcez"    , m_plt_forcez);
        pp.get("chemip"    , m_plt_chemip);
        pp.get("tau"       , m_plt_tau);
        pp.get("vort"      , m_plt_vort);
        pp.get("strainrate", m_plt_strainrate);
        pp.get("divu"      , m_plt_divu);
        pp.get("volfrac"   , m_plt_volfrac);
        pp.get("hDVDF"     , m_plt_hDVDF);
        pp.get("fDVDF"     , m_plt_fDVDF);
        pp.get("gDVDF"     , m_plt_gDVDF);
    }
    }
}

void DUGKSAmr::SetMutableTraceTimeStep(const double hTime)
{
    trace_time_step[max_level] = hTime;

    for(int lev = max_level-1;lev >= 0;--lev)
    {
        trace_time_step[lev] = ref_ratio[lev].min()*trace_time_step[lev + 1];
    }
    for(int lev = max_level;lev >=0;--lev)
    {
        amrex::Print()  <<"tracing time for level " <<lev<<" : "
                        << trace_time_step[lev] << std::endl;
    }
}
// void DUGKSAmr::SetFluxScheme()
// {
//     if(FluxScheme() == "CD" or FluxScheme() == "cd")
//     {
//         Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_CD;
//     }
//     else if(FluxScheme() == "UW" or FluxScheme() == "uw")
//     {
//         Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW;
//     }
//     else
//     {
//         amrex::Abort("FulxScheme must be either Central Scheme(CD) or Upwind Scheme(UW)");
//     }
// }

void DUGKSAmr::SetFluxScheme()
{
    if(FluxScheme() == "UW")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW_Default;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Default;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 2, "Fatal Error : invalid n_ghost!!!");
    }
    else if(FluxScheme() == "UW-D-VanLeer")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW_Default_Second;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Default_VanLeer;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 2, "Fatal Error : invalid n_ghost!!!");
    }
    else if(FluxScheme() == "UW-D-2")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW_Default_Second;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Default_Second;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 2, "Fatal Error : invalid n_ghost!!!");
    }
    else if(FluxScheme() == "UW-1")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW_First;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Null;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 2, "Fatal Error : invalid n_ghost!!!");
    }
    else if(FluxScheme() == "UW-3")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW_Third;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Third;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 3, "Fatal Error : invalid n_ghost!!!");
    }
    else if (FluxScheme() == "UW-4F")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW_Fourth;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Fourth;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 3, "Fatal Error : invalid n_ghost!!!");
    }
    else if (FluxScheme() == "UW-4C")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW_Fourth_Cross;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Fourth;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 3, "Fatal Error : invalid n_ghost!!!");
    }
    else if(FluxScheme() == "UW-4C-3")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW_Fourth_Cross_Third;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Fourth_Third;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 3, "Fatal Error : invalid n_ghost!!!");
    }
    else if (FluxScheme() == "UW-5")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW_Fifth;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Fifth;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 3, "Fatal Error : invalid n_ghost!!!");
    }
    else if (FluxScheme() == "CD")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_CD;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Null;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 2, "Fatal Error : invalid n_ghost!!!");
    }
    else if (FluxScheme() == "WENO-3")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_WENO_Third;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Null;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 2, "Fatal Error : invalid n_ghost!!!");
    }
    else if (FluxScheme() == "WENO-5")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_WENO_Fifth;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Null;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 3, "Fatal Error : invalid n_ghost!!!");
    }
    else if (FluxScheme() == "FD-WENO-3")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_FD_DVDF_Barh_WENO_Third;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Null;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 2, "Fatal Error : invalid n_ghost!!!");
    }
    else if (FluxScheme() == "FD-WENO-5")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_FD_DVDF_Barh_WENO_Fifth;
        Update_DVDF_Grad_BarPlus = &DUGKSAmr::Update_DVDF_Grad_BarPlus_Null;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_ghost == 3, "Fatal Error : invalid n_ghost!!!");
    }
    else
    {
        std::ostringstream ss;
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  <<"FulxScheme must be one of :\n"
            << "0. UW\n"
            << "1. UW-3\n"
            << "2. UW-4F\n"
            << "3. UW-4C\n"
            << "4. UW-4C-3\n"
            << "5. UW-5\n"
            << "6. CD\n"
            << "7. WENO-3\n"
            << "8. WENO-5\n"
            << "9. FD-WENO-3\n"
            << "10. FD-WENO-5\n";
        amrex::Abort(ss.str());
    }
}

void DUGKSAmr::SetAdvanceScheme()
{
    if(AdvanceScheme() == "CN")
    {
        Advance = &DUGKSAmr::Advance_CN;
    }
    else if(AdvanceScheme() == "NEQ")
    {
        Advance = &DUGKSAmr::Advance_NEQ;
    }
    else if(AdvanceScheme() == "IMEX")
    {
        Advance = &DUGKSAmr::Advance_IMEX;
    }
    else if(AdvanceScheme() == "BE")
    {
        Advance = &DUGKSAmr::Advance_BackwardEuler;
    }
    else if(AdvanceScheme() == "SS")
    {
        Advance = &DUGKSAmr::Advance_SS;
    }
    else if(AdvanceScheme() == "FD")
    {
        Advance = &DUGKSAmr::Advance_FD;
    }
    else if(AdvanceScheme() == "MRTSS") {
        Advance = &DUGKSAmr::Advance_Moment_SS;
    }
    else if(AdvanceScheme() == "DVM")
    {
        Advance = &DUGKSAmr::Advance_DVM;
    }
    else
    {
        std::ostringstream ss;
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  <<"AdvanceScheme must be one of :\n"
            << "0. CN\n"
            << "1. BE\n"
            << "2. IMEX\n"
            << "3. SS\n"
            << "4. FD\n"
            << "5. DVM\n";
        amrex::Abort(ss.str());
    }
}
void DUGKSAmr::SetL2NormCase()
{
    if(CaseIndex() == 1 and CaseName() == "UnsteadyTaylorGreenVortex") {
        MultiLevelL2Norm = &DUGKSAmr::L2NormUnsteadyTaylorGreenVortex;
    }
    else if(CaseIndex() == 222 and CaseName() == "ZalesakDisk")
    {
        MultiLevelL2Norm = &DUGKSAmr::L2NormZalesakDisk;
    }
    else if(CaseIndex() == 221 and CaseName() == "TransDrop" or
            CaseIndex() == 223 and CaseName() == "ElongatedInterface" or
            CaseIndex() == 224 and CaseName() == "SmoothDeformation")
    {
        MultiLevelL2Norm = &DUGKSAmr::L2NormPhase;
    }
    else {
        std::cout << "L2-Norm inactivated" << std::endl;
    }
}