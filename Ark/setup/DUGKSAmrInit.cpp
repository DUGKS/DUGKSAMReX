#include <DUGKSAmr.H>

// Resize all arrays when instance of DUGKSAmr class is constructed.
// This is only done at the very start of the simulation. 
void DUGKSAmr::ResizeArrays ()
{
    // Time holders for fillpatch stuff
    m_t_new.resize(max_level + 1, 0.);
    m_t_old.resize(max_level + 1, 0.);

    trace_time_step.resize(max_level + 1, 0.);;

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
    ppTime.get("fixed_trace", use_fixed_trace);
    ppTime.get("CFL", Courant_Friedrichs_Lewy);
    ppTime.get("advance_time_step", advance_time_step);

    double hTime = 0.5*advance_time_step;
    ppTime.query("trace_time_step", hTime);

    if(use_fixed_trace) {
        for(auto &d : trace_time_step) { d = hTime; }
    }
    else {
        SetMutableTraceTimeStep(hTime);
    }

    }

    { // Prefix DUGKSAmr
	ParmParse pp("DUGKSAmr");
    
    pp.query("verbose", m_verbose);
    pp.query("tag_region",tag_region);

    pp.queryarr("phi_error", phi_error);
    if (phi_error.size() > 0) {
        Real last = phi_error.back();
        phi_error.resize(max_level+1, last);
    }

    pp.queryarr("grad_phi_error", grad_phi_error);
    if (grad_phi_error.size() > 0) {
        Real last = grad_phi_error.back();
        grad_phi_error.resize(max_level+1, last);
    }

    tag_region_lo.resize(3);
    tag_region_hi.resize(3);

    pp.queryarr("tag_region_lo", tag_region_lo);
    pp.queryarr("tag_region_hi", tag_region_hi);
    }// end prefix DUGKSamr
}

void DUGKSAmr::ReadIOParameters()
{
    // Prefix control
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

   // The plt_ccse_regtest resets the defaults,
    //     but we can over-ride those below
    int plt_ccse_regtest = 0;
    pp.query("plt_ccse_regtest", plt_ccse_regtest);

    // Which variables to write to plotfile
    if (plt_ccse_regtest != 0)
    {
        m_plt_phi         =   0;
        m_plt_rho         =   1;
        m_plt_velx        =   1;
        m_plt_vely        =   1;
        m_plt_velz        =   1;
        m_plt_pres        =   1;
        m_plt_temp        =   0;
        m_plt_gphix       =   0;
        m_plt_gphiy       =   0;
        m_plt_gphiz       =   0;
        m_plt_forcex      =   0;
        m_plt_forcey      =   0;
        m_plt_forcez      =   0;
        m_plt_tau         =   0;
        m_plt_vort        =   0;
        m_plt_strainrate  =   0;
        m_plt_divu        =   0;
        m_plt_volfrac     =   0;
    }    
}
void DUGKSAmr::InitData()
{
    BL_PROFILE("DUGKSAmr::InitData()");

    int restart_flag = 0;
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
        restart_flag = 1;
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
void DUGKSAmr::SetMutableTraceTimeStep(const double hTime)
{
    trace_time_step[max_level] = hTime;

    for(int lev = max_level-1;lev >= 0;--lev)
    {
        trace_time_step[lev] = ref_ratio[lev].min()*trace_time_step[lev + 1];
        amrex::Print() << trace_time_step[lev] << std::endl;
    }
}
void DUGKSAmr::SetFluxScheme()
{
    if(FluxScheme() == "CD" or FluxScheme() == "cd")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_CD;
    }
    else if(FluxScheme() == "UW" or FluxScheme() == "uw")
    {
        Update_DVDF_Barh = &DUGKSAmr::Update_DVDF_Barh_UW;
    }
    else
    {
        amrex::Abort("FulxScheme must be either Central Scheme(CD) or Upwind Scheme(UW)");
    }
}