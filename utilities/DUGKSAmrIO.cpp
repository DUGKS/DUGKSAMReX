#include <DUGKSAmr.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_buildInfo.H>

using namespace::amrex;
using std::ofstream;

namespace { const std::string level_prefix{"Level_"}; }

void Ark::PrintSplitLine(const char c) 
{ 
    amrex::Print()<<"----------------------------------"<<c<<std::endl;
}
void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

void DUGKSAmr::WriteHeader (const std::string& name, bool is_checkpoint) const
{
	if(ParallelDescriptor::IOProcessor())
	{
	    std::string HeaderFileName(name + "/Header");
	    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
	    std::ofstream HeaderFile;

	    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

	    HeaderFile.open(HeaderFileName.c_str(),
	                    std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

	    if(!HeaderFile.good()) {
	        amrex::FileOpenFailed(HeaderFileName);
	    }

	    HeaderFile.precision(17);

	    if(is_checkpoint) {
	        HeaderFile << "Checkpoint version: 1\n";
	    } else {
	        HeaderFile << "HyperCLaw-V1.1\n";
	    }

	    HeaderFile << finest_level << "\n";

	    // Time stepping controls
	    HeaderFile << nstep << "\n";
	    HeaderFile << cur_time << "\n";

	    // Geometry
	    for(int i = 0; i < BL_SPACEDIM; ++i) {
	        HeaderFile << Geom(0).ProbLo(i) << ' ';
	    }
	    HeaderFile << '\n';

	    for(int i = 0; i < BL_SPACEDIM; ++i)
	        HeaderFile << Geom(0).ProbHi(i) << ' ';
	    HeaderFile << '\n';

	    // BoxArray
	    for(int lev = 0; lev <= finest_level; ++lev)
	    {
	        boxArray(lev).writeOn(HeaderFile);
	        HeaderFile << '\n';
	    }
	}
}
void DUGKSAmr::WriteCheckPointFile() const
{
    BL_PROFILE("DUGKSAmr::WriteCheckPointFile()");

    VisMF::SetNOutFiles(128);

    const std::string& checkpointname = amrex::Concatenate(name_check_file, nstep);

    amrex::Print() << "\n\t Writing checkpoint " << checkpointname << std::endl;

    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, finest_level + 1, true);

    bool is_checkpoint = true;
    WriteHeader(checkpointname, is_checkpoint);
    WriteJobInfo(checkpointname);

    for(int lev = 0; lev <= finest_level; ++lev)
    {
        VisMF::Write(m_leveldata[lev]->msv[Ark::PHI],
            amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Phase"));

        VisMF::Write(m_leveldata[lev]->msv[Ark::RHO],
            amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Density"));

        VisMF::Write(m_leveldata[lev]->msv[Ark::VEL],
            amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Velocity"));

        VisMF::Write(m_leveldata[lev]->msv[Ark::P],
            amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "P"));

        VisMF::Write(m_leveldata[lev]->msv[Ark::T],
            amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "T"));
    }
}

void DUGKSAmr::ReadCheckpointFile()
{
    BL_PROFILE("DUGKSAmr::ReadCheckpointFile()");

    amrex::Print() << "Restarting from checkpoint " << name_restart_file << std::endl;

    Real prob_lo[BL_SPACEDIM];
    Real prob_hi[BL_SPACEDIM];

    /***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              allocate DUGKAmr memory (DUGKSAmr::AllocateArrays)            *
     *              (by calling MakeNewLevelFromScratch)
     ***************************************************************************/

    std::string File(name_restart_file + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // Start reading from checkpoint file 
    
    // Title line
    std::getline(is, line);

    // Finest level
    is >> finest_level;
    GotoNextLine(is);

    // Step count
    is >> nstep;
    GotoNextLine(is);

    // Current time
    is >> cur_time;
    GotoNextLine(is);

    // // Time step size
    // is >> m_dt;
    // GotoNextLine(is);

    // is >> m_prev_dt;
    // GotoNextLine(is);

    // is >> m_prev_prev_dt;
    // GotoNextLine(is);

    // Low coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while(lis >> word)
        {
            prob_lo[i++] = std::stod(word);
        }
    }

    // High coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while(lis >> word)
        {
            prob_hi[i++] = std::stod(word);
        }
    }

    // Set up problem domain
    RealBox rb(prob_lo, prob_hi);
    Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        SetGeometry(lev, Geometry(Geom(lev).Domain(), rb, Geom(lev).CoordInt(),
                                  Geom(lev).isPeriodic()));
    }

    for(int lev = 0; lev <= finest_level; ++lev)
    {
        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // Create distribution mapping
        DistributionMapping dm{ba, ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(lev, cur_time, ba, dm);
    }

    /***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/

    // Load the field data
    for(int lev = 0; lev <= finest_level; ++lev)
    {
        VisMF::Read(m_leveldata[lev]->msv[Ark::PHI],
            amrex::MultiFabFileFullPrefix(lev, name_restart_file, level_prefix, "Phase"));

        VisMF::Read(m_leveldata[lev]->msv[Ark::RHO],
            amrex::MultiFabFileFullPrefix(lev, name_restart_file, level_prefix, "Density"));

        VisMF::Read(m_leveldata[lev]->msv[Ark::VEL],
            amrex::MultiFabFileFullPrefix(lev, name_restart_file, level_prefix, "Velocity"));

        VisMF::Read(m_leveldata[lev]->msv[Ark::P],
            amrex::MultiFabFileFullPrefix(lev, name_restart_file, level_prefix, "P"));

        VisMF::Read(m_leveldata[lev]->msv[Ark::T],
            amrex::MultiFabFileFullPrefix(lev, name_restart_file, level_prefix, "T"));
    }

    amrex::Print() << "Restart complete" << std::endl;
}

void DUGKSAmr::WriteJobInfo(const std::string& path) const
{
    if(ParallelDescriptor::IOProcessor())
    {
        // job_info file with details about the run
        std::ofstream jobInfoFile;
        std::string FullPathJobInfoFile = path;
        std::string PrettyLine =
            "===============================================================================\n";

        FullPathJobInfoFile += "/DUGKSAmr_job_info";
        jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

        // job information
        jobInfoFile << PrettyLine;
        jobInfoFile << " DUGKSAmr Job Information\n";
        jobInfoFile << PrettyLine;

#ifdef ARK_USE_PHASE
        jobInfoFile << "USE_PHASE------->";
        jobInfoFile << "Velocity Model : D"<<AMREX_SPACEDIM<<"V"<<Ark::DVM_PHASE<<"\n";
#endif
#ifdef ARK_USE_MOMENTUM
        jobInfoFile << "USE_MOMENTUM---->";
        jobInfoFile << "Velocity Model : D"<<AMREX_SPACEDIM<<"V"<<Ark::DVM_MOMENTUM<<"\n";
#endif
#ifdef ARK_USE_ENERGY
        jobInfoFile << "USE_ENERGY------>";
        jobInfoFile << "Velocity Model : D"<<AMREX_SPACEDIM<<"V"<<Ark::DVM_ENERGY<<"\n";
#endif

        jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
        jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

        jobInfoFile << "\n\n";

        // build information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Build Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
        jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
        jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
        jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
        jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
        jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << "\n";
        jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

        jobInfoFile << "\n";

        const char* githash1 = buildInfoGetGitHash(1);
        const char* githash2 = buildInfoGetGitHash(2);
        if(std::strlen(githash1) > 0)
        {
            jobInfoFile << "DUGKSAmr git hash: " << githash1 << "\n";
        }
        if(std::strlen(githash2) > 0)
        {
            jobInfoFile << "AMReX git hash: " << githash2 << "\n";
        }

        jobInfoFile << "\n\n";

        // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for(int i = 0; i <= finest_level; i++)
        {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for(int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                jobInfoFile << geom[i].Domain().length(dir) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << "\n\n";

        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();
    }
}

void DUGKSAmr::WritePlotFile() 
{
    BL_PROFILE("DUGKSAmr::WritePlotFile()");
    VisMF::SetNOutFiles(128);

//     if (m_plt_vort or m_plt_divu or m_plt_strainrate) {
//         for (int lev = 0; lev <= finest_level; ++lev) {
// #ifdef AMREX_USE_EB
//             const int ng = (EBFactory(0).isAllRegular()) ? 1 : 2;
// #else
//             const int ng = 1;
// #endif
//             fillpatch_velocity(lev, cur_time, m_leveldata[lev]->msv[Ark::VEL], ng);
//             fillpatch_density(lev, cur_time, m_leveldata[lev]->msv[Ark::RHO], ng);
//         }
//     }

    const std::string& plotfilename = amrex::Concatenate(name_plot_file, nstep);

    amrex::Print()  << "  Writing plotfile " << plotfilename 
                    << " at time " << cur_time << std::endl;

    int ncomp = 0;

    // Phase
    if (m_plt_phi) ++ncomp;

    // Density
    if (m_plt_rho) ++ncomp;

    // Velocity components
    if (m_plt_velx) ++ncomp;
    if (m_plt_vely) ++ncomp;
#if (AMREX_SPACEDIM == 3)
    if (m_plt_velz) ++ncomp;
#endif

    // Pressure
    if(m_plt_pres) ++ncomp;

    // Temperature
    if(m_plt_temp) ++ncomp;

    // Phase gradient components
    if (m_plt_gphix) ++ncomp;
    if (m_plt_gphiy) ++ncomp;
#if (AMREX_SPACEDIM == 3)
    if (m_plt_gphiz) ++ncomp;
#endif

    // Force components
    if (m_plt_forcex) ++ncomp;
    if (m_plt_forcey) ++ncomp;
#if (AMREX_SPACEDIM == 3)
    if (m_plt_forcez) ++ncomp;
#endif

    // Apparent viscosity
    if(m_plt_tau) ++ncomp;

    // Vorticity
    if(m_plt_vort) ++ncomp;

    // Magnitude of the rate-of-strain tensor 
    if(m_plt_strainrate) ++ncomp;

    // Divergence of velocity field
    if(m_plt_divu) ++ncomp;

#ifdef AMREX_USE_EB
    // Cut cell volume fraction
    if(m_plt_volfrac) ++ncomp;
#endif

    Vector<MultiFab> mf(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
    }

    Vector<std::string> pltscaVarsName;
    int icomp = 0;
    if (m_plt_phi) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::PHI], 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("phase");
        ++icomp;
    }
    if (m_plt_rho) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::RHO], 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("density");
        ++icomp;
    }
    if (m_plt_velx) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::VEL], 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("velx");
        ++icomp;
    }
    if (m_plt_vely) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::VEL], 1, icomp, 1, 0);
        }
        pltscaVarsName.push_back("vely");
        ++icomp;
    }
#if (AMREX_SPACEDIM == 3)
    if (m_plt_velz) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::VEL], 2, icomp, 1, 0);
        }
        pltscaVarsName.push_back("velz");
        ++icomp;
    }
#endif
    if (m_plt_pres) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::P], 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("pressure");
        ++icomp;
    }
    if (m_plt_temp) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::T], 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("temperature");
        ++icomp;
    }
    if (m_plt_gphix) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::GradPHI], 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gphix");
        ++icomp;
    }
    if (m_plt_gphiy) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::GradPHI], 1, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gphiy");
        ++icomp;
    }
#if (AMREX_SPACEDIM == 3)
    if (m_plt_gphiz) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::GradPHI], 2, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gphiz");
        ++icomp;
    }
#endif
    if (m_plt_forcex) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::FORCE], 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("forcex");
        ++icomp;
    }
    if (m_plt_forcey) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::FORCE], 1, icomp, 1, 0);
        }
        pltscaVarsName.push_back("forcey");
        ++icomp;
    }
#if (AMREX_SPACEDIM == 3)
    if (m_plt_forcez) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->msv[Ark::FORCE], 2, icomp, 1, 0);
        }
        pltscaVarsName.push_back("forcez");
        ++icomp;
    }
#endif
    if (m_plt_tau) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->fDVDF[Ark::TAU], 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("tauMOM");
        ++icomp;
    }
    // if (m_plt_p) {
    //     for (int lev = 0; lev <= finest_level; ++lev) {
    //         amrex::average_node_to_cellcenter(mf[lev], icomp, m_leveldata[lev]->p, 0, 1);
    //     }
    //     pltscaVarsName.push_back("p");
    //     ++icomp;
    // }
    // if (m_plt_eta) {
    //     for (int lev = 0; lev <= finest_level; ++lev) {
    //         MultiFab vel_eta(mf[lev], amrex::make_alias, icomp, 1);
    //         compute_viscosity_at_level(lev,
    //                                    &vel_eta,
    //                                    &m_leveldata[lev]->msv[Ark::RHO],
    //                                    &m_leveldata[lev]->msv[Ark::VEL],
    //                                    Geom(lev),
    //                                    cur_time, 0);
    //     }
    //     pltscaVarsName.push_back("eta");
    //     ++icomp;
    // }
    if (m_plt_vort) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab vort(mf[lev], amrex::make_alias, icomp, 1);
            ComputeVorticity(lev, cur_time, vort, m_leveldata[lev]->msv[Ark::VEL]);
        }
        pltscaVarsName.push_back("vort");
        ++icomp;
    }
    if (m_plt_strainrate) {
        
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab strainrate(mf[lev], amrex::make_alias, icomp, 1);
            ComputeStrainrateAtLevel(lev,
                                        &strainrate,
                                        &m_leveldata[lev]->msv[Ark::VEL],
                                        Geom(lev),
                                        cur_time, 0);
        }
        pltscaVarsName.push_back("strainrate");
        ++icomp;
    }
    if (m_plt_divu) {
        amrex::Abort("plt_divu: xxxxx TODO");
        pltscaVarsName.push_back("divu");
        ++icomp;
    }
#ifdef AMREX_USE_EB
    if (m_plt_volfrac) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], EBFactory(lev).getVolFrac(), 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("volfrac");
        ++icomp;
    }
#endif

#ifdef AMREX_USE_EB
    for (int lev = 0; lev <= finest_level; ++lev) {
        EB_set_covered(mf[lev], 0.0);
    }
#endif

    AMREX_ALWAYS_ASSERT(ncomp == static_cast<int>(pltscaVarsName.size()));

    // This needs to be defined in order to use amrex::WriteMultiLevelPlotfile, 
    // but will never change unless we use subcycling. 
    // If we do use subcycling, this should be a incflo class member. 
    Vector<int> istep(finest_level + 1, nstep);

    // Write the plotfile
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf), 
                                   pltscaVarsName, Geom(), cur_time, istep, refRatio());
    WriteJobInfo(plotfilename);
}
void DUGKSAmr::WriteL2NormFile(const Vector<double> &cur_L2_norm)
{
    if(ParallelDescriptor::IOProcessor()) 
    {
        std::ofstream L2NormFile;
        std::ostringstream oss_L2;   
        oss_L2<<"L2_norm_Pe"<<Pe()
            <<"_TauM"<<TauPhase()
            <<"_NL"<<ChNum()<<"_"
            <<CaseName()<<".dat";
        L2NormFile.open(oss_L2.str().c_str(), std::ios::app);
        L2NormFile  << std::setiosflags(std::ios::scientific) 
                    << std::setprecision(Ark::digitsNum);
        L2NormFile << cur_time <<"    ";
        for(auto d : cur_L2_norm) 
        {
            L2NormFile << d <<"    ";
        }
        L2NormFile<<'\n';

        L2NormFile.close();

        // std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(16);
        // std::cout << cur_time <<"    "<<cur_L2_norm[0]<<'\n';
    }
}
void DUGKSAmr::WriteResidualFile(const Vector<double> &cur_residual)
{
    if(ParallelDescriptor::IOProcessor()) 
    {
        std::ofstream ResidualFile;
        std::ostringstream oss_res;  
        oss_res <<"Residual_Pe"<<Pe()
                <<"_Cn"<<Cn()
                <<"_NL"<<ChNum()<<"_"
                <<CaseName()<<".dat";
        ResidualFile.open(oss_res.str().c_str(), std::ios::app);
        ResidualFile << std::setiosflags(std::ios::scientific) 
                     << std::setprecision(Ark::digitsNum);
        ResidualFile << cur_time <<"    ";
        for(auto d : cur_residual) 
        {
            ResidualFile << d <<"    ";
        }
        ResidualFile<<'\n';

        ResidualFile.close();
    }
}
void DUGKSAmr::WriteL2NormVelocity(int const lev, double const vel_diff, double const vel_norm)
{
    if(ParallelDescriptor::IOProcessor()) 
    {
        std::ofstream L2NormVelocity;
        std::ostringstream oss_velocity;
        oss_velocity <<"Level"<<lev<<"_"
                     <<"L2norm_velocity"
                     << "_Pe"<<Pe()
                     << "_Cn"<<Cn()
                     << "_NL"<<ChNum()<<"_"
                     << CaseName()<<".dat";
        L2NormVelocity.open(oss_velocity.str().c_str(), std::ios::app);
        L2NormVelocity  << std::setiosflags(std::ios::scientific) 
                        << std::setprecision(Ark::digitsNum);
        L2NormVelocity  << cur_time <<"    "<<vel_diff<<"    "<<vel_norm<<"\n";
        L2NormVelocity.close();
    }
}
void DUGKSAmr::WriteSumPhaseFile(const Vector<double> &cur_sum)
{
    if(ParallelDescriptor::IOProcessor()) 
    {
        std::ofstream SumFile;
        std::ostringstream oss_sum;  
        oss_sum <<"SumPhase_Pe"<<Pe()
                <<"_Cn"<<Cn()
                <<"_NL"<<ChNum()<<"_"
                <<CaseName()<<".dat";
        SumFile.open(oss_sum.str().c_str(), std::ios::app);
        SumFile << std::setiosflags(std::ios::scientific) 
                     << std::setprecision(Ark::digitsNum);
        SumFile << cur_time <<"    ";
        for(auto d : cur_sum) 
        {
            SumFile << d <<"    ";
        }
        SumFile<<'\n';

        SumFile.close();
    }
}
// void DUGKSAmr::WritePressureFile(int i = 0, j = 0)
// {
//     if(ParallelDescriptor::IOProcessor()) 
//     {
//         std::ofstream PresFile;
//         std::ostringstream oss_pres;  
//         oss_pres <<"Pressure_Pe"<<Pe()
//                 <<"_Cn"<<Cn()
//                 <<"_NL"<<ChNum()<<"_"
//                 <<CaseName()<<".dat";
//         PresFile.open(oss_pres.str().c_str(), std::ios::app);
//         PresFile << std::setiosflags(std::ios::scientific) 
//                      << std::setprecision(Ark::digitsNum);
//         PresFile << cur_time <<"    "<<Data(lev);
//         for(auto d : cur_sum) 
//         {
//             PresFile << d <<"    ";
//         }
//         PresFile<<'\n';

//         PresFile.close();
//     }
// }