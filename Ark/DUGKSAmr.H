#ifndef _ARK_DUGKSAMR_H
#define _ARK_DUGKSAMR_H

#include <DUGKSCore.H>
#include <DUGKSFlowBase.H>
#include <AMReX_Math.H>
#include <AMReX_BCRec.H>
#include <AMReX_AmrCore.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MultiFabUtil.H>
// #include <AMReX_ParallelDescriptor.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBSupport.H>
#include <AMReX_EBInterpolater.H>
using amrex::EBSupport;
#endif

using amrex::Dim3;
using amrex::XDim3;
using amrex::Array;
using amrex::BCRec;
using amrex::MFIter;
using amrex::MFInfo;
using amrex::Vector;
using amrex::IntVect;
using amrex::RealBox;
using amrex::BoxArray;
using amrex::Geometry;
using amrex::FArrayBox;
using amrex::ParmParse;
using amrex::FabFactory;
using amrex::Orientation;
using amrex::Periodicity;
using amrex::TagBoxArray;
using amrex::GeometryData;
using amrex::FArrayBoxFactory;
using amrex::DistributionMapping;

namespace Ark {

    void FillExternalDirichlet (const IntVect& iv, Array4<Real> const& q,
                                const int dcomp, const int numcomp,
                                GeometryData const& geom, const Real time,
                                const BCRec* bcr, const int bcomp,
                                const int orig_comp);

    void FillExternalNeumann (const IntVect& iv, Array4<Real> const& q,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp);

    void FillExternalRobin (const IntVect& iv, Array4<Real> const& q,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp);
}

class DUGKSAmr final: public Ark::FlowParm, public amrex::AmrCore
{

public:
    static const std::string DmVnModel;

public:
	DUGKSAmr();

	DUGKSAmr (DUGKSAmr&& rhs) = default;
    DUGKSAmr& operator= (DUGKSAmr&& rhs) = default;

    DUGKSAmr (const DUGKSAmr& rhs) = delete;
    DUGKSAmr& operator= (const DUGKSAmr& rhs) = delete;

	virtual ~DUGKSAmr();

    // Initialize multilevel AMR data
    void InitData ();

    // Tag cells for refinement
	virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, 
                           amrex::Real time, int ngrow) override;

    //! Manually tag.  Note that tags is built on level lev grids coarsened by bf_lev[lev].
    virtual void ManualTagsPlacement (int lev, TagBoxArray& tags, 
                                      const Vector<IntVect>& bf_lev) override;

	// Make a new level from scratch using provided BoxArray and DistributionMapping
	// Only used during initialization
	virtual void MakeNewLevelFromScratch (int lev, amrex::Real time,
	                                      const amrex::BoxArray& new_grids,
	                                      const amrex::DistributionMapping& new_dmap) override;

	// Make a new level using provided BoxArray and DistributionMapping,
	// and fill with interpolated coarse level data
	virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time,
	                                     const amrex::BoxArray& ba,
	                                     const amrex::DistributionMapping& dm) override;

	// Remake an existing level using provided BoxArray and DistributionMapping,
	// and fill with existing fine and coarse data
	virtual void RemakeLevel (int lev, amrex::Real time,
	                          const amrex::BoxArray& ba,
	                          const amrex::DistributionMapping& dm) override;

	// Delete level data
	virtual void ClearLevel (int lev) override;

    void (DUGKSAmr::*Advance)() = nullptr;

    //trapezium rule
    void Advance_CN();

    //NonEquilibrium rule
    void Advance_NEQ();

    //Implicit Euler
    void Advance_BackwardEuler();

    //Mixed
    void Advance_IMEX();

    //Strang-Splitting 
    void Advance_SS();

    // Finite difference scheme
    void Advance_FD(); 

    //DVM
    void Advance_DVM();

    void Evolve();

    //
    double CFL() const { return Courant_Friedrichs_Lewy; } 
    double Dt() const { return advance_time_step; }
    double hDt(int const lev) const { return trace_time_step[lev]; }

public:
    ///////////////////////////////////////////////////////////////////////////
    //
    // boundary_conditions
    //
    ///////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////////
    //
    // case
    //
    ///////////////////////////////////////////////////////////////////////////
    enum class LV { drop = -1, bubble =1 };

    void CaseInitFluid (int const lev);

#ifdef AMREX_USE_EB
    static int m_eb_basic_grow_cells;
    static int m_eb_volume_grow_cells;
    static int m_eb_full_grow_cells;

    static EBSupport      m_eb_support_level;

    static void SetEBMaxGrowCells (int nbasic, int nvolume, int nfull) noexcept {
        m_eb_basic_grow_cells = nbasic;
        m_eb_volume_grow_cells = nvolume;
        m_eb_full_grow_cells = nfull;
    }

    static void SetEBSupportLevel (EBSupport ebs) { m_eb_support_level = ebs; }
#endif

private:
    // Be verbose?
    int m_verbose = 0;

    //ghost layer
    int n_ghost = 0;

    amrex::Vector<amrex::Real> m_t_old;
    amrex::Vector<amrex::Real> m_t_new;

    // current time
    amrex::Real cur_time     = 0.0;

    // Time step
    bool use_mutable_trace = true;
    int nstep         = -1;
    amrex::Real Courant_Friedrichs_Lewy;
    amrex::Real advance_time_step;
    amrex::Vector<amrex::Real> trace_time_step;

    // AMR / refinement settings
    int m_refine_cutcells = -1;
    int m_regrid_int = -1;

    amrex::Vector<amrex::Real> phi_error;
    amrex::Vector<amrex::Real> grad_phi_error;
    amrex::Vector<amrex::Real> vel_error;

    //tag a specified region
    amrex::Vector<int> tag_region;
    amrex::Vector<amrex::Real> tag_region_lo;
    amrex::Vector<amrex::Real> tag_region_hi;

    //control parameters
    int end_step = 0;
    bool steady_state = false;
    Real residual_tolerance = 1;
    Real adapt_tolerance = 0;

    //
    Real residual = -1;
    //plot control
    std::string name_plot_file{"plt"};
    std::string name_check_file{"chk"};
    std::string name_restart_file{""};
    std::string name_tag_file{""};

    int write_conv_interval = -1;
    int write_residual_interval = -1;
    int write_check_interval = -1;
    int write_plot_interval = -1;
    bool write_plot_on_restart = false;

    int m_plt_phi         = 1;
    int m_plt_rho         = 0;
    int m_plt_velx        = 0;
    int m_plt_vely        = 0;
    int m_plt_velz        = 0;
    int m_plt_pres        = 0;
    int m_plt_temp        = 0;
    int m_plt_gphix       = 0;
    int m_plt_gphiy       = 0;
    int m_plt_gphiz       = 0;
    int m_plt_forcex      = 0;
    int m_plt_forcey      = 0;
    int m_plt_forcez      = 0;

    int m_plt_chemip      = 0;
    int m_plt_tau         = 0;
    int m_plt_vort        = 0;
    int m_plt_strainrate  = 0;
    int m_plt_divu        = 0;
    int m_plt_volfrac     = 0; //volume fraction for cut cells

    int m_plt_hDVDF       = 0;
    int m_plt_fDVDF       = 0;
    int m_plt_gDVDF       = 0;

    struct LevelData {
        LevelData () = default;
        LevelData (amrex::BoxArray const& ba,
                   amrex::DistributionMapping const& dm,
                   amrex::FabFactory<amrex::FArrayBox> const& fact, 
                   int nghost);

        LevelData(const LevelData& ld) = delete;
        LevelData(LevelData &&ld) = delete;

        Vector<MultiFab*> get_hDVDFnF_flux() noexcept;
        Vector<MultiFab*> get_fDVDFnF_flux() noexcept;
        Vector<MultiFab*> get_gDVDFnF_flux() noexcept;

        Vector<MultiFab const*> get_hDVDFnF_flux_const() const noexcept;
        Vector<MultiFab const*> get_fDVDFnF_flux_const() const noexcept;
        Vector<MultiFab const*> get_gDVDFnF_flux_const() const noexcept;

        Vector<MultiFab const*> get_hDVDFnF_eq_const() const noexcept;
        Vector<MultiFab const*> get_fDVDFnF_eq_const() const noexcept;
        Vector<MultiFab const*> get_gDVDFnF_eq_const() const noexcept;

        Vector<MultiFab const*> get_face_phase_const() const noexcept;
        Vector<MultiFab const*> get_face_density_const() const noexcept;
        Vector<MultiFab const*> get_face_velocity_const() const noexcept;
        Vector<MultiFab const*> get_face_pressure_const() const noexcept;
        Vector<MultiFab const*> get_face_temperature_const() const noexcept;

        // macroscopic multifabs : cell-centered 
       
        // MultiFab pressure;
        // MultiFab temperature;
        // MultiFab force;
        // MultiFab grad_phase; // phase gradient
        // MultiFab grad_density;

        //macroscopic multifabs
        Array<MultiFab, Ark::MSVNUM> msv;
        Array<MultiFab, Ark::MSVNUM> msv_xF;
        Array<MultiFab, Ark::MSVNUM> msv_yF;
#if(AMREX_SPACEDIM == 3)
        Array<MultiFab, Ark::MSVNUM> msv_zF;
#endif
        MultiFab phase_prev;
        MultiFab velocity_prev;

        MultiFab phase_old;
        MultiFab density_old;
        MultiFab velocity_old;

        //mesoscopic multifabs for MASS

        //cell centered	
        Array<MultiFab, Ark::DVDFNUM> hDVDF;
        //face centered
        Array<MultiFab, Ark::FACEDVDFNUM> hDVDFxF;
        Array<MultiFab, Ark::FACEDVDFNUM> hDVDFyF;
#if(AMREX_SPACEDIM == 3)
        Array<MultiFab, Ark::FACEDVDFNUM> hDVDFzF;
#endif
        //mesoscopic multifabs for MOMENTUM

        //cell centered 
        Array<MultiFab, Ark::DVDFNUM> fDVDF;
        //face centered
        Array<MultiFab, Ark::FACEDVDFNUM> fDVDFxF;
        Array<MultiFab, Ark::FACEDVDFNUM> fDVDFyF;
#if(AMREX_SPACEDIM == 3)
        Array<MultiFab, Ark::FACEDVDFNUM> fDVDFzF;
#endif
        //mesoscopic multifabs for ENERGY
        
        //cell centered 
        Array<MultiFab, Ark::DVDFNUM> gDVDF;
        //face centered
        Array<MultiFab, Ark::FACEDVDFNUM> gDVDFxF;
        Array<MultiFab, Ark::FACEDVDFNUM> gDVDFyF;  
#if(AMREX_SPACEDIM == 3)
        Array<MultiFab, Ark::FACEDVDFNUM> gDVDFzF;
#endif
    private:
        void ConstructMSVMultiFabAtFace (Array<MultiFab, Ark::MSVNUM>& msv_nF,
                                        amrex::BoxArray const& ba_nF,
                                        amrex::DistributionMapping const& dm,
                                        amrex::FabFactory<FArrayBox> const& fact);

        void ConstructDVDFMultiFabAtCell (Array<MultiFab, Ark::DVDFNUM>& xDVDF,
                                        amrex::BoxArray const& ba,
                                        amrex::DistributionMapping const& dm,
                                        amrex::FabFactory<FArrayBox> const& fact,
                                        int const Q, int const ng);

        void ConstructDVDFMultiFabAtFace (Array<MultiFab, Ark::FACEDVDFNUM>& xDVDFnF,
                                        amrex::BoxArray const& ba_nF,
                                        amrex::DistributionMapping const& dm,
                                        amrex::FabFactory<FArrayBox> const& fact,
                                        int const Q, int const ng);
      
    };

    amrex::Vector<std::unique_ptr<LevelData> > m_leveldata;

    amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_factory;

    amrex::FabFactory<amrex::FArrayBox> const&
    Factory (int lev) const noexcept { return *m_factory[lev]; }

    LevelData& Data(int lev) noexcept {return *m_leveldata[lev]; }

    const LevelData& Data(int lev) const noexcept {return *m_leveldata[lev]; }

#ifdef AMREX_USE_EB
    amrex::EBFArrayBoxFactory const&
    EBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*m_factory[lev]);
    }
#endif

    //boundary data
    amrex::Vector<amrex::BCRec> m_bcrec_physical;
    amrex::Vector<amrex::BCRec> m_bcrec_phase;
    amrex::Vector<amrex::BCRec> m_bcrec_density;
    amrex::Vector<amrex::BCRec> m_bcrec_velocity;
    amrex::Vector<amrex::BCRec> m_bcrec_pressure;
    amrex::Vector<amrex::BCRec> m_bcrec_temperature;
    amrex::Vector<amrex::BCRec> m_bcrec_chemicalP;
    amrex::Vector<amrex::BCRec> m_bcrec_force;

    // amrex::Vector<amrex::BCRec> m_bcrec_gphi;
    // amrex::Vector<amrex::BCRec> m_bcrec_gdensity;

    amrex::Vector<amrex::BCRec> m_bcrec_DVDF_wave;
    amrex::Vector<amrex::BCRec> m_bcrec_DVDF_barplus;

    using BCValue = amrex::Array<amrex::Real, 2*AMREX_SPACEDIM>;

    BCValue m_dirichlet_phase;
    BCValue m_dirichlet_density;
    BCValue m_dirichlet_velocity_x;
    BCValue m_dirichlet_velocity_y;
    BCValue m_dirichlet_velocity_z;
    BCValue m_dirichlet_pressure;
    BCValue m_dirichlet_temperature;

    Vector<BCRec> const& get_physical_bcrec () const noexcept { return m_bcrec_physical;} 
    Vector<BCRec> const& get_phase_bcrec () const noexcept { return m_bcrec_phase;}
    Vector<BCRec> const& get_density_bcrec () const noexcept { return m_bcrec_density; }
    Vector<BCRec> const& get_velocity_bcrec () const noexcept { return m_bcrec_velocity; }
    Vector<BCRec> const& get_pressure_bcrec () const noexcept { return m_bcrec_pressure; }
    Vector<BCRec> const& get_temperature_bcrec () const noexcept { return m_bcrec_temperature;}
    Vector<BCRec> const& get_chemicalP_bcrec () const noexcept { return m_bcrec_chemicalP;}
    Vector<BCRec> const& get_force_bcrec () const noexcept { return m_bcrec_force;}
    Vector<BCRec> const& get_DVDF_wave_bcrec () const noexcept { return m_bcrec_DVDF_wave;}
    Vector<BCRec> const& get_DVDF_barplus_bcrec () const noexcept { return m_bcrec_DVDF_barplus;}

    BCValue const& get_dirichlet_phase() const noexcept 
    {
        return m_dirichlet_phase; 
    }
    BCValue const& get_dirichlet_density() const noexcept 
    {
        return m_dirichlet_density; 
    }
    BCValue const& get_dirichlet_velocity_x() const noexcept 
    {
        return m_dirichlet_velocity_x; 
    }
    BCValue const& get_dirichlet_velocity_y() const noexcept 
    {
        return m_dirichlet_velocity_y; 
    }
    BCValue const& get_dirichlet_velocity_z() const noexcept 
    {
        return m_dirichlet_velocity_z; 
    }
    BCValue const& get_dirichlet_pressure() const noexcept 
    {
        return m_dirichlet_pressure; 
    }
    BCValue const& get_dirichlet_temperature() const noexcept 
    {
        return m_dirichlet_temperature; 
    }
//
    int nghost_state () const 
    {
    #ifdef AMREX_USE_EB
            return n_ghost + 1;
    #endif
            return n_ghost;
    }

    int zero_ghost () const { return 0;}
    int unit_ghost () const { return 1;}
    int double_ghost() const { return 2;}
    int triple_ghost() const { return 3;}

private:
    ///////////////////////////////////////////////////////////////////////////
    //
    // tag region
    //
    ///////////////////////////////////////////////////////////////////////////

    //DUGKSAmrTag.cpp

    void TagScalarField (TagBoxArray& Tab, const MultiFab &ScalarFab, const Real &tagError);

    void TagVectorField (TagBoxArray& Tab, const MultiFab &VectorFab, const Real &tagError);

    void TagMultiVectorField (TagBoxArray& Tab, Vector<MultiFab const*> &VectorFabPtr, 
                              const Vector<double> &tagError);
    
    void TagDesignatedRegion(TagBoxArray& Tab, const Real *prob_lo, const Real *dx);

    ///////////////////////////////////////////////////////////////////////////
    //
    // boundary_conditions
    //
    ///////////////////////////////////////////////////////////////////////////

    //DUGKSSetBC.cpp
    void InitBoundaryConditions();

    void ReadPhysicalBoundaryCondition();
    
    void ReadBoundaryCondition(std::string const &var, Vector<BCRec> &m_bcrec, 
                               int const ncomp = 1);

    void SetPeriodicityOverwrite(int const dir);
    void SetPeriodicityOverwrite(int const dir, Vector<BCRec> &m_bcrec);
    void SetPeriodicBoundaryCondition();

    void ReadDirichletBoundaryCondition(std::string const &var, BCValue &m_dirichlet);


    //DUGKSPhysBC.cpp
    void FillPhysicalBoundaryMacroVars(int const lev);
    void FillPhysicalBoundaryMacroVars(int const lev, int const bcid, int const dir);
    void FillPhysicalBoundaryMacroVarsNoSlipWall (int const lev, int const dir);
    void FillPhysicalBoundaryMacroVarsSymmetry (int const lev, int const dir);
    void FillPhysicalBoundaryMacroVarsOutFlow (int const lev, int const dir);
    void FillInteriorBoundary (MultiFab &xFab, const Periodicity& period);

    //DUGKSDVDFBC.cpp
    void FillPhysicalBoundaryDVDF(int const lev);
    void FillPhysicalBoundaryDVDF(int const lev, int const bcid, int const dir);

    
    //family of average down
    void AverageDownAll();
    void AverageDownMacroVars();
    void AverageDownPhase();
    void AverageDownDensity();
    void AverageDownVelocity();
    void AverageDownPressure();
    void AverageDownTemperature();
    
    void AverageDownDVDFWave();
    void AverageDownDVDFBarPlus();
    void AverageDownDVDFFlux();

    void AverageDownLevelAll(int const crse_lev);
    void AverageDownLevelMacroVars(int const crse_lev);
    void AverageDownLevelPhase(int const crse_lev);
    void AverageDownLevelDensity(int const crse_lev);
    void AverageDownLevelVelocity(int const crse_lev);
    void AverageDownLevelPressure(int const crse_lev);
    void AverageDownLevelPotential(int const crse_lev);
    void AverageDownLevelTemperature(int const crse_lev);
    void AverageDownLevelGradPhase(int const crse_lev);
    void AverageDownLevelGradDensity(int const crse_lev);
    void AverageDownLevelForce(int const crse_lev);

    void AverageDownLevelDVDFWave(int const crse_lev);
    void AverageDownLevelDVDFBarPlus(int const crse_lev);
    void AverageDownLevelDVDFFlux(int const crse_lev);

    //family of fill patch
    void FillPatchNew (int const lev, Real const time, DUGKSAmr::LevelData &newld);
    void FillPatchMacroVars ();
    void PostFillPatchNew (int const lev, DUGKSAmr::LevelData &newld);
    
    
    void FillPatchPhase       (int const lev, Real const time, MultiFab& phase);
    void FillPatchPhasePrev   (int const lev, Real const time, MultiFab& phasePrev);
    void FillPatchDensity     (int const lev, Real const time, MultiFab& density);
    void FillPatchVelocity    (int const lev, Real const time, MultiFab& velocity);
    void FillPatchVelocityPrev(int const lev, Real const time, MultiFab& velocityPrev);
    void FillPatchPressure    (int const lev, Real const time, MultiFab& pres);
    void FillPatchTemperature (int const lev, Real const time, MultiFab& temp);
    void FillPatchChemicalP   (int const lev, Real const time, MultiFab& chemiP);
    void FillPatchGradPhase   (int const lev, Real const time, MultiFab& grad_phi);
    void FillPatchGradDensity (int const lev, Real const time, MultiFab& grad_rho);
    void FillPatchForce       (int const lev, Real const time, MultiFab& force);

    void FillPatchDVDFWave (int const lev, Real const time, MultiFab& xWave, 
                            Array<MultiFab, Ark::DVDFNUM> LevelData::*xDVDF);

    void FillPatchDVDFBarPlus (int const lev, Real const time, MultiFab& xBP, 
                               Array<MultiFab, Ark::DVDFNUM> LevelData::*xDVDF);
    template<typename F>
    void FillPatchCore (int const lev, Real const time, MultiFab &mf, 
                        const Vector<MultiFab*>& cmf, const Vector<Real>& ctime,
                        const Vector<MultiFab*>& fmf, const Vector<Real>& ftime,
                        const Vector<BCRec>& bcs);

    //family of fill coarse patch
    void FillCoarsePatchNew(int const lev, Real const time, DUGKSAmr::LevelData &newld);

    void FillCoarsePatchPhase       (int const lev, Real const time, MultiFab& phase);
    void FillCoarsePatchDensity     (int const lev, Real const time, MultiFab& density);
    void FillCoarsePatchVelocity    (int const lev, Real const time, MultiFab& velocity);
    void FillCoarsePatchPressure    (int const lev, Real const time, MultiFab& pres);
    void FillCoarsePatchTemperature (int const lev, Real const time, MultiFab& temp);

    void FillCoarsePatchDVDFWave (int const lev, Real const time, MultiFab& xWave, 
                                  Array<MultiFab, Ark::DVDFNUM> LevelData::*xDVDF);

    void FillCoarsePatchDVDFBarPlus (int const lev, Real const time, MultiFab& xBP, 
                                     Array<MultiFab, Ark::DVDFNUM> LevelData::*xDVDF);
    template<typename F>
    void FillCoarsePatchCore(int const lev, Real const time, 
                             MultiFab &mf, MultiFab const &cmf,
                             const Vector<BCRec> &bcrec);

    //DUGKSSetBC.cpp
    void average_scalar_x (Array4<Real> const &dst, Array4<Real const> const &src, 
                           const Dim3 &lo, const Dim3 &hi, int const /*fixed index*/, int const ncomp = 0);

    void average_scalar_y (Array4<Real> const &dst, Array4<Real const> const &src, 
                           const Dim3 &lo, const Dim3 &hi, int const /*fixed index*/, int const ncomp = 0);

    void average_scalar_z (Array4<Real> const &dst, Array4<Real const> const &src, 
                           const Dim3 &lo, const Dim3 &hi, int const /*fixed index*/, int const ncomp = 0);

    void average_vector_y (Array4<Real> const &dst, Array4<Real const> const &src, 
                           const Dim3 &lo, const Dim3 &hi, int const /*fixed index*/);

    void set_dirichlet_velocity_xF (Array4<Real> const &vel, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_velocity_yF (Array4<Real> const &vel, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_velocity_zF (Array4<Real> const &vel, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_pressure_xF (Array4<Real> const &pres, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_pressure_yF (Array4<Real> const &pres, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_pressure_zF (Array4<Real> const &pres, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_phase_xF (Array4<Real> const &phi, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_phase_yF (Array4<Real> const &phi, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_phase_zF (Array4<Real> const &phi, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_density_xF (Array4<Real> const &rho, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_density_yF (Array4<Real> const &rho, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    void set_dirichlet_density_zF (Array4<Real> const &rho, const Dim3 &lo, const Dim3 &hi, 
                                    int const /*fixed index*/, Orientation face);

    ///////////////////////////////////////////////////////////////////////////
    //
    // case
    //
    ///////////////////////////////////////////////////////////////////////////
    //initial the macroscopic variables in level data
    void case_init_macroscopic_variables(int const lev);

    //initial the distribution function to equilibrium
    void case_init_mesoscopic_equilibrium(int const lev);

    void post_case_init(int const lev);

    void BasicCase(int const lev);

    void UniformFlow(int const lev);

    void UnsteadyTaylorGreenVortex(int const lev);

    void ForceTaylorGreenVortex(int const lev);

    void QuiescentDrop(int const lev);

    void TransDrop(int const lev);

    void ZalesakDisk(int const lev);

    void ElongatedInterface (int const lev);

    void SmoothDeformation (int const lev);

    void ReversedSingleVortex(int const lev);

    void LayeredPoiseuille (int const lev);

    void SpinodalDecomposition(int const lev);

    void RayleighTaylorInstability(int const lev);

    void FlatInterface(int const lev);

    void DropletImpact(int const lev);
    
    void TaylorBubble(int const lev);

    void RisingBubble(int const lev);

    void FallingDrop(int const lev);

    void DropletCoalescence(int const lev);


    ///////////////////////////////////////////////////////////////////////////
    //
    // discrete_velocity_model
    //
    ///////////////////////////////////////////////////////////////////////////

    void Update_MacroVars(int const lev);
    void Update_MacroVars_hDt(int const lev, double const hTime);

    void Update_MacroBaseVars (Array<MultiFab, Ark::MSVNUM> &msv);
    void Update_MacroSourceVars(Array<MultiFab, Ark::MSVNUM> &msv, const Real interval);

    void Update_Tau_Phase(MultiFab &tauFab, Array<MultiFab, Ark::MSVNUM> &msv);
    void Update_Tau_Momentum(MultiFab &tauFab, Array<MultiFab, Ark::MSVNUM> &msv);
    void Update_Tau_Energy(MultiFab &tauFab, Array<MultiFab, Ark::MSVNUM> &msv);

    void Update_DVDF_Eq(int const lev);
    void Update_DVDF_EqSo(int const lev);

    void Update_DVDF_Eq_hDt(int const lev);
    void Update_DVDF_EqSo_hDt(int const lev);

    ///////////////////////////////////////////////////////////////////////////
    //
    // Evolve
    //
    ///////////////////////////////////////////////////////////////////////////

    void Nullify_Source (int const lev);

    template<bool ENABLESOURCE=true>
    void Construct_DVDF_Wave(int const lev, double const dTime);

    template<bool ENABLESOURCE=true>
    void Construct_DVDF_Prime (int const lev, double const dTime);

    template<bool ENABLESOURCE=true>
    void Construct_DVDF_NonEq (int const lev, double const dTime);

    void Collision_DVM_BarPlus(int const lev);

    void Collision_DVDF_BarPlus (int const lev, double const hTime);

    void Collision_DVDF_BarPlus_NonEq(int const lev, double const hTime);

    void Update_FD_DVDF_BarPlus (int const lev, double const hTime);

    void Update_DVDF_Grad_BarPlus_Null (int const lev);

    void Update_DVDF_Grad_BarPlus_Default (int const lev);

    void Update_DVDF_Grad_BarPlus_Default_Second (int const lev);

    void Update_DVDF_Grad_BarPlus_Default_VanLeer (int const lev);

    void Update_DVDF_Grad_BarPlus_Third (int const lev);

    void Update_DVDF_Grad_BarPlus_Fourth (int const lev);

    void Update_DVDF_Grad_BarPlus_Fourth_Third (int const lev);

    void Update_DVDF_Grad_BarPlus_Fifth (int const lev);

    void Update_DVDF_Barh_UW_Default (int const lev);

    void Update_DVDF_Barh_UW_First (int const lev);
    
    void Update_DVDF_Barh_UW_Third (int const lev);

    void Update_DVDF_Barh_UW_Fourth (int const lev);

    void Update_DVDF_Barh_UW_Fourth_Cross (int const lev);

    void Update_DVDF_Barh_UW_Fourth_Cross_Third (int const lev);

    void Update_DVDF_Barh_UW_Fifth (int const lev);

    void Update_DVDF_Barh_UW_Default_Second (int const lev);

    void Update_DVDF_Barh_CD (int const lev);

    void Update_DVDF_Barh_WENO_Third (int const lev);

    void Update_DVDF_Barh_WENO_Fifth (int const lev);

    void Update_FD_DVDF_Barh_WENO_Third (int const lev);

    void Update_FD_DVDF_Barh_WENO_Fifth (int const lev);

    template <bool ENABLESOURCE=true>
    void Collision_DVDF_Flux (int const lev, double const hTime);

    template <bool ENABLESOURCE=true>
    void Update_DVDF_MesoFlux (int const lev, double const hTime);

    void Update_DVDF_Wave(int const lev, double const hTime);

    void Update_DVDF_Wave_NonEq(int const lev);

    void Update_DVM_Wave(int const lev);

    void Update_DVDF_Wave_SS (int const lev, double const ssTime);

    void Update_DVDF_Tau(int const lev);

    void Update_DVDFnF_Tau (int const lev);

    void (DUGKSAmr::*Update_DVDF_Barh) (int const lev) = nullptr;

    void (DUGKSAmr::*Update_DVDF_Grad_BarPlus) (int const lev) = nullptr;


    ///////////////////////////////////////////////////////////////////////////
    //
    // embedded_boundaries
    //
    ///////////////////////////////////////////////////////////////////////////

#ifdef AMREX_USE_EB
    void MakeEBGeometry ();
    void make_eb_annulus ();
    void make_eb_box ();
    void make_eb_cylinder ();
    void make_eb_twocylinders ();
    void make_eb_regular ();
    void make_eb_sphere ();
    void make_eb_spherecube ();
    void make_eb_cyl_tuscan ();
    void make_eb_tuscan ();
    void WriteMyEBSurface ();
#endif

    ///////////////////////////////////////////////////////////////////////////
    //
    // field operation
    //
    ///////////////////////////////////////////////////////////////////////////
    double HyperbolicTangentProfile (double const distance, double const heavy, double const light, LV token);

    double SourcePhi (double const phi);

    void InterpolateFaceSource (Array<MultiFab, Ark::MSVNUM> &msv_nF,
                                Array<MultiFab, Ark::MSVNUM> const &msv,
                                const IntVect &faceIndex);

    void AllenCahnSource (int const lev);

    void AllenCahnTemporalPhiU(LevelData &ld);

    void CahnHilliardSource (int const lev);

    void FreeEnergySource (int const lev);

    void PseudoPotential (int const lev);

    void ComputeVorticity (int lev, amrex::Real time, amrex::MultiFab& vort,
                           amrex::MultiFab const& vel);
    void ComputeDivU (amrex::Real time);
    void ComputeDrag ();
    amrex::Real ComputeKineticEnergy () const;

    void ComputeStrainrateAtLevel (int lev,
                                   amrex::MultiFab* sr,
                                   amrex::MultiFab* vel,
                                   amrex::Geometry& lev_geom,
                                   amrex::Real time, int nghost);

    void StaircasePhase(MultiFab &PhaseFab);

    void BoundPhase(MultiFab &PhaseFab);

    void BoundGradient(MultiFab &GradPhaseFab);

    void ElongationVelocityField (MultiFab &VelFab, const Real *dx, 
                                  const RealBox &rb, int const token);

    void ElongationVelocityField3D (MultiFab &VelFab, const Real *dx, 
                                    const RealBox &rb, int const token);

    void DeformationVelocityField (MultiFab &VelFab, const Real *dx, 
                                   const Real Time, int const token);

    void ReversedVortexVelocityField (MultiFab &VelFab, const Real *dx, 
                                      const Real time, int const token);

    ///////////////////////////////////////////////////////////////////////////
    //
    // setup
    //
    ///////////////////////////////////////////////////////////////////////////

    void ReadParametersDUGKSAmr ();
    void ReadIOParameters ();
    void ResizeArrays (); // Resize arrays to fit (up to) max_level + 1 AMR levels
    void SetMutableTraceTimeStep(const double hTime);
    void SetFluxScheme();
    void SetAdvanceScheme();

    ///////////////////////////////////////////////////////////////////////////
    //
    // utilities
    //
    ///////////////////////////////////////////////////////////////////////////

    //input and output
    void SteadyStateReached(const Vector<double> & residualv);
    void WriteHeader (const std::string& name, bool is_checkpoint) const;
    void WriteJobInfo (const std::string& dir) const;
    void WritePlotFile ();
    void WriteCheckPointFile () const;
    void ReadCheckpointFile ();

    Vector<double> MultiLevelResidual();
    Vector<double> MultiLevelL2Norm();
    Vector<double> MultiLevelSum();
    Vector<double> MultiLevelSumDensity();

    std::string FileNameSuffix();
    void WriteL2NormFile(const Vector<double> &);
    void WriteResidualFile(const Vector<double> &);
    void WriteSumPhaseFile(const Vector<double> &);
    void WriteSumDensityFile(const Vector<double> &);
    void WriteL2NormVelocity(int const lev, double const vel_diff, double const vel_norm);
    void WriteInterfacePosition(const Vector<double> &interfacePosition);

    void WriteComputingTime(const double interval);

    void WriteTecplotFlowField(int const lev);
    void WriteTecplotGrid(int const lev);
    //diagnose
    void PrintMaxValues (amrex::Real time);
    void PrintMaxVel (int lev);
    void PrintMaxGp (int lev);
    void CheckForNans (int lev);
    void PublicFillPatchNew(int const lev);

    double L2NormSubtract (const MultiFab& xFabA, const MultiFab& xFabB, int const comp = 0);
    double L2NormSubtractZalesakDisk (const MultiFab& xFabA, const MultiFab& xFabB);
    double L2NormVelocity(int lev);
    double L2NormPhase(int lev);
    double L2NormDensity(int lev);
    double SumPhase(int lev);
    double InterfacePosition (const MultiFab& phi, const double rb_lo_y, const Real* dl,
                              int const indexI, int const indexJ, int const indexK);

    double InterfacePositionHighest (const MultiFab& PhaseFab, const double probDomain_lo_y,
                                     const Real* dx, int const index);

    double InterfacePositionLowest (const MultiFab& PhaseFab, const double probDomain_lo_y,
                                    const Real* dx, int const index);

    Vector<double> InterfacePositionRayleighTaylor(int lev);
    Vector<double> InterfacePositionTaylorBubble(int lev);
    Vector<double> InterfacePositionRisingBubble(int lev);
    Vector<double> InterfacePositionLiquidBridge(int lev);

    void FillLevelPhaseZalesakDisk(MultiFab &phase, const Real *dl);
    void FillLevelScalarFieldSphere (MultiFab &scalarFab, const double heavy, const double light,
                                    const Real *dl, LV token);
    void FillLevelUnsteadyTaylorGreenVortex(int const lev, MultiFab &VelFab, MultiFab &PresFab);
    void FillLevelForceTaylorGreenVortex (int const lev, MultiFab &VelFab, MultiFab &PresFab,
                                          MultiFab &ForceFab, int const token);

    void CheckDUGKSAmr();
//
public:
    void CheckForPhysicalBoundary (int const lev);
#ifdef AMREX_DEBUG
    void Verify(int const lev);
    void SetUnit(int const lev);
    void Test();
    void TestRegrid(int const crse_lev);
    void TestFillCoarsePatch(int const fine_lev);
    double BiquadraticValue(double x, double y, double z) { return x*y + y*y + z*x; }
    double BilinearValue(double x, double y, double z) { return 2*x + x*y + y*z; }
#endif

};
#endif