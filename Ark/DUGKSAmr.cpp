#include <DUGKSAmr.H>

// Need this for TagCutCells
#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#endif

#ifdef AMREX_USE_EB
	int DUGKSAmr::m_eb_basic_grow_cells = 5;
	int DUGKSAmr::m_eb_volume_grow_cells = 4;
	int DUGKSAmr::m_eb_full_grow_cells = 2;
	EBSupport DUGKSAmr::m_eb_support_level = EBSupport::volume;
#endif

DUGKSAmr::DUGKSAmr():FlowParm(),AmrCore()
{
    // NOTE: Geometry on all levels has just been defined in the AmrCore
    // constructor. No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    // Initialize memory for data-array internals
    // must be executed at the very start
    ResizeArrays();

    // Read inputs file using ParmParse
    ReadParametersDUGKSAmr();

#ifdef AMREX_USE_EB
    // This is needed before initializing level MultiFab
    MakeEBGeometry();
#endif

    InitBoundaryConditions();

    SetFluxScheme();

    SetAdvanceScheme();

    CheckDUGKSAmr();
}

DUGKSAmr::~DUGKSAmr()
{

}
// Make a new level from scratch using provided BoxArray and DistributionMapping.
// overrides the pure virtual function in AmrCore
//! Only used during initialization
void DUGKSAmr::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& new_grids,
                                        const DistributionMapping& new_dmap)
{
    BL_PROFILE("DUGKSAmr::MakeNewLevelFromScratch()");

    if (m_verbose > 0)
    {
        amrex::Print() << "Making new level " << lev << " from scratch" << std::endl;
        if (m_verbose > 2) {
            amrex::Print() << "with BoxArray " << new_grids << std::endl;
        }
    }

    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

#ifdef AMREX_USE_EB
    m_factory[lev] = makeEBFabFactory(geom[lev], new_grids, new_dmap,
                                      {m_eb_basic_grow_cells, 
                                       m_eb_volume_grow_cells, 
                                       m_eb_full_grow_cells},
                                 	    m_eb_support_level);
#else
    m_factory[lev].reset(new FArrayBoxFactory());
#endif

    m_leveldata[lev].reset(new LevelData(new_grids, new_dmap, *m_factory[lev], nghost_state()));

        m_t_new[lev] = time;
        m_t_old[lev] = time;

        // if (name_restart_file.empty()) {
            CaseInitFluid(lev);
        // }
}

Vector<MultiFab*> DUGKSAmr::LevelData::get_hDVDFnF_flux() noexcept
{
    Vector<MultiFab*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&hDVDFxF[Ark::FLUX]);
    r.push_back(&hDVDFyF[Ark::FLUX]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&hDVDFzF[Ark::FLUX]);
#endif
    return r;
}

Vector<MultiFab*> DUGKSAmr::LevelData::get_fDVDFnF_flux() noexcept
{
    Vector<MultiFab*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&fDVDFxF[Ark::FLUX]);
    r.push_back(&fDVDFyF[Ark::FLUX]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&fDVDFzF[Ark::FLUX]);
#endif
    return r;
}

Vector<MultiFab*> DUGKSAmr::LevelData::get_gDVDFnF_flux() noexcept
{
    Vector<MultiFab*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&gDVDFxF[Ark::FLUX]);
    r.push_back(&gDVDFyF[Ark::FLUX]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&gDVDFzF[Ark::FLUX]);
#endif
    return r;
}

Vector<MultiFab const*> DUGKSAmr::LevelData::get_hDVDFnF_eq_const() const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&hDVDFxF[Ark::EQF]);
    r.push_back(&hDVDFyF[Ark::EQF]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&hDVDFzF[Ark::EQF]);
#endif
    return r;
}
Vector<MultiFab const*> DUGKSAmr::LevelData::get_fDVDFnF_eq_const() const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&fDVDFxF[Ark::EQF]);
    r.push_back(&fDVDFyF[Ark::EQF]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&fDVDFzF[Ark::EQF]);
#endif
    return r;
}

Vector<MultiFab const*> DUGKSAmr::LevelData::get_gDVDFnF_eq_const() const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&gDVDFxF[Ark::EQF]);
    r.push_back(&gDVDFyF[Ark::EQF]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&gDVDFzF[Ark::EQF]);
#endif
    return r;
}
Vector<MultiFab const*> DUGKSAmr::LevelData::get_hDVDFnF_flux_const() const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&hDVDFxF[Ark::FLUX]);
    r.push_back(&hDVDFyF[Ark::FLUX]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&hDVDFzF[Ark::FLUX]);
#endif
    return r;
}
Vector<MultiFab const*> DUGKSAmr::LevelData::get_fDVDFnF_flux_const() const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&fDVDFxF[Ark::FLUX]);
    r.push_back(&fDVDFyF[Ark::FLUX]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&fDVDFzF[Ark::FLUX]);
#endif
    return r;
}
Vector<MultiFab const*> DUGKSAmr::LevelData::get_gDVDFnF_flux_const() const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&gDVDFxF[Ark::FLUX]);
    r.push_back(&gDVDFyF[Ark::FLUX]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&gDVDFzF[Ark::FLUX]);
#endif
    return r;
}
Vector<MultiFab const*> DUGKSAmr::LevelData::get_face_phase_const() const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&msv_xF[Ark::PHI]);
    r.push_back(&msv_yF[Ark::PHI]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&msv_zF[Ark::PHI]);
#endif
    return r;
}
Vector<MultiFab const*> DUGKSAmr::LevelData::get_face_density_const() const noexcept
{
    Vector<const MultiFab*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&msv_xF[Ark::RHO]);
    r.push_back(&msv_yF[Ark::RHO]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&msv_zF[Ark::RHO]);
#endif
    return r;
}
Vector<MultiFab const*> DUGKSAmr::LevelData::get_face_velocity_const() const noexcept
{
    Vector<const MultiFab*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&msv_xF[Ark::VEL]);
    r.push_back(&msv_yF[Ark::VEL]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&msv_zF[Ark::VEL]);
#endif
    return r;
}
Vector<MultiFab const*> DUGKSAmr::LevelData::get_face_pressure_const() const noexcept
{
    Vector<const MultiFab*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&msv_xF[Ark::P]);
    r.push_back(&msv_yF[Ark::P]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&msv_zF[Ark::P]);
#endif
    return r;
}
Vector<MultiFab const*> DUGKSAmr::LevelData::get_face_temperature_const() const noexcept
{
    Vector<const MultiFab*> r;
    r.reserve(AMREX_SPACEDIM);
//
    r.push_back(&msv_xF[Ark::T]);
    r.push_back(&msv_yF[Ark::T]);
#if(AMREX_SPACEDIM == 3)
    r.push_back(&msv_zF[Ark::T]);
#endif
    return r;
}