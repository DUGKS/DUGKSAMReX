#include <DUGKSAmr.H>

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void DUGKSAmr::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
	BL_PROFILE("DUGKSAmr::MakeNewLevelFromCoarse()");

	if (m_verbose > 10) {
        amrex::Print() << "Making new level " << lev << " from coarse" << std::endl;
    }

#ifdef AMREX_USE_EB
    std::unique_ptr<FabFactory<FArrayBox> > new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                        {m_eb_basic_grow_cells, 
                                                                         m_eb_volume_grow_cells, 
                                                                         m_eb_full_grow_cells},
                                                                        m_eb_support_level);
#else
    std::unique_ptr<FabFactory<FArrayBox>> new_fact(new FArrayBoxFactory());
#endif
    std::unique_ptr<LevelData> new_leveldata(new LevelData(ba, dm, *new_fact, nghost_state()));

    FillCoarsePatchNew(lev, time, *new_leveldata);

    m_leveldata[lev] = std::move(new_leveldata);
    m_factory[lev] = std::move(new_fact);
}


// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void DUGKSAmr::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    BL_PROFILE("DUGKSAmr::RemakeLevel()");

    if (m_verbose > 10) {
        amrex::Print() << "Remaking level " << lev << std::endl;
    }

#ifdef AMREX_USE_EB
    std::unique_ptr<FabFactory<FArrayBox> > new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                        {m_eb_basic_grow_cells, 
                                                                         m_eb_volume_grow_cells, 
                                                                         m_eb_full_grow_cells},
                                                                        m_eb_support_level);
#else
    std::unique_ptr<FabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
#endif
    std::unique_ptr<LevelData> new_leveldata(new LevelData(ba, dm, *new_fact, nghost_state()));

    FillPatchNew(lev, time, *new_leveldata);
    PostFillPatchNew(lev,*new_leveldata);

    m_leveldata[lev] = std::move(new_leveldata);
    m_factory[lev] = std::move(new_fact);
}

// Delete level data
// overrides the pure virtual function in AmrCore
void DUGKSAmr::ClearLevel (int lev)
{
    BL_PROFILE("DUGKSAmr::ClearLevel()");
    m_leveldata[lev].reset();
    m_factory[lev].reset();
}