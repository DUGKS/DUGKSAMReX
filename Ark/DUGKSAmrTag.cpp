#include <DUGKSAmr.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#endif

// tag cells for refinement
// overrides the pure virtual function in AmrCore
void DUGKSAmr::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
	BL_PROFILE("DUGKSAmr::ErrorEst()");

#ifdef AMREX_USE_EB
    auto const& factory = EBFactory(lev);
    auto const& flags = factory.getMultiEBCellFlagFab();
#endif

    auto &ld = Data(lev);
    const Real *prob_lo = Geom(lev).ProbLo(); 
    const double *dl = Geom(lev).CellSize();

    const MultiFab &PhiFab     = ld.msv[Ark::PHI];
    const MultiFab &GradPhiFab = ld.msv[Ark::GradPHI];

    // const auto tagval = amrex::TagBox::SET;
    // const auto clearval = amrex::TagBox::CLEAR;

    const bool tag_scalar = lev < phi_error.size();
    const bool tag_vector = lev < grad_phi_error.size();
    const bool tag_designated_region = lev < tag_region.size();

    if(tag_scalar) {
        TagScalarField(tags, PhiFab, dl, phi_error[lev]);
    }

    if(tag_vector) {
        TagVectorField(tags, GradPhiFab, dl, grad_phi_error[lev]);
    }
    if(tag_designated_region and tag_region[lev]) {
        TagDesignatedRegion(tags, prob_lo, dl);
    }

#ifdef AMREX_USE_EB
    m_refine_cutcells = true;
    // Refine on cut cells
    if (m_refine_cutcells)
    {
        amrex::TagCutCells(tags, m_leveldata[lev]->velocity);
    }
#endif
}

void DUGKSAmr::TagScalarField (TagBoxArray& Tab, const MultiFab &ScalarFab, const Real *dx, 
                               const Real &tagError)
{
    const char tagval = static_cast<char>(amrex::TagBox::SET);
//    const auto clearval = amrex::TagBox::CLEAR;

    const Real dlx = dx[Ark::X];
    const Real dly = dx[Ark::Y];
#if (AMREX_SPACEDIM == 3)
    const Real dlz = dx[Ark::Z];
#else
    const Real dlz = 0.;
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(ScalarFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box &bx = mfi.tilebox();

        const auto &tagArray = Tab.array(mfi);
        Array4<Real const> const &scalar = ScalarFab.const_array(mfi);   

        amrex::ParallelFor(bx, 
        [scalar, tagError, tagval, tagArray] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if(scalar(i,j,k) > tagError) {
                tagArray(i,j,k) = tagval;
            }
        });
    }
}
void DUGKSAmr::TagVectorField (TagBoxArray& Tab, const MultiFab &VectorFab, const Real *dx, 
                               const Real &tagError)
{
    const auto tagval = amrex::TagBox::SET;
//    const auto clearval = amrex::TagBox::CLEAR;

    const Real dlx = dx[Ark::X];
    const Real dly = dx[Ark::Y];
#if (AMREX_SPACEDIM == 3)
    const Real dlz = dx[Ark::Z];
#else
    const Real dlz = 0.;
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(VectorFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box &bx = mfi.tilebox();

        const auto &tagArray = Tab.array(mfi);
        Array4<Real const> const &vec = VectorFab.const_array(mfi);   

        amrex::ParallelFor(bx, 
        [vec, tagError, tagval, tagArray] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if(AMREX_SPACEDIM == 3)
                Real const mod_vector = Ark::Hypot (vec(i,j,k,Ark::X),
                                                    vec(i,j,k,Ark::Y), 
                                                    vec(i,j,k,Ark::Z));
                                                          
#else
                Real const mod_vector = std::hypot(vec(i,j,k,Ark::X),vec(i,j,k,Ark::Y)); 
#endif
                if (mod_vector > tagError) {
                    tagArray(i,j,k) = tagval;
                }
        });
    }
}

void DUGKSAmr::TagDesignatedRegion(TagBoxArray& Tab, const Real *real_lo, const Real *dx)
{
    const auto tagval = amrex::TagBox::SET;
//    const auto clearval = amrex::TagBox::CLEAR;

    const Real dlx = dx[Ark::X];
    const Real dly = dx[Ark::Y];
    const Real real_lo_x = real_lo[Ark::X];
    const Real real_lo_y = real_lo[Ark::Y];
    
#if (AMREX_SPACEDIM == 3)
    const Real dlz = dx[Ark::Z];
    const Real real_lo_z = real_lo[Ark::Z];
#else
    const Real dlz = 0.;
    const Real real_lo_z = 0;
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Tab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box &bx = mfi.tilebox();

        const auto &tagArray = Tab.array(mfi);

        Real xlo = tag_region_lo[0];
        Real xhi = tag_region_hi[0];

        Real ylo = tag_region_lo[1];
        Real yhi = tag_region_hi[1];

        Real zlo = tag_region_lo[2];
        Real zhi = tag_region_hi[2];

        amrex::ParallelFor(bx,
        [xlo, xhi, ylo, yhi, zlo, zhi, real_lo_x, real_lo_y, real_lo_z, dlx, dly, dlz, tagArray, tagval]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+0.5)*dlx + real_lo_x;
            Real y = (j+0.5)*dly + real_lo_y;
            Real z = (k+0.5)*dlz + real_lo_z;

            // Tag if we are inside the specified box
            if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
            {
               tagArray(i,j,k) = tagval;
            }
        });
    }
}
void DUGKSAmr::ManualTagsPlacement(int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev)
{

}