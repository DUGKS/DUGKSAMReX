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

    const auto tagval = amrex::TagBox::SET;
//    const auto clearval = TagBox::CLEAR;

    const bool tag_gradphi = lev < grad_phi_error.size();

    const Real dlx = geom[lev].CellSize(0);
    const Real dly = geom[lev].CellSize(1);
#if (AMREX_SPACEDIM == 3)
    const Real dlz = geom[lev].CellSize(2);
#else
    const Real dlz = 0.;
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_leveldata[lev]->msv[Ark::PHI], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& tag = tags.array(mfi);

        if (tag_gradphi) 
        {
            Array4<Real const> const& grad_phi = 
            m_leveldata[lev]->msv[Ark::GradPHI].const_array(mfi);

            // Real phierr = tag_phi ? phi_error[lev]: std::numeric_limits<Real>::max();
            Real gradphierr = tag_gradphi 
                              ? grad_phi_error[lev] : std::numeric_limits<Real>::max();
            
            amrex::ParallelFor(bx,
            [tag_gradphi, grad_phi, gradphierr, tagval, tag]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // if (tag_phi and phi(i,j,k) > phierr) {
                //     tag(i,j,k) = tagval;
                // }
#if(AMREX_SPACEDIM == 3)
                Real const mod_grad_phi = std::sqrt( Ark::Square (grad_phi(i,j,k,Ark::X), 
                                                                  grad_phi(i,j,k,Ark::Y), 
                                                                  grad_phi(i,j,k,Ark::Z)));
                                                                  
                                                    
#else
                Real const mod_grad_phi = std::sqrt( Ark::Square (grad_phi(i,j,k,Ark::X), 
                                                                  grad_phi(i,j,k,Ark::Y))); 
#endif
                if (mod_grad_phi > gradphierr) {
                    tag(i,j,k) = tagval;
                }
            });
        } 

        if(tag_region)
        {
           Real xlo = tag_region_lo[0];
           Real xhi = tag_region_hi[0];

           Real ylo = tag_region_lo[1];
           Real yhi = tag_region_hi[1];

           Real zlo = tag_region_lo[2];
           Real zhi = tag_region_hi[2];

           amrex::ParallelFor(bx,
           [xlo, xhi, ylo, yhi, zlo, zhi, dlx, dly, dlz, tagval, tag]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
                Real x = (i+0.5)*dlx;
                Real y = (j+0.5)*dly;
                Real z = (k+0.5)*dlz;

                // Tag if we are inside the specified box
                if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
                {
                   tag(i,j,k) = tagval;
                }
           }); 
        }
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

void DUGKSAmr::ManualTagsPlacement(int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev)
{

}