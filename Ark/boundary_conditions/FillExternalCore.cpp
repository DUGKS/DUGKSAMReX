#include <DUGKSAmr.H>

inline 
Real FillExternalContactAngle(const Real tangent_grad_hypotenuse,
                              const Real inner, const Real theta = 60)
{
    return inner + std::tan((0.5 - theta/180.0)*Ark::PI)*tangent_grad_hypotenuse;
    // return 1;
}

namespace amrex {

void FillExternalCore (const IntVect& iv, Array4<Real> const& q,
                       const int dcomp, const int numcomp,
                       GeometryData const& geom, const Real time,
                       const BCRec* bcr, const int bcomp,
                       const int orig_comp)
{
    const Box& domain_box = geom.Domain(); 
    const auto& domain_lo = domain_box.loVect();
    const auto& domain_hi = domain_box.hiVect();

    const Box qbx(q);
    const auto& qlo = qbx.loVect();
    const auto& qhi = qbx.hiVect();

    const int i = iv[0];
    const int iminus = i == qlo[0] ? qlo[0] : i-1;
    const int iplus  = i == qhi[0] ? qhi[0] : i+1;
    const int ilo = domain_lo[0];
    const int ihi = domain_hi[0];


#if AMREX_SPACEDIM >= 2
    const int j = iv[1];
    const int jminus = j == qlo[1] ? qlo[1] : j-1;
    const int jplus  = j == qhi[1] ? qhi[1] : j+1;
    const int jlo = domain_lo[1];
    const int jhi = domain_hi[1];
#endif

#if AMREX_SPACEDIM == 3
    const int k = iv[2];
    const int kminus = k == qlo[2] ? qlo[2] : k-1;
    const int kplus  = k == qhi[2] ? qhi[2] : k+1;
    const int klo = domain_lo[2];
    const int khi = domain_hi[2];
#else
    const int k = 0;
#endif

    for (int n = dcomp; n < numcomp+dcomp; ++n)
        {
            const BCRec& bc = bcr[bcomp+n-dcomp];

            if (i < ilo)
            {
                switch (bc.lo(0)) {
                case (amrex::BCType::ext_dir):
                {
#if AMREX_SPACEDIM == 2
                    const Real hypotenuse = std::fabs(q(ilo,jplus,k,n) - q(ilo,jminus,k,n));
#elif AMREX_SPACEDIM == 3
                    const Real hypotenuse = 
                    std::hypot(q(ilo,jplus,k    ,n) - q(ilo,jminus,k     ,n), 
                               q(ilo,j    ,kplus,n) - q(ilo,j     ,kminus,n));
#endif
                    q(i,j,k) = FillExternalContactAngle(hypotenuse,q(ilo+1,j,k,n));
                    break;
                }
                }
            }
            else if (i > ihi)
            {
                switch (bc.hi(0)) {
                case (amrex::BCType::ext_dir):
                {
#if AMREX_SPACEDIM == 2
                    const Real hypotenuse = std::fabs(q(ihi,jplus,k,n) - q(ihi,jminus,k,n));
#elif AMREX_SPACEDIM == 3
                    const Real hypotenuse = 
                    std::hypot(q(ihi,jplus,k    ,n) - q(ihi,jminus,k     ,n),
                               q(ihi,j    ,kplus,n) - q(ihi,j     ,kminus,n));
#endif
                    q(i,j,k) = FillExternalContactAngle(hypotenuse,q(ihi-1,j,k,n));
                    break;
                }
                }
            }

#if AMREX_SPACEDIM >= 2

            if (j < jlo)
            {
                switch (bc.lo(1)) {
                case (amrex::BCType::ext_dir):
                {
#if AMREX_SPACEDIM == 2
                    const Real hypotenuse = std::fabs(q(iplus,jlo,k,n) - q(iminus,jlo,k,n));
#elif AMREX_SPACEDIM == 3
                    const Real hypotenuse = 
                    std::hypot(q(iplus,jlo,k    ,n) - q(iminus,jlo,k     ,n), 
                               q(i    ,jlo,kplus,n) - q(i     ,jlo,kminus,n));
#endif
                    q(i,j,k) = FillExternalContactAngle(hypotenuse,q(i,jlo+1,k,n));
                    break;
                }
                }
            }
            else if (j > jhi)
            {
                switch (bc.hi(1)) {
                case (amrex::BCType::ext_dir):
                {
#if AMREX_SPACEDIM == 2
                    const Real hypotenuse = std::fabs(q(iplus,jhi,k,n) - q(iminus,jhi,k,n));
#elif AMREX_SPACEDIM == 3
                    const Real hypotenuse = 
                    std::hypot(q(iplus,jhi,k    ,n) - q(iminus,jhi,k     ,n),
                               q(i    ,jhi,kplus,n) - q(i     ,jhi,kminus,n));
#endif
                    q(i,j,k) = FillExternalContactAngle(hypotenuse,q(i,jhi-1,k,n));
                    break;
                }                
                }
            }
#endif

#if AMREX_SPACEDIM == 3

            if (k < klo)
            {
                switch (bc.lo(2)) {
                case (amrex::BCType::ext_dir):
                {
                    const Real hypotenuse = 
                    std::hypot(q(iplus,j    ,klo,n) - q(iminus,j     ,klo,n),
                               q(i    ,jplus,klo,n) - q(i     ,jminus,klo,n));
                    q(i,j,k) = FillExternalContactAngle(hypotenuse,q(i,j,klo+1,n));
                    break;
                }
                
                }
            }
            else if (k > khi)
            {
                switch (bc.hi(2)) {
                case (amrex::BCType::ext_dir):
                {
                    const Real hypotenuse = 
                    std::hypot(q(iplus,j,khi,n) - q(iminus,j,khi,n),
                               q(i,jplus,khi,n) - q(i,jminus,khi,n));

                    q(i,j,k) = FillExternalContactAngle(hypotenuse,q(i,j,khi-1,n));
                    break;
                }            
                }
            }
#endif
        }
}
}

/*void FillExternalPhase (Box const& bx, Array4<Real> const& dest,
                        const int dcomp, const int numcomp,
                        GeometryData const& geom, const Real time,
                        const BCRec* bcr, const int bcomp,
                        const int orig_comp)
{
    // do something for external Dirichlet (amrex::BCType::ext_dir) if there are
    
}
void FillExternalDensity (Box const& bx, Array4<Real> const& dest,
                          const int dcomp, const int numcomp,
                          GeometryData const& geom, const Real time,
                          const BCRec* bcr, const int bcomp,
                          const int orig_comp)
{
    // do something for external Dirichlet (amrex::BCType::ext_dir) if there are
}
void FillExternalVelocity (Box const& bx, Array4<Real> const& dest,
                           const int dcomp, const int numcomp,
                           GeometryData const& geom, const Real time,
                           const BCRec* bcr, const int bcomp,
                           const int orig_comp)
{
    // do something for external Dirichlet (amrex::BCType::ext_dir) if there are
}
void FillExternalPressure (Box const& bx, Array4<Real> const& dest,
                           const int dcomp, const int numcomp,
                           GeometryData const& geom, const Real time,
                           const BCRec* bcr, const int bcomp,
                           const int orig_comp)
{
    // do something for external Dirichlet (amrex::BCType::ext_dir) if there are
}
void FillExternalTemperature (Box const& bx, Array4<Real> const& dest,
                              const int dcomp, const int numcomp,
                              GeometryData const& geom, const Real time,
                              const BCRec* bcr, const int bcomp,
                              const int orig_comp)
{
    // do something for external Dirichlet (amrex::BCType::ext_dir) if there are
}
void FillExternalForce (Box const& bx, Array4<Real> const& dest,
                              const int dcomp, const int numcomp,
                              GeometryData const& geom, const Real time,
                              const BCRec* bcr, const int bcomp,
                              const int orig_comp)
{
  // do something for external Dirichlet (amrex::BCType::ext_dir) if there are
}                              
void FillExternalDVDFWave (Box const& bx, Array4<Real> const& dest,
                           const int dcomp, const int numcomp,
                           GeometryData const& geom, const Real time,
                           const BCRec* bcr, const int bcomp,
                           const int orig_comp)
{
    // do something for external Dirichlet (amrex::BCType::ext_dir) if there are
}
void FillExternalDVDFBarPlus (Box const& bx, Array4<Real> const& dest,
                              const int dcomp, const int numcomp,
                              GeometryData const& geom, const Real time,
                              const BCRec* bcr, const int bcomp,
                              const int orig_comp)
{
    // do something for external Dirichlet (amrex::BCType::ext_dir) if there are
}
*/