#ifndef DUGKS_EXTERNAL_BCFUNC_H
#define DUGKS_EXTERNAL_BCFUNC_H

#include <DUGKSAmr.H>

struct FillExternalPhase;
struct FillExternalDensity;
struct FillExternalVelocity;
struct FillExternalPressure;
struct FillExternalTemperature;
struct FillExternalForce;
struct FillExternalChemiP;
struct FillExternalDVDFWave;
struct FillExternalDVDFBarPlus;

struct FillExternalEmpty;

struct FillExternalPhase {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        Ark::FillExternalNeumann(iv,dest,dcomp,numcomp,geom,time,bcr,bcomp,orig_comp);
    }
};

struct FillExternalDensity {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        Ark::FillExternalNeumann(iv,dest,dcomp,numcomp,geom,time,bcr,bcomp,orig_comp);
    }
};

struct FillExternalVelocity {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        
    }
};

struct FillExternalPressure {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        
    }
};

struct FillExternalTemperature {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        
    }
};

struct FillExternalForce {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        
    }
};

struct FillExternalChemicalP {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        // Ark::FillExternalDirichlet(iv,dest,dcomp,numcomp,geom,time,bcr,bcomp,orig_comp);
        // Ark::FillExternalNeumann(iv,dest,dcomp,numcomp,geom,time,bcr,bcomp,orig_comp);
    }
};

struct FillExternalDVDFWave {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        
    }
};

struct FillExternalDVDFBarPlus {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        
    }
};

struct FillExternalEmpty {
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        
    }
};
#endif

/* used for cpu boundaries, deprecated.
extern void FillExternalPhase (Box const& bx, Array4<Real> const& dest,
                        const int dcomp, const int numcomp,
                        GeometryData const& geom, const Real time,
                        const BCRec* bcr, const int bcomp,
                        const int orig_comp);

extern void FillExternalDensity (Box const& bx, Array4<Real> const& dest,
                          const int dcomp, const int numcomp,
                          GeometryData const& geom, const Real time,
                          const BCRec* bcr, const int bcomp,
                          const int orig_comp);

extern void FillExternalVelocity (Box const& bx, Array4<Real> const& dest,
                           const int dcomp, const int numcomp,
                           GeometryData const& geom, const Real time,
                           const BCRec* bcr, const int bcomp,
                           const int orig_comp);

extern void FillExternalPressure (Box const& bx, Array4<Real> const& dest,
                           const int dcomp, const int numcomp,
                           GeometryData const& geom, const Real time,
                           const BCRec* bcr, const int bcomp,
                           const int orig_comp);

extern void FillExternalTemperature (Box const& bx, Array4<Real> const& dest,
                              const int dcomp, const int numcomp,
                              GeometryData const& geom, const Real time,
                              const BCRec* bcr, const int bcomp,
                              const int orig_comp);

extern void FillExternalForce (Box const& bx, Array4<Real> const& dest,
                              const int dcomp, const int numcomp,
                              GeometryData const& geom, const Real time,
                              const BCRec* bcr, const int bcomp,
                              const int orig_comp);

extern void FillExternalDVDFWave (Box const& bx, Array4<Real> const& dest,
                           const int dcomp, const int numcomp,
                           GeometryData const& geom, const Real time,
                           const BCRec* bcr, const int bcomp,
                           const int orig_comp);

extern void FillExternalDVDFBarPlus (Box const& bx, Array4<Real> const& dest,
                              const int dcomp, const int numcomp,
                              GeometryData const& geom, const Real time,
                              const BCRec* bcr, const int bcomp,
                              const int orig_comp);
#endif
*/