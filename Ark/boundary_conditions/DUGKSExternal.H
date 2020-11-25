#ifndef DUGKS_EXTERNAL_BCFUNC_H
#define DUGKS_EXTERNAL_BCFUNC_H

#include <DUGKSAmr.H>

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