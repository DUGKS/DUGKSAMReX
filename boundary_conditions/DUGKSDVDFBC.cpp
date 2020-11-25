#include <DUGKSAmr.H>

//non-equilibrium extrapolation
template <typename DVM>
void DUGKSBoundary<DVM>::FillPhysicalBoundaryFlux_NEE
(
	const Vector<MultiFab*> &xFluxFabPtr,
	const Vector<const MultiFab*> &xEqFabPtr,
	const Array<MultiFab, Ark::DVDFNUM> &xDVDF,
	int const lo_min, int const hi_max,
	Real const dTime, int const dir
)
{

AMREX_ASSERT(dir < AMREX_SPACEDIM);

AMREX_ASSERT(xFluxFabPtr.size() == AMREX_SPACEDIM);
AMREX_ASSERT(xEqFabPtr.size() == AMREX_SPACEDIM);

MultiFab &xFluxFabnF = *xFluxFabPtr[dir];
MultiFab const &xEqFabnF = *xEqFabPtr[dir];

Real const hTime = 0.5*dTime;

if(dir == Ark::X)
{

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for(amrex::MFIter mfi(xFluxFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Box &bx = mfi.tilebox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		//target
		Array4<Real 	 > const &xFlux_xF = xFluxFabnF.array(mfi);

		//source
		Array4<Real const> const &xEq_xF = xEqFabnF.const_array(mfi);
		Array4<Real const> const &xtau  = xDVDF[Ark::TAU].const_array(mfi);
		Array4<Real const> const &xWave = xDVDF[Ark::WAVE].const_array(mfi);
		Array4<Real const> const &xEq   = xDVDF[Ark::EQ].const_array(mfi);
		Array4<Real const> const &xSo   = xDVDF[Ark::SO].const_array(mfi);

		if(lo.x == lo_min)
		{
			int const ixC = lo_min; //cell index

			LoopVS(DVM::Q)
			LoopPS_YZ(lo,hi)
			{
				Real const xneq = xtau(ixC,j,k)/(xtau(ixC,j,k) + 0.5*hTime);

				xFlux_xF(lo_min,j,k,n) = xEq_xF(lo_min,j,k,n)
				+ xneq*(xWave(ixC,j,k,n) - xEq(ixC,j,k,n) + hTime*xSo(ixC,j,k,n));
			}
		}
		if(hi.x == hi_max)
		{
			int const ixC = hi_max - 1; //cell index

			LoopVS(DVM::Q)
			LoopPS_YZ(lo,hi)
			{
				Real const xneq = xtau(ixC,j,k)/(xtau(ixC,j,k) + hTime);

				xFlux_xF(hi_max,j,k,n) = xEq_xF(hi_max,j,k,n)
				+ xneq*(xWave(ixC,j,k,n) - xEq(ixC,j,k,n) + hTime*xSo(ixC,j,k,n));
			}
		}
	}
}
else if(dir == Ark::Y)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for(amrex::MFIter mfi(xFluxFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Box &bx = mfi.tilebox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		//target
		Array4<Real 	 > const &xFlux_yF = xFluxFabnF.array(mfi);

		//source
		Array4<Real const> const &xEq_yF = xEqFabnF.const_array(mfi);
		Array4<Real const> const &xtau  = xDVDF[Ark::TAU].const_array(mfi);
		Array4<Real const> const &xWave = xDVDF[Ark::WAVE].const_array(mfi);
		Array4<Real const> const &xEq   = xDVDF[Ark::EQ].const_array(mfi);
		Array4<Real const> const &xSo   = xDVDF[Ark::SO].const_array(mfi);

		if(lo.y == lo_min)
		{
			int const jyC = lo_min; //cell index

			LoopVS(DVM::Q)
			LoopPS_XZ(lo,hi)
			{
				Real const xneq = xtau(i,jyC,k)/(xtau(i,jyC,k) + hTime);

				xFlux_yF(i,lo_min,k,n) = xEq_yF(i,lo_min,k,n)
				+ xneq*(xWave(i,jyC,k,n) - xEq(i,jyC,k,n) + hTime*xSo(i,jyC,k,n));
				;
			}
		}
		if(hi.y == hi_max)
		{
			int const jyC = hi_max - 1; //cell index

			LoopVS(DVM::Q)
			LoopPS_XZ(lo,hi)
			{
				Real const xneq = xtau(i,jyC,k)/(xtau(i,jyC,k) + hTime);

				xFlux_yF(i,hi_max,k,n) = xEq_yF(i,hi_max,k,n)
				+ xneq*(xWave(i,jyC,k,n) - xEq(i,jyC,k,n) + hTime*xSo(i,jyC,k,n));
				;
			}
		}
	}	
}
#if(AMREX_SPACEDIM == 3)
else if(dir == Ark::Z)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for(amrex::MFIter mfi(xFluxFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Box &bx = mfi.tilebox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		//target
		Array4<Real 	 > const &xFlux_zF = xFluxFabnF.array(mfi);

		//source
		Array4<Real const> const &xEq_zF = xEqFabnF.const_array(mfi);
		Array4<Real const> const &xtau  = xDVDF[Ark::TAU].const_array(mfi);
		Array4<Real const> const &xWave = xDVDF[Ark::WAVE].const_array(mfi);
		Array4<Real const> const &xEq   = xDVDF[Ark::EQ].const_array(mfi);
		Array4<Real const> const &xSo   = xDVDF[Ark::SO].const_array(mfi);

		if(lo.z == lo_min)
		{
			int const kzC = lo_min; //cell index

			LoopVS(DVM::Q)
			LoopPS_XY(lo,hi)
			{
				Real const xneq = xtau(i,j,kzC)/(xtau(i,j,kzC) + hTime);

				xFlux_zF(i,j,lo_min,n) = xEq_zF(i,j,lo_min,n)
				+ xneq*(xWave(i,j,kzC,n) - xEq(i,j,kzC,n) + hTime*xSo(i,j,kzC,n));
			}
		}
		if(hi.z == hi_max)
		{
			int const kzC = hi_max - 1; //cell index

			LoopVS(DVM::Q)
			LoopPS_XY(lo,hi)
			{
				Real const xneq = xtau(i,j,kzC)/(xtau(i,j,kzC) + hTime);

				xFlux_zF(i,j,hi_max,n) = xEq_zF(i,j,hi_max,n)
				+ xneq*(xWave(i,j,kzC,n) - xEq(i,j,kzC,n) + hTime*xSo(i,j,kzC,n));
			}
		}
	}	
}
#endif
	else {
		std::ostringstream ss;   
        ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
        ss  <<"Invalid face direction\n";
        amrex::Abort(ss.str());	
	}
}
template <typename DVM>
void DUGKSBoundary<DVM>::BounceBack (Array4<Real> const &xFlux,
                 				 	 Array4<const Real> const &var, 
                 				 	 Array4<const Real> const &vel,
                 				 	 int const i, int const j, int const k, int const n)
{
	xFlux(i,j,k,n) = xFlux(i,j,k,DVM::BB[n]) 
			   		+ 2.0*var(i,j,k)*DVM::omega[n]
			   		*(
				   		DVM::xi_u[n]*vel(i,j,k,Ark::X)
				   	 +	DVM::xi_v[n]*vel(i,j,k,Ark::Y)	
				   	 #if(AMREX_SPACEDIM == 3)
				   	 +	DVM::xi_w[n]*vel(i,j,k,Ark::Z)
				   	 #endif	
				   	 )/DVM::RT;
#ifdef AMREX_DEBUG
	if(vel(i,j,k,Ark::X) != 0 or vel(i,j,k,Ark::Y) != 0)
	amrex::Print()  << vel(i,j,k,Ark::X) <<"    "
					<< vel(i,j,k,Ark::Y) <<"\n";
#endif
}
//
template <typename DVM>
void DUGKSBoundary<DVM>::FillPhysicalBoundaryFlux_BB (const Vector<MultiFab*> &xFluxFabPtr, 
									  				  const Vector<const MultiFab*> &xVarFabPtr,
									  				  const Vector<const MultiFab*> &VelFabPtr,
									  				  int const lo_min, int const hi_max,
									  				  int const dir)
{

AMREX_ASSERT(dir < AMREX_SPACEDIM);

AMREX_ASSERT(xFluxFabPtr.size() == AMREX_SPACEDIM);
AMREX_ASSERT(xVarFabPtr.size() == AMREX_SPACEDIM);
AMREX_ASSERT(VelFabPtr.size() == AMREX_SPACEDIM);

MultiFab& xFluxFab      = *xFluxFabPtr[dir];
const MultiFab& xVarFab = *xVarFabPtr[dir];
const MultiFab& VelFab  = *VelFabPtr[dir];

AMREX_ASSERT(xFluxFab.boxArray() == xVarFab.boxArray());
AMREX_ASSERT(xFluxFab.boxArray() == VelFab.boxArray());
//
AMREX_ASSERT(xFluxFab.DistributionMap() == xVarFab.DistributionMap());
AMREX_ASSERT(xFluxFab.DistributionMap() == VelFab.DistributionMap());


if(dir == Ark::X)
{

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for(amrex::MFIter mfi(xFluxFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Box &bx = mfi.tilebox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		//target
		Array4<Real 	 > const &xFlux_xF = xFluxFab.array(mfi);

		//source
		Array4<Real const> const &xvar_xF = xVarFab.const_array(mfi);
		Array4<Real const> const &vel_xF  = VelFab.const_array(mfi);

		if(lo.x == lo_min)
		{
			LoopVS(DVM::Q)
			{
				if(DVM::xi_u[n] > 0)
				{
					LoopPS_YZ(lo,hi)
					{
						BounceBack(xFlux_xF,xvar_xF,vel_xF,lo_min,j,k,n);
					}
				}
			}
		}
		if(hi.x == hi_max)
		{
			LoopVS(DVM::Q)
			{
				if(DVM::xi_u[n] < 0)
				{
					LoopPS_YZ(lo,hi)
					{
						BounceBack(xFlux_xF,xvar_xF,vel_xF,hi_max,j,k,n);
					}
				}
			}
		}
	}
}
else if(dir == Ark::Y)
{

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for(amrex::MFIter mfi(xFluxFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Box &bx = mfi.tilebox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		//target
		Array4<Real 	 > const &xFlux_yF = xFluxFab.array(mfi);

		//source
		Array4<Real const> const &xvar_yF = xVarFab.const_array(mfi);
		Array4<Real const> const &vel_yF  = VelFab.const_array(mfi);

		if(lo.y == lo_min)
		{
			LoopVS(DVM::Q)
			{
				if(DVM::xi_v[n] > 0)
				{
					LoopPS_XZ(lo,hi)
					{
						BounceBack(xFlux_yF,xvar_yF,vel_yF,i,lo_min,k,n);
					}
				}
			}
		}
		if(hi.y == hi_max)
		{
			LoopVS(DVM::Q)
			{
				if(DVM::xi_v[n] < 0)
				{
					LoopPS_XZ(lo,hi)
					{
						BounceBack(xFlux_yF,xvar_yF,vel_yF,i,hi_max,k,n);
					}
				}
			}
		}
	}
}
#if(AMREX_SPACEDIM == 3)
else if(dir == Ark::Z)
{

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for(amrex::MFIter mfi(xFluxFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Box &bx = mfi.tilebox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		//target
		Array4<Real 	 > const &xFlux_zF = xFluxFab.array(mfi);

		//source
		Array4<Real const> const &xvar_zF = xVarFab.const_array(mfi);
		Array4<Real const> const &vel_zF  = VelFab.const_array(mfi);

		if(lo.z == lo_min)
		{
			LoopVS(DVM::Q)
			{
				if(DVM::xi_w[n] > 0)
				{
					LoopPS_XY(lo,hi)
					{
						BounceBack(xFlux_zF,xvar_zF,vel_zF,i,j,lo_min,n);
					}
				}
			}
		}
		if(hi.z == hi_max)
		{
			LoopVS(DVM::Q)
			{
				if(DVM::xi_w[n] < 0)
				{
					LoopPS_XY(lo,hi)
					{
						BounceBack(xFlux_zF,xvar_zF,vel_zF,i,j,hi_max,n);
					}
				}
			}
		}
	}
}
#endif
else{
	std::ostringstream ss;   
    ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
    ss  <<"Invalid face direction\n";
    amrex::Abort(ss.str());	
}
}
void DUGKSAmr::FillPhysicalBoundaryDVDF(int const lev)
{
	for(int dir = 0; dir < AMREX_SPACEDIM;++dir)
	{
		FillPhysicalBoundaryDVDF(lev, get_physical_bcrec()[0].lo(dir), dir);

		FillPhysicalBoundaryDVDF(lev, get_physical_bcrec()[0].hi(dir), dir);
	}
}
void DUGKSAmr::FillPhysicalBoundaryDVDF(int const lev, int const bcid, int const dir)
{
  	auto &ld = Data(lev);
  	
	int const lo_min_face = Geom(lev).Domain().smallEnd(dir);
	int const hi_max_face = Geom(lev).Domain().bigEnd(dir) + 1;

	switch(bcid)
	{
		case amrex::PhysBCType::interior :
		{
			//do nothing
			break;
		}
		case amrex::PhysBCType::inflow :
		{
			//to be done...
			amrex::Abort("inflow boundary has not finished");
			break;
		}
		case amrex::PhysBCType::outflow :
		{
			//HOEXTRAPCC
			break;
		}
		case amrex::PhysBCType::symmetry :
		{
			//to be done... may be never
			amrex::Abort("symmetry boundary has not finished");
			break;
		}
		case amrex::PhysBCType::slipwall :
		{
			//to be done... may be never
			amrex::Abort("slipwall boundary has not finished");
			break;
		}
		case amrex::PhysBCType::noslipwall :
		{
#ifdef ARK_USE_PHASE
DUGKSBoundary<Ark::DmVnPhase>::FillPhysicalBoundaryFlux_NEE (ld.get_hDVDFnF_flux(),
													   	 	 ld.get_hDVDFnF_eq_const(),
													   		 ld.hDVDF, 
													   		 lo_min_face, hi_max_face,
													   		 Dt(), dir);

DUGKSBoundary<Ark::DmVnPhase>::FillPhysicalBoundaryFlux_BB (ld.get_hDVDFnF_flux(), 
									  				  		ld.get_face_phase_const(),
									  				  		ld.get_face_velocity_const(),
									  				  		lo_min_face, hi_max_face, dir);
#endif

#ifdef ARK_USE_MOMENTUM
DUGKSBoundary<Ark::DmVnMomentum>::FillPhysicalBoundaryFlux_NEE (ld.get_fDVDFnF_flux(),
													   	 	    ld.get_fDVDFnF_eq_const(),
													   		 	ld.fDVDF, 
													   		 	lo_min_face, hi_max_face,
													   		 	Dt(), dir);

DUGKSBoundary<Ark::DmVnMomentum>::FillPhysicalBoundaryFlux_BB (ld.get_fDVDFnF_flux(), 
									  				  		   ld.get_face_density_const(),
									  				  		   ld.get_face_velocity_const(),
									  				  		   lo_min_face, hi_max_face, dir);
#endif

#ifdef ARK_USE_ENERGY
DUGKSBoundary<Ark::DmVnEnergy>::FillPhysicalBoundaryFlux_NEE (ld.get_gDVDFnF_flux(),
													   	 	  ld.get_gDVDFnF_eq_const(),
													   		  ld.gDVDF, 
													   		  lo_min_face, hi_max_face,
													   		  Dt(), dir);

DUGKSBoundary<Ark::DmVnEnergy>::FillPhysicalBoundaryFlux_BB (ld.get_gDVDFnF_flux(), 
									  				  		 ld.get_face_temperature_const(),
									  				  		 ld.get_face_velocity_const(),
									  				  		 lo_min_face, hi_max_face, dir);
#endif
			break;
		}	
		default :
    	{
    		std::ostringstream ss;   
       	    ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
       	    ss  <<"Invalid physical boundary\n";
       	    amrex::Abort(ss.str());
    	}
	}
}
void DUGKSAmr::FillInteriorBoundaryGradBarPlus (Array<MultiFab, Ark::DVDFNUM> &xDVDF,
											   const Periodicity& period)
{
	for(int n = Ark::BPx;n < Ark::DVDFNUM;++n)
	{
		xDVDF[n].FillBoundary(period);
	}
// 	xDVDF[Ark::BPx].FillBoundary(period);
// 	xDVDF[Ark::BPy].FillBoundary(period);
// 	xDVDF[Ark::BPxx].FillBoundary(period);
// 	xDVDF[Ark::BPyy].FillBoundary(period);
// 	xDVDF[Ark::BPxy].FillBoundary(period);
// #if(AMREX_SPACEDIM == 3)
// 	xDVDF[Ark::BPz].FillBoundary(period);
// 	xDVDF[Ark::BPzz].FillBoundary(period);
// 	xDVDF[Ark::BPyz].FillBoundary(period);
// 	xDVDF[Ark::BPzx].FillBoundary(period);
// #endif
}