#include <DUGKSAmr.H>
//
void DUGKSAmr::FillPhysicalBoundaryMacroVars(int const lev)
{
	for(int dir = 0; dir < AMREX_SPACEDIM;++dir)
	{
		FillPhysicalBoundaryMacroVars(lev, get_physical_bcrec()[0].lo(dir), dir);

		FillPhysicalBoundaryMacroVars(lev, get_physical_bcrec()[0].hi(dir), dir);
	}
}
void DUGKSAmr::FillPhysicalBoundaryMacroVars(int const lev, int const bcid, int const dir)
{
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
			FillPhysicalBoundaryMacroVarsNoSlipWall(lev,dir);
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
void DUGKSAmr::FillPhysicalBoundaryMacroVarsNoSlipWall(int const lev, int const dir)
{
	auto &ld = Data(lev);

	int const lo_min = Geom(lev).Domain().smallEnd(dir);
	int const hi_max = Geom(lev).Domain().bigEnd(dir) + 1;

  if(dir == Ark::X) 
  {

	#ifdef _OPENMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif
	for(amrex::MFIter mfi(ld.msv_xF[Ark::RHO], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Box &bx = mfi.tilebox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		Array4<Real const> const &phi  = ld.msv[Ark::PHI].const_array(mfi);
		Array4<Real const> const &rho  = ld.msv[Ark::RHO].const_array(mfi);
		Array4<Real const> const &pres = ld.msv[Ark::P].const_array(mfi);

		Array4<Real> const &phi_xF  = ld.msv_xF[Ark::PHI].array(mfi);
		Array4<Real> const &rho_xF  = ld.msv_xF[Ark::RHO].array(mfi);
		Array4<Real> const &pres_xF = ld.msv_xF[Ark::P].array(mfi);
		Array4<Real> const &vel_xF  = ld.msv_xF[Ark::VEL].array(mfi);

		if(lo.x == lo_min) 
		{				
			average_scalar_x (phi_xF,  phi,  lo, hi, lo_min);
			average_scalar_x (rho_xF,  rho,  lo, hi, lo_min);
			average_scalar_x (pres_xF, pres, lo, hi, lo_min);
			// set_dirichlet_pressure_xF (pres_xF, lo, hi, lo_min,
									  // Orientation(dir, Orientation::Side::low));

			set_dirichlet_velocity_xF (vel_xF, lo, hi, lo_min, 
									  Orientation(dir, Orientation::Side::low));
		}

		if(hi.x == hi_max) 
		{
			average_scalar_x (phi_xF,  phi,  lo, hi, hi_max);
			average_scalar_x (rho_xF,  rho,  lo, hi, hi_max);
			average_scalar_x (pres_xF, pres, lo, hi, hi_max);

			// set_dirichlet_pressure_xF (pres_xF, lo, hi, hi_max,
									  // Orientation(dir, Orientation::Side::high));

			set_dirichlet_velocity_xF (vel_xF, lo, hi, hi_max, 
									  Orientation(dir, Orientation::Side::high));
		}
	}
  }
  else if(dir == Ark::Y)
  {
  	#ifdef _OPENMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif
	for(amrex::MFIter mfi(ld.msv_yF[Ark::RHO], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Box &bx = mfi.tilebox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		Array4<Real const> const &phi  = ld.msv[Ark::PHI].const_array(mfi);
		Array4<Real const> const &rho  = ld.msv[Ark::RHO].const_array(mfi);
		Array4<Real const> const &pres = ld.msv[Ark::P].const_array(mfi);

		Array4<Real> const &phi_yF  = ld.msv_yF[Ark::PHI].array(mfi);
		Array4<Real> const &rho_yF  = ld.msv_yF[Ark::RHO].array(mfi);
		Array4<Real> const &pres_yF = ld.msv_yF[Ark::P].array(mfi);
		Array4<Real> const &vel_yF  = ld.msv_yF[Ark::VEL].array(mfi);

		if(lo.y == lo_min) 
		{				
			average_scalar_y(phi_yF,  phi,  lo, hi, lo_min);
			average_scalar_y(rho_yF,  rho,  lo, hi, lo_min);
			average_scalar_y(pres_yF, pres, lo, hi, lo_min);

			// set_dirichlet_pressure_yF (pres_yF, lo, hi, lo_min,
									  // Orientation(dir, Orientation::Side::low));

			set_dirichlet_velocity_yF (vel_yF, lo, hi, lo_min, 
									  Orientation(dir, Orientation::Side::low));
		}
		if(hi.y == hi_max)
		{
			average_scalar_y(phi_yF,  phi,  lo, hi, hi_max);
			average_scalar_y(rho_yF,  rho,  lo, hi, hi_max);
			average_scalar_y(pres_yF, pres, lo, hi, hi_max);

			// set_dirichlet_pressure_yF (pres_yF, lo, hi, hi_max, 
									  // Orientation(dir, Orientation::Side::high));

			set_dirichlet_velocity_yF (vel_yF, lo, hi, hi_max, 
									  Orientation(dir, Orientation::Side::high));
		}
	}
  }
#if(AMREX_SPACEDIM == 3)
  else if(dir == Ark::Z)
  {
  	#ifdef _OPENMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif
	for(amrex::MFIter mfi(ld.msv_zF[Ark::RHO], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Box &bx = mfi.tilebox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		Array4<Real const> const &phi  = ld.msv[Ark::PHI].const_array(mfi);
		Array4<Real const> const &rho  = ld.msv[Ark::RHO].const_array(mfi);
		Array4<Real const> const &pres = ld.msv[Ark::P].const_array(mfi);

		Array4<Real> const &phi_zF  = ld.msv_zF[Ark::PHI].array(mfi);
		Array4<Real> const &rho_zF  = ld.msv_zF[Ark::RHO].array(mfi);
		Array4<Real> const &pres_zF = ld.msv_zF[Ark::P].array(mfi);
		Array4<Real> const &vel_zF  = ld.msv_zF[Ark::VEL].array(mfi);

		if(lo.z == lo_min) 
		{				
			average_scalar_z(phi_zF,  phi,  lo, hi, lo_min);
			average_scalar_z(rho_zF,  rho,  lo, hi, lo_min);
			average_scalar_z(pres_zF, pres, lo, hi, lo_min);

			// set_dirichlet_pressure_zF (pres_zF, lo, hi, lo_min, 
			// 						  Orientation(dir, Orientation::Side::low));

			set_dirichlet_velocity_zF (vel_zF, lo, hi, lo_min, 
									  Orientation(dir, Orientation::Side::low));
		}
		if(hi.z == hi_max)
		{
			average_scalar_z( phi_zF, phi,  lo, hi, hi_max);
			average_scalar_z( rho_zF, rho,  lo, hi, hi_max);
			average_scalar_z(pres_zF, pres, lo, hi, hi_max);

			// set_dirichlet_pressure_zF (pres_zF, lo, hi, hi_max, 
			// 						  Orientation(dir, Orientation::Side::high));

			set_dirichlet_velocity_zF (vel_zF, lo, hi, hi_max, 
									  Orientation(dir, Orientation::Side::high));
		}
	}
  }
#endif
}
