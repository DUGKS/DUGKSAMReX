#ifndef ARK_DUGKS_CORE_IMPLEMENTATION_H
#define ARK_DUGKS_CORE_IMPLEMENTATION_H

template<typename DVM>
void DUGKSCore<DVM>::Update_BarPlus(Array<MultiFab, Ark::DVDFNUM> &xDVDF, 
							   		Real const hTime, Real const dTime)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xDVDF[Ark::BP], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xTau  = xDVDF[Ark::TAU].const_array(mfi);
        Array4<Real const> const &xWave = xDVDF[Ark::WAVE].const_array(mfi);
        Array4<Real const> const &xEq   = xDVDF[Ark::EQ].const_array(mfi);
        Array4<Real const> const &xSo   = xDVDF[Ark::SO].const_array(mfi);
        Array4<Real>       const &xBP   = xDVDF[Ark::BP].array(mfi);

        amrex::ParallelFor(bx, DVM::Q, 
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	const Real xtau = xTau(i,j,k);
        	const Real b = (hTime + dTime)/(2.0*xtau + dTime);
        	// const Real c = b*xTau(i,j,k);

        	xBP(i,j,k,n) = xWave(i,j,k,n)
        				 + b*(xEq(i,j,k,n) - xWave(i,j,k,n) + xtau*xSo(i,j,k,n));
        });

  //       const auto lo = lbound(bx);
		// const auto hi = ubound(bx);
  //       LoopPS(lo,hi)
  //       {
  //       	const Real a = (2.0*xtau(i,j,k) - htime)/(2.0*xtau(i,j,k) + dtime);
  //       	const Real b = 1.0 - a;
  //       	const Real c = b*xtau(i,j,k);
  //       	LoopVS(DVM::Q)
  //       	{
  //       		xBP(i,j,k,n) = a*xWave(i,j,k,n) + b*xEq(i,j,k,n) + c*xSo(i,j,k,n);
  //       	}
  //       }
    }
}

template<typename DVM>
void DUGKSCore<DVM>::Update_Grad_BarPlus(Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real* dx)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xDVDF[Ark::BP], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.growntilebox(1);

		Array4<Real> const &xBPx    = xDVDF[Ark::BPx].array(mfi);
		Array4<Real> const &xBPy    = xDVDF[Ark::BPy].array(mfi);
		Array4<Real> const &xBPxx   = xDVDF[Ark::BPxx].array(mfi);
		Array4<Real> const &xBPyy   = xDVDF[Ark::BPyy].array(mfi);
		Array4<Real> const &xBPxy   = xDVDF[Ark::BPxy].array(mfi);

#if(AMREX_SPACEDIM == 3)
		Array4<Real> const &xBPz  = xDVDF[Ark::BPz].array(mfi);
		Array4<Real> const &xBPzz = xDVDF[Ark::BPzz].array(mfi);
		Array4<Real> const &xBPyz = xDVDF[Ark::BPyz].array(mfi);
		Array4<Real> const &xBPzx = xDVDF[Ark::BPzx].array(mfi);
#endif
		
		Array4<Real  const> const &xBP = xDVDF[Ark::BP].const_array(mfi);

		GradFunc::Update_Grad_x(   bx, DVM::Q, dx, xBP, xBPx);
		GradFunc::Update_Grad_y(   bx, DVM::Q, dx, xBP, xBPy);
		//
		GradFunc::Update_Grad_xx(  bx, DVM::Q, dx, xBP, xBPxx);
		GradFunc::Update_Grad_yy(  bx, DVM::Q, dx, xBP, xBPyy);
		GradFunc::Update_Grad_xy(  bx, DVM::Q, dx, xBP, xBPxy);
#if(AMREX_SPACEDIM == 3)
		GradFunc::Update_Grad_z(   bx, DVM::Q, dx, xBP, xBPz);
		GradFunc::Update_Grad_zz(  bx, DVM::Q, dx, xBP, xBPzz);
		GradFunc::Update_Grad_yz(  bx, DVM::Q, dx, xBP, xBPyz);
		GradFunc::Update_Grad_zx(  bx, DVM::Q, dx, xBP, xBPzx);
#endif
	}
}
template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_UW (MultiFab &xBarhFabnF, 
									 Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
									 const Real* dx, const Real hTime, 
									 amrex::Direction FaceID)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xBarhFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real> const &xBh = xBarhFabnF.array(mfi);

		Array4<Real const> const &xBP   = xDVDF[Ark::BP].const_array(mfi);
		Array4<Real const> const &xBPx  = xDVDF[Ark::BPx].const_array(mfi);
		Array4<Real const> const &xBPy  = xDVDF[Ark::BPy].const_array(mfi);
		Array4<Real const> const &xBPxx = xDVDF[Ark::BPxx].const_array(mfi);
		Array4<Real const> const &xBPyy = xDVDF[Ark::BPyy].const_array(mfi);
		Array4<Real const> const &xBPxy = xDVDF[Ark::BPxy].const_array(mfi);

#if(AMREX_SPACEDIM == 3)
		Array4<Real const> const &xBPz  = xDVDF[Ark::BPz].const_array(mfi);
		Array4<Real const> const &xBPzz = xDVDF[Ark::BPzz].const_array(mfi);
		Array4<Real const> const &xBPyz = xDVDF[Ark::BPyz].const_array(mfi);
		Array4<Real const> const &xBPzx = xDVDF[Ark::BPzx].const_array(mfi);
#endif

		if(FaceID == amrex::Direction::x){

		const Real delta = dx[Ark::X];

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_u[n]);
			int const multiplier = token < 0 ? -1 : 1;
			int const di = token < 0 ? 0 : -1;
			//
			const Real dlx  = 0.5*delta*multiplier - hTime*DVM::xi_u[n];
			const Real dlx2 = dlx*dlx;	
			//
			const Real dly  = -hTime*DVM::xi_v[n];
			const Real dly2 = dly*dly;	

			const Real dlxy = dlx*dly; 

#if(AMREX_SPACEDIM == 3)
			const Real dlz  = -hTime*DVM::xi_w[n];
			const Real dlz2 = dlz*dlz;
			const Real dlyz = dly*dlz;
			const Real dlzx = dlz*dlx;
#endif
			
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const ii = i + di;

#if(AMREX_SPACEDIM == 2)
				xBh(i,j,k,n) = xBP(ii,j,k,n) 
							 + xBPx(ii,j,k,n)*dlx + xBPy(ii,j,k,n)*dly
							 + (xBPxx(ii,j,k,n)*dlx2 + xBPyy(ii,j,k,n)*dly2)*0.5
							 + xBPxy(ii,j,k,n)*dlxy;
#endif

#if(AMREX_SPACEDIM == 3)
				xBh(i,j,k,n) = xBP(ii,j,k,n) 
							 + xBPx(ii,j,k,n)*dlx
							 + xBPy(ii,j,k,n)*dly
							 + xBPz(ii,j,k,n)*dlz
							 + (xBPxx(ii,j,k,n)*dlx2 
							   +xBPyy(ii,j,k,n)*dly2 
							   +xBPzz(ii,j,k,n)*dlz2)*0.5							 
							 + xBPxy(ii,j,k,n)*dlxy
							 + xBPyz(ii,j,k,n)*dlyz
							 + xBPzx(ii,j,k,n)*dlzx;			 
#endif							 

			});
			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
#if(AMREX_SPACEDIM == 2)
				Real tmp = xBP(i,j,k,n) 
						 - xBPx(i,j,k,n)*dlx + xBPy(i,j,k,n)*dly
						 + (xBPxx(i,j,k,n)*dlx2 + xBPyy(i,j,k,n)*dly2)*0.5
						 - xBPxy(i,j,k,n)*dlxy;
#endif

#if(AMREX_SPACEDIM == 3)
				Real tmp = xBP(i,j,k,n) 
						 - xBPx(i,j,k,n)*dlx 
						 + xBPy(i,j,k,n)*dly 
						 + xBPz(i,j,k,n)*dlz
						 + (xBPxx(i,j,k,n)*dlx2
						   +xBPyy(i,j,k,n)*dly2
						   +xBPzz(i,j,k,n)*dlz2)*0.5
						 - xBPxy(i,j,k,n)*dlxy 
						 + xBPyz(i,j,k,n)*dlyz 
						 - xBPzx(i,j,k,n)*dlzx;		 
#endif						 
				xBh(i,j,k,n) += tmp;
				xBh(i,j,k,n) /= 2;
			});
			}		
		}
		}
		else if(FaceID == amrex::Direction::y){

		const Real delta = dx[Ark::Y];

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_v[n]);
			int const multiplier = token < 0 ? -1 : 1;
			int const dj = token < 0 ? 0 : -1;
			//
			const Real dly  = 0.5*delta*multiplier - hTime*DVM::xi_v[n];
			const Real dly2 = dly*dly;
			//
			const Real dlx  = -hTime*DVM::xi_u[n];
			const Real dlx2 = dlx*dlx;	
			//
			const Real dlxy = dlx*dly;

#if(AMREX_SPACEDIM == 3)
			const Real dlz  = -hTime*DVM::xi_w[n];
			const Real dlz2 = dlz*dlz;
			const Real dlyz = dly*dlz;
			const Real dlzx = dlz*dlx;
#endif

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const jj = j + dj;
#if(AMREX_SPACEDIM == 2)
				xBh(i,j,k,n) = xBP(i,jj,k,n) 
							 + xBPx(i,jj,k,n)*dlx + xBPy(i,jj,k,n)*dly
							 + (xBPxx(i,jj,k,n)*dlx2 + xBPyy(i,jj,k,n)*dly2)*0.5
							 + xBPxy(i,jj,k,n)*dlxy;
#endif

#if(AMREX_SPACEDIM == 3)
				xBh(i,j,k,n) = xBP(i,jj,k,n) 
							 + xBPx(i,jj,k,n)*dlx 
							 + xBPy(i,jj,k,n)*dly 
							 + xBPz(i,jj,k,n)*dlz
							 + (xBPxx(i,jj,k,n)*dlx2 
							   +xBPyy(i,jj,k,n)*dly2
							   +xBPzz(i,jj,k,n)*dlz2)*0.5
							 + xBPxy(i,jj,k,n)*dlxy
							 + xBPyz(i,jj,k,n)*dlyz
							 + xBPzx(i,jj,k,n)*dlzx;			 
#endif
			});

			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
#if(AMREX_SPACEDIM == 2)					
				Real tmp = xBP(i,j,k,n) 
						 + xBPx(i,j,k,n)*dlx - xBPy(i,j,k,n)*dly
						 + (xBPxx(i,j,k,n)*dlx2 + xBPyy(i,j,k,n)*dly2)*0.5
						 - xBPxy(i,j,k,n)*dlxy;
#endif

#if(AMREX_SPACEDIM == 3)
				Real tmp = xBP(i,j,k,n) 
						 + xBPx(i,j,k,n)*dlx 
						 - xBPy(i,j,k,n)*dly
						 + xBPz(i,j,k,n)*dlz
						 + (xBPxx(i,j,k,n)*dlx2 
						   +xBPyy(i,j,k,n)*dly2
						   +xBPzz(i,j,k,n)*dlz2)*0.5
						 - xBPxy(i,j,k,n)*dlxy
						 - xBPyz(i,j,k,n)*dlyz
						 + xBPzx(i,j,k,n)*dlzx;		 
#endif							 
				xBh(i,j,k,n) += tmp;
				xBh(i,j,k,n) /= 2.0;
			});
			}
		}
		}
#if(AMREX_SPACEDIM == 3)
		else if(FaceID == amrex::Direction::z){

		const Real delta = dx[Ark::Z];

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_w[n]);
			int const multiplier = token < 0 ? -1 : 1;
			int const dk = token < 0 ? 0 : -1;
			//
			const Real dlz  = 0.5*delta*multiplier - hTime*DVM::xi_w[n];
			const Real dlz2 = dlz*dlz;
			//
			const Real dlx  = -hTime*DVM::xi_u[n];
			const Real dlx2 = dlx*dlx;
			//
			const Real dly  = -hTime*DVM::xi_v[n];
			const Real dly2 = dly*dly;

			const Real dlxy = dlx*dly;
			const Real dlyz = dly*dlz;
			const Real dlzx = dlz*dlx;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const kk = k + dk;

				xBh(i,j,k,n) = xBP(i,j,kk,n) 
							 + xBPx(i,j,kk,n)*dlx 
							 + xBPy(i,j,kk,n)*dly 
							 + xBPz(i,j,kk,n)*dlz
							 + (xBPxx(i,j,kk,n)*dlx2 
							   +xBPyy(i,j,kk,n)*dly2
							   +xBPzz(i,j,kk,n)*dlz2)*0.5
							 + xBPxy(i,j,kk,n)*dlxy
							 + xBPyz(i,j,kk,n)*dlyz
							 + xBPzx(i,j,kk,n)*dlzx;
			});

			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				Real tmp = xBP(i,j,k,n) 
						 + xBPx(i,j,k,n)*dlx 
						 + xBPy(i,j,k,n)*dly
						 - xBPz(i,j,k,n)*dlz
						 + (xBPxx(i,j,k,n)*dlx2 
						   +xBPyy(i,j,k,n)*dly2
						   +xBPzz(i,j,k,n)*dlz2)*0.5
						 + xBPxy(i,j,k,n)*dlxy
						 - xBPyz(i,j,k,n)*dlyz
						 - xBPzx(i,j,k,n)*dlzx;		 

				xBh(i,j,k,n) += tmp;
				xBh(i,j,k,n) /= 2.0;
			});
			}
		}
		}			
#endif
		else {
			amrex::Abort("UW : How did we get here?");
		}
		
	}
}
template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_CD (MultiFab &xBarhFab, const MultiFab &xBarPlusFab,
									 const Real* dx, const Real hTime, 
									 amrex::Direction FaceID)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xBarhFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Real dlx = dx[Ark::X];
		const Real dly = dx[Ark::Y];
#if(AMREX_SPACEDIM == 3)
		const Real dlz = dx[Ark::Z];
#endif
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xBP = xBarPlusFab.const_array(mfi);
		Array4<Real      > const &xBh = xBarhFab.array(mfi);

		if(FaceID == amrex::Direction::x) {

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
			Real xBP_x = (xBP(i,j,k,n) - xBP(i-1,j,k,n))/dlx;

			Real xBP_y = (xBP(i,  j+1,k,n) - xBP(i,  j-1,k,n) 
					     +xBP(i-1,j+1,k,n) - xBP(i-1,j-1,k,n)
					     )*0.25/dly;

#if(AMREX_SPACEDIM == 2)
			xBh(i,j,k,n) = 0.5*(xBP(i,j,k,n) + xBP(i-1,j,k,n)) 
        			     - hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y);
#endif
#if(AMREX_SPACEDIM == 3)
			Real xBP_z = (xBP(i,  j,k+1,n) - xBP(i,  j,k-1,n)
        				 +xBP(i-1,j,k+1,n) - xBP(i-1,j,k-1,n)
        				 )*0.25/dlz;

			xBh(i,j,k,n) = 0.5*(xBP(i,j,k,n) + xBP(i-1,j,k,n))
        			     - hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y + DVM::xi_w[n]*xBP_z);
#endif
		});
		}

		else if(FaceID == amrex::Direction::y){

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	Real xBP_x = (xBP(i+1,j,  k,n) - xBP(i-1,j,  k,n)
        				 +xBP(i+1,j-1,k,n) - xBP(i-1,j-1,k,n)
        				 )*0.25/dlx;

        	Real xBP_y = (xBP(i,j,k,n) - xBP(i,j-1,k,n))/dly;

#if(AMREX_SPACEDIM == 2)
        	xBh(i,j,k,n) = 0.5*(xBP(i,j,k,n) + xBP(i,j-1,k,n)) 
        			     - hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y);
#endif
#if(AMREX_SPACEDIM == 3)
        	Real xBP_z = (xBP(i,j,  k+1,n) - xBP(i,j,  k-1,n)
        				 +xBP(i,j-1,k+1,n) - xBP(i,j-1,k-1,n)
        				 )*0.25/dlz;

        	xBh(i,j,k,n) = 0.5*(xBP(i,j,k,n) + xBP(i,j-1,k,n)) 
        			     - hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y + DVM::xi_w[n]*xBP_z);
#endif        			    
        });
		}
#if(AMREX_SPACEDIM == 3)
		else if(FaceID == amrex::Direction::z){

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
			Real xBP_x = (xBP(i+1,j,k,  n) - xBP(i-1,j,k,  n) 
					     +xBP(i+1,j,k-1,n) - xBP(i-1,j,k-1,n)
					     )*0.25/dlx;

			Real xBP_y = (xBP(i,j+1,k,  n) - xBP(i,j-1,k,  n) 
					     +xBP(i,j+1,k-1,n) - xBP(i,j-1,k-1,n)
					     )*0.25/dly;

			Real xBP_z = (xBP(i,j,k,n) - xBP(i,j,k-1,n))/dlz;

			xBh(i,j,k,n) = 0.5*(xBP(i,j,k,n) + xBP(i,j,k-1,n)) 
        			     - hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y + DVM::xi_w[n]*xBP_z);
		});
		}
#endif
		else {
			amrex::Abort("CD : How did we get here?");
		}
	}
}

template<typename DVM>
void DUGKSCore<DVM>::Update_xFluxX(Array<MultiFab, Ark::FACEDVDFNUM> &xDVDFnF, const Real hTime)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xDVDFnF[Ark::FLUX], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xTauh = xDVDFnF[Ark::TAUF].const_array(mfi);
        Array4<Real const> const &xEqh  = xDVDFnF[Ark::EQF].const_array(mfi);
        Array4<Real const> const &xSoh  = xDVDFnF[Ark::SOF].const_array(mfi);

        Array4<Real 	 > const &xFlux = xDVDFnF[Ark::FLUX].array(mfi);

        amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	const Real xtauh = xTauh(i,j,k);
        	const Real bh = hTime/(2.0*xtauh + hTime);
        	// const Real ch = bh*xtauh;
       	//
        	xFlux(i,j,k,n) = xFlux(i,j,k,n) 
        				   + bh*(xEqh(i,j,k,n) - xFlux(i,j,k,n) + xtauh*xSoh(i,j,k,n));
        	// xFlux(i,j,k,n) = xBh(i,j,k,n);
        });

  //       const auto lo = lbound(bx);
		// const auto hi = ubound(bx);
  //       LoopPS(lo,hi)
  //       {
  //       	const Real ah = 2.0*xtauh(i,j,k)/(2.0*xtauh(i,j,k) + hTime);
  //       	const Real bh = 1.0 - ah;
  //       	const Real ch = bh*xtauh(i,j,k);

  //       	LoopVS(DVM::Q)
  //       	{
  //       		xFlux(i,j,k,n) = ah*xBh(i,j,k,n) + bh*xEqh(i,j,k,n) + ch*xSoh(i,j,k,n);
  //       	}
  //       }
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Update_Wave (MultiFab &xWaveFab, const MultiFab &xBarPlusFab, 
							 	  const Vector<MultiFab const*> &xFluxFabnF,
							 	  const Real *dx, const Real hTime, const Real dTime)
{
	const Real dlx = dx[Ark::X];
	const Real dly = dx[Ark::Y];

#if(AMREX_SPACEDIM == 3)
	const Real dlz = dx[Ark::Z];
	const Real dSz = dlx*dly;
#else
	const Real dlz = 1.0;
#endif

	const Real dSx = dly*dlz;
	const Real dSy = dlz*dlx;

	const Real volume = dlx*dly*dlz;

	const MultiFab& xFluxFabxF = *xFluxFabnF[Ark::X];
	const MultiFab& xFluxFabyF = *xFluxFabnF[Ark::Y];

#if(AMREX_SPACEDIM == 3)
	const MultiFab& xFluxFabzF = *xFluxFabnF[Ark::Z];
#endif

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xWaveFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xFlux_xF = xFluxFabxF.const_array(mfi);
		Array4<Real const> const &xFlux_yF = xFluxFabyF.const_array(mfi);
#if(AMREX_SPACEDIM == 3)
		Array4<Real const> const &xFlux_zF = xFluxFabzF.const_array(mfi);
#endif
		Array4<Real const> const &xBP       = xBarPlusFab.const_array(mfi);
		Array4<Real 	 > const &xWave     = xWaveFab.array(mfi);

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	Real const xFluxSum
        	= (xFlux_xF(i,j,k,n) - xFlux_xF(i + 1,j,    k,    n))*dSx*DVM::xi_u[n]
			+ (xFlux_yF(i,j,k,n) - xFlux_yF(i,    j + 1,k,    n))*dSy*DVM::xi_v[n]
#if(AMREX_SPACEDIM == 3)
			+ (xFlux_zF(i,j,k,n) - xFlux_zF(i,    j    ,k + 1,n))*dSz*DVM::xi_w[n]
#endif			
			;

			const Real b = (hTime - dTime)/(hTime + dTime);
			Real const xWavePlus = xBP(i,j,k,n) + b*(xWave(i,j,k,n) - xBP(i,j,k,n));

			xWave(i,j,k,n) = xWavePlus + xFluxSum*dTime/volume;
        });
	}
}

//
#endif