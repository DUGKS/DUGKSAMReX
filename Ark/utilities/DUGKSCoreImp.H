#ifndef ARK_DUGKS_CORE_IMPLEMENTATION_H
#define ARK_DUGKS_CORE_IMPLEMENTATION_H

template<typename DVM>
void DUGKSCore<DVM>::Update_BarPlus(Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hTime, Real const dTime)
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

        amrex::ParallelFor(bx, DVM::Q, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
void DUGKSCore<DVM>::Update_BarPlus(Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hhTime)
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

        amrex::ParallelFor(bx, DVM::Q, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	const Real xtau = xTau(i,j,k);
        	const Real xneq = (xtau - hhTime)/xtau;

        	xBP(i,j,k,n) = xEq(i,j,k,n) + xneq*xWave(i,j,k,n) + hhTime*xSo(i,j,k,n);
        });
    }
}

// - Construct Prime from Wave 
template<typename DVM>
template<bool ENABLESOURCE>
void DUGKSCore<DVM>::Construct_Prime (Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hdTime)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xDVDF[Ark::WAVE], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xTau  = xDVDF[Ark::TAU].const_array(mfi);
        Array4<Real const> const &xEq   = xDVDF[Ark::EQ].const_array(mfi);
        Array4<Real const> const &xSo   = xDVDF[Ark::SO].const_array(mfi);

        Array4<Real>       const &xWave   = xDVDF[Ark::BP].array(mfi);

        amrex::ParallelFor(bx, DVM::Q, 
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	const Real xtau = xTau(i,j,k);
        	const Real b = hdTime/(xtau + hdTime);

        	// xWave(i,j,k,n) += b*(xEq(i,j,k,n) + xtau*xSo(i,j,k,n)*ENABLESOURCE - xWave(i,j,k,n));

        	// c++ 17
        	if constexpr (ENABLESOURCE) {
        		xWave(i,j,k,n) += b*(xEq(i,j,k,n) + xtau*xSo(i,j,k,n) - xWave(i,j,k,n));
        	}
        	else {
        		xWave(i,j,k,n) += b*(xEq(i,j,k,n)- xWave(i,j,k,n));
        	}
        });
    }
}

// - Construct NonEq from Wave or Prime
/*******************************************************************************\
f^{N} = \frac{2\tau}{2\tau+{\Delta}t}(\tilde{f}-f^{E}+\frac{{\Delta}t}{2}f^{S})
\*******************************************************************************/
template<typename DVM>
template<bool ENABLESOURCE>
void DUGKSCore<DVM>::Construct_NonEq (Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hdTime)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xDVDF[Ark::WAVE], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xTau  = xDVDF[Ark::TAU].const_array(mfi);
        Array4<Real const> const &xEq   = xDVDF[Ark::EQ].const_array(mfi);
        Array4<Real const> const &xSo   = xDVDF[Ark::SO].const_array(mfi);

        Array4<Real>       const &xWave    = xDVDF[Ark::WAVE].array(mfi);

        amrex::ParallelFor(bx, DVM::Q, 
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	const Real xtau = xTau(i,j,k);
        	const Real xneq = xtau/(xtau + hdTime);

        	// xWave(i,j,k,n) = xneq*(xWave(i,j,k,n) - xEq(i,j,k,n) + ENABLESOURCE*hdTime*xSo(i,j,k,n));

        	// c++ 17
        	if constexpr (ENABLESOURCE) {
        		xWave(i,j,k,n) = xneq*(xWave(i,j,k,n) - xEq(i,j,k,n) + hdTime*xSo(i,j,k,n));
        	}
        	else {
        		xWave(i,j,k,n) = xneq*(xWave(i,j,k,n) - xEq(i,j,k,n));
        	}
        });
    }
}

// - Construct Wave from NonEq
/*******************************************************************************\
\tilde{f} = f^{E} + \frac{2\tau+{\Delta}t}{2\tau}f^{N} - \frac{{\Delta}t}{2}f^{S}
\*******************************************************************************/
template<typename DVM>
template<bool ENABLESOURCE>
void DUGKSCore<DVM>::Construct_Wave (Array<MultiFab, Ark::DVDFNUM> &xDVDF, Real const hdTime)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xDVDF[Ark::WAVE], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xTau  = xDVDF[Ark::TAU].const_array(mfi);
        Array4<Real const> const &xEq   = xDVDF[Ark::EQ].const_array(mfi);
        Array4<Real const> const &xSo   = xDVDF[Ark::SO].const_array(mfi);

        Array4<Real>       const &xBP   = xDVDF[Ark::BP].array(mfi);
        Array4<Real>       const &xWave = xDVDF[Ark::WAVE].array(mfi);

        amrex::ParallelFor(bx, DVM::Q, 
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	const Real xtau = xTau(i,j,k);
        	const Real xneq = (xtau + hdTime)/xtau;

        	// xWave(i,j,k,n) = xEq(i,j,k,n) + xneq*xWave(i,j,k,n) - ENABLESOURCE*hdTime*xSo(i,j,k,n);

        	// c++ 17
        	if constexpr (ENABLESOURCE) {
        		xWave(i,j,k,n) = xEq(i,j,k,n) + xneq*xWave(i,j,k,n) - hdTime*xSo(i,j,k,n);
        	}
        	else {
        		xWave(i,j,k,n) = xEq(i,j,k,n) + xneq*xWave(i,j,k,n);
        	}
        });
    }
}

template<typename DVM>
void DUGKSCore<DVM>::Update_Grad_BarPlus (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real* dx)
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

template<typename DVM>
void DUGKSCore<DVM>::Update_Grad_BarPlus_VanLeer (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real* dx)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xDVDF[Ark::BP], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.growntilebox(1);

		Array4<Real> const &xBPx = xDVDF[Ark::BPx].array(mfi);
		Array4<Real> const &xBPy = xDVDF[Ark::BPy].array(mfi);

		Array4<Real  const> const &xBP = xDVDF[Ark::BP].const_array(mfi);

		GradFunc::Update_Grad_VanLeer_x(bx, DVM::Q, dx, xBP, xBPx);
		GradFunc::Update_Grad_VanLeer_y(bx, DVM::Q, dx, xBP, xBPy);

#if(AMREX_SPACEDIM == 3)
		Array4<Real> const &xBPz  = xDVDF[Ark::BPz].array(mfi);
		GradFunc::Update_Grad_VanLeer_z(bx, DVM::Q, dx, xBP, xBPz);
#endif
	}
}

template<typename DVM>
void DUGKSCore<DVM>::Update_Grad_BarPlus_2nd (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real* dx)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xDVDF[Ark::BP], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.growntilebox(1);

		Array4<Real> const &xBPx    = xDVDF[Ark::BPx].array(mfi);
		Array4<Real> const &xBPy    = xDVDF[Ark::BPy].array(mfi);

		Array4<Real  const> const &xBP = xDVDF[Ark::BP].const_array(mfi);

		GradFunc::Update_Grad_x(bx, DVM::Q, dx, xBP, xBPx);
		GradFunc::Update_Grad_y(bx, DVM::Q, dx, xBP, xBPy);

#if(AMREX_SPACEDIM == 3)
		Array4<Real> const &xBPz = xDVDF[Ark::BPz].array(mfi);
		GradFunc::Update_Grad_z(bx, DVM::Q, dx, xBP, xBPz);
#endif
	}
}

template<typename DVM>
void DUGKSCore<DVM>::Update_Grad_BarPlus_3rd (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real* dx)
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

		GradFunc::Update_Grad_x_5th(  bx, DVM::Q, dx, xBP, xBPx);
		GradFunc::Update_Grad_y_5th(  bx, DVM::Q, dx, xBP, xBPy);
		//
		GradFunc::Update_Grad_xx_5th( bx, DVM::Q, dx, xBP, xBPxx);
		GradFunc::Update_Grad_yy_5th( bx, DVM::Q, dx, xBP, xBPyy);
		GradFunc::Update_Grad_xy(     bx, DVM::Q, dx, xBP, xBPxy);
#if(AMREX_SPACEDIM == 3)
		GradFunc::Update_Grad_z_5th(  bx, DVM::Q, dx, xBP, xBPz);
		GradFunc::Update_Grad_zz_5th( bx, DVM::Q, dx, xBP, xBPzz);
		GradFunc::Update_Grad_yz(     bx, DVM::Q, dx, xBP, xBPyz);
		GradFunc::Update_Grad_zx(     bx, DVM::Q, dx, xBP, xBPzx);
#endif
	}
}

template<typename DVM>
void DUGKSCore<DVM>::Update_Grad_BarPlus_4th (Array<MultiFab, Ark::DVDFNUM> &xDVDF, 
											  const Real* dx)
{

#if(AMREX_SPACEDIM == 2)

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
		Array4<Real> const &xBPxxx  = xDVDF[Ark::BPxxx].array(mfi);
		Array4<Real> const &xBPxxy  = xDVDF[Ark::BPxxy].array(mfi);
		Array4<Real> const &xBPxyy  = xDVDF[Ark::BPxyy].array(mfi);
		Array4<Real> const &xBPyyy  = xDVDF[Ark::BPyyy].array(mfi);
		Array4<Real> const &xBPxxyy = xDVDF[Ark::BPxxyy].array(mfi);
		
		Array4<Real  const> const &xBP = xDVDF[Ark::BP].const_array(mfi);

		GradFunc::Update_Grad_x_5th(   bx, DVM::Q, dx, xBP, xBPx);
		GradFunc::Update_Grad_y_5th(   bx, DVM::Q, dx, xBP, xBPy);
		//
		GradFunc::Update_Grad_xx_5th(  bx, DVM::Q, dx, xBP, xBPxx);
		GradFunc::Update_Grad_yy_5th(  bx, DVM::Q, dx, xBP, xBPyy);
		GradFunc::Update_Grad_xy(  bx, DVM::Q, dx, xBP, xBPxy);
		GradFunc::Update_Grad_xxx( bx, DVM::Q, dx, xBP, xBPxxx);
		GradFunc::Update_Grad_xxy( bx, DVM::Q, dx, xBP, xBPxxy);
		GradFunc::Update_Grad_xyy( bx, DVM::Q, dx, xBP, xBPxyy);
		GradFunc::Update_Grad_yyy( bx, DVM::Q, dx, xBP, xBPyyy);
		GradFunc::Update_Grad_xxyy(bx, DVM::Q, dx, xBP, xBPxxyy);
	}

#endif
}

template<typename DVM>
void DUGKSCore<DVM>::Update_Grad_BarPlus_5th (Array<MultiFab, Ark::DVDFNUM> &xDVDF, 
											  const Real* dx)
{

#if(AMREX_SPACEDIM == 2)

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

		Array4<Real> const &xBPxxx  = xDVDF[Ark::BPxxx].array(mfi);
		Array4<Real> const &xBPyyy  = xDVDF[Ark::BPyyy].array(mfi);
		Array4<Real> const &xBPxxxx = xDVDF[Ark::BPxxxx].array(mfi);
		Array4<Real> const &xBPyyyy = xDVDF[Ark::BPyyyy].array(mfi);
		
		Array4<Real  const> const &xBP = xDVDF[Ark::BP].const_array(mfi);

		GradFunc::Update_Grad_x_5th(   bx, DVM::Q, dx, xBP, xBPx);
		GradFunc::Update_Grad_y_5th(   bx, DVM::Q, dx, xBP, xBPy);
		//
		GradFunc::Update_Grad_xx_5th(  bx, DVM::Q, dx, xBP, xBPxx);
		GradFunc::Update_Grad_yy_5th(  bx, DVM::Q, dx, xBP, xBPyy);
		GradFunc::Update_Grad_xy(  bx, DVM::Q, dx, xBP, xBPxy);

		GradFunc::Update_Grad_xxx( bx, DVM::Q, dx, xBP, xBPxxx);
		GradFunc::Update_Grad_yyy( bx, DVM::Q, dx, xBP, xBPyyy);
		GradFunc::Update_Grad_xxxx(bx, DVM::Q, dx, xBP, xBPxxxx);
		GradFunc::Update_Grad_yyyy(bx, DVM::Q, dx, xBP, xBPyyyy);
	}

#endif
}

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_UW_1st (MultiFab &xBarhFabnF, 
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

		if(FaceID == amrex::Direction::x) {

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_u[n]);
			int const di = token < 0 ? 0 : -1;
			//
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const ii = i + di;

				xBh(i,j,k,n) = xBP(ii,j,k,n);						 
			});
			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{					 
				xBh(i,j,k,n) += xBP(i,j,k,n);
				xBh(i,j,k,n) /= 2;
			});
			}		
		}
		}
		else if(FaceID == amrex::Direction::y){

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_v[n]);
			int const dj = token < 0 ? 0 : -1;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const jj = j + dj;
				xBh(i,j,k,n) = xBP(i,jj,k,n);
			});

			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{					 
				xBh(i,j,k,n) += xBP(i,j,k,n);
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
			int const dk = token < 0 ? 0 : -1;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const kk = k + dk;

				xBh(i,j,k,n) = xBP(i,j,kk,n); 
			});

			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) += xBP(i,j,k,n);
				xBh(i,j,k,n) /= 2.0;
			});
			}
		}
		}			
#endif
		else {
			amrex::Abort("UW-1st : How did we get here?");
		}
		
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_UW_2nd (MultiFab &xBarhFabnF, 
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

#if(AMREX_SPACEDIM == 3)
		Array4<Real const> const &xBPz  = xDVDF[Ark::BPz].const_array(mfi);
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
			const Real dly  = -hTime*DVM::xi_v[n];

#if(AMREX_SPACEDIM == 3)
			const Real dlz  = -hTime*DVM::xi_w[n];
#endif
			
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const ii = i + di;

#if(AMREX_SPACEDIM == 2)

				xBh(i,j,k,n) = xBP(ii,j,k,n) + xBPx(ii,j,k,n)*dlx + xBPy(ii,j,k,n)*dly;

#elif(AMREX_SPACEDIM == 3)

				xBh(i,j,k,n) = xBP(ii,j,k,n) 
							 + xBPx(ii,j,k,n)*dlx
							 + xBPy(ii,j,k,n)*dly
							 + xBPz(ii,j,k,n)*dlz;		

#endif							 

			});
			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
#if(AMREX_SPACEDIM == 2)

				Real tmp = xBP(i,j,k,n) - xBPx(i,j,k,n)*dlx + xBPy(i,j,k,n)*dly;

#elif(AMREX_SPACEDIM == 3)

				Real tmp = xBP(i,j,k,n) 
						 - xBPx(i,j,k,n)*dlx 
						 + xBPy(i,j,k,n)*dly 
						 + xBPz(i,j,k,n)*dlz;

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
			const Real dlx  = -hTime*DVM::xi_u[n];

#if(AMREX_SPACEDIM == 3)
			const Real dlz  = -hTime*DVM::xi_w[n];
#endif

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const jj = j + dj;

#if(AMREX_SPACEDIM == 2)

				xBh(i,j,k,n) = xBP(i,jj,k,n) + xBPx(i,jj,k,n)*dlx + xBPy(i,jj,k,n)*dly;

#elif(AMREX_SPACEDIM == 3)

				xBh(i,j,k,n) = xBP(i,jj,k,n) 
							 + xBPx(i,jj,k,n)*dlx 
							 + xBPy(i,jj,k,n)*dly 
							 + xBPz(i,jj,k,n)*dlz;			 
#endif
			});

			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
#if(AMREX_SPACEDIM == 2)	

				Real tmp = xBP(i,j,k,n) + xBPx(i,j,k,n)*dlx - xBPy(i,j,k,n)*dly;

#elif(AMREX_SPACEDIM == 3)

				Real tmp = xBP(i,j,k,n) 
						 + xBPx(i,j,k,n)*dlx 
						 - xBPy(i,j,k,n)*dly
						 + xBPz(i,j,k,n)*dlz;	
						 	 
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
			const Real dlx  = -hTime*DVM::xi_u[n];
			const Real dly  = -hTime*DVM::xi_v[n];

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const kk = k + dk;

				xBh(i,j,k,n) = xBP(i,j,kk,n) 
							 + xBPx(i,j,kk,n)*dlx 
							 + xBPy(i,j,kk,n)*dly 
							 + xBPz(i,j,kk,n)*dlz;
			});

			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				Real tmp = xBP(i,j,k,n) 
						 + xBPx(i,j,k,n)*dlx 
						 + xBPy(i,j,k,n)*dly
						 - xBPz(i,j,k,n)*dlz;		 

				xBh(i,j,k,n) += tmp;
				xBh(i,j,k,n) /= 2.0;
			});
			}
		}
		}			
#endif
		else {
			amrex::Abort("UW-2 : How did we get here?");
		}
		
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_UW_3rd (MultiFab &xBarhFabnF, 
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
#elif(AMREX_SPACEDIM == 3)

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

#elif(AMREX_SPACEDIM == 3)

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

#elif(AMREX_SPACEDIM == 3)

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

#elif(AMREX_SPACEDIM == 3)

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
			amrex::Abort("UW-3 : How did we get here?");
		}
		
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_UW_4th_cross (MultiFab &xBarhFabnF, 
									     		Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
									   	 		const Real* dx, const Real hTime, 
									 	 		amrex::Direction FaceID)
{

#if(AMREX_SPACEDIM == 2)

	AMREX_ASSERT(dx[Ark::X] == dx[Ark::Y]);
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
		Array4<Real const> const &xBPxxy  = xDVDF[Ark::BPxxy].const_array(mfi);
		Array4<Real const> const &xBPxyy  = xDVDF[Ark::BPxyy].const_array(mfi);
		Array4<Real const> const &xBPxxyy = xDVDF[Ark::BPxxyy].const_array(mfi);

		if(FaceID == amrex::Direction::x){

		const Real delta   = dx[Ark::X];
		const Real deltaxy = dx[Ark::X]*dx[Ark::Y]/24.;
		const Real deltaxxyy = deltaxy*deltaxy;

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_u[n]);
			int const multiplier = token < 0 ? -1 : 1;
			int const di = token < 0 ? 0 : -1;
			//
			const Real dlx  = 0.5*delta*multiplier - hTime*DVM::xi_u[n];
			const Real dlx2 = dlx*dlx;
			const Real dlx3 = dlx*dlx2;				
			//
			const Real dly  = -hTime*DVM::xi_v[n];
			const Real dly2 = dly*dly;
			const Real dly3 = dly*dly2;	

			const Real dlxy   = dlx*dly; 
			const Real dlxxy  = dlx2*dly;
			const Real dlxyy  = dlx*dly2;
			const Real dlxxyy = dlx2*dly2;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const ii = i + di;

				xBh(i,j,k,n) = xBP(ii,j,k,n) 
							 + xBPx(ii,j,k,n)*dlx + xBPy(ii,j,k,n)*dly
							 + xBPxy(ii,j,k,n)*dlxy
							 + xBPxx(ii,j,k,n)*(0.5*dlx2 - deltaxy)
							 + xBPyy(ii,j,k,n)*(0.5*dly2 - deltaxy)						 
							 + xBPxxy(ii,j,k,n)*(0.5*dlx2 - deltaxy)*dly
							 + xBPxyy(ii,j,k,n)*(0.5*dly2 - deltaxy)*dlx
							 + xBPxxyy(ii,j,k,n)*(0.25*dlxxyy + deltaxxyy 
							 					 - 0.5*(dlx2 + dly2)*deltaxy);

				// xBh(i,j,k,n) = xBP(ii,j,k,n) 
				// 			 + xBPx(ii,j,k,n)*dlx + xBPy(ii,j,k,n)*dly
				// 			 + xBPxy(ii,j,k,n)*dlxy
				// 			 + (xBPxx(ii,j,k,n)*dlx2 + xBPyy(ii,j,k,n)*dly2)*0.5
				// 			 + (xBPxxy(ii,j,k,n)*dlxxy + xBPxyy(ii,j,k,n)*dlxyy)*0.5
				// 			 + xBPxxyy(ii,j,k,n)*dlxxyy*0.25;
							 					 
			});
			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{

				Real tmp = xBP(i,j,k,n) 
						 - xBPx(i,j,k,n)*dlx + xBPy(i,j,k,n)*dly
						 - xBPxy(i,j,k,n)*dlxy
						 + xBPxx(i,j,k,n)*(0.5*dlx2 - deltaxy)
						 + xBPyy(i,j,k,n)*(0.5*dly2 - deltaxy)						 
						 + xBPxxy(i,j,k,n)*(0.5*dlx2 - deltaxy)*dly
						 - xBPxyy(i,j,k,n)*(0.5*dly2 - deltaxy)*dlx
						 + xBPxxyy(i,j,k,n)*(0.25*dlxxyy + deltaxxyy 
							 					 - 0.5*(dlx2 + dly2)*deltaxy);

				// Real tmp = xBP(i,j,k,n) 
				// 		 - xBPx(i,j,k,n)*dlx + xBPy(i,j,k,n)*dly
				// 		 - xBPxy(i,j,k,n)*dlxy
				// 		 + (xBPxx(i,j,k,n)*dlx2 + xBPyy(i,j,k,n)*dly2)*0.5					 
				// 		 + (xBPxxy(i,j,k,n)*dlxxy - xBPxyy(i,j,k,n)*dlxyy)*0.5
				// 		 + xBPxxyy(i,j,k,n)*dlxxyy*0.25;
					 
				xBh(i,j,k,n) += tmp;
				xBh(i,j,k,n) /= 2;
			});
			}		
		}
		}
		else if(FaceID == amrex::Direction::y){

		const Real delta = dx[Ark::Y];
		const Real deltaxy = dx[Ark::X]*dx[Ark::Y]/24.;
		const Real deltaxxyy = deltaxy*deltaxy;

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_v[n]);
			int const multiplier = token < 0 ? -1 : 1;
			int const dj = token < 0 ? 0 : -1;
			//
			const Real dly  = 0.5*delta*multiplier - hTime*DVM::xi_v[n];
			const Real dly2 = dly*dly;
			const Real dly3 = dly*dly2;
			//
			const Real dlx  = -hTime*DVM::xi_u[n];
			const Real dlx2 = dlx*dlx;	
			const Real dlx3 = dlx*dlx2;
			//
			const Real dlxy   = dlx*dly; 
			const Real dlxxy  = dlx2*dly;
			const Real dlxyy  = dlx*dly2;
			const Real dlxxyy = dlx2*dly2;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const jj = j + dj;

				xBh(i,j,k,n) = xBP(i,jj,k,n) 
							 + xBPx(i,jj,k,n)*dlx + xBPy(i,jj,k,n)*dly
							 + xBPxy(i,jj,k,n)*dlxy
							 + xBPxx(i,jj,k,n)*(0.5*dlx2 - deltaxy)
							 + xBPyy(i,jj,k,n)*(0.5*dly2 - deltaxy)						 
							 + xBPxxy(i,jj,k,n)*(0.5*dlx2 - deltaxy)*dly
							 + xBPxyy(i,jj,k,n)*(0.5*dly2 - deltaxy)*dlx
							 + xBPxxyy(i,jj,k,n)*(0.25*dlxxyy + deltaxxyy 
							 					 - 0.5*(dlx2 + dly2)*deltaxy);

				// xBh(i,j,k,n) = xBP(i,jj,k,n) 
				// 			 + xBPx(i,jj,k,n)*dlx + xBPy(i,jj,k,n)*dly
				// 			 + xBPxy(i,jj,k,n)*dlxy
				// 			 + (xBPxx(i,jj,k,n)*dlx2 + xBPyy(i,jj,k,n)*dly2)*0.5					 
				// 			 + (xBPxxy(i,jj,k,n)*dlxxy + xBPxyy(i,jj,k,n)*dlxyy)*0.5
				// 			 + xBPxxyy(i,jj,k,n)*dlxxyy*0.25;

			});

			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				
				Real tmp = xBP(i,j,k,n) 
						 + xBPx(i,j,k,n)*dlx - xBPy(i,j,k,n)*dly
						 - xBPxy(i,j,k,n)*dlxy
						 + xBPxx(i,j,k,n)*(0.5*dlx2 - deltaxy)
						 + xBPyy(i,j,k,n)*(0.5*dly2 - deltaxy)						 
						 - xBPxxy(i,j,k,n)*(0.5*dlx2 - deltaxy)*dly
						 + xBPxyy(i,j,k,n)*(0.5*dly2 - deltaxy)*dlx
						 + xBPxxyy(i,j,k,n)*(0.25*dlxxyy + deltaxxyy 
							 					 - 0.5*(dlx2 + dly2)*deltaxy);

				// Real tmp = xBP(i,j,k,n) 
				// 		 + xBPx(i,j,k,n)*dlx - xBPy(i,j,k,n)*dly - xBPxy(i,j,k,n)*dlxy						 
				// 		 + (xBPxx(i,j,k,n)*dlx2 + xBPyy(i,j,k,n)*dly2)*0.5					 
				// 		 + (xBPxyy(i,j,k,n)*dlxyy - xBPxxy(i,j,k,n)*dlxxy)*0.5
				// 		 + xBPxxyy(i,j,k,n)*dlxxyy*0.25;
						 
				xBh(i,j,k,n) += tmp;
				xBh(i,j,k,n) /= 2.0;
			});
			}
		}
		}
		else {
			amrex::Abort("UW-4th-cross : How did we get here?");
		}
		
	}
#endif

}

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_UW_4th (MultiFab &xBarhFabnF, 
									     Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
									   	 const Real* dx, const Real hTime, 
									 	 amrex::Direction FaceID)
{
#if(AMREX_SPACEDIM == 2)

	AMREX_ASSERT(dx[Ark::X] == dx[Ark::Y]);
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
		Array4<Real const> const &xBPxxx  = xDVDF[Ark::BPxxx].const_array(mfi);
		Array4<Real const> const &xBPxxy  = xDVDF[Ark::BPxxy].const_array(mfi);
		Array4<Real const> const &xBPxyy  = xDVDF[Ark::BPxyy].const_array(mfi);
		Array4<Real const> const &xBPyyy  = xDVDF[Ark::BPyyy].const_array(mfi);

		if(FaceID == amrex::Direction::x){

		const Real delta   = dx[Ark::X];
		const Real deltaxy = dx[Ark::X]*dx[Ark::Y]/24.;
		const Real deltaxxyy = deltaxy*deltaxy;

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_u[n]);
			int const multiplier = token < 0 ? -1 : 1;
			int const di = token < 0 ? 0 : -1;
			//
			const Real dlx  = 0.5*delta*multiplier - hTime*DVM::xi_u[n];
			const Real dlx2 = dlx*dlx;
			const Real dlx3 = dlx*dlx2;				
			//
			const Real dly  = -hTime*DVM::xi_v[n];
			const Real dly2 = dly*dly;
			const Real dly3 = dly*dly2;	

			const Real dlxy   = dlx*dly; 
			const Real dlxxy  = dlx2*dly;
			const Real dlxyy  = dlx*dly2;
			const Real dlxxyy = dlx2*dly2;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const ii = i + di;

				xBh(i,j,k,n) = xBP(ii,j,k,n) 
							 + xBPx(ii,j,k,n)*dlx + xBPy(ii,j,k,n)*dly + xBPxy(ii,j,k,n)*dlxy
							 + (xBPxx(ii,j,k,n)*dlx2 + xBPyy(ii,j,k,n)*dly2)*0.5							 
							 + (xBPxxx(ii,j,k,n)*dlx3 + xBPyyy(ii,j,k,n)*dly3)/6.0						 
							 + (xBPxxy(ii,j,k,n)*dlxxy + xBPxyy(ii,j,k,n)*dlxyy)*0.5
							 ;
			});
			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				Real tmp = xBP(i,j,k,n) 
						 - xBPx(i,j,k,n)*dlx + xBPy(i,j,k,n)*dly - xBPxy(i,j,k,n)*dlxy
						 + (xBPxx(i,j,k,n)*dlx2 + xBPyy(i,j,k,n)*dly2)*0.5						 
						 + (xBPyyy(i,j,k,n)*dly3  - xBPxxx(i,j,k,n)*dlx3)/6.0
						 + (xBPxxy(i,j,k,n)*dlxxy - xBPxyy(i,j,k,n)*dlxyy)*0.5						 
						 ;

				xBh(i,j,k,n) += tmp;
				xBh(i,j,k,n) /= 2;
			});
			}		
		}
		}
		else if(FaceID == amrex::Direction::y){

		const Real delta = dx[Ark::Y];
		const Real deltaxy = dx[Ark::X]*dx[Ark::Y]/24.;
		const Real deltaxxyy = deltaxy*deltaxy;

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_v[n]);
			int const multiplier = token < 0 ? -1 : 1;
			int const dj = token < 0 ? 0 : -1;
			//
			const Real dly  = 0.5*delta*multiplier - hTime*DVM::xi_v[n];
			const Real dly2 = dly*dly;
			const Real dly3 = dly*dly2;
			//
			const Real dlx  = -hTime*DVM::xi_u[n];
			const Real dlx2 = dlx*dlx;	
			const Real dlx3 = dlx*dlx2;
			//
			const Real dlxy   = dlx*dly; 
			const Real dlxxy  = dlx2*dly;
			const Real dlxyy  = dlx*dly2;
			const Real dlxxyy = dlx2*dly2;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const jj = j + dj;

				xBh(i,j,k,n) = xBP(i,jj,k,n) 
							 + xBPx(i,jj,k,n)*dlx + xBPy(i,jj,k,n)*dly + xBPxy(i,jj,k,n)*dlxy
							 + (xBPxx(i,jj,k,n)*dlx2 + xBPyy(i,jj,k,n)*dly2)*0.5
							 + (xBPxxx(i,jj,k,n)*dlx3 + xBPyyy(i,jj,k,n)*dly3)/6.0
							 + (xBPxxy(i,jj,k,n)*dlxxy + xBPxyy(i,jj,k,n)*dlxyy)*0.5							 
							;

			});

			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
					
				Real tmp = xBP(i,j,k,n) 
						 + xBPx(i,j,k,n)*dlx - xBPy(i,j,k,n)*dly - xBPxy(i,j,k,n)*dlxy
						 + (xBPxx(i,j,k,n)*dlx2 + xBPyy(i,j,k,n)*dly2)*0.5
						 + (xBPxxx(i,j,k,n)*dlx3 - xBPyyy(i,j,k,n)*dly3)/6.0
						 + (xBPxyy(i,j,k,n)*dlxyy - xBPxxy(i,j,k,n)*dlxxy)*0.5						 
						;
						 
				xBh(i,j,k,n) += tmp;
				xBh(i,j,k,n) /= 2.0;
			});
			}
		}
		}
		else {
			amrex::Abort("UW-4th-Full : How did we get here?");
		}
		
	}
#endif
}

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_UW_5th (MultiFab &xBarhFabnF, 
									 	 Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
									 	 const Real* dx, const Real hTime, 
									 	 amrex::Direction FaceID)
{
#if(AMREX_SPACEDIM == 2)

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
		Array4<Real const> const &xBPxxx  = xDVDF[Ark::BPxxx].const_array(mfi);
		Array4<Real const> const &xBPyyy  = xDVDF[Ark::BPyyy].const_array(mfi);
		Array4<Real const> const &xBPxxxx = xDVDF[Ark::BPxxxx].const_array(mfi);
		Array4<Real const> const &xBPyyyy = xDVDF[Ark::BPyyyy].const_array(mfi);

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
			const Real dlx3 = dlx2*dlx;
			const Real dlx4 = dlx3*dlx;	
			//
			const Real dly  = -hTime*DVM::xi_v[n];
			const Real dly2 = dly*dly;
			const Real dly3 = dly2*dly;
			const Real dly4 = dly3*dly;	

			const Real dlxy = dlx*dly; 
			
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const ii = i + di;

				xBh(i,j,k,n) = xBP(ii,j,k,n) 
							 + xBPx(ii,j,k,n)*dlx + xBPy(ii,j,k,n)*dly + xBPxy(ii,j,k,n)*dlxy
							 + (xBPxx(ii,j,k,n)*dlx2 + xBPyy(ii,j,k,n)*dly2)*0.5							 
							 + (xBPxxx(ii,j,k,n)*dlx3 + xBPyyy(ii,j,k,n)*dly3)/6.0
							 + (xBPxxxx(ii,j,k,n)*dlx4 + xBPyyyy(ii,j,k,n)*dly4)/24.0;

			});
			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				Real tmp = xBP(i,j,k,n) 
						 - xBPx(i,j,k,n)*dlx + xBPy(i,j,k,n)*dly - xBPxy(i,j,k,n)*dlxy
						 + (xBPxx(i,j,k,n)*dlx2 + xBPyy(i,j,k,n)*dly2)*0.5						 
						 + (xBPyyy(i,j,k,n)*dly3  - xBPxxx(i,j,k,n)*dlx3)/6.0
						 + (xBPxxxx(i,j,k,n)*dlx4 + xBPyyyy(i,j,k,n)*dly4)/24.0;
						 
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
			const Real dly3 = dly*dly2;
			const Real dly4 = dly*dly3;
			//
			const Real dlx  = -hTime*DVM::xi_u[n];
			const Real dlx2 = dlx*dlx;
			const Real dlx3 = dlx*dlx2;
			const Real dlx4 = dlx*dlx3;	
			//
			const Real dlxy = dlx*dly;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const jj = j + dj;

				xBh(i,j,k,n) = xBP(i,jj,k,n) 
							 + xBPx(i,jj,k,n)*dlx + xBPy(i,jj,k,n)*dly + xBPxy(i,jj,k,n)*dlxy
							 + (xBPxx(i,jj,k,n)*dlx2 + xBPyy(i,jj,k,n)*dly2)*0.5
							 + (xBPxxx(i,jj,k,n)*dlx3 + xBPyyy(i,jj,k,n)*dly3)/6.0
							 + (xBPxxxx(i,jj,k,n)*dlx4 + xBPyyyy(i,jj,k,n)*dly4)/24.0;

			});

			if(token == 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{

				Real tmp = xBP(i,j,k,n) 
						 + xBPx(i,j,k,n)*dlx - xBPy(i,j,k,n)*dly - xBPxy(i,j,k,n)*dlxy
						 + (xBPxx(i,j,k,n)*dlx2 + xBPyy(i,j,k,n)*dly2)*0.5
						 + (xBPxxx(i,j,k,n)*dlx3 - xBPyyy(i,j,k,n)*dly3)/6.0
						 + (xBPxxxx(i,j,k,n)*dlx4 + xBPyyyy(i,j,k,n)*dly4)/24.0;

// #if(AMREX_SPACEDIM == 3)
// 				Real tmp = xBP(i,j,k,n) 
// 						 + xBPx(i,j,k,n)*dlx 
// 						 - xBPy(i,j,k,n)*dly
// 						 + xBPz(i,j,k,n)*dlz
// 						 - xBPxy(i,j,k,n)*dlxy
// 						 - xBPyz(i,j,k,n)*dlyz
// 						 + xBPzx(i,j,k,n)*dlzx
// 						 + (xBPxx(i,j,k,n)*dlx2 
// 						   +xBPyy(i,j,k,n)*dly2
// 						   +xBPzz(i,j,k,n)*dlz2)*0.5
// 						 + (xBPxxx(i,j,k,n)*dlx3
// 						   -xBPyyy(i,j,k,n)*dly3
// 						   +xBPzzz(i,j,k,n)*dlz3)/6.0
// 						 + (xBPxxxx(i,j,k,n)*dlx4
// 						   +xBPyyyy(i,j,k,n)*dly4
// 						   +xBPzzzz(i,j,k,n)*dlz4)/24.0;		 
// #endif							 
				xBh(i,j,k,n) += tmp;
				xBh(i,j,k,n) /= 2.0;
			});
			}
		}
		}
// #if(AMREX_SPACEDIM == 3)
// 		else if(FaceID == amrex::Direction::z){

// 		const Real delta = dx[Ark::Z];

// 		LoopVS(DVM::Q)
// 		{
// 			int const token = Ark::sign(DVM::xi_w[n]);
// 			int const multiplier = token < 0 ? -1 : 1;
// 			int const dk = token < 0 ? 0 : -1;
// 			//
// 			const Real dlz  = 0.5*delta*multiplier - hTime*DVM::xi_w[n];
// 			const Real dlz2 = dlz*dlz;
// 			const Real dlz3 = dlz*dlz2;
// 			const Real dlz4 = dlz*dlz3;
// 			//
// 			const Real dlx  = -hTime*DVM::xi_u[n];
// 			const Real dlx2 = dlx*dlx;
// 			const Real dlx3 = dlx*dlx2;
// 			const Real dlx4 = dlx*dlx3;
// 			//
// 			const Real dly  = -hTime*DVM::xi_v[n];
// 			const Real dly2 = dly*dly;
// 			const Real dly3 = dly*dly2;
// 			const Real dly4 = dly*dly3;

// 			const Real dlxy = dlx*dly;
// 			const Real dlyz = dly*dlz;
// 			const Real dlzx = dlz*dlx;

// 			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
// 			{
// 				int const kk = k + dk;

// 				xBh(i,j,k,n) = xBP(     i,j,kk,n) 
// 							 + xBPx(    i,j,kk,n)*dlx 
// 							 + xBPy(    i,j,kk,n)*dly 
// 							 + xBPz(    i,j,kk,n)*dlz
// 							 + xBPxy(   i,j,kk,n)*dlxy
// 							 + xBPyz(   i,j,kk,n)*dlyz
// 							 + xBPzx(   i,j,kk,n)*dlzx
// 							 + (xBPxx(  i,j,kk,n)*dlx2 
// 							   +xBPyy(  i,j,kk,n)*dly2
// 							   +xBPzz(  i,j,kk,n)*dlz2)*0.5
// 							 + (xBPxxx( i,j,kk,n)*dlx3
// 							   +xBPyyy( i,j,kk,n)*dly3
// 							   +xBPzzz( i,j,kk,n)*dlz3)/6.0
// 							 + (xBPxxxx(i,j,kk,n)*dlx4
// 							   +xBPyyyy(i,j,kk,n)*dly4
// 							   +xBPzzzz(i,j,kk,n)*dlz4)/24.0;
// 			});

// 			if(token == 0) {
// 			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
// 			{
// 				Real tmp = xBP(i,j,k,n) 
// 						 + xBPx(i,j,k,n)*dlx 
// 						 + xBPy(i,j,k,n)*dly
// 						 - xBPz(i,j,k,n)*dlz
// 						 + xBPxy(i,j,k,n)*dlxy
// 						 - xBPyz(i,j,k,n)*dlyz
// 						 - xBPzx(i,j,k,n)*dlzx
// 						 + (xBPxx(i,j,k,n)*dlx2 
// 						   +xBPyy(i,j,k,n)*dly2
// 						   +xBPzz(i,j,k,n)*dlz2)*0.5
// 						 + (xBPxxx(i,j,k,n)*dlx3
// 						   +xBPyyy(i,j,k,n)*dly3
// 						   -xBPzzz(i,j,k,n)*dlz3)/6.0
// 						 + (xBPxxxx(i,j,k,n)*dlx4
// 						   +xBPyyyy(i,j,k,n)*dly4
// 						   +xBPzzzz(i,j,k,n)*dlz4)/24.0;

// 				xBh(i,j,k,n) += tmp;
// 				xBh(i,j,k,n) /= 2.0;
// 			});
// 			}
// 		}
// 		}			
// #endif
		else {
			amrex::Abort("UW-5th : How did we get here?");
		}
		
	}
#endif
}

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_UW_Limiter (MultiFab &xBarhFabnF, 
											Array<MultiFab, Ark::DVDFNUM> const &xDVDF,
						 					amrex::Direction FaceID)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xBarhFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real> const &xBh = xBarhFabnF.array(mfi);

		Array4<Real const> const &xBP = xDVDF[Ark::BP].const_array(mfi);

		if(FaceID == amrex::Direction::x) {

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_u[n]);
			int const di = token < 0 ? 0 : -1;
			//
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const ii = i + di;

				const double vL =  xBP(ii,  j,k,n) - xBP(ii-1,j,k,n);
				const double vR =  xBP(ii+1,j,k,n) - xBP(ii,  j,k,n);

				// const double vanLeer = vanC > 0 ? 2*vanL/(vanL + vanR) : 0;
				// const double MinMod = vanC > 0 ? std::min(vL/vR,1) : 0;
				const double limiter = Limiter::VanAlbada2(vL,vR);

				if(token == 0) {

					const double vL0 =  xBP(i,  j,k,n) - xBP(i-1,j,k,n);
					const double vR0 =  xBP(i+1,j,k,n) - xBP(i,  j,k,n);
					const double limiter0 = std::min(Limiter::VanAlbada2(vL0,vR0), limiter);
					
					const double xtmp = (xBP(ii,j,k,n) + xBP(i,j,k,n))*0.5;

					xBh(i,j,k,n) = xtmp + limiter0*(xBh(i,j,k,n) - xtmp);								 
				}
				else {
					xBh(i,j,k,n) = xBP(ii,j,k,n) + limiter*(xBh(i,j,k,n) - xBP(ii,j,k,n));
				}					 

			});	
		}
		}
		else if(FaceID == amrex::Direction::y) {

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_v[n]);
			int const dj = token < 0 ? 0 : -1;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const jj = j + dj;

				const double vF =  xBP(i,jj,  k,n) - xBP(i,jj-1,k,n);
				const double vB =  xBP(i,jj+1,k,n) - xBP(i,jj,  k,n);

				const double limiter = Limiter::VanAlbada2(vF,vB);

				if(token == 0) {

					const double vF0 =  xBP(i,j,  k,n) - xBP(i,j-1,k,n);
					const double vB0 =  xBP(i,j+1,k,n) - xBP(i,j,  k,n);
					const double limiter0 = std::min(Limiter::VanAlbada2(vF0,vB0), limiter);

					const double xtmp = (xBP(i,jj,k,n) + xBP(i,j,k,n))*0.5;

					xBh(i,j,k,n) = xtmp + limiter0*(xBh(i,j,k,n) - xtmp);
				}
				else {
					xBh(i,j,k,n) = xBP(i,jj,k,n) + limiter*(xBh(i,j,k,n) - xBP(i,jj,k,n));
				}
			});
		}
		}
#if(AMREX_SPACEDIM == 3)
		else if(FaceID == amrex::Direction::z) {

		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_w[n]);
			int const dk = token < 0 ? 0 : -1;

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const kk = k + dk;

				const double vF =  xBP(i,j,kk  ,n) - xBP(i,j,kk-1,n);
				const double vB =  xBP(i,j,kk+1,n) - xBP(i,j,kk  ,n);

				const double limiter = Limiter::VanAlbada2(vF,vB);

				if(token == 0) {

					const double vF0 =  xBP(i,j,k  ,n) - xBP(i,j,k-1,n);
					const double vB0 =  xBP(i,j,k+1,n) - xBP(i,j,k  ,n);
					const double limiter0 = std::min(Limiter::VanAlbada2(vF0,vB0), limiter);

					const double xtmp = (xBP(i,j,kk,n) + xBP(i,j,k,n))*0.5;

					xBh(i,j,k,n) = xtmp + limiter0*(xBh(i,j,k,n) - xtmp);
				}
				else {
					xBh(i,j,k,n) = xBP(i,j,kk,n) + limiter*(xBh(i,j,k,n) - xBP(i,j,kk,n));
				}
			});
		}
		}			
#endif
		else {
			amrex::Abort("Limiter : How did we get here?");
		}
		
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_CD (MultiFab &xBarhFabnF, const MultiFab &xBarPlusFab,
									 const Real* dx, const Real hTime, 
									 amrex::Direction FaceID)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xBarhFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Real dlx = dx[Ark::X];
		const Real dly = dx[Ark::Y];
#if(AMREX_SPACEDIM == 3)
		const Real dlz = dx[Ark::Z];
#endif
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xBP = xBarPlusFab.const_array(mfi);
		Array4<Real      > const &xBh = xBarhFabnF.array(mfi);

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

#elif(AMREX_SPACEDIM == 3)

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

#elif(AMREX_SPACEDIM == 3)

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

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_WENO_3rd (MultiFab &xBarhFabnF, const MultiFab &xBarPlusFab,
									 	   const Real* dx, const Real hTime, 
									 	   amrex::Direction FaceID)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xBarhFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Real dlx = dx[Ark::X];
		const Real dly = dx[Ark::Y];
#if(AMREX_SPACEDIM == 3)
		const Real dlz = dx[Ark::Z];
#endif
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xBP = xBarPlusFab.const_array(mfi);
		Array4<Real      > const &xBh = xBarhFabnF.array(mfi);

		if(FaceID == amrex::Direction::x) {

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
			int const token = Ark::sign(DVM::xi_u[n]);

			if(token > 0) {				
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i-3,j,k,n), xBP(i-2,j,k,n), xBP(i-1,j,k,n),
												xBP(i  ,j,k,n), xBP(i+1,j,k,n));
			}
			else if(token < 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i+2,j,k,n), xBP(i+1,j,k,n), xBP(i,j,k,n),
												xBP(i-1,j,k,n), xBP(i-2,j,k,n));
			}
			else {
				xBh(i,j,k,n) = GradFunc::WENO_C(xBP(i-3,j,k,n), xBP(i-2,j,k,n), xBP(i-1,j,k,n),
								  				xBP(i  ,j,k,n), xBP(i+1,j,k,n), xBP(i+2,j,k,n));
			}

			Real xBP_x = (xBP(i,j,k,n) - xBP(i-1,j,k,n))/dlx;

			Real xBP_y = (xBP(i,  j+1,k,n) - xBP(i,  j-1,k,n) 
					     +xBP(i-1,j+1,k,n) - xBP(i-1,j-1,k,n)
					     )*0.25/dly;

#if(AMREX_SPACEDIM == 2)

			xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y);

#elif(AMREX_SPACEDIM == 3)

			Real xBP_z = (xBP(i,  j,k+1,n) - xBP(i,  j,k-1,n)
        				 +xBP(i-1,j,k+1,n) - xBP(i-1,j,k-1,n)
        				 )*0.25/dlz;

			xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y + DVM::xi_w[n]*xBP_z);
#endif
		});
		}

		else if(FaceID == amrex::Direction::y){

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	int const token = Ark::sign(DVM::xi_v[n]);
			
			if(token > 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z (xBP(i,j-3,k,n), xBP(i,j-2,k,n), xBP(i,j-1,k,n),
												 xBP(i,j  ,k,n), xBP(i,j+1,k,n));				
			}
			else if (token < 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z (xBP(i,j+2,k,n), xBP(i,j+1,k,n), xBP(i,j,k,n),
												 xBP(i,j-1,k,n), xBP(i,j-2,k,n));
			}
			else {
				xBh(i,j,k,n) = GradFunc::WENO_C (xBP(i,j-3,k,n), xBP(i,j-2,k,n), xBP(i,j-1,k,n),
								  				 xBP(i,j  ,k,n), xBP(i,j+1,k,n), xBP(i,j+2,k,n));
			}

        	Real xBP_x = (xBP(i+1,j,  k,n) - xBP(i-1,j,  k,n)
        				 +xBP(i+1,j-1,k,n) - xBP(i-1,j-1,k,n)
        				 )*0.25/dlx;

        	Real xBP_y = (xBP(i,j,k,n) - xBP(i,j-1,k,n))/dly;

#if(AMREX_SPACEDIM == 2)

        	xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y);

#elif(AMREX_SPACEDIM == 3)

        	Real xBP_z = (xBP(i,j,  k+1,n) - xBP(i,j,  k-1,n)
        				 +xBP(i,j-1,k+1,n) - xBP(i,j-1,k-1,n)
        				 )*0.25/dlz;

        	xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y + DVM::xi_w[n]*xBP_z);
#endif        			    
        });
		}
#if(AMREX_SPACEDIM == 3)
		else if(FaceID == amrex::Direction::z){

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
			int const token = Ark::sign(DVM::xi_w[n]);
			
			if(token > 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i,j,k-3,n), xBP(i,j,k-2,n), xBP(i,j,k-1,n),
												xBP(i,j,k  ,n), xBP(i,j,k+1,n));							 
			}
			else if (token < 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i,j,k+2,n), xBP(i,j,k+1,n), xBP(i,j,k,n),
												xBP(i,j,k-1,n), xBP(i,j,k-2,n));
			}
			else {
				xBh(i,j,k,n) = GradFunc::WENO_C (xBP(i,j,k-3,n), xBP(i,j,k-2,n), xBP(i,j,k-1,n),
								  				 xBP(i,j,k  ,n), xBP(i,j,k+1,n), xBP(i,j,k+2,n));
			}

			Real xBP_x = (xBP(i+1,j,k,  n) - xBP(i-1,j,k,  n) 
					     +xBP(i+1,j,k-1,n) - xBP(i-1,j,k-1,n)
					     )*0.25/dlx;

			Real xBP_y = (xBP(i,j+1,k,  n) - xBP(i,j-1,k,  n) 
					     +xBP(i,j+1,k-1,n) - xBP(i,j-1,k-1,n)
					     )*0.25/dly;

			Real xBP_z = (xBP(i,j,k,n) - xBP(i,j,k-1,n))/dlz;

			xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y + DVM::xi_w[n]*xBP_z);
		});
		}
#endif
		else {
			amrex::Abort("WENO-3 : How did we get here?");
		}
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Update_Barh_WENO_5th (MultiFab &xBarhFabnF, const MultiFab &xBarPlusFab,
									 	   const Real* dx, const Real hTime, 
									 	   amrex::Direction FaceID)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xBarhFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const Real dlx = dx[Ark::X];
		const Real dly = dx[Ark::Y];

#if(AMREX_SPACEDIM == 3)
		const Real dlz = dx[Ark::Z];
#endif
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real const> const &xBP = xBarPlusFab.const_array(mfi);
		Array4<Real      > const &xBh = xBarhFabnF.array(mfi);

		if(FaceID == amrex::Direction::x) {

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
			int const token = Ark::sign(DVM::xi_u[n]);

			if(token > 0) {				
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i-3,j,k,n), xBP(i-2,j,k,n), xBP(i-1,j,k,n),
												xBP(i  ,j,k,n), xBP(i+1,j,k,n));
			}
			else if(token < 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i+2,j,k,n), xBP(i+1,j,k,n), xBP(i,j,k,n),
												xBP(i-1,j,k,n), xBP(i-2,j,k,n));
			}
			else {
				xBh(i,j,k,n) = GradFunc::WENO_C(xBP(i-3,j,k,n), xBP(i-2,j,k,n), xBP(i-1,j,k,n),
								  				xBP(i  ,j,k,n), xBP(i+1,j,k,n), xBP(i+2,j,k,n));
			}

			Real xBP_x = (xBP(i,j,k,n) - xBP(i-1,j,k,n))/dlx;

			Real xBP_y = (xBP(i,  j+1,k,n) - xBP(i,  j-1,k,n) 
					     +xBP(i-1,j+1,k,n) - xBP(i-1,j-1,k,n)
					     )*0.25/dly;

#if(AMREX_SPACEDIM == 2)

			xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y);

#elif(AMREX_SPACEDIM == 3)

			Real xBP_z = (xBP(i,  j,k+1,n) - xBP(i,  j,k-1,n)
        				 +xBP(i-1,j,k+1,n) - xBP(i-1,j,k-1,n)
        				 )*0.25/dlz;

			xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y + DVM::xi_w[n]*xBP_z);
#endif
		});
		}

		else if(FaceID == amrex::Direction::y){

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
        	int const token = Ark::sign(DVM::xi_v[n]);
			
			if(token > 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z (xBP(i,j-3,k,n), xBP(i,j-2,k,n), xBP(i,j-1,k,n),
												 xBP(i,j  ,k,n), xBP(i,j+1,k,n));				
			}
			else if (token < 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z (xBP(i,j+2,k,n), xBP(i,j+1,k,n), xBP(i,j,k,n),
												 xBP(i,j-1,k,n), xBP(i,j-2,k,n));
			}
			else {
				xBh(i,j,k,n) = GradFunc::WENO_C (xBP(i,j-3,k,n), xBP(i,j-2,k,n), xBP(i,j-1,k,n),
								  				 xBP(i,j  ,k,n), xBP(i,j+1,k,n), xBP(i,j+2,k,n));
			}

        	Real xBP_x = (xBP(i+1,j,  k,n) - xBP(i-1,j,  k,n)
        				 +xBP(i+1,j-1,k,n) - xBP(i-1,j-1,k,n)
        				 )*0.25/dlx;

        	Real xBP_y = (xBP(i,j,k,n) - xBP(i,j-1,k,n))/dly;

#if(AMREX_SPACEDIM == 2)

        	xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y);

#elif(AMREX_SPACEDIM == 3)

        	Real xBP_z = (xBP(i,j,  k+1,n) - xBP(i,j,  k-1,n)
        				 +xBP(i,j-1,k+1,n) - xBP(i,j-1,k-1,n)
        				 )*0.25/dlz;

        	xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y + DVM::xi_w[n]*xBP_z);
#endif        			    
        });
		}
#if(AMREX_SPACEDIM == 3)
		else if(FaceID == amrex::Direction::z){

		amrex::ParallelFor(bx, DVM::Q,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
			int const token = Ark::sign(DVM::xi_w[n]);
			
			if(token > 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i,j,k-3,n), xBP(i,j,k-2,n), xBP(i,j,k-1,n),
												xBP(i,j,k  ,n), xBP(i,j,k+1,n));							 
			}
			else if (token < 0) {
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i,j,k+2,n), xBP(i,j,k+1,n), xBP(i,j,k,n),
												xBP(i,j,k-1,n), xBP(i,j,k-2,n));
			}
			else {
				xBh(i,j,k,n) = GradFunc::WENO_C (xBP(i,j,k-3,n), xBP(i,j,k-2,n), xBP(i,j,k-1,n),
								  				 xBP(i,j,k  ,n), xBP(i,j,k+1,n), xBP(i,j,k+2,n));
			}

			Real xBP_x = (xBP(i+1,j,k,  n) - xBP(i-1,j,k,  n) 
					     +xBP(i+1,j,k-1,n) - xBP(i-1,j,k-1,n)
					     )*0.25/dlx;

			Real xBP_y = (xBP(i,j+1,k,  n) - xBP(i,j-1,k,  n) 
					     +xBP(i,j+1,k-1,n) - xBP(i,j-1,k-1,n)
					     )*0.25/dly;

			Real xBP_z = (xBP(i,j,k,n) - xBP(i,j,k-1,n))/dlz;

			xBh(i,j,k,n) -= hTime*(DVM::xi_u[n]*xBP_x + DVM::xi_v[n]*xBP_y + DVM::xi_w[n]*xBP_z);
		});
		}
#endif
		else {
			amrex::Abort("WENO-5 : How did we get here?");
		}
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Update_FD_Barh_WENO_3rd (MultiFab &xBarhFabnF, const MultiFab &xBarPlusFab,
								  	   		  amrex::Direction FaceID)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xBarhFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real> const &xBh = xBarhFabnF.array(mfi);
		Array4<Real const> const &xBP = xBarPlusFab.const_array(mfi);

		if(FaceID == amrex::Direction::x){
		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_u[n]);
			
			if(token > 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i-3,j,k,n), xBP(i-2,j,k,n), xBP(i-1,j,k,n),
												xBP(i  ,j,k,n), xBP(i+1,j,k,n));							 
			});
			}
			else if (token < 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i+2,j,k,n), xBP(i+1,j,k,n), xBP(i,j,k,n),
												xBP(i-1,j,k,n), xBP(i-2,j,k,n));
			});
			}
			else {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = 
				GradFunc::WENO_C (xBP(i-3,j,k,n), xBP(i-2,j,k,n), xBP(i-1,j,k,n),
								  xBP(i  ,j,k,n), xBP(i+1,j,k,n), xBP(i+2,j,k,n));
			});
			}		
		}
		}
		else if(FaceID == amrex::Direction::y){
		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_v[n]);
			
			if(token > 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z (xBP(i,j-3,k,n), xBP(i,j-2,k,n), xBP(i,j-1,k,n),
												xBP(i,j  ,k,n), xBP(i,j+1,k,n));							 
			});
			}
			else if (token < 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z (xBP(i,j+2,k,n), xBP(i,j+1,k,n), xBP(i,j,k,n),
												 xBP(i,j-1,k,n), xBP(i,j-2,k,n));
			});
			}
			else {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = 
				GradFunc::WENO_C (xBP(i,j-3,k,n), xBP(i,j-2,k,n), xBP(i,j-1,k,n),
								  xBP(i,j  ,k,n), xBP(i,j+1,k,n), xBP(i,j+2,k,n));
			});
			}		
		}
		}
#if(AMREX_SPACEDIM == 3)
		else if(FaceID == amrex::Direction::z){
		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_w[n]);
			
			if(token > 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i,j,k-3,n), xBP(i,j,k-2,n), xBP(i,j,k-1,n),
												xBP(i,j,k  ,n), xBP(i,j,k+1,n));							 
			});
			}
			else if (token < 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i,j,k+2,n), xBP(i,j,k+1,n), xBP(i,j,k,n),
												xBP(i,j,k-1,n), xBP(i,j,k-2,n));
			});
			}
			else {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = 
				GradFunc::WENO_C (xBP(i,j,k-3,n), xBP(i,j,k-2,n), xBP(i,j,k-1,n),
								  xBP(i,j,k  ,n), xBP(i,j,k+1,n), xBP(i,j,k+2,n));
			});
			}		
		}
		}
#endif
		else {
			amrex::Abort("FD-WENO-3 : How did we get here?");
		}
	}
}

template <typename DVM>
void DUGKSCore<DVM>::Update_FD_Barh_WENO_5th (MultiFab &xBarhFabnF, const MultiFab &xBarPlusFab,
								  	   		  amrex::Direction FaceID)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xBarhFabnF, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real> const &xBh = xBarhFabnF.array(mfi);
		Array4<Real const> const &xBP = xBarPlusFab.const_array(mfi);

		if(FaceID == amrex::Direction::x){
		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_u[n]);
			
			if(token > 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i-3,j,k,n), xBP(i-2,j,k,n), xBP(i-1,j,k,n),
												xBP(i  ,j,k,n), xBP(i+1,j,k,n));							 
			});
			}
			else if (token < 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i+2,j,k,n), xBP(i+1,j,k,n), xBP(i,j,k,n),
												xBP(i-1,j,k,n), xBP(i-2,j,k,n));
			});
			}
			else {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = 
				GradFunc::WENO_C (xBP(i-3,j,k,n), xBP(i-2,j,k,n), xBP(i-1,j,k,n),
								  xBP(i  ,j,k,n), xBP(i+1,j,k,n), xBP(i+2,j,k,n));
			});
			}		
		}
		}
		else if(FaceID == amrex::Direction::y){
		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_v[n]);
			
			if(token > 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z (xBP(i,j-3,k,n), xBP(i,j-2,k,n), xBP(i,j-1,k,n),
												xBP(i,j  ,k,n), xBP(i,j+1,k,n));							 
			});
			}
			else if (token < 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z (xBP(i,j+2,k,n), xBP(i,j+1,k,n), xBP(i,j,k,n),
												 xBP(i,j-1,k,n), xBP(i,j-2,k,n));
			});
			}
			else {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = 
				GradFunc::WENO_C (xBP(i,j-3,k,n), xBP(i,j-2,k,n), xBP(i,j-1,k,n),
								  xBP(i,j  ,k,n), xBP(i,j+1,k,n), xBP(i,j+2,k,n));
			});
			}		
		}
		}
#if(AMREX_SPACEDIM == 3)
		else if(FaceID == amrex::Direction::z){
		LoopVS(DVM::Q)
		{
			int const token = Ark::sign(DVM::xi_w[n]);
			
			if(token > 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i,j,k-3,n), xBP(i,j,k-2,n), xBP(i,j,k-1,n),
												xBP(i,j,k  ,n), xBP(i,j,k+1,n));							 
			});
			}
			else if (token < 0) {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = GradFunc::WENO_Z(xBP(i,j,k+2,n), xBP(i,j,k+1,n), xBP(i,j,k,n),
												xBP(i,j,k-1,n), xBP(i,j,k-2,n));
			});
			}
			else {
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				xBh(i,j,k,n) = 
				GradFunc::WENO_C (xBP(i,j,k-3,n), xBP(i,j,k-2,n), xBP(i,j,k-1,n),
								  xBP(i,j,k  ,n), xBP(i,j,k+1,n), xBP(i,j,k+2,n));
			});
			}		
		}
		}
#endif
		else {
			amrex::Abort("FD_WENO_5 : How did we get here?");
		}
	}
}	

template<typename DVM>
template<bool ENABLESOURCE>
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

        	xFlux(i,j,k,n) = xFlux(i,j,k,n) 
        				   + bh*(xEqh(i,j,k,n) - xFlux(i,j,k,n) + xtauh*xSoh(i,j,k,n)*ENABLESOURCE);
        	// xFlux(i,j,k,n) = xBh(i,j,k,n);
        });
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

// template <typename DVM>
// void DUGKSCore<DVM>::Update_Wave (Array<MultiFab, Ark::DVDFNUM> &xDVDF, 
// 							 	  const Vector<MultiFab const*> &xFluxFabnF,
// 							 	  const Real *dx, const Real dTime)
// {
// 	const Real hdTime = dTime/2.0;

// 	const Real dlx = dx[Ark::X];
// 	const Real dly = dx[Ark::Y];

// #if(AMREX_SPACEDIM == 3)
// 	const Real dlz = dx[Ark::Z];
// 	const Real dSz = dlx*dly;
// #else
// 	const Real dlz = 1.0;
// #endif

// 	const Real dSx = dly*dlz;
// 	const Real dSy = dlz*dlx;

// 	const Real volume = dlx*dly*dlz;

// 	MultiFab& xWaveFab = xDVDF[Ark::WAVE];

// 	const MultiFab& TauFab = xDVDF[Ark::TAU];
// 	const MultiFab& xEqFab = xDVDF[Ark::EQ];
// 	const MultiFab& xSoFab = xDVDF[Ark::SO];

// 	const MultiFab& xFluxFabxF = *xFluxFabnF[Ark::X];
// 	const MultiFab& xFluxFabyF = *xFluxFabnF[Ark::Y];

// #if(AMREX_SPACEDIM == 3)
// 	const MultiFab& xFluxFabzF = *xFluxFabnF[Ark::Z];
// #endif

// 	#ifdef _OPENMP
// 	#pragma omp parallel
// 	#endif
// 	for(amrex::MFIter mfi(xWaveFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
// 	{
// 		const amrex::Box& bx = mfi.tilebox();

// 		Array4<Real const> const &xFlux_xF = xFluxFabxF.const_array(mfi);
// 		Array4<Real const> const &xFlux_yF = xFluxFabyF.const_array(mfi);
// #if(AMREX_SPACEDIM == 3)
// 		Array4<Real const> const &xFlux_zF = xFluxFabzF.const_array(mfi);
// #endif
// 		Array4<Real const> const &xEq       = xEqFab.const_array(mfi);
// 		Array4<Real const> const &xSo       = xSoFab.const_array(mfi);
// 		Array4<Real const> const &xtau      = TauFab.const_array(mfi);
// 		Array4<Real 	 > const &xWave     = xWaveFab.array(mfi);

// 		amrex::ParallelFor(bx, DVM::Q,
//         [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//         {
//         	Real const xFluxSum
//         	= (xFlux_xF(i,j,k,n) - xFlux_xF(i + 1,j,    k,    n))*dSx*DVM::xi_u[n]
// 			+ (xFlux_yF(i,j,k,n) - xFlux_yF(i,    j + 1,k,    n))*dSy*DVM::xi_v[n]
// #if(AMREX_SPACEDIM == 3)
// 			+ (xFlux_zF(i,j,k,n) - xFlux_zF(i,    j    ,k + 1,n))*dSz*DVM::xi_w[n]
// #endif			
// 			;

// 			const Real xneq = (xtau(i,j,k)-hdTime)/xtau(i,j,k);
// 			const Real xWavePlus = xEq(i,j,k,n) + xneq*xWave(i,j,k,n) + hdTime*xSo(i,j,k,n);

// 			xWave(i,j,k,n) = xWavePlus + xFluxSum*dTime/volume;
//         });
// 	}
// }

template <typename DVM>
void DUGKSCore<DVM>::Update_Wave_SS (MultiFab &xWaveFab, const MultiFab &xSoFab, const Real ssTime)
{
	amrex::MultiFab::Saxpy(xWaveFab, ssTime, xSoFab, 0, 0, DVM::Q, 0);
}

template<typename DVM>
void DVMCore<DVM>::Update_Wave (Array<MultiFab, Ark::DVDFNUM> &xDVDF, const Real dTime, 
								const Real minL)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(xDVDF[Ark::WAVE], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		const amrex::Box& bx = mfi.tilebox();

		Array4<Real 	 > const &xWave = xDVDF[Ark::WAVE].array(mfi);

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

		LoopVS(DVM::Q)
		{
			const Real dlx  = (minL - dTime)*DVM::xi_u[n];
			const Real dlx2 = dlx*dlx;

			const Real dly  = (minL - dTime)*DVM::xi_v[n];
			const Real dly2 = dly*dly;

			const Real dlxy = dlx*dly;

#if(AMREX_SPACEDIM == 3)
			const Real dlz = (minL - dTime)*DVM::xi_w[n];
			const Real dlz2 = dlz*dlz;
			const Real dlyz = dly*dlz;
			const Real dlzx = dlz*dlx;
#endif			
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				int const ii = i - static_cast<int>(DVM::xi_u[n]);
				int const jj = j - static_cast<int>(DVM::xi_v[n]);
#if(AMREX_SPACEDIM == 3)
				int const kk = k - static_cast<int>(DVM::xi_w[n]);

				xWave(i,j,k,n) = xBP(ii,jj,kk,n)
							   + xBPx(ii,jj,kk,n)*dlx 
							   + xBPy(ii,jj,kk,n)*dly
							   + xBPz(ii,jj,kk,n)*dlz
							   + (xBPxx(ii,jj,kk,n)*dlx2 
							     +xBPyy(ii,jj,kk,n)*dly2
							     +xBPzz(ii,jj,kk,n)*dlz2)*0.5
							   + xBPxy(ii,jj,kk,n)*dlxy
							   + xBPyz(ii,jj,kk,n)*dlyz
							   + xBPzx(ii,jj,kk,n)*dlzx;
#elif(AMREX_SPACEDIM == 2)
				xWave(i,j,k,n) = xBP(ii,jj,k,n)
							   + xBPx(ii,jj,k,n)*dlx 
							   + xBPy(ii,jj,k,n)*dly
							   + (xBPxx(ii,jj,k,n)*dlx2 
							   	 +xBPyy(ii,jj,k,n)*dly2)*0.5
							   + xBPxy(ii,jj,k,n)*dlxy;				   
#endif				
			});	
		}
	}
}
//
#endif