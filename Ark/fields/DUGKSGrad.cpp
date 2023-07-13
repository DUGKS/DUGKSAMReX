#include <DUGKSAmr.H>

void GradFunc::Update_Grad_x (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i+1,j,k,n) - src(i-1,j,k,n))/(2.0*dlx);
	});
}

void GradFunc::Update_Grad_VanLeer_x (const Box &bx, int const ncomp, Real const *dx,
							  	  	  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		const double vF = src(i  ,j,k,n) - src(i-1,j,k,n);
		const double vB = src(i+1,j,k,n) - src(i  ,j,k,n);

		int const tokenF = Ark::sign(vF);
		int const tokenB = Ark::sign(vB);

		dst(i,j,k,n) = tokenF*tokenB > 0 ? 2*vF*vB/((vF + vB)*dlx) : 0;
	});
}
void GradFunc::Update_Grad_x_5th (const Box &bx, int const ncomp, Real const *dx,
							  	  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i+1,j,k,n) - src(i-1,j,k,n))*8.
		-	(src(i+2,j,k,n) - src(i-2,j,k,n))
		)/(12.0*dlx);
	});
}
void GradFunc::Update_Grad_y (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i,j+1,k,n) - src(i,j-1,k,n))/(2.0*dly);
	});
}
void GradFunc::Update_Grad_VanLeer_y (const Box &bx, int const ncomp, Real const *dx,
							  	  	  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		const double vF = src(i,j  ,k,n) - src(i,j-1,k,n);
		const double vB = src(i,j+1,k,n) - src(i,j  ,k,n);

		int const tokenF = Ark::sign(vF);
		int const tokenB = Ark::sign(vB);

		dst(i,j,k,n) = tokenF*tokenB > 0 ? 2*vF*vB/((vF + vB)*dly) : 0;
	});
}
void GradFunc::Update_Grad_y_5th (const Box &bx, int const ncomp, Real const *dx,
							  	  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i,j+1,k,n) - src(i,j-1,k,n))*8.
		-   (src(i,j+2,k,n) - src(i,j-2,k,n))
		)/(12.0*dly);
	});
}
void GradFunc::Update_Grad_xx (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];
	double const dlx2 = dlx*dlx;

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i-1,j,k,n) - 2.0*src(i,j,k,n) + src(i+1,j,k,n)))/(dlx2);
	});
}
void GradFunc::Update_Grad_xx_5th (const Box &bx, int const ncomp, Real const *dx,
							   	   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];
	double const dlx2 = dlx*dlx;

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i-1,j,k,n) + src(i+1,j,k,n))*16
		-	(src(i-2,j,k,n) + src(i+2,j,k,n))
		-	30.0*src(i,j,k,n)
		)/(12.*dlx2);
	});
}

void GradFunc::Update_Grad_yy (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];
	double const dly2 = dly*dly;

	ParallelFor(bx,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i,j-1,k,n) - 2.0*src(i,j,k,n) + src(i,j+1,k,n)))/(dly2);
	});
}

void GradFunc::Update_Grad_yy_5th(const Box &bx, int const ncomp,Real const *dx,
							   	  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];
	double const dly2 = dly*dly;

	ParallelFor(bx,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i,j-1,k,n) + src(i,j+1,k,n))*16
		-	(src(i,j-2,k,n) + src(i,j+2,k,n))
		-	30.0*src(i,j,k,n)
		)/(12.*dly2);
	});
}

void GradFunc::Update_Grad_xy (const Box &bx, int const ncomp,Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dlx = dx[Ark::X];
	const Real dly = dx[Ark::Y];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			src(i+1,j+1,k,n) + src(i-1,j-1,k,n)
		-	src(i+1,j-1,k,n) - src(i-1,j+1,k,n)
		)/(4.0*dlx*dly);
	});
}

void GradFunc::Update_Grad_xxx (const Box &bx, int const ncomp, Real const *dx,
							    Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx  = dx[Ark::X];
	double const dlx3 = dlx*dlx*dlx;

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
		-	(src(i+1,j,k,n) - src(i-1,j,k,n))
		+	(src(i+2,j,k,n) - src(i-2,j,k,n))*0.5
		)/(dlx3);
	});
}
void GradFunc::Update_Grad_yyy (const Box &bx, int const ncomp, Real const *dx,
							    Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];
	double const dly3 = dly*dly*dly;

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
		-	(src(i,j+1,k,n) - src(i,j-1,k,n))
		+	(src(i,j+2,k,n) - src(i,j-2,k,n))*0.5
		)/(dly3);
	});
}

void GradFunc::Update_Grad_xxy (const Box &bx, int const ncomp,Real const *dx,
							    Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dlx = dx[Ark::X];
	const Real dly = dx[Ark::Y];
	const Real dlxxy = dlx*dlx*dly;

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i-1,j+1,k,n) - 2.0*src(i,j+1,k,n) + src(i+1,j+1,k,n))
		-	(src(i-1,j-1,k,n) - 2.0*src(i,j-1,k,n) + src(i+1,j-1,k,n))
		)/(2.0*dlxxy);
	});
}

void GradFunc::Update_Grad_xyy (const Box &bx, int const ncomp,Real const *dx,
							    Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dlx = dx[Ark::X];
	const Real dly = dx[Ark::Y];
	const Real dlxyy = dlx*dly*dly;
	
	ParallelFor(bx,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i+1,j-1,k,n) - 2.0*src(i+1,j,k,n) + src(i+1,j+1,k,n))
		-	(src(i-1,j-1,k,n) - 2.0*src(i-1,j,k,n) + src(i-1,j+1,k,n))
		)/(2.0*dlxyy);
	});
}

void GradFunc::Update_Grad_xxxx (const Box &bx, int const ncomp, Real const *dx,
							     Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];
	double const dlx4 = dlx*dlx*dlx*dlx;

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			src(i,j,k,n)*6.0
		-	(src(i+1,j,k,n) + src(i-1,j,k,n))*4.0
		+	(src(i+2,j,k,n) + src(i-2,j,k,n))
		)/dlx4;
	});
}
void GradFunc::Update_Grad_yyyy (const Box &bx, int const ncomp, Real const *dx,
							     Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];
	double const dly4 = dly*dly*dly*dly;

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			src(i,j,k,n)*6.0
		-	(src(i,j+1,k,n) + src(i,j-1,k,n))*4.0
		+	(src(i,j+2,k,n) + src(i,j-2,k,n))
		)/dly4;
	});
}
void GradFunc::Update_Grad_xxyy (const Box &bx, int const ncomp,Real const *dx,
							     Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dlx = dx[Ark::X];
	const Real dly = dx[Ark::Y];
	const Real dlxxyy = dlx*dlx*dly*dly;

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i+1,j+1,k,n) - 2.0*src(i+1,j,k,n) + src(i+1,j-1,k,n))
		-	(src(i  ,j+1,k,n) - 2.0*src(i  ,j,k,n) + src(i  ,j-1,k,n))*2
		+   (src(i-1,j+1,k,n) - 2.0*src(i-1,j,k,n) + src(i-1,j-1,k,n))
		)/dlxxyy;
	});
}

#if(AMREX_SPACEDIM == 3)

void GradFunc::Update_Grad_z (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::Z];

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i,j,k+1,n) - src(i,j,k-1,n))/(2.0*dlz);
	});
}
void GradFunc::Update_Grad_z_5th (const Box &bx, int const ncomp, Real const *dx,
							  	  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::Z];

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i,j,k+1,n) - src(i,j,k-1,n))*8.
		-   (src(i,j,k+2,n) - src(i,j,k-2,n))
		)/(12.0*dlz);
	});
}
void GradFunc::Update_Grad_VanLeer_z (const Box &bx, int const ncomp, Real const *dx,
							  	  	  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::Z];

	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		const double vF = src(i,j,k  ,n) - src(i,j,k-1,n);
		const double vB = src(i,j,k+1,n) - src(i,j,k  ,n);

		int const tokenF = Ark::sign(vF);
		int const tokenB = Ark::sign(vB);

		dst(i,j,k,n) = tokenF*tokenB > 0 ? 2*vF*vB/((vF + vB)*dlz) : 0;
	});
}
void GradFunc::Update_Grad_zz (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::Z];

	ParallelFor(bx,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i,j,k-1,n) - 2.0*src(i,j,k,n) + src(i,j,k+1,n)))/(dlz*dlz);
	});
}

void GradFunc::Update_Grad_zz_5th(const Box &bx, int const ncomp,Real const *dx,
							   	  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::Z];
	double const dlz2 = dlz*dlz;

	ParallelFor(bx,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i,j,k-1,n) + src(i,j,k+1,n))*16
		-	(src(i,j,k-2,n) + src(i,j,k+2,n))
		-	30.0*src(i,j,k,n)
		)/(12.*dlz2);
	});
}
void GradFunc::Update_Grad_yz (const Box &bx, int const ncomp,Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dly = dx[Ark::Y];
	const Real dlz = dx[Ark::Z];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			src(i,j+1,k+1,n) + src(i,j-1,k-1,n)
		-	src(i,j+1,k-1,n) - src(i,j-1,k+1,n)
		)/(4.0*dly*dlz);
	});
}
void GradFunc::Update_Grad_zx (const Box &bx, int const ncomp,Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dlz = dx[Ark::Z];
	const Real dlx = dx[Ark::X];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			src(i+1,j,k+1,n) + src(i-1,j,k-1,n)
		-	src(i-1,j,k+1,n) - src(i+1,j,k-1,n)
		)/(4.0*dlz*dlx);
	});
}

// void GradFunc::Update_Grad_zzz (const Box &bx, int const ncomp, Real const *dx,
// 							    Array4<Real const> const &src, Array4<Real> const &dst)
// {
// 	double const dlz  = dx[Ark::Z];
// 	double const dlz3 = dlz*dlz*dlz;

// 	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
// 	{
// 		dst(i,j,k,n) = 
// 		(
// 		-	(src(i,j,k+1,n) - src(i,j,k-1,n))
// 		+	(src(i,j,k+2,n) - src(i,j,k-2,n))*0.5
// 		)/(dlz3);
// 	});
// }

// void GradFunc::Update_Grad_zzzz (const Box &bx, int const ncomp, Real const *dx,
// 							     Array4<Real const> const &src, Array4<Real> const &dst)
// {
// 	double const dlz  = dx[Ark::Z];
// 	double const dlz4 = dlz*dlz*dlz*dlz;

// 	ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
// 	{
// 		dst(i,j,k,n) = 
// 		(
// 			src(i,j,k,n)*6.0
// 		-	(src(i,j,k+1,n) + src(i,j,k-1,n))*4.0
// 		+	(src(i,j,k+2,n) + src(i,j,k-2,n))
// 		)/dlz4;
// 	});
// }

void GradFunc::gradient_operator_isotropy (const Box &bx, Real const* dl, 
								 		   Array4<Real const> const &src,
								 		   Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];

	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dlx == dlz);

	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::X) = 
		(
			(src(i+1,j,k) - src(i-1,j,k))*16

		+	(
				src(i+1,j+1,k) + src(i+1,j-1,k) + src(i+1,j,k+1) + src(i+1,j,k-1) 
			-	src(i-1,j+1,k) - src(i-1,j-1,k) - src(i-1,j,k+1) - src(i-1,j,k-1) 
			)*4

		+	(
			  	src(i+1,j+1,k+1) + src(i+1,j-1,k+1) + src(i+1,j+1,k-1) + src(i+1,j-1,k-1)
			- 	src(i-1,j+1,k+1) - src(i-1,j-1,k+1) - src(i-1,j+1,k-1) - src(i-1,j-1,k-1)
			)
		)/(72*dlx);
	});
	//
	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Y) =	
		(
			(src(i,j+1,k) - src(i,j-1,k))*16

		+	(
				src(i+1,j+1,k) + src(i-1,j+1,k) + src(i,j+1,k+1) + src(i,j+1,k-1) 
			-	src(i+1,j-1,k) - src(i-1,j-1,k) - src(i,j-1,k+1) - src(i,j-1,k-1) 
			)*4

		+	(
			  	src(i+1,j+1,k+1) + src(i-1,j+1,k+1) + src(i+1,j+1,k-1) + src(i-1,j+1,k-1)
			- 	src(i+1,j-1,k+1) - src(i-1,j-1,k+1) - src(i+1,j-1,k-1) - src(i-1,j-1,k-1)
			)
		)/(72*dly);
	});
	//
	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Z) =
		(
			(src(i,j,k+1) - src(i,j,k-1))*16

		+	(
				src(i,j+1,k+1) + src(i,j-1,k+1) + src(i+1,j,k+1) + src(i-1,j,k+1) 
			-	src(i,j+1,k-1) - src(i,j-1,k-1) - src(i+1,j,k-1) - src(i-1,j,k-1) 
			)*4

		+	(
			  	src(i+1,j+1,k+1) + src(i+1,j-1,k+1) + src(i-1,j+1,k+1) + src(i-1,j-1,k+1)
			- 	src(i+1,j+1,k-1) - src(i+1,j-1,k-1) - src(i-1,j+1,k-1) - src(i-1,j-1,k-1)
			)
		)/(72*dlz);
	});
}

void GradFunc::gradient_operator_5th (const Box &bx, Real const* dl, 
								 	  Array4<Real const> const &src,
								 	  Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];

	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dlx == dlz);

	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx,	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::X) = 
		(
			(src(i+1,j,k) - src(i-1,j,k))*8.
		-	(src(i+2,j,k) - src(i-2,j,k))
		)/(12.0*dlx);
	});
	//
	ParallelFor(bx,	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Y) = 
		(
			(src(i,j+1,k) - src(i,j-1,k))*8.
		-   (src(i,j+2,k) - src(i,j-2,k))
		)/(12.0*dly);
	});
	//
	ParallelFor(bx,	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Z) = 
		(
			(src(i,j,k+1) - src(i,j,k-1))*8.
		-   (src(i,j,k+2) - src(i,j,k-2))
		)/(12.0*dlz);
	});
}

//doi:10.1016/j.jcp.2004.05.005	  		   
void GradFunc::laplacian_operator_isotropy (const Box &bx, Real const *dl, 
										    Array4<Real const> const &src, 
					  	   					Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];
	
	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dly == dlz);

	const Real dln = dlx;

	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k) = 
		(
			(
			  	src(i+1,j  ,k)   + src(i-1,j  ,k)
			+ 	src(i,  j+1,k)   + src(i,  j-1,k) 
			+ 	src(i,  j  ,k+1) + src(i,  j  ,k-1)
			)*20.
		+	(
			  	src(i+1,j+1,k) + src(i-1,j+1,k) + src(i+1,j-1,k) + src(i-1,j-1,k)
			+ 	src(i+1,j,k+1) + src(i-1,j,k+1) + src(i+1,j,k-1) + src(i-1,j,k-1)
			+ 	src(i,j+1,k+1) + src(i,j-1,k+1) + src(i,j+1,k-1) + src(i,j-1,k-1)
			)*6.
		+	(
				src(i+1,j+1,k+1) + src(i+1,j-1,k+1) + src(i+1,j+1,k-1) + src(i+1,j-1,k-1)
			+	src(i-1,j+1,k+1) + src(i-1,j-1,k+1) + src(i-1,j+1,k-1) + src(i-1,j-1,k-1)
			)
		-   src(i,j,k)*200
		)/(48*dln*dln);
	});
}
//
void GradFunc::laplacian_operator_5th (const Box &bx, Real const *dl, 
									   Array4<Real const> const &src, 
					  	   			   Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];
	
	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dly == dlz);

	const Real dln  = dlx;
	const Real dln2 = dln*dln;

	ParallelFor(bx,	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k) =  
		(
			(
				src(i-1,j  ,k  ) + src(i+1,j,  k  ) 
			+ 	src(i  ,j-1,k  ) + src(i,  j+1,k  )
			+ 	src(i  ,j  ,k-1) + src(i,  j  ,k+1)
			)*16
		-	
			(
				src(i-2,j  ,k  ) + src(i+2,j  ,k  ) 
			+	src(i  ,j-2,k  ) + src(i  ,j+2,k  )
			+	src(i  ,j  ,k-2) + src(i  ,j  ,k+2)
			)
		-	src(i,j,k)*90.
		)/(12.*dln2);
	});
}

void GradFunc::laplacian_normal_operator_isotropy (const Box &bx, Real const* dl,
												   Array4<Real const> const &src,	 
								 				   Array4<Real const> const &gradsrc,
								 				   Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];
	
	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dly == dlz);

	const Real dln  = dlx;
	const Real dln2 = dln*dln;

	amrex::Abort("laplacian_normal_operator_isotropy to be done...");
}
void GradFunc::laplacian_normal_operator_5th (const Box &bx, Real const* dl,
											  Array4<Real const> const &src,	 
								 			  Array4<Real const> const &gradsrc,
								 			  Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];
	
	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dly == dlz);

	const Real dln  = dlx;
	const Real dln2 = dln*dln;

	amrex::Abort("laplacian_normal_operator_5th to be done...");
}

#elif(AMREX_SPACEDIM == 2)

void GradFunc::gradient_operator_central (const Box &bx, Real const* dl, 
								 		   Array4<Real const> const &src,
								 		   Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::X) = ((src(i+1,j,k) - src(i-1,j,k)))/(2.0*dlx);
	});
	//
	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Y) = ((src(i,j+1,k) - src(i,j-1,k)))/(2.0*dly);
	});
}

void GradFunc::gradient_operator_isotropy (const Box &bx, Real const* dl, 
								 		   Array4<Real const> const &src,
								 		   Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::X) = 
		(
			(src(i+1,j  ,k) - src(i-1,j  ,k))*4.0
		+	(src(i+1,j-1,k) - src(i-1,j-1,k))
		+	(src(i+1,j+1,k) - src(i-1,j+1,k))
		)/(12.0*dlx);
	});
	//
	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Y) =
		(
			(src(i,  j+1,k) - src(i,  j-1,k))*4.0
		+	(src(i-1,j+1,k) - src(i-1,j-1,k))
		+	(src(i+1,j+1,k) - src(i+1,j-1,k))
		)/(12.0*dly);
	});
}

void GradFunc::gradient_operator_isotropy_blend (const Box &bx, Real const* dl, 
								 		 		 Array4<Real const> const &src,
								 		 		 Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		Real biased_x = 
		(
			(src(i+1,j  ,k) - src(i-1,j  ,k))*16.0 - (src(i+2,j  ,k) - src(i-2,j  ,k))*4.0
		+	(src(i+1,j-1,k) - src(i-1,j-1,k))*4.0  - (src(i+2,j-2,k) - src(i-2,j-2,k))
		+	(src(i+1,j+1,k) - src(i-1,j+1,k))*4.0  - (src(i+2,j+2,k) - src(i-2,j+2,k))
		)/(24.0*dlx);

		Real central_x = 
		(
			(src(i+1,j  ,k) - src(i-1,j  ,k))*4.0
		+	(src(i+1,j-1,k) - src(i-1,j-1,k))
		+	(src(i+1,j+1,k) - src(i-1,j+1,k))
		)/(12.0*dlx);

		dst(i,j,k,Ark::X) = 0.5*(biased_x + central_x);
	});
	//
	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		Real biased_y =
		(
			(src(i,  j+1,k) - src(i,  j-1,k))*16.0 - (src(i,  j+2,k) - src(i,  j-2,k))*4.0
		+	(src(i-1,j+1,k) - src(i-1,j-1,k))*4.0  - (src(i-2,j+2,k) - src(i-2,j-2,k))
		+	(src(i+1,j+1,k) - src(i+1,j-1,k))*4.0  - (src(i+2,j+2,k) - src(i+2,j-2,k))
		)/(24.0*dly);

		Real central_y =
		(
			(src(i,  j+1,k) - src(i,  j-1,k))*4.0
		+	(src(i-1,j+1,k) - src(i-1,j-1,k))
		+	(src(i+1,j+1,k) - src(i+1,j-1,k))
		)/(12.0*dly);

		dst(i,j,k,Ark::Y) = 0.5*(biased_y + central_y);
	});
}

void GradFunc::gradient_operator_isotropy_biased (const Box &bx, Real const* dl, 
								 		 		 Array4<Real const> const &src,
								 		 		 Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::X) = 
		(
			(src(i+1,j  ,k) - src(i-1,j  ,k))*16.0 - (src(i+2,j  ,k) - src(i-2,j  ,k))*4.0
		+	(src(i+1,j-1,k) - src(i-1,j-1,k))*4.0  - (src(i+2,j-2,k) - src(i-2,j-2,k))
		+	(src(i+1,j+1,k) - src(i-1,j+1,k))*4.0  - (src(i+2,j+2,k) - src(i-2,j+2,k))
		)/(24.0*dlx);
	});
	//
	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Y) =
		(
			(src(i,  j+1,k) - src(i,  j-1,k))*16.0 - (src(i,  j+2,k) - src(i,  j-2,k))*4.0
		+	(src(i-1,j+1,k) - src(i-1,j-1,k))*4.0  - (src(i-2,j+2,k) - src(i-2,j-2,k))
		+	(src(i+1,j+1,k) - src(i+1,j-1,k))*4.0  - (src(i+2,j+2,k) - src(i+2,j-2,k))
		)/(24.0*dly);
	});
}

void GradFunc::gradient_operator_5th (const Box &bx, Real const* dl, 
								 	  Array4<Real const> const &src,
								 	  Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::X) = 
		(
			(src(i+1,j,k) - src(i-1,j,k))*8.
		-	(src(i+2,j,k) - src(i-2,j,k))
		)/(12.0*dlx);
	});
	//
	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Y) = 
		(
			(src(i,j+1,k) - src(i,j-1,k))*8.
		-   (src(i,j+2,k) - src(i,j-2,k))
		)/(12.0*dly);
	});
}
void GradFunc::laplacian_operator_isotropy (const Box &bx, Real const *dl, 
										    Array4<Real const> const &src, 
					  	   					Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlxy = dlx*dly;

	AMREX_ASSERT(dlx == dly);

	ParallelFor(bx,	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k) = 
		(
			(src(i+1,j,  k) + src(i-1,j,  k) + src(i,  j+1,k) + src(i,  j-1,k))*4.0
		+	(src(i+1,j+1,k) + src(i-1,j+1,k) + src(i-1,j-1,k) + src(i+1,j-1,k))
		-	src(i,j,k)*20
		)/(6.0*dlxy);
	});
}
double GradFunc::laplacian_operator_isotropy (Real const *dl, Array4<Real const> const &var, 
					  	   					  int const i,int const j,int const k)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlxy = dlx*dly;
	
	AMREX_ASSERT(dlx == dly);
//
	return
	(	
		(var(i+1,j,  k) + var(i-1,j,  k) + var(i,j+1,  k) + var(i,j-1,  k))*4.0
	+	(var(i+1,j+1,k) + var(i-1,j+1,k) + var(i-1,j-1,k) + var(i+1,j-1,k))
	-	var(i,j,k)*20
	)/6.0*dlxy;

}
void GradFunc::laplacian_operator_5th (const Box &bx, Real const *dl, 
									   Array4<Real const> const &src, 
					  	   			   Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlxy = dlx*dly;

	AMREX_ASSERT(dlx == dly);

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k) =  
		(
			(src(i-1,j,k) + src(i+1,j,k) + src(i,j-1,k) + src(i,j+1,k))*16
		-	(src(i-2,j,k) + src(i+2,j,k) + src(i,j-2,k) + src(i,j+2,k))
		-	src(i,j,k)*60.
		)/(12.*dlxy);
	});
}
// nabla^2{_n}{\phi}
void GradFunc::laplacian_normal_operator_isotropy (const Box &bx, Real const* dl,
												Array4<Real const> const &src,	 
								 				Array4<Real const> const &gradsrc,
								 				Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	const Real dlx2 = dlx*dlx;
	const Real dly2 = dly*dly;
	const Real dlxy = dlx*dly;

	ParallelFor(bx,	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		const Real src_x = gradsrc(i,j,k,Ark::X);
		const Real src_y = gradsrc(i,j,k,Ark::Y);

		const double denominator = Ark::Square(src_x,src_y);

		if(Ark::tinyEqual(denominator,0.0,1E-30))
		{
			dst(i,j,k) = 0;
		}
		else
		{
			const Real src_xx = 
			(
				(src(i-1,j,  k) - 2.0*src(i,j,  k) + src(i+1,j,  k))*10.0
			+	(src(i-1,j+1,k) - 2.0*src(i,j+1,k) + src(i+1,j+1,k))
			+	(src(i-1,j-1,k) - 2.0*src(i,j-1,k) + src(i+1,j-1,k))
			)/(12.0*dlx2);

			const Real src_yy =
			(
				(src(i,  j-1,k) - 2.0*src(i,  j,k) + src(i,  j+1,k))*10.0
			+	(src(i+1,j-1,k) - 2.0*src(i+1,j,k) + src(i+1,j+1,k))
			+	(src(i-1,j-1,k) - 2.0*src(i-1,j,k) + src(i-1,j+1,k))
			)/(12.0*dly2);

			const Real src_xy = 
			(
				src(i+1,j+1,k) + src(i-1,j-1,k)
			-	src(i+1,j-1,k) - src(i-1,j+1,k)
			)/(4.0*dlxy);

			const Real numerator = src_x*src_x*src_xx + src_y*src_y*src_yy + 2.*src_x*src_y*src_xy;

			dst(i,j,k) = numerator/denominator;
		}
	});
}

void GradFunc::laplacian_normal_operator_5th (const Box &bx, Real const* dl,
											  Array4<Real const> const &src,	 
								 			  Array4<Real const> const &gradsrc,
								 			  Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	const Real dlx2 = dlx*dlx;
	const Real dly2 = dly*dly;
	const Real dlxy = dlx*dly;

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		const Real src_x = gradsrc(i,j,k,Ark::X);
		const Real src_y = gradsrc(i,j,k,Ark::Y);

		const double denominator = Ark::Square(src_x,src_y);

		if(Ark::tinyEqual(denominator,0.0,1E-30))
		{
			dst(i,j,k) = 0;
		}
		else
		{
			const Real src_xx = 
			(
				(src(i-1,j,k) + src(i+1,j,k))*16
			-	(src(i-2,j,k) + src(i+2,j,k))
			-	src(i,j,k)*30.0
			)/(12.*dlx2);

			const Real src_yy =
			(
				(src(i,j-1,k) + src(i,j+1,k))*16
			-	(src(i,j-2,k) + src(i,j+2,k))
			-	src(i,j,k)*30.
			)/(12.*dly2);

			const Real src_xy = 
			(
				src(i+1,j+1,k) + src(i-1,j-1,k)
			-	src(i+1,j-1,k) - src(i-1,j+1,k)
			)/(4.0*dlxy);

			const Real numerator = 
			(src_y*src_xy - src_x*src_yy)*src_x + (src_x*src_xy - src_y*src_xx)*src_y;

			dst(i,j,k) = numerator/denominator;
		}
	});
}

void GradFunc::laplacian_normal_chemip_isotropy (const Box &bx, Real const* dl,
												   Array4<Real const> const &src,	 
								 				   Array4<Real const> const &gradsrc,
								 				   Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	const Real dlx2 = dlx*dlx;
	const Real dly2 = dly*dly;
	const Real dlxy = dlx*dly;

	ParallelFor(bx,	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		const Real src_x = gradsrc(i,j,k,Ark::X);
		const Real src_y = gradsrc(i,j,k,Ark::Y);

		const double denominator = Ark::Square(src_x,src_y);

		if(Ark::tinyEqual(denominator,0.0,1E-30))
		{
			dst(i,j,k) = 0;
		}
		else
		{
			const Real src_xx = 
			(
				(src(i-1,j,  k) - 2.0*src(i,j,  k) + src(i+1,j,  k))*10.0
			+	(src(i-1,j+1,k) - 2.0*src(i,j+1,k) + src(i+1,j+1,k))
			+	(src(i-1,j-1,k) - 2.0*src(i,j-1,k) + src(i+1,j-1,k))
			)/(12.0*dlx2);

			const Real src_yy =
			(
				(src(i,  j-1,k) - 2.0*src(i,  j,k) + src(i,  j+1,k))*10.0
			+	(src(i+1,j-1,k) - 2.0*src(i+1,j,k) + src(i+1,j+1,k))
			+	(src(i-1,j-1,k) - 2.0*src(i-1,j,k) + src(i-1,j+1,k))
			)/(12.0*dly2);

			const Real src_xy = 
			(
				src(i+1,j+1,k) + src(i-1,j-1,k)
			-	src(i+1,j-1,k) - src(i-1,j+1,k)
			)/(4.0*dlxy);

			const Real numerator = (src_y*src_xy - src_x*src_yy)*src_x + (src_x*src_xy - src_y*src_xx)*src_y;

			dst(i,j,k) = numerator/denominator;
		}
	});
}
#endif

