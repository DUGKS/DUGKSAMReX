#include <DUGKSAmr.H>

#if(AMREX_SPACEDIM == 3)
void GradFunc::Update_Grad_x (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i+1,j,k,n) - src(i-1,j,k,n))/(2.0*dlx);
	});
}
void GradFunc::Update_Grad_y (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i,j+1,k,n) - src(i,j-1,k,n))/(2.0*dly);
	});
}
void GradFunc::Update_Grad_z (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::Z];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i,j,k+1,n) - src(i,j,k-1,n))/(2.0*dlz);
	});
}
void GradFunc::Update_Grad_xx (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i-1,j,k,n) - 2.0*src(i,j,k,n) + src(i+1,j,k,n)))/(dlx*dlx);
	});
}
void GradFunc::Update_Grad_yy (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i,j-1,k,n) - 2.0*src(i,j,k,n) + src(i,j+1,k,n)))/(dly*dly);
	});
}
void GradFunc::Update_Grad_zz (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::Z];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i,j,k-1,n) - 2.0*src(i,j,k,n) + src(i,j,k+1,n)))/(dlz*dlz);
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

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
double GradFunc::laplacian_operator_isotropy (Real const *dl, Array4<Real const> const &var, 
					  	   					  int const i,int const j,int const k)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];
	
	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dly == dlz);

	const Real dln = dlx;
//	
	Real tmp =
			(
				(
				  	var(i+1,j  ,k)   + var(i-1,j  ,k)
				+ 	var(i,  j+1,k)   + var(i,  j-1,k) 
				+ 	var(i,  j  ,k+1) + var(i,  j  ,k-1)
				)*20.
	//
			+	(
				  	var(i+1,j+1,k) + var(i-1,j+1,k) + var(i+1,j-1,k) + var(i-1,j-1,k)
				+ 	var(i+1,j,k+1) + var(i-1,j,k+1) + var(i+1,j,k-1) + var(i-1,j,k-1)
				+ 	var(i,j+1,k+1) + var(i,j-1,k+1) + var(i,j+1,k-1) + var(i,j-1,k-1)
				)*6.
	//
			+	(
					var(i+1,j+1,k+1) + var(i+1,j-1,k+1) + var(i+1,j+1,k-1) + var(i+1,j-1,k-1)
				+	var(i-1,j+1,k+1) + var(i-1,j-1,k+1) + var(i-1,j+1,k-1) + var(i-1,j-1,k-1)
				)
			-   var(i,j,k)*200
			)/(48*dln*dln);

	return tmp;
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

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k) = 
		(
			(
			  	src(i+1,j  ,k)   + src(i-1,j  ,k)
			+ 	src(i,  j+1,k)   + src(i,  j-1,k) 
			+ 	src(i,  j  ,k+1) + src(i,  j  ,k-1)
			)*20.
//
		+	(
			  	src(i+1,j+1,k) + src(i-1,j+1,k) + src(i+1,j-1,k) + src(i-1,j-1,k)
			+ 	src(i+1,j,k+1) + src(i-1,j,k+1) + src(i+1,j,k-1) + src(i-1,j,k-1)
			+ 	src(i,j+1,k+1) + src(i,j-1,k+1) + src(i,j+1,k-1) + src(i,j-1,k-1)
			)*6.
//
		+	(
				src(i+1,j+1,k+1) + src(i+1,j-1,k+1) + src(i+1,j+1,k-1) + src(i+1,j-1,k-1)
			+	src(i-1,j+1,k+1) + src(i-1,j-1,k+1) + src(i-1,j+1,k-1) + src(i-1,j-1,k-1)
			)
		-   src(i,j,k)*200
		)/(48*dln*dln);
	});
}
#endif

double GradFunc::central_interpolation (const double upstream1, const double upstream0,
                            		    const double downstream0, const double downstream1,
										int const token)
{
	// -----------------------------------------
	// |up2 | up1 | up0 | down0 | down1 | down2|
	// -----------------------------------------
	// Interface value = o*(up0 + down0) + n*(up1 + down1)
	// o + n == 1/2; o + 9n == -token/64

	switch (token)
	{
		case 6040:
        {
            return (0.6*upstream0 + 0.4*downstream0);
            break;
        }
		case -32:
        {
            return (upstream0 + downstream0)/2.0;
            break;
        }
		case 0:
        {
            return ( (upstream0 + downstream0)*9.0 - (upstream1 + downstream1)*1)/16.0;
            break;
        }
		case 16:
        {
            return ( (upstream0 + downstream0)*19.0 - (upstream1 + downstream1)*3)/32.0;
            break;
        }
		case 24:
		{
			return ( (upstream0 + downstream0)*39.0 - (upstream1 + downstream1)*7)/64.0;
            break;
		}
		case 32:
        {
            return ( (upstream0 + downstream0)*5.0 - (upstream1 + downstream1)*1)/8.0;
            break;
        }
        case 40:
        {
        	return ( (upstream0 + downstream0)*41.0 - (upstream1 + downstream1)*9)/64.0;
            break;
        }
        case 42:
        {
        	return ( (upstream0 + downstream0)*165.0 - (upstream1 + downstream1)*37)/256.0;
            break;
        }
        case 44:
        {
        	return ( (upstream0 + downstream0)*83.0 - (upstream1 + downstream1)*19)/128.0;
            break;
        }
		case 46:
        {
            return ( (upstream0 + downstream0)*167.0 - (upstream1 + downstream1)*39)/256.0;
            break;
        }
		case 48:
        {
            return ( (upstream0 + downstream0)*21.0 - (upstream1 + downstream1)*5)/32.0;
            break;
        }
		case 50:
        {
            return ( (upstream0 + downstream0)*169.0 - (upstream1 + downstream1)*41)/256.0;
            break;
        }
		case 52:
        {
            return ( (upstream0 + downstream0)*85.0 - (upstream1 + downstream1)*21)/128.0;
            break;
        }
        case 56:
        {
            return ( (upstream0 + downstream0)*43.0 - (upstream1 + downstream1)*11)/64.0;
            break;
        }
        case 60:
        {
            return ( (upstream0 + downstream0)*87.0 - (upstream1 + downstream1)*23)/128.0;
            break;
        }
		case 64:
        {
            return ( (upstream0 + downstream0)*11.0 - (upstream1 + downstream1)*3)/16.0;
            break;
        }
        case 68:
        {
            return ( (upstream0 + downstream0)*89.0 - (upstream1 + downstream1)*25)/128.0;
            break;
        }
        case 72:
        {
            return ( (upstream0 + downstream0)*45.0 - (upstream1 + downstream1)*13)/64.0;
            break;
        }
        case 76:
        {
            return ( (upstream0 + downstream0)*91.0 - (upstream1 + downstream1)*27)/128.0;
            break;
        }
        case 80:
        {
            return ( (upstream0 + downstream0)*23.0 - (upstream1 + downstream1)*7)/32.0;
            break;
        }
        case 84:
        {
        	return ( (upstream0 + downstream0)*93.0 - (upstream1 + downstream1)*29)/128.0;
        	break;
        }
        case 88:
        {
        	return ( (upstream0 + downstream0)*47.0 - (upstream1 + downstream1)*15)/64.0;
        	break;
        }
        case 92:
        {
        	return ( (upstream0 + downstream0)*95.0 - (upstream1 + downstream1)*31)/128.0;
        	break;
        }
        case 96:
        {
        	return ( (upstream0 + downstream0)*3.0 - (upstream1 + downstream1)*1)/4.0;
        	break;
        }
        case 100:
        {
        	return ( (upstream0 + downstream0)*97.0 - (upstream1 + downstream1)*33)/128.0;
        	break;
        }
        case 104:
        {
        	return ( (upstream0 + downstream0)*49.0 - (upstream1 + downstream1)*17)/64.0;
        	break;
        }
        case 108:
        {
        	return ( (upstream0 + downstream0)*99.0 - (upstream1 + downstream1)*35)/128.0;
        	break;
        }
        case 112:
        {
        	return ( (upstream0 + downstream0)*25.0 - (upstream1 + downstream1)*9)/32.0;
        	break;
        }
        case 116:
        {
        	return ( (upstream0 + downstream0)*101.0 - (upstream1 + downstream1)*37)/128.0;
        	break;
        }
        case 120:
        {
        	return ( (upstream0 + downstream0)*51.0 - (upstream1 + downstream1)*19)/64.0;
        	break;
        }
        case 124:
        {
        	return ( (upstream0 + downstream0)*103.0 - (upstream1 + downstream1)*39)/128.0;
        	break;
        }
        case 128:
        {
        	return ( (upstream0 + downstream0)*13.0 - (upstream1 + downstream1)*5)/16.0;
        	break;
        }
		default :
		{
			std::ostringstream ss;   
            ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
            ss  <<"invalid scheme token\n";
            amrex::Abort(ss.str());
		}
	}

	return -999999;
	// return ( (upstream2 + downstream2)*3 
    //        - (upstream1 + downstream1)*25.0
    //        + (upstream0 + downstream0)*150.0)/256.0;
}