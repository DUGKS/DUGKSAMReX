#include <DUGKSAmr.H>

//symmetric
double Limiter::VanAlbada2(const double vL, const double vR)
{
	const int Sign = Ark::sign(vL)*Ark::sign(vR);

	return Sign > 0 ? 2.0/(vL/vR + vR/vL): 0;
}
double Limiter::VanLeer(const double vL, const double vR)
{
	const int Sign = Ark::sign(vL)*Ark::sign(vR);

	return Sign > 0 ? 2.0*vL/(vR + vL): 0;
}