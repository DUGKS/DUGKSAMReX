#ifndef ARK_DVM_DmVn_H
#define ARK_DVM_DmVn_H

namespace Ark{

struct D2V9
{

static const int Q = 9;

static const double MaxU;

static const double MaSpan;//sqrt(3.0*RT)

static const double RT;

static const double Lambda;

static const double A;
static const double B;
static const double C;
static const double D;
static const double E;
static const double F;

static const double omega[Q];

static const double xi_u[Q];

static const double xi_v[Q];

static const int BB[Q];

};

using DmVnPhase = D2V9;

using DmVnMomentum = D2V9;

}

#endif