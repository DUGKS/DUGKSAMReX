#ifndef ARK_DVM_DmVn_H
#define ARK_DVM_DmVn_H

namespace Ark{

struct D2V5
{
static const unsigned D = 2;

static const unsigned Q = 5;

static const double MaxU;

static const double MaSpan;//sqrt(3.0*RT)

static const double RT;

static const double Lambda;

static const double A;
static const double B;

static const double omega[Q];

static const double xi_u[Q];

static const double xi_v[Q];

static const double xi_w[Q];

static const unsigned BB[Q];

};

struct D2V9
{

static const unsigned D = 2;

static const unsigned Q = 9;

static const double MaxU;

static const double MaSpan;//sqrt(3.0*RT)

static const double RT;

static const double Lambda;

static const double A;
static const double B;
static const double C;

static const double E;
static const double F;

static const double omega[Q];

static const double xi_u[Q];

static const double xi_v[Q];

static const double xi_w[Q];

static const unsigned BB[Q];

};

struct D3V7
{

static const unsigned D = 3;

static const unsigned Q = 7;

static const double MaxU;

static const double MaSpan;//sqrt(3.0*RT)

static const double RT;

static const double Lambda;

static const double A;
static const double B;

static const double omega[Q];

static const double xi_u[Q];

static const double xi_v[Q];

static const double xi_w[Q];

static const unsigned BB[Q];

};

struct D3V15
{
static const unsigned D = 3;

static const unsigned Q = 15;

static const double MaxU;

static const double MaSpan;//sqrt(3.0*RT)

static const double RT;

static const double Lambda;

static const double A;
static const double B;
static const double C;

static const double omega[Q];

static const double xi_u[Q];

static const double xi_v[Q];

static const double xi_w[Q];

static const unsigned BB[Q];

};

}

#endif