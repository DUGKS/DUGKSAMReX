#include <DUGKSAmr.H>

inline 
double unnormalizedWeight (const double tau5, const double SI, const double epsilon)
{
    return 1 + tau5/(SI + epsilon);
}

double GradFunc::WENO_Z (const double upstream2,
                         const double upstream1,
                         const double upstream0,
                         const double downstream0,
                         const double downstream1)
{
    constexpr double A = 13.0/12.0;
    constexpr double B = 0.25;
    constexpr double epsilonWENO = 1e-16;
//
    double const targetU = (2.*upstream2 - 7.*upstream1 + 11.*upstream0)/6.;
    double const targetC = (-upstream1   + 5.*upstream0 + 2.*downstream0)/6.;
    double const targetD = (2.*upstream0 + 5.*downstream0 - downstream1)/6.;
//
    double const SIU = A*Ark::Square(upstream2 - 2.*upstream1 + upstream0) 
                     + B*Ark::Square(upstream2 - 4.*upstream1 + 3.*upstream0);
                     
    double const SIC = A*Ark::Square(upstream1 - 2.*upstream0 + downstream0) 
                     + B*Ark::Square(upstream1 - downstream0);
    
    double const SID = A*Ark::Square(upstream0    - 2.*downstream0 + downstream1)
                     + B*Ark::Square(3.*upstream0 - 4.*downstream0 + downstream1);

    double const tau5 = std::fabs(SIU - SID);

    double const uweightU = 0.1*unnormalizedWeight(tau5,SIU,epsilonWENO);
    double const uweightC = 0.6*unnormalizedWeight(tau5,SIC,epsilonWENO);
    double const uweightD = 0.3*unnormalizedWeight(tau5,SID,epsilonWENO);

    double const sumWeight = uweightU + uweightC + uweightD;
    double const weightU = uweightU/sumWeight;
    double const weightC = uweightC/sumWeight;
    double const weightD = uweightD/sumWeight;
    
    return weightU*targetU + weightC*targetC + weightD*targetD;
}

double GradFunc::WENO_C (const double upstream2, const double upstream1,
                         const double upstream0, const double downstream0,
                         const double downstream1, const double downstream2)
{
    return ( (upstream0 + downstream0)*37 
           - (upstream1 + downstream1)*8.0
           + (upstream2 + downstream2))/60.0;
}



double GradFunc::double_interpolation  (const double upstream1, const double upstream0,
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

double GradFunc::triple_interpolation  (const double upstream2, 
                                        const double upstream1,
                                        const double upstream0, 
                                        const double downstream0,
                                        const double downstream1, 
                                        const double downstream2,
                                        const int token)
{    
    // -----------------------------------------
    // |up2 | up1 | up0 | down0 | down1 | down2|
    // -----------------------------------------
    // Interface value = o*(up0 + down0) + n*(up1 + down1) + m*(up2 + down2)
    // o + n + m == 1/2; o + 9n + 25m == -token/64; o + 81n + 625m == 0
    //
    switch (token)
    {
        case -32:
        {
            return ((upstream0 + downstream0))/2.0;
            break;
        }
        case 0:
        {
            return ( (upstream0 + downstream0)*150 
                   - (upstream1 + downstream1)*25.0
                   + (upstream2 + downstream2)*3.0)/256.0;
            break;
        }
        case 16:
        {
            return ( (upstream0 + downstream0)*242 
                   - (upstream1 + downstream1)*57.0
                   + (upstream2 + downstream2)*7.0)/384.0;
            break;
        }
        case 24:
        {
            return ( (upstream0 + downstream0)*334 
                   - (upstream1 + downstream1)*89.0
                   + (upstream2 + downstream2)*11.0)/512.0;
            break;
        }
        case 28:
        {
            return ( (upstream0 + downstream0)*2038. 
                   - (upstream1 + downstream1)*573.
                   + (upstream2 + downstream2)*71.)/3072.;
            break;
        }
        case 30:
        {
            return ( (upstream0 + downstream0)*1370. 
                   - (upstream1 + downstream1)*395.
                   + (upstream2 + downstream2)*49.)/2048.;
            break;
        }
        case 32:
        {
            return ( (upstream0 + downstream0)*518. 
                   - (upstream1 + downstream1)*153.
                   + (upstream2 + downstream2)*19.)/768.0;
            break;
        }
        case 34:
        {
            return ( (upstream0 + downstream0)*4178. 
                   - (upstream1 + downstream1)*1263.
                   + (upstream2 + downstream2)*157.)/6144.;
            break;
        }
        case 36:
        {
            return ( (upstream0 + downstream0)*702. 
                   - (upstream1 + downstream1)*217.
                   + (upstream2 + downstream2)*27.)/1024.0;
            break;
        }
        case 40:
        {
            return ( (upstream0 + downstream0)*1070. 
                   - (upstream1 + downstream1)*345.
                   + (upstream2 + downstream2)*43.)/1536.0;
            break;
        }
        case 42:
        {
            return ( (upstream0 + downstream0)*1438 
                   - (upstream1 + downstream1)*473.0
                   + (upstream2 + downstream2)*59.0)/2048.0;
            break;
        }
        case 44:
        {
            return ( (upstream0 + downstream0)*2174. 
                   - (upstream1 + downstream1)*729.0
                   + (upstream2 + downstream2)*91.0)/3072.0;
            break;
        }
        case 46:
        {
            return ( (upstream0 + downstream0)*4382. 
                   - (upstream1 + downstream1)*1497.0
                   + (upstream2 + downstream2)*187.0)/6144.0;
            break;
        }
        case 48:
        {
            return ( (upstream0 + downstream0)*23 
                   - (upstream1 + downstream1)*8.0
                   + (upstream2 + downstream2)*1.0)/32.0;
            break;
        }
        case 50:
        {
            return ( (upstream0 + downstream0)*4450. 
                   - (upstream1 + downstream1)*1575.0
                   + (upstream2 + downstream2)*197.0)/6144.0;
            break;
        }
        case 52:
        {
            return ( (upstream0 + downstream0)*2242. 
                   - (upstream1 + downstream1)*807.
                   + (upstream2 + downstream2)*101.)/3072.;
            break;
        }
        case 56:
        {
            return ( (upstream0 + downstream0)*1138 
                   - (upstream1 + downstream1)*423.0
                   + (upstream2 + downstream2)*53.0)/1536.0;
            break;
        }
        case 60:
        {
            return ( (upstream0 + downstream0)*770.
                   - (upstream1 + downstream1)*295.0
                   + (upstream2 + downstream2)*37.0)/1024.;
            break;
        }
        case 64:
        {
            return ( (upstream0 + downstream0)*586 
                   - (upstream1 + downstream1)*231.0
                   + (upstream2 + downstream2)*29.0)/768.0;
            break;
        }        
        case 68:
        {
            return ( (upstream0 + downstream0)*2378. 
                   - (upstream1 + downstream1)*963.
                   + (upstream2 + downstream2)*121.)/3072.;
            break;
        }
        case 72:
        {
            return ( (upstream0 + downstream0)*402 
                   - (upstream1 + downstream1)*167.0
                   + (upstream2 + downstream2)*21.0)/512.0;
            break;
        }
        case 76:
        {
            return ( (upstream0 + downstream0)*2446 
                   - (upstream1 + downstream1)*1041.0
                   + (upstream2 + downstream2)*131.0)/3072.0;
            break;
        }
        case 80:
        {
            return ( (upstream0 + downstream0)*310 
                   - (upstream1 + downstream1)*135.0
                   + (upstream2 + downstream2)*17.0)/384.0;
            break;
        }
        case 84:
        {
            return ( (upstream0 + downstream0)*838 
                   - (upstream1 + downstream1)*373.0
                   + (upstream2 + downstream2)*47.0)/1024.0;
            break;
        }
        case 88:
        {
            return ( (upstream0 + downstream0)*1274 
                   - (upstream1 + downstream1)*579.0
                   + (upstream2 + downstream2)*73.0)/1536.0;
        }
        case 92:
        {
            return ( (upstream0 + downstream0)*2582 
                   - (upstream1 + downstream1)*1197.0
                   + (upstream2 + downstream2)*151.0)/3072.0;
        }
        case 96:
        {
            return ( (upstream0 + downstream0)*218
                   - (upstream1 + downstream1)*103.0
                   + (upstream2 + downstream2)*13.0)/256.0;
        }
        case 100:
        {
            return ( (upstream0 + downstream0)*2650 
                   - (upstream1 + downstream1)*1275.0
                   + (upstream2 + downstream2)*161.0)/3072.0;
        }
        case 104:
        {
            return ( (upstream0 + downstream0)*1342 
                   - (upstream1 + downstream1)*657.0
                   + (upstream2 + downstream2)*83.0)/1536.0;
        }
        case 108:
        {
            return ( (upstream0 + downstream0)*906 
                   - (upstream1 + downstream1)*451.0
                   + (upstream2 + downstream2)*57.0)/1024.0;
        }
        case 112:
        {
            return ( (upstream0 + downstream0)*172 
                   - (upstream1 + downstream1)*87.0
                   + (upstream2 + downstream2)*11.0)/192.0;
        }
        case 116:
        {
            return ( (upstream0 + downstream0)*2786 
                   - (upstream1 + downstream1)*1431.0
                   + (upstream2 + downstream2)*181.0)/3072.0;
        }
        case 120:
        {
            return ( (upstream0 + downstream0)*470 
                   - (upstream1 + downstream1)*245.0
                   + (upstream2 + downstream2)*31.0)/512.0;
        }
        case 124:
        {
            return ( (upstream0 + downstream0)*2854 
                   - (upstream1 + downstream1)*1509.0
                   + (upstream2 + downstream2)*191.0)/3072.0;
        }
        case 128:
        {
            return ( (upstream0 + downstream0)*722 
                   - (upstream1 + downstream1)*387.0
                   + (upstream2 + downstream2)*49.0)/768.0;
        }
        default:
        {
            std::ostringstream ss;   
            ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
            ss  <<"invalid scheme token\n";
            amrex::Abort(ss.str());
        }
    }
    // return ( (upstream0 + downstream0)*3 
    //        - (upstream1 + downstream1)*25.0
    //        + (upstream2 + downstream2)*150.0)/256.0;

    // return ( (upstream0 + downstream0)*23.0
    //        - (upstream1 + downstream1)*8.0
    //        + (upstream2 + downstream0))/32.0;
    return 0;
}