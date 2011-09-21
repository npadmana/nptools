#include "bao_forecast.h"

using namespace std;
using namespace Eigen;

// Define the fiducial cosmology
detf _fidcosmo = detf_fiducial();
const double sigma0 = 12.4 * 0.8/0.9;


double Sigma_perp(double a) {
    return sigma0 * growth(a, _fidcosmo) ;
}

double Sigma_par(double a) {
    return Sigma_perp(a)*(1+fgrowth(a, _fidcosmo));
}

#define WMAP_THREE
/* #define WMAP_ONE */

#define KSTEP 0.01
#define NUM_KSTEP 50

double Pbao_list[NUM_KSTEP] =
  { 14.10, 20.19, 16.17, 11.49, 8.853, 7.641, 6.631, 5.352, 4.146, 3.384,
    3.028, 2.799, 2.479, 2.082, 1.749, 1.551, 1.446, 1.349, 1.214, 1.065,
    0.9455, 0.8686, 0.8163, 0.7630, 0.6995, 0.6351, 0.5821, 0.5433, 0.5120, 0.4808,
    0.4477, 0.4156, 0.3880, 0.3655, 0.3458, 0.3267, 0.3076, 0.2896, 0.2734, 0.2593,
    0.2464, 0.2342, 0.2224, 0.2112, 0.2010, 0.1916, 0.1830, 0.1748, 0.1670, 0.1596};
    /* This is the power spectrum of WMAP-3, normalized to 1 at k=0.2 */
#define BAO_POWER 2710.0    /* The power spectrum at k=0.2h Mpc^-1 for sigma8=1 */
#define BAO_SILK 8.38
#define BAO_AMP 0.05169


/* Return the covariance matrix between D/s, H*s */
Matrix2d bao_forecast (
        const double number_density,   /* The number density in h^3 Mpc^-3 */
        const double sigma8,           /* The real-space, linear clustering amplitude */
        const double Sigma_perp,       /* The transverse rms Lagrangian displacement */
        const double Sigma_par,        /* The line of sight rms Lagrangian displacement */
        const double Sigma_z,          /* The line of sight rms comoving distance error due to redshift uncertainties */
                /* Note that Sigma_perp and Sigma_par are for pairwise differences,
                   while Sigma_z is for each individual object */
        const double beta,             /* The redshift distortion parameter */
        const double volume           /* The survey volume in h^-3 Gpc^3 */
    )
{
    /* This routine takes about 300 microseconds to run with mustep=0.05
     * on my Intel workstation. */
    double mustep = 0.05;
    int ik;
    double mu, mu2, k;
    double Sigma_perp2, Sigma_par2, Sigma_z2, nP, redshift_distort, tmp, tmpz, Sigma2_tot, sum;
    Matrix2d fish, cov;
    fish.setZero();

    double Silk_list[NUM_KSTEP];

    Sigma_perp2 = Sigma_perp*Sigma_perp;  /* We only use the squares of these */
    Sigma_par2 = Sigma_par*Sigma_par;
    Sigma_z2 = Sigma_z*Sigma_z;
    nP = number_density*sigma8*sigma8*BAO_POWER;   /* At k=0.2 h Mpc^-1 */
    if (sqrt(Sigma_par2+Sigma_z2)>3*Sigma_perp) mustep /= 10.0;
        /* Take finer steps if integrand is anisotropic. */
        /* One might need to adjust this further */

    for (ik=0, k=0.5*KSTEP; ik<NUM_KSTEP; ik++,k+=KSTEP)
        Silk_list[ik] = exp(-2.0*pow(k*BAO_SILK,1.40))*k*k;
        /* Pre-compute this for speed.  However, if you need an extra 10% speed bump,
        move this outside of the function and reuse the computation in each new call. */

    for (mu=0.5*mustep; mu<1; mu+=mustep) {
        mu2 = mu*mu;
        redshift_distort = (1+beta*mu2)*(1+beta*mu2);
        tmp = 1.0/(nP*redshift_distort);
        Sigma2_tot = Sigma_perp2*(1-mu2)+Sigma_par2*mu2;

        for (sum=0.0, ik=0, k=0.5*KSTEP; ik<NUM_KSTEP; ik++,k+=KSTEP) {
            tmpz = Pbao_list[ik]+tmp*exp(k*k*Sigma_z2*mu2);
            sum += Silk_list[ik]*exp(-k*k*Sigma2_tot)/tmpz/tmpz;
            /* These two exp() take nearly all of the run time */
        }
        fish(0,0) += sum*(1-mu2)*(1-mu2);
        fish(0,1) += sum*(1-mu2)*mu2;
        fish(1,1) += sum*mu2*mu2;
    }
    // Flip the sign, since the fisher matrix parameter is actually H^-1
    fish(0,1) = -fish(0,1);
    // Symmetrize the matrix
    fish(1,0) = fish(0,1);

    // Normalization quantities
    fish *= (BAO_AMP*BAO_AMP/8.0/M_PI/M_PI*1.0e9*KSTEP*mustep*volume);

    // Compute the covariance matrix
    cov = fish.inverse();

    return cov;
}

// Some convenience routines
Matrix2d bao_forecast_shell (
        const double number_density,   /* The number density in h^3 Mpc^-3 */
        const double sigma8,           /* The real-space, linear clustering amplitude */
        const double Sigma_z,          /* The line of sight rms comoving distance error due to redshift uncertainties */
                /* Note that Sigma_perp and Sigma_par are for pairwise differences,
                   while Sigma_z is for each individual object */
        const double beta,             /* The redshift distortion parameter */
        const double zmin, const double zmax, const double area, /* zmin, zmax, area in deg^2 */
        const double recon,
        const bool isnumber /* False (default) if number_density is a true density, otherwise True if number density is number/dz/area */
    )
{
    double vol, nbar; // Volume of the shell in Gpc^3 h^-3
    vol = shellVol_Gpc_h(z2a(zmin), z2a(zmax), area, _fidcosmo);
    if (isnumber) {
       nbar = number_density * (zmax-zmin) * area/vol * 1.e9;
    } else {
       nbar = number_density;
    }
    double zmid = (zmin + zmax)/2.0;
    double amid = z2a(zmid);

    double sperp = Sigma_perp(amid)*recon;
    double spar = Sigma_par(amid)*recon;

    return bao_forecast(number_density, sigma8, sperp, spar, Sigma_z, beta, vol);
}

