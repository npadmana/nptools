/* Simple implementation of the fitting formulae of Seo & Eisenstein (2006) */
/* This file can be compiled with gcc -O -o bao_forecast bao_forecast.c -lm
   to get a standalone code, or one can strip off main.c and use it as 
   a function call. */

/* In using the results, we remind the reader to include the covariance between
   the transverse and line of sight distance measurements and to treat
   the constraints here as those on D/s and H*s rather than D and H. */

#include <math.h>

/* Uncomment one of the following two lines */
#define WMAP_THREE
/* #define WMAP_ONE */

#define KSTEP 0.01
#define NUM_KSTEP 50

#ifdef WMAP_THREE    /* Omega_m = 0.24 */
float Pbao_list[NUM_KSTEP] = 
  { 14.10, 20.19, 16.17, 11.49, 8.853, 7.641, 6.631, 5.352, 4.146, 3.384, 
    3.028, 2.799, 2.479, 2.082, 1.749, 1.551, 1.446, 1.349, 1.214, 1.065, 
    0.9455, 0.8686, 0.8163, 0.7630, 0.6995, 0.6351, 0.5821, 0.5433, 0.5120, 0.4808, 
    0.4477, 0.4156, 0.3880, 0.3655, 0.3458, 0.3267, 0.3076, 0.2896, 0.2734, 0.2593, 
    0.2464, 0.2342, 0.2224, 0.2112, 0.2010, 0.1916, 0.1830, 0.1748, 0.1670, 0.1596};
    /* This is the power spectrum of WMAP-3, normalized to 1 at k=0.2 */
#define BAO_POWER 2710.0    /* The power spectrum at k=0.2h Mpc^-1 for sigma8=1 */
#define BAO_SILK 8.38
#define BAO_AMP 0.05169
#endif   /* WMAP_ONE */

#ifdef WMAP_ONE     /* Omega_m = 0.27 */
float Pbao_list[NUM_KSTEP] = 
  { 9.034, 14.52, 12.63, 9.481, 7.409, 6.397, 5.688, 4.804, 3.841, 3.108,
    2.707, 2.503, 2.300, 2.014, 1.707, 1.473, 1.338, 1.259, 1.174, 1.061,
    0.9409, 0.8435, 0.7792, 0.7351, 0.6915, 0.6398, 0.5851, 0.5376, 0.5018, 0.4741,
    0.4484, 0.4210, 0.3929, 0.3671, 0.3456, 0.3276, 0.3112, 0.2950, 0.2788, 0.2635,
    0.2499, 0.2379, 0.2270, 0.2165, 0.2062, 0.1965, 0.1876, 0.1794, 0.1718, 0.1646};
    /* This is the power spectrum of WMAP-1, normalized to 1 at k=0.2 */
#define BAO_POWER 2875.0    /* The power spectrum at k=0.2h Mpc^-1 for sigma8=1 */
#define BAO_SILK 7.76
#define BAO_AMP 0.04024
#endif   /* WMAP_ONE */


void bao_forecast ( 
	float number_density,	/* The number density in h^3 Mpc^-3 */
	float sigma8,		/* The real-space, linear clustering amplitude */
	float Sigma_perp,	/* The transverse rms Lagrangian displacement */
	float Sigma_par,	/* The line of sight rms Lagrangian displacement */
	float Sigma_z,		/* The line of sight rms comoving distance error due to redshift uncertainties */
		/* Note that Sigma_perp and Sigma_par are for pairwise differences,
		   while Sigma_z is for each individual object */
	float beta, 		/* The redshift distortion parameter */
	float volume,		/* The survey volume in h^-3 Gpc^3, set to 1 if input <=0 */
	float *Drms,		/* The rms error forecast for D/s, in percent */
	float *Hrms,		/* The rms error forecast for H*s, in percent */
	float *r,		/* The correlation coefficient between D and H */
				/* The covariance matrix for D/s and H*s is hence
					Drms**2    Drms*Hrms*r
					Drms*Hrms*r   Hrms**2     */
	float *Rrms		/* The rms error forecast for D/s and H*s, in 
				   percent, if one requires that the radial 
				   and transverse scale changes are the same. */
    )
{
    /* This routine takes about 300 microseconds to run with mustep=0.05 
     * on my Intel workstation. */
    float mustep = 0.05;
    int ik;
    float mu, mu2, k;
    float Fdd, Fdh, Fhh, sum;
    float Sigma_perp2, Sigma_par2, Sigma_z2, nP, redshift_distort, tmp, tmpz, Sigma2_tot;
    float Silk_list[NUM_KSTEP];

    Fdd = Fdh = Fhh = 0.0;
    if (volume<=0) volume=1.0;   /* Note: even setting 0 volume to 1, to avoid divide by zero below */
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
	Fdd += sum*(1-mu2)*(1-mu2);
	Fdh += sum*(1-mu2)*mu2;
	Fhh += sum*mu2*mu2;
    }
    *r = Fdh/sqrt(Fhh*Fdd);   /* This is the correlation coeff between D and H */
    	/* Recall that the Fisher matrix parameter is actually H^-1, hence
	   there is one extra minus sign to cancel that of the F to C inversion */
    Fdd *= BAO_AMP*BAO_AMP/8.0/M_PI/M_PI*1.0e9*KSTEP*mustep*volume;
    Fhh *= BAO_AMP*BAO_AMP/8.0/M_PI/M_PI*1.0e9*KSTEP*mustep*volume;
    /* Invert the Fisher matrix and quote the diagonal elements */
    *Drms = 1.0/sqrt(Fdd*(1.0-(*r)*(*r)));
    *Hrms = 1.0/sqrt(Fhh*(1.0-(*r)*(*r)));
    *Rrms = (*Drms)*sqrt((1-(*r)*(*r))/(1+(*Drms)/(*Hrms)*(2*(*r)+(*Drms)/(*Hrms))));
    return;
}


/* ========== Remove everything below this line if one wants ============ */
/* ============ to call bao_forecast() from other programs   ============ */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    float D,H,r,R;
    int j;
    if (argc!=7) {
        fprintf(stderr, "Call with n, sigma8, Sigma_perp, Sigma_par, Sigma_z, and beta.\n");
        fprintf(stderr, "Output is sigma(D/s), sigma(H*s), r, and sigma(spherical)\n    for a 1 h^-3 Gpc^3 survey.\n");
        exit(1);
    }
    /* Uncomment the following line to invoke 1e5 iterations for timing */
    /* for (j=0;j<1e5;j++)  */
    bao_forecast( 
       atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]),
       1.0, &D, &H, &r, &R);
    printf("%f %f %f %f\n", D, H, r, R);
    return 0;
}
