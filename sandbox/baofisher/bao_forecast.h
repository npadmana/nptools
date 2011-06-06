#ifndef BAO_FISHER_H
#define BAO_FISHER_H

#ifdef __cpluscplus
extern "C" {
#endif


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
    );

#ifdef __cpluscplus
}
#endif
#endif // BAO_FISHER_H
