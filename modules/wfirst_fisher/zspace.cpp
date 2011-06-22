// ***********************************************************
// Code to estimate Fisher matrix constraints on b*sig8 and f*sig8, for 
// multiple populations. If used in publications, please cite
// White, Song & Percival (2008)
//
// Code calls GSL libraries (www.gnu.org/software/gsl/) available under the 
// GNU public license. It should compile using 
// 
// g++ -LGSL_DIR -lm -lgsl -lgslcblas -o fisher_gsl fisher_gsl.c
// 
// where GSL_DIR is the location of the libraries.
//
// Please report any bugs or problems to will.percival@port.ac.uk
//
// This version 8/10/2008, written by Will Percival
// ***********************************************************



#include "zspace.h"
#include <math.h>

using namespace std;
using namespace Eigen;

double tk_eh98(double);

void zspace_mbias_pk  ( 
                      const VectorXd & nbar,    // The number density in h^3 Mpc^-3
		      double sigma8,   // The real-space, linear clustering amplitude
                      const VectorXd & bias,    // The real-space, linear bias
		      double f,	       // f ~ Omega_m^(0.6) 
		      double Sigma_z,  // z error translated into comoving distance
		      double vol_mpc,  // The survey volume in h^-3 Mpc^3
		      double kmax,     // Maximum k-value to integrate to
                      MatrixXd &invfish  // inverse covariance matrix for bs and fs
    ) {

  // variables
  double k, mu;
  int is, i, j, s, l, m;

  // make sure things have the same size
  const int NSAMP = nbar.size();
  if (bias.size() != NSAMP) { fprintf(stderr, "bias and nbar should be compatible"); exit(0);}

  // catcher to avoid problems caused by odd inputs
  if(vol_mpc<=0.0) { fprintf(stderr,"volume<0"); exit(0); }
  for(is=0;is<NSAMP;is++) 
    if(bias[is]<=0.0 || nbar[is]<=0.0) 
      { fprintf(stderr,"bias or nbar <=0"); exit(0); }
  
  double kstep     = 0.001;                     // step in k integration
  double mustep    = 0.001;                     // step in mu intergration
  double Sigma_z2  = Sigma_z*Sigma_z;           // (redshift error)^2

  int NPOW = NSAMP*(NSAMP+1)/2;

  MatrixXd bigfish(NSAMP+1, NSAMP+1);
  bigfish.setZero(NSAMP+1, NSAMP+1);
  VectorXd Ps(NSAMP), err(NSAMP);
  MatrixXd cov(NPOW, NPOW), icov(NPOW, NPOW), dPdp(NPOW, NSAMP+1);

  // integral over kmax
  for(k=0.5*kstep; k<kmax; k+=kstep) {
    
    double tf   = tk_eh98(k);
    double pk = 7.03563e+06*k*tf*tf*sigma8*sigma8; 
    
    // integral over mu
    for(mu=0.5*mustep; mu<1.; mu+=mustep) {
      double mu2 = mu*mu;
      
      // has been written for simplicity, but could be rewritten for speed
      
      double zdamp  = exp(-k*k*Sigma_z2*mu2); // power spectrum damping due to z-error (inc photo-z)

      // calculate power spectra and multiplicative error terms
      for(i=0;i<NSAMP;i++) {
	Ps[i]      = (bias[i]+mu2*f)*(bias[i]+mu2*f)*pk*zdamp; // damped redshift-space P(k)
	err[i]     = (1.+1./(nbar[i]*Ps[i]));                  // multiplicative shot noise
      }

      // need covariance < P_ij P_lm >. Loop through these, assuming j>=i, m>=l
      int ip=-1;
      for(i=0;i<NSAMP;i++)
	for(j=i;j<NSAMP;j++) {
	  ip++;

	  int jp=-1;
	  for(l=0;l<NSAMP;l++) 
	    for(m=l;m<NSAMP;m++) {
	      jp++;


	      if(ip==jp) {
		// diagonal elements
                if(i==j) cov(ip, jp) = 2.*Ps[i]*Ps[j]*err[i]*err[j];
                if(i!=j) cov(ip, jp) = Ps[i]*Ps[j] + Ps[i]*Ps[j]*err[i]*err[j];
              } else {
		// off-diagonal elements
                cov(ip, jp) = 2.*sqrt(Ps[i]*Ps[j]*Ps[l]*Ps[m]);
                if(i==j && (i==l || i==m)) cov(ip, jp) *= err[i];
                if(l==m && (l==i || l==j)) cov(ip, jp) *= err[l];
                if(i!=j && l!=m && (i==l || i==m)) cov(ip, jp) = 0.5*cov(ip, jp)*(1.+err[i]);
                if(i!=j && l!=m && (j==l || j==m)) cov(ip, jp) = 0.5*cov(ip, jp)*(1.+err[j]);
              }

	    } // end first loop through power spectra P_lm

	  // set up derivatives of power spectra
	  for(s=0;s<NSAMP;s++) {
            //int index = ip*(NSAMP+1) + s;
            if(i==j && i==s) dPdp(ip, s) = Ps[i]*2./(sigma8*(bias[i]+mu2*f));
            if(i!=j && i==s) dPdp(ip, s) = Ps[j]   /(sigma8*(bias[j]+mu2*f));
            if(i!=j && j==s) dPdp(ip, s) = Ps[i]   /(sigma8*(bias[i]+mu2*f));
	  }
          dPdp(ip, NSAMP) = (Ps[i]*mu2/(sigma8*(bias[i]+mu2*f))+
					    Ps[j]*mu2/(sigma8*(bias[j]+mu2*f)));
 
	} // end second loop through power spectra P_ij
      
      // invert covariance matrix
      icov = cov.inverse();
      // now calculate Fisher matrix					
      bigfish += (dPdp.transpose() * icov * dPdp)*(k*k*vol_mpc);
      
    }
  }
  bigfish *= 2.0/4.0/M_PI/M_PI*kstep*mustep;
  invfish = bigfish.inverse();
}

// Transfer function of Eisenstein & Hu 1998 
// (Equation numbers refer to this paper)
double tk_eh98(double k)
{
  double rk,e,thet,thetsq,thetpf,b1,b2,zd,ze,rd,re,rke,s,rks,q,y,g;
  double ab,a1,a2,ac,bc,f,c1,c2,tc,bb,bn,ss,tb,tk_eh;
  double h,hsq,om_mhsq,om_b,om_m;

  // set up cosmology
  h    = 0.72;
  om_m = 0.25;
  om_b = 0.15*om_m;

  // convert k to Mpc^-1 rather than hMpc^-1
  rk=k*h;
  hsq=h*h;
  om_mhsq=om_m*hsq;

  // constants
  e=exp(1.);      
  thet=2.728/2.7;
  thetsq=thet*thet;
  thetpf=thetsq*thetsq;

  // Equation 4 - redshift of drag epoch
  b1=0.313*pow(om_mhsq,-0.419)*(1.+0.607*pow(om_mhsq,0.674));
  b2=0.238*pow(om_mhsq,0.223);
  zd=1291.*(1.+b1*pow(om_b*hsq,b2))*pow(om_mhsq,0.251)
    /(1.+0.659*pow(om_mhsq,0.828));

  // Equation 2 - redshift of matter-radiation equality
  ze=2.50e4*om_mhsq/thetpf;

  // value of R=(ratio of baryon-photon momentum density) at drag epoch
  rd=31500.*om_b*hsq/(thetpf*zd);

  // value of R=(ratio of baryon-photon momentum density) at epoch of
  // matter-radiation equality
  re=31500.*om_b*hsq/(thetpf*ze);

  // Equation 3 - scale of ptcle horizon at matter-radiation equality
  rke=7.46e-2*om_mhsq/(thetsq);

  // Equation 6 - sound horizon at drag epoch
  s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1.+sqrt(re)));

  // Equation 7 - silk damping scale
  rks=1.6*pow(om_b*hsq,0.52)*pow(om_mhsq,0.73)*(1.+pow(10.4*om_mhsq,-0.95));

  // Equation 10  - define q
  q=rk/13.41/rke;
      
  // Equations 11 - CDM transfer function fits
  a1=pow(46.9*om_mhsq,0.670)*(1.+pow(32.1*om_mhsq,-0.532));
  a2=pow(12.0*om_mhsq,0.424)*(1.+pow(45.0*om_mhsq,-0.582));
  ac=pow(a1,(-om_b/om_m))*pow(a2,pow(-(om_b/om_m),3.));

  // Equations 12 - CDM transfer function fits
  b1=0.944/(1.+pow(458.*om_mhsq,-0.708));
  b2=pow(0.395*om_mhsq,-0.0266);
  bc=1./(1.+b1*(pow(1.-om_b/om_m,b2)-1.));

  // Equation 18
  f=1./(1.+pow(rk*s/5.4,4.));

  // Equation 20
  c1=14.2 + 386./(1.+69.9*pow(q,1.08));
  c2=14.2/ac + 386./(1.+69.9*pow(q,1.08));

  // Equation 17 - CDM transfer function
  tc=f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q*q) +
    (1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q*q);

  // Equation 15
  y=(1.+ze)/(1.+zd);
  g=y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)));

  // Equation 14
  ab=g*2.07*rke*s/pow(1.+rd,0.75);

  // Equation 23
  bn=8.41*pow(om_mhsq,0.435);

  // Equation 22
  ss=s/pow(1.+pow(bn/rk/s,3.),1./3.);

  // Equation 24
  bb=0.5+(om_b/om_m) + (3.-2.*om_b/om_m)*sqrt(pow(17.2*om_mhsq,2.)+1.);

  // Equations 19 & 21
  tb=log(e+1.8*q)/(log(e+1.8*q)+c1*q*q)/(1+pow(rk*s/5.2,2.));
  tb=(tb+ab*exp(-pow(rk/rks,1.4))/(1.+pow(bb/rk/s,3.)))*sin(rk*ss)/rk/ss;
    
  // Equation 8
  tk_eh=(om_b/om_m)*tb+(1.-om_b/om_m)*tc;
  
  return tk_eh;
}
