#include <iostream>
#include <cstdio>
#include <string>
#include <complex>

#include <boost/format.hpp>
#include <drfftw.h>

#include "np_functional.h"
#include "multi_array_wrap.h"
#include "gslwrap.h"
#include "eigen_utils.h"

using namespace std;
using namespace boost;
using namespace ma;


const int Ndim = 3;		// particles live in 3 dimensions
const double pi = atan(1.0)*4.0;
int Ngrid = 128; // Keep these non-constant so that we might change them.
double Lbox = 2000.0;

// Useful typedefs
typedef complex<double> cdouble;
typedef multi_array_ref<complex<double>, 4> ref4c;
typedef MA<array3d>::IndexArr indices3d;

// Convenience functions
double _ki_impl(i_ ndx) { return (2.0*pi/Lbox)*(ndx > Ngrid/2  ? ndx - Ngrid : ndx);}
boost::function<double(i_)> ki = _ki_impl;

template <class A>
void fftw_forward(MA<A> grid) {
    const int Ng = grid.shape()[1];
    const int narr = grid.shape()[0];

    // Set up the plan
    rfftwnd_plan plan = rfftw3d_create_plan(Ng, Ng, Ng, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
    if (plan == NULL) throw runtime_error("Error creating the FFTW plan\n");

    // Do it
    double *pos;
    for (int ii = 0; ii < narr; ++ii) {
        pos = grid.sub(ii).ref();
        rfftwnd_one_real_to_complex(plan, pos, reinterpret_cast< fftw_complex* >(pos));
    }

    // Normalize
    grid /= pow(static_cast<double>(Ng), 3);

    // Destroy
    rfftwnd_destroy_plan(plan);
}

template <class A>
void fftw_reverse(MA<A> grid) {
    const int Ng = grid.shape()[1];
    const int narr = grid.shape()[0];

    // Set up the plan
    rfftwnd_plan plan = rfftw3d_create_plan(Ng, Ng, Ng, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
    if (plan == NULL) throw runtime_error("Error creating the FFTW plan\n");

    // Do it
    double *pos;
    for (int ii = 0; ii < narr; ++ii) {
        pos = grid.sub(ii).ref();
        rfftwnd_one_complex_to_real(plan, reinterpret_cast< fftw_complex* > (pos), pos);
    }

    // Destroy
    rfftwnd_destroy_plan(plan);
}

double index2k(const indices3d& ilist) {
    Vector3d tmp;
    transform(ilist.begin(), ilist.end(), tmp.data(), ki);
    return tmp.norm();
}

// Functor to compute k_{i} v_{i}
struct kdot_impl {
    int idim;
    void operator()(cdouble& val, const indices3d& ilist) {
        val *= ki(ilist[idim]);
    }
};

// Functor to compute x* y and store it in z
struct complex_mult_impl {
    void operator()(cdouble& x, cdouble& y, cdouble& z) {
        z = conj(x) * y;
    }
};

class PkStruct : public GSLHist {
    public :
       void operator()(const cdouble& val, const indices3d& ilist);
       MatrixXd pk(double norm=1.0);
};

void PkStruct::operator()(const cdouble& val, const indices3d& ilist) {
    double kk = index2k(ilist);
    add(log(kk), real(val));
}

MatrixXd PkStruct::pk(double norm) {
    // Get the histogram outputs
    MatrixXd retval = val();
    int nrows = retval.rows();

    // I hate log outputs, convert to a human readable format
    retval.col(0) = retval.col(0).array().exp();
    retval.col(1) = retval.col(1).array().exp();

    //Normalize
    retval.col(2) *= pow(Lbox, 3) * norm;

    // Divide by n, taking care not to divide by zero
    boost::function< double(double, double) > nicedivide = if_then_else_return(_2 > 0, _1/_2, 0.0);
    transform(retval.col(2).data(), retval.col(2).data() + nrows,
              retval.col(3).data(), retval.col(2).data(),
              nicedivide);

    return retval;
}






// Define the position and velocity data
class Coords {
    public :
        MA<array2f> pos, vel;
        int npart;
        Coords(string fn);
};

Coords::Coords(string fn) {
    // Header structure
    struct FileHeader
    {
       int   npart;		// Total number of particles
       int   nsph;		// Number of gas particles
       int   nstar;		// Number of star particles
       float aa;			// Scale factor
       float softlen;		// Gravitational softening
    } header;

    // Declarations yanked from Martin's code
    int   ngot;
    int    eflag,hsize;
    FILE   *fpr = NULL;


    // Can we open the file?
    if ((fpr = fopen(fn.c_str(),"r")) == NULL)  throw "Failed to open file";

    // Read endian flag.
    ngot = fread (&eflag,sizeof(int),1,fpr);
    if (ngot != 1) throw "Can't read file";
    if (eflag != 1) throw "Endian flag is not 1";

    // Read header size.
    ngot = fread(&hsize,sizeof(int),1,fpr);
    if (ngot != 1) throw "Error reading file header";
    if (hsize != sizeof(struct FileHeader)) throw "Incorrect header size";

    // Read and unpack header.
    ngot = fread(&header,sizeof(struct FileHeader),1,fpr);
    if (ngot != 1) throw "Unable to read header";
    npart = header.npart;

    pos.resize(extents[npart][3]);
    vel.resize(extents[npart][3]);

    // The usual code mucks around here, but just suck it all in,
    // since we had bloody well have enough memory for these data.
    int nget = Ndim * npart;
    ngot = fread(pos.ref(), sizeof(float), nget, fpr);
    if (ngot != nget) {
        cout << boost::format("Only got %1% of %2% particle positions....\n") % ngot % npart;
        throw "Error reading....";
    }
    ngot = fread(vel.ref(), sizeof(float), nget, fpr);
    if (ngot != nget) {
        cout << boost::format("Only got %1% of %2% particle velocities....\n") % ngot % npart;
        throw "Error reading....";
    }

    fclose(fpr);
}



/* CIC assign to a grid.
 * If idim<0, then just count particles, otherwise add in
 * velocities
 */
template <class Arr>
void cic (MA<Arr> grid, Coords& posvel, int idim=-1) {

   int Ng = grid.shape()[0]; // Assume cube
   BOOST_STATIC_ASSERT(Ndim == 3);
   Vector3i im, ip;
   Vector3f dx, x1;

   // It will be useful to consider the positions as an Eigen array
   Map< MatrixXf > epos = posvel.pos.matrix();

   double val;
   for (int ii = 0; ii < posvel.npart; ++ii) {
       x1 = Ng*epos.col(ii);
       im = x1.cast<int>(); // Slight ugliness due to casting
       ip = im.array() + 1; ip = ip.unaryExpr(_1%Ng); // % not defined in Eigen
       dx = x1 - im.cast<float>(); //slight ugliness due to casting
       if (idim >= 0) {val = posvel.vel[ii][idim];} else {val = 1;}
       grid[im[0]][im[1]][im[2]] += (1.0-dx[0])*(1.0-dx[1])*(1.0-dx[2])*val;
       grid[ip[0]][im[1]][im[2]] +=      dx[0] *(1.0-dx[1])*(1.0-dx[2])*val;
       grid[im[0]][ip[1]][im[2]] += (1.0-dx[0])*     dx[1] *(1.0-dx[2])*val;
       grid[im[0]][im[1]][ip[2]] += (1.0-dx[0])*(1.0-dx[1])*     dx[2] *val;
       grid[ip[0]][ip[1]][im[2]] +=      dx[0] *     dx[1] *(1.0-dx[2])*val;
       grid[ip[0]][im[1]][ip[2]] +=      dx[0] *(1.0-dx[1])*     dx[2] *val;
       grid[im[0]][ip[1]][ip[2]] += (1.0-dx[0])*     dx[1] *     dx[2] *val;
       grid[ip[0]][ip[1]][ip[2]] +=      dx[0] *     dx[1] *     dx[2] *val;
   }

}

int main() {

    // read in the file
    Coords posvel("dm_1.0000.bin");
    cout << "Read in " << posvel.npart << " particles\n";

    // Exercise the interface and code
    {
           Vector3f av;
           Map< MatrixXf > epos = posvel.pos.matrix();
           av = epos.rowwise().sum()/ posvel.npart;
           cout << "Average position : " << av.transpose() << endl;

           Map< MatrixXf > evel = posvel.vel.matrix();
           av = evel.rowwise().sum()/ posvel.npart;
           cout << "Average velocity : " << av.transpose() << endl;
    }


    // Define the grid
    MA<array4d> grid(extents[4][Ngrid][Ngrid][Ngrid+2]); // Pad for FFT, 1 component for overdensity, 3 for velocity
    // Define views on this grid
    MA<view4_4d> rgrid = grid.view<view4_4d>(indices[r_()][r_()][r_()][r_(0, Ngrid)]);
    MA<ref4c> cgrid ( new ref4c(reinterpret_cast< complex<double>* >(grid.ref()), extents[4][Ngrid][Ngrid][Ngrid/2 + 1]));
    MA<array4d::reference> dense = rgrid.sub(0);


    // CIC grid
    cout << "Integrated density, density-weighted velocity -->\n";
    for (int ii = 0; ii < 4; ++ii) {
        cic(rgrid.sub(ii), posvel, ii-1);
        cout << boost::format("%1$6i %2$20.10f\n") % ii % (rgrid.sub(ii).sum());
    }



    // Normalize by density
    int nzeros;
    cout << "Integrated velocity field -->\n";
    for (int ii=1; ii < 4; ++ii) {
        nzeros = 0;
        boost::function < void(double&, double&) > ff = if_(_2 > 0)[_1 /= _2].else_[var(nzeros)++];
        multi_for(rgrid.sub(ii), dense,ff);
        cout << boost::format("%1$6i %2$20.10f %3$10i\n") % ii % (rgrid.sub(ii).sum()*Lbox) % nzeros;
       }
    // Normalize the density by rho_mean
    cout << dense.sum() << endl;
    dense /= static_cast<double>(posvel.npart)/pow( static_cast<double>(Ngrid), 3);
    cout << dense.sum() << endl;


    // FFT
    fftw_forward(grid);


    // Compute k_i v_i for each grid separately
    kdot_impl kdot;
    for (int ii=1; ii < 4; ++ii) {
        kdot.idim = ii-1;
        multi_for_native_indices(cgrid.sub(ii), kdot);
    }

    // Now add the vectors together --- things are contiguous, so go ahead and use the
    // easy route
    cgrid.sub(1)() += cgrid.sub(2)() + cgrid.sub(3)();
    // This should be what we need modulo normalization factors

    // Compute P(k)
    PkStruct pk;
    double lkmin = log(0.008); double lkmax = log(0.3);
    pk.setbins(lkmin, lkmax, 20);

    // We want to compute delta* delta, delta* theta, and theta* theta
    complex_mult_impl cmult;
    // First do delta* theta -- and store this in cgrid.sub(2)
    multi_for(cgrid.sub(0), cgrid.sub(1), cgrid.sub(2), cmult);
    // delta* delta and theta* theta -- store in place
    multi_for(cgrid.sub(0), cgrid.sub(0), cgrid.sub(0), cmult); // We do it this way to make the connection clear with delta* theta
    multi_for(cgrid.sub(1), cgrid.sub(1), cgrid.sub(1), cmult); // We do it this way to make the connection clear with delta* theta


    // Do the matter power spectrum
    multi_for_native_indices(cgrid.sub(0), pk);
    cout << pk.pk() << endl;
    writeMatrix("pk_matter.dat", "%1$10.4e ", pk.pk());

    // Do the velocity power spectrum
    pk.reset();
    multi_for_native_indices(cgrid.sub(1), pk);
    double fac; fac = pow(Lbox, 2);
    cout << pk.pk(fac) << endl; // The true adds in the correct velocity normalization
    writeMatrix("pk_theta.dat", "%1$10.4e ", pk.pk(fac));


    // Do the cross power spectrum
    pk.reset();
    multi_for_native_indices(cgrid.sub(2), pk);
    fac = Lbox;
    cout << pk.pk(fac) << endl; // The true adds in the correct velocity normalization
    writeMatrix("pk_delta_theta.dat", "%1$10.4e ", pk.pk(fac));

}
