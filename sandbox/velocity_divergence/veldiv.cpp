#include <iostream>
#include <cstdio>
#include <string>
#include <complex>

#include <boost/format.hpp>
#include <drfftw.h>

#include "np_functional.h"
#include "multi_array_wrap.h"

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
    grid /= (1./pow(static_cast<double>(Ng), 3));

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
   MatrixXf epos = posvel.pos.eig2<MatrixXf>();

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
           MatrixXf epos = posvel.pos.eig2<MatrixXf>();
           av = epos.rowwise().sum()/ posvel.npart;
           cout << "Average position : " << av.transpose() << endl;

           MatrixXf evel = posvel.vel.eig2<MatrixXf>();
           av = evel.rowwise().sum()/ posvel.npart;
           cout << "Average velocity : " << av.transpose() << endl;
    }


    // Now define the grid
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
    dense /= static_cast<double>(posvel.npart)/pow( static_cast<double>(Ngrid), 3);


    // FFT
    fftw_forward(grid);


}
