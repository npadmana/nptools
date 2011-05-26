#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include "multi_array_wrap.h"

using namespace std;
using namespace boost;
using namespace ma;
using namespace boost::lambda;

const int Ndim = 3;		// particles live in 3 dimensions


// Define the position and velocity data
class Coords {
    public :
        array2f pos, vel;
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
    ngot = fread(pos.data(), sizeof(float), nget, fpr);
    if (ngot != nget) {
        cout << boost::format("Only got %1% of %2% particle positions....\n") % ngot % npart;
        throw "Error reading....";
    }
    ngot = fread(vel.data(), sizeof(float), nget, fpr);
    if (ngot != nget) {
        cout << boost::format("Only got %1% of %2% particle velocities....\n") % ngot % npart;
        throw "Error reading....";
    }

    fclose(fpr);
}



/* CIC assign to a grid.
 * If idim=1, then just count particles, otherwise add in
 * velocities
 */
template <class MultiArray>
void cic(MultiArray& grid, const Coords& posvel) {
    vector<int> N(Ndim), ix(Ndim), ix1(Ndim);

    N.assign(grid.shape(), grid.shape()+3);
    cout << "The dimensions of the grid are : ";
    for_each(N.begin(), N.end(), cout << _1 << " ");
    cout << endl;

    array2f::const_reference part1 = posvel.pos[0];
    transform(part1.begin(), part1.end(), N.begin(), ix.begin(), _1 * _2); // Get lower grid position
    transform(ix.begin(), ix.end(), N.begin(), ix1.begin(), (_1+1)%(_2)); // Get upper grid position
    for (int ii =0; ii < 3; ++ii) {
        cout << format("%1% %2% %3% %4% \n") % part1[ii] % N[ii] % ix[ii] % ix1[ii];
    }

}




int main() {

    // Ngrid
    const int Ngrid = 128;

    // read in the file
    Coords posvel("dm_1.0000.bin");
    cout << "Read in " << posvel.npart << " particles\n";

    // define the grid -- and appropriate views of it
    array4d grid_full(extents[4][Ngrid][Ngrid][Ngrid+2]);

    // Remove FFT buffer
    view4_4d rgrid_full = grid_full[indices[r_()][r_()][r_()][r_(0, Ngrid)]];
    multi_for(rgrid_full, _1=0.0);

    array4d::reference ri = rgrid_full[0];
    cic(ri, posvel);


}
