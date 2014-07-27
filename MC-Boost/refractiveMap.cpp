/*
 * refractiveMap.cpp
 *
 *  Created on: Apr 23, 2012
 *      Author: StaleyJW
 */


#include "refractiveMap.h"
#include "vector3D.h"


#include <boost/lexical_cast.hpp>





RefractiveMap::RefractiveMap(const int Nx, const int Ny, const int Nz)
: X_PML_OFFSET(25),
  Y_PML_OFFSET(10),
  Z_PML_OFFSET(10)
{
    total_medium_size.X = Nx;
    total_medium_size.Y = Ny;
    total_medium_size.Z = Nz;
    
    refractive_total = refractive_x = refractive_y = refractive_z = NULL;
}


RefractiveMap::~RefractiveMap()
{

    /// Clean up memory.
    if (refractive_total)
    {
        delete refractive_total;
        refractive_total = NULL;
    }
    
    if (refractive_x)
    {
        delete refractive_x;
        refractive_x = NULL;
    }
    
    if (refractive_y)
    {
        delete refractive_y;
        refractive_y = NULL;
    }
    
    if (refractive_z)
    {
        delete refractive_z;
        refractive_z = NULL;
    }
}






float
RefractiveMap::getRefractiveIndexFromGradientGrid(const char axis, const int x_photon, const int y_photon, const int z_photon)
{
    /// Refractive index changes induced from pressure propagating along x-axis.
    if (axis == 'x')
    {
        assert(refractive_x != NULL);
        return refractive_x->GetElementFrom3D(x_photon,
                                              y_photon,
                                              z_photon);
    }

    /// Refractive index changes induced from pressure propagating along y-axis.
    if (axis == 'y')
    {
        assert(refractive_y != NULL);
        return refractive_y->GetElementFrom3D(x_photon,
                                              y_photon,
                                              z_photon);
    }
    
    /// Refractive index changes induced from pressure propagating along z-axis.    
    if (axis == 'z')
    {
        assert(refractive_z != NULL);
        return refractive_z->GetElementFrom3D(x_photon,
                                              y_photon,
                                              z_photon);
    }
    
    /// Error
    return -1.0;
}



// Returns the spatially located refractive index value from the supplied location in the grid.
/// x_photon, y_photon and z_photon are the voxel coordinates obtained by translating the spatial coordinate
/// of the current scattering event.
float
RefractiveMap::getRefractiveIndexFromGrid(const int x_photon, const int y_photon, const int z_photon)
{
    /// NEW: Used to access k-Wave data structure.
//    return Get_refractive_index_TRealMatrix(x_photon + X_PML_OFFSET,
//                                            y_photon + Y_PML_OFFSET,
//                                            z_photon + Z_PML_OFFSET);
    
    return Get_refractive_index_TRealMatrix(x_photon,
                                            y_photon,
                                            z_photon);

}


/// Return the refractive index from the TRealMatrix.
float
RefractiveMap::Get_refractive_index_TRealMatrix(const int x_photon, const int y_photon, const int z_photon)
{
    assert(refractive_total != NULL);
    
    /// NOTE: 'GetElementFrom3D' is a method of class 'TRealMatrix' provided by kWave.
    ///       To have the ultrasound z-axis orthogonal to the light z-axis we need to take
    ///       into account the way in which kWave stores data.  3-D grids are accessed as,
    ///       GetElementFrom3D(x-axis, y-axis, z-axis).  To accomodate this, Monte-carlo
    ///       x-coordinate is used to obtain kWave's z-coordinate data.  This means the
    ///       x-y plane in Monte-carlo forms the y-z plane in kWave.  Simply put, a point
    ///       in x-y from Monte-carlo retrieves data from z in kWave.
    ///       Below is Monte-carlo coordinates used to access kWave data so that they propagate
    ///       orthogonally.
    
    /// ON FURTHER INSPECTION I BELIEVE THE ABOVE IS NOT TRUE.  I THINK IT IS ALREADY MADE BY LOOKING
    /// AT THE TRANSDUCER PLOT IN KWAVE.  VERIFY!
    ///
	///              !!!!   Verified by Jiri (k-Wave developer) !!!!
    return refractive_total->GetElementFrom3D(x_photon, y_photon, z_photon);

}


/// Invert the phase of the refractive index data 180 degrees by multiplying through the matrix by -1.
void
RefractiveMap::Invert_phase(void)
{
    float * raw_data        = refractive_total->GetRawData();
    const size_t sensor_size  = refractive_total->GetTotalElementCount();

    for (size_t i = 0; i < sensor_size; i++)
    {
        /// Perform the inversion.
        raw_data[i] *= -1.0f;
    }
}



void
RefractiveMap::Update_refractive_map_from_sensor(TRealMatrix * refractive_total_sensor, const long * sensor_index)
{
//    if (refractive_total == NULL) refractive_total = new TRealMatrix(total_medium_size);
//    
//    const size_t sensor_size = refractive_total_sensor->GetTotalElementCount();
    
    /// Loop over the entire 'refractive_total' matrix and only update over the region where the sensors are stared in the medium.
    /// FIXME:
    /// - Assignment is protected in the 'TRealMatrix' class. Maybe 'RefractiveMap' should inherent from 'TRealMatrix'???
//    for (size_t i = 0; i < sensor_size; i++)
//    {
//        refractive_total[sensor_index[i]] = refractive_total_sensor[i];
//    }
    
    
}








