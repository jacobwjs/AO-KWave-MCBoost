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
  Z_PML_OFFSET(10),
  refractive_total(NULL),
  background_refractive_total(NULL),
  refractive_x(NULL),
  refractive_y(NULL),
  refractive_z(NULL)
{
    total_medium_size.X = Nx;
    total_medium_size.Y = Ny;
    total_medium_size.Z = Nz;
    
}


RefractiveMap::~RefractiveMap()
{

    /// Clean up memory.
    if (refractive_total)
    {
        delete refractive_total;
        refractive_total = NULL;
    }
    
    if (background_refractive_total)
    {
        delete background_refractive_total;
        background_refractive_total = NULL;
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


/// x_photon, y_photon and z_photon are the voxel coordinates obtained by translating the spatial coordinate
/// of the current scattering event.
float
RefractiveMap::getBackgroundRefractiveIndexFromGrid(const int x_photon, const int y_photon, const int z_photon)
{
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
    return background_refractive_total->GetElementFrom3D(x_photon, y_photon, z_photon);
    
}


/// Invert the phase of the refractive index data 180 degrees by calculating the
/// modulated value, inverting its sign (+ or -) and placing the new value back.
void
RefractiveMap::Invert_phase(TLongMatrix * sensor_mask_index)
{
    float * nmap = refractive_total->GetRawData();
    const float * nmap_background = background_refractive_total->GetRawData();
    
    const size_t sensor_size = sensor_mask_index->GetTotalElementCount();
    const long *  index = sensor_mask_index->GetRawData();
    
    float modulated_value = 0.0f;
    
    /// Only invert the medium where the sensor locations exist.
    for (size_t i = 0; i < sensor_size; i++)
    {
        /// The change in refractive index due to ultrasound.
        modulated_value = nmap[index[i]] - nmap_background[index[i]];
        
        
        /// Perform the inversion.
        /// If the modulated value was initially positive, then we flip its sign (-) and add it to the background medium.
        /// If initially negative the opposite occurs. Essentially the same as sending in another pulse of ultrasound with
        /// an initial phase shift of 180 degrees.
        nmap[index[i]] = nmap_background[index[i]] + (-1*modulated_value);
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








