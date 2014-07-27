#ifndef COMMONSTRUCTS_H
#define COMMONSTRUCTS_H


#include <MC-Boost/pressureMap.h>
#include <MC-Boost/refractiveMap.h>
#include <MC-Boost/displacementMap.h>

#include "multikey.h"


/// ---------------------------------------------------------------------
//typedef struct {
//    // Create a pointer to a PressureMap object.  For use with
//	// modeling acousto-optics.
//	PressureMap * pmap;
//    
//	// Create a pointer to a RefractiveMap object.  For use with
//	// modeling acousto-optics.
//	RefractiveMap * nmap;
//    
//    // Create a pointer to a Displacement object.  For use with
//    // modeling acousto-optics.
//    DisplacementMap * dmap;
//    
//    // Frequency of the ultrasound used.
//    float  US_freq;
//    
//    // The wavenumber of the ultrasound.
//    double  waveNumber;
//    
//    /// The density of the medium.
//    ///double  density;
//    
//    /// The speed-of-sound of the medium
//    //double  speed_of_sound;
//    
//    /// The total number of elements in the sensor mask.
//    long    sensor_mask_index_size;
//    
//    // Number of time steps in the simulation.
//    int totalTimeSteps;
//    
//    /// Time step (i.e. dt used in k-Wave for updating calculations).
//    float dt;
//    
//    
//} kWaveSim;



///// ---------------------------------------------------------------------
//// Structure for containing cartesian coordinates.
//struct _coords {
//	double x;
//	double y;
//	double z;
//};
//
//
//
//// Structure for containing spherical coordinates.
//struct _sphereCoords{
//    double r;
//    double theta;
//    double phi;
//};
//
//
//
//// Structure for containing the direction cosines of the photon.
//struct _directionCos{
//    double x;
//    double y;
//    double z;
//};
//
//
//
///// ---------------------------------------------------------------------
//struct VoxelAttributes {
//    float dx;
//    float dy;
//    float dz;
//    
//    size_t Nx;
//    size_t Ny;
//    size_t Nz;
//};
//
//
//
///// ---------------------------------------------------------------------
//struct _Aperture_Properties{
//    double radius;
//    
//    struct _coords coordinates;
//    
//	bool xy_plane;
//	bool xz_plane;
//	bool yz_plane;
//    
//};



/// ---------------------------------------------------------------------
typedef struct {
    bool DISPLACE;
    bool REFRACTIVE_TOTAL;
    bool REFRACTIVE_GRADIENT;
    bool MODULATION_DEPTH;
    bool SAVE_SEEDS;
    bool USE_SEEDS;
    bool STORE_FLUENCE;
} MC_Parameters;


/// ---------------------------------------------------------------------
typedef struct {
    double refractive_index_contribution; /// \delta{n}
    double displacement_contribution;     /// \delta{d}
    double combined_contribution;         /// \delta{n} + \delta{d}
} opticalPathLengths;


typedef struct {
    double refractive_index_contribution; /// \delta{n}
    double displacement_contribution;     /// \delta{d}
    double combined_contribution;         /// \delta{n} + \delta{d}
    double weight;
    double exit_x_coord;
    double exit_y_coord;
    double exit_z_coord;
} exitInfo;

/// Multikey map for storing the optical path lengths (OPL's)
/// for comparison of modulation depth.
typedef std::map<MultiKey, std::vector<opticalPathLengths> > MultiKeyMap_OPLs;

/// Multikey map for storing the exit information of the photon
/// when it had left the medium and was detected.
typedef std::map<MultiKey, std::vector<exitInfo> > MultiKeyMap_ExitInfo;






#endif   /// COMMONSTRUCTS_H