// Class defines the properties of a photon.
#ifndef PHOTON_H
#define PHOTON_H

#include <MatrixClasses/RealMatrix.h>

#include <boost/math/constants/constants.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include <MC-Boost/RNG.h>
#include <MC-Boost/coordinates.h>
#include <MC-Boost/vectorMath.h>
#include "common_structs.h"

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;
//#include <boost/random/uniform_real.hpp>
//#include <boost/random/variate_generator.hpp>
//#include <boost/random/mersenne_twister.hpp>

#define ALIVE 1		// Value depicting Photon should continue propagation.
#define DEAD  0	    // Photon has lost all energy and failed roulette.
#define ONE_MINUS_COSZERO 1.0e-12
#define THRESHOLD	0.01		// Threshold for determining if we should perform roulette
#define CHANCE      0.1  		// Used in roulette
#define SIGN(x)           ((x)>=0 ? 1:-1)



// Forward decleration of objects.
class Medium;
class Vector3d;
class Layer;
class InjectionAperture;
class Aperture;



//typedef struct coords InjectionCoords;

class Photon
{
public:

	friend class LoggerBase;

	// Constructors
	Photon(void);
	
	// Destructor
	~Photon(void);
    
    // Common function to initialize basic values of the photon object.
    void    initCommon(void);
	
	// Set the number of iterations this Photon (i.e. thread) will run.
	void	setIterations(const int n);

	// Move photon to new position
	void	hop(void);

	// Drop absorbed energy from photon weight due to absorption in medium.
	void	drop(void);

	// Change the trajectory of the photon due to scattering in the medium.
	void	spin(void);
	
	// Set the step size of the photon.
	void 	setStepSize(void);

	// Decide whether the photon should be transmitted to another layer
	// or internally reflected.
	void	transmitOrReflect(const char *);

	// Reset the Photon attributes so it can be propogated again.
	void	reset(void);
		
	// Give the photon a probabilistic chance of surviving or terminating
	// if its weight has dropped below a specified threshold.
	void	performRoulette(void);

	// Return the cartesian coordinates
	//double	getX(void) {return photonVect->location.x;}
	//double	getY(void) {return photonVect->location.y;}
	//double	getZ(void) {return photonVect->location.z;}

	// Return the direction cosines
	//double	getDirX(void) {return photonVect->direction.x;}
	//double	getDirY(void) {return photonVect->direction.y;}
	//double	getDirZ(void) {return photonVect->direction.z;}

	// Return the current weight of the photon
	double	getWeight(void) {return weight;}
	
	// Returns a random number 'n': 0 < n < 1
	double	getRandNum(void);
	
	// Return the calculated reflectance.
	double	getLayerReflectance(void);

	// Return the calculated medium reflectance (the boundary of the tissue).
	double	getMediumReflectance(void);
    
    // Return the photon's current location in the medium.
    boost::shared_ptr<Vector3d> getPhotonCoords(void) {return currLocation;}

	// Return the status of the photon.
	bool	isAlive(void) {return status;}

	// Terminate the current photon.
	void	kill(void) {status = DEAD;}

	// Calculate the new location of the photon and
	// move it to those coordinates.
	void	updateLocation(void);

	// Update weight based on specular reflectance.
	void	specularReflectance(double n1, double n2);
	
	// Update the direction cosine when internal reflection occurs on z-axis.
	void	internallyReflectZ(void);

	// Update the direction cosine when internal reflection occurs on y-axis.
	void	internallyReflectY(void);                              
    
	// Update the direction cosine when internal reflection occurs on z-axis.
	void	internallyReflectX(void);
    
	// Transmit the photon.
	void	transmit(const char *type);

	// Plot the photon's path.
	void	plotPath(void);
	
	/// Simulate photons in the medium.  This object is responsible for running 'num_iterations' of
    /// photons.
    /// NOTE:
    /// - This is the function call given to the boost thread object.
    void    Inject_photon_through_aperture(Medium *m, const int num_iterations, RNG_seed_vector *rng, Aperture *input_aperture,
                                   MC_Parameters &Params);
    
    // Hop, Drop, Spin, Roulette and everything in between.
    // NOTE: 'iterations' are the number of photons simulated by this 'Photon' object.
    void    propagatePhoton(const int iterations);
	
	// Sets initial trajectory values.
	void	initTrajectory(void);
	
	// Zero's out the local detection array.
	void	initAbsorptionArray(void);

	// Initialize the RNG.
	void	initRNG(unsigned int s1, unsigned int s2, unsigned int s3, unsigned int s4);

    // Tests if the photon will come into contact with a layer boundary
    // after setting the new step size.  If so the process of transmitting or
    // reflecting the photon begins.
    bool    checkLayerBoundary(void);
    
	// Check if photon has come into contact with a layer boundary.
	bool 	hitLayerBoundary(void);

	// Displace (i.e. update the location) the photon some distance
	// based on the pressure at that location.
	void 	displacePhotonFromPressure(void);
    
    // Displace (i.e. update the location) the photon on an arc shaped path
    // due to the changes in the refractive index changes due to changes in 
    // in pressure.
    void    displacePhotonFromRefractiveGradient(const double n1, const double n2);
    
    // Alter the optical path length due to the changes in refractive
    // index of the medium from the changes in pressure.
    void    alterOPLFromAverageRefractiveChanges(void);

    // Combination of displacement and average refractive changes with respect to influencing
    // OPL.
    void	displacePhotonAndAlterOPLFromAverageRefractiveChanges(void);

	// Write the coordinates of each scattering event to file for use
	// with plotting in matlab.
	//void 	writephoton_dataToFile(void);

	// Add the coordinates of the photon at it's current position
	// to the 'photon_data' vector.  Used for tracking the positiion of
	// scattering events in the medium.
	void	captureLocationCoords(void);

	// Add the coordinates and path length to the vector.
	void	captureExitCoordsAndLength(void);

	// Add the coordinates, path length, and weight of photon to the vector.
	void 	captureExitCoordsLengthWeight(void);

	// Check if exit location is through the aperture that will fall
	// on the detector.
	bool	didExitThroughDetectorAperture(void);

	// Write the coordinates of this photon to file.
	void	writeCoordsToFile(void);
    
    // Tests if the photon will come into contact with a medium boundary
    // after setting the new step size.  If so the process of transmitting or
    // reflecting the photon begins.
    bool    checkMediumBoundary(void);
    
	// Check if photon has left the bounds of the medium.
	bool	hitMediumBoundary(void);
    
    // Tests if the photon has crossed the plane defined by the detector.  Since
    // the detector (at this stage) only is concerned with photons that make their
    // way to the medium boundary, and would exit through the detector, we only
    // make this check in the case where the photon has hit the medium boundary.
    bool    checkDetector(void);
    
    // Check if photon has hit the detector during it's step.
    bool    hitDetector(void);
    
    // Store the energy lost into a local array that will be written to a global array
    // for all photons once they are DEAD.
    // This relieves contention between threads trying to update a single global data
    // structure and improves speed.
    void    updateLocalWeightArray(const double absorbed);

	// Write the x-y coordinates of the exit location when the photon left the medium, path length
	// and also the weight of the photon when it exited the medium.
	void	writeExitLocationsLengthWeight(void);
    
    

	
private:
	// Number of times this Photon (i.e., thread) will execute; where one execution
	// is the full cycle of photon propagation.
	int iterations;
	
	// Location of the photon with displacement from ultrasound source taken into account.
	double x_disp, y_disp, z_disp;

	
    
    // A vector object that contains the photon's location and direction.
    //boost::shared_ptr<Vector3d> photonVect;
    boost::shared_ptr<Vector3d> currLocation;
    boost::shared_ptr<Vector3d> prevLocation;
    

	// Weight of the photon.
	double	weight;
	
	// Step size for the photon.
	double	step;
	
	// Step size to boundary.  Used when calculating distance from layer
	// boundary to current position of the photon.  Specifically, it is
	// the remainder of the step size after calculating the distance to
	// the layer boundary.
	// i.e. (step_size - distance_to_boundary)/mu_t
	double step_remainder;

	// status for current photon - dead (false) or alive (true).
	bool	status;
	
	// cosine and sine theta.  Used for trajectory calculations.
	double	cos_theta, sin_theta;
	
	// The azimuthal angle
	double	psi;
	
	// The value of internal reflectance that is compared to a random
	// number (uniform between (0,1]) to determine if the photon should
	// be transmitted or reflected on a stochastic basis.
	double	reflectance;
    
    // Transmission angle for a photon when it hits a layer boundary.
    double transmission_angle;

	
	// The number of steps this photon has taken while propagating through
	// the medium.
	int num_steps;
	
	// Pointer to the medium which this photon will propagate through.
	Medium *m_medium;
    
    /// Pointer to the input aperture which photons enter the medium from.
    Aperture *m_input_aperture;

	// The thread id associated with this photon object.  The value is passed
	// in from the loop that creates the threads in main.cpp.
	int thread_id;


    // A Pseudo-random number generator with a large period.
    RNG * RNG_generator;


	RNG_seed_vector *iteration_seeds;
	
    
    // Seeds that this photon used to seed the RNG that allowed it to produce
    // steps that eventually led it to exiting the medium through the exit-aperture.
    RNGSeeds exit_seeds;

	// Tracks the path length of the photon through the medium.
    ///
	double displaced_optical_path_length;			/// OPL from only the displacement of scattering events.
    double refractiveIndex_optical_path_length;		/// OPL from only the refractive index changes.
	double combined_OPL;                            /// OPL from displacement and refractive index changes.

	// Tracks whether or not a photon has hit a medium boundary.
	bool hit_x_bound, hit_y_bound, hit_z_bound;
    
    
    // Pointer to the current layer the photon is in.
    Layer *currLayer;


    // Count through the layer.
    double cnt_through_aperture;

    // The time-of-flight of this photon bundle through the medium.
    double time_of_flight;
    
    // flags that define what (and what not) should take place during the simulation.
    bool SIM_MODULATION_DEPTH;
    bool SIM_DISPLACEMENT;
    bool SIM_REFRACTIVE_TOTAL;
    bool SIM_REFRACTIVE_GRADIENT;
    bool SAVE_RNG_SEEDS;
    bool STORE_FLUENCE;
    
    /// Fluence map for absorbed energy for this 'Photon' object. This gets put into a 'total_fluence_map' in 'Medium'
    /// after execution that this 'Photon' thread simulates all of it's photons.
    /// NOTE:
    /// - Innacuracies arise when the step size (i.e. Photon::Hop()) is larger than the voxel size.
    TRealMatrix *partial_fluence_map;

}; 		

#endif // end PHOTON_H
