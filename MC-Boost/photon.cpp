
#include "debug.h"
#include <MC-Boost/logger.h>
#include <MC-Boost/absorber.h>
#include <MC-Boost/vector3D.h>
#include <MC-Boost/layer.h>
#include <MC-Boost/medium.h>
#include <MC-Boost/photon.h>
#include <MC-Boost/displacementMap.h>
#include <MC-Boost/pressureMap.h>
#include <MC-Boost/refractiveMap.h>
#include <MC-Boost/injectionAperture.h>



#undef DEBUG




Photon::Photon(void)
{
#ifdef DEBUG
	cout << "Creating Photon...\n";
#endif

    m_medium = NULL;
    m_input_aperture = NULL;
    
	// The current location of the photon.
	currLocation = boost::shared_ptr<Vector3d> (new Vector3d);
	currLocation->withDirection();  // Enable direction for this vector.

	// The location before the photon hopped.
	prevLocation = boost::shared_ptr<Vector3d> (new Vector3d);
	prevLocation->withDirection();

	this->initCommon();
}


Photon::~Photon(void)
{
#ifdef DEBUG
	cout << "Destructing Photon...\n";
#endif

}


// Common initialization function.
void Photon::initCommon(void)
{
	// Photon just created, so it is alive.
	status = ALIVE;

	// Set to initial weight values.
	weight = 1;

	step = 0;
	step_remainder = 0;

	// Set the number of interactions back to zero.
	num_steps = 0;


	// Set the flags for hitting a layer boundary.
	hit_x_bound = hit_y_bound = hit_z_bound = false;

	// Set the transmission angle for a photon.
	transmission_angle = 0.0f;

	// Initialize the time-of-flight value.
	time_of_flight = 0.0f;

	// Set the path lengths during initialization.
	unmodulated_optical_path_length 	= 0.0f;
	displaced_optical_path_length 		= 0.0f;
	refractiveIndex_optical_path_length = 0.0f;
	combined_OPL 		= 0.0f;

    /// Everything defaults to false.
    SIM_MODULATION_DEPTH = SIM_DISPLACEMENT = SIM_REFRACTIVE_TOTAL = SIM_REFRACTIVE_GRADIENT = SAVE_RNG_SEEDS = false;

}


// Set the number of iterations this thread will run.
void Photon::setIterations(const int num)
{
	iterations = num;
}


void Photon::initTrajectory(void)
{
    double injection_point;
    injection_point = m_input_aperture->Get_radius() * sqrt(getRandNum());

	// Randomly set photon trajectory to yield anisotropic source.
	cos_theta = (2.0 * getRandNum()) - 1;
	sin_theta = sqrt(1.0 - cos_theta*cos_theta);
	psi = 2.0 * PI * getRandNum();

    /// Set the injection points to spread across the beam diameter (Guassian beam profile).
    currLocation->location.x = m_input_aperture->Get_x_coord() + injection_point*cos(psi);
    currLocation->location.y = m_input_aperture->Get_y_coord() + injection_point*sin(psi);
    currLocation->location.z = m_input_aperture->Get_z_coord();
    

	// Set the initial direction cosines for this photon.
    /// Initial direction is ballistic light as it enters the medium.  This means
    /// the photon enters the medium fully aligned with the optical axis as no scattering has
    /// taken place and it is assumed that the laser and plane of the medium where light is injected
    /// are normal to each other.
	currLocation->setDirX(0.0f);
	currLocation->setDirY(0.0f);
	currLocation->setDirZ(1.0f);

}




// BOOST thread library starts execution here.
// 1) Hop - move the photon
// 2) Drop - drop weight due to absorption
// 3) Spin - update trajectory accordingly
// 4) Roulette - test to see if photon should live or die.
void Photon::Inject_photon_through_aperture(Medium *medium, const int num_iterations, RNG_seed_vector *rng_seeds, Aperture *input_aperture,
                                            MC_Parameters &Params)
{
    
	// Before propagation we set the medium which will be used by the photon.
	this->m_medium = medium;
    
    /// We set the input aperture in which photons enter the medium.
    this->m_input_aperture = input_aperture;
    
    
	// Set the GLOBAL values that dictate what mechanisms of AO are simulated and/or if the seeds
	// for the RNG should be saved.
    SIM_MODULATION_DEPTH        = Params.MODULATION_DEPTH;
	SIM_DISPLACEMENT 			= Params.DISPLACE;
	SIM_REFRACTIVE_TOTAL        = Params.REFRACTIVE_TOTAL;
    SIM_REFRACTIVE_GRADIENT     = Params.REFRACTIVE_GRADIENT;
	SAVE_RNG_SEEDS				= Params.SAVE_SEEDS;
    
    
    // Assign this photon object a random number generator, which is passed in from main().
    RNG_generator = new RNG();
    iteration_seeds = rng_seeds;
    
    
	// Set the current layer the photon starts propagating through.  This will
	// be updated as the photon moves through layers by checking 'hitLayerBoundary'.
	currLayer = m_medium->getLayerFromDepth(currLocation->location.z);
    
	// Move the photon through the medium. 'iterations' represents the number of photons this
	// object (which is a thread) will execute.
	propagatePhoton(num_iterations);
    
}


void Photon::propagatePhoton(const int iterations)
{

	// Inject 'iterations' number of photons into the medium.
	int i;
	for (i = 0; i < iterations; i++)
	{

		/// If we are saving seeds, then we are propagating photons until they
		/// exit, using the continous stream of the RNG in the process.  Otherwise
		/// the exit seeds have already been generated, which means we have the seeds
		/// and each iteration should use seeds that detected photons (i.e. photons that
		/// exited the medium through the aperture).  So, the RNG is updated each iteration
		/// with the seed values that caused it to hop() through the medium, eventually
		/// finding its way through the exit aperture.
		if ((iteration_seeds->size() == iterations))
		{
			/// Case where we want each iteration to contain a RNG from new seeds.
			/// That is, each iteration will contain, from the start of photon propagation,
			/// seeds that had caused this photon to hop() until it exited the detection aperture.
            //RNG_seed_vector temp = (*iteration_seeds);
            //cout << temp[i].s1 << " " << temp[i].s2 << " " << temp[i].s3 << " " << temp[i].s4 << "\n";
			RNG_generator->initRNG((*iteration_seeds)[i].s1,
									(*iteration_seeds)[i].s2,
									(*iteration_seeds)[i].s3,
									(*iteration_seeds)[i].s4);
		}
		else if ((iteration_seeds->size() == 1) && (i == 0))
		{
			/// Case where only one set of seeds has been given.
			/// This is when the RNG needs to produce a continuous stream,
			/// typically when seeds are going to be saved.
            //RNG_seed_vector temp = (*iteration_seeds);
            //cout << temp[1].s1 << " " << temp[1].s2 << " " << temp[1].s3 << " " << temp[1].s4 << "\n";
            //cout.flush();
			RNG_generator->initRNG((*iteration_seeds)[i].s1,
									(*iteration_seeds)[i].s2,
									(*iteration_seeds)[i].s3,
									(*iteration_seeds)[i].s4);
		}



        if (SAVE_RNG_SEEDS || SIM_MODULATION_DEPTH)
        {
            exit_seeds = RNG_generator->getState();
        }

        // Initialize the photon's trajectory before any scattering event has taken place.
        initTrajectory();


		// While the photon has not been terminated by absorption or leaving
		// the medium we propagate it through he medium.
		while (isAlive())
		{

			// Calculate and set the step size for the photon.
			setStepSize();


			// Make various checks on the photon to see if layer or medium boundaries
			// are hit and whether the photon should be transmitted or reflected.

			// Flags for testing if a photon hit/passed through a layer
			// or medium boundary.
			//bool hitLayer = checkLayerBoundary();
			bool hitMedium = checkMediumBoundary();



			//if (!hitLayer && !hitMedium)
			if (!hitMedium)
			{
				// sanity check.
				assert(this->status == ALIVE);

				// Move the photon in the medium.
				hop();

				// Now displace the photon at its new location some distance depending on
				// how the pressure has moved scattering particles and/or due to the change
				// in the path of the photon due to refractive index gradient.
				if (SIM_REFRACTIVE_TOTAL && !SIM_DISPLACEMENT)
				{
					alterOPLFromAverageRefractiveChanges();
				}

				if (SIM_DISPLACEMENT && !SIM_REFRACTIVE_TOTAL)
				{
					displacePhotonFromPressure();
				}

				if (SIM_DISPLACEMENT && SIM_REFRACTIVE_TOTAL)
				{
					displacePhotonAndAlterOPLFromAverageRefractiveChanges();
				}



				// Drop weight of the photon due to an interaction with the medium.
				drop();

				// Calculate the new coordinates of photon propagation.
				spin();

				// Test whether the photon should continue propagation from the
				// Roulette rule.
				performRoulette();

			}

		} // end while() loop


		// Photon is no longer ALIVE. Reset the photon and start propogation over from the beginning.
		reset();

	} // end for() loop


}




void Photon::reset()
{
#ifdef DEBUG
	cout << "Resetting Photon...\n";
#endif

	// Photon just created, so it is alive.
	status = ALIVE;

	// Set back to initial weight values.
	weight = 1;

    /// Zero out the location of the photon for peace of mind. This will
    /// be updated in 'initTrajectory' before propagation begins again.
	currLocation->location.x = 0.0f;
	currLocation->location.y = 0.0f;
	currLocation->location.z = 0.0f;

    /// Zero out the step size.
	step = 0;
	step_remainder = 0;

	// Reset the number of interactions back to zero.
	num_steps = 0;

	// Reset the path lengths of the photon.
	unmodulated_optical_path_length 	= 0.0f;
	displaced_optical_path_length 		= 0.0f;
	refractiveIndex_optical_path_length = 0.0f;
	combined_OPL				 		= 0.0f;


	// Reset the flags for hitting a layer boundary.
	hit_x_bound = false;
	hit_y_bound = false;
	hit_z_bound = false;

	// Reset the transmission angle for a photon.
	transmission_angle = 0;

	// Set photon trajectory to yield a uniform distribution of light over a given beam radius.
	initTrajectory();

	// Reset the current layer from the injection coordinates of the photon.
	currLayer = m_medium->getLayerFromDepth(m_input_aperture->Get_z_coord());

}


// Update the step size for the photon.
void Photon::setStepSize()
{
	// Update the current values of the absorption and scattering coefficients
	// based on the depth in the medium (i.e. which layer the photon is in).
	double mu_a = currLayer->getAbsorpCoeff(currLocation);
	double mu_s = currLayer->getScatterCoeff(currLocation);


	// If last step put the photon on the layer boundary
	// then we need a new step size.  Otherwise, the step
	// is set to remainder size and update step_remainder to zero.
	if (step_remainder == 0.0) {
		double rnd = getRandNum();
		// Calculate the new step length of the photon.
		step = -log(rnd)/(mu_a	+ mu_s);
	}
	else
	{
		step = step_remainder/(mu_a + mu_s);
		step_remainder = 0.0;
	}
}




// Tests if the photon will come into contact with a layer boundary
// after setting the new step size.  If so the process of transmitting or
// reflecting the photon begins.
bool Photon::checkLayerBoundary(void)
{
	if (hitLayerBoundary())
	{
#ifdef DEBUG
		cout << "Hit layer boundary\n";
#endif

		hop();	// Move the photon to the layer boundary.
		transmitOrReflect("layer");
		return true;
	}

	return false;
}


bool Photon::checkMediumBoundary(void)
{
	if (hitMediumBoundary())
	{
#ifdef DEBUG
		cout << "Hit medium boundary\n";
#endif
		hop();  // Move the photon to the medium boundary.
		transmitOrReflect("medium");
		return true;
	}

	return false;
}


// Check if the photon passed through the detection window.
// Returns the number of detectors that were crossed in this
// hop in the case of multiple detectors present.
bool Photon::checkDetector(void)
{
	int cnt =  m_medium->photonHitDetectorPlane(currLocation);
	// If cnt > 0 the photon exited through the bounds of the detector.
	if (cnt > 0)
	{
		return true;
	}
	else
		return false;
}


// Step photon to new position.
void Photon::hop()
{
#ifdef DEBUG
	cout << "Hopping...\n";
#endif

	num_steps++;


	// Save the location before making the hop.
	prevLocation->location.x = currLocation->location.x;
	prevLocation->location.y = currLocation->location.y;
	prevLocation->location.z = currLocation->location.z;
	prevLocation->setDirX(currLocation->getDirX());
	prevLocation->setDirY(currLocation->getDirY());
	prevLocation->setDirZ(currLocation->getDirZ());


	// Update the location
	currLocation->location.x += step * currLocation->getDirX();
	currLocation->location.y += step * currLocation->getDirY();
	currLocation->location.z += step * currLocation->getDirZ();


}


// Return absorbed energy from photon weight at this location.
void Photon::drop()
{
#ifdef DEBUG
	cout << "Dropping...\n";
#endif

	if (this->status == DEAD) return;




	double mu_a = 0.0f;
	double mu_s = 0.0f;
	double albedo = 0.0f;
	double absorbed = 0.0f;


	Absorber * absorber = currLayer->getAbsorber(currLocation);
	// If an absorber was returned, then we get the absorption and
	// scattering coefficients from it.  Otherwise we use the values
	// from the background layer.
	if (absorber != NULL)
	{
		mu_a = absorber->getAbsorberAbsorptionCoeff();
		mu_s = absorber->getAbsorberScatteringCoeff();

		// Calculate the albedo and remove a portion of the photon's weight for this
		// interaction.
		albedo = mu_s / (mu_a + mu_s);
		absorbed = weight * (1 - albedo);

		// Update the absorbed weight in this absorber.
		absorber->updateAbsorbedWeight(absorbed);


	}
	else
	{
		// Update the current values of the absorption and scattering coefficients
		// based on the depth in the medium (i.e. which layer the photon is in).
		// NOTE:
		// - No need to index into the layer and see if absorption and scattering coefficients
		//   should be pulled from absorber, because we verified above in the if() that this was
		//   not the case.  Saves a small amount of time searching through the absorbers.
		mu_a = currLayer->getAbsorpCoeff();
		mu_s = currLayer->getScatterCoeff();

		// Calculate the albedo and remove a portion of the photon's weight for this
		// interaction.
		albedo = mu_s / (mu_a + mu_s);
		absorbed = weight * (1 - albedo);
	}



	// Remove the portion of energy lost due to absorption at this location.
	weight -= absorbed;

	// Deposit lost energy in the grid of the medium.
	//m_medium->absorbEnergy(z, absorbed);


}




// Calculate the new trajectory of the photon.
void Photon::spin()
{
#ifdef DEBUG
	cout << "Spinning...\n";
#endif

	if (this->status == DEAD) return;

    // Get the anisotropy factor from the current location
    // of the photon in the medium.
    double g = currLayer->getAnisotropy(currLocation);

	double rnd = getRandNum();

	double dirZ = currLocation->getDirZ();
	double dirY = currLocation->getDirY();
	double dirX = currLocation->getDirX();

	if (g == 0.0) {
		cos_theta = (2.0 * rnd) - 1.0;
	}
	else {
		double temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd);
		cos_theta = (1.0 + g*g - temp*temp)/(2.0*g);
	}
	sin_theta = sqrt(1.0 - cos_theta*cos_theta); /* sqrt() is faster than sin(). */

	// Sample 'psi'.
	psi = 2.0 * PI * getRandNum();
	double cos_psi = cos(psi);
	double sin_psi = 0.0;
	if (psi < PI) {
		sin_psi = sqrt(1.0 - cos_psi*cos_psi);     /* sqrt() is faster than sin(). */
	}
	else {
		sin_psi = -sqrt(1.0 - cos_psi*cos_psi);
	}

	double uxx, uyy, uzz;
	double temp;
	/* New trajectory. */
	if (1 - fabs(dirZ) <= ONE_MINUS_COSZERO) {      /* close to perpendicular. */
		uxx = sin_theta * cos_psi;
		uyy = sin_theta * sin_psi;
		uzz = cos_theta * SIGN(dirZ);   /* SIGN() is faster than division. */
	}
	else {					/* usually use this option */


		temp = sqrt(1.0 - dirZ * dirZ);
		uxx = sin_theta * (dirX * dirZ * cos_psi - dirY * sin_psi) / temp + dirX * cos_theta;
		uyy = sin_theta * (dirY * dirZ * cos_psi + dirX * sin_psi) / temp + dirY * cos_theta;
		uzz = -sin_theta * cos_psi * temp + dirZ * cos_theta;
	}

	// Update trajectory.
	currLocation->setDirX(uxx);
	currLocation->setDirY(uyy);
	currLocation->setDirZ(uzz);
}



void Photon::performRoulette(void)
{
	// Photon has already been killed, presumably by leaving the medium.
	if (this->status == DEAD) return;


	if (weight < THRESHOLD) {
		if (getRandNum() <= CHANCE) {
			weight /= CHANCE;
		}
		else {
#ifdef DEBUG
			cout << "Photon died in Roulette\n";
#endif
			status = DEAD;
		}
	}
}



void Photon::displacePhotonFromPressure(void)
{

	// Photon does not get displaced on boundaries of medium.
	if (hit_x_bound || hit_y_bound || hit_z_bound) return;

#ifdef DEBUG
	int Nx = m_medium->kwave.dmap->getnumVoxelsXaxis();
	int Ny = m_medium->kwave.dmap->getnumVoxelsYaxis();
	int Nz = m_medium->kwave.dmap->getnumVoxelsZaxis();

	assert((_x < Nx && _x >= 0) &&
			(_y < Ny && _y >= 0) &&
			(_z < Nz && _z >= 0) ||
			assert_msg("_x=" << _x << " _y=" << _y << " _z=" << _z << "\n"
					<< currLocation->location.x << " "
					<< currLocation->location.y << " "
					<< currLocation->location.z));
#endif



	// Transform the location of the photon in the medium to discrete locations in the grid.
	//
	static double dx = m_medium->dx;
	static double Nx = m_medium->Nx;

	static double dy = m_medium->dy;
	static double Ny = m_medium->Ny;

	static double dz = m_medium->dz;
	static double Nz = m_medium->Nz;

    /// Get the appropriate voxel number for the 3-D grid.
	int _x = currLocation->location.x/dx - (currLocation->location.x/dx)/Nx;
	int _y = currLocation->location.y/dy - (currLocation->location.y/dy)/Ny;
	int _z = currLocation->location.z/dz - (currLocation->location.z/dz)/Nz;




	// Subtle case where index into grid is negative because of rounding errors above.
#ifdef DEBUG
	if (_z < 0 || _x < 0 || _y < 0)
	{
		// FIXME:
		// - This should be removed when detection of line-plane intersection
		//   is implemented.  For now, it is a necessary evil.
		cout << "Error in array index calculation: Photon::displacePhotonFromPressure()\n";
		this->status = DEAD;
		return;
	}
#endif






    // Get the displacement caused by ultrasound pressure at this location of the scattering event.
	//
    // Index into the displacement grids to retrieve the value of how much to
    // displace the photon from it's current location along each axis (x,y,z).
    double x_displacement = m_medium->kwave.dmap->getDisplacementFromGridX(_x, _y, _z);
    double y_displacement = m_medium->kwave.dmap->getDisplacementFromGridY(_x, _y, _z);
    double z_displacement = m_medium->kwave.dmap->getDisplacementFromGridZ(_x, _y, _z);

    /// To mitigate any numerical instability in the displacement values that is
    /// compounded over the time ultrasound propagates through the medium, we threshold
    /// the displacement such that we don't displace this scattering event until the
    /// displacement value exceeds the DISPLACEMENT_THRESHOLD value, which can only be caused
    /// when we are in (or near) the ultrasound focus.
    /// FIXME:
    /// - This needs to be taken care of in the computation of the displacement values, so that
    ///   this check never need take place.
    ///
    ///const float DISPLACEMENT_THRESHOLD = 0.0f; /// (meters)
    ///if (abs(x_displacement) >= DISPLACEMENT_THRESHOLD) currLocation->location.x += x_displacement;
    ///if (abs(y_displacement) >= DISPLACEMENT_THRESHOLD) currLocation->location.y += y_displacement;
    ///if (abs(z_displacement) >= DISPLACEMENT_THRESHOLD) currLocation->location.z += z_displacement;
    currLocation->location.x += x_displacement;
    currLocation->location.y += y_displacement;
    currLocation->location.z += z_displacement;


	// Get the local refractive index based on the coordinates of the displaced photon.
	/// NOTE:
	/// - Not needed since we assume the mechanism of modulation is strictly the displacement
	///   and nothing to do with spatially varying refractive index values.  Therefore, use
	///   the background refractive index.
    /// FIXME:
	/// - How does this change with non-uniform density???
    /// - Need to have a refractive index map that contains the variations (without US modulation induced changes),
    ///   which requires the creation of an unmodulated nmap (i.e. make the computation at t=0 and store it without
    ///   further pertubation).
	//double local_refractive_index = m_medium->kwave.nmap->getRefractiveIndexFromGrid(_x, _y, _z);
	double local_refractive_index = 1.33;

	// Get the distance between the previous location and the current location.
	double distance_traveled = VectorMath::Distance(prevLocation, currLocation);


	// Update the optical path length of the photon through the medium by
	// calculating the distance between the two points and multiplying by the refractive index.
	//displaced_optical_path_length += VectorMath::Distance(prevLocation, currLocation) * currLayer->getRefractiveIndex();
	displaced_optical_path_length += distance_traveled * local_refractive_index;
}



// XXX: After simulation of optical path lengths (OPL) it has been verified
//      that the use of a piecewise integration of refractive indices, and
//      their respective distances (i.e. their distance through the corresponding voxel)
//      is innaccurate.  This is because the distance calculated through each voxel is never
//      completely exact, and is in truth quite erroneous.  It has been found it is better to find the refractive index
//      values the photon bundle has moved through during this step, average those values,
//      and multiply by the accurate distance.
// ------------------------------------------------------------------------------------------
/*
// Alter the optical path length of the photon through the medium as it encounters
// refractive index changes along it's step.
// Basic idea is that the line segment from the photon step encounters changes
// in the mediums refractive indeces.  From Fermat's theorem we know that the
// arc (S) that the line segment will make is an extremum based on the values
// of refractive index it encounters on the physical path length (L).  The change
// in the optical length is given by L*refractive_index.
// By using the parametric equation of a line from the previous and current locations
// of the photon it is possible to move along this line, index into the K-Wave supplied
// grid of pre-calculated refractive index changes, and alter the pathlength accordingly.
void Photon::alterPathLengthFromRefractiveChanges(void)
{
	// Photon does not get its optical path length changed on boundaries of medium.
	if (hit_x_bound || hit_y_bound || hit_z_bound) return;

	// Current value of the change in refractive index (i.e. n(x,y,z) coordinates).
	double n_curr = 0.0f;
	// Previous value of the change in refractive index.
	double n_prev = 0.0f;
	// The value of the background index of refraction (i.e. no change due to pressure).
	double n_background = 0.0f;

	// Scales the position along the line segment from previous to current location.
	// i.e. P(t) = prevLocation + t*currLocation
	double t_curr = 0.0f;
	double t_prev = 0.0f;

	// The adiabatic piezo-optical coefficient of the medium.
	double eta = 0.3211;

	// The wave number of the acoustic wave.
	double k = 2*PI/m_medium->kwave.transducerFreq;

	// Create the new vector that is the line segment from the previous position of the photon (p0) to
	// it's current position (p1).
	boost::shared_ptr<Vector3d> lineSegment = (*currLocation) - (*prevLocation);

	// These track the distance on the line segment the photon travels that contains it's respective
	// index of refraction.  That is, by tracking where a given point in space starts and ends with
	// a given value of refractive index, we know which portion (i.e. the distance from
	// 'previousPointOnSegment' to 'pointOnSegment') should be multiplied by it's corresponding index
	// of refraction, thus yielding an optical path length dependent on the gradient changes of the
	// refraction of index due to variations in pressure.
	boost::shared_ptr<Vector3d> pointOnSegment;
	boost::shared_ptr<Vector3d> previousPointOnSegment = prevLocation;  // Set to last end point from hopping in medium.


	// Prime the check below by setting the previous and current value of the refractive index changes to the same values.
	//
	// The grid indices that are calculated form the photon's position.
	// Transform the location of the photon in the medium to discrete locations in the grid.
	//
	double dx = m_medium->kwave.pmap->getDx();
	double Nx = m_medium->kwave.pmap->getNumVoxelsXaxis();

	double dy = m_medium->kwave.pmap->getDy();
	double Ny = m_medium->kwave.pmap->getNumVoxelsYaxis();

	double dz = m_medium->kwave.pmap->getDz();
	double Nz = m_medium->kwave.pmap->getNumVoxelsZaxis();



	int _x = previousPointOnSegment->location.x/dx - (previousPointOnSegment->location.x/dx)/Nx;
	int _y = previousPointOnSegment->location.y/dy - (previousPointOnSegment->location.y/dy)/Ny;
	int _z = previousPointOnSegment->location.z/dz - (previousPointOnSegment->location.z/dz)/Nz;


	// Prime the refractive index values so they can be compared as 'n_curr' is updated based on position.
	n_background = currLayer->getRefractiveIndex();
	n_curr = n_prev = n_background + n_background*eta*k*m_medium->kwave.pmap->getPressureFromGrid(_x, _y, _z);




	// Because we must make the integration of the refractive index changes discrete, we
	// have a step size (ds) along the arc (S).
	// More specifically, value of (steps) dictates how many portions (i.e. the resolution) we take along the arc (S).
	int steps = 8;
	double ds = 1/(double)steps;
	for (int i = 0; i < steps; i++, t_curr += ds)
	{
		// Update the (t_curr) value to move along the line segment.
		// 'prevLocation' is where the last hop moved the photon. 'lineSegment'
		// is the line formed form 'prevLocation' to the new location of the photon 'currLocation'.
		// 'pointOnSegment' is a point on that line ('lineSegment') based on the value of 't'.
		//t_curr += ds;
		pointOnSegment = (*prevLocation) + (*((*lineSegment)*t_curr));

		// Transform the location of the photon in the medium to discrete locations in the grid.
		//
		_x = pointOnSegment->location.x/dx - (pointOnSegment->location.x/dx)/Nx;
		_y = pointOnSegment->location.y/dy - (pointOnSegment->location.y/dy)/Ny;
		_z = pointOnSegment->location.z/dz - (pointOnSegment->location.z/dz)/Nz;


#ifdef DEBUG
		int Nx = m_medium->kwave.pmap->getnumVoxelsXaxis();
		int Ny = m_medium->kwave.pmap->getnumVoxelsYaxis();
		int Nz = m_medium->kwave.pmap->getnumVoxelsZaxis();

		assert((_x < Nx && _x >= 0) &&
				(_y < Ny && _y >= 0) &&
				(_z < Nz && _z >= 0) ||
				assert_msg("_x=" << _x << " _y=" << _y << " _z=" << _z << "\n"
						<< currLocation->location.x << " "
						<< currLocation->location.y << " "
						<< currLocation->location.z));
#endif

		// Update the value of the refractive index based on coordinates that retrieve the
		// current discrete location of the pressure from the pressure map grid.
		//
		// background index of refraction.  Might be different within some occlusion, so we pass
		// in coordinates to get the value at a specific location.  If simulating a homogenous medium
		// a constant can be placed here, which should speed up the simulation since no searching will occur.
		// n_background = currLayer->getRefractiveIndex(pointOnSegment);
		n_curr = n_background + n_background*eta*k*m_medium->kwave.pmap->getPressureFromGrid(_x, _y, _z);

		// If there was a change in refractive index value, due to a change in pressure found in the pressure
		// grid, we accumulate this portion of the step, which is multiplied by this portion's refractive index
		// value, and update the optical path length accordingly.
		//
		if (n_curr != n_prev)
		{

			// Displace the photon on its arced path due to the refractive index changes along the step.
			//displacePhotonFromRefractiveGradient(n_prev, n_curr);

			// The point on the line segment where the this refractive index change started.
			previousPointOnSegment = (*prevLocation) + (*((*lineSegment)*t_prev));

			// The current point on the segment has moved to another voxel in the K-Wave simulation, therefore
			// we need to step back a portion of 'ds' that was moved so we have the length over a given
			// refractive index change, not just beyond.  Therefore we subtract off the portion that was added
			// to 'ds' when updating the current point.
			pointOnSegment = (*prevLocation) + (*((*lineSegment)*(t_curr)));

			refractiveIndex_optical_path_length += VectorMath::Distance(pointOnSegment, previousPointOnSegment) * n_curr;

			// Update the values for the next iterations.
			n_prev = n_curr;
			t_prev = t_curr;

		}


	} // end for() loop


	// If the last portion of the line segment (pointOnSegment - previousPointOnSegment) didn't see an index of
	// refraction change, its optical path length won't be updated.  Also, if the entire step fit into a voxel, this takes
	// that into account as well.
	refractiveIndex_optical_path_length += VectorMath::Distance(pointOnSegment, previousPointOnSegment) * n_curr;

}
 */

void Photon::alterOPLFromAverageRefractiveChanges(void)
{
	// Photon does not get its optical path length changed on boundaries of medium.
	if (hit_x_bound || hit_y_bound || hit_z_bound) return;


	static double dx = m_medium->dx;
	static double Nx = m_medium->Nx;

	static double dy = m_medium->dy;
	static double Ny = m_medium->Ny;

	static double dz = m_medium->dz;
	static double Nz = m_medium->Nz;


	// Used for discrete values that will index into the appropriate voxel in the refractive map grid.
	//
	int _x = 0;
	int _y = 0;
	int _z = 0;

	// Current value of the change in refractive index (i.e. n(x,y,z) coordinates).
	double n_cumulative = 0.0f;

	// Scales the position along the line segment from previous to current location.
	// i.e. P(t) = prevLocation + t*currLocation
	//double t_curr = 0.0f;

	// NOTE:
	// - Because the step size is much smaller than the size of a voxel, which means any given step
	//   at best will at one end lie in a voxel adjecent to the other end, we only have a chance of
	//   having different values of refractive index at the far ends of the segment connecting the previous location
	//   and the current location.  Therefore, only check those to simplify things.
	//   -----------------------------------------------------------------------------------------------------------
	// Transform the location of the photon in the medium to discrete locations in the grid based
    // on the photon's previous location in the medium.
	//
	_x = prevLocation->location.x/dx - (prevLocation->location.x/dx)/Nx;
	_y = prevLocation->location.y/dy - (prevLocation->location.y/dy)/Ny;
	_z = prevLocation->location.z/dz - (prevLocation->location.z/dz)/Nz;

	// Accumulate the refractive index encountered along the segment between scattering events for this
	// step through the medium.
	n_cumulative = m_medium->kwave.nmap->getRefractiveIndexFromGrid(_x, _y, _z);

	// Transform the location of the photon in the medium to discrete locations in the grid based
    // on the photon's current location in the medium after 'hopping'.
	//
    _x = currLocation->location.x/dx - (currLocation->location.x/dx)/Nx;
	_y = currLocation->location.y/dy - (currLocation->location.y/dy)/Ny;
	_z = currLocation->location.z/dz - (currLocation->location.z/dz)/Nz;

	// Accumulate the refractive index encountered along the segment between scattering events for this
	// step through the medium.
	n_cumulative += m_medium->kwave.nmap->getRefractiveIndexFromGrid(_x, _y, _z);

	// Create the average refractive index.
	double n_avg = n_cumulative / 2.0;

	// Get the distance between the previous location and the current location.
	double distance_traveled = VectorMath::Distance(prevLocation, currLocation);
	//cout << "distance = " << distance_traveled << "\n";
	//cout << "n_avg = " << n_avg << "\n";

	// Make the final calculation of the OPL.
	refractiveIndex_optical_path_length += distance_traveled * n_avg;

/*

    // Create the new vector that is the line segment from the previous position of the photon (p0) to
	// it's current position (p1).
	boost::shared_ptr<Vector3d> lineSegment = (*currLocation) - (*prevLocation);

	// These track the distance on the line segment the photon travels that contains it's respective
	// index of refraction.  That is, by tracking where a given point in space starts and ends with
	// a given value of refractive index, we know which portion (i.e. the distance from
	// 'previousPointOnSegment' to 'pointOnSegment') should be multiplied by it's corresponding index
	// of refraction, thus yielding an optical path length dependent on the gradient changes of the
	// refraction of index due to variations in pressure.
	boost::shared_ptr<Vector3d> pointOnSegment;

    // Because we must make the integration of the refractive index changes discrete, we
	// have a step size (ds) along the arc (S).
	// More specifically, value of (num_differential_steps) dictates how many portions (i.e. the resolution) we take along the arc (S).

	int num_differential_steps = 1;
	double ds = 1/(double)num_differential_steps;
	for (int i = 0; i <= num_differential_steps; i++, t_curr += ds)
	{
		// Update the (t_curr) value to move along the line segment.
		// 'prevLocation' is where the last hop moved the photon. 'lineSegment'
		// is the line formed form 'prevLocation' to the new location of the photon 'currLocation'.
		// 'pointOnSegment' is a point on that line ('lineSegment') based on the value of 't'.
		//t_curr += ds;
		pointOnSegment = (*prevLocation) + (*((*lineSegment)*t_curr));

		// Transform the location of the photon in the medium to discrete locations in the grid.
		//
		_x = pointOnSegment->location.x/dx - (pointOnSegment->location.x/dx)/Nx;
		_y = pointOnSegment->location.y/dy - (pointOnSegment->location.y/dy)/Ny;
		_z = pointOnSegment->location.z/dz - (pointOnSegment->location.z/dz)/Nz;


		// Accumulate the refractive index encountered along the segment between scattering events for this
		// step through the medium.
		n_cumulative += m_medium->kwave.nmap->getRefractiveIndexFromGrid(_x, _y, _z);
	} // end for() loop


	// Calculate the OPL through the average index of refraction encountered.
	double distance_traveled = VectorMath::Distance(prevLocation, pointOnSegment);  // The distance between the two scattering locations (i.e. length of line segment)
	double n_avg = n_cumulative / num_differential_steps; // Average refractive index.
	refractiveIndex_optical_path_length +=  distance_traveled * n_avg;
*/


	// Update the time-of-flight value for this portion of propagation.
	time_of_flight += distance_traveled / (3e8/n_avg);

}

// XXX: FINISH ME
// - Combine both displacement and average refractive index changes to OPL contribution.
void Photon::displacePhotonAndAlterOPLFromAverageRefractiveChanges(void)
{
	displacePhotonFromPressure();
	alterOPLFromAverageRefractiveChanges();

	/// Due to the order in which the above functions are called, the OPL of
	/// 'refractiveIndex_optical_path_length' contains the OPL under the influence of
	/// both AO mechanisms (displacement, pressure induced refractive index).  However,
	/// it's possible to run the simulation with any combination (all mechanisms on, one, etc.).
	/// If we make it here, both have been turned on, and to remove any disambiguity we assign
	/// 'combined_OPL' the value of 'refractiveIndex_optical_path_length' to
	/// reduce the overhead of having to redo what is essentially embedded in 'refractiveIndex_optical_path_length'.
    /// FIXME:
    /// - This is very confusing and should be updated someday to track non-displaced coordinates for separating
    ///   the mechanisms clearly.
    /// - Or another way?
	combined_OPL = refractiveIndex_optical_path_length;

}


// FIXME:
// - This currently does NOT work correctly.
// - Need to calculate the change in direction as the photon moves into a refractive
//   index change.
// - Need to figure out the orientation of the plane that the photon passes through
//   which causes refraction.  Is it normal to the direction of photon propagation?
//   Is it normal to the direction of ultrasound.
void Photon::displacePhotonFromRefractiveGradient(const double n1, const double n2)
{


}



// FIXME:
// - The photon, if moving back up to the 'Air' layer, should be able to transmit out of the medium.
//   Currently this does not happen (maybe because of reflection on the layer boundary?).
// - If the photon hits the layer at the bottom of the medium, it gets stuck there until
//   it probabilistically is transmitted through the medium.  It should hit the layer, check
//   the medium, and determine if it should transmit or reflect.  This bug only arises occasionally.
void Photon::transmit(const char *type)
{


	// 'tempLayer' is used
	Layer *tempLayer = NULL;

	if (strcmp("layer", type) == 0)
	{
#ifdef DEBUG
		cout << "Transmitting through layer\n";
		cout << currLocation;
#endif


		// If we transmit through the layer to another we must update
		// the current layer pointer of the photon so it will correctly
		// calculate the next step size.
		if (currLocation->getDirZ() > 0)
			tempLayer = m_medium->getLayerBelowCurrent(currLocation->location.z);
		else
			tempLayer = m_medium->getLayerAboveCurrent(currLayer);




		// If 'tempLayer' is NULL we are at the edge of the medium since
		// the layer above or below does not exist.  Therefore we decide if the
		// photon should reflect or transmit from the medium boundary.
		if (tempLayer == NULL)
		{
			// Since 'tempLayer' is NULL we have tried to retrieve a layer outside of the
			// medium's bounds.  Therefore 'currLayer' still is valid (i.e. not NULL) and we check
			// if the photon should reflect or transmit through the medium.
			if (hitMediumBoundary())
			{
				hop(); // Move the photon to the medium boundary.
				transmitOrReflect("medium");
				return;
			}
		}
		else // Transmitted to new layer.
		{
            
            /// FIXME:
            /// - If multiple layers are ever used, with different spatial refractive index values,
            ///   this all becomes error prone. Especially under insonification at the medium boundaries
            ///   and use with 'AO_sim'.
			// Store index of refraction from layer the photon is coming from.
			double ni = currLayer->getRefractiveIndex();

			// Get the new index of refraction from the layer the photon is moving to.
			double nt = tempLayer->getRefractiveIndex();

			// Since the photon is being transmitted, update to the pointer to the new layer that
			// was found above.
			currLayer = tempLayer;

			// NOTE: Assumes layers are normal to initial direction, which is dirZ.
			// Set the direction cosines.
			currLocation->setDirZ(cos(this->transmission_angle));
			currLocation->setDirY(currLocation->getDirY()*ni/nt);
			currLocation->setDirX(currLocation->getDirX()*ni/nt);
		}



	}
	else if (strcmp("medium", type) == 0)
	{
#ifdef DEBUG
		cout << "Transmitting through medium boundary\n";
		cout << currLocation;
#endif
		// If we transmit through the medium, meaning the photon exits the medium,
		// we see if the exit location passed through the detector.  If so, the exit
		// data is written out to file.
		if (checkDetector())
		{

            // Update the direction cosines upon leaving the medium so calculations can be mode
			// to see if this photon makes it's way to the detector (i.e. CCD camera).
			// NOTE:
			// - We assume outside of the medium is nothing but air with an index of refraction 1.0
			//   so no division for x & y direction cosines because nt = 1.0
			// - It is also assumed that the photon is transmitted through the x-y plane.

            if (SIM_DISPLACEMENT || SIM_REFRACTIVE_TOTAL || SIM_REFRACTIVE_GRADIENT)
         	{

                /// FIXME:
                /// - This will produce an error in the exit angles if insonification occurs
                ///   near the exit aperture. This is because the background medium (i.e. unmodulated)
                ///   is being probed for the refractive index value.
                double ni = currLayer->getRefractiveIndex();
                currLocation->setDirZ(cos(this->transmission_angle));
                currLocation->setDirY(currLocation->getDirY()*ni);
                currLocation->setDirX(currLocation->getDirX()*ni);

                // Write exit data via the logger.
                Logger::getInstance()->Write_weight_OPLs_coords(*this);

                /// Is simulation of modulation depth enabled? And if so store the values for later use for
                /// calculating modulation depth.
                if (SIM_MODULATION_DEPTH)
                {
                    Logger::getInstance()->Store_OPL(exit_seeds, displaced_optical_path_length, refractiveIndex_optical_path_length);
                }
         	}

			// Write out the seeds that caused this photon to hop, drop and spin its way out the
			// exit-aperture.
			//
			if (SAVE_RNG_SEEDS)
			{
				Logger::getInstance()->writeRNGSeeds(exit_seeds.s1, exit_seeds.s2, exit_seeds.s3, exit_seeds.s4);
			}

		}

		// The photon has left the medium, so kill it.
		this->status = DEAD;
	}
	else
	{
		cout << "Error:  Wrong paramater passed to Photon::transmit()\n";
	}

}


void Photon::transmitOrReflect(const char *type)
{
	// Test whether to transmit the photon or reflect it.  If the reflectance is
	// greater than the random number generated, we internally reflect, otherwise
	// transmit.
	// Case when interaction is with a layer boundary.
	if (strcmp("layer", type) == 0)
	{
		// Stochastically determine if the photon should be transmitted or reflected.
		if (getLayerReflectance() > getRandNum())
		{
#ifdef DEBUG
			cout << "Internally reflecting on layer boundary\n";
#endif
			internallyReflectZ();

			// Since the photon has interacted with the tissue we deposit weight.
			drop();

			// Perform roulette to determine if the photon should continue propagation.
			// NOTE: spin() is not performed since the photon has been internally reflected
			//	   	and the new direction has been taken care of by internallyReflect().
			performRoulette();

		}
		else
		{
#ifdef DEBUG
			cout << "Transmitting through layer\n";
#endif
			transmit("layer");
		}
	}
	// Case when interaction is with a medium boundary.
	else if (strcmp("medium", type) == 0)
	{
		// Stochastically determine if the photon should be transmitted or reflected.
		if (getMediumReflectance() > getRandNum())
		{
#ifdef DEBUG
			cout << "Reflecting photon on medium boundary\n";
#endif
			// Depending on which medium boundary was hit, we reflect on that axis,
			// change the direction of the direction cosign, and reset the boolean flag.
			if (hit_x_bound)
			{
				internallyReflectX();
			}
			else if (hit_y_bound)
			{
				internallyReflectY();
			}
			else if (hit_z_bound)
			{
				internallyReflectZ();
			}
			else
			{
				cout << "Error, no medium boundary hit\n";
			}


			// Since the photon has interacted with the tissue we deposit weight.
			drop();

			// Perform roulette to determine if the photon should continue propagation.
			// NOTE: spin() is not performed since the photon has been internally reflected
			//	   	and the new direction has been taken care of by internallyReflect().
			performRoulette();
		}
		else
		{
			transmit("medium");
		}
	}
	else
	{
		cout << "Error: transmitOrReflect()\n";
	}
}



// XXX: *** Need to verify the logic below is correct ***
double Photon::getMediumReflectance(void)
{
	// Sanity check.
	assert((hit_x_bound == true) ||
			(hit_y_bound == true) ||
			(hit_z_bound == true));


	//Layer *currLayer = m_medium->getLayerFromDepth(currLocation->location.z);
	double refract_index_n1 = currLayer->getRefractiveIndex();	// Current layer's refractive index.
	double refract_index_n2 = 1.0;	// Outside of the medium is only air.


	double axis_direction = 0.0;
	if (hit_x_bound)
		axis_direction = currLocation->getDirX();
	else if (hit_y_bound)
		axis_direction = currLocation->getDirY();
	else
		axis_direction = currLocation->getDirZ();

	// Calculate the incident angle based on the axis in which the photon hit the medium.
	double incident_angle = acos(abs(axis_direction));

	// Calculate the critical angle.
	double critical_angle = asin(refract_index_n2 / refract_index_n1);

	// If the incident angle is larger than the critical angle, the reflectance
	// is set to 1, otherwise the reflectance is calculated from Fresnel's equation.
	if (incident_angle > critical_angle)
	{
		reflectance = 1;
	}
	else
	{
		// Calculate the transmission angle through the medium boundary.
		this->transmission_angle = asin(refract_index_n1/refract_index_n2 * sin(incident_angle));

		// Calcualte the Fresnel reflection
		reflectance = 0.5 * (pow(sin(incident_angle-transmission_angle), 2)/pow(sin(incident_angle+transmission_angle), 2)
				+ pow(tan(incident_angle-transmission_angle), 2)/pow(tan(incident_angle+transmission_angle), 2));
	}


	return reflectance;

}


double Photon::getLayerReflectance(void)
{
	double refract_index_n1 = 0.0;	// Current layer's refractive index.
	double refract_index_n2 = 0.0;	// Next layer's refractive index.
	Layer *nextLayer = NULL;

	double incident_angle = acos(abs(currLocation->getDirZ()));
	refract_index_n1 = currLayer->getRefractiveIndex();

	// If the photon is moving towards a deeper layer.
	if (currLocation->getDirZ() > 0)
	{
		nextLayer = m_medium->getLayerBelowCurrent(currLocation->location.z);
	}
	// If the photon is moving towards a more shallow layer.
	else if (currLocation->getDirZ() < 0)
	{
		nextLayer = m_medium->getLayerAboveCurrent(currLayer);
	}
	// Perpendicular propagation.
	else
	{
		// FIXME:
		// This is where propagation is normal to the boundary/layer
		// and should be transmitted regardless.
		cout << "Photon should transmit...\n";
	}


	// If the layer above or below is outside of the medium
	// we assign the refractive index to be air, otherwise
	// use the value from the layer.
	if (nextLayer == NULL)
		refract_index_n2 = 1.0;
	else
		refract_index_n2 = nextLayer->getRefractiveIndex();


	// Calculate the critical angle.
	if (refract_index_n2 > refract_index_n1)
	{
		// For specular reflection we always remove some portion of the weight and transmit the photon
		// to the next layer.
		transmission_angle = asin(refract_index_n1/refract_index_n2 * sin(incident_angle));
		specularReflectance(refract_index_n1, refract_index_n2);

		// Since we know we transmit, set reflectance to zero.
		reflectance = 0;
	}
	else
	{
		double critical_angle = asin(refract_index_n2 / refract_index_n1);
		// If the incident angle is larger than the critical angle, the reflectance
		// is set to 1, otherwise the reflectance is calculated from Fresnel's equation.
		if (incident_angle > critical_angle)
		{
			reflectance = 1;
		}
		else
		{
			transmission_angle = asin(refract_index_n1/refract_index_n2 * sin(incident_angle));

			reflectance = 0.5 * (pow(sin(incident_angle-transmission_angle), 2)/pow(sin(incident_angle+transmission_angle), 2)
					+ pow(tan(incident_angle-transmission_angle), 2)/pow(tan(incident_angle+transmission_angle), 2));
		}
	}

	return reflectance;
}




// FIXME:
// - Switch to planar bounds detection.  Below is so convoluted it's scary.
//   There is also an edge case of injection the photon, at its first step,
//   has directon cosignZ == 1.  If diffuse reflection occurs and it leaves
//   the medium boundary on its first run, this case is not accounted for.
//   THIS IS BAD WHEN INDEXING INTO A DISPLACEMENT ARRAY AND HAVING A NEGATIVE
//   VALUE CAUSES THINGS TO BLOW UP AND SPEND AN ENTIRE DAY FINDING A SUBTLE BUG.
//   THIS IS A REMINDER TO SWITCH TO PLANE INTERSECTION DETECTION....
// Note: We take the absolute value in the case where the direction
//       cosine is negative, since the distance to boundary would be
//       negative otherwise, which is untrue.  This arises due to assuming
//       the lower axis bound in each dimension (x, y, z) begins at zero.
//       This could also be achieved by simply subtracting the current location
//       from zero (e.g. 0-y/diry), which would change the sign as well.
bool Photon::hitMediumBoundary(void)
{
	double distance_to_boundary = 0.0;
	double distance_to_boundary_X = 0.0;
	double distance_to_boundary_Y = 0.0;
	double distance_to_boundary_Z = 0.0;

	//Layer *layer = m_medium->getLayerFromDepth(currLocation->location.z);
	double mu_t = currLayer->getTotalAttenuationCoeff(currLocation);
	double x_bound = m_medium->Get_X_bound();
	double y_bound = m_medium->Get_Y_bound();
	double z_bound = m_medium->Get_Z_bound();

	// Check if the photon has been given a step size past the outer edges of the medium.
	// If so we calculate the distance to the boundary.

	// The step that puts the photon in a new location based on its current location.
	double x_step = step*currLocation->getDirX() + currLocation->location.x;
	double y_step = step*currLocation->getDirY() + currLocation->location.y;
	double z_step = step*currLocation->getDirZ() + currLocation->location.z;




	if (x_step >= x_bound || x_step <= 0.0f)
	{
		hit_x_bound = true;
		if (currLocation->getDirX() > 0.0f) // Moving towards positive x_bound
			distance_to_boundary_X = (x_bound - currLocation->location.x) / currLocation->getDirX();
		else
			distance_to_boundary_X = abs(currLocation->location.x / currLocation->getDirX());
	}

	if (y_step >= y_bound || y_step <= 0.0f)
	{
		hit_y_bound = true;
		if (currLocation->getDirY() > 0.0f) // Moving towards positive y_bound
			distance_to_boundary_Y = (y_bound - currLocation->location.y) / currLocation->getDirY();
		else
			distance_to_boundary_Y = abs(currLocation->location.y / currLocation->getDirY());
	}

	if (z_step >= z_bound || z_step <= 0.0f)
	{
		hit_z_bound = true;
		if (currLocation->getDirZ() > 0.0f) // Moving towards positive z_bound
			distance_to_boundary_Z = (z_bound - currLocation->location.z) / currLocation->getDirZ();
		else
			distance_to_boundary_Z = abs(currLocation->location.z / currLocation->getDirZ());
	}

	// If none were hit, we can return false and forego any further checking.
	if (!hit_x_bound && !hit_y_bound && !hit_z_bound)
		return false;


	// Need to take care of the rare case that photons can be at a corner edge and step
	// past the boundary of the medium in 2 or more dimensions.  If it is found that
	// the photon has passed through the bounds of the grid in more than 2 bounds we always
	// want to move it the smallest distance (i.e. to the closest bounds), therefore we should
	// never cross 2 bounds (or more!) simultaneously.
	if (hit_x_bound && hit_y_bound)
	{
		if (distance_to_boundary_X < distance_to_boundary_Y)
		{
			distance_to_boundary = distance_to_boundary_X;
			hit_x_bound = true;
			hit_y_bound = false;
			hit_z_bound = false;
		}
		else
		{
			distance_to_boundary = distance_to_boundary_Y;
			hit_x_bound = false;
			hit_y_bound = true;
			hit_z_bound = false;
		}
	}
	else if (hit_x_bound && hit_z_bound)
	{
		if (distance_to_boundary_X < distance_to_boundary_Z)
		{
			distance_to_boundary = distance_to_boundary_X;
			hit_x_bound = true;
			hit_y_bound = false;
			hit_z_bound = false;
		}
		else
		{
			distance_to_boundary = distance_to_boundary_Z;
			hit_x_bound = false;
			hit_y_bound = false;
			hit_z_bound = true;
		}
	}
	else if (hit_y_bound && hit_z_bound)
	{
		if (distance_to_boundary_Y < distance_to_boundary_Z)
		{
			distance_to_boundary = distance_to_boundary_Y;
			hit_x_bound = false;
			hit_y_bound = true;
			hit_z_bound = false;
		}
		else
		{
			distance_to_boundary = distance_to_boundary_Z;
			hit_x_bound = false;
			hit_y_bound = false;
			hit_z_bound = true;
		}

	}
	else if (hit_x_bound && hit_y_bound && hit_z_bound)
	{
		cout << "ERROR: hit 3 boundaries\n";
		cout.flush();
		distance_to_boundary = distance_to_boundary_X < distance_to_boundary_Y ?
				distance_to_boundary_X : distance_to_boundary_Y;
		distance_to_boundary = distance_to_boundary < distance_to_boundary_Z ?
				distance_to_boundary : distance_to_boundary_Z;
		hit_x_bound = true;
		hit_y_bound = hit_z_bound = false;
	}
	else
	{
		// If we make it here, we have hit only one boundary, which means the other 2 have values
		// of zero for their distances.  Therefore, without checking we can simply add them all together
		// to get the single distance to boundary.
		distance_to_boundary = distance_to_boundary_X + distance_to_boundary_Y + distance_to_boundary_Z;
	}



	// If the step size of the photon is larger than the distance
	// to the boundary and we are moving in some direction along
	// the axis we calculate the left over step size and then step
	// the photon to the boundary with the next call to 'hop()'.
	//if (step > distance_to_boundary && distance_to_boundary != 0)
	if (step > distance_to_boundary)
	{
		step_remainder = (step - distance_to_boundary)*mu_t;
		step = distance_to_boundary;
		return true;
	}
	else
	{
		return false;
	}
}


bool Photon::hitLayerBoundary(void)
{
	// 1)Find distance to layer boundary where there could potentially be
	//   refractive index mismatches.
	// 2) If distance to boundary is less than the current step_size of the
	//   photon, we move the photon to the boundary and keep track of how
	//   much distance is left over from step_size - distance_to_boundary.

	double distance_to_boundary = 0.0;
	//Layer *layer = m_medium->getLayerFromDepth(currLocation->location.z);
	double mu_t = currLayer->getTotalAttenuationCoeff(currLocation);



	// If the direction the photon is traveling is towards the deeper boundary
	// of the layer, we execute the first clause.  Otherwise we are moving to
	// the more superficial boundary of the layer.
	if (currLocation->getDirZ() > 0.0)
	{
		distance_to_boundary = (currLayer->getDepthEnd() - currLocation->location.z) / currLocation->getDirZ();
	}
	else if (currLocation->getDirZ() < 0.0)
	{
		distance_to_boundary = (currLayer->getDepthStart() - currLocation->location.z) / currLocation->getDirZ();
	}


	// If the step size of the photon is larger than the distance
	// to the boundary and we are moving in some direction along
	// the z-axis (i.e. not parallel to the layer boundary) we calculate
	// the left over step size and then step the photon to the boundary.
	if (currLocation->getDirZ() != 0.0 && step > distance_to_boundary)
	{
		step_remainder = (step - distance_to_boundary)*mu_t;
		step = distance_to_boundary;
		return true;
	}
	else
	{
		return false;
	}
}


void Photon::specularReflectance(double n1, double n2)
{
	// update the weight after specular reflectance has occurred.
	weight = weight - (pow((n1 - n2), 2) / pow((n1 + n2), 2)) * weight;
}


// Update the direction cosine when internal reflection occurs on z-axis.
void Photon::internallyReflectZ(void)
{
	currLocation->setDirZ(-1*currLocation->getDirZ());

	// Reset the flag.
	hit_z_bound = false;
}

// Update the direction cosine when internal reflection occurs on y-axis.
void Photon::internallyReflectY(void)
{
	currLocation->setDirY(-1*currLocation->getDirY());

	// Reset the flag.
	hit_y_bound = false;
}


// Update the direction cosine when internal reflection occurs on z-axis.
void Photon::internallyReflectX(void)
{
	currLocation->setDirX(-1*currLocation->getDirX());

	// Reset the flag.
	hit_x_bound = false;
}



double Photon::getRandNum(void)
{
	// Thread safe RNG.
	return RNG_generator->getRandNum();

}


