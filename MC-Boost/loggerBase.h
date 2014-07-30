//
//  loggerBase.h
//  Xcode
//
//  Created by Jacob Staleu on 7/18/11.
//  Copyright 2014 BMPI, University of Twente. All rights reserved.
//

#ifndef LOGGERBASE_H
#define LOGGERBASE_H

#include <vector>
#include <map>
#include <ctime>
#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <boost/thread/mutex.hpp>
#include <boost/lexical_cast.hpp>

//#include "RNG.h"
#include "common_structs.h"

// Forward decleration of objects.
class Photon;
class Vector3d;



class LoggerBase
{
public:
    
    LoggerBase();
    virtual ~LoggerBase();
    
    virtual void Destroy();
    virtual void Destroy(const std::string &filename);
    
    virtual void openExitFile(const std::string &filename);
    virtual void createRNGSeedFile(const std::string &filename);
    virtual void openAbsorberFile(const std::string &filename);
    virtual void openTOFFile(const std::string &filename);
	virtual void Open_vel_disp_file(const std::string &filename);
    virtual void Open_modulation_depth_file(const std::string &filename);
    
    
	/// STUB
    virtual void	Write_absorber_data(const double absorbed_weight) {cout << "Logger::Write_absorber_data ... STUB\n";};
    
	/// Writes the weight, optical path lengths (displaced, refractive changes, combined) and exit coords.
	virtual void	Write_weight_OPLs_coords(Photon &p);                /// Write straight to disk as they are received during execution.
    virtual void    Write_weight_OPLs_seeds_coordinates_from_MAP();     /// Write data that has been stored in the MAP during execution.
    
    /// Stores the weight, optical path lengths (displaced, refractive changes, combined) and exit coords.
    virtual void    Store_weight_OPLs_seeds_coordinates(RNGSeeds &seeds, Photon &p);
    
    
	/// Writes velocity and displacements obtained from k-Wave/MC-Boost.
	virtual void	Write_velocity_displacement(float ux, float uy, float uz,
                                        float disp_x, float disp_y, float disp_z);
    
    /// Stores the OPL of a photon as it exits through the detector.
    virtual void    Store_OPLs(RNGSeeds &seeds, double refractive_OPL, double displacment_OPL, double combined_OPL);
    
    /// Writes all the stored OPL data to disk.
    virtual void    Write_OPL_data(void);
    
    
    // Writes the seed that generated the random events that lead this photon to escape through
    // the aperture.
    //
    virtual void    writeRNGSeeds(const unsigned int s1, const unsigned int s2,
    					  const unsigned int s3, const unsigned int s4);
    
    
    // Returns the number of photons that were detected through the exit-aperture.
    //
    virtual size_t     	Get_num_exited_photons(void) 		{return exit_cnt;};
	virtual size_t 		Get_num_detected_seed_photons(void)	{return seed_cnt;};
    
    
    // Writes the time-of-flight value for the photon bundle when it exits the medium.
    virtual void    writeTOFData(const double tof);
    
    
    // Starts the clock for timing.
    virtual void    startClock() {start = clock();}
    
    // Ends the clock for timing.
    virtual void    endClock() {end = ((double)clock() - start) / CLOCKS_PER_SEC;}
    
    // Returns the current time.
    virtual std::string     getCurrTime(void);
    
    /// Set the logger's name.
    void    Set_name(const std::string &name)   {logger_name = name;};
    
    
    
protected:
    
    // The output streams associated with data for the photon and data for
    // the absorbers.
    ofstream exit_data_stream;              // Photon exit from medium data stream.
    ofstream absorber_data_stream;          // Absorber stream.
    ofstream rng_seed_stream;               // Random-number-generator stream.
    ofstream tof_stream;                    // Time-of-flight stream.
	ofstream velocity_displacement_stream;	// Velocity and displacement stream.
    ofstream modulation_depth_stream;       // OPL stream.
    
    
    // Tracks how many photons were detected through the aperture.
    size_t exit_cnt;
	size_t seed_cnt;
    
    boost::mutex m_mutex;
    boost::mutex m_tof_mutex;
    
    
    /// Map with multiple keys that point to a vector of OPL's to
    /// compare tagged to untagged portions of light.
    MultiKeyMap_OPLs OPL_Map;
    
    /// Map with multiple keys that point to a vector of 'exitInfo' to
    /// eventually write out to disk.
    MultiKeyMap_ExitInfo Exit_Map;
    
    
    // Used to append time/date to saved data files and update execution time.
    time_t epoch;
    struct tm *ptr_ts;
    clock_t start, end;
    
    std::string logger_name;
    
    
private:
    LoggerBase(const LoggerBase&){};             // copy constructor is private
    LoggerBase& operator=(const LoggerBase& rhs) {};  // assignment operator is private
};

#endif  /// LOGGERBASE_H
