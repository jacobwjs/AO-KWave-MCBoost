//
//  logger.cpp
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 BMPI. All rights reserved.
//

#include "vector3D.h"
#include "photon.h"
#include "logger.h"
#include <cmath>
using std::cos;



Logger * Logger::pInstance = 0;



Logger::Logger()
{
    exit_cnt = seed_cnt = 0;
}


Logger::~Logger()
{
    exit_data_stream.close();
    absorber_data_stream.close();
    rng_seed_stream.close();
    tof_stream.close();
}

Logger * Logger::getInstance(void)
{
    if (!pInstance)
    {
        pInstance = new Logger();		
   }
    
    return pInstance;
}


void Logger::openExitFile(const std::string &filename)
{
    exit_cnt = 0;
    // Ensure file stream is not already open.
    if (exit_data_stream.is_open())
        exit_data_stream.close();
    
    exit_data_stream.open(filename.c_str());
    if (!exit_data_stream)
    {
    	cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
    	exit(1);
    }
}


void Logger::openTOFFile(const std::string &filename)
{
	// Ensure file stream is not already open.
	if (tof_stream.is_open())
		tof_stream.close();

	tof_stream.open(filename.c_str());
	if (!tof_stream)
	{
		cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
		exit(1);
	}
}


void Logger::Open_vel_disp_file(const std::string &filename)
{
	// Ensure file stream is not already open.
	if (velocity_displacement_stream.is_open())
		velocity_displacement_stream.close();

	velocity_displacement_stream.open(filename.c_str());
	if (!velocity_displacement_stream)
	{
		cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
		exit(1);
	}
}




void Logger::createRNGSeedFile(const std::string &filename)
{
    seed_cnt = 0;
    rng_seed_stream.open(filename.c_str());
    if (!rng_seed_stream)
    {
        cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
        exit(1);
    }
}

void Logger::openAbsorberFile(const std::string &filename)
{
    // Ensure file stream is not already open.
    if (absorber_data_stream.is_open())
        absorber_data_stream.close();
    
    absorber_data_stream.open(filename.c_str());
    if (!absorber_data_stream)
    {
    	cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
    	exit(1);
    }
}




void Logger::Write_weight_OPLs_coords(Photon &p)
{
	boost::mutex::scoped_lock lock(m_mutex);

    exit_cnt++;
	exit_data_stream << p.weight << " "
					 << p.displaced_optical_path_length << " "
					 << p.refractiveIndex_optical_path_length << " "
					 << p.combined_OPL << " "
					 << p.currLocation->location.x << " "
					 << p.currLocation->location.y << " " 
					 << p.currLocation->location.z << " "
					 << "\n";
	exit_data_stream.flush();
}
                                          
                                  

void Logger::writeRNGSeeds(const unsigned int s1, const unsigned int s2,
							const unsigned int s3, const unsigned int s4)
{
    boost::mutex::scoped_lock lock(m_mutex);

    seed_cnt++;
    rng_seed_stream << s1 << " " <<
                       s2 << " " <<
                       s3 << " " <<
                       s4 << " " << "\n";
    rng_seed_stream.flush();
}


void Logger::writeTOFData(const double tof)
{

	boost::mutex::scoped_lock lock(m_tof_mutex);

	tof_stream << tof << " \n";
	tof_stream.flush();
}



/// For 3-D axis.
void Logger::Write_velocity_displacement(float ux, float uy, float uz,
                                         float disp_x, float disp_y, float disp_z)
{
    boost::mutex::scoped_lock lock(m_mutex);

    velocity_displacement_stream << ux << " "
                                 << uy << " "
                                 << uz << " "
                                 << disp_x << " "
                                 << disp_y << " "
                                 << disp_z << "\n";
    velocity_displacement_stream.flush();
    
}




// Returns the current time.
std::string
Logger::getCurrTime(void)
{
    
	// Set current time variable to be used with naming data files that are saved from the simulations.
	epoch = time(NULL);
	ptr_ts = localtime(&epoch);
    
	return (boost::lexical_cast<std::string>(ptr_ts->tm_hour) + "_" +
			boost::lexical_cast<std::string>(ptr_ts->tm_min) + "_" +
			boost::lexical_cast<std::string>(ptr_ts->tm_sec));
}

 


