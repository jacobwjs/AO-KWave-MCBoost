//
//  LoggerBase.cpp
//  Xcode
//
//  Created by Jacob Staley on 7/18/14.
//  Copyright 2014 BMPI, University of Twente. All rights reserved.
//

#include "multikey.h"
#include "vector3D.h"
#include "photon.h"
#include "LoggerBase.h"
#include <cmath>
using std::cos;




LoggerBase::LoggerBase()
{
    exit_cnt = seed_cnt = 0;
    startClock();
}


void LoggerBase::Destroy()
{
    endClock();
    this->~LoggerBase();
}

void LoggerBase::Destroy(const std::string &detector_name)
{
    //Write_weight_OPLs_coordinates(detector_name);
}

LoggerBase::~LoggerBase()
{
    if (exit_data_stream.is_open())
        exit_data_stream.close();
    if (absorber_data_stream.is_open())
        absorber_data_stream.close();
    if (rng_seed_stream.is_open())
        rng_seed_stream.close();
    if (tof_stream.is_open())
        tof_stream.close();
    if (velocity_displacement_stream.is_open())
        velocity_displacement_stream.close();
    if (modulation_depth_stream.is_open())
        modulation_depth_stream.close();
    
}


void LoggerBase::openExitFile(const std::string &filename)
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


void LoggerBase::openTOFFile(const std::string &filename)
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


void LoggerBase::Open_vel_disp_file(const std::string &filename)
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


void LoggerBase::Open_modulation_depth_file(const std::string &filename)
{
    if (modulation_depth_stream.is_open())
        modulation_depth_stream.close();
    
    modulation_depth_stream.open(filename.c_str());
    if (!modulation_depth_stream)
    {
        cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
        exit(1);
    }
}



void LoggerBase::createRNGSeedFile(const std::string &filename)
{
    seed_cnt = 0;
    rng_seed_stream.open(filename.c_str());
    if (!rng_seed_stream)
    {
        cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
        exit(1);
    }
}

void LoggerBase::openAbsorberFile(const std::string &filename)
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




void LoggerBase::Write_weight_OPLs_coords(Photon &p)
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



void LoggerBase::writeRNGSeeds(const unsigned int s1, const unsigned int s2,
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


void LoggerBase::writeTOFData(const double tof)
{
    
	boost::mutex::scoped_lock lock(m_tof_mutex);
    
	tof_stream << tof << " \n";
	tof_stream.flush();
}



/// For 3-D axis.
void LoggerBase::Write_velocity_displacement(float ux, float uy, float uz,
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



/// Store the OPL based on the initial seeds of the photon.
void LoggerBase::Store_weight_OPLs_coordinates(RNGSeeds &seeds, Photon &p)
{
    boost::mutex::scoped_lock lock(m_mutex);
    exit_cnt++;
    
    /// Create a new key for the map based on this detected photon's initial seeds.
    MultiKey key(seeds.s1, seeds.s2, seeds.s3, seeds.s4);
    
    //OPL new_vals;
    exitInfo new_vals;
    new_vals.weight                         = p.weight;
    new_vals.refractive_index_contribution  = p.displaced_optical_path_length;
    new_vals.displacement_contribution      = p.refractiveIndex_optical_path_length;
    new_vals.combined_contribution          = p.combined_OPL;
    new_vals.exit_x_coord                   = p.currLocation->location.x;
    new_vals.exit_y_coord                   = p.currLocation->location.y;
    new_vals.exit_z_coord                   = p.currLocation->location.z;
    
    
    /// Check if the key already exists.
    if (Exit_Map.count(key) != 0)
    {
        /// get the vector of OPLs for this key (i.e. seed values)
        std::vector<exitInfo> &temp = Exit_Map[key];
        
        /// and store it back
        temp.push_back(new_vals);
        Exit_Map[key].swap(temp);
    }
    else
    {
        /// Insert into the map.
        std::vector<exitInfo> temp;
        
        temp.push_back(new_vals);
        Exit_Map[key] = temp;
    }
    
    
}


void LoggerBase::Write_weight_OPLs_coordinates_from_MAP()
{

    std::string output_file = "./Data/Detected_photons/" + logger_name + "-" + getCurrTime() + "_exit_data.dat";
    openExitFile(output_file);
    
    cout << " Detected" << " (" << logger_name << "): " << Exit_Map.size() << " photons\n";
    
    if (!exit_data_stream)
    {
        cout << "!!! ERROR: Output stream is not open. LoggerBase::Write_weight_OPLs_coords_from_MAP()\n";
        exit(1);
    }
    
    
    std::map<MultiKey, std::vector<exitInfo> >::const_iterator map_iter;
    for (map_iter = Exit_Map.begin(); map_iter != Exit_Map.end(); map_iter++)
    {
        /// NOTE:
        /// - map_iter->first  = key
        /// - map_iter->second = value
        std::vector<exitInfo> temp = map_iter->second;
        for (std::vector<exitInfo>::const_iterator vec_iter = temp.begin(); vec_iter != temp.end(); vec_iter++)
        {
            exit_data_stream << std::fixed << std::setprecision(15)
                             << (*vec_iter).weight << " "
                             << (*vec_iter).refractive_index_contribution << " "
                             << (*vec_iter).displacement_contribution << " "
                             << (*vec_iter).combined_contribution << " "
                             << (*vec_iter).exit_x_coord << " "
                             << (*vec_iter).exit_y_coord << " "
                             << (*vec_iter).exit_z_coord;
        }
        exit_data_stream << endl;
    }
    
    exit_data_stream.flush();
    
    /// Once we write out the MAP we want to clear it.
    Exit_Map.clear();
}




/// Store the OPL based on the initial seeds of the photon.
void LoggerBase::Store_OPL(RNGSeeds &seeds, double n_OPL, double d_OPL)
{
    boost::mutex::scoped_lock lock(m_mutex);
    
    /// Create a new key for the map based on this detected photon's initial seeds.
    MultiKey key(seeds.s1, seeds.s2, seeds.s3, seeds.s4);
    
    //OPL new_vals;
    opticalPathLengths new_vals;
    new_vals.refractive_index_contribution  = n_OPL;
    new_vals.displacement_contribution      = d_OPL;
    new_vals.combined_contribution          = 0.0f;
    
    /// Check if the key already exists.
    if (OPL_Map.count(key) != 0)
    {
        /// get the vector of OPLs for this key (i.e. seed values)
        std::vector<opticalPathLengths> &temp = OPL_Map[key];
        
        /// and store it back
        temp.push_back(new_vals);
        OPL_Map[key].swap(temp);
    }
    else
    {
        /// Insert into the map.
        std::vector<opticalPathLengths> temp;
        
        temp.push_back(new_vals);
        OPL_Map[key] = temp;
    }
    
    
}


void LoggerBase::Write_OPL_data()
{
    if (!modulation_depth_stream)
    {
        cout << "!!! ERROR: Output stream is not open. LoggerBase::Write_OPL_data()\n";
        exit(1);
    }
    
    
    std::map<MultiKey, std::vector<opticalPathLengths> >::const_iterator map_iter;
    for (map_iter = OPL_Map.begin(); map_iter != OPL_Map.end(); map_iter++)
    {
        /// NOTE:
        /// - map_iter->first  = key
        /// - map_iter->second = value
        std::vector<opticalPathLengths> temp = map_iter->second;
        for (std::vector<opticalPathLengths>::const_iterator vec_iter = temp.begin(); vec_iter != temp.end(); vec_iter++)
        {
            modulation_depth_stream << std::fixed << std::setprecision(15)
            << (*vec_iter).refractive_index_contribution<< ','
            << (*vec_iter).displacement_contribution << ' ';
        }
        modulation_depth_stream << endl;
    }
    
    modulation_depth_stream.flush();
}



// Returns the current time.
std::string
LoggerBase::getCurrTime(void)
{
    
	// Set current time variable to be used with naming data files that are saved from the simulations.
	epoch = time(NULL);
	ptr_ts = localtime(&epoch);
    
    return (boost::lexical_cast<std::string>(ptr_ts->tm_yday) + "_" +
            boost::lexical_cast<std::string>(ptr_ts->tm_hour) + "_" +
			boost::lexical_cast<std::string>(ptr_ts->tm_min) + "_" +
			boost::lexical_cast<std::string>(ptr_ts->tm_sec));
}




