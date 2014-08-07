#include "ExitData.h"
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <iostream>
using std::cout;



ExitData::ExitData()
{
    
}


ExitData::ExitData(const int num_detected_photons)
{
	values.reserve(num_detected_photons);
}


ExitData::~ExitData()
{

}


/// The number of lines written to the exit data file is the number
/// of photons that were detected.
size_t ExitData::Get_num_exit_data_lines(const std::string &filename)
{
    
	int i = 0;
	std::string line;
    
	// Input stream.
	std::ifstream temp_stream;
	temp_stream.open(filename.c_str());
    
	do
	{
    	getline(temp_stream,line);
		if (temp_stream.fail())
		{
			break;
		}
 		++i;
	}
	while (temp_stream.good());
    
	temp_stream.close();
    
    
    /// Update the internal value.    
    m_num_detected_photons = i;
    
	return i;
}

size_t ExitData::Get_num_exit_data_entries_per_line()
{
    /// To know how large to make each vector we need to know how many data points
	/// were written out to file while collecting data on exit photons.  Data is
	/// written line by line, and based on what is chosen to be saved in the simulation
	/// (e.g. weight, coords, etc.) this value changes.  Here we read in one line and
	/// find out how many data points are on a single line.
	std::istringstream stream1;
	std::string line;
	getline(exit_file_stream, line);
	stream1.str(line);
	double temp_num;
	size_t cnt = 0;
	while (stream1 >> temp_num) cnt++;
    
    /// Reset the ifstream back to the beginning of the file.
    exit_file_stream.clear();
    exit_file_stream.seekg(0, std::ios::beg);
    
    /// Update the internal value.
    m_num_exit_data_entries_per_line = cnt;
    
    return cnt;
}


// Allocate a 2D vector to hold all of the exit aperture data.
void ExitData::loadExitData(const std::string &filename)
{
	
    values.reserve(Get_num_exit_data_lines(filename));
    
    // Open the file that contains the exit data from the medium's aperture.
	//
	if (exit_file_stream.is_open())
	{
        exit_file_stream.close();
    }
   
	exit_file_stream.open(filename.c_str());
    if (!exit_file_stream)
    {
        cout << "!!! ERROR: Could not open file (" << filename << ")\n";
        exit(1);
    }

	
    /// Find out how many columns of data were written to the exit data file.
	size_t COLS = Get_num_exit_data_entries_per_line();
    

	/// Capacity has been reserved upon creation, but to remove any issues
	/// with old data we clear the vector before we populate it.
	values.clear();
	assert(values.size() == 0);
	
    

    // Read and store the exit data to the 2D array.
    double temp = 0.0;
	size_t i = 0;
	do
	{	 
	        values.push_back(std::vector<double>(COLS));
	        
	        for (size_t j = 0; j < COLS; j++)
	        {
	            exit_file_stream >> temp;
				if (exit_file_stream.fail())
				{
					break;
				}

	            values[i][j] = temp; 
	        }
 		++i;
	}
	while (exit_file_stream.good());
}
