#ifndef EXIT_DATA_H
#define EXIT_DATA_H

#include <vector>
#include <string>
#include <fstream>
using std::ifstream;

class ExitData {
public:
    ExitData();
	ExitData(const int num_detected_photons);
	~ExitData();

    
	
	// Load in the data from the exit file.
	void loadExitData(const std::string &filename);
    
    size_t  Get_num_exit_data_lines(const std::string &filename);
    
    size_t  Get_num_exit_data_entries_per_line();
    
    size_t  getNumPhotonAttributes()    const {return m_num_exit_data_entries_per_line;};
    
    size_t  getNumPhotons()             const {return m_num_detected_photons;};


   	std::vector<std::vector<double> > values;
private:
    size_t m_num_detected_photons;
	
    size_t m_num_exit_data_entries_per_line;

	// Input stream.
	std::ifstream exit_file_stream;
};

#endif //EXIT_DATA_H  
