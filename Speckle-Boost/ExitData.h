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
    
    int  Get_num_exit_data_entries(const std::string &filename);
    int  getNumPhotons()    const {return m_num_detected_photons;};


   	std::vector<std::vector<double> > values;
private:
    int m_num_detected_photons;
	//std::vector<std::vector<double> > values;

	// Input stream.
	std::ifstream exit_file_stream;
};

#endif //EXIT_DATA_H  
