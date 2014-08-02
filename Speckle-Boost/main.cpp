/// Simulates a CCD camera and forms speckle patterns by transforming photons into waves and interfering them on the virtual CCD.
#include <cstdlib>
#include <sys/stat.h> 
#include <boost/filesystem.hpp>
using namespace boost::filesystem;
#include <boost/thread/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <ctime>
#include <iostream>
#include <complex>
#include <algorithm>    
#include <vector>

using namespace std;
#include "CCDGrid.h"

/// Holds the filename and timestamp of the file.
/// Used for sorting and loading data.
struct filename_tstamp
{
	std::string filename;
	size_t tstamp;
};


bool SortFunction (struct filename_tstamp a, struct filename_tstamp b) { return (a.tstamp < b.tstamp); };
void printFiles(std::vector<filename_tstamp> files);
int  Get_num_detected_photons(std::string &filename);

const std::string separator = "---------------------------------------------------\n";

int main()
{


	

	// Let boost decide how many threads to run on this architecture.
	//
    //const int NUM_THREADS = boost::thread::hardware_concurrency();
    const int NUM_THREADS = 2;

	/// Pixel count of the CCD.
	const size_t PIXEL_CNT = 512;


    // Display the possible number of "concurrent" processes possible on this architecture.
	//
    cout << "CPU threads supported: = " << NUM_THREADS << endl;
	cout << separator;
    
    path p_exit_data = "../Data/Detected_photons";
	path p_speckle_data = "../Data/Speckles";
    if (is_directory(p_exit_data))
	{
		cout << "Data directory found: " << p_exit_data << '\n';
		
		/// Check if we have a directory to store the generated speckle data.
		if (is_directory(p_speckle_data))
		{
			cout << "Storing speckle data: " << p_speckle_data << '\n';
		}
		else
		{
			cout << "!!!ERROR: Directory for storing speckle data to location [" << p_speckle_data << "] does not exist.\n";
			exit(1);
		}
		
	}
	else
	{
		cout << "!!!ERROR: Data directory does not exist.  Given the following path: " << p_exit_data << '\n';
		exit(1);
	}



	/// Get the timestamps of all the files and add them to the vector for sorting later.
	struct stat st;
	std::vector<filename_tstamp> files;
	for (directory_iterator itr(p_exit_data); itr!=directory_iterator(); ++itr)
    {
		std::string f = itr->path().string(); // + itr->path().filename().string();
		if (stat(f.c_str(), &st) != 0)
		{
			cout << "!!!ERROR: Unable to read time stamp of " << f << '\n';
			cout << "st.st_mtime = " << st.st_mtime << '\n';
			exit(1);
		}

		/// Ignore the seeds file used for generating exit photons (i.e. through exit aperture) and directories.
		if ((itr->path().filename().string() != "seeds_for_exit.dat") && is_regular_file(itr->path()))
		{
			filename_tstamp temp;
			temp.filename = f;
			temp.tstamp = st.st_mtime;
			files.push_back(temp);
		}
    }

	//printFiles(files);
	/// Sort the files based on their timestamp.
	std::sort (files.begin(), files.end(), SortFunction);
	//printFiles(files);

	
    CCDGrid *ccd[NUM_THREADS];
    boost::thread_group threads;
    
	// The number of exit-data files to read in and operate on.
	//
    const int NUM_FILES = files.size();
	int num_detected_photons = Get_num_detected_photons((files.at(0)).filename);
	cout << separator;
	cout << "Processing " << NUM_FILES << " exit data files.\n";
	cout << "Detected photons: " << num_detected_photons << '\n';
	cout << separator;


	

	// Capture the time before interfering photons onto the CCD.
    //
    clock_t start, end;
    start = clock();

	
	/// Allocate a CCD for each thread, allowing multiple sets of data (time-steps)
	/// to be constructed in parallel.
	///
	for (size_t i = 0; i < NUM_THREADS; i++)
	{
		//CCDGrid *grid = new CCDGrid(128, 128, 10e-6);
        ccd[i] = new CCDGrid(PIXEL_CNT,
							   PIXEL_CNT,
							   6e-6,
							   num_detected_photons);
	}

    /// Track how many times the threads have run, which allows displaying of which
    /// file (i.e. time step of ultrasound propagation) of exit data has been opened.
    size_t thread_cnt = 0;

	// Spawn threads to open and process the data contained in the
	// exit aperture files.
	//
    for (size_t i = thread_cnt; i <= NUM_FILES-1; i += NUM_THREADS)
	{
		
		// Launch threads.
        for (size_t j = 0; j < NUM_THREADS; j++)
		{ 
            
        	// The file name that holds the exit aperture data.  Within the loop below the time step
        	// is appended to this name so each thread operates on a separate data file.
        	//
        	//std::string filename = "exit-aperture-" + boost::lexical_cast<std::string>(i+j) + ".dat";
			std::string exit_data_filename = (files.at(i+j)).filename;            
			std::string speckle_data_filename = p_speckle_data.string() + "/speckle_t" + boost::lexical_cast<std::string>(i+j) + ".dat";

            /// Load the data into the virtual ccd for processing before the thread is spawned.
            ///
            ccd[j]->Load_exit_data(exit_data_filename);

        	cout << "Launching thread (" << thread_cnt++ << " of " << NUM_FILES << ").  Processing " << exit_data_filename << endl;
            boost::thread *t = new boost::thread(&CCDGrid::makeSpeckle,
                                                 ccd[j],
                                                 exit_data_filename,
                                                 speckle_data_filename);
			threads.add_thread(t);
			
		}

		
        // Join threads.
		threads.join_all();
        
        // Clean up memory.
        for (int k = 0; k < NUM_THREADS; k++)
        {
            //delete ccd[k];
            (ccd[k])->Zero_grid();
        }
        
	}


	/// Free CCD grid memory.
	///
	for (size_t i = 0; i < NUM_THREADS; i++)
	{
        if (ccd[i] != NULL)
		{		
            delete ccd[i];
            ccd[i] = NULL;
		}
	}


	// Print out the elapsed time it took from beginning to end.
    //
    end = (((double)clock() - start) / CLOCKS_PER_SEC)*10;
    cout << "\n\nTotal time elapsed: " << end << " seconds" << endl;
     
     
	return 0;
}



void printFiles(std::vector<filename_tstamp> f)
{
	std::vector<filename_tstamp>::iterator it = f.begin();
	for(; it != f.end(); ++it) 
	{
    	cout << it->filename << '\n';
	}



}


int Get_num_detected_photons(std::string &filename)
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

	return i;
}



