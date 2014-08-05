/// Simulates a CCD camera and forms speckle patterns by transforming photons into waves and interfering them on the virtual CCD.
#include <cstdlib>
#include <sys/stat.h> 
//#include <boost/filesystem/operations.hpp>
//#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;
#include <boost/thread/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <ctime>
#include <iostream>
#include <complex>
#include <algorithm>    
#include <vector>
#include <glob.h>

using namespace std;
#include "CCDGrid.h"

#include "CommandLineParametersSpeckleBoost.h"

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
void Print_CCD_attributes(const int pix_cnt_x,
                          const int pix_cnt_y,
                          const double pix_size,
                          const double center_x,
                          const double center_y,
                          const std::string mechanism);


/// Work around for compiler errors with boost
vector<string> globVector(const string& pattern){
    glob_t glob_result;
    glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> files;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        files.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return files;
}


const std::string separator = "---------------------------------------------------\n";

int main(int argc, char** argv)
{

    TCommandLineParametersSpeckleBoost CommandLineParams;
	CommandLineParams.ParseCommandLine(argc, argv);
    

	// Let boost decide how many threads to run on this architecture.
	//
    //const int NUM_THREADS = boost::thread::hardware_concurrency();
    const size_t NUM_THREADS = CommandLineParams.GetNumberOfThreads();

    /// What AO mechanism we are using to form the speckle pattern.
    std::string mechanism;
    if (CommandLineParams.IsCompute_displacement_OPL()) mechanism = "d_OPL";
    if (CommandLineParams.IsCompute_refractive_OPL())   mechanism = "n_OPL";
    if (CommandLineParams.IsCompute_combined_OPL())     mechanism = "combined_OPL";
    
	/// Pixel attributes of the CCD.
    double center_ccd_x_coord = CommandLineParams.GetCCDxcoord();
    double center_ccd_y_coord = CommandLineParams.GetCCDycoord();
	size_t tmp_CNT_X = CommandLineParams.GetNumberOfPixelsXdim();
    size_t tmp_CNT_Y = CommandLineParams.GetNumberOfPixelsYdim();
    float tmp_SIZE   = CommandLineParams.GetPixelSize();

    /// If nothing has been specified via the commandline, supply defaults for CCD attributes.
    const size_t PIXEL_CNT_X = (tmp_CNT_X == 0) ? 512:tmp_CNT_X;
    const size_t PIXEL_CNT_Y = (tmp_CNT_Y == 0) ? 512:tmp_CNT_Y;
    const double PIXEL_SIZE  = (tmp_SIZE  == 0.0f) ? 6e-6:tmp_SIZE;
    


    // Display the possible number of "concurrent" processes possible on this architecture.
	//
    cout << "CPU threads supported: = " << NUM_THREADS << endl;
	cout << separator;
    
    std::string exit_data_dir = CommandLineParams.GetInputDirectoryName(); //"../Data/Detected_photons";
    path p_exit_data = exit_data_dir;
	path p_speckle_data = CommandLineParams.GetOutputDirectoryName(); //"../Data/Speckles";
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
			cout << "Directory for storing speckle data to location [" << p_speckle_data << "] does not exist.\n"
                 << "Creating directory\n";
            
            create_directories(p_speckle_data); 
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
    exit_data_dir+= "*";
    std::vector<std::string> filenames = globVector(exit_data_dir);
    for (std::vector<std::string>::iterator it = filenames.begin(); it != filenames.end(); it++)
    {
		if (stat(it->c_str(), &st) != 0)
		{
			cout << "!!!ERROR: Unable to read time stamp of " << *it << '\n';
			cout << "st.st_mtime = " << st.st_mtime << '\n';
			exit(EXIT_FAILURE);
		}
        
        /// Remove any directories, only interested in files.
        path temp_path = *it;
        if (is_regular_file(temp_path))
        {
            filename_tstamp temp;
            temp.filename = *it;
            temp.tstamp = st.st_mtime;
            files.push_back(temp);
        }
       
    }
    
    
    /// THIS IS CAUSING AN ERROR DUE TO COMPILER VERSIONS BEING DIFFERENT FOR BUILDING BOOST LIBRARIES AND BUILDING THIS EXECUTABLE
//	for (directory_iterator itr(p_exit_data); itr!=directory_iterator(); ++itr)
//    {
//		std::string f = itr->path().string(); // + itr->path().filename().string();
//		if (stat(f.c_str(), &st) != 0)
//		{
//			cout << "!!!ERROR: Unable to read time stamp of " << f << '\n';
//			cout << "st.st_mtime = " << st.st_mtime << '\n';
//			exit(1);
//		}
//
//		/// Ignore the seeds file used for generating exit photons (i.e. through exit aperture) and directories.
//		if ((itr->path().filename().string() != "seeds_for_exit.dat") && is_regular_file(itr->path()))
//		{
//			filename_tstamp temp;
//			temp.filename = f;
//			temp.tstamp = st.st_mtime;
//			files.push_back(temp);
//		}
//    }

	//printFiles(files);
	/// Sort the files based on their timestamp.
	std::sort (files.begin(), files.end(), SortFunction);
	//printFiles(files);

	
    CCDGrid *ccd[NUM_THREADS];
    boost::thread_group threads;
    
	// The number of exit-data files to read in and operate on.
	//
    const size_t NUM_FILES = files.size();
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
		Print_CCD_attributes(PIXEL_CNT_X,
                             PIXEL_CNT_Y,
                             PIXEL_SIZE,
                             center_ccd_x_coord,
                             center_ccd_y_coord,
                             mechanism);
        
        ccd[i] = new CCDGrid(PIXEL_CNT_X,
                             PIXEL_CNT_Y,
                             PIXEL_SIZE,
                             center_ccd_x_coord,
                             center_ccd_y_coord,
                             mechanism,
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
			std::string exit_data_filename = (files.at(i+j)).filename;
            
            /// Prepend 'SPECKLE' to the name of the exit data that is being processed so we can later match any speckle data file with it's exit data file.
			std::string speckle_data_filename = p_speckle_data.string() + "SPECKLE_" + exit_data_filename.substr(exit_data_filename.find_last_of("/\\")+1);
            

            /// Load the data into the virtual ccd for processing before the thread is spawned.
            ///
            ccd[j]->Load_exit_data(exit_data_filename);

        	cout << "Launching thread (" << ++thread_cnt << " of " << NUM_FILES << ").  Processing " << exit_data_filename << endl;
            boost::thread *t = new boost::thread(&CCDGrid::makeSpeckle,
                                                 ccd[j],
                                                 exit_data_filename,
                                                 speckle_data_filename);
			threads.add_thread(t);
			
		}

		
        // Join threads.
		threads.join_all();
        
        // Clean up memory.
        for (size_t k = 0; k < NUM_THREADS; k++)
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



void Print_CCD_attributes(const int pix_cnt_x,
                          const int pix_cnt_y,
                          const double pix_size,
                          const double center_x,
                          const double center_y,
                          const std::string mechanism)
{
    cout << separator
         << "CCD attributes /" << endl
         << "---------------\n"
         << " - pixels: " << pix_cnt_x << "x" << pix_cnt_y << endl
         << " - pixel size: " << pix_size << " [meters]" << endl
         << " - dimensions: " << pix_size*pix_cnt_x << "x" << pix_size*pix_cnt_y << " [meters]" << endl
         << " - coordinates (x,y): " << "(" << center_x << "," << center_y << ") [meters]" << endl
         << " - detecting AO mechanism: " << mechanism << endl
         << separator << endl;
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



