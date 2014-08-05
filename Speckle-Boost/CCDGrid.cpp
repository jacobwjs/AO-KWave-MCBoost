#include "CCDGrid.h"
#include "ExitData.h"
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <iostream>
using std::cout;
using std::endl;



struct less_than_key
{
    inline bool operator() (const std::vector<double> vec1, const std::vector<double> vec2)
    {
        return (vec1[1] < vec2[1]);
    }
};





CCDGrid::CCDGrid(int x_pixels, 
				 int y_pixels, 
				 double pixel_size,
                 double center_x,
                 double center_y,
                 std::string mechanism)
{
    //initCommon();
    m_center_x = center_x;
    m_center_y = center_y;
    
    m_pixel_size = pixel_size;
    
    m_OPL_mechanism = mechanism;
    
	exit_data = new ExitData();
	setGrid(x_pixels, y_pixels, pixel_size);
}


CCDGrid::~CCDGrid()
{
    //cout << "Destructor: ~CCDGrid()\n";
    if (exit_data)
    {
        delete exit_data;
        exit_data = NULL;
    }
}


void CCDGrid::initCommon(void)
{
    
}

int CCDGrid::Get_num_detected_photons() const
{
    assert(exit_data != NULL);
    return exit_data->getNumPhotons();
}


void CCDGrid::setGrid(int x_pixels, int y_pixels, double pixel_size)
{
	// Set the pixel dimensions.
	dx = dy = pixel_size;
    num_x_pixels = x_pixels;
    num_y_pixels = y_pixels;
    
	// Allocate the 2D array of pixels.
	for (int i = 0; i < x_pixels; i++)
	{
		grid.push_back(std::vector<std::complex<double> >(y_pixels));
	}
    
    Zero_grid();
}


void CCDGrid::Zero_grid(void)
{
   // Ensure the grid is zero'd out.
    for (int x = 0; x < num_x_pixels; x++)
        for (int y = 0; y < num_y_pixels; y++)
            grid[x][y] = 0.0;
}


void CCDGrid::Load_exit_data(const std::string &input_filename)
{

    if (exit_data)
    {
        exit_data->loadExitData(input_filename);
    }
    else
    {
        cout << "!!! Error: 'ExitData' object does not exist !!!\n"
             << "CCDGrid::Load_exit_data()\n";

        exit(EXIT_FAILURE);
    }
}

void CCDGrid::makeSpeckle(const std::string &input_filename,
						  const std::string &output_filename)
{

    
    // Load the exit data, from file, for this thread.
    //
    //exit_data->loadExitData(input_filename);




	/// Sort the exit data based on the OPL.
	//std::sort(exit_data->values.begin(), exit_data->values.end(), less_than_key());
 	//cout << "Sorted...\n";

    
    double distance_to_pixel = 0.0;
    double x_pixel = 0.0;
    double y_pixel = 0.0;
    
    // Starting leftmost location of the CCD.
    //
    double start_x = m_center_x - (num_x_pixels/2 * dx);
    double start_y = m_center_y - (num_y_pixels/2 * dy);

    

    
    // Exit weight of the photon when it left the medium.
    //
    double weight = 0.0;    
    
    // The optical path length used in the calculation.
	double displaced_OPL 		= 0.0;
    double refraction_OPL       = 0.0;
	double combined_OPL			= 0.0;
    double OPL		 			= 0.0;
    
    // Total path length (medium + distance to CCD pixel).
    double L = 0.0;
    
    // Holds the x,y exit locations of the photon from the medium.
    double exit_location_x = 0.0;
    double exit_location_y = 0.0;
    double exit_location_z = 0.0;   


    bool DISPLACEMENT_ENABLED = false;
    bool REFRACTION_ENABLED   = false;
    bool COMBINED_ENABLED     = false;
    if (m_OPL_mechanism == "d_OPL") DISPLACEMENT_ENABLED      = true;   // Displacement AO mechanism
    if (m_OPL_mechanism == "n_OPL") REFRACTION_ENABLED        = true;  // Refraction AO mechanism
    if (m_OPL_mechanism == "combined_OPL") COMBINED_ENABLED   = true; // Displacement and Refraction

    if (DISPLACEMENT_ENABLED)
	{
        cout << " - Displacement: ON\n";
		cout.flush();
	}
    if (REFRACTION_ENABLED)
	{
        cout << " - Refraction: ON\n";
		cout.flush();
	}
    if (COMBINED_ENABLED)
    {
        cout << " - Combined (Refraction + Displacement): ON\n";
        cout.flush();
    }
 
    // For each detected photon.
    //
    for (size_t i = 0; i < exit_data->values.size(); i++)
    {


        // Store this photon's weight and path length through the medium.
        //
        weight 				= exit_data->values[i][0];
		displaced_OPL 		= exit_data->values[i][1];
        refraction_OPL      = exit_data->values[i][2];
        combined_OPL		= exit_data->values[i][3];
        exit_location_x 	= exit_data->values[i][4];        
        exit_location_y		= exit_data->values[i][5];
        exit_location_z		= exit_data->values[i][6];


		/// Assign the OPL used in the calculation.
        ///
        if (DISPLACEMENT_ENABLED)
        {
            OPL = displaced_OPL;		/// Form speckle from the displacement contribution of AO.
        }
        else if (REFRACTION_ENABLED)
        {
            OPL = refraction_OPL;
        }
        else if (COMBINED_ENABLED)
        {
            OPL = combined_OPL;
        }
        else
        {
            cout << "\n\n!!! No AO mechanism assigned for the OPL. Exiting !!!\n";
            exit(-1);
        }

        // For each pixel in the x-axis.
        //
        for (size_t x = 0; x < grid.size(); x++)
        {
            x_pixel = start_x + (dx * (x+1));
            
            // For each pixel in the y-axis.
            //
            for (size_t y = 0; y < grid[0].size(); y++)
            {
                y_pixel = start_y + (dy * (y+1));
                
                // Calculate the distance from the medium's exit aperture to
                // the pixel (x,y) of the CCD.
                distance_to_pixel = sqrt(pow(DISTANCE,2) + 
                                         pow((x_pixel - exit_location_x),2) + 
                                         pow((y_pixel - exit_location_y),2));
                
                
                // The total path length (medium + distance to CCD pixel).
                L = OPL + distance_to_pixel;
                
                // Define an imaginary number.
                std::complex<double> im(0.0, 1.0);
                
                // Compute the interference pattern on this pixel (produces complex value).
                //
                //std::complex<double> temp = distance_to_pixel * weight * exp(-1.0 * im * (L*PI*2/LAMBDA));

                // XXX:
                // - Note that each of the below equations gives a different result and it has been found that the
                //   last is the most accurate, which was verified by looking at the mean intensity provided and comparing
                //   it to changes in the number of photons simulated.  That is, simulating 10K photons and 50K photons
                //   should result in roughly a 5x increase in the mean intensity.
                //grid[x][y] = grid[x][y] + 1.0/((distance_to_pixel) * (1)) * exp(-1.0 * im * (L*PI*2/LAMBDA));
                grid[x][y] = grid[x][y] + (1.0/(distance_to_pixel*distance_to_pixel)) * sqrt(weight) * exp(-1.0 * im * (L*PI*2/LAMBDA));
                
            }
        }
    }
    
    
#ifdef DEBUG
    printGrid();
#endif
    
    writeGridToFile(output_filename);
    
}



void CCDGrid::writeGridToFile(const std::string &output_filename)
{
    //std::string output_filename = "speckle-" + boost::lexical_cast<std::string>(time_step) + ".txt";
	if (speckle_data_stream.is_open())
		speckle_data_stream.close();

    speckle_data_stream.open(output_filename.c_str());
    if (!speckle_data_stream)
    {
        cout << "!!! ERROR: Output file ('" << output_filename << "') could not be opened !!!\n";
        exit(EXIT_FAILURE);
    }
    
    // Set the precision and width of the data written to file.
    //speckle_data_stream.width(10);
    //speckle_data_stream.setf(std::ios::showpoint | std::ios::fixed);
    speckle_data_stream.precision(6);
    
    for (size_t x = 0; x < grid.size(); x++)
    {
        for (size_t y = 0; y < grid[0].size(); y++)
        {
            //double temp = abs(pow(grid[x][y],2));
        	speckle_data_stream << abs(pow(grid[x][y],2)) << "\t";

        }
        speckle_data_stream << "\n";
        speckle_data_stream.flush();
    }
    speckle_data_stream.flush();
}

void CCDGrid::printGrid(void)
{
    // For each pixel in the x-axis.
    //
    for (size_t x = 0; x < grid.size(); x++)
    {

        // For each pixel in the y-axis.
        //
        for (size_t y = 0; y < grid[0].size(); y++)
        {
            cout << "Grid[" << x << "][" << y << "] = " << grid[x][y] << "\n";
        }
    }
}


void CCDGrid::Print_CCD_attributes()
{
    cout << "------------------------------------------" << endl
        << "CCD attributes /" << endl
        << "---------------\n"
        << " - pixels: " << num_x_pixels << "x" << num_y_pixels << endl
        << " - pixel size: " << m_pixel_size << " [meters]" << endl
        << " - dimensions: " << m_pixel_size*num_x_pixels << "x" << m_pixel_size*num_y_pixels << " [meters]" << endl
        << " - coordinates (x,y): " << "(" << m_center_x << "," << m_center_y << ") [meters]" << endl
        << " - detecting AO mechanism: " << m_OPL_mechanism << endl
        << "------------------------------------------" << endl;
}
