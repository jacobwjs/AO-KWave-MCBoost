#ifndef CCD_GRID_H
#define CCD_GRID_H

#include <vector>
#include <complex>
#include <fstream>
using std::ofstream;

// Forward declaration of class.
class ExitData;

// Distance of exit-aperture of medium to CCD grid. [meters]
#define DISTANCE 0.50

// The center location of the CCD.  That
// is, if the medium exit aperture is 
// centered at (1.0,1.0), then so should
// the CCD.
#define CENTER_X 0.02250
#define CENTER_Y 0.01145

#define PI 3.14159265358979
#define LAMBDA 532e-9



class CCDGrid {
public:
    CCDGrid(void);
	CCDGrid(int x_pixels, 
			int y_pixels,
			double pixel_size,
			int num_detected_photons);
	~CCDGrid();

    // Place for various constructors to call a common initiliazation routine.
    void initCommon(void);
    
    // Sets the dimensions of the grid.
    void setGrid(int x_pixels, int y_pixels, double pixel_size);
    
    // Computes the speckle data which can be used to form the speckle image.
	void makeSpeckle(const std::string &input_filename, const std::string &output_filename);
    
    // Write CCD grid data (interference intensities) to file.
    void writeGridToFile(const std::string &output_filename);

    /// Load data into the CCD
    void Load_exit_data(const std::string &exit_data_filename);
    
    // Debugging.
    void printGrid(void);

	/// Zero the CCD grid, effectively resetting it.
	void Zero_grid(void);
    
private:
	ExitData *exit_data;

	// The 2D grid that holds the intensity values.
	std::vector<std::vector<std::complex<double> > > grid;

	// The pixel dimensions [meters];
	double dx;
	double dy;
    
    // Number of pixels in each axis direction.
    int num_x_pixels;
    int num_y_pixels;
    
    // Distance from exit aperture of medium to CCD.
    double distance_to_CCD;
    
    
    // The output file stream.
    std::ofstream speckle_data_stream;

};

#endif //CCD_GRID_H
