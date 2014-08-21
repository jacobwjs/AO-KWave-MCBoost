/**
 * @file        CommandLineParameters.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * @brief       The header file containing the command line parameters
 * 
 * @version     kspaceFirstOrder3D 2.13
 * @date        29 August 2012, 11:25 (created) \n        
 *              11 October 2012, 17:05 (revised) 
 *              03 August 2014, 16:49 (revised - Jacob Staley)
 * 
 * @section Params Command Line Parameters
 * 
 * The  C++ code requires two mandatory parameters and accepts a few optional parameters and flags. 
 * The mandatory parameters \c -i and \c -o specify the input and output file. The file names respect the path conventions for particular operating system.  
 * If any of the files is not specified, cannot be found or created, an error message is shown.
 * 
 * The \c -t parameter sets the number of threads used, which defaults the system maximum. On CPUs with Intel Hyper-Threading (HT), 
 * performance will generally be better if HT is disabled in the BIOS settings. If HT is switched on, the default will be to create 
 * twice as many threads as there are physical processor cores, which might slow down the code execution. Therefore, if the HT is on,
 *  try specifying the number of threads manually for best performance (e.g., 4 for Intel Quad-Core). We recommend experimenting 
 * with this parameter to find the best configuration. Note, if there are other tasks being executed on the system, 
 * it might be useful to further limit the number of threads to prevent system overload.
 * 
 *
\verbatim
---------------------------------- Usage ---------------------------------
 Simulation flags (which simulation to run):
  --AO_sim                         : Run the Acousto-Optic simulation with data provided at runtime
 
Mandatory parameters:
  -i <input_file_name>            : HDF5 input file 
  -o <output_file_name>           : HDF5 output file
 
Optional parameters: 
  -t <num_threads>                : Number of CPU threads (default = MAX)
  -r <interval_in_%>              : Progress print interval (default = 5%)
  -c <comp_level>                 : Output file compression level <0,9>
                                      (default = 3)
  --benchmark <steps>             : Run a specified number of time steps
   
  -h                              : Print help
  --help                          : Print help
  --version                       : Print version
   
Output flags:   
  -p                              : Store acoustic pressure 
                                      (default if nothing else is on)
                                      (the same as --p_raw)
 
----------------------------------------  
\endverbatim
 * 
 * 
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
 * Copyright (C) 2012 Jiri Jaros and Bradley Treeby
 * 
 * This file is part of k-Wave. k-Wave is free software: you can redistribute it 
 * and/or modify it under the terms of the GNU Lesser General Public License as 
 * published by the Free Software Foundation, either version 3 of the License, 
 * or (at your option) any later version.
 * 
 * k-Wave is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the GNU Lesser General Public License for more details. 
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with k-Wave. If not, see <http://www.gnu.org/licenses/>.
 */ 
 



#ifndef TCOMMANDLINESPARAMETERS_SPECKLE_H
#define	TCOMMANDLINESPARAMETERS_SPECKLE_H

#include <cstdlib>
#include <string>


/**
 * @class TCommandLineParametersSpeckleBoost
 * @brief The class to parse and store command line parameters
 */
class TCommandLineParametersSpeckleBoost {
public:

    /// Constructor
    TCommandLineParametersSpeckleBoost();
    /// Destructor
    virtual ~TCommandLineParametersSpeckleBoost() {};
    
    /// Get input file name
    std::string GetInputDirectoryName()      const {return InputDirectoryName;};
    /// Get output file name
    std::string GetOutputDirectoryName()     const {return OutputDirectoryName;};
    
    
    /// Is --version flag set
    bool IsVersion()                    const {return PrintVersion; };
        /// Get number of threads
    int  GetNumberOfThreads()           const {return NumberOfThreads;};
    
    /// Get start time index when sensor data collection begins
    int GetSpeckleStartTimeIndex()             const {return speckle_start_time_step;};
    
    /// Get end time index when sensor data collection ends
    int GetSpeckleEndTimeIndex()        const {return speckle_end_time_step;};
    
    /// Return the size of an individual pixel (meters)
    double GetPixelSize()                const {return pixel_size;};
    
    /// Return the number of pixels along each dimension.
    int GetNumberOfPixelsXdim()         const {return x_pixel_cnt;};
    int GetNumberOfPixelsYdim()         const {return y_pixel_cnt;};
    
    /// Return the location of the center of the CCD.
    double GetCCDxcoord()               const {return center_ccd_x;};
    double GetCCDycoord()               const {return center_ccd_y;};

    /// Is --refractive_OPL set
    bool IsCompute_refractive_OPL()     const {return compute_refractive_OPL;};
    /// Is --displacement_OPL set
    bool IsCompute_displacement_OPL()   const {return compute_displacement_OPL;};
    /// is --combination_OPL set
    bool IsCompute_combined_OPL()       const {return compute_combined_OPL;};
    
    /// Is --AO_sim set
    bool IsRun_AO_sim()                 const {return Run_AO_sim;};
    
    /// Is --complex_data set
    bool IsWrite_complex_data()         const {return write_complex_data;};
    

    /// Print usage and exit
    void PrintUsageAndExit();   
    /// Print setup
    void PrintSetup();          
    /// Parse command line
    void ParseCommandLine(int argc, char** argv);    
    
    
protected:
    /// Copy constructor not allowed for public
    TCommandLineParametersSpeckleBoost(const TCommandLineParametersSpeckleBoost& src);
    
    /// operator = not allowed for public
    TCommandLineParametersSpeckleBoost& operator = (const TCommandLineParametersSpeckleBoost& src);

private:
    /// Input file name
    std::string InputDirectoryName;
    /// Output file name
    std::string OutputDirectoryName;
    
    /// NumberOfThreads value
    int         NumberOfThreads;
  
    /// PrintVersion value
    bool        PrintVersion;
    
    /// Boolean to decide if data should be left in complex form when it is written out to disk.
    bool        write_complex_data;
    
    /// Create the interference pattern based on the contribution from displacements
    bool        compute_displacement_OPL;
    /// Create the interference pattern based on the contribution from refractive index changes
    bool        compute_refractive_OPL;
    /// Create the interference pattern based on the contribution from displacements + refractive index changes
    bool        compute_combined_OPL;
    
    /// Run AO_sim that uses data as it's produced from kWave (i.e. US simulation runs)
    bool        Run_AO_sim;

    /// StartTimeStep value (what time step to start processing photon exit data)
    int         speckle_start_time_step;
    
    /// EndTimeStep value (what time step to stop processing photon exit data)
    int         speckle_end_time_step;
    
    /// Center location of the CCD in cartesian coordinates.
    double       center_ccd_x;
    double       center_ccd_y;
    
    /// Pixel count along the x and y axis.
    int         x_pixel_cnt;
    int         y_pixel_cnt;
    
    /// The size of the individual pixels.
    double       pixel_size;

};// end of class TCommandLineParametersSpeckleBoost

#endif	/* TCOMMANDLINESPARAMETERS_SPECKLE_H */

