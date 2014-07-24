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
 * The \c -r parameter specifies how often information about the simulation progress is printed to the command line. 
 * By default, the C++ code prints out the progress of the simulation, the elapsed time, and the estimated time of 
 * completion in intervals corresponding to 5% of the total number of times steps.
 * 
 * The \c -c parameter specifies the compression level used by the ZIP library to reduce the 
 * size of the output file. The actual compression rate is highly dependent on the shape of 
 * the sensor mask and the range of stored quantities. In general, the output data is very 
 * hard to compress, and using higher compression levels can greatly increase the time to 
 * save data while not having a large impact on the final file size. The default compression 
 * level of 3 represents a balance between compression ratio and performance that is suitable in most cases.
 * 
 * The \c --benchmark parameter enables the total length of simulation (i.e., the number of time steps) 
 * to be overwritten by setting a new number of time steps to simulate. This is particularly useful 
 * for performance evaluation and benchmarking. As the code performance is relatively stable, 50-100 time steps is 
 * usually enough to predict the simulation duration. This parameter can also be used to quickly find the ideal number of CPU threads to use.
 *
 * The \c -h and \c --help parameters print all the parameters of the C++ code, while the \c --version parameter reports the code version and internal build number.
 * 
 * 
 * The remaining flags specify the output quantities to be recorded during the simulation and stored on disk analogous to 
 * the sensor.record input. If the \c -p or \c --p raw flags are set (these are equivalent), a time series of the acoustic pressure at the 
 * grid points specified by the sensor mask is recorded. If the \c --p rms and/or \c --p max flags are set, the root mean square and/or maximum values 
 * of the pressure at the grid points specified by the sensor mask are recorded. Finally, if the \c --p final flag is set, the values for the entire acoustic pressure 
 * field in the final time step of the simulation is stored (this will always include the PML, regardless of the setting for <tt> `PMLInside' </tt>). Flags to record the acoustic 
 * particular velocity are defined in an analogous fashion.
 * 
 * In addition to returning the acoustic pressure and particle velocity, the acoustic intensity at the grid points specified by 
 * the sensor mask can also be calculated and stored. To account for the staggered grid scheme, approximate values for the particle 
 * velocity at the unstaggered grid points are automatically calculated using linear interpolation before multiplication by the acoustic pressure. 
 * Two means of aggregation are possible: \c -I or \c --I_avg calculates and stores the average acoustic intensity, while \c --I max calculates the maximum acoustic intensity.
 * 
 * Any combination of \c p, \c u and \c I fags is admissible. If no output flag is set, a time-series for
 * the acoustic pressure is recorded. If it is not necessary to collect the output quantities over 
 * the entire simulation, the starting time step when the collection begins can be specified 
 * using the -s parameter. Note, the index for the first time step is 1 (this follows the MATLAB indexing convention).
 *
 *
\verbatim
---------------------------------- Usage ---------------------------------
 /// -------------------------------------- JWJS ----------------------------------------
Simulation flags (which simulation to run): 
  --AO_sim                         : Run the Acousto-Optic simulation with data provided at runtime
  --AO_sim_loadData                : Run the Acousto-Optic simulation with precomputed data
  --MC_sim                         : Run the Monte-Carlo simulation
  --kWave_sim                      : Run the kWave simulation
 /// --------------------------------------------
 
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
  --p_raw                         : Store raw time series of p (default)
  --p_rms                         : Store rms of p
  --p_max                         : Store max of p
  --p_final                       : Store final pressure field 
   
  -u                              : Store ux, uy, uz
                                      (the same as --u_raw)
  --u_raw                         : Store raw time series of ux, uy, uz
  --u_rms                         : Store rms of ux, uy, uz
  --u_max                         : Store max of ux, uy, uz
  --u_final                       : Store final acoustic velocity
   
  -I                              : Store intensity
                                      (the same as --I_avg) 
  --I_avg                         : Store avg of intensity
  --I_max                         : Store max of intensity
   
  -s <timestep>                   : Time step when data collection begins
                                      (default = 1)
 
 /// ----------------------------- JWJS -------------------------------------------------
 --US_freq                        : Ultrasound frequency used in the simulation (used with 'phase_inversion' option)
 --phase_inversion                : Run the acousto-optic simulation at time 't' with ultrasound phase (phi) and (phi+180)
 
 --save_seeds                     : Save the RNG seeds that created photon paths that were detected
 --load_seeds <input_file_name>   : Load RNG seeds to use for photon propagation
 
 --modulation_depth               : Save the optical path lengths to disk for comparison
 
 --n                               : Store index of refraction\n");
                                        (all axial components nx, ny, nz)
 --refractive_total               : Store the norm of the index of refraction
 --refractive_x                   : Store the x-component of the index of refraction
 --refractive_y                   : Store the y-component of the index of refraction
 --refractive_z                   : Store the z-component of the index of refraction
 
 --d                               : Store displacements
                                        (all axial components disp_x, disp_y, disp_z)
 --disp_x                         : Store displacement along x-axis
 --disp_y                         : Store displacement along y-axis
 --disp_z                         : Store displacement along z-axis
 
 --fluence_map                    : Store the fluence in the medium

 -e <timestep>                    : Time step when data collection ends
 
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
 



#ifndef TCOMMANDLINESPARAMETERS_H
#define	TCOMMANDLINESPARAMETERS_H

#include <cstdlib>
#include <string>


/**
 * @class TCommandLineParameters
 * @brief The class to parse and store command line parameters
 */
class TCommandLineParameters {
public:

    /// Constructor
    TCommandLineParameters();
    /// Destructor
    virtual ~TCommandLineParameters() {};
    
    /// Get input file name
    std::string GetInputFileName()      const {return InputFileName;};    
    /// Get output file name
    std::string GetOutputFileName()     const {return OutputFileName;};

    /// ------------------------------- JWJS --------------------------------------------------
    std::string GetPlaneWaveAxis()      const {return Plane_wave_axis;};
    
    /// Get input file name for RNG seeds. Used for loading seeds in monte-carlo simulation.
    std::string GetRNGSeedsFileName()   const {return InputRNGSeedsFileName;};
    /// ------------------------------------/
    
    /// Is --benchmark flag set?
    bool IsBenchmarkFlag()              const {return BenchmarkFlag;};
    /// Is --version flag set
    bool IsVersion()                    const {return PrintVersion; };
    /// Get benchmark time step count
    int  GetBenchmarkTimeStepsCount()   const {return BenchmarkTimeStepsCount;};
    
    /// Get compression level
    int  GetCompressionLevel()          const {return CompressionLevel;};
    /// Get number of threads 
    int  GetNumberOfThreads()           const {return NumberOfThreads;};
    /// Get verbose interval
    int  GetVerboseInterval()           const {return VerboseInterval;};
    /// Get start time index when sensor data collection begins
    int GetStartTimeIndex()             const {return StartTimeStep;};
    
    
    /// ---------------------- JWJS -----------------------------------
    float Get_US_freq()                 const {return US_freq;};
    
    /// Get end time index when sensor data collection ends
    int GetEndTimeIndex()               const {return EndTimeStep;};
    
    /// Set end time index when sensor data collection ends
    void SetEndTimeIndex(const int end) 	  {EndTimeStep = end;};
    /// ----------------------------/
   
    /// Is --p_raw set?
    bool IsStore_p_raw()                const {return Store_p_raw;};
    /// Is --p_rms set?
    bool IsStore_p_rms()                const {return Store_p_rms;};
    /// Is --p_max set?
    bool IsStore_p_max()                const {return Store_p_max;};
    /// Is --p_final set?
    bool IsStore_p_final()              const {return Store_p_final;};
    
    /// Is --u_raw set?
    bool IsStore_u_raw()                const {return Store_u_raw;};
    /// Is --u_rms set?
    bool IsStore_u_rms()                const {return Store_u_rms;};
    /// Is --u_max set?
    bool IsStore_u_max()                const {return Store_u_max;};    
    /// Is --u_final set?
    bool IsStore_u_final()              const {return Store_u_final;};
    
    /// Is --I_avg set
    bool IsStore_I_avg()                const {return Store_I_avg;};
    /// Is --I_max set
    bool IsStore_I_max()                const {return Store_I_max;};
    
    
    /// ---------------------------------------------- JWJS ------------------------------------
    /// Is --Plane_wave
    bool IsPlane_wave()                 const {return Plane_wave;};
    
    /// Is --US_freq set
    bool IsUS_freq_known()              const {return US_freq_known;};
    
    /// Is --save_seeds set
    bool IsStore_RNG_seeds()            const {return Store_seeds;};

    /// Is --load_seeds set
    bool IsLoad_seeds()                 const {return Load_seeds;};

    /// Is --phase_inversion set
    bool IsPhase_inversion()            const {return Phase_inversion;};

    /// Is --modulation_depth set
    bool IsStore_modulation_depth()     const {return Store_modulation_depth;};
    
    /// Is --refractive_total set
    bool IsStore_refractive_total()     const {return Store_refractive_total;};
    /// Is --refractive_x set
    bool IsStore_refractive_x()         const {return Store_refractive_x;};
    /// Is --refractive_y set
    bool IsStore_refractive_y()         const {return Store_refractive_y;};
    /// Is --refractive_z set
    bool IsStore_refractive_z()         const {return Store_refractive_z;};
    
    /// Is --disp_x set
    bool IsStore_disp_x()               const {return Store_disp_x;};
    /// Is --disp_y set
    bool IsStore_disp_y()               const {return Store_disp_y;};
    /// Is --disp_z set
    bool IsStore_disp_z()               const {return Store_disp_z;};
    
    /// Is --fluence_map set
    bool IsStore_fluence_map()          const {return Store_fluence_map;};
    
    /// Is --AO_sim set
    bool IsRun_AO_sim()                 const {return Run_AO_sim;};
    /// Is --AO_sim_loadData set
    bool IsRun_AO_sim_loadData()        const {return Run_AO_sim_loadData;};
    /// Is --MC_sim set
    bool IsRun_MC_sim()                 const {return Run_MC_sim;};
    /// Is --kWave_sim set
    bool IsRun_kWave_sim()              const {return Run_kWave_sim;};
    /// Is --AO_sim_sphere set
    bool IsRun_AO_sim_sphere()          const {return Run_AO_sim_sphere;};
    
    /// ---------------------------------------------------/
    

    /// Print usage and exit
    void PrintUsageAndExit();   
    /// Print setup
    void PrintSetup();          
    /// Parse command line
    void ParseCommandLine(int argc, char** argv);    
    
    
protected:
    /// Copy constructor not allowed for public
    TCommandLineParameters(const TCommandLineParameters& src);
    
    /// operator = not allowed for public
    TCommandLineParameters& operator = (const TCommandLineParameters& src);

private:
    /// Input file name
    std::string InputFileName;
    /// Output file name
    std::string OutputFileName;
    
    /// NumberOfThreads value
    int         NumberOfThreads;
    /// VerboseInterval value
    int         VerboseInterval;
    /// CompressionLevel value
    int         CompressionLevel;
    
    /// BenchmarkFlag value
    bool        BenchmarkFlag;    
    /// BenchmarkTimeStepsCount value
    int         BenchmarkTimeStepsCount;
    /// PrintVersion value
    bool        PrintVersion;    
    
    /// Store_p_raw value
    bool        Store_p_raw;
    /// Store_p_rms value
    bool        Store_p_rms;
    /// Store_p_max value
    bool        Store_p_max;
    /// Store_p_final value
    bool        Store_p_final;
    
    /// Store_u_raw value
    bool        Store_u_raw;
    /// Store_u_rms value
    bool        Store_u_rms;
    /// Store_u_max value
    bool        Store_u_max;
    /// Store_u_final value
    bool        Store_u_final;
    
    
    /// -------------------------------------------- JWJS --------------------------------------------
    ///
    ///
    ///             ---- Options to decide what is run, via the commandline ----         ///
    ///
    /// Used to signify that a 'perfect' plane wave is propagated. That is,
    /// the wave equation is only propagated along 1-dimension.
    bool        Plane_wave;
    /// Input name for the axis to propagate the plane wave along (e.g. x-axis)
    std::string Plane_wave_axis;
    
    /// Used to signal if the US frequency is input via commandline option.
    bool        US_freq_known;
    /// The ultrasound frequency used in the simulation.
    float       US_freq;
    
    /// Store the RNG seeds used in monte-carlo simulation.
    bool        Store_seeds;

    /// Use RNG seeds from a previous run of the monte-carlo simulation.
    bool        Load_seeds;
    /// Input file name for RNG seeds
    std::string InputRNGSeedsFileName;

    /// Perform a phase inversion acousto-optic simulation, which
    /// simply runs the acousto-optic simulation on displacement
    /// and/or refractive index data, at a recorded time step,
    /// and inverts the data such that it was as if the ultrasound
    /// had made it to the same location and time, but with a 180
    /// phase shift.
    bool        Phase_inversion;

    /// Store the modulation depth (tagged vs. untagged photons)
    bool        Store_modulation_depth;
    
    /// Store index of refraction values (total)
    bool        Store_refractive_total;
    /// Store index of refraction values of the x-component
    bool        Store_refractive_x;
    /// Store index of refraction values of the y-component
    bool        Store_refractive_y;
    /// Store index of refraction values of the z-component
    bool        Store_refractive_z;
    
    /// Store displacement along x-axis
    bool        Store_disp_x;
    /// Store displacement along y-axis
    bool        Store_disp_y;
    /// Store displacement along z-axis
    bool        Store_disp_z;
    
    /// Store the fluence map
    bool        Store_fluence_map;
    
    
    /// Run AO_sim that uses data as it's produced from kWave (i.e. US simulation runs)
    bool        Run_AO_sim;
    /// Run AO_sim that loads precomputed displacement and refractive index data
    bool        Run_AO_sim_loadData;
    /// Run only the monte-carlo simulation
    bool        Run_MC_sim;
    /// Run only the kWave simulation
    bool        Run_kWave_sim;
    /// Run AO_sim that approximates the tagging volume by a sphere filled with given values.
    bool        Run_AO_sim_sphere;
    
    
    
    /// -------------------------------------------------/
    
    
    /// Store_I_avg value
    bool        Store_I_avg;
    /// Store_I_max value
    bool        Store_I_max;
    /// StartTimeStep value
    int         StartTimeStep;
    
    
    /// ------------------------- JWJS -----------------------------------------------
    /// EndTimeStep value (when to stop collecting data)
    int         EndTimeStep;
    /// -------------------------------
    

    /// Default compression level 
    static const int DefaultCompressionLevel = 3;
    /// Default verbose interval
    static const int DefaultVerboseInterval  = 5;
    
    
};// end of class TCommandLineParameters

#endif	/* TCOMMANDLINESPARAMETERS_H */

