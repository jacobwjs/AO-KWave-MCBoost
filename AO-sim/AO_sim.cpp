//
//  AO_sim.cpp
//  K_Wave_C
//
//  Created by jacob on 12/3/12.
//  Copyright (c) 2012 BMPI. All rights reserved.
//

#include <exception>

#include <AO-sim/AO_sim.h>
#include <MC-Boost/logger.h>
#include <HDF5/HDF5_File.h>

#include <string>
#include <sstream>
#include <fstream>
using std::ofstream;

#include <MC-Boost/vectorMath.h>
using namespace VectorMath;


void PrintMatrix(TRealMatrix &matrix, TParameters *params);
void PrintMatrix(Medium * m_medium, TParameters *params);


AO_Sim::AO_Sim()
{
    // Responsible for running the monte-carlo simulation.
    da_boost = new MC_Boost();

    // Responsible for running the k-Wave simulation.
    KSpaceSolver = new TKSpaceFirstOrder3DSolver();

    /// Set default values for private members.
    m_Laser_injection_coords.x = -1;
    m_Laser_injection_coords.y = -1;
    m_Laser_injection_coords.z = -1;

	/// Set default value.
	MC_time_step = 0.0f;

    /// For use when data is loaded from a previous simulation of ultrasound.
    refractive_total_InputStream = NULL;
    refractive_x_InputStream = refractive_y_InputStream = refractive_z_InputStream = NULL;
    displacement_x_InputStream = displacement_y_InputStream = displacement_z_InputStream = NULL;
}



AO_Sim::~AO_Sim()
{
    cout << "\n\nAO_Sim:: destructor\n";



    if (da_boost)
        delete da_boost;

    if (KSpaceSolver)
        delete KSpaceSolver;

//    if (m_medium)
//        delete m_medium;
}



/// Run the kWave simulation.
void
AO_Sim::Run_kWave_sim(TParameters * Parameters)
{
    /// Load data from disk
    cout << "Loading k-Wave data ........";
    try {
        KSpaceSolver->LoadInputData();
    }catch (ios::failure e) {
        cout << "Failed!\nK-Wave panic: Data loading was not successful!\n%s\n"
        << e.what();
        cerr << "K-Wave panic: Data loading was not successful! \n%s\n"
        << e.what();
        exit(EXIT_FAILURE);
    }
    cout << ".... done\n";

    //fprintf(stdout,"Elapsed time:          %8.2fs\n",KSpaceSolver.GetDataLoadTime());

    /// start computation of k-Wave simulation.
    cout << ".......... k-Wave Computation ...........\n";
    KSpaceSolver->PreCompute();

    
    /// It should still be possible from a kWave_sim to compute and save refractive index changes,
    /// as well as displacements for later use. We see if any of those are enabled and if so we make
    /// those computations in the loop below.
    bool sim_refractive_total = Parameters->IsSim_refractive_total();
    bool sim_refractive_grad  = Parameters->IsSim_refractive_grad();
    bool sim_displacement     = Parameters->IsSim_displacement();
    /// Decide what to simulate (refractive gradient, refractive_total, displacement).
    if (sim_refractive_grad)
    {
        da_boost->Simulate_refractive_gradient(true);
        cout << "Refractive gradient: ON\n";
    }
    if (sim_refractive_total)
    {
        da_boost->Simulate_refractive_total(true);
        cout << "Refraction: ON\n";
    }
    if (sim_displacement)
    {
        da_boost->Simulate_displacement(true);
        cout << "Displacement: ON\n";
    }
    

    //    ActPercent = 0;
    KSpaceSolver->FromAO_sim_PrintOutputHeader();
    KSpaceSolver->IterationTimeStart();
	size_t k_wave_Nt = Parameters->Get_Nt();
    k_wave_Nt = 285;
    for (KSpaceSolver->SetTimeIndex(0); KSpaceSolver->GetTimeIndex() < k_wave_Nt; KSpaceSolver->IncrementTimeIndex()){

        cout << ".......... Running k-Wave ........... ("
        << KSpaceSolver->GetTimeIndex() << " of "
        << Parameters->Get_Nt() << ")\n";

        KSpaceSolver->FromAO_sim_Compute_uxyz();

        // add in the velocity u source term
        KSpaceSolver->FromAO_sim_Add_u_source();

        // add in the transducer source term (t = t1) to ux
        if (Parameters->Get_transducer_source_flag() > KSpaceSolver->GetTimeIndex())
            KSpaceSolver->FromAO_sim_AddTransducerSource();


        KSpaceSolver->FromAO_sim_Compute_duxyz();


        if (Parameters->Get_nonlinear_flag())
            KSpaceSolver->FromAO_sim_Compute_rhoxyz_nonlinear();
        else
            KSpaceSolver->FromAO_sim_Compute_rhoxyz_linear();


        // add in the source pressure term
        KSpaceSolver->FromAO_sim_Add_p_source();

        if (Parameters->Get_nonlinear_flag()){
            KSpaceSolver->FromAO_sim_Compute_new_p_nonlinear();
        }else {
            KSpaceSolver->FromAO_sim_Compute_new_p_linear();

        }


        //-- calculate initial pressure
        if ((KSpaceSolver->GetTimeIndex() == 0) && (Parameters->Get_p0_source_flag() == 1))
            KSpaceSolver->FromAO_sim_Calculate_p0_source();

        
        
        /// If enabled, compute the refractive index for the medium at this time step.
        if (sim_refractive_total)
        {
            KSpaceSolver->FromAO_sim_compute_refractive_index_total();
        }

        
        /// If enabled, compute the displacement for the medium at this time step.
        if (sim_displacement)
        {
            KSpaceSolver->FromAO_sim_compute_displacement();
        }
        
        
        
        /// FIXME:
        /// - There needs to be a calculation done based on the attributes and the medium (SOS) and ultrasound
        ///   frequency, when to save data if --phase_inversion is enabled (i.e. only save data during PI phase shifts
        ///   of the ultrasound propagation). For now this is just a place holder until that is done, and data is always
        ///   saved at every time step over the interval specified.
        ///
        /// If the 'phase_inversion' option is enabled, then we only save data during PI phase shifts.
        if (Parameters->IsPhase_inversion())
        {
            /// Here is where data is stored in their respective sensor data structure if enabled via the commandline.
            /// The stored data is updated over the run of the simulation and written out to disk at the completion
            /// of the simulation.
            KSpaceSolver->FromAO_sim_StoreSensorData();
        }
        else
        {
            /// Here is where data is stored in their respective sensor data structure if enabled via the commandline.
            /// The stored data is updated over the run of the simulation and written out to disk at the completion
            /// of the simulation.
            KSpaceSolver->FromAO_sim_StoreSensorData();
        }


        KSpaceSolver->FromAO_sim_PrintStatistics();

    }

    cout << " .......... Done\n";
    
    KSpaceSolver->IterationTimeStop();
    
    cout << "\n\n-------------------------------------------------------------\n";
    cout << " Timing / \n";
    cout << " -------";
    cout << "\nComputational loop elapsed time: " << KSpaceSolver->GetIterationTime();
    
    cout << "\nPost-processing phase......."; cout.flush();
    KSpaceSolver->PostProcessingTimeStart();
    KSpaceSolver->FromAO_sim_PostProcessing();
    KSpaceSolver->PostProcessingTimeStop();
    cout << "Done \n";
    cout << "Post-processing elapsed time: " << KSpaceSolver->GetPostProcessingTime();
    
    KSpaceSolver->FromAO_sim_WriteOutputDataInfo();
    Parameters->HDF5_OutputFile.Close();
    
    cout << "\nTotal elapsed simulation time: " << KSpaceSolver->GetTotalTime();
    
    
    
    /// Display attributes about the simulation (e.g. max pressure, intensity, displacement, etc.).
    Print_runtime_ultrasound_statistics(Parameters);

}

/// Run the monte-carlo simulation.
void
AO_Sim::Run_monte_carlo_sim(TParameters * Parameters)
{
    /// Explicitly set that we are not running anything other than the monte-carlo simulation.
    da_boost->Simulate_refractive_gradient(false);
    da_boost->Simulate_refractive_total(false);
    da_boost->Simulate_displacement(false);
  
    
    /// Set the monte-carlo simulation to use, or save, RNG seeds based on command line args.
    //Parameters->IsStore_seeds()   ?   da_boost->Save_RNG_seeds(true) : da_boost->Save_RNG_seeds(false);
    //Parameters->IsLoad_seeds()    ?   da_boost->Use_RNG_seeds(true)  : da_boost->Use_RNG_seeds(false);

    
    size_t time = 1;
    da_boost->Run_MC_sim_timestep(m_medium,
                                  m_Laser_injection_coords,
                                  time);
    
}





/// Run the acousto-optic simulation.
void
AO_Sim::Run_acousto_optics_sim(TParameters * Parameters)
{
    /// ------------------ Begin MC-Boost Configuration ------------------------------------------
    /// Configure the monte-carlo simulation to run what was specified via the commandline.
    ///
    bool sim_refractive_total = Parameters->IsSim_refractive_total();
    bool sim_refractive_grad  = Parameters->IsSim_refractive_grad();
    bool sim_displacement     = Parameters->IsSim_displacement();

    /// Decide what to simulate (refractive gradient, refractive_total, displacement).
    if (sim_refractive_grad)
    {
       da_boost->Simulate_refractive_gradient(true);
       cout << "Refractive gradient: ON\n";
    }
    if (sim_refractive_total)
    {
       da_boost->Simulate_refractive_total(true);
       cout << "Refraction: ON\n";
    }
    if (sim_displacement)
    {
       da_boost->Simulate_displacement(true);
       cout << "Displacement: ON\n";
    }
    
    /// FIXME:
    /// - Need to die gracefully with an appropriate error message.
    if (!sim_displacement || !sim_refractive_total || !sim_refractive_grad)
    {
        /// Should fail here, since nothing has been enabled for the AO_sim.
    }

    
    /// If any of the above are turned on, we need the sensor mask indices to map the sensor
    /// to the respective medium (refractive index medium, displacement medium, etc.) so that
    /// it can be updated each run.
    const long * sensor_index = KSpaceSolver->FromAO_sim_Get_sensor_mask_ind().GetRawData();

    /// Load data from disk
    cout << "Loading k-Wave data ........";
    try {
        KSpaceSolver->LoadInputData();
    }catch (ios::failure e) {
        cout << "Failed!\nK-Wave panic: Data loading was not successful!\n%s\n"
             << e.what();
        cerr << "K-Wave panic: Data loading was not successful! \n%s\n"
             << e.what();
        exit(EXIT_FAILURE);
    }
    cout << ".... done\n";



    /// start computation of k-Wave simulation.
    cout << ".......... k-Wave Computation ...........\n";
    KSpaceSolver->PreCompute();


//    ActPercent = 0;
    KSpaceSolver->FromAO_sim_PrintOutputHeader();
    KSpaceSolver->IterationTimeStart();
	size_t k_wave_Nt = Parameters->Get_Nt();
	//k_wave_Nt = 952;
    for (KSpaceSolver->SetTimeIndex(0); KSpaceSolver->GetTimeIndex() < k_wave_Nt; KSpaceSolver->IncrementTimeIndex()){

        cout << ".......... Running k-Wave ........... ("
             << KSpaceSolver->GetTimeIndex() << " of "
             << Parameters->Get_Nt() << ")\n";

        KSpaceSolver->FromAO_sim_Compute_uxyz();

	// add in the velocity u source term
        KSpaceSolver->FromAO_sim_Add_u_source();

	// add in the transducer source term (t = t1) to ux
        if (Parameters->Get_transducer_source_flag() > KSpaceSolver->GetTimeIndex())
            KSpaceSolver->FromAO_sim_AddTransducerSource();


        KSpaceSolver->FromAO_sim_Compute_duxyz();


        if (Parameters->Get_nonlinear_flag())
            KSpaceSolver->FromAO_sim_Compute_rhoxyz_nonlinear();
        else
            KSpaceSolver->FromAO_sim_Compute_rhoxyz_linear();


	// add in the source pressure term
        KSpaceSolver->FromAO_sim_Add_p_source();

        if (Parameters->Get_nonlinear_flag()){
            KSpaceSolver->FromAO_sim_Compute_new_p_nonlinear();
        }else {
            KSpaceSolver->FromAO_sim_Compute_new_p_linear();

        }


	//-- calculate initial pressure
        if ((KSpaceSolver->GetTimeIndex() == 0) && (Parameters->Get_p0_source_flag() == 1))
            KSpaceSolver->FromAO_sim_Calculate_p0_source();



        



        /// --------------------------- Begin Monte-Carlo Simulation ------------------------------------------------------
        ///
		/// If the flag for simulating the influence of refractive index changes or displacements is set
		/// we need to grab the appropriate data from k-Wave and build grids for
		/// photon propagation.
		if (sim_refractive_grad || sim_displacement || sim_refractive_total)
		{
            TRealMatrix *refractive_total_full_medium = NULL;

            TRealMatrix *refractive_x = NULL;
            TRealMatrix *refractive_y = NULL;
            TRealMatrix *refractive_z = NULL;

            TRealMatrix *disp_x_full_medium = NULL;
            TRealMatrix *disp_y_full_medium = NULL;
            TRealMatrix *disp_z_full_medium = NULL;


            if (sim_refractive_grad)
            {
                /// FIXME:
                /// - Not currently functional.
//                KSpaceSolver->FromAO_sim_compute_refractive_index();
//                refractive_x = &(KSpaceSolver->FromAO_sim_Get_refractive_x());
//                refractive_y = &(KSpaceSolver->FromAO_sim_Get_refractive_y());
//                refractive_z = &(KSpaceSolver->FromAO_sim_Get_refractive_z());
            }

            if (sim_refractive_total)
            {
                KSpaceSolver->FromAO_sim_compute_refractive_index_total();
                refractive_total_full_medium = &(KSpaceSolver->FromAO_sim_Get_refractive_total_full_medium());
            }

            if (sim_displacement)
            {
                KSpaceSolver->FromAO_sim_compute_displacement();
                disp_x_full_medium = &(KSpaceSolver->FromAO_sim_Get_disp_x_full_medium());
                disp_y_full_medium = &(KSpaceSolver->FromAO_sim_Get_disp_y_full_medium());
                disp_z_full_medium = &(KSpaceSolver->FromAO_sim_Get_disp_z_full_medium());
            }

			/// NOTE:
			/// - The refractive index and displacement values need to be updated every time step
            ///   within kWave, otherwise what is assigned here is invalid.
            if (sim_refractive_total && sim_displacement)
            {
                /// Create a refractive map (total) based upon the pressure and density changes at this time step.
                m_medium->Create_refractive_map_from_full_medium(refractive_total_full_medium);

                /// Create a displacment map based upon the pressure at this time step.
        		m_medium->Create_displacement_map(disp_x_full_medium,
                                                  disp_y_full_medium,
                                                  disp_z_full_medium);
            }
			else if (sim_refractive_grad && sim_displacement)
			{
				/// Create a refractive map (gradient) based upon the pressure and density changes at this time step.
                /// FIXME:
                /// - Not currently functional.
//        		m_medium->Create_refractive_map(refractive_x,
//                                                refractive_y,
//                                                refractive_z);

				/// Create a displacment map based upon the pressure at this time step.
        		m_medium->Create_displacement_map(disp_x_full_medium,
                                                  disp_y_full_medium,
                                                  disp_z_full_medium);
			}
            else if (sim_refractive_total)
            {
                /// Create a refractive map (total) based upon the pressure and density changes at this time step.
                m_medium->Create_refractive_map_from_full_medium(refractive_total_full_medium);
            }
			else if (sim_refractive_grad)
			{
                /// Create a refractive map (gradient) based upon the pressure and density changes at this time step.
                /// FIXME:
                /// - Not currently functional.
//        		m_medium->Create_refractive_map(refractive_x,
//                                                refractive_y,
//                                                refractive_z);
			}
			else if (sim_displacement)
			{
                /// Create a displacment map based upon the pressure at this time step.
        		m_medium->Create_displacement_map(disp_x_full_medium,
                                                  disp_y_full_medium,
                                                  disp_z_full_medium);

        	}




        	/// Only run the MC-sim after ultrasound has propagated a certain distance (or time).
			/// Similar to the stroboscopic experiments. In essence we only run the the MC_sim
            /// (i.e. light injection) every PI phase step.
            /// If the 'phase_inversion' option was enabled via the commandline we only save
            /// data (i.e. store sensor data) during these phase steps. It limits the total amount
            /// of data saved (pressure values, displacement values, whatever was enabled on the commandline)
            /// to the PI phase shifts of the US wave propagation. This is useful for later simulations where
            /// the medium optical properties change, but the acoustic properties and US propagation is the same.
            /// FIXME:
            /// - Currently 'MC_time_step' requires calculation in 'main.cpp', but this should done
            ///   automatically when reading in the INPUT h5 file. Something along the line of,
            ///   MC_time_step = (US_wavelength/2) * (1/medium_speed_of_sound);
            ///   NOTE:
            ///   - For what is below to work 'MC_time_step' must be evenly divisible by 'dt', which is something
            ///     that needs to be ensured when generating the INPUT h5 file. In fact this entire calculation
            ///     should be made over there and loaded from the INPUT h5 file at runtime. This ensures no mistakes
            ///     arise. The problem can occur when the US wave is advanced beyond a PI phase shift and light is never
            ///     injected at the proper time. This will reduce the detected tagged fraction.
        	static size_t cnt = MC_time_step/Parameters->Get_dt();
            size_t curr_time = KSpaceSolver->GetTimeIndex();
        	if (((curr_time % cnt) == 0) && (curr_time >= 0))
        	{
                
            	da_boost->Run_MC_sim_timestep(m_medium,
                                              m_Laser_injection_coords,
                                              KSpaceSolver->GetTimeIndex());
                
                /// If the 'phase_inversion' option is enabled, then we only save data during PI phase shifts.
                if (Parameters->IsPhase_inversion())
                {
                    /// Here is where data is stored in their respective sensor data structure if enabled via the commandline.
                    /// The stored data is updated over the run of the simulation and written out to disk at the completion
                    /// of the simulation.
                    KSpaceSolver->FromAO_sim_StoreSensorData();
                }
        	}
		} // END if (sim_refractive_grad || sim_displacement )

        //PrintMatrix(KSpaceSolver->FromAO_sim_Get_refractive_total());
        ///
        /// --------------------------- End Monte-Carlo Simulation ------------------------------------------------------

        
        //-- store the initial pressure at the first time step --//
        /// If 'phase_inversion' is NOT enabled, then we save data every time step over the time period specified on the commandline.
        if (! Parameters->IsPhase_inversion())
        {
            /// Here is where data is stored in their respective sensor data structure if enabled via the commandline.
            /// The stored data is updated over the run of the simulation and written out to disk at the completion
            /// of the simulation.
            KSpaceSolver->FromAO_sim_StoreSensorData();
        }
        
        KSpaceSolver->FromAO_sim_PrintStatistics();
        
    }

    cout << " .......... Done\n";

    KSpaceSolver->IterationTimeStop();
    
    cout << "\n\n-------------------------------------------------------------\n";
    cout << " Timing / \n";
    cout << " -------";
    cout << "\nComputational loop elapsed time: " << KSpaceSolver->GetIterationTime();
    
    cout << "\nPost-processing phase......."; cout.flush();
    KSpaceSolver->PostProcessingTimeStart();
    KSpaceSolver->FromAO_sim_PostProcessing();
    KSpaceSolver->PostProcessingTimeStop();
    cout << "Done \n";
    cout << "Post-processing elapsed time: " << KSpaceSolver->GetPostProcessingTime();
    
    KSpaceSolver->FromAO_sim_WriteOutputDataInfo();
    Parameters->HDF5_OutputFile.Close();
    
    cout << "\nTotal elapsed simulation time: " << KSpaceSolver->GetTotalTime();
    
    
    /// Display attributes about the simulation (e.g. max pressure, intensity, displacement, etc.).
    Print_runtime_ultrasound_statistics(Parameters);
        


}



void
AO_Sim::Run_acousto_optics_sim_loadData(TParameters * Parameters)
{
    
    /// FIXME:
    /// - Needs to be an update to this section. Sensor data is read in, but that might not
    ///   be the size of the full medium. This needs to be taken into account by reading
    ///   in sensor data from the h5 file and populating a larger matrix that accounts for
    ///   the full medium.
    ///
    cout << "\n\n\n**************** Implement me *********************\n";
    cout << "AO_Sim::Run_acousto_optics_sim_loadData(TParameters * Parameters)\n\n\n";
    /*
    cout << ".............. Loading precomputed kWave data .............. \n";
    cout << "Opened: " << Parameters->GetOutputFileName() << '\n';


    cout << "Domain dims: [" << Parameters->GetFullDimensionSizes().X << ", "
                             << Parameters->GetFullDimensionSizes().Y << ", "
                             << Parameters->GetFullDimensionSizes().Z << "]"
                             << '\n';
    cout << "Simulation time steps (total): " << Parameters->Get_Nt() << '\n';

    /// How many time steps in the kWave simulation was the refractive index
    /// and/or displacement data recorded.
    size_t recorded_time_steps = -1;

    TRealMatrix *refractive_total_sensor = NULL;
    //TRealMatrix *refractive_x = NULL;
    //TRealMatrix *refractive_y = NULL;
    //TRealMatrix *refractive_z = NULL;

    TRealMatrix *disp_x = NULL;
    TRealMatrix *disp_y = NULL;
    TRealMatrix *disp_z = NULL;


    THDF5_File& HDF5_OutputFile = Parameters->HDF5_OutputFile;
    //THDF5_File& HDF5_InputFile  = Parameters->HDF5_InputFile;


    long int Nx, Ny, Nz;
    TDimensionSizes ScalarSizes(1,1,1);
    HDF5_OutputFile.ReadCompleteDataset(Nx_Name,  ScalarSizes, &Nx);
    HDF5_OutputFile.ReadCompleteDataset(Ny_Name,  ScalarSizes, &Ny);
    HDF5_OutputFile.ReadCompleteDataset(Nz_Name,  ScalarSizes, &Nz);

    /// The dimensions of the recorded data from kWave.
    /// FIXME:
    /// - Should take into account the PML sizes if the sensor.mask
    ///   size did also during the generation of the INPUT h5 file from matlab.
    ///   This might take changing the matlab code for
    ///   saving these values to the INPUT h5 file, which will be written
    ///   out to the OUTPUT h5 file.
    /// NOTE:
    /// - This is not correct if the sensor.mask was not the full medium.
    TDimensionSizes sensor_size(Nx, Ny, Nz);

    /// Is simulation of refractive index changes (total) enabled.  If so create the
    /// input stream to read the data from the HDF5 file and create the refracitve_map
    /// in the medium in which the monte-carlo simulation will propagate through.
    if (Parameters->IsSim_refractive_total())
    {
        /// Notify the monte-carlo simulation that the optical path length should be altered based on spatially varying
        /// refractive index values created from ultrasound via the kWave simulation.
        da_boost->Simulate_refractive_total(true);

        refractive_total_InputStream = new TInputHDF5Stream();
        if (!refractive_total_InputStream)  throw bad_alloc();

        /// Set the HDF5 file to read from, which was the output file from the previous run.
        refractive_total_InputStream->SetHDF5File(HDF5_OutputFile);

        /// FIXME: Currently not used, but the number of time steps to run the simulation are found
        ///      from how many time steps of the ultrasound propagation were recoreded (e.g. refract_total_size.Y).
        TDimensionSizes refract_total_size = HDF5_OutputFile.GetDatasetDimensionSizes(refractive_total_sensor_Name);


        refractive_total_sensor = new TRealMatrix(sensor_size);
        if (!refractive_total_full_medium) throw bad_alloc();

        /// Set the number of iterations to load in data from the HDF5 file.
        /// Due to the way data is written out to the file, the 'Y' dimension of the
        /// data structure contains the number of iterations that were actually written to the HDF5.
        recorded_time_steps = refract_total_size.Y;

    }

    /// Is simulation of displacements enabled.  If so create the
    /// input stream to read the data from the HDF5 file and create the displacement_map
    /// in the medium in which the monte-carlo simulation will propagate through.
    if (Parameters->IsSim_displacement())
    {
        /// Notify the monte-carlo simulation that scattering events should be displaced based on ultrasound from
        /// the kWave simulation.
        da_boost->Simulate_displacement(true);

        displacement_x_InputStream = new TInputHDF5Stream();
        displacement_y_InputStream = new TInputHDF5Stream();
        displacement_z_InputStream = new TInputHDF5Stream();

        if ((!displacement_x_InputStream) || (!displacement_y_InputStream) || (!displacement_z_InputStream))
            throw bad_alloc();

        /// Set the HDF5 file to read from, which was the output file from the previous kWave/AO_sim run,
        /// which is when the data was stored.
        displacement_x_InputStream->SetHDF5File(HDF5_OutputFile);
        displacement_y_InputStream->SetHDF5File(HDF5_OutputFile);
        displacement_z_InputStream->SetHDF5File(HDF5_OutputFile);

        /// FIXME: Currently not used, but the number of time steps to run the simulation are found
        ///      from how many time steps of the ultrasound propagation were recoreded (e.g. disp_x_size.Y).
        TDimensionSizes disp_x_size = HDF5_OutputFile.GetDatasetDimensionSizes(disp_x_Name);
        TDimensionSizes disp_y_size = HDF5_OutputFile.GetDatasetDimensionSizes(disp_y_Name);
        TDimensionSizes disp_z_size = HDF5_OutputFile.GetDatasetDimensionSizes(disp_z_Name);

        disp_x = new TRealMatrix(sensor_size);
        disp_y = new TRealMatrix(sensor_size);
        disp_z = new TRealMatrix(sensor_size);
        if ((!disp_x) || (!disp_y) || (!disp_z)) throw bad_alloc();

        /// Set the number of iterations to load in data from the HDF5 file.
        /// Due to the way data is written out to the file, the 'Y' dimension of the
        /// data structure contains the number of iterations that were actually written to the HDF5.
        /// NOTE:
        /// - The number of iterations of recording data for displacement data will always be the same,
        ///   so we simply use the x-axis displacement number of iterations recorded.
        recorded_time_steps = disp_x_size.Y;
    }

    /// Print the number of time steps that were recorded and stored out of the total
    /// number of time steps simulated.
    cout << "Simulation time steps (recorded): " << recorded_time_steps << '\n';

    /// Is storage of the modulation depth enabled?  If so inform the monte-carlo simulation
    /// to log the optical path lengths via the logger.
    if (Parameters->IsStore_modulation_depth())
    {
        da_boost->Simulate_modulation_depth(true);
    }


    /// Typically when we want to look at the modulation depth (taggged ratio to untagged) we only want
    /// to look at the resulting acousto-optic signal at the ultrasound focus and with the ultrasound
    /// at with a 180 degree phase shift at the same location (i.e. a given propagation time 't'). To achieve
    /// this we simply need to load in the data that was specified (i.e. displacement and/or refractive index vals)
    /// and run the monte-carlo simulation twice, with the recorded data at time 't' and the 180 degree phase
    /// shifted data, also at time 't'. To create the "phase inversion" we simply multiply the recorded data
    /// by -1, which is done with the calls below.
    if (Parameters->IsPhase_inversion())
    {
        /// Typically when data is saved to the HDF5 file, the first recorded time step is the displacement
        /// values at t=1, when no ultrasound is present. Therefore we read past it and use the next
        /// chunk of data, which if the data has been saved correctly, should be at time t=focus.
        /// FIXME:
        /// - Should be able to choose which set of data to invert in cases of when many time steps
        ///   of data were saved. For the time being this will need to be altered here.
        const size_t INVERT_RECORDED_TIME_STEP = 2;
        for (size_t i = 1; i <= INVERT_RECORDED_TIME_STEP; i++)
        {
            if (Parameters->IsSim_refractive_total())
            {
                /// Read refractive total data in from the HDF5 file that holds the precomputed values for
                /// a previously run kWave simulation.
                refractive_total_InputStream->ReadData(refractive_total_sensor_Name, refractive_total_sensor->GetRawData());
            }

            if (Parameters->IsSim_displacement())
            {
                /// Read displacment data in from the HDF5 file that holds the precomputed values for
                /// a previously run kWave simulation.
                displacement_x_InputStream->ReadData(disp_x_Name, disp_x->GetRawData());
                displacement_y_InputStream->ReadData(disp_y_Name, disp_y->GetRawData());
                displacement_z_InputStream->ReadData(disp_z_Name, disp_y->GetRawData());
            }
        }

        /// Now that we have read past the time t=1 and have the appropriate time step read from the
        /// HDF5 file (i.e. time t=focus), we run the monte-carlo simulation twice - once with the
        /// data in the original recorded phase, and again with a 180 degree phase shift.
        ///
        /// First simulation with original phase.
        if (Parameters->IsSim_refractive_total())
        {
            /// Create a refractive map (total) based on the values read in from the HDF5 file.
            //m_medium->Create_refractive_map(refractive_total_sensor, sensor_size);
            cout << "\n\n\n ********* Implement me *********  \n\n\n Parameters->IsSim_refractive_total()\n\n\n";
        }
        if (Parameters->IsSim_displacement())
        {
            /// Create a displacement map based on the values read in from the HDF5 file.
            m_medium->Create_displacement_map(disp_x, disp_y, disp_z);
        }
        cout << "-------------------- Running Phase Inversion (original phase) --------------------\n";
        da_boost->Run_MC_sim_timestep(m_medium,
                                      m_Laser_injection_coords,
                                      INVERT_RECORDED_TIME_STEP);

        /// Now perform the inversion and run the simulation as if the ultrasound focus had reached the
        /// same location and same time, but with a phase that had been shifted 180 degrees.
        if (Parameters->IsSim_refractive_total())
        {
            /// Notify the medium to invert the data for creating the phase shift.
            m_medium->Invert_refractive_map_phase();
        }
        if (Parameters->IsSim_displacement())
        {
            /// Notify the medium to invert the data for creating the phase shift.
            m_medium->Invert_displacement_map_phase();
        }
        cout << "-------------------- Running Phase Inversion (180 deg shifted) ------------------\n";
        da_boost->Run_MC_sim_timestep(m_medium,
                                      m_Laser_injection_coords,
                                      INVERT_RECORDED_TIME_STEP);

    }
    else {
        /// Load data from the HDF5 for each time step of the kWave simulation, and run the monte-carlo simulation
        /// to use the refractive index and displacement data loaded in, which is the acousto-optic simulation at
        /// that point in time.
        for (size_t i = 0; i < recorded_time_steps; i++)
        {
            if (Parameters->IsSim_refractive_total())
            {
                /// Read refractive total data in from the HDF5 file that holds the precomputed values for
                /// a previously run kWave simulation.
                refractive_total_InputStream->ReadData(refractive_total_Name, refractive_total_sensor->GetRawData());


                /// Create a refractive map (total) based on the values read in from the HDF5 file.
                //m_medium->Create_refractive_map(refractive_total_sensor, sensor_size);
                cout << "\n\n\n ********* Implement me *********  \n\n\n Parameters->IsSim_refractive_total()\n\n\n";
                
                ///PrintMatrix((*refractive_total_sensor), Parameters);

                /// Zero out the matrix for the next read in from the HDF5 file.
                ///refractive_total_sensor->ZeroMatrix();
            }

            if (Parameters->IsSim_displacement())
            {
                /// Read displacment data in from the HDF5 file that holds the precomputed values for
                /// a previously run kWave simulation.
                displacement_x_InputStream->ReadData(disp_x_Name, disp_x->GetRawData());
                displacement_y_InputStream->ReadData(disp_y_Name, disp_y->GetRawData());
                displacement_z_InputStream->ReadData(disp_z_Name, disp_y->GetRawData());


                /// Create a displacement map based on the values read in from the HDF5 file.
                m_medium->Create_displacement_map(disp_x, disp_y, disp_z);

                ///PrintMatrix((*disp_x), Parameters);
                ///PrintMatrix(m_medium, Parameters);

                /// Zero out the matrices for the next read in from the HDF5 file.
                //disp_x->ZeroMatrix();
                //disp_y->ZeroMatrix();
                //disp_z->ZeroMatrix();


            }

            /// Set the monte-carlo simulation to use, or save, RNG seeds based on command line args.
            //Parameters->IsStore_seeds()   ?   da_boost->Save_RNG_seeds(true) : da_boost->Save_RNG_seeds(false);
            //Parameters->IsLoad_seeds()    ?   da_boost->Use_RNG_seeds(true)  : da_boost->Use_RNG_seeds(false);

            int time_step = i;
            /// Run the monte-carlo simulation with the loaded in data (displacements, refractive index vals).
            /// XXX:
            /// - Needs testing!!!
            ///cout << "............. Running MC-Boost ........... ";
            ///cout << "(time step: " << time_step << ")\n";
            ///cout.flush();
            da_boost->Run_MC_sim_timestep(m_medium,
                                          m_Laser_injection_coords,
                                          time_step);

        }/// end FOR LOOP

    } /// end IF (PHASE_INVERSION)


    /// Clean up memory.
    if (Parameters->IsSim_refractive_total())
    {
        delete refractive_total_sensor;
        refractive_total_sensor = NULL;
    }
    if (Parameters->IsSim_displacement())
    {
        delete disp_x;
        delete disp_y;
        delete disp_z;
        disp_x = disp_y = disp_z = NULL;
    }

     */
     
}




void
AO_Sim::kWave_allocate_memory()
{
    assert(KSpaceSolver != NULL);
    try {
        KSpaceSolver->AllocateMemory();
    } catch (exception e){
        cout << "Failed!\nK-Wave panic: Not enough memory to run this simulation!\n";
        cout << e.what() << "\n";
        cerr << "K-Wave panic: Not enough memory to run this simulation!\n";
        cerr << e.what() << "\n";

        exit(EXIT_FAILURE);
    }

    cout << ".... done\n";

}


// Creates the grid for the monte-carlo simulation based upon the dimensions
// used in the k-Wave simulation.
//
void
AO_Sim::Create_MC_grid(TParameters * Parameters)
{
    // Number of voxels in each axis, with the PML taken into account.
    // NOTE:
    // - It is assumed that the PML is always included INSIDE the computational
    //   medium (i.e. it is not expanded), as defined in the matlab script, and that there is a small
    //   cushion (5 voxels along US transmission axis) to ensure the probe sits away from the PML.
    //   These things are taken into account below in the monte-carlo simulation grid.
    // - As an example, if the x-axis is 512 elements and the X_PML is 20 elements,
    //   the US probe should sit at X_PML+OFFSET in the k-Wave computation medium, which
    //   means the monte-carlo medium should reflect these changes because US data in the PML
    //   is heavily distorted, as it should be, and we don't want photons moving through those
    //   regions.
    size_t x_pml_offset = Parameters->Get_pml_x_size();
    size_t y_pml_offset = Parameters->Get_pml_y_size();
    size_t z_pml_offset = Parameters->Get_pml_z_size();
    
    
    size_t Nx, Ny, Nz;
    /// Dimension size of the MC_sim grid.
    Nx = Parameters->GetFullDimensionSizes().X;
    Ny = Parameters->GetFullDimensionSizes().Y;
    Nz = Parameters->GetFullDimensionSizes().Z;
    
        

    // Size of each voxel along its respective axis.
    double dx = Parameters->Get_dx();
    double dy = Parameters->Get_dy();
    double dz = Parameters->Get_dz();


    // Create the medium object that represents the grid for the
    // monte-carlo simulation based upon the dimensions for simulating
    // ultrasound.
    //
    m_medium = new Medium(Nx*dx,
                          Ny*dy,
                          Nz*dz);



    // Set the voxel sizes of the 'Medium'.
    this->m_medium->Set_dx(dx);
    this->m_medium->Set_dy(dy);
    this->m_medium->Set_dz(dz);

    // Set the number of voxels of the 'Medium'.
    this->m_medium->Set_Nx(Nx);
    this->m_medium->Set_Ny(Ny);
    this->m_medium->Set_Nz(Nz);

    /// Set the size of the PML's used in the k-Wave simulation.
    this->m_medium->Set_X_PML(x_pml_offset);
    this->m_medium->Set_Y_PML(y_pml_offset);
    this->m_medium->Set_Z_PML(z_pml_offset);


    /// Set the size of the sensor mask used in the kWave simulation.  This is the number
    /// of elements in the TRealMatrix that recorded data.
    // const size_t  sensor_size = Get_sensor_mask_ind().GetTotalElementCount();
    this->m_medium->kwave.sensor_mask_index_size    = Parameters->Get_sensor_mask_index_size();

    /// Set the time-step used in k-Wave.
    this->m_medium->kwave.dt = Parameters->Get_dt();
    //this->m_medium->kwave.speed_of_sound            = Parameters->Get_c0_scalar();


}







void
AO_Sim::Print_MC_sim_params(TParameters * Parameters)
{
    assert (m_medium != NULL);
    assert (da_boost != NULL);
    cout << "-----------------------------------------------------\n"
         << "MC-Boost parameters /\n"
         << "--------------------\n";
    cout << " Number of CPU threads: " << da_boost->Get_CPU_threads() << '\n';
    cout << " Medium size: [x=" << m_medium->Get_X_bound() << ", y=" << m_medium->Get_Y_bound() << ", z=" << m_medium->Get_Z_bound() << "] (meters)\n";
    cout << " Medium dims: [Nx=" << m_medium->Get_Nx() << ", Ny=" << m_medium->Get_Ny() << ", z=" << m_medium->Get_Nz() << "] (voxels)\n";
    cout << " Time step: " << MC_time_step << '\n';
    /// If we are NOT loading seeds, display how many photon energy packets are going to be simulated.
    if (!Parameters->IsLoad_seeds()) cout << " Photons: " << da_boost->Get_num_photons_to_sim() << '\n';
}


void
AO_Sim::Add_circular_detector_MC_medium(Detector_Properties &props)
{


    
    assert(m_medium != NULL);

    /*
    /// Remove error from significant digits.  Just truncate value at the defined decimal place.
    double bottom_x_axis = VectorMath::truncate_to_places(m_medium->Get_X_bound(), 6);
    double top_x_axis = 0.0f;

    double bottom_y_axis = VectorMath::truncate_to_places(m_medium->Get_Y_bound(), 6);
    double top_y_axis = 0.0f;

    double bottom_z_axis = VectorMath::truncate_to_places(m_medium->Get_Z_bound(), 6);
    double top_z_axis = 0.0f;


    double x_coord_plus_radius = VectorMath::truncate_to_places((props.x_coord + props.radius), 6);
    double x_coord_minus_radius = VectorMath::truncate_to_places((props.x_coord - props.radius), 6);
    double x_coord = VectorMath::truncate_to_places(props.x_coord, 6);

    double y_coord_plus_radius = VectorMath::truncate_to_places((props.y_coord + props.radius), 6);
    double y_coord_minus_radius = VectorMath::truncate_to_places((props.y_coord - props.radius), 6);
    double y_coord = VectorMath::truncate_to_places(props.y_coord, 6);

    double z_coord_plus_radius = VectorMath::truncate_to_places((props.z_coord + props.radius), 6);
    double z_coord_minus_radius = VectorMath::truncate_to_places((props.z_coord - props.radius), 6);
    double z_coord = VectorMath::truncate_to_places(props.z_coord, 6);

	

    /// Ensure the detector fits in the medium.
    /// Note: Because the detector is a circular 2-D plane, we need to
    /// check above and below the origin point for the plane (xy, xz, yz) that the detector
    /// sits on (to take into acount its radius).  Otherwise it just needs to
    /// fit between, since it is flat and has no depth along that dimension.
    if (props.xy_plane)
    {
        if ((bottom_x_axis < x_coord_plus_radius) || (x_coord_minus_radius < top_x_axis)) // xy_plane
        {
            cerr << "!!! Error: Circular detector does not fit in medium (x-axis)\n";
            cerr << "AO_Sim::Add_circular_detector_MC_medium()\n";
            exit(EXIT_FAILURE);
        }
        if ((bottom_y_axis < y_coord_plus_radius) || (y_coord_minus_radius < top_y_axis)) // xy_plane
        {
            cerr << "!!! Error: Circular detector does not fit in medium (y-axis)\n";
            cerr << "AO_Sim::Add_circular_detector_MC_medium()\n";
            exit(EXIT_FAILURE);
        }
        if ((bottom_z_axis < z_coord) || (top_z_axis-threshold > z_coord))    // xz_plane
        {
            cerr << "!!! Error: Circular detector does not fit in medium (z-axis)\n";
            cerr << "AO_Sim::Add_circular_detector_MC_medium()\n";
            exit(EXIT_FAILURE);
        }
    }

    if (props.yz_plane)
    {
        if ((bottom_y_axis < y_coord_plus_radius) || (y_coord_minus_radius < top_y_axis))  // xz_plane
        {
            cerr << "!!! Error: Circular detector does not fit in medium (y-axis)\n";
            cerr << "AO_Sim::Add_circular_detector_MC_medium()\n";
            exit(EXIT_FAILURE);
        }
        if ((bottom_z_axis < z_coord_plus_radius) || (z_coord_minus_radius < top_z_axis))  // xz_plane
        {
            cerr << "!!! Error: Circular detector does not fit in medium (z-axis)\n";
            cerr << "AO_Sim::Add_circular_detector_MC_medium()\n";
            exit(EXIT_FAILURE);
        }
        if ((bottom_x_axis < x_coord) || (top_x_axis > x_coord))
        {
            cerr << "!!! Error: Circular detector does not fit in medium (z-axis)\n";
            cerr << "AO_Sim::Add_circular_detector_MC_medium()\n";
            exit(EXIT_FAILURE);
        }

    }

    if (props.xz_plane)
    {
        if ((bottom_z_axis < z_coord_plus_radius) || (z_coord_minus_radius < top_z_axis))
        {
            cerr << "!!! Error: Circular detector does not fit in medium (z-axis)\n";
            cerr << "AO_Sim::Add_circular_detector_MC_medium()\n";
            exit(EXIT_FAILURE);
        }
        if ((bottom_x_axis < x_coord_plus_radius) || (x_coord_minus_radius < top_x_axis)) // xy_plane
        {
            cerr << "!!! Error: Circular detector does not fit in medium (x-axis)\n";
            cerr << "AO_Sim::Add_circular_detector_MC_medium()\n";
            exit(EXIT_FAILURE);
        }
        if ((bottom_y_axis < y_coord) || (top_y_axis > y_coord))    // yz_plane
        {
            cerr << "!!! Error: Circular detector does not fit in medium (y-axis)\n";
            cerr << "AO_Sim::Add_circular_detector_MC_medium()\n";
            exit(EXIT_FAILURE);
        }
    }

    */

    m_medium->addDetector(new CircularDetector(props));

}




/// TEST CASES
///------------------------------------------------------------------------------------------------------------------------------
void
AO_Sim::Test_Read_HDF5_File(TParameters * Parameters)
{
    cout << "------------------- In AO_Sim::Test_Read_HDF5_File() ------------------\n";
    cout << "Opened: " << Parameters->GetOutputFileName() << '\n';


    cout << "Domain dims: [" << Parameters->GetFullDimensionSizes().X << ", "
                             << Parameters->GetFullDimensionSizes().Y << ", "
                             << Parameters->GetFullDimensionSizes().Z << "]"
                             << '\n';

    cout << "Simulation time steps (total): " << Parameters->Get_Nt() << '\n';




//    if (Parameters->IsSim_refractive_index() || Parameters->IsSim_displacement())
//    {
//        precomputed_InputStream = new TInputHDF5Stream();
//        if (!precomputed_InputStream) throw bad_alloc();
//
//        cout << "Loading precomputed data\n";
//        precomputed_InputStream.SetHDF5File(Parameters->HDF5_OutputFile);
//    }

//    if (Parameters->IsSim_refractive_index() || Parameters->IsSim_displacement())
//    {
//        refractive_x_InputStream = new TInputHDF5Stream();
//        if (!refractive_x_InputStream)  throw bad_alloc();
//
//        refractive_y_InputStream = new TInputHDF5Stream();
//        if (!refractive_y_InputStream)  throw bad_alloc();
//
//        refractive_z_InputStream = new TInputHDF5Stream();
//        if (!refractive_z_InputStream)  throw bad_alloc();
//
//        refractive_x_InputStream->SetHDF5File(Parameters->HDF5_OutputFile);
//        refractive_y_InputStream->SetHDF5File(Parameters->HDF5_OutputFile);
//        refractive_z_InputStream->SetHDF5File(Parameters->HDF5_OutputFile);
//
//
//        cout << "Loading refractive index data\n";
//    }
//    if (Parameters->IsSim_displacement())
//    {
//        displacement_x_InputStream = new TInputHDF5Stream();
//        if (!displacement_x_InputStream)  throw bad_alloc();
//
//        displacement_y_InputStream = new TInputHDF5Stream();
//        if (!displacement_y_InputStream)  throw bad_alloc();
//
//        displacement_z_InputStream = new TInputHDF5Stream();
//        if (!displacement_z_InputStream)  throw bad_alloc();
//
//        displacement_x_InputStream->SetHDF5File(Parameters->HDF5_OutputFile);
//        displacement_y_InputStream->SetHDF5File(Parameters->HDF5_OutputFile);
//        displacement_z_InputStream->SetHDF5File(Parameters->HDF5_OutputFile);
//
//
//        cout << "Loading displacement data\n";
//    }

    if (Parameters->IsSim_refractive_total())
    {
        refractive_total_InputStream = new TInputHDF5Stream();
        if (!refractive_total_InputStream)  throw bad_alloc();

        refractive_total_InputStream->SetHDF5File(Parameters->HDF5_OutputFile);
    }


/// --------------------------------------- BEGIN WORKING VERS ---------------------------
//    /// The data is read into a 1D array.
//    const int RANK_OUT = 1;
//
//    herr_t  status;
//    hid_t   dataspace, memspace, dataset, file;
//    hsize_t dims_out[2];
//
//    //float Data[75][384][288];
//
//    THDF5_File& HDF5_OutputFile = Parameters->HDF5_OutputFile;
//    /// Close, so that we can open in a known state.
//    HDF5_OutputFile.Close();
//
//    /// Open the hdf5 file, with flags.
//    file = H5Fopen("FF_debug_OUTPUT_refract_total.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
//    dataset = H5Dopen(file, refractive_total_Name, H5P_DEFAULT);
//
//    /// Select hyperslab in the file.
//    dataspace       = H5Dget_space(dataset);
//    int rank        = H5Sget_simple_extent_ndims(dataspace);
//    int status_n    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
//    cout << "Rank: " << rank << "\nDimensions: " << dims_out[0] << "x" << dims_out[1] << '\n';
//
//
//    for (int i = 0; i < 4; i++)
//    {
//    /// Define hyperslab in the dataset.
//    hsize_t offset[rank];
//    offset[0] = 0;
//    offset[1] = i;
//    offset[2] = 0;
//    hsize_t count[rank];
//    count[0] = 1;
//    count[1] = 1;
//    count[2] = 8294400;
//    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
//    if (status < 0) cout << "Error: H5Sselect_hyperslab() defining hyperslab in dataset\n";
//
//    /// Define the memory dataspace.
//    hsize_t dimsm[1];
//    //dimsm[0] = 1;
//    dimsm[0] = 8294400;
//    memspace = H5Screate_simple(RANK_OUT, dimsm, NULL);
//
//    /// Define the memory hyperslab.
//    hsize_t count_out[1];
//    hsize_t offset_out[1];
//    offset_out[0] = 0;
//    count_out[0] = 8294400;
//    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
//    if (status < 0) cout << "Error: H5Sselect_hyperslab() defining memory in hyperslab\n";
//
//    /// Read data from hyperslab in the file into the hyperslab in memory and display
//    float *data_out = new float[8294400];
//    status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, data_out);
//
//    /// Assign the data to a TRealMatrix
//    TDimensionSizes sensor_size(75, 384, 288);
//    TRealMatrix refract_vals(sensor_size);
//    float *raw = refract_vals.GetRawData();
//    for (int i = 0; i < 8294400; i++)
//    {
//        raw[i] = data_out[i];
//    }
//    delete data_out;
//
//    /// Print data out to file for imagesc() in matlab to verify.
//    PrintMatrix(refract_vals);
//
//    } // end for LOOP
//
/// --------------------------------------- END WORKING VERS ---------------------------




/// --------------------------------------- BEGIN InputHDF5Stream WORKING VERS ---------
    THDF5_File& HDF5_OutputFile = Parameters->HDF5_OutputFile;
    THDF5_File& HDF5_InputFile  = Parameters->HDF5_InputFile;



    TDimensionSizes temp = HDF5_OutputFile.GetDatasetDimensionSizes(refractive_total_full_medium_Name);

    long int Nx, Ny, Nz;
    TDimensionSizes ScalarSizes(1,1,1);
    HDF5_OutputFile.ReadCompleteDataset(Nx_Name,  ScalarSizes, &Nx);
    HDF5_OutputFile.ReadCompleteDataset(Ny_Name,  ScalarSizes, &Ny);
    HDF5_OutputFile.ReadCompleteDataset(Nz_Name,  ScalarSizes, &Nz);

    
    /// FIXME:
    /// - This is incorrect when sensor.mask does not match the full size of the medium.
    TDimensionSizes sensor_size(Nx, Ny, Nz);
    TRealMatrix refract_total_full_medium(sensor_size);

    refractive_total_InputStream->ReadData(refractive_total_full_medium_Name, refract_total_full_medium.GetRawData());
    PrintMatrix(refract_total_full_medium, Parameters);
/// --------------------------------------- END InputHDF5Stream WORKING VERS ---------


}
/// end Test_Read_HDF5_File()




void
AO_Sim::Print_runtime_ultrasound_statistics(TParameters * Parameters)
{
    /// Print statistics about what was specified via the commandline.
    ///
    bool sim_refractive_total = Parameters->IsSim_refractive_total();
    bool sim_refractive_grad  = Parameters->IsSim_refractive_grad();
    bool sim_displacement     = Parameters->IsSim_displacement();
    bool store_p_max          = Parameters->IsStore_p_max();
    bool store_I_max          = Parameters->IsStore_I_max();
    
    
    /// Display statistics about the simulation if options were enabled.
    if (sim_refractive_grad || sim_displacement || sim_refractive_total ||
        Parameters->IsStore_p_max() ||
        Parameters->IsStore_I_max())
    {
        cout << "\n\n-------------------------------------------------------------\n";
        cout << " Statistics / \n";
        cout << " -----------";
        
        if (store_p_max)
        {
            cout << "\n Max pressure: " << KSpaceSolver->stats.max_pressure / 1e6 << " [MPa]";
            cout << "\n - Time index: " << KSpaceSolver->stats.pressure_t_index;
            cout << "\n - Simulation time: " << KSpaceSolver->stats.pressure_t_index * Parameters->Get_dt();
        }
        
        if (store_I_max)
        {
            cout << "\n Max intensity [x-axis]: " << KSpaceSolver->stats.max_intensity_x;
            cout << "\n - Time index: " << KSpaceSolver->stats.intensity_t_index_xaxis;
            cout << "\n - Simulation time: " << KSpaceSolver->stats.intensity_t_index_xaxis * Parameters->Get_dt() << " [sec]";
            
            cout << "\n Max intensity [y-axis]: " << KSpaceSolver->stats.max_intensity_y;
            cout << "\n - Time index: " << KSpaceSolver->stats.intensity_t_index_yaxis;
            cout << "\n - Simulation time: " << KSpaceSolver->stats.intensity_t_index_yaxis * Parameters->Get_dt() << " [sec]";
            
            cout << "\n Max intensity [z-axis]: " << KSpaceSolver->stats.max_intensity_z;
            cout << "\n - Time index: " << KSpaceSolver->stats.intensity_t_index_zaxis;
            cout << "\n - Simulation time: " << KSpaceSolver->stats.intensity_t_index_zaxis * Parameters->Get_dt() << " [sec]";
        }
        if (sim_displacement)
        {
            cout << "\n Max displacement: Implement me!!!";
            cout << "\n - Time index:  Implement me!!!";
            cout << "\n - Simulation time:  Implement me!!!";
        }
        if (sim_refractive_total)
        {
            cout << "\n Max refractive index value: Implement me!!!";
            cout << "\n - Time index:  Implement me!!!";
            cout << "\n - Simulation time:  Implement me!!!";
        }
    } /// end if
    
}


void PrintMatrix(TRealMatrix &data, TParameters *Parameters)
{
    TDimensionSizes Dims;

    Dims.X = Parameters->GetFullDimensionSizes().X;
    Dims.Y = Parameters->GetFullDimensionSizes().Y;
    Dims.Z = Parameters->GetFullDimensionSizes().Z;



    static int time_step;
    std::string filename = "matrix_data";
    std::stringstream ss;
    ss << time_step++;
    filename = filename + "_" + ss.str() + ".txt";

    cout << "Writing Matrix to file: " << filename << '\n';

    std::ofstream data_file_stream;	// Velocity and displacement stream;
    data_file_stream.open(filename.c_str());
	if (!data_file_stream)
	{
		cout << "!!! ERROR: Could not open file for writing.  Check directory structure.\n";
		exit(1);
	}


    int x, y, z;
    z = 151;
    for (x = 0; x < Dims.X; x++)
    {
        for (y = 0; y < Dims.Y; y++)
        {
            data_file_stream << data.GetElementFrom3D(x, y, z) << ' ';
        }
        data_file_stream << '\n';
    }

    data_file_stream.close();
}



void PrintMatrix(Medium * m_medium, TParameters *Parameters)
{
    TDimensionSizes Dims;

    Dims.X = Parameters->GetFullDimensionSizes().X;
    Dims.Y = Parameters->GetFullDimensionSizes().Y;
    Dims.Z = Parameters->GetFullDimensionSizes().Z;



    static int time_step;
    std::string filename = "matrix_data";
    std::stringstream ss;
    ss << time_step++;
    filename = filename + "_" + ss.str() + ".txt";

    cout << "Writing Matrix to file: " << filename << '\n';

    std::ofstream data_file_stream;	// Velocity and displacement stream;
    data_file_stream.open(filename.c_str());
	if (!data_file_stream)
	{
		cout << "!!! ERROR: Could not open file for writing.  Check directory structure.\n";
		exit(1);
	}


    int x, y, z;
    z = 151;
    for (x = 0; x < Dims.X; x++)
    {
        for (y = 0; y < Dims.Y; y++)
        {
            data_file_stream << m_medium->kwave.dmap->getDisplacementFromGridX(x, y, z) << ' ';
        }
        data_file_stream << '\n';
    }

    data_file_stream.close();
}
