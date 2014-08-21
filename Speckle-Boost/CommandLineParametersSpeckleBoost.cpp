/**
 * @file        CommandLineParameters.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * @brief       The implementation file containing the command line parameters
 *
 * @version     kspaceFirstOrder3D 2.13
 * @date        29 August 2012, 11:25 (created) \n
 *              11 October 2012, 17:05 (revised)
 *
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

#include <stdio.h>
#include <getopt.h>
#include <string.h>

#include "CommandLineParametersSpeckleBoost.h"
#include "../Utils/ErrorMessages.h"


/**
 * Constructor
 */
TCommandLineParametersSpeckleBoost::TCommandLineParametersSpeckleBoost() :
        PrintVersion(false),
        InputDirectoryName(""), OutputDirectoryName (""),
        NumberOfThreads(1),
        center_ccd_x(0.0f), center_ccd_y(0.0f),
        x_pixel_cnt(0), y_pixel_cnt(0), pixel_size(0.0f),
        speckle_start_time_step(0), speckle_end_time_step(0),
        compute_displacement_OPL(false), compute_refractive_OPL(false), compute_combined_OPL(false)
{

}// end of TCommandLineParametersSpeckleBoost
//------------------------------------------------------------------------------

/**
 * Print usage and exit.
 */
void TCommandLineParametersSpeckleBoost::PrintUsageAndExit()
{

 printf("---------------------------------- Usage ---------------------------------\n");
 printf("  -i <dir_path>                    : Directory path to exit data that to process\n");
 printf("  -o <dir_path>                    : Directory path where results will be stored\n");
 printf("                                      (creates directory if it doesn't already exist)\n");
 printf("  -t <num_threads>                 : Number of CPU threads to run\n");
 printf("                                      (default = 1\n");
 printf("  --complex_data               : Write data to file in complex valued format\n ");
 printf("\n");
 printf("  --x_pixels <pixels>              : Number of pixels along the x-dimension\n");
 printf("  --y_pixels <pixels>              : Number of pixels along the y-dimension\n");
 printf("  --pixel_size <size>              : Pixel size (meters)\n");
 printf("\n");
 printf("  --x_center <x-coord>             : Center location of the CCD in the x-axis (in meters).\n");
 printf("  --y_center <y-coord>             : Center location of the CCD in the y-axis (in meters).\n");
 printf("\n");
 printf("  --d_OPL                          : Compute speckle pattern based on displacements\n");
 printf("  --n_OPL                          : Compute speckle pattern based on refractive index changes\n");
 printf("  --combined_OPL                   : Compute speckle pattern based on displacements and refractive index changes\n");
 printf("\n");
 printf("  -s <timestep>                    : Time step when speckle pattern processing begins\n");
 printf("                                      (default = 1)\n");
 printf("  -e <timestep>                    : Time step when speckle pattern processing ends\n");
 printf("--------------------------------------------------------------------------\n");
 printf("\n");


 exit(EXIT_FAILURE);

}// end of PrintUsageAndExit
//------------------------------------------------------------------------------

/**
 * Print setup.
 */
void TCommandLineParametersSpeckleBoost::PrintSetup(){

    printf("List of enabled parameters:\n");
    printf("  Input directory           %s\n",InputDirectoryName.c_str());
    printf("  Output directory           %s\n",OutputDirectoryName.c_str());
    printf("\n");
    printf("  Number of threads     %d\n", NumberOfThreads);
    printf("\n");
    printf("  x-pixel count         %d\n", x_pixel_cnt);
    printf("  y-pixel count         %d\n", y_pixel_cnt);
    printf("  pixel size (meters)   %f\n", pixel_size);
    printf("\n");
    printf("  Compute displacement OPL          %d\n", compute_displacement_OPL);
    printf("  Compute refractive OPL            %d\n", compute_refractive_OPL);
    printf("  Compute combination OPL          %d\n", compute_combined_OPL);
    printf("\n");
    printf("  Speckle pattern generation begins at  %d\n", speckle_start_time_step);
    printf("\n");
    printf("  Speckle pattern generation ends at %d\n", speckle_end_time_step);
    /// -----------------------


}// end of PrintSetup
//------------------------------------------------------------------------------

/**
 * Parse command line.
 * @param [in, out] argc
 * @param [in, out] argv
 */
void TCommandLineParametersSpeckleBoost::ParseCommandLine(int argc, char** argv){

   char c;
   int longIndex;
   const char * shortOpts = "i:o:v:c:t:s:e:puIhnd:";

   const struct option longOpts[] = {
        { "version", no_argument, NULL, 0 },

        { "x_pixels", required_argument, NULL, 0 },
        { "y_pixels", required_argument, NULL, 0 },
        { "pixel_size", required_argument, NULL, 0 },
       
        { "x_center", required_argument, NULL, 0 },
        { "y_center", required_argument, NULL, 0 },

        { "d_OPL", no_argument, NULL, 0 },
        { "n_OPL", no_argument, NULL, 0 },
        { "combined_OPL", no_argument, NULL, 0 },

        { "s", required_argument, NULL, 's'},
        { "e", required_argument, NULL, 'e'},

        { "complex_data", no_argument, NULL, 0},
       
        { NULL, no_argument, NULL, 0 }
    };


   // Short parameters //
   while ((c = getopt_long (argc, argv, shortOpts, longOpts, &longIndex )) != -1){
       switch (c){

          case 'i':{
             InputDirectoryName = optarg;
             break;
          }
          case 'o':{
             OutputDirectoryName = optarg;
             break;
          }

          case 't':{
              if ((optarg == NULL) || (atoi(optarg) <= 0)) {
                  fprintf(stderr,"%s", CommandlineParameters_ERR_FMT_NoThreadNumbers);
                  PrintUsageAndExit();
              }else {
                NumberOfThreads = atoi(optarg);
              }

              break;
          }
          case 'h':{

             PrintUsageAndExit();
             break;
          }

          case 's':{
               if ((optarg == NULL) || (atoi(optarg) < 1)) {
                  fprintf(stderr,"%s", CommandlineParameters_ERR_FMT_NoStartTimestep);
                  PrintUsageAndExit();
               }
               speckle_start_time_step = atoi(optarg);

             break;
          }
        
           case 'e':{
               if ((optarg == NULL) || (atoi(optarg) < 1)) {
                   fprintf(stderr,"%s", CommandlineParameters_ERR_FMT_NoEndTimestep);
                   PrintUsageAndExit();
               }
               speckle_end_time_step = atoi(optarg);

               break;

           }
         

           case 0:{   /* long option without a short arg */
                if( strcmp( "version", longOpts[longIndex].name ) == 0 ) {
                    PrintVersion = true;
                    return;
                } else
                if( strcmp( "x_pixels", longOpts[longIndex].name ) == 0) {
                    x_pixel_cnt = atoi(optarg);
                } else
                if( strcmp( "y_pixels", longOpts[longIndex].name ) == 0) {
                    y_pixel_cnt = atoi(optarg);
                } else
                if( strcmp( "pixel_size", longOpts[longIndex].name ) == 0) {
                    x_pixel_cnt = atof(optarg);
                } else
                if( strcmp( "x_center", longOpts[longIndex].name ) == 0) {
                    center_ccd_x = atof(optarg);
                } else
                if( strcmp( "y_center", longOpts[longIndex].name ) == 0) {
                    center_ccd_y = atof(optarg);
                } else
                if( strcmp( "n_OPL", longOpts[longIndex].name ) == 0) {
                    compute_refractive_OPL = true;
                } else
                if( strcmp( "d_OPL", longOpts[longIndex].name ) == 0) {
                    compute_displacement_OPL = true;
                } else
                if( strcmp( "combined_OPL", longOpts[longIndex].name ) == 0) {
                    compute_combined_OPL = true;
                } else
                if( strcmp( "complex_data", longOpts[longIndex].name ) == 0) {
                    write_complex_data = true;
                }
                else {
                    PrintUsageAndExit();
                }
                break;
           }
          default:{
               PrintUsageAndExit();
          }
       }
   }


    //-- Post checks --//



   
   if (InputDirectoryName == "") {
       fprintf(stderr,"%s",CommandlineParameters_ERR_FMT_NoSpeckleInputDirectory);
       PrintUsageAndExit();
   }
    
   if (OutputDirectoryName == "") {

       fprintf(stderr,"%s",CommandlineParameters_ERR_FMT_NoSpeckleOutputDirectory);
       PrintUsageAndExit();
   }


   if (!(compute_displacement_OPL || compute_refractive_OPL || compute_combined_OPL))
       /// -------------------------------
   {
       fprintf(stderr, "%s", "!!! ERROR: Nothing has been specified to simulate. Exiting!\n");
       PrintUsageAndExit();
   }
    
   if ((center_ccd_x == 0.0f) || (center_ccd_y == 0.0f))
   {
       fprintf(stderr, "%s", "!!! ERROR: Center X-Y coordinates for the CCD have not been set. Exiting!\n");
       PrintUsageAndExit();
   }


}// end of ParseCommandLine
//------------------------------------------------------------------------------
