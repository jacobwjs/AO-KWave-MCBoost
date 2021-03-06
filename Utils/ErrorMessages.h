/**
 * @file        ErrorMessages.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia     \n
 *              jiri.jaros@anu.edu.au
 *
 * @brief       The header file containing all error messages of the project
 *
 * @version     kspaceFirstOrder3D 2.13
 * @date        9 August 2011, 12:34          (created) \n
 *              14 September 2012, 14:20      (revised)
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


#ifndef ERRORMESSAGES_H
#define	ERRORMESSAGES_H

//----------------------------- HDF5 error messages --------------------------//

/// HDF5 error message
const char * const HDF5_ERR_FMT_FileNotCreated          = "Error: File \"%s\" could not be created!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_FileCannotRecreated      = "Error: Cannot recreate an opened file \"%s\"!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_FileCannotReopen        = "Error: Cannot reopen an opened file \"%s\"!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_FileNotClosed           = "Error: File \"%s\" could not be closed!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_CouldNotWriteTo         = "Error: Could not write into \"%s\" dataset!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_CouldNotReadFrom        = "Error: Could not read from \"%s\" dataset!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_WrongDimensionSizes     = "Error: Dataset \"%s\"  has wrong dimension sizes!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_FileNotOpened           = "Error: File \"%s\" could not be opened!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_NotHDF5File             = "Error: File \"%s\" is not a valid HDF5 file!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_DatasetNotOpened        = "Error: File \"%s\" could not open dataset \"%s\"!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_CouldNotSetCompression  = "Error: File \"%s\", dataset \"%s\" could set compression level [%d]!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_BadAttributeValue       = "Error: Bad attribute value: [%s,%s] = %s";
/// HDF5 error message
const char * const HDF5_ERR_FMT_CouldNotWriteToAttribute  = "Error: Could not write into \"%s\" attribute of \"%s\" dataset!\n";
/// HDF5 error message
const char * const HDF5_ERR_FMT_CouldNotReadFromAttribute = "Error: Could not read from \"%s\" attribute of \"%s\" dataset!\n";



//---------------------------------- Matrix Classes  -------------------------//

/// Matrix class error message
const char * const  Matrix_ERR_FMT_NotEnoughMemory  = "Error: Class %s: Memory allocation failed: Not Enough Memory\n";
/// Matrix class error message
const char * const  Matrix_ERR_FMT_MatrixNotFloat     = "Error: Matrix [%s] data type is not of single precision floating point!\n";
/// Matrix class error message
const char * const  Matrix_ERR_FMT_MatrixNotReal      = "Error: Matrix [%s] domain is not real!\n";
/// Matrix class error message
const char * const  Matrix_ERR_FMT_MatrixNotComplex   = "Error: Matrix [%s] domain is not complex!\n";
/// Matrix class error message
const char * const  Matrix_ERR_FMT_MatrixNotLong      = "Error: Matrix [%s] data type is not unsigned long!\n";

//--------------------------------- Matrix Container  ------------------------//

/// Matrix container error message
const char * const  MatrixContainer_ERR_FMT_RecordUnknownDistributionType =
    "K-Space panic: Matrix [%s] has unknown distribution type in the C++ code!\n\
    [File, line] : [%s,%d]!\n";

/// Matrix container error message
const char * const  MatrixContainer_ERR_FMT_ReloactaionError =
    "K-Space panic: Matrix [%s] is being reallocated!\n\
    [File, line] : [%s,%d]!\n";


//-------------------------- Command line Parameters  ------------------------//

/// Command line parameters error message
const char * const CommandlineParameters_ERR_FMT_NoVerboseIntreval        = "Command line parsing error: No or invalid verbose interval provided!\n";
/// Command line parameters error message
const char * const CommandlineParameters_ERR_FMT_NoThreadNumbers          = "Command line parsing error: No or invalid number of CPU threads!\n";
/// Command line parameters error message
const char * const CommandlineParameters_ERR_FMT_NoCompressionLevel       = "Command line parsing error: No or invalid compression level!\n";
/// Command line parameters error message
const char * const CommandlineParameters_ERR_FMT_NoStartTimestep          = "Command line parsing error: No or invalid collection start time step!\n";

/// ------------------------------------------------- JWJS ---------------------------------------------------------
/// Command line parameters error message
const char * const CommandlineParameters_ERR_FMT_NoEndTimestep          = "Command line parsing error: No or invalid collection end time step!\n";
/// -----------------------------------------------------------------------------

/// Command line parameters error message
const char * const CommandlineParameters_ERR_FMT_NoBenchmarkTimeStepCount = "Command line parsing error: No or invalid benchmark time step count!\n";
/// Command line parameters error message
const char * const CommandlineParameters_ERR_FMT_NoInputFile =  "Error: The input file was not specified!\n";
/// Command line parameters error message
const char * const CommandlineParameters_ERR_FMT_NoOutputFile=  "Error: The output file was not specified!\n";
/// Command line parameters error message
const char * const Parameters_ERR_FMT_Illegal_alpha_power_value = "Error: Illegal value of alpha_power!";
/// Command line parameters error message
const char * const Parameters_ERR_FMT_Illegal_StartTime_value   = "Error: The start index is out of the simulation span <%ld, %ld>!\n";
/// Command line parameters error message
const char * const Parameters_ERR_FMT_IncorrectInputFileFormat = "Error: Incorrect input file\"%s\" format!\n";
/// Command line parameters error message
const char * const Parameters_ERR_FMT_IncorrectMajorHDF5FileVersion = "Error: Incorrect major version of the HDF5 file %s (expected is %s)!\n";
/// Command line parameters error message
const char * const Parameters_ERR_FMT_IncorrectMinorHDF5FileVersion = "Error: Incorrect minor version of the HDF5 file %s (expected is %s)!\n";

#endif	/* ERRORMESSAGES_H */

