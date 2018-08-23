/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL                               
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This file implements functions for getting SUNDIALS version
 * information.
 * -----------------------------------------------------------------*/

#include <string.h>

#include <sundials/sundials_version.h>

/* fill string with SUNDIALS version information */
int SUNDIALSGetVersion(char *version, int len)
{
  if (strlen(SUNDIALS_VERSION) > len) {
    return(-1);
  }
  
  strncpy(version, SUNDIALS_VERSION, len);
  return(0);
}

/* fill integers with SUNDIALS major, minor, and patch release 
   numbers and fill a string with the release label */
int SUNDIALSGetVersionNumber(int *major, int *minor, int *patch, 
                             char *label, int len)
{
  if (strlen(SUNDIALS_VERSION_LABEL) > len) {
    return(-1);
  }
  
  *major = SUNDIALS_VERSION_MAJOR;
  *minor = SUNDIALS_VERSION_MINOR;
  *patch = SUNDIALS_VERSION_PATCH;
  strncpy(label, SUNDIALS_VERSION_LABEL, len);

  return(0);
}
