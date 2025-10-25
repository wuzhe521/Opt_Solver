#ifndef GONS_CORE_CONFIG_H_
#define GONS_CORE_CONFIG_H_

// ==== Version Information ====
#define GONS_VERSION "1.0.0"

// ==== Floating Point Precision ====
#define GONS_FLT_EPSILON 1.0e-6f

// ==== Feature Toggles ====
#define ENABLE_PRINT_MATRIX 1

#define EMBEDDED_MODE 0
// ==========================

#if (EMBEDDED_MODE == 0)

// ==== Type Definitions ====
#define GONS_INT int64_t

#define GONS_UINT uint64_t

#define GONS_SIZE size_t

#define GONS_FLOAT double
// ==========================
#elif (EMBEDDED_MODE == 1)

// ==== Type Definitions ====
#define GONS_INT int

#define GONS_UINT unsigned int

#define GONS_SIZE long int

#define GONS_FLOAT float
// ==========================
#endif

#define LOG_LVL                                                                \
  1 // 0 Print Nothing
    // 1 Print MSG
    // 2 Pint Location and MSG

const char *SOLVER_HEADER = 
" =================================================== \n\
|  GONS: General Opensource Numerical Solver        |\n\
|  Version: 1.0.0                                   |\n\
|  Author: Wu                                       |\n\
|  Website: https://github.com/nullptr/GONS         |\n\
|  GONS is licensed under the MIT License.          |\n\
|  Copyright (c) 2025 <Wu>. All rights reserved.    |\n\
|  [Warnning]:                                      |\n\
|             Dont expect me to be responsible for  |\n\
|             any damages caused by this software.  |\n\
|             So, use it at your own risk.          |\n\
 =================================================== \n ";



#endif

