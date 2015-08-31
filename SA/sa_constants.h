#ifndef SA_CONSTANTS_H
#define SA_CONSTANTS_H


/**
 * initial temperature
 */
#define FIXEDINITTEMP "FIXEDINITTEMP"
#define INITFROMSOL   "INITFROMSOL"


/**
 * acceptance criteria
 */
#define METROPOLIS "METROPOLIS"
#define BASICACC   "BASICACC"
#define GEOMACC    "GEOMACC"


/**
 * TempLength
 */
#define CONSTTEMPLEN     "CONSTTEMPLEN"
#define NEIGHSIZETEMPLEN "NEIGHSIZETEMPLEN"

/**
 * cooling schemes
 */
#define GEOM          "GEOM"
#define MARKOV        "MARKOV"
#define LOGCOOLING    "LOGCOOLING"
#define CONSTCOOLING  "CONSTCOOLING"
#define LUNDYMEES     "LUNDYMEES"
#define LINEARCOOLING "LINEARCOOLING"
#define NOCOOLING     "NOCOOLING"

/**
 * termination criteria
 */
#define MAXBADITERS "MAXBADITERS" 
#define MAXITERS    "MAXITERS"

#endif
