#ifndef SA_CONSTANTS_H
#define SA_CONSTANTS_H


/**
 * initial temperature
 */
#define FIXEDINITTEMP "FIXEDINITTEMP"
#define INITFROMSOL   "INITFROMSOL"
#define RANDOMWALKINITTEMP "RANDOMWALKINITTEMP"


/**
 * acceptance criteria
 */
#define METROPOLIS "METROPOLIS"
#define BASICACC   "BASICACC"
#define GEOMACC    "GEOMACC"
#define DETERMINISTICACC "DETERMINISTICACC"


/**
 * exploration method
 */
#define SARANDOMEXPLORATION "SARANDOMEXPLORATION"
#define SASEQUENTIALEXPLORATION "SASEQUENTIALEXPLORATION"
#define SAFIRSTIMPROVEMENTEXPLORATION "SAFIRSTIMPROVEMENTEXPLORATION"
#define SABESTIMPROVEMENTEXPLORATION "SABESTIMPROVEMENTEXPLORATION"

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
 * temperature restart
 */
#define SANOTEMPRESTART "SANOTEMPRESTART"
#define SAMINTEMPRESTART "SAMINTEMPRESTART"
#define SAPERCTEMPRESTART "SAPERCTEMPRESTART"


/**
 * termination criteria
 */
#define MAXBADITERS     "MAXBADITERS" 
#define MAXITERS        "MAXITERS"
#define NEVERTERM       "NEVERTERM"
#define ACCRATETERM     "ACCRATETERM"
#define LASTACCRATETERM "LASTACCRATETERM"

#endif
