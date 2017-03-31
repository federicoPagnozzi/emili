#ifndef SA_CONSTANTS_H
#define SA_CONSTANTS_H

#ifndef ALBERTOSA
#define ALBERTOSA "SA"
#endif


/**
 * initial temperature
 */
#define FIXEDINITTEMP "FIXEDINITTEMP"
#define INITFROMSOL   "INITFROMSOL"
#define RANDOMWALKINITTEMP "RANDOMWALKINITTEMP"
#define RANDOMWALKAVGINITTEMP "RANDOMWALKAVGINITTEMP"
#define RANDOMWALKINITPROB "RANDOMWALKINITPROB"
#define CONNOLLYRWIT "CONNOLLYRWIT"
#define MISEVICIUSINITTEMP "MISEVICIUSINITTEMP"
#define SIMPLEMISEVICIUSINITTEMP "SIMPLEMISEVICIUSINITTEMP"
#define OSMANPOTTSINITTEMP "OSMANPOTTSINITTEMP"


/**
 * acceptance criteria
 */
#define METROPOLIS        "METROPOLIS"
#define METROPOLISWFORCED "METROPOLISWFORCED"
#define BASICACC          "BASICACC"
#define GEOMACC           "GEOMACC"
#define DETERMINISTICACC  "DETERMINISTICACC"
#define LAHCACC           "LAHCACC"
#define GDAACC            "GDAACC"
#define RTRACC            "RTRACC"
#define APPROXEXPACC      "APPROXEXPACC"
#define GENSAACC          "GENSAACC"
#define PRECOMPUTEDMETROPOLIS "PRECOMPUTEDMETROPOLIS"
#define PRECOMPUTEDMETROPOLISWFORCED "PRECOMPUTEDMETROPOLISWFORCED"


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
#define PROBSIZETEMPLEN "PROBSIZETEMPLEN"
#define SQUAREDPROBSIZETEMPLEN "SQUAREDPROBSIZETEMPLEN"
#define BRNEIGHSIZETEMPLEN "BRNEIGHSIZETEMPLEN"
#define BRGEOMTEMPLEN "BRGEOMTEMPLEN"
#define MAXACCEPTEDTEMPLEN "MAXACCEPTEDTEMPLEN"
#define CAPPEDMAXACCEPTEDTEMPLEN "CAPPEDMAXACCEPTEDTEMPLEN"
#define NEIGHCAPPEDMAXACCEPTEDTEMPLEN "NEIGHCAPPEDMAXACCEPTEDTEMPLEN"
#define ARITMTEMPLEN "ARITMTEMPLEN"
#define GEOMTEMPLEN "GEOMTEMPLEN"
#define LOGTEMPLEN "LOGTEMPLEN"
#define EXPTEMPLEN "EXPTEMPLEN"
#define NOTEMPLEN "NOTEMPLEN"

/**
 * cooling schemes
 */
#define GEOM          "GEOM"
#define MARKOV        "MARKOV"
#define LOGCOOLING    "LOGCOOLING"
#define CONSTCOOLING  "CONSTCOOLING"
#define LUNDYMEES     "LUNDYMEES"
#define LUNDYMEESCONNOLLY "LUNDYMEESCONNOLLY"
#define LINEARCOOLING "LINEARCOOLING"
#define TEMPBANDCOOLING "TEMPBANDCOOLING"
#define QUADRATICCOOLING "QUADRATICCOOLING"
#define Q87COOLING    "Q87COOLING"
#define CONNOLLYQ87COOLING    "CONNOLLYQ87COOLING"
#define NOCOOLING     "NOCOOLING"
#define OSMANPOTTSPFSPCOOLING "OSMANPOTTSPFSPCOOLING"

/**
 * temperature restart
 */
#define SANOTEMPRESTART "SANOTEMPRESTART"
#define SAMINTEMPRESTART "SAMINTEMPRESTART"
#define SAPERCTEMPRESTART "SAPERCTEMPRESTART"
#define SALOWRATERESTART "SALOWRATERESTART"
#define SALASTRATERESTART "SALASTRATERESTART"
#define SALOWRATEREHEAT "SALOWRATEREHEAT"
#define SALASTRATEREHEAT "SALASTRATEREHEAT"
#define SALOCALMINREHEAT "SALOCALMINREHEAT"
#define SALOCALMINTEMPRESTART "SALOCALMINTEMPRESTART"
#define SALOCALMINENHANCEDREHEAT "SALOCALMINENHANCEDREHEAT"
#define SALOWRATERESTARTBEST "SALOWRATERESTARTBEST"
#define SALOCALMINRESTARTBEST "SALOCALMINRESTARTBEST"
#define SAMAXITERSTEMPRESTART "SAMAXITERSTEMPRESTART"
#define SAMAXITERSREHEAT "SAMAXITERSREHEAT"
#define SANEIGHSIZEMAXITERSTEMPRESTART "SANEIGHSIZEMAXITERSTEMPRESTART"
#define SASQUAREDNEIGHSIZEMAXITERSTEMPRESTART "SASQUAREDNEIGHSIZEMAXITERSTEMPRESTART"
#define SANEIGHSIZEMAXITERSREHEAT "SANEIGHSIZEMAXITERSREHEAT"
#define SAMAXSTEPSTEMPRESTART "SAMAXSTEPSTEMPRESTART"
#define SAMAXSTEPSREHEAT "SAMAXSTEPSREHEAT"
#define SANEIGHSIZEMAXSTEPSTEMPRESTART "SANEIGHSIZEMAXSTEPSTEMPRESTART"
#define SANEIGHSIZEMAXSTEPSREHEAT "SANEIGHSIZEMAXSTEPSREHEAT"

/**
 * termination criteria
 */
#define MAXBADITERS     "MAXBADITERS" 
#define MAXITERS        "MAXITERS"
#define NEVERTERM       "NEVERTERM"
#define ACCRATETERM     "ACCRATETERM"
#define LASTACCRATETERM "LASTACCRATETERM"
#define NEIGHSIZEITERTERM "NEIGHSIZEITERTERM"
#define SQUAREDNSITERTERM "SQUAREDNSITERTERM"
#define MAXTEMPRESTARTSTERM "MAXTEMPRESTARTSTERM"
#define LOCALMINTERM "LOCALMINTERM"
#define NEIGHSIZELOCALMINTERM "NEIGHSIZELOCALMINTERM"
#define MAXSTEPSTERM "MAXSTEPSTERM"

#endif
