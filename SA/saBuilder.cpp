#include "SABuilder.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>

#define ALBERTOSAI "SAI"

bool prs::SABuilder::isCompatibleWith(char* problem_definition)
{
    return true;
}

prs::Component prs::SABuilder::buildComponent(int type)
{
    Component c;
    void* raw_pointer=nullptr;
    switch(type)
    {
    case COMPONENT_COOLING: raw_pointer = buildCooling();break;
    case COMPONENT_TEMP_RESTART: raw_pointer = buildTempRestart();break;
    case COMPONENT_TEMP_LENGTH: raw_pointer =  buildTempLength();break;
    case COMPONENT_EXPLORATION: raw_pointer = buildExploration();break;
    case COMPONENT_INIT_TEMP: raw_pointer = buildInitTemp();break;
    case COMPONENT_ACCEPTANCE: raw_pointer = buildAcceptance();break;
    }

    if(raw_pointer != nullptr)
    {
        c.setType(type);
        c.setRawComponent(raw_pointer);
    }
    else{
        return prs::Builder::buildComponent(type);
    }
    return c;
}

emili::LocalSearch* prs::SABuilder::buildAlgo()
{
   prs::incrementTabLevel();
   emili::LocalSearch* ls = nullptr;
   if(tm.checkToken(ALBERTOSA))
   {
       printTab("Alberto's handmade Simulated Annealing! ");

       emili::InitialSolution* initsol = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();

       emili::Neighborhood* nei = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();

       emili::sa::SAInitTemp* inittemp = retrieveComponent(COMPONENT_INIT_TEMP).get<emili::sa::SAInitTemp>();
       inittemp->setInitialSolution(initsol);
       inittemp->setNeighborhood(nei);
       inittemp->setInstance(gp.getInstance());
       inittemp->setup();

       emili::sa::SAAcceptance* acceptance = retrieveComponent(COMPONENT_ACCEPTANCE).get<emili::sa::SAAcceptance>();
       acceptance->setStartTemperature(inittemp->get());

       emili::sa::SACooling* cooling = retrieveComponent(COMPONENT_COOLING).get< emili::sa::SACooling>();
       cooling->setInitTemp(inittemp->get());
       cooling->setNeighborhood(nei->size());

       emili::sa::SATempRestart* temprestart = retrieveComponent(COMPONENT_TEMP_RESTART).get<emili::sa::SATempRestart>();
       temprestart->setInitTemp(inittemp->get());
       temprestart->setNeighborhoodSize(nei->size());
       cooling->setTempRestart(temprestart);

       emili::sa::SATermination* term = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::sa::SATermination>();
       term->setMinTemp(inittemp->get());
       term->setNeighborhoodsize(nei->size());

       emili::sa::SATempLength* templ = retrieveComponent(COMPONENT_TEMP_LENGTH).get<emili::sa::SATempLength>();
       templ->setNeighborhoodSize(nei->size());
       cooling->setTempLength(templ);

       emili::sa::SAExploration* explo = retrieveComponent(COMPONENT_EXPLORATION).get<emili::sa::SAExploration>();
       explo->setNeighborhood(nei);
       explo->setCooling(cooling);
       explo->setTermination(term);
       explo->setAcceptance(acceptance);
       explo->setTenure();

       ls = new emili::sa::SimulatedAnnealing(initsol,
                                     inittemp,
                                     acceptance,
                                     cooling,
                                     temprestart,
                                     term,
                                     templ,
                                     explo,
                                     nei,
                                     NULL);
   } else if(tm.checkToken(ALBERTOSAI)) {
       printTab("Alberto's handmade incumbent-returning Simulated Annealing! ");

       emili::InitialSolution* initsol = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();

       emili::Neighborhood* nei = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();

       emili::sa::SAInitTemp* inittemp = retrieveComponent(COMPONENT_INIT_TEMP).get<emili::sa::SAInitTemp>();
       inittemp->setInitialSolution(initsol);
       inittemp->setNeighborhood(nei);
       inittemp->setInstance(gp.getInstance());
       inittemp->setup();

       emili::sa::SAAcceptance* acceptance = retrieveComponent(COMPONENT_ACCEPTANCE).get<emili::sa::SAAcceptance>();
       acceptance->setStartTemperature(inittemp->get());

       emili::sa::SACooling* cooling = retrieveComponent(COMPONENT_COOLING).get< emili::sa::SACooling>();
       cooling->setInitTemp(inittemp->get());
       cooling->setNeighborhood(nei->size());

       emili::sa::SATempRestart* temprestart = retrieveComponent(COMPONENT_TEMP_RESTART).get<emili::sa::SATempRestart>();
       temprestart->setInitTemp(inittemp->get());
       temprestart->setNeighborhoodSize(nei->size());
       cooling->setTempRestart(temprestart);

       emili::sa::SATermination* term = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::sa::SATermination>();
       term->setMinTemp(inittemp->get());
       term->setNeighborhoodsize(nei->size());

       emili::sa::SATempLength* templ = retrieveComponent(COMPONENT_TEMP_LENGTH).get<emili::sa::SATempLength>();
       templ->setNeighborhoodSize(nei->size());
       cooling->setTempLength(templ);

       emili::sa::SAExploration* explo = retrieveComponent(COMPONENT_EXPLORATION).get<emili::sa::SAExploration>();
       explo->setNeighborhood(nei);
       explo->setCooling(cooling);
       explo->setTermination(term);
       explo->setAcceptance(acceptance);
       explo->setTenure();

       ls = new emili::sa::SimulatedAnnealingIncumbent(initsol,
                                     inittemp,
                                     acceptance,
                                     cooling,
                                     temprestart,
                                     term,
                                     templ,
                                     explo,
                                     nei,
                                     NULL);
   }
   prs::decrementTabLevel();
   return ls;
}

emili::Termination* prs::SABuilder::buildTermination()
{
    prs::incrementTabLevel();
    emili::Termination* term = nullptr;
    if (tm.checkToken(MAXBADITERS)) {
        prs::printTab("SAMaxBadIterTermination");
        int    mb = tm.getInteger();
        printTabPlusOne("mb",mb);
        term = new emili::sa::SAMaxBadIterTermination(mb);
    } else if (tm.checkToken(MAXITERS)) {
        prs::printTab("SAMaxIterTermination");
        int    mi = tm.getInteger();
        printTabPlusOne("mi",mi);
        term = new emili::sa::SAMaxIterTermination(mi);
    } else if (tm.checkToken(MAXLOCALITERS)) {
        prs::printTab("SAMaxLocalIterTermination");
        int    mi = tm.getInteger();
        printTabPlusOne("mi",mi);
        term = new emili::sa::SAMaxLocalIterTermination(mi);
    } else if (tm.checkToken(NEVERTERM)) {
        prs::printTab("SAWhileTrueTermination");
        term = new emili::sa::SAWhileTrueTermination();
    } else if (tm.checkToken(ACCRATETERM)) {
        prs::printTab("SAAcceptanceRateTermination");
        float rate = tm.getDecimal();
        printTabPlusOne("rate",rate);
        term = new emili::sa::SAAcceptanceRateTermination(rate);
    } else if (tm.checkToken(LASTACCRATETERM)) {
        prs::printTab("SALastAcceptanceRateTermination");
        int te = tm.getInteger();
        float rate = tm.getDecimal();
        printTabPlusOne("te",te);
        printTabPlusOne("rate",rate);
        term = new emili::sa::SALastAcceptanceRateTermination(te, rate);
    } else if (tm.checkToken(MAXTEMPRESTARTSTERM)) {
        prs::printTab("SANeighSizeIterTermination");
        int tr = tm.getInteger();
        printTabPlusOne("tr",tr);
        term = new emili::sa::SaMaxTempRestartsTermination(tr);
    } else if (tm.checkToken(NEIGHSIZEITERTERM)) {
        prs::printTab("SANeighSizeIterTermination");
        float co = tm.getDecimal();
        printTabPlusOne("co",co);
        term = new emili::sa::SANeighSizeIterTermination(co);
    } else if (tm.checkToken(SQUAREDNSITERTERM)) {
        prs::printTab("SASquaredNeighSizeIterTermination");
        float co = tm.getDecimal();
        printTabPlusOne("co",co);
        term = new emili::sa::SASquaredNeighSizeIterTermination(co);
    } else if (tm.checkToken(LOCALMINTERM)) {
        prs::printTab("SALocalMinTermination");
        int te = tm.getInteger();
        printTabPlusOne("te",te);
        term = new emili::sa::SALocalMinTermination(te);
    } else if (tm.checkToken(NEIGHSIZELOCALMINTERM)) {
        prs::printTab("SANeighSizeLocalMinTermination");
        float co = tm.getDecimal();
        printTabPlusOne("co",co);
        term = new emili::sa::SANeighSizeLocalMinTermination(co);
    } else if (tm.checkToken(MAXSTEPSTERM)) {
        prs::printTab("SAMaxStepsTermination");
        int ms = tm.getInteger();
        printTabPlusOne("ms",ms);
        term = new emili::sa::SAMaxStepsTermination(ms);
    } else if (tm.checkToken(MINTEMPTERM)) {
        prs::printTab("SAMinTempTermination");
   /*     char* p = tm.peek();
          double v;//NO NO NO
          try {
           v = std::stod(p);
          } catch(...) {
            return new SAMinTempTermination(inittemp->getMinTemp());
        }*/
        double v = tm.getDecimal();
        printTabPlusOne("v",v);
        term = new emili::sa::SAMinTempTermination(v);
    } else if (tm.checkToken(MINTEMPTERM2)) {
            prs::printTab("SAMinTempTermination without v");
            term = new emili::sa::SAMinTempTermination();
    }
    prs::decrementTabLevel();
    return term;
}

emili::Acceptance* prs::SABuilder::buildAcceptance()
{
    prs::incrementTabLevel();
    emili::Acceptance* acc = nullptr;
    if (tm.checkToken(METROPOLIS)) {
        printTab("SAMetropolisAcceptance");
        acc = new emili::sa::SAMetropolisAcceptance();
    } else if (tm.checkToken(METROPOLISWFORCED)) {
        printTab("SAPrecomputedMetropolisWithForcedAcceptance");
        acc = new emili::sa::SAMetropolisWithForcedAcceptance();
    } else if (tm.checkToken(BASICACC)) {
        printTab("SABasicAcceptance");
        acc = new emili::sa::SABasicAcceptance();
    } else if (tm.checkToken(APPROXEXPACC)) {
        printTab("SAApproxExpAcceptance");
        acc = new emili::sa::SAApproxExpAcceptance();
    } else if (tm.checkToken(GEOMACC)) {
        printTab("SAGeometricAcceptance");
        double rf = tm.getDecimal();
        printTabPlusOne("rf",rf);
        acc = new emili::sa::SAGeometricAcceptance( rf);
    } else if (tm.checkToken(GENSAACC)) {
        printTab("GeneralizedSAAcceptance");
        double beta = tm.getDecimal();
        double g = tm.getDecimal();
        printTabPlusOne("beta",beta);
        printTabPlusOne("g",g);
        acc = new emili::sa::GeneralizedSAAcceptance( beta, g);
    } else if (tm.checkToken(DETERMINISTICACC)) {
        printTab("SADeterministicAcceptance");
        acc = new emili::sa::SADeterministicAcceptance();
    } else if (tm.checkToken(GDAACC)) {
        printTab("GreatDelugeAcceptance");
        acc = new emili::sa::GreatDelugeAcceptance();
    } else if (tm.checkToken(RTRACC)) {
        printTab("RecordToRecordAcceptance");
        double de = tm.getDecimal();
        printTabPlusOne("de",de);
        acc = new emili::sa::RecordToRecordAcceptance(de);
    } else if (tm.checkToken(LAHCACC)) {
        printTab("LAHCAcceptance");
        int te = tm.getInteger();
        printTabPlusOne("te",te);
        acc = new emili::sa::LAHCAcceptance(te);
    } else if (tm.checkToken(PRECOMPUTEDMETROPOLIS)) {
        printTab("SAPrecomputedMetropolisAcceptance");
        int te = tm.getInteger();
        printTabPlusOne("te",te);
        acc = new emili::sa::SAPrecomputedMetropolisAcceptance( te);
    } else if (tm.checkToken(PRECOMPUTEDMETROPOLISWFORCED)) {
        printTab("SAPrecomputedMetropolisWithForcedAcceptance");
        int te = tm.getInteger();
        printTabPlusOne("te",te);
        acc = new emili::sa::SAPrecomputedMetropolisWithForcedAcceptance( te);
    } else if (tm.checkToken(BOUNDEDMETROPOLIS)) {
        printTab("SABoundedMetropolisAcceptance");
        double rd = tm.getDecimal();
        printTabPlusOne("rd",rd);
        acc = new emili::sa::SABoundedMetropolisAcceptance(rd);
    } else if (tm.checkToken(ALLACC)) {
        printTab("SAAcceptanceAll");
        acc = new emili::sa::SAAcceptanceAll();
    }
    prs::decrementTabLevel();
    return acc;
}

emili::sa::SACooling* prs::SABuilder::buildCooling()
{
    prs::incrementTabLevel();
    emili::sa::SACooling* co = nullptr;
    emili::Problem* instance = gp.getInstance();
    if (tm.checkToken(GEOM)) {
        prs::printTab("GeomCooling");
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        printTabPlusOne("a",a);
        printTabPlusOne("b",b);
        co =  new emili::sa::GeomCooling(a,b);
    } else if (tm.checkToken(MARKOV)) {
        printTab("MarkovCooling");
        float b = tm.getDecimal();
        printTabPlusOne("b",b);
        co =  new emili::sa::MarkovCooling(b);
    } else if (tm.checkToken(LOGCOOLING)) {
        printTab("LogCooling");
        float b = tm.getDecimal();
        printTabPlusOne("b",b);
        co =  new emili::sa::LogCooling(b);
    } else if (tm.checkToken(CONSTCOOLING)) {
        printTab("ConstantCooling");
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        printTabPlusOne("a",a);
        printTabPlusOne("b",b);
        co =  new emili::sa::ConstantCooling(a,b);
    } else if (tm.checkToken(LUNDYMEES)) {
        printTab("LundyMeesCooling");
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        printTabPlusOne("a",a);
        printTabPlusOne("b",b);
        co =  new emili::sa::LundyMeesCooling(a,b);
    } else if (tm.checkToken(LUNDYMEESCONNOLLY)) {
        printTab("LundyMeesConnollyCooling");
        co =  new emili::sa::LundyMeesConnollyCooling();
    } else if (tm.checkToken(Q87COOLING)) {
        printTab("Q87Cooling");
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        int   m = tm.getInteger();
        printTabPlusOne("a",a);
        printTabPlusOne("b",b);
        printTabPlusOne("m",m);
        co =  new emili::sa::Q87Cooling(a, b, m);
    } else if (tm.checkToken(CONNOLLYQ87COOLING)) {
        printTab("ConnollyQ87Cooling");
        co =  new emili::sa::ConnollyQ87Cooling();
    } else if (tm.checkToken(LINEARCOOLING)) {
        printTab("LinearCooling");
        float a = tm.getDecimal();
        printTabPlusOne("a",a);
        co =  new emili::sa::LinearCooling(a);
    } else if (tm.checkToken(NOCOOLING)) {
        printTab("NoCooling");
        co =  new emili::sa::NoCooling();
    } else if (tm.checkToken(TEMPBANDCOOLING)) {
        printTab("SATemperatureBandCooling");
        float a = tm.getDecimal();
        printTabPlusOne("a",a);
        co =  new emili::sa::SATemperatureBandCooling(a);
    } else if (tm.checkToken(QUADRATICCOOLING)) {
        printTab("SAQuadraticCooling");
        co =  new emili::sa::SAQuadraticCooling();
    } else if (tm.checkToken(ARITHMETICCOOLING)) {
        printTab("ArithmeticCooling");
        double a = tm.getDecimal();
        printTabPlusOne("a",a);
        co =  new emili::sa::ArithmeticCooling(a);
    } else if (tm.checkToken(OBA1)) {
        printTab("OldBachelor1");
        long   M = tm.getInteger();
        float   delta = tm.getDecimal();
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        float c = tm.getDecimal();
        printTabPlusOne("M",M);
        printTabPlusOne("delta",delta);
        printTabPlusOne("a",a);
        printTabPlusOne("b",b);
        printTabPlusOne("c",c);
        co =  new emili::sa::OldBachelor1(M, delta, a, b, c, instance);
    } else if (tm.checkToken(OBA2)) {
        printTab("OldBachelor2");
        long   M = tm.getInteger();
        float   delta = tm.getDecimal();
        float  d = tm.getDecimal();
        printTabPlusOne("M",M);
        printTabPlusOne("delta",delta);
        printTabPlusOne("d",d);
        co =  new emili::sa::OldBachelor2(M, delta, d);
    } else if (tm.checkToken(OBA3)) {
        printTab("OldBachelor3");
        long   ac = tm.getInteger();
        float  ad = tm.getDecimal();  
        long   dc = tm.getInteger();
        float  dd = tm.getDecimal();  
        long   cc = tm.getInteger();
        printTabPlusOne("accepted_cap",ac);
        printTabPlusOne("accepted_delta",ad);
        printTabPlusOne("discarded_cap",dc);
        printTabPlusOne("discarded_delta",dd);
        printTabPlusOne("updated cap",cc);
        co =  new emili::sa::OldBachelor3(ac, ad, dc, dd, cc);
    } else if (tm.checkToken(OBA4)) {
        printTab("OldBachelor4");
        long   ac = tm.getInteger();
        float  ad = tm.getDecimal();
        long   dc = tm.getInteger();
        float  dd = tm.getDecimal();
        long   cc = tm.getInteger();
        printTabPlusOne("accepted_cap",ac);
        printTabPlusOne("accepted_delta",ad);
        printTabPlusOne("discarded_cap",dc);
        printTabPlusOne("discarded_delta",dd);
        printTabPlusOne("updated cap",cc);
        co =  new emili::sa::OldBachelor4(ac, ad, dc, dd, cc);
    } else if (tm.checkToken(OBA5)) {
        printTab("OldBachelor4");
        long   ac = tm.getInteger();
        float  ad = tm.getDecimal();
        long   dc = tm.getInteger();
        float  dd = tm.getDecimal();
        long   cc = tm.getInteger();
        float  cd = tm.getDecimal();
        printTabPlusOne("accepted_cap",ac);
        printTabPlusOne("accepted_delta",ad);
        printTabPlusOne("discarded_cap",dc);
        printTabPlusOne("discarded_delta",dd);
        printTabPlusOne("updated cap",cc);
        printTabPlusOne("ratio",cd);
        co =  new emili::sa::OldBachelor5(ac, ad, dc, dd, cc, cd);
    } else if (tm.checkToken(OSMANPOTTSPFSPCOOLING)) {
        printTab("OsmanPottsPFSPCooling");
        co = new emili::sa::OsmanPottsPFSPCooling(instance);
    }
    prs::decrementTabLevel();
    return co;
}

emili::sa::SATempRestart* prs::SABuilder::buildTempRestart()
{
    prs::incrementTabLevel();
    emili::sa::SATempRestart* tr = nullptr;
    if (tm.checkToken(SANOTEMPRESTART)) {
        printTab("SANoRestart");
        tr = new emili::sa::SANoRestart();
    } else if (tm.checkToken(SAMINTEMPRESTART)) {
        printTab("SAMinRestart");
        float va = tm.getDecimal();
        printTabPlusOne("va",va);
        tr = new emili::sa::SAMinRestart(va);
    } else if (tm.checkToken(SAPERCTEMPRESTART)) {
        printTab("SAPercRestart");
        float va = tm.getDecimal();
        printTabPlusOne("va",va);
        tr = new emili::sa::SAPercRestart( va);
    } else if (tm.checkToken(SALOWRATERESTART)) {
        printTab("SALowRateRestart");
        float va = tm.getDecimal();
        printTabPlusOne("va",va);
        tr = new emili::sa::SALowRateRestart(va);
    } else if (tm.checkToken(SALASTRATERESTART)) {
        printTab("SALastRateRestart");
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        printTabPlusOne("te",te);
        printTabPlusOne("va",va);
        tr = new emili::sa::SALastRateRestart(te, va);
    } else if (tm.checkToken(SALOWRATEREHEAT)) {
        printTab("SALowRateReheat");
        float th = tm.getDecimal();
        float va = tm.getDecimal();
        printTabPlusOne("th",th);
        printTabPlusOne("va",va);
        tr = new emili::sa::SALowRateReheat(  th, va);
    } else if (tm.checkToken(SALASTRATEREHEAT)) {
        printTab("SALastRateReheat");
        int   te = tm.getInteger();
        float th = tm.getDecimal();
        float va = tm.getDecimal();
        printTabPlusOne("te",te);
        printTabPlusOne("th",th);
        printTabPlusOne("va",va);
        tr = new emili::sa::SALastRateReheat(  te, th, va);
    } else if (tm.checkToken(SALOCALMINREHEAT)) {
        printTab("SALocalMinReheat");
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        printTabPlusOne("te",te);
        printTabPlusOne("va",va);
        tr = new emili::sa::SALocalMinReheat(  te, va);
    } else if (tm.checkToken(SALOWRATERESTARTBEST)) {
        printTab("SALowRateRestartToBest");
        float va = tm.getDecimal();
        printTabPlusOne("va",va);
        tr = new emili::sa::SALowRateRestartToBest(  va);
    } else if (tm.checkToken(SALOCALMINRESTARTBEST)) {
        printTab("SALocalMinRestartToBest");
        int   te = tm.getInteger();
        printTabPlusOne("te",te);
        tr = new emili::sa::SALocalMinRestartToBest(  te);
    } else if (tm.checkToken(SALOCALMINTEMPRESTART)) {
        printTab("SALocalMinTempRestart");
        int   te = tm.getInteger();
        printTabPlusOne("te",te);
        tr = new emili::sa:: SALocalMinTempRestart(  te);
    } else if (tm.checkToken(SALOCALMINENHANCEDREHEAT)) {
        printTab("SALocalMinEnhancedReheat");
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        float ep = tm.getDecimal();
        printTabPlusOne("te",te);
        printTabPlusOne("va",va);
        printTabPlusOne("ep",ep);
        tr = new emili::sa::SALocalMinEnhancedReheat(  te, va, ep);
    } else if (tm.checkToken(SAMAXITERSTEMPRESTART)) {
        printTab("SAMaxItersTempRestart");
        int   te = tm.getInteger();
        printTabPlusOne("te",te);
        tr = new emili::sa::SAMaxItersTempRestart(te);
    } else if (tm.checkToken(SAMAXITERSREHEAT)) {
        printTab("SAMaxItersReheat");
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        printTabPlusOne("te",te);
        printTabPlusOne("va",va);
        tr = new emili::sa::SAMaxItersReheat(  te, va);
    } else if (tm.checkToken(SANEIGHSIZEMAXITERSTEMPRESTART)) {
        printTab("SANeighSizeMaxItersTempRestart");
        float   te = tm.getDecimal();
        printTabPlusOne("te",te);
        tr = new emili::sa::SANeighSizeMaxItersTempRestart(te);
    } else if (tm.checkToken(SASQUAREDNEIGHSIZEMAXITERSTEMPRESTART)) {
        printTab("SASquaredNeighSizeMaxItersTempRestart");
        float   te = tm.getDecimal();
        printTabPlusOne("te",te);
        tr = new emili::sa::SASquaredNeighSizeMaxItersTempRestart(te);
    } else if (tm.checkToken(SANEIGHSIZEMAXITERSREHEAT)) {
        printTab("SANeighSizeMaxItersReheat");
        float   te = tm.getDecimal();
        float va = tm.getDecimal();
        printTabPlusOne("te",te);
        printTabPlusOne("va",va);
        tr = new emili::sa::SANeighSizeMaxItersReheat(te, va);
    } else if (tm.checkToken(SAMAXSTEPSTEMPRESTART)) {
        printTab("SAMaxStepsTempRestart");
        int   te = tm.getInteger();
        printTabPlusOne("te",te);
        tr = new emili::sa::SAMaxStepsTempRestart(  te);
    } else if (tm.checkToken(SAMAXSTEPSREHEAT)) {
        printTab("SAMaxStepsReheat");
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        printTabPlusOne("te",te);
        printTabPlusOne("va",va);
        tr = new emili::sa::SAMaxStepsReheat(  te, va);
    } else if (tm.checkToken(SANEIGHSIZEMAXSTEPSTEMPRESTART)) {
        printTab("SANeighSizeMaxStepsTempRestart");
        float   te = tm.getDecimal();
        printTabPlusOne("te",te);
        tr = new emili::sa::SANeighSizeMaxStepsTempRestart(  te);
    } else if (tm.checkToken(SANEIGHSIZEMAXSTEPSREHEAT)) {
        printTab("SANeighSizeMaxStepsReheat");
        float   te = tm.getDecimal();
        float va = tm.getDecimal();
        printTabPlusOne("te",te);
        printTabPlusOne("va",va);
        tr = new emili::sa::SANeighSizeMaxStepsReheat(te, va);
    }
    prs::decrementTabLevel();
    return tr;
}

emili::sa::SATempLength* prs::SABuilder::buildTempLength()
{
    prs::incrementTabLevel();
    emili::sa::SATempLength* sat = nullptr;
    emili::Problem* instance = gp.getInstance();
    if (tm.checkToken(CONSTTEMPLEN)) {
        printTab("ConstantTempLength");
        int l = tm.getInteger();
        printTabPlusOne("l",l);
        sat = new emili::sa::ConstantTempLength(l);
    } else if (tm.checkToken(RANDOMTEMPLEN)) {
        printTab("RandomTempLength");
        int lb = tm.getInteger();
        printTabPlusOne("lb",lb);
        int ub = tm.getInteger();
        printTabPlusOne("ub",ub);
        sat = new emili::sa::RandomTempLength(lb, ub);
    } else if (tm.checkToken(NEIGHSIZETEMPLEN)) {
        printTab("NeighSizeTempLength");
        float a = tm.getDecimal();
        printTabPlusOne("a",a);
        sat = new emili::sa::NeighSizeTempLength(  a);
    } else if (tm.checkToken(PROBSIZETEMPLEN)) {
        printTab("ProblemSizeTempLength");
        float a = tm.getDecimal();
        printTabPlusOne("a",a);
        sat = new emili::sa::ProblemSizeTempLength(instance, a);
    } else if (tm.checkToken(SQUAREDPROBSIZETEMPLEN)) {
        printTab("SquaredProblemSizeTempLength");
        float a = tm.getDecimal();
        printTabPlusOne("a",a);
        sat = new emili::sa::SquaredProblemSizeTempLength(instance, a);
    } else if (tm.checkToken(BRNEIGHSIZETEMPLEN)) {
        printTab("BurkardRendlNeighSizeTempLength");
        float a = tm.getDecimal();
        printTabPlusOne("a",a);
        sat = new emili::sa::BurkardRendlNeighSizeTempLength(  a);
    } else if (tm.checkToken(MAXACCEPTEDTEMPLEN)) {
        printTab("MaxAcceptedTempLength");
        int a = tm.getInteger();
        printTabPlusOne("a",a);
        sat = new emili::sa::MaxAcceptedTempLength(a);
    } else if (tm.checkToken(CAPPEDMAXACCEPTEDTEMPLEN)) {
        printTab("CappedMaxAcceptedTempLength");
        int a = tm.getInteger();
        printTabPlusOne("a",a);
        int c = tm.getInteger();
        printTabPlusOne("c",c);
        sat = new emili::sa::CappedMaxAcceptedTempLength(a, c);
    } else if (tm.checkToken(NEIGHCAPPEDMAXACCEPTEDTEMPLEN)) {
        printTab("NeighSizeCappedMaxAcceptedTempLength");
        float a = tm.getDecimal();
        printTabPlusOne("a",a);
        int c = tm.getInteger();
        printTabPlusOne("c",c);
        sat = new emili::sa::NeighSizeCappedMaxAcceptedTempLength(  a, c);
    } else if (tm.checkToken(ARITMTEMPLEN)) {
        printTab("ArithmeticTempLength");
        int a = tm.getInteger();
        printTabPlusOne("a",a);
        int b = tm.getInteger();
        printTabPlusOne("b",b);
        sat = new emili::sa::ArithmeticTempLength(a, b);
    } else if (tm.checkToken(GEOMTEMPLEN)) {
        printTab("GeomTempLength");
        int a = tm.getInteger();
        printTabPlusOne("a",a);
        float b = tm.getDecimal();
        printTabPlusOne("b",b);
        sat = new emili::sa::GeomTempLength(a, b);
    } else if (tm.checkToken(LOGTEMPLEN)) {
        printTab("LogTempLength");
        int a = tm.getInteger();
        printTabPlusOne("a",a);
        int b = tm.getInteger();
        printTabPlusOne("b",b);
        sat = new emili::sa::LogTempLength(a, b);
    } else if (tm.checkToken(EXPTEMPLEN)) {
        printTab("ExpTempLength");
        int a = tm.getInteger();
        printTabPlusOne("a",a);
        float b = tm.getDecimal();
        printTabPlusOne("b",b);
        sat = new emili::sa::ExpTempLength(a, b);
    } else if (tm.checkToken(NOTEMPLEN)) {
        printTab("NoTempLength");
        sat = new emili::sa::NoTempLength();
    } else if (tm.checkToken(BRGEOMTEMPLEN)) {
        printTab("BurkardRendlGeomTempLength");
        float b = tm.getDecimal();
        printTabPlusOne("b",b);
        sat = new emili::sa::BurkardRendlGeomTempLength(  b);
    }
    prs::decrementTabLevel();
    return sat;
}

emili::sa::SAExploration* prs::SABuilder::buildExploration()
{
    prs::incrementTabLevel();
    emili::sa::SAExploration* sae = nullptr;
    if (tm.checkToken(SARANDOMEXPLORATION)) {
        printTab("SARandomExploration");
        sae = new emili::sa::SARandomExploration();
    } else if (tm.checkToken(SASEQUENTIALEXPLORATION)) {
        printTab("SASequentialExploration");
        sae = new emili::sa::SASequentialExploration();
    } else if (tm.checkToken(SABESTOFKEXPLORATION)) {
        printTab("SABestOfKExploration");
        long k = tm.getInteger();
        printTabPlusOne("k",k);
        sae = new emili::sa::SABestOfKExploration(k);
    } else if (tm.checkToken(SABESTOFKSEQUENTIALEXPLORATION)) {
        printTab("SABestOfKSequentialExploration");
        long k = tm.getInteger();
        printTabPlusOne("k",k);
        sae = new emili::sa::SABestOfKSequentialExploration(k);
    } else if (tm.checkToken(SANSBESTOFKSEQUENTIALEXPLORATION)) {
        printTab("SANSBestOfKSequentialExploration");
        double k = tm.getDecimal();
        printTabPlusOne("k",k);
        sae = new emili::sa::SANSBestOfKSequentialExploration(k);
    } else if (tm.checkToken(SANSBESTOFKRANDOMEXPLORATION)) {
        printTab("SANSBestOfKRandomExploration");
        double k = tm.getDecimal();
        printTabPlusOne("k",k);
        sae = new emili::sa::SANSBestOfKRandomExploration(k);
    } else if (tm.checkToken(SAFIRSTBESTOFKEXPLORATION)) {
        printTab("SAFirstBestOfKExploration");
        long k = tm.getInteger();
        printTabPlusOne("k",k);
        sae = new emili::sa::SAFirstBestOfKExploration(k);
    } else if (tm.checkToken(SANSFIRSTBESTOFKEXPLORATION)) {
        printTab("SANSFirstBestOfKExploration");
        double k = tm.getDecimal();
        printTabPlusOne("k",k);
        sae = new emili::sa::SANSFirstBestOfKExploration(k);
    }
    prs::decrementTabLevel();
    return sae;
}

emili::sa::SAInitTemp* prs::SABuilder::buildInitTemp()
{
    prs::incrementTabLevel();
    emili::sa::SAInitTemp* it = nullptr;
    if (tm.checkToken(FIXEDINITTEMP)) {
        printTab("Set the initial temperature at a fixed value");
        double value = tm.getDecimal();
        printTabPlusOne("start value : ", value);
        it = new emili::sa::FixedInitTemp(value);
    } else if (tm.checkToken(INITFROMSOL)) {
        printTab("Set the initial temperature starting from an initial solution");
        double value = tm.getDecimal();
        printTabPlusOne("start value : ", value);
        it = new emili::sa::InitTempFromSolution(value);
    } else if (tm.checkToken(RANDOMWALKINITTEMP)) {
        printTab("Do a random walk and take the highest gap as initial temperature");
        int length = tm.getInteger();
        double value = tm.getDecimal();
        printTabPlusOne("start value : ", value);
        printTabPlusOne("length : ", length);
        it = new emili::sa::RandomWalkInitTemp(length, value);
    } else if (tm.checkToken(CONNOLLYRWIT)){
        printTab("ConnollyRandomWalkInitTemp");
        int length = tm.getInteger();
        double value = tm.getDecimal();
        printTabPlusOne("start value : ", value);
        printTabPlusOne("length : ", length);
        it = new emili::sa::ConnollyRandomWalkInitTemp(length,value);
    } else if (tm.checkToken(RANDOMWALKAVGINITTEMP)) {
        printTab("RandomWalkAvgInitTemp");
        int length = tm.getInteger();
        double value = tm.getDecimal();
        printTabPlusOne("start value : ", value);
        printTabPlusOne("length : ", length);
        it = new emili::sa::RandomWalkAvgInitTemp(length,value);
    } else if (tm.checkToken(RANDOMWALKINITPROB)) {
        printTab("RandomWalkInitProb");
        float ip = tm.getDecimal();
        int length = tm.getInteger();
        double value = tm.getDecimal();
        printTabPlusOne("start value : ", value);
        printTabPlusOne("length : ", length);
        it = new emili::sa::RandomWalkInitProb(ip, length, value);
    } else if (tm.checkToken(MISEVICIUSINITTEMP)) {
        printTab("MiseviciusInitTemp");
        int length = tm.getInteger();
        float l11 = tm.getDecimal();
        float l12 = tm.getDecimal();
        float l21 = tm.getDecimal();
        float l22 = tm.getDecimal();
        printTabPlusOne("length : ", length);
        printTabPlusOne("l11 : ", l11);
        printTabPlusOne("l12 : ", l12);
        printTabPlusOne("l21 : ", l21);
        printTabPlusOne("l22 : ", l22);
        it = new emili::sa::MiseviciusInitTemp(length, l11, l12, l21, l22, 1);
    } else if (tm.checkToken(SIMPLEMISEVICIUSINITTEMP)) {
        printTab("SimplifiedMiseviciusInitTemp");
        int length = tm.getInteger();
        float l1 = tm.getDecimal();
        float l2 = tm.getDecimal();
        printTabPlusOne("length : ", length);
        printTabPlusOne("l1 : ", l1);
        printTabPlusOne("l2 : ", l2);
        it = new emili::sa::SimplifiedMiseviciusInitTemp(length, l1, l2, 1);
    } else if (tm.checkToken(OSMANPOTTSINITTEMP)) {
        printTab("OsmanPottsInitTemp");
        float dc = tm.getDecimal();
        float tf = tm.getDecimal();
        float coeff = tm.getDecimal();
        printTabPlusOne("dc : ", dc);
        printTabPlusOne("tf : ", tf);
        printTabPlusOne("coeff : ", coeff);
        it = new emili::sa::OsmanPottsInitTemp(dc, tf, coeff);
    } else if (tm.checkToken(RANDOMWALKSTATSINITTEMP)) {
        printTab("RandomWalkStatsInitTemp");
        int length = tm.getInteger();
        double value = tm.getDecimal();
        printTabPlusOne("length : ", length);
        printTabPlusOne("value : ", value);
        it = new emili::sa::RandomWalkStatsInitTemp(length, value);
    }

    prs::decrementTabLevel();
    return it;
}

/*
 *
 * std::ostringstream oss;
        oss << "number of max iterations "<< ti;
        printTabPlusOne(oss.str().c_str());
 */



emili::LocalSearch* prs::MABuilder::buildAlgo()
{
   prs::incrementTabLevel();
   emili::LocalSearch* ls = nullptr;
   if(tm.checkToken(METROPOLISALGO))
   {
       printTab("Metropolis Algorithm ");

       emili::InitialSolution* initsol = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();

       emili::sa::SAAcceptance* acceptance = retrieveComponent(COMPONENT_ACCEPTANCE).get<emili::sa::SAAcceptance>();
       //acceptance->setStartTemperature(inittemp->get());

       emili::sa::SATermination* term = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::sa::SATermination>();
       //termination->setMinTemp(0);
       //term->setNeighborhoodsize(nei->size());
       // Initial temperature for acceptance, minimum temperature and neighbourhood size for termination are not
       // considered here, to not go crazy...

       SAStatus* sastatus = new SAStatus();

       emili::LocalSearch* ls1 = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
       emili::LocalSearch* ls2 = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();

       //ls1->setSearchStatus(sastatus);
       //ls2->setSearchStatus(sastatus);
       acceptance->set_status(sastatus);

       ls = new emili::metropolis::MetropolisAlgorithm(*initsol,
                                     ls1,
                                     ls2,
                                     acceptance,
                                     term);

       ls->setSearchStatus(sastatus);
   }
   prs::decrementTabLevel();
   return ls;
}

