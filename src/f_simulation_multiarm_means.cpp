#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

#include "f_utilities.h"
#include "rpact_types.h"

using namespace Rcpp;

double matrixMult(NumericVector x, NumericVector y) {
  double result = 0;
  for(int i = 0; i < x.length(); i++) {
    result += x[i] * y[i];
  }
  return result;
}


// [[Rcpp::export(name = ".getSimulatedStageMeansMultiArmCpp")]]
List getSimulatedStageMeansMultiArmCpp(Environment design,
                                       NumericVector muVector,
                                       double stDev,
                                       NumericVector plannedSubjects,
                                       std::string typeOfSelection,
                                       std::string effectMeasure,
                                       LogicalVector adaptations,
                                       double epsilonValue,
                                       double rValue,
                                       double threshold,
                                       NumericVector allocationRatioPlanned,
                                       NumericVector minNumberOfSubjectsPerStage,
                                       NumericVector maxNumberOfSubjectsPerStage,
                                       double conditionalPower,
                                       double thetaH1,
                                       double stDevH1,
                                       Function calcSubjectsFunction,
                                       bool calcSubjectsFunctionIsUserDefined,
                                       Nullable<Function> selectArmsFunction) {
  
  int kMax = plannedSubjects.length();
  int gMax = muVector.length();
  NumericMatrix simMeans(gMax + 1, kMax, ((NumericVector) rep(NA_REAL,(gMax + 1) * kMax)).begin());
  NumericMatrix overallEffects(gMax, kMax, ((NumericVector) rep(NA_REAL,gMax * kMax)).begin());
  NumericMatrix subjectsPerStage(gMax + 1, kMax, ((NumericVector) rep(NA_REAL,(gMax + 1) * kMax)).begin());
  NumericMatrix testStatistics(gMax, kMax, ((NumericVector) rep(NA_REAL,gMax * kMax)).begin());
  NumericMatrix overallTestStatistics(gMax, kMax, ((NumericVector) rep(NA_REAL,gMax * kMax)).begin());
  NumericMatrix separatePValues(gMax, kMax, ((NumericVector) rep(NA_REAL,gMax * kMax)).begin());
  NumericVector conditionalCriticalValue(kMax - 1, NA_REAL);
  NumericVector conditionalPowerPerStage(kMax, NA_REAL);
  LogicalMatrix selectedArms(gMax + 1, kMax, ((LogicalVector) rep(false,(gMax + 1) * kMax)).begin());
  std::fill(selectedArms.column(0).begin(), selectedArms.column(0).end(), true);
  NumericVector adjustedPValues(kMax, NA_REAL);
  
  std::string typeOfDesign = Rcpp::as<std::string>(((CharacterVector) design.attr("class"))[0]);
  
  Function getWeightsFisher(".getWeightsFisher");
  Function getWeightsInverseNormal(".getWeightsInverseNormal");
        
  NumericVector weights(kMax);
  if (typeOfDesign == "TrialDesignFisher") {
    weights = getWeightsFisher(design);
  } else if (typeOfDesign == "TrialDesignInverseNormal") {
    weights = getWeightsInverseNormal(design);
  }
  for(int k = 0; k <= kMax - 1; k++) {
    if (k == 0) {
        subjectsPerStage(gMax, k) = plannedSubjects[k] / allocationRatioPlanned[k];
    } else {
        subjectsPerStage(gMax, k) = (plannedSubjects[k] - plannedSubjects[k - 1]) / allocationRatioPlanned[k];
    }
    if (subjectsPerStage(gMax, k) > 0) {
        simMeans(gMax, k) = R::rnorm(0, stDev / sqrt(subjectsPerStage(gMax, k)));
    }
    
      for (int treatmentArm = 0; treatmentArm < gMax; treatmentArm++) {
        if (selectedArms(treatmentArm, k)) {
          if (k == 0) {
            subjectsPerStage(treatmentArm, k) = plannedSubjects[k];
          } else {
            subjectsPerStage(treatmentArm, k) = plannedSubjects[k] - plannedSubjects[k - 1];
          }
          if (subjectsPerStage(treatmentArm, k) > 0) {
            simMeans(treatmentArm, k) = R::rnorm(muVector[treatmentArm], stDev / sqrt(subjectsPerStage(treatmentArm, k)));
            testStatistics(treatmentArm, k) = (simMeans(treatmentArm, k) - simMeans(gMax, k)) / (stDev * sqrt(1 / subjectsPerStage(treatmentArm, k) + 1 / subjectsPerStage(gMax, k)));
          }
          double vector1 = sum(((NumericVector) subjectsPerStage(treatmentArm, _ ))[Range(0, k)] * ((NumericVector) simMeans(treatmentArm, _ ))[Range(0, k)]) / sum(((NumericVector) subjectsPerStage(treatmentArm, _ ))[Range(0, k)]);
          double vector2 = sum(((NumericVector) subjectsPerStage(gMax, _ ))[Range(0, k)] * ((NumericVector) simMeans(gMax, _ ))[Range(0, k)]) / sum(((NumericVector) subjectsPerStage(gMax, _ ))[Range(0, k)]);
          overallEffects(treatmentArm, k) = vector1 - vector2;
          overallTestStatistics(treatmentArm, k) = overallEffects(treatmentArm, k) / (stDev * sqrt(1 / sum(((NumericVector) subjectsPerStage(treatmentArm, _ ))[Range(0, k)]) + 1 / sum(((NumericVector) subjectsPerStage(gMax, _ ))[Range(0, k)])));
          separatePValues(treatmentArm, k) = 1 - R::pnorm(testStatistics(treatmentArm, k), 0, 1, 1, 0);
        }
      }
      if ((k + 1) < kMax) {
        if (colSums(selectedArms)[k] == 1) {
            break;
        }
        // Bonferroni adjustment
        if(na_omit(separatePValues( _ , k)).length() == 0) {
          adjustedPValues[k] = 1 - 1e-07;
        } else {
          adjustedPValues[k] = std::min(min(na_omit(separatePValues( _ , k))) * (colSums(selectedArms)[k] - 1), 1 - 1e-07);
        }
          
          
        // conditional critical value to reject the null hypotheses at the next stage of the trial
        if (typeOfDesign == "TrialDesignConditionalDunnett") {
          conditionalCriticalValue[k] = (getOneMinusQNorm(design["alpha"]) - getOneMinusQNorm(adjustedPValues[k]) * sqrt((double) design["informationAtInterim"])) / sqrt(1 - ((double) design["informationAtInterim"]));
        } else {
          if (typeOfDesign == "TrialDesignFisher") {
            conditionalCriticalValue[k] = getOneMinusQNorm(std::min(pow(((NumericVector) design["criticalValues"])[k + 1] / ((NumericVector) cumprod((NumericVector) vectorPow(adjustedPValues[Range(0, k)],  weights[Range(0, k)])))[k], (int) 1 / weights[k + 1]), 1 - 1e-07));
          } else {
            conditionalCriticalValue[k] = (((NumericVector) design["criticalValues"])[k + 1] * sqrt(((NumericVector) design["informationRates"])[k + 1]) - sum(sapply(adjustedPValues[Range(0, k)], [&](double x) {return getOneMinusQNorm(x);}) * weights[Range(0, k)])) / sqrt((double) ((NumericVector) design["informationRates"])[k + 1] - ((NumericVector) design["informationRates"])[k]);
          }
        }
        Function selectTreatmentArms(".selectTreatmentArms");
        if (adaptations[k]) {
          if (effectMeasure == "testStatistic") {
            selectedArms( _ , k + 1) = selectedArms( _ , k) & (LogicalVector) selectTreatmentArms(k, overallTestStatistics( _ , k), typeOfSelection, epsilonValue, rValue, threshold, selectArmsFunction);
          } else if (effectMeasure == "effectEstimate") {
            selectedArms( _ , k + 1) = selectedArms( _ , k) & (LogicalVector) selectTreatmentArms(k, overallEffects( _ , k), typeOfSelection, epsilonValue, rValue, threshold, selectArmsFunction);
          }
          double newSubjects = Rf_asReal(calcSubjectsFunction(
            _["stage"] = k + 2, // to be consistent with non-multiarm situation, cf. line 37
            _["conditionalPower"] = conditionalPower,
            _["conditionalCriticalValue"] = conditionalCriticalValue,
            _["plannedSubjects"] = plannedSubjects,
            _["allocationRatioPlanned"] = allocationRatioPlanned,
            _["selectedArms"] = selectedArms,
            _["thetaH1"] = thetaH1,
            _["stDevH1"] = stDevH1,
            _["overallEffects"] = overallEffects,
            _["minNumberOfSubjectsPerStage"] = minNumberOfSubjectsPerStage,
            _["maxNumberOfSubjectsPerStage"] = maxNumberOfSubjectsPerStage
          ));
  
          if (!NumericVector::is_na(conditionalPower) || calcSubjectsFunctionIsUserDefined) {
            plannedSubjects[Range(k + 1, kMax - 1)] = sum(((NumericVector) subjectsPerStage(gMax, _ ))[Range(0, k)] * allocationRatioPlanned[Range(0, k)]) + ((NumericVector) cumsum(rep(newSubjects, kMax - (k + 1))));
          }
        } else {
            selectedArms( _ , k + 1) = selectedArms( _ , k);
        }
        double thetaStandardized;
        if (NumericVector::is_na(thetaH1)) {
          NumericVector vector1 = ((NumericVector) overallEffects( _ , k))[((LogicalVector) selectedArms( _ , k))[Range(0, gMax - 1)]];
          if(na_omit(vector1).length() == 0) {
            thetaStandardized = R_PosInf;
          } else {
            thetaStandardized = std::max(min(na_omit(vector1 / stDevH1)), 1e-12);
          }
              
        } else {
          thetaStandardized = thetaH1 / stDevH1;
        }
        conditionalPowerPerStage[k] = 1 - R::pnorm(conditionalCriticalValue[k] - thetaStandardized * sqrt(plannedSubjects[k + 1] - plannedSubjects[k]) * sqrt(1 / (1 + allocationRatioPlanned[k])), 0, 1, 1, 0);
    }
  }
  return List::create(
    _["subjectsPerStage"] = subjectsPerStage,
    _["allocationRatioPlanned"] = allocationRatioPlanned,
    _["overallEffects"] = overallEffects,
    _["testStatistics"] = testStatistics,
    _["overallTestStatistics"] = overallTestStatistics,
    _["separatePValues"] = separatePValues,
    _["conditionalCriticalValue"] = conditionalCriticalValue,
    _["conditionalPowerPerStage"] = conditionalPowerPerStage,
    _["selectedArms"] = selectedArms
  );
}


// [[Rcpp::export(name = ".getSimulationMultiArmMeansLoopsCpp")]]
List getSimulationMultiArmMeansLoopsCpp(int maxNumberOfIterations,
                                        Environment design,
                                        NumericMatrix effectMatrix,
                                        double stDev,
                                        NumericVector plannedSubjects,
                                        std::string typeOfSelection,
                                        std::string effectMeasure,
                                        LogicalVector adaptations,
                                        double epsilonValue,
                                        double rValue,
                                        double threshold,
                                        NumericVector allocationRatioPlanned,
                                        NumericVector minNumberOfSubjectsPerStage,
                                        NumericVector maxNumberOfSubjectsPerStage,
                                        double conditionalPower,
                                        double thetaH1,
                                        double stDevH1,
                                        Function calcSubjectsFunction,
                                        bool calcSubjectsFunctionIsUserDefined,
                                        Nullable<Function> selectArmsFunction,
                                        NumericMatrix indices,
                                        std::string intersectionTest,
                                        std::string successCriterion,
                                        int gMax,
                                        int kMax,
                                        NumericVector criticalValuesDunnett,
                                        NumericVector muMaxVector) {
  std::string typeOfDesign = Rcpp::as<std::string>(((CharacterVector) design.attr("class"))[0]);
  
  int cols = muMaxVector.length();

  NumericVector simulatedSelections(Dimension(kMax, cols, gMax + 1));
  NumericVector simulatedRejections(Dimension(kMax, cols, gMax));
  NumericVector simulatedSubjectsPerStage(Dimension(kMax, cols, gMax + 1));
  
  NumericMatrix simulatedNumberOfActiveArms(kMax, cols);
  NumericMatrix simulatedSuccessStopping(kMax, cols);
  NumericMatrix simulatedFutilityStopping(kMax - 1, cols);
  NumericMatrix simulatedConditionalPower(kMax, cols);
  NumericVector simulatedRejectAtLeastOne(cols, 0);
  NumericVector expectedNumberOfSubjects(cols, 0);
  NumericMatrix iterations(kMax, cols);
        
  int len = maxNumberOfIterations * kMax * gMax * cols;
      
  NumericVector dataIterationNumber(len, NA_REAL);
  NumericVector dataStageNumber(len, NA_REAL);
  NumericVector dataArmNumber(len, NA_REAL);
  NumericVector dataAlternative(len, NA_REAL);
  NumericVector dataEffect(len, NA_REAL);
  NumericVector dataSubjectsControlArm(len, NA_REAL);
  NumericVector dataSubjectsActiveArm(len, NA_REAL);
  NumericVector dataNumberOfSubjects(len, NA_REAL);
  NumericVector dataNumberOfCumulatedSubjects(len, NA_REAL);
  NumericVector dataRejectPerStage(len, NA_REAL);
  NumericVector dataSuccessStop(len, NA_REAL);
  LogicalVector dataFutilityStop(len, NA_REAL);
  NumericVector dataTestStatistics(len, NA_REAL);
  NumericVector dataConditionalCriticalValue(len, NA_REAL);
  NumericVector dataConditionalPowerAchieved(len, NA_REAL);
  NumericVector dataEffectEstimate(len, NA_REAL);
  NumericVector dataPValuesSeparate(len, NA_REAL);
  
  if(is_na(stDevH1)) {
    stDevH1 = stDev;
  }
  
  int index = 0;
  for(int i = 0; i < cols; i++) {
    for(int j = 0; i < maxNumberOfIterations; j++) {
      List stageResults = getSimulatedStageMeansMultiArmCpp(
                              design = design,
                              muVector = effectMatrix[i, ],
                              stDev = stDev,
                              plannedSubjects = plannedSubjects,
                              typeOfSelection = typeOfSelection,
                              effectMeasure = effectMeasure,
                              adaptations = adaptations,
                              epsilonValue = epsilonValue,
                              rValue = rValue,
                              threshold = threshold,
                              allocationRatioPlanned = allocationRatioPlanned,
                              minNumberOfSubjectsPerStage = minNumberOfSubjectsPerStage,
                              maxNumberOfSubjectsPerStage = maxNumberOfSubjectsPerStage,
                              conditionalPower = conditionalPower,
                              thetaH1 = thetaH1,
                              stDevH1 = stDevH1,
                              calcSubjectsFunction = calcSubjectsFunction,
                              calcSubjectsFunctionIsUserDefined = calcSubjectsFunctionIsUserDefined,
                              selectArmsFunction = selectArmsFunction
      );
      
      List closedTest;
      Function performClosedConditionalDunnettTestForSimulation(".performClosedConditionalDunnettTestForSimulation");
      Function performClosedCombinationTestForSimulationMultiArm(".performClosedCombinationTestForSimulationMultiArm");
        
      if(typeOfDesign == "TrialDesignConditionalDunnett") {
        closedTest = performClosedConditionalDunnettTestForSimulation(
                    stageResults = stageResults,
                    design = design, indices = indices,
                    criticalValuesDunnett = criticalValuesDunnett, successCriterion = successCriterion);
      } else {
        closedTest = performClosedCombinationTestForSimulationMultiArm(
                    stageResults = stageResults,
                    design = design, indices = indices,
                    intersectionTest = intersectionTest, successCriterion = successCriterion);
      }
      
      bool rejectAtSomeStage = false;
      LogicalVector rejectedArmsBefore(gMax, false);
    
      for(int k = 0; k < kMax; k++) {
      
        //simulatedRejections[k, i, ] = simulatedRejections[k, i, ] + (as<LogicalMatrix>(closedTest["rejected"])(_, k) & as<NumericMatrix>(closedTest["selectedArms"])(Range(0, gMax - 1), k) | rejectedArmsBefore); //TODO
        //simulatedSelections[k, i, ] = simulatedSelections[k, i, ] + as<NumericMatrix>(closedTest["selectedArms"])(_, k); //TODO
        for(int idx = 0; idx < gMax; idx++) {//TODO sapply?
          simulatedRejections[k + i * kMax + idx * kMax * cols] = simulatedRejections[k + i * kMax + idx * kMax * cols] + (as<LogicalMatrix>(closedTest["rejected"])(idx, k) & as<NumericMatrix>(closedTest["selectedArms"])(idx, k) | rejectedArmsBefore);
          simulatedRejections[k + i * kMax + idx * kMax * cols] = simulatedRejections[k + i * kMax + idx * kMax * cols] + as<NumericMatrix>(closedTest["selectedArms"])(idx, k);
        }

        simulatedNumberOfActiveArms(k, i) = simulatedNumberOfActiveArms(k, i) + sum(as<NumericMatrix>(closedTest["selectedArms"])(_, k));
        if(!any(is_na(as<LogicalVector>(closedTest["successStop"])))) {
          simulatedSuccessStopping(k, i) = simulatedSuccessStopping(k, i) + as<LogicalVector>(closedTest["successStop"])[k];
        }
      
        if((kMax > 1) && ((k + 1) < kMax)) {
          if(!any(is_na(as<LogicalVector>(closedTest["futilityStop"])))) {
            simulatedFutilityStopping(k, i) = simulatedFutilityStopping(k, i) + (as<LogicalVector>(closedTest["futilityStop"])[k] && !as<LogicalVector>(closedTest["successStop"])[k]); //TODO r sugar +
          }
          if(!as<LogicalVector>(closedTest["successStop"])[k] && !as<LogicalVector>(closedTest["futilityStop"])[k]) {
            simulatedConditionalPower(k + 1, i) = simulatedConditionalPower(k + 1, i) + as<NumericVector>(stageResults["conditionalPowerPerStage"])[k];
          }
        }
      
        iterations(k, i) = iterations(k, i) + 1;

        for(int g = 0; g < gMax + 1; g++) {
          if(!is_na(as<NumericMatrix>(stageResults["subjectsPerStage"])(g, k))) {
            simulatedSubjectsPerStage[k + i * kMax + g * kMax * cols] = simulatedSubjectsPerStage[k + i * kMax + g * kMax * cols] + as<NumericMatrix>(stageResults["subjectsPerStage"])(g, k);
          }
        }
      
        for(int g = 0; g < gMax; g++) {
          dataIterationNumber[index] = j;
          dataStageNumber[index] = k;
          dataArmNumber[index] = g;
          dataAlternative[index] = muMaxVector[i];
          dataEffect[index] = effectMatrix(i, g);
          dataSubjectsControlArm[index] = round(as<NumericMatrix>(stageResults["subjectsPerStage"])(gMax + 1, k), 1);
          dataSubjectsActiveArm[index] = round(as<NumericMatrix>(stageResults["subjectsPerStage"])(g, k), 1);
          dataNumberOfSubjects[index] = round(sum(na_omit(as<NumericMatrix>(stageResults["subjectsPerStage"])(_ , k))), 1);
          dataNumberOfCumulatedSubjects[index] = round(sum(na_omit(as<NumericMatrix>(stageResults["subjectsPerStage"])(_, Range(0, k)))), 1);
          dataRejectPerStage[index] = as<LogicalMatrix>(closedTest["rejected"])(g, k);
          dataTestStatistics[index] = as<NumericMatrix>(stageResults["testStatistics"])(g, k);
          dataSuccessStop[index] = as<LogicalVector>(closedTest["successStop"])[k];
          if((k + 1) < kMax) {
            dataFutilityStop[index] = as<LogicalVector>(closedTest["futilityStop"])[k];
            dataConditionalCriticalValue[index] = as<NumericVector>(stageResults["conditionalCriticalValue"])[k];
            dataConditionalPowerAchieved[index + 1] = as<NumericVector>(stageResults["conditionalPowerPerStage"])[k];
          }
          dataEffectEstimate[index] = as<NumericMatrix>(stageResults["overallEffects"])(g, k);
          dataPValuesSeparate[index] = as<NumericMatrix>(closedTest["separatePValues"])(g, k);
          index++;
        }
      
        if(!rejectAtSomeStage && any(as<LogicalMatrix>(closedTest["rejected"])(_, k) & as<LogicalMatrix>(closedTest["selectedArms"])(Range(0, gMax - 1), k) | rejectedArmsBefore)) {
          simulatedRejectAtLeastOne[i] = simulatedRejectAtLeastOne[i] + 1;
          rejectAtSomeStage = true;
        }
      
        if(((k + 1) < kMax) && (as<LogicalVector>(closedTest["successStop"])[k] || as<LogicalVector>(closedTest["futilityStop"])[k])) {
          // rejected hypotheses remain rejected also in case of early stopping
          //TODO simulatedRejections[(k + 1):kMax, i, ] <- simulatedRejections[(k + 1):kMax, i, ] + matrix((closedTest$rejected[, k] & closedTest$selectedArms[1:gMax, k] | rejectedArmsBefore), kMax - k, gMax, byrow = TRUE)
          for(int idx = 0; idx < gMax; idx++) {//TODO sapply?
            simulatedRejections[Range(k + 1 + (i * kMax + idx * kMax * cols), kMax - 1 + (i * kMax + idx * kMax * cols))] = simulatedRejections[Range(k + 1 + (i * kMax + idx * kMax * cols), kMax - 1 + (i * kMax + idx * kMax * cols))] + ;
          }
          break;
        }
      
        rejectedArmsBefore = as<LogicalMatrix>(closedTest["rejected"])(_, k) & as<LogicalMatrix>(closedTest["selectedArms"])(Range(0, gMax - 1), k) | rejectedArmsBefore;
      }
    }
    
    simulatedSubjectsPerStage = sapply(simulatedSubjectsPerStage,[&](double x){return NumericVector::is_na(x) ? 0 : x;});

    //TODO simulatedSubjectsPerStage[, i, ] <- simulatedSubjectsPerStage[, i, ] / iterations[, i]
    //for(int idx = 0; idx < gMax; idx++) {//TODO sapply?
    //  simulatedSubjectsPerStage[, i, ] <- simulatedSubjectsPerStage[, i, ] / iterations[, i];
    //}
    
    if(kMax > 1) {
      //TODO simulatedRejections[2:kMax, i, ] <- simulatedRejections[2:kMax, i, ] - simulatedRejections[1:(kMax - 1), i, ]
      NumericVector stopping = cumsum(simulatedSuccessStopping(Range(0, kMax - 2), i) + simulatedFutilityStopping(_, i)) / maxNumberOfIterations;
      //TODO  expectedNumberOfSubjects[i] <- sum(simulatedSubjectsPerStage[1, i, ] + t(1 - stopping) %*% simulatedSubjectsPerStage[2:kMax, i, ])
    } else {
      //TODO expectedNumberOfSubjects[i] <- sum(simulatedSubjectsPerStage[1, i, ])
      //expectedNumberOfSubjects[i] = sum(simulatedSubjectsPerStage[Range(0, gMax) * kMax * cols + i * kMax]);
    }
  }
  return List::create();
}