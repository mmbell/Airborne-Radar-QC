/*
 *  QCscript.h
 *
 *  Copyright 2011 Michael Bell and Cory Wolff
 *  All rights reserved.
 *
 */

#ifndef QCSCRIPT_H
#define QCSCRIPT_H

#include "AirborneRadarQC.h"
#include "rice/Array.hpp"

class QCscript : public AirborneRadarQC
{

public:
     QCscript();
     ~QCscript();

     // Ruby wrappers
     bool rb_setInputPath(const std::string& in);
     bool rb_setOutputPath(const std::string& out);
     std::string rb_getInputPath();
     std::string rb_getOutputPath();
     bool rb_load(int swpIndex);
     bool rb_saveQCedSwp(int swpIndex);
     void rb_thresholdData(const std::string& threshfield, const std::string& fldname,
                           float threshold, const std::string& direction = "below");
     void rb_probGroundGates(const std::string& oriFieldName, const std::string& newFieldName,
			     float eff_beamwidth, const std::string& demFilename);
     void rb_despeckleRadial(const std::string& fieldName, int speckle);
     void rb_despeckleAzimuthal(const std::string& fieldName, int speckle);
     void rb_copyEdits(const std::string& oriFieldName, const std::string& newFieldName);
     void rb_calcRatio(const std::string& topFieldName, const std::string& bottomFieldName,
                       const std::string& newFieldName, bool zflag);
     void rb_setNavigationCorrections(const std::string& cfacFileName, const std::string& radarName);
     void rb_removeAircraftMotion(const std::string& vrFieldName, const std::string& vgFieldName);
     void rb_copyField(const std::string& oriFieldName, const std::string& newFieldName);
	 void rb_wxProbability(const std::string& vgFieldName, const std::string& wxFieldName, const Rice::Array weight);

private:

	float* c_weight;

};

#endif
