/*
 *  QCscript.cpp
 *
 *  Copyright 2011 Michael Bell and Cory Wolff
 *  All rights reserved.
 *
 */

#include "rice/Data_Type.hpp"
#include "rice/Constructor.hpp"
#include "rice/Array.hpp"
#include "QCscript.h"
#include "AirborneRadarQC.h"
#include <QString>
#include <iostream>

using namespace Rice;

extern "C"
void Init_QCscript(void)
{
  RUBY_TRY
  {

    define_class<AirborneRadarQC>("AirborneRadarQC");

    define_class<QCscript, AirborneRadarQC>("QCscript")
      .define_constructor(Constructor<QCscript>())
      .define_method("processSweeps", &QCscript::processSweeps)
      .define_method("setInputPath", &QCscript::rb_setInputPath)
      .define_method("setOutputPath", &QCscript::rb_setOutputPath)
      .define_method("getInputPath", &QCscript::rb_getInputPath)
      .define_method("getOutputPath", &QCscript::rb_getOutputPath)
      .define_method("getFileListSize", &QCscript::getfileListsize)
      .define_method("load", &QCscript::rb_load)
      .define_method("save", &QCscript::rb_saveQCedSwp)
      .define_method("setNavigationCorrections", &QCscript::rb_setNavigationCorrections)
      .define_method("removeAircraftMotion", &QCscript::rb_removeAircraftMotion)
      .define_method("thresholdData", &QCscript::rb_thresholdData)
      .define_method("probGroundGates", &QCscript::rb_probGroundGates)
      .define_method("calcRatio", &QCscript::rb_calcRatio)
      .define_method("despeckleRadial", &QCscript::rb_despeckleRadial)
      .define_method("despeckleAzimuthal", &QCscript::rb_despeckleAzimuthal)
      .define_method("copyEdits", &QCscript::rb_copyEdits)
      .define_method("copyField", &QCscript::rb_copyField)
	  .define_method("wxProbability", &QCscript::rb_wxProbability);
  }
  RUBY_CATCH
} 

QCscript::QCscript()
   : AirborneRadarQC()
{

}

QCscript::~QCscript()
{
}


bool QCscript::rb_setInputPath(const std::string& in)
{
	QString path = QString::fromStdString(in);
	return setInputPath(path);
}

bool QCscript::rb_setOutputPath(const std::string& out)
{
        QString path = QString::fromStdString(out);
        return setOutputPath(path);
}

std::string QCscript::rb_getInputPath()
{
        return getInputPath().toStdString();
}

std::string QCscript::rb_getOutputPath()
{
        return getOutputPath().toStdString();
}

bool QCscript::rb_load(int swpIndex)
{
	if(load(swpIndex)) return true;
	return false;
}

bool QCscript::rb_saveQCedSwp(int swpIndex)
{
        if(saveQCedSwp(swpIndex)) return true;
        return false;
}

void QCscript::rb_thresholdData(const std::string& threshfield, const std::string& fldname,
                                float threshold, const std::string& direction)
{
        QString thr = QString::fromStdString(threshfield);
	QString fld = QString::fromStdString(fldname);
        QString dir = QString::fromStdString(direction);
        thresholdData(thr, fld, threshold, dir);
}

void QCscript::rb_probGroundGates(const std::string& oriFieldName, const std::string& newFieldName,
                                  float eff_beamwidth, const std::string& demFilename)
{
        QString ori = QString::fromStdString(oriFieldName);
        QString fld = QString::fromStdString(newFieldName);
		QString dem = QString::fromStdString(demFilename);
        probGroundGates(ori,fld,eff_beamwidth,dem);
}

void QCscript::rb_despeckleRadial(const std::string& fieldName, int speckle)
{
        QString fld = QString::fromStdString(fieldName);
        despeckleRadial(fld,speckle);
}

void QCscript::rb_despeckleAzimuthal(const std::string& fieldName, int speckle)
{
        QString fld = QString::fromStdString(fieldName);
        despeckleAzimuthal(fld,speckle);
}

void QCscript::rb_copyEdits(const std::string& oriFieldName, const std::string& newFieldName)
{
        QString ori = QString::fromStdString(oriFieldName);
        QString fld = QString::fromStdString(newFieldName);
        copyEdits(ori,fld);
}

void QCscript::rb_calcRatio(const std::string& topFieldName, const std::string& bottomFieldName,
			    const std::string& newFieldName, bool zflag)
{
        QString top = QString::fromStdString(topFieldName);
        QString bottom = QString::fromStdString(bottomFieldName);
        QString fld = QString::fromStdString(newFieldName);
        calcRatio(top, bottom, fld, zflag);
}

void QCscript::rb_setNavigationCorrections(const std::string& cfacFileName, const std::string& radarName)
{
        QString cfac = QString::fromStdString(cfacFileName);
        QString radar = QString::fromStdString(radarName);
        setNavigationCorrections(cfac, radar);
}

void QCscript::rb_removeAircraftMotion(const std::string& vrFieldName, const std::string& vgFieldName)
{
        QString vr = QString::fromStdString(vrFieldName);
		QString vg = QString::fromStdString(vgFieldName);
		removeAircraftMotion(vr,vg);
}

void QCscript::rb_copyField(const std::string& oriFieldName, const std::string& newFieldName)
{
        QString ori = QString::fromStdString(oriFieldName);
        QString fld = QString::fromStdString(newFieldName);
        copyField(ori,fld);
}

void QCscript::rb_wxProbability(const std::string& vgFieldName, const std::string& wxFieldName, const Rice::Array weight)
{
        QString vg = QString::fromStdString(vgFieldName);
		QString wx = QString::fromStdString(wxFieldName);
		c_weight = new float[7];
		for (int i=0; i< 7; i++) c_weight[i] = from_ruby<float>(weight[i]);
		wxProbability(vg,wx,c_weight);
}
