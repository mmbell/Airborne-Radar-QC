/*
 *  AirborneRadarQC.h
 *
 *  Copyright 2011 Michael Bell and Cory Wolff
 *  All rights reserved.
 *
 */

#ifndef AIRBORNEQC_H
#define AIRBORNEQC_H

#include "Dorade.h"
#include <QList>
#include <QDir>

class AirborneRadarQC
{
	
public:
	// Conor / Deor
	AirborneRadarQC(const QString& in, const QString& out, const QString& suffix);
	~AirborneRadarQC();

	// I/O
	bool readSwpDir();
	bool load(const int& swpIndex);
	bool loadAuxSwp(const int& swpIndex);
	bool saveQCedSwp(const int& swpIndex);
	bool saveDorade(const QString& doradeFilename);
	int getfileListsize() { return swpfileList.size(); }
	QString getswpfileName(int n) { return swpfileList[n]; }
	
	// QC
	bool processSweeps();
	bool newField(const QString& oldFieldName,const QString& newFieldName, 
				  const QString& newFieldDesc,const QString& newFieldUnits);
	bool copyAuxField(const QString& oldFieldName,const QString& newFieldName, 
					  const QString& newFieldDesc,const QString& newFieldUnits);
	void thresholdData(const QString& threshfield, const QString& fldname, 
				const float& threshold, const QString& direction = "below");
	void despeckleRadial(const QString& fldname, const int& speckle);
	void despeckleAzimuthal(const QString& fldname, const int& speckle);
	void GaussianSmooth(const QString& oriFieldName, const QString& newFieldName, const int& scale);
	void GaussianSmooth(const QString& oriFieldName, float** field, const int& scale);
	void swpField2array(const QString& oriFieldName, float** field);
	
	// REC Fields
	void calcTexture(const QString& oriFieldName, const QString& fldname);
	void calcSpinSteiner(const QString& oriFieldName, const QString& fldname);
	void calcSpinKessinger(const QString& oriFieldName, const QString& fldname);
	void calcStdDev(const QString& oldFieldName, const QString& fldname);
	void calcMeanRef(const QString& fldname);
	void calcSpatialMean(const QString& oriFieldName, const QString& newFieldName, const int& gateWindow, const int& rayWindow);
	void calcTemporalMean(const QString& oriFieldName, const QString& newFieldName);
	void calcGate2GateRefGrad(const QString& fldname);
	void calcAzimuthRefGrad(const QString& fldname);
	
	void calcRatio(const QString& topFieldName, const QString& bottomFieldName,
				   const QString& newFieldName, const bool& zflag);
	void calcGradientMagnitude(const QString& oriFieldName, const QString& newFieldName, const int& order);
	void calcGradientMagnitude(const QString& oriFieldName, float** field, const int& order);
	void calcGradientMagnitude(float** orifield, float** field, const int& order);

	// Generic derivatives 
	// Create new swp field for display
	void calc1stAzimuthalDerivative(const QString& oriFieldName, const QString& newFieldName, const int& order);
	void calc1stRadialDerivative(const QString& oriFieldName, const QString& newFieldName, const int& order);
	void calc2ndAzimuthalDerivative(const QString& oriFieldName, const QString& newFieldName, const int& order);
	void calc2ndRadialDerivative(const QString& oriFieldName, const QString& newFieldName, const int& order);
	void calcLaplacian(const QString& oriFieldName, const QString& newFieldName);
	void calcMixedPartial(const QString& oriFieldName, const QString& newFieldName);
	
	// Create from swp field to temporary
	void calc1stAzimuthalDerivative(const QString& oriFieldName, float** field, const int& order);
	void calc1stRadialDerivative(const QString& oriFieldName, float** field, const int& order);
	void calc2ndAzimuthalDerivative(const QString& oriFieldName, float** field, const int& order);
	void calc2ndRadialDerivative(const QString& oriFieldName, float** field, const int& order);
	void calcLaplacian(const QString& oriFieldName, float** field);
	void calcMixedPartial(const QString& oriFieldName, float** field);
	
	// Create from temporary to temporary
	void calc1stAzimuthalDerivative(float** orifield, float** field, const int& order);
	void calc1stRadialDerivative(float** orifield, float** field, const int& order);
	void calc2ndAzimuthalDerivative(float** orifield, float** field, const int& order);
	void calc2ndRadialDerivative(float** orifield, float** field, const int& order);
	void calcLaplacian(float** orifield, float** field);
	void calcMixedPartial(float** orifield, float** field);
	
	// Flag sensitive areas
	void flagGroundGates(const QString& fldname, const float& eff_beamwidth);
	
	// Probabilities
	void probGroundGates(const QString& oriFieldName, const QString& newFieldName, const float& eff_beamwidth);
	void calcWeatherProb(const QString& mdbzt_name, const QString& mdbzs_name, const QString& mdbzl_name, const QString& mvgs_name, const QString& mncp_name);
	void wxProbability2();
	void mapRefTexture(const QString& fldname);
	void mapMeanRef(const QString& fldname);
	void mapRefSpin(const QString& fldname);
	void mapVelStd(const QString& fldname);
	void mapNCP(const QString& fldname);
	void mapRefLaplacian(const QString& fldname);

	// Other
	void compareForeAftRef();
	void dumpFLwind();
	void compareFLwind();
	void removeAircraftMotion(const QString& vrFieldName, const QString& vgFieldName);
	void setNavigationCorrections(const QString& cfacFileName, const QString& radarName);

	// Verification
	void BrierSkillScore();
	void RelativeOperatingCharacteristic();
	void ReliabilityDiagram();
	void verify();
	void soloiiScriptROC();
	void soloiiScriptVerification();
	
private:
	QList<QString> swpfileList;
	QDir dataPath;
	QDir outPath;
	QString swpSuffix;
	
	Dorade swpfile;
	Dorade auxSwpfile; // Used to merge or thin sweeps

	int getRayIndex(int ri, int nrays);
	float calcRefTextInterestMap(float texture);
	float calcRefSpinInterestMap(float spin);
	float calcVelStdInterestMap(float std);
	float calcMeanRefInterestMap(float ref);
	float calcNCPInterestMap(float ncp);
	float calcRefLaplacianInterestMap(float lap);
	
};

#endif
