/* AirborneRadarQC */
/* Copyright 2011 Michael Bell and Cory Wolff */
/* All rights reserved */

#include <iostream>
#include <QApplication>
#include "radarqc/ext/QCscript/Dorade.h"
#include "radarqc/ext/QCscript/AirborneRadarQC.h"

using namespace std;

int main (int argc, char *argv[]) {

	// Get the arguments
	if (argc < 3) {
		// Eventually, no arguments would start an interactive GUI mode
		//QApplication app(argc, argv);
		cout << "Usage: eldoraqc /path/to/sweepfiles /path/to/output\n";
		exit(1);
	}
	
	// The qc object will read from one directory and write to another
	QString inpath = argv[1];
	QString outpath = argv[2];	
	QString suffix = "QC";
	AirborneRadarQC QC(inpath, outpath, suffix);
	
	// Process the data
	QC.processSweeps();
	
	return 0;
	
}

