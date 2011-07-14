require './QCscript'

input = "/Users/mbell/Science/current/eldora/autoqc/test/ori"
output = "/Users/mbell/Science/current/eldora/autoqc/test/qced"
ta_cfac = "/Users/mbell/Science/current/eldora/autoqc/test/rf12.cfac.aft"
tf_cfac = "/Users/mbell/Science/current/eldora/autoqc/test/rf12.cfac.fore"

qc = QCscript.new
qc.setInputPath(input)
qc.setOutputPath(output)
puts 'Processing data from ' + qc.getInputPath + ' to ' + qc.getOutputPath
0.upto(qc.getFileListSize - 1) { |file|
	puts "Processing #{file}"
	qc.load(file)

	qc.setNavigationCorrections(ta_cfac, "TA-ELDR")
        qc.setNavigationCorrections(tf_cfac, "TF-ELDR")
	qc.removeAircraftMotion("VR", "VQC");

	qc.thresholdData("NCP", "ZZ", 0.2, "below")
	qc.thresholdData("NCP","VQC", 0.2, "below")

	qc.probGroundGates("ZZ", "GG", 2.0);
	qc.thresholdData("GG","ZZ", 0.7, "above");
	qc.thresholdData("GG","VQC", 0.7, "above");
		
	qc.calcRatio("SW", "ZZ", "SWZ", true);
	qc.thresholdData("SWZ","ZZ", 0.6, "above");
	qc.thresholdData("SWZ","VQC", 0.6, "above");

	qc.despeckleRadial("ZZ", 2);
	qc.despeckleRadial("VQC", 2);
	qc.despeckleAzimuthal("ZZ", 2);
	qc.despeckleAzimuthal("VQC", 2);

	qc.save(file)
}
puts 'QC complete'

