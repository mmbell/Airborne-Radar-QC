require './QCscript'

input = "/Users/mbell/Science/current/eldora/autoqc/test/ori"
output = "/Users/mbell/Science/current/eldora/autoqc/test/qced"
ta_cfac = "/Users/mbell/Science/current/eldora/autoqc/test/rf12.cfac.aft"
tf_cfac = "/Users/mbell/Science/current/eldora/autoqc/test/rf12.cfac.fore"

qc = QCscript.new
qc.setInputPath(input)
qc.setOutputPath(output)
puts "Processing data from #{qc.getInputPath} to #{qc.getOutputPath}"
qc.getFileListSize.times { |file|
	puts "Processing file #{file}"
	qc.instance_eval do
	   load file

	   setNavigationCorrections(ta_cfac, "TA-ELDR")
           setNavigationCorrections(tf_cfac, "TF-ELDR")
	   removeAircraftMotion("VR", "VQC");

	   thresholdData("NCP", "ZZ", 0.2, "below")
	   thresholdData("NCP","VQC", 0.2, "below")

	   probGroundGates("ZZ", "GG", 2.0);
	   thresholdData("GG","ZZ", 0.7, "above");
	   thresholdData("GG","VQC", 0.7, "above");
		
	   calcRatio("SW", "ZZ", "SWZ", true);
	   thresholdData("SWZ","ZZ", 0.6, "above");
	   thresholdData("SWZ","VQC", 0.6, "above");

	   despeckleRadial("ZZ", 2);
	   despeckleRadial("VQC", 2);
	   despeckleAzimuthal("ZZ", 2);
	   despeckleAzimuthal("VQC", 2);

	   save file
	end
}
puts 'QC complete'

