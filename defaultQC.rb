require 'radarqc'

input = "./raw"
output = "./qced"
ta_cfac = "./corrections/cfac.aft"
tf_cfac = "./corrections/cfac.fore"

qc = RadarQC.new
qc.setInputPath(input)
qc.setOutputPath(output)
puts "Processing data from #{qc.getInputPath} to #{qc.getOutputPath}"
qc.getFileListSize.times { |file|
	puts "Processing file #{file}"
	qc.instance_eval do
	   load file
	   applyNavigationCorrections :to => "TA-ELDR", :using => ta_cfac
           applyNavigationCorrections :to => "TF-ELDR", :using => tf_cfac
	   removeAircraftMotion :from => "VR", :to => "VG"
           threshold :field => "VG", :on => "NCP", :below => "0.3"
	   flagSurfaceGates :with_beamwidth => "2.0", :in => "GG" #:using => 'ASTGTM2_N46E008_dem.tif'
	   threshold :field => "VG", :on => "GG", :above => "0.7"
           calcRatio :of => "SW", :over => "ZZ", :in => "SWZ", :using => "linear_z"
	   threshold :field => "VG", :on => "SWZ", :above => "0.6"
	   despeckle :gates => '5', :in => "VG", :along => "radial";
	   despeckle :gates => '5', :in => "VG", :along => "azimuth";
	   copyEdits :from => "VG", :to => "DBZ"
	   save file
	end
}
puts 'QC complete'

