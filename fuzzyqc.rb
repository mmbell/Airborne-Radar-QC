require 'radarqc'

input = "./raw"
output = "./qced"

qc = RadarQC.new
qc.setInputPath(input)
qc.setOutputPath(output)
puts "Processing data from #{qc.getInputPath} to #{qc.getOutputPath}"
qc.getFileListSize.times { |file|
	puts "Processing file #{file}"
	qc.instance_eval do
	   load file
     #removeAircraftMotion :from => "VR", :to => "VG"
     wxProbability :of => "VV", :in => "WXP", :gradient => 1.0, :laplacian => 10.0, :mixedpartial => 1.0, \
       :ncp => 1.0, :stddev => 1.0, :ground => 10.0, :swz => 1.0
	   #copyEdits :from => "VV", :to => "DBZ"
	   save file
	end
}
puts 'QC complete'

