#!/usr/bin/env ruby

###########################################################################
# This script controls the end to end flow of the navigation correction and
# quality control programs as part of the Airborne Radar QC package.
# The input should be a directory that contains dorade files or a file that
# contains the list of files.
# The output will be dorade, sweep, or netcdf files that are corrected and
# quality controlled.
# Copyright 2011, Michael Bell and Cory Wolff
###########################################################################

if (ARGV.size < 1)
  puts "Usage: RunAirborneRadarQC.rb </path/to/dorade/files>\n"
  exit
end

indir = ARGV[0]

#------------------------------------------------------------------------
# 1. Run RadxConvert to generate cfRadial files for navigation correction
#------------------------------------------------------------------------
`rm ./cfradial/*`
command = `./bin/RadxConvert -d -params RadxConvert.params -f #{indir}/*`
puts command

#-------------------------------------
# 2. Generate the list of cfRadial files
#-------------------------------------
numfiles = `ls -1 ./cfradial | wc -l`
files = `ls -1 ./cfradial`
numfiles.to_i
filelist = File.new('./cfradial/filelist', 'w')
filelist.puts(numfiles)
filelist.puts(files)
filelist.close

#--------------------------------------------
# 3. Run the preprocessor on the netcdf files
#--------------------------------------------
command = `cd ./cfradial; ../bin/readnetcdf_DBZ_VR filelist; cd ../`
puts command
`rm ./cfradial/*.txt`
`rm ./cfradial/filelist`

#--------------------------------------------
# 4. Modify the navigation input file
#--------------------------------------------
template = File.open('DATA_cns_template', 'r')
data = File.new('DATA_cns_run', 'w')
template.each { |line|
   if line =~ /NUMFILES/
	data.puts(numfiles.to_s)
   else
	data.puts(line)
   end
}
data.close
template.close

#---------------------------------
# 5. Run the Navigation Correction
#---------------------------------
command = `./bin/cns_eldo_cai DATA_cns_run`
puts command

#-------------------
# 5. Run the QC code
# (must 1st apply navigation corrections to the files before running; should this be done as
#  a separate program or in the QC code?)
#-------------------
`rm ./qced/*`
load 'defaultQC.rb'

#-------------------
# 6. Run RadxConvert
#-------------------
command = `./bin/RadxConvert -d -params RadxConvert.params -f ./qced/*`
puts command

#--------------------------------------
# Output is cfRadial files that are
# ready for Dual Doppler analysis.
#--------------------------------------
