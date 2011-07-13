#!/usr/bin/ruby

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

@indir = ARGV[0]

#-------------------------------------
# 1. Generate the list of dorade files
#-------------------------------------
command = `ls -1 #@indir > .filelist`
puts command

#-------------------
# 2. Run RadxConvert
#-------------------
command = `./RadxConvert -params RadxConvert.params -paths #@indir`
puts command

#--------------------------------------------
# 3. Run the preprocessor on the netcdf files
#--------------------------------------------



#---------------------------------
# 4. Run the Navigation Correction
#---------------------------------


#-------------------
# 5. Run the QC code
# (must 1st apply navigation corrections to the files before running; should this be done as
#  a separate program or in the QC code?)
#-------------------
#command = `./AirborneRadarQC #@indir ./qced`
#puts command
QC = AirborneRadarQC.new
QC.processSweeps

#-------------------
# 6. Run RadxConvert
#-------------------
command = `./RadxConvert -params RadxConvert.params -paths ./qced`
puts command

#--------------------------------------
# Output is cfRadial files that are
# ready for Dual Doppler analysis.
#--------------------------------------
