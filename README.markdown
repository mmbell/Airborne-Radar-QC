# Airborne-Radar-QC

This code is part of a project for meteorological data quality control of airborne tail Doppler radars.
The project is sponsored by the National Science Foundation and the National Center for Atmospheric Research.
It currently consists of:

  * radarqc -- Tools for batch editing of radar data using Ruby and C++
  * navigation -- Calculation of INS navigation corrections using Fortran

A complete quality control procedure consists of creating navigation corrections, applying them to the raw data, and removing non-meteorological radar echoes from the data.

A sample Ruby script `RunAirborneRadarQC.rb` that performs these procedures in order to a set of Dorade sweepfiles is included.

To run the script, the `RadxConvert` tool is required, which is available at the [Radx] (http://www.ral.ucar.edu/projects/titan/docs/radial_formats/radx.html) website

## RadarQC

Use rake to build and install the `radarqc` Ruby gem onto your local machine. 

     $ rake install

You will also need the `bundler` gem to install the required dependencies. Once the `radarqc` gem is installed, you can include it in an editing script. A sample `defaultQC.rb` script is included in the top level directory. The script is based on the procedure described in:

Bell, M. M., W.-C. Lee, C. Wolff, and H. Cai, 2011: An automated procedure for quality control of airborne tail Doppler radar data. *In preparation*

### Qt

The Ruby gem requires the [Qt] (http://http://qt.nokia.com/) library to build the native extension.
A command line C++ application can also be built that does not require the Ruby interface. From the top-level directory run

     $ qmake

to create a Makefile or Xcode project for your machine.

## Navigation
 
The navigation correction code is based on the algorithm described in:

Georgis, J.-F., F. Roux, and P. H. Hildebrand, 2000: Observation of precipitating systems over complex orography with meteorological Doppler radars: A feasibility study. *Meteor. Atmos. Phys.*, **72**, 185-202.

To compile, a Makefile is included

     $ make

More complete documentation on how to run the software will be described here shortly.

## Contributing to Airborne-Radar-QC

* Check out the latest master to make sure the feature hasn't been implemented or the bug hasn't been fixed yet
* Check out the [issue tracker](http://github.com/mmbell/Airborne-Radar-QC/issues) to make sure someone already hasn't requested it and/or contributed it
* Fork the project
* Start a feature/bugfix branch
* Commit and push until you are happy with your contribution
* Make sure to add tests for the feature/bugfix. This is important so I don't break it in a future version unintentionally.
* Please try not to mess with the Rakefile, version, or history. If you want to have your own version, or is otherwise necessary, that is fine, but please isolate it to its own commit so I can cherry-pick around it.

## Copyright

Radarqc Copyright (c) 2011 Michael Bell, Cory Wolff

Navigation Copyright (c) 2011 Michael Bell, Huaqing Cai, and Frank Roux.

See LICENSE for details.

