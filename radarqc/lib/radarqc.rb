#
# RadarQC class
#
# This is a Ruby extension of the QCscript and AirborneRadarQC c++ classes.
# Pure Ruby methods can be declared here that make the script idioms more natural
# or additional computational or logic methods the build upon the c++
# methods that interface directly with the radar data
#
# Author::	Michael Bell
# Copyright::	Copyright (c) 2011, Michael Bell and Cory Wolff
# License::	MIT

require 'QCscript' 

# This is a Ruby extension of the QCscript and AirborneRadarQC c++ classes
# Pure Ruby methods can be declared here that make the script idioms more natural
# or additional computational or logic methods the build upon the c++
# methods that interface directly with the radar data

class RadarQC < QCscript

  # Syntax: threshold :field => 'VG' :on => 'NCP' :below => '0.2'
  def threshold(args)
     if args.has_key?(:above)
	thresholdData(args[:on], args[:field], args[:above].to_f, "above")
     else
	# Assume its below
        thresholdData(args[:on], args[:field], args[:below].to_f, "below")
     end
  end

  # Syntax: applyNavCorrections :to => 'radar', :using => 'file'
  def applyNavigationCorrections(args)
     setNavigationCorrections(args[:using], args[:to])
  end

  # Syntax: removeAircraftMotion :from => 'VR', :to => 'VG'
  def removeAircraftMotion(args)
     super(args[:from],args[:to])
  end

  # Syntax: flagSurfaceGates :with_beamwidth => '2.0' :in => 'GG'
  def flagSurfaceGates(args)
     probGroundGates("ZZ", args[:in], args[:with_beamwidth].to_f);
  end

  # Syntax: calcRatio :of => 'SW', :over => 'ZZ', :in => 'SWZ', [optional] :using => 'linear_z'
  def calcRatio(args)
     if (args[:using] == 'linear_z')
	super(args[:of], args[:over], args[:in], true);
     else
	super(args[:of], args[:over], args[:in], false);
     end
  end

 # Syntax: despeckle :gates => '5', :in => 'VG', :along => 'radial' (or 'azimuthal')
  def despeckle(args)
     if (args[:along] == 'radial')
	despeckleRadial(args[:in], args[:gates].to_i)
     else
	despeckleAzimuthal(args[:in], args[:gates].to_i)
     end
  end

  # Syntax: copyEdits :from => 'VG', :to => 'DBZ'
  def copyEdits(args)
     super(args[:from], args[:to])
  end

end
