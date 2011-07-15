###########################################################################
#
# RadarQC class
#
# This is a Ruby extension of the QCscript and AirborneRadarQC c++ classes
# Pure Ruby methods can be declared here that make the script idioms more natural
# or additional computational or logic methods the build upon the c++
# methods that interface directly with the radar data
#
# Copyright 2011, Michael Bell and Cory Wolff
###########################################################################

require './lib/radarqc/QCscript' 

class RadarQC < QCscript

  def threshold(args)
     # Syntax: threshold :field => 'VG' :on => 'NCP' :below => '0.2'
     if args.has_key?(:above)
	thresholdData(args[:on], args[:field], args[:above].to_f, "above")
     else
	# Assume its below
        thresholdData(args[:on], args[:field], args[:below].to_f, "below")
     end
  end

  def applyNavigationCorrections(args)
     # Syntax: applyNavCorrections :to => 'radar', :using => 'file'
     setNavigationCorrections(args[:using], args[:to])
  end

  def removeAircraftMotion(args)
     # Syntax: removeAircraftMotion :from => 'VR', :to => 'VG'
     super(args[:from],args[:to])
  end

  def flagSurfaceGates(args)
     # Syntax: flagSurfaceGates :with_beamwidth => '2.0' :in => 'GG'
     probGroundGates("ZZ", args[:in], args[:with_beamwidth].to_f);
  end

  def calcRatio(args)
     # Syntax: calcRatio :of => 'SW', :over => 'ZZ', :in => 'SWZ', :using => 'linear_z'
     if (args[:using] == 'linear_z')
	super(args[:of], args[:over], args[:in], true);
     else
	super(args[:of], args[:over], args[:in], false);
     end
  end

  def despeckle(args)
     # Syntax: despeckle :gates => '5', :in => 'VG', :along => 'radial' (or 'azimuthal')
     if (args[:along] == 'radial')
	despeckleRadial(args[:in], args[:gates].to_i)
     else
	despeckleAzimuthal(args[:in], args[:gates].to_i)
     end
  end

  def copyEdits(args)
     # Syntax: copyEdits :from => 'VG', :to => 'DBZ'
     super(args[:from], args[:to])
  end

end
