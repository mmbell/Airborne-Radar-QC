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

class RadarQCDynamicMapperMatch
  attr_accessor :attribute
  def initialize(method_sym)
    if method_sym.to_s =~ /^map_(.*)$/
      @attribute = $1.to_sym
    end
  end

  def match?
    @attribute != nil
  end
end

class RadarQC < QCscript

  def self.method_missing(method_sym, *arguments, &block)
    match = RadarQCDynamicMapperMatch.new(method_sym)
    if match.match?
      map(match.attribute => arguments.first)
    else
      super
    end
  end

  def self.respond_to?(method_sym, include_private = false)
    if RadarQCDynamicMapperMatch.new(method_sym).match?
      true
    else
      super
    end
  end

  # Syntax: defineMembership :in => 'VG', :from => '0, :to => '1', :as => 'linear', :slope => '1', :intercept => '0' 
  def map(args)
    if (args[:field] == "NCP")
      ncp = args[:value].to_f
      if (ncp < 0.2)
      	0.0;
      elsif (ncp > 0.4)
      	1.0;
      else
      	((ncp-0.2)/0.2);
      end
    end
  end

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

  # Syntax: flagSurfaceGates :with_beamwidth => '2.0' :in => 'GG', [optional] :using => 'ASTGTM2_N46E008_dem.tif'
  def flagSurfaceGates(args)
    if (args[:using].nil?)
      probGroundGates("ZZ", args[:in], args[:with_beamwidth].to_f, '');
    else 
      probGroundGates("ZZ", args[:in], args[:with_beamwidth].to_f, args[:using]);
    end
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

  # Syntax: newField :from => 'VG', :to => 'DBZ'
  def copyField(args)
     super(args[:from], args[:to])
  end

  # Syntax: wxProbability :of => "VG", in => "WXP", :gradient => 1.0, :laplacian => 1.0, :mixedpartial => 1.0,
  # :ncp => 1.0, :stddev => 1.0, :ground => 1.0, :swz => 1.0
  def wxProbability(args) 
    weight = Array.new(7)
    weight[0] = args[:gradient].to_f
    weight[1] = args[:laplacian].to_f
    weight[2] = args[:mixedpartial].to_f
    weight[3] = args[:ncp].to_f
    weight[4] = args[:stddev].to_f
    weight[5] = args[:ground].to_f
    weight[6] = args[:swz].to_f
    super(args[:of],args[:in], weight)
  end

end
