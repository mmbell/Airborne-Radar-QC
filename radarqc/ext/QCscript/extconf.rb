require 'rubygems'
require 'mkmf-rice'
require 'fileutils'

$: << File.expand_path("../../../lib", __FILE__)

require 'qt/qmake'

Qt::Qmake.config('radarqc')

create_makefile('QCscript')
