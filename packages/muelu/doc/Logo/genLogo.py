#! /usr/bin/env python
import math
import optparse
import os
import subprocess
import sys

def getstatusoutput(cmd):
    """Return (status, output) of executing cmd in a shell."""
    pipe = os.popen(cmd + ' 2>&1', 'r')
    text = pipe.read()
    sts = pipe.close()
    if sts is None: sts = 0
    if text[-1:] == '\n': text = text[:-1]
    return sts, text


def deleteDir(path):
    """deletes the path entirely"""
    cmd = "rm -rf "+path
    result = getstatusoutput(cmd)
    if(result[0]!=0):
        raise RuntimeError(result[1])

def createDir(path):
    """deletes the path entirely"""
    cmd = "mkdir "+path
    result = getstatusoutput(cmd)
    if(result[0]!=0):
        raise RuntimeError(result[1])

def runCommand(cmd):
    """deletes the path entirely"""
    result = getstatusoutput(cmd)
    #if(result[0]!=0):
    #    raise RuntimeError(result[1])
    return result[1]

def clearWindow():
    os.system('cls' if os.name == 'nt' else 'clear')

def waitForKey():
    os.system('read -s -n 1 -p "Press any key to continue..."')
    print

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# some colors
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    OKDARKGREEN = '\033[32m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def disable(self):
    self.HEADER = ''
    self.OKBLUE = ''
    self.OKGREEN = ''
    self.OKDARKGREEN = ''
    self.WARNING = ''
    self.FAIL = ''
    self.ENDC = ''

def print_part_output(filenum,filenameformat,line):
    filename = filenameformat.replace('#',str(filenum) )
    fout=open(filename,'w')
    fout.write(line)
    fout.close

def splitfileintofragments(inputfile,splittingtxt):

  filenameformat = inputfile + '_#.fragment'
  file = open( inputfile )
  lines=file.read().split(splittingtxt)

  bLeadingLines = True

  for i in range(0,len( lines ) ):
    line = lines[i]
    if not line.strip():
      continue
    else:
      print_part_output( i+1 , filenameformat, line )

  # move fragment files to this folder
  runCommand("mv ../src/*.fragment .")


################################# MAIN routine
if __name__ == '__main__':

  if not os.path.isfile("logo.tmpl") or not os.access("logo.tmpl", os.R_OK):
    print bcolors.FAIL+" Failure: Could not find logo.tmpl"+bcolors.ENDC
    sys.exit(0)

  p = optparse.OptionParser()

  p.add_option('-w', '--width', dest="pxwidth", default=1024, type='int', \
          help='Desired logo width in pixels (note that the actual width will be slightly different/smaller)')

  p.add_option('-c', '--color', dest="color", default="Blue",
          help='Color scheme. Predefined values: "Blue", "Orange" or "Gray"')

  p.add_option('-b', '--background', dest="background", default="light",
          help='Parameter controlling the color of "Multigrid Framework" text options.' +
          'Possible values: "light", "dark". "light" means that the logo is supposed to be used with a light background')

  # parse
  options, arguments = p.parse_args()

  color = options.color

  if os.path.isfile("logo_tmpl.tex"):
    os.remove("logo_tmpl.tex")
  o = open("logo_tmpl.tex","a")
  #for line in open("tmpl/twoblocks.head_TMPL"):
  for line in open("logo.tmpl"):
    line = line.replace("$$COLOR$$", color)
    if options.background == "dark":
      line = line.replace("$$FONTCOLOR$$", "colorL!50!colorEnd")
    else:
      line = line.replace("$$FONTCOLOR$$", "colorM")
    o.write(line)
  o.close()

  print bcolors.OKGREEN + "Run lualatex... " + bcolors.ENDC
  runCommand("lualatex logo_tmpl.tex")
  runCommand("mv logo_tmpl.pdf logo.pdf")

  # calculate density resolution
  dens = options.pxwidth / 2.45

  print bcolors.OKGREEN + "Create logo on white background... " + bcolors.ENDC
  runCommand("convert -density " + str(dens) + "x" + str(dens) + " logo.pdf logo_" + color + ".png")
  runCommand("convert logo_" + color + ".png -trim +repage -bordercolor white    -border 10x10 logo_" + color + ".png")
  print bcolors.OKGREEN + "Create transparent logo... " + bcolors.ENDC
  runCommand("convert logo_" + color + ".png -fuzz 10% -transparent white logo_" + color + "_transparent.png")

  print bcolors.OKDARKGREEN + "Finished. " + bcolors.ENDC
