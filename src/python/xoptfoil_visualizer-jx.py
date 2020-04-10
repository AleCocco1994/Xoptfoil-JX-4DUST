#!/usr/bin/env python

#  This file is part of XOPTFOIL.

#  XOPTFOIL is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  XOPTFOIL is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with XOPTFOIL.  If not, see <http://www.gnu.org/licenses/>.

#  Copyright (C) 2014 -- 2016 Daniel Prosser

import argparse
from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np
from math import log10, floor
from sys import version_info

# Default plottiong options

plotoptions = dict(show_seed_airfoil = True,
                   show_seed_polar = True,
                   show_seed_airfoil_only = False,
                   show_seed_polar_only = False,
                   show_airfoil_info = True,
                   plot_airfoils = True,
                   plot_polars = True,
                   plot_optimization_history = True,
                   drag_plot_type = "vs. lift",
                   save_animation_frames = False,
                   color_for_seed = "blue",
                   color_for_new_designs = "red",
                   monitor_update_interval = 10)

################################################################################
#
# Airfoil class
#
################################################################################
class Airfoil:

  def __init__(self):
    self.x = np.zeros((0))
    self.y = np.zeros((0))
    self.maxt = 0.
    self.xmaxt = 0.
    self.maxc = 0.
    self.xmaxc = 0.
    self.alpha = np.zeros((0))
    self.cl = np.zeros((0))
    self.cd = np.zeros((0))
    self.cm = np.zeros((0))
    self.xtrt = np.zeros((0))
    self.xtrb = np.zeros((0))
    self.npt = 0
    self.noper = 0
    # jx-mod additional 2nd derivation
    self.deriv2 = np.zeros((0))
    self.deriv3 = np.zeros((0))
    # jx-mod additional glide and climb ratio, falpangle
    self.glide = np.zeros((0))
    self.climb = np.zeros((0))
    self.flapangle = np.zeros((0))
    
  def setCoordinates(self, x, y):
    self.x = x
    self.y = y
    self.npt = x.shape[0]

  # jx-mod set  2nd, 3rd derivation
  def setDerivatives(self, deriv2, deriv3):
    self.deriv2 = deriv2
    self.deriv3 = deriv3
    iLE = np.argmin(self.x)
    self.deriv2[0:iLE] = np.multiply(self.deriv2[0:iLE], -1)   
    self.deriv3[0:iLE] = np.multiply(self.deriv3[0:iLE], -1) 

  def setGeometryInfo(self, maxt, xmaxt, maxc, xmaxc):
    self.maxt = maxt
    self.xmaxt = xmaxt
    self.maxc = maxc
    self.xmaxc = xmaxc

  def setPolars(self, alpha, cl, cd, cm, xtrt, xtrb, flapangle):
    self.alpha = alpha
    self.cl = cl
    self.cd = cd
    self.cm = cm
    self.xtrt = xtrt
    self.xtrb = xtrb
    self.noper = alpha.shape[0]
    # jx-mod additional glide and climb ratio
    self.glide = self.cl / self.cd
    self.climb = np.zeros((len(cl)))
    for i in range(len(cl)): 
      if cl[i] > 0.0: 
        self.climb[i] = (cl[i]**1.5)/cd[i] 
    self.flapangle = flapangle

################################################################################
# Reads airfoil coordinates from file
def read_airfoil_coordinates(filename, zonetitle, designnum):

  ioerror = 0
  x = []
  y = []
  maxt = 0.
  xmaxt = 0.
  maxc = 0.
  xmaxc = 0.
  # jx-mod additionally 2nd and 3rd derivative
  deriv2 = []
  deriv3 = []


  # Try to open the file

  try:
    f = open(filename) 
  except IOError:
    ioerror = 1
    return x, y, maxt, xmaxt, maxc, xmaxc, ioerror

  # Read lines until we get to the correct zone

  zonefound = False
  zonelen = len(zonetitle)

  for textline in f:

    if (not zonefound):

      # Check for the zone we are looking for, and read geometry info

      if (textline[0:zonelen] == zonetitle):
        if (designnum != 0):
          checkline = textline.split("SOLUTIONTIME=")
          checkdesign = int(checkline[1])
          if (checkdesign == designnum): zonefound = True
        else: zonefound = True
    
      if zonefound:
        splitline = textline.split(",")
        if len(splitline) > 2:
          maxt = float((splitline[1].split("="))[1])
          xmaxt = float((splitline[2].split("="))[1])
          maxc = float((splitline[3].split("="))[1])
          xmaxcline = splitline[4].split('"')[0]
          xmaxc = float((xmaxcline.split("="))[1])

    else:

      # Check to see if we've read all the coordinates

      if (textline[0:4] == "zone"): break

      # Otherwise read coordinates

      else:

        line = textline.split()
        x.append(float(line[0]))
        y.append(float(line[1]))

        # jx-mod additionally 2nd and 3rd derivative
        if (len(line) > 2):
          deriv2.append(float(line[2]))
          deriv3.append(float(line[3]))

  # Error if zone has not been found after reading the file

  if (not zonefound):
    ioerror = 2

  # Close the file

  f.close()

  # Return coordinate data

  # jx-mod additionally 2nd and 3rd derivative
  return x, y, maxt, xmaxt, maxc, xmaxc, ioerror, deriv2, deriv3

################################################################################
# Reads airfoil polars from file
def read_airfoil_polars(filename, zonetitle):

  ioerror = 0
  alpha = []
  cl = []
  cd = []
  cm = []
  xtrt = []
  xtrb = []
  #jx-mod Read also flap angle from polar file  
  flapangle = []

  # Try to open the file

  try:
    f = open(filename) 
  except IOError:
    ioerror = 1
    return ioerror

  # Read lines until we get to the correct zone

  zonefound = False
  zonelen = len(zonetitle)

  for textline in f:

    if (not zonefound):

      # Check for the zone we are looking for

      if (textline[0:zonelen] == zonetitle):
    
        zonefound = True

    else:

      # Check to see if we've read all the polars

      if (textline[0:4] == "zone"): break

      # Otherwise read polars

      else:

        line = textline.split()
        alpha.append(float(line[0]))
        cl.append(float(line[1]))
        cd.append(float(line[2]))
        cm.append(float(line[3]))
        xtrt.append(float(line[4]))
        xtrb.append(float(line[5]))

        #jx-mod Read flap-angle if exists
        if (len(line) > 6):
          flapangle.append(float(line[6]))

  # Error if zone has not been found after reading the file

  if (not zonefound):
    ioerror = 2

  # Close the file

  f.close()

  # Return polar data

  return alpha, cl, cd, cm, xtrt, xtrb, flapangle, ioerror

################################################################################
# Reads optimization history
def read_optimization_history(step):

  ioerror = 0
  fmin = 0.
  relfmin = 0.
  rad = 0.

  # Try to open the file

  try:
    f = open('optimization_history.dat') 
  except IOError:
    ioerror = 1
    return fmin, relfmin, rad, ioerror

  # Read lines until we get to the step

  stepfound = False
  for textline in f:

    if (not stepfound):

      # Check for the step we are looking for

      splitline = textline.split()
      try:
        linestep = int(splitline[0])
      except ValueError:
        continue

      if (linestep == step):
        stepfound = True
        fmin = float(splitline[1])
        relfmin = float(splitline[2])
        rad = float(splitline[3])
    
  # Error if step has not been found after reading the file

  if (not stepfound):
    ioerror = 2

  # Close the file

  f.close()

  # Return optiimzation history data

  return fmin, relfmin, rad, ioerror

################################################################################
# Loads airfoil coordinates and polars from files
def load_airfoils_from_file(coordfilename, polarfilename):

  # Initialize output data

  seedfoil = Airfoil()
  designfoils = []
  ioerror = 0

  # Check for seed airfoil coordinates

  print("Checking for airfoil coordinates file " + coordfilename + "...")

  zonetitle = 'zone t="Seed airfoil'

  # jx-mod additional 2nd and 3rd derivation
  x, y, maxt, xmaxt, maxc, xmaxc, ioerror, deriv2, deriv3 = read_airfoil_coordinates(
                                                    coordfilename, zonetitle, 0)
  if (ioerror == 1):
    print("Warning: file " + coordfilename + " not found.")
    return seedfoil, designfoils, ioerror
  elif (ioerror == 2):
    print("Error: zone labeled " + zonetitle + " not found in " + coordfilename
          + ".")
    return seedfoil, designfoils, ioerror

  seedfoil.setCoordinates(np.array(x), np.array(y))
  seedfoil.setGeometryInfo(maxt, xmaxt, maxc, xmaxc)
  # jx-mod additional 2nd and 3rd derivation
  seedfoil.setDerivatives(deriv2, deriv3)

  # Read coordinate data for designs produced by optimizer

  print("Reading airfoil coordinates from file " + coordfilename + "...")

  read_finished = False
  counter = 1
  while (not read_finished):

    zonetitle = 'zone t="Airfoil'
    x, y, maxt, xmaxt, maxc, xmaxc, ioerror, deriv2, deriv3 = read_airfoil_coordinates(
                                              coordfilename, zonetitle, counter)
    if (ioerror == 2):
      read_finished = True
      numfoils = counter - 1
      ioerror = 0
    else:
      currfoil = Airfoil()
      currfoil.setCoordinates(np.array(x), np.array(y))
      currfoil.setGeometryInfo(maxt, xmaxt, maxc, xmaxc)
      
      # jx-mod additional 2nd and 3rd derivation
      currfoil.setDerivatives(deriv2, deriv3)

      designfoils.append(currfoil)
      counter += 1

  print("Found " + str(numfoils) + " airfoil coordinates plus seed airfoil.")

  # Read seed airfoil polars (note: negative error code means coordinates were
  # read but not polars)

  print("Checking for airfoil polars file " + coordfilename + "...")

  zonetitle = 'zone t="Seed airfoil polar"'
  alpha, cl, cd, cm, xtrt, xtrb, flapangle, ioerror = read_airfoil_polars(polarfilename,
                                                                      zonetitle)
  if (ioerror == 1):
    print("Warning: file " + polarfilename + " not found.")
    return seedfoil, designfoils, 0 - ioerror
  elif (ioerror == 2):
    print("Error: zone labeled " + zonetitle + " not found in " + polarfilename
          + ".")
    return seedfoil, designfoils, 0 - ioerror

  seedfoil.setPolars(np.array(alpha), np.array(cl), np.array(cd), np.array(cm),
                     np.array(xtrt), np.array(xtrb), np.array(flapangle))

  # Read polar data for designs produced by optimizer

  print("Reading airfoil polars from file " + polarfilename + "...")

  read_finished = False
  counter = 1
  while (not read_finished):

    zonetitle = 'zone t="Polars", SOLUTIONTIME=' + str(counter)
    alpha, cl, cd, cm, xtrt, xtrb, flapangle, ioerror = read_airfoil_polars(polarfilename,
                                                                      zonetitle)
    if (ioerror == 2):
      read_finished = True
      ioerror = 0
    else:
      designfoils[counter-1].setPolars(np.array(alpha), np.array(cl), 
                                       np.array(cd), np.array(cm), 
                                       np.array(xtrt), np.array(xtrb),
                                       np.array(flapangle))
      counter += 1

  print("Found " + str(counter-1) + " airfoil polars plus seed airfoil.")
  if (counter != numfoils+1):
    print("Error: number of airfoil coordinates and polars does not match.")
    ioerror = 3 
    return seedfoil, designfoils, ioerror

  return seedfoil, designfoils, ioerror

################################################################################
# Plots airfoil coordinates
def plot_airfoil_coordinates(seedfoil, designfoils, plotnum, firsttime=True, 
                             animation=False, prefix=None):
  global plotoptions

  # Select requested airfoil

  if (plotnum > 0): foil = designfoils[plotnum-1]

  # Aliases for colors

  sc = plotoptions["color_for_seed"]
  nc = plotoptions["color_for_new_designs"]

  # Set up coordinates plot

  window_name = "Geometry  " + str(prefix)

  if (firsttime): plt.close(window_name)
  cfig = plt.figure(num= window_name)
  # jx-mod  New size and location of Windows - only hacked for res 1920x1024
  if (firsttime): plt.get_current_fig_manager().window.setGeometry(850,600,1000,450)

  ax = plt.subplot(111)
  plt.cla()

  # Auto plotting bounds

  if ( (plotoptions["show_seed_airfoil_only"]) or (plotnum == 0) ):
    xmax = np.max(seedfoil.x)
    xmin = np.min(seedfoil.x)
    ymax = np.max(seedfoil.y)
    ymin = np.min(seedfoil.y)
  elif (plotoptions["show_seed_airfoil"]):
    xmax = max([np.max(seedfoil.x), np.max(foil.x)])
    xmin = min([np.min(seedfoil.x), np.min(foil.x)])
    ymax = max([np.max(seedfoil.y), np.max(foil.y)])
    ymin = min([np.min(seedfoil.y), np.min(foil.y)])
  else:
    xmax = np.max(foil.x)
    xmin = np.min(foil.x)
    ymax = np.max(foil.y)
    ymin = np.min(foil.y)
  xrng = xmax - xmin
  # jx-mod changed from xmax= 0.1
  xmax= xmax + 0.05*xrng
  xmin= xmin - 0.05*xrng
  yrng = ymax - ymin
  ymax= ymax + 0.05*yrng
  ymin= ymin - 0.05*yrng


  # Plot airfoil coordinates

  if ( (plotoptions["show_seed_airfoil_only"]) or (plotnum == 0) ):
    ax.plot(seedfoil.x, seedfoil.y, color=sc)
  elif (plotoptions["show_seed_airfoil"]):
    ax.plot(seedfoil.x, seedfoil.y, color=sc)
    ax.plot(foil.x, foil.y, color=nc)
  else:
    ax.plot(foil.x, foil.y, color=nc)
  ax.set_aspect('equal', 'datalim')
  ax.set_xlabel('x')
  ax.set_ylabel('z')
  ax.set_xlim([xmin,xmax])


  # jx-mod Start - Plot of 2nd, 3rd derivative and delta y
  plot_2nd_deriv = False
  plot_3rd_deriv = False
  plot_delta_y   = True
  iLE = np.argmin(seedfoil.x)

  if plot_2nd_deriv:
    ax.set_ylim([-ymax,ymax])       # achieve plot of 2nd derivative is aligned in x-axis
    mirrorax0 = ax.twinx()
    mirrorax0.set_ylabel('2nd derivation', color='grey')
    mirrorax0.set_ylim(0.8, -0.8)
    mirrorax0.axhline(0, color='grey', linewidth=0.5)
    mirrorax0.plot(seedfoil.x[0:(iLE-15)],  seedfoil.deriv2[0:(iLE-15)],  color='blue', linewidth=0.8, linestyle='--') #top
    mirrorax0.plot(seedfoil.x[(-iLE+15):],  seedfoil.deriv2[(-iLE+15):],  color='grey', linewidth=0.8, linestyle='--') #bottom
    if (plotnum > 0): 
      mirrorax0.plot(foil.x[0:(iLE-15)], foil.deriv2[0:(iLE-15)], color='red', linewidth=0.8, linestyle='--')
      mirrorax0.plot(foil.x[(-iLE+15):], foil.deriv2[(-iLE+15):], color='red', linewidth=0.8, linestyle='--')

  if plot_3rd_deriv:
    ax.set_ylim([-ymax,ymax])       # achieve plot of 2nd derivative is aligned in x-axis
    mirrorax0 = ax.twinx()
    mirrorax0.set_ylabel('3rd derivation', color='grey')
    mirrorax0.set_ylim(20, -20)
    mirrorax0.axhline(0, color='grey', linewidth=0.5)
    mirrorax0.plot(seedfoil.x[0:(iLE-15)],  seedfoil.deriv3[0:(iLE-15)],  color='grey', linewidth=0.8, linestyle=':')
    mirrorax0.plot(seedfoil.x[(-iLE+15):],  seedfoil.deriv3[(-iLE+15):],  color='grey', linewidth=0.8, linestyle=':')
    if (plotnum > 0): 
      mirrorax0.plot(foil.x[0:(iLE-15)],  foil.deriv3[0:(iLE-15)],  color='magenta', linewidth=0.8, linestyle='--')
      mirrorax0.plot(foil.x[(-iLE+15):],  foil.deriv3[(-iLE+15):],  color='magenta', linewidth=0.8, linestyle='-.')

  # jx-mod Plot of delta between seed and current
  if plot_delta_y and (plotnum > 0):
    ax.set_ylim([-ymax,ymax])       # achieve plot  is aligned in x-axis
    mirrorax0 = ax.twinx()
    # mirrorax0.set_ylabel('delta z to seed ', color='green')
    mirrorax0.set_ylim(-ymax/5, ymax/5)
    mirrorax0.axes.set_yticks([])      # hide y-ticks
    mirrorax0.axhline(0, color='grey', linewidth=0.5)
    mirrorax0.plot(foil.x, (foil.y - seedfoil.y), color='green', linewidth=0.8, linestyle=':')

  # jx-mod End - Plot of 2nd, 3rd derivative and delta y

  # Display geometry info

  if plotoptions["show_airfoil_info"]:
    if ( (plotoptions["show_seed_airfoil_only"]) or (plotnum == 0) ):
      mytext = ("Thickness: " + str(seedfoil.maxt) + '\n' +
                "   at x/c: " + str(seedfoil.xmaxt) + '\n' +
                "Camber: " + str(seedfoil.maxc) + '\n' +
                "   at x/c: " + str(seedfoil.xmaxc))
      #ax.text(-0.05, -0.2, mytext, color=sc, verticalalignment='top')
      ax.text(0.0, -0.1, mytext, color=sc, verticalalignment='top', horizontalalignment='left')
    elif (plotoptions["show_seed_airfoil"]):
      mytext = ("Thickness: " + str(seedfoil.maxt) + '\n' +
                "   at x/c: " + str(seedfoil.xmaxt) + '\n' +
                "Camber: " + str(seedfoil.maxc) + '\n' +
                "   at x/c: " + str(seedfoil.xmaxc))
      # jx-mod
      # ax.text(-0.05, -0.2, mytext, color=sc, verticalalignment='top')
      ax.text(0.0, -0.1, mytext, color=sc, verticalalignment='top', horizontalalignment='left')
      mytext = ("Thickness: " + str(foil.maxt) + '\n' +
                "   at x/c: " + str(foil.xmaxt) + '\n' +
                "Camber: " + str(foil.maxc) + '\n' +
                "   at x/c: " + str(foil.xmaxc))
      # jx-mod
      # ax.text(0.70, -0.2, mytext, color=nc, verticalalignment='top')
      ax.text(1.0, -0.1, mytext, color=nc, verticalalignment='top', horizontalalignment='right')
    else:
      mytext = ("Thickness: " + str(foil.maxt) + '\n' +
                "   at x/c: " + str(foil.xmaxt) + '\n' +
                "Camber: " + str(foil.maxc) + '\n' +
                "   at x/c: " + str(foil.xmaxc))
      ax.text(0.70, -0.2, mytext, color=nc, verticalalignment='top')

  # Legend for coordinates plot

  # bbox_loc = (0.5, 1.1)

  # Fake lines for legend

  lines = []
  if ( (plotoptions["show_seed_airfoil"]) or 
       (plotoptions["show_seed_airfoil_only"]) or (plotnum == 0) ):
    fakeline = plt.Line2D((0,1),(0,0), color=sc, label="Seed airfoil")
    lines.append(fakeline)
  if ( (not plotoptions["show_seed_airfoil_only"]) and (plotnum != 0) ):
    fakeline = plt.Line2D((0,1),(0,0), color=nc, 
                          label="Design number " + str(plotnum))
    lines.append(fakeline)

  # Create legend

  labels = [l.get_label() for l in lines]
  ax.legend(lines, labels, loc="upper right", numpoints=1)
    
  # Update plot for animation only (for others, plt.show() must be called
  # separately)

  if animation:
    if (firsttime): cfig.show()
    else: plt.pause(0.0001)
    cfig.canvas.draw()

    # Save animation frames if requested

    if plotoptions["save_animation_frames"]:
      if (prefix == None):
        print("Error: no file prefix specified - cannot save animation frames.")
      else:
        imagefname = prefix + '_coordinates.png'
        print("Saving image frame to file " + imagefname + ' ...')
        plt.savefig(imagefname)

################################################################################
# Plots polars
def plot_polars(seedfoil, designfoils, plotnum, firsttime=True, animation=False,
                prefix=None, pfig=None, axarr=None, legend=None):
  global plotoptions

  # Select requested airfoil

  if (plotnum > 0): foil = designfoils[plotnum-1]

  # Aliases for colors

  sc = plotoptions["color_for_seed"]
  nc = plotoptions["color_for_new_designs"]

  window_name = "Polars  " + str(prefix)

  # Set up polars plot. Note: for monitoring, must pass pfig, axarr, and legend
  # after the initial plotting, because currently plt.subplots always creates a
  # new figure and the only other options would be to save these as global
  # variables or to destroy and recreate the figure each time. (Or: make the
  # polar plot its own class?)

  if firsttime: 
    plt.close(window_name)
    # jx-mod new subplots for glide ratio and climb
    pfig, axarr = plt.subplots(2, 3, num= window_name)
    pfig.subplots_adjust(hspace=0.3, wspace=0.3)
    pfig.set_size_inches(11, 8, forward=True)
    # jx-mod new sizing
    plt.get_current_fig_manager().window.setGeometry(100,30,1300,550)

  else:
    plt.figure(num= window_name)
    axarr[0,0].clear()
    axarr[0,1].clear()
    axarr[1,0].clear()
    axarr[1,1].clear()
    #jx-mod glide ration & climb
    axarr[0,2].clear()
    axarr[1,2].clear()
  plt.cla()

  # Auto plotting bounds

  if ( (plotoptions["show_seed_polar_only"]) or (plotnum == 0) ):
    almax = np.max(seedfoil.alpha)
    almin = np.min(seedfoil.alpha)
    clmax = np.max(seedfoil.cl)
    clmin = np.min(seedfoil.cl)
    cdmax = np.max(seedfoil.cd)
    cdmin = np.min(seedfoil.cd)
    cmmax = np.max(seedfoil.cm)
    cmmin = np.min(seedfoil.cm)
    xtrmax = max([np.max(seedfoil.xtrt), np.max(seedfoil.xtrb)])
    xtrmin = min([np.min(seedfoil.xtrt), np.min(seedfoil.xtrb)])
    # jx-mod
    glidemax = np.max(seedfoil.glide)
    glidemin = np.min(seedfoil.glide)
    climbmax = np.max(seedfoil.climb)
    climbmin = np.min(seedfoil.climb)
  elif (plotoptions["show_seed_polar"]):
    almax = max([np.max(seedfoil.alpha), np.max(foil.alpha)])
    almin = min([np.min(seedfoil.alpha), np.min(foil.alpha)])
    clmax = max([np.max(seedfoil.cl), np.max(foil.cl)])
    clmin = min([np.min(seedfoil.cl), np.min(foil.cl)])
    cdmax = max([np.max(seedfoil.cd), np.max(foil.cd)])
    cdmin = min([np.min(seedfoil.cd), np.min(foil.cd)])
    cmmax = max([np.max(seedfoil.cm), np.max(foil.cm)])
    cmmin = min([np.min(seedfoil.cm), np.min(foil.cm)])
    xtrmax = max([np.max(seedfoil.xtrt), np.max(seedfoil.xtrb),
                  np.max(foil.xtrt), np.max(foil.xtrb)])
    xtrmin = min([np.min(seedfoil.xtrt), np.min(seedfoil.xtrb),
                  np.min(foil.xtrt), np.min(foil.xtrb)])
    # jx-mod
    glidemax = max([np.max(seedfoil.glide), np.max(foil.glide)])
    glidemin = min([np.min(seedfoil.glide), np.min(foil.glide)])
    climbmax = max([np.max(seedfoil.climb), np.max(foil.climb)])
    climbmin = min([np.min(seedfoil.climb), np.min(foil.climb)])
  else:
    almax = np.max(foil.alpha)
    almin = np.min(foil.alpha)
    clmax = np.max(foil.cl)
    clmin = np.min(foil.cl)
    cdmax = np.max(foil.cd)
    cdmin = np.min(foil.cd)
    cmmax = np.max(foil.cm)
    cmmin = np.min(foil.cm)
    xtrmax = max([np.max(foil.xtrt), np.max(foil.xtrb)])
    xtrmin = min([np.min(foil.xtrt), np.min(foil.xtrb)])
    # jx-mod
    glidemax = np.max(foil.glide)
    glidemin = np.min(foil.glide)
    climbmax = np.max(foil.climb)
    climbmin = np.min(foil.climb)
  alrng = almax - almin
  almax = almax + 0.1*alrng
  almin = almin - 0.1*alrng
  cdrng = cdmax - cdmin
  cdmax = cdmax + 0.1*cdrng
  cdmin = cdmin - 0.1*cdrng
  clrng = clmax - clmin
  clmax = clmax + 0.1*clrng
  clmin = clmin - 0.1*clrng
  cmrng = cmmax - cmmin
  cmmax = cmmax + 0.1*cmrng
  cmmin = cmmin - 0.1*cmrng
  xtrrng = xtrmax - xtrmin
  xtrmax = xtrmax + 0.1*xtrrng
  xtrmin = xtrmin - 0.1*xtrrng
  # jx-mod
  gliderng = glidemax - glidemin
  glidemax = glidemax + 0.1*gliderng
  glidemin = glidemin - 0.1*gliderng
  if (glidemin < 0.0):
    glidemin = 0.0
  climbrng = climbmax - climbmin
  climbmax = climbmax + 0.1*climbrng
  climbmin = climbmin - 0.1*climbrng
  if (climbmin < 0.0):
    climbmin = 0.0
  glideclmin = clmin
  if (glideclmin < 0.0):
    glideclmin = 0.0 


  # Plot polars

  if ( (plotoptions["show_seed_polar_only"]) or (plotnum == 0) ):
    axarr[0,0].plot(seedfoil.alpha, seedfoil.cl, linestyle='-', color=sc,
                    marker='o')
    if plotoptions["drag_plot_type"] == "vs. lift":
      axarr[0,1].plot(seedfoil.cd, seedfoil.cl, linestyle='-', color=sc, 
                      marker='o')
    else:
      axarr[0,1].plot(seedfoil.alpha, seedfoil.cd, linestyle='-', color=sc,
                      marker='o')
    axarr[1,0].plot(seedfoil.alpha, seedfoil.cm, linestyle='-', color=sc,
                    marker='o')
    axarr[1,1].plot(seedfoil.xtrt, seedfoil.alpha, linestyle='-', color=sc, 
                    marker='o')
    axarr[1,1].plot(seedfoil.xtrb, seedfoil.alpha, linestyle='--', color=sc, 
                    marker='o')
    # jx-mod show glide and climb plot
    if plotoptions["drag_plot_type"] == "vs. lift":
      axarr[0,2].plot(seedfoil.cl, seedfoil.glide, linestyle='-', color=sc, 
                      marker='o')
      axarr[1,2].plot(seedfoil.cl, seedfoil.climb, linestyle='-', color=sc, 
                      marker='o')
    else:
      axarr[0,2].plot(seedfoil.alpha, seedfoil.glide, linestyle='-', color=sc, 
                      marker='o')
      axarr[1,2].plot(seedfoil.alpha, seedfoil.climb, linestyle='-', color=sc, 
                      marker='o')
  elif (plotoptions["show_seed_polar"]):
    axarr[0,0].plot(seedfoil.alpha, seedfoil.cl, linestyle='-', color=sc, 
                    marker='o')
    axarr[0,0].plot(foil.alpha, foil.cl, linestyle='-', color=nc, marker='s')
    if plotoptions["drag_plot_type"] == "vs. lift":
      axarr[0,1].plot(seedfoil.cd, seedfoil.cl, linestyle='-', color=sc, 
                      marker='o')
      axarr[0,1].plot(foil.cd, foil.cl, linestyle='-', color=nc, marker='s')

      # jx-mod show cd-value in graph
      for i in range(len(foil.cl)): 
        if ((len(foil.flapangle) > 0) and (foil.flapangle[i] != 0)):
          axarr[0,1].annotate(('f {:5.2f}'.format(foil.flapangle[i])), (cdmax - 0.29*cdrng, foil.cl[i]), fontsize = 9)
        else:
          axarr[0,1].annotate(('   cd {:5.4f}'.format(foil.cd[i])), (foil.cd[i], foil.cl[i]), fontsize = 9)
    else:
      axarr[0,1].plot(seedfoil.alpha, seedfoil.cd, linestyle='-', color=sc, 
                      marker='o')
      axarr[0,1].plot(foil.alpha, foil.cd, linestyle='-', color=nc, marker='s')
    axarr[1,0].plot(seedfoil.alpha, seedfoil.cm, linestyle='-', color=sc, 
                    marker='o')
    axarr[1,0].plot(foil.alpha, foil.cm, linestyle='-', color=nc, marker='s')
    axarr[1,1].plot(seedfoil.xtrt, seedfoil.alpha, linestyle='-', color=sc, 
                    marker='o')
    axarr[1,1].plot(foil.xtrt, foil.alpha, linestyle='-', color=nc, marker='s')
    axarr[1,1].plot(seedfoil.xtrb, seedfoil.alpha, linestyle='--', color=sc, 
                    marker='o')
    axarr[1,1].plot(foil.xtrb, foil.alpha, linestyle='--', color=nc, marker='s')
    # jx-mod
    if plotoptions["drag_plot_type"] == "vs. lift":
      axarr[0,2].plot(seedfoil.cl, seedfoil.glide, linestyle='-', color=sc, marker='o')
      axarr[0,2].plot(foil.cl, foil.glide, linestyle='-', color=nc, marker='s')
      axarr[1,2].plot(seedfoil.cl, seedfoil.climb, linestyle='-', color=sc, marker='o')
      axarr[1,2].plot(foil.cl, foil.climb, linestyle='-', color=nc, marker='s')
    else:
      axarr[0,2].plot(seedfoil.alpha, seedfoil.glide, linestyle='-', color=sc, marker='o')
      axarr[0,2].plot(foil.alpha, foil.glide, linestyle='-', color=nc, marker='s')
      axarr[1,2].plot(seedfoil.alpha, seedfoil.climb, linestyle='-', color=sc, marker='o')
      axarr[1,2].plot(foil.alpha, foil.climb, linestyle='-', color=nc, marker='s')
  else: 
    axarr[0,0].plot(foil.alpha, foil.cl, linestyle='-', color=nc, marker='s') 
    if plotoptions["drag_plot_type"] == "vs. lift":
      axarr[0,1].plot(foil.cd, foil.cl, linestyle='-', color=nc, marker='s') 
      # jx-mod
      axarr[0,2].plot(foil.cl, foil.glide, linestyle='-', color=nc, marker='s')
      axarr[1,2].plot(foil.cl, foil.climb, linestyle='-', color=nc, marker='s')
    else:
      axarr[0,1].plot(foil.alpha, foil.cd, linestyle='-', color=nc, marker='s') 
      # jx-mod
      axarr[0,2].plot(foil.alpha, foil.glide, linestyle='-', color=nc, marker='s')
      axarr[1,2].plot(foil.alpha, foil.climb, linestyle='-', color=nc, marker='s')
    axarr[1,0].plot(foil.alpha, foil.cm, linestyle='-', color=nc, marker='s') 
    axarr[1,1].plot(foil.xtrt, foil.alpha, linestyle='-', color=nc, marker='s')
    axarr[1,1].plot(foil.xtrb, foil.alpha, linestyle='--', color=nc, marker='s')
  axarr[0,0].set_xlabel('Angle of attack')
  axarr[0,0].set_ylabel('Lift coefficient')
  axarr[0,0].set_xlim([almin,almax])
  axarr[0,0].set_ylim([clmin,clmax])
  axarr[0,0].grid()
  if plotoptions["drag_plot_type"] == "vs. lift":
    axarr[0,1].set_xlabel('Drag coefficient')
    axarr[0,1].set_ylabel('Lift coefficient')
    axarr[0,1].set_xlim([cdmin,cdmax])
    axarr[0,1].set_ylim([clmin,clmax])
  else:
    axarr[0,1].set_xlabel('Angle of attack')
    axarr[0,1].set_ylabel('Drag coefficient')
    axarr[0,1].set_xlim([almin,almax])
    axarr[0,1].set_ylim([cdmin,cdmax])
  axarr[0,1].grid()
  axarr[1,0].set_xlabel('Angle of attack')
  axarr[1,0].set_ylabel('Pitching moment coefficient')
  axarr[1,0].set_xlim([almin,almax])
  axarr[1,0].set_ylim([cmmin,cmmax])
  axarr[1,0].grid()
  axarr[1,1].set_xlabel('Transition x/c\n(top: solid, bottom: dashed)')
  axarr[1,1].set_ylabel('Angle of attack')
  axarr[1,1].set_xlim([xtrmin,xtrmax])
  axarr[1,1].set_ylim([almin,almax])
  axarr[1,1].grid()
  #jx-mod
  if plotoptions["drag_plot_type"] == "vs. lift":
    axarr[0,2].set_xlabel('Lift coefficient')
    axarr[0,2].set_ylabel('Glide ratio')
    axarr[0,2].set_xlim([glideclmin,clmax])
    axarr[0,2].set_ylim([glidemin,glidemax])
  else:
    axarr[0,2].set_xlabel('Angle of attack')
    axarr[0,2].set_ylabel('Glide ratio')
    axarr[0,2].set_xlim([almin,almax])
    axarr[0,2].set_ylim([glidemin,glidemax])
  axarr[0,2].grid(True)
  if plotoptions["drag_plot_type"] == "vs. lift":
    axarr[1,2].set_xlabel('Lift coefficient')
    axarr[1,2].set_ylabel('Climb ratio')
    axarr[1,2].set_xlim([glideclmin,clmax])
    axarr[1,2].set_ylim([climbmin,climbmax])
  else:
    axarr[1,2].set_xlabel('Angle of attack')
    axarr[1,2].set_ylabel('Climb ratio')
    axarr[1,2].set_xlim([almin,almax])
    axarr[1,2].set_ylim([climbmin,climbmax])
  axarr[1,2].grid(True)

  # Draw legend

  lines = []
  if ( (plotoptions["show_seed_polar"]) or 
       (plotoptions["show_seed_polar_only"]) or (plotnum == 0) ):
    fakeline = plt.Line2D((0,1),(0,0), linestyle='-', color=sc, marker='o',
                          label="Seed airfoil")
    lines.append(fakeline)
  if ( (not plotoptions["show_seed_polar_only"]) and (plotnum != 0) ):
    fakeline = plt.Line2D((0,1),(0,0), linestyle='-', color=nc, marker='s', 
                          label="Design number " + str(plotnum))
    lines.append(fakeline)
  
  bbox_loc = (0.5, 1.00)
  labels = [l.get_label() for l in lines]
  if legend: legend.remove()
  legend = pfig.legend(lines, labels, loc="upper center", 
                       bbox_to_anchor=bbox_loc, numpoints=1)

  # Update plot for animation only (for others, plt.show() must be called
  # separately)

  if animation:
    if (firsttime): pfig.show()
    else: plt.pause(0.0001)
    pfig.canvas.draw()

    # Save animation frames if requested
  
    if plotoptions["save_animation_frames"]:
      if (prefix == None):
        print("Error: no file prefix specified - cannot save animation frames.")
      else:
        imagefname = prefix + '_polars.png'
        print("Saving image frame to file " + imagefname + ' ...')
        plt.savefig(imagefname)

  return pfig, axarr, legend

################################################################################
# Plots optimization history
def plot_optimization_history(steps, fmins, relfmins, rads, firsttime=True,
                              animation=False, ofig=None, axarr=None, 
                              prefix=None, mirrorax0=None):
  global plotoptions

  # Set up optimization history plot. Note: for monitoring, must pass ofig, 
  # axarr, and mirrorax0 after the initial plotting, because currently
  # plt.subplots always creates a new figure, and the only other option would be
  # to save these as global variables or to destry and recreate the figure each
  # time. (Or: make the optimization history plot its own class?)

  window_name = "Optimization History  " + str(prefix)

  if firsttime: 
    plt.close(window_name)
    ofig, axarr = plt.subplots(2, 1, num= window_name)
  else:
    plt.figure(num= window_name)
    axarr[0].clear()
    mirrorax0.clear()
    axarr[1].clear()
  plt.cla()

  # jx-mod
  if firsttime: plt.get_current_fig_manager().window.setGeometry(1420,70,480,380)

  # Plot optimization history

  axarr[0].plot(steps, fmins, color='blue')
  for t1 in axarr[0].get_yticklabels(): t1.set_color('blue')
  if (firsttime): mirrorax0 = axarr[0].twinx()
  mirrorax0.plot(steps, relfmins, color='red')
  for t2 in mirrorax0.get_yticklabels(): t2.set_color('red')
  axarr[1].plot(steps, rads)

  axarr[0].set_xlabel('Iteration')
  axarr[0].set_ylabel('Objective function', color='blue')
  mirrorax0.set_ylabel('% Improvement over seed', color='red')
  axarr[1].set_xlabel('Iteration')
  axarr[1].set_ylabel('Design radius')
  axarr[1].set_yscale("log")
  axarr[1].grid()

  # Update plot for animation only (for others, plt.show() must be called
  # separately)

  if animation:
    if (firsttime): ofig.show()
    else: plt.pause(0.0001)
    ofig.canvas.draw()

  return ofig, axarr, mirrorax0

################################################################################
# Input function that checks python version
def my_input(message):

  # Check python version

  python_version = version_info[0]

  # Issue correct input command

  if (python_version == 2):
    return raw_input(message)
  else:
    return input(message)

################################################################################
# Plotting menu
def plotting_menu(seedfoil, designfoils):
  global plotoptions

  # Load optimization history data if it's available

  steps, fmins, relfmins, rads, ioerror = read_new_optimization_history()

  numfoils = len(designfoils)
  plotting_complete = False
  validchoice = False
  while (not validchoice):
    print("")
    print("There are " + str(numfoils) + " designs.")
    plotnum = int(my_input("Enter design to plot (or 0 to return): "))

    # Return to main menu

    if (plotnum == 0): 
      validchoice = True
      plotting_complete = True

    # Check for bad index

    elif ( (plotnum < 1) or (plotnum > numfoils) ):
      validchoice = False
      print("Error: index out of bounds.")

    # Plot design

    else:
      validchoice = True
      plt.close()
      if plotoptions["plot_airfoils"]:
        plot_airfoil_coordinates(seedfoil, designfoils, plotnum, firsttime=True)
      if plotoptions["plot_polars"]:
        plot_polars(seedfoil, designfoils, plotnum, firsttime=True)
      if (plotoptions["plot_optimization_history"] and steps.shape[0] > 0):
        plot_optimization_history(steps, fmins, relfmins, rads, firsttime=True)
      plt.show()
      plotting_complete = True 

  return plotting_complete

################################################################################
# Reads new airfoil coordinates and polar files for updates during optimization
def read_new_airfoil_data(seedfoil, designfoils, prefix):

  # Temporary airfoil struct

  foil = Airfoil()

  # Set up file names to monitor

  coordfilename = prefix + '_design_coordinates.dat'
  polarfilename = prefix + '_design_polars.dat'

  # Loop through files until we reach latest available design

  reading = True
  while reading:

    if (seedfoil.npt == 0):
      zonetitle = 'zone t="Seed airfoil"'
      foilstr = 'seed'
      nextdesign = 0
    else:
      nextdesign = len(designfoils) + 1
      zonetitle = 'zone t="Airfoil'
      foilstr = 'design number ' + str(nextdesign)
  
    # Read data from coordinate file
  
    x, y, maxt, xmaxt, maxc, xmaxc, ioerror, deriv2, deriv3 = read_airfoil_coordinates(
                                           coordfilename, zonetitle, nextdesign)
    if (ioerror == 1):
      print("Airfoil coordinates file " + coordfilename + " not available yet.")
      reading = False
      break
    elif (ioerror == 2):
      reading = False
      break
    else:
      print("Read coordinates for " + foilstr + ".")
      foil.setCoordinates(np.array(x), np.array(y))
      # jx-mod additional 2nd and 3rd deriv
      foil.setDerivatives (deriv2, deriv3)
      # jx-mod additional 2nd and 3rd deriv
      foil.setGeometryInfo(maxt, xmaxt, maxc, xmaxc)
  
    # Set zone title for polars
  
    if (foilstr == 'seed'):
      zonetitle = 'zone t="Seed airfoil polar"'
    else:
      zonetitle = 'zone t="Polars", SOLUTIONTIME=' + str(nextdesign)
  
    # Read data from polar file (not: negative error code means coordinates were
    # read but not polars)
  
    alpha, cl, cd, cm, xtrt, xtrb, flapangle, ioerror = read_airfoil_polars(polarfilename,
                                                                      zonetitle)
    if ( (ioerror == 1) or (ioerror == 2) ):
      print("Warning: polars will not be available for this design.")
      ioerror = 3
      reading = False
    else:
      print("Read polars for " + foilstr + ".")
      foil.setPolars(np.array(alpha), np.array(cl), np.array(cd), np.array(cm),
                     np.array(xtrt), np.array(xtrb), np.array(flapangle))
  
    # Copy data to output objects
  
    if (foilstr == 'seed'): seedfoil = foil
    else: designfoils.append(foil)

  return seedfoil, designfoils, ioerror

################################################################################
# Reads new optimization history data for updates during optimization
def read_new_optimization_history(steps=None, fmins=None, relfmins=None, 
                                  rads=None):

  if ((steps is None) or steps.shape[0] == 0):
    steps = np.zeros((0), dtype=int)
    fmins = np.zeros((0))
    relfmins = np.zeros((0))
    rads = np.zeros((0))
    currstep = 0
  else: 
    numsteps = steps.shape[0]
    currstep = steps[numsteps-1]

  # Loop through file until we reach latest available step

  reading = True
  while reading:

    if (steps.shape[0] == 0): nextstep = 1
    else:
      numsteps = steps.shape[0]
      nextstep = steps[numsteps-1] + 1
  
    # Read data from optimization history file
  
    fmin, relfmin, rad, ioerror = read_optimization_history(nextstep)
    if (ioerror == 1):
      print("optimization_history.dat not available yet.")
      reading = False
    elif (ioerror == 2):
      reading = False
      if (nextstep - 1 > currstep):
        print("Read optimization data to step " + str(nextstep-1) + ".")
    else: 
      # Copy data to output objects
  
      steps = np.append(steps, nextstep)
      fmins = np.append(fmins, fmin)
      relfmins = np.append(relfmins, relfmin)
      rads = np.append(rads, rad)

  return steps, fmins, relfmins, rads, ioerror

################################################################################
# Gets boolean input from user
def get_boolean_input(key, keyval):

  validchoice = False
  while (not validchoice):
    print("Current value for " + key + ": " + str(keyval))
    print("Available choices: True, False\n")
    sel = my_input("Enter new value: ")
    if ( (sel == "True") or (sel == "true")):
      retval = True
      validchoice = True
    elif ( (sel == "False") or (sel == "false")):
      retval = False
      validchoice = True
    else:
      print("Please enter True or False.")
      validchoice = False

  return retval

################################################################################
# Gets color input from user
def get_color_input(key, keyval):

  colors = ["blue", "green", "red", "cyan", "magenta", "yellow", "black"]

  validchoice = False
  while (not validchoice):
    print("Current value for " + key + ": " + str(keyval))
    print("Available choices: blue, green, red, cyan, magenta, yellow, black\n")
    sel = my_input("Enter new value: ")

    # Check for valid color

    for c in colors:
      if (sel == c):
        validchoice = True
        retval = sel
        break
    if (not validchoice):
      print("Invalid color specified.  Please enter a valid color.")
      validchoice = False

  return retval

################################################################################
# Gets float input from user, subject to user-supplied min and max values
def get_float_input(key, keyval, minallow=None, maxallow=None):

  validchoice = False
  while (not validchoice):
    print("Current value for " + key + ": " + str(keyval) + '\n')
    sel = my_input("Enter new value: ")

    # Check for bad format

    try:
      val = float(sel)
    except ValueError:
      print("Error: " + key + " must be a floating point number.")
      validchoice = False
      continue

    # Check for out-of-bounds selection

    if (minallow != None):
      if (val <= minallow):
        print("Error: " + key + " must be greater than " + str(minallow) + ".")
        validchoice = False
        continue
    if (maxallow != None):
      if (val >= maxallow):
        print("Error: " + key + " must be less than " + str(maxallow) + ".")
        validchoice = False
        continue

    # If it passed all these checks, it's an acceptable input

    validchoice = True
    retval = val

  return retval

################################################################################
# Gets drag plot type from user input
def get_drag_plot_type(key, keyval):

  validchoice = False
  while (not validchoice):
    print("Current value for " + key + ": " + str(keyval))
    print("Available choices: vs. lift, vs. alpha\n")
    sel = my_input("Enter new value: ")
    if ( (sel == "vs. lift") or (sel == "vs lift") ):
      retval = "vs. lift"
      validchoice = True
    elif ( (sel == "vs. alpha") or (sel == "vs alpha") ):
      retval = "vs. alpha"
      validchoice = True
    else:
      print("Please enter vs. lift or vs. alpha.")
      validchoice = False

  return retval

################################################################################
# Options menu: allows user to change plot options
def options_menu():
  global plotoptions

  # Status variable

  options_complete = False
 
  # Print list of plotting options

  print("")
  print("Available plotting options:")
  print("")
  for key in sorted(plotoptions):
    print(key + " [" + str(plotoptions[key]) + "]")
  print("")

  # Get user input

  key = my_input("Enter option to change (or 0 to return): ")

  # True/False settings

  if ( (key == "show_seed_airfoil") or (key == "show_seed_airfoil_only") or
       (key == "show_seed_polar") or (key == "show_seed_polar_only") or
       (key == "save_animation_frames") or (key == "plot_airfoils") or 
       (key == "plot_polars") or (key == "show_airfoil_info") or
       (key == "plot_optimization_history") ):
    options_complete = False
    plotoptions[key] = get_boolean_input(key, plotoptions[key])

  # Change colors

  elif ( (key == "color_for_seed") or (key == "color_for_new_designs") ):
    options_complete = False
    plotoptions[key] = get_color_input(key, plotoptions[key])

  # Change drag plot type

  elif key == "drag_plot_type":
    options_complete = False
    plotoptions[key] = get_drag_plot_type(key, plotoptions[key])

  # Change monitor update interval

  elif (key == "monitor_update_interval"):
    options_complete = False
    plotoptions[key] = get_float_input(key, plotoptions[key], minallow=0.0)

  # Exit options menu

  elif (key == "0"): 
    options_complete = True

  # Error for invalid input

  else: 
    options_complete = False
    print("Unrecognized plot option.")

  return options_complete

################################################################################
# Main menu
################################################################################
#
# jx-mod initialchoice to autostart operation
#
def main_menu(initialchoice, seedfoil, designfoils, prefix):
  global plotoptions

  exitchoice = False
  while (not exitchoice):
    
    if initialchoice:
      choice = initialchoice
      initialchoice = ""
    else:
      print("")
      print("Options:")
      print("[0] Exit")
      print("[1] Plot a specific design")
      print("[2] Animate all designs")
      print("[3] Monitor an ongoing optimization")
      print("[4] Change plotting options")
      print("")

      choice = my_input("Enter a choice [0-4]: ")

    # Exit design_visualizer

    if (choice == "0"):
      exitchoice = True

    # Plot a single design

    elif (choice == "1"):
      exitchoice = False

      # Turn on matplotlib toolbar

      # rcParams['toolbar'] = 'toolbar2'
      rcParams['toolbar'] = 'None'


      # Go to plotting menu

      plotting_complete = False
      while (not plotting_complete): plotting_complete = plotting_menu(
                                                          seedfoil, designfoils)

    # Animate all designs

    elif (choice == "2"):
      exitchoice = False

      # Number of digits in design counter string

      numfoils = len(designfoils)
      if (numfoils == 0):
        print("There are no designs to animate.  Run xoptfoil first.")
        continue
      width = int(floor(log10(float(numfoils)))) - 1

      # Turn off matplotlib toolbar

      rcParams['toolbar'] = 'None'

      # Loop through designs, updating plot

      pfig = None
      axarr = None
      leg = None
      plt.close()
      for i in range(0, numfoils):
        if (i == 0): init = True
        else: init = False

        if (plotoptions["save_animation_frames"]):

          # Determine number of zeroes to pad with and image file prefix
  
          currwidth = int(floor(log10(float(i+1)))) - 1
          numzeroes = width - currwidth
          imagepref = prefix + numzeroes*'0' + str(i+1)

        else: imagepref = None

        # Update plots

        if plotoptions["plot_airfoils"]:
          plot_airfoil_coordinates(seedfoil, designfoils, i+1, firsttime=init,
                                   animation=True, prefix=imagepref)
        if plotoptions["plot_polars"]:
          pfig, axarr, leg = plot_polars(seedfoil, designfoils, i+1, 
                                         firsttime=init, animation=True, 
                                         prefix=imagepref, pfig=pfig, 
                                         axarr=axarr, legend=leg)

    # Monitor optimization progress

    elif (choice == "3"):
      exitchoice = False
      print('Monitoring optimization progress. To stop, enter the command ' +
            '"stop_monitoring" in run_control.')

      # Turn off matplotlib toolbar and temporarily disable saving images

      rcParams['toolbar'] = 'None'
      temp_save_frames = plotoptions["save_animation_frames"]
      plotoptions["save_animation_frames"] = False

      # Read airfoil coordinates, polars, and optimization history
      # (clears any data from previous run)

      seedfoil, designfoils, ioerror = load_airfoils_from_file(
                                                   coordfilename, polarfilename)
      steps, fmins, relfmins, rads, ioerror = read_new_optimization_history()

      # Periodically read data and update plot

      init = True
      monitoring = True
      pfig = None
      axarr = None
      leg = None
      ofig = None
      oaxarr = None
      mirrorax0 = None
      plt.close()
      while (monitoring):

        # Update plot

        if (ioerror != 1): 
          numfoils = len(designfoils)
          if plotoptions["plot_airfoils"]:
            plot_airfoil_coordinates(seedfoil, designfoils, numfoils, 
                                     firsttime=init, animation=True, prefix = prefix)
          if plotoptions["plot_polars"]:
            pfig, axarr, leg = plot_polars(seedfoil, designfoils, numfoils,
                                           firsttime=init, animation=True, prefix = prefix,
                                           pfig=pfig, axarr=axarr, legend=leg)
          if plotoptions["plot_optimization_history"]:
            ofig, oaxarr, mirrorax0 = plot_optimization_history(steps, fmins, 
                                      relfmins, rads, firsttime=init, prefix = prefix,
                                      animation=True, ofig=ofig, axarr=oaxarr,
                                      mirrorax0=mirrorax0)
          init = False

        # Pause for requested update interval

        plt.pause(plotoptions["monitor_update_interval"])

        # Update airfoil and optimization data

        seedfoil, designfoils, ioerror = read_new_airfoil_data(seedfoil,
                                                            designfoils, prefix)
        steps, fmins, relfmins, rads, ioerror = read_new_optimization_history(
                                                   steps, fmins, relfmins, rads)
        
        # Check for stop_monitoring in run_control file

        try:
          f = open('run_control')
        except IOError:
          continue

        commands = []
        for line in f:
          commands += [line.strip()]
          if (line.strip() == "stop_monitoring"):
            print("stop_monitoring command found. Returning to main menu.")
            monitoring = False
# jx-mod  stop also if "stop" for xoptfoil found ...
          if (line.strip() == "stop"):
            print("stop command found. Returning to main menu.")
            monitoring = False
        f.close()
        if len(commands) > 0:
          f = open('run_control', 'w')
          for command in commands:
            if command != "stop_monitoring":
              f.write(command + '\n') 
          f.close()
 
      # Change save_animation_frames back to original setting when done

      plotoptions["save_animation_frames"] = temp_save_frames

    # Change plotting options

    elif (choice == "4"):
      exitchoice = False
      options_complete = False
      while (not options_complete): options_complete = options_menu()

    # Invalid choice

    else:
      print("Error: please enter a choice 0-4.")

################################################################################
# Main design_visualizer program
if __name__ == "__main__":

  # Read command line arguments

# jx-mod Implementation of command line arguments - start
# 
#  Following registry keys must be set to handle command line arguments to python.exe
#        HKEY_CLASSES_ROOT\Applications\python.exe\shell\open\command --> "xx:\Anaconda3\python.exe" "%1" %*
#        HKEY_CLASSES_ROOT\py_auto_file\shell\open\command            --> "xx:\Anaconda3\python.exe" "%1" %*
#
  # initiate the parser
  parser = argparse.ArgumentParser('')
  parser.add_argument("--option", "-o", help="set initial action option", choices=['1','2','3','4'])
  parser.add_argument("--case", "-c", help="the case name for the optimization (e.g., optfoil)")

  # read arguments from the command line
  args = parser.parse_args()

  if args.case:
    prefix = args.case
  else:
    print("Enter the case name for the optimization (e.g., optfoil, which ")
    prefix = my_input("is the default case name): ")
    print("")

  coordfilename = prefix + '_design_coordinates.dat'
  polarfilename = prefix + '_design_polars.dat'

  # Read airfoil coordinates and polars

  seedfoil, designfoils, ioerror = load_airfoils_from_file(
                                                   coordfilename, polarfilename)

  # Warning if file is not found

  if (ioerror == 1):
    print("You will not be able to create plots until coordinate data is read.")
  elif (ioerror < 0):
    print("Only airfoils are available for plotting (no polars).")

  # Call main menu

  # jx-mod additional inital option choice (for autostart)
  if (abs(ioerror) <= 1): main_menu(args.option, seedfoil, designfoils, prefix)
