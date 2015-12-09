#!/usr/bin/python
"""
Author: Wong Matthew Lun
Date: 2015-12-08 6:21PM
"""

from PolyDataHandler import CenterLineHandler, ArmSurfaceHandler
import vtk
import optparse
import sys

def main(args):
    parser = optparse.OptionParser()
    parser.add_option("-s", "--surface",action="store", dest="surface", default=True,help="Input surface filename.")
    parser.add_option("-q", "--quiet",action="store_true", dest="quiet", default=False,help="Suppress console outputs")
    parser.add_option("-c", "--centerline",action="store", dest="centerline", default=True,help="Input centerline filename.")
    parser.add_option("-o", "--output", action="store", dest="outFileName", type=str, default="drilled.stl", help="Set output file name")
    parser.add_option("-m", "--holesPerSlice", action="store", dest="holesPerSlice", type=int, default=5, help="Set number of holes per slice")
    parser.add_option("-n", "--numOfSlice", action="store", dest="numOfSlice", type=int ,default=5, help="Set number of slices")
    parser.add_option("-r", "--radius", action="store", dest="radius", type=float, default=5, help="Set hole radius required")
    parser.add_option("-p", "--padding", action="store", dest="padding", type=str, default="1,0", help="Set padding level where no holes are drilled")
    parser.add_option("-e", "--errorTorlerance", action="store", dest="error", type=float, default=1, help="Set maximum error tolerance from idea grid in degrees")

    (options, args) = parser.parse_args()
    surfaceFileName = options.surface
    centerlineFileName = options.centerline
    outFileName = options.outFileName
    [startPadding, endPadding] = [int(options.padding.split(",")[i]) for i in xrange(2)]

    # create center line object
    cl = CenterLineHandler(centerlineFileName)
    cl.Read()

    # careate arm object
    arm = ArmSurfaceHandler(surfaceFileName, cl)
    arm.Read()

    # Get a list of holes and then drill
    holelist = arm.GetSemiUniDistnaceGrid(options.holesPerSlice, options.numOfSlice, options.error, startPadding, endPadding)
    arm.Drill(holelist, options.radius, options.quiet)

    clippermapper = vtk.vtkPolyDataMapper()
    clippermapper.SetInputData(arm._data)

    writer = vtk.vtkSTLWriter()
    writer.SetFileName(outFileName)
    writer.SetInputData(arm._data)
    writer.Write()

if __name__ == '__main__':
    main(sys.argv)



