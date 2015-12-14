#!/usr/bin/python
"""
Author: Wong Matthew Lun
Date: 2015-12-08 6:21PM

Return error list:
0   Success
1   Most likely Write Failed
2   Cannot find surface/centerline files
3   Error tolerance errors

"""

from PolyDataHandler import CenterLineHandler, ArmSurfaceHandler
import vtk
import optparse
import sys
import os

def main(args):
    parser = optparse.OptionParser()
    parser.add_option("-s", "--surface",action="store", dest="surface", default=True,help="Input surface filename.")
    parser.add_option("-q", "--quiet",action="store_true", dest="quiet", default=False,help="Suppress console outputs")
    parser.add_option("-c", "--centerline",action="store", dest="centerline", default=True,help="Input centerline filename.")
    parser.add_option("-o", "--output", action="store", dest="outFileName", type=str, default="drilled.stl", help="Set output file name")
    parser.add_option("-m", "--holesPerSlice", action="store", dest="holesPerSlice", type=int, default=5, help="Set number of holes per slice")
    parser.add_option("-n", "--numOfSlice", action="store", dest="numOfSlice", type=int ,default=5, help="Set number of slices")
    parser.add_option("-r", "--radius", action="store", dest="radius", type=float, default=5, help="Set hole radius required")
    parser.add_option("-p", "--padding", action="store", dest="padding", type=str, default="20,20", help="Set padding level where no holes are drilled")
    parser.add_option("-e", "--errorTorlerance", action="store", dest="error", type=float, default=1, help="Set maximum error tolerance from idea grid in degrees")
    parser.add_option("-a", "--auto", action="store_true", dest="auto", default=False, help="Automatically determine parameters")

    (options, args) = parser.parse_args()
    surfaceFileName = options.surface
    centerlineFileName = options.centerline
    outFileName = options.outFileName
    [startPadding, endPadding] = [int(options.padding.split(",")[i]) for i in xrange(2)]

    try:
        if not os.path.isfile(surfaceFileName):
            raise IOError("Surface file %s dosen't exist!"%surfaceFileName)
        if not os.path.isfile(centerlineFileName):
            raise IOError("Centerline file %s dosen't exist!"%centerlineFileName)

        # create center line object
        cl = CenterLineHandler(centerlineFileName)
        cl.Read()

        # careate arm object
        arm = ArmSurfaceHandler(surfaceFileName, cl)
        arm.Read()

        # Get a list of holes and then drill
        holelist = arm.GetSemiUniDistnaceGrid(options.holesPerSlice, options.numOfSlice - 1, options.error, startPadding, endPadding)
        arm.SphereDrill(holelist, options.radius, options.quiet)

        clippermapper = vtk.vtkPolyDataMapper()
        if vtk.vtkVersion().GetVTKVersion < 6:
            clippermapper.SetInput(arm._data)
        else:
            clippermapper.SetInputData(arm._data)

        writer = vtk.vtkSTLWriter()
        writer.SetFileName(outFileName)
        writer.SetInputData(arm._data)
        if writer.Write() == 1:
            return 0
        else:
            return 1
    except IOError:
        return 2
    except RuntimeError:
        return 3

if __name__ == '__main__':
    main(sys.argv)



