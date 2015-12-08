#!/usr/bin/python
"""
Author: Wong Matthew Lun
Date: 2015-12-08 6:21PM
"""

from PolyDataHandler import CenterLineHandler, ArmSurfaceHandler
import vtk
import math

def main():
    vtkmath = vtk.vtkMath()

    cl = CenterLineHandler("./Data/centerline.vtp")
    cl.Read()

    # get all points roughly on the given surface
    arm = ArmSurfaceHandler("./Data/clipped_surface.stl", cl)
    arm.Read()
    holelist = arm.GetSemiUniDistnaceGrid(10, 7, 1)

    for j in xrange(len(holelist)): # TODO: Replace 7 with arbitary value
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(holelist[j])
        sphere.SetRadius(5)
        sphereMapper = vtk.vtkPolyDataMapper()
        sphereMapper.SetInputConnection(sphere.GetOutputPort())

        sphereActor = vtk.vtkActor()
        sphereActor.SetMapper(sphereMapper)
        if j % 10 == 0:
            sphereActor.GetProperty().SetColor(0,0.5,0.5)
        else:
            sphereActor.GetProperty().SetColor(1,0,0)

        arm._renderer.AddActor(sphereActor)


    print len(holelist)
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(arm._renderer)


    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    renWin.Render()
    iren.Initialize()
    iren.Start()

if __name__ == '__main__':
    main()



