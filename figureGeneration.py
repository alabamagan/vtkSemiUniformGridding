#!/bin/bash

import PolyDataHandler
import vtk
import numpy as np
import numpy.linalg as linalg

def CreateArrowActor(start, end):
    if len(start) != len(end):
        raise IndexError
        return

    start = np.array(start)
    end = np.array(end)

    arbitrary = np.random.rand(3)
    arrowSource = vtk.vtkArrowSource()
    arrowSource.SetTipLength(0.05)
    arrowSource.SetShaftRadius(0.01)
    arrowSource.SetTipRadius(0.025)

    normalizedX = end - start
    length = linalg.norm(normalizedX)
    normalizedX /= linalg.norm(normalizedX)

    normalizedZ = np.cross(normalizedX, arbitrary)
    normalizedZ /= linalg.norm(normalizedZ)

    normalizedY = np.cross(normalizedZ, normalizedX)

    matrix = vtk.vtkMatrix4x4()
    for i in xrange(3):
        matrix.SetElement(i, 0, normalizedX[i])
        matrix.SetElement(i, 1, normalizedY[i])
        matrix.SetElement(i, 2, normalizedZ[i])

    transform = vtk.vtkTransform()
    transform.Translate(start)
    transform.Concatenate(matrix)
    transform.Scale(length, length, length)

    transformPD = vtk.vtkTransformPolyDataFilter()
    transformPD.SetTransform(transform)
    transformPD.SetInputConnection(arrowSource.GetOutputPort())

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(transformPD.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    return actor

def FindClosestPoint(point, dataset):
    m_kdtree = vtk.vtkKdTreePointLocator()
    m_kdtree.SetDataSet(dataset)
    m_kdtree.BuildLocator()
    m_closestCenterlinePointId = m_kdtree.FindClosestPoint(point)
    m_closestCenterlinePoint = dataset.GetPoint(m_closestCenterlinePointId)
    return m_closestCenterlinePoint, m_closestCenterlinePointId

def figureX():
    openingmarker = [-24.877422332763672, 2.268731117248535, 7.306049823760986]
    userSelectedPoint = [42,42,51]
    lightgreen = [0.4,0.9,0.4]
    lightred = [0.9,0.3,0.3]

    cl = PolyDataHandler.CenterLineHandler("data/Centerline_extended.vtp")
    cl.Read()
    print "Reading centerline finished..."
    cl._actor.GetProperty().SetColor(0.1,0.3,1)

    arm = PolyDataHandler.ArmSurfaceHandler("data/Flare.stl", cl, openingmarker)
    arm.Read()
    print "Reading arm surface finished..."

    userSelectedPointSnap = FindClosestPoint(userSelectedPoint, arm._data)
    userSelectedPointSnapActor = arm.GetPointActor(userSelectedPointSnap[1])
    userSelectedCenterPoint = FindClosestPoint(userSelectedPointSnap[0], cl._data)

    userSelectedPointArrow = CreateArrowActor(userSelectedCenterPoint[0], userSelectedPointSnap[0])
    userSelectedPointSnapActor.GetProperty().SetColor(lightgreen)
    userSelectedPointArrow.GetProperty().SetColor(lightgreen)

    alphaVectorCenter = FindClosestPoint(openingmarker, cl._data)
    alphaVectorArrow = CreateArrowActor(alphaVectorCenter[0], openingmarker)
    alphaVectorArrow.GetProperty().SetColor(lightred)

    arm._actor.GetProperty().SetColor(0.1,0.1,0.1)
    arm._actor.GetProperty().SetOpacity(0.1)
    arm._actor.GetProperty().SetShading(1)
    arm._actor.GetProperty().SetAmbient(0)
    arm._actor.GetProperty().SetDiffuse(0.3)
    arm._actor.GetProperty().SetSpecular(1)
    arm._actor.GetProperty().SetRepresentationToWireframe()
    # arm._renderer.AddActor(userSelectedPointArrow)
    # arm._renderer.AddActor(alphaVectorArrow)
    arm._renderer.AddActor(cl._actor)
    # arm._renderer.AddActor(userSelectedPointSnapActor)
    arm._renderer.SetBackground(1, 1, 1)

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(arm._renderer)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renwin)
    iren.Initialize()
    renwin.Render()
    iren.Start()

def figureSlice():
    openingmarker = [-24.877422332763672, 2.268731117248535, 7.306049823760986]
    userSelectedPoint = [42,42,51]
    lightgreen = [0.4,0.9,0.4]
    lightred = [0.9,0.3,0.3]

    cl = PolyDataHandler.CenterLineHandler("data/Centerline_extended.vtp")
    cl.Read()
    print "Reading centerline finished..."
    cl._actor.GetProperty().SetColor(0.1,0.3,1)

    arm = PolyDataHandler.ArmSurfaceHandler("data/Flare.stl", cl, openingmarker)
    arm.Read()
    print "Reading arm surface finished..."

    sliceactors = []
    arm.GetSemiUniDistnaceGrid(10, 9, 5, 25, 25, 40)
    for clinterval in arm._centerLineIntervals:
        cutter = arm.SliceSurfaceCutter(cl.GetPoint(clinterval), arm._averageTangent)
        sm = vtk.vtkPolyDataMapper()
        sm.SetInputConnection(cutter.GetOutputPort())
        planeActor=vtk.vtkActor()
        planeActor.GetProperty().SetColor(0.7, 0, 0)
        planeActor.GetProperty().SetLineWidth(2)
        planeActor.SetMapper(sm)
        sliceactors.append(planeActor)

    arm._actor.GetProperty().SetColor(0.1,0.1,0.1)
    arm._actor.GetProperty().SetOpacity(0.1)
    arm._actor.GetProperty().SetShading(1)
    arm._actor.GetProperty().SetAmbient(0)
    arm._actor.GetProperty().SetDiffuse(0.3)
    arm._actor.GetProperty().SetSpecular(1)
    arm._actor.GetProperty().SetRepresentationToWireframe()
    for actors in sliceactors:
        arm._renderer.AddActor(actors)
    # arm._renderer.AddActor(userSelectedPointArrow)
    # arm._renderer.AddActor(alphaVectorArrow)
    arm._renderer.AddActor(cl._actor)
    # arm._renderer.AddActor(userSelectedPointSnapActor)
    arm._renderer.SetBackground(1, 1, 1)

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(arm._renderer)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renwin)
    iren.Initialize()
    renwin.Render()
    iren.Start()

def figureSphere():
    openingmarker = [-24.877422332763672, 2.268731117248535, 7.306049823760986]
    userSelectedPoint = [42,42,51]
    lightgreen = [0.4,0.9,0.4]
    lightred = [0.9,0.3,0.3]

    cl = PolyDataHandler.CenterLineHandler("data/Centerline_extended.vtp")
    cl.Read()
    print "Reading centerline finished..."
    cl._actor.GetProperty().SetColor(0.1,0.3,1)

    arm = PolyDataHandler.ArmSurfaceHandler("data/Flare.stl", cl, openingmarker)
    arm.Read()
    print "Reading arm surface finished..."

    sphereActors = []
    arm.GetSemiUniDistnaceGrid(18, 9, 5, 25, 25, 40)
    for i in xrange(17):
        ssource = vtk.vtkSphereSource()
        ssource.SetCenter(arm._holeList[i])
        ssource.SetRadius(3.)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(ssource.GetOutputPort())
        sactor = vtk.vtkActor()
        sactor.SetMapper(mapper)
        sactor.GetProperty().SetColor(0.8,0,0)
        sphereActors.append(sactor)


    arm._actor.GetProperty().SetColor(0.1,0.1,0.1)
    arm._actor.GetProperty().SetOpacity(0.1)
    arm._actor.GetProperty().SetShading(1)
    arm._actor.GetProperty().SetAmbient(0)
    arm._actor.GetProperty().SetDiffuse(0.3)
    arm._actor.GetProperty().SetSpecular(1)
    arm._actor.GetProperty().SetRepresentationToWireframe()
    for actors in sphereActors:
        arm._renderer.AddActor(actors)
    # arm._renderer.AddActor(userSelectedPointArrow)
    # arm._renderer.AddActor(alphaVectorArrow)
    arm._renderer.AddActor(cl._actor)
    # arm._renderer.AddActor(userSelectedPointSnapActor)
    arm._renderer.SetBackground(1, 1, 1)

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(arm._renderer)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renwin)
    iren.Initialize()
    renwin.Render()
    iren.Start()


def figureDrilled():
    openingmarker = [-24.877422332763672, 2.268731117248535, 7.306049823760986]
    userSelectedPoint = [42,42,51]
    lightgreen = [0.4,0.9,0.4]
    lightred = [0.9,0.3,0.3]

    cl = PolyDataHandler.CenterLineHandler("data/Centerline_extended.vtp")
    cl.Read()
    print "Reading centerline finished..."
    cl._actor.GetProperty().SetColor(0.1,0.3,1)

    arm = PolyDataHandler.ArmSurfaceHandler("data/Flare.stl", cl, openingmarker)
    arm.Read()
    print "Reading arm surface finished..."

    arm.GetSemiUniDistnaceGrid(18, 9, 5, 25, 25, 40)


    arm.SphereDrill(arm._holeList, 3)

    arm._actor.GetProperty().SetColor(0.1,0.1,0.1)
    arm._actor.GetProperty().SetOpacity(1)
    arm._actor.GetProperty().SetShading(1)
    arm._actor.GetProperty().SetAmbient(0)
    arm._actor.GetProperty().SetDiffuse(0.3)
    arm._actor.GetProperty().SetSpecular(1)
    # arm._actor.GetProperty().SetRepresentationToWireframe()

    # arm._renderer.AddActor(userSelectedPointArrow)
    # arm._renderer.AddActor(alphaVectorArrow)
    arm._renderer.AddActor(cl._actor)
    # arm._renderer.AddActor(userSelectedPointSnapActor)
    arm._renderer.SetBackground(1, 1, 1)

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(arm._renderer)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renwin)
    iren.Initialize()
    renwin.Render()
    iren.Start()

def main():
    openingmarker = [-24.877422332763672, 2.268731117248535, 7.306049823760986]
    userSelectedPoint = [42,42,51]
    lightgreen = [0.4,0.9,0.4]
    lightred = [0.9,0.3,0.3]

    cl = PolyDataHandler.CenterLineHandler("data/Centerline_extended.vtp")
    cl.Read()
    print "Reading centerline finished..."
    cl._actor.GetProperty().SetColor(0.1,0.3,1)

    arm = PolyDataHandler.ArmSurfaceHandler("data/Flare.stl", cl, openingmarker)
    arm.Read()
    print "Reading arm surface finished..."

    holelist = arm.GetSemiUniDistnaceGrid(22, 9, 5, 25, 25)
    print holelist[0]

    alphaVectorCenter = FindClosestPoint(openingmarker, cl._data)
    alphaVectorArrow = CreateArrowActor(alphaVectorCenter[0], openingmarker)
    alphaVectorArrow.GetProperty().SetColor(lightgreen)

    holesVectorActors = []
    for i in xrange(19):
        if i == 0:
            subVectorCenter = FindClosestPoint(holelist[i], cl._data)

        subVectorArros = CreateArrowActor(subVectorCenter[0], holelist[i+1])
        subVectorArros.GetProperty().SetColor(lightred)
        holesVectorActors.append(subVectorArros)

    arm._actor.GetProperty().SetColor(0.1,0.1,0.1)
    arm._actor.GetProperty().SetOpacity(0.1)
    arm._actor.GetProperty().SetShading(1)
    arm._actor.GetProperty().SetAmbient(0)
    arm._actor.GetProperty().SetDiffuse(0.3)
    arm._actor.GetProperty().SetSpecular(1)
    arm._actor.GetProperty().SetRepresentationToWireframe()

    arm._renderer.AddActor(alphaVectorArrow)
    arm._renderer.AddActor(cl._actor)
    for i in holesVectorActors:
        arm._renderer.AddActor(i)
    arm._renderer.SetBackground(1, 1, 1)

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(arm._renderer)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renwin)
    iren.Initialize()
    renwin.Render()
    iren.Start()


if __name__ == '__main__':
    figureDrilled()



