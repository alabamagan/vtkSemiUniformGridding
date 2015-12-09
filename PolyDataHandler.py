#!/usr/bin/python
"""
Author: Wong Matthew Lun
Date: 2015-12-08 6:21PM
"""
import vtk
import math

class CenterLineHandler(vtk.vtkPolyData):
    def __init__(self, filename):
        self.filename = filename
        self._reader = None
        self._renderer = vtk.vtkRenderer()
        self._renderWindow = vtk.vtkRenderWindow()
        self._renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self._IS_READ_FLAG=False

    def Read(self):
        """
        Initial step, should be performed before everything else starts, considering to add
        this into the class constructor so that it is read automatically.
        :return:
        """
        if self.filename.split('.')[-1] == "vtp":
            m_reader = vtk.vtkXMLPolyDataReader()
        m_reader.SetFileName(self.filename)
        m_reader.Update()

        m_mapper = vtk.vtkPolyDataMapper()
        m_mapper.SetInputConnection(m_reader.GetOutputPort())

        m_actor = vtk.vtkActor()
        m_actor.SetMapper(m_mapper)

        self._reader = m_reader
        self._renderer.AddActor(m_actor)
        self._data = m_reader.GetOutput()
        self._IS_READ_FLAG=True
        pass

    def ShowInteractor(self):
        """
        Useless function, For Debug, will delete
        :return:
        """
        self._renderWindow.AddRenderer(self._renderer)
        self._renderWindowInteractor.SetRenderWindow(self._renderWindow)
        self._renderer.SetBackground(0,0,0)
        self._renderer.Render()
        self._renderWindowInteractor.Initialize()
        self._renderWindowInteractor.Start()
        pass

    def WriteImage(self, m_outFileName="./Dump/tmp.png",m_dimension=[400,400]):
        """
        Write current renderer to a png file. For Debug

        :param m_outFileName:   [str] Output name of the file, can be directory name. Default="./Dump/tmp.png"
        :param m_dimension:     [x, y]. Dimension, i.e. width and height of the image file.
        :return:
        """
        if self._renderer.GetRenderWindow() == None:
            self._renderWindow.AddRenderer(self._renderer)

        elif self._renderer.GetRenderWindow() != self._renderWindow:
            self._renderer.GetRenderWindow().Finalize()
            self._renderWindow = vtk.vtkRenderWindow()
            self._renderWindow.AddRenderer(self._renderer)
        else:
            self._renderWindow = vtk.vtkRenderWindow()
            self._renderWindow.AddRenderer(self._renderer)

        self._renderWindow.SetOffScreenRendering(1)
        self._renderWindow.SetSize(m_dimension)
        self._renderWindow.Render()
        self._renderWindow.SetAAFrames(10)

        m_writer = vtk.vtkPNGWriter()
        m_wintoim = vtk.vtkWindowToImageFilter()
        m_wintoim.SetInput(self._renderWindow)
        m_wintoim.Update()

        m_writer.SetInputConnection(m_wintoim.GetOutputPort())
        m_writer.SetFileName(m_outFileName)
        m_writer.Write()

        pass

    def GetData(self):
        """
        Return the data output from the reader.

        Require sequence: Read()

        :return:
        """
        return self._data

    def GetDistance(self, int1, int2):
        """
        Return the distance between two points of the data set, specified by its ID
        Can be replaced by vtkMath() methods
        Require sequence: Read()

        :param int1: [int] ID of the first point
        :param int2: [int] ID of the second point
        :return:
        """
        m_p1, m_p2 = [self._data.GetPoint(i) for i in [int1, int2]]
        d = math.sqrt(sum([(m_p1[i] - m_p2[i])**2 for i in xrange(3)]))
        return d

    def GetEqualDistanceIntervalsIndex(self, m_distance, m_startPadding=0, m_endPadding=0):
        """
        Return a list of index indicating the id of a roughly uni-distance points.

        :param m_distance:      [float] Distance between each segments
        :param m_startingPoint: [int] Index of the user specified starting point. Default=0
        :return:
        """
        m_intevals = []
        m_loopValue = 0
        for i in xrange(m_startPadding+1, self._data.GetNumberOfPoints() - m_endPadding):
            m_loopValue += self.GetDistance(i, i-1)
            if m_loopValue > m_distance:
                m_intevals.append(i)
                m_loopValue = 0
        return m_intevals

    def PrintPoints(self):
        """
        Useless, for DEBUG

        Require sequence: Read()
        :return:
        """
        m_data = self._data
        m_numberOfPoints = m_data.GetNumberOfPoints()
        for i in xrange(m_numberOfPoints):
            print m_data.GetPoint(i)
        pass

    def LoopPoints(self):
        """
        Useless, for DEBUG

        Require sequence: Read()
        :return:
        """
        outList = []
        m_data = self._data
        m_numberOfPoints = m_data.GetNumberOfPoints()
        for i in xrange(m_numberOfPoints):
            outList.append(m_data.GetPoint(i))

        return outList

    def GetNormalizedTangent(self, m_pointID, range=3, step=1):
        m_indexlist = [m_pointID + i for i in xrange(-range, range+1, step)]
        m_preAverage = [0,0,0]
        for i in m_indexlist:
            l_n = [self._data.GetPoint(i)[k] - self._data.GetPoint(i - 1)[k] for k in xrange(3)]
            l_magnitude = math.sqrt(sum([l_n[k]**2 for k in xrange(3)]))
            m_preAverage = [m_preAverage[k] + l_n[k]/l_magnitude for k in xrange(3)]
        m_preAverage = [m_preAverage[k]*step/float(2*range+1) for k in xrange(3)]
        return m_preAverage

    def GetPoint(self, m_int):
        return self._data.GetPoint(m_int)

class ArmSurfaceHandler(vtk.vtkPolyData):
    def __init__(self, filename, centerline):
        """
        Create an ArmSurface object

        :param filename:    STL file of the casting
        :param centerline:  Centerline Object of the casting
        :return:
        """
        self.filename = filename
        self._reader = None
        self._renderer = vtk.vtkRenderer()
        self._renderWindow = vtk.vtkRenderWindow()
        self._renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self._IS_READ_FLAG=False
        self._centerLine = centerline

    def IsRead(self):
        return self._IS_READ_FLAG

    def Read(self):
        """
        Initial step, should be performed before everything else starts, considering to add
        this into the class constructor so that it is read automatically.
        :return:
        """
        if self.filename.split('.')[-1] == "vtp":
            m_reader = vtk.vtkXMLPolyDataReader()
        elif self.filename.split('.')[-1] == "stl":
            m_reader = vtk.vtkSTLReader()
        m_reader.SetFileName(self.filename)
        m_reader.Update()

        m_mapper = vtk.vtkPolyDataMapper()
        m_mapper.SetInputConnection(m_reader.GetOutputPort())

        m_actor = vtk.vtkActor()
        m_actor.SetMapper(m_mapper)

        self._reader = m_reader
        self._renderer.AddActor(m_actor)
        self._data = m_reader.GetOutput()
        self._IS_READ_FLAG = True
        pass

    def GetPoint(self, m_int):
        return self._data.GetPoint(m_int)

    def SliceSurfaceOld(self, m_pt, m_normalVector, m_thickness=0.1):
        """
        Cut surface by checking dot product of the relative position vector from center line m_pt to
        all points and the normal vector of the plane. If the dot product is smaller than a
        given value, that point is accepted

        Replaced by vtkcutter function.

        :param m_pt:            [vtkID] a point of the centerline on the desired cutting plane
        :param m_normalVector:  [x, y, z] normal vector of the desired cutting plane
        :return:
        """
        m_slicepoints = []
        for i in xrange(self._data.GetNumberOfPoints()):
            m_tmp = self._data.GetPoint(i)
            m_ringVector = [m_tmp[k] - m_pt[k] for k in xrange(3)]
            m_dotProduct = sum([m_ringVector[k] * m_normalVector[k] for k in xrange(3)])
            if math.fabs(m_dotProduct) < m_thickness :
                m_slicepoints.append(i)
            # if i%100==0:
            #     print math.fabs(m_dotProduct)
        m_vtkslicepoints = vtk.vtkPolyData()
        m_vtkpts = vtk.vtkPoints()
        for i in m_slicepoints:
            m_vtkpts.InsertNextPoint(self._data.GetPoint(i))
        m_vtkslicepoints.SetPoints(m_vtkpts)
        return m_vtkslicepoints

    def SliceSurface(self, m_pt, m_normalVector):
        """
        Use vtkCutter to obtain cutted.
        :param m_pt:
        :param m_normalVector:
        :return:
        """
        m_plane=vtk.vtkPlane()
        m_plane.SetOrigin(m_pt)
        m_plane.SetNormal(m_normalVector)

        #create cutter
        m_cutter=vtk.vtkCutter()
        m_cutter.SetCutFunction(m_plane)
        m_cutter.SetInputConnection(self._reader.GetOutputPort())
        m_cutter.Update()

        m_vtkpoints = m_cutter.GetOutput().GetPoints()
        return m_vtkpoints

    def GetSemiUniDistnaceGrid(self, m_holePerSlice, m_numberOfSlice, m_errorTolerance=1, m_startPadding = 0, m_endPadding=0):
        """
        Obtain a set of coordinates roughly equal to a projection of periodic square grid vertex on the arm
        surface.

        :param m_holePerSlice:      [int]   Desired number of holes per slice
        :param m_numberOfSlice:     [int]   Desired number of slices
        :param m_errorTolerance:    [float] The maximum allowed deviation of hole coordinate from idea grid
                                            Note that large value causes offsets to inter-hole spacing
        :param m_startPadding:      [int]   Starting side padding where no holes will be drilled
        :param m_endPadding:        [int]   Ending side padding where no holes will be drilled
        :return: [list] List of hole coordinates
        """
        vtkmath = vtk.vtkMath()

        if not self._centerLine._IS_READ_FLAG:
            self._centerLine.Read()

        m_totalDistance = 0
        for i in xrange(1 + m_startPadding, self._centerLine._data.GetNumberOfPoints() - m_endPadding):
            m_totalDistance += self._centerLine.GetDistance(i, i-1)

        m_sliceSpacing = m_totalDistance/(1. + m_numberOfSlice)
        m_intervalIndexes = self._centerLine.GetEqualDistanceIntervalsIndex(m_sliceSpacing)

        m_tangents = []
        for k in xrange(len(m_intervalIndexes)):
            m_tmp = self._centerLine.GetNormalizedTangent(m_intervalIndexes[k], range=30, step=3)
            m_tangents.append(m_tmp)

        m_average = [sum([m_tangents[i][j] for i in xrange(3)])/float(len(m_tangents)) for j in xrange(3)]

        m_holeList = []
        m_alphaNormal = None
        m_masterPt = self._centerLine.GetPoint(m_intervalIndexes[0])
        for i in xrange(len(m_intervalIndexes)):
            l_sliceCenter = self._centerLine.GetPoint(m_intervalIndexes[i])
            l_slice = self.SliceSurface(l_sliceCenter, m_average)
            if i == 0:
                l_ringAlphaPt = l_slice.GetPoint(i)
                l_ringAlphaVect = [l_ringAlphaPt[j] - l_sliceCenter[j] for j in xrange(3)]
                m_alphaNormal = [0,0,0]
                vtkmath.Cross(m_average, l_ringAlphaVect, m_alphaNormal)

            for j in xrange(l_slice.GetNumberOfPoints()):
                l_ringSliceAlphaVect = [l_slice.GetPoint(j)[k] - l_sliceCenter[k] for k in xrange(3)]
                l_ringSliceMasterVect = [l_slice.GetPoint(j)[k] - m_masterPt[k] for k in xrange(3)]
                if math.fabs(vtkmath.Dot(l_ringSliceMasterVect, m_alphaNormal)) < 10 and vtkmath.Dot(l_ringSliceAlphaVect, l_ringAlphaVect) > 0:
                    break

            l_loopbreak = 0
            l_holeList = [[l_ringSliceAlphaVect[k] + l_sliceCenter[k] for k in xrange(3) ]]
            while(len(l_holeList) < m_holePerSlice):
                for j in xrange(l_slice.GetNumberOfPoints()):
                    l_p1 = [0.,0.,0.]
                    l_ringVect = [l_slice.GetPoint(j)[k] - l_sliceCenter[k] for k in xrange(3)]
                    vtkmath.Cross(l_ringSliceAlphaVect, l_ringVect,l_p1)
                    l_p2 = vtkmath.Dot(l_p1, m_average)
                    l_angleBetweenRunningAndInitialVector = vtkmath.AngleBetweenVectors(l_ringSliceAlphaVect, l_ringVect)
                    if l_angleBetweenRunningAndInitialVector > vtkmath.RadiansFromDegrees((360.- m_holePerSlice*m_errorTolerance)/m_holePerSlice) and l_angleBetweenRunningAndInitialVector < vtkmath.RadiansFromDegrees((360. + m_holePerSlice*m_errorTolerance)/m_holePerSlice)and l_p2 > 0:
                        l_ringSliceAlphaVect = l_ringVect
                        l_holeList.append([l_ringVect[k] + l_sliceCenter[k] for k in xrange(3)])
                        break
                if l_loopbreak == m_holePerSlice:
                    raise RuntimeError("Current error tolerence setting is to low to produce anything.")
                l_loopbreak += 1
            m_holeList.extend(l_holeList)

        return m_holeList

    def WriteImage(self, m_outFileName="./Dump/tmp.png",m_dimension=[400,400]):
        """
        Write current renderer to a png file. For Debug

        :param m_outFileName:   [str] Output name of the file, can be directory name. Default="./Dump/tmp.png"
        :param m_dimension:     [x, y]. Dimension, i.e. width and height of the image file.
        :return:
        """
        if self._renderer.GetRenderWindow() == None:
            self._renderWindow.AddRenderer(self._renderer)

        elif self._renderer.GetRenderWindow() != self._renderWindow:
            self._renderer.GetRenderWindow().Finalize()
            self._renderWindow = vtk.vtkRenderWindow()
            self._renderWindow.AddRenderer(self._renderer)
        else:
            self._renderWindow = vtk.vtkRenderWindow()
            self._renderWindow.AddRenderer(self._renderer)

        self._renderWindow.SetOffScreenRendering(1)
        self._renderWindow.SetSize(m_dimension)
        self._renderWindow.Render()
        self._renderWindow.SetAAFrames(10)

        m_writer = vtk.vtkPNGWriter()
        m_wintoim = vtk.vtkWindowToImageFilter()
        m_wintoim.SetInput(self._renderWindow)
        m_wintoim.Update()

        m_writer.SetInputConnection(m_wintoim.GetOutputPort())
        m_writer.SetFileName(m_outFileName)
        m_writer.Write()

        pass

    def Drill(self, m_holelist, m_holeRadius, m_quiet=False):
        m_totalNumOfHoles = len(m_holelist)
        if not m_quiet:
            print "Drilling"

        for i in xrange(m_totalNumOfHoles):
            m_sphere = vtk.vtkSphere()
            m_sphere.SetCenter(m_holelist[i])
            m_sphere.SetRadius(m_holeRadius)


            clipper = vtk.vtkClipPolyData()
            clipper.SetInputData(self._data)
            clipper.SetClipFunction(m_sphere)
            clipper.Update()

            clipped = clipper.GetOutput()
            self._data.DeepCopy(clipped)

            if not m_quiet:
                print "%s/%s -- %.2f %%"%(i+1, m_totalNumOfHoles, (i+1)*100/float(m_totalNumOfHoles))

        pass