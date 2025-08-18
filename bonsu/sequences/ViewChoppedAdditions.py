"""
Created on Fri Aug 15 10:37:56 2025.

@author: davids_work
"""
################################################################################
#                                 Imports                                      #
################################################################################
import wx
import vtk
import math
import numpy
from vtk.util import numpy_support
from ..interface.common import CNTR_CLIP
from ..operations.loadarray import LoadArray, LoadCoordsArray
from ..interface.ViewChoppedCommon import rotate_zhat_by_theta_psi


################################################################################
#                                Constants                                     #
################################################################################


################################################################################
#                                 Methods                                      #
################################################################################
def Sequence_ViewChopped(self, ancestor):
    """Parse and view the Chopped Slices."""
    RSVC: Refactored_Sequence_ViewChopped = \
        Refactored_Sequence_ViewChopped(self, ancestor)
    RSVC.obtain_data_and_coords_filenames()
    RSVC.obtain_background()
    RSVC.load_data()
    RSVC.obtain_slicing_parameters()
    RSVC.report_progress()
    RSVC.setup_vtk_data()
    RSVC.draw_colourbar()
    RSVC.get_source_geometry()
    RSVC.create_slabs()
    RSVC.finalise_rendering()
    return None


################################################################################
#                                 Classes                                      #
################################################################################
class Refactored_Sequence_ViewChopped():
    """Sequence_ViewChopped but refactored into a class."""

    def __init__(self, sequence, ancestor) -> None:
        """Initialise the class."""
        self._panelphase = ancestor.GetPage(0)
        self._panelvisual = ancestor.GetPage(1)
        self._sequence = sequence
        return None

    def obtain_data_and_coords_filenames(self) -> None:
        """Access and retrieve paths to data and coords."""
        self._data_filename: str = \
            self._sequence.input_filename.objectpath.GetValue()
        self._coords_filename: str = \
            self._sequence.coords_filename.objectpath.GetValue()
        return None

    def obtain_background(self) -> None:
        """Retrieve the background colour from PanelVisual."""
        self._background: tuple[float] = \
            (float(self._panelvisual.r)/255.0, float(self._panelvisual.g)/255.0,
             float(self._panelvisual.b)/255.0)
        return None

    def load_data(self) -> None:
        """Load the data and coords into PanelVisual."""
        try:
            data: numpy.ndarray = \
                LoadArray(self._panelphase, self._data_filename)
            data_max: float = numpy.abs(data).max()
            coords: numpy.ndarray = \
                LoadCoordsArray(self._panelphase, self._coords_filename, data)
        except Exception as e:
            print(f'RSVC.load_data() - {type(e)} - {e}')
            data = numpy.random.random((5, 4, 3)).astype(numpy.complex128)
            data_max = numpy.abs(data).max()
            coords = numpy.random.random((data.size, 3))
            dlg: wx.MessageBox = \
                wx.MessageDialog(self._sequence, "Could not load array.",
                                 "Sequence View Object", wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
        finally:
            self._panelvisual.data = data
            self._panelvisual.data_max = data_max
            self._panelvisual.coords = coords
        return None

    def obtain_slicing_parameters(self) -> None:
        """Calculate the slicing parameters and store in memory."""
        theta: float = float(self._sequence.theta.value.GetValue())
        psi: float = float(self._sequence.psi.value.GetValue())
        self._slicing_plane_normal: numpy.ndarray = \
            rotate_zhat_by_theta_psi(theta_deg=theta, psi_deg=psi)
        self._origins: list[tuple[float]] = []
        return None

    def draw_axes_or_not(self) -> None:
        """If ticked box, draw axes in PanelVisual. Else don't."""
        if not self._sequence.chkbox_axes.GetValue():
            return None
        self._panelvisual.axis.SetBounds(self._panelvisual.object_amp.
                                         GetBounds())
        self._panelvisual.axis.SetInputData(self._panelvisual.object_amp)
        self._panelvisual.axis.SetCamera(self._panelvisual.renderer_amp_real.
                                         GetActiveCamera())
        self._panelvisual.axis.SetLabelFormat("%6.1f")
        self._panelvisual.axis.SetFlyModeToOuterEdges()
        self._panelvisual.axis.ScalingOff()
        self._panelvisual.axis. \
            SetFontFactor(float(self._sequence.axes_fontfactor.value.
                                GetValue()))
        self._panelvisual.axis.SetXLabel("X")
        self._panelvisual.axis.SetYLabel("Y")
        self._panelvisual.axis.SetZLabel("Z")
        self._panelvisual.axis.Modified()
        self._panelvisual.renderer_amp_real.AddViewProp(self._panelvisual.axis)
        return None

    def report_progress(self) -> None:
        """Log information."""
        self._panelphase.queue_info.put("Preparing object chopping.")
        self._panelphase.queue_info. \
            put(''.join(("Array Size: ", str(self._panelvisual.data.shape))))
        return None

    def setup_vtk_data(self) -> None:
        """Set up the vtk data only once."""
        # Extract the phase info and flatten (VTK preferred)
        self._panelvisual.flat_data_phase: numpy.ndarray = \
            (numpy.angle(self._panelvisual.data)).transpose(2, 1, 0).flatten()
        # Convert to VTK's native data type
        self._panelvisual.vtk_data_array_phase: numpy.ndarray = \
            numpy_support.numpy_to_vtk(self._panelvisual.flat_data_phase)
        # Label for quick access (?)
        self._panelvisual.vtk_data_array_phase.SetName("mapscalar")
        # Repeat for amplitude.
        self._panelvisual.flat_data: numpy.ndarray = \
            (numpy.abs(self._panelvisual.data)).transpose(2, 1, 0).flatten()
        self._panelvisual.vtk_data_array: numpy.ndarray = \
            numpy_support.numpy_to_vtk(self._panelvisual.flat_data)
        # Similar operation for the coords.
        self._panelvisual.vtk_coordarray: numpy.ndarray = \
            numpy_support.numpy_to_vtk(self._panelvisual.coords)
        # Prepare the coordinate list to store highly accurate numbers.
        self._panelvisual.object_amp_points.SetDataTypeToDouble()
        # Memory-optimised coordinate list preparation.
        self._panelvisual.object_amp_points. \
            SetNumberOfPoints(self._panelvisual.data.size)
        # Assign the coordinate array to the vtk.
        self._panelvisual.object_amp_points. \
            SetData(self._panelvisual.vtk_coordarray)
        # Just store the shape locally and make more accessible.
        self._data_shp: numpy.ndarray = \
            numpy.array(self._panelvisual.data.shape, dtype=numpy.int64)
        # Assign 3D coordinates to 3D data grid.
        self._panelvisual.object_amp. \
            SetPoints(self._panelvisual.object_amp_points)
        # Assign the amplitude values as the main data.
        self._panelvisual.object_amp.GetPointData(). \
            SetScalars(self._panelvisual.vtk_data_array)
        # Assign phase data as extra data.
        self._panelvisual.object_amp.GetPointData().\
            AddArray(self._panelvisual.vtk_data_array_phase)
        # Structure the flat data points into a 3D grid of specific dimensions.
        self._panelvisual.object_amp.SetDimensions(self._data_shp)
        # Inform VTK of the changed data object.
        self._panelvisual.object_amp.Modified()
        # Repeat for phase.
        self._panelvisual.object_phase. \
            SetPoints(self._panelvisual.object_amp_points)
        self._panelvisual.object_phase.GetPointData(). \
            SetScalars(self._panelvisual.vtk_data_array_phase)
        self._panelvisual.object_phase.SetDimensions(self._data_shp)
        self._panelvisual.object_phase.Modified()
        objectnpoints = self._panelvisual.flat_data_phase.size
        objectbounds = self._panelvisual.object_amp.GetBounds()
        density = (objectbounds[1] - objectbounds[0]) * \
            (objectbounds[3] - objectbounds[2]) * \
            (objectbounds[5] - objectbounds[4])/objectnpoints
        linedensity = math.pow(density, 1.0/3.0)
        self._mesh_maxedge = 20.0*linedensity
        return None

    def draw_colourbar(self) -> None:
        """Draw the colourbar."""
        self._panelvisual.image_probe: numpy.ndarray = \
            self._panelvisual.object_phase
        self._panelvisual.lut_phase_real.SetNumberOfTableValues(256)
        self._panelvisual.lut_phase_real.SetScaleToLinear()
        self._panelvisual.lut_phase_real. \
            SetTableRange([float(self._sequence.phasemin.value.GetValue()),
                           float(self._sequence.phasemax.value.GetValue())])
        lutsource: numpy.ndarray = \
            self._panelphase.cms[self._panelphase.cmls[1][0]][1]
        if self._panelphase.cmls[1][1] == 0:
            for k in range(256):
                self._panelvisual.lut_phase_real. \
                    SetTableValue(k, lutsource[k][0], lutsource[k][1],
                                  lutsource[k][2], 1)
        else:
            for k in range(256):
                self._panelvisual.lut_phase_real. \
                    SetTableValue(255-k, lutsource[k][0], lutsource[k][1],
                                  lutsource[k][2], 1)
        self._panelvisual.lut_phase_real.SetRamp(0)
        self._panelvisual.lut_phase_real.Build()
        self._panelvisual.lut_amp_real.SetNumberOfTableValues(256)
        self._panelvisual.lut_amp_real.SetScaleToLinear()
        self._panelvisual.lut_amp_real. \
            SetTableRange(self._panelvisual.object_amp.GetPointData().
                          GetScalars().GetRange())
        lutsource: numpy.ndarray = \
            self._panelphase.cms[self._panelphase.cmls[0][0]][1]
        if self._panelphase.cmls[0][1] == 0:
            for k in range(256):
                self._panelvisual.lut_amp_real. \
                    SetTableValue(k, lutsource[k][0], lutsource[k][1],
                                  lutsource[k][2], 1)
        else:
            for k in range(256):
                self._panelvisual.lut_amp_real. \
                    SetTableValue(255-k, lutsource[k][0], lutsource[k][1],
                                  lutsource[k][2], 1)
        self._panelvisual.lut_amp_real.SetRamp(0)
        self._panelvisual.lut_amp_real.Build()
        self._panelvisual.scalebar_amp_real.SetTitle("")
        self._panelvisual.scalebar_amp_real. \
            SetLookupTable(self._panelvisual.lut_phase_real)
        return None

    def get_source_geometry(self) -> None:
        """Return the final, unclipped, renderable geometry."""
        # Contour!
        self._panelvisual.filter_amp_real. \
            SetInputData(self._panelvisual.object_amp)
        self._panelvisual.filter_amp_real.ComputeNormalsOn()
        self._panelvisual.filter_amp_real.ComputeScalarsOn()
        self._panelvisual.filter_amp_real.SetNumberOfContours(1)
        self._contour: float = float(self._sequence.contour.value.GetValue())
        if self._contour > self._panelvisual.data_max:
            self._contour = CNTR_CLIP*self._panelvisual.data_max
        self._panelvisual.filter_amp_real. \
            SetValue(0, float(self._contour))
        self._panelvisual.filter_amp_real.Modified()
        self._panelvisual.filter_amp_real.Update()
        # Smooth!
        self._panelvisual.smooth_filter_real. \
            SetInputConnection(self._panelvisual.filter_amp_real.
                               GetOutputPort())
        self._panelvisual.smooth_filter_real.SetNumberOfIterations(15)
        self._panelvisual.smooth_filter_real.SetRelaxationFactor(0.1)
        self._panelvisual.smooth_filter_real.FeatureEdgeSmoothingOff()
        self._panelvisual.smooth_filter_real.BoundarySmoothingOn()
        self._panelvisual.smooth_filter_real.Update()
        # Calculate the normals for each point to define reactions to light
        self._panelvisual.normals_amp_real.\
            SetInputConnection(self._panelvisual.smooth_filter_real.
                               GetOutputPort())
        self._panelvisual.normals_amp_real. \
            SetFeatureAngle(float(self._sequence.feature_angle.value.
                                  GetValue()))
        self._panelvisual.normals_amp_real.ConsistencyOff()
        self._panelvisual.normals_amp_real.SplittingOff()
        self._panelvisual.normals_amp_real.AutoOrientNormalsOff()
        self._panelvisual.normals_amp_real.ComputePointNormalsOn()
        self._panelvisual.normals_amp_real.ComputeCellNormalsOff()
        self._panelvisual.normals_amp_real.NonManifoldTraversalOff()
        # Optimise
        self._panelvisual.triangles_amp_real. \
            SetInputConnection(self._panelvisual.normals_amp_real.
                               GetOutputPort())
        self._panelvisual.strips_amp_real. \
            SetInputConnection(self._panelvisual.triangles_amp_real.
                               GetOutputPort())
        self._panelvisual.strips_amp_real.Update()
        return None

    def create_slab_actors(self, origin_plane_01: numpy.ndarray) -> None:
        """Return the actors of a singular slab."""
        self._contour: float = float(self._sequence.contour.value.GetValue())
        if self._contour > self._panelvisual.data_max:
            self._contour = CNTR_CLIP*self._panelvisual.data_max
        # Determine the origin of the second cut plane.
        origin_plane_02: numpy.ndarray = origin_plane_01 + \
            self._slicing_plane_normal * \
            float(self._sequence.thickness.value.GetValue())
        # Put slab-cutting planes' origins into a singular variable.
        slab_cutting_origins: tuple[numpy.ndarray] = \
            (origin_plane_01, origin_plane_02)
        # Initialise the cutting planes.
        slicing_planes: tuple[vtk.vtkPlane] = (vtk.vtkPlane(), vtk.vtkPlane())
        # Iterate through the cutting planes and attribute origins and normals.
        for plane_index, plane in enumerate(slicing_planes):
            plane.SetOrigin(*(slab_cutting_origins[plane_index]))
            plane.SetNormal(*(self._slicing_plane_normal *
                              numpy.power(-1, plane_index+1)))
        # Calculate the space between the planes.
        boolean_operation: vtk.vtkImplicitBoolean = vtk.vtkImplicitBoolean()
        boolean_operation.SetOperationTypeToIntersection()
        [boolean_operation.AddFunction(plane) for plane in slicing_planes]
        # Clip the triangles in the space between the planes.
        clipper: vtk.vtkClipPolyData = vtk.vtkClipPolyData()
        clipper.SetInputConnection(self._panelvisual.strips_amp_real.
                                   GetOutputPort())
        clipper.SetClipFunction(boolean_operation)
        clipper.GenerateClippedOutputOn()
        clipper.SetInsideOut(True)
        clipper.SetValue(0)
        clipper.Update()
        # Initialise the amplitude mappers between and outsite de planes.
        amp_mappers: list[vtk.vtkPolyDataMapper] = \
            [vtk.vtkPolyDataMapper(), vtk.vtkPolyDataMapper()]
        amp_actors: list[vtk.vtkActor] = []
        # Iterate through and amplitude mappers to draw them.
        for amp_mapper_index, amp_mapper in enumerate(amp_mappers):
            amp_mapper.SetInputConnection(clipper.GetOutputPort())
            amp_mapper.SetLookupTable(self._panelvisual.lut_amp_real)
            amp_mapper.SetScalarRange(self._panelvisual.object_amp.
                                      GetPointData().GetScalars().GetRange())
            amp_mapper.SetScalarModeToUsePointData()
            amp_mapper.Modified()
            amp_mapper.Update()
            # ... and initialise their respective actors.
            amp_actors.append(vtk.vtkActor())
        # Iterate through the amplitude actors and attribute their mappers.
        for amp_actor_index, amp_actor in enumerate(amp_actors):
            amp_actor.GetProperty().SetOpacity(amp_mapper_index)
            amp_actor.SetMapper(amp_mappers[amp_actor_index])
        # Something similar for phase.
        filter_plane: vtk.vtkContourFilter = vtk.vtkContourFilter()
        filter_plane.SetInputData(self._panelvisual.object_amp)
        filter_plane.ComputeNormalsOn()
        filter_plane.ComputeScalarsOn()
        filter_plane.SetNumberOfContours(1)
        filter_plane.SetValue(0, float(self._contour))
        filter_plane.Modified()
        filter_plane.Update()
        smooth_plane: vtk.vtkSmoothPolyDataFilter = \
            vtk.vtkSmoothPolyDataFilter()
        smooth_plane.SetInputConnection(filter_plane.GetOutputPort())
        smooth_plane.SetNumberOfIterations(15)
        smooth_plane.SetRelaxationFactor(0.1)
        smooth_plane.FeatureEdgeSmoothingOff()
        smooth_plane.BoundarySmoothingOn()
        smooth_plane.Modified()
        smooth_plane.Update()
        cutters: list[vtk.vtkCutter] = [vtk.vtkCutter(), vtk.vtkCutter()]
        filter_tris: list[vtk.vtkContourTriangulator] = []
        for cutter_index, cutter in enumerate(cutters):
            cutter.SetInputConnection(smooth_plane.GetOutputPort())
            cutter.SetCutFunction(slicing_planes[cutter_index])
            cutter.GenerateTrianglesOn()
            cutter.Update()
            cutter.GenerateCutScalarsOn()
            filter_tris.append(vtk.vtkContourTriangulator())
        meshsubs: list[vtk.vtkAdaptiveSubdivisionFilter] = []
        for filter_tri_index, filter_tri in enumerate(filter_tris):
            filter_tri.SetInputConnection(cutters[filter_tri_index].
                                          GetOutputPort())
            meshsubs.append(vtk.vtkAdaptiveSubdivisionFilter())
        probefilters: list[vtk.vtkProbeFilter] = []
        for meshsub_index, meshsub in enumerate(meshsubs):
            meshsub.SetInputConnection(filter_tris[meshsub_index].
                                       GetOutputPort())
            meshsub.SetMaximumEdgeLength(self._mesh_maxedge)
            meshsub. \
                SetMaximumNumberOfPasses(int(float(self._sequence.meshsubiter.
                                                   value.GetValue())))
            meshsub.SetOutputPointsPrecision(vtk.vtkAlgorithm.SINGLE_PRECISION)
            meshsub.Update()
            probefilters.append(vtk.vtkProbeFilter())
        triangles_planes: list[vtk.vtkTriangleFilter] = []
        for probefilter_index, probefilter in enumerate(probefilters):
            probefilter. \
                SetInputConnection(meshsubs[probefilter_index].GetOutputPort())
            probefilter.SetSourceData(self._panelvisual.object_phase)
            probefilter.Update()
            triangles_planes.append(vtk.vtkTriangleFilter())
        normals_phases: list[vtk.vtkPolyDataNormals] = []
        for tri_plane_index, tri_plane in enumerate(triangles_planes):
            tri_plane.SetInputConnection(probefilters[tri_plane_index].
                                         GetOutputPort())
            normals_phases.append(vtk.vtkPolyDataNormals())
        triangles_phases: list[vtk.vtkTriangleFilter] = []
        for normal_phase_index, normal_phase in enumerate(normals_phases):
            normal_phase. \
                SetInputConnection(triangles_planes[normal_phase_index].
                                   GetOutputPort())
            normal_phase.SetFeatureAngle(float(self._sequence.feature_angle.
                                               value.GetValue()))
            normal_phase.ConsistencyOff()
            normal_phase.SplittingOff()
            normal_phase.AutoOrientNormalsOff()
            normal_phase.ComputePointNormalsOn()
            normal_phase.ComputeCellNormalsOn()
            normal_phase.NonManifoldTraversalOff()
            triangles_phases.append(vtk.vtkTriangleFilter())
        strips_phases: list[vtk.vtkStripper] = []
        for tri_phase_index, tri_phase in enumerate(triangles_phases):
            tri_phase.SetInputConnection(normals_phases[tri_phase_index].
                                         GetOutputPort())
            strips_phases.append(vtk.vtkStripper())
        phase_mappers: list[vtk.vtkPolyDataMapper] = []
        for strip_phase_index, strip_phase in enumerate(strips_phases):
            strip_phase.SetInputConnection(triangles_phases[strip_phase_index].
                                           GetOutputPort())
            phase_mappers.append(vtk.vtkPolyDataMapper())
        phase_actors: list[vtk.vtkActor] = []
        for phase_mapper_index, phase_mapper in enumerate(phase_mappers):
            phase_mapper.SetInputConnection(strips_phases[phase_mapper_index].
                                            GetOutputPort())
            phase_mapper.SetLookupTable(self._panelvisual.lut_phase_real)
            phase_mapper. \
                SetScalarRange([float(self._sequence.phasemin.value.GetValue()),
                                float(self._sequence.phasemax.value.
                                      GetValue())])
            phase_mapper.SetScalarModeToUsePointData()
            phase_mapper.Update()
            phase_mapper.Modified()
            phase_actors.append(vtk.vtkActor())
        for actor_index, actor in enumerate(phase_actors):
            actor.GetProperty().SetOpacity(1.0)
            actor.SetMapper(phase_mappers[actor_index])
        return (amp_actors, phase_actors)

    def create_slabs(self) -> None:
        number_of_slices: int = \
            int(self._sequence.number_of_slices.value.GetValue())
        thickness: float = \
            float(self._sequence.thickness.value.GetValue())
        separation: float = float(self._sequence.separation.value.GetValue())
        center = self._panelvisual.object_amp.GetCenter()
        total_slice_coverage = number_of_slices * thickness
        start_point: numpy.ndarray = numpy.array(center) - \
            total_slice_coverage * self._slicing_plane_normal / 2
        # --- Main Loop to Create Slabs ---
        for slice_index in range(number_of_slices):
            # Calculate the origin for this specific slab
            slab_origin: numpy.ndarray = start_point + \
                (slice_index * thickness) * self._slicing_plane_normal
            # Create the actors for this one slab
            amp_actors, phase_actors = self.create_slab_actors(slab_origin)
            # --- Apply Separation Transform ---
            transform = vtk.vtkTransform()
            offset = (slice_index - (number_of_slices - 1) / 2.0) * \
                separation * self._slicing_plane_normal
            transform.Translate(offset)
            # Apply the transform to ALL actors for this slab
            for actor in amp_actors + phase_actors:
                actor.SetUserTransform(transform)
            # --- Add the finished, transformed actors to the main list ---
            for actor in amp_actors + phase_actors:
                self._panelvisual.actors.append(actor)
        return None

    def finalise_rendering(self) -> None:
        """Finalise the rendering and release buttons."""
        self._panelvisual.SetPicker()
        self._panelvisual.renderer_amp_real.RemoveAllViewProps()
        self._panelvisual.renderer_phase_real.RemoveAllViewProps()
        self._panelvisual.renderer_amp_recip.RemoveAllViewProps()
        self._panelvisual.renderer_phase_recip.RemoveAllViewProps()
        self._panelvisual.renWin.GetRenderWindow(). \
            RemoveRenderer(self._panelvisual.renderer_amp_real)
        self._panelvisual.renWin.GetRenderWindow(). \
            RemoveRenderer(self._panelvisual.renderer_phase_real)
        self._panelvisual.renWin.GetRenderWindow(). \
            RemoveRenderer(self._panelvisual.renderer_amp_recip)
        self._panelvisual.renWin.GetRenderWindow(). \
            RemoveRenderer(self._panelvisual.renderer_phase_recip)
        self._panelvisual.renWin.GetRenderWindow().Modified()
        for actor in self._panelvisual.actors:
            self._panelvisual.renderer_amp_real.AddActor(actor)
        self._panelvisual.renderer_amp_real. \
            AddActor2D(self._panelvisual.scalebar_amp_real)
        self._panelvisual.renderer_amp_real. \
            AddActor2D(self._panelvisual.textActor)
        self._panelvisual.renderer_amp_real.SetBackground(*(self._background))
        self._panelvisual.renWin.GetRenderWindow(). \
            AddRenderer(self._panelvisual.renderer_phase_real)
        self._panelvisual.renWin.GetRenderWindow(). \
            AddRenderer(self._panelvisual.renderer_amp_real)
        self._panelvisual.renWin.GetRenderWindow().SetMultiSamples(0)
        self._panelvisual.renWin.SetInteractorStyle(self._panelvisual.style3D)
        self._panelvisual.renderer_amp_real.ResetCamera()
        self._panelvisual.renWin.SetPicker(self._panelvisual.picker_amp_real)
        self._panelvisual.renderer_amp_real.SetViewport(0, 0, 1, 1)
        self._panelvisual.renderer_phase_real.SetViewport(1, 1, 1, 1)
        self.draw_axes_or_not()
        self._panelvisual.Layout()
        self._panelvisual.Show()
        self._panelvisual.datarangelist = []
        self._panelvisual.ReleaseVisualButtons(gotovisual=True)
        self._panelvisual.button_vremove.Enable(False)
        self._panelvisual.RefreshSceneFull()
        return None


################################################################################
#                         Dunder Name Dunder Main                              #
################################################################################
if __name__ == '__main__':
    pass
