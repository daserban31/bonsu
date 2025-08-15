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
    RSVC.old_ViewDataAmpClippedPhase()
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

    def old_ViewDataAmpClippedPhase(self) -> None:
        self._panelphase.queue_info.put("Preparing object chopping.")
        self._panelphase.queue_info. \
            put(''.join(("Array Size: ", str(self._panelvisual.data.shape))))
        number_of_slices: int = \
            int(self._sequence.number_of_slices.value.GetValue())
        separation: float = float(self._sequence.separation.value.GetValue())
        meshsubiter: int = int(float(self._sequence.meshsubiter.value.GetValue()))
        opacity: float = float(self._sequence.opacity.value.GetValue())
        contour: float = float(self._sequence.contour.value.GetValue())
        if contour > self._panelvisual.data_max:
            contour = CNTR_CLIP*self._panelvisual.data_max
        feature_angle: float = float(self._sequence.feature_angle.value.GetValue())
        phasemax: float = float(self._sequence.phasemax.value.GetValue())
        phasemin: float = float(self._sequence.phasemin.value.GetValue())

        self._panelvisual.flat_data_phase: numpy.ndarray = \
            (numpy.angle(self._panelvisual.data)).transpose(2, 1, 0).flatten()
        self._panelvisual.vtk_data_array_phase: numpy.ndarray = \
            numpy_support.numpy_to_vtk(self._panelvisual.flat_data_phase)
        self._panelvisual.vtk_data_array_phase.SetName("mapscalar")
        shp: numpy.ndarray = numpy.array(self._panelvisual.data.shape, dtype=numpy.int64)
        self._panelvisual.flat_data: numpy.ndarray = \
            (numpy.abs(self._panelvisual.data)).transpose(2, 1, 0).flatten()
        self._panelvisual.vtk_data_array: numpy.ndarray = \
            numpy_support.numpy_to_vtk(self._panelvisual.flat_data)
        self._panelvisual.vtk_coordarray: numpy.ndarray = \
            numpy_support.numpy_to_vtk(self._panelvisual.coords)
        self._panelvisual.object_amp_points.SetDataTypeToDouble()
        self._panelvisual.object_amp_points.SetNumberOfPoints(self._panelvisual.data.size)
        self._panelvisual.object_amp_points.SetData(self._panelvisual.vtk_coordarray)
        self._panelvisual.object_amp.SetPoints(self._panelvisual.object_amp_points)
        self._panelvisual.object_amp.GetPointData(). \
            SetScalars(self._panelvisual.vtk_data_array)
        self._panelvisual.object_amp.GetPointData().\
            AddArray(self._panelvisual.vtk_data_array_phase)
        self._panelvisual.object_amp.SetDimensions(shp)
        self._panelvisual.object_amp.Modified()
        self._panelvisual.object_phase.SetPoints(self._panelvisual.object_amp_points)
        self._panelvisual.object_phase.GetPointData(). \
            SetScalars(self._panelvisual.vtk_data_array_phase)
        self._panelvisual.object_phase.SetDimensions(shp)
        self._panelvisual.object_phase.Modified()
        self._panelvisual.image_probe: numpy.ndarray = self._panelvisual.object_phase
        self._panelvisual.lut_phase_real.SetNumberOfTableValues(256)
        self._panelvisual.lut_phase_real.SetScaleToLinear()
        self._panelvisual.lut_phase_real.SetTableRange([phasemin, phasemax])
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
            SetTableRange(self._panelvisual.object_amp.GetPointData().GetScalars().
                          GetRange())
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
        self._panelvisual.scalebar_amp_real.SetLookupTable(self._panelvisual.lut_phase_real)
        # --- Plane drawing starts here
        self._panelvisual.plane.SetNormal(*self._slicing_plane_normal)
        self._panelvisual.plane.SetOrigin(1528, -2450, -4670)
        self._panelvisual.filter_amp_real.SetInputData(self._panelvisual.object_amp)
        self._panelvisual.filter_amp_real.ComputeNormalsOn()
        self._panelvisual.filter_amp_real.ComputeScalarsOn()
        self._panelvisual.filter_amp_real.SetNumberOfContours(1)
        self._panelvisual.filter_amp_real.SetValue(0, contour)
        self._panelvisual.filter_amp_real.Modified()
        self._panelvisual.filter_amp_real.Update()
        self._panelvisual.smooth_filter_real. \
            SetInputConnection(self._panelvisual.filter_amp_real.GetOutputPort())
        self._panelvisual.smooth_filter_real.SetNumberOfIterations(15)
        self._panelvisual.smooth_filter_real.SetRelaxationFactor(0.1)
        self._panelvisual.smooth_filter_real.FeatureEdgeSmoothingOff()
        self._panelvisual.smooth_filter_real.BoundarySmoothingOn()
        self._panelvisual.smooth_filter_real.Update()
        self._panelvisual.normals_amp_real.\
            SetInputConnection(self._panelvisual.smooth_filter_real.GetOutputPort())
        self._panelvisual.normals_amp_real.SetFeatureAngle(feature_angle)
        self._panelvisual.normals_amp_real.ConsistencyOff()
        self._panelvisual.normals_amp_real.SplittingOff()
        self._panelvisual.normals_amp_real.AutoOrientNormalsOff()
        self._panelvisual.normals_amp_real.ComputePointNormalsOn()
        self._panelvisual.normals_amp_real.ComputeCellNormalsOff()
        self._panelvisual.normals_amp_real.NonManifoldTraversalOff()
        self._panelvisual.triangles_amp_real. \
            SetInputConnection(self._panelvisual.normals_amp_real.GetOutputPort())
        self._panelvisual.strips_amp_real. \
            SetInputConnection(self._panelvisual.triangles_amp_real.GetOutputPort())
        self._panelvisual.clipper.SetInputConnection(self._panelvisual.strips_amp_real.
                                               GetOutputPort())
        self._panelvisual.clipper.SetClipFunction(self._panelvisual.plane)
        self._panelvisual.clipper.GenerateClippedOutputOn()
        self._panelvisual.clipper.SetValue(0)
        self._panelvisual.clipper.Update()
        self._panelvisual.mapper_amp_real.SetInputConnection(self._panelvisual.clipper.
                                                       GetOutputPort())
        self._panelvisual.mapper_amp_real.SetLookupTable(self._panelvisual.lut_amp_real)
        self._panelvisual.mapper_amp_real.SetScalarRange(self._panelvisual.object_amp.
                                                   GetPointData().GetScalars().
                                                   GetRange())
        self._panelvisual.mapper_amp_real.SetScalarModeToUsePointData()
        self._panelvisual.mapper_amp_real.Modified()
        self._panelvisual.mapper_amp_real.Update()
        self._panelvisual.actor_amp_real.GetProperty().SetOpacity(1.0)
        self._panelvisual.actor_amp_real.SetMapper(self._panelvisual.mapper_amp_real)
        self._panelvisual.mapper_amp_real2.SetInputConnection(self._panelvisual.clipper.
                                                        GetClippedOutputPort())
        self._panelvisual.mapper_amp_real2.SetLookupTable(self._panelvisual.lut_amp_real)
        self._panelvisual.mapper_amp_real2.SetScalarRange(self._panelvisual.object_amp.
                                                    GetPointData().GetScalars().
                                                    GetRange())
        self._panelvisual.mapper_amp_real2.SetScalarModeToUsePointData()
        self._panelvisual.mapper_amp_real2.Modified()
        self._panelvisual.mapper_amp_real2.Update()
        self._panelvisual.actor_amp_real2.GetProperty().SetOpacity(opacity)
        self._panelvisual.actor_amp_real2.SetMapper(self._panelvisual.mapper_amp_real2)
        self._panelvisual.filter_plane.SetInputData(self._panelvisual.object_amp)
        self._panelvisual.filter_plane.ComputeNormalsOn()
        self._panelvisual.filter_plane.ComputeScalarsOn()
        self._panelvisual.filter_plane.SetNumberOfContours(1)
        self._panelvisual.filter_plane.SetValue(0, contour)
        self._panelvisual.filter_plane.Modified()
        self._panelvisual.filter_plane.Update()
        self._panelvisual.smooth_plane.SetInputConnection(self._panelvisual.filter_plane.
                                                    GetOutputPort())
        self._panelvisual.smooth_plane.SetNumberOfIterations(15)
        self._panelvisual.smooth_plane.SetRelaxationFactor(0.1)
        self._panelvisual.smooth_plane.FeatureEdgeSmoothingOff()
        self._panelvisual.smooth_plane.BoundarySmoothingOn()
        self._panelvisual.smooth_plane.Update()
        self._panelvisual.cutter.SetInputConnection(self._panelvisual.smooth_plane.
                                              GetOutputPort())
        self._panelvisual.cutter.SetCutFunction(self._panelvisual.plane)
        self._panelvisual.cutter.GenerateTrianglesOn()
        self._panelvisual.cutter.Update()
        self._panelvisual.cutter.GenerateCutScalarsOn()
        self._panelvisual.filter_tri.SetInputConnection(self._panelvisual.cutter.
                                                  GetOutputPort())
        objectnpoints = shp[0]*shp[1]*shp[2]
        objectbounds = self._panelvisual.object_amp.GetBounds()
        density = (objectbounds[1] - objectbounds[0]) * \
            (objectbounds[3] - objectbounds[2]) * \
            (objectbounds[5] - objectbounds[4])/objectnpoints
        linedensity = math.pow(density, 1.0/3.0)
        mesh_maxedge = 20.0*linedensity
        self._panelvisual.meshsub.SetInputConnection(self._panelvisual.filter_tri.
                                               GetOutputPort())
        self._panelvisual.meshsub.SetMaximumEdgeLength(mesh_maxedge)
        self._panelvisual.meshsub.SetMaximumNumberOfPasses(meshsubiter)
        self._panelvisual.meshsub.SetOutputPointsPrecision(vtk.vtkAlgorithm.
                                                     SINGLE_PRECISION)
        self._panelvisual.meshsub.Update()
        self._panelvisual.probefilter.SetInputConnection(self._panelvisual.meshsub.
                                                   GetOutputPort())
        self._panelvisual.probefilter.SetSourceData(self._panelvisual.object_phase)
        self._panelvisual.probefilter.Update()
        self._panelvisual.triangles_plane.SetInputConnection(self._panelvisual.probefilter.
                                                       GetOutputPort())
        self._panelvisual.normals_phase_real. \
            SetInputConnection(self._panelvisual.triangles_plane.GetOutputPort())
        self._panelvisual.normals_phase_real.SetFeatureAngle(feature_angle)
        self._panelvisual.normals_phase_real.ConsistencyOff()
        self._panelvisual.normals_phase_real.SplittingOff()
        self._panelvisual.normals_phase_real.AutoOrientNormalsOff()
        self._panelvisual.normals_phase_real.ComputePointNormalsOn()
        self._panelvisual.normals_phase_real.ComputeCellNormalsOn()
        self._panelvisual.normals_phase_real.NonManifoldTraversalOff()
        self._panelvisual.triangles_phase_real. \
            SetInputConnection(self._panelvisual.normals_phase_real.GetOutputPort())
        self._panelvisual.strips_phase_real. \
            SetInputConnection(self._panelvisual.triangles_phase_real.GetOutputPort())
        self._panelvisual.mapper_phase_real. \
            SetInputConnection(self._panelvisual.strips_phase_real.GetOutputPort())
        self._panelvisual.mapper_phase_real. \
            SetLookupTable(self._panelvisual.lut_phase_real)
        self._panelvisual.mapper_phase_real.SetScalarRange([phasemin, phasemax])
        self._panelvisual.mapper_phase_real.SetScalarModeToUsePointData()
        self._panelvisual.mapper_phase_real.Update()
        self._panelvisual.mapper_phase_real.Modified()
        self._panelvisual.actor_phase_real.GetProperty().SetOpacity(1.0)
        self._panelvisual.actor_phase_real.SetMapper(self._panelvisual.mapper_phase_real)
        self._panelvisual.SetPicker()
        self._panelvisual.renderer_amp_real.RemoveAllViewProps()
        self._panelvisual.renderer_phase_real.RemoveAllViewProps()
        self._panelvisual.renderer_amp_recip.RemoveAllViewProps()
        self._panelvisual.renderer_phase_recip.RemoveAllViewProps()
        self._panelvisual.renWin.GetRenderWindow().RemoveRenderer(self._panelvisual.
                                                            renderer_amp_real)
        self._panelvisual.renWin.GetRenderWindow().RemoveRenderer(self._panelvisual.
                                                            renderer_phase_real)
        self._panelvisual.renWin.GetRenderWindow().RemoveRenderer(self._panelvisual.
                                                            renderer_amp_recip)
        self._panelvisual.renWin.GetRenderWindow(). \
            RemoveRenderer(self._panelvisual.renderer_phase_recip)
        self._panelvisual.renWin.GetRenderWindow().Modified()
        self._panelvisual.renderer_amp_real.AddActor(self._panelvisual.actor_amp_real)
        self._panelvisual.renderer_amp_real.AddActor(self._panelvisual.actor_amp_real2)
        self._panelvisual.renderer_amp_real.AddActor(self._panelvisual.actor_phase_real)
        self._panelvisual.renderer_amp_real.AddActor2D(self._panelvisual.scalebar_amp_real)
        self._panelvisual.renderer_amp_real.AddActor2D(self._panelvisual.textActor)
        self._panelvisual.renderer_amp_real.SetBackground(*(self._background))
        self._panelvisual.renWin.GetRenderWindow().AddRenderer(self._panelvisual.
                                                         renderer_phase_real)
        self._panelvisual.renWin.GetRenderWindow().AddRenderer(self._panelvisual.
                                                         renderer_amp_real)
        self._panelvisual.renWin.GetRenderWindow().SetMultiSamples(0)
        self._panelvisual.renWin.SetInteractorStyle(self._panelvisual.style3D)
        self._panelvisual.renderer_amp_real.ResetCamera()
        self._panelvisual.renWin.SetPicker(self._panelvisual.picker_amp_real)
        self._panelvisual.renderer_amp_real.SetViewport(0, 0, 1, 1)
        self._panelvisual.renderer_phase_real.SetViewport(1, 1, 1, 1)
        self.draw_axes_or_not()
        self._panelvisual.Layout()
        self._panelvisual.Show()
        
        return None

    def finalise_rendering(self) -> None:
        """Finalise the rendering and release buttons."""
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
