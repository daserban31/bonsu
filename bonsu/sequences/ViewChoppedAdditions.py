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

    def ViewDataAmpClippedPhase(self, ancestor, data, r, g, b):
        self.ancestor.GetPage(0).queue_info. \
            put("Preparing object visualisation")
        self.ancestor.GetPage(0).queue_info.put(''.join(("Array Size: ",
                                                         str(data.shape))))
        panelvisual = ancestor.GetPage(1)
        number_of_slices: int = int(self.number_of_slices.value.GetValue())
        separation: float = float(self.separation.value.GetValue())
        theta: float = float(self.theta.value.GetValue())
        psi = float(self.psi.value.GetValue())
        # --- Calculate the slicing coordinates
        slicing_plane_normal: numpy.ndarray = \
            rotate_zhat_by_theta_psi(theta_deg=theta, psi_deg=psi)
        projections_along_normal = numpy.matmul(panelvisual.coords,
                                                slicing_plane_normal)
        line_min = projections_along_normal.min()
        line_max = projections_along_normal.max()
        bin_edges = numpy.linspace(line_min, line_max, number_of_slices)
        coordindxs = numpy.digitize(projections_along_normal, bin_edges,
                                    right=True)
        # panelvisual.coords += separation * \
        #     numpy.multiply(coordindxs.reshape((coordindxs.size, 1)),
        #                    slicing_plane_normal)
        # --- Continue the usual script
        meshsubiter: int = int(float(self.meshsubiter.value.GetValue()))
        opacity: float = float(self.opacity.value.GetValue())
        contour: float = float(self.contour.value.GetValue())
        maxval: numpy.float32 = numpy.abs(data).max()
        if contour > maxval:
            contour = CNTR_CLIP*maxval
        feature_angle: float = float(self.feature_angle.value.GetValue())
        phasemax: float = float(self.phasemax.value.GetValue())
        phasemin: float = float(self.phasemin.value.GetValue())
        panelvisual.flat_data_phase: numpy.ndarray = \
            (numpy.angle(data)).transpose(2, 1, 0).flatten()
        panelvisual.vtk_data_array_phase: numpy.ndarray = \
            numpy_support.numpy_to_vtk(panelvisual.flat_data_phase)
        panelvisual.vtk_data_array_phase.SetName("mapscalar")
        shp: numpy.ndarray = numpy.array(data.shape, dtype=numpy.int64)
        panelvisual.flat_data: numpy.ndarray = \
            (numpy.abs(data)).transpose(2, 1, 0).flatten()
        panelvisual.vtk_data_array: numpy.ndarray = \
            numpy_support.numpy_to_vtk(panelvisual.flat_data)
        panelvisual.vtk_coordarray: numpy.ndarray = \
            numpy_support.numpy_to_vtk(panelvisual.coords)
        panelvisual.object_amp_points.SetDataTypeToDouble()
        panelvisual.object_amp_points.SetNumberOfPoints(data.size)
        panelvisual.object_amp_points.SetData(panelvisual.vtk_coordarray)
        panelvisual.object_amp.SetPoints(panelvisual.object_amp_points)
        panelvisual.object_amp.GetPointData(). \
            SetScalars(panelvisual.vtk_data_array)
        panelvisual.object_amp.GetPointData().\
            AddArray(panelvisual.vtk_data_array_phase)
        panelvisual.object_amp.SetDimensions(shp)
        panelvisual.object_amp.Modified()
        panelvisual.object_phase.SetPoints(panelvisual.object_amp_points)
        panelvisual.object_phase.GetPointData(). \
            SetScalars(panelvisual.vtk_data_array_phase)
        panelvisual.object_phase.SetDimensions(shp)
        panelvisual.object_phase.Modified()
        panelvisual.image_probe: numpy.ndarray = panelvisual.object_phase
        panelvisual.lut_phase_real.SetNumberOfTableValues(256)
        panelvisual.lut_phase_real.SetScaleToLinear()
        panelvisual.lut_phase_real.SetTableRange([phasemin, phasemax])
        lutsource: numpy.ndarray = \
            self.ancestor.GetPage(0).cms[self.ancestor.GetPage(0).cmls[1][0]][1]
        if self.ancestor.GetPage(0).cmls[1][1] == 0:
            for k in range(256):
                panelvisual.lut_phase_real. \
                    SetTableValue(k, lutsource[k][0], lutsource[k][1],
                                  lutsource[k][2], 1)
        else:
            for k in range(256):
                panelvisual.lut_phase_real. \
                    SetTableValue(255-k, lutsource[k][0], lutsource[k][1],
                                  lutsource[k][2], 1)
        panelvisual.lut_phase_real.SetRamp(0)
        panelvisual.lut_phase_real.Build()
        panelvisual.lut_amp_real.SetNumberOfTableValues(256)
        panelvisual.lut_amp_real.SetScaleToLinear()
        panelvisual.lut_amp_real. \
            SetTableRange(panelvisual.object_amp.GetPointData().GetScalars().
                          GetRange())
        lutsource: numpy.ndarray = \
            self.ancestor.GetPage(0).cms[self.ancestor.GetPage(0).cmls[0][0]][1]
        if self.ancestor.GetPage(0).cmls[0][1] == 0:
            for k in range(256):
                panelvisual.lut_amp_real. \
                    SetTableValue(k, lutsource[k][0], lutsource[k][1],
                                  lutsource[k][2], 1)
        else:
            for k in range(256):
                panelvisual.lut_amp_real. \
                    SetTableValue(255-k, lutsource[k][0], lutsource[k][1],
                                  lutsource[k][2], 1)
        panelvisual.lut_amp_real.SetRamp(0)
        panelvisual.lut_amp_real.Build()
        panelvisual.scalebar_amp_real.SetTitle("")
        panelvisual.scalebar_amp_real.SetLookupTable(panelvisual.lut_phase_real)
        # --- Plane drawing starts here
        panelvisual.plane.SetNormal(*slicing_plane_normal)
        panelvisual.plane.SetOrigin(1528, -2450, -4670)
        panelvisual.filter_amp_real.SetInputData(panelvisual.object_amp)
        panelvisual.filter_amp_real.ComputeNormalsOn()
        panelvisual.filter_amp_real.ComputeScalarsOn()
        panelvisual.filter_amp_real.SetNumberOfContours(1)
        panelvisual.filter_amp_real.SetValue(0, contour)
        panelvisual.filter_amp_real.Modified()
        panelvisual.filter_amp_real.Update()
        panelvisual.smooth_filter_real. \
            SetInputConnection(panelvisual.filter_amp_real.GetOutputPort())
        panelvisual.smooth_filter_real.SetNumberOfIterations(15)
        panelvisual.smooth_filter_real.SetRelaxationFactor(0.1)
        panelvisual.smooth_filter_real.FeatureEdgeSmoothingOff()
        panelvisual.smooth_filter_real.BoundarySmoothingOn()
        panelvisual.smooth_filter_real.Update()
        panelvisual.normals_amp_real.\
            SetInputConnection(panelvisual.smooth_filter_real.GetOutputPort())
        panelvisual.normals_amp_real.SetFeatureAngle(feature_angle)
        panelvisual.normals_amp_real.ConsistencyOff()
        panelvisual.normals_amp_real.SplittingOff()
        panelvisual.normals_amp_real.AutoOrientNormalsOff()
        panelvisual.normals_amp_real.ComputePointNormalsOn()
        panelvisual.normals_amp_real.ComputeCellNormalsOff()
        panelvisual.normals_amp_real.NonManifoldTraversalOff()
        panelvisual.triangles_amp_real. \
            SetInputConnection(panelvisual.normals_amp_real.GetOutputPort())
        panelvisual.strips_amp_real. \
            SetInputConnection(panelvisual.triangles_amp_real.GetOutputPort())
        panelvisual.clipper.SetInputConnection(panelvisual.strips_amp_real.
                                               GetOutputPort())
        panelvisual.clipper.SetClipFunction(panelvisual.plane)
        panelvisual.clipper.GenerateClippedOutputOn()
        panelvisual.clipper.SetValue(0)
        panelvisual.clipper.Update()
        panelvisual.mapper_amp_real.SetInputConnection(panelvisual.clipper.
                                                       GetOutputPort())
        panelvisual.mapper_amp_real.SetLookupTable(panelvisual.lut_amp_real)
        panelvisual.mapper_amp_real.SetScalarRange(panelvisual.object_amp.
                                                   GetPointData().GetScalars().
                                                   GetRange())
        panelvisual.mapper_amp_real.SetScalarModeToUsePointData()
        panelvisual.mapper_amp_real.Modified()
        panelvisual.mapper_amp_real.Update()
        panelvisual.actor_amp_real.GetProperty().SetOpacity(1.0)
        panelvisual.actor_amp_real.SetMapper(panelvisual.mapper_amp_real)
        panelvisual.mapper_amp_real2.SetInputConnection(panelvisual.clipper.
                                                        GetClippedOutputPort())
        panelvisual.mapper_amp_real2.SetLookupTable(panelvisual.lut_amp_real)
        panelvisual.mapper_amp_real2.SetScalarRange(panelvisual.object_amp.
                                                    GetPointData().GetScalars().
                                                    GetRange())
        panelvisual.mapper_amp_real2.SetScalarModeToUsePointData()
        panelvisual.mapper_amp_real2.Modified()
        panelvisual.mapper_amp_real2.Update()
        panelvisual.actor_amp_real2.GetProperty().SetOpacity(opacity)
        panelvisual.actor_amp_real2.SetMapper(panelvisual.mapper_amp_real2)
        panelvisual.filter_plane.SetInputData(panelvisual.object_amp)
        panelvisual.filter_plane.ComputeNormalsOn()
        panelvisual.filter_plane.ComputeScalarsOn()
        panelvisual.filter_plane.SetNumberOfContours(1)
        panelvisual.filter_plane.SetValue(0, contour)
        panelvisual.filter_plane.Modified()
        panelvisual.filter_plane.Update()
        panelvisual.smooth_plane.SetInputConnection(panelvisual.filter_plane.
                                                    GetOutputPort())
        panelvisual.smooth_plane.SetNumberOfIterations(15)
        panelvisual.smooth_plane.SetRelaxationFactor(0.1)
        panelvisual.smooth_plane.FeatureEdgeSmoothingOff()
        panelvisual.smooth_plane.BoundarySmoothingOn()
        panelvisual.smooth_plane.Update()
        panelvisual.cutter.SetInputConnection(panelvisual.smooth_plane.
                                              GetOutputPort())
        panelvisual.cutter.SetCutFunction(panelvisual.plane)
        panelvisual.cutter.GenerateTrianglesOn()
        panelvisual.cutter.Update()
        panelvisual.cutter.GenerateCutScalarsOn()
        panelvisual.filter_tri.SetInputConnection(panelvisual.cutter.
                                                  GetOutputPort())
        objectnpoints = shp[0]*shp[1]*shp[2]
        objectbounds = panelvisual.object_amp.GetBounds()
        density = (objectbounds[1] - objectbounds[0]) * \
            (objectbounds[3] - objectbounds[2]) * \
            (objectbounds[5] - objectbounds[4])/objectnpoints
        linedensity = math.pow(density, 1.0/3.0)
        mesh_maxedge = 20.0*linedensity
        panelvisual.meshsub.SetInputConnection(panelvisual.filter_tri.
                                               GetOutputPort())
        panelvisual.meshsub.SetMaximumEdgeLength(mesh_maxedge)
        panelvisual.meshsub.SetMaximumNumberOfPasses(meshsubiter)
        panelvisual.meshsub.SetOutputPointsPrecision(vtk.vtkAlgorithm.
                                                     SINGLE_PRECISION)
        panelvisual.meshsub.Update()
        panelvisual.probefilter.SetInputConnection(panelvisual.meshsub.
                                                   GetOutputPort())
        panelvisual.probefilter.SetSourceData(panelvisual.object_phase)
        panelvisual.probefilter.Update()
        panelvisual.triangles_plane.SetInputConnection(panelvisual.probefilter.
                                                       GetOutputPort())
        panelvisual.normals_phase_real. \
            SetInputConnection(panelvisual.triangles_plane.GetOutputPort())
        panelvisual.normals_phase_real.SetFeatureAngle(feature_angle)
        panelvisual.normals_phase_real.ConsistencyOff()
        panelvisual.normals_phase_real.SplittingOff()
        panelvisual.normals_phase_real.AutoOrientNormalsOff()
        panelvisual.normals_phase_real.ComputePointNormalsOn()
        panelvisual.normals_phase_real.ComputeCellNormalsOn()
        panelvisual.normals_phase_real.NonManifoldTraversalOff()
        panelvisual.triangles_phase_real. \
            SetInputConnection(panelvisual.normals_phase_real.GetOutputPort())
        panelvisual.strips_phase_real. \
            SetInputConnection(panelvisual.triangles_phase_real.GetOutputPort())
        panelvisual.mapper_phase_real. \
            SetInputConnection(panelvisual.strips_phase_real.GetOutputPort())
        panelvisual.mapper_phase_real. \
            SetLookupTable(panelvisual.lut_phase_real)
        panelvisual.mapper_phase_real.SetScalarRange([phasemin, phasemax])
        panelvisual.mapper_phase_real.SetScalarModeToUsePointData()
        panelvisual.mapper_phase_real.Update()
        panelvisual.mapper_phase_real.Modified()
        panelvisual.actor_phase_real.GetProperty().SetOpacity(1.0)
        panelvisual.actor_phase_real.SetMapper(panelvisual.mapper_phase_real)
        panelvisual.SetPicker()
        panelvisual.renderer_amp_real.RemoveAllViewProps()
        panelvisual.renderer_phase_real.RemoveAllViewProps()
        panelvisual.renderer_amp_recip.RemoveAllViewProps()
        panelvisual.renderer_phase_recip.RemoveAllViewProps()
        panelvisual.renWin.GetRenderWindow().RemoveRenderer(panelvisual.
                                                            renderer_amp_real)
        panelvisual.renWin.GetRenderWindow().RemoveRenderer(panelvisual.
                                                            renderer_phase_real)
        panelvisual.renWin.GetRenderWindow().RemoveRenderer(panelvisual.
                                                            renderer_amp_recip)
        panelvisual.renWin.GetRenderWindow(). \
            RemoveRenderer(panelvisual.renderer_phase_recip)
        panelvisual.renWin.GetRenderWindow().Modified()
        panelvisual.renderer_amp_real.AddActor(panelvisual.actor_amp_real)
        panelvisual.renderer_amp_real.AddActor(panelvisual.actor_amp_real2)
        panelvisual.renderer_amp_real.AddActor(panelvisual.actor_phase_real)
        panelvisual.renderer_amp_real.AddActor2D(panelvisual.scalebar_amp_real)
        panelvisual.renderer_amp_real.AddActor2D(panelvisual.textActor)
        panelvisual.renderer_amp_real.SetBackground(r, g, b)
        panelvisual.renWin.GetRenderWindow().AddRenderer(panelvisual.
                                                         renderer_phase_real)
        panelvisual.renWin.GetRenderWindow().AddRenderer(panelvisual.
                                                         renderer_amp_real)
        panelvisual.renWin.GetRenderWindow().SetMultiSamples(0)
        panelvisual.renWin.SetInteractorStyle(panelvisual.style3D)
        panelvisual.renderer_amp_real.ResetCamera()
        panelvisual.renWin.SetPicker(panelvisual.picker_amp_real)
        panelvisual.renderer_amp_real.SetViewport(0, 0, 1, 1)
        panelvisual.renderer_phase_real.SetViewport(1, 1, 1, 1)
        panelvisual.Layout()
        panelvisual.Show()
        if self.chkbox_axes.GetValue():
            panelvisual.axis.SetBounds(panelvisual.object_amp.GetBounds())
            panelvisual.axis.SetInputData(panelvisual.object_amp)
            panelvisual.axis.SetCamera(panelvisual.renderer_amp_real.
                                       GetActiveCamera())
            panelvisual.axis.SetLabelFormat("%6.1f")
            panelvisual.axis.SetFlyModeToOuterEdges()
            panelvisual.axis.ScalingOff()
            panelvisual.axis.SetFontFactor(float(self.axes_fontfactor.value.
                                                 GetValue()))
            panelvisual.axis.SetXLabel("X")
            panelvisual.axis.SetYLabel("Y")
            panelvisual.axis.SetZLabel("Z")
            panelvisual.axis.Modified()
            panelvisual.renderer_amp_real.AddViewProp(panelvisual.axis)
    data_file = self.input_filename.objectpath.GetValue()
    coords_file = self.coords_filename.objectpath.GetValue()
    panelvisual = self.ancestor.GetPage(1)
    r = float(panelvisual.r)/255.0
    g = float(panelvisual.g)/255.0
    b = float(panelvisual.b)/255.0
    try:
        panelvisual.data = LoadArray(self.ancestor.GetPage(0), data_file)
        panelvisual.data_max = numpy.abs(panelvisual.data).max()
        panelvisual.coords = LoadCoordsArray(self.ancestor.GetPage(0),
                                             coords_file, panelvisual.data)
    except Exception as e:
        print(f'View Chopped slices errored {type(e)} - {e}.')
        msg = "Could not load array."
        dlg = wx.MessageDialog(self, msg, "Sequence View Object", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
        return
    else:
        ViewDataAmpClippedPhase(self, ancestor, panelvisual.data, r, g, b)
        panelvisual.datarangelist = []
        panelvisual.ReleaseVisualButtons(gotovisual=True)
        panelvisual.button_vremove.Enable(False)
        panelvisual.RefreshSceneFull()
    return None


################################################################################
#                                 Classes                                      #
################################################################################


################################################################################
#                         Dunder Name Dunder Main                              #
################################################################################
if __name__ == '__main__':
    pass
