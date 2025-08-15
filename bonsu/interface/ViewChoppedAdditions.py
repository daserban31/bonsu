"""
Created on Fri Aug 15 10:24:59 2025.

@author: davids_work
"""
################################################################################
#                                 Imports                                      #
################################################################################
import wx
from .common import TextPanelObject, SpinnerObject, MAX_INT, MIN_INT, \
    MAX_INT_16, MIN_INT_16
from .ViewChoppedCommon import rotate_zhat_by_theta_psi
from ..sequences.functions import Sequence_ViewChopped


################################################################################
#                                Constants                                     #
################################################################################


################################################################################
#                                 Methods                                      #
################################################################################


################################################################################
#                                 Classes                                      #
################################################################################
class SubPanel_ViewChopped(wx.ScrolledWindow):
    """
    Similar to View Object, but have the object sliced across a direction.

    Addable to interface/subpanel.py.
    """

    treeitem = {'name':  'View Chopped', 'type': 'operpreview'}

    def sequence(self, selff, pipelineitem):
        """View the sequence - I suppose just in case other code uses this."""
        pass

    def __init__(self, parent, ancestor) -> None:
        """Initialise the Chopper Subpanel."""
        pi: float = 3.141593
        self.ancestor = ancestor
        self.panelvisual = self.ancestor.GetPage(1)
        wx.ScrolledWindow.__init__(self, parent, style=wx.SUNKEN_BORDER)
        vbox: wx.BoxSizer = wx.BoxSizer(wx.VERTICAL)

        # --- Title / Description
        title: wx.StaticText = \
            wx.StaticText(self,
                          label=' '.join(("View Numpy array with coordinate",
                                          "correction chopped and sliced",
                                          "across a direction.")))
        vbox.Add(title, 0, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                 border=2)

        # --- Input File line
        self.input_filename: TextPanelObject = \
            TextPanelObject(self, "Input file: ", "", 100, '*.npy')
        vbox.Add(self.input_filename, 0,
                 flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=2)

        # --- Coords File line
        self.coords_filename: TextPanelObject = \
            TextPanelObject(self, "Co-ord's file: ", "", 100, '*.npy')
        vbox.Add(self.coords_filename, 0,
                 flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=2)

        # --- Amplitude Parameters
        vbox.Add((-1, 10))
        self.sbox1: wx.StaticBox = \
            wx.StaticBox(self, label="Amplitude", style=wx.BORDER_DEFAULT)
        self.sboxs1: wx.StaticBoxSizer = \
            wx.StaticBoxSizer(self.sbox1, wx.VERTICAL)
        self.hbox1: wx.BoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.contour: SpinnerObject = \
            SpinnerObject(self, "Isosurf.: ", MAX_INT, MIN_INT, 1, 100, 75, 80)
        self.contour.SetToolTip("Isosurface value")
        self.hbox1.Add(self.contour, 0,
                       flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=2)
        self.hbox1.Add((5, -1))
        self.opacity: SpinnerObject = \
            SpinnerObject(self, "Opacity: ", 1.0, 0.0, 0.1, 0.5, 75, 80)
        self.opacity.SetToolTip("Opacity of isosurface")
        self.hbox1.Add(self.opacity, 0,
                       flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=2)
        self.sboxs1.Add(self.hbox1, 0,
                        flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=2)
        self.feature_angle: SpinnerObject = \
            SpinnerObject(self, "FA: ", 180, 0, 1, 90, 75, 40)
        self.feature_angle.SetToolTip("Feature Angle of rendered objects.")
        self.sboxs1.Add(self.feature_angle, 0,
                        flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP |
                        wx.BOTTOM, border=2)
        vbox.Add(self.sboxs1, 0, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                 border=2)
        vbox.Add((-1, 5))

        # --- Phase Parameters
        self.sbox2: wx.StaticBox = \
            wx.StaticBox(self, label="Phase", style=wx.BORDER_DEFAULT)
        self.sboxs2: wx.StaticBoxSizer = \
            wx.StaticBoxSizer(self.sbox2, wx.VERTICAL)
        self.hbox2: wx.BoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.phasemax: SpinnerObject = \
            SpinnerObject(self, "Max: ", pi, 0.0, 0.01, pi, 80, 80)
        self.hbox2.Add(self.phasemax, 0,
                       flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM,
                       border=2)
        self.hbox2.Add((5, -1))
        self.phasemin: SpinnerObject = \
            SpinnerObject(self, "Min: ", 0.0, -pi, 0.01, -pi, 80, 80)
        self.hbox2.Add(self.phasemin, 0,
                       flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM,
                       border=2)
        self.sboxs2.Add(self.hbox2, 0,
                        flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=2)
        vbox.Add(self.sboxs2, 0,
                 flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=2)

        # --- Cut plane params
        vbox.Add((-1, 5))
        self.sbox4: wx.StaticBox = \
            wx.StaticBox(self, label="Slicing Parameters ",
                         style=wx.BORDER_DEFAULT)
        self.sboxs4: wx.StaticBoxSizer = \
            wx.StaticBoxSizer(self.sbox4, wx.VERTICAL)
        self.hbox4 = wx.BoxSizer(wx.HORIZONTAL)
        self.number_of_slices: SpinnerObject = \
            SpinnerObject(self, "No. of Slices: ", MAX_INT_16, MIN_INT_16, 1,
                          100, 15, 60)
        self.theta: SpinnerObject = \
            SpinnerObject(self, "Polar \u0398 (deg): ", MAX_INT_16,
                          MIN_INT_16, 1, 0, 15, 60)
        self.psi: SpinnerObject = \
            SpinnerObject(self, "Azimuthal \u03A8 (deg): ", MAX_INT_16,
                          MIN_INT_16, 1, 0, 15, 60)
        self.separation: SpinnerObject = \
            SpinnerObject(self, "Separation (nm): ", MAX_INT_16, MIN_INT_16,
                          50, 100, 15, 60)
        self.thickness: SpinnerObject = \
            SpinnerObject(self, 'Thickness (nm)', MAX_INT_16, MIN_INT_16, 50,
                          200, 15, 60)
        self.hbox4.Add(self.number_of_slices, 0, flag=wx.EXPAND | wx.RIGHT,
                       border=5)
        self.hbox4.Add(self.theta, 0, flag=wx.EXPAND | wx.LEFT | wx.RIGHT,
                       border=5)
        self.hbox4.Add(self.psi, 0, flag=wx.EXPAND | wx.LEFT | wx.RIGHT,
                       border=5)
        self.hbox4.Add(self.separation, 0, flag=wx.EXPAND | wx.LEFT | wx.Right,
                       border=5)
        self.hbox4.Add(self.thickness, 0, flag=wx.EXPAND | wx.LEFT | wx.Right,
                       border=5)
        self.sboxs4.Add(self.hbox4, 0,
                        flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP |
                        wx.BOTTOM,
                        border=2)
        self.sboxs4.Add((-1, 5))
        vbox.Add((-1, 5))
        vbox.Add(self.sboxs4, 0, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                 border=2)
        self.number_of_slices.spin.SetEventFunc(self.OnPlaneSpin)
        self.separation.spin.SetEventFunc(self.OnPlaneSpin)
        self.thickness.spin.SetEventFunc(self.OnPlaneSpin)
        self.theta.spin.SetEventFunc(self.OnPlaneSpin)
        self.psi.spin.SetEventFunc(self.OnPlaneSpin)
        self.number_of_slices.value.Bind(wx.EVT_KEY_DOWN, self.OnPlaneKey)
        self.separation.value.Bind(wx.EVT_KEY_DOWN, self.OnPlaneKey)
        self.thickness.value.Bind(wx.EVT_KEY_DOWN, self.OnPlaneKey)
        self.theta.value.Bind(wx.EVT_KEY_DOWN, self.OnPlaneKey)
        self.psi.value.Bind(wx.EVT_KEY_DOWN, self.OnPlaneKey)

        # --- Clipped Mesh Iterations Params
        vbox.Add((-1, 5))
        self.hbox5: wx.BoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.meshsubiter: SpinnerObject = \
            SpinnerObject(self, "Clipped mesh iterations: ", MAX_INT_16, 1, 1,
                          5, 120, 120)
        self.hbox5.Add(self.meshsubiter, 0, flag=wx.EXPAND | wx.LEFT, border=10)
        vbox.Add(self.hbox5, 0,
                 flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM,
                 border=2)

        # --- View Axes Tick Box
        vbox.Add((-1, 5))
        self.hbox6: wx.BoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.chkbox_axes: wx.CheckBox = \
            wx.CheckBox(self, -1, 'View axes', size=(200, 20))
        self.chkbox_axes.SetValue(False)
        self.hbox6.Add(self.chkbox_axes, 1,
                       flag=wx.EXPAND | wx.LEFT | wx.RIGHT, border=2)
        self.hbox6.Add((-1, 5))
        self.axes_fontfactor: SpinnerObject = \
            SpinnerObject(self, "Font Factor:", MAX_INT, 1, 1, 2, 100, 100)
        self.hbox6.Add(self.axes_fontfactor, 0,
                       flag=wx.EXPAND | wx.LEFT | wx.RIGHT, border=2)
        vbox.Add(self.hbox6, 0,
                 flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM,
                 border=2)

        # --- View Button
        vbox.Add((-1, 5))
        button_view = wx.Button(self, label="View", size=(70, 30))
        button_view.Bind(wx.EVT_BUTTON, self.SeqParser)
        vbox.Add(button_view, 0, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                 border=2)

        # Finishing the whole panel layout
        self.SetAutoLayout(True)
        self.SetSizer(vbox)
        self.FitInside()
        self.SetScrollRate(5, 5)
        return None

    def SeqParser(self, event) -> None:
        """Parse the Chopped Slices."""
        Sequence_ViewChopped(self, self.ancestor)
        self.ancestor.GetPage(4).data_poll_timer.Start(1000)
        return None

    def OnPlaneSpin(self, event) -> None:
        """Update the visualisation when params change."""
        number_of_slices: float = float(self.number_of_slices.value.GetValue())
        thickness: float = float(self.thickness.value.GetValue())
        separation: float = float(self.separation.value.GetValue())
        theta: float = float(self.theta.value.GetValue())
        psi: float = float(self.psi.value.GetValue())
        self.panelvisual.plane.SetOrigin(1528, -2450, -4670)
        self.panelvisual.plane.SetNormal(*rotate_zhat_by_theta_psi(theta, psi))
        self.panelvisual.plane.Modified()
        self.panelvisual.RefreshScene()
        return None

    def OnPlaneKey(self, event) -> None:
        """Update visualisation when change in panel."""
        if event.GetKeyCode() == wx.WXK_RETURN:
            self.OnPlaneSpin(None)
        else:
            event.Skip()
        return None


################################################################################
#                         Dunder Name Dunder Main                              #
################################################################################
if __name__ == '__main__':
    pass
