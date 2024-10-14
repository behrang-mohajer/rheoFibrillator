# state file generated using paraview version 5.6.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1357, 795]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [20.0, 0.0, 0.0]
renderView1.StereoType = 0
renderView1.CameraPosition = [14.895099214394772, 28.815515022552244, 71.72686930659135]
renderView1.CameraFocalPoint = [19.999999999999993, 2.1413501698662904e-15, 2.5124606574725998e-15]
renderView1.CameraViewUp = [0.21575678375030302, 0.9112791507608282, -0.3507411006060579]
renderView1.CameraParallelScale = 20.049937655763422
renderView1.Background = [0.32, 0.34, 0.43]

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVFoamReader'
channelOpenFOAM = PVFoamReader(FileName='/home/behrang/OpenFOAM/behrang-9/run/channel/channel.OpenFOAM')
channelOpenFOAM.MeshParts = ['internalMesh', 'inlet - patch', 'outlet - patch', 'walls - wall']
channelOpenFOAM.Fields = ['U', 'tau']

# create a new 'Calculator'
n_11 = Calculator(Input=channelOpenFOAM)
n_11.ResultArrayName = 'N_11'
n_11.Function = 'tau_XX - tau_YY'

# create a new 'Calculator'
shear = Calculator(Input=n_11)
shear.ResultArrayName = 'Shear'
shear.Function = 'tau_XY'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from shear
shearDisplay = Show(shear, renderView1)

# get color transfer function/color map for 'Shear'
shearLUT = GetColorTransferFunction('Shear')
shearLUT.RGBPoints = [-3.263613224029541, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 3.263613224029541, 0.705882, 0.0156863, 0.14902]
shearLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
shearDisplay.Representation = 'Surface'
shearDisplay.ColorArrayName = ['POINTS', 'Shear']
shearDisplay.LookupTable = shearLUT
shearDisplay.OSPRayScaleArray = 'Shear'
shearDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
shearDisplay.SelectOrientationVectors = 'N_11'
shearDisplay.ScaleFactor = 4.0
shearDisplay.SelectScaleArray = 'Shear'
shearDisplay.GlyphType = 'Arrow'
shearDisplay.GlyphTableIndexArray = 'Shear'
shearDisplay.GaussianRadius = 0.2
shearDisplay.SetScaleArray = ['POINTS', 'Shear']
shearDisplay.ScaleTransferFunction = 'PiecewiseFunction'
shearDisplay.OpacityArray = ['POINTS', 'Shear']
shearDisplay.OpacityTransferFunction = 'PiecewiseFunction'
shearDisplay.DataAxesGrid = 'GridAxesRepresentation'
shearDisplay.SelectionCellLabelFontFile = ''
shearDisplay.SelectionPointLabelFontFile = ''
shearDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
shearDisplay.DataAxesGrid.XTitleFontFile = ''
shearDisplay.DataAxesGrid.YTitleFontFile = ''
shearDisplay.DataAxesGrid.ZTitleFontFile = ''
shearDisplay.DataAxesGrid.XLabelFontFile = ''
shearDisplay.DataAxesGrid.YLabelFontFile = ''
shearDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
shearDisplay.PolarAxes.PolarAxisTitleFontFile = ''
shearDisplay.PolarAxes.PolarAxisLabelFontFile = ''
shearDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
shearDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for shearLUT in view renderView1
shearLUTColorBar = GetScalarBar(shearLUT, renderView1)
shearLUTColorBar.Title = 'Shear'
shearLUTColorBar.ComponentTitle = ''
shearLUTColorBar.TitleFontFile = ''
shearLUTColorBar.LabelFontFile = ''

# set color bar visibility
shearLUTColorBar.Visibility = 1

# show color legend
shearDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Shear'
shearPWF = GetOpacityTransferFunction('Shear')
shearPWF.Points = [-3.263613224029541, 0.0, 0.5, 0.0, 3.263613224029541, 1.0, 0.5, 0.0]
shearPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(shear)
# ----------------------------------------------------------------