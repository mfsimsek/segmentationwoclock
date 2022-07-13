# segmentationwoclock
Simseketal_2022_Manuscript

This code (clockSimulation.m) incorporates segmentation clock as an inhibitor into FGF signaling pathway simulation. No clock, no diffusion and pulsatile inhibitor drug treatment cases can be simulated by changing relevant parameters. Clock is defined according to its experimentally determined period change from posterior to anterior PSM, with a simplified sine function. Amplitude of the clock is constant along the PSM.

Createfigure.m converts mesh kymographs into publication style 2-D plots.

borderliner.m finds the border trend of kymographs to interpret discrete dynamics of wavefront.

FIJI macros beginning with IHC are for analyzing immunostaining data (generic pipeline code can be used for any number of staining channels). Analysis gives intensity quantifications along an LOI of choice as nuclear masked, non-nuclear (cytoplasmic) masked and unmasked data of z-layers or maximum intensity projection of 3 consecutive z-layers.

Xirp macros in two parts (before and after SNT algorithm) are for quantifying somite boundary occupancy/intactness from xirp2a ISH staining. 
