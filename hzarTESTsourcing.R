## Files involved:

## Core code files
## hzarClasses.R
## +-> modelEquationPrototypeGrouping.R


## ignoring
## plotClineFunc

#### prep work functions
source("hzarPrep.R");
## setup transect with molecular data.
## sampleLikelihoodMolecularPop
## doMolecularData1DPops

## helper functions for model setup
## mkParam
## setupMoleCenterClineParameters
## buildCline1D

## model templates

## model setup
## makeSimpleCline1D
## +->buildCline1D
## makeTailedCline1D
## +->buildCline1D
## makeCline1D
## +->buildCline1D

#### fitting functions and their helpers

## helper functions for fitCline
## getCredibleLLspace
## getCredibleCutG
## getCredibleCut
## cline.unscale
## cline.logit
## splitParameters
## cssp
## cline.samp.dist

## model fitting
## fitClineModel
## +->cssp
## +->splitParameters
## reFitClineFunc
## +->splitParameters


#### Stuff I am ignoring for now.

## Extensions
## hzarPlotting.R
## hzarMultiModel.R

## Test data
## molecular-data-robb2.csv
## RailGenData.csv

## Logging / testing 
## oldTests/
