\name{simulateSpectra}

\alias{simulateSpectra}
\alias{simulateImage}

\alias{addShape}
\alias{presetImageDef}

\title{Simulate a mass spectrum or MS imaging experiment}

\description{
	Simulate mass spectra or complete MS imaging experiments, including a possible baseline, spatial and spectral noise, mass drift, mass resolution, and multiplicative variation, etc.

    A number of preset imaging designs are available for quick-and-dirty simulation of images.

    These functions are designed for small proof-of-concept examples and testing, and may not scale well to simulating larger datasets.
}

\usage{
simulateSpectra(n = 1L, npeaks = 50L,
    mz = rlnorm(npeaks, 7, 0.3), intensity = rlnorm(npeaks, 1, 0.9),
    from = 0.9 * min(mz), to = 1.1 * max(mz), by = 400,
    sdpeaks = sdpeakmult * log1p(intensity), sdpeakmult = 0.2,
    sdnoise = 0.1, sdmz = 10, resolution = 1000, fmax = 0.5,
    baseline = 0, decay = 10, units=c("ppm", "mz"),
    centroided = FALSE, \dots)

simulateImage(pixelData, featureData, preset,
    from = 0.9 * min(mz), to = 1.1 * max(mz), by = 400,
    sdrun = 1, sdpixel = 1, spcorr = 0.3, SAR = TRUE,
    resolution = 1000, fmax = 0.5, units=c("ppm", "mz"),
    centroided = FALSE, continuous = TRUE,
    verbose = getCardinalVerbose(), chunkopts = list(),
    BPPARAM = getCardinalBPPARAM(), \dots)

addShape(pixelData, center, size, shape=c("circle", "square"), name=shape)

presetImageDef(preset = 1L, nrun = 1, npeaks = 30L,
    dim = c(20L, 20L), peakheight = exp(1), peakdiff = exp(1),
    sdsample = 0.2, jitter = TRUE, \dots)
}

\arguments{
    \item{n}{The number of spectra to simulate.}

    \item{npeaks}{The number of peaks to simulate. Not used if \code{mz} and \code{intensity} are provided.}

    \item{mz}{The theoretical m/z values of the simulated peaks.}

    \item{intensity}{The mean intensities of the simulated peaks.}

    \item{from}{The minimum m/z value used for the mass range.}

    \item{to}{The maximum m/z value used for the mass range.}

    \item{by}{The step-size used for the observed m/z-values of the profile spectrum.}

    \item{sdpeaks}{The standard deviation(s) for the distributions of observed peak intensities on the log scale.}

    \item{sdpeakmult}{A multiplier used to calculate \code{sdpeaks} based on the mean intensities of peaks; used to simulate multiplicative variance. Not used if \code{sdpeaks} is provided.}

    \item{sdnoise}{The standard deviation of the random noise in the spectrum on the log scale.}

    \item{sdmz}{The standard deviation of the mass error in the observed m/z values of peaks, in units indicated by \code{units}.}

    \item{resolution}{The mass resolution as defined by \code{m / dm}, where \code{m} is the observed mass and \code{dm} is the width of the peak at a proportion of its maximum height defined by \code{fmax} (defaults to full-width-at-half-maximum -- FWHM -- definition). Note that this is NOT the same as the definition of \code{resolution} in the \code{\link{readImzML}} function.}

    \item{fmax}{The fraction of the maximum peak height to use when defining the mass resolution.}

    \item{baseline}{The maximum intensity of the baseline. Note that \code{baseline=0} means there is no baseline.}

    \item{decay}{A constant used to calculate the exponential decay of the baseline. Larger values mean the baseline decays more sharply at the lower mass range of the spectrum.}

    \item{units}{The units for \code{by} and \code{sdmz}. Either parts-per-million or absolute m/z units.}

    \item{centroided}{Should the simulated spectrum representation be centroided (\code{TRUE}) or profile (\code{FALSE})?}

    \item{continuous}{Should the simulated spectrum storage type be continuous (\code{TRUE}) or processed (\code{FALSE}), where "continuous" means a dense representation and "processed" means a sparse representation?}

    \item{verbose}{Should progress messages be printed?}

    \item{chunkopts}{Chunk processing options. See \code{\link{chunkApply}} for details.}

    \item{BPPARAM}{An optional instance of \code{BiocParallelParam}. See documentation for \code{\link{bplapply}}.}

    \item{pixelData}{A \code{\linkS4class{PositionDataFrame}} giving the pixel design of the experiment. The names of the columns should match the names of columns in \code{featureData}. Each column should be a logical vector corresponding to a morphological substructure, indicate which pixels belong to that substructure.}

    \item{featureData}{A \code{\linkS4class{MassDataFrame}} giving the feature design of the experiment. Each row should correspond to an expected peak. The names of the columns should match the names of columns in \code{pixelData}. Each column should be a numeric vector corresponding to a morphological substructure, giving the mean intensity of that peak for that substructure.}

    \item{preset}{A number indicating a preset image definition to use.}

    \item{nrun}{The number of runs to simulate for each condition.}

    \item{sdrun}{A standard deviation giving the run-to-run variance.}

    \item{sdpixel}{A standard deviation giving the pixel-to-pixel variance.}

    \item{spcorr}{The spatial autocorrelation. Must be between 0 and 1, where \code{spcorr=0} indicates no spatial autocorrelation.}

    \item{SAR}{Should a spatial autoregressive (SAR) model be used for simulating spatially-correlated noise (\code{TRUE}) versus a simpler model that uses spatially-smoothed Gaussian noise (\code{FALSE})? The calculation of the SAR matrix for large images can be very time-consuming, so if the simpler model is adequate, then setting this to \code{FALSE} can result in significantly faster simulation.}

    \item{\dots}{Additional arguments to pass to \code{\link{simulateSpectra}} or \code{\link{presetImageDef}}.}

    \item{dim}{The dimensions of the preset image.}

    \item{peakheight}{Reference intensities used for peak heights by the preset.}

    \item{peakdiff}{A reference intensity difference used for the mean peak height difference between conditions, for presets that simulate multiple conditions.}

    \item{sdsample}{A standard deviation giving the amount of variation from the true peak heights for this simulated sample.}

    \item{jitter}{Should random noise be added to the location and size of the shapes?}

    \item{center}{The center of the shape.}

    \item{size}{The size of the shape (from the center).}

    \item{shape}{What type of shape to add.}

    \item{name}{The name of the added column.}
}

\details{
    The \code{simulateSpectra()} and \code{simulateImage()} functions are used to simulate mass spectra and MS imaging experiments. They provide a great deal of control over the parameters of the simulation, including all sources of variation.

    For \code{simulateImage()}, the user should provide the design of the simulated experiment via matching columns in \code{pixelData} and \code{featureData}, where each column corresponds to different morphological substructures or differing conditions. These design data frames are returned in the \code{metadata()} of the returned object for later reference.

    A number of presets are defined by \code{presetImageDef()}, which returns only the \code{pixelData} and \code{featureData} necessary to define the experiment for \code{simulateImage()}. These can be referenced for help in understanding how to define experiments for \code{simulateImage()}.

    The preset images are:

    \itemize{
        \item{1: a centered circle}

        \item{2: a topleft circle and a bottomright square}

        \item{3: two corner squares and a centered circle}

        \item{4: a centered circle with conditions A and B in different runs}

        \item{5: a topleft circle and a bottomright square with conditions A and B in different runs}

        \item{6: two corner squares and a centered circle; the circle has conditions A and B in different runs}

        \item{7: matched pairs of circles with conditions A and B within the same runs; includes reference peaks}

        \item{8: matched pairs of circles inside squares with conditions A and B within the same runs; includes reference peaks}

        \item{9: a small sphere inside a larger sphere (3D)}
    }

    The \code{addShape()} function is provided for convenience when generating the \code{pixelData} for \code{simulateImage()}, as a simple way of adding morphological substructures using basic shapes such as squares and circles.
}

\value{
	For \code{simulateSpectra}, a \code{MassDataFrame} with elements:

	\itemize{
        \item{\code{mz}: a numeric vector of the observed m/z values}

        \item{\code{intensity}: a numeric vector or matrix of the intensities}
    }

    For \code{simulateImage}, a \code{\linkS4class{MSImagingExperiment}} object.

    For \code{addShape}, a new \code{\linkS4class{PositionDataFrame}} with a logical column added for the corresponding shape.

    For \code{presetImageDef}, a list with two elements: the \code{pixelData} and \code{featureData} to be used as input to \code{simulateImage()}.
}

\author{
	Kylie A. Bemis
}

\seealso{
	\code{\link[matter]{simspec}}
}

\examples{
set.seed(1, kind="L'Ecuyer-CMRG")

# generate a spectrum
s <- simulateSpectra(1)
plot(s$intensity ~ s$mz, type="l")

# generate a noisy low-resolution spectrum with a baseline
s <- simulateSpectra(1, baseline=2, sdnoise=0.3, resolution=100)
plot(s$intensity ~ s$mz, type="l")

# generate a high-resolution spectrum
s <- simulateSpectra(1, npeaks=100, resolution=10000)
plot(s$intensity ~ s$mz, type="l")

# generate an image
mse <- simulateImage(preset=1, npeaks=10, dim=c(10,10))
peaks <- mz(metadata(mse)$design$featureData)

image(mse, mz=peaks[c(1,4,5,6)])
plot(mse, coord=c(x=3,y=3))
}

\keyword{datagen}
