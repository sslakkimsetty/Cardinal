
#### Peak processing ####
## ----------------------

setMethod("peakProcess", "MSImagingExperiment_OR_Arrays",
	function(object, ref,
		spectra = "intensity", index = "mz",
		method = c("diff", "sd", "mad", "quantile", "filter", "cwt"),
		SNR = 2, type = c("height", "area"),
		tolerance = NA, units = c("ppm", "mz"),
		sampleSize = NA, filterFreq = TRUE, outfile = NULL,
		verbose = getCardinalVerbose(), chunkopts = list(),
		BPPARAM = getCardinalBPPARAM(), ...)
{
	if ( !is.na(sampleSize) ||
		(!isCentroided(object) && !missing(ref) && !is.null(ref)) )
	{
		# need to do peak picking
		if ( missing(ref) || is.null(ref) )
		{
			# create reference peaks from sample spectra
			if ( sampleSize < 1 ) {
				# sample size is a proportion
				n <- ceiling(sampleSize * length(object))
				perc <- 100 * sampleSize
			} else if ( sampleSize > 0 ) {
				# sample size is a count
				n <- min(sampleSize, length(object))
				perc <- round(100 * n / length(object))
			} else {
				.Error("'sampleSize' must be positive")
			}
			label <- if (n != 1L) "spectra" else "spectrum"
			.Log("processing peaks for ", n, " ", label, " ",
				"(~", perc, "% of data)",
				message=verbose)
			i <- seq.default(1L, length(object), length.out=n)
			ref <- peakProcess(object[i],
				spectra=spectra, index=index,
				method=method, SNR=SNR, type=type,
				tolerance=tolerance, units=units, filterFreq=filterFreq,
				verbose=verbose, chunkopts=chunkopts,
				BPPARAM=BPPARAM, ...)
			domain <- mz(ref)
		} else {
			if ( is(ref, "MSImagingExperiment") || is(ref, "MassDataFrame") )
				ref <- mz(ref)
			domain <- as.numeric(ref)
		}
		# extract the peaks based on reference
		.Log("extracting reference peaks from all spectra",
			message=verbose)
		object <- peakPick(object, ref=ref,
			tolerance=tolerance, units=units, type=type)
		object <- process(object,
			spectra=spectra, index=index,
			domain=domain, outfile=outfile,
			verbose=verbose, chunkopts=chunkopts,
			BPPARAM=BPPARAM, ...)
		if ( is(ref, "MSImagingExperiment") )
			featureData(object) <- featureData(ref)
	} else {
		# check for peak picking
		if ( isCentroided(object) ) {
			if ( length(processingData(object)) == 0L &&
				!is.sparse(spectra(object, spectra)) &&
				!is(object, "MSImagingArrays") )
			{
				.Log("peaks are already processed",
					message=verbose)
				return(object)
			} else {
				.Log("peaks are already detected",
					message=verbose)
			}
		} else {
			# pick peaks on all spectra
			object <- peakPick(object,
				method=method, SNR=SNR, type=type)
		}
		# align peaks across all spectra
		object <- peakAlign(object, ref=ref,
			spectra=spectra, index=index, outfile=outfile,
			tolerance=tolerance, units=units,
			verbose=verbose, chunkopts=chunkopts,
			BPPARAM=BPPARAM, ...)
		# filter peaks
		if ( !is.null(featureData(object)[["count"]]) &&
			(isTRUE(filterFreq) || filterFreq > 0) )
		{
			if ( is.numeric(filterFreq) ) {
				if ( filterFreq < 1 ) {
					# filterFeq is a proportion
					n <- ceiling(filterFreq * length(object))
				} else if ( filterFreq > 0) {
					# filterFeq is a count
					n <- as.integer(filterFreq)
				} else {
					.Error("'filterFreq' must be positive")
				}
			} else {
				# remove singleton peaks
				n <- 1L
			}
			label <- if (n / length(object) < 0.01) "<" else "~"
			.Log("filtering to keep only peaks with counts > ", n, " ",
				"(", label, round(100 * n / length(object), digits=2L),
				"% of considered spectra)",
				message=verbose)
			object <- object[featureData(object)[["count"]] > n,]
		}
	}
	# return object
	.Log("processed to ", nrow(object), " peaks",
		message=verbose)
	if ( validObject(object) )
		object
})


#### Peak alignment ####
## ---------------------

setMethod("peakAlign", "MSImagingExperiment",
	function(object, ref,
		spectra = "intensity", index = "mz",
		binfun = "min", binratio = 2,
		tolerance = NA, units = c("ppm", "mz"), ...)
{
	if ( !missing(ref) ) {
		if ( is(ref, "MSImagingExperiment") || is(ref, "MassDataFrame") )
			ref <- mz(ref)
	}
	if ( missing(units) && !missing(tolerance) )
		units <- get_units_from_names(tolerance, units)
	units <- switch(match.arg(units), ppm="relative", mz="absolute")
	if ( !is.na(tolerance) )
		tolerance <- switch(units,
			relative=1e-6 * tolerance,
			absolute=tolerance)
	ans <- callNextMethod(object, ref=ref,
		spectra=spectra, index=index,
		binfun=binfun, binratio=binratio,
		tolerance=tolerance, units=units, ...)
	spectraData <- spectraData(ans)
	featureData <- as(featureData(ans), "MassDataFrame")
	new("MSImagingExperiment", spectraData=spectraData,
		featureData=featureData, elementMetadata=pixelData(ans),
		experimentData=experimentData(object), centroided=TRUE,
		metadata=metadata(ans), processing=list())
})

setMethod("peakAlign", "MSImagingArrays",
	function(object, ref,
		spectra = "intensity", index = "mz",
		binfun = "min", binratio = 2,
		tolerance = NA, units = c("ppm", "mz"), ...)
{
	if ( !missing(ref) ) {
		if ( is(ref, "MSImagingExperiment") || is(ref, "MassDataFrame") )
			ref <- mz(ref)
	}
	if ( missing(units) && !missing(tolerance) )
		units <- get_units_from_names(tolerance, units)
	units <- switch(match.arg(units), ppm="relative", mz="absolute")
	if ( !is.na(tolerance) )
		tolerance <- switch(units,
			relative=1e-6 * tolerance,
			absolute=tolerance)
	ans <- callNextMethod(object, ref=ref,
		spectra=spectra, index=index,
		binfun=binfun, binratio=binratio,
		tolerance=tolerance, units=units, ...)
	spectraData <- spectraData(ans)
	featureData <- as(featureData(ans), "MassDataFrame")
	new("MSImagingExperiment", spectraData=spectraData,
		featureData=featureData, elementMetadata=pixelData(ans),
		experimentData=experimentData(object), centroided=TRUE,
		metadata=metadata(ans), processing=list())
})

setMethod("peakAlign", "SpectralImagingExperiment",
	function(object, ref,
		spectra = "intensity", index = NULL,
		binfun = "min", binratio = 2,
		tolerance = NA, units = c("relative", "absolute"),
		verbose = getCardinalVerbose(), chunkopts = list(),
		BPPARAM = getCardinalBPPARAM(), ...)
{
	if ( length(processingData(object)) > 0L )
		object <- process(object, spectra=spectra, index=index,
			verbose=verbose, chunkopts=chunkopts,
			BPPARAM=BPPARAM, ...)
	if ( missing(units) && !missing(tolerance) )
		units <- get_units_from_names(tolerance, units)
	units <- match.arg(units)
	if ( length(index) > 1L )
		.Error("more than 1 'index' array not supported")
	snm <- spectra
	inm <- index
	spectra <- spectra(object, snm)
	if ( is.null(inm) ) {
		inm <- "index"
		domain <- seq_len(nrow(object))
	} else {
		domain <- featureData(object)[[inm]]
		if ( is.null(domain) )
			.Error("index ", sQuote(inm), " not found")
	}
	if ( !missing(ref) && !is.null(ref) ) {
		if ( isTRUE(all.equal(domain, ref)) ) {
			.Log("peaks are already aligned",
				message=verbose)
			return(object)
		}
	}
	if ( is.sparse(spectra) ) {
		index <- atomindex(spectra)
		spectra <- atomdata(spectra)
	} else {
		.Error("nothing to align for spectra ", sQuote(snm), "; ",
			"has peakPick() been used?")
	}
	.Log("detected ~", round(mean(lengths(index)), digits=1L),
		" peaks per spectrum",
		message=verbose)
	.peakAlign(object, ref=ref, spectra=spectra, index=index,
		domain=domain, spectraname=snm, indexname=inm,
		binfun=binfun, binratio=binratio,
		tolerance=tolerance, units=units,
		verbose=verbose, chunkopts=chunkopts,
		BPPARAM=BPPARAM)
})

setMethod("peakAlign", "SpectralImagingArrays",
	function(object, ref,
		spectra = "intensity", index = NULL,
		binfun = "min", binratio = 2,
		tolerance = NA, units = c("relative", "absolute"),
		verbose = getCardinalVerbose(), chunkopts = list(),
		BPPARAM = getCardinalBPPARAM(), ...)
{
	if ( length(processingData(object)) > 0L )
		object <- process(object, spectra=spectra, index=index,
			verbose=verbose, chunkopts=chunkopts,
			BPPARAM=BPPARAM, ...)
	if ( missing(units) && !missing(tolerance) )
		units <- get_units_from_names(tolerance, units)
	units <- match.arg(units)
	if ( length(index) > 1L )
		.Error("more than 1 'index' array not supported")
	snm <- spectra
	inm <- index
	spectra <- spectra(object, snm)
	if ( is.null(inm) ) {
		inm <- "index"
		index <- lapply(lengths(spectra), seq_len)
	} else {
		index <- spectra(object, inm)
		if ( is.null(index) )
			.Error("index ", sQuote(inm), " not found")
	}
	.Log("detected ~", round(mean(lengths(index)), digits=1L),
		" peaks per spectrum",
		message=verbose)
	if ( missing(ref) || is.null(ref) ) {
		domain <- NULL
	} else {
		domain <- ref
	}
	.peakAlign(object, ref=ref, spectra=spectra, index=index,
		domain=domain, spectraname=snm, indexname=inm,
		binfun=binfun, binratio=binratio,
		tolerance=tolerance, units=units,
		verbose=verbose, chunkopts=chunkopts,
		BPPARAM=BPPARAM)
})

.peakAlign <- function(object, ref, spectra, index, domain,
	spectraname, indexname, binfun, binratio, tolerance, units,
	verbose, chunkopts, BPPARAM)
{
	tol.ref <- switch(units, relative="x", absolute="abs")
	if ( is.null(domain) || is.na(tolerance) ) {
		width <- match.arg(binfun, c("median", "min", "max", "mean"))
		.Log("summarizing peak gaps for alignment",
			message=verbose)
		.Log("using bin function ", sQuote(width),
			" to summarize peak gaps across spectra",
			message=verbose)
		indexbins <- estimateDomain(index, width=width, units=units,
			verbose=verbose, chunkopts=chunkopts, BPPARAM=BPPARAM)
	}
	if ( is.na(tolerance) ) {
		# estimate tolerance as (binratio x min peak-to-peak gap)
		tol <- binratio * estres(indexbins, ref=tol.ref)
		tol <- switch(units,
			relative=round(2 * tol, digits=6L) * 0.5,
			absolute=round(tol, digits=4L))
		.Log("estimated ", units, " tolerance of ", tol,
			message=verbose)
	} else {
		# validate user-specified tolerance
		tol <- setNames(unname(tolerance), units)
	}
	if ( is.null(domain) ) {
		# set peak bins estimated from index
		res <- estres(indexbins, ref=tol.ref)
		res <- switch(units,
			relative=round(2 * res, digits=6L) * 0.5,
			absolute=round(res, digits=4L))
		domain <- indexbins
	} else {
		# set peak bins to (tolerance / binratio)
		res <- tol / binratio
		domain <- switch(units,
			relative=seq_rel(min(domain), max(domain), by=res),
			absolute=seq(min(domain), max(domain), by=res))
	}
	if ( missing(ref) || is.null(ref) ) {
		.Log("using bin ratio of ", binratio,
			" to create peak bins (per tolerance half-window)",
			message=verbose)
		.Log("using peak bins with ", units,
			" resolution of ", res,
			message=verbose)
		.Log("binning peaks to create shared reference",
			message=verbose)
		FUN <- isofun(function(x, domain, tol, tol.ref) {
			matter::binpeaks(x, domain=domain, tol=tol, tol.ref=tol.ref,
				merge=FALSE, na.drop=FALSE)
		}, CardinalEnv())
		peaks <- chunk_lapply(index, FUN,
			domain=domain, tol=tol, tol.ref=tol.ref,
			simplify=matter::stat_c,
			verbose=verbose, chunkopts=chunkopts,
			BPPARAM=BPPARAM)
		.Log("merging peak bins with ", units,
			" centroid differences <= ", tol,
			message=verbose)
		peaks <- peaks[!is.na(peaks)]
		peaks <- mergepeaks(peaks, tol=tol, tol.ref=tol.ref)
		n <- nobs(peaks)
		ref <- structure(as.vector(peaks), n=n)
	} else {
		n <- NULL
	}
	if ( verbose ) {
		ppm <- switch(units,
			relative=paste0("(", 1e6 * tol, " ppm)"),
			absolute="")
		.Log("aligned to ", length(ref),
			" reference peaks with ", units,
			" tolerance ", tol, " ", ppm,
			message=verbose)
	}
	spectra <- sparse_mat(index=index,
		data=spectra, domain=ref,
		nrow=length(ref), ncol=length(object),
		tolerance=tol, sampler="none")
	spectraData <- SpectraArrays(setNames(list(spectra), spectraname))
	featureData <- DataFrame(setNames(list(ref), indexname))
	if ( !is.null(n) ) {
		featureData[["count"]] <- n
		featureData[["freq"]] <- n / length(object)
	}
	label <- "peak alignment"
	metadata <- list(
		tolerance=unname(tol), units=units,
		binfun=binfun, binratio=binratio)
	metadata <- setNames(list(metadata), label)
	metadata <- setNames(list(metadata), .processing_id())
	metadata <- c(metadata(object), metadata)
	names(metadata) <- make.unique(names(metadata))
	new("SpectralImagingExperiment", spectraData=spectraData,
		featureData=featureData, elementMetadata=pixelData(object),
		metadata=metadata, processing=list())
}


#### Peak picking ####
## -------------------

setMethod("peakPick", "MSImagingExperiment",
	function(object, ref,
		method = c("diff", "sd", "mad", "quantile", "filter", "cwt"),
		SNR = 2, type = c("height", "area"),
		tolerance = NA, units = c("ppm", "mz"), ...)
{
	if ( !missing(ref) ) {
		if ( is(ref, "MSImagingExperiment") || is(ref, "MassDataFrame") )
			ref <- mz(ref)
	}
	if ( isCentroided(object) ) {
		.Warn("object is already centroided")
	} else {
		centroided(object) <- TRUE
	}
	if ( missing(units) && !missing(tolerance) )
		units <- get_units_from_names(tolerance, units)
	units <- switch(match.arg(units), ppm="relative", mz="absolute")
	if ( !is.na(tolerance) )
		tolerance <- switch(units,
			relative=1e-6 * tolerance,
			absolute=tolerance)
	callNextMethod(object, ref=ref, method=method, SNR=SNR,
		type=type, tolerance=tolerance, units=units, ...)
})

setMethod("peakPick", "MSImagingArrays",
	function(object, ref,
		method = c("diff", "sd", "mad", "quantile", "filter", "cwt"),
		SNR = 2, type = c("height", "area"),
		tolerance = NA, units = c("ppm", "mz"), ...)
{
	if ( !missing(ref) ) {
		if ( is(ref, "MSImagingExperiment") || is(ref, "MassDataFrame") )
			ref <- mz(ref)
	}
	if ( isCentroided(object) ) {
		.Warn("object is already centroided")
	} else {
		centroided(object) <- TRUE
	}
	if ( missing(units) && !missing(tolerance) )
		units <- get_units_from_names(tolerance, units)
	units <- switch(match.arg(units), ppm="relative", mz="absolute")
	if ( !is.na(tolerance) )
		tolerance <- switch(units,
			relative=1e-6 * tolerance,
			absolute=tolerance)
	callNextMethod(object, ref=ref, method=method, SNR=SNR,
		type=type, tolerance=tolerance, units=units, ...)
})

setMethod("peakPick", "SpectralImagingData",
	function(object, ref,
		method = c("diff", "sd", "mad", "quantile", "filter", "cwt"),
		SNR = 2, type = c("height", "area"),
		tolerance = NA, units = c("relative", "absolute"), ...)
{
	method <- match.arg(method)
	type <- match.arg(type)
	if ( missing(ref) || is.null(ref) ) {
		if ( !is.na(tolerance) )
			.Warn("no 'ref' given so 'tolerance' will be ignored")
		if ( method == "cwt" ) {
			FUN <- .peakPick_cwt
		} else {
			FUN <- .peakPick
		}
		addProcessing(object, FUN,
			label=paste0(type, " peak picking"),
			method=method, SNR=SNR, type=type, ...)
	} else {
		if ( missing(units) && !missing(tolerance) )
			units <- get_units_from_names(tolerance, units)
		units <- match.arg(units)
		if ( is.unsorted(ref) )
			ref <- sort(ref)
		tol.ref <- switch(units, relative="x", absolute="abs")
		if ( is.na(tolerance) ) {
			tol <- 0.5 * estres(ref, ref=tol.ref)
		} else {
			tol <- tolerance
		}
		FUN <- .peakPick_ref
		addProcessing(object, FUN,
			label=paste0(type, " peak picking"),
			ref=ref, tol=tol, tol.ref=tol.ref, type=type, ...)
	}
})

.peakPick <- function(x, t, method, ..., SNR = 2, type = "height")
{
	peaks <- matter::findpeaks(x, noise=method, snr=SNR, relheight=NULL, ...)
	if ( type == "height" ) {
		cbind(t[peaks], x[peaks])
	} else if ( type == "area" ) {
		cbind(t[peaks], matter::peakareas(x, peaks, domain=t))
	} else {
		.Error("invalid peak type: ", sQuote(type))
	}
}

.peakPick_cwt <- function(x, t, method, ..., SNR = 2, type = "height")
{
	peaks <- matter::findpeaks_cwt(x, snr=SNR, ...)
	if ( type == "height" ) {
		cbind(t[peaks], x[peaks])
	} else if ( type == "area" ) {
		cbind(t[peaks], matter::peakareas(x, peaks, domain=t))
	} else {
		.Error("invalid peak type: ", sQuote(type))
	}
}

.peakPick_ref <- function(x, t, ref, tol, tol.ref,..., type = "height")
{
	peaks <- matter::findpeaks(x, relheight=NULL, bounds=FALSE, ...)
	hits <- bsearch(ref, t[peaks], tol=tol, tol.ref=tol.ref)
	nz <- !is.na(hits)
	values <- numeric(length(ref))
	if ( type == "height" ) {
		values[nz] <- matter::peakheights(x, peaks[hits[nz]])
	} else if ( type == "area" ) {
		values[nz] <- matter::peakareas(x, peaks[hits[nz]], domain=t)
	} else {
		.Error("invalid peak type: ", sQuote(type))
	}
	cbind(ref, values)
}

