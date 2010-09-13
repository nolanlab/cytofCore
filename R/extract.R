cytofCore.extract.cells <- function(data, cols=NULL, thresh=10.0, sigma=3, num_sigma=3, noise.min.length=30, freq=77000) {
    if (!is.matrix(data)) {
		stop("'data' must be a matrix");
    }

    if (is.null(cols)) {
		cols <- 1:ncol(data)
    } else {
		cols <- match(cols,colnames(data))
		if (any(is.na(cols))) {
			stop("Invalid column specifier")
		}
    }
    
    # Apply gaussian filter to total intensity
    coeff <- dnorm((-num_sigma*sigma):(num_sigma*sigma),sd=sigma)
    d <- stats::filter(ts(rowSums(data[,cols]),frequency=freq), coeff, method="convolution", sides=2)
     
    # Cells are runs of data >= threshold, while noise is preceding region
    runs  <- rle(as.vector(d >= thresh))
    cells <- which(runs$value)
    
    found <- matrix(nrow=length(cells),ncol=length(cols)+3)
    found[,1] <- (cumsum(runs$lengths))[cells-1] + 1 # Leading push
    found[,2] <- found[,1]/freq*1000
    found[,3] <- runs$lengths[cells]  # Cell length 
    for (i in 1:nrow(found)) {
		cell_range <- seq(found[i,1],length.out=found[i,3])
		found[i,4:ncol(found)] <- colSums(data[cell_range,cols])  # Integrate raw data for cell
    }
	colnames(found) <- c("Leading_Push", "Time", "Cell_length",colnames(data)[cols])

    noise <- matrix(nrow=length(cells),ncol=length(cols)+2)
    noise[,1] <- (cumsum(runs$lengths))[cells-2] + 1 # Leading push for region preceding cell
    noise[,2] <- runs$lengths[cells-1]  # Noise region duration
    for (i in 1:nrow(noise)) {
		if (noise[i,2] > 1) {	
			noise_range <- seq(noise[i,1],length.out=noise[i,2])
			noise[i,3:ncol(noise)] <- colSums(data[noise_range,cols]) / noise[i,2]  # Produce per push noise estimate
		} else {
			noise[i,3:ncol(noise)] <- data[noise[i,1],cols]
		}
    }
	noise <- subset(noise, noise[,2] > noise.min.length)
    colnames(noise) <- c("Leading_Push", "Duration",colnames(data)[cols])
    
	list(cells=found, noise=noise, intensity=d)
}

cytofCore.filter.cells <- function(cells) {
    fail <- c()
    for (i in 1:nrow(cells$cells)) {
		beg <- cells$cells[i,1]
		len <- cells$cells[i,3]
	
		slope <- cells$intensity[seq(beg+1,length.out=len)] - cells$intensity[seq(beg,length.out=len)]
		runs  <- rle(slope > 0)
		if (length(runs$values) > 2) {  # More than one zero crossing ...
			fail <- c(fail, i)
		}
    }
    return(fail) 
}

cytofCore.extract.R <- function(imd, conf,
	pulse_thresh=3.0, num_pushes=NULL, thresh=10.0, sigma=3, num_sigma=3, min_length=10, max_length=75, noise.subtraction=TRUE,
	noise.min.length=30, slope.filter=TRUE, freq=77000) {

	dual <- cytofCore.read.imd(imd, conf=conf, pulse_thresh=pulse_thresh, num_pushes=num_pushes)
	cell <- cytofCore.extract.cells(dual$dual, thresh=thresh, sigma=sigma, num_sigma=num_sigma, noise.min.length=noise.min.length, freq=freq)

	cell$cells <- subset(cell$cells, cell$cells[,"Cell_length"] >= min_length & cell$cells[,"Cell_length"] < max_length)

	quality <- rep(1., nrow(cell$cells))
	if (slope.filter) {
		quality[cytofCore.filter.cells(cell)] <- 0.0
	}
	
	if (noise.subtraction) {
		cell$noise <- rbind(matrix(0., nrow=1, ncol=ncol(cell$noise)),cell$noise)
		noise.idxs <- findInterval(cell$cells[,"Leading_Push"],cell$noise[,"Leading_Push"])
		# 3 is start of analyte columns for noise, 4 is start for cells
		cell$cells[,4:ncol(cell$cells)] <- 
			cell$cells[,4:ncol(cell$cells)] - cell$cells[,"Cell_length"]*cell$noise[noise.idxs,3:ncol(cell$noise)]  
	}

	list(cells.found=cell$cells,quality=quality)
}

cytofCore.extract.native <- function(imd, conf, 
	pulse_thresh=3.0, num_pushes=.Machine$integer.max, thresh=10.0, sigma=3, num_sigma=3, min_length=10, max_length=75, noise.subtraction=TRUE,
	noise.min.length=30, slope.filter=FALSE, freq=77000) {

	# Load and format "conf" values for compute dual counts
	if (is.character(conf) && file.exists(conf)) {
		conf <- cytofCore.read.conf(conf)
	}
	toDual <- as.matrix(conf[,c("Intercept","Slope")])
	storage.mode(toDual) <- "double"

	# Generate smoothing filter
	smooth <-  dnorm((-num_sigma*sigma):(num_sigma*sigma),sd=sigma)

	# Extract cells
	r <- .Call(
			"CC_extract", 
			imd, toDual, 
			list(
			smooth=smooth, pulse.threshold=pulse_thresh, cell.threshold=thresh, 
			cell.min.length=as.integer(min_length), cell.max.length=as.integer(max_length),
			noise.subtraction=noise.subtraction, noise.min.length=as.integer(noise.min.length),
			slope.filter=slope.filter,
			num.pushes=as.integer(num_pushes)
			))
	r$cells.extent <- cbind(r$cells.extent[,1], r$cells.extent[,1]/freq*1000, r$cells.extent[,2]) 

	# Build return data structure to mimic output from CyToF software
	analytes <- paste(round(conf$Mass),conf$Symbol,sep="")
	colnames(r$cells.obs) <- analytes
	colnames(r$cells.extent) <- c("Leading_Push","Time","Cell_length")
	list(cells.found=cbind(r$cells.extent,r$cells.obs),quality=r$cells.quality)
}


