cytofCore.extract.cells <- function(data, cols=NULL, thresh=10.0, sigma=3, num_sigma=3, min_length=10, max_length=75, freq=77000) {
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
    cells <- which(runs$value & runs$lengths >= min_length & runs$lengths <= max_length)
    
    found <- matrix(nrow=length(cells),ncol=length(cols)+3)
    found[,1] <- (cumsum(runs$lengths))[cells-1] + 1 # Leading push
    found[,2] <- found[,1]/freq
    found[,3] <- runs$lengths[cells]  # Cell length 
    for (i in 1:nrow(found)) {
	cell_range <- seq(found[i,1],length.out=found[i,3])
	found[i,4:ncol(found)] <- colSums(data[cell_range,cols])  # Integrate raw data for cell
    }
    colnames(found) <- c("Leading_Push", "Time", "Cell_length",colnames(data)[cols])

    noise <- matrix(nrow=length(cells),ncol=length(cols)+2)
    noise[,1] <- (cumsum(runs$lengths))[cells-2]  # Leading push for region preceding cell
    noise[,2] <- runs$lengths[cells-1]  # Noise region duration
    for (i in 1:nrow(noise)) {
	if (noise[i,2] > 1) {	
	    noise_range <- seq(noise[i,1],length.out=noise[i,2])
	    noise[i,3:ncol(noise)] <- colSums(data[noise_range,cols]) / noise[i,2]  # Produce per push noise estimate
	} else {
	    noise[i,3:ncol(noise)] <- data[noise[i,1],cols]
	}
    }
    colnames(noise) <- c("Leading_Push", "Duration",colnames(data)[cols])
    list(cells=found, noise=noise, intensity=d)
}

cytofCore.extract.native <- function(file, conf, 
	pulse_thresh=3.0, num_pushes=.Machine$integer.max, thresh=10.0, sigma=3, num_sigma=3, min_length=10, max_length=75, noise.subtraction=TRUE,
	freq=77000) {

	# Load and format "conf" values for compute dual counts
	if (is.character(conf) && file.exists(conf)) {
		conf <- read.table(conf,header=TRUE)
	}
	toDual <- as.matrix(conf[,c("Intercept","Slope")])
	storage.mode(toDual) <- "double"

	# Generate smoothing filter
	smooth <-  dnorm((-num_sigma*sigma):(num_sigma*sigma),sd=sigma)

	# Extract cells
	r <- .Call(
			"CC_extract", 
			file, toDual, 
			list(
			smooth=smooth, pulse.threshold=pulse_thresh, cell.threshold=thresh, 
			cell.min.length=as.integer(min_length), cell.max.length=as.integer(max_length),
			noise.subtraction=noise.subtraction,
			num.pushes=as.integer(num_pushes)
			))
	r$cells.extent <- cbind(r$cells.extent[,1], r$cells.extent[,1]/freq*1000, r$cells.extent[,2]) 

	# Build return data structure to mimic output from CyToF software
	analytes <- paste(round(conf$Mass),conf$Symbol,sep="")
	colnames(r$cells.obs) <- analytes
	colnames(r$cells.extent) <- c("Leading_Push","Time","Cell_length")
	list(cells.found=cbind(r$cells.extent,r$cells.obs),quality=r$cells.quality)

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
