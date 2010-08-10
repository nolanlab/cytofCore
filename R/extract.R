cytofCore.extract.cells <- function(data, cols=NULL, thresh=10.0, sigma=3, num_sigma=3, min_length=10, max_length=75, freq=90000) {
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
    found[,1] <- (cumsum(runs$lengths))[cells-1]  # Leading push
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
