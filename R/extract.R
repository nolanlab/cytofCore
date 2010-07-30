cytofCore.extract.cells <- function(data, cols=NULL, thresh=10.0, sigma=3, num_sigma=3, min_length=10, max_length=75) {
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
    d <- filter(rowSums(data[,cols]), coeff, method="convolution", sides=2)
     
    # Cells are runs of data >= threshold, while noise is preceding region
    runs  <- rle(as.vector(d >= thresh))
    cells <- which(runs$value & runs$lengths >= min_length & runs$lengths <= max_length)
    
    found <- matrix(nrow=length(cells),ncol=length(cols)+2)
    found[,1] <- (cumsum(runs$lengths))[cells-1]  # Leading push
    found[,2] <- runs$lengths[cells]  # Cell length 
    for (i in 1:nrow(found)) {
	cell_range <- seq(found[i,1],length.out=found[i,2])
	found[i,3:ncol(found)] <- colSums(data[cell_range,cols])  # Integrate raw data for cell
    }
    colnames(found) <- c("Leading_Push", "Cell_length",colnames(data)[cols])

    noise <- matrix(nrow=length(cells),ncol=length(cols)+2)
    noise[,1] <- (cumsum(runs$lengths))[cells-2]  # Leading push for region preceding cell
    noise[,2] <- runs$lengths[cells-1]  # Noise region duration
    for (i in 1:nrow(noise)) {
	noise_range <- seq(noise[i,1],length.out=noise[i,2])
	noise[i,3:ncol(noise)] <- colSums(data[noise_range,cols]) / noise[i,2]  # Produce per push noise estimate
    }
    colnames(noise) <- c("Leading_Push", "Duration",colnames(data)[cols])
    list(cells=found, noise=noise)
}
