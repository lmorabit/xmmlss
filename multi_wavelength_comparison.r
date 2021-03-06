## multi-wavelength cross matching for radio multi-wavelength busy week
library('raster')
library( 'gplots' )
library( 'viridis' )
library( 'minpack.lm' )
library( 'mixtools' )
mycols <- viridis( 5 )
source('/vardy/leah/xmmlss/plotting_libraries.r')


read_in_catalogue <- function( infile ){
    ## read in a csv table with column headers
    qmc <- read.table( infile, stringsAsFactors=FALSE, sep=",", header=TRUE )
    return( qmc )
}

find_nearest_neighbours <- function( df ){

    if ( 'ALPHA_J2000' %in% colnames( df ) ){
	racol <- 'ALPHA_J2000'
	decol <- 'DELTA_J2000'
    } else {
	racol <- 'RA'
	decol <- 'DEC'
    }

    nearest_neighbours <- c()
    for ( ii in 1:dim(df)[1] ){
        distances <- cenang( df[ii,racol], df[ii,decol], df[,racol], df[,decol] )
        nearest_neighbours <- c( nearest_neighbours, min( distances[which(distances >0)] ) )
    }
    nearest_neighbours <- nearest_neighbours * 60. * 60. ## convert to arcsec
    return( nearest_neighbours )

}


plot_detection_fractions_vs_radio_size <- function( qmc, plotfile='Match_fraction_vs_size.pdf', radio_limit=0 ){

    ## get a subset of radio sources
    qmc1 <- qmc[ which( qmc$Total_flux > radio_limit ), ]

    ## flag if detected in W1 or i band
    W1_det <- which( qmc1$e_W1mag < 0.2 )
    i_det <- which( qmc1$iMeanPSFMag > 0 )
    detected <- intersect( W1_det, i_det )
    
    ## get the sizes, convert to arcsec, and log
    log_radio_dc_maj <- log10( qmc1$DC_Maj * 60 * 60 )
    my_breaks <- pretty( log_radio_dc_maj, 20 )
    
    det_frac <- c()
    for ( ii in 1:(length(my_breaks)) ) {
        ## number of sources
        bin_index <- which( log_radio_dc_maj >= my_breaks[ii] & log_radio_dc_maj < my_breaks[ii+1] )
        n_bin <- length( bin_index )
        n_detected <- length( intersect( bin_index, detected ) )
        bin_frac <- n_detected / n_bin 
        det_frac <- c( det_frac, bin_frac )
    }
      
    pdf( plotfile )
    par( mar=c(5,5,2,2) )
    plot( 10.^my_breaks, det_frac, type='l', lwd=3, axes=FALSE, xlab="", ylab="" )
    axis( 1 )
    axis( 2 )
    box( which='plot' )
    mtext( 'DC Maj size [arcsec]', side=1, line=3, cex=1.25 )
    mtext( 'Match fraction (W1 or i)', side=2, line=3, cex=1.25 )
    dev.off()

}

get_band_index <- function( qmc, band_names ){

    my_cols <- colnames( qmc )
    band_cols <- c()
    for ( bn in band_names ) band_cols <- c( band_cols, which( my_cols == bn ) )
    return( band_cols )

}

get_band_error_index <- function( qmc, band_error_names ){

    my_cols <- colnames( qmc )
    band_error_cols <- c()
    for ( ben in band_error_names ) band_error_cols <- c( band_error_cols, which( my_cols == ben ) )
    return( band_error_cols )

}




## also want to plot the detection fraction vs. radio size (DC_Maj)

plot_radio_detection_fractions <- function( qmc, bands, band_names, band_error_names, plotfile='Fraction_of_matches.pdf', radio_limit=0 ){

    ## print out some info
    cat( 'Processing catalogue with a radio flux limit of', radio_limit, 'Jy.\n' )
    
    radio_index <- which( qmc$Total_flux > radio_limit )
    n_radio_sources <- length( radio_index )
    qmc_radio_detected <- qmc[ radio_index, ]

    ## print out some info
    cat( 'There are', n_radio_sources, 'radio sources.\n' )

    ## get subtables with just the bands we're interested in
    qmc_colnames <- colnames( qmc )

    band_mag_index <- get_band_index( qmc, band_names )
    band_mag_errors_index <- get_band_error_index( qmc, band_error_names )
    col_vec <- viridis( length( bands )+1 )[1:length(bands)]
    
    ## find the number of sources
    n_radio_matches <- c()      ## all
    n_radio_matches_S <- c()    ## single
    n_radio_matches_M <- c()    ## multiple
    n_radio_matches_C <- c()    ## complex  

    ## find 5-sigma detections
    for ( ii in band_mag_errors_index ) {
        n_radio_matches <- c( n_radio_matches, length( which( qmc_radio_detected[ ,ii] < 0.2 ) ) )
        n_radio_matches_S <- c( n_radio_matches_S, length( intersect( which( qmc_radio_detected[ ,ii] < 0.2 ), which( qmc_radio_detected$S_Code == 'S' ) ) ) )
        n_radio_matches_M <- c( n_radio_matches_M, length( intersect( which( qmc_radio_detected[ ,ii] < 0.2 ), which( qmc_radio_detected$S_Code == 'M' ) ) ) )
        n_radio_matches_C <- c( n_radio_matches_C, length( intersect( which( qmc_radio_detected[ ,ii] < 0.2 ), which( qmc_radio_detected$S_Code == 'C' ) ) ) )
    }

    percentage_of_matches <- n_radio_matches / n_radio_sources
    percentage_of_matches_S <- n_radio_matches_S / length( which( qmc_radio_detected$S_Code == 'S' ) )
    percentage_of_matches_M <- n_radio_matches_M / length( which( qmc_radio_detected$S_Code == 'M' ) )
    percentage_of_matches_C <- n_radio_matches_C / length( which( qmc_radio_detected$S_Code == 'C' ) )

    ## stack to make a barplot
    my_matches <- rbind( percentage_of_matches, percentage_of_matches_S, percentage_of_matches_M, percentage_of_matches_C ) 
    my_cols <- rbind( col_vec, tcol( col_vec, trans=150 ), tcol( col_vec, trans=100 ), tcol( col_vec, trans=100 ) )
    my_dens <- rbind( rep( 200, length( percentage_of_matches ) ), rep( 75, length( percentage_of_matches ) ), rep( 50, length( percentage_of_matches ) ), rep( 25, length( percentage_of_matches ) ) )  
    my_angle <- rbind( rep( 45, length( percentage_of_matches ) ), rep( 45, length( percentage_of_matches ) ), rep( -45, length( percentage_of_matches ) ), rep( 45, length( percentage_of_matches ) ) )  

    pdf( plotfile )
    par( mar=c(5,5,2,2) )
    barplot( my_matches, names.arg = bands, col=my_cols, density=my_dens, angle=my_angle, ylim=c(0,1), beside=TRUE )
    #text( 10, 0.65, "WISE", col=mycols[1], cex=2 )
    #text( 27, 0.12, "2MASS", col=mycols[2], cex=2 )
    #text( 46, 0.475, "Pan-STARRS", col=mycols[3], cex=2 )
    #text( 62, 0.4, "UHS", col=mycols[4], cex=2 )
    box( which="plot" )
    mtext( "Fraction of Radio matches", side=2, line=3, cex=1.25 )
    legend( 'topright', c('All', 'S', 'M', 'C' ), col=c('lightgray', tcol('lightgray',trans=150), tcol('lightgray',trans=100), tcol('lightgray',trans=100) ), density=c(200,75,50,25), angle=c(45,45,-45,45), bty='n', cex=1.5 )
    dev.off()

}


plot_radio_detection_fractions_gaul <- function( qmc, plotfile='Gaussian_Fraction_of_matches.pdf', radio_limit=0 ){

    ## print out some info
    cat( 'Processing catalogue with a radio flux limit of', radio_limit, 'Jy.\n' )
        
    radio_index <- which( qmc$Total_flux > radio_limit )
    n_radio_sources <- length( radio_index )
    qmc_radio_detected <- qmc[ radio_index, ]

    ## print out some info
    cat( 'There are', n_radio_sources, 'radio sources.\n' )

    ## Get things which are duplicated
    radio_source_ids <- qmc_radio_detected$Source_id
    ## initialise a vector with number of gaussians
    n_gaussians <- rep( 1, length( radio_source_ids ) )
    ## the values in this index will be 0 for 1 gaussian (unique entry) and >0 for things with more than one entry
    duplicated_index <- duplicated( radio_source_ids ) + duplicated( radio_source_ids, fromLast=TRUE )

    ## things with more than one gaussian
    duplicated_ids <- unique( radio_source_ids[ which( duplicated_index > 0) ] )
    for ( di in duplicated_ids ){
        di_index <- which( radio_source_ids == di )
        n_gaussians[ di_index ] <- length( di_index ) 
    }

    ## divide into bins
    my_breaks <- c( 0, 1, 2, 5, 60 )


    ## get subtables with just the bands we're interested in
    band_mag_index <- get_colnames_of_bands( qmc )
    band_mag_errors_index <- get_colnames_of_band_errors( qmc )

    ## what data structure do i need to save this to?
    
    ## wise band detections
    for ( ii in band_mag_index[1:7] ) {
        ## loop through bins ...
        for ( jj in 1:(length(my_breaks)-1) ){
            bin_index <- which( n_gaussians > my_breaks[jj] & n_gaussians <= my_breaks[jj+1] )
            n_total <- length( bin_index )
            n_bin <- length( which( qmc_radio_detected[bin_index, ii] < 0.2 ) )
        }
    }
    ## panstarrs
    for ( ii in seq( 50, 58, 2 ) ) {
        n_radio_matches <- c( n_radio_matches, length( which( qmc_radio_detected[ ,ii] > 0 ) ) )
        n_radio_matches_S <- c( n_radio_matches_S, length( intersect( which( qmc_radio_detected[ ,ii] > 0 ), which( qmc_radio_detected$S_Code == 'S' ) ) ) )
        n_radio_matches_M <- c( n_radio_matches_M, length( intersect( which( qmc_radio_detected[ ,ii] > 0 ), which( qmc_radio_detected$S_Code == 'M' ) ) ) )
        n_radio_matches_C <- c( n_radio_matches_C, length( intersect( which( qmc_radio_detected[ ,ii] > 0 ), which( qmc_radio_detected$S_Code == 'C' ) ) ) )
    }
    ## UHS 
    n_radio_matches <- c( n_radio_matches, length( which( qmc$Total_flux > radio_limit & qmc$JAPERMAG3 > 0 ) ) )
    n_radio_matches_S <- c( n_radio_matches_S, length( intersect( which( qmc$Total_flux > radio_limit & qmc$JAPERMAG3 > 0 ), which( qmc_radio_detected$S_Code == 'S' ) ) ) )
    n_radio_matches_M <- c( n_radio_matches_M, length( intersect( which( qmc$Total_flux > radio_limit & qmc$JAPERMAG3 > 0 ), which( qmc_radio_detected$S_Code == 'M' ) ) ) )
    n_radio_matches_C <- c( n_radio_matches_C, length( intersect( which( qmc$Total_flux > radio_limit & qmc$JAPERMAG3 > 0 ), which( qmc_radio_detected$S_Code == 'C' ) ) ) )

    bands <- c( 'w1', 'w2', 'w3', 'w4', 'J', 'H', 'K', 'g', 'r', 'i', 'z', 'y', 'J' )
    col_vec <- c( rep( mycols[1], 4 ), rep( mycols[2], 3 ), rep( mycols[3], 5 ), mycols[4] )    
    

    percentage_of_matches <- n_radio_matches / n_radio_sources
    percentage_of_matches_S <- n_radio_matches_S / length( which( qmc_radio_detected$S_Code == 'S' ) )
    percentage_of_matches_M <- n_radio_matches_M / length( which( qmc_radio_detected$S_Code == 'M' ) )
    percentage_of_matches_C <- n_radio_matches_C / length( which( qmc_radio_detected$S_Code == 'C' ) )


    ## stack to make a barplot
    my_matches <- rbind( percentage_of_matches, percentage_of_matches_S, percentage_of_matches_M, percentage_of_matches_C ) 
    my_cols <- rbind( col_vec, tcol( col_vec, trans=150 ), tcol( col_vec, trans=100 ), tcol( col_vec, trans=100 ) )
    my_dens <- rbind( rep( 200, length( percentage_of_matches ) ), rep( 75, length( percentage_of_matches ) ), rep( 50, length( percentage_of_matches ) ), rep( 25, length( percentage_of_matches ) ) )  
    my_angle <- rbind( rep( 45, length( percentage_of_matches ) ), rep( 45, length( percentage_of_matches ) ), rep( -45, length( percentage_of_matches ) ), rep( 45, length( percentage_of_matches ) ) )  

    pdf( plotfile )
    par( mar=c(5,5,2,2) )
    barplot( my_matches, names.arg = bands, col=my_cols, density=my_dens, angle=my_angle, ylim=c(0,0.7), beside=TRUE )
    text( 10, 0.65, "WISE", col=mycols[1], cex=2 )
    text( 27, 0.12, "2MASS", col=mycols[2], cex=2 )
    text( 46, 0.475, "Pan-STARRS", col=mycols[3], cex=2 )
    text( 62, 0.4, "UHS", col=mycols[4], cex=2 )
    box( which="plot" )
    mtext( "Fraction of radio matches", side=2, line=3, cex=1.25 )
    legend( 'topright', c('All', 'S', 'M', 'C' ), col=c('lightgray', tcol('lightgray',trans=150), tcol('lightgray',trans=100), tcol('lightgray',trans=100) ), density=c(200,75,50,25), angle=c(45,45,-45,45), bty='n', cex=1.5 )
    dev.off()

}


compare_w1_with_panstarrs <- function( qmc, radio_limit=0 ){

    ## get panstarrs data
    panstarrs <- qmc[ ,seq(50,58,2)]

    p_bands <- gsub( 'MeanPSFMag', '', colnames( panstarrs ) )

    ## define things that don't depend on the bands
    my_breaks <- seq( 1, 9, 0.5 )
    valid_index <- which( qmc$e_W1mag < 0.2 )
    radio_index <- which( qmc$Total_flux > radio_limit )
    mycols <- viridis( 6 )
    
    ## start a plot
    pdf( 'W1_PanSTARRS_colours.pdf' )
    par( mar=c(5,5,2,2) )
    plot( 1, type='n', axes=FALSE, xlab="", ylab="", xlim=c(1,10), ylim=c(0,0.3) )

    ## loop through the bands
    for ( ii in 1:length(p_bands) ){
        cat( ' ... comparing W1 and ', p_bands[ii], '\n' )
        ## get colours
        panstarrs_minus_w1mag <- panstarrs[,ii] - qmc$W1mag
        
        lof_frac <- c()
        for ( jj in 1:(length(my_breaks)) ){ 
            ## get the number of matches in the bin
            bin_index <- which( panstarrs_minus_w1mag >= my_breaks[jj] & panstarrs_minus_w1mag < my_breaks[jj+1] )
            my_indices <- list( bin_index, valid_index, radio_index )
            bin_frac <- length( Reduce( intersect, my_indices ) ) / length( intersect( bin_index, valid_index ) )
            lof_frac <- c( lof_frac, bin_frac )
        }
        lines( my_breaks, lof_frac, lwd=2, col=mycols[ii] )
    }
    axis( 1 )
    axis( 2 )
    mtext( "PanSTARRS - W1 [mag]", side=1, line=3, cex=1.25 )
    mtext( "Fraction of radio matches", side=2, line=3, cex=1.25 )
    box( which="plot" )
    legend( 'topright', p_bands, col=mycols[1:5], lwd=2, bty="n" )
    dev.off()

}


create_footprint_mask <- function( df, cellsize=60, tolerance=1e-7, outfile='mask', ra_hours=FALSE, method='2dhist', exclude_halo=TRUE, use_minmax=FALSE ){

    ## allowable methods are '2dhist' (default) and '2dkde'

    #df <- read.table( infile, header=TRUE, stringsAsFactors=FALSE, sep=',' )

    if ( 'RA' %in% colnames( df ) ){
        ra_name <- 'RA'
        dec_name <- 'DEC'
    } else {
        ra_name <- 'ALPHA_J2000'
        dec_name <- 'DELTA_J2000'
    }

    ra <- df[ , ra_name ]
    ## if h
    if ( ra_hours ) ra <- ra*15
    dec <- df[ , dec_name ]
    ## cartesian projection
    ra <- ra * cos( dec * pi / 180 )
    if ( exclude_halo ) hf <- df[ ,'HALOFLAG' ]

    ## first define the area to use 
    xmin <- min( ra ) - ( cellsize / 60. / 60. )
    xmax <- max( ra ) + ( cellsize / 60. / 60. )

    ymin <- min( dec ) - ( cellsize / 60. / 60. )
    ymax <- max( dec ) + ( cellsize / 60. / 60. )

    if ( use_minmax ){
        cat( 'Using min/max RA/DEC instead of histogram\n' )
        xnsteps <- ceiling( ( xmax - xmin ) / ( cellsize / 60. / 60. ) )
        xstepsize <- ( xmax - xmin ) / xnsteps 
        xbr <- seq( xmin, xmax, xstepsize )
        ynsteps <- ceiling( ( ymax - ymin ) / ( cellsize / 60. / 60. ) )
        ystepsize <- ( ymax - ymin ) / ynsteps 
        ybr <- seq( ymin, ymax, ystepsize )
        
        counts <- matrix( 1, nrow=length(xbr), ncol=length(ybr) )

        k <- list( x.breaks=xbr, y.breaks=ybr, counts=counts )        
        

    } else {
        if ( method == '2dkde' ){
        
            cat('Using kernal density estimation.\n' )
            xnsteps <- ceiling( ( xmax - xmin ) / ( cellsize / 60. / 60. ) ) + 10
            xrange_padding <- 0.5 * ( xnsteps * ( cellsize / 60. / 60. ) - ( xmax - xmin ) )
            xl <- xmin - xrange_padding
            xu <- xmax + xrange_padding 

            ynsteps <- ceiling( ( ymax - ymin ) / ( cellsize / 60. / 60. ) ) + 10
            yrange_padding <- 0.5 * ( ynsteps * ( cellsize / 60. / 60. ) - ( ymax - ymin ) )
            yl <- ymin - yrange_padding
            yu <- ymax + yrange_padding

            ## using kernal density estimate (2d) from MASS package
            k <- kde2d( ra, dec, n=c(xnsteps,ynsteps), lims=c(xl,xu,yl,yu) )
            k$z[ which( k$z >= tolerance ) ] <- 1
            k$z[ which( k$z < tolerance ) ] <- 0
            image( k,col=gray.colors(2,start=0,end=1) )
        
        }

        if ( method == '2dhist' ){
            ## using 2d histogram from gplots package
            
            xnsteps <- ceiling( ( xmax - xmin ) / ( cellsize / 60. / 60. ) )
            ynsteps <- ceiling( ( ymax - ymin ) / ( cellsize / 60. / 60. ) )

            k <- hist2d( ra, dec, nbins=c(xnsteps,ynsteps), col=viridis(20), show=FALSE )
            k$counts[ which( k$counts > 0 ) ] <- 1

            ## exclude haloflags
            if ( exclude_halo ){
                hf_ra <- ra[ which( hf == 1 ) ]
                hf_dec <- dec[ which( hf == 1 ) ]

                ## create a mask using k$x.breaks and k$y.breaks
                k_halo <- k$counts
                k_halo[1:dim(k_halo)[1],1:dim(k_halo)[2]] <- 0
                for ( ii in 1:length(hf_ra) ){
                    ## first check if it's in the mask area
                    ra_check <- ( hf_ra[ii] >= min( k$x.breaks ) & hf_ra[ii] < max( k$x.breaks ) )
                    dec_check <- ( hf_dec[ii] >= min( k$y.breaks ) & hf_dec[ii] < max( k$y.breaks ) )
                    if ( ra_check+dec_check == 2 ){
                        bin_x <- max( which( k$x.breaks <= hf_ra[ii] ) )
                        bin_y <- max( which( k$y.breaks <= hf_dec[ii] ) )
                        k_halo[bin_x,bin_y] <- 1
                    } # end if
                } # end for
                k$counts <- k$counts - k_halo
            } # end if

            kim <- raster( list( x=k$x, y=k$y, z=k$counts ) )

            #plot( 1, type='n', xlim=c(37.5,33), ylim=c(-6,-3.5),xlab='RA [deg]', ylab='Dec [deg]')
            #rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
            #image( kim, col=gray.colors(2,start=0,end=1), add=TRUE )
            #my_decision <- readline('Is this map acceptable? (y/n): ')
	    my_decision <- 'y'
            if ( my_decision == 'y' ){
                plotfile <- paste( strsplit( outfile, '.csv' )[[1]], '_mask_',cellsize,'arcsec.pdf', sep='' )
                pdf( plotfile )
                par( mar=c(5,5,2,2) )
                plot( 1, type='n', xlim=c(37.5,33), ylim=c(-6,-3.5),xlab='RA [deg]', ylab='Dec [deg]')
                rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
                image( kim, col=gray.colors(2,start=0,end=1), add=TRUE )
                dev.off()
                ## save the mask somewhere ....
                outfile <- paste( strsplit( outfile, '.csv' )[[1]], '_mask_',cellsize,'arcsec.txt', sep='' )
                writeLines( c( paste( k$x.breaks, collapse="," ), paste( k$y.breaks, collapse="," )), con=outfile )
                write.table( k$counts, file=outfile, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE )

            } else {
                cat( 'Not saving plot. Try re-running with a different cell size.\n' )
                cat( '   (current cell size', cellsize, 'arcsec)\n' )
            }
	    #dev.off()

        }
    }

    return( k )

}

apply_mask <- function( df, my_mask, filestem='' ){

    ## check if it's a radio catalogue or not
    if ( 'RA' %in% colnames( df ) ){
        ra_name <- 'RA'
        dec_name <- 'DEC'
    } else {
        ra_name <- 'ALPHA_J2000'
        dec_name <- 'DELTA_J2000'
    }

    ## cartesian projection
    ra <- df[ , ra_name ]
    ## if h
    #if ( ra_hours ) ra <- ra*15
    dec <- df[ , dec_name ]
    ## cartesian projection
    ra <- ra * cos( dec * pi / 180 )

    mask_vec <- rep( 0, dim(df)[1] )
    ## loop through data frame
    for ( ii in 1:dim(df)[1] ){
        ## first check if it's in the mask area
        ra_check <- ( ra[ ii ] >= min( my_mask$x.breaks ) & ra[ ii ] < max( my_mask$x.breaks ) )
        dec_check <- ( dec[ ii ] >= min( my_mask$y.breaks ) & dec[ ii ] < max( my_mask$y.breaks ) )
        if ( ra_check+dec_check == 2 ){
            bin_x <- max( which( my_mask$x.breaks <= ra[ ii ] ) )
            bin_y <- max( which( my_mask$y.breaks <= dec[ ii ] ) ) 
            if ( my_mask$counts[bin_x,bin_y] == 1 ) mask_vec[ii] <- 1           
        }
    }
    masked_df <- df[ which( mask_vec == 1 ), ]
    outfile <- paste( filestem, '_masked.csv', sep='' )
    write.table( masked_df, file=outfile, row.names=FALSE, quote=FALSE, sep="," )
    return( masked_df )

}

apply_masked_regions <- function( df ){

    df_cols <- colnames( df )

    if ( 'RA' %in% df_cols ){
        ra_name <- 'RA'
        dec_name <- 'DEC'
    } else {
        ra_name <- 'ALPHA_J2000'
        dec_name <- 'DELTA_J2000'
    }

    ## get the mask definition information
    area_def_file <- '/vardy/leah/data/xmmlss/mask_definition/Kband_area_definition.dat'
    halo_regions <- '/vardy/leah/data/xmmlss/mask_definition/ds9.reg'

    area_def <- read.table( area_def_file, stringsAsFactors=FALSE, header=TRUE )

    include_area_x <- area_def$include_x
    include_area_y <- area_def$include_y 

    ## trim down to only the 'include' area
    include_index_x <- which( df[ , ra_name ] > min( include_area_x ) & df[ , ra_name ] < max( include_area_x ) )
    include_index_y <- which( df[ , dec_name ] > min( include_area_y ) & df[ , dec_name ] < max( include_area_y ) )
    include_index <- intersect( include_index_x, include_index_y )

    ## loop through the exclude regions
    area_def_cols <- colnames( area_def )
    excols <- area_def_cols[ which( grepl( 'exclude_x', area_def_cols ) )]

    exclude_index <- c()
    for ( ii in 1:length( excols ) ){
        excl_x <- area_def[ , paste( 'exclude_x', ii, sep='' ) ]
        excl_y <- area_def[ , paste( 'exclude_y', ii, sep='' ) ]
        excl_index_x <- which( df[ , ra_name ] > min( excl_x ) & df[ , ra_name ] < max( excl_x ) )
        excl_index_y <- which( df[ , dec_name ] > min( excl_y ) & df[ , dec_name ] < max( excl_y ) )
        excl_index <- intersect( excl_index_x, excl_index_y )
        exclude_index <- c( exclude_index, excl_index )
    }
    exclude_index <- unique( exclude_index )

    ## find where things are in include but not exclude.
    include_in_exclude <- ( include_index %in% exclude_index )

    final_index <- include_index[ which( !include_in_exclude ) ]

    new_df <- df[ final_index, ]

    ## now exclude the halo regions
    ## first read the file
    ds9_reg <- readLines( halo_regions )
    ds9_reg <- ds9_reg[3:length(ds9_reg)]

    exclude_index <- c()
    for ( ds9r in ds9_reg ){
        tmp <- strsplit( ds9r, '(', fixed=TRUE )[[1]]
        if ( tmp[1] == 'circle' ){
            x_cen <- as.numeric( strsplit( tmp[2], ',' )[[1]][1] )
            y_cen <- as.numeric( strsplit( tmp[2], ',' )[[1]][2] )
            radius <- as.numeric( substr( strsplit( tmp[2], ',' )[[1]][3], 1, 6 ) ) 
            distances <- cenang( x_cen, y_cen, new_df[,ra_name], new_df[,dec_name] ) * 60 * 60 ## convert to arcsec
            excl_index <- which( distances < radius )
            exclude_index <- c( exclude_index, excl_index )

        } #endif
    } #endfor ds9r
    exclude_index <- unique( exclude_index )

    tmp_index <- seq( 1, dim( new_df )[1] )

    circle_index <- tmp_index[ which( !( tmp_index %in% exclude_index ) ) ]

    ## deal with polygons --which are actually only squares
    exclude_index <- c()
    for ( ds9r in ds9_reg ){
        tmp <- strsplit( ds9r, '(', fixed=TRUE )[[1]]
        if ( tmp[1] == 'polygon' ){
            tmp1 <- strsplit( tmp[2], ',' )[[1]]
            tmp1[length(tmp1)] <- substr( tmp1[length(tmp1)], 1, nchar(tmp1[length(tmp1)])-1 )
            tmp1 <- as.numeric( tmp1 )
            point1 <- list( x=tmp1[1], y=tmp1[2] )
            point2 <- list( x=tmp1[3], y=tmp1[4] )
            point3 <- list( x=tmp1[5], y=tmp1[6] )
            point4 <- list( x=tmp1[7], y=tmp1[8] )

            ## line between point 1 and point 2
            p12_m <- ( point2$y - point1$y ) / ( point2$x - point1$x )
            p12_b <- point1$y - p12_m*point1$x 
                
            my_y <- new_df[,ra_name] * p12_m + p12_b 
            excl_index <- which( new_df[,dec_name] < my_y )

            ## line between point 3 and point 4
            p43_m <- ( point3$y - point4$y ) / ( point3$x - point4$x )
            p43_b <- point4$y - p43_m*point4$x 

            my_y <- new_df[,ra_name] * p43_m + p43_b 
            excl_index <- intersect( excl_index, which( new_df[,dec_name] > my_y ) )

            ## line between point 1 and point 4
            p41_m <- ( point1$y - point4$y ) / ( point1$x - point4$x )
            p41_b <- point4$y - p41_m*point4$x 

            my_y <- new_df[,ra_name] * p41_m + p41_b 
            excl_index <- intersect( excl_index, which( new_df[,dec_name] < my_y ) )

            ## line between point 2 and point 3
            p32_m <- ( point2$y - point3$y ) / ( point2$x - point3$x )
            p32_b <- point3$y - p32_m*point3$x 
                
            my_y <- new_df[,ra_name] * p32_m + p32_b 
            excl_index <- intersect( excl_index, which( new_df[,dec_name] > my_y ) )

            exclude_index <- c( exclude_index, excl_index )

        }
    }
    exclude_index <- unique( exclude_index )

    polygon_index <- tmp_index[ which( !(tmp_index %in% exclude_index ) ) ]

    final_include_index <- intersect( circle_index, polygon_index )
        
    df_out <- new_df[ final_include_index, ]

    return( df_out )

}


find_ra_dec_plot_limits <- function( ra, dec, padding=30 ){

    xmin <- min( ra ) - ( padding / 60. / 60. )
    xmax <- max( ra ) + ( padding / 60. / 60. )
    xnsteps <- ceiling( ( xmax - xmin ) / ( padding / 60. / 60. ) ) 
    xrange_padding <- 0.5 * ( xnsteps * ( padding / 60. / 60. ) - ( xmax - xmin ) )
    xl <- xmin - xrange_padding
    xu <- xmax + xrange_padding 

    ymin <- min( dec ) - ( padding / 60. / 60. )
    ymax <- max( dec ) + ( padding / 60. / 60. )
    ynsteps <- ceiling( ( ymax - ymin ) / ( padding / 60. / 60. ) ) 
    yrange_padding <- 0.5 * ( ynsteps * ( padding / 60. / 60. ) - ( ymax - ymin ) )
    yl <- ymin - yrange_padding
    yu <- ymax + yrange_padding

    my_limits <- c( xl, xu, yl, yu )
    return( my_limits )

}

plot_sky_distribution <- function( parm1='', ra1=0, de1=0, parm2='', ra2=0, de2=0, parm3='', ra3=0, de3=0 ){

    mycols=viridis(4)
    ## set up the plotting
    xlab <- "RA [hh:mm:ss]"
    ylab <- "Dec [dd:mm:ss]"
    comb_ra <- c( ra1, ra2, ra3 )
    comb_dec <- c( de1, de2, de3 )
    plot_limits <- find_ra_dec_plot_limits( comb_ra[ which( comb_ra != 0 )], comb_dec[ which( comb_dec != 0 )] )

    ## open file
    filestem <- 'sky_distribution'
    myparms <- c( parm1, parm2, parm3 )
    parm_names <- c()
    for ( p in myparms ) if ( nchar( p ) > 0 ) {
        filestem <- paste( p, '_', filestem, sep='' )
        parm_names <- c( parm_names, gsub( '_', ' ', p ) )
    }
    pdf( paste( filestem, '.pdf', sep='' ) )
    par( mar=c(5,5,2,2) )
    ## don't forget to reverse RA
    plot( ra1, de1, pch=16, cex=0.2, xlab=xlab, ylab=ylab, xlim=rev(plot_limits[1:2]), ylim=plot_limits[3:4], axes=FALSE, col=mycols[1] )
    if ( length(ra2) > 1 ) points( ra2, de2, pch=16, cex=0.6, col=mycols[2] )
    if ( length(ra3) > 1 ) points( ra3, de3, pch=16, cex=0.6, col=mycols[3] )
    nticks <- 5
    ## format x-axis
    #atvec <- seq(plot_limits[1],plot_limits[2],(plot_limits[2]-plot_limits[1])/(nticks-1))
    #stepsize <- ( plot_limits[2] - plot_limits[1] ) / ( length( atvec ) - 1 )
    #plothms <- deg2hms( rev(seq(plot_limits[1],plot_limits[2],stepsize)), type='cat', sep=':' )
    atvec <- seq(par("xaxp")[1],par("xaxp")[2],(par("xaxp")[2]-par("xaxp")[1])/(nticks-1))
    stepsize <- ( par("xaxp")[2] - par("xaxp")[1] ) / ( length(atvec)-1 )
    plothms <- deg2hms( rev(seq(par("xaxp")[1],par("xaxp")[2],stepsize)), type='cat', sep=':' )
    axis( side=1, at=rev(atvec), labels=plothms, cex=2 )
    ## format y-axis
    atvec <- seq(par("yaxp")[1],par("yaxp")[2],(par("yaxp")[2]-par("yaxp")[1])/(nticks-1))
    stepsize <- ( par("yaxp")[2] - par("yaxp")[1] ) / ( length(atvec)-1 )
    plotdms <- deg2dms( seq(par("yaxp")[1],par("yaxp")[2],stepsize), type='cat', sep=':' )
    axis( side=2, at=atvec, labels=plotdms, cex=2 )
    if ( length(parm_names) > 1 ) legend( 'topright', parm_names, col=mycols[1:length(parm_names)], pch=rep(16,length(parm_names)), pt.cex=c(0.4,rep(0.7,(length(parm_names)-1))), bty='n' )
    box( which='plot' )
    dev.off()
        
}

find_number_no_counterparts <- function( radio_df, video_df, radii ){

    ## what is the number of sources with no possible counterparts
    n_counterparts <- matrix( nrow=dim(radio_df)[1], ncol=length(radii) )
    pb <- txtProgressBar( min=0, max=dim(radio_df)[1], style=3 )
    for ( ii in 1:dim(radio_df)[1] ){
        ## calculate the distance from the source to all other sources, convert to arcsec 
        distances <- cenang( radio_df$RA[ii], radio_df$DEC[ii], video_df$ALPHA_J2000, video_df$DELTA_J2000 ) * 60. * 60.
        ## loop through the radii to find the number of counterparts
        for ( jj in 1:length(radii) ) n_counterparts[ii,jj] <- length( which( distances <= radii[jj] ) )
        setTxtProgressBar( pb, ii )        
    }
    close( pb )
    cat( '\n' )

    n_blanks <- c()
    for ( jj in 1:length(radii) ) n_blanks <- c( n_blanks, length( which( n_counterparts[,jj] == 0 ) ) )
    return( n_blanks )

}

find_Q0_fleuren <- function( radio_df, video_df, radii ){

    ## first check if this has already been done
    ## get the band first for the filename
    tmp <- colnames( video_df )
    tmp1 <- tmp[which( grepl( 'MAG_', tmp ) )]
    my_band <- strsplit( tmp1, '_' )[[1]][1]

    f_file <- paste( 'Fleuren_no_counterparts_', my_band, '.csv', sep='' )

    if ( !file.exists( f_file ) ){

	    ## for the real catalogue, find the number of sources that don't have counterparts as a function of r
	    cat( 'Finding blanks for REAL catalogue.\n' )
	    REAL_no_counterparts <- find_number_no_counterparts( radio_df, video_df, radii )
    
	    ## make a random radio catalogue
	    n_srcs <- dim( radio_df )[1]
	    min_RA <- min( video_df$ALPHA_J2000 )
	    max_RA <- max( video_df$ALPHA_J2000 )
	    min_DEC <- min( video_df$DELTA_J2000 )
	    max_DEC <- max( video_df$DELTA_J2000 )

	    ## randomly generate RA/DEC, scaling to the appropriate area in the sky
	    set.seed(30)
	    rand_DEC <- runif( n_srcs * 4 )
	    rand_DEC <- rand_DEC * ( max_DEC-min_DEC ) + min_DEC

	    set.seed(20)
	    rand_RA <- runif( n_srcs * 4 )
	    rand_RA <- rand_RA * ( max_RA-min_RA ) + min_RA 
	    ## do i need a correction for cos( DEC ) ?

	    random_RA_DEC <- data.frame( rand_RA, rand_DEC, stringsAsFactors=FALSE )
	    colnames( random_RA_DEC ) <- c( 'RA', 'DEC' )

        ## apply the mask
        masked_random <- apply_masked_regions( random_RA_DEC )

	    ## trim down the array so it's the right size 
	    if ( dim( masked_random )[1] >= n_srcs ) random_RA_DEC <- masked_random[ 1:n_srcs, ] else cat( 'NOT ENOUGH DATA POINTS IN THE RANDOM CATALOGUE.\n' )

	    ## find the number of no counterparts as a function of radius
	    cat( 'Finding blanks for RANDOM catalogue.\n' )
	    RANDOM_no_counterparts <- find_number_no_counterparts( random_RA_DEC, video_df, radii )

	    ## take the ratio
	    no_counterpart_ratio <- REAL_no_counterparts / RANDOM_no_counterparts
    
	    cp_ratio <- data.frame( REAL_no_counterparts, RANDOM_no_counterparts, no_counterpart_ratio )
	    write.table( cp_ratio, file=f_file, quote=FALSE, row.names=FALSE, sep=',' )

    } else {
	cat( 'Finding blanks already done for this band, reading file.\n' )
	cp_ratio <- read.table( f_file, stringsAsFactors=FALSE, header=TRUE, sep=',' )
	REAL_no_counterparts <- cp_ratio$REAL_no_counterparts
	RANDOM_no_counterparts <- cp_ratio$RANDOM_no_counterparts
	no_counterpart_ratio <- cp_ratio$no_counterpart_ratio
    }

    ## fit the model to the data
    y <- no_counterpart_ratio
    result <- nlsLM( y ~ 1 - q_0 * ( 1 - exp( -radii^2 / ( 2 * sig^2 ) ) ), start=list( q_0 = 1.3, sig=1 ) )
    coeffs <- coef( summary( result ) )
    q_0 <- coeffs[1,1]
    q_0_err <- coeffs[1,2]

    cat( '... Q_0 found:', q_0, '+/-', q_0_err, '\n' )

    q0_df <- data.frame( q_0, q_0_err )

    ## make a plot 
    ## get the band first for the filename

    mycols <- viridis( 6 )
    pdf( paste( 'Q_0_estimate_', my_band, '.pdf', sep='' ) )
    lplot( radii, RANDOM_no_counterparts/dim(radio_df)[1], type='l', col=mycols[2], x_lab='Radius [arcsec]', y_lab='1-Q(r)', lwd=3, lty=3, ylim=c(0,1) )
    lines( radii, REAL_no_counterparts/dim(radio_df)[1], col=mycols[4], lwd=3 )
    lines( radii, predict( result ), lwd=3, col='gray' )
    points( radii, no_counterpart_ratio, pch=18, cex=1.5 )
    legend( 'topright', c( 'Random', 'Real', 'Data', 'Fit' ), col=c(mycols[2],mycols[4],'black','gray'), lty=c(3,1,NA,1), pch=c(NA,NA,18,NA), lwd=3, bty='n', cex=1.25 )
    text( max(radii)/2-0.5, 0.6, bquote( Q[0] == .(q_0) %+-% .(q_0_err) ) )
    dev.off()

    return( q0_df )    

}

calculate_hist <- function( mydata, mybreaks ){

    mycounts <- c()
    for ( ii in 2:length(mybreaks) ) mycounts <- c( mycounts, length( which( mydata >= mybreaks[ii-1] & mydata < mybreaks[ii] ) ) )
    break_step <- mybreaks[2] - mybreaks[1]
    mymids <- mybreaks[1:(length(mybreaks)-1)]+0.5*break_step

    df <- list( counts=mycounts, mids=mymids, breaks=mybreaks )
    return( df )

}

calculate_positional_errors <- function( df, beam_size=5 ){

	## CALCULATE POSITIONAL INFORMATION FOR THE RADIO
	## -- positional errors: condon
	#pos_noise <- sqrt( my_radio_cat$e_RA^2. + my_radio_cat$e_DEC^2. )
	## -- positional errors: ivison
	pos_noise <- 0.6 * beam_size / ( df$Total_flux / df$e_Total_flux )
	## add in a calibration errors	
	sigma_pos <- sqrt( 0.1^2. + pos_noise^2. )
	r_max <- 5 * max( sqrt( sigma_pos^2 ) )
	cat( 'r_max is', r_max, 'arcsec.\n' )
    return( list( r_max=r_max, sigma_pos=sigma_pos ) )
}

calculate_nm <- function( df, mag_cols, mag_bins, area_arcsec ){

    ## CALCULATE nm -- density of background sources in VIDEO band
    ## histogram the sources and find the source density
    nmhist <- calculate_hist( df[,mag_cols[1]], mag_bins )
    nm <- nmhist$counts / area_arcsec  
    return( nm )

}

make_mag_bins <- function( df, mag_cols ){

    mag_bins <- seq( floor(min( df[ , mag_cols[1]] )), ceiling(max( df[,mag_cols[1]] )), 0.4 )
    ## check that it spans the entire range
    if ( ceiling(max( df[ ,mag_cols[1]] )) > max( mag_bins ) ) mag_bins <- c( mag_bins, max(mag_bins)+0.4 )
    if ( floor(min( df[ ,mag_cols[1]] )) < min( mag_bins ) ) mag_bins <- c( min(mag_bins)-0.4, mag_bins )
    return( mag_bins )

}

calculate_matched_mags <- function( my_band, radio_df, video_df, r_max, filename='' ){

	n_radio_sources <- dim(radio_df)[1]
	if ( !file.exists( filename ) ){
		match_magnitudes <- c()
        radio_ids <- c()
		## this is the bit that takes a long time ... start a progress bar
		cat( 'Finding magnitudes of all sources within r_max.\n' )
		pb <- txtProgressBar(min = 0, max = n_radio_sources, style = 3)
		for ( ii in 1:n_radio_sources ){
			## find distances
			distances <- cenang( radio_df$RA[ii], radio_df$DEC[ii], video_df$ALPHA_J2000, video_df$DELTA_J2000 )*60*60
			candidate_index <- which( distances <= r_max ) 
			match_magnitudes <- c( match_magnitudes, video_df[candidate_index,mag_cols[1]] )
            radio_ids <- c( radio_ids, rep( radio_df$Source_id[ii], length( candidate_index ) ) )
		## update the progress bar
		setTxtProgressBar( pb, ii )
		}
		cat( '\n' )
		close( pb )
        out_df <- data.frame( radio_ids, match_magnitudes )
        write.table( out_df, file=filename, row.names=FALSE, sep=',', quote=FALSE )
	} else out_df <- read.table( filename, header=TRUE, sep=',' )
    return( out_df )
}

read_VLA_data <- function( vla_file, mask=NULL ){

    ## read in the data
    VLA <- read.table( vla_file,  skip=30, stringsAsFactors=FALSE, header=FALSE, sep="," )
    ## label the columns
    colnames( VLA ) <- c('ID', 'RA', 'DEC', 'e_RA', 'e_DEC', 'Total_flux', 'e_Total_flux', 'Peak_flux', 'e_Peak_flux', 'rms', 'DC_Maj', 'e_DC_Maj', 'DC_Min', 'e_DC_Min', 'DC_PA', 'e_DC_PA', 'resolved', 'Gaussian_ID', 'Source_id', 'Island_id' )

    outfile <- '/vardy/leah/data/xmmlss/VLA/13B-308_080916.csv'

    ## write a file that can be read in by topcat
    write.table( VLA, file=outfile, row.names=FALSE, sep=",", quote=FALSE )
    ## return the table    
    return( VLA )

}

prepare_radio_data <- function( starting_file, snr_cutoff=5, outfile='sum_VLA.csv' ){

    VLA <- read_VLA_data( starting_file )

    ## get rid of duplicated rows
    VLA <- VLA[ which( !duplicated( VLA ) ), ]

    ## and find signal to noise
    snr <- VLA$Peak_flux / ( VLA$rms / 1e3 ) 
    good_index <- which( snr >= snr_cutoff )

    my_VLA <- VLA[ good_index, ]
    
    ## count unique sources
    n_unique <- length( unique( my_VLA$Source_id ) )
    
    ## sum the multiple-component sources 
    my_sum_VLA <- sum_components( my_VLA )

    nn <- find_nearest_neighbours( my_sum_VLA )
    my_sum_VLA[ ,'nearest_neighbour' ] <- nn

    write.table( my_sum_VLA, file=outfile, row.names=FALSE, quote=FALSE, sep=',' )

    return( my_sum_VLA )

}

sum_components <- function( df, thresh=5 ){ #, outfile='sum_VLA.csv' ){

    ##--- CREATE AN EMPTY DATA FRAME
    sum_df <- data.frame( ID=character(), RA=double(), DEC=double(), e_RA=double(), e_DEC=double(), Total_flux=double(), e_Total_flux=double(), Peak_flux=double(), e_Peak_flux=double(), rms=integer(), DC_Maj=double(), e_DC_Maj=double(), DC_Min=double(), e_DC_Min=double(), DC_PA=double(), e_DC_PA=double(), resolved=integer(), Gaussian_ID=integer(), Source_id=integer(), Island_id=integer(), stringsAsFactors=FALSE )

    ## and an empty id string
    source_type <- c()

    ## get the source id's
    source_ids <- df$Source_id

    ##--- DUPLICATED SOURCE ID's (i.e., multiple gaussian components)    
    ## get a duplicated index
    duplicated_flags <- duplicated( source_ids ) + duplicated( source_ids, fromLast=TRUE )
    duplicated_ids <- unique( source_ids[ which( duplicated_flags > 0 ) ] )

    ## loop through the duplicated sources
    low_score_ids <- c()
    for ( d_id in duplicated_ids ){

        tmp_df <- df[ which( df$Source_id == d_id ), ]

        ## identify artefacts ?
        ## do this using e_RA, e_DEC, e_DC_Maj, e_DC_Min, e_DC_PA
        ## everything gets a score of 1 if it's more than the standard deviation of 
        ## all the components
        scores <- rep( 0, dim(tmp_df)[1] )
        scores[ which( tmp_df$e_RA > sd( tmp_df$e_RA ) ) ] <- scores[ which( tmp_df$e_RA > sd( tmp_df$e_RA ) ) ] + 1.
        scores[ which( tmp_df$e_DEC > sd( tmp_df$e_DEC ) ) ] <- scores[ which( tmp_df$e_DEC > sd( tmp_df$e_DEC ) ) ] + 1.
        scores[ which( tmp_df$e_DC_Maj > sd( tmp_df$e_DC_Maj ) ) ] <- scores[ which( tmp_df$e_DC_Maj > sd( tmp_df$e_DC_Maj ) ) ] + 1.
        scores[ which( tmp_df$e_DC_Min > sd( tmp_df$e_DC_Min ) ) ] <- scores[ which( tmp_df$e_DC_Min > sd( tmp_df$e_DC_Min ) ) ] + 1.
        scores[ which( tmp_df$e_DC_PA > sd( tmp_df$e_DC_PA ) ) ] <- scores[ which( tmp_df$e_DC_PA > sd( tmp_df$e_DC_PA ) ) ] + 1.
        ## check for low scores
        good_scores <- which( scores < 2 )
        if ( length( good_scores ) < 1 ){
            low_score_ids <- c( low_score_ids, d_id )
            source_type <- c( source_type, 'C' )
        } else { 
            tmp_df <- tmp_df[ good_scores, ] 
            source_type <- c( source_type, 'M' )
        }

        ## find geometric mean of RA, DEC (and errors)
        mean_RA <- mean( tmp_df$RA )
        mean_DEC <- mean( tmp_df$DEC )
        e_mean_RA <- sqrt( sum( tmp_df$e_RA^2. ) )
        e_mean_DEC <- sqrt( sum( tmp_df$e_DEC^2. ) )
        ## convert RA, DEC to ID
        mean_ID <- paste( 'J', substr( deg2hms( mean_RA, type='cat', sep='' ), 1, 8 ), substr( deg2dms( mean_DEC, type='cat', sep='' ), 1, 7 ), sep="" )
        ## sum the total flux & errors
        sum_total_flux <- sum( tmp_df$Total_flux )
        e_sum_total_flux <- sqrt( sum( tmp_df$e_Total_flux^2. ) )
        ## take the max peak flux & error
        max_peak_flux <- max( tmp_df$Peak_flux )
        e_max_peak_flux <- tmp_df$e_Peak_flux[ which( tmp_df$Peak_flux == max_peak_flux ) ]
        if ( length( e_max_peak_flux ) > 0 ) e_max_peak_flux <- max( e_max_peak_flux )
        ## average the rms
        avg_rms <- mean( tmp_df$rms[ which( tmp_df$rms > 0 ) ] )
        ## use DC_Maj, DC_Min, and DC_PA to find largest size ... ?
        my_sizes <- find_largest_size( tmp_df$RA, tmp_df$DEC, tmp_df$DC_PA, tmp_df$DC_Maj, tmp_df$DC_Min )
        my_DC_Maj <- my_sizes[1]
        e_my_DC_Maj <- my_sizes[2]
        my_DC_Min <- my_sizes[3]
        e_my_DC_Min <- my_sizes[4]
        my_DC_PA <- my_sizes[5]
        e_my_DC_PA <- my_sizes[6]

        ## set resolved = sum( resolved )
        sum_resolved <- sum( tmp_df$resolved )
        ## set Gaussian_ID = which( Peak_flux == max( Peak_flux ) )
        gauss_id <- tmp_df$Gaussian_ID[ which( tmp_df$Peak_flux == max( tmp_df$Peak_flux ) ) ]
        if ( length( gauss_id ) > 0 ) gauss_id <- min( gauss_id )
        ## set Island_ID = which( Peak_flux == max( Peak_flux ) )
        isl_id <- tmp_df$Island_id[ which( tmp_df$Peak_flux == max( tmp_df$Peak_flux ) ) ]
        if ( length( isl_id ) > 0 ) isl_id <- min( isl_id )
        ## Source_id = d_id

        ##--- ADD TO DATA FRAME
        my_row <- data.frame( mean_ID, mean_RA, mean_DEC, e_mean_RA, e_mean_DEC, sum_total_flux, e_sum_total_flux, max_peak_flux, e_max_peak_flux, avg_rms, my_DC_Maj, e_my_DC_Maj, my_DC_Min, e_my_DC_Min, my_DC_PA, e_my_DC_PA, sum_resolved, gauss_id, d_id, isl_id, stringsAsFactors=FALSE )
        colnames( my_row ) <- colnames( sum_df )
        sum_df <- merge( sum_df, my_row, all=TRUE )

    }
    
    cat( '... done processing duplicated Source_ids.\n' )

    ## that means everything else is a single gaussian
    single_ids <- source_ids[ which( duplicated_flags == 0 ) ]
    for ( s_id in single_ids ){
        sum_df <- merge( sum_df, df[which( df$Source_id == s_id ),], all=TRUE )
        source_type <- c( source_type, 'S' )
    }

    sum_df[,"S_Code"] <- source_type

    cat( ' -- There are', length( single_ids ), 'single sources. (', 100*length( single_ids )/length( source_ids ), '%)\n' )
    n_dup <- length( duplicated_ids )
    cat( ' -- There are', n_dup, 'duplicated Source_ids. (', 100*n_dup/length( source_ids ), '%)\n' )
    cat( ' -- There are', n_dup-length(low_score_ids), 'source groups. (', 100*(n_dup-length(low_score_ids))/length( source_ids ), '%)\n' )

    #write.table( sum_df, file=outfile, row.names=FALSE, quote=FALSE, sep="," )
    return( sum_df )

}

find_largest_size <- function( ra, dec, pa, l_maj, l_min, do_plot=FALSE ){

    ## convert major and minor axes to sigma and degrees
    sigma_maj <- l_maj / 60. / 60. / (2*sqrt(2*log(2)))
    sigma_min <- l_min / 60. / 60. / (2*sqrt(2*log(2)))

    ## and PA to radians
    pa_rad <- ( 180 - pa ) * pi / 180.

    if ( do_plot ){
        ## start a plot 
        x_lim <- rev( range( ra ) )
        x_lim[1] <- x_lim[1] + 0.5*( max( x_lim ) - min( x_lim ) )
        x_lim[2] <- x_lim[2] - 0.5*( max( x_lim ) - min( x_lim ) )
        y_lim <- range( dec )
        y_lim[1] <- y_lim[1] - 0.5*( max( y_lim ) - min( y_lim ) )
        y_lim[2] <- y_lim[2] + 0.5*( max( y_lim ) - min( y_lim ) )
        plot( ra, dec, xlim=x_lim, ylim=y_lim )
    }
    
    ## get the coordinates for each gaussian (default is 95% limit)
    coordinates_x <- c()
    coordinates_y <- c()
    for ( ii in 1:length(ra) ){
        rho_sigmax_sigmay <- cos( pa_rad[ii] ) * sigma_maj[ii] * sigma_min[ii]
        my_mat <- matrix( c( sigma_maj[ii]^2., rho_sigmax_sigmay, rho_sigmax_sigmay, sigma_min[ii]^2 ), 2, 2 )
        if ( do_plot ) values <- ellipse( mu=c( ra[ii], dec[ii] ), sigma=my_mat ) else values <- ellipse( mu=c( ra[ii], dec[ii] ), sigma=my_mat, draw=FALSE )
        coordinates_x <- c( coordinates_x, values[,1] )
        coordinates_y <- c( coordinates_y, values[,2] )
    }

    ## for now, just take the rectangle that everything fits in.
    ## later, do more awesome stuff. 
    xvec <- c( min( coordinates_x ), min(coordinates_x ), max( coordinates_x ), max( coordinates_x ) )
    yvec <- c( min( coordinates_y ), max( coordinates_y ), max( coordinates_y), min( coordinates_y ) )
    if ( do_plot ) polygon( xvec, yvec, border='green' )
    
    x_size <- ( max( coordinates_x ) - min( coordinates_x ) ) * 60. * 60.
    y_size <- ( max( coordinates_y ) - min( coordinates_y ) ) * 60. * 60.
    pa <- 0.
    x_cen <- ( max( coordinates_x ) - min( coordinates_x ) )/2 + min( coordinates_x )
    y_cen <- ( max( coordinates_y ) - min( coordinates_y ) )/2 + min( coordinates_y )

    ## CASE 1: x > y 
    if ( x_size > y_size ){
        maj_size <- sqrt( x_size^2. + y_size^2. )
        min_size <- y_size / cos( asin( y_size / maj_size ) )
        ## determine which direction the major axis points
        ## can do this by making a mask with the gaussians and summing in quad_1 and quad_4 ...
        ## for now, do it by the number of points (and remember that RA is backwards)
        quad_1 <- length( which( coordinates_x > x_cen & coordinates_y > y_cen ) )
        quad_4 <- length( which( coordinates_x < x_cen & coordinates_y > y_cen ) )
        pa <- asin( y_size / maj_size ) * 180 / pi
        if ( quad_1 > quad_4 ) pa <- pa + 90.
    }            

    ## CASE 2: y > x
    if ( x_size < y_size ){
        maj_size <- sqrt( x_size^2. + y_size^2. )
        min_size <- x_size / sin( acos( x_size / maj_size ) )
        ## determine which direction the major axis points
        ## can do this by making a mask with the gaussians and summing in quad_1 and quad_4 ...
        ## for now, do it by the number of points (and remember that RA is backwards)
        quad_1 <- length( which( coordinates_x > x_cen & coordinates_y > y_cen ) )
        quad_4 <- length( which( coordinates_x < x_cen & coordinates_y > y_cen ) )
        pa <- asin( x_size / maj_size ) * 180 / pi
        if ( quad_4 < quad_1 ) pa <- pa + 90
    }

    ## CASE 3: y == x
    if ( x_size == y_size ){
        maj_size <- sqrt( x_size^2. + y_size^2. )   
        min_size <- maj_size
        quad_1 <- length( which( coordinates_x > x_cen & coordinates_y > y_cen ) )
        quad_4 <- length( which( coordinates_x < x_cen & coordinates_y > y_cen ) )
        if ( quad_1 > quad_4 ) pa <- 135 else pa <- 45
    }


    maj_size <- max( x_size, y_size )
    min_size <- min( x_size, y_size )
    if ( x_size > y_size ) pa <- 0 else pa <- 90    

    ## for now, set the errors to zero
    my_size <- c( maj_size, 0., min_size, 0., pa, 0. )
    if ( do_plot ) dev.off()
    return( my_size )

}

box_area <- function( xvec, yvec ){

    area <- ( max( xvec ) - min( xvec ) ) * ( max( yvec ) - min( yvec ) )
    return( area )
}

get_video_area_arcsec <- function(){

    df_cols <- colnames( df )

    if ( 'RA' %in% df_cols ){
        ra_name <- 'RA'
        dec_name <- 'DEC'
    } else {
        ra_name <- 'ALPHA_J2000'
        dec_name <- 'DELTA_J2000'
    }

    ## get the mask definition information
    area_def_file <- '/vardy/leah/data/xmmlss/mask_definition/Kband_area_definition.dat'
    halo_regions <- '/vardy/leah/data/xmmlss/mask_definition/ds9.reg'

    area_def <- read.table( area_def_file, stringsAsFactors=FALSE, header=TRUE )

    ## convert to arcsec
    ad <- area_def * 60 * 60

    main_area <- box_area( ad$include_x, ad$include_y )

    ## subtract off bottom area
    bottom_area <- box_area( ad$exclude_x1, ad$exclude_y1 ) + box_area( c(max(ad$exclude_x2)-max(ad$exclude_x1)), ad$exclude_y2 ) + box_area( ad$exclude_x3, c(max(ad$exclude_y3),max(ad$exclude_y1)) ) + box_area( ad$exclude_x4, c(max(ad$exclude_y4),max(ad$exclude_y2)) )

    ## subtract off top area
    top_area <- box_area( ad$exclude_x5, ad$exclude_y5 ) + box_area( ad$exclude_x6, c(min(ad$exclude_y6),min(ad$exclude_y5)) ) + box_area( ad$exclude_x7, c(min(ad$exclude_y6),min(ad$exclude_y7)) )

    ## interior areas
    interior_areas <- box_area( ad$exclude_x8, ad$exclude_y8 ) + box_area( ad$exclude_x9, ad$exclude_y9 ) + box_area( ad$exclude_x10, ad$exclude_y10 )
    
    main_area <- main_area - bottom_area - top_area - interior_areas
    
    ## now exclude the halo regions
    ## first read the file
    ds9_reg <- readLines( halo_regions )
    ds9_reg <- ds9_reg[3:length(ds9_reg)]

    circle_areas <- c()
    for ( ds9r in ds9_reg ){
        tmp <- strsplit( ds9r, '(', fixed=TRUE )[[1]]
        if ( tmp[1] == 'circle' ){
            radius <- as.numeric( substr( strsplit( tmp[2], ',' )[[1]][3], 1, 6 ) ) ## this is in arcsec
            circle_areas <- c( circle_areas, pi*radius^2 )
        } #endif
    } #endfor ds9r
    circle_area <- sum( circle_areas )

    main_area <- main_area - circle_area

    ## deal with polygons --which are actually only squares
    polygon_areas <- c()
    for ( ds9r in ds9_reg ){
        tmp <- strsplit( ds9r, '(', fixed=TRUE )[[1]]
        if ( tmp[1] == 'polygon' ){
            tmp1 <- strsplit( tmp[2], ',' )[[1]]
            tmp1[length(tmp1)] <- substr( tmp1[length(tmp1)], 1, nchar(tmp1[length(tmp1)])-1 )
            tmp1 <- as.numeric( tmp1 ) * 60 * 60
            point1 <- list( x=tmp1[1], y=tmp1[2] )
            point2 <- list( x=tmp1[3], y=tmp1[4] )
            point3 <- list( x=tmp1[5], y=tmp1[6] )
            point4 <- list( x=tmp1[7], y=tmp1[8] )

            mat1 <- as.matrix( cbind(c(point1$x, point1$y ), c(point2$x, point2$y)) )
            mat2 <- as.matrix( cbind(c(point2$x, point2$y ), c(point3$x, point3$y)) )
            mat3 <- as.matrix( cbind(c(point3$x, point3$y ), c(point4$x, point4$y)) )
            mat4 <- as.matrix( cbind(c(point4$x, point4$y ), c(point1$x, point1$y)) )
    
            sum_dets <- 0.5 * ( det( mat1 ) + det( mat2 ) + det( mat3 ) + det( mat4 ) )

            polygon_areas <- c( polygon_areas, sum_dets )

        }
    }

    polygon_area <- sum( polygon_areas )
    
    total_area <- main_area - polygon_area

    return( total_area )

}

