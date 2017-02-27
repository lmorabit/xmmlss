## multi-wavelength cross matching for radio multi-wavelength busy week
library('raster')
library( 'gplots' )
library( 'viridis' )
mycols <- viridis( 5 )
source('~/Dropbox/scripts/plotting_libraries.r')


read_in_catalogue <- function( infile ){
    ## read in a csv table with column headers
    qmc <- read.table( infile, stringsAsFactors=FALSE, sep=",", header=TRUE )
    return( qmc )
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


create_footprint_mask <- function( df, ra_col_name, dec_col_name, cellsize=60, tolerance=1e-7, outfile='mask', ra_hours=FALSE, method='2dhist', exclude_halo=TRUE ){

    ## allowable methods are '2dhist' (default) and '2dkde'

    #df <- read.table( infile, header=TRUE, stringsAsFactors=FALSE, sep=',' )

    ra <- df[ , ra_col_name ]
    ## if h
    if ( ra_hours ) ra <- ra*15
    dec <- df[ , dec_col_name ]
    ## cartesian projection
    ra <- ra * cos( dec * pi / 180 )
    hf <- df[ ,'HALOFLAG' ]

    ## first define the area to use 
    xmin <- min( ra ) - ( cellsize / 60. / 60. )
    xmax <- max( ra ) + ( cellsize / 60. / 60. )

    ymin <- min( dec ) - ( cellsize / 60. / 60. )
    ymax <- max( dec ) + ( cellsize / 60. / 60. )

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

        k <- hist2d( ra, dec, nbins=c(xnsteps,ynsteps), col=viridis(20) )
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
                ra_check <- ( hf_ra[ii] >= min( k$x.breaks ) & hf_ra[ii] <= max( k$x.breaks ) )
                dec_check <- ( hf_dec[ii] >= min( k$y.breaks ) & hf_dec[ii] <= max( k$y.breaks ) )
                if ( ra_check+dec_check == 2 ){
                    bin_x <- max( which( k$x.breaks < hf_ra[ii] ) )
                    bin_y <- max( which( hf_dec[ii] > k$y.breaks ) )
                    k_halo[bin_x,bin_y] <- 1
                } # end if
            } # end for
        } # end if
            
        k$counts <- k$counts - k_halo

        kim <- raster( list( x=k$x, y=k$y, z=k$counts ) )

        plot( 1, type='n', xlim=c(37.5,33), ylim=c(-6,-3.5),xlab='RA [deg]', ylab='Dec [deg]')
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
        image( kim, col=gray.colors(2,start=0,end=1), add=TRUE )
        my_decision <- readline('Is this map acceptable? (y/n): ')
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

    mask_vec <- rep( 0, dim(df)[1] )
    ## loop through data frame
    for ( ii in 1:dim(df)[1] ){
        ## first check if it's in the mask area
        ra_check <- ( df[ ii, ra_name ] >= min( my_mask$x.breaks ) & df[ ii, ra_name ] <= max( my_mask$x.breaks ) )
        dec_check <- ( df[ ii, dec_name ] >= min( my_mask$y.breaks ) & df[ ii, dec_name ] <= max( my_mask$y.breaks ) )
        if ( ra_check+dec_check == 2 ){
            ## ra will always be positive 
            bin_x <- max( which( my_mask$x.breaks < df[ ii, ra_name ] ) )
            ## dec can be negative
            bin_y <- max( which( df[ ii, dec_name ] >= my_mask$y.breaks ) ) ## all values are greater than zero
            ## case 3: if they straddle zero ... ?
            if ( my_mask$counts[bin_x,bin_y] == 1 ) mask_vec[ii] <- 1            
        }
    }
    masked_df <- df[ which( mask_vec == 1 ), ]
    outfile <- paste( filestem, '_masked.csv', sep='' )
    write.table( masked_df, file=outfile, row.names=FALSE, quote=FALSE, sep="," )
    return( masked_df )

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



























