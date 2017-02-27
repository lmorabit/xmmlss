## PROCESS VIDEO RADIO DATA
##  written by Leah K. Morabito
##  started on 25 Jan 2017
##  last updated ...
##  contact: lmorabit@gmail.com

## read in the VLA data

library( 'astro' )
library( 'celestial' )
library( 'mixtools' )

## band information
band_wavelengths <<- c( 2.149, 1.646, 1.254, 1.021, 1.020, 0.9256, 0.878, 0.7764, 0.6404, 0.472, 0.3538 )
band_order <<- c( 'K', 'H', 'J', 'Y', 'YC', 'ZC', 'ZV', 'I', 'R', 'G', 'U' )

## STILL TO DO:
#- source counts

read_VLA_data <- function( vla_file, mask=NULL ){

    ## read in the data
    VLA <- read.table( vla_file,  skip=30, stringsAsFactors=FALSE, header=FALSE, sep="," )
    ## label the columns
    colnames( VLA ) <- c('ID', 'RA', 'DEC', 'e_RA', 'e_DEC', 'Total_flux', 'e_Total_flux', 'Peak_flux', 'e_Peak_flux', 'rms', 'DC_Maj', 'e_DC_Maj', 'DC_Min', 'e_DC_Min', 'DC_PA', 'e_DC_PA', 'resolved', 'Gaussian_ID', 'Source_id', 'Island_id' )

    outfile <- 'VLA/13B-308_080916.csv'

    ## write a file that can be read in by topcat
    write.table( VLA, file=outfile, row.names=FALSE, sep=",", quote=FALSE )
    ## return the table    
    return( VLA )

}

## plot some distributions
plot_distributions <- function( df, plotfile='distributions.pdf' ){

    ## HISTOGRAMS
    #-- DC_Maj
    #-- DC_Min
    #-- DC_PA
    #-- Total_flux
    #-- Peak_flux
    ## PLOTS
    #-- Total vs. Peak flux

    ## do this in landscape pdf
    pdf( plotfile, paper='a4r' )
    #par( mar=c(5,5,2,2) )
    par( mar=c(5,5,2,2), oma = c(0,0,0,0) )
    par( mfrow=c(2,2) )
    par( pty="s" )

    hist( df$DC_Maj, axes=FALSE, xlab="", ylab="", main="" ) 
    axis( 1 )
    axis( 2 )
    box( which="plot" )
    mtext( 'DC Major axis', side=1, line=3 )
    mtext( 'Frequency', side=2, line=3 )
    hist( df$DC_Min, axes=FALSE, xlab="", ylab="", main="" )
    axis( 1 )
    axis( 2 )
    box( which="plot" )
    mtext( 'DC Minor axis', side=1, line=3 )
    mtext( 'Frequency', side=2, line=3 )
    hist( df$DC_PA, axes=FALSE, xlab="", ylab="", main="" )
    axis( 1 )
    axis( 2 )
    box( which="plot" )
    mtext( 'DC Position Angle', side=1, line=3 )
    mtext( 'Frequency', side=2, line=3 )
    plot( 1, type='n', axes=FALSE, xlab="", ylab="", main="" )
    hist( log10( df$Total_flux ), axes=FALSE, xlab="", ylab="", main="" )
    axis( 1 )
    axis( 2 )
    box( which="plot" )
    mtext( 'Total Flux [log(mJy)]', side=1, line=3 )
    mtext( 'Frequency', side=2, line=3 )
    hist( log10( df$Peak_flux ), axes=FALSE, xlab="", ylab="", main="" )
    axis( 1 )
    axis( 2 )
    box( which="plot" )
    mtext( 'Peak Flux [log(mJy/beam)]', side=1, line=3 )
    mtext( 'Frequency', side=2, line=3 )
    plot( df$Total_flux, df$Peak_flux, axes=FALSE, xlab="", ylab="", main="", cex=0.8 )
    lines( seq(0,1000), seq(0,1000), lty=2 )
    axis( 1 )
    axis( 2 )
    box( which="plot" )
    mtext( 'Total flux [mJy]', side=1, line=3 )
    mtext( 'Peak flux [mJy/beam]', side=2, line=3 )
    
    dev.off()

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

## sum up components
sum_components <- function( df, thresh=5 ){

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

        ## sigma threshold for flux
        #tmp_df <- tmp_df[ which( tmp_df$Total_flux >= thresh*1e-3*tmp_df$rms ), ]

        ## identify artefacts
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

    write.table( sum_df, file='sum_VLA.csv', row.names=FALSE, quote=FALSE, sep="," )
    return( sum_df )
}

cross_match_positional_errors <- function( df, bands, plotfile='Cross-matched-positional-errors.pdf' ){

    col_vec <- viridis( length( bands )+1 )[1:length(bands)]

    band_ra <- paste( bands, '_ALPHA_J2000', sep='' )
    band_dec <- paste( bands, '_DELTA_J2000', sep='' )

    pdf( plotfile )
    par( mar=c(5,5,2,2) )
    plot( 1, type='n', xlim=c( -4, 4 ), ylim=c( -4, 4 ), axes=FALSE, xlab='RA offset [arcsec]', ylab='Dec offset [arcsec]', main='' )
    grid()
    
    for ( ii in 1:length(bands) ){
        ra_diff <- ( df$RA - df[ ,band_ra[ii] ] ) * 60. * 60. 
        dec_diff <- ( df$DEC - df[ , band_dec[ii] ] ) * 60. * 60. 
        points( ra_diff, dec_diff, pch='.', col=col_vec[ii] )
    }
    axis( 1 )
    axis( 2 )
    box( which='plot' )
    dev.off()

}

plot_radio_nearest_neighbours <- function( df, plotfile='Nearest_neighbours.pdf' ){

    df <- df[ which( is.finite( df$RA ) ), ]

    nearest_neighbours <- c()
    for ( ii in 1:dim(df)[1] ){
        distances <- cenang( df$RA[ii], df$DEC[ii], df$RA, df$DEC )
        nearest_neighbours <- c( nearest_neighbours, min( distances[which(distances >0)] ) )
    }
    nearest_neighbours <- nearest_neighbours * 60. * 60. 

    pdf( plotfile )
    par( mar=c(5,5,2,2) )
    hist( nearest_neighbours, main='', xlab='Nearest neighbour distance [arcsec]', breaks=40 )
    box( which='plot' )
    dev.off()

}



plot_photometry <- function( source_id, catalogue ){

    system( paste('grep',source_id,catalogue,'> tmp.out' ) )
    my_cols <- strsplit( readLines( catalogue, n=1 ), ',' )[[1]]
    A <- read.table( 'tmp.out', stringsAsFactors=FALSE, sep=',', header=FALSE ) 
    colnames( A ) <- my_cols
    my_source <- A[ which( A$ID == source_id ), ]
    
    ## flux .... but for now, magnitudes.
    yval <- c()
    yval_err <- c()
    for ( bb in band_order ){
        yval <- c( yval, my_source[1, which( my_cols == paste(bb,'_MAG_AUTO',sep='') ) ]  )
        yval_err <- c( yval_err, my_source[1, which( my_cols == paste(bb,'_MAGERR_AUTO',sep='') ) ] )
    }

    pdf( paste( source_id, '.pdf', sep='' ) )
    par( mar=c(5,5,2,2) )
    plot( band_wavelengths, yval, pch=15, xlim=c(2.5,0), ylim=c(0.9*min(yval),1.1*max(yval)), axes=FALSE, xlab='', ylab=''  )
    arrows(band_wavelengths, yval-yval_err, band_wavelengths, yval+yval_err, length=0.05, angle=90, code=3)
    legend( 'topleft', c(source_id,''), col='white', bty='n', cex=2 )
    axis( 1 )
    axis( 2 )
    mtext( 'Wavelength (microns)', side=1, line=3, cex=1.25 )
    mtext( 'Apparent magnitude', side=2, line=3, cex=1.25 )
    box( which='plot' )
    dev.off()

}

find_snr <- function( infile, snr=5, band_aper='MAG_AUTO', band_aper_err='MAGERR_AUTO' ){

    df <- read.table( infile, sep=',', header=TRUE, stringsAsFactors=FALSE )
    ## get the bands
    tmp_bands <- strsplit( colnames( df )[which( grepl( 'DELTA', colnames( df ) ) ) ], '_' )
    my_bands <- c()
    for ( ii in 1:length( tmp_bands ) ) if ( length( tmp_bands[[ii]] ) == 3 ) my_bands <- c( my_bands, tmp_bands[[ii]][1] )

    ## loop through the bands to find the SNR
    for ( my_b in my_bands ){
        my_snr <- 2.5 / log( 10 ) / df[ ,paste(my_b,band_aper_err, sep='_') ]        
        df[ , paste(my_b,'SNR',sep='_')] <- my_snr
    }
    
    ## write out the file
    write.table( df, file=gsub( '.csv', '_snr.csv', infile ), quote=FALSE, sep=',', row.names=FALSE )

    cat( 'File with SNR written to',gsub( '.csv', '_snr.csv', infile ), '\n' )

    return( df )

}




make_snr_cut_catalogues <- function( df, snr=5 ){

    tmp_bands <- strsplit( colnames( df )[which( grepl( 'DELTA', colnames( df ) ) ) ], '_' )
    my_bands <- c()
    for ( ii in 1:length( tmp_bands ) ) if ( length( tmp_bands[[ii]] ) == 3 ) my_bands <- c( my_bands, tmp_bands[[ii]][1] )

}



















