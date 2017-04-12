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
