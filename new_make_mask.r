create_footprint_mask <- function( df, cellsize=60, tolerance=1e-7, outfile='mask', ra_hours=FALSE, exclude_halo=TRUE ){

    df_cols <- colnames( df )

    if ( 'RA' %in% df_cols ){
        ra_name <- 'RA'
        dec_name <- 'DEC'
    } else {
        ra_name <- 'ALPHA_J2000'
        dec_name <- 'DELTA_J2000'
    }

    ## get RA and DEC of everything
    ra <- df[ , ra_name ]
    if ( ra_hours ) ra <- ra*15
    dec <- df[ , dec_name ]

    ## get SNR information
    snr_col <- df_cols[ which( grepl( 'SNR', df_cols ) ) ]
    snr_index <- ( df[,snr_col] >= 5 )

    ## get haloflag information
    if ( exclude_halo ) hf_index <- ( df[ ,'HALOFLAG' ] == 1 )

    ## first define the area to use CFHTLS-W1_2016-04-14_fullcat_errfix_MAGAUTO_Galaxies_K.csv
    xmin <- min( ra ) - ( cellsize / 60. / 60. )
    xmax <- max( ra ) + ( cellsize / 60. / 60. )

    ymin <- min( dec ) - ( cellsize / 60. / 60. )
    ymax <- max( dec ) + ( cellsize / 60. / 60. )

    xnsteps <- ceiling( ( xmax - xmin ) / ( cellsize / 60. / 60. ) )
    xsteps <- ( xmax - xmin ) / xnsteps 
    xbins <- seq( xmin, xmax, xsteps )

    ynsteps <- ceiling( ( ymax - ymin ) / ( cellsize / 60. / 60. ) )
    ysteps <- ( ymax - ymin ) / ynsteps
    ybins <- seq( ymin, ymax, ysteps )

    ## make a matrix that is of size c(xbins-1,ybins-1)
    counts <- matrix( nrow=length(xbins)-1, ncol=length(ybins)-1 )

    ## find the counts in each bin
	pb <- txtProgressBar( min=2, max=dim(counts)[1], style=3 )
    if ( exclude_halo ){
        for ( ii in 1:dim(counts)[1] ){
            for ( jj in 1:dim(counts)[2] ){
                ra_index <- ( ra >= xbins[ii] & ra < xbins[ii+1] )
                dec_index <- ( dec >= ybins[jj] & dec < ybins[jj+1] )
                counts[ii,jj] <- sum( (ra_index + dec_index + snr_index - hf_index )/3 == 1 ) 
            } #endfor jj
        	setTxtProgressBar( pb, ii )   
        } #endfor ii
    } else {
        for ( ii in 1:dim(counts)[1] ){
            for ( jj in 1:dim(counts)[2] ){
                ra_index <- ( ra >= xbins[ii] & ra < xbins[ii+1] )
                dec_index <- ( dec >= ybins[jj] & dec < ybins[jj+1] )
                counts[ii,jj] <- sum( (ra_index + dec_index + snr_index )/3 == 1 ) 
            } #endfor jj
        	setTxtProgressBar( pb, ii )   
        } #endfor ii
    }
	close( pb )
	cat( '\n' )

    counts[ which( counts > 0 ) ] <- 1

    k <- list( x.breaks=xbins, y.breaks=ybins, counts=counts )
    kim <- raster( list( x=xbins, y=ybins, z=counts ) )

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

    return( k )

}
