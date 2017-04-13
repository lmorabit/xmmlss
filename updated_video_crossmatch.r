## import necessary libraries
source('/vardy/leah/xmmlss/plotting_libraries.r')
source('/vardy/leah/xmmlss/multi_wavelength_comparison.r')
source('/vardy/leah/xmmlss/leah_helper_functions.r')

################# HOUSEKEEPING
mycols <- viridis( 5 )
stilts_exec <- '/vardy/leah/stilts/stilts'
master_cat <- '/vardy/leah/data/xmmlss/cats/CFHTLS-W1_2016-04-14_fullcat_errfix_MAGAUTO.fits'
snr_cat <- gsub( '.fits', '_SNR.csv', master_cat )
gs_class_cat <- 'GSCLASS.csv'
## which magnitude values to use
mag_col_to_use <- 'MAG_AUTO'
mag_err_col_to_use <- gsub( 'MAG', 'MAGERR', mag_col_to_use ) 
my_bands <- c('J','H','K','U','G','R','I','ZC','ZV','Y','YC')

make_footprint_maps <- FALSE  ## set to true if you want to make maps to see the coverage of each band. takes a while to run.

################# MAKE FOOTPRINT MAPS OF ALL BANDS
## make a catalogue with just the SNR columns, ID, and RA/DEC for mask-making purposes
if ( make_footprint_maps ){
    cat( 'Making footprint maps ... ' )
    snr_cols <- paste( c( 'ID', 'ALPHA_J2000', 'DELTA_J2000', paste( my_bands, '_SNR', sep='' ) ), collapse=' ' )
    ss <- paste( stilts_exec, " tpipe in=", master_cat, " out=", snr_cat, " omode=out ofmt=csv", ' cmd=\'keepcols "', snr_cols, '"', "\' ", sep='' )
    system( ss )
    ## make masks for areas with SNR >= 5
    snr_dat <- read.table( snr_cat, header=TRUE, stringsAsFactors=FALSE, sep=',' )
    snr_dat[ ,'HALOFLAG'] <- rep( 0, dim( snr_dat )[1] )
    ## get band names
    snr_cols <- colnames( snr_dat )[4:14]
    for ( snrc in snr_cols ){
        tmp_dat <- snr_dat[ which( snr_dat[ , snrc ] >= 5 ), ]
        my_mask <- create_footprint_mask( tmp_dat, cellsize=40, tolerance=1e-7, outfile=snrc, ra_hours=FALSE, method='2dhist', exclude_halo=FALSE, use_minmax=FALSE )
    } # end for
    cat( 'done.\n' )
} # end if
##################################


################# STAR GALAXY SEPARATION
## match the master cat with the star ID catalogue and remove where GS_CLASS==1 (i.e., stars)
## if files already exist, this uses the output
cat( 'Performing the star/galaxy sepration to remove stars ... \n' )

if ( !file.exists( 'stars.fits' ) ){
    ## make a star catalogue from the star/galaxy separation
    ss <- paste( stilts_exec, ' tpipe in=', gs_class_cat, ' ifmt=csv out=stars.fits cmd=\'select GS_CLASS==1\'', sep='' )
    system( ss )
}

## remove the stars
galaxy_cat <- gsub( '.fits', '_Galaxies.fits', master_cat )
if ( !file.exists( galaxy_cat ) ){
    ss <- paste( stilts_exec, ' tmatch2 in1=', master_cat, ' in2=stars.fits out=', galaxy_cat, ' omode=out matcher=exact values1=ID values2=SID join=1not2', sep='' )
    system( ss )
}

cat( ' ... done.\n' )
##################################

################# READ IN RADIO DATA AND PREPARE

cat( 'Preparing the VLA data!\n' )
## read in VLA data
## master radio file
vla_file <- '/vardy/leah/data/xmmlss/VLA/13B-308_080916_csv.txt'

snr_cutoff <- 5 ## self-explanatory
myoutfile <- paste( 'sum_VLA_SNR', snr_cutoff, '.csv', sep='' )
if ( !file.exists( myoutfile ) ) VLA_dat <- prepare_radio_data( vla_file, snr_cutoff=snr_cutoff, outfile=myoutfile ) else VLA_dat <- read.table( myoutfile, stringsAsFactors=FALSE, header=TRUE, sep=',' )

## apply the mask
VLA_dat <- apply_masked_regions( VLA_dat )
## write out the file
write.table( VLA_dat, file=gsub( 'sum', 'masked_sum', myoutfile ), row.names=FALSE, quote=FALSE, sep=',' )

## complete data will be VLA_dat
cat( 'There are', dim( VLA_dat )[1], 'radio sources within the VIDEO coverage.\n' )

## flag sources that are appropriate for doing LR matching
cat( 'Finding isolated, unresolved sources for estimating LR values ... \n' )

##-- use things that pybdsm have identified as 'S' type sources
single_index <- which( VLA_dat$S_Code == 'S' )
cat( 'There are', length( single_index ), 'S type sources.\n' )

##-- find unresolved sources
unresolved_index <- find_unresolved_index( VLA_dat$Total_flux, VLA_dat$Peak_flux, VLA_dat$e_Total_flux, VLA_dat$e_Peak_flux )
cat( 'There are', length( unresolved_index ), 'unresolved sources.\n' )

##-- find isolated sources
hdens <- density( VLA_dat$nearest_neighbour )
peak_hdens <- hdens$x[ which( hdens$y == max( hdens$y ) ) ]
cat( 'Peak in nearest neighbour distribution at', peak_hdens, 'arcsec.\n' )
my_cutoff <- round( peak_hdens/10 ) * 10
isolated_index <- which( VLA_dat$nearest_neighbour >= my_cutoff )
cat( 'There are', length( isolated_index ), 'isolated sources.\n' )

##-- combine the indices
i_list <- list( single_index, unresolved_index, isolated_index )
combined_index <- Reduce( intersect, i_list )

cat( 'There are', length( intersect( single_index, unresolved_index ) ), 'single, unresolved sources.\n' )
cat( 'There are', length( combined_index ), 'single, unresolved, isolated sources.\n' )

## test this
combined_index <- intersect( single_index, unresolved_index )


lr_flag <- rep( 0, dim( VLA_dat )[1] )
lr_flag[ combined_index ] <- 1 

##-- add to data frame
VLA_dat[ , 'LR_flag' ] <- lr_flag

##-- make a data frame with just those things for LR matching
my_radio_cat <- VLA_dat[ which( VLA_dat$LR_flag == 1 ), ]

################# CROSS-MATCHING!

##---- PRE-LOOP HOUSEKEEPING
## LR relevant things
LR_threshold <- 0.8
match_band <- c()
radio_ID <- c()
video_ID <- c()
lr_value <- c()
lr_rel <- c()
n_cont <- c()
separation <- c()
## band catalogue relevant things
basic_cols <- c( 'ID', 'ALPHA_J2000', 'DELTA_J2000', 'HALOFLAG', 'CONTEXT' )
info_cols <- c( 'FLAGS', 'CLASS_STAR' )


#my_bands <- c( 'K' )
my_bands <- c('J','H','K','U','G','R','ZC','Y')

cat( 'Beginning LR matching ... \n' )
for ( my_band in my_bands ){
    
	cat( ' ... Processing band', my_band, '\n' )
    ################# MAKE A BAND CATALOGUE AND APPLY MASK
    keep_bands <- paste( my_band, mag_col_to_use, sep='_' )
    keep_band_errs <- paste( my_band, mag_err_col_to_use, sep='_' )
    keep_info_cols <- paste( my_band, info_cols, sep='_' )
    keep_cols <- c( basic_cols, keep_bands, keep_band_errs, paste( my_band, '_SNR', sep='' ), paste( my_band, info_cols, sep='_' ) )

    band_cat <- gsub( '.fits', paste( '_', my_band, '.csv', sep='' ), galaxy_cat )
    if ( !file.exists( band_cat ) ){
        cat( my_band, 'band catalogue does not exist yet, creating ...\n' )
        ## get rid of things with MAG_AUTO = 99
        ## get rid of things with FLAGS > 4
        ss <- paste( stilts_exec, ' tpipe in=', galaxy_cat, ' out=', band_cat, ' omode=out ofmt=csv cmd=\'select ', my_band, '_FLAGS<5\' cmd=\'select ', my_band, '_MAG_AUTO<99\' cmd=\'keepcols "', paste( keep_cols, collapse=' ' ), '"', "\' ", sep='' )
        system( ss )
    } # end if

    ## read in data   
    cat( 'Reading in', my_band, 'band data ... ' )
    band_dat <- read.table( band_cat, header=TRUE, stringsAsFactors=FALSE, sep=',' )
    ## get the right magnitude columns
    mag_cols <- c( paste( my_band, mag_col_to_use, sep='_' ), paste( my_band, mag_err_col_to_use, sep='_' ) )
    cat( 'done.\n' )

    ## APPLY THE MASK
    masked_band_dat <- apply_masked_regions( band_dat )    
    ## get rid of lingering HALOFLAG (there's only 139 of them)
    masked_band_dat <- masked_band_dat[ which( masked_band_dat$HALOFLAG == 0 ), ]
    ## write out the table
    write.table( masked_band_dat, file=gsub( paste( '_', my_band, sep='' ), paste( '_masked_', my_band, sep='' ), band_cat ), row.names=FALSE, quote=FALSE, sep=',' )
    
    snr_cutoff <- 5
    valid_index <- which( masked_band_dat[ , paste( my_band, '_SNR', sep='' ) ] >= snr_cutoff )

    ##-- CATALOGUE FOR USE IN LR MATCHING
    video_dat <- masked_band_dat[ valid_index, ]

    ## plot this for later
    pdf( 'video_cat_coverage.pdf' )
    lplot( video_dat$ALPHA_J2000, video_dat$DELTA_J2000, pch='.', col='gray48', xlim=c(max(video_dat$ALPHA_J2000),min(video_dat$ALPHA_J2000)), x_lab='RA [deg]', y_lab='Dec [deg]' )
    points( VLA_dat$RA, VLA_dat$DEC, pch=23, col=mycols[2], bg=mycols[4], cex=0.4 )#, lwd=0.5 )
    dev.off()


    ###################################

    #VLA_dat[ , 'LR_flag' ] <- lr_flag

    n_radio_sources <- dim( my_radio_cat )[1]
    cat( 'Starting LR matching for', n_radio_sources, 'radio sources.\n' )

    ## POSITIONAL ERROR
    poserr <- calculate_positional_errors( my_radio_cat, beam_size=5 )
    r_max <- poserr$r_max
    sigma_pos <- poserr$sigma_pos 
    ## where sigma pos < 0.5 arcsec, make it = 0.5 arcsec
    sigma_pos[ which( sigma_pos < 0.5 ) ] <- 0.5

    ## MAKE THE MAGNITUDE BINS
    mag_bins <- make_mag_bins( video_dat, mag_cols )

    ## FIND n(m)
    video_area_arcsec <- get_video_area_arcsec()
    nm <- calculate_nm( video_dat, mag_cols, mag_bins, video_area_arcsec )

    ## FIND MATCHED MAGNITUDES
	mag_file <- paste( 'matched_magnitudes_', my_band, '_r', format( r_max, digits=3, nmall=2 ), '_SNR', snr_cutoff, '.txt', sep='' )
    mag_matches <- calculate_matched_mags( my_band, my_radio_cat, video_dat, r_max, filename=mag_file )
    match_magnitudes <- mag_matches$match_magnitudes
    match_radio_sources <- mag_matches$radio_ids

    ## FIND q0
    ## --- Ciliegi+ (2003)
	n_matches <- length( match_magnitudes )
	Q0 <- ( n_matches - sum( nm ) * pi * r_max^2. * n_radio_sources ) / n_radio_sources 
	cat( 'Ciliegi+ Q_0:',Q0,'\n' )

    radii <- seq( 1, 5, 0.2 ) ## don't go below the resoluion limit of the K-band data
	fleuren_Q0 <- find_Q0_fleuren( my_radio_cat, video_dat, radii )
	fQ0 <- fleuren_Q0$q_0 ## the second element is the error
    cat( 'Fleuren+ Q_0:', fQ0,'\n' )

    ## FIND q(m) -- the expected distribution of the true counterparts
    ## -- find total(m)
	tmhist <- calculate_hist( match_magnitudes, mag_bins )
	total_m <- c()
	for ( ii in 1:length( tmhist$counts ) ) total_m <- c( total_m, sum( tmhist$counts[1:ii] ) )

	background <- nm * n_radio_sources * pi * r_max^2
	real_m <- total_m - background

	qm <- real_m / sum( real_m ) * fQ0

	## FIND RATIO OF q(m)/n(m)
	qm_nm <- qm / nm

	## plot some things so far
	pdf( paste( 'LR_values_', my_band, '.pdf',sep='' ) )
	m <- rbind( c(1), c(2), c(3) )
	layout(m, heights=c(1.15,1,1.3))
    ## plot 1
	lplot( tmhist$mids, log10( total_m ), xlim=c(14,22.2), y_lab='log(N(counterparts))', xaxon=FALSE, type='s', lwd=2, col='gray', ylim=c(-1,3.5), margins=c(0,5,2,2) )
	lines( tmhist$mids, log10( real_m ), lty=2, lwd=2, type='s' )
	lines( tmhist$mids, log10( background ), lty=3, lwd=2, type='s' )
	legend( 'topleft', c('Total','Real','Background'), col=c('gray','black','black'), lwd=2, lty=c(1,2,3), bty='n' )
    ## plot 2
	lplot( tmhist$mids, log10( qm_nm ), xlim=c(14,22.2), y_lab='log(P(m)=q(m)/n(m))', xaxon=FALSE, type='s', lwd=2, ylim=c(1,3), margins=c(0,5,0,2) )
    ## plot 3
	lplot( tmhist$mids, log10( qm ), xlim=c(14,22.2), x_lab='Ks mag', y_lab='log(q(m))', type='s', lwd=2, margins=c(5,5,0,2) )
	dev.off()

    ## FIND LR AND RELIABILITY VALUES
    ## do this for the entire catalogue
    my_radio_cat <- VLA_dat

    n_radio_sources <- dim(my_radio_cat)[1]
	pb <- txtProgressBar( min=0, max=n_radio_sources, style=3 )
	count <- 0 
	for ( ii in 1:n_radio_sources ){

		## CALCULATE f(r)
		## -- find distances 
		distances <- cenang( my_radio_cat$RA[ii], my_radio_cat$DEC[ii], video_dat$ALPHA_J2000, video_dat$DELTA_J2000 ) * 60. * 60.
		candidate_index <- which( distances <= r_max ) 
	    if ( length( candidate_index ) > 0 ){
            ## get the data of the candidates
		    tmp_dat <- video_dat[candidate_index,]
            ## get their magnitudes
			band_magnitudes <- tmp_dat[,mag_cols[1]]
            ## calculate the radial probability distribution for the candidates
			f_r <- 1 / ( 2 * pi * sigma_pos[ii]^2 ) * exp( - distances[candidate_index]^2 / ( 2 * sigma_pos[ii]^2 ) )
			## loop through candidates to calculate the LR
			LR <- c()
			for ( jj in 1:length(candidate_index) ){	
				## find the qm and nm for the magnitude bin of the source
				qm_nm_jj <- qm_nm[ max( which( mag_bins <= band_magnitudes[jj]) ) ]
				LR <- c( LR, qm_nm_jj * f_r[jj] )
			}

		LR_reliability <- LR / ( sum( LR ) + 1 - fQ0 )
	
		radio_ID <- c( radio_ID, rep( my_radio_cat$Source_id[ii], length(candidate_index)) )
		video_ID <- c( video_ID, tmp_dat$ID )
		lr_value <- c( lr_value, LR )
		lr_rel <- c( lr_rel, LR_reliability )
		lr_cont <- sum( 1 - LR_reliability[ which( LR_reliability > LR_threshold ) ] )
		n_cont <- c( n_cont, rep( lr_cont, length(candidate_index) ) )
        separation <- c( separation, distances[candidate_index] )
		count <- count + length( candidate_index )

	}
	match_band <- c( match_band, rep( my_band, length( candidate_index ) ) )
	setTxtProgressBar( pb, ii )   
	}
	close( pb )
	cat( '\n' )

    cat( ' ... done processing band.\n' )

} ## end for (my_bands)

cat( 'Writing final matches to a file.\n' )
final_matches <- data.frame( radio_ID, match_band, video_ID, lr_value, lr_rel, n_cont, separation, stringsAsFactors=FALSE )
write.table( final_matches, file='final_matches.csv', quote=FALSE, row.names=FALSE, sep=',' )
cat( 'done.\n' )







