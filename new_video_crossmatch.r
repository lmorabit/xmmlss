## import necessary libraries
source('/vardy/leah/xmmlss/plotting_libraries.r')
source('/vardy/leah/xmmlss/multi_wavelength_comparison.r')
source('/vardy/leah/xmmlss/video_data.r')
source('/vardy/leah/xmmlss/leah_helper_functions.r')

################# HOUSEKEEPING
mycols <- viridis( 5 )
stilts_exec <- '/vardy/leah/stilts/stilts'
master_cat <- 'cats/CFHTLS-W1_2016-04-14_fullcat_errfix_MAGAUTO.fits'
snr_cat <- gsub( '.fits', '_SNR.csv', master_cat )
gs_class_cat <- 'GSCLASS.csv'
## which magnitude values to use
mag_col_to_use <- 'MAG_AUTO'
mag_err_col_to_use <- gsub( 'MAG', 'MAGERR', mag_col_to_use ) 
my_bands <- c('J','H','K','U','G','R','I','ZC','ZV','Y','YC')

################# MAKE FOOTPRINT MAPS OF ALL BANDS
## make a catalogue with just the SNR columns, ID, and RA/DEC for mask-making purposes
make_footprint_maps <- FALSE
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
## match the master cat with the star ID catalogue and remove where GS_CLASS==1

cat( 'Using the star/galaxy sepration to remove stars ... \n' )

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

## OLD STUFF
## keep only things with HALOFLAG==0 
#galaxy_nohalo_cat <- gsub( '.fits', '_HALOFLAG0.fits', galaxy_cat )
#ss <- paste( stilts_exec, ' tpipe in=', galaxy_cat, ' out=', galaxy_nohalo_cat, ' omode=out cmd=\'select HALOFLAG==0\'', sep='' )
#system( ss )

##################################

################# MAKE A K-BAND CATALOGUE

cat( 'Making a K-band catalogue to use as a mask ... \n' )

my_bands <- c( 'K' )

## column information
basic_cols <- c( 'ID', 'ALPHA_J2000', 'DELTA_J2000', 'HALOFLAG', 'CONTEXT' )
info_cols <- c( 'FLAGS', 'CLASS_STAR' )

for ( my_band in my_bands ){
    keep_bands <- paste( my_band, mag_col_to_use, sep='_' )
    keep_band_errs <- paste( my_band, mag_err_col_to_use, sep='_' )
    keep_info_cols <- paste( my_band, info_cols, sep='_' )
    keep_cols <- c( basic_cols, keep_bands, keep_band_errs, paste( my_band, '_SNR', sep='' ), paste( my_band, info_cols, sep='_' ) )

    band_cat <- gsub( '.fits', paste( '_', my_band, '.csv', sep='' ), galaxy_cat )
    if ( !file.exists( band_cat ) ){
        ## get rid of things with MAG_AUTO = 99
        ## get rid of things with FLAGS > 4
        ss <- paste( stilts_exec, ' tpipe in=', galaxy_cat, ' out=', band_cat, ' omode=out ofmt=csv cmd=\'select ', my_band, '_FLAGS<5\' cmd=\'select ', my_band, '_MAG_AUTO<99\' cmd=\'keepcols "', paste( keep_cols, collapse=' ' ), '"', "\' ", sep='' )
        system( ss )
    } # end if
} # end for

cat( ' ... done. Reading in data ... ' )
## read in data -- just K-band at the moment!
kband_dat <- read.table( band_cat, header=TRUE, stringsAsFactors=FALSE, sep=',' )
cat( 'done.\n' )

##################################

################# MAKE A MASK FROM THE K-BAND CATALOGUE

## make a mask
cat( 'Making a mask to use for the radio data ... ' )
radio_mask <- create_footprint_mask( kband_dat, cellsize=40, tolerance=1e-7, outfile='Kband', ra_hours=FALSE, method='2dhist', exclude_halo=TRUE, use_minmax=FALSE )
fivesig_index <- which( kband_dat$K_SNR >= 5 )
radio_mask_fivesig <- create_footprint_mask( kband_dat[fivesig_index,], cellsize=40, tolerance=1e-7, outfile='Kband_5sig', ra_hours=FALSE, method='2dhist', exclude_halo=TRUE, use_minmax=FALSE )

video_area <- data.frame( min( radio_mask$x.breaks ), min( radio_mask$y.breaks ), max( radio_mask$x.breaks ), max( radio_mask$y.breaks ) )
colnames( video_area ) <- c( 'xl', 'yl', 'xr', 'yr' )
cat( 'done.\n' )

## apply the mask to the video data (as a check)
cat( 'Applying mask to K-band data ... ' )
masked_kband <- apply_mask( kband_dat, radio_mask_fivesig, filestem='kband_5sig' )
cat( 'done.\n' )

##################################

################# READ IN RADIO DATA AND PREPARE

cat( 'Preparing the VLA data!\n' )
## read in VLA data
## master radio file
vla_file <- 'VLA/13B-308_080916_csv.txt'

snr_cutoff <- 5
myoutfile <- paste( 'sum_VLA_',paste( format(video_area,digits=3,nsmall=3), collapse='_' ), '_SNR', snr_cutoff, '.csv', sep='' )
if ( !file.exists( myoutfile ) ) VLA_dat <- prepare_radio_data( vla_file, video_area, snr_cutoff=snr_cutoff, outfile=myoutfile ) else VLA_dat <- read.table( myoutfile, stringsAsFactors=FALSE, header=TRUE, sep=',' )

cat( 'Applying the mask to the radio data ... ' )
## use fivesig mask???
if ( !file.exists('VLA_fivesig_masked.csv') ) masked_VLA_fivesig <- apply_mask( VLA_dat, radio_mask_fivesig, filestem='VLA_fivesig' ) else masked_VLA_fivesig <- read.table( 'VLA_fivesig_masked.csv', stringsAsFactors=FALSE, header=TRUE, sep=',' )
## or normal mask???
if ( !file.exists('VLA_masked.csv') ) masked_VLA <- apply_mask( VLA_dat, radio_mask, filestem='VLA' ) else masked_VLA <- read.table( 'VLA_masked.csv', stringsAsFactors=FALSE, header=TRUE, sep=',' )
cat( 'done.\n' )

###################################

################# FIND ISOLATED, UNRESOLVED SOURCES FOR CROSS-MATCHING

cat( 'Finding isolated, unresolved sources for estimating LR values ... \n' )
my_radio_cat <- masked_VLA_fivesig

#h <- hist( my_radio_cat$nearest_neighbour, breaks=60, prob=TRUE, plot=FALSE )
h <- hist( my_radio_cat$nearest_neighbour, breaks=60, plot=FALSE )
hdens <- density( my_radio_cat$nearest_neighbour )
peak_hdens <- hdens$x[ which( hdens$y == max( hdens$y ) ) ]
cat( 'Peak in nearest neighbour distribution at', peak_hdens, 'arcsec.\n' )
my_cutoff <- round( peak_hdens/10 ) * 10
h_index <- which( my_radio_cat$nearest_neighbour >= my_cutoff )
##-- find unresolved
unresolved_index <- find_unresolved_index( my_radio_cat$Total_flux, my_radio_cat$Peak_flux, my_radio_cat$e_Total_flux, my_radio_cat$e_Peak_flux )
##-- and take only 'S' type sources
single_index <- which( my_radio_cat$S_Code == 'S' )
i_list <- list( h_index, unresolved_index, single_index )
isolated_index <- Reduce( intersect, i_list )

cat( 'There are', length( isolated_index ), 'isolated, unresolved sources within the K-band mask.\n' )

radio_cat_isolated <- my_radio_cat[ isolated_index, ]

##################################

################# CROSS-MATCHING!

##---- PRE-LOOP HOUSEKEEPING
LR_threshold <- 0.8
match_band <- c()
radio_ID <- c()
video_ID <- c()
lr_value <- c()
lr_rel <- c()
n_cont <- c()

my_bands <- c( 'K' )

cat( 'Beginning LR matching ... \n' )
for ( my_band in my_bands ){

	cat( ' ... Processing band', my_band, '\n' )
    n_radio_sources <- dim( my_radio_cat )[1]
    n_radio_sources_isolated <- length( isolated_index )
    cat( 100*n_radio_sources_isolated/n_radio_sources, 'percent of sources are isolated/unresolved\n' )

    ## read the catalogue
    band_cat <- Sys.glob( paste('cats/*Galaxies_',my_band,'.csv',sep='' ) )
    band_dat <- read.table( band_cat, header=TRUE, stringsAsFactors=FALSE, sep=',' )
    ## and use only things where HALOFLAG==0
    band_dat <- band_dat[ which( band_dat$HALOFLAG == 0 ), ]

    ## get the right magnitude columns
    mag_cols <- c( paste( my_band, mag_col_to_use, sep='_' ), paste( my_band, mag_err_col_to_use, sep='_' ) )

    ## apply five-sigma mask
    masked_band_dat <- apply_mask( band_dat, radio_mask_fivesig, filestem=my_band )
    ## and take only five-sigma values
    masked_band_dat_fivesig <- masked_band_dat[ which( masked_band_dat[ , mag_cols[2] ] <= 0.217 ), ]


    ## POSITIONAL ERROR
    poserr <- calculate_positional_errors( my_radio_cat, beam_size=5 )
    r_max <- poserr$r_max
    sigma_pos <- poserr$sigma_pos 

    ## MAKE THE MAGNITUDE BINS
    mag_bins <- make_mag_bins( masked_band_dat_fivesig, mag_cols )

    ## FIND n(m)
    video_area_arcsec <- length( which( radio_mask_fivesig$counts == 1 ) ) * 40 * 40 
    nm <- calculate_nm( masked_band_dat_fivesig, mag_cols, mag_bins, video_area_arcsec )

    ## FIND MATCHED MAGNITUDES
    mag_matches <- calculate_matched_mags( my_band, my_radio_cat, masked_band_dat_fivesig, r_max )
    match_magnitudes <- mag_matches$match_magnitudes 
    match_radio_sources <- mag_matches$radio_id
    ## match_magnitudes_isolated 
    match_magnitudes_isolated <- match_magnitudes[ which( match_radio_sources %in% my_radio_cat$Source_id[isolated_index] ) ]

    ## FIND q0
    ## --- Ciliegi+ (2003)
	n_matches <- length( match_magnitudes )
	n_matches_isolated <- length( match_magnitudes_isolated )
	Q0 <- ( n_matches - sum( nm * pi * r_max^2. * n_radio_sources ) ) / n_radio_sources 
	Q0_isolated <- ( n_matches_isolated - sum( nm * pi * r_max^2. * n_radio_sources_isolated ) ) / n_radio_sources_isolated
	cat( 'Ciliegi+ Q_0:',Q0_isolated,'\n' )

    radii <- seq( 1, 5, 0.2 ) ## don't go below the resoluion limit of the K-band data
	fleuren_Q0 <- find_Q0_fleuren( my_radio_cat[isolated_index,], radio_mask_fivesig, masked_band_dat_fivesig, radii )
	fQ0 <- fleuren_Q0$q_0 ## the second element is the error
    cat( 'Fleuren+ Q_0:', fQ0,'\n' )

    ## FIND q(m) -- the expected distribution of the true counterparts
    ## -- find total(m)
	tmhist <- calculate_hist( match_magnitudes_isolated, mag_bins )
	total_m <- c()
	for ( ii in 1:length( tmhist$counts ) ) total_m <- c( total_m, sum( tmhist$counts[1:ii] ) )

	background <- nm * n_radio_sources_isolated * pi * r_max^2
	real_m <- total_m - background

	qm <- real_m / sum( real_m ) * fQ0

	## FIND RATIO OF q(m)/n(m)
	qm_nm <- qm / nm

	## plot some things so far
	pdf( paste( 'LR_values_', my_band, '.pdf',sep='' ) )
	m <- rbind( c(1), c(2), c(3) )
	layout(m)
    ## plot 1
	lplot( tmhist$mids, log10( total_m ), xlim=c(13,22.2), y_lab='log(N(counterparts))', type='s', lwd=2, col='gray', ylim=c(-1,3.5) )
	lines( tmhist$mids, log10( real_m ), lty=2, lwd=2, type='s' )
	lines( tmhist$mids, log10( background ), lty=3, lwd=2, type='s' )
	legend( 'topleft', c('Total','Real','Background'), col=c('gray','black','black'), lwd=2, lty=c(1,2,3), bty='n' )
    ## plot 2
	lplot( tmhist$mids, log10( qm_nm ), xlim=c(13,22.2), y_lab='log(P(m)=q(m)/n(m))', type='s', lwd=2, ylim=c(1,3) )
    ## plot 3
	lplot( tmhist$mids, log10( qm ), xlim=c(13,22.2), x_lab='Ks mag', y_lab='log(q(m))', type='s', lwd=2 )
	dev.off()

    ## FIND RELIABILITY VALUES
    n_radio_sources <- dim(my_radio_cat)[1]
	pb <- txtProgressBar( min=0, max=n_radio_sources, style=3 )
	count <- 0 
	for ( ii in 1:n_radio_sources ){

		## CALCULATE f(r)
		## -- find distances 
		distances <- cenang( my_radio_cat$RA[ii], my_radio_cat$DEC[ii], masked_band_dat_fivesig$ALPHA_J2000, masked_band_dat_fivesig$DELTA_J2000 ) * 60. * 60.
		candidate_index <- which( distances <= r_max ) 
	    if ( length( candidate_index ) > 0 ){
            ## get the data of the candidates
		    tmp_dat <- masked_band_dat_fivesig[candidate_index,]
            ## get their magnitudes
			band_magnitudes <- tmp_dat[,mag_cols[1]]
            ## calculate the radial probability distribution for the candidates
			f_r <- 1 / ( 2 * pi * sigma_pos[ii]^2 ) * exp( distances[candidate_index]^2 / ( 2 * sigma_pos[ii]^2 ) )

			## loop through candidates to calculate the LR
			LR <- c()
			for ( jj in 1:length(candidate_index) ){	
				## find the qm and nm for the magnitude bin of the source
				#qm_nm_jj <- qm_nm[ min( which( band_magnitudes[jj] > mag_bins ) ) ]
				qm_nm_jj <- qm_nm[ max( which( mag_bins <= band_magnitudes[jj] ) ) ]
				LR <- c( LR, qm_nm_jj * f_r[jj] )
			}

		LR_reliability <- LR / ( sum( LR ) + 1 - fQ0 )
	
		radio_ID <- c( radio_ID, rep( my_radio_cat$Source_id[ii], length(candidate_index)) )
		video_ID <- c( video_ID, tmp_dat$ID )
		lr_value <- c( lr_value, LR )
		lr_rel <- c( lr_rel, LR_reliability )
		lr_cont <- sum( 1 - LR_reliability[ which( LR_reliability > LR_threshold ) ] )
		n_cont <- c( n_cont, rep( lr_cont, length(candidate_index) ) )
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
final_matches <- data.frame( radio_ID, match_band, video_ID, lr_value, lr_rel, n_cont, stringsAsFactors=FALSE )
write.table( final_matches, file='final_matches.csv', quote=FALSE, row.names=FALSE, sep=',' )
cat( 'done.\n' )







