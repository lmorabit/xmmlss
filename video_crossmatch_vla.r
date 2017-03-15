## PROCESS VIDEO RADIO DATA
##  written by Leah K. Morabito
##  started on 26 Jan 2017
##  last updated ...
##  contact: lmorabit@gmail.com

## starting from /data/xmmlss/
script_dir <- '/vardy/leah/xmmlss/'

source( paste( script_dir, 'video_data.r', sep='' ))
source( paste( script_dir, 'multi_wavelength_comparison.r', sep='' ))
source( paste( script_dir, 'leah_helper_functions.r', sep='' ))

stilts_exec <- '/vardy/leah/stilts/stilts'

## VIDEO master catalogue
master_cat <- '/vardy/videouser/v1.3-20160414/cats/CFHTLS-W1_2016-04-14_fullcat_errfix.fits'

## STAR-GALAXY SEPARATION -----------------------------------------------------
if ( !file.exists( 'sg_separation_values.txt' ) ) source( paste( script_dir, 'star_galaxy_separation.r', sep='' ))
## read in the values for separating stars and galaxies
sg_vals <- read.table( 'sg_separation_values.txt', stringsAsFactors=FALSE, header=TRUE )

## find the g-i and J-K values for the whole catalogue
master_cat_sg <- gsub( '.fits', '_SGsep.csv', master_cat )
## if on vardy change the name
if ( strsplit( master_cat_sg, '/' )[[1]][2] == 'vardy' ) master_cat_sg <- gsub( 'videouser/v1.3-20160414/cats', 'leah/data/xmmlss/cats', master_cat_sg )

if ( !file.exists( master_cat_sg ) ){
	ss <- paste( stilts_exec, ' tpipe in=', master_cat, ' out=', master_cat_sg, ' omode=out ofmt=csv cmd=\'addcol G_I G_MAG_APER_2-I_MAG_APER_2\' cmd=\'addcol J_K J_MAG_APER_2-K_MAG_APER_2\' cmd=\'keepcols "ID G_I J_K"\'', sep='' )
	system( ss )

	## read in the values and locate the portion of space classified as 'galaxies'
	sg_dat <- read.table( master_cat_sg, stringsAsFactors=FALSE, header=TRUE, sep=',' )
	## calculate quad values
	locus_values <- sg_vals$c2 + sg_vals$b2 * sg_dat$G_I + sg_vals$a2 * sg_dat$G_I^2
	locus_values[ which( sg_dat$G_I < 0.4 ) ] <- sg_vals$c1
	locus_values[ which( sg_dat$G_I > 1.9 ) ] <- sg_vals$c3
	galaxy_index <- which( sg_dat$J_K > locus_values )

	## the final IDs of 'galaxies'
	GS_CLASS <- rep( 1, length( sg_dat$ID ) )
	GS_CLASS[ galaxy_index ] <- 0
	df <- data.frame( sg_dat$ID, GS_CLASS )
	colnames( df ) <- c( 'GAL_ID', 'GS_CLASS' )

	## write out a catalogue with these IDs that can be cross-matched with the original catalogue
	write.table( df, file='GAL_IDs.csv', row.names=FALSE, quote=FALSE, sep=',' )

	master_cat_out <-  gsub( '.fits', '_GSclass.fits', master_cat )
	## if on vardy
	if ( strsplit( my_cat_name, '/' )[[1]][2] == 'vardy' ) master_cat_out <- gsub( 'videouser/v1.3-20160414/cats', 'leah/data/xmmlss/cats', master_cat_out )

	ss <- paste( stilts_exec, ' tmatch2 in1=GAL_IDs.csv ifmt1=csv in2=', master_cat, ' out=', master_cat_out, ' matcher=exact values1=GAL_ID values2=ID join=1and2', sep='' )
	system( ss )
}

master_cat_out <-  gsub( '.fits', '_GSclass.fits', master_cat )
## if on vardy
if ( strsplit( master_cat_out, '/' )[[1]][2] == 'vardy' ) master_cat_out <- gsub( 'videouser/v1.3-20160414/cats', 'leah/data/xmmlss/cats', master_cat_out )
master_cat <- master_cat_out
 


## MAKE A MASK FROM K-BAND -----------------------------------------------------
## STEP 2: MATCH BAND WITH VLA DATA
## 2.1: make a mask from K-band data
## 2.2: apply to radio data
## 2.3: look for cross matches

my_band <- 'K'
mag_cut <- 23.5
#mag_col <- 'APER_2' 
mag_col <- 'AUTO'

## set up output files, etc.
basic_cols <- c( 'ID', 'ALPHA_J2000', 'DELTA_J2000', 'HALOFLAG', 'CONTEXT', 'GS_CLASS' )
tmp_keep <- c( '_FLAGS', '_CLASS_STAR', '_DET_FLAG' )
mag_cols <- paste( my_band, c( '_MAG_', '_MAGERR_' ), mag_col, sep='' )
keep_cols <- c( basic_cols, mag_cols, 'SNR', paste( my_band, tmp_keep, sep='') )
my_cat_name <- gsub( '.fits', paste( '_', my_band, '_', mag_col, '_', mag_cut, '.csv', sep='' ), master_cat )

## this stilts command reads in the master catalogue and outputs a catalogue keeping only keep_cols
## where the magnitude cut has been applied, and snr is calculated and added to the table. 
ss <- paste( stilts_exec, " tpipe in=", master_cat, " out=", my_cat_name, " omode=out ofmt=csv cmd=\'select ", mag_cols[1], "<=", mag_cut, "\' cmd=\'addcol SNR 1.08574/", mag_cols[2], '\' cmd=\'keepcols "', paste( keep_cols, collapse=' ' ), '"', "\' ", sep='' )
system( ss )

## read in the data to make a mask
cv <- read.table( my_cat_name, stringsAsFactors=FALSE, header=TRUE, sep=',' )
## remove stars
galaxy_index <- which( cv$GS_CLASS == 0 )
cv <- cv[ galaxy_index, ]

## make a mask for the radio data avoiding haloflag areas
radio_mask <- create_footprint_mask( cv, 'ALPHA_J2000', 'DELTA_J2000', cellsize=60, outfile=my_cat_name, method='2dhist', exclude_halo=TRUE )

## APPLY MASK TO RADIO DATA -----------------------------------------------------

## VLA 
## start by reading in the catalogue and plotting some things
vla_file <- 'VLA/13B-308_080916_csv.txt'
VLA <- read_VLA_data( vla_file )
plot_sky_distribution( parm1='VLA', ra1=VLA$RA, de1=VLA$DEC )
plot_distributions( VLA, plotfile='VLA_distributions.pdf' )

## sum components by Source_id and replot things
if ( !file.exists('sum_VLA.csv') ) sum_VLA <- sum_components( VLA ) else sum_VLA <- read.table( 'sum_VLA.csv', stringsAsFactors=FALSE, header=TRUE, sep=',' )
plot_sky_distribution( parm1='sum_VLA', ra1=VLA$RA, de1=VLA$DEC )
plot_distributions( sum_VLA, plotfile='sum_VLA_distributions.pdf' )

## find nearest neighbours
nn <- find_nearest_neighbours( sum_VLA )
sum_VLA[ , 'nearest_neighbour' ] <- nn

## apply the mask from the VIDEO data
if ( !file.exists('VLA_masked.csv') ) masked_VLA <- apply_mask( VLA, radio_mask, filestem='VLA' ) else masked_VLA <- read.table( 'VLA_masked.csv', stringsAsFactors=FALSE, header=TRUE, sep=',' )
if ( !file.exists( 'sum_VLA_masked.csv' ) ) masked_sum_VLA <- apply_mask( sum_VLA, radio_mask, filestem='sum_VLA' ) else masked_sum_VLA <- read.table( 'sum_VLA_masked.csv', stringsAsFactors=FALSE, header=TRUE, sep=',' )
plot_sky_distribution( parm1='sum_VLA', ra1=sum_VLA$RA, de1=sum_VLA$DEC, parm2='Masked', ra2=masked_sum_VLA$RA, de2=masked_sum_VLA$DEC )

plot_radio_nearest_neighbours( VLA, plotfile='Nearest_neighbours_VLA.pdf' )
#plot_radio_nearest_neighbours( sum_VLA, plotfile='Nearest_neighbours_sum_VLA.pdf' )

## CROSS-MATCH WITH BANDS -----------------------------------------------------
LR_threshold <- 0.8

## select only the galaxies in the master catalogue
ss <- paste( stilts_exec, ' tpipe in=', master_cat, ' out=', gsub( 'GSclass', 'galaxies', master_cat ), ' cmd=\'select GS_CLASS==0\'',sep='' )
if ( !file.exists(gsub( 'GSclass', 'galaxies', master_cat )) ) system( ss )
galaxies_cat <- gsub( 'GSclass', 'galaxies', master_cat )


#my_bands <- c('J','H','K','U','G','R','I','ZC','ZV','Y','YC')
my_bands <- c('J','H','K','Y')
#my_bands <- c('K')

match_band <- c()
radio_ID <- c()
video_ID <- c()
lr_value <- c()
lr_rel <- c()
n_cont <- c()


## find isolated sources
my_radio_cat <- masked_sum_VLA
h <- hist( my_radio_cat$nearest_neighbour, breaks=60, prob=TRUE, plot=FALSE )
hdens <- density( my_radio_cat$nearest_neighbour )
peak_hdens <- hdens$x[ which( hdens$y == max( hdens$y ) ) ]
cat( 'Peak in nearest neighbour distribution at', peak_hdens, 'arcsec.\n' )
my_cutoff <- round( peak_hdens/10 ) * 10

my_radio_cat <- my_radio_cat[ which( my_radio_cat$nearest_neighbour >= my_cutoff ), ]


## select only 'single' sources?
#my_radio_cat <- masked_sum_VLA[ which( masked_sum_VLA$S_Code == 'S' ), ]
## get the unresolved things
unresolved_index <- find_unresolved_index( my_radio_cat$Total_flux, my_radio_cat$Peak_flux, my_radio_cat$e_Total_flux, my_radio_cat$e_Peak_flux )
my_radio_cat <- my_radio_cat[ unresolved_index, ]

#my_radio_cat <- masked_sum_VLA


for ( my_band in my_bands ){

	cat( 'Processing band', my_band, '\n' )

	## make a catalogue for the band
	band_cat <- gsub( '.fits', paste( '_', my_band, '.csv', sep='' ), galaxies_cat )
	mag_cols <- paste( my_band, c( '_MAG_', '_MAGERR_' ), mag_col, sep='' )
	keep_cols <- c( basic_cols, mag_cols, 'SNR', paste( my_band, tmp_keep, sep='') )
	ss <- paste( stilts_exec, ' tpipe in=', galaxies_cat, ' omode=out out=', band_cat, ' ofmt=csv cmd=\'addcol SNR 1.08574/', mag_cols[2], '\' cmd=\'select SNR>=5\' cmd=\'select HALOFLAG==0\' cmd=\'keepcols "', paste( keep_cols, collapse=' ' ), '"', "\' ", sep='' )
	system( ss )

	## read in the catalogue
	band_dat <- read.table( band_cat, stringsAsFactors=FALSE, header=TRUE, sep=',' )

	## get area -- use the same for every band
	area_degrees <- ( max( band_dat$ALPHA_J2000 ) - min( band_dat$ALPHA_J2000 ) ) * ( max( band_dat$DELTA_J2000 ) - min( band_dat$DELTA_J2000 ) )  
	## 5.9 square degrees! woo - nice check
	area_arcsec <- area_degrees * 3600.^2
	## that's a lot of arcseconds ...

	## CALCULATE POSITIONAL INFORMATION FOR THE RADIO
	## -- positional errors: condon
	#pos_noise <- sqrt( my_radio_cat$e_RA^2. + my_radio_cat$e_DEC^2. )
	## -- positional errors: ivison
	beam_size <- 5.
	pos_noise <- 0.6 * beam_size / ( my_radio_cat$Total_flux / my_radio_cat$e_Total_flux )
	## add in a calibration errors	
	sigma_pos <- sqrt( 0.1^2. + pos_noise^2. )
	r_max <- 5 * max( sqrt( sigma_pos^2 ) )
	cat( 'r_max is', r_max, 'arcsec.\n' )

	## CALCULATE nm -- density of background sources
	## -- find source counts in bins
	## numbers ~ 1e-4
	mag_bins <- seq( floor(min( band_dat[ , mag_cols[1]] )), ceiling(max( band_dat[,mag_cols[1]] )), 0.4 )
	## check that it spans the entire range
	if ( ceiling(max( band_dat[,mag_cols[1]] )) > max( mag_bins ) ) mag_bins <- c( mag_bins, max(mag_bins)+0.4 )
	if ( floor(min( band_dat[ , mag_cols[1]] )) < min( mag_bins ) ) mag_bins <- c( min(mag_bins)-0.4, mag_bins )
	nmhist <- calculate_hist( band_dat[,mag_cols[1]], mag_bins )
	nm <- nmhist$counts / area_arcsec
	
	## CALCULATE qm -- expected distribution of counter parts

	##-- first find magnitudes of all sources within r_max of all radio positions
	n_radio_sources <- dim(my_radio_cat)[1]
	matched_mags <- paste( 'matched_magnitudes_', my_band, '_r', r_max, '.txt', sep='' )
	if ( !file.exists( matched_mags ) ){
		match_magnitudes <- c()
		## this is the bit that takes a long time ... start a progress bar
		cat( 'Finding magnitudes of all sources within r_max.\n' )
		pb <- txtProgressBar(min = 0, max = n_radio_sources, style = 3)
		for ( ii in 1:n_radio_sources ){
			## find distances
			distances <- cenang( my_radio_cat$RA[ii], my_radio_cat$DEC[ii], band_dat$ALPHA_J2000, band_dat$DELTA_J2000 )*60*60
			candidate_index <- which( distances <= r_max ) 
			match_magnitudes <- c( match_magnitudes, band_dat[candidate_index,mag_cols[1]] )
		## update the progress bar
		setTxtProgressBar( pb, ii )
		}
		cat( '\n' )
		close( pb )
	} else match_magnitudes <- as.numeric( readLines( matched_mags ) )

	## -- find total(m)
	tmhist <- calculate_hist( match_magnitudes, mag_bins )
	total_m <- c()
	for ( ii in 1:length( tmhist$counts ) ) total_m <- c( total_m, sum( tmhist$counts[1:ii] ) )

	## -- find real(m)
	background <- nm * n_radio_sources * pi * r_max^2
	real_m <- total_m - background
	n_matches <- length( match_magnitudes )

	## -- find Q_0
	## --- Ciliegi+ (2003)
	Q0 <- ( n_matches - sum( nm * pi*r_max^2. * n_radio_sources ) ) / n_radio_sources 
	cat( 'Ciliegi+ Q_0:',Q0,'\n' )

	## --- fleuren+ (2012)
	## first calculate number of sources in the radio catalogue with NO counterpart within a radius of r, as a function of r
	## then do the same for a random catalogue
	## take the ratio. this will provide values that depend on r.
	## also calculate F(r) = 1 - exp( r^2/2sig^2 ) (eqn 11 in McAlpine)
	## then determine Q_0 via a fit to the ratio of U_obs/U_random
	## I think I'll do that in a separate function.
	radii <- seq( 1, 30 )
	fleuren_Q0 <- find_Q0_fleuren( my_radio_cat, radio_mask, band_dat, radii )
	fQ0 <- fleuren_Q0$q_0 ## the second element is the error

	## -- find qm: the expected distribution of the true counterparts
	qm <- real_m / sum( real_m ) * fQ0

	## find ratio of q(m)/n(m)
	qm_nm <- qm / nm

	## plot some things so far
	pdf( paste( 'LR_values_', my_band, '.pdf',sep='' ) )
	m <- rbind( c(1), c(2), c(3) )
	layout(m)
	lplot( tmhist$mids, log10( total_m ), xlim=c(13,24), y_lab='log(N(counterparts))', type='s', lwd=2, col='gray' )
	lines( tmhist$mids, log10( real_m ), lty=2, lwd=2, type='s' )
	lines( tmhist$mids, log10( background ), lty=3, lwd=2, type='s' )
	legend( 'topleft', c('Total','Real','Background'), col=c('gray','black','black'), lwd=2, lty=c(1,2,3), bty='n' )

	lplot( nmhist$mids, log10( qm_nm ), xlim=c(13,24), y_lab='log(P(m)=q(m)/n(m))', type='s', lwd=2 )

	lplot( nmhist$mids, log10( qm ), xlim=c(13,24), x_lab='Ks mag', y_lab='log(q(m))', type='s', lwd=2 )
	dev.off()

	## find the reliability values

	## loop through the radio sources   
	n_radio_sources <- dim(my_radio_cat)[1]
	pb <- txtProgressBar( min=0, max=n_radio_sources, style=3 )
	count <- 0 
	for ( ii in 1:n_radio_sources ){

		## CALCULATE f(r)
		## -- find distances 
		distances <- cenang( my_radio_cat$RA[ii], my_radio_cat$DEC[ii], band_dat$ALPHA_J2000, band_dat$DELTA_J2000 ) * 60. * 60.
		candidate_index <- which( distances <= r_max ) 
	if ( length( candidate_index ) > 0 ){
		tmp_dat <- band_dat[candidate_index,]
			band_magnitudes <- tmp_dat[,mag_cols[1]]
			f_r <- 1 / ( 2 * pi * sigma_pos[ii]^2 ) * exp( distances[candidate_index]^2 / ( 2 * sigma_pos[ii]^2 ) )

			## loop through candidates to calculate the LR
			LR <- c()
			for ( jj in 1:length(candidate_index) ){	
				## find the qm and nm for the magnitude bin of the source
				qm_nm_jj <- qm[ max( which( band_magnitudes[jj] > mag_bins ) ) ]
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

}

## plot completeness!  record information!
final_matches <- data.frame( radio_ID, match_band, video_ID, lr_value, lr_rel, n_cont, stringsAsFactors=FALSE )
write.table( final_matches, file='final_matches.csv', quote=FALSE, row.names=FALSE, sep=',' )



