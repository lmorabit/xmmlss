## PROCESS VIDEO RADIO DATA
##  written by Leah K. Morabito
##  started on 26 Jan 2017
##  last updated ...
##  contact: lmorabit@gmail.com

## starting from /data/xmmlss/
source('~/Dropbox/scripts/xmmlss/video_data.r')
source('~/Dropbox/scripts/xmmlss/multi_wavelength_comparison.r')

stilts_exec <- '/home/morabito/software/stilts/stilts'

## VIDEO master catalogue
master_cat <- 'cats/CFHTLS-W1_2016-04-14_fullcat_errfix.fits'

## STAR-GALAXY SEPARATION -----------------------------------------------------
if ( !file.exists( 'sg_separation_values.txt' ) ) source('~/Dropbox/scripts/xmmlss/star_galaxy_separation.r' )
## read in the values for separating stars and galaxies
sg_vals <- read.table( 'sg_separation_values.txt', stringsAsFactors=FALSE, header=TRUE )

## find the g-i and J-K values for the whole catalogue
master_cat_sg <- gsub( '.fits', '_SGsep.csv', master_cat )
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

    ss <- paste( stilts_exec, ' tmatch2 in1=GAL_IDs.csv ifmt1=csv in2=', master_cat, ' out=', gsub( '.fits', '_GSclass.fits', master_cat ), ' matcher=exact values1=GAL_ID values2=ID join=1and2', sep='' )
    system( ss )
}

master_cat <- gsub( '.fits', '_GSclass.fits', master_cat )
 


## MAKE A MASK FROM K-BAND -----------------------------------------------------
## STEP 2: MATCH BAND WITH VLA DATA
## 2.1: make a mask from K-band data
## 2.2: apply to radio data
## 2.3: look for cross matches

my_band <- 'K'
mag_cut <- 23.5
mag_col <- 'APER_2' 

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

## apply the mask from the VIDEO data
if ( !file.exists('VLA_masked.csv') ) masked_VLA <- apply_mask( VLA, radio_mask, filestem='VLA' ) else masked_VLA <- read.table( 'VLA_masked.csv', stringsAsFactors=FALSE, header=TRUE, sep=',' )
if ( !file.exists( 'sum_VLA_masked.csv' ) ) masked_sum_VLA <- apply_mask( sum_VLA, radio_mask, filestem='sum_VLA' ) else masked_sum_VLA <- read.table( 'sum_VLA_masked.csv', stringsAsFactors=FALSE, header=TRUE, sep=',' )
plot_sky_distribution( parm1='sum_VLA', ra1=sum_VLA$RA, de1=sum_VLA$DEC, parm2='Masked', ra2=masked_sum_VLA$RA, de2=masked_sum_VLA$DEC )

plot_radio_nearest_neighbours( VLA, plotfile='Nearest_neighbours_VLA.pdf' )
plot_radio_nearest_neighbours( sum_VLA, plotfile='Nearest_neighbours_sum_VLA.pdf' )

## CROSS-MATCH WITH BANDS -----------------------------------------------------

## select only the galaxies in the master catalogue
ss <- paste( stilts_exec, ' tpipe in=', master_cat, ' out=', gsub( 'GSclass', 'galaxies', master_cat ), ' cmd=\'select GS_CLASS==0\'',sep='' )
if ( !file.exists(gsub( 'GSclass', 'galaxies', master_cat )) ) system( ss )
galaxies_cat <- gsub( 'GSclass', 'galaxies', master_cat )


my_bands <- c('J','H','K','U','G','R','I','ZC','ZV','Y','YC')

for ( my_band in my_bands ){

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

    ## CALCULATE SOME POSITIONAL INFORMATION
    ## -- positional errors
    pos_noise <- sqrt( masked_sum_VLA$e_RA^2. + masked_sum_VLA$e_DEC^2. )
    sigma_pos <- sqrt( 0.1^2. + pos_noise^2. )
    r_max <- 5 * max( sqrt( sigma_pos^2 ) )
    cat( 'r_max is', r_max, 'arcsec.\n' )

    ## CALCULATE nm -- density of background sources
    ## -- find source counts in bins
    ## numbers ~ 1e-4
    mag_bins <- seq( floor(min( band_dat$K_MAG_APER_2 )), ceiling(max( band_dat$K_MAG_APER_2 )), 0.5 )
    nmhist <- hist( band_dat$K_MAG_APER_2, breaks=mag_bins )
    nm <- nmhist$counts / area_arcsec
    
    ## CALCULATE qm -- expected distribution of counter parts

    ##-- first find magnitudes of all sources within r_max of all radio positions
    n_radio_sources <- dim(masked_sum_VLA)[1]
    match_magnitudes <- c()
    for ( ii in 1:n_radio_sources ){
        ## find distances
        distances <- cenang( masked_sum_VLA$RA[ii], masked_sum_VLA$DEC[ii], band_dat$ALPHA_J2000, band_dat$DELTA_J2000 )*60*60
        candidate_index <- which( distances <= r_max ) 
        match_magnitudes <- c( match_magnitudes, band_dat$K_MAG_APER_2[candidate_index] )
    }
    ## -- find total(m)
    tmhist <- hist( match_magnitudes, breaks=mag_bins )
    total_m <- c()
    for ( ii in 1:length( tmhist$counts ) total_m <- c( total_m, sum( tmhist$counts[1:ii] ) )

    ## -- find real(m)
    real_m <- total_m - nm * n_radio_sources * pi * r_max^2 
    n_matches <- length( match_magnitudes )

    ## -- find Q_0
    ## --- Ciliegi+ (2003)
    Q_0 <- ( n_matches - sum( nm * pi*r_max^2. * n_radio_sources ) ) / n_radio_sources 

    ## --- fleuren+ (2012)
    ## first calculate number of sources in the radio catalogue with NO counterpart within a radius of r, as a function of r
    ## then do the same for a random catalogue
    ## take the ratio. this will provide values that depend on r.
    ## also calculate F(r) = 1 - exp( r^2/2sig^2 ) (eqn 11 in McAlpine)
    ## then determine Q_0 via a fit to the ratio of U_obs/U_random

    ## -- find qm
    qm <- real_m / sum( real_m )

    ## find ratio of q(m)/n(m)
    qm_nm <- qm / nm

    ## loop through the radio sources   
    n_radio_sources <- dim(masked_sum_VLA)[1]
    for ( ii in 1:n_radio_sources ){

        ## CALCULATE f(r)
        ## -- find distances 
        distances <- cenang( masked_sum_VLA$RA[ii], masked_sum_VLA$DEC[ii], band_dat$ALPHA_J2000, band_dat$DELTA_J2000 )*60*60
        candidate_index <- which( distances <= r_max ) 
        band_magnitudes <- band_dat$K_MAG_APER_2[candidate_index]
        f_r <- 1 / ( 2 * pi * sigma_pos[ii]^2 ) * exp( distances[candidate_index]^2 / ( 2 * sigma_pos[ii]^2 ) )
        ## convert to radians? -- blargh
        #f_r <- 1 / ( 2 * pi * (sigma_pos[ii]/206265)^2 ) * exp( (distances[candidate_index]/206265)^2. / ( 2 * (sigma_pos[ii]/206265)^2 ) )
        ## these are very big numbers!
 
        n_matches <- length( candidate_index )
        Q_0 <- n_matches - sum( nm * pi*r_max^2 * n_radio_sources ) / n_radio_sources
        qm <- almost_qm * Q_0

        ## loop through candidates to calculate the LR
        LR <- c()
        for ( jj in 1:length(candidate_index) ){    
            ## find the qm and nm for the magnitude bin of the source
            qm_jj <- qm[ max( which( band_magnitudes[jj] > mag_bins ) ) ]
            nm_jj <- nm[ max( which( band_magnitudes[jj] > mag_bins ) ) ]
            LR <- c( LR, qm_jj * f_r[jj] / nm_jj )
        }

        LR_reliability <- c()
        rel_factor <- sum( LR ) + ( 1 - Q_0 )
        for ( jj in 1:length(candidate_index ) ) LR_reliability <- c( LR_reliability, LR[jj] / rel_factor )
        
    }

}



### BELOW THIS IS OLD STUFF ------------------------------------------------

## start with isolated sources
## find likelihood ratios

## read in the cross-matched file
#cross_matched_cat <- 'cats/VLA_CFHTLS-W1_2016-04-14_fullcat_errfix_2arcsec_matched.csv'
#cross_matched_cat <- 'CFHTLS-W1_2016-04-14_fullcat_errfix_slim_Ks_23.5_VLA_sum_masked.csv'
cm_dat <- read_in_catalogue( 'xmatch.csv' )
## fix the column names
mycols <- colnames( cm_dat )
for ( ii in 1:length(mycols) ) {
    mycols[ii] <- gsub( '_1', '', mycols[ii] )
    tmp <- strsplit( mycols[ii], '_' )[[1]]
    if ( tmp[length(tmp)] == '2' ) mycols[ii] <- paste('TMPCOL_', ii, sep='' )
}
colnames( cm_dat ) <- mycols

my_bands <- c('J','H','K','U','G','R','I','ZC','ZV','Y','YC')
band_names <- paste( my_bands, '_MAG_AUTO', sep='' )
band_error_names <- paste( my_bands, '_MAGERR_AUTO', sep='' )
filename='Fraction_of_matches.pdf'
plot_radio_detection_fractions( cm_dat, my_bands, band_names, band_error_names, plotfile=filename )

## likelihood matching


cross_match_positional_errors( cm_dat, my_bands )
plot_radio_nearest_neighbours( VLA, plotfile='Nearest_neighbours_VLA.pdf' )
plot_radio_nearest_neighbours( sum_VLA, plotfile='Nearest_neighbours_sum_VLA.pdf' )
plot_radio_nearest_neighbours( cm_dat )

