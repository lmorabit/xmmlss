## PROCESS VIDEO RADIO DATA
##  written by Leah K. Morabito
##  started on 26 Jan 2017
##  last updated ...
##  contact: lmorabit@gmail.com

## starting from /data/xmmlss/
source('~/Dropbox/scripts/xmmlss/video_data.r')
source('~/Dropbox/scripts/xmmlss/multi_wavelength_comparison.r')

## flags for processing steps
start_vla_from_scratch <- FALSE
make_a_footprint_mask <- FALSE


##-- VLA CATALOGUE
if ( start_vla_from_scratch ){
    vla_file <- 'VLA/13B-308_080916_csv.txt'
    VLA <- read_VLA_data( vla_file )

    ## plot the sky distribution
    plot_sky_distribution( parm1='VLA', ra1=VLA$RA, de1=VLA$DEC )
    ## using the entire catalogue, which is uncombined
    plot_distributions( VLA, plotfile='VLA_distributions.pdf' )

    ## sum things by Source_id
    sum_VLA <- sum_components( VLA )
    plot_sky_distribution( parm1='sum_VLA', ra1=VLA$RA, de1=VLA$DEC )
    plot_distributions( sum_VLA, plotfile='sum_VLA_distributions.pdf' )

} else {
    ## read in sum_VLA data
    VLA <- read.table( 'VLA/13B-308_080916.csv', header=TRUE, sep=',', stringsAsFactors=FALSE )
    sum_VLA <- read.table( 'sum_VLA.csv', header=TRUE, sep=',', stringsAsFactors=FALSE )

}

## MAKE A MASK FROM THE IR DATA
## To do this I used the fullcat but with Ks_band_only so the file size is smaller
infile <- 'CFHTLS-W1_2016-04-14_fullcat_errfix_Ks_band_only.csv'
if ( make_a_footprint_mask ) {
    ## MAKE A MASK
    my_mask <- create_footprint_mask( infile, 'ALPHA_J2000', 'DELTA_J2000', cellsize=60, method='2dhist' )
} else { 
    ## READ IN A MASK ALREADY MADE
    infile <- Sys.glob( 'mask*txt' )
    ## open a file connection
    zz <- file( infile, 'r' )
    x_vals <- readLines( con=zz, n=1 )
    my_x <- as.numeric( strsplit( x_vals, ',' )[[1]] ) 
    y_vals <- readLines( con=zz, n=1 )
    my_y <- as.numeric( strsplit( y_vals, ',' )[[1]] ) 
    close( zz )
    mask_counts <- read.table( infile, skip=2, stringsAsFactors=FALSE, header=FALSE )
    my_mask <- data.frame( x.breaks=my_x, y.breaks=my_y, counts=mask_counts )
}


## APPLY THE MASK
masked_VLA <- apply_mask( VLA, my_mask, filestem='VLA' )
masked_sum_VLA <- apply_mask( sum_VLA, my_mask, filestem='sum_VLA' )

## plot the difference
plot_sky_distribution( parm1='sum_VLA', ra1=sum_VLA$RA, de1=sum_VLA$DEC, parm2='Masked', ra2=masked_sum_VLA$RA, de2=masked_sum_VLA$DEC )


## CROSS-MATCH IN TOPCAT!!!

## read in the cross-matched file
#cross_matched_cat <- 'cats/VLA_CFHTLS-W1_2016-04-14_fullcat_errfix_2arcsec_matched.csv'
cross_matched_cat <- 'CFHTLS-W1_2016-04-14_fullcat_errfix_sum_VLA_mask_2arcsec_matched.csv'
cm_dat <- read_in_catalogue( cross_matched_cat )

my_bands <- c('J','H','K','U','G','R','I','ZC','ZV','Y','YC')
band_names <- paste( my_bands, '_MAG_AUTO', sep='' )
band_error_names <- paste( my_bands, '_MAGERR_AUTO', sep='' )
filename='Fraction_of_matches.pdf'
plot_radio_detection_fractions( cm_dat, my_bands, band_names, band_error_names, plotfile=filename )

cross_match_positional_errors( cm_dat, my_bands )
plot_radio_nearest_neighbours( VLA, plotfile='Nearest_neighbours_VLA.pdf' )
plot_radio_nearest_neighbours( sum_VLA, plotfile='Nearest_neighbours_sum_VLA.pdf' )
plot_radio_nearest_neighbours( cm_dat )

