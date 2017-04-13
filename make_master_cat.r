stilts_exec <- '/vardy/leah/stilts/stilts'
master_cat <- '/vardy/videouser/v1.3-20160414/cats/CFHTLS-W1_2016-04-14_fullcat_errfix.fits'

## basic columns to keep
basic_cols <- c( 'ID', 'ALPHA_J2000', 'DELTA_J2000', 'HALOFLAG', 'CONTEXT' )

## which magnitude values to use
mag_col_to_use <- 'MAG_AUTO'
mag_err_col_to_use <- gsub( 'MAG', 'MAGERR', mag_col_to_use ) 

## these are the bands
my_bands <- c('J','H','K','U','G','R','I','ZC','ZV','Y','YC')

## magnitudes
keep_bands <- paste( my_bands, mag_col_to_use, sep='_' )
keep_band_errs <- paste( my_bands, mag_err_col_to_use, sep='_' )

## other band info
tmp_keep <- c( '_FLAGS', '_CLASS_STAR', '_DET_FLAG' )
other_cols <- c()
for ( my_band in my_bands ) other_cols <- c( other_cols, paste( my_band, tmp_keep, sep='' ) )

## calculate SNR
add_snr <- c()
for ( keep_band_err in keep_band_errs ) add_snr <- c( add_snr, paste( " cmd=\'addcol ", strsplit( keep_band_err, '_' )[[1]][1], "_SNR 1.08574/", keep_band_err, "\'",  sep='' ) )

## also calculate G-I and J-K (and G-R)
gijkgr <- paste( " cmd=\'addcol G_MINUS_I ", keep_bands[which( my_bands=='G' )], "-", keep_bands[which( my_bands=='I' )], "\' cmd=\'addcol J_MINUS_K ", keep_bands[which( my_bands=='J' )], "-", keep_bands[which( my_bands=='K' )], "\' cmd=\'addcol G_MINUS_R ", keep_bands[which( my_bands=='G' )], "-", keep_bands[which( my_bands=='R' )], "\'", sep='' )

keep_cols <- c( basic_cols, keep_bands, keep_band_errs, other_cols, paste( my_bands, '_SNR', sep='' ), c( 'G_MINUS_I', 'G_MINUS_R', 'J_MINUS_K' ) )


outcat <- gsub( '/vardy/videouser/v1.3-20160414/', '', master_cat )
outcat <- gsub( '.fits', '_MAGAUTO.fits', outcat )

ss <- paste( stilts_exec, " tpipe in=", master_cat, " out=", outcat, " omode=out ", paste( add_snr, collapse='' ), gijkgr, ' cmd=\'keepcols "', paste( keep_cols, collapse=' ' ), '"', "\' ", sep='' )
system( ss )


