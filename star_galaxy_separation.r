my_bands <- c( 'K', 'J', 'I', 'G' )
mag_cut <- 30
mag_col <- 'APER_2' 
basic_cols <- c( 'ID', 'ALPHA_J2000', 'DELTA_J2000', 'HALOFLAG' )

for ( my_band in my_bands ){

    ## set up output files, etc.
    mag_cols <- paste( my_band, c( '_MAG_', '_MAGERR_' ), mag_col, sep='' )
    keep_cols <- c( basic_cols, mag_cols, 'SNR' )
    my_cat_name <- gsub( '.fits', paste( '_', my_band, '_', mag_col, '_', mag_cut, '.csv', sep='' ), master_cat )

    ## this stilts command reads in the master catalogue and outputs a catalogue keeping only keep_cols
    ## where the magnitude cut has been applied, and snr is calculated and added to the table. 
    ss <- paste( stilts_exec, " tpipe in=", master_cat, " out=", my_cat_name, " omode=out ofmt=csv cmd=\'select ", mag_cols[1], "<=", mag_cut, "\' cmd=\'select HALOFLAG==0\' cmd=\'addcol SNR 1.08574/", mag_cols[2], '\' cmd=\'select SNR>=5\' cmd=\'keepcols "', paste( keep_cols, collapse=' ' ), '"', "\' ", sep='' )
    system( ss )

}

## combine the catalogues together
my_files <- Sys.glob('cats/*30.csv') ## these will be sorted alphabetically
G_dat <- read.table( my_files[1], stringsAsFactors=FALSE, header=TRUE, sep=',' )
I_dat <- read.table( my_files[2], stringsAsFactors=FALSE, header=TRUE, sep=',' )
J_dat <- read.table( my_files[3], stringsAsFactors=FALSE, header=TRUE, sep=',' )
K_dat <- read.table( my_files[4], stringsAsFactors=FALSE, header=TRUE, sep=',' )

## match the IDs
good_IDs <- Reduce( intersect, list( G_dat$ID, I_dat$ID, J_dat$ID, K_dat$ID ) )
cat( 'There are', length( good_IDs ), 'with a five-sigma detection in G,I,J,K.\n' )

g_mag <- G_dat$G_MAG_APER_2[ which(G_dat$ID %in% good_IDs) ]
i_mag <- I_dat$I_MAG_APER_2[ which(I_dat$ID %in% good_IDs) ]
J_mag <- J_dat$J_MAG_APER_2[ which(J_dat$ID %in% good_IDs) ]
K_mag <- K_dat$K_MAG_APER_2[ which(K_dat$ID %in% good_IDs) ]

df <- data.frame( g_mag, i_mag, J_mag, K_mag )
good_index <- which( df$K_mag <= 23.5 )

df <- df[ good_index, ]

g_minus_i <- df$g_mag - df$i_mag
J_minus_K <- df$J_mag - df$K_mag

x_lim <- c(0,4)
y_lim <- c(-1,2.5)
pdf( 'sg_separation.pdf' )
lplot( g_minus_i, J_minus_K , pch='.', xlim=x_lim, ylim=y_lim, x_lab='g-i [mag]', y_lab='J-K [mag]', col=mycols[3] )
# col='gray48'
## make a nice density plot?
#stuff <- kde2d( g_minus_i, J_minus_K, n=c(100,100), lims=c(x_lim,y_lim) )
#contour( stuff )

## now fit the stellar locus ...
## first limit the outliers
x_index <- which( g_minus_i >= 0.4 & g_minus_i <= 1.9 )
y_index <- which( J_minus_K >= y_lim[1] & J_minus_K <= y_lim[2] )
good_data <- intersect( x_index, y_index )

x <- g_minus_i[good_data]
y <- J_minus_K[good_data]

## also limit the galaxy contributin
x_div <- seq( floor(min(x)), ceiling(max(x)), 0.1 )
y_div <- 0.4 * x_div - 0.5
#lines( x_div, y_div, col='red' )

## find where points are below this line ...
y_vals <- 0.4 * x - 0.5 
sl_vals <- which( y <= y_vals )

x <- x[sl_vals]
x2 <- x^2.
y <- y[sl_vals]
quad_model <- lm( y ~ x + x2 )

x_loc <- seq( 0.4, 1.9, 0.1 )
pred_vals <- predict( quad_model, list( x=x_loc, x2=x_loc^2 ) )
lines( x_loc, pred_vals , col='black', lwd=2.5 )
x_up <- seq( 1.9, 10, 0.1 )
lines( x_up, rep(pred_vals[which( x_loc == max(x_loc)) ], length( x_up ) ), col='black', lwd=2.5 )
x_down <- seq( -1, 0.4, 0.1 )
lines( x_down, rep( pred_vals[ which( x_loc == 0.4 ) ], length( x_down ) ), col='black', lwd=2.5 )

## shift it up
shift_val <- 0.12
lines( x_loc, pred_vals+shift_val , col='black', lty=2, lwd=2.5 )
lines( x_up, rep(pred_vals[which( x_loc == max(x_loc)) ], length( x_up ) )+shift_val, col='black', lty=2, lwd=2.5 )
lines( x_down, rep( pred_vals[ which( x_loc == 0.4 ) ], length( x_down ) )+shift_val, col='black', lty=2, lwd=2.5 )
legend( 'topright', c(expression(0.12+f[locus]),expression(f[locus])), lty=c(2,1), lwd=2.5, bty='n', cex=1.5 )
dev.off()

## write a file with the results
c1 <- pred_vals[ which( x_loc == 0.4 ) ]
c2 <- quad_model$coefficients[[1]]
b2 <- quad_model$coefficients[[2]]
a2 <- quad_model$coefficients[[3]]
c3 <- pred_vals[which( x_loc == max(x_loc)) ]

df <- data.frame( c1, c2, b2, a2, c3 )
write.table( df, file='sg_separation_values.txt', row.names=FALSE, quote=FALSE )

cat( 'f_locus =', c1, ' ( x < 0.4 )\n' )
cat( 'f_locus =', c2,'+', b2,'x',a2,'x^2  (0.4 < x < 1.9 )\n' )
cat( 'f_locus =', c3, ' ( x > 1.9 )\n' )

