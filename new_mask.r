## read in the K-band catalogue
band_dat <- read.table( 'cats/CFHTLS-W1_2016-04-14_fullcat_errfix_MAGAUTO_Galaxies_K.csv', sep=',', stringsAsFactors=FALSE, header=TRUE )

## get K-band detected values
kband_det <- which( substr( as.character( band_dat$ID ), 1, 1 ) == '3' )
kband_dat <- band_dat[ kband_det, ]

## get five-sigma values
kband_dat_fivesig <- kband_dat[ which(kband_dat$K_SNR>=5), ]

## and get things with no flags
kband_dat_fivesig_noflags <- kband_dat_fivesig[ which( kband_dat_fivesig$K_FLAGS < 1 ), ]

## define an overall rectangle
xmin <- min( kband_dat_fivesig_noflags$ALPHA_J2000 )
xmax <- max( kband_dat_fivesig_noflags$ALPHA_J2000 )
ymin <- min( kband_dat_fivesig_noflags$DELTA_J2000 )
ymax <- max( kband_dat_fivesig_noflags$DELTA_J2000 )

pdf( 'area_definition.pdf' )
plot( kband_dat_fivesig_noflags$ALPHA_J2000, kband_dat_fivesig_noflags$DELTA_J2000, pch=16, cex=0.1, col='black', xlab='RA [deg]', ylab='Dec [deg]' )
#abline( v=c(33.83938,37.18655), col='purple' )
#abline( h=c(-5.641043,-3.985834), col='purple' )

## field division lines to get the bottom border
div_lines = c( 35.15, 36.2 )
#abline( v=div_lines, col='green' )

new_ymin1 = min( kband_dat_fivesig_noflags$DELTA_J2000[which( kband_dat_fivesig_noflags$ALPHA_J2000 > div_lines[1] )] )
new_ymin2 = min( kband_dat_fivesig_noflags$DELTA_J2000[which( kband_dat_fivesig_noflags$ALPHA_J2000 > div_lines[2] )] )
#abline( h=c(new_ymin1, new_ymin2), col='purple' )

## field division lines to get the top border
div_lines = c( 34.85, 35.9)
#abline( v=div_lines, col='green' )
new_ymax1 = max( kband_dat_fivesig_noflags$DELTA_J2000[which( kband_dat_fivesig_noflags$ALPHA_J2000 < div_lines[1] )] )
new_ymax2 = max( kband_dat_fivesig_noflags$DELTA_J2000[which( kband_dat_fivesig_noflags$ALPHA_J2000 < div_lines[2] )] )
#abline( h=c(new_ymax1, new_ymax2), col='purple' )

## cutouts bottom border
c_index_1 <- which( kband_dat_fivesig_noflags$ALPHA_J2000 < 35.15 & kband_dat_fivesig_noflags$DELTA_J2000 < new_ymin1 )
c_xmax1 = max( kband_dat_fivesig_noflags$ALPHA_J2000[ c_index_1 ] )
c_index_2 <- which( kband_dat_fivesig_noflags$ALPHA_J2000 > 35.15 & kband_dat_fivesig_noflags$DELTA_J2000 < new_ymin2 )
c_xmin1 = min( kband_dat_fivesig_noflags$ALPHA_J2000[ c_index_2 ] )
#abline( v=c(c_xmax1,c_xmin1), col='purple' )

c_ymax1 = min( kband_dat_fivesig_noflags$DELTA_J2000[which(kband_dat_fivesig_noflags$ALPHA_J2000 > c_xmax1+0.01 & kband_dat_fivesig_noflags$ALPHA_J2000 < c_xmin1-0.01 )] )
#abline( h=c(c_ymax1 ), col='purple' )

c_index_3 <- which( kband_dat_fivesig_noflags$ALPHA_J2000 < 36.2 & kband_dat_fivesig_noflags$DELTA_J2000 < new_ymin2 )
c_xmax2 = max( kband_dat_fivesig_noflags$ALPHA_J2000[ c_index_3 ] )
c_index_4 <- which( kband_dat_fivesig_noflags$ALPHA_J2000 > 36.2 & kband_dat_fivesig_noflags$DELTA_J2000 < -5.42 )
c_xmin2 = min( kband_dat_fivesig_noflags$ALPHA_J2000[ c_index_4 ] )
#abline( v=c(c_xmax2,c_xmin2), col='purple' )

c_ymax2 = min( kband_dat_fivesig_noflags$DELTA_J2000[which(kband_dat_fivesig_noflags$ALPHA_J2000 > c_xmax2+0.01 & kband_dat_fivesig_noflags$ALPHA_J2000 < c_xmin2-0.01 )] )
#abline( h=c(c_ymax2 ), col='purple' )

## cutouts top border
c_index_5 <- which( kband_dat_fivesig_noflags$ALPHA_J2000 < 34.85 & kband_dat_fivesig_noflags$DELTA_J2000 > -4.23 )
c_xmax3 = max( kband_dat_fivesig_noflags$ALPHA_J2000[c_index_5] )
c_index_6 <- which( kband_dat_fivesig_noflags$ALPHA_J2000 > 34.85 & kband_dat_fivesig_noflags$DELTA_J2000 > -4.23 )
c_xmin3 = min( kband_dat_fivesig_noflags$ALPHA_J2000[c_index_6] )
c_ymin3 = max( kband_dat_fivesig_noflags$DELTA_J2000[ which(kband_dat_fivesig_noflags$ALPHA_J2000 > c_xmax3 & kband_dat_fivesig_noflags$ALPHA_J2000 < c_xmin3 )] )
#abline( v=c(c_xmax3,c_xmin3), col='purple' )
#abline( h=c(c_ymin3 ), col='purple' )

c_xmin4 = min( kband_dat_fivesig_noflags$ALPHA_J2000[which(kband_dat_fivesig_noflags$DELTA_J2000 > new_ymax2)] )
#abline( v=c(c_xmin4), col='purple' )

## polygon vectors
xvec <- c( xmin, xmin, c_xmax3, c_xmax3, c_xmin3, c_xmin3, c_xmin4, c_xmin4, xmax, xmax, c_xmin2, c_xmin2, c_xmax2, c_xmax2, c_xmin1, c_xmin1, c_xmax1, c_xmax1 )
yvec <- c( ymin, new_ymax1, new_ymax1, c_ymin3, c_ymin3, new_ymax2, new_ymax2, ymax, ymax, new_ymin2, new_ymin2, c_ymax2, c_ymax2, new_ymin1, new_ymin1, c_ymax1, c_ymax1, ymin )

polygon( xvec, yvec, border='purple' )


## interior cutouts

div_y <- c(-5.3,-5.35)
#abline( h=div_y, col='green' )

ic_index_1 <- which( kband_dat_fivesig_noflags$DELTA_J2000 > div_y[2] & kband_dat_fivesig_noflags$DELTA_J2000 < div_y[1] )
ic_xmin1 <- max( kband_dat_fivesig_noflags$ALPHA_J2000[ intersect( ic_index_1, which( kband_dat_fivesig_noflags$ALPHA_J2000 < c_xmin1-0.01 ) ) ] )
ic_xmax1 <- min( kband_dat_fivesig_noflags$ALPHA_J2000[ intersect( ic_index_1, which( kband_dat_fivesig_noflags$ALPHA_J2000 > c_xmax1 ) ) ] )

#abline( v=c( ic_xmin1, ic_xmax1 ), col='purple' )

ic_index_2 <- which( kband_dat_fivesig_noflags$ALPHA_J2000 > ic_xmin1+0.01 & kband_dat_fivesig_noflags$ALPHA_J2000 < ic_xmax1-0.01 )
ic_ymin1 <- max( kband_dat_fivesig_noflags$DELTA_J2000[ intersect( ic_index_2, which( kband_dat_fivesig_noflags$DELTA_J2000 < -5.3) ) ] )
ic_ymax1 <- min( kband_dat_fivesig_noflags$DELTA_J2000[ intersect( ic_index_2, which( kband_dat_fivesig_noflags$DELTA_J2000 > -5.35) ) ] )

#abline( h=c( ic_ymin1, ic_ymax1 ), col='purple' )

ic1_xvec <- c( ic_xmin1, ic_xmin1, ic_xmax1, ic_xmax1 )
ic1_yvec <- c( ic_ymin1, ic_ymax1, ic_ymax1, ic_ymin1 )

polygon( ic1_xvec, ic1_yvec, border='purple' )

div_y <- c( -5.26, -5.22 )

ic_index_3 <- which( kband_dat_fivesig_noflags$DELTA_J2000 > div_y[1] & kband_dat_fivesig_noflags$DELTA_J2000 < div_y[2] )
ic_xmin2 <- max( kband_dat_fivesig_noflags$ALPHA_J2000[ intersect( ic_index_3, which( kband_dat_fivesig_noflags$ALPHA_J2000 < c_xmin2-0.01 ) ) ] )
ic_xmax2 <- min( kband_dat_fivesig_noflags$ALPHA_J2000[ intersect( ic_index_3, which( kband_dat_fivesig_noflags$ALPHA_J2000 > c_xmax2+0.01 ) ) ] )

#abline( v=c( ic_xmin2, ic_xmax2 ), col='purple' )

ic_index_4 <- which( kband_dat_fivesig_noflags$ALPHA_J2000 > ic_xmin2+0.01 & kband_dat_fivesig_noflags$ALPHA_J2000 < ic_xmax2-0.01 )
ic_ymin2 <- max( kband_dat_fivesig_noflags$DELTA_J2000[ intersect( ic_index_4, which( kband_dat_fivesig_noflags$DELTA_J2000 < -5.22) ) ] )
ic_ymax2 <- min( kband_dat_fivesig_noflags$DELTA_J2000[ intersect( ic_index_4, which( kband_dat_fivesig_noflags$DELTA_J2000 > -5.26) ) ] )

#abline( h=c( ic_ymin2, ic_ymax2 ), col='purple' )

ic2_xvec <- c( ic_xmin2, ic_xmin2, ic_xmax2, ic_xmax2 )
ic2_yvec <- c( ic_ymin2, ic_ymax2, ic_ymax2, ic_ymin2 )

polygon( ic2_xvec, ic2_yvec, border='purple' )

div_y <- c(-4.41, -4.37)

ic_index_5 <- which( kband_dat_fivesig_noflags$DELTA_J2000 > div_y[1] & kband_dat_fivesig_noflags$DELTA_J2000 < div_y[2] )
ic_xmin3 <- max( kband_dat_fivesig_noflags$ALPHA_J2000[ intersect( ic_index_5, which( kband_dat_fivesig_noflags$ALPHA_J2000 < c_xmin3-0.01 ) ) ] )
ic_xmax3 <- min( kband_dat_fivesig_noflags$ALPHA_J2000[ intersect( ic_index_5, which( kband_dat_fivesig_noflags$ALPHA_J2000 > c_xmax3 ) ) ] )

#abline( v=c( ic_xmin3, ic_xmax3 ), col='purple' )

ic_index_6 <- which( kband_dat_fivesig_noflags$ALPHA_J2000 > ic_xmin3+0.01 & kband_dat_fivesig_noflags$ALPHA_J2000 < ic_xmax3-0.01 )
ic_ymin3 <- max( kband_dat_fivesig_noflags$DELTA_J2000[ intersect( ic_index_6, which( kband_dat_fivesig_noflags$DELTA_J2000 < -4.37) ) ] )
ic_ymax3 <- min( kband_dat_fivesig_noflags$DELTA_J2000[ intersect( ic_index_6, which( kband_dat_fivesig_noflags$DELTA_J2000 > -4.41) ) ] )

#abline( h=c( ic_ymin3, ic_ymax3 ), col='purple' )

ic3_xvec <- c( ic_xmin3, ic_xmin3, ic_xmax3, ic_xmax3 )
ic3_yvec <- c( ic_ymin3, ic_ymax3, ic_ymax3, ic_ymin3 )

polygon( ic3_xvec, ic3_yvec, border='purple' )
dev.off()

include_x <- c( xmin, xmax )
include_y <- c( ymin, ymax )

exclude_x1 <- c( c_xmax1, xmax )
exclude_y1 <- c( ymin, new_ymin1 )

exclude_x2 <- c( c_xmax2, xmax )
exclude_y2 <- c( ymin, new_ymin2 )

exclude_x3 <- c( c_xmax1, c_xmin1 )
exclude_y3 <- c( ymin, c_ymax1 )

exclude_x4 <- c( c_xmax2, c_xmin2 )
exclude_y4 <- c( ymin, c_ymax2 )

exclude_x5 <- c( xmin, c_xmin4 )
exclude_y5 <- c( new_ymax2, ymax )

exclude_x6 <- c( xmin, c_xmin3 )
exclude_y6 <- c( new_ymax1, ymax )

exclude_x7 <- c( c_xmax3, c_xmin3 )
exclude_y7 <- c( c_ymin3, ymax )

exclude_x8 <- c( ic_xmin1, ic_xmax1 )
exclude_y8 <- c( ic_ymin1, ic_ymax1 )

exclude_x9 <- c( ic_xmin2, ic_xmax2 )
exclude_y9 <- c( ic_ymin2, ic_ymax2 )

exclude_x10 <- c( ic_xmin3, ic_xmax3 )
exclude_y10 <- c( ic_ymin3, ic_ymax3 )


## write these to a file
my_info <- data.frame( include_x, include_y, exclude_x1, exclude_y1, exclude_x2, exclude_y2, exclude_x3, exclude_y3, exclude_x4, exclude_y4, exclude_x5, exclude_y5, exclude_x6, exclude_y6, exclude_x7, exclude_y7, exclude_x8, exclude_y8, exclude_x9, exclude_y9, exclude_x10, exclude_y10 )

write.table( my_info, file='Kband_area_definition.dat', row.names=FALSE, quote=FALSE )

