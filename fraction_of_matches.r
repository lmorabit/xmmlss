source('../../xmmlss/plotting_libraries.r')
source('../../xmmlss/leah_helper_functions.r')

############ HOUSEKEEPING

rel_threshold <- 0.8

#########################

## read in the data
final_matches <- read.table( 'final_matches.csv', stringsAsFactors=FALSE, sep=',', header=TRUE )
my_radio_cat <- read.table( 'VLA_fivesig_masked.csv', stringsAsFactors=FALSE, header=TRUE, sep=',' )

## get some information on the results of the cross-matching
n_radio_sources <- dim( my_radio_cat )[1]
n_radio_matched <- length( unique( final_matches$radio_ID ) )
cat( n_radio_matched/n_radio_sources*100, 'percent of sources matched.\n' )

rel_index <- which( final_matches$lr_rel > rel_threshold )
final_rel <- final_matches[ rel_index, ]
n_rel_matched <- length( unique( final_rel$radio_ID ) )
cat( n_rel_matched/n_radio_sources*100, 'percent of sources reliably matched.\n' )

## find unresolved sources
unresolved_index <- find_unresolved_index( my_radio_cat$Total_flux, my_radio_cat$Peak_flux, my_radio_cat$e_Total_flux, my_radio_cat$e_Peak_flux )

## split things by s, m, c
s_index <- which( my_radio_cat$S_Code == 'S' )
m_index <- which( my_radio_cat$S_Code == 'M' )
c_index <- which( my_radio_cat$S_Code == 'C' )

## get a list of the bands
matched_bands <- unique( final_rel$match_band )

## find the number of matches for all
all_matches <- c()
for ( mb in matched_bands ) all_matches <- c( all_matches, length( unique( final_rel[ which( final_rel$match_band == mb ), 'radio_ID' ] ) ) )

## divide into single, multiple, complex
s_radio_id <- my_radio_cat$Source_id[ s_index ]
m_radio_id <- my_radio_cat$Source_id[ m_index ]
c_radio_id <- my_radio_cat$Source_id[ c_index ]

s_matches <- c()
m_matches <- c()
c_matches <- c()
for ( mb in matched_bands ) { 
    tmp_vec <- final_rel$radio_ID[ which( final_rel$match_band == mb ) ] 
    s_matches <- c( s_matches, length( unique( tmp_vec[ which( tmp_vec %in% s_radio_id ) ] ) ) )
    m_matches <- c( m_matches, length( unique( tmp_vec[ which( tmp_vec %in% m_radio_id ) ] ) ) )
    c_matches <- c( c_matches, length( unique( tmp_vec[ which( tmp_vec %in% c_radio_id ) ] ) ) )
}

percentage_of_matches <- all_matches / n_radio_sources
percentage_of_s_matches <- s_matches / length( s_radio_id )
percentage_of_m_matches <- m_matches / length( m_radio_id )
percentage_of_c_matches <- c_matches / length( c_radio_id )

plotcols <- viridis( length( all_matches )+2 )

my_matches <- rbind( percentage_of_matches, percentage_of_s_matches, percentage_of_m_matches, percentage_of_c_matches )
col_vec <- plotcols[2:(length(plotcols)-1)] 
my_cols <- rbind( col_vec, tcol( col_vec, trans=150 ), tcol( col_vec, trans=100 ), tcol( col_vec, trans=100 ) )
my_dens <- rbind( rep( 200, length( percentage_of_matches ) ), rep( 75, length( percentage_of_matches ) ), rep( 50, length( percentage_of_matches ) ), rep( 25, length( percentage_of_matches ) ) )
my_angle <- rbind( rep( 45, length( percentage_of_matches ) ), rep( 45, length( percentage_of_matches ) ), rep( -45, length( percentage_of_matches ) ), rep( 45, length( percentage_of_matches ) ) )


pdf( 'Matches_by_band.pdf')
par( mar=c(5,5,2,2) )
barplot( my_matches, names.arg = matched_bands, col=my_cols, density=my_dens, angle=my_angle, ylim=c(0,1), beside=TRUE )
box( which="plot" )
mtext( "Fraction of Radio matches", side=2, line=3, cex=1.25 )
legend( 'topright', c('All', 'S', 'M', 'C' ), col=c('lightgray', tcol('lightgray',trans=150), tcol('lightgray',trans=100), tcol('lightgray',trans=100) ), density=c(200,75,50,25), angle=c(45,45,-45,45), bty='n', cex=1.5 )
#legend( 'topright', c('All', 'S' ), col=c('lightgray', tcol('lightgray',trans=150) ), density=c(200,75), angle=c(45,45), bty='n', cex=1.5 )
dev.off()



match_cols <- viridis( length( matched_bands )+2 )[2:(length(matched_bands)+1)]


all_matches <- c()
all_colors <- c()
## find out the number of sources detected in more than one band
for ( ii in 1:length(matched_bands) ){
    band_ids <- final_rel$radio_ID[ which( final_rel$match_band == matched_bands[ii] ) ]
    other_bands <- c()
    for ( jj in 1:length(matched_bands) ) {
        tmp_dat <- final_rel[ which( final_rel$match_band == matched_bands[jj] ), ]
        other_bands <- c( other_bands, length( which( tmp_dat$radio_ID %in% band_ids ) )/length(band_ids) )
    }
    no_match_index <- which( matched_bands != matched_bands[ii] )
    all_matches <- cbind( all_matches, other_bands )
    all_colors <- cbind( all_colors, match_cols )

}

pdf( 'Multi-band_fractions.pdf' )
par( mar=c(5,5,5,2) )
barplot( all_matches, col=all_colors, names.arg=matched_bands, beside=TRUE, ylim=c(0,1) )
#legend( 0, -0.05, rep('',length(matched_bands)), col=match_cols, pch=15, xpd=TRUE, horiz=TRUE, cex=2.2, bty='n' )
box( which='plot', lwd=2 )
dev.off()

pdf( 'Multi-band_fractions_2D.pdf' )
plotcols <- viridis( 256 )
par( mar=c( 5,5,7,7) )
image( all_matches, col=plotcols[1:250], axes=FALSE )
grid( nx=length(matched_bands), ny=length(matched_bands), col='gray48' )
xstep <- 1/(length(matched_bands)-1)
atvec <- seq( 0,1,xstep )
axis( 1, at=atvec, labels=matched_bands )
axis( 2, at=atvec, labels=matched_bands )
rgn <- par('usr')
color.legend( rgn[2]+0.2*xstep, rgn[3], rgn[2]+xstep, rgn[4], legend=format(atvec,digits=2), rect.col=plotcols[seq(1,251,25)], gradient="n", align="rb", lwd=3 )
box( which='plot', lwd=3 )
mtext( 'PRIMARY', side=2, line=3, cex=1.5 )
dev.off()

## find the number with no matches at all ...
matched_ids <- unique( final_rel$radio_ID )
total_percentage_with_matches <- length( matched_ids ) / n_radio_sources * 100.

unmatched <- my_radio_cat[ which( ! (my_radio_cat$Source_id %in% matched_ids ) ),  ]
matched <- my_radio_cat[ which( (my_radio_cat$Source_id %in% matched_ids ) ),  ]

write.table( unmatched, file='unmatched_sources.csv', row.names=FALSE, quote=FALSE, sep=',' )
write.table( matched, file='matched_sources.csv', row.names=FALSE, quote=FALSE, sep=',' )
