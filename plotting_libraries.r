##################################
##	Plot a fits image	##
##				##
## written by Leah K. Morabito  ##
##    lmorabit@gmail.com	##
##################################


require("viridis")
require('FITSio')
require('astro')
require('celestial')
require('fields')
library('RColorBrewer')

flux_density_Jy_per_beam <<- expression('Flux density [Jy bm'^-1*']')
lum_erg_per_sec_hz <<- expression('Luminosity [erg s'^-1~'Hz'^-1*']')
log_nuLnu_erg_per_sec <<- expression('log('*nu*'Luminosity/erg s'^-1*')')
power_watts_per_hz <<- expression('Power [W'~'Hz'^-1*']')
RA_hms_label <<- "RA [hh:mm:ss]"
Dec_dms_label <<- "Dec [dd:mm:ss]"


tcol <- function(color, trans = 100) {

  if (length(color) != length(trans) & 
        !any(c(length(color), length(trans)) == 1)) 
    stop('Vector lengths not correct')
  if (length(color) == 1 & length(trans) > 1) 
    color <- rep(color, length(trans))
  if (length(trans) == 1 & length(color) > 1) 
    trans <- rep(trans, length(color))

  res <- paste0('#', apply(apply(rbind(col2rgb(color)), 2, function(x) 
    format(as.hexmode(x), 2)), 2, paste, collapse = ''))
  res <- unlist(unname(Map(paste0, res, as.character(as.hexmode(trans)))))
  res[is.na(color)] <- NA
  return(res)
}

format_log_axis_labels <- function( vals ){

	my_labels <- rep( "", length( vals ) )
	for ( rr in 1:length(vals) ){
		if ( vals[rr] < 1 ){
			my_labels[rr] <- format( vals[rr], nsmall=1, digits=1 )
		} else {
			my_labels[rr] <- as.character( round( vals[rr] ) )
		}
	}
	return( my_labels )
}

get_legend_position <- function( xaxp, yaxp, l_pos='br', xlog=FALSE, ylog=FALSE ) {

	x_range <- xaxp[2] - xaxp[1]
	y_range <- yaxp[2] - yaxp[1]
	x_min <- xaxp[1]
	y_min <- yaxp[1]
	if (l_pos == 'br'){
		if ( xlog ) x_pos <- 10^( 0.5 * log10(x_range) + log10(x_min) ) else x_pos <- 0.7 * x_range + x_min
		if ( ylog ) y_pos <- 0.1 * y_range + y_min else y_pos <- 0.1 * y_range + y_min
	}
	if (l_pos == 'tl'){
	        if ( xlog ) x_pos <- 10.^( 0.02 * log10(x_range) + log10(x_min) ) else x_pos <- 0.8 * x_range + x_min
                if ( ylog ) y_pos <- 0.1 * y_range + y_min else y_pos <- 0.8 * y_range + y_min
	}
	if (l_pos == 'tr'){
                if ( xlog ) x_pos <- 10.^( 0.02 * log10(x_range) + log10(x_min) ) else x_pos <- 0.8 * x_range + x_min
                if ( ylog ) y_pos <- 0.1 * y_range + y_min else y_pos <- 0.8 * y_range + y_min
	}
	if (l_pos == 'bl'){
                if ( xlog ) x_pos <- 10.^( 0.02 * log10(x_range) + log10(x_min) ) else x_pos <- 0.8 * x_range + x_min
                if ( ylog ) y_pos <- 0.1 * y_range + y_min else y_pos <- 0.8 * y_range + y_min
	}
	return( data.frame( x=x_pos, y=y_pos ) )

}

get_aips_beam_size <- function( fitsheader ){

        infoline <- fitsheader[ grep( "BMAJ", fitsheader ) ]
        tmp <- strsplit( infoline, " " )
        bmaj <- as.numeric(tmp[[1]][8])
        bmin <- as.numeric(tmp[[1]][11])
        bpa <- as.numeric(tmp[[1]][13])

        beam <- data.frame( bmaj, bmin, bpa )
        return(beam)
}

get_color_scale <- function ( colorscale ){

	if (colorscale == "linear"){
		mycols <- viridis(256)
	}
	if (colorscale == "log"){
		tmpcols <- viridis(2560)
		a <- 100.
		xseq <- seq(0,1,1/256)
		yvals <- log10( a * xseq + 1. )/log10( a )
		tmpind <- (yvals - min(yvals)) / max( yvals - min(yvals)) * 2560 
		mycols <- tmpcols[tmpind]
	}
	if (colorscale == "pow"){
		tmpcols <- viridis(2560)
		a <- 100.
		xseq <- seq(0,1,1/256)
		yvals <- (a^xseq - 1.)/a
		tmpind <- (yvals - min(yvals)) / max( yvals - min(yvals)) * 2560
		mycols <- tmpcols[tmpind]
	}
	if (colorscale == "sqrt"){
		tmpcols <- viridis(2560)
		xseq <- seq(0,1,1/256)
		yvals <- sqrt(xseq)
                tmpind <- (yvals - min(yvals)) / max( yvals - min(yvals)) * 2560
                mycols <- tmpcols[tmpind]
        }
	if (colorscale == "square"){
                tmpcols <- viridis(2560)
                xseq <- seq(0,1,1/256)
                yvals <- xseq^2.
                tmpind <- (yvals - min(yvals)) / max( yvals - min(yvals)) * 2560
                mycols <- tmpcols[tmpind]
        }
        if (colorscale == "asinh"){
                tmpcols <- viridis(2560)
                xseq <- seq(0,1,1/256)
                yvals <- asinh( 10 * xseq ) / 3.
                tmpind <- (yvals - min(yvals)) / max( yvals - min(yvals)) * 2560
                mycols <- tmpcols[tmpind]
        }
	        if (colorscale == "sinh"){
                tmpcols <- viridis(2560)
                xseq <- seq(0,1,1/256)
                yvals <- sinh( 3. * xseq ) / 10.
                tmpind <- (yvals - min(yvals)) / max( yvals - min(yvals)) * 2560
                mycols <- tmpcols[tmpind]
        }
	

	return(mycols)
}

get_log_tick_vals <- function( mypar ){

	if ( is.finite( log10(mypar[1]) ) ) parseq <- seq(log10(mypar[1]), log10(mypar[2]) ) else parseq <- seq( 0, log10(mypar[2] ) )
	tickseq <- seq(0,9)
	tick_vals <- numeric()
	for (pars in parseq) tick_vals <- c( tick_vals, 10^(pars) * tickseq )
	return( tick_vals )
}

get_scale_range <- function( my_imdat, zscale, maxvalue, mycols ){

	if (zscale){	
		print( "Setting zscale." )
		med_val <- median( my_imdat, na.rm=TRUE )
		maxval <- 0.0005 * max( my_imdat, na.rm=TRUE )
		minval <- -0.00025 * max( my_imdat, na.rm=TRUE )
		#print( minval )
		#print( maxval )
		my_index <- which( my_imdat < minval , arr.ind=TRUE )
		my_imdat[ my_index ] <- minval
		my_index <- which( my_imdat > maxval , arr.ind=TRUE )
		my_imdat[ my_index ] <- maxval
	}

	if (maxvalue > 1){
		mymaxval <<- max(my_imdat, na.rm=TRUE)
		mycols <- mycols[1:floor(length(mycols)/maxvalue)]
	} else{
		mymaxval <<- maxvalue*max(my_imdat, na.rm=TRUE)
	}
	#print( min(my_imdat,na.rm=TRUE) )
	#print( mymaxval )

	my_index <- which( my_imdat > mymaxval, arr.ind=TRUE )
	my_imdat[ my_index ] <- mymaxval

	return( my_imdat )
}

plot_two_histograms <- function( outfile, thing1, parm1, col1, thing2, parm2, col2, mybreaks, my_xlab="Parameter", my_ylab="Number", log=FALSE ){

	## histogram of sample 1
	hist1 <- hist( thing1, breaks=mybreaks, plot=FALSE )
	hist1$counts[ hist1$counts == 0 ] <- NA
	## histogram of sample 2
	hist2 <- hist( thing2, breaks=mybreaks, plot=FALSE )
	hist2$counts[ hist2$counts == 0 ] <- NA
	

	## final values
	hist1_values <- hist1$counts 
	hist2_values <- hist2$counts 
	xvals <- hist1$mids

	## plot
	pdf(outfile)
	par(mar=c(5,5,2,2),lwd=1.5)
	## x and y limits
	my_xlim <- c( min(mybreaks), max(mybreaks) )
	minval <- min( c(min(hist1_values, na.rm=TRUE), min(hist2_values, na.rm=TRUE) ))
	maxval <- max( c(max(hist1_values, na.rm=TRUE), max(hist2_values, na.rm=TRUE) ))
	my_ylim <- c(minval, maxval)
	if ( log ){
		## transform for log plot
		scl_fac <- 1e6
		my_ylim <- log10( my_ylim * scl_fac ) 
	} 

	plot( 1, type="n", axes=FALSE, xlim=my_xlim, ylim=my_ylim, xlab="", ylab="" )
	## plot optical

	## get box width
	boxw <- ( hist1$mids[2] - hist1$mids[1] )/2.

	for ( rr in 1:length( hist1_values ) ){
		xvec <- c( xvals[rr]-boxw, xvals[rr]-boxw, xvals[rr]+boxw, xvals[rr]+boxw, xvals[rr]-boxw  )
		## plot hist1
		if ( log ) plotvals <- log10( hist1_values[rr] * scl_fac ) else plotvals <- hist1_values[rr]
		yvec <- c( 0, plotvals, plotvals, 0, 0 )
		polygon( xvec, yvec, col=col1, density=15, angle=-45, border="black" )
		## plot hist2
		if ( log ) plotvals <- log10( hist2_values[rr] * scl_fac ) else plotvals <- hist2_values[rr]
		yvec <- c( 0,  plotvals, plotvals, 0, 0 )
		polygon( xvec, yvec, col=tcol(col2,120) )		
	}


	box( which="plot" )
	yvals <- pretty( seq( par("yaxp")[1], par("yaxp")[2] ) )
	if ( log ){
		tmpy <- 10.^( yvals - log10( scl_fac ) )
		ylabs <- format_log_axis_labels( tmpy )
	} else ylabs <- format( yvals )

	axis( 1, cex.lab=1.2 )
	axis( 2, at=yvals , labels=ylabs, cex.lab=1.2 )
	mtext( my_xlab, side=1, line=3, cex=1.5 )
	mtext( my_ylab, side=2, line=3, cex=1.5 )

	legend( 'topleft', c( parm1, parm2 ), col=c("black","black"), fill=c(col1,tcol(col2,120)), angle=c(-45,45), density=c(20,100), pt.cex=2, cex=1.5, bty="n" )	

        dev.off()
	
}

plot_three_histograms <- function( outfile, thing1, parm1, col1, thing2, parm2, col2, thing3, parm3, col3, mybreaks, my_xlab="Parameter", my_ylab="Number", plog=FALSE ){

	## histogram of sample 1
	hist1 <- hist( thing1, breaks=mybreaks, plot=FALSE )
	hist1$counts[ hist1$counts == 0 ] <- NA
	## histogram of sample 2
	hist2 <- hist( thing2, breaks=mybreaks, plot=FALSE )
	hist2$counts[ hist2$counts == 0 ] <- NA
	## histogram of sample 3
	hist3 <- hist( thing3, breaks=mybreaks, plot=FALSE )
	hist3$counts[ hist3$counts == 0 ] <- NA
	

	## final values
	hist1_values <- hist1$counts 
	hist2_values <- hist2$counts 
	hist3_values <- hist3$counts
	xvals <- hist1$mids

	## plot
	pdf(outfile)
	par(mar=c(5,5,2,2),lwd=1.5)
	## x and y limits
	my_xlim <- c( min(mybreaks), max(mybreaks) )
	minval <- min( c(min(hist1_values, na.rm=TRUE), min(hist2_values, na.rm=TRUE) ))
	maxval <- max( c(max(hist1_values, na.rm=TRUE), max(hist2_values, na.rm=TRUE) ))
	my_ylim <- c(minval, maxval)
	if ( plog ){
		## transform for log plot
		scl_fac <- 1e6
		my_ylim <- log10( my_ylim * scl_fac ) 
	} 

	plot( 1, type="n", axes=FALSE, xlim=my_xlim, ylim=my_ylim, xlab="", ylab="" )
	## plot optical

	## get box width
	boxw <- ( hist1$mids[2] - hist1$mids[1] )/2.

	for ( rr in 1:length( hist1_values ) ){
		xvec <- c( xvals[rr]-boxw, xvals[rr]-boxw, xvals[rr]+boxw, xvals[rr]+boxw, xvals[rr]-boxw  )
		## plot hist1
		if ( plog ) plotvals <- log10( hist1_values[rr] * scl_fac ) else plotvals <- hist1_values[rr]
		yvec <- c( 0, plotvals, plotvals, 0, 0 )
		polygon( xvec, yvec, col=col1, density=15, angle=-45, border="black" )
		## plot hist2
		if ( plog ) plotvals <- log10( hist2_values[rr] * scl_fac ) else plotvals <- hist2_values[rr]
		yvec <- c( 0,  plotvals, plotvals, 0, 0 )
		polygon( xvec, yvec, col=tcol(col2,120) )
		## plot hist3
		if ( plog ) plotvals <- log10( hist3_values[rr] * scl_fac ) else plotvals <- hist3_values[rr]
		yvec <- c( 0,  plotvals, plotvals, 0, 0 )
		polygon( xvec, yvec, col=tcol(col3,120) )		
	}


	box( which="plot" )
	yvals <- pretty( seq( par("yaxp")[1], par("yaxp")[2] ) )
	if ( plog ){
		tmpy <- 10.^( yvals - log10( scl_fac ) )
		ylabs <- format_log_axis_labels( tmpy )
	} else ylabs <- format( yvals )

	axis( 1, cex.lab=1.2 )
	axis( 2, at=yvals , labels=ylabs, cex.lab=1.2 )
	mtext( my_xlab, side=1, line=3, cex=1.5 )
	mtext( my_ylab, side=2, line=3, cex=1.5 )

	legend( 'topleft', c( parm1, parm2, parm3 ), col=c("black","black","black"), fill=c(col1,tcol(col2,120),tcol(col3,120)), angle=c(-45,45,45), density=c(20,100,100), pt.cex=2, cex=1.5, bty="n" )	

        dev.off()
	
}

lplot <- function( x, y, x_lab='', y_lab='', leg=TRUE, ... ){

	par( mar=c(5,5,2,2) )
	plot( x, y, ..., axes=FALSE, xlab='', ylab='' )
	axis( 1 )
	axis( 2 )
	mtext( x_lab, side=1, line=3, cex=1.25 )
	mtext( y_lab, side=2, line=3, cex=1.25 )
	box( which='plot', lwd=1.5 )
	if ( !leg ) par( mar=c(5.1,4.1,4.1,2.1) )
	
}
