source('/vardy/leah/xmmlss/plotting_libraries.r')
mycols <- viridis( 5 )
stilts_exec <- '/vardy/leah/stilts/stilts'
master_cat <- 'cats/CFHTLS-W1_2016-04-14_fullcat_errfix_MAGAUTO.fits'


star_galaxy_separation <- function( master_cat ){

    ## --- housekeeping: index appropriately

    ## get the subset of the master catalogue with good SNR and no HALOFLAG
    ## and FLAGS < 5
    ss <- paste( stilts_exec, " tpipe in=", master_cat, " out=cats/GRJK_star_galaxy.csv omode=out ofmt=csv cmd=\'select J_SNR>=5\' cmd=\'select K_SNR>=5\' cmd=\'select G_SNR>=5\' cmd=\'select R_SNR>=5\' cmd=\'select J_FLAGS<5\' cmd=\'select K_FLAGS<5\' cmd=\'select G_FLAGS<5\' cmd=\'select R_FLAGS<5\' cmd=\'select HALOFLAG==0\'", sep='' )
    system( ss )

    ## read in the catalogue
    my_dat <- read.table( 'cats/GRJK_star_galaxy.csv', header=TRUE, stringsAsFactors=FALSE, sep=',' )
  
    ## find the variables
    g_minus_r <- my_dat$G_MINUS_R
    J_minus_K <- my_dat$J_MINUS_K

    ## make a plot
    x_lim <- c(-0.25,2)
#    y_lim <- c(-1,2.5)
    y_lim <- c(-1,1.5)
    pdf( 'sg_separation.pdf' )
    #lplot( g_minus_r, J_minus_K , pch='.', xlim=x_lim, ylim=y_lim, x_lab='g-i [mag]', y_lab='J-K [mag]', col=tcol(mycols[3],trans=100) )
    lplot( g_minus_r, J_minus_K , pch=16, xlim=x_lim, ylim=y_lim, x_lab='g-r [mag]', y_lab='J-K [mag]', col=tcol(mycols[3],trans=100), cex=0.2 )
    # col='gray48'
    ## make a nice density plot
    stuff <- kde2d( g_minus_r, J_minus_K, n=c(200,200), lims=c(x_lim,y_lim) )
    mylevels <- pretty( stuff$z, n=20 )
    mylevels <- mylevels[ c(1,2,3,4,5,6,7,8,12,16,20) ]
    contour( stuff, add=TRUE, drawlabels=FALSE, levels=mylevels )

    ## now fit the stellar locus ...
    l1 <- 0.25  ## 0.4 g-i 
    l2 <- 1.3  ## 1.9 g-i

    ## first limit the outliers
    x_index <- which( g_minus_r >= l1 & g_minus_r <= l2 )  ## for g-i
    y_index <- which( J_minus_K >= y_lim[1] & J_minus_K <= y_lim[2] )
    good_data <- intersect( x_index, y_index )

    x <- g_minus_r[good_data]
    y <- J_minus_K[good_data]

    ## also limit the galaxy contributin
    x_div <- seq( floor(min(x)), ceiling(max(x)), 0.1 )
#    y_div <- 0.4 * x_div - 0.5  ## for g-i
    y_div <- 0.5 * x_div - 0.475
    #lines( x_div, y_div, col='red' )

    ## find where points are below this line ...
#    y_vals <- 0.4 * x - 0.5  ## for g-i
    y_vals <- 0.5 * x - 0.475
    sl_vals <- which( y <= y_vals )

    x <- x[sl_vals]
    x2 <- x^2.
    y <- y[sl_vals]
    quad_model <- lm( y ~ x + x2 )

    x_loc <- seq( l1, l2, 0.01 )
    pred_vals <- predict( quad_model, list( x=x_loc, x2=x_loc^2 ) )
    lines( x_loc, pred_vals , col='black', lwd=2.5 )
    x_up <- seq( l2, 10, 0.01 )
    lines( x_up, rep(pred_vals[which( x_loc == max(x_loc)) ], length( x_up ) ), col='black', lwd=2.5 )
    x_down <- seq( -1, l1, 0.01 )
    lines( x_down, rep( pred_vals[ which( x_loc == l1 ) ], length( x_down ) ), col='black', lwd=2.5 )

    ## shift it up
    shift_val <- 0.14
    lines( x_loc, pred_vals+shift_val , col='black', lty=2, lwd=2.5 )
    lines( x_up, rep(pred_vals[which( x_loc == max(x_loc)) ], length( x_up ) )+shift_val, col='black', lty=2, lwd=2.5 )
    lines( x_down, rep( pred_vals[ which( x_loc == l1 ) ], length( x_down ) )+shift_val, col='black', lty=2, lwd=2.5 )
    legend( 'topright', c(expression(0.14+f[locus]),expression(f[locus])), lty=c(2,1), lwd=2.5, bty='n', cex=1.5 )
    dev.off()

    ## write a file with the results
    c1 <- pred_vals[ which( x_loc == l1 ) ]
    c2 <- quad_model$coefficients[[1]]
    b2 <- quad_model$coefficients[[2]]
    a2 <- quad_model$coefficients[[3]]
    c3 <- pred_vals[which( x_loc == max(x_loc)) ]

    df <- data.frame( c1, c2, b2, a2, c3 )
    write.table( df, file='sg_separation_values.txt', row.names=FALSE, quote=FALSE )

    cat( 'f_locus =', c1, ' ( x <', l1, ')\n' )
    cat( 'f_locus =', c2,'+', b2,'x',a2,'x^2  (0.4 < x < 1.9 )\n' )
    cat( 'f_locus =', c3, ' ( x >', l2, ')\n' )

    ## and output a catalogue of galaxies
    locus_values <- c2 + b2 * g_minus_r + a2 * g_minus_r^2 
    locus_values[ which( g_minus_r < l1 ) ] <- c1 
    locus_values[ which( g_minus_r > l2 ) ] <- c3 
    
    add_offset <- 0.14
    star_index <- which( J_minus_K < (locus_values+add_offset) )keepc
    gs_class <- rep( 0, length( J_minus_K ) )
    gs_class[ star_index ] <- 1
    
    gs_df <- data.frame( my_dat$ID, gs_class )
    colnames( gs_df ) <- c( 'SID', 'GS_CLASS' )
    write.table( gs_df, file='GSCLASS.csv', row.names=FALSE, sep=',', quote=FALSE )

    ## also get things with CLASS_STAR > 0.97 and K_MAG_AUTO < 20 
    #ss <- paste( stilts_exec, " tpipe in=", master_cat, " out=cats/CLASS_STAR_star_galaxy.csv omode=out ofmt=csv cmd=\'select K_CLASS_STAR>0.97\' cmd=\'keepcols ", '"ID"\'', sep='' )  ## 26980
    ss <- paste( stilts_exec, " tpipe in=", master_cat, " out=cats/CLASS_STAR_star_galaxy.csv omode=out ofmt=csv cmd=\'select K_CLASS_STAR>0.97\' cmd=\'select K_MAG_AUTO<20\' cmd=\'keepcols ", '"ID"\'', sep='' )  ## 12691
    system( ss )
    

}
