library('astro')

apply_masked_regions <- function( df ){

    df_cols <- colnames( df )

    if ( 'RA' %in% df_cols ){
        ra_name <- 'RA'
        dec_name <- 'DEC'
    } else {
        ra_name <- 'ALPHA_J2000'
        dec_name <- 'DELTA_J2000'
    }

    ## get the mask definition information
    area_def_file <- '/vardy/leah/data/xmmlss/mask_definition/Kband_area_definition.dat'
    halo_regions <- '/vardy/leah/data/xmmlss/mask_definition/ds9.reg'

    area_def <- read.table( area_def_file, stringsAsFactors=FALSE, header=TRUE )

    include_area_x <- area_def$include_x
    include_area_y <- area_def$include_y 

    ## trim down to only the 'include' area
    include_index_x <- which( df[ , ra_name ] > min( include_area_x ) & df[ , ra_name ] < max( include_area_x ) )
    include_index_y <- which( df[ , dec_name ] > min( include_area_y ) & df[ , dec_name ] < max( include_area_y ) )
    include_index <- intersect( include_index_x, include_index_y )

    ## loop through the exclude regions
    area_def_cols <- colnames( area_def )
    excols <- area_def_cols[ which( grepl( 'exclude_x', area_def_cols ) )]

    exclude_index <- c()
    for ( ii in 1:length( excols ) ){
        excl_x <- area_def[ , paste( 'exclude_x', ii, sep='' ) ]
        excl_y <- area_def[ , paste( 'exclude_y', ii, sep='' ) ]
        excl_index_x <- which( df[ , ra_name ] > min( excl_x ) & df[ , ra_name ] < max( excl_x ) )
        excl_index_y <- which( df[ , dec_name ] > min( excl_y ) & df[ , dec_name ] < max( excl_y ) )
        excl_index <- intersect( excl_index_x, excl_index_y )
        exclude_index <- c( exclude_index, excl_index )
    }
    exclude_index <- unique( exclude_index )

    ## find where things are in include but not exclude.
    include_in_exclude <- ( include_index %in% exclude_index )

    final_index <- include_index[ which( !include_in_exclude ) ]

    new_df <- df[ final_index, ]

    ## now exclude the halo regions
    ## first read the file
    ds9_reg <- readLines( halo_regions )
    ds9_reg <- ds9_reg[3:length(ds9_reg)]

    exclude_index <- c()
    for ( ds9r in ds9_reg ){
        tmp <- strsplit( ds9r, '(', fixed=TRUE )[[1]]
        if ( tmp[1] == 'circle' ){
            x_cen <- as.numeric( strsplit( tmp[2], ',' )[[1]][1] )
            y_cen <- as.numeric( strsplit( tmp[2], ',' )[[1]][2] )
            radius <- as.numeric( substr( strsplit( tmp[2], ',' )[[1]][3], 1, 6 ) ) 
            distances <- cenang( x_cen, y_cen, new_df[,ra_name], new_df[,dec_name] ) * 60 * 60 ## convert to arcsec
            excl_index <- which( distances < radius )
            exclude_index <- c( exclude_index, excl_index )

        } #endif
    } #endfor ds9r
    exclude_index <- unique( exclude_index )

    tmp_index <- seq( 1, dim( new_df )[1] )

    circle_index <- tmp_index[ which( !( tmp_index %in% exclude_index ) ) ]

    ## deal with polygons --which are actually only squares
    exclude_index <- c()
    for ( ds9r in ds9_reg ){
        tmp <- strsplit( ds9r, '(', fixed=TRUE )[[1]]
        if ( tmp[1] == 'polygon' ){
            tmp1 <- strsplit( tmp[2], ',' )[[1]]
            tmp1[length(tmp1)] <- substr( tmp1[length(tmp1)], 1, nchar(tmp1[length(tmp1)])-1 )
            tmp1 <- as.numeric( tmp1 )
            point1 <- list( x=tmp1[1], y=tmp1[2] )
            point2 <- list( x=tmp1[3], y=tmp1[4] )
            point3 <- list( x=tmp1[5], y=tmp1[6] )
            point4 <- list( x=tmp1[7], y=tmp1[8] )

            ## line between point 1 and point 2
            p12_m <- ( point2$y - point1$y ) / ( point2$x - point1$x )
            p12_b <- point1$y - p12_m*point1$x 
                
            my_y <- new_df[,ra_name] * p12_m + p12_b 
            excl_index <- which( new_df[,dec_name] < my_y )

            ## line between point 3 and point 4
            p43_m <- ( point3$y - point4$y ) / ( point3$x - point4$x )
            p43_b <- point4$y - p43_m*point4$x 

            my_y <- new_df[,ra_name] * p43_m + p43_b 
            excl_index <- intersect( excl_index, which( new_df[,dec_name] > my_y ) )

            ## line between point 1 and point 4
            p41_m <- ( point1$y - point4$y ) / ( point1$x - point4$x )
            p41_b <- point4$y - p41_m*point4$x 

            my_y <- new_df[,ra_name] * p41_m + p41_b 
            excl_index <- intersect( excl_index, which( new_df[,dec_name] < my_y ) )

            ## line between point 2 and point 3
            p32_m <- ( point2$y - point3$y ) / ( point2$x - point3$x )
            p32_b <- point3$y - p32_m*point3$x 
                
            my_y <- new_df[,ra_name] * p32_m + p32_b 
            excl_index <- intersect( excl_index, which( new_df[,dec_name] > my_y ) )

            exclude_index <- c( exclude_index, excl_index )

        }
    }
    exclude_index <- unique( exclude_index )

    polygon_index <- tmp_index[ which( !(tmp_index %in% exclude_index ) ) ]

    final_include_index <- intersect( circle_index, polygon_index )
        
    df_out <- new_df[ final_include_index, ]

    return( df_out )

}


















