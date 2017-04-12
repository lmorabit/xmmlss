library('astro')

box_area <- function( xvec, yvec ){

    area <- ( max( xvec ) - min( xvec ) ) * ( max( yvec ) - min( yvec ) )
    return( area )
}

get_video_area_arcsec <- function(){

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

    ## convert to arcsec
    ad <- area_def * 60 * 60

    main_area <- box_area( ad$include_x, ad$include_y )

    ## subtract off bottom area
    bottom_area <- box_area( ad$exclude_x1, ad$exclude_y1 ) + box_area( c(max(ad$exclude_x2)-max(ad$exclude_x1)), ad$exclude_y2 ) + box_area( ad$exclude_x3, c(max(ad$exclude_y3),max(ad$exclude_y1)) ) + box_area( ad$exclude_x4, c(max(ad$exclude_y4),max(ad$exclude_y2)) )

    ## subtract off top area
    top_area <- box_area( ad$exclude_x5, ad$exclude_y5 ) + box_area( ad$exclude_x6, c(min(ad$exclude_y6),min(ad$exclude_y5)) ) + box_area( ad$exclude_x7, c(min(ad$exclude_y6),min(ad$exclude_y7)) )

    ## interior areas
    interior_areas <- box_area( ad$exclude_x8, ad$exclude_y8 ) + box_area( ad$exclude_x9, ad$exclude_y9 ) + box_area( ad$exclude_x10, ad$exclude_y10 )
    
    main_area <- main_area - bottom_area - top_area - interior_areas
    
    ## now exclude the halo regions
    ## first read the file
    ds9_reg <- readLines( halo_regions )
    ds9_reg <- ds9_reg[3:length(ds9_reg)]

    circle_areas <- c()
    for ( ds9r in ds9_reg ){
        tmp <- strsplit( ds9r, '(', fixed=TRUE )[[1]]
        if ( tmp[1] == 'circle' ){
            radius <- as.numeric( substr( strsplit( tmp[2], ',' )[[1]][3], 1, 6 ) ) ## this is in arcsec
            circle_areas <- c( circle_areas, pi*radius^2 )
        } #endif
    } #endfor ds9r
    circle_area <- sum( circle_areas )

    main_area <- main_area - circle_area

    ## deal with polygons --which are actually only squares
    polygon_areas <- c()
    for ( ds9r in ds9_reg ){
        tmp <- strsplit( ds9r, '(', fixed=TRUE )[[1]]
        if ( tmp[1] == 'polygon' ){
            tmp1 <- strsplit( tmp[2], ',' )[[1]]
            tmp1[length(tmp1)] <- substr( tmp1[length(tmp1)], 1, nchar(tmp1[length(tmp1)])-1 )
            tmp1 <- as.numeric( tmp1 ) * 60 * 60
            point1 <- list( x=tmp1[1], y=tmp1[2] )
            point2 <- list( x=tmp1[3], y=tmp1[4] )
            point3 <- list( x=tmp1[5], y=tmp1[6] )
            point4 <- list( x=tmp1[7], y=tmp1[8] )

            mat1 <- as.matrix( cbind(c(point1$x, point1$y ), c(point2$x, point2$y)) )
            mat2 <- as.matrix( cbind(c(point2$x, point2$y ), c(point3$x, point3$y)) )
            mat3 <- as.matrix( cbind(c(point3$x, point3$y ), c(point4$x, point4$y)) )
            mat4 <- as.matrix( cbind(c(point4$x, point4$y ), c(point1$x, point1$y)) )
    
            sum_dets <- 0.5 * ( det( mat1 ) + det( mat2 ) + det( mat3 ) + det( mat4 ) )

            polygon_areas <- c( polygon_areas, sum_dets )

        }
    }

    polygon_area <- sum( polygon_areas )
    
    total_area <- main_area - polygon_area

    return( total_area )

}


















