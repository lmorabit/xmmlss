############## READ IN VISUAL MATCH DATA
visual_matches_cat <- '/vardy/leah/data/xmmlss/visual_matches/eyeball_out.dat'
vm_widths <- c( 16, 3, 8, 14, 3, 7, 200 )
visual_matches <- read.fwf( visual_matches_cat, widths=vm_widths, skip=1, stringsAsFactors=FALSE )
tmp <- readLines( visual_matches_cat, n=1 )
colnames( visual_matches ) <- strsplit( tmp, '\t' )[[1]]

############## READ IN VLA DATA
## add the corresponding source id column
## read in the un-combined VLA data
vla_file <- '/vardy/leah/data/xmmlss/VLA/13B-308_080916_csv.txt'
vla_data <- read_VLA_data( vla_file )

## add the Source_id to the visual match data
tmp <- c()
for ( jid in visual_matches$name ) tmp <- c( tmp, unique( vla_data$Source_id[ which( vla_data$ID == jid ) ] ) )
visual_matches[ ,'Source_id'] <- tmp

## NO VISUAL MATCHES
no_visual_matches <- visual_matches[ which( visual_matches$match_ID == 0 ), ]

## weird things?
weird_visual_matches <- visual_matches[ which( visual_matches$match_ID < 300000 ), ]

## visual matches
visual_matches <- visual_matches[ which( visual_matches$match_ID > 300000 ), ]

## unresolved sources
visual_unresolved_index <- which( trimws(visual_matches$FR_class) == '0' )

visual_matches <- visual_matches[ visual_unresolved_index, ] 

## there are still a few Source_ids that are repeated .... remove these
repeated_ids <- unique( visual_matches$Source_id[ which( duplicated( visual_matches$Source_id ) )] )

no_repeat_index <- which( ! visual_matches$Source_id %in% repeated_ids )

visual_matches <- visual_matches[ no_repeat_index, ]


############ TANGENTIAL:  check that the 'match_ID' and 'Source_id' groupings are the same
check_cols <- paste( visual_matches$match_ID, visual_matches$Source_id, sep=' ' )
mismatch_index <- c()
for ( ii in 1:dim( visual_matches )[1] ){

    jname <- visual_matches$name[ii]
    match_name_ind <- which( visual_matches$match_ID == visual_matches$match_ID[ii] ) 
    source_id_ind <- which( visual_matches$Source_id == visual_matches$Source_id[ii] )
    ## check if the indices are the same length
    if ( length( match_name_ind ) == length( source_id_ind ) ){
        ind_check <- match_name_ind - source_id_ind
        if ( max( ind_check ) > 0 ) mismatch_index <- c( mismatch_index, ii )
    } else mismatch_index <- c( mismatch_index, ii )
}
#valid_index <- which( ! seq(1,dim(visual_matches)[1]) %in% mismatch_index ) 
#visual_matches <- visual_matches[ valid_index, ]

## NOW all the Source_ids are unique.  compare with the final matches.

############### READ IN/PREPARE FINAL LR MATCH DATA
final_matches <- read.table( 'final_matches.csv', stringsAsFactors=FALSE, sep=',', header=TRUE )



rel_threshold <- 0.8
rel_index <- which( final_matches$lr_rel > rel_threshold )
final_rel <- final_matches[ rel_index, ]
## pick the best matches so every pair is unique
good_index <- c()
for ( ii in 1:dim(final_rel)[1] ){
    r_id_index <- which( final_rel$radio_ID == final_rel$radio_ID[ii] )
    if ( length( r_id_index ) == 1 ){
        good_index <- c( good_index, ii )
    } else {
        ## find the best match
        best_match_index <- which( final_rel$lr_rel[r_id_index] == max( final_rel$lr_rel[r_id_index] ) )
        good_index <- c( good_index, r_id_index[ best_match_index ] )
    }
}
good_index <- good_index[ which( !duplicated( good_index ) ) ]
final_lr_matches <- final_rel[ good_index, ]


## FIRST CONSIDER ONLY SOURCE IDS THAT ARE IN THE VISUAL MATCH CATALOGUE
visual_index <- which( final_matches$radio_ID %in% visual_matches$Source_id )
## this narrows it down to 571 sources
## 2007 sources

final_matches_visual <- final_matches[ visual_index, ]
rel_threshold <- 0.8
rel_index <- which( final_matches_visual$lr_rel > rel_threshold )
## this narrows it down to 190 sources
final_matches_visual_rel <- final_matches_visual[ rel_index, ]

## get rid of visual matches that are not in this final thing
final_visual_index <- which( visual_matches$Source_id %in% final_matches_visual_rel$radio_ID )
visual_matches <- visual_matches[ final_visual_index, ]

## compare!
v_rad <- visual_matches$Source_id
v_nir <- visual_matches$match_ID

r_rad <- final_matches_visual_rel$radio_ID
r_nir <- final_matches_visual_rel$video_ID

## sort and subtract
v_sort_order <- order( v_rad )
v_rad <- v_rad[ v_sort_order ]
v_nir <- v_nir[ v_sort_order ]

r_sort_order <- order( r_rad )
r_rad <- r_rad[ r_sort_order ]
r_nir <- r_nir[ r_sort_order ]

## check that the radio ids are the same
radio_check <- v_rad - r_rad
if ( max( radio_check ) > 0 ) print('there is a problem!' )

mismatch_index <- which( v_nir - r_nir != 0 )
cat( 'Mismatches:\n' )
cat( v_rad[ mismatch_index ], '\n' )



############### READ IN K-BAND DATA TO GET POSITIONAL INFO
band_cat <- "cats/CFHTLS-W1_2016-04-14_fullcat_errfix_MAGAUTO_Galaxies_K.csv"
kband_dat <- read.table( band_cat, header=TRUE, stringsAsFactors=FALSE, sep=',' )

############### READ IN/PREPARE FINAL LR MATCH DATA
final_matches <- read.table( 'final_matches.csv', stringsAsFactors=FALSE, sep=',', header=TRUE )
rel_threshold <- 0.8
rel_index <- which( final_matches$lr_rel > rel_threshold )
final_rel <- final_matches[ rel_index, ]
## pick the best matches so every pair is unique
good_index <- c()
for ( ii in 1:dim(final_rel)[1] ){
    r_id_index <- which( final_rel$radio_ID == final_rel$radio_ID[ii] )
    if ( length( r_id_index ) == 1 ){
        good_index <- c( good_index, ii )
    } else {
        ## find the best match
        best_match_index <- which( final_rel$lr_rel[r_id_index] == max( final_rel$lr_rel[r_id_index] ) )
        good_index <- c( good_index, r_id_index[ best_match_index ] )
    }
}
good_index <- good_index[ which( !duplicated( good_index ) ) ]
final_lr_matches <- final_rel[ good_index, ]



## VIDEO --> RADIO
## match visual_matches$match_ID to final_matches$video_ID
common_matches_vr <- which( visual_matches$match_ID %in% final_lr_matches$video_ID )
n_common_matches_vr <- length( common_matches_vr )
cat( 'There are', n_common_matches_vr, 'VIDEO IDs in common.\n' ) 

## RADIO --> VIDEO
## match visual_matches$Source_id to final_rel$radio_ID
common_matches_rv <- which( visual_matches$Source_id %in% final_lr_matches$radio_ID )
n_common_matches_rv <- length( common_matches_rv )
cat( 'There are', n_common_matches_rv, 'Radio IDs in common.\n' ) 

## SAME IDs
same_id_index <- intersect( common_matches_vr, common_matches_rv )

## MISMATCHED (i.e., there is both an LR match and a visual match and they disagree)
mismatched_id_index <- which( !( common_matches_rv %in% common_matches_vr ) )

make_output_cat <- function( my_index, visual_dat, lr_dat, kband_dat, filename='' ){
    
    ## get the source_id and match_id from the visual catalogue
    s_id <- visual_dat$Source_id[ my_index ]
    m_id <- visual_dat$match_ID[ my_index ]

    ## get the subset of lr_dat with this source_id
    lr_index <- which( lr_dat$radio_ID %in% s_id )
    lr_dat_out <- lr_dat[ lr_index, ]
    lr_dat_out[ ,'Visual_ID'] <- m_id 

    ## find the info on the video sources in the k-band data
 
    ## write a table
    write.table( lr_dat_out, file=filename, quote=FALSE, sep=',', row.names=FALSE )


}
