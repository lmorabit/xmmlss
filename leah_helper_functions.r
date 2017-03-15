find_resolved_index <- function( total_flux, peak_flux, total_flux_error, peak_flux_error, p2c=1 ){

        peak_flux <- peak_flux * p2c
        peak_flux_error <- peak_flux_error * p2c

        peak_to_total_ratio <- peak_flux / total_flux
        peak_to_total_error <- peak_to_total_ratio * sqrt( (peak_flux_error / peak_flux)^2. + (total_flux_error / total_flux)^2. )
        resolved_index <- which( abs(peak_to_total_ratio - 1) > peak_to_total_error )
        return( resolved_index )
}

find_unresolved_index <- function( total_flux, peak_flux, total_flux_error, peak_flux_error, p2c=1 ){

        peak_flux <- peak_flux * p2c
        peak_flux_error <- peak_flux_error * p2c

        peak_to_total_ratio <- peak_flux / total_flux
        peak_to_total_error <- peak_to_total_ratio * sqrt( (peak_flux_error / peak_flux)^2. + (total_flux_error / total_flux)^2. )
        unresolved_index <- which( abs(peak_to_total_ratio - 1 ) <= peak_to_total_error )
        return( unresolved_index )
}


erf <- function(x) 2 * pnorm( x * sqrt(2) ) - 1
erfc <- function(x) 2 * pnorm( x * sqrt(2), lower = FALSE )

