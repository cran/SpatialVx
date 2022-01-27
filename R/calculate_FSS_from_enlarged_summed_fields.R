calculate_FSS_from_enlarged_summed_fields <-
function (fsum1enl, fsum2enl, n, d) 
{
    ffrac1 = calculate_fractions_from_enlarged_summed_field(fsum1enl, 
        n, d)
    ffrac2 = calculate_fractions_from_enlarged_summed_field(fsum2enl, 
        n, d)
    top = sum((ffrac1 - ffrac2)^2)
    bottom = sum(ffrac1^2 + ffrac2^2)
    if (bottom > 0) {
        return(1 - top/bottom)
    }
    else {
        return(0)
    }
}
