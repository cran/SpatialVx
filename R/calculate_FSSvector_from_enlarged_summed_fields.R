calculate_FSSvector_from_enlarged_summed_fields <-
function (fsum1enl, fsum2enl, nvector, d) 
{
    FSSvector = nvector
    for (inn in 1:length(nvector)) FSSvector[[inn]] = calculate_FSS_from_enlarged_summed_fields(fsum1enl, 
        fsum2enl, nvector[inn], d)
    return(FSSvector)
}
