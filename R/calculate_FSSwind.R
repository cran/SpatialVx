calculate_FSSwind <-
function (findex1, findex2, nvector) 
{
    if (ncol(findex1) != ncol(findex2) || nrow(findex1) != nrow(findex2)) {
        stop("The two input wind class index fields don't have the same dimensions !!!")
    }
    max_index = max(max(findex1), max(findex2))
    min_value = min(min(findex1), min(findex2))
    if (identical(findex1, round(findex1)) == FALSE || identical(findex2, 
        round(findex2)) == FALSE || min_value < 1) {
        stop("At least one of the two index fields contains non-integer or non-positive values. Only positive integer values that represent wind classes are allowed in the input index fields! The smallest allowed windclass index is 1 and not 0.")
    }
    if (identical(nvector, round(nvector)) == FALSE || length(which(nvector%%2 != 
        1)) > 0) {
        stop("nvector can only contain odd positive integer values representing the neigberhood size!")
    }
    fsum1enl_list = list()
    fsum2enl_list = list()
    d = max(ncol(findex1), nrow(findex1))
    for (iclass in 1:max_index) {
        ftemp = findex1
        ftemp[ftemp != iclass] <- 0
        ftemp[ftemp == iclass] <- 1
        ftemp2 = enlarge_matrix(ftemp, d)
        fsum1enl_list[[iclass]] <- calculate_summed_field(ftemp2)
        ftemp = findex2
        ftemp[ftemp != iclass] <- 0
        ftemp[ftemp == iclass] <- 1
        ftemp2 = enlarge_matrix(ftemp, d)
        fsum2enl_list[[iclass]] <- calculate_summed_field(ftemp2)
    }
    FSSwind = nvector
    for (inn in 1:length(nvector)) {
        numerator = 0
        denominator = 0
        for (iclass in 1:max_index) {
            ffrac1 <- calculate_fractions_from_enlarged_summed_field(fsum1enl_list[[iclass]], 
                nvector[inn], d)
            ffrac2 <- calculate_fractions_from_enlarged_summed_field(fsum2enl_list[[iclass]], 
                nvector[inn], d)
            numerator = numerator + sum((ffrac1 - ffrac2)^2)
            denominator = denominator + sum(ffrac1^2 + ffrac2^2)
        }
        FSSwind[inn] = 1 - numerator/denominator
    }
    return(FSSwind)
}
