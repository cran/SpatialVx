calculate_FSSvector_from_binary_fields <-
function (fbin1, fbin2, nvector) 
{
    if (ncol(fbin1) != ncol(fbin2) || nrow(fbin1) != nrow(fbin2)) {
        stop("The two input binary fields don't have the same dimensions !!!")
    }
    if (sum(fbin1 == 0) + sum(fbin1 == 1) != length(fbin1) || 
        sum(fbin2 == 0) + sum(fbin2 == 1) != length(fbin2)) {
        stop("At least one of the two fields contains a non-binary value. Only values 0 and 1 are allowed in the binary input fields.")
    }
    if (identical(nvector, round(nvector)) == FALSE || length(which(nvector%%2 != 
        1)) > 0) {
        stop("nvector can only contain odd positive integer values representing the neigberhood size!")
    }
    d = max(ncol(fbin1), nrow(fbin1))
    fbin1enl = enlarge_matrix(fbin1, d)
    fsum1enl = calculate_summed_field(fbin1enl)
    fbin2enl = enlarge_matrix(fbin2, d)
    fsum2enl = calculate_summed_field(fbin2enl)
    return(calculate_FSSvector_from_enlarged_summed_fields(fsum1enl, 
        fsum2enl, nvector, d))
}
