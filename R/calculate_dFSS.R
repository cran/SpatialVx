calculate_dFSS <-
function (fbin1, fbin2) 
{
    if (ncol(fbin1) != ncol(fbin2) || nrow(fbin1) != nrow(fbin2)) {
        stop("The two input binary fields don't have the same dimensions !!!")
    }
    if (sum(fbin1 == 0) + sum(fbin1 == 1) != length(fbin1) || 
        sum(fbin2 == 0) + sum(fbin2 == 1) != length(fbin2)) {
        stop("At least one of the two fields contains a non-binary value. Only values 0 and 1 are allowed in the binary input fields.")
    }
    d = max(ncol(fbin1), nrow(fbin1))
    fbin1enl = enlarge_matrix(fbin1, d)
    fbin2enl = enlarge_matrix(fbin2, d)
    fsum1 = calculate_summed_field(fbin1enl)
    fsum2 = calculate_summed_field(fbin2enl)
    dimx = ncol(fsum1)
    dimy = nrow(fsum1)
    sum1 = fsum1[dimy, dimx]
    sum2 = fsum2[dimy, dimx]
    if (sum1 == 0 || sum2 == 0) {
        stop("One of the two input binary fields has no precipation (all values are 0) !!!")
    }
    bias = max(sum1/sum2, sum2/sum1)
    if (bias > 2) {
        stop("The freqeuncy bias for the two input binary field is larger than 2. This is not allowed since it will produce unrealistic results for dFSS. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
    }
    if (bias > 1.5) {
        warning("The freqeuncy bias for the two input binary field is larger than 1.5. In some cases this might result in unrealistic results for the dFSS. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
    }
    FSSn1 = calculate_FSS_from_enlarged_summed_fields(fsum1, 
        fsum2, 1, d)
    dFSS = 0
    if (FSSn1 == 1) {
        dFSS = 0
    }
    else {
        funion = fbin1 * fbin2
        fbin1or = fbin1 - funion
        fbin2or = fbin2 - funion
        fbin1orenl = enlarge_matrix(fbin1or, d)
        fbin2orenl = enlarge_matrix(fbin2or, d)
        fsum1or = calculate_summed_field(fbin1orenl)
        fsum2or = calculate_summed_field(fbin2orenl)
        sum1or = fsum1or[dimy, dimx]
        sum2or = fsum2or[dimy, dimx]
        if (sum1or == 0 || sum2or == 0) {
            stop("One of the two fields with overlap removed has no precipation (all values are 0) !!!. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
        }
        biasor = max(sum1or/sum2or, sum2or/sum1or)
        if (biasor > 2) {
            stop("The freqeuncy bias for the two fields with overlap removed is larger than 2. This is not allowed since it will produce unrealistic results for dFSS. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
        }
        if (biasor > 1.5) {
            warning("The freqeuncy bias for the two fields with overlap removed is larger than 1.5. In some cases this might result in unrealistic results for the dFSS. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
        }
        n_FSS05 = calculate_n05_using_bisection_from_the_summed_fields(fsum1or, 
            fsum2or, d)
        dFSS = (1 - FSSn1) * floor(n_FSS05/2)
    }
    return(dFSS)
}
