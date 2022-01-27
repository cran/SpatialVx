calculate_n05_using_bisection_from_the_summed_fields <-
function (fsum1or, fsum2or, d) 
{
    dimx = ncol(fsum1or)
    dimy = nrow(fsum1or)
    rows = nrow(fsum1or) - 2 * d
    cols = ncol(fsum1or) - 2 * d
    x1 = 1
    x2 = 2 * max(rows, cols) - 1
    if (x2 < 5) {
        stop("Domain size needs to be at least 2 grid points")
    }
    FSS1 = calculate_FSS_from_enlarged_summed_fields(fsum1or, 
        fsum2or, x1, d)
    FSS2 = calculate_FSS_from_enlarged_summed_fields(fsum1or, 
        fsum2or, x2, d)
    if (FSS1 > 0.5) 
        return(1)
    if (FSS2 <= 0.5) {
        stop("FSS does never reach value 0.5. There is something wrong.")
    }
    repeat {
        xnew = (x1 + x2)/2
        if (xnew%%2 == 0) {
            xnew = xnew + 1
        }
        FSSnew = calculate_FSS_from_enlarged_summed_fields(fsum1or, 
            fsum2or, xnew, d)
        if (FSSnew > 0.5) {
            x2 = xnew
            FSS2 = FSSnew
        }
        else {
            x1 = xnew
            FSS1 = FSSnew
        }
        if (x2 - x1 <= 2) {
            break
        }
    }
    return(x2)
}
