calculate_summed_field <-
function (f) 
{
    fsum1 <- f
    fsum2 <- f
    for (k in 1:nrow(f)) {
        fsum1[k, ] <- cumsum(f[k, ])
    }
    for (k in 1:ncol(f)) {
        fsum2[, k] <- cumsum(fsum1[, k])
    }
    return(fsum2)
}
