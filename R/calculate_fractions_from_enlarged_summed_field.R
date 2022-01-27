calculate_fractions_from_enlarged_summed_field <-
function (fsumenl, n, d) 
{
    if (n%%2 != 1) {
        stop("Neighbourhood size can only be an odd positive integer")
    }
    dr = (n - 1)/2
    rows = nrow(fsumenl) - 2 * d
    cols = ncol(fsumenl) - 2 * d
    if (n > 2 * max(rows, cols) - 1) {
        return(calculate_fractions_from_enlarged_summed_field(fsumenl, 
            2 * max(rows, cols) - 1, d))
    }
    f_topright = fsumenl[(d + 1 + dr):(d + rows + dr), (d + 1 + 
        dr):(d + cols + dr)]
    f_bottomleft = fsumenl[(d + 1 - dr - 1):(d + rows - dr - 
        1), (d + 1 - dr - 1):(d + cols - dr - 1)]
    f_topleft = fsumenl[(d + 1 + dr):(d + rows + dr), (d + 1 - 
        dr - 1):(d + cols - dr - 1)]
    f_bottomright = fsumenl[(d + 1 - dr - 1):(d + rows - dr - 
        1), (d + 1 + dr):(d + cols + dr)]
    ffrac = f_topright - f_topleft - f_bottomright + f_bottomleft
    ffrac = ffrac/(n * n)
    return(ffrac)
}
