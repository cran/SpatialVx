enlarge_matrix <-
function (f, d) 
{
    fenl = matrix(0, nrow(f) + 2 * d, ncol(f) + 2 * d, byrow = FALSE)
    fenl[(1 + d):(nrow(f) + d), (1 + d):(ncol(f) + d)] = f
    return(fenl)
}
