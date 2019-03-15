# function to sample from a defined uniform probability distribution
"multsample" <- function(P) {
    o <- sample(1:length(P), size = 1, replace = TRUE, prob = P)
    o
}
