"checking.RT" <- function(RT, i, RT.win) {
    out <- c(i, which(abs(RT - RT[i]) > RT.win))
    out
}
