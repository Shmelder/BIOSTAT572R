
test.once <- function(data, method = "lr"){
  if ( method == "lr"){
    ful.lm <- lm(data[, 1] ~ data[, -1])
    p.val <- anova(ful.lm)$`Pr(>F)`[1]
  }
  return(as.numeric(p.val < 0.05))
}
