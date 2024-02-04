pred_comb_mediator <- function(cont, BM, MCdata_trt, MCdata_control, IT){
  if(cont == FALSE) {
    x_hat <- pnorm(pwbart_it(MCdata_control, BM$treedraws, mu = BM$binaryOffset, it=IT)) #make predictions conditional on counterfactual data
    new_MCdata <- cbind(MCdata_trt, rbinom(length(x_hat), 1, prob = x_hat))
  } else if(cont == TRUE) {
    x_hat <- pwbart_it(MCdata_control, BM$treedraws, mu = BM$mu, it = IT)
    new_MCdata <- cbind(MCdata_trt, rnorm(length(x_hat), mean = x_hat, sd = mean(BM$sigma)))
  }
  return(new_MCdata)
}