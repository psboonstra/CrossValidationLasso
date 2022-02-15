## ----setup, include=FALSE-----------------------------------------------------------
library(hexbin); 
library(gganimate);
library(tidyverse); 
library(broom);
library(transformr)
library(glue);
library(gifski);
library(glmnet)
library(glmnetUtils);
library(broom);

knitr::opts_knit$set(root.dir = "~/Desktop/Work/Teaching/Bios699_Winter2022/Presentations/CrossValidation");
knitr::opts_chunk$set(echo = T, warning = F, message = F);
knitr::knit_hooks$set(mysize = function(before, options, envir) {
  if (before) 
    return(options$size);
})

#knitr::opts_chunk$set(width = 10);
#knitr::opts_chunk$set(tidy.opts=list(width.cutoff=40));
recache = F;
options(digits = 3);
figure_scaler = 1/2;#1/2 for ioslides; ~1/3 for word, pdf
text_scaler = 3/3;#1 for ioslides; 2/3 for word, pdf
fig.x = 16 * figure_scaler;
fig.y = 9 * figure_scaler;
theme_set(theme_bw())



## ---- out.width = "95%",echo = F----------------------------------------------------
knitr::include_graphics("cv_example.png")


## ---- include = T, echo = F, cache = T, fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
set.seed(1);
eigenvalues = sort(0.25 + rexp(50, rate = 1));
true_beta = c(rep(0.5, 25), numeric(25));
lambda_seq = seq(0, 1.5, length = 201);
bias_sq = trace_variance = numeric(length(lambda_seq));
for(k in seq_along(lambda_seq)) {
  bias_sq[k] = lambda_seq[k] ^ 2 * sum(true_beta^2 / (eigenvalues + lambda_seq[k])^2)
  trace_variance[k] = 0.05 * sum(eigenvalues / (eigenvalues + lambda_seq[k])^2);
}
par(mar = c(3.1, 7.1, 0.1, 0.1));
plot.new();
plot.window(xlim = range(lambda_seq), ylim = c(0, max(trace_variance + bias_sq)));
axis(1, labels = NA);
axis(2, las = 1, at = c(0, trace_variance[1]), labels = c(0, expression(var(beta[MLE]))), cex.axis = 1.5);
mtext(side = 1, text = expression(lambda),  line = 1, cex = 1.5);
text(lambda_seq[150], 0.8 * bias_sq[150], labels = expression(bias(beta[R](lambda))^2), cex = 1.5)
points(lambda_seq, bias_sq, type = 'l',  lty = 2, col = "red", lwd = 5);
text(lambda_seq[150], 1.5 * trace_variance[150], labels = expression(var(beta[R](lambda))), cex = 1.5)
points(lambda_seq, trace_variance, type = 'l',  lty = 3, col = "blue", lwd = 5);
text(lambda_seq[50], 2 * trace_variance[50], labels = expression(MSE(beta[R](lambda))), cex = 1.5)
points(lambda_seq, trace_variance + bias_sq, type = 'l',  lty = 1, col = "black", lwd = 5);


## ---- include = T, echo = T, cache = T----------------------------------------------
library(glmnet);
args(cv.glmnet);


## ---- echo = F----------------------------------------------------------------------
options(digits = 3);


## ---- include = T, echo = T, cache = T----------------------------------------------
library(glmnet);
set.seed(100);
n = 200;#number obs
p = 1e3;#number covariates
rho = 0.025;#compound symmetric correlation
true_beta = c(0.3, 0.3, numeric(p - 2));
chol_x = chol(matrix(rho, nrow = p, ncol = p) + 
                diag(1 - rho, nrow = p));
sigma_y = sqrt(0.20);
x = matrix(rnorm(n * p), nrow = n) %*% chol_x;
y = x%*%true_beta + rnorm(n, sd = sigma_y);


## -----------------------------------------------------------------------------------
1 / (1 + sigma_y^2 / drop(crossprod(chol_x%*%true_beta)))


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
plot(glmnet(x, y, alpha = 0), xvar = "lambda")


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_ridge_cv = cv.glmnet(x, y, nfolds = 5, alpha = 0);
plot(ex1_ridge_cv);


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_ridge <- 
  glmnet(x, y, lambda = ex1_ridge_cv$lambda.1se, alpha = 0) %>% 
  tidy() %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "ridge", 
            partition = 1)

ggplot(ex1_ridge) + 
  geom_point(aes(x = truth,  
                 y = estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 24))


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
plot(glmnet(x, y, alpha = 1), xvar = "lambda")


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_lasso_cv = cv.glmnet(x, y, nfolds = 5, alpha = 1);
plot(ex1_lasso_cv);


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_lasso <- 
  glmnet(x, y, lambda = ex1_lasso_cv$lambda.1se, alpha = 1) %>%
  tidy(return_zeros = TRUE) %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "lasso", 
            partition = 1)

ggplot(ex1_lasso) + 
  geom_point(aes(x = truth,  
                 y = estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 24))


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
plot(glmnet(x, y, alpha = 0.5), xvar = "lambda")


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_enet_cv = cv.glmnet(x, y, nfolds = 5, alpha = 0.5);
plot(ex1_enet_cv);


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_enet <- 
  glmnet(x, y, lambda = ex1_enet_cv$lambda.1se, alpha = 0.1) %>% 
  tidy(return_zeros = TRUE) %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "enet", 
            partition = 1)

ggplot(ex1_enet) + 
  geom_point(aes(x = truth,  
                 y = estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 24))


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_ridge_cv_p2 = cv.glmnet(x, y, nfolds = 5, alpha = 0);
plot(ex1_ridge_cv_p2);


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_ridge_p2 <- 
  glmnet(x, y, lambda = ex1_ridge_cv_p2$lambda.1se, alpha = 0) %>% 
  tidy() %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "ridge", 
            partition = 2)

ggplot(ex1_ridge_p2) + 
  geom_point(aes(x = truth,  
                 y = estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 24))


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_lasso_cv_p2 = cv.glmnet(x, y, nfolds = 5, alpha = 1);
plot(ex1_lasso_cv_p2);


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_lasso_p2 <- 
  glmnet(x, y, lambda = ex1_lasso_cv_p2$lambda.1se, alpha = 1) %>%
  tidy(return_zeros = TRUE) %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "lasso", 
            partition = 2)

ggplot(ex1_lasso_p2) + 
  geom_point(aes(x = truth,  
                 y = estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 24))


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_enet_cv_p2 = cv.glmnet(x, y, nfolds = 5, alpha = 0.5);
plot(ex1_enet_cv_p2);


## ---- include = T, echo = T, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex1_enet_p2 <- 
  glmnet(x, y, lambda = ex1_enet_cv_p2$lambda.1se, alpha = 0.1) %>% 
  tidy(return_zeros = TRUE) %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "enet", 
            partition = 2)

ggplot(ex1_enet_p2) + 
  geom_point(aes(x = truth,  
                 y = estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 24))


## -----------------------------------------------------------------------------------
ex1_all_results <- 
  bind_rows(
    ex1_ridge,
    ex1_ridge_p2,
    ex1_lasso,
    ex1_lasso_p2,
    ex1_enet,
    ex1_enet_p2) %>%
  mutate(method = factor(method, levels = c("ridge", "lasso", "enet")), 
         term = fct_inorder(term) %>% as.numeric(), 
         partition = factor(partition))

ex1_all_results %>%
  group_by(method, partition) %>%
  summarize(rmse = sqrt(mean((estimate - truth)^2)))

ggplot(ex1_all_results) + 
  geom_point(aes(x = term, 
                 y = estimate, 
                 color = partition), 
             size = 0.5) + 
  facet_grid(~method) + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none", x = "none") + 
  theme(text = element_text(size = 24))


## ---- include = T, echo = T, cache = T----------------------------------------------
set.seed(100);
n = 200;#number obs
p = 1e3;#number covariates
rho = 0.4;#compound symmetric correlation
true_beta = 0.6/p + numeric(p);
chol_x = chol(matrix(rho, nrow = p, ncol = p) + 
                diag(1 - rho, nrow = p));
sigma_y = sqrt(0.20);
x = matrix(rnorm(n * p), nrow = n) %*% chol_x;
y = x%*%true_beta + rnorm(n, sd = sigma_y);


## -----------------------------------------------------------------------------------
1 / (1 + sigma_y^2 / drop(crossprod(chol_x%*%true_beta)))


## ---- include = T, echo = F, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex2_ridge_cv = cv.glmnet(x, y, nfolds = 5, alpha = 0);
plot(ex2_ridge_cv);


## ---- include = T, echo = F, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex2_ridge <- 
  glmnet(x, y, lambda = ex2_ridge_cv$lambda.1se, alpha = 0) %>% 
  tidy() %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "ridge", 
            partition = 1)

ggplot(ex2_ridge) + 
  geom_point(aes(x = truth,  
                 y = estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 24))


## ---- include = T, echo = F, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex2_lasso_cv = cv.glmnet(x, y, nfolds = 5, alpha = 1);
plot(ex2_lasso_cv);


## ---- include = T, echo = F, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex2_lasso <- 
  glmnet(x, y, lambda = ex2_lasso_cv$lambda.1se, alpha = 1) %>%
  tidy(return_zeros = TRUE) %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "lasso", 
            partition = 1)

ggplot(ex2_lasso) + 
  geom_point(aes(x = truth,  
                 y = estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 24))


## ---- include = T, echo = F, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex2_enet_cv = cv.glmnet(x, y, nfolds = 5, alpha = 0.5);
plot(ex2_enet_cv);


## ---- include = T, echo = F, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
ex2_enet <- 
  glmnet(x, y, lambda = ex2_enet_cv$lambda.1se, alpha = 0.1) %>% 
  tidy(return_zeros = TRUE) %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "enet", 
            partition = 1)

ggplot(ex2_enet) + 
  geom_point(aes(x = truth,  
                 y = estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 24))


## ---- include = T, echo = F, cache = T,  fig.width = fig.x, fig.height = fig.y, fig.align = 'center'----
# not shown 
ex2_ridge_cv_p2 = cv.glmnet(x, y, nfolds = 5, alpha = 0);

ex2_ridge_p2 <- 
  glmnet(x, y, lambda = ex2_ridge_cv_p2$lambda.1se, alpha = 0) %>% 
  tidy() %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "ridge", 
            partition = 2)

ex2_lasso_cv_p2 = cv.glmnet(x, y, nfolds = 5, alpha = 1);

ex2_lasso_p2 <- 
  glmnet(x, y, lambda = ex2_lasso_cv_p2$lambda.1se, alpha = 1) %>%
  tidy(return_zeros = TRUE) %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "lasso", 
            partition = 2)


ex2_enet_cv_p2 = cv.glmnet(x, y, nfolds = 5, alpha = 0.5);

ex2_enet_p2 <- 
  glmnet(x, y, lambda = ex2_enet_cv_p2$lambda.1se, alpha = 0.1) %>% 
  tidy(return_zeros = TRUE) %>%
  filter(term != "(Intercept)") %>%
  bind_cols(truth = true_beta, 
            method = "enet", 
            partition = 2)



## ---- echo = F----------------------------------------------------------------------
ex2_all_results <- 
  bind_rows(
    ex2_ridge,
    ex2_ridge_p2,
    ex2_lasso,
    ex2_lasso_p2,
    ex2_enet,
    ex2_enet_p2) %>%
  mutate(method = factor(method, levels = c("ridge", "lasso", "enet")), 
         term = fct_inorder(term) %>% as.numeric(), 
         partition = factor(partition))

ex2_all_results %>%
  group_by(method, partition) %>%
  summarize(rmse = sqrt(mean((estimate - truth)^2))) %>%
  knitr::kable(digits = 6)

ggplot(ex2_all_results) + 
  geom_point(aes(x = term, 
                 y = estimate, 
                 color = partition), 
             size = 0.5) + 
  facet_grid(~method) + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none", x = "none") + 
  theme(text = element_text(size = 24))


## ---- echo = F----------------------------------------------------------------------
ex2_lasso_nonzero <- 
  ex2_lasso %>%
  mutate(term = fct_inorder(term) %>% as.numeric()) %>%
  filter(estimate != 0)
x_reduced <- 
  x[,ex2_lasso_nonzero$term]
colnames(x_reduced) = ex2_lasso_nonzero$term

ex2_relaxed_lasso <- 
  lm(y ~ x_reduced) %>%
  tidy() %>%
  filter(term!= "(Intercept)") %>%
  mutate(term = str_replace(term, "x_reduced","") %>% as.numeric())

ggplot() + 
  geom_point(aes(x = ex2_lasso_nonzero$estimate,
                 y = ex2_relaxed_lasso$estimate)) + 
  geom_abline(intercept = 0, slope = 1) + 
  xlab("Original Lasso") + 
  ylab("MLE using\nLasso-selected variables") +
  theme(text = element_text(size = 24))


