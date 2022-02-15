require(glue);
require(tidyverse);
require(ggplot);
require(gganimate);
require(gifski);

length_lambda_seq = 25
n_circle_points = 200;

# These are the number of points (equidistant on the quantile unit of measurement)
# to calculate the likelihood surface
n_quants = length_lambda_seq^2
p_seq_x = plogis(seq(-3, 4, length = n_quants));
p_seq_y = plogis(seq(-8, 8, length = n_quants));

# These are the parameters to construct the fake likelihood surface
# It's a bivariate normal, so 2 mean parameters and 3 covariance parameters
x_means = c(0.5, 0.5, 0.5);
x_sds = c(0.5, 1, 0.55);
y_means = c(1.5, 1.5, 1.5);
y_sds = c(0.4, 0.15, 0.65);
yx_regs = c(-0.1, -0.8, -0.8);


for(j in 1:3) {
  
  dist_data <-
    expand_grid(
      x_coord = qnorm(p_seq_x, mean = x_means[j], sd = x_sds[j]), 
      y_coord = qnorm(p_seq_y, mean = yx_regs[j] * x_coord + y_means[j], sd = y_sds[j]))
  
  dist_data <- 
    bind_rows(
      dist_data,
      expand_grid(
        x_coord = c(-0.005, -0.0025, 0, 0.0025, 0.005), 
        y_coord = seq(min(dist_data$y_coord), max(dist_data$y_coord), length = 100))
    ) %>%
    arrange(x_coord, y_coord) %>%
    mutate(count = dnorm(x_coord, mean = x_means[j], sd = x_sds[j], log = T) + 
             dnorm(y_coord, mean = yx_regs[j] * x_coord + y_means[j], sd = y_sds[j], log = T),
           count = exp(1 + count - min(count)),
           count_color = grey(.8*(1-count^1.5/max(count)^1.5))
    ) 
  
  xrange = range(dist_data$x_coord) + c(-1,1) * 0.1 * diff(range(dist_data$x_coord))
  yrange = range(dist_data$y_coord) + c(-1,1) * 0.1 * diff(range(dist_data$y_coord))
  
  foo1 = 
    with(dist_data,
         abs(x_coord[which.max(count)]) + abs(y_coord[which.max(count)]));
  foo2 = 
    with(dist_data, 
         sqrt(x_coord[which.max(count)]^2 + y_coord[which.max(count)]^2));
  classo_seq = seq(foo1,(2/5) * foo1,length = length_lambda_seq); 
  cridge_seq = seq(foo2,(2/5) * foo2,length = length_lambda_seq); 
  
  
  polygon_data = NULL;
  for(i in 1:length_lambda_seq) {
    
    curr_lasso_results <- 
      dist_data %>%
      mutate(lasso_criterion = (abs(x_coord) + abs(y_coord) <= classo_seq[i]) * count) %>%
      slice_max(lasso_criterion, n = 1, with_ties = F)
    
    curr_ridge_results <- 
      dist_data %>%
      mutate(ridge_criterion = (sqrt((x_coord)^2 + (y_coord)^2) <= cridge_seq[i]) * count) %>%
      slice_max(ridge_criterion, n = 1, with_ties = F)
    
    
    polygon_data = 
      bind_rows(
        polygon_data,
        tibble(
          lambda_index = i,
          label = "Lasso",
          lambda_x = curr_lasso_results %>% pull(x_coord),
          lambda_y = curr_lasso_results %>% pull(y_coord),
          x = c(-classo_seq[i],0,classo_seq[i],0,-classo_seq[i]),
          y = c(0,classo_seq[i],0,-classo_seq[i], 0)),
        tibble(
          lambda_index = i,
          label = "Ridge",
          lambda_x = curr_ridge_results %>% pull(x_coord),
          lambda_y = curr_ridge_results %>% pull(y_coord),
          x = cridge_seq[i]*c(seq(-1, 1,length = n_circle_points),
                              seq(1, -1,length = n_circle_points)),
          y = sqrt(cridge_seq[i]^2-x^2)*c(rep(-1, n_circle_points),rep(1, n_circle_points))))
  }
  polygon_data = 
    polygon_data %>%
    mutate(lambda_index = factor(lambda_index, levels = seq_len(length_lambda_seq), ordered = T))
  rm(foo1, foo2, curr_lasso_results, curr_ridge_results, classo_seq, cridge_seq, i)
  
  
  p <- 
    ggplot(data = polygon_data) +
    geom_point(data = dist_data, 
               mapping = aes(x = x_coord, 
                             y = y_coord),
               color = dist_data$count_color,
               size = 0.175, 
               alpha = 1) + 
    geom_hline(yintercept = 0, linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_path(mapping = aes(x = x, y = y, group = label)) + 
    geom_point(mapping = aes(x = lambda_x, 
                             y = lambda_y, 
                             color = label,
                             shape = label),
               size = 3) + 
    scale_color_discrete(name = "") + 
    scale_shape_discrete(name = "") +
    scale_x_continuous(name = expression(beta[1])) + 
    scale_y_continuous(name = expression(beta[2])) +
    coord_fixed(
      xlim = xrange, 
      ylim = yrange, 
      expand = FALSE) + 
    theme(text = element_text(size = 14),
          legend.position = c(0.7, 0.75),
          legend.background = element_blank(),
          legend.box.background = element_blank()) +
    transition_states(lambda_index,
                      transition_length = 1,
                      state_length = 1,
                      wrap = FALSE) +
    shadow_trail(alpha = 0.5, size = size / 2, exclude_layer = 4)
    
  
  anim_save(glue("plot{j}.gif"), p, 
            nframes = 4 * length_lambda_seq, 
            duration = 12,  
            units = "in", 
            res = 96,
            height = 6, 
            width = 6,
            end_pause = 3,
            rewind = T)
  
}



if(0) {
  ggplot() +
    geom_point(data = dist_data, 
               mapping = aes(x = x_coord, 
                             y = y_coord),
               color = dist_data$count_color,
               size = 0.175)  + 
    scale_color_discrete(name = "") + 
    scale_shape_discrete(name = "") +
    scale_x_continuous(name = expression(beta[1])) + 
    scale_y_continuous(name = expression(beta[2])) +
    coord_fixed(ratio = 1, 
                xlim = c(-xmax, xmax), 
                ylim = c(-xmax, xmax), expand = 0) + 
    theme(text = element_text(size = 14),
          legend.position = c(0.1, 0.33),
          legend.background = element_blank(),
          legend.box.background = element_blank()) 
}
