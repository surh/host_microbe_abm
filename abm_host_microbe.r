library(tidyverse)
library(gganimate)


world_size <- 25
nutA_0 <- 10
nutB_0 <- 10
nutA_prod <- nutB_prod <- 0.1
nutA_K <- nutB_K <- 10
host_occupancy <- 0.3
prop_microbe <- 0.1
D_H <- D_Hm <- 0.01
R_H <- R_Hm <- 0.1
transmission_threshold <- 2
S_M <- 0.01
n_gens <- 100

set.seed(1291)
W_0 <- expand_grid(x = 1:world_size,
            y = 1:world_size,
            nut_A = nutA_0,
            nut_B = nutB_0) %>%
  mutate(occupant = sample(x = c("H", "Hm", "E"),
                           size = world_size ^ 2,
                           replace = TRUE,
                           prob = c(host_occupancy * (1 - prop_microbe),
                                    host_occupancy * prop_microbe,
                                    1 - host_occupancy)))
# W_0 %>%
#   ggplot(aes(x = x, y = y)) +
#   geom_tile(aes(fill = occupant)) +
#   theme_classic()

# W_prev <- W_0
# source_id <- 0

W_curr <- W_0
Res <- W_curr %>%
  mutate(gen = 0)
for(gen in 1:10){
  cat("====== Generation ", gen, " =====\n")
  # 1. Deaths
  ii_H <- W_curr$occupant == "H"
  ii_D <- ( runif(n = sum(ii_H), min = 0, max = 1) < D_H ) | ( W_curr$nut_A[ii_H] < 1 )
  W_curr$occupant[ which(ii_H)[which(ii_D)] ] <- "E"
  
  ii_H <- W_curr$occupant == "Hm"
  ii_D <- ( runif(n = sum(ii_H), min = 0, max = 1) < D_Hm ) | ( W_curr$nut_A[ii_H] + W_curr$nut_B[ii_H] < 1 )
  W_curr$occupant[ which(ii_H)[which(ii_D)] ] <- "E"
  
  # 2. Reproduce
  # NOTE: I should join the 2 reproduce section in just one for loop
  ii_H <- W_curr$occupant == "H" & W_curr$nut_A > 0
  ii_R <- runif(n = sum(ii_H), min = 0, max = 1) < R_H
  
  To_reproduce <- W_curr[ which(ii_H)[which(ii_R)], ]
  # Maybe I can randomize the order of reproductions to ensure
  # no bias in the position
  for(i in 1:nrow(To_reproduce)){
    neighbors <- expand_grid(xn = (To_reproduce$x[i]-1):(To_reproduce$x[i]+1),
                             yn = (To_reproduce$y[i]-1):(To_reproduce$y[i]+1)) %>%
      filter(xn > 0) %>%
      filter(yn > 0) %>%
      filter(xn <= world_size) %>%
      filter(yn <= world_size) %>%
      filter(!(xn == To_reproduce$x[i] & yn == To_reproduce$y[i])) %>%
      left_join(W_curr, by = join_by(xn == x, yn == y))
    
    ii_E <- neighbors$occupant == "E"
    if(sum(ii_E) == 0){
      next
    }else{
      ii_new <- sample(x = which(ii_E), size = 1)
      
      W_curr$occupant[ W_curr$x == neighbors$xn[ii_new] & W_curr$y == neighbors$yn[ii_new] ] <- "H"
      ii_nut <- W_curr$x == To_reproduce$x[i] & W_curr$y == To_reproduce$y[i]
      W_curr$nut_A[ ii_nut ] <- W_curr$nut_A[ ii_nut ] - 1
      
    }
  }
  
  ii_H <- W_curr$occupant == "Hm" & (W_curr$nut_A > 0 | W_curr$nut_B > 0) 
  ii_R <- runif(n = sum(ii_H), min = 0, max = 1) < R_Hm
  
  To_reproduce <- W_curr[ which(ii_H)[which(ii_R)], ]
  # Maybe I can randomize the order of reproductions to ensure
  # no bias in the position
  for(i in 1:nrow(To_reproduce)){
    neighbors <- expand_grid(xn = (To_reproduce$x[i]-1):(To_reproduce$x[i]+1),
                             yn = (To_reproduce$y[i]-1):(To_reproduce$y[i]+1)) %>%
      filter(xn > 0) %>%
      filter(yn > 0) %>%
      filter(xn <= world_size) %>%
      filter(yn <= world_size) %>%
      filter(!(xn == To_reproduce$x[i] & yn == To_reproduce$y[i])) %>%
      left_join(W_curr, by = join_by(xn == x, yn == y))
    
    ii_E <- neighbors$occupant == "E"
    if(sum(ii_E) == 0){
      next
    }else{
      ii_new <- sample(x = which(ii_E), size = 1)
      
      W_curr$occupant[ W_curr$x == neighbors$xn[ii_new] & W_curr$y == neighbors$yn[ii_new] ] <- "Hm"
      ii_nut <- W_curr$x == To_reproduce$x[i] & W_curr$y == To_reproduce$y[i]
      nut_A_curr <- W_curr$nut_A[ ii_nut ]
      nut_B_curr <- W_curr$nut_B[ ii_nut ]
      nut_consume <- sample(c("A", "B"), size = 1, prob = c(nut_A_curr, nut_B_curr))
      if(nut_consume == "A"){
        W_curr$nut_A[ ii_nut ] <- nut_A_curr - 1
      }else if(nut_consume == "B"){
        W_curr$nut_B[ ii_nut ] <- nut_B_curr - 1
      }else{
        stop("ERROR: nut_consume", call. = TRUE)
      }
    }
  }
  
  
  # 3. Transmit microbe
  ii_H <- W_curr$occupant == "H"
  for(i in which(ii_H)){
    i <- which(ii_H)[1]
    neighbors <- expand_grid(xn = (W_curr$x[i]-1):(W_curr$x[i]+1),
                             yn = (W_curr$y[i]-1):(W_curr$y[i]+1)) %>%
      filter(xn > 0) %>%
      filter(yn > 0) %>%
      filter(xn <= world_size) %>%
      filter(yn <= world_size) %>%
      filter(!(xn == W_curr$x[i] & yn == W_curr$y[i])) %>%
      left_join(W_curr, by = join_by(xn == x, yn == y))
    
    
    if(sum(neighbors$occupant == "Hm") >= transmission_threshold){
      W_curr$occupant[i] <- "Hm" 
    }
  }
  
  # 4. Shed microbe
  ii_H <- W_curr$occupant == "Hm"
  ii_S <- runif(n = sum(ii_H), min = 0, max = 1) < S_M
  W_curr$occupant[ which(ii_H)[which(ii_S)] ] <- "H"
  
  # 5. Replenish nutrient
  ii_nut <- runif(n = world_size ^ 2, min = 0, max = 1) < nutA_prod & W_curr$nut_A < nutA_K
  W_curr$nut_A[ ii_nut ] <- W_curr$nut_A[ ii_nut ] + 1
  
  ii_nut <- runif(n = world_size ^ 2, min = 0, max = 1) < nutB_prod & W_curr$nut_B < nutB_K
  W_curr$nut_B[ ii_nut ] <- W_curr$nut_B[ ii_nut ] + 1
  
  Res <- Res %>%
    bind_rows(W_curr %>%
                mutate(gen = gen))
}

write_tsv(Res, "sim.tsv")


a1 <- Res %>%
  ggplot(aes(x = x, y = y)) +
  geom_tile(aes(fill = occupant)) +
  theme_classic() +
  transition_manual(gen) +
  labs(title = "Generation: {current_frame}")
anim_save("sim.gif", a1, fps = 5, nframes = max(Res$gen))


# Res %>%
#   filter(gen == 0) %>%
#   ggplot(aes(x = x, y = y)) +
#   geom_tile(aes(fill = occupant)) +
#   theme_classic()
# 
# Res %>%
#   filter(gen == 10) %>%
#   ggplot(aes(x = x, y = y)) +
#   geom_tile(aes(fill = occupant)) +
#   theme_classic()
# 
# 
# 
# Res %>%
#   filter(gen == 10) %>%
#   ggplot(aes(x = x, y = y)) +
#   geom_tile(aes(fill = nut_A)) +
#   theme_classic()
# 
# Res %>%
#   filter(gen == 10) %>%
#   ggplot(aes(x = x, y = y)) +
#   geom_tile(aes(fill = nut_B)) +
#   theme_classic()

# 
# 
# W_temp <- W_curr %>%
#   pmap_dfr(function(x, y, nut_A, nut_B, occupant, World, world_size = 100){
#     neighbors <- expand_grid(xn = (x-1):(x+1),
#                              yn = (y-1):(y+1)) %>%
#       filter(xn > 0) %>%
#       filter(yn > 0) %>%
#       filter(xn <= world_size) %>%
#       filter(yn <= world_size) %>%
#       filter(!(xn == x & yn == y)) %>%
#       left_join(World, by = join_by(xn == x, yn == y)) %>%
#       mutate(x = x, y = y)
#     
#     to_grow <- NULL
#     
#     if(occupant == "H"){
#       if(nut_A == 0){
#         # Die if no nutrient
#         occupant <- "E"
#       }else if(nut_A > 0){
#         # If nutrient
#         prob <- runif(n = 1, min = 0, max = 1)
#         if(prob < D_H){
#           # Die with D_H prob
#           occupant <- "E"
#         }else if(prob < D_H + R_H){
#           # Use nutrient to grow if empty space
#           ii_E <- which(neighbors$occupant == "E")
#           to_grow <- neighbors[sample(ii_E, size = 1),]
#         }
#       }else{
#         stop("ERROR: nutA", call. = TRUE)
#       }
#     }else if(occupant == "Hm"){
#       if(nut_A == 0 && nut_B == 0){
#         # Die if no nutrient
#         occupant <- "E"
#       }else{
#         # If nutrient
#         prob <- runif(n = 1, min = 0, max = 1)
#         if(prob < D_Hm){
#           # Die with prob D_Hm
#           occupant <- "E"
#         }else if(prob < D_Hm + R_Hm){
#           # Use nutrient to grow if empty space
#           ii_E <- which(neighbors$occupant == "E")
#           to_grow <- neighbors[sample(ii_E, size = 1),]
#         }
#       }
#     }
#     
#     tibble(x = x,
#            y = y,
#            nut_A = nut_A,
#            nut_B = nut_B,
#            occupant = occupant,
#            xn = NA,
#            yn = NA) %>%
#       bind_rows(to_grow)
# 
#   },World = W_curr, world_size = world_size)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 


# for(i in 1:nrow(W_curr)){
#   prob_i <- runif(n = 1, min = 0, max = 1)
#   # Die if no nutrient or with prob D_H
#   if(W_curr$occupant[i] == "H"){
#     if(prob_i < D_H || W_curr$nut_A[i] == 0){
#       W_curr$occupant[i] <- "E"
#     }else{
#       # If not died, we check neighbors
#       neighbors <- expand_grid(xn = (W_prev$x[i]-1):(W_prev$x[i]+1),
#                                yn = (W_prev$y[i]-1):(W_prev$y[i]+1)) %>%
#         filter(xn > 0) %>%
#         filter(yn > 0) %>%
#         filter(xn <= world_size) %>%
#         filter(yn <= world_size) %>%
#         filter(!(xn == W_curr$x[i] & yn == W_curr$y[i])) %>%
#         left_join(W_prev, by = join_by(xn == x, yn == y)) %>%
#         mutate(xs = W_prev$x[i], ys = W_prev$y[i])
#       
#       if(prob_i < D_H + R_H && any(neighbors$occupant == "E")){
#         # If nutrient and and space, reproduce with prob_i using nut_A
#         n_ii <- sample(nrow(neighbors),size = 1, replace = FALSE, prob = 1*(neighbors$occupant == "E"))
#         new_ii <- W_curr$x == neighbors$xn[n_ii] & W_curr$y == neighbors$yn[n_ii]
#         W_curr$occupant[new_ii] <- "H"
#         W_curr$nut_A[i] <- W_curr$nut_A[i]
#       }
#     }
#   }
#   
# }

