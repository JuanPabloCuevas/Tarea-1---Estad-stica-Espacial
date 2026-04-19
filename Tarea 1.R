
library(lattice); library(grid)
set.seed(1012)

# Pregunta 1

grid <- expand.grid(
  x = seq(0, 1, length.out = 50),
  y = seq(0, 1, length.out = 50)
)

n <- nrow(grid)

# simulación
simular_tlc <- function(alpha, m, tau2){
  Z.sum <- rep(0, n)
  
  for (k in 1:m) {
    Wk <- c(rnorm(1, 0, sqrt(2 * alpha)),
            rnorm(1, 0, sqrt(2 * alpha)))
    Uk <- runif(1, 0, 2 * pi)
    
    arg <- grid$x * Wk[1] + grid$y * Wk[2] + Uk
    Zk <- sqrt(2) * cos(arg)
    
    Z.sum <- Z.sum + Zk
  }
  
  Z <- Z.sum / sqrt(m)
  
  if (tau2 > 0) {
    Z <- Z + rnorm(n, 0, sqrt(tau2))
  }
  
  return(Z)
}

# Ejemplo
Z1 <- simular_tlc(alpha = 5,  m = 100, tau2 = 0)
Z2 <- simular_tlc(alpha = 10, m = 100, tau2 = 0.1)
Z3 <- simular_tlc(alpha = 20, m = 100, tau2 = 0.5)

# levelplot(Z1 ~ grid$x * grid$y, main = "alpha = 5, tau2 = 0")
# levelplot(Z2 ~ grid$x * grid$y, main = "alpha = 10, tau2 = 0.1")
# levelplot(Z3 ~ grid$x * grid$y, main = "alpha = 20, tau2 = 0.5")

print(levelplot(Z1 ~ grid$x * grid$y, main="alpha = 5, tau2 = 0",
                aspect="iso"), split=c(1,1,2,2), more=TRUE)

print(levelplot(Z2 ~ grid$x * grid$y, main="alpha = 10, tau2 = 0.1",
                aspect="iso"), split=c(2,1,2,2), more=TRUE)

print(levelplot(Z3 ~ grid$x * grid$y, main="alpha = 20, tau2 = 0.5",
                aspect="iso"), split=c(1,2,2,2))

grid.text("Figura 1: Simulación de campos aleatorios variando alpha y efecto pepita",
          x = 0.5, y = 0.02, gp = gpar(fontsize = 10))

set.seed(1012)
# Tiempo
tiempo1 <- system.time({
  
  grid = expand.grid(1:50,1:50)
  n = nrow(grid)
  dista = as.matrix(dist(grid))
  covmat = exp(-3*dista/30)
  cc = chol(covmat)
  z = t(cc) %*% rnorm(n)
  
})

tiempo2 <- system.time({
  
  Z <- simular_tlc(alpha = 20, m = 100, tau2 = 0.5)
  
})

tabla <- data.frame(
  Metodo = c("Cholesky", "TLC"),
  User = c(tiempo1[1], tiempo2[1]),
  System = c(tiempo1[2], tiempo2[2]),
  Elapsed = c(tiempo1[3], tiempo2[3])
)

library(kableExtra)

kable(tabla, caption = "Tiempos de ejecución") %>%
  kable_styling(latex_options = "hold_position")

# Pregunta 2

library(geoR); library(ggplot2); library(dplyr); library(purrr); library(fields)

data <- read.csv("precipitaciones.txt", sep = "")

geo.datos <- as.geodata(data, coords.col = 1:2, data.col = 3)

par(mar = c(5, 4, 4, 2) + 0.1)

hist(geo.datos$data, probability = TRUE,
     main = "Histograma de precipitaciones")

lines(density(geo.datos$data, adjust = 2), col = "red", lwd = 2)

mtext("Figura 2: Histograma de precipitaciones", side = 1, line = 4)

z <- geo.datos$data
coords <- geo.datos$coords

breaks <- quantile(z, probs = seq(0, 1, length.out = 5), na.rm = TRUE)
colores <- c("lightblue", "deepskyblue", "dodgerblue3", "navy")
z.col <- colores[cut(z, breaks = breaks, include.lowest = TRUE)]

par(mar = c(5, 4, 4, 6))

plot(coords, col = z.col, pch = 19,
     main = "Mapa de precipitaciones")

image.plot(legend.only = TRUE, zlim = range(z), col = colores,
           smallplot = c(0.85, 0.88, 0.2, 0.8))

mtext("Figura 3: Mapa de precipitaciones", side = 1, line = 4)

compute_variograms_geor <- function(geo, widths, cutoff = 150) {
  setNames(widths, widths) |>
    map(\(w) {
      vg <- variog(geo,
                   max.dist = cutoff,
                   uvec = seq(0, cutoff, by = w))
      
      data.frame(
        dist = vg$u,
        gamma = vg$v,
        np = vg$n,
        lag_width = paste0("width = ", w)
      )
    })
}

plot_variogram_sensitivity_geor <- function(v_list) {
  bind_rows(v_list) |>
    ggplot(aes(x = dist, y = gamma, color = lag_width)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    theme_minimal() +
    labs(
      title = "Sensibilidad del variograma al ancho de lag",
      x = "Distancia",
      y = expression(hat(gamma)(h)),
      color = NULL
    )
}

par(mar = c(5, 4, 4, 2) + 0.1)

v_list <- compute_variograms_geor(geo.datos, widths = c(10, 30, 60), cutoff = 450)
plot_variogram_sensitivity_geor(v_list) +
  labs(caption = "Figura 4: Variograma experimental con cutoff = 150 y distintas widths")

vg <- variog(geo.datos,
             max.dist = 450,
             uvec = seq(0, 450, by = 30))

ajuste.esf <- variofit(vg,
                       cov.model = "spherical",
                       ini.cov.pars = c(0.003, 100),
                       nugget = 0)

# Ver ajustes sobre el variograma
plot(vg, main = "Modelo esférico")
lines(ajuste.esf, col = "blue")
legend("bottomright",
       legend = c("Esférico"),
       col = c("blue"), lty = 1, bty = "n")
mtext("Figura 5: Ajuste del modelo esférico en el variograma con width = 30", side = 1, line = 4)


library(gstat); library(sp)

data_res <- data.frame(
  x = geo.datos$coords[,1],
  y = geo.datos$coords[,2],
  z = geo.datos$data
)

coordinates(data_res) <- ~x+y

x_seq <- seq(min(geo.datos$coords[,1]),
             max(geo.datos$coords[,1]),
             length.out = 50)

y_seq <- seq(min(geo.datos$coords[,2]),
             max(geo.datos$coords[,2]),
             length.out = 50)

grid_df <- expand.grid(x = x_seq, y = y_seq)
coordinates(grid_df) <- ~x+y
gridded(grid_df) <- TRUE

vgm_model <- vgm(
  psill = ajuste.esf$cov.pars[1],   # varianza parcial
  model = "Sph",                    # esférico
  range = ajuste.esf$cov.pars[2],   # alcance
  nugget = ajuste.esf$nugget
)

krige_ok <- krige(
  z ~ 1,
  locations = data_res,
  newdata = grid_df,
  model = vgm_model
)

krige_df <- as.data.frame(krige_ok)

ggplot(krige_df, aes(x = x, y = y, fill = var1.pred)) +
  geom_raster() +
  geom_point(
    data = as.data.frame(data_res),
    aes(x = x, y = y),
    size = 0.8, color = "white", shape = 21, stroke = 0.2,
    inherit.aes = FALSE
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Mapa de interpolaciones por krigging",
    x = "coordx", y = "coordy",
    fill = "Predicción",
    caption = "Figura 6: Predicciones realizadas con Kriggin ordinario"
  )


ggplot(krige_df, aes(x = x, y = y, fill = var1.var)) +
  geom_raster() +
  geom_point(
    data = as.data.frame(data_res),
    aes(x = x, y = y),
    size = 0.8, color = "white", shape = 21, stroke = 0.2,
    inherit.aes = FALSE
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Varianza de predicción",
    x = "coordx", y = "coordy",
    fill = "Varianza",
    caption = "Figura 7: Varianza utilizando Kriggin ordinario"
  )

