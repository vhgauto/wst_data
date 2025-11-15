library(terra)
library(gt)
library(tidyterra)
library(tidyverse)

gain_b10 <- 0.0003342
bias_b10 <- .1

gain_b11 <- 0.0003342
bias_b11 <- .1

k1_b10 <- 774.8853
k2_b10 <- 1321.0789

k1_b11 <- 480.8883
k2_b11 <- 1201.1442

e_b10 <- 0.992
e_b11 <- 0.998
e_prom <- 0.995
e_delta <- -0.006

b0 <- 2.29250
b1 <- 0.99290
b2 <- 0.15450
b3 <- -0.31220
b4 <- 3.71860
b5 <- 0.35020
b6 <- -3.58890
b7 <- 0.18250

l <- list.files(
  path = "Matias_2025_RGB/",
  pattern = ".tif$",
  full.names = TRUE
)

l <- list.files(
  path = "Matias_2025_RGB/",
  pattern = ".tif$",
  full.names = TRUE
)

r <- map(l, rast)

f_mask <- function(x) {
  m <- ifel(x$QA_PIXEL == 21952, 1, NA)
  if (is.na(terra::global(m, na.rm = TRUE)$mean)) {
    r <- NULL
  } else {
    r <- x * m
  }
  return(r)
}

r_mask <- map(r, f_mask, .progress = TRUE)

r_agua <- tibble(
  x = r_mask
) |>
  mutate(raster = map_lgl(x, is.null)) |>
  filter(!raster) |>
  pull(x)

nubes <- tibble(
  x = r_mask
) |>
  mutate(nro = row_number()) |>
  mutate(raster = map_lgl(x, is.null)) |>
  filter(raster) |>
  pull(nro)

png(
  width = 30,
  height = 30 * 481 / 524,
  units = "cm",
  res = 300,
  bg = "transparent",
  filename = "fig/nubes.png"
)

par(mfcol = c(3, 3), mar = c(1, 1, 0, 0))

walk(
  nubes[1:9],
  ~ terra::plotRGB(
    r[[.x]],
    r = 3,
    g = 2,
    b = 1,
    smooth = FALSE,
    stretch = "lin",
    mar = c(1, 1, 1, 1)
  )
)

dev.off()

f_area <- function(x) {
  print(varnames(x))
  tibble(
    area = expanse(x$QA_PIXEL)$area,
    fecha = str_sub(varnames(x), 13, 20) |> ymd()
  )
}

area_embalse <- map(r_agua, f_area) |>
  list_rbind()

write_csv(wst_area, "datos/area_embalse.csv")

area_embalse <- read_csv("datos/area_embalse.csv", show_col_types = FALSE)
max_area <- max(area_embalse$area)

fechas_nubes_parciales <- area_embalse |>
  slice_min(order_by = area, n = 9) |>
  pull(fecha)

no_areas <- tibble(
  x = r
) |>
  mutate(nro = row_number()) |>
  mutate(
    fecha = map_chr(x, varnames)
  ) |>
  mutate(
    fecha = ymd(str_sub(fecha, -8))
  ) |>
  filter(fecha %in% fechas_nubes_parciales) |>
  pull(nro)

png(
  width = 30,
  height = 30 * 481 / 524,
  units = "cm",
  res = 300,
  bg = "transparent",
  filename = "fig/nubes_parcial.png"
)

par(mfcol = c(3, 3), mar = c(1, 1, 0, 0))

walk(
  no_areas[1:9],
  ~ terra::plotRGB(
    r[[.x]],
    r = 3,
    g = 2,
    b = 1,
    smooth = FALSE,
    stretch = "lin",
    mar = c(1, 1, 1, 1)
  )
)

dev.off()

area_embalse <- read_csv("datos/area_embalse.csv", show_col_types = FALSE)
max_area <- max(area_embalse$area)

areas <- area_embalse |>
  mutate(nro = row_number()) |>
  filter(area >= max_area * .7) |>
  pull(nro)

r_clean <- tibble(
  x = r_agua
) |>
  mutate(nro = row_number()) |>
  filter(nro %in% areas) |>
  pull(x)

f_radiancia_10 <- function(x) {
  x$B10 * gain_b10 + bias_b10
}

f_radiancia_11 <- function(x) {
  x$B11 * gain_b11 + bias_b11
}

r_l_10 <- map(r_clean, f_radiancia_10)
r_l_11 <- map(r_clean, f_radiancia_11)

f_t_10 <- function(x) {
  k2_b10 / (log(k1_b10 / x + 1))
}

f_t_11 <- function(x) {
  k2_b11 / (log(k1_b11 / x + 1))
}

r_t_10 <- map(r_l_10, f_t_10)
r_t_11 <- map(r_l_11, f_t_11)

f_wst <- function(t10, t11) {
  term1 <- b1 + b2 * ((1 - e_prom) / (e_prom)) + b3 * ((e_delta) / (e_prom^2))
  term2 <- b4 + b5 * ((1 - e_prom) / (e_prom)) + b6 * ((e_delta) / (e_prom^2))

  b0 + term1 * (t10 + t11) / 2 + term2 * (t10 - t11) / 2 + b7 * (t10 - t11)^2
}

wst <- map2(r_t_10, r_t_11, ~ f_wst(.x, .y))

f_promedio <- function(x) {
  tibble(
    temp = terra::global(x, fun = "mean", na.rm = TRUE)$mean,
    fecha = str_sub(varnames(x), 13, 20) |> ymd()
  )
}

wst_prom <- map(wst, f_promedio) |>
  list_rbind()

write_csv(wst_prom, "datos/wst_prom.csv")
wst_prom <- read_csv("datos/wst_prom.csv", show_col_types = FALSE)


g1 <- ggplot(wst_prom, aes(fecha, temp - 273.15)) +
  geom_line(linewidth = .3, color = "grey40") +
  geom_point(
    size = 1,
    shape = 21,
    color = "grey40",
    fill = "violetred",
    alpha = 1,
    stroke = .1
  ) +
  geom_smooth(aes(color = "a"), se = FALSE, formula = y ~ x, method = "loess") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(
    breaks = scales::breaks_width(2),
    limits = c(9, 31),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = "violetred",
    labels = "Tendencia",
    name = NULL
  ) +
  labs(x = NULL, y = "wst (°C)") +
  theme_classic() +
  theme(
    text = element_text(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_line(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    legend.background = element_rect(
      color = "black",
      linewidth = .3,
      fill = "#FAFAFA"
    ),
    legend.margin = margin(2, 2, 2, 2),
    legend.position = "inside",
    legend.key = element_blank(),
    legend.position.inside = c(1, 1),
    legend.justification.inside = c(1, 1)
  )

wst_año <- tibble(
  x = wst
) |>
  mutate(fecha = map_chr(x, ~ str_sub(varnames(.x), 13, 20))) |>
  mutate(fecha = ymd(fecha)) |>
  mutate(año = year(fecha)) |>
  slice_head(n = 1, by = año) |>
  pull(x)

f_mapa <- function(r) {
  titulo <- ymd(str_sub(varnames(r), 13, 20))

  ggplot() +
    geom_spatraster(
      data = r - 273.15
    ) +
    scale_fill_whitebox_c(
      palette = "muted",
      name = "wst (°C)",
      limits = c(16, 38),
      breaks = seq(16, 38, 2)
    ) +
    coord_sf(expand = FALSE) +
    labs(title = titulo) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = .5, face = "bold", size = 8),
    )
}

wst_mapa <- map(wst_año, f_mapa)

wst_mapa_comp <- patchwork::wrap_plots(wst_mapa, nrow = 3, guides = "collect") &
  theme(
    plot.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(margin = margin(b = 10, r = 5)),
    legend.key.height = unit(5, "pt"),
    legend.key.width = unit(45, "pt")
  )

f_fecha <- function(x) {
  varnames(x) |>
    str_sub(-8) |>
    ymd()
}

f_rgb <- function(x, carpeta) {
  tit <- f_fecha(r[[x]])

  png(
    width = 1000,
    height = round(1000 * 481 / 524),
    units = "px",
    bg = "transparent",
    filename = paste0("fig/", carpeta, "/", tit, ".png")
  )

  terra::plotRGB(
    r[[x]],
    r = 3,
    g = 2,
    b = 1,
    smooth = FALSE,
    stretch = "lin",
    mar = c(1, 1, 4, 1),
    main = as.character(tit),
    cex.main =
  )

  dev.off()
}

r_tbl <- tibble(
  x = r
) |>
  mutate(
    fecha = map(x, f_fecha) |> list_c()
  ) |>
  mutate(nro = row_number())

r_agua_tbl <- tibble(
  x_agua = r_agua
) |>
  mutate(
    fecha = map(x_agua, f_fecha) |> list_c()
  )

todas_nubes <- anti_join(r_tbl, r_agua_tbl, by = join_by(fecha)) |>
  pull(nro)

walk(todas_nubes, ~ f_rgb(.x, "nubes_png"), .progress = TRUE)

l_nubes <- list.files(
  path = "fig/nubes_png/",
  pattern = "png",
  full.names = TRUE
)

r_clean_tbl <- tibble(
  x_clean = r_clean
) |>
  mutate(
    fecha = map(x_clean, f_fecha) |> list_c()
  )

todas_nubes_parciales <- anti_join(r_tbl, r_clean_tbl, by = join_by(fecha)) |>
  inner_join(r_agua_tbl, by = join_by(fecha)) |>
  pull(nro)

walk(
  todas_nubes_parciales,
  ~ f_rgb(.x, "nubes_parciales_png"),
  .progress = TRUE
)

l_nubes_parciales <- list.files(
  path = "fig/nubes_parciales_png/",
  pattern = "png",
  full.names = TRUE
)

r_clean_tbl <- tibble(
  x_clean = r_clean
) |>
  mutate(
    fecha = map(x_clean, f_fecha) |> list_c()
  )

todas_despejado <- inner_join(r_tbl, r_clean_tbl, by = join_by(fecha)) |>
  pull(nro)

walk(todas_despejado, ~ f_rgb(.x, "despejado_png"), .progress = TRUE)

l_despejado <- list.files(
  path = "fig/despejado_png/",
  pattern = "png",
  full.names = TRUE
)
