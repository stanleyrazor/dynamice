
pacman::p_load(dplyr, ggplot2, haven, tidyr, purrr, readxl,
               lubridate, stringr, naniar, INLA)

# data files
load("~/Documents/GitHub/KDHS Model Updates/data/2022 DHS/KEBR8AFL_2022.RData")
birth22 <- birth; rm(birth)
birth14 <- read_dta("~/Documents/GitHub/KDHS Model Updates/data/2014 BIRTH DHS/KEBR72FL.DTA")

# loading national MCV1 & MCV2
g1 <- readxl::read_xlsx("own-reports/Measles vaccination coverage 2024-04-12 10-40 UTC-2.xlsx")
g2 <- g1 |>
  filter(COVERAGE_CATEGORY == "ADMIN" & ANTIGEN == "MCV1") |>
  select(year = YEAR, target = TARGET_NUMBER, gcov = COVERAGE) |>
  mutate(gcov = gcov / 100)

# loading county MCV1 & MCV2
c0 <- read.csv("own-reports/county-cases-measles.csv") |>
  filter(Data != "Measles") |>
  select(county = Organisation.unit, antigen = Data, year = Period, value = Value) |>
  mutate(value = value / 100,
         value = ifelse(value > 1, .99, value),

         antigen = ifelse(antigen == "Proportion of under 1 year receiving vaccine against Measles and Rubella 1",
                          "mcv1",
                          "mcv2"),

         county = str_replace_all(county, " County", ""),
         county = str_to_lower(county),
         county = case_when(
           county == "elgeyo marakwet" ~ "elgeyo-marakwet",
           county == "muranga" ~ "murang'a",
           county == "tharaka nithi" ~ "tharaka-nithi",
           T ~ as.character(county)
         )) |>
  pivot_wider(names_from = antigen, values_from = value)

c1 <- read.csv("own-reports/county-cases-measles.csv") |>
  filter(Data == "Proportion of under 1 year receiving vaccine against Measles and Rubella 1") |>
  select(county = Organisation.unit, year = Period, value = Value) |>
  mutate(value = value / 100,
         value = ifelse(value > 1, .99, value),

         county = str_replace_all(county, " County", ""),
         county = str_to_lower(county),
         county = case_when(
           county == "elgeyo marakwet" ~ "elgeyo-marakwet",
           county == "muranga" ~ "murang'a",
           county == "tharaka nithi" ~ "tharaka-nithi",
           T ~ as.character(county)
         ))

c2 <- c1 |>
  pivot_wider(names_from = year, values_from = value)

cor(c2[, -1], method = "pearson")


# -------------------------------------------------------------------------

b14 <- birth14 |>
  select(caseid, v005, county = sregion, midx, mcv1 = h9,
         date = h9y, yob = b2) |>
  as_factor() |>
  filter(!is.na(midx)) |>
  filter(mcv1 != "don't know") |>
  mutate(v005 = v005 / 1e6,
         date = as.character(date),
         date = ifelse(!is.na(mcv1) & is.na(date), yob + 1, date),

         indicator = data.table::fifelse(test = mcv1 %in% c("vaccination date on card",
                                        "reported by mother",
                                        "vaccination marked on card"),
                            yes = 1,
                            no = 0,
                            na = NA)) |>
  select(-yob) |>
  filter(date %in% 2010:2014) |>

  mutate(
    county = case_when(
      county == "muranga" ~ "murang'a",
      county == "tharaka" ~ "tharaka-nithi",
      county == "trans-nzoia" ~ "trans nzoia",
      county == "elgeyo marak" ~ "elgeyo-marakwet",
      TRUE ~ as.character(county)
    )

  )
b14


# -------------------------------------------------------------------------

b22 <- birth22 |>
  select(caseid, v005, county = v024, midx, mcv1 = h9,
         date = h9y, yob = b2) |>
  as_factor() |>
  filter(!is.na(midx)) |>
  filter(mcv1 != "don't know") |>
  mutate(v005 = v005 / 1e6,
         date = as.character(date),
         date = ifelse(!is.na(mcv1) & is.na(date), yob, date),

         indicator = data.table::fifelse(test = mcv1 %in% c("vaccination date on card",
                                                            "reported by mother",
                                                            "vaccination marked on card"),
                                         yes = 1,
                                         no = 0,
                                         na = NA)) |>
  select(-yob) |>
  filter(date %in% 2018:2022)
b22


# -------------------------------------------------------------------------

d1 <- bind_rows(b14 |>
            select(v005, county, date, indicator) |>
            mutate(county = as.character(county),
                   year = 2013),
          b22 |>
            select(v005, county, date, indicator) |>
            mutate(county = as.character(county),
                   year = 2020)) |>
  group_by(county, year) |>
  reframe(mcv1 = weighted.mean(indicator, v005)) |>
  mutate(year = as.numeric(year))

d2 <- bind_rows(d1,
                c1 |> select(county, year, mcv1 = value),
                expand.grid(county = unique(d1 |> pull(county)),
                            year = c(2014,2015:2019),
                            mcv1 = NA)) |>
  mutate(county = factor(county),
         countynum = as.numeric(county))

ggplot(d2) +
  geom_point(aes(x = year, y = mcv1)) +
  facet_wrap(~county) +
  theme_bw()


# -------------------------------------------------------------------------

d3 <- d2 |>
  mutate(l_mcv1 = qlogis(mcv1),
         g_year = year) |>
  merge(g2 |> select(year, gcov),
        by = "year", all = T)


# National-county-coverage ------------------------------------------------

dhs <- bind_rows(b14 |> select(v005, indicator) |>
                   mutate(year = 2013),
                 b22 |> select(v005, indicator) |>
                   mutate(year = 2020)) |>
  group_by(year) |>
  reframe(gcov = weighted.mean(indicator, v005))


f1 <- merge(bind_rows(c1, d1 |> select(county, year, value = mcv1)),
            g2 |> select(-target), by = "year", all.x = T)
f2 <- f1 |>
  mutate(multiplier = value / gcov)

f2 |>
  group_by(county) |>
  mutate(med = median(multiplier)) |>
  ungroup() |>

  ggplot(aes(x = reorder(county, med), y = multiplier)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 1), col = "red") +
  theme(axis.text.x = element_blank())

mult_vals <- f2 |>
  group_by(county) |>
  reframe(mult = median(multiplier)) |>
  ungroup()

f3 <- expand.grid(
  county = mult_vals |> pull(county),
  year = g2 |> pull(year)
) |>
  merge(g2, by = "year") |>
  merge(mult_vals, by = "county") |>
  mutate(mcv1 = mult * gcov)

f3 |>
  ggplot() +
  geom_line(aes(x = year, y = mcv1)) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~county) +
  theme_classic()


# Population backprojection -----------------------------------------------

p1 <- rKenyaCensus::V3_T2.3 |>
  filter(Age == 0 & SubCounty == "ALL") |>
  select(county = County, pop = Total) |>
  mutate(prop = pop / sum(pop),
         county = str_to_lower(county),
         county = ifelse(county == "taita/ taveta", "taita taveta", county)) |>
  select(county, prop)

load("data/data_pop.rda")
data_pop <- data_pop |>
  filter(country == "Kenya" & age_from == 0 & year %in% 2013:2030) |>
  select(year, n_pop = value) |>
  data.frame()

p2 <- expand.grid(county = p1 |> pull(county),
            year = data_pop |> pull(year)) |>
  merge(p1, by = "county") |>
  merge(data_pop, by = "year") |>
  mutate(c_pop = n_pop * prop)

p3 <- p2 |>
  select(county, year, population = c_pop)

## merging the mcv1 coverage and the population datasets
f4 <- merge(p3, f3, by = c("county", "year")) |>
  group_by(year) |>
  mutate(target = sum(population),
         national_vax = gcov * target) |>
  ungroup()


### trying out a single year
# tval <- f4 |> filter(year == 2013) |> pull(national_vax) |> unique()
# popvec <- f4 |> filter(year == 2013) |> pull(population)
# init_par <- f4 |> filter(year == 2013) |> pull(mcv1)
#
# fn <- \(par, tval, popvec) {
#   cval <- sum(popvec * par)
#   abs(tval - cval)
# }
#
# (opar <- optim(par = init_par,
#                fn,
#                tval = tval, popvec = popvec,
#                method = "L-BFGS-B",
#                lower = c(0, 0), upper = c(1, 1),
#                control = list(maxit = 100)))

# function for benchmarking in the dasymetric redistribution approach
optim_bench <- \(mcv, pop, target) {

  fn <- \(par, tval, popvec) {
    cval <- sum(popvec * par)
    abs(tval - cval)
  }

  opar <- optim(par = mcv,
                fn,

                tval = first(target),
                popvec = pop,

                method = "L-BFGS-B",
                lower = c(0, 0), upper = c(1, 1),
                control = list(maxit = 100))

  opar$par
}

f5 <- f4 |>
  group_by(year) |>
  mutate(mcv1_dasred = optim_bench(mcv1 = mcv1, pop = population, target = national_vax)) |>
  ungroup()

ggplot(f5) +
  geom_point(aes(x = mcv1, y = mcv1_dasred),
             shape = '*', size = 1.6) +
  geom_abline(aes(intercept = 0, slope = 1), col = "red") +
  facet_wrap(~year) +
  theme_bw(base_line_size = 0)

f5 |>
  ggplot() +
  geom_line(aes(x = year, y = mcv1), col = "red") +
  geom_line(aes(x = year, y = mcv1_dasred)) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~county) +
  theme_bw(base_line_size = 0)

f6 <- f5 |>
  select(county, year, population, mcv1 = mcv1_dasred)

# MCV2 at the county level ------------------------------------------------

c2 <- g1 |>
  filter(COVERAGE_CATEGORY == "ADMIN") |>
  select(antigen = ANTIGEN, year = YEAR, coverage = COVERAGE) |>
  mutate(antigen = str_to_lower(antigen),
         coverage = coverage / 1e2) |>
  pivot_wider(names_from = antigen, values_from = coverage)

c2_train <- c2 |> na.omit()
c2_test <- c2 |> filter(is.na(mcv2))

with(c2, plot(mcv1, mcv2, pch = '*'))

c2_m1 <- mgcv::gam(mcv2 ~ s(year, k = 8, bs = "ts") + mcv1,
                   data = c2_train,
                   family = binomial)
c2 |>
  mutate(pred = predict(c2_m1, c2, type = "response")) |>
  ggplot() +
  geom_point(aes(x = mcv1, y = mcv2)) +
  geom_line(aes(x = mcv1, y = pred))

# now using the national MCV1-MCV2 relationship to predict county MCV2
f7 <- f6 |>
  merge(c0 |>
          mutate(mult = mcv2 / mcv1) |>
          group_by(county) |>
          reframe(mult = mean(mult)),
        by = "county") |>
  mutate(mcv2 = .5 * predict(c2_m1, f6, type = "response") |> as.numeric() +
           .5 * (mcv1 * mult))

# benchmarking the MCV2 to match nationals
bench1 <- g1 |>
  filter(ANTIGEN == "MCV2" & COVERAGE_CATEGORY == "ADMIN") |>
  select(year = YEAR, mcv2 = COVERAGE) |>
  mutate(mcv2 = mcv2 / 100,
         mcv2 = ifelse(year < 2015, 0, mcv2))

# imputing the 2020 value, with regression model of 2015-2023 data
# with(bench1, plot(year, mcv2, pch = '*'))
# lm(mcv2~year, data=bench1 |> filter(year >= 2015))
bench1 <- bench1 |>
  mutate(mcv2 = ifelse(year == 2020, .51789, mcv2))

load("data/data_pop.rda")
data_pop <- data_pop |>
  filter(country == "Kenya" & age_from == 1 & year %in% 2013:2030) |>
  select(year, n_pop = value) |>
  data.frame()

p2 <- expand.grid(county = p1 |> pull(county),
                  year = data_pop |> pull(year)) |>
  merge(p1, by = "county") |>
  merge(data_pop, by = "year") |>
  mutate(c_pop = n_pop * prop)

p3 <- p2 |>
  select(county, year, population = c_pop) |>
  merge(bench1 |> select(year, n_mcv2 = mcv2), by = "year") |>
  merge(data_pop, by = "year") |>
  mutate(national_vax = n_mcv2 * n_pop)

f8 <- p3 |>
  merge(f7 |> select(county, year, mcv2),
        by = c("county", "year")) |>
  group_by(year) |>
  mutate(mcv2_dasred = optim_bench(mcv = mcv2, pop = population, target = national_vax)) |>
  ungroup()

### aggregated results
f9 <- f7 |> select(county, year, mcv1) |>
  merge(f8 |> select(county, year, mcv2 = mcv2_dasred),
        by = c("county", "year"))

# -------------------------------------------------------------------------

# very important plot:
a1 <- f9 |>
  filter(county %in% c("migori", "homa bay", "west pokot", "nairobi", "makueni", "vihiga", "samburu", "kilifi"))|>
  ggplot() +
  geom_line(aes(x = year, y = mcv1, col = "MCV1")) +
  geom_line(aes(x = year, y = mcv2, col = "MCV2")) +

  geom_point(data = d1 |>
               filter(county %in% c("migori", "homa bay", "west pokot", "nairobi", "makueni", "vihiga", "samburu", "kilifi")),
             aes(x = year, y = mcv1, shape = "DHS"), col = "black") +
  geom_point(data = c0 |>
               filter(county %in% c("migori", "homa bay", "west pokot", "nairobi", "makueni", "vihiga", "samburu", "kilifi")),
             aes(x = year, y = mcv1, shape = "KHIS"), col = "black") +
  geom_point(data = c0 |>
               filter(county %in% c("migori", "homa bay", "west pokot", "nairobi", "makueni", "vihiga", "samburu", "kilifi")),
             aes(x = year, y = mcv2, shape = "KHIS"), col = "blue") +

  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(from = 2013, to = 2023, by = 3)) +
  scale_color_manual(values = c("MCV1" = "black", "MCV2" = "blue")) +

  labs(title = "Estimates for historical MCV coverage",
       shape = "Data source:",
       col = "Estimate:",
       x = "Year",
       y = "Coverage (%)") +
  facet_wrap(~county) +
  theme_bw(base_line_size = 0,
           base_family = "Times New Roman")
a1

ggsave("img/dasred_county_plot.png", plot = a1,
width = 8, height = 8, units = "in", dpi = 300)


# creating the results data -----------------------------------------------

county_codes <- read.csv("own-reports/county_codes.csv")

# routine vaccination data
routine <- f9 |>
  select(county, year, mcv1, mcv2) |>
  bind_rows(
    expand.grid(county = f9 |> pull(county) |> unique(),
                year = 2024:2030) |>
      merge(f9 |> filter(year == 2023) |> select(county, mcv1, mcv2),
            by = "county") |>
      group_by(county) |>
      arrange(year) |>
      mutate(id = row_number()) |>
      ungroup() |>

      # projecting to grow by 1% yearly
      mutate(mcv1 = mcv1 + (id * .01),
             mcv2 = mcv2 + (id * .01)
             ) |>
      select(county, year, mcv1, mcv2)
  ) |>
  mutate(mcv1 = ifelse(mcv1 > .99, .99, mcv1),
         mcv2 = ifelse(mcv2 > .99, .99, mcv2))

routine_mcv1 <- routine |>
  select(county, year, coverage = mcv1) |>
  mutate(vaccine = "MCV1") |>
  merge(county_codes |> select(county = name, country_code = code),
        by = "county") |>
  select(vaccine, country_code, country = county, year, coverage)


routine_mcv2 <- routine |>
    merge(county_codes |> select(county = name, country_code = code),
          by = "county") |>
    select(country = county, country_code, year, MCV1 = mcv1, MCV2 = mcv2) |>
    pivot_longer(c(MCV1, MCV2)) |>
    select(vaccine = name, country_code, country, year, coverage = value)

routine_no_vacc <- routine_mcv2 |>
  mutate(coverage = 0)

write.csv(routine_mcv1, "vac_coverage_top10/base case/scenarios/county/temp_routine_mcv1.csv", row.names = F)
write.csv(routine_mcv2, "vac_coverage_top10/base case/scenarios/county/temp_routine_mcv2.csv", row.names = F)
write.csv(routine_no_vacc, "vac_coverage_top10/base case/scenarios/county/temp_routine_no_vacc.csv", row.names = F)
