
# pkgs
devtools::load_all()
pacman::p_load(dplyr, stringr, tidyr)

# assume perfect timeliness: all receive MCV1 at 9-month
adj.timeliness = TRUE  # FALSE
if (adj.timeliness){
  data_timeliness [!is.na(age) & age < 52*(9/12), timeliness := 0]
  data_timeliness [!is.na(age) & age >= 52*(9/12), timeliness := 1]
}
# assume a fixed R0
adj.fixR0 = 14.25  # based on Kenya: 14.25
if (!is.na(adj.fixR0)){
  data_r0 [ , r0 := adj.fixR0]
}


# -------------------------------------------------------------------------

fortran_names <- c("vaccine2019_sia_adjBeta_mcv2_rSIAdpd.exe",          # most updated version
                   "vaccine2019_sia_adjBeta_mcv2_nodpdSIA.exe",
                   "vaccine2019_sia_adjBeta_mcv2_stepVE_rSIAdpd.exe",
                   "vaccine2019_sia_adjBeta_mcv2_nodpdSIA_stepVE.exe",
                   "vaccine2019_sia_adjBeta_mcv2_stepVE_rSIAdpdSA.exe", # sensitivity analysis for SIA delivery
                   "vaccine2019_sia_adjBeta_mcv2_stepVE_rSIAnever.exe"  # alternative assumption for SIA delivery
)

contact_mats <- c("syn", "polymod", "prpmix", "unimix")

main <- list (
  scenarios    = c("Base",          # (1) Base: step change VE + random SIA delivery to zero-dose children
                   "Efficacy",      # (2) Base + linear VE
                   "SIAdeliver",    # (3) Base + non-random SIA delivery to zero-dose children
                   "ContactSyn",    # (4) Base + synthetic contact patterns
                   "ContactPrp",    # (5) Base + proportional mixing
                   "ContactUni",    # (6) Base + uniform mixing
                   "Update",        # (7) Update: linear VE + synthetic matrix + non-random SIA delivery
                   "SIAdeliverSA",  # (8) Base + sensitivity analysis for non-random SIA delivery
                   "SIAnever"       # (9) Base + SIA delivery assuming 7.7% never reached
  ),
  fortran_model = fortran_names [c(4, 2, 3, 4, 4, 4, 1, 5, 6)],
  contact_mat   = contact_mats  [c(2, 2, 2, 1, 3, 4, 1, 2, 2)],
  step_ve       =                c(T, F, T, T, T, T, F, T, T)
)


# -------------------------------------------------------------------------

var <- list (
  # vaccine coverage
  vaccine_coverage_folder           = "vac_coverage_top10/base case/",
  coverage_prefix                   = "coverage",
  touchstone                        = "_",
  antigen                           = NULL,
  vaccine_coverage_subfolder        = "scenarios/",

  # disease burden
  central_burden_estimate_folder    = "central_burden_estimate/",
  stochastic_burden_estimate_folder = "stochastic_burden_estimate/",

  # diagnostic plots folder
  plot_folder                       = "plots/",

  # modelling group name
  group_name                        = NULL,

  # log file name
  log_name                          = "test_log",

  # countries - specify iso3 codes to analyse only these countries
  #             or set it to "all" to analyse all included countries
  countries                         = c("KEN"),

  cluster_cores                     = 1,    # number of cores
  psa                               = 0     # psa runs; 0 for single central run, was previously 5
)

# for central run: set number of runs to 0 (var$psa = 0)
# for stochastic runs: set number of runs to greater than 0 (eg: var$psa = 200)
# burden_estimate_folder (central or stochastic)
if (var$psa == 0) {
  burden_estimate_folder <- var$central_burden_estimate_folder
} else {
  burden_estimate_folder <- var$stochastic_burden_estimate_folder
}

# folders for burden estimate results
dir.create (file.path (paste0 (getwd(), "/", burden_estimate_folder, "Wolfson/")), recursive = T)
dir.create (file.path (paste0 (getwd(), "/", burden_estimate_folder, "Portnoy/")), recursive = T)


vaccine_strategies <- c("no-vacc",    # (1) no vaccination (set vaccination and using_sia to 0)
                        "mcv1",       # (2) MCV1 only
                        "mcv2",       # (3) MCV1 & MCV2
                        "sia-conti",  # (4) MCV1 & MCV2 and SIAs
                        "sia-stop"    # (5) MCV1 & MCV2 and SIAs (stop SIAs after 2019)
)

# corresponding vaccination strategies used for Portnoy's CFR
portnoy_scenarios <- c("no-vaccination",    # (1) no vaccination (set vaccination and using_sia to 0)
                       "mcv1-default",      # (2) MCV1 only
                       "mcv2-default",      # (3) MCV1 & MCV2
                       "campaign-default",  # (4) MCV1 & MCV2 and SIAs
                       "campaign-default"   # (5) MCV1 & MCV2 and SIAs (stop SIAs after 2019)
)

# set SIAs and vaccination parameters for each scenario to minimize errors for running
set_sia         <- c (0, 0, 0, 1, 1)
set_vaccination <- c (0, 1, 2, 2, 2)


# -------------------------------------------------------------------------

update_coverage_inputs = FALSE
if (update_coverage_inputs) {
  # (1) create vaccine coverage file (0% coverage) for no vaccination scenario
  create_no_vaccination_coverage_file (
    no_vaccination_coverage_file = paste0 (var$vaccine_coverage_folder,
                                           "coverage_", scenarios[1],".csv"),
    vaccination_coverage_file    = paste0 (var$vaccine_coverage_folder,
                                           "coverage_", scenarios[2],".csv")
  )

  # (2) create a vaccine coverage file for no future SIA scenario
  continuous_sia_coverage <- fread (paste0(var$vaccine_coverage_folder,
                                           "coverage_", scenarios[4],".csv"))
  stop_sia_coverage <- continuous_sia_coverage [!(year > 2019 & activity_type == "campaign")]
  fwrite (stop_sia_coverage, file = paste0(var$vaccine_coverage_folder,
                                           "coverage_", scenarios[5],".csv"))

  # (3) generate 2 vaccine coverage files per scenario for routine and SIA
  for (index in 1:5) {
    create_vaccine_coverage_routine_sia (
      vaccine_coverage_folder    = var$vaccine_coverage_folder,
      vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
      coverage_prefix            = var$coverage_prefix,
      touchstone                 = var$touchstone,
      antigen                    = var$antigen,
      scenario_name              = scenarios [index],
      rev_cov                    = FALSE # coverage data have been pre-processed
    )
  }
}

imain <- 7

# -------------------------------------------------------------------------

# grid data:
grd <- expand.grid(
  scenario_index = 1:3,
  contact_mat = c("polymod", "syn", "unimix"), # prpmix
  step_ve = c(T, F)
) |>
  dplyr::mutate(contact_mat = as.character(contact_mat),
                scenario_name = vaccine_strategies[scenario_index],
                scenario_number = sprintf ("scenario%02d", scenario_index),
                filename = paste0("cbe_", scenario_name, "_", contact_mat, "_", c("stepVE", "linearVE")[step_ve + 1], '.csv'),
                plotname = paste0("cbe_", scenario_name, "_", contact_mat, "_", c("stepVE", "linearVE")[step_ve + 1]))

for (i in 1:nrow(grd)) {

  # deleting the fortran generated files:
  folder_path <- "outcome/"
  txt_files <- list.files(path = folder_path, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
  file.remove(txt_files)

  index <- grd[i, "scenario_index"]
  message(paste0(
    "contact: ", grd[i, 2], ' | ',
    "Step VE: ", grd[i, 3], ' | ',
    "Vaccination: ", grd[i, 4], "\n"
  ))

  # run model and estimate cases
  burden_estimate_file <- runScenario (
    vaccine_coverage_folder    = var$vaccine_coverage_folder,
    vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
    coverage_prefix            = var$coverage_prefix,
    antigen                    = var$antigen,
    scenario_name              = grd[i, "scenario_name"],
    save_scenario              = grd[i, "scenario_number"],
    burden_estimate_folder     = burden_estimate_folder,
    group_name                 = var$group_name,
    log_name                   = var$log_name,
    countries                  = var$countries,
    cluster_cores              = var$cluster_cores,
    psa                        = var$psa,
    vaccination                = set_vaccination [index],
    using_sia                  = set_sia         [index],
    measles_model              = main$fortran_model [imain],
    debug_model                = FALSE,
    contact_mat                = grd[i, "contact_mat"], #main$contact_mat [imain],
    step_ve                    = grd[i, "step_ve"], #main$step_ve     [imain],
    sim_years                  = 2000:2050,

    filename = grd[i, "filename"]
  )


  # estimate deaths and DALYs by Wolfson and Portnoy's methods
  estimate_deaths_dalys (cfr_option             = "Wolfson",
                         burden_estimate_file   = burden_estimate_file,
                         burden_estimate_folder = burden_estimate_folder,
                         psa                    = var$psa)

  estimate_deaths_dalys (cfr_option             = "Portnoy",
                         burden_estimate_file   = burden_estimate_file,
                         burden_estimate_folder = burden_estimate_folder,
                         vimc_scenario          = portnoy_scenarios[index],
                         portnoy_scenario       = "s6",  # constant CFR after 2018
                         psa                    = var$psa)


}




# -------------------------------------------------------------------------

## graph 1
dg1 <- read.csv("vac_coverage_top10/base case/scenarios/routine_mcv2.csv")
g1 <- dg1 |>
  filter(country == "Kenya" & year <= 2023 & year >= 2013) |>
  ggplot() +
  geom_point(data = data.frame(
    vaccine = rep("SIA", 11),
    country_code = rep("KEN", 11),
    country = rep("Kenya", 11),
    year = c(1994, 1999, 2002, 2004, 2005, 2006, 2009, 2012, 2016, 2018, 2021),
    coverage = c(0.01137331, 0.07129982, 0.96028941, 0.03450731, 0.0165366, 0.46507228, 0.84041985, 0.8664371, 1, 0.00917591, 0.5374028)
  ) |> filter(year >= 2013),
  aes(x = as.character(year), y = coverage, shape = "SIA", group = year),
  size = 2.5
  ) +
  geom_line(aes(x = as.character(year), y = coverage, col = vaccine, group = vaccine)) +
  # geom_vline(aes(xintercept = 2023), col = 'black', lty = "dashed") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Measles vaccine coverage by dose (2013-2023)",
       x = "Year", y = "Vaccine coverage",
       col = NULL,
       shape = NULL,
       caption = "Historical coverage data was obtained from WHO databases, \nwhile future trends were assumed to increase
       absolutely by 1% every year from 2023.") +
  theme_bw(base_line_size = 0) +
  theme(legend.position = "bottom")
g1

ggsave("img/vac_coverage.png", plot = g1,
width = 6, height = 6, units = "in", dpi = 300)

g2 <- read.csv("own-reports/country-cases-measles.csv") |>
  transmute(date = lubridate::my(periodname),
         value = Measles) |>
  ggplot(aes(x = date, y = value)) +
  geom_point(alpha = .2) +
  geom_smooth(method = "loess", se = F, span = .1, lwd = .5) +
  labs(title = "Measles cases over time", x = NULL, y = "Count") +
  theme_bw(base_line_size = 0)
g2

ggsave("img/cases_over_time.png", plot = g2,
       width = 6, height = 6, units = "in", dpi = 300)

ggsave("img/kenya_pic.png", plot = g1 | g2,
       width = 12, height = 6, units = "in", dpi = 300)


# -------------------------------------------------------------------------

load("data/data_cfr_portnoy.rda")
load("data/data_cfr_wolfson.rda")

cfr_wolfson <- data_cfr_wolfson |> filter(country == "Kenya") |> pull(CFR)

scenarios <- c("cfr_over5_no-vaccination_s6", "cfr_over5_mcv1-default_s6", "cfr_over5_mcv2-default_s6",
               "cfr_under5_no-vaccination_s6", "cfr_under5_mcv1-default_s6", "cfr_under5_mcv2-default_s6")

cfr_portnoy <- data_cfr_portnoy |>
  filter(country_name == "Kenya") |>
  select(country_name, year, all_of(scenarios)) |>
  pivot_longer(-c(country_name, year)) |>
  separate_wider_delim(name, delim = '_',
                       names = c("type", "age", "vaccine", "index")) |>
  select(-c("country_name", "type", "index")) |>

  mutate(vaccine = str_replace_all(vaccine, "-default", ""),
         vaccine = str_to_upper(vaccine),
         age = ifelse(age == "under5", "Under 5", "Over 5"),
         age = factor(age, levels = c("Under 5", "Over 5")),
         wolfson = ifelse(age == "Under 5", cfr_wolfson / 2, cfr_wolfson))

g2 <- ggplot(cfr_portnoy) +
  geom_line(aes(x = year, y = value, col = vaccine)) +
  geom_hline(aes(yintercept = wolfson), col = 'black')+
  facet_wrap(~age, scales = "free") +
  labs(title = "Case fatality risk", x = "Year", y = "CFR",
       col = NULL,
       caption = "The black solid lines are the time invariant CFR estimates from Wolfson, \nthe rest are the time-varying estimates from Portnoy et al.") +
  theme_bw(base_line_size = 0) +
  theme(legend.position = "bottom")
g2

ggsave("img/cfr.png", plot = g2,
       width = 10, height = 5, units = "in", dpi = 300)


# -------------------------------------------------------------------------

# social contact patterns
fn_cmat <- \(contact_mat, countries = "KEN") {

  amat <- switch (contact_mat,
          "syn"     = sapply (countries,
                              function(cty){data_contact_syn[[cty]]},
                              simplify = FALSE, USE.NAMES = TRUE),
          "polymod" = sapply (countries,
                              function(cty){data_contact_polymod[[cty]]},
                              simplify = FALSE, USE.NAMES = TRUE),
          "prpmix"  = sapply (countries,
                              function(x = NULL){
                                pop <- data_pop |> filter(country == "Kenya" & year == 2020) |> arrange(age_from) |> pull(value);
                                amat <- matrix(NA, nrow = 101, ncol = 101)
                                for (i in 1:length(pop)) {
                                  amat[, i] <- pop[i] / sum(pop)
                                }
                                amat
                              },
                              simplify = FALSE, USE.NAMES = TRUE), # not actually used but to fit in the model structure
          "unimix"  = sapply (countries,
                              function(x = NULL){matrix (1/101, nrow = 101, ncol = 101)},
                              simplify = FALSE, USE.NAMES = TRUE))

  ret <- amat[[1]] |>
    cbind(age=1:101) |>
    data.frame() |>
    pivot_longer(-age) |>
    mutate(name = str_replace_all(name, "V", "") |> as.numeric()) |>
    setNames(c("from", "to", "contact"))
  return(ret)

}

c1 <- fn_cmat(contact_mat = "prpmix", countries = "KEN") |>
  ggplot() +
  geom_tile(aes(x = from, y = to, fill = contact),
            show.legend = F) +
  labs(x = NULL, y = NULL, title = "Proportional (Homogeonous) mixing") +
  theme_minimal(base_line_size = 0)

c2 <- fn_cmat(contact_mat = "polymod", countries = "KEN") |>
  ggplot() +
  geom_tile(aes(x = from, y = to, fill = contact),
            show.legend = F) +
  labs(x = NULL, y = NULL, title = "POLYMOD (Great Britain)") +
  theme_minimal(base_line_size = 0)

c3 <- fn_cmat(contact_mat = "syn", countries = "KEN") |>
  ggplot() +
  geom_tile(aes(x = from, y = to, fill = contact),
            show.legend = F) +
  labs(x = NULL, y = NULL, title = "Synthetic (Prem et al.)") +
  theme_minimal(base_line_size = 0)

c4 <- fn_cmat(contact_mat = "unimix", countries = "KEN") |>
  ggplot() +
  geom_tile(aes(x = from, y = to, fill = contact),
            show.legend = F) +
  labs(x = NULL, y = NULL, title = "Uniform (No age-dependency)") +
  theme_minimal(base_line_size = 0)


ggsave("img/contact_matrices.png", plot = gridExtra::grid.arrange(c2, c3, c4, nrow = 2),
       width = 8, height = 8, units = "in", dpi = 300)


# -------------------------------------------------------------------------


v1 <- data.frame(age = seq(from = 6/12, to = 3, length = 100)) |>
  mutate(`Step function` = ifelse(age < 1, .85, .95),
         `Linear function` = .68 * (1.0149^(0:99)),
         `Linear function` = ifelse(`Linear function` > .98, .98, `Linear function`)) |>
  pivot_longer(-age) |>

  ggplot() +
  geom_line(aes(x = age, y = value, col = name)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "First dose efficacy",
       x = "Age", y = "Vaccine efficacy",
       col = NULL) +
  theme_bw(base_line_size = 0) +
  theme(legend.position = "bottom")

ggsave("img/vaxx_eff.png", plot = v1,
       width = 4, height = 4, units = "in", dpi = 300)


# -------------------------------------------------------------------------

fls.p <- list.files("central_burden_estimate/Portnoy")
fls.w <- list.files("central_burden_estimate/Wolfson")

d1.p <- fls.p |>
  purrr::map(\(fl) {
    read.csv(paste0("central_burden_estimate/Portnoy/", fl)) |>
      dplyr::mutate(file = fl)
  })
d1.w <- fls.w |>
  purrr::map(\(fl) {
    read.csv(paste0("central_burden_estimate/Wolfson/", fl)) |>
      dplyr::mutate(file = fl)
  })


d2 <- dplyr::bind_rows(d1.p, d1.w)
d3 <- d2 |>
  separate_wider_delim(file, delim = '_',
                       names = c("model", "vaccine", "contactmat", "vaxeff", "cfr_method")) |>
  mutate(cfr_method = str_replace_all(cfr_method, ".csv", ""),
         vaccine = str_to_upper(vaccine),
         contactmat = str_to_upper(contactmat)) |>
  select(-model)


d4 <- d3 |>
  group_by(year, vaccine, contactmat, vaxeff, cfr_method) |>
  reframe(across(c(cases, deaths, dalys), sum))

d5 <- d4 |>
  pivot_longer(c(cases, deaths, dalys)) |>
  mutate(name = ifelse(name == "dalys", "DALYs", str_to_title(name)),
         value = value / 1e6,
         vaccine = ifelse(vaccine == "NO-VACC", "No vaccination", vaccine),
         vaccine = factor(vaccine, levels = c("No vaccination", "MCV1", "MCV2"))
         )

m1 <- ggplot(d5 |> filter(vaxeff == "linearVE" &
                            cfr_method == "Portnoy" &
                            contactmat == "SYN")) +
  geom_line(aes(x = year, y = value, col = vaccine)) +
  facet_wrap(name ~ ., scales = "free", nrow = 2) +
  labs(x = "Year", y = "Health burden (millions)",
       title = "Measles vaccination impact",
       col = NULL,
       caption = "Measles cases, deaths and disability-adjusted life years by vaccination strategies and countries.") +
  theme_bw(base_line_size = 0) +
  theme(legend.position = "bottom")
m1

ggsave("img/vaxx_impact.png", plot = m1,
       width = 8, height = 8, units = "in", dpi = 300)


# -------------------------------------------------------------------------

## How the age distribution of cases is affected by contact matrices
struct <- data.frame(age = 0:100,
           # group = c(rep(
           #   c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
           #     "45-49","50-54","55-59","60-64","65-69"),
           #   each = 5
           # ),
           # rep("70+", 31))
           group = c("0", "1", "2", rep("3+", 98))
           )

a1 <- d3 |>
  filter(cfr_method == "Portnoy" & vaxeff == "linearVE") |>
  select(year, age, cases, vaccine, contactmat) |>
  merge(struct, by = "age") |>
  mutate(group = factor(group,
                        levels = c("0", "1", "2", "3+"))) |>
  group_by(year, contactmat, vaccine, group) |>
  reframe(cases = sum(cases)) |>
  group_by(year, contactmat, vaccine) |>
  mutate(dist = cases / sum(cases)) |>
  ungroup() |>
  filter(contactmat != "PRPMIX") |>
  mutate(contactmat = case_when(
    contactmat == "POLYMOD" ~ "POLYMOD (Great Britain)",
    contactmat == "SYN" ~ "Synthetic (Prem et al.)",
    contactmat == "UNIMIX" ~ "Uniform mixing",
    TRUE ~ "Proportional mixing"
  ),
  vaccine = ifelse(vaccine == "NO-VACC", "No vaccinations", vaccine),
  vaccine = factor(vaccine, levels = c("No vaccinations", "MCV1", "MCV2")))

a2 <- a1 |>
  ggplot(aes(x = year, y = dist, fill = group)) +
  geom_area(position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(vaccine ~ contactmat) +
  labs(
    title = "Age distribution of cases",
    x = "Year",
    y = "Percentage distribution",
    fill = NULL
  ) +
  scale_fill_manual(values = c("#27aae1", "#fe7501", "#084887", "#efca08")) +
  theme_minimal(base_line_size = 0) +
  theme(legend.position = "bottom")

a2
ggsave("img/age_dist.png", plot = a2,
       width = 8, height = 8, units = "in", dpi = 300)


# -------------------------------------------------------------------------

ggsave("img/joint_plot.png", plot = (g2 | v1) / (c2 | c3 |  c1) / (c4 | ggplot() + theme_void() | ggplot() + theme_void()),
       width = 8, height = 8, units = "in", dpi = 300)


# Tables ------------------------------------------------------------------

t1 = d4 |>
  select(-year) |>
  group_by(contactmat, vaccine, cfr_method, vaxeff) |>
  reframe(across(everything(), sum)) |>
  mutate(across(where(is.numeric), \(x) round(x/1e6, 1)),
         vaccine = ifelse(vaccine == "NO-VACC", "No vaccination", vaccine),
         vaccine = factor(vaccine, levels = c("No vaccination", "MCV1", "MCV2")))

t2 <- t1 |>
  filter(vaxeff == "linearVE" & contactmat != "PRPMIX" & cfr_method == "Portnoy") |>
  select(Contact = contactmat, `Vaccine strategy` = vaccine,
         Cases = cases, Deaths = deaths, DALYs = dalys)


# Assuming `y` is your dataset
# Pivot the data for desired format
t3 <- t2 |>
  mutate(`Vaccine strategy` = as.character(`Vaccine strategy`),
         `Vaccine strategy` = ifelse(`Vaccine strategy` == "MCV2", "MCV1 & MCV2", `Vaccine strategy`)) |>
  pivot_longer(cols = c(Cases, Deaths, DALYs),
               names_to = "Metric",
               values_to = "Value") |>
  pivot_wider(names_from = c(`Vaccine strategy`, Metric),
              values_from = Value)

# Create the gt table
t4 <- t3 |>
  select(all_of(
    c("Contact", "No vaccination_Cases", "No vaccination_Deaths",
      "No vaccination_DALYs", "MCV1_Cases", "MCV1_Deaths", "MCV1_DALYs", "MCV1 & MCV2_Cases",
      "MCV1 & MCV2_Deaths", "MCV1 & MCV2_DALYs")
  )) |>
  gt(rowname_col = "Contact") |>
  tab_spanner_delim(delim = "_") |>
  tab_options(
    table.font.names = "Times New Roman",
    table.font.size = 12
  )
t4

# -------------------------------------------------------------------------

t2 |>
  filter(`Vaccine strategy` == "No vaccination")

v1 <- t1 |>
  filter(cfr_method == "Portnoy" & contactmat == "SYN" & vaccine != "No vaccination") |>
  pivot_longer(c(cases, deaths, dalys)) |>
  mutate(name = ifelse(name == "dalys", "DALYs", str_to_title(name)),
         vaccine = as.character(vaccine),
         vaccine = ifelse(vaccine == "MCV2", "MCV1 & MCV2", vaccine),
         vaxeff = ifelse(vaxeff == "linearVE", "Linear VE", "Step VE")) |>
  filter(name == "Cases") |>
  ggplot() +
  geom_point(aes(x = vaxeff, y = value,
                 group = vaxeff),
             shape = 15, size = 3) +
  geom_line(aes(x = vaxeff, y = value, group = 1)) +
  labs(title = "Sensitivity to vaccine efficacy",
       y = "Health burden (millions)", x = "Vaccine efficcy function") +
  facet_grid(name ~ vaccine, scales = "free") +
  theme_bw(base_line_size = 0) +
  theme(legend.position = "bottom")
v1
ggsave("img/sens_vax_eff.png", plot = v1,
       width = 10, height = 5, units = "in", dpi = 300)


# County plots ------------------------------------------------------------

require(sf)

cc_mcv <- read.csv("own-reports/county-cases-measles.csv") |>
  select(variable = Data, county = Organisation.unit, year = Period, value = Value) |>
  mutate(county = str_replace_all(county, " County", ""),
         variable = ifelse(variable == "Measles", "Measles cases",
                           ifelse(variable == "Proportion of under 1 year receiving vaccine against Measles and Rubella 1",
                                  "Proportion of under 1 year receiving vaccine\nagainst Measles and Rubella 1",
                                  "Proportion of under two years receiving vaccine\nagainst Measles and Rubella 2")))

shp <- rKenyaCensus::KenyaCounties_SHP |> st_as_sf() |>
  select(county = County) |>
  mutate(county = str_to_title(county)) |>
  mutate(county = case_when(
    county == "Elgeyo/Marakwet" ~ "Elgeyo Marakwet",
    county == "Murang'a" ~ "Muranga",
    county == "Nairobi City" ~ "Nairobi",
    county == "Taita/Taveta" ~ "Taita Taveta",
    county == "Tharaka-Nithi" ~ "Tharaka Nithi",
    TRUE ~ county
  ))

cc_mcv <- merge(shp, cc_mcv, by = "county")

uvals <- cc_mcv |> pull(variable) |> unique()
plts <- list()

for (i in 1:length(uvals)) {
  plts[[i]] <- ggplot(cc_mcv |> filter(variable == uvals[i])) +
    geom_sf(aes(fill = value)) +
    facet_wrap(.~year) +
    labs(title = uvals[i], fill = "Value") +
    theme_void() +
    scale_fill_viridis_c(option = "viridis")  # Use "viridis" color scale
}

ggsave("img/county_plots.png", plot = (plts[[1]] / plts[[2]] / plts[[3]]),
       width = 8, height = 8, units = "in", dpi = 300)


