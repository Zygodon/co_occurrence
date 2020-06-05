library("RMySQL")
library(tidyverse)
library(cooccur)
library(igraph)
library(tidygraph)
library(lsmeans)

# Functions
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}

# General SQL query
query <- function(q)
{
  # Remote DB with password
  con <- dbConnect(MySQL(), 
                   user  = "guest",
                   password    = "guest",
                   dbname="meadows",
                   port = 3306,
                   host   = "sxouse.ddns.net")
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

# Load the database
GetTheData <-  function()
{
  # GET DATA FROM DB
  # Remote DB with password
  con <- dbConnect(MySQL(), 
                   user  = "guest",
                   password    = "guest",
                   dbname="meadows",
                   port = 3306,
                   host   = "sxouse.ddns.net")
  
  
  q <- sprintf('select assembly_id, assembly_name, quadrat_count, community, quadrat_id, quadrat_size, visit_date, records_id, species.species_id, 
    species.species_name from assemblies
      join quadrats on quadrats.assembly_id = assemblies_id
      join visit_dates on quadrats.vd_id = visit_dates.vds_id
      join records on records.quadrat_id = quadrats_id
      join species on species.species_id = records.species_id
      # Two assemblies have 0 quadrat count; exclude A.capillaris_stolonifera;
      # exclude some odd assemblies with no assigned community
    where quadrat_count = 5 and species.species_id != 4 and community is not null
    and quadrat_size = "2x2";') 
  # NOTE: this extract includes "MG5", i.e. some MG5 communities where 
  # the team have not decided
  # on a sub-group.
  
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

the_data <- GetTheData()
# Make the stand (assembly) based occupancy matrix site_occ. Start with d2:
d2 <- (the_data %>% select(assembly_id, species_name)
       %>% group_by(assembly_id, species_name)
       %>% summarise(n=n())
       %>% ungroup()
       %>% pivot_wider(names_from = species_name, values_from = n))
# At this point, d2 has the number of hits for each assembly and species.
# Replace anything numeric with 1, and any NA with 0
d3 <- (d2 %>% select(-assembly_id) %>% replace(., !is.na(.), 1)
       %>% replace(., is.na(.), 0)) # Replace NAs with 0)
# Insert a column with the stand IDs. I know that it is usual with occupancy
# matrices to have the species as rows and sites as columns but as I intend
# analysis by species (R analysis as opposed to Q) I find it easier for
# the moment to think of the species in columns.
cn <- colnames(d3)
d3 <- (mutate(d3, stand = paste("S", as.character(d2$assembly_id), sep = "_"))
       %>% select(stand, cn)) #Re order columns so stand 1st
# Rename and delete d3
site_occ <- d3
rm(the_data,d2, d3)
species <- colnames(select(site_occ, -1))[which(colSums(select(site_occ, -1)) > 20)]
m <- t(site_occ %>% select(-1))

analysis <- cooccur(mat = m, type = "spp_site", thresh = F, spp_names = T)
tab <- as_tibble(prob.table(analysis))
tab <- tab %>% mutate(negative = p_lt < 0.05, positive = p_gt < 0.05, value = "random")
tab$value <- ifelse(tab$positive, "positive", tab$value)
tab$value <- ifelse(tab$negative, "negative", tab$value)

p0 <- ggplot(data = tab, aes(x = exp_cooccur, y = obs_cooccur)) +
  geom_point(aes(colour = value)) +
  labs(title = "Co-occurrence; observed vs expected")
plot(p0)

# At the moment, everything is "connected" to everything else just because it is there in the table
# Remove disconnected pairs to make a non-densely connected graph.
g_data <- tab %>% filter(obs_cooccur > 0)
vertices <- (g_data %>% select(sp1_name, sp2_name)
             %>% pivot_longer(cols = c(sp1_name, sp2_name),values_to = "species")
             %>% distinct(species)) # vertices: just a list of species

# Add hits to vertices
sp1 <- g_data %>% distinct(sp1_name, sp1_inc) %>% rename(species = sp1_name, hits = sp1_inc)
sp2 <- g_data %>% distinct(sp2_name, sp2_inc) %>% rename(species = sp2_name, hits = sp2_inc)
hits <- full_join(sp1, sp2)
vertices <- vertices %>% left_join(hits, by = "species")
rm(sp1, sp2, hits)

# Make a graph in order to calculate the degree of each speices
# (degre: how many other species each one co-occurs with)
G0 <- tbl_graph(nodes = vertices, 
                edges = g_data %>% select("sp1_name", "sp2_name", everything()))
G0 <- G0 %>% activate(nodes) %>% mutate(deg = centrality_degree())
# Re-do vertices to include degree
rm(vertices)
vertices <- G0 %>% activate(nodes) %>% as_tibble(.)

# Exponential plot, degree ~ hits
p1 <- ggplot(vertices, aes(hits, deg)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Degree vs hits")
plot(p1)

# Transform log(degree) - looks like capacitor charging d = Dmax(1-exp(h0*h))
# Where h0 is an inverse "hit rate parameter"
logplot <- tibble(h = vertices$hits, deg = vertices$deg)
# logplot <- logplot %>% mutate(d = -log(1 - (vertices$deg/length(vertices$species))))
logplot <- logplot %>% mutate(d = log(length(vertices$species) - deg))

# One parameter model - h0 is the slope
p2 <- ggplot(logplot, aes(h, d)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, colour = "grey50") +
  labs(title = "One parameter model")
plot(p2)

# But it looks like the slope breaks at about h = 30. 
# Is there evidence for a "fast" hitrate constant h1 and a "slow" hitrate constant h2?

# model 1: 1 independent variable (only h). The one parameter model
fit_one <- lm(d~h, data=logplot)
h0 <- fit_one$coefficients[2] # the one-parameter hitrate constant

# Assume after h = 30, the slope is dominated by h2, th slow hitrate constant
# Add a factor to divide the two parts of the relationship
logplot <- logplot %>% mutate(h_const = ifelse(h > 30, "slow", "fast"))
p3 <- ggplot(logplot, aes(h, d)) +
  geom_point(aes(colour = h_const)) +
  geom_smooth(data = subset(logplot, h > 30), method = "lm", se = F) +
  labs(title = "Points dominated by 'slow' coefficient")
plot(p3)

# Extract the "slow" hitrate constant using points hit rate > 30
fit_slow <- lm(d~h, data=logplot, subset=(h > 30))
h2 <- coefficients(fit_slow)[2]

# Are the slopes of the two parts of the graph different
# fit_fast <- lm(d~h, data=logplot, subset=(h <= 30))
hhc_interaction <- lm(d ~ h * h_const, data = logplot)
an <- anova(hhc_interaction) # an suggests the h_h_const interaction is significant

# Obtain slopes
hhc_interaction$coefficients
f_s_lst <- lstrends(hhc_interaction, "h_const", var="h")
# Compare slopes
pairs(f_s_lst) # Yes the slopes are different, p.value 0.0011

# Is the 2-parameter ( h and h_const) model better?
# Model 2: 2 IVs (h and h_const)
fit_two <- lm(d ~ h * h_const, data = logplot)
# Compare model 1 to model 2
anova(fit_one, fit_two) # yes they are different, but which is best?
# anova table suggests fit_2 beter than fit_one,
# Pr(>F) 0.002148,
# Resdidual SS 5.35 for fit_one cf 5.03 for fit_two
# see https://bookdown.org/ndphillips/YaRrr/comparing-regression-models-with-anova.html

# So it's worth stripping out the slow constant h2 to leave the fast constant h1
# slow_const <- summary(f_s_lst)$h.trend[2] # Checking: same as h2 above
logplot <- logplot %>% mutate(fast_d = d + (h*h2))
fit_fast <- lm(fast_d~h, data = subset(logplot, h <= 30))
h1 <- fit_fast$coefficients[2] # the "fast" hit rate constant

# And plotting the points dominated by the fast constant separately
p5 <- ggplot(logplot, aes(h, d)) +
  geom_point(aes(colour = h_const)) +
  geom_point(data = subset(logplot, h <= 30), aes(h, fast_d), shape = 1, colour = "grey50") +
  geom_smooth(data = subset(logplot, h > 30), aes(colour = h_const), method = "lm", se = F) +
  geom_smooth(data = subset(logplot, h <= 30), aes(h, fast_d), method = "lm", se = F, colour = "grey50") +
  labs(title = "Fast and slow components separated")
plot(p5)

# But the really important question is why we see this 
# exponential relationship in the first place.

