# R is a programming language, the R console is the primary interface for executing R code

#RStudio is an IDE or Integrated Development Environment.
####It's more user-friendly and allows users to edit code, look at console output, visualize data and plots, manage packages

#Tour RStudio, 4 panes--Console (where code is executed), Script (where you write code), Environment/History,  Files/Plots/Packages/Help/Viewer

#Basic code execution, click "Run", use Ctrl+enter
1+1 # notice a few things here--symbols appear in different colors, annotations start with "#", answer displays in console

#To store objects, we use vectors. We establish a vector with a "<-"
calculation <- 1+1 #appears in environment and is stored in this session, we can highlight and run "calculation"

# At the beginning of a session/project, we set a working directory. This will house the files we need to work with in R and output we create in R.
# Can do this manually, or with code.If following along, go ahead and create this folder on your desktop.
setwd("C:/Users/Brittany.Bankovich/Desktop/RWorkshop") #anything in this folder will display in the files/plots pane
#setwd is what we call a "function." There are thousands of R functions, some exist in base R (what we call what comes standard with the download)
# We can look at function details by typing ?function
?setwd #This is one of the most helpful pieces of code, it allows you to see everything that needs to go into your function (arguments/parameters)

# Package basics and loading packages
# We'll start with tidyverse, which is a newish package (2016 is new for me, back off); it's useful for importing, cleaning, transforming, visualizing, and modeling data in R. 

#First we need to install tidyverse
install.packages('tidyverse', dep=TRUE) #dep=TRUE tells R to install all dependencies, or other packages that integrate with our package of interest. This will usually save you some frustration later (BUT NOT ALWAYS)
#Let's try a function in this package.
tibble(calculation)

# Why didn't this work? We need to tell R we want to use this package in this session, we do this with library(). We need to do this every session.
library(tidyverse) #Can also click it on under "packages" tab
tibble(calculation)
?tidyverse
browseVignettes("tidyverse")
# We can write this to a file, by default it's going to store in our working directory
write.csv(calculation, "calculation.csv") #refresh files, it shows up
# We can read it back in
read.csv("calculation.csv")
# A lot of people like read_csv() instead as it will read as a tibble, for this you'll need to install 'readr' and load the 'readr' library
install.packages('readr', dep=TRUE)
library(readr)
read_csv("calculation.csv")

# Most packages contain some sample data so you can run through example script. We're going to use sample data in the palmerpenguins package. 
install.packages('palmerpenguins', dep=TRUE)
library(palmerpenguins)
# Tell R to load and store sample dataset
data(penguins)
?palmerpenguins
# head() is a good function to look at a few rows of data
head(penguins)
# with View() you can look at the entire dataset
View(penguins)
# and if you need to explore field types, use str(); this will often let you know if R imported the fields as expected or if you'll need to change them
str(penguins) #for example, maybe we need year to be a factor instead of an integer
# We can alter field types; we reference data columns with '$'
penguins$year <- as.factor(penguins$year)
str(penguins)
# the summary() function is useful for quickly identifying outliers/typos; maybe there shouldn't be any NAs, or maybe we know bill depth shouldn't exceed 20 mm
summary(penguins)
# Tidy basics--pipes. Think of pipes like a physical pipe, data passes through a pipe and undergoes a process and a new dataset emerges
penguins |>  # Ctrl + shift + m shortcut; this is the native pipe, there is also a pipe associated with the 'magrittr' package that looks like "%>%". You can choose your pipe in Tools → Global Options → Code → Use native pipe operator
filter(species == "Adelie") # filter() chooses rows; when referencing text we use "==", if we're referencing a number we use "="

#We can save tidy operations to new vectors
species_mass <- penguins |> 
  select(body_mass_g, species) # select chooses columns, you can also use it to reorder based on the order in which you type them 
  
#Add new columns with mutate()
penguins |> 
  mutate(body_mass_kg = body_mass_g / 1000)

#Calculate some summary statistics with group_by() and summarise()
penguins |> 
  group_by(species) |> 
  summarise(mean_bill_length = mean(bill_length_mm, na.rm = TRUE)) #na.rm tells R to remove NAs from the calculation; now we are starting to see how parentheses come into play

#And do a bunch of tidy operations in sequence
penguins_summary <- penguins |> 
  # Select only the columns we care about
  select(species, island, body_mass_g, sex) |> 
  # Keep only female penguins on "Biscoe" island
  filter(sex == "female", island == "Biscoe") |> 
  # Convert body mass to kg
  mutate(body_mass_kg = body_mass_g / 1000) |> 
  # Group by species
  group_by(species) |> 
  # Calculate summary statistics
  summarise(
    mean_mass_kg = mean(body_mass_kg, na.rm = TRUE), #mean mass in kg
    max_mass_kg = max(body_mass_kg, na.rm = TRUE), #maximum mass in kg
    n = n() #give us a count
  )
# You can start to see how this could be very useful after initial setup when bringing in new data over time

# Let's look at a dataset that contains some common inconsistencies re: data entry
# The 'messy' package lets us take a clean dataset and create inconsistencies for practice
install.packages("messy")
library(messy)

set.seed(123) # this is a useful function that allows for reproducibility with random processes (e.g., you want to generate the same random data each time)
data(penguins, package = "palmerpenguins")
messy_penguins <- messy(penguins)
View(messy_penguins) # Can you identify different issues with the data? List a few.

# We'll use 'stringr', which should already be loaded with tidyverse, but no harm in reloading it

library(stringr)

# stringr is for: 
#1 Character manipulation: these functions allow you to manipulate individual characters within the strings in character vectors.
#2 Whitespace tools to add, remove, and manipulate whitespace.
#3 Locale sensitive operations whose operations will vary from locale to locale (in the world).
#4 Pattern matching functions. These recognise four engines of pattern description. The most common is regular expressions (regex), but there are three other tools.

# Work in sequence
clean_penguins_text <- messy_penguins |> 
  mutate(
    # Step 1: convert everything to lowercase
    species = str_to_lower(species),
    island  = str_to_lower(island),
    
    # Step 2: remove all non-letter characters
    species = str_replace_all(species, "[^a-z]", ""), #[^a-z] → any character that is NOT a lowercase letter, ^=not, replace with nothing
    island  = str_replace_all(island, "[^a-z]", ""),
    
    # Step 3: convert to proper title case
    species = str_to_title(species),
    island  = str_to_title(island),
    
    # Step 4: clean sex column
    sex = str_to_lower(sex),
    sex = str_replace_all(sex, "[^a-z]", ""),
    sex = case_when(
      sex %in% c("male", "m") ~ "male", # after removing non-letters and converting to lowercase, checks to see if value is "male" or "m", and if so, forces value to be "male"
      sex %in% c("female", "f") ~ "female",
      TRUE ~ NA_character_ #if value isn't "m","f","male","female" convert to NA
    )
  )

# Now we'll fix issues in the numeric columns 
clean_penguins_text_numbers <- clean_penguins_text |>  #passing through the dataset with the cleaned text from above
  mutate(
    body_mass_kg = as.numeric(body_mass_g) / 1000,  # convert to kg
    flipper_length_mm = as.numeric(flipper_length_mm),
    bill_length_mm = as.numeric(bill_length_mm),
    bill_depth_mm = as.numeric(bill_depth_mm)
  )

View(clean_penguins_text_numbers)# looks pretty good now

############# Short break? #############

# Let's move on to some visuals
# The package we're going to be using is 'ggplot2'. It's based on a book called The Grammar of Graphics, which treats data graphics like a structured language. Figures are created using systematically combined elements. Here is a great cheatsheet: https://github.com/rstudio/cheatsheets/blob/main/data-visualization.pdf
install.packages('ggplot2', dep=TRUE)
library(ggplot2)

# We'll go back to the penguins dataset and create a very basic ggplot

ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g)) + #aes() tells R which columns to plot; notice we build with "+"
  geom_point() 

# We can save our plots as vectors
simple_plot <- ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g)) + #aes() tells R which columns to plot; notice we build with "+"
  geom_point() 
simple_plot

# This is pretty basic and would definitely look better with some labels 
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g)) +  # map columns
  geom_point() +                                                  # add points
  labs(
    title = "Penguins: Flipper Length vs. Body Mass",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)"
  )

# Say we want to color the points by a variable (sex), we can simply add this to the aes() argument
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color= sex)) +  # sex informs number of colors to plot
  geom_point() +                                                
  labs(
    title = "Penguins: Flipper Length vs. Body Mass",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)"
  )

# One step further, we can change point shapes by species
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color= sex, shape=species)) +  # species informs number of shapes to plot
  geom_point() +                                                 
  labs(
    title = "Penguins: Flipper Length vs. Body Mass",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)"
  )

# What if I don't like the default colors? We can adjust with a new element
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color= sex, shape=species)) +  
  geom_point() +                                                 
  labs(
    title = "Penguins: Flipper Length vs. Body Mass",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)"
  )+
  scale_color_manual(values = c("male" = "blue", "female" = "hotpink")) #you can write colors by name or use hex codes (see list of colors here: https://r-charts.com/colors/)

# What about the point shapes, can I change those too?
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color= sex, shape=species)) +  
  geom_point() +                                                 
  labs(
    title = "Penguins: Flipper Length vs. Body Mass",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)")+
  scale_color_manual(values = c("male" = "blue", "female" = "hotpink"))+
  scale_shape_manual(values = c("Adelie" = 18,  "Chinstrap" = 8, "Gentoo" = 0)) #https://r-charts.com/base-r/pch-symbols/ 

# We can even add a trendline 
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color= sex, shape=species)) +  
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black")+ 
  labs(
    title = "Penguins: Flipper Length vs. Body Mass",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)")+
  scale_color_manual(values = c("male" = "blue", "female" = "hotpink"))+
  scale_shape_manual(values = c("Adelie" = 18,  "Chinstrap" = 8, "Gentoo" = 0))
#And adjust axis scales
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color= sex, shape=species)) +  
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black")+ 
  labs(
    title = "Penguins: Flipper Length vs. Body Mass",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)")+
  scale_color_manual(values = c("male" = "blue", "female" = "hotpink"))+
  scale_shape_manual(values = c("Adelie" = 18,  "Chinstrap" = 8, "Gentoo" = 0))+
  scale_x_continuous(limits = c(0, 250)) +
  scale_y_continuous(limits = c(0, 6500))
# We can begin adjusting text, background, grid, legend, axes (anything non-data) with theme()
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color= sex, shape=species)) +  
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black")+ #linear model with confidence intervals
  labs(
    title = "Penguins: Flipper Length vs. Body Mass",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)")+
  scale_color_manual(values = c("male" = "blue", "female" = "hotpink"))+
  scale_shape_manual(values = c("Adelie" = 18,  "Chinstrap" = 8, "Gentoo" = 0))+
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(size = 20, face = "bold"),       # title size and bold
    axis.title = element_text(size = 16),                      # axis label size
    axis.text = element_text(size = 14),                       # axis tick labels
    legend.title = element_text(size = 14),                    # legend title
    legend.text = element_text(size = 12)                      # legend items
  )


## Let's explore different plot types 
# Histogram
ggplot(penguins, aes(x = body_mass_g)) +
  geom_histogram()
# Histogram by group
ggplot(penguins, aes(x = body_mass_g, fill = sex)) +
  geom_histogram(alpha = 0.6) # alpha controls transparency

# Box plot
ggplot(penguins, aes(x = species, y = body_mass_g)) +
  geom_boxplot()

# Violin plot
ggplot(penguins, aes(x = species, y = body_mass_g)) +
  geom_violin()
# Bar plot
ggplot(penguins, aes(x = species)) +
  geom_bar()

#Sometimes it's useful to create multiple plots and visualize side-by-side, we control this with facet_wrap() in ggplot2
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color = species)) +
  geom_point() +
  facet_wrap(~ species)

# Journals always have figure formatting requirements, so it's useful to export our figures to see how they will display at set resolution and set
# dimensions
# First, let's store the plot as a vector
penguin_figure_for_ms <- ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color= sex, shape=species)) +  
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black")+ #linear model with confidence intervals
  labs(
    title = "Penguins: Flipper Length vs. Body Mass",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)")+
  scale_color_manual(values = c("male" = "blue", "female" = "hotpink"))+
  scale_shape_manual(values = c("Adelie" = 18,  "Chinstrap" = 8, "Gentoo" = 0))+
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(size = 20, face = "bold"),       # title size and bold
    axis.title = element_text(size = 16),                      # axis label size
    axis.text = element_text(size = 14),                       # axis tick labels
    legend.title = element_text(size = 14),                    # legend title
    legend.text = element_text(size = 12)                      # legend items
  )

# Say the journal you're submitting to wants figures to be .jpeg, 180 mm wide
ggsave("figure1.png", plot = penguin_figure_for_ms,  width = 7.09, height = 5, dpi = 300) #units are in inches, so convert mm to in


##### Basic statistics ##### 
## Putting it ALL together ###
# Base R has some basic statistical tests and models
# lm = linear model, cor.test = correlation, t.test = t-test

# 'car' package extends this to ANOVA and includes tools for testing for multicollinearity and exploring regression diagnostics
# 'MASS' includes stepwise regression tools
# 'lme4' extends regression to mixed-effects models

install.packages('car', dep=TRUE)
install.packages('MASS', dep=TRUE)
install.packages('lme4', dep=TRUE)

library(car)
library(MASS)
library(lme4)

# Test for correlation
cor.test(penguins$flipper_length_mm, penguins$body_mass_g, use = "complete.obs")
# Two-sample t-test
# Isolate one species of penguin
adelie <- penguins |> 
  filter(species == "Adelie")
# Does body mass differ significantly between the sexes?
t.test(body_mass_g ~ sex, data = adelie, var.equal = TRUE, na.rm = TRUE)

# Type II ANOVA with 'car'
# Remove NAs first*
penguins_clean <- penguins |> 
  filter(!is.na(body_mass_g), !is.na(species))
# 'car' takes lm in Anova() function
mod <- lm(body_mass_g ~ species, data = penguins_clean) 
Anova(mod, type = 2)
# Post-hoc test; Tukey test
tukey <- TukeyHSD(aov(mod))
# Convert TukeyHSD results to data frame
tukey_df <- as.data.frame(tukey$species)
tukey_df$comparison <- rownames(tukey_df)
# Keep only columns we need
tukey_df <- tukey_df %>% 
  select(comparison, diff, lwr, upr, `p adj`)
#Make a plot to visualize comparisons
ggplot(tukey_df, aes(x = comparison, y = diff)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Tukey HSD: Pairwise Differences in Penguin Body Mass",
    x = "Species Comparison",
    y = "Difference in Body Mass (g)"
  ) +
  theme_minimal() #we can see Gentoo-Adelie and Gentoo-Chinstrap do not cross zero, thus are statistically different (Gentoo penguins are heavier than Adelie and Chinstrap penguins)

# Stepwise with MASS
# We need to remove all missing values first with tidy/dplyr
penguins_stepwise <- penguins |> 
  select(body_mass_g,
         flipper_length_mm,
         bill_length_mm,
         bill_depth_mm,
         sex,
         species) |> 
  na.omit() ## Why is this throwing an error? This is a frustrating R quirk--sometimes, different packages contain functions with the same names and it overrides the function you are trying to use. In this case, 'MASS' has a function called select(). We need to tell R we want to use the select() function in 'dplyr'. We do this using the format package::function

#You can check to see which package/function is being used by default:
select
# Here we'll specify we want to use dplyr::select. You could also just reload 'dplyr'/'tidyverse' library AFTER 'MASS'
penguins_stepwise <- penguins |> 
  dplyr::select(body_mass_g,
         flipper_length_mm,
         bill_length_mm,
         bill_depth_mm,
         sex,
         species) |> 
  na.omit()
# Start with the full model (every possible variable)
full_mod <- lm(
  body_mass_g ~ flipper_length_mm +
    bill_length_mm +
    bill_depth_mm +
    sex +
    species,
  data = penguins_stepwise
)
# Now define the null/intercept model
null_mod <- lm(body_mass_g ~ 1, data = penguins_stepwise)
# Perform the stepwise regression
step_mod <- stepAIC(
  null_mod,
  scope = list(lower = null_mod, upper = full_mod),
  direction = "both",
  trace = TRUE
)
# See results
summary(step_mod)# kept all variables
# It can be a little challenging to interpret this output given that there are multiple factorial variables, sometimes it's good to visualize it in a plot to help understand magnitude and direction of relationships
#'broom' can help us take messy statistical output and turn them into tibbles to make visualization clearer
install.packages('broom', dep=TRUE)
library(broom)

coef_df <- tidy(step_mod, conf.int = TRUE) |> 
  filter(term != "(Intercept)") |>    # drop intercept
  mutate(term = reorder(term, estimate))

#plot with ggplot2
ggplot(coef_df, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    width = 0.2,
    orientation = "y"
  ) +
  labs(
    title = "Effects on Penguin Body Mass",
    x = "Change in body mass (g)",
    y = ""
  ) +
  theme_minimal(base_size = 14) # This is great, but because Adelie is the intercept/reference level, we're only seeing how other species deviate from Adelie
# We can generate predicted means by species to visualize all of them together with the 'emmeans' package (estimated marginal means)
install.packages('emmeans', dep=TRUE)
library(emmeans)
# Estimated marginal means for species
emm_species <- emmeans(step_mod, ~ species)
# Convert to data frame
emm_df <- as.data.frame(emm_species)
# Plot
ggplot(emm_df, aes(x = species, y = emmean)) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15
  ) +
  labs(
    title = "Predicted Body Mass by Species",
    y = "Body mass (g)",
    x = "Species"
  ) +
  theme_minimal(base_size = 14)
# Or look at predictions by variable (can change it up)
ggplot(penguins_stepwise,
       aes(flipper_length_mm, body_mass_g, color = species)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Model-Predicted Body Mass vs Flipper Length",
    x = "Flipper length (mm)",
    y = "Body mass (g)"
  ) +
  theme_minimal(base_size = 14)

## You'll likely want to explore some model diagnostics to make sure your inferences are valid
# basic residual plots
par(mfrow = c(2,2)) #this changes the layout of the plot window, right now I am telling it to go from a single plot to a 2x2 plot
lm_model <- lm(body_mass_g ~ flipper_length_mm + sex + species + bill_depth_mm + bill_length_mm, 
               data = penguins_stepwise)
plot(lm_model)
#1 residuals vs. fitted, checks for non-linearity and unequal variance (heteroscedasticity)
#2 normal q-q plot, checks for normality of residuals
#3 scale-location another heteroscedasticity check
#4 residuals vs. leverage, helps to identify influential points (possible outliers/data entry issues), problem points are labeled
# Diagnostics are generally fine, a few points may be worth investigating or removing from the dataset (282, 35, 164, 76, 314)

### Mixed-effects models with lme4 ###
# Linear model with island as random effect (how does flipper length predict body mass, each island can have its own baseline body mass)
mod <- lmer(body_mass_g ~ flipper_length_mm + (1 | island), data = penguins)
summary(mod)# between island variation (16817) lower than within/Residual (146471); for every 1 mm increase in flipper length we expect a +45g increase in body mass after accounting for differences by island

# Generalized linear model with logit link (what's the probability a penguin is male given flipper length, with a random intercept for species)
mod2 <- glmer(sex ~ flipper_length_mm + (1 | species), 
              data = penguins, 
              family = binomial)
summary(mod2) #random effects: species matters, fixed effects: longer flippers = greater odds of being a male penguin
exp(.20689) # odds ratio = 1.23, for every 1mm increase in flipper length the odds of being male increases by ~23%