##Setup

library(tidyverse); library(ggbeeswarm); library(lmerTest); library(janitor); library(here); library(pwr); library(tidyr); library(influence.ME); library(performance); library(see); library(emmeans)

theme_Publication <- function(base_size=18, base_family="Helvetica") {
  library(grid); library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
            text = element_text(),panel.background = element_rect(colour = NA), plot.background =
              element_rect(colour = NA),
            panel.border = element_rect(colour = NA),axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), axis.line = element_line(colour="black"),axis.ticks = element_line(),
            panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            legend.position = "right", legend.direction = "vertical", legend.key.size= unit(0.8, "cm"),
            legend.title = element_text(face="bold", size = rel(0.8)), legend.key = element_rect(colour = NA), 
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),strip.text = element_text(face="bold")))
}

rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}

format_p <- function(pval){
  if(pval < 0.001){"< .001"} else if(pval < 0.009){str_remove(round(pval,3), "^0+")} else {str_remove(round(pval,2), "^0+")}
}
 

#Load sample data

file_path <- #file_path
  file_list <- list.files(file_path, pattern = "\\.csv$", full.names = TRUE)
#previous df hidden for privacy
file_list <- file_list[grep(paste(participants$Participant_ID, collapse = "|"), file_list)]

# Some participants have the left eye recorded and some have the right eye recorded (EyeLink switches when infant moves), putting in filters so that only areas of interest are in the dataframe.

raw_data <- file_list %>%
  map(function(file) {
    data <- read.csv(file, header = TRUE)
    filtered_data <- data %>%
      mutate(across(c(LEFT_GAZE_X, LEFT_GAZE_Y, RIGHT_GAZE_X, RIGHT_GAZE_Y), as.numeric),
             GAZE_X = ifelse(is.na(LEFT_GAZE_X), RIGHT_GAZE_X, LEFT_GAZE_X),
             GAZE_Y = ifelse(is.na(LEFT_GAZE_Y), RIGHT_GAZE_Y, LEFT_GAZE_Y),
             pupil_size = ifelse(is.na(LEFT_pupil_size), RIGHT_pupil_size, LEFT_pupil_size)) %>%
      filter(GAZE_X > 0 & GAZE_Y > 0 & GAZE_X < 1024 & GAZE_Y < 768)
    
    # Print summary of filtered data
    print(paste("File:", file))
    return(filtered_data)
  }) %>%
  reduce(rbind)


print(unique(raw_data$subject_id))

#selecting variables we want
raw_data_select = raw_data %>%
  select(subject_id, TRIAL_NUMBER, condition, phase, version, GAZE_X, GAZE_Y, pupil_size, TIME_ACCUMULATED, timestamp, sample_message) #here we still have the trails for which the attention getter was skipped.
#keeping time accumulated here so I can do a sanity check later/see if what I calculate is close. 

#make empty df
new_time_clean <- data.frame()

# Reset so all time stamps start at 0 when the video starts (after attention grabber)
for (this_sub in unique(raw_data_select$subject_id)) {
  for (this_trial in unique(raw_data_select$TRIAL_NUMBER)) {
    curr_trial <- raw_data_select %>%
      filter(subject_id == this_sub & TRIAL_NUMBER == this_trial)
    
    # Find the first occurrence of "DISPLAY_VIDEO"in sample message
    start_index <- which(curr_trial$sample_message == "DISPLAY_VIDEO")[1]
    
    # If "DISPLAY_VIDEO" not found, look for "INVISIBLE_BOUNDARY_ON_VIDEO;NULL_ACTION_ON"
    if (is.na(start_index)) {
      start_index <- which(curr_trial$sample_message == "INVISIBLE_BOUNDARY_ON_VIDEO;NULL_ACTION_ON")[1]
    }
    
    # Skip if neither message is found in the trial
    if (is.na(start_index)) {
      next
    }
    
    start_time <- curr_trial$timestamp[start_index]  # What time does the trial start
    curr_trial$new_time <- curr_trial$timestamp - start_time  #Calculate time stamps
    
    if (this_trial == unique(raw_data_select$TRIAL_NUMBER)[1] & 
        this_sub == unique(raw_data_select$subject_id)[1]) {
      new_time_clean <- curr_trial
    } else {
      new_time_clean <- bind_rows(new_time_clean, curr_trial)  # Merge to create new dataframe
    }
  }
}

new_time_clean <- new_time_clean %>%
  mutate(new_time_s = new_time/1000,
         Participant_ID = sub("^[^0-9]*([0-9]+)$", "\\1", subject_id))


#The number of non-negative values in each trial (so we can calculate the time since each sample is taken for 0.002 seconds) to get a more accurate trial looking time (EyeLink can project looks outside of screen coordinates so we need to recalculate here)
looking_by_trial <- new_time_clean %>%
  group_by(Participant_ID, TRIAL_NUMBER) %>%
  summarize(
    non_negative_count = sum(new_time_s >= 0),
    phase = first(phase),
    condition = first(condition),
    version = first(version),
    TRIAL_NUMBER = first(TRIAL_NUMBER)
  ) %>%
  mutate(
    look_time = non_negative_count * 0.002
  )
 

###Test trial looking

# Test trials only
test_trials <- looking_by_trial %>%
  filter(phase == "test") %>%
  group_by(Participant_ID, TRIAL_NUMBER) %>%
  summarise(version = first(version),
            look_time = first(look_time))

# Marking habituation condition
test_trials <- test_trials %>%
  mutate(condition = case_when(
    version %in% c("A1", "A2") ~ "LR",
    version %in% c("B1","B2") ~ "random"))

# Marking test trial condition based on habituation condition
test_trials <- test_trials %>%
  arrange(Participant_ID, TRIAL_NUMBER) %>%
  group_by(Participant_ID) %>%
  mutate(test_condition = case_when(
    version[1] == 'A1' ~ rep(c('congruent', 'incongruent', 'congruent', 'incongruent'), length.out = n())[1:n()],
    version[1] == 'A2' ~ rep(c('incongruent', 'congruent', 'incongruent', 'congruent'), length.out = n())[1:n()],
    version[1] == 'B1' ~ rep(c('congruent', 'incongruent', 'congruent', 'incongruent'), length.out = n())[1:n()],
    version[1] == 'B2' ~ rep(c('incongruent', 'congruent', 'incongruent', 'congruent'), length.out = n())[1:n()]
  ),
  pair = case_when(
    row_number() %in% 1:2 ~ 1, #first pair
    row_number() %in% 3:4 ~ 2  #second pair
  )) %>%
  ungroup()

test_trials <- rename(test_trials, hab_condition = condition)
test_trials$pair <- as.factor(test_trials$pair)
 
#Exclusion
###Data quality: remove participants with low quality

# Identify participants with 50% or more trials with less than 5 seconds of looking to the test trials
exclusion_data_quality <- test_trials %>%
  group_by(Participant_ID, pair) %>%
  filter(look_time < 5) %>%
  summarize(low_trial_count = n(),
            TRIAL_NUMBER = TRIAL_NUMBER, 
            look_time = look_time) %>%
  ungroup()

testPair_included <- test_trials %>% 
  anti_join(exclusion_data_quality, by = c("Participant_ID", "pair"))
 

###Outliers, using Cook's Distance
#using the influence.ME package, cooks distance method:
lmermodelInitial <- lmer(look_time ~ hab_condition*test_condition + pair + (1|Participant_ID), data = testPair_included)
summary(lmermodelInitial)

alt.est <- influence(lmermodelInitial, group = "Participant_ID")
cooks_values <- cooks.distance(alt.est)

cooks_values.df <- data.frame(
  Participant_ID = unique(testPair_included$Participant_ID),
  cooks_values = cooks_values
)
# An observation with Cook’s distance larger than three times the mean will be marked an outlier.
mean_cooks_d <- mean(cooks_values)

# Flag and remove outliers
cooks_influential_points <- cooks_values > (3 * mean_cooks_d)

# Visualize all points
plot(influence(lmermodelInitial, "Participant_ID"), which = "cook",
     cutoff = (3 * mean_cooks_d), sort = TRUE, 
     xlab = "Cook's Distance",
     ylab = "Subject ID")

cooks_values.df$Outlier <- cooks_influential_points 
cooks_included<- subset(cooks_values.df, Outlier == FALSE) # remove  outliers

#remove all outliers from all dataframes
testPair_included <- testPair_included  %>% 
  filter(Participant_ID %in% cooks_included$Participant_ID)

participants <- participants %>% 
  filter(Participant_ID %in% testPair_included$Participant_ID)

# check how many participants left
length(unique(test_trials$Participant_ID))
length(unique(participants$Participant_ID))

# write final df of included participants for later demographics analyses
write.csv(participants, "includedParticipants.csv")
 
#Data exploration
###Habituation trials

#Do an initial look:
# All habituation trials
hab_trials <- looking_by_trial %>%
  filter(phase == "habituation") %>%
  group_by(Participant_ID, TRIAL_NUMBER) %>%
  summarise(version = first(version),
            look_time = first(look_time))

# What do habituation trials look like?
hab_looking <- hab_trials %>%
  group_by(Participant_ID) %>%
  summarize(mean_look = mean(look_time), se = sd(look_time)/sqrt(n())) %>%
  mutate(se_lo = mean_look - se,
         se_hi = mean_look + se)

# Are babies completing a different number of trials based on habituation condition?
hab_mean <- hab_condition %>%
 mutate(condition = case_when(
    version %in% c("A1", "A2") ~ "LR",
    version %in% c("B1","B2") ~ "random")) %>%
  group_by(condition) %>%
  summarize(mean_hab = mean(as.numeric(num_hab)))

# Are more babies habituating in some versions?
hab_count <- hab_trials %>%
  group_by(version) %>%
  count(num_hab < 14) 

 

###Test trials

# Looking times to incongruent and congruent test trials
looks_to_trialType <- testPair_included %>%
  select(Participant_ID, look_time, hab_condition, test_condition, pair) %>%
  group_by(Participant_ID, pair) %>%
  pivot_wider(names_from = test_condition, values_from = look_time) %>%
  mutate(inc_longer = incongruent > congruent)

# How many babies in each condition were looking longer at incongruent test trials
group_cond_looking <- looks_to_trialType %>%
  group_by(hab_condition, pair) %>%
  summarize(
    total_participants = n_distinct(Participant_ID),
    num_looking_longer = sum(inc_longer, na.rm = TRUE))

 
#Stats and plots

# Looking time to test trials, by habituation condition. 
#group-level: mean look duration to congruent trials (avg by participant), grouped by habituation condition
all_test_looking <- testPair_included %>%
  na.exclude(testPair_included) %>%
  group_by(hab_condition, test_condition) %>%
  summarise(looking = mean(look_time), sd = sd(look_time), se = sd(look_time)/sqrt(n())) %>%
  mutate(se_lo = looking - se,
         se_hi = looking + se)

# Make sure data types are what we want them to be
testPair_included <- testPair_included %>%
  mutate(across(c(hab_condition, test_condition, Participant_ID), as.factor))

#1. Logistic mixed-effects model:
modelS <- lmer(look_time ~ hab_condition*test_condition + pair + (1|Participant_ID), data = testPair_included)
summary(modelS)
reflmermodelS <- summary(modelS)

# Visualize:
# Test pairs averaged across test condition:
testLooking_plot <- testPair_included %>%
  group_by(Participant_ID, test_condition) %>%
  mutate(avglook_time =  mean(look_time)) %>%
  slice_head() %>%
  select(-c("TRIAL_NUMBER","version","pair", "look_time")) %>%
  pivot_wider(names_from = test_condition, values_from = avglook_time) %>%
  mutate(inc_longer = incongruent > congruent) %>%
  mutate(inc_longer = case_when( 
           inc_longer == TRUE ~ 1,
         inc_longer == FALSE ~ 0)) %>%
  pivot_longer(3:4, names_to = "test_condition", values_to = "avglook_time")

ggplot(testLooking_plot, aes(x = test_condition, y = avglook_time, fill = hab_condition)) +
  geom_boxplot(width = 0.5, outlier.colour = NA, alpha = 0.6, color="black", fill="blue") +
  #geom_quasirandom(width = .15, size = 1, alpha = 1) + geom_line(aes(group = Participant_ID)) + 
  geom_point(size = 1, alpha = 1) + geom_line(aes(group = interaction(Participant_ID)), alpha=.6) +
  geom_errorbar(data=all_test_looking, aes(ymin = se_lo, ymax = se_hi, y=NULL), width=0) +
  geom_point(data=all_test_looking, aes(y=looking), fill='white', shape=23, size=3) + 
  guides(fill = F) + 
  xlab('Test Trial Type') +
  ylab('Looking Time (seconds)') + ggtitle('Looking Time by Habituation Condition') +
  facet_wrap(~hab_condition, labeller = labeller(hab_condition = c("LR" ="Left-to-right" , "random" = "Nonlinear"))) + scale_x_discrete(breaks = c('congruent', 'incongruent'), labels = c('Congruent', 'Incongruent')) +
  theme_Publication()

#2. Contrasts and follow up tests using emmeans 
# Start with the interaction:
modelS.emm <- emmeans(modelS, ~ hab_condition | test_condition)

# Look at the contrast
contrast(modelSeattle.emm, 'tukey') %>%
  broom::tidy() %>%
  head()

# It might be nice to extract the estimates and plot them:
modelS.emm.df <-
  modelS.emm %>%
  broom::tidy()

modelS.emm.df %>%
  ggplot(aes(hab_condition, estimate, ymin=estimate - std.error, ymax=estimate + std.error, color = test_condition)) + 
  scale_color_manual(labels = c("Congruent", "Incongruent"), 
                     values = c("congruent" = "darkgoldenrod3", "incongruent" = "blue3"), 
                     name = "Test Condition") + 
  labs(x = "Habituation Condition",
    y = "Looking Time (seconds)") + 
  scale_x_discrete(labels=c("Left-to-right", "Nonlinear")) +
  geom_pointrange() + theme_Publication()

# 3. Are there more babies looking at incongruent vs congruent trials based on habituation condition?
# Average across the pairs because there is no effect of pair
testLooking_model <- testLooking_plot %>%
  pivot_wider(names_from = test_condition, values_from = avglook_time)

testLooking_model$hab_condition <- as.factor(testLooking_model$hab_condition)

cond_glmer <- glmer(inc_longer ~ hab_condition + (1 | Participant_ID), family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"), data = testLooking_model)
summary(cond_glmer)
# No effect

# Because there is no effect, let's histogram differences in looking times so we can visualize it
testLooking_model <- testLooking_model %>%
  mutate(look_difference = incongruent - congruent)
  
# Plot it
ggplot(testLooking_model, aes(x = look_difference, fill = hab_condition)) +
  geom_histogram(binwidth = 0.6, position = "dodge", color = "black", alpha = 0.7) +
  labs(
    title = "Difference in Looking at Test Trials by Condition",
    x = "Difference (seconds)",
    y = "Number of participants") +
  theme_Publication() +
  scale_fill_manual(labels = c("Left-to-right", "Nonlinear"), values = c("LR" = "purple", "random" = "yellow"), name = "Habituation Condition")
 

##Demographics
#avg age
mean(seattle_participants$Age)
#average and median of max habituation trials, across conditions.
mean(hab_participant$max_hab)
median(hab_participant$max_hab)

#how many participants in each condition
cbl_count <- seattle_participants %>%
  group_by(Parent_Counterbalance) %>%
  count(Counterbalance)

cbl_noside <- seattle_participants %>%
  group_by(Counterbalance) %>%
  count(Counterbalance)

#sex distrihution
seattle_participants %>%
  count(Sex)

countAgain <- testPair_included %>%
  group_by(Participant_ID) %>%
  slice(1) %>%
  group_by(version) %>%
  count(version)
 

# For additional quality check:
#Draw the AOI and plot each baby for each trial, save the plots, look at them.
raw_data_select <- read.csv("")
raw_data_select <- raw_data_select %>% filter(phase =="test")

#there are a lot of data points in sample data so I will bin these in 50 ms timebins"
bin_timestamps <- function(timestamps)#{
  return(floor(timestamps / 50) * 50)  #round down to nearest 50 ms
}

raw_data_binned <- raw_data_select %>%
  mutate(binned_timestamp = bin_timestamps(timestamp))

#get the first row per bin, per trial per participant
raw_data_binned <- raw_data_binned %>%
  group_by(Participant_ID, TRIAL_NUMBER ,binned_timestamp) %>%
  slice(1) #so now we should have gaze coordinates for every 50 ms instead of every 2 ms

#Screen dimensions
screen_width <- 1024
screen_height <- 768

#AOI coordinates:
aoi_left <- 315
aoi_top <- 768 - 410 #to adjust for inverted y-axis
aoi_right <- 685
aoi_bottom <- 768 - 670 #to adjust for inverted y-axis
output_folder <- "gazePlots"

# Plot for each participant and trial
raw_data_binned %>%
  group_by(Participant_ID, TRIAL_NUMBER) %>%
  do({ 
    p <- ggplot(., aes(x = GAZE_X, y = GAZE_Y)) +
      geom_rect(aes(xmin = aoi_left, xmax = aoi_right, ymin = aoi_bottom, ymax = aoi_top),
                fill = "orange", alpha = 0.3) +
      #plot gaze points, this is from each sample, not from fix report. 
      geom_point(color = "black", size = 1, alpha = 0.6) +
      coord_fixed(xlim = c(0, screen_width), ylim = c(0, screen_height)) +
      scale_y_reverse(limits = c(768, 0)) +
      labs(title = paste("Participant:", unique(.$Participant_ID), 
                         "Trial:", unique(.$TRIAL_NUMBER)),
           x = "GAZE_X", y = "GAZE_Y") +
      theme_Publication()
    
    #    #Save each plot as a file
    filename <- file.path(output_folder, paste0("Participant_", unique(.$Participant_ID), 
                                                "_Trial_", unique(.$TRIAL_NUMBER), ".png"))
    ggsave(filename, plot = p, width = 7, height = 7, dpi = 300)
    NULL
  })
 

