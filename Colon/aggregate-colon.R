# Using the Yale High Performance Clusters (HPC), runs would not always finish
# if I set the number of trials too high and so I would run smaller batches of
# trials at a time. This script (with slight modifications for the various runs)
# was used to aggregate the colon runs into one RData file

# Merge three different runs of d = 1000, p = 2000 for colon data
load("successes_1000_2000_1.RData")
combined <- successes
load("successes_1000_2000_2.RData")
combined$success <- rbind(combined$success, successes$success)
combined$repaired <- rbind(combined$repaired, successes$repaired)
load("successes_1000_2000_3.RData")
combined$success <- rbind(combined$success, successes$success)
combined$repaired <- rbind(combined$repaired, successes$repaired)
successes <- combined
save(successes, file = "successes_1000_2000.RData")
